using Printf
# ===========================================================================
# RAY FAN: tomography feasibility under InSight. Traces teleseismic P rays
# through TAYAK for every LOCATED event in results/tomo/event_geometry.csv
# (MQS catalog v14, events with a P pick; unlocated events sit at distance
# 0 and are excluded).
#
# Method: exact earth-flattening transform, z = R ln(R/r), v_f = (R/r) v,
# which preserves travel time, ray parameter and surface incidence.
# TAYAK becomes piecewise-linear layers with velocity jumps (1, 10, 77.368
# km, CMB); per-layer ray integrals are closed-form, so X(p) and T(p) are
# exact. TAYAK's Moho jump (5.85 -> 7.40 km/s) makes X(p) multi-branched:
# crustal rays, the Moho head wave (Pn) and mantle-turning P all compete.
# For each event ALL branches are evaluated (ray roots of X(p) = R*delta
# plus head waves at the jumps) and the first arrival (min T) is kept.
# Observed P travel times from the catalog (observed_P_times.csv) are
# printed for comparison.
#
# Writes results/tomo/ray_fan.csv, results/tomo/ray_map.png (polar map)
# and results/tomo/ray_fan.png (E-W / N-S cross-sections of the fan).
#
# Run from the repo root:  julia tools/ray_fan.jl
# Stdlib only. TAYAK.mod: data/ (Khan et al. 2016, via the IPGP
# dataverse, doi:10.18715/IPGP.2021.kqwc4yaf); catalog: MQS v14.
# ===========================================================================

const R_MARS = 3389.5
const MOHO_TAYAK = 77.368
const TAYAK = get(ENV, "MARS_TAYAK_MOD",
                  normpath(joinpath(@__DIR__, "..", "data", "TAYAK.mod")))
const GEOM = normpath(joinpath(@__DIR__, "..", "results", "tomo", "event_geometry.csv"))
const OBS = normpath(joinpath(@__DIR__, "..", "results", "tomo", "observed_P_times.csv"))
const OUTDIR = normpath(joinpath(@__DIR__, "..", "results", "tomo"))

# ---- TAYAK -> flattened layer list ---------------------------------------------
# raw rows (depth, vp), repeated depth = velocity jump
const ROWS = Tuple{Float64,Float64}[]
let seen_header = false
    for line in eachline(TAYAK)
        seen_header || (occursin("depth", line) && (seen_header = true); continue)
        f = split(line)
        length(f) < 2 && continue
        d = parse(Float64, f[1]); vp = parse(Float64, f[2])
        d > 1597.0 && continue          # mantle only: P never crosses the CMB
        push!(ROWS, (d, vp))
    end
end
@assert ROWS[1][1] == 0.0
zf(d) = R_MARS * log(R_MARS / (R_MARS - d))          # flattened depth
vf(d, v) = R_MARS / (R_MARS - d) * v                 # flattened velocity
realdepth(z) = R_MARS * (1.0 - exp(-z / R_MARS))

# ITEMS in depth order: ("seg", z0, v0, z1, v1) or ("jmp", z, v_below)
const ITEMS = Tuple{String,Float64,Float64,Float64,Float64}[]
let d0 = ROWS[1][1], v0 = ROWS[1][2]
    for (d, v) in ROWS[2:end]
        if d == d0                                       # jump
            v > v0 && push!(ITEMS, ("jmp", zf(d), vf(d, v), 0.0, 0.0))
        else
            push!(ITEMS, ("seg", zf(d0), vf(d0, v0), zf(d), vf(d, v)))
        end
        d0, v0 = d, v
    end
end
const V_SURF = vf(0.0, ROWS[1][2])                      # 3.67771 km/s
asech(u) = log((1.0 + sqrt(max(0.0, 1.0 - u * u))) / u)

# ---- exact flat-layer ray (descending from the surface) -------------------------
# returns (X, T, z_turn, kind, v_hw): X/T are surface-to-surface (2x the
# receiver-side half), kind in ("turn", "hw", "core", "stop").
function ray(p, z_stop = Inf)
    x = 0.0; t = 0.0
    for (kind0, z0, v0, z1, v1) in ITEMS
        if kind0 == "jmp"
            if z0 > z_stop
                return (x, t, z0, "stop", 0.0)
            elseif p * v0 >= 1.0
                return (2x, 2t, z0, "hw", v0)
            end
        else
            if z1 > z_stop
                z1 = z_stop
            end
            u0 = p * v0; u1 = p * v1
            g = (v1 - v0) / (z1 - z0)
            if u0 >= 1.0
                return (2x, 2t, z0, "turn", 0.0)
            elseif g > 0.0 && u1 >= 1.0                  # turning inside layer
                zt = z0 + (1.0 / p - v0) / g
                x += sqrt(1.0 - u0 * u0) / (p * g)
                t += asech(u0) / g
                return (2x, 2t, zt, "turn", 0.0)
            elseif abs(g) < 1e-12
                x += p * v0 / sqrt(1.0 - u0 * u0) * (z1 - z0)
                t += (z1 - z0) / (v0 * sqrt(1.0 - u0 * u0))
            else
                x += (sqrt(1.0 - u0 * u0) - sqrt(1.0 - u1 * u1)) / (p * g)
                t += (asech(u0) - asech(u1)) / g
            end
            if z1 >= z_stop
                return (x, t, z1, "stop", 0.0)
            end
        end
    end
    (2x, 2t, Inf, "core", 0.0)                          # reached the CMB
end
X_of_p(p) = ray(p)[1]

function bisect_X(target, lo, hi)
    for _ in 1:80
        m = 0.5 * (lo + hi)
        X_of_p(m) < target ? (lo = m) : (hi = m)
    end
    0.5 * (lo + hi)
end

# all candidate first arrivals for a surface distance X_t (km)
function first_arrival(X_t)
    cands = NamedTuple[]                                # (T, p, turn, kind, vhw)
    # ray roots: scan p, bisect sign changes; verify against discontinuities
    plo, phi = 1e-5, 0.185
    npg = 3000
    xprev = nothing
    for i in 0:npg
        p = plo + (phi - plo) * i / npg
        X, T, zt, kind, vhw = ray(p)
        if kind == "core"
            xprev = nothing
            continue
        end
        if xprev !== nothing
            if (xprev[1] - X_t) * (X - X_t) < 0
                pr = bisect_X(X_t, xprev[2], p)
                Xr, Tr, ztr, kindr, _ = ray(pr)
                if kindr != "core" && abs(Xr - X_t) < 0.5
                    push!(cands, (T = Tr, p = pr, turn = ztr, kind = kindr, vhw = 0.0))
                end
            end
        end
        xprev = (X, p)
    end
    # head waves at velocity jumps: T = T_crit + (X - X_crit)/v_below
    for (kind0, z, vbelow, _, _) in ITEMS
        kind0 == "jmp" || continue
        pc = (1.0 + 1e-12) / vbelow
        Xc, Tc, _, kindc, _ = ray(pc)
        kindc == "hw" || continue
        X_t >= Xc - 0.5 && push!(cands, (T = Tc + (X_t - Xc) / vbelow, p = pc,
                                          turn = z, kind = "hw", vhw = vbelow))
    end
    isempty(cands) && return nothing
    cands[argmin([c.T for c in cands])]
end
# receiver-side horizontal offset at flattened depth z_stop (km)
function xdown(p, z_stop)
    X, _, _, kind, _ = ray(p, z_stop)
    return kind == "stop" ? X : 0.5 * X
end

# ---- events -----------------------------------------------------------------------
events = []
let lines = readlines(GEOM)
    for line in lines[2:end]
        f = split(line, ',')
        dist = parse(Float64, f[9]); baz = parse(Float64, f[10])
        0.0 < dist <= 90.0 || continue     # beyond 90 deg: PKP/PP territory
        push!(events, (name = f[1], cls = f[2], dist = dist, baz = baz))
    end
    global NUNLOCATED = count(l -> parse(Float64, split(l, ',')[9]) == 0.0, lines[2:end])
end
obs = Dict{String,Float64}()
isfile(OBS) && for line in readlines(OBS)[2:end]
    f = split(line, ','); obs[f[1]] = parse(Float64, f[3])
end

# ---- trace --------------------------------------------------------------------------
mkpath(OUTDIR)
rows = NamedTuple[]
for ev in events
    X_t = R_MARS * deg2rad(ev.dist)
    ray0 = first_arrival(X_t)
    if ray0 === nothing
        @printf("%-8s: %5.1f deg: no P branch reached (core branch), skipped\n", ev.name, ev.dist)
        continue
    end
    inc = asind(clamp(ray0.p * V_SURF, 0.0, 1.0))
    xmo = xdown(ray0.p, zf(MOHO_TAYAK))
    Tob = get(obs, ev.name, NaN)
    push!(rows, (name = ev.name, cls = ev.cls, dist = ev.dist, baz = ev.baz,
                 p = ray0.p, inc = inc, turn = ray0.turn, T = ray0.T,
                 kind = ray0.kind, xmoho = xmo, Tobs = Tob))
    @printf("%-8s %-13s %5.1f deg  baz %5.1f  p %6.4f s/km  inc %4.1f  %s  T %5.0f s (obs %s)\n",
            ev.name, ev.cls, ev.dist, ev.baz, ray0.p, inc,
            ray0.kind == "hw" ? "Pn head wave" : "turn $(round(Int, realdepth(ray0.turn))) km",
            ray0.T, isnan(Tob) ? "  -  " : @sprintf("%4.0f s", Tob))
end

open(joinpath(OUTDIR, "ray_fan.csv"), "w") do fh
    println(fh, "event,class,dist_deg,baz_deg,p_s_per_km,inc_deg,kind,turn_depth_km,T_pred_s,T_obs_s,x_at_moho_km")
    for r in rows
        @printf(fh, "%s,%s,%.3f,%.2f,%.4f,%.2f,%s,%.1f,%.1f,%s,%.2f\n",
                r.name, r.cls, r.dist, r.baz, r.p, r.inc, r.kind,
                realdepth(r.turn), r.T, isnan(r.Tobs) ? "" : string(round(Int, r.Tobs)), r.xmoho)
    end
end

# fan paths for the figure: horizontal offset vs real depth, 0..300 km
function fanpath(p, hmax = 300.0)
    hs = collect(0.0:2.0:hmax)
    (hs, [h <= 0.0 ? 0.0 : xdown(p, zf(h)) for h in hs])
end
PATHS = [fanpath(r.p) for r in rows]

# ---- PNG writer -------------------------------------------------------------------
const CRC_TAB = let tab = Vector{UInt32}(undef, 256)
    for n in 0:255
        cc = UInt32(n)
        for _ in 1:8
            cc = (cc & 0x00000001) == 0 ? cc >> 1 : (cc >> 1) ⊻ 0xedb88320
        end
        tab[n+1] = cc
    end
    tab
end
crc32(data) = let c = 0xffffffff
    for b in data; c = CRC_TAB[((c ⊻ b) & 0xff) + 1] ⊻ (c >> 8); end
    c ⊻ 0xffffffff
end
adler32(data) = let a = 1, b = 0
    for bb in data; a = (a + bb) % 65521; b = (b + a) % 65521; end
    b * 65536 + a
end
be32(n) = UInt8[(n >> 24) & 255, (n >> 16) & 255, (n >> 8) & 255, n & 255]
chunk(tag, data) = vcat(be32(length(data)), Vector{UInt8}(tag), data,
                        be32(crc32(vcat(Vector{UInt8}(tag), data))))
function zlib_stored(raw)
    out = UInt8[0x78, 0x01]
    i = 1; n = length(raw)
    while true
        m = min(65535, n - i + 1)
        final = (i + m - 1 >= n)
        push!(out, final ? 0x01 : 0x00)
        push!(out, UInt8(m & 0xff), UInt8((m >> 8) & 0xff))
        push!(out, UInt8((m ⊻ 0xffff) & 0xff), UInt8(((m ⊻ 0xffff) >> 8) & 0xff))
        append!(out, raw[i:i+m-1])
        i += m
        final && break
    end
    a = adler32(raw)
    push!(out, UInt8((a >> 24) & 0xff), UInt8((a >> 16) & 0xff),
                 UInt8((a >> 8) & 0xff), UInt8(a & 0xff))
    out
end
function write_png(path, W, H, pix)
    raw = UInt8[]
    for y in 1:H
        push!(raw, 0x00)
        for x in 1:W; append!(raw, @view pix[y, x, :]); end
    end
    open(path, "w") do fh
        write(fh, UInt8[137, 80, 78, 71, 13, 10, 26, 10])
        write(fh, chunk("IHDR", vcat(be32(W), be32(H), UInt8[0x08, 0x02, 0, 0, 0])))
        write(fh, chunk("IDAT", zlib_stored(raw)))
        write(fh, chunk("IEND", UInt8[]))
    end
end

const BLACK = (0x00, 0x00, 0x00); const WHITE = (0xff, 0xff, 0xff)
const RED = (0xc8, 0x1e, 0x1e); const BLUE = (0x32, 0x82, 0xc8)
const GREY = (0x60, 0x60, 0x60)
const FS = 2
const DIG = Dict{Char,Vector{Tuple{Int,Int}}}(
    '0' => [(1,1),(2,1),(3,1),(1,2),(3,2),(1,3),(3,3),(1,4),(3,4),(2,5),(1,5),(3,5),(1,6),(3,6)],
    '1' => [(3,1),(2,1),(3,2),(3,3),(3,4),(3,5),(3,6),(2,4),(1,3)],
    '2' => [(1,1),(2,1),(3,2),(3,3),(2,4),(1,5),(1,6),(2,6),(3,6)],
    '3' => [(1,1),(2,1),(3,1),(3,2),(2,3),(3,4),(3,5),(1,5),(2,6),(1,6),(3,6)],
    '4' => [(1,1),(1,2),(1,3),(2,3),(3,1),(3,2),(3,3),(3,4),(3,5),(3,6)],
    '5' => [(1,1),(2,1),(3,1),(1,2),(1,3),(2,3),(3,4),(3,5),(1,5),(2,6),(3,6)],
    '6' => [(2,1),(3,2),(1,2),(1,3),(1,4),(2,4),(3,4),(3,5),(3,6),(1,5),(1,6),(2,6)],
    '7' => [(1,1),(2,1),(3,1),(3,2),(3,3),(2,4),(2,5),(1,6)],
    '8' => [(2,1),(1,1),(3,1),(1,2),(3,2),(2,3),(1,4),(3,4),(1,5),(3,5),(2,6),(1,6),(3,6)],
    '9' => [(2,1),(1,1),(3,1),(1,2),(3,2),(1,3),(2,3),(3,4),(3,5),(1,4),(3,6),(2,6)],
    '.' => [(2,5),(2,6)],
    '-' => [(1,4),(2,4),(3,4)],
    'E' => [(1,1),(2,1),(3,1),(1,2),(1,3),(2,3),(1,4),(1,5),(2,5),(1,6),(3,6)],
    'W' => [(1,1),(3,1),(1,2),(3,2),(1,3),(3,3),(2,3),(1,4),(3,4),(1,5),(3,5),(1,6),(3,6)],
    'N' => [(1,1),(1,2),(1,3),(1,4),(1,5),(1,6),(2,2),(3,3),(3,1),(3,2),(3,4),(3,5),(3,6)],
    'S' => [(2,1),(3,2),(1,2),(2,3),(3,4),(3,5),(1,5),(2,6),(1,6),(3,6)],
    'T' => [(1,1),(2,1),(3,1),(2,2),(2,3),(2,4),(2,5),(2,6)],
    'A' => [(2,1),(1,2),(3,2),(1,3),(2,3),(3,3),(1,4),(3,4),(1,5),(3,5),(1,6),(3,6)],
    'Y' => [(1,1),(3,1),(1,2),(3,2),(2,3),(2,4),(2,5),(2,6)],
    'K' => [(1,1),(1,2),(1,3),(1,4),(1,5),(1,6),(3,1),(2,2),(2,4),(3,5),(3,6)],
    'M' => [(1,1),(2,1),(3,2),(1,2),(1,3),(2,3),(1,4),(1,5),(3,5),(1,6),(3,6)],
    'I' => [(1,1),(2,1),(3,1),(2,2),(2,3),(2,4),(2,5),(1,6),(2,6),(3,6)],
    'o' => [(2,2),(1,3),(3,3),(1,4),(3,4),(1,5),(3,5),(2,6)],
    'h' => [(1,1),(1,2),(1,3),(1,4),(1,5),(1,6),(2,2),(3,3),(3,4),(3,5),(3,6)],
    'k' => [(1,1),(1,2),(1,3),(1,4),(1,5),(1,6),(3,3),(2,4),(3,5)],
    'm' => [(1,2),(2,2),(1,3),(1,4),(2,4),(1,5),(3,3),(3,5),(1,6),(3,6)],
    'p' => [(1,2),(1,3),(1,4),(1,5),(1,6),(2,2),(3,3),(2,4)],
    'e' => [(2,2),(1,3),(3,3),(1,4),(2,4),(3,4),(1,5),(2,6),(3,6)],
    'g' => [(2,2),(1,3),(3,3),(1,4),(2,4),(3,5),(1,5),(2,6),(3,6)],
    'd' => [(3,1),(3,2),(3,3),(3,4),(3,5),(3,6),(2,2),(1,3),(1,4),(1,5),(1,6)],
    't' => [(2,1),(2,2),(2,3),(2,4),(2,5),(2,6),(1,2),(3,2)],
    'c' => [(2,2),(3,3),(1,3),(1,4),(1,5),(2,6),(3,6)],
    'u' => [(1,2),(3,2),(1,3),(3,3),(1,4),(3,4),(1,5),(3,5),(1,6),(2,6),(3,6)],
    'n' => [(1,2),(1,3),(1,4),(1,5),(1,6),(2,2),(3,3),(3,4),(3,5),(3,6)],
    'v' => [(1,2),(3,2),(1,3),(3,3),(2,4),(2,5)],
    'r' => [(1,2),(1,3),(1,4),(1,5),(1,6),(2,2),(3,3)],
    'l' => [(2,1),(2,2),(2,3),(2,4),(2,5),(2,6)],
    'i' => [(2,1),(2,3),(2,4),(2,5),(2,6)],
    'a' => [(2,2),(1,3),(3,3),(2,4),(1,5),(3,5),(2,6),(1,6),(3,6)],
    'b' => [(1,1),(1,2),(1,3),(1,4),(1,5),(1,6),(2,2),(3,3),(3,4),(2,5)],
    'f' => [(2,1),(2,2),(2,3),(2,4),(2,5),(2,6),(1,1),(3,2)],
    's' => [(2,2),(1,3),(2,3),(3,4),(1,5),(2,6),(3,6),(3,5)],
    'q' => [(2,2),(1,3),(3,3),(1,4),(2,4),(3,4),(1,5),(3,5),(2,6),(3,6)],
    'P' => [(1,1),(2,1),(1,2),(3,2),(1,3),(2,3),(1,4),(1,5),(1,6)],
)
function draw_in(pix, s, x0, y0, col)
    for (ci, ch) in enumerate(s)
        for (dx, dy) in get(DIG, ch, Tuple{Int,Int}[])
            for sy in 0:FS-1, sx in 0:FS-1
                x = x0 + (ci - 1) * 4 * FS + (dx - 1) * FS + sx
                y = y0 + (dy - 1) * FS + sy
                (1 <= x <= size(pix, 2) && 1 <= y <= size(pix, 1)) &&
                    (pix[y, x, 1] = col[1]; pix[y, x, 2] = col[2]; pix[y, x, 3] = col[3])
            end
        end
    end
end
text_w(s) = length(s) * 4 * FS
function putc_in(pix, x, y, c)
    (1 <= x <= size(pix, 2) && 1 <= y <= size(pix, 1)) &&
        (pix[y, x, 1] = c[1]; pix[y, x, 2] = c[2]; pix[y, x, 3] = c[3])
end
hline_in(pix, y, x0, x1, c) = (for x in min(x0, x1):max(x0, x1); putc_in(pix, x, y, c); end)

# ---- figure 1: polar map ------------------------------------------------------------
let W = 560, H = 620
    pix = fill(UInt8(255), H, W, 3)
    draw(s, x, y, c) = draw_in(pix, s, x, y, c)
    putc(x, y, c) = putc_in(pix, x, y, c)
    cx, cy, RPX = 280, 340, 240
    for ring in (30, 60, 90)
        rr = RPX * ring / 90
        for a in 0:0.01:2pi
            putc(round(Int, cx + rr * sin(a)), round(Int, cy - rr * cos(a)), GREY)
        end
        draw(string(ring), cx + 3, round(Int, cy - rr) + 2, BLACK)
    end
    draw("N", cx - 4, cy - RPX - 24, BLACK)
    draw("E", cx + RPX + 8, cy - 8, BLACK)
    draw("MQS v14: located events with P picks", 20, 12, BLACK)
    draw("rings: epicentral distance, deg", 20, H - 40, BLACK)
    draw("red: broadband, blue: low frequency", 20, H - 22, BLACK)
    lab = ["S1222a", "S0173a", "S0235b", "S1102a", "S1094b", "S0325a"]
    for r in rows
        rr = RPX * min(r.dist, 90) / 90
        xx = round(Int, cx + rr * sind(r.baz)); yy = round(Int, cy - rr * cosd(r.baz))
        col = r.cls == "BROADBAND" ? RED : BLUE
        for dy in -2:2, dx in -2:2; putc(xx + dx, yy + dy, col); end
        r.name in lab && draw(r.name, xx + 5, yy - 6, BLACK)
    end
    write_png(joinpath(OUTDIR, "ray_map.png"), W, H, pix)
end

# ---- figure 2: cross-section fan ------------------------------------------------------
let W = 1000, H = 580
    pix = fill(UInt8(255), H, W, 3)
    draw(s, x, y, c) = draw_in(pix, s, x, y, c)
    putc(x, y, c) = putc_in(pix, x, y, c)
    panels = [("E-W section", 70, 470), ("N-S section", 530, 930)]
    HMAX = 300.0; PM = 300.0
    for (title, x0, x1) in panels
        top, bot = 60, H - 60
        pw, ph = x1 - x0, bot - top
        for y in top:bot; putc(x0, y, BLACK); putc(x1, y, BLACK); end
        for x in x0:x1; putc(x, top, BLACK); putc(x, bot, BLACK); end
        ym = round(Int, top + MOHO_TAYAK / HMAX * ph)
        hline_in(pix, ym, x0, x1, GREY)
        draw("TAYAK Moho", x0 + 8, ym + 6, GREY)
        for v in 0:100:300
            y = round(Int, top + v / HMAX * ph)
            hline_in(pix, y, x0, x1, BLACK)
            draw(string(v), x0 - 40, y - 6, BLACK)
        end
        for v in -200:200:200
            x = round(Int, x0 + pw / 2 + v / PM * (pw / 2))
            for y in top:bot; putc(x, y, BLACK); end
            draw(string(v), x - 12, bot + 8, BLACK)
        end
        draw(title, x0 + 8, 16, BLACK)
        draw("km", x1 - text_w("km") - 4, bot + 26, BLACK)
        xs = round(Int, x0 + pw / 2)
        for d in 0:6; for xx in (xs - d):(xs + d); putc(xx, top + 1 + d, RED); end; end
        for (i, r) in enumerate(rows)
            t = clamp((r.dist - 27.0) / (73.0 - 27.0), 0.0, 1.0)
            col = (round(UInt8, 50 + t * 170), round(UInt8, 130 - t * 80), round(UInt8, 200 - t * 150))
            (hs, xs_) = PATHS[i]
            prevx, prevy = xs, top
            for (j, h) in enumerate(hs)
                off = xs_[j] * (title[1] == 'E' ? sind(r.baz) : cosd(r.baz))
                x = round(Int, x0 + pw / 2 + off / PM * (pw / 2))
                y = round(Int, top + h / HMAX * ph)
                n = max(abs(x - prevx), abs(y - prevy))
                for k in 0:n
                    putc(round(Int, prevx + (x - prevx) * k / max(n, 1)),
                         round(Int, prevy + (y - prevy) * k / max(n, 1)), col)
                end
                prevx, prevy = x, y
            end
        end
    end
    draw("P ray fan beneath InSight, colored 27-73 deg distance", 70, H - 16, BLACK)
    write_png(joinpath(OUTDIR, "ray_fan.png"), W, H, pix)
end

# ---- summary ------------------------------------------------------------------------
npn = count(r -> r.kind == "hw", rows)
nmantle = count(r -> r.kind == "turn", rows)
@printf("\n%d located events traced: %d Pn head wave, %d mantle-turning P; %d unlocated (dist 0) excluded\n",
        length(rows), npn, nmantle, NUNLOCATED)
incs = [r.inc for r in rows]
@printf("incidence at receiver: %.1f-%.1f deg (median %.1f)\n",
        minimum(incs), maximum(incs), incs[div(end + 1, 2)] + 0.0 * 0)
sorted_inc = sort(incs)
@printf("incidence sorted: min %.1f, median %.1f, max %.1f\n",
        sorted_inc[1], sorted_inc[div(length(sorted_inc), 2) + 1], sorted_inc[end])
xm = [r.xmoho for r in rows]
@printf("TAYAK Moho piercing offsets: %.1f-%.1f km\n", minimum(xm), maximum(xm))
bazs = [r.baz for r in rows if 44 <= r.baz <= 137]
@printf("main wedge: %d events, baz %.0f-%.0f deg\n", length(bazs), minimum(bazs), maximum(bazs))
res = [r.Tobs - r.T for r in rows if !isnan(r.Tobs)]
if !isempty(res)
    @printf("T_pred vs T_obs: mean resid %+.0f s (rms %.0f s) -> TAYAK %.0f%% slow\n",
            sum(res) / length(res), sqrt(sum(x -> x * x, res) / length(res)),
            100 * sum(res) / sum(r.Tobs for r in rows if !isnan(r.Tobs)))
end
