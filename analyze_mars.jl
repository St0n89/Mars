#!/usr/bin/env julia
# ===========================================================================
# ANALYSIS: Mars FDTD — plane P wave through TAYAK crust  (pure stdlib)
#
# Compares arrivals in the synthetic traces against 1-D ray predictions
# for vertical incidence: direct P, Ps conversions at 77.4 / 10 / 1 km,
# PPmP Moho multiple. Primary metric: relative delays measured by cross-
# correlating each phase window against the direct-P envelope (all phases
# share the source wavelet, so the estimator is bias-free).
#
# No packages, no zlib: the PNG writer is hand-rolled (CRC32 + Adler32 +
# stored DEFLATE blocks). For a compressed DEFLATE encoder see bmp2png.jl.
#
# Usage: julia analyze_mars.jl        (data + output in ./mars_run)
# ===========================================================================

using Printf

# ---- TAYAK layers, grid-quantised (identical to config/env_mars.jl) ------
const LAYERS = [  # (d_top_km, d_bot_km, vp, vs)
    (0.0,   1.0,   3677.71, 1739.80),
    (1.0,   10.0,  4952.25, 2780.97),
    (10.0,  77.4,  5846.66, 3281.16),
    (77.4,  80.0,  7400.90, 4246.33),
    (80.0,  100.0, 7413.64, 4247.35),
    (100.0, 110.0, 7436.12, 4180.19),
    (110.0, 120.0, 7446.65, 4181.17),
]
const Z_REC = 0.1      # receiver depth in km (iz = 2)
const Z_SRC = 100.0     # source depth in km

function travel(z0, z1, mode)
    """1-D travel time between depths z0 (deeper) and z1 (shallower)."""
    t = 0.0
    for (top, bot, vp, vs) in LAYERS
        a = max(min(z1, bot), top)      # shallow end, clipped to layer
        b = min(max(z0, top), bot)      # deep end, clipped to layer
        b > a && (t += (b - a) * 1000.0 / (mode == 'P' ? vp : vs))
    end
    t
end

# ---- predictions ----------------------------------------------------------
t_P = travel(Z_SRC, Z_REC, 'P')
ps_delay(d)  = travel(d, Z_REC, 'S') - travel(d, Z_REC, 'P')  # S above, P below
ppm_delay(d) = 2.0 * travel(d, Z_REC, 'P')                    # surface + interface bounce

const PHASES = [
    ("P",      t_P,                    "direct P from 100 km"),
    ("Ps10",   t_P + ps_delay(10.0),   "Ps conversion at 10 km"),
    ("Ps77",   t_P + ps_delay(77.4),   "Ps conversion at Moho (77.4 km)"),
    ("PPmP10", t_P + ppm_delay(10.0),  "multiple off 10 km interface"),
    ("PPmP77", t_P + ppm_delay(77.4),  "PPmP Moho multiple"),
]

# ---- read traces -----------------------------------------------------------
base = joinpath(@__DIR__, "mars_run")
data = [split(ln, ',') for ln in readlines(joinpath(base, "traces.csv")) if !isempty(ln)][2:end]
t = [parse(Float64, r[1]) for r in data]
nt = length(t)
nrec = length(data[1]) - 1
traces = [[parse(Float64, r[j + 1]) for r in data] for j in 1:nrec]

meta = Dict{String,Float64}()
for ln in readlines(joinpath(base, "meta.jl"))
    p = findfirst('=', ln)
    p === nothing && continue
    k = strip(replace(strip(ln[1:p-1]), "const " => ""))
    meta[String(k)] = parse(Float64, strip(ln[p+1:end]))
end
dt = meta["dt"]
f0 = meta["f0"]

# ---- envelope --------------------------------------------------------------
PULSE_SPAN = 4.0 + 1.5     # source wavelet length + ringdown (s)

function envelope(tr)
    """Moving-average of |trace| over ~0.5 s."""
    w = max(1, round(Int, 0.5 / dt))
    n = length(tr)
    out = Vector{Float64}(undef, n)
    s = 0.0
    for i in 1:n
        s += abs(tr[i])
        i > w && (s -= abs(tr[i - w]))
        out[i] = s / min(i, w)
    end
    out
end

# ---- cross-correlation picking ---------------------------------------------
# All phases share the source wavelet, so aligning each phase window against
# the direct-P envelope estimates the onset delay without the envelope-peak
# / wavelet-shape bias.
idxf(tv) = clamp(round(Int, tv / dt), 0, nt - 1)   # 0-based index, like the array layout

function xcorr_pick(env, t_ref0, t_search0, t_search1)
    """Lag maximising normalised correlation of env with the P reference."""
    k0 = idxf(t_ref0); k1 = idxf(t_ref0 + 5.5)     # P reference wavelet
    K = k1 - k0
    K <= 0 && return nothing, 0.0
    r = env[k0+1:k1]
    rr = 0.0
    for v in r; rr += v * v; end
    rr = sqrt(rr)
    rr == 0.0 && return nothing, 0.0
    s0 = idxf(t_search0); s1 = idxf(t_search1) - K
    best = -1; best_sc = -1.0
    for L in s0:max(s0, s1)
        (L < 0 || L + K > nt) && continue
        seg = env[L+1:L+K]
        ss = 0.0
        for v in seg; ss += v * v; end
        ss = sqrt(ss)
        ss == 0.0 && continue
        sc = 0.0
        for i in 1:K
            sc += r[i] * seg[i]
        end
        sc /= (rr * ss)
        sc > best_sc && (best_sc = sc; best = L)
    end
    best < 0 && return nothing, 0.0
    pk = maximum(env[best+1:best+K])
    # onset = search start + lag, corrected for the reference offset
    t[best+1] - (t_ref0 - t_P), pk
end

@printf("receivers: %d, nt: %d, dt = %.5f s, f0 = %.1f Hz\n\n", nrec, nt, dt, f0)
envs = envelope.(traces)

fmt2(v) = v === nothing ? "  -  " : @sprintf("%.2f", v)
@printf("%-8s %7s %8s %7s %10s %10s %8s  %s\n",
        "phase", "pred s", "onset s", "err s", "delay pred", "delay meas",
        "peakenv", "note")
picks_rows = ["phase,pred_onset_s,meas_onset_s,onset_err_s,pred_delay_s,meas_delay_s,peak_env,note"]
for (name, tp, note) in PHASES
    ons = Float64[]; pks = Float64[]
    for env in envs
        to, pk = xcorr_pick(env, t_P - 1.0, tp - 2.5, tp + 8.0)
        to !== nothing && (push!(ons, to); push!(pks, pk))
    end
    on = isempty(ons) ? nothing : begin
        m = 0.0
        for v in ons; m += v; end
        m / length(ons)
    end
    pkv = isempty(pks) ? 0.0 : maximum(pks)
    pred_d = tp - t_P
    meas_d = on === nothing ? nothing : on - t_P
    if name != "P" && tp < t_P + PULSE_SPAN
        flag = "window contaminated by direct-P coda (unresolved)"
    elseif meas_d !== nothing
        flag = note * (name == "P" ? "" : @sprintf(", err %+.2f s", meas_d - pred_d))
    else
        flag = note
    end
    @printf("%-8s %7.2f %8s %7s %10.2f %10s %8.1e  %s\n",
            name, tp, fmt2(on), fmt2(on === nothing ? nothing : on - tp),
            pred_d, fmt2(meas_d), pkv, flag)
    cells = [name,
             @sprintf("%.3f", tp),
             on === nothing ? "" : @sprintf("%.3f", on),
             on === nothing ? "" : @sprintf("%.3f", on - tp),
             @sprintf("%.3f", pred_d),
             meas_d === nothing ? "" : @sprintf("%.3f", meas_d),
             @sprintf("%.3e", pkv),
             flag]
    push!(picks_rows, join(((',' in c) ? "\"" * c * "\"" : c for c in cells), ","))
end
open(joinpath(base, "picks.csv"), "w") do io
    println(io, join(picks_rows, "\n"))
end

# ===========================================================================
# RECORD-SECTION PNG — hand-rolled writer (CRC32, Adler32, stored DEFLATE)
# ===========================================================================
SPP = 16           # px per second
RH  = 34           # px per trace row
LM, RM = 6, 0      # left/right margin
W = LM + Int(floor(t[end] * SPP)) + RM
H = nrec * RH + 40

const BG  = (0x18, 0x18, 0x1e)
const POS = (0x3c, 0x5a, 0xc8)   # positive lobes: blue
const NEG = (0xc8, 0x5a, 0x3c)   # negative lobes: red
const GRN = (0x00, 0xaa, 0x3c)   # predicted arrivals
const AXS = (0x3c, 0x3c, 0x46)
const TXT = (0x96, 0x96, 0xa0)
const WIG = (0xeb, 0xeb, 0xeb)

img = fill(BG, W, H)
put0(x, y, c) = (0 <= x < W && 0 <= y < H) && (img[x+1, y+1] = c)  # 0-based

gmax = maximum(maximum(abs.(tr)) for tr in traces)
clip = gmax / 3.0

# stacked, individually-normalised wiggles + variable-density fill
for (j, tr) in enumerate(traces)
    y0 = (j - 1) * RH + RH ÷ 2            # 0-based row centre
    amp = maximum(abs.(tr))
    amp == 0.0 && continue
    half = RH ÷ 2 - 1
    prev_y = y0
    for (i, tv) in enumerate(t)
        x = LM + floor(Int, tv * SPP)   # 0-based column
        if abs(tr[i]) > amp * 0.15
            c = tr[i] > 0 ? POS : NEG
            sgn = tr[i] > 0 ? 1 : -1
            h = min(half, floor(Int, abs(tr[i]) / amp * half))
            for yy in min(y0, y0 + sgn * h):max(y0, y0 + sgn * h)
                put0(x, yy, c)
            end
        end
        y = round(Int, y0 - tr[i] / clip * half)
        y = max(y0 - half, min(y0 + half, y))
        for yy in min(prev_y, y):max(prev_y, y)
            put0(x, yy, abs(tr[i]) > amp * 0.02 ? WIG : BG)
        end
        prev_y = y
    end
end

# predicted arrival curves (green, dashed)
for (_, tp, _) in PHASES
    x = LM + floor(Int, tp * SPP)
    for y in 0:(nrec * RH - 1)
        (y ÷ 3) % 2 == 0 && put0(x, y, GRN)
    end
end

# receiver labels + time axis (tiny hand-drawn bitmap font)
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
    'k' => [(1,1),(2,1),(1,2),(2,2),(1,3),(2,3),(2,4),(3,5),(3,6)],
    'm' => [(1,1),(2,1),(1,2),(1,3),(2,3),(1,4),(1,5),(3,2),(3,5),(1,6),(3,6)],
    's' => [(2,1),(3,2),(1,2),(2,3),(3,4),(3,5),(1,5),(2,6),(1,6),(3,6)],
)
function draw_text0(s, x0, y0, col)
    for (ci, ch) in enumerate(s)
        for (dx, dy) in get(DIG, ch, Tuple{Int,Int}[])
            put0(x0 + (ci - 1) * 4 + dx - 1, y0 + dy - 1, col)
        end
    end
end

for j in 0:nrec-1
    draw_text0(string(14 + j) * "km", 0, j * RH + RH ÷ 2 - 4, TXT)
end
for sec in 0:10:Int(floor(t[end]))
    draw_text0(string(sec) * "s", LM + sec * SPP - 4, H - 14, TXT)
    for y in 0:(nrec * RH - 1)
        put0(LM + sec * SPP, y, AXS)
    end
end

# ---- PNG assembly ----------------------------------------------------------
const CRC_TAB = let tab = Vector{UInt32}(undef, 256)
    for n in 0:255
        c = UInt32(n)
        for _ in 1:8
            c = (c & 0x00000001) == 0 ? c >> 1 : (c >> 1) ⊻ 0xedb88320
        end
        tab[n+1] = c
    end
    tab
end
function crc32(data)
    c = 0xffffffff
    for b in data
        c = CRC_TAB[((c ⊻ b) & 0xff) + 1] ⊻ (c >> 8)
    end
    c ⊻ 0xffffffff
end

function adler32(data)
    a = UInt32(1); b = UInt32(0)
    for byte in data
        a = (a + byte) % 65521
        b = (b + a) % 65521
    end
    (b << 16) | a
end

"""zlib stream with STORED deflate blocks (BTYPE=00): valid, uncompressed."""
function zlib_stored(raw)
    out = UInt8[0x78, 0x01]
    i = 1
    n = length(raw)
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

be32(x) = UInt8[(x >> 24) & 0xff, (x >> 16) & 0xff, (x >> 8) & 0xff, x & 0xff]
function chunk(tag, data)
    b = collect(codeunits(tag))
    vcat(be32(length(data)), b, data, be32(crc32(vcat(b, data))))
end

function pack_rows(img, W, H)
    raw = Vector{UInt8}(undef, H * (1 + 3W))
    o = 0
    for y in 1:H
        raw[o += 1] = 0x00            # filter type: None
        for x in 1:W
            c = img[x, y]
            raw[o += 1] = c[1]
            raw[o += 1] = c[2]
            raw[o += 1] = c[3]
        end
    end
    raw
end
raw = pack_rows(img, W, H)

ihdr = vcat(be32(W), be32(H), UInt8[8, 2, 0, 0, 0])   # 8-bit RGB
png = vcat(UInt8[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a],
           chunk("IHDR", ihdr),
           chunk("IDAT", zlib_stored(raw)),
           chunk("IEND", UInt8[]))
open(joinpath(base, "section.png"), "w") do io
    write(io, png)
end
@printf("section.png: %dx%dpx (%.1f MB, stored deflate), picks.csv written\n",
        W, H, length(png) / 1e6)
