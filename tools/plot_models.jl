using Printf
# ===========================================================================
# PLOT MODELS: range/depth slice of the FDTD environment for a given
# model — where the inversion runs. Draws the 1-D layering (colored by
# vp), crust interfaces, the Moho, the 13-receiver surface array
# (x 53-65 km), the tapered source line at z=100 km (x 10-82 km), the
# 10-deg slant arrow and the schematic conversion ray from the array
# centre down to its Moho piercing point.
#
# Writes results/models/<model>_slice.png.
#
# Run from the repo root:
#   julia tools/plot_models.jl                 # the 4 inversion scenarios
#   julia tools/plot_models.jl TAYAK invA_r2   # any models
# Stdlib only; models come from config/models_gen.jl.
# ===========================================================================

include(joinpath(@__DIR__, "..", "config", "models.jl"))

const MODELS = isempty(ARGS) ? ["invA_r8", "invB_r8", "invB_r5", "invA_r1"] : ARGS
const SLANT_DEG = 10.0

# ---- PNG writer (CRC32 + Adler32 + stored DEFLATE, same as analyze_mars.jl)
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
    for b in data
        c = CRC_TAB[((c ⊻ b) & 0xff) + 1] ⊻ (c >> 8)
    end
    c ⊻ 0xffffffff
end
adler32(data) = let a = 1, b = 0
    for bb in data
        a = (a + bb) % 65521
        b = (b + a) % 65521
    end
    b * 65536 + a
end
be32(n) = UInt8[(n >> 24) & 255, (n >> 16) & 255, (n >> 8) & 255, n & 255]
chunk(tag, data) = vcat(be32(length(data)), Vector{UInt8}(tag), data,
                        be32(crc32(vcat(Vector{UInt8}(tag), data))))
function zlib_stored(raw)
    out = UInt8[0x78, 0x01]                       # zlib header (CM=8, FLEVEL=0)
    i = 1
    n = length(raw)
    while true
        m = min(65535, n - i + 1)
        final = (i + m - 1 >= n)
        push!(out, final ? 0x01 : 0x00)
        push!(out, UInt8(m & 0xff), UInt8((m >> 8) & 0xff))          # LEN little-endian
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
        for x in 1:W
            append!(raw, @view pix[y, x, :])
        end
    end
    open(path, "w") do fh
        write(fh, UInt8[137, 80, 78, 71, 13, 10, 26, 10])   # PNG signature
        write(fh, chunk("IHDR", vcat(be32(W), be32(H), UInt8[0x08, 0x02, 0, 0, 0])))
        write(fh, chunk("IDAT", zlib_stored(raw)))
        write(fh, chunk("IEND", UInt8[]))
    end
end

# ---- canvas -----------------------------------------------------------------
const W, H = 880, 560
const LM, RM, TM, BM = 56, 24, 40, 100
const PW, PH = W - LM - RM, H - TM - BM
const XMAX_KM, ZMAX_KM = 120.0, 120.0
px(x) = LM + round(Int, x / XMAX_KM * PW)
py(z) = TM + round(Int, z / ZMAX_KM * PH)

const BLACK = (0x00, 0x00, 0x00)
const WHITE = (0xff, 0xff, 0xff)
const RED   = (0xc8, 0x1e, 0x1e)
const GREY  = (0x60, 0x60, 0x60)
pix = fill(UInt8(255), H, W, 3)
putc(x, y, c) = (1 <= x <= W && 1 <= y <= H &&
                (pix[y, x, 1] = c[1]; pix[y, x, 2] = c[2]; pix[y, x, 3] = c[3]))
hlinec(y, x0, x1, c) = (for x in min(x0, x1):max(x0, x1); putc(x, y, c); end)
vlinec(x, y0, y1, c) = (for y in min(y0, y1):max(y0, y1); putc(x, y, c); end)
function linec(x0, y0, x1, y1, c, thick = 1)
    n = max(abs(x1 - x0), abs(y1 - y0))
    for i in 0:n
        x = round(Int, x0 + (x1 - x0) * i / max(n, 1))
        y = round(Int, y0 + (y1 - y0) * i / max(n, 1))
        for t in 0:thick-1
            putc(x + t, y, c)
        end
    end
end

# vp colormap: 5 stops, 2000..9000 m/s
const CSTOPS = [(0.00, (60, 20, 120)), (0.25, (50, 130, 200)),
                (0.50, (40, 190, 140)), (0.75, (250, 210, 70)),
                (1.00, (220, 50, 50))]
function vpcolor(vp)
    t = clamp((vp - 2000.0) / 7000.0, 0.0, 1.0)
    k = 1
    while k < length(CSTOPS) && t > CSTOPS[k + 1][1]; k += 1; end
    (k == length(CSTOPS)) && return map(UInt8, CSTOPS[end][2])
    t0, c0 = CSTOPS[k]
    t1, c1 = CSTOPS[k + 1]
    f = (t - t0) / (t1 - t0)
    (round(UInt8, c0[1] + f * (c1[1] - c0[1])),
     round(UInt8, c0[2] + f * (c1[2] - c0[2])),
     round(UInt8, c0[3] + f * (c1[3] - c0[3])))
end

# ---- 3x6 pixel font ----------------------------------------------------------
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
    '/' => [(3,1),(3,2),(2,3),(1,4),(1,5)],
    '_' => [(1,6),(2,6),(3,6)],
    'A' => [(2,1),(1,2),(3,2),(1,3),(2,3),(3,3),(1,4),(3,4),(1,5),(3,5),(1,6),(3,6)],
    'B' => [(1,1),(2,1),(1,2),(3,2),(1,3),(2,3),(1,4),(3,4),(1,5),(3,5),(1,6),(2,6)],
    'C' => [(2,1),(3,2),(1,2),(1,3),(1,4),(1,5),(3,6),(2,6)],
    'D' => [(1,1),(2,1),(1,2),(3,2),(1,3),(3,3),(1,4),(3,4),(1,5),(3,5),(1,6),(2,6)],
    'E' => [(1,1),(2,1),(3,1),(1,2),(1,3),(2,3),(1,4),(1,5),(2,5),(1,6),(3,6)],
    'H' => [(1,1),(3,1),(1,2),(3,2),(1,3),(2,3),(3,3),(1,4),(3,4),(1,5),(3,5),(1,6),(3,6)],
    'I' => [(2,1),(2,2),(2,3),(2,4),(2,5),(2,6)],
    'K' => [(1,1),(1,2),(1,3),(1,4),(1,5),(1,6),(3,1),(2,2),(2,4),(3,5),(3,6)],
    'L' => [(1,1),(1,2),(1,3),(1,4),(1,5),(1,6),(2,6),(3,6)],
    'M' => [(1,1),(2,1),(3,2),(1,2),(1,3),(2,3),(1,4),(1,5),(3,5),(1,6),(3,6)],
    'N' => [(1,1),(1,2),(1,3),(1,4),(1,5),(1,6),(2,2),(3,3),(3,1),(3,2),(3,4),(3,5),(3,6)],
    'O' => [(2,1),(1,2),(3,2),(1,3),(3,3),(1,4),(3,4),(1,5),(3,5),(2,6),(1,6),(3,6)],
    'P' => [(1,1),(2,1),(1,2),(3,2),(1,3),(2,3),(1,4),(1,5),(1,6)],
    'R' => [(1,1),(2,1),(1,2),(3,2),(1,3),(2,3),(1,4),(3,4),(1,5),(3,5),(1,6),(3,6)],
    'S' => [(2,1),(3,2),(1,2),(2,3),(3,4),(3,5),(1,5),(2,6),(1,6),(3,6)],
    'T' => [(1,1),(2,1),(3,1),(2,2),(2,3),(2,4),(2,5),(2,6)],
    'U' => [(1,1),(3,1),(1,2),(3,2),(1,3),(3,3),(1,4),(3,4),(1,5),(3,5),(1,6),(2,6),(3,6)],
    'V' => [(1,1),(3,1),(1,2),(3,2),(1,3),(3,3),(2,4),(2,5),(2,6)],
    'Y' => [(1,1),(3,1),(1,2),(3,2),(2,3),(2,4),(2,5),(2,6)],
    'a' => [(2,2),(1,3),(3,3),(2,4),(1,5),(3,5),(2,6),(1,6),(3,6)],
    'b' => [(1,1),(1,2),(1,3),(1,4),(1,5),(1,6),(2,2),(3,3),(3,4),(2,5)],
    'c' => [(2,2),(3,3),(1,3),(1,4),(1,5),(2,6),(3,6)],
    'd' => [(3,1),(3,2),(3,3),(3,4),(3,5),(3,6),(2,2),(1,3),(1,4),(1,5),(1,6)],
    'e' => [(2,2),(1,3),(3,3),(1,4),(2,4),(3,4),(1,5),(2,6),(3,6)],
    'g' => [(2,2),(1,3),(3,3),(1,4),(2,4),(3,5),(1,5),(2,6),(3,6)],
    'h' => [(1,1),(1,2),(1,3),(1,4),(1,5),(1,6),(2,2),(3,3),(3,4),(3,5),(3,6)],
    'i' => [(2,1),(2,3),(2,4),(2,5),(2,6)],
    'k' => [(1,1),(1,2),(1,3),(1,4),(1,5),(1,6),(3,3),(2,4),(3,5)],
    'l' => [(2,1),(2,2),(2,3),(2,4),(2,5),(2,6)],
    'm' => [(1,2),(2,2),(1,3),(1,4),(2,4),(1,5),(3,3),(3,5),(1,6),(3,6)],
    'n' => [(1,2),(1,3),(1,4),(1,5),(1,6),(2,2),(3,3),(3,4),(3,5),(3,6)],
    'o' => [(2,2),(1,3),(3,3),(1,4),(3,4),(1,5),(3,5),(2,6)],
    'p' => [(1,2),(1,3),(1,4),(1,5),(1,6),(2,2),(3,3),(2,4)],
    'r' => [(1,2),(1,3),(1,4),(1,5),(1,6),(2,2),(3,3)],
    's' => [(2,2),(1,3),(2,3),(3,4),(1,5),(2,6),(3,6),(3,5)],
    't' => [(2,1),(2,2),(2,3),(2,4),(2,5),(2,6),(1,2),(3,2)],
    'u' => [(1,2),(3,2),(1,3),(3,3),(1,4),(3,4),(1,5),(3,5),(1,6),(2,6),(3,6)],
    'v' => [(1,2),(3,2),(1,3),(3,3),(2,4),(2,5)],
    'y' => [(1,2),(3,2),(1,3),(3,3),(2,4),(2,5),(3,6)],
)
const FS = 2                                   # font scale: each font pixel is FS x FS
draw_text(s, x0, y0, col) = (for (ci, ch) in enumerate(s)
    for (dx, dy) in get(DIG, ch, Tuple{Int,Int}[])
        for sy in 0:FS-1, sx in 0:FS-1
            putc(x0 + (ci - 1) * 4 * FS + (dx - 1) * FS + sx, y0 + (dy - 1) * FS + sy, col)
        end
    end
end)
text_w(s) = length(s) * 4 * FS

# ---- one model slice ---------------------------------------------------------
function plot_model(model)
    haskey(MARS_MODELS, model) || error("unknown model $model")
    layers = MARS_MODELS[model]
    moho = MARS_MOHO[model]
    fill!(pix, UInt8(255))

    # layer fill, colored by vp (per pixel row: no gaps at thin layers)
    for y in (TM + 1):(TM + PH - 1)
        z_km = (y - TM) / PH * ZMAX_KM
        col = let c = (0xff, 0xff, 0xff)
            for L in layers
                if L[1] <= z_km < L[2]
                    c = vpcolor(L[3])
                    break
                end
            end
            c
        end
        hlinec(y, LM + 1, LM + PW - 1, col)
    end

    # crust interfaces (thin) and Moho (thick)
    for b in mars_boundaries(model)
        hlinec(py(b), LM, LM + PW, GREY)
    end
    for t in 0:2
        hlinec(py(moho) + t, LM, LM + PW, BLACK)
    end
    draw_text(@sprintf("Moho %.1f km", moho), LM + 8, py(moho) - 14, BLACK)

    # receiver array: 13 triangles at the free surface, x 53-65 km
    for xr in 53.0:1.0:65.0
        x0 = px(xr)
        for d in 0:5
            for xx in (x0 - d):(x0 + d)
                putc(xx, TM + 1 + d, RED)
            end
        end
    end
    draw_text("array", px(67.0), TM + 2, BLACK)

    # source line at z=100 km, x 10-82 km
    for t in 0:2
        hlinec(py(100.0) + t, px(10.0), px(82.0), BLACK)
    end
    draw_text("src line", px(10.0), py(100.0) - 16, BLACK)

    # slant arrow: 10 deg off vertical, rising toward smaller x
    let x0 = px(105.0), y0 = py(85.0), L = 60.0
        x1 = round(Int, x0 - L * sind(SLANT_DEG))      # rising toward -x
        y1 = round(Int, y0 - L * cosd(SLANT_DEG))
        linec(x0, y0, x1, y1, BLACK, 2)
        # arrowhead
        for a in (-25.0, 25.0)
            hx = x1 + round(Int, 14 * sind(SLANT_DEG + a))
            hy = y1 + round(Int, 14 * cosd(SLANT_DEG + a))
            linec(x1, y1, hx, hy, BLACK, 2)
        end
        draw_text("slant 10 deg", x1 + 6, y1 - 16, BLACK)
    end

    # schematic conversion ray: array centre down to its Moho piercing
    # point; the rising plane wave is 10 deg off vertical, so downward
    # from the surface the ray drifts to larger x
    tan_s = tand(SLANT_DEG)
    for y in (TM + 1):(py(moho) - 1)
        z_km = (y - TM) / PH * ZMAX_KM
        x_km = 59.0 + z_km * tan_s
        (y % 6) < 3 && putc(px(x_km), y, WHITE)
    end
    for dy in -2:2, dx in -2:2
        putc(px(59.0 + moho * tan_s) + dx, py(moho) + dy, WHITE)
    end
    draw_text("Ps conv", px(59.0 + moho * tan_s) + 8, py(moho) - 14, WHITE)

    # frame + ticks
    hlinec(TM, LM, LM + PW, BLACK)
    hlinec(TM + PH, LM, LM + PW, BLACK)
    vlinec(LM, TM, TM + PH, BLACK)
    vlinec(LM + PW, TM, TM + PH, BLACK)
    for v in 0:20:120
        vlinec(px(v), TM + PH, TM + PH + 5, BLACK)
        draw_text(string(v), px(v) - 12, TM + PH + 8, BLACK)
        hlinec(py(v), LM - 5, LM, BLACK)
        draw_text(string(v), LM - 44, py(v) - 6, BLACK)
    end
    draw_text("range km", LM + PW - text_w("range km") - 4, TM + PH + 22, BLACK)
    # vertical "depth km" label
    let s = "depth km", y0 = TM + div(PH - length(s) * 10, 2)
        for (ci, ch) in enumerate(s)
            draw_text(string(ch), LM - 54, y0 + (ci - 1) * 10, BLACK)
        end
    end

    # vp colorbar
    let cb_y = H - 34, cb_h = 12
        for xx in 0:(PW - 1)
            c = vpcolor(2000.0 + xx / (PW - 1) * 7000.0)
            for yy in 0:(cb_h - 1)
                putc(LM + 1 + xx, cb_y + yy, c)
            end
        end
        hlinec(cb_y - 1, LM, LM + PW, BLACK)
        hlinec(cb_y + cb_h, LM, LM + PW, BLACK)
        vlinec(LM, cb_y - 1, cb_y + cb_h, BLACK)
        vlinec(LM + PW, cb_y - 1, cb_y + cb_h, BLACK)
        for v in 3000:2000:9000
            draw_text(string(v), LM + round(Int, (v - 2000) / 7000 * PW) - 12, cb_y + cb_h + 6, BLACK)
        end
        draw_text("vp m/s", W - text_w("vp m/s") - 4, cb_y - 16, BLACK)
    end

    # title
    draw_text(model * " slice", LM, 8, BLACK)

    out = joinpath(@__DIR__, "..", "results", "models", "$(model)_slice.png")
    mkpath(dirname(out))
    write_png(out, W, H, pix)
    out
end

for model in MODELS
    out = plot_model(model)
    @printf("wrote %s\n", out)
end
