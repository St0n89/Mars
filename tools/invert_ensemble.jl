# ===========================================================================
# ENSEMBLE-SPACE INVERSION — synthetic recovery over the published InSight
# RF crust posterior (Knapmeyer-Endrun et al. 2021).
#
# Model space: all 10000 invA/invB posterior models (the published
# posterior is our prior; uniform weight per model).
# Forward operator: results/ray_operator.csv (tools/ray_operator.jl,
# validated against 18 FDTD runs).
# Data: the picks from the FDTD sweep (runs/<truth>/summary.csv) — an
# end-to-end test: solver -> matched filter -> inversion.
#
# Likelihood per model i:
#   detection:  P_pick_i = 0.95 if pkg_onset_i >= PULSE_SPAN else 0.05
#     data = package pick at d_obs (matched filter, sigma 1.2 s absorbs
#            the 0.0-1.4 s pick bias seen in the sweep):
#       L_i = P_pick_i * exp(-0.5((pkg_onset_i-d_obs)/sigma)^2) * Lp_i
#     data = NO pick (package buried in the P coda):
#       L_i = (1 - P_pick_i) * Lp_i
#   Lp_i = exp(-0.5((p_i-p_obs)/sigma_p)^2), sigma_p = 2% (moveout fit)
#
# The no-pick case is informative: it EXCLUDES the resolvable (thick-crust)
# models. The pick case localises Moho within the resolvable subset.
#
# Usage:   julia tools/invert_ensemble.jl        (all scenarios)
# Writes:  results/inversion/<truth>_posterior.csv, posterior.png,
#          results/inversion/inversion_summary.csv
# ===========================================================================
using Printf

const PULSE_SPAN = 4.0 + 1.5          # s, 0.5 Hz coda span
const SIG_PKG    = 1.2                # s, matched-filter pick scatter
const SIG_P_FRAC = 0.02               # moveout slowness scatter
const P_DETECT   = 0.95               # detection probability if resolvable

const TRUTHS = ["invA_r8", "invB_r8", "invB_r5", "invA_r1"]

# ---- load operator ----------------------------------------------------------
rows = []
open(joinpath(@__DIR__, "..", "results", "ray_operator.csv")) do fh
    hdr = split(readline(fh), ',')
    idx = Dict(k => i for (i, k) in enumerate(hdr))
    for ln in eachline(fh)
        v = split(ln, ',')
        push!(rows, (ens = v[idx["ensemble"]],
                     moho = parse(Float64, v[idx["moho_km"]]),
                     p = parse(Float64, v[idx["p_s_m"]]),
                     pkg = parse(Float64, v[idx["pkg_onset_s"]]),
                     psm = parse(Float64, v[idx["psm_s"]]),
                     res = parse(Int, v[idx["resolvable"]]) == 1))
    end
end
@printf("operator: %d models (%d resolvable, %.1f%%)\n",
        length(rows), count(r -> r.res, rows), 100 * count(r -> r.res, rows) / length(rows))

# ---- load one truth run's summary.csv ---------------------------------------
function read_summary(run)
    path = joinpath(@__DIR__, "..", "runs", run, "summary.csv")
    lines = readlines(path)
    k = split(lines[1], ',')
    v = split(lines[2], ',')
    d = Dict{String,String}()
    for (kk, vv) in zip(k, v)
        d[kk] = vv
    end
    d
end

# second observable: the deepest-interface Ps delay (Ps@b2), picked only
# when its predicted delay clears the P coda (delay_pred >= PULSE_SPAN).
# Its model prediction is ps_del(b_deep) = psm - pkg_onset.
function read_ps_b2(run)
    path = joinpath(@__DIR__, "..", "runs", run, "picks.csv")
    best = nothing
    for ln in readlines(path)[2:end]
        f = split(ln, ',')
        startswith(f[1], "Ps") && f[1] != "PsPkg" || continue
        dp = parse(Float64, f[2])
        dp >= PULSE_SPAN || continue
        best === nothing || dp > best[2] || continue
        best = (f[1], dp, parse(Float64, f[3]))
    end
    best === nothing ? (hasb2 = false, d2 = 0.0) : (hasb2 = true, d2 = best[3])
end

# ---- posterior over Moho ----------------------------------------------------
const NBIN = 90                        # 1 km bins, Moho 0-90 km
binof(m) = clamp(floor(Int, m) + 1, 1, NBIN)

function posterior(obs)
    # obs = (haspick, d_obs, p_obs, hasb2, d2_obs)
    sp = SIG_P_FRAC * obs.p_obs
    w = zeros(length(rows))
    for (i, r) in enumerate(rows)
        lp = exp(-0.5 * ((r.p - obs.p_obs) / sp)^2)
        obs.hasb2 && (lp *= exp(-0.5 * ((r.psm - r.pkg - obs.d2_obs) / SIG_PKG)^2))
        pdet = r.res ? P_DETECT : (1.0 - P_DETECT)
        if obs.haspick
            w[i] = pdet * exp(-0.5 * ((r.pkg - obs.d_obs) / SIG_PKG)^2) * lp
        else
            w[i] = (1.0 - pdet) * lp
        end
    end
    s = sum(w)
    s > 0 || error("zero evidence")
    w ./= s
    post = zeros(NBIN)
    pres = 0.0
    for (i, r) in enumerate(rows)
        post[binof(r.moho)] += w[i]
        r.res && (pres += w[i])
    end
    prior = zeros(NBIN)
    for r in rows
        prior[binof(r.moho)] += 1.0 / length(rows)
    end
    (prior = prior, post = post, p_resolvable = pres)
end

function cred(post)
    # discrete percentiles over bin centres (km)
    c = cumsum(post)
    tot = c[end]
    q(f) = begin
        t = f * tot
        for i in eachindex(c)
            c[i] >= t && return i - 0.5
        end
        NBIN - 0.5
    end
    (lo95 = q(0.025), lo68 = q(0.16), med = q(0.5),
     hi68 = q(0.84), hi95 = q(0.975))
end

# ---- PNG writer (CRC32 + Adler32 + stored DEFLATE, same as analyze_mars.jl) --
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
    a = adler32(raw)
    i = 1
    blocks = UInt8[]
    while true
        n = length(raw) - i + 1
        fin = i + n - 1 >= length(raw)
        n > 65535 && (n = 65535; fin = false)
        push!(blocks, fin ? 0x01 : 0x00)
        append!(blocks, be32(n)[1:2])
        append!(blocks, be32(n ⊻ 0xffffffff)[1:2])
        append!(blocks, @view raw[i:i+n-1])
        i += n
        fin && break
    end
    vcat(blocks, be32(a))
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

# ---- plot --------------------------------------------------------------------
const W, H = 880, 460
const LM, RM, TM, BM = 56, 24, 30, 46
const PW, PH = W - LM - RM, H - TM - BM
const WHITE = UInt8[255, 255, 255]
const BLACK = UInt8[0, 0, 0]
const GREY  = UInt8[170, 170, 170]
const BLUE  = UInt8[30, 60, 200]
const RED   = UInt8[200, 30, 30]
const AXS   = UInt8[200, 200, 200]
pix = fill(UInt8(255), H, W, 3)
putc(x, y, c) = (1 <= x <= W && 1 <= y <= H && (pix[y, x, 1] = c[1]; pix[y, x, 2] = c[2]; pix[y, x, 3] = c[3]))
vlinec(x, y0, y1, c) = (for y in min(y0, y1):max(y0, y1); putc(x, y, c); end)
hlinec(y, x0, x1, c) = (for x in min(x0, x1):max(x0, x1); putc(x, y, c); end)
function bar(x0, x1, y0, y1, c)      # filled bar clipped at y1 (bottom)
    for x in x0:x1
        for y in y1:y0                # y grows downward; y1 = base
            putc(x, y, c)
        end
    end
end
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
draw_text(s, x0, y0, col) = (for (ci, ch) in enumerate(s)
    for (dx, dy) in get(DIG, ch, Tuple{Int,Int}[])
        putc(x0 + (ci - 1) * 4 + dx - 1, y0 + dy - 1, col)
    end
end)

function plot_posterior(path, prior, post, truth_moho, haspick)
    fill!(pix, UInt8(255))
    xm(m) = LM + round(Int, m / 90.0 * PW)
    ymax = maximum(vcat(prior, post)) * 1.1
    yh(v) = TM + PH - round(Int, v / ymax * PH)
    # axes
    hlinec(TM + PH, LM, LM + PW, BLACK)
    vlinec(LM, TM, TM + PH, BLACK)
    for m in 0:10:90
        vlinec(xm(m), TM + PH, TM + PH + 5, BLACK)
        m > 0 && draw_text(string(m), xm(m) - 6, TM + PH + 10, BLACK)
        vlinec(xm(m), TM, TM + PH, AXS)
    end
    draw_text("km", LM + PW - 10, TM + PH + 18, BLACK)
    base = TM + PH
    for i in 1:NBIN
        x0 = LM + round(Int, (i - 1) / NBIN * PW)
        x1 = LM + round(Int, i / NBIN * PW) - 1
        prior[i] > 0 && bar(max(x0, LM) + 1, max(x1, LM), base, max(yh(prior[i]), TM + 1), GREY)
        post[i]   > 0 && bar(max(x0, LM) + 1, max(x1, LM), base, max(yh(post[i]), TM + 1), BLUE)
    end
    xt = xm(truth_moho)
    for y in TM:base
        putc(xt, y, RED)
        putc(xt + 1, y, RED)
    end
    draw_text(haspick ? "pkg pick" : "no pick", LM + 8, TM + 8, haspick ? BLUE : GREY)
    write_png(path, W, H, pix)
end

# ---- run all scenarios ---------------------------------------------------------
mkpath(joinpath(@__DIR__, "..", "results", "inversion"))
summ_path = joinpath(@__DIR__, "..", "results", "inversion", "inversion_summary.csv")
open(summ_path, "w") do sh
    println(sh, "truth,truth_moho_km,haspick,has_ps_b2,d_obs_s,ps_b2_obs_s,p_obs_s_m,mode_km,median_km,ci68_lo,ci68_hi,ci95_lo,ci95_hi,truth_in_ci68,truth_in_ci95,p_resolvable_model")
    for truth in TRUTHS
        s = read_summary(truth)
        onset_pred = parse(Float64, s["pkg_onset_pred"])
        p_obs = parse(Float64, s["p_meas"])
        haspick = onset_pred >= PULSE_SPAN
        d_obs = haspick ? parse(Float64, s["psm_delay_meas"]) : 0.0
        truth_moho = parse(Float64, s["moho_km"])
        b2 = read_ps_b2(truth)
        pr = posterior((haspick = haspick, d_obs = d_obs, p_obs = p_obs,
                        hasb2 = b2.hasb2, d2_obs = b2.d2))
        ci = cred(pr.post)
        mode = argmax(pr.post) - 0.5
        in68 = ci.lo68 <= truth_moho <= ci.hi68
        in95 = ci.lo95 <= truth_moho <= ci.hi95
        @printf("%-9s truth Moho %.1f km  %s%s  ->  mode %.1f, median %.1f, 68%% [%.1f, %.1f], 95%% [%.1f, %.1f], P(resolvable)=%.3f\n",
                truth, truth_moho,
                haspick ? "pick @ $(round(d_obs, digits=2))s" : "NO pick",
                b2.hasb2 ? " + Ps@b2 @ $(round(b2.d2, digits=2))s" : "",
                mode, ci.med, ci.lo68, ci.hi68, ci.lo95, ci.hi95, pr.p_resolvable)
        csv = joinpath(@__DIR__, "..", "results", "inversion", "$(truth)_posterior.csv")
        open(csv, "w") do fh
            println(fh, "moho_bin_km,n_models_prior,prior_dens,post_dens")
            for i in 1:NBIN
                n = count(r -> binof(r.moho) == i, rows)
                @printf(fh, "%.1f,%d,%.6f,%.6f\n", i - 0.5, n, pr.prior[i], pr.post[i])
            end
        end
        png = joinpath(@__DIR__, "..", "results", "inversion", "$(truth)_posterior.png")
        plot_posterior(png, pr.prior, pr.post, truth_moho, haspick)
        @printf(sh, "%s,%.1f,%d,%d,%.4f,%.4f,%.4e,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%d,%d,%.4f\n",
                truth, truth_moho, haspick ? 1 : 0, b2.hasb2 ? 1 : 0,
                d_obs, b2.d2, p_obs,
                mode, ci.med, ci.lo68, ci.hi68, ci.lo95, ci.hi95,
                in68 ? 1 : 0, in95 ? 1 : 0, pr.p_resolvable)
    end
    println("wrote results/inversion/")
end
