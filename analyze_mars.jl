#!/usr/bin/env julia
# ===========================================================================
# ANALYSIS: Mars FDTD — slanted plane P wave, receiver-function picks
#
# Predicts arrivals with fixed horizontal slowness p (conserved through
# the stack): vertical slowness eta = sqrt(1/v^2 - p^2), so
#   t_P(x_r)   = sum eta_p dz + p (x_r - x_line_center)
#   Ps delay   = sum above interface (eta_s - eta_p) dz
#   PPmP delay = 2 sum above interface eta_p dz
# The Moho P->S conversion is measured as a PACKAGE: the converted S
# re-converts S->P at every crustal interface it crosses on the way up,
# so the observable at the surface is the superposition of PsP@b
# (delay = sum_{b..moho}(eta_s-eta_p) dz) for each crust interface b
# plus PsM (S all the way to the surface). At vertical incidence all
# P<->S conversions vanish (control run); at slant the package onset is
# the S->P re-conversion at the DEEPEST crust interface. Picking a
# single "PsM" arrival against the P reference biases the delay early
# by the package spread (verified on TAYAK 15 deg: measured 9.34 s vs
# PsM 10.82 s vs package onset PsP@10km 9.12 s).
# Onsets are measured per receiver by cross-correlation against that
# receiver's own direct-P envelope, then averaged into delays. Also
# fits the P moveout slope across the array as an independent check on
# p, and records a late-trace envelope ratio as a resonance/QA metric:
# slant runs carry a numerically growing surface-wave train (phase
# speed ~2.5 km/s, e-folding ~1 s, marching outward from the line
# centre) that reaches the array from ~t=38 s; late_env_ratio << 1
# means the run is clean (see config/sr_mars.jl for the full story).
#
# ENV: MARS_RUN (default "mars_run"). Model, Moho, slant and p are read
# from <MARS_RUN>/meta.jl; the layer table from config/models.jl.
# Stdlib only. Writes picks.csv, summary.csv and section.png (hand-rolled
# PNG: CRC32 + Adler32 + stored DEFLATE, no zlib).
# ===========================================================================

using Printf

const RUN_DIR = get(ENV, "MARS_RUN", "mars_run")
base = joinpath(@__DIR__, RUN_DIR)

# ---- meta ------------------------------------------------------------------
meta = Dict{String,Any}()
for ln in readlines(joinpath(base, "meta.jl"))
    p = findfirst('=', ln)
    p === nothing && continue
    k = String(strip(replace(strip(ln[1:p-1]), "const " => "")))
    v = String(strip(ln[p+1:end]))
    meta[k] = startswith(v, "\"") ? String(strip(v, '"')) : parse(Float64, v)
end
dt    = meta["dt"];  f0 = meta["f0"]
p     = meta["p_slant"]; slant = meta["slant_deg"]
MODEL = meta["model"]

include(joinpath(@__DIR__, "config", "models.jl"))
layers = MARS_MODELS[MODEL]
moho   = MARS_MOHO[MODEL]
bnd    = mars_boundaries(MODEL)          # crust interfaces above the Moho
b1     = isempty(bnd) ? 0.0 : bnd[1]
b2     = length(bnd) >= 2 ? bnd[2] : moho

const Z_REC  = 0.1      # receiver depth, km
const Z_SRC  = 100.0     # source depth, km
X_CENTER    = get(meta, "x_center", 30000.0)  # line centre, m (legacy runs: 30 km)

# ---- slowness ray theory ----------------------------------------------------
function eta(v, p_)                       # vertical slowness (s/m)
    s = p_ * v
    abs(s) < 1 || error("post-critical slowness at v=$v")
    sqrt(1.0 / (v * v) - p_ * p_)
end
function travel_slant(z0, z1, mode, p_)    # km -> s, fixed horizontal slowness
    t = 0.0
    for L in layers
        a = max(min(z1, L[2]), L[1])
        b = min(max(z0, L[1]), L[2])
        b > a && (t += (b - a) * 1000.0 * eta(mode == 'P' ? L[3] : L[4], p_))
    end
    t
end

ps_del(d)  = travel_slant(d, Z_REC, 'S', p) - travel_slant(d, Z_REC, 'P', p)
ppm_del(d) = 2.0 * travel_slant(d, Z_REC, 'P', p)
T_P_col    = travel_slant(Z_SRC, Z_REC, 'P', p)

lbl(d) = replace(@sprintf("%.1f", d), "." => "p")
# deepest crustal interface above the Moho; the conversion package onset
# is the S->P re-conversion there (equals PsM delay for 1-layer crusts)
b_deep = isempty(bnd) ? 0.0 : bnd[end]
pkg_onset = isempty(bnd) ? ps_del(moho) : ps_del(moho) - ps_del(b_deep)
psm_pred  = ps_del(moho)

PHASES = [("P", 0.0, "direct P from 100 km")]
isempty(bnd) || push!(PHASES, ("Ps$(lbl(b1))", ps_del(b1), "Ps at $(b1) km"))
length(bnd) >= 2 && push!(PHASES, ("Ps$(lbl(b2))", ps_del(b2), "Ps at $(b2) km"))
push!(PHASES, ("PsPkg", pkg_onset,
    isempty(bnd) ? "Moho conversion package (PsM only, 1-layer crust)" :
    "Moho conversion package: S->P at $(b_deep) km (onset) .. PsM $(@sprintf("%.2f", psm_pred)) s"))
length(bnd) >= 2 && push!(PHASES, ("PPmP$(lbl(b2))", ppm_del(b2), "multiple off $(b2) km"))
push!(PHASES, ("PPmPM", ppm_del(moho), "PPmP Moho multiple"))

# ---- traces -----------------------------------------------------------------
lines = readlines(joinpath(base, "traces.csv"))
data  = [split(ln, ',') for ln in lines if !isempty(ln)][2:end]
t     = [parse(Float64, r[1]) for r in data]
nt    = length(t)
hdr   = split(lines[1], ',')
xrec  = [parse(Float64, match(r"recv_(\d+)_m", h).captures[1]) for h in hdr[2:end]]
nrec  = length(xrec)
traces = [[parse(Float64, r[j + 1]) for r in data] for j in 1:nrec]

@printf("model=%s, Moho=%.1f km, slant=%.1f deg, p=%.3e s/m\n", MODEL, moho, slant, p)
p == 0 && @printf("slant=0: P<->S conversions are not excited at vertical incidence; Ps/PsPkg picks are coda noise\n")
@printf("receivers: %d (x=%.0f-%.0f m), nt=%d, dt=%.5f s, f0=%.1f Hz\n\n",
        nrec, xrec[1], xrec[end], nt, dt, f0)

# ---- envelope + per-receiver cross-correlation -------------------------------
PULSE_SPAN = 4.0 + 1.5

function envelope(tr)
    w = max(1, round(Int, 0.5 / dt))
    out = Vector{Float64}(undef, nt)
    s = 0.0
    for i in 1:nt
        s += abs(tr[i])
        i > w && (s -= abs(tr[i - w]))
        out[i] = s / min(i, w)
    end
    out
end

idxf(tv) = clamp(round(Int, tv / dt), 0, nt - 1)   # 0-based

# lag maximising normalised correlation of env with its own P reference;
# returns (onset, peak env in the aligned window). t_ref0 = start of the P
# reference wavelet (tP_j - 1); window starts range +/- half around
# t_search; tP_j = this receiver's predicted direct-P onset.
function xcorr_pick(env, t_ref0, t_search, tP_j; half = 2.5)
    k0 = idxf(t_ref0); k1 = idxf(t_ref0 + 5.5)
    K = k1 - k0
    K <= 0 && return nothing, 0.0
    r = env[k0+1:k1]
    rr = 0.0
    for v in r; rr += v * v; end
    rr = sqrt(rr)
    rr == 0.0 && return nothing, 0.0
    s0 = idxf(t_search - half); s1 = idxf(t_search + half + 5.5) - K
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
    t[best+1] + (tP_j - t_ref0), maximum(env[best+1:best+K])
end

envs = envelope.(traces)
meanv(v) = let m = 0.0
    for u in v; m += u; end
    m / length(v)
end

# matched filter for the Moho conversion package: template = the measured
# P-wavelet envelope, placed at every predicted member onset (PsP at each
# crust interface + PsM), cross-correlated over +-1.5 s shift. The bare
# P-wavelet correlation fails here because the PPmP coda tail overlaps the
# package onset and biases the pick late (verified on TAYAK slant 10 deg).
function matched_pick(env, tP_j, mdelays, onset0)
    Lp = round(Int, 8.0 / dt)              # Hann pulse length
    i0 = idxf(tP_j)                        # P onset, 0-based
    (i0 + Lp > nt) && return nothing, 0.0
    pshape = env[i0+1:i0+Lp]
    maximum(pshape) <= 0 && return nothing, 0.0
    m0, m1 = minimum(mdelays), maximum(mdelays)
    Ls = round(Int, (m1 - m0) / dt) + Lp   # template support
    Td = zeros(Ls)
    for m in mdelays
        off = round(Int, (m - m0) / dt)
        for i in 1:Lp
            Td[off+i] += pshape[i]
        end
    end
    nT = sqrt(sum(v -> v * v, Td))
    nT == 0 && return nothing, 0.0
    best_sc = -1.0; best_d = 0.0
    for k in round(Int, -1.5 / dt):round(Int, 1.5 / dt)
        s = idxf(tP_j + m0) + k            # 0-based segment start
        (s < 0 || s + Ls > nt) && continue
        seg = env[s+1:s+Ls]
        ns = sqrt(sum(v -> v * v, seg))
        ns == 0 && continue
        sc = 0.0
        for i in 1:Ls
            sc += Td[i] * seg[i]
        end
        sc /= (nT * ns)
        sc > best_sc && (best_sc = sc; best_d = k * dt)
    end
    (tP_j + onset0 + best_d,
     maximum(env[idxf(tP_j + m0 + best_d)+1:min(idxf(tP_j + m1 + best_d + 8.0)+1, nt)]))
end

# per-receiver picks: onset_j(phase), peak_j(phase)
mdelays = sort(unique(vcat([psm_pred - ps_del(b) for b in bnd], [psm_pred])))
onsets = Dict{String,Vector{Float64}}()
peaks  = Dict{String,Vector{Float64}}()
for (name, delay, _) in PHASES
    on = Float64[]; pk = Float64[]
    for (j, env) in enumerate(envs)
        tP_j = T_P_col + p * (xrec[j] - X_CENTER)
        if name == "PsPkg"
            to, pv = matched_pick(env, tP_j, mdelays, pkg_onset)
        else
            to, pv = xcorr_pick(env, tP_j - 1.0, tP_j + delay, tP_j)
        end
        to !== nothing && (push!(on, to); push!(pk, pv))
    end
    onsets[name] = on; peaks[name] = pk
end

# ---- P: absolute onset, moveout ----------------------------------------------
onP = onsets["P"]
tP_meas = let
    m = 0.0
    for (j, v) in enumerate(onP); m += v - p * (xrec[j] - X_CENTER); end
    m / length(onP)
end
p_meas = let
    mx = meanv(xrec); my = meanv(onP)
    sxx = sum(v -> (v - mx)^2, xrec)
    sxy = 0.0
    for j in eachindex(xrec); sxy += (xrec[j] - mx) * (onP[j] - my); end
    sxy / sxx
end

# ---- table -------------------------------------------------------------------
fmt2(v) = v === nothing ? "  -  " : @sprintf("%.2f", v)
@printf("%-9s %9s %9s %8s %10s %10s %9s %8s  %s\n",
        "phase", "delay pd", "delay ms", "err s", "onset pd", "onset ms",
        "amp/P", "peakenv", "note")
picks_rows = ["phase,delay_pred_s,delay_meas_s,delay_err_s,onset_pred_s,onset_meas_s,amp_over_P,peak_env,note"]
summary = Dict{String,Any}("model" => MODEL,
                          "ensemble" => startswith(MODEL, "invA") ? "invA" :
                          startswith(MODEL, "invB") ? "invB" : "TAYAK",
                          "b1_km" => b1, "b2_km" => b2, "moho_km" => moho,
                          "slant_deg" => slant, "p_pred" => p, "p_meas" => p_meas,
                          "tP_pred" => T_P_col, "tP_meas" => tP_meas,
                          "tP_err" => tP_meas - T_P_col)
pkP = isempty(peaks["P"]) ? 0.0 : meanv(peaks["P"])
for (name, delay, note) in PHASES
    on = onsets[name]
    pkv = isempty(peaks[name]) ? 0.0 : meanv(peaks[name])
    # receiver-local delay: phase onset minus that receiver's P onset
    if name == "P" || isempty(on) || isempty(onP) || length(on) != length(onP)
        d_meas = nothing
    else
        d_meas = let
            m = 0.0
            for i in eachindex(on); m += on[i] - onP[i]; end
            m / length(on)
        end
    end
    amp = name == "P" ? 1.0 : (pkP > 0 ? pkv / pkP : 0.0)
    onset_pred = name == "P" ? T_P_col : T_P_col + delay
    onset_meas = isempty(on) ? nothing : meanv(on)
    if name != "P" && delay < PULSE_SPAN
        flag = "inside direct-P coda (unresolved)"
    elseif d_meas !== nothing
        flag = note * @sprintf(", err %+.2f s", d_meas - delay)
    else
        flag = note
    end
    @printf("%-9s %9.2f %9s %8s %10.2f %10s %9.3f %8.1e  %s\n",
            name, delay, fmt2(d_meas),
            fmt2(d_meas === nothing ? nothing : d_meas - delay),
            onset_pred, fmt2(onset_meas), amp, pkv, flag)
    push!(picks_rows, join((
        name, @sprintf("%.3f", delay),
        d_meas === nothing ? "" : @sprintf("%.3f", d_meas),
        d_meas === nothing ? "" : @sprintf("%.3f", d_meas - delay),
        @sprintf("%.3f", onset_pred),
        onset_meas === nothing ? "" : @sprintf("%.3f", onset_meas),
        @sprintf("%.3f", amp), @sprintf("%.3e", pkv), flag), ","))
    if name == "PsPkg"
        summary["pkg_onset_pred"] = pkg_onset
        summary["psm_delay_pred"] = psm_pred
        summary["psm_delay_meas"] = d_meas === nothing ? "" : d_meas
        summary["pkg_err"] = d_meas === nothing ? "" : d_meas - pkg_onset
        summary["psm_err"] = d_meas === nothing ? "" : d_meas - psm_pred
        summary["psm_over_P"] = amp
    end
    if name == "PPmPM"
        summary["ppmpm_delay_pred"] = delay
        summary["ppmpm_delay_meas"] = d_meas === nothing ? "" : d_meas
        summary["ppmpm_err"] = d_meas === nothing ? "" : d_meas - delay
        summary["ppmpm_over_P"] = amp
    end
end
open(joinpath(base, "picks.csv"), "w") do io
    println(io, join(picks_rows, "\n"))
end

# late-trace envelope ratio: growing-surface-wave QA metric
late_ratio = let
    m = 0.0; n = 0
    for env in envs
        pk = maximum(env[idxf(tP_meas - 1)+1:idxf(tP_meas + 4)+1])
        pk > 0 || continue
        m += env[nt] / pk; n += 1
    end
    n > 0 ? m / n : 0.0
end
summary["late_env_ratio"] = late_ratio

# summary.csv (one row; sweep aggregation = cat across run dirs)
skeys = ["model", "ensemble", "b1_km", "b2_km", "moho_km", "slant_deg",
         "p_pred", "p_meas", "tP_pred", "tP_meas", "tP_err",
         "pkg_onset_pred", "psm_delay_pred", "psm_delay_meas",
         "pkg_err", "psm_err", "psm_over_P",
         "ppmpm_delay_pred", "ppmpm_delay_meas", "ppmpm_err", "ppmpm_over_P",
         "late_env_ratio"]
open(joinpath(base, "summary.csv"), "w") do io
    println(io, join(skeys, ","))
    fmtv(v) = v isa AbstractString ? v :
              abs(v) < 1e-3 ? @sprintf("%.4e", v) : @sprintf("%.4f", v)
    println(io, join([fmtv(summary[k]) for k in skeys], ","))
end
@printf("moveout slope: p_meas = %.3e s/m (pred %.3e), tP err %+.2f s, late_env_ratio = %.2f\n",
        p_meas, p, tP_meas - T_P_col, late_ratio)

# ===========================================================================
# RECORD-SECTION PNG — hand-rolled writer (CRC32, Adler32, stored DEFLATE)
# ===========================================================================
SPP = 16           # px per second
RH  = 34           # px per trace row
LM  = 6            # left margin
W = LM + Int(floor(t[end] * SPP))
H = nrec * RH + 40

const BG  = (0x18, 0x18, 0x1e)
const POS = (0x3c, 0x5a, 0xc8)
const NEG = (0xc8, 0x5a, 0x3c)
const GRN = (0x00, 0xaa, 0x3c)
const AXS = (0x3c, 0x3c, 0x46)
const TXT = (0x96, 0x96, 0xa0)
const WIG = (0xeb, 0xeb, 0xeb)

img = fill(BG, W, H)
put0(x, y, c) = (0 <= x < W && 0 <= y < H) && (img[x+1, y+1] = c)

gmax = maximum(maximum(abs.(tr)) for tr in traces)
clip = gmax / 3.0
for (j, tr) in enumerate(traces)
    y0 = (j - 1) * RH + RH ÷ 2
    amp = maximum(abs.(tr))
    amp == 0.0 && continue
    half = RH ÷ 2 - 1
    prev_y = y0
    for (i, tv) in enumerate(t)
        x = LM + floor(Int, tv * SPP)
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

# predicted arrivals: slanted dashed green lines (moveout across the
# array); the conversion package is drawn as a band between its onset
# (S->P at the deepest crust interface) and its PsM member
for (name, delay, _) in PHASES
    for (j, xv) in enumerate(xrec)
        tP_j = T_P_col + p * (xv - X_CENTER)
        xd = LM + floor(Int, (tP_j + delay) * SPP)
        for y in (j - 1) * RH:(j * RH - 1)
            (y ÷ 3) % 2 == 0 && put0(xd, y, GRN)
        end
        if name == "PsPkg"
            xd2 = LM + floor(Int, (tP_j + psm_pred) * SPP)
            for x in min(xd, xd2):max(xd, xd2)
                (x ÷ 3) % 2 == 0 && put0(x, (j - 1) * RH + 2, GRN)
            end
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
function draw_text0(s, x0, y0, col)
    for (ci, ch) in enumerate(s)
        for (dx, dy) in get(DIG, ch, Tuple{Int,Int}[])
            put0(x0 + (ci - 1) * 4 + dx - 1, y0 + dy - 1, col)
        end
    end
end

for (j, xv) in enumerate(xrec)
    draw_text0(string(round(Int, xv / 1000)) * "km", 0, (j - 1) * RH + RH ÷ 2 - 4, TXT)
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
        raw[o += 1] = 0x00
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
ihdr = vcat(be32(W), be32(H), UInt8[8, 2, 0, 0, 0])
png = vcat(UInt8[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a],
           chunk("IHDR", ihdr),
           chunk("IDAT", zlib_stored(raw)),
           chunk("IEND", UInt8[]))
open(joinpath(base, "section.png"), "w") do io
    write(io, png)
end
@printf("section.png: %dx%dpx, picks.csv + summary.csv written\n", W, H)
