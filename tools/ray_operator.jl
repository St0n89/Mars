# ===========================================================================
# RAY OPERATOR over the full InSight RF crust ensembles (forward step of the
# ensemble-space inversion).
#
# For every invA/invB posterior model (5000 each, Knapmeyer-Endrun et al.
# 2021) compute, at a fixed incidence slant, the SAME observables that
# analyze_mars.jl picks from the FDTD runs:
#
#   p          = sin(SLANT)/vp(100 km)   horizontal slowness (per-model,
#                                        same convention as run_mars.jl)
#   tP         = T_P(100 km -> receiver)
#   pkg_onset  = ps_del(moho) - ps_del(b_deep)   conversion-package onset
#                                        (PsP at the deepest crust interface)
#   psm        = ps_del(moho)            PsM member delay
#   ppmp_b2    = 2*T_P(b_deep)           P multiple off the deepest interface
#   ppmpm      = 2*T_P(moho)            P multiple off the Moho
#
# Layers are quantised exactly like config/models_gen.jl (200 m, clipped at
# 120 km), so operator predictions and FDTD runs see the same model. The
# 18-run sweep validated these formulas against the solver: slowness within
# 0.6%, resolvable-package onsets within 0.0-1.4 s (matched-filter bias).
#
# Usage:  julia tools/ray_operator.jl [MARS_SLANT=10]
# Writes: results/ray_operator.csv (one row per model, pooled ensembles)
# ===========================================================================
using Printf

const SLANT = parse(Float64, get(ENV, "MARS_SLANT", "10"))
const Z_REC  = 0.1                            # receiver depth, km
const PULSE_SPAN = 4.0 + 1.5                   # 0.5 Hz coda span, s (f0=0.5)

include(joinpath(@__DIR__, "ensemble_pick.jl"))  # parse_invA/parse_invB, q02

eta(v, p_) = sqrt(1.0 / (v * v) - p_ * p_)     # vertical slowness, s/m

# quantise + clip one model into env-format layers + crust interfaces
function prep(m)
    moho_q = q02(m.moho)
    layers = Tuple{Float64,Float64,Float64,Float64,Float64}[]
    bnds   = Float64[]
    for L in m.layers
        dtop = q02(max(L[1], 0.0))
        dbot = q02(min(L[2], 120.0))
        dbot <= dtop && continue
        push!(layers, (dtop, dbot, L[3], L[4], L[5]))
        dtop > 0 && dtop < moho_q && push!(bnds, dtop)
    end
    (moho_q, layers, bnds)
end

# fixed-slowness travel time, km -> s (same loop as analyze_mars.jl)
function travel(layers, z0, z1, mode, p_)
    t = 0.0
    for L in layers
        a = max(min(z1, L[2]), L[1])
        b = min(max(z0, L[1]), L[2])
        b > a && (t += (b - a) * 1000.0 *
                  eta(mode == 'P' ? L[3] : L[4], p_))
    end
    t
end

function operator_row(m)
    moho, layers, bnds = prep(m)
    # vp at 100 km (source depth): the layer that contains it
    vp100 = 0.0
    for L in layers
        L[1] <= 100.0 < L[2] && (vp100 = L[3])
    end
    vp100 > 0 || return nothing
    p = sind(SLANT) / vp100
    tP(d)   = travel(layers, d, Z_REC, 'P', p)
    tS(d)   = travel(layers, d, Z_REC, 'S', p)
    ps_del(d)  = tS(d) - tP(d)
    b_deep = isempty(bnds) ? moho : maximum(bnds)
    pkg    = ps_del(moho) - ps_del(b_deep)
    (moho        = moho,
     b1          = isempty(bnds) ? 0.0 : minimum(bnds),
     b2          = b_deep,
     p           = p,
     tP          = tP(100.0),
     pkg_onset   = pkg,
     psm         = ps_del(moho),
     ppmp_b2     = 2.0 * tP(b_deep),
     ppmpm       = 2.0 * tP(moho),
     resolvable  = pkg >= PULSE_SPAN)
end

function main()
    A = parse_invA(INV_A)
    B = parse_invB(INV_B)
    @printf("parsed: %d invA + %d invB models, slant %.1f deg\n",
            length(A), length(B), SLANT)
    out = joinpath(@__DIR__, "..", "results", "ray_operator.csv")
    open(out, "w") do fh
        println(fh, "ensemble,idx,moho_km,b1_km,b2_km,p_s_m,tP_s," *
                    "pkg_onset_s,psm_s,ppmp_b2_s,ppmpm_s,resolvable")
        n = 0
        for (tag, models) in (("invA", A), ("invB", B))
            for (i, m) in enumerate(models)
                r = operator_row(m)
                r === nothing && continue
                @printf(fh, "%s,%d,%.2f,%.2f,%.2f,%.4e,%.4f,%.4f,%.4f,%.4f,%.4f,%d\n",
                        tag, i, r.moho, r.b1, r.b2, r.p, r.tP,
                        r.pkg_onset, r.psm, r.ppmp_b2, r.ppmpm,
                        r.resolvable ? 1 : 0)
                n += 1
            end
        end
        @printf("wrote %s: %d rows\n", out, n)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
