using Printf
# ===========================================================================
# ENSEMBLE PICK: turn the InSight joint-Ps-RF crust ensembles
# (Knapmeyer-Endrun et al. 2021, Science 373, doi:10.18715/IPGP.2021.kqwc4yaf)
# into 16 FDTD env models: 8 members per ensemble, evenly spaced in
# Moho-depth rank across the full 5000-member posterior.
#
#   invA (3lay jointPs HF, N=5000): per model "% model N" block:
#     misfit / crustal VP/VS, then rows "depth  VS  VP  rho" (km, km/s,
#     g/cc). Row j's values apply [d_j, d_{j+1}); crust layers are
#     listed twice (top + bottom), giving zero-width rows that carry
#     the interface depths. 3-layer crust -> Moho = depth of row 6.
#     Values ARE (vs, vp, rho); mantle rows are the model's own mantle.
#   invB (3lay jointPs HF, N=5000): one line per model:
#     d1 vs1 r1 d2 vs2 r2 d3 vs3 r3 vs4 r4 misfit   (r = VP/VS)
#     4 layers: [0,d1) [d1,d2) [d2,d3) + mantle half-space below d3;
#     Moho = d3. Densities via the Birch fit to invA:
#     rho (g/cc) = 0.0262 + 0.4389 * vp (km/s).
#
# Output: config/models_gen.jl with ENSEMBLE_MODELS, ENSEMBLE_MOHO,
# ENSEMBLE_BOUNDS — layers (dtop_km, dbot_km, vp_m_s, vs_m_s, rho_kg_m3),
# boundaries 200 m quantised, clipped at 120 km depth.
#
# Run from the repo root:  julia tools/ensemble_pick.jl
# Stdlib only. Data: /home/user/knowledge/insight-mars-fdtd/data/ (durable
# copies) or the dataverse API (see file header comments).
# ===========================================================================

const DATA_DIR = get(ENV, "MARS_ENSEMBLE_DATA",
                     "/home/user/knowledge/insight-mars-fdtd/data")
const INV_A = joinpath(DATA_DIR, "mars_invA_3lay_jointPs_ensemble_N5000_HF.dat")
const INV_B = joinpath(DATA_DIR, "mars_invB_3lay_jointPs_ensemble_N5000_HF.dat")

q02(d) = round(d / 0.2) * 0.2        # 200 m quantisation, km

# ---- invA -------------------------------------------------------------------
# returns per-model NamedTuple: (moho, layers[(dtop,dbot,vs,vp,rho)...])
function parse_invA(path)
    models = []
    cur_rows = Vector{Tuple{Float64,Float64,Float64,Float64}}()  # d, vs, vp, rho
    flush_model!() = begin
        isempty(cur_rows) && return
        # rows: values apply [d_j, d_{j+1}); zero-width rows mark interfaces
        layers = Tuple{Float64,Float64,Float64,Float64,Float64}[]
        interfaces = Float64[]
        for j in 1:length(cur_rows)
            d_top = cur_rows[j][1]
            d_bot = j < length(cur_rows) ? cur_rows[j+1][1] : Inf
            d_bot <= d_top && (push!(interfaces, d_bot); continue)
            push!(layers, (d_top, min(d_bot, 120.0), cur_rows[j][3] * 1000.0,
                           cur_rows[j][2] * 1000.0, cur_rows[j][4] * 1000.0))
        end
        # 3-layer crust: 3 interface rows -> last one is the Moho
        length(interfaces) >= 3 ||
            error("invA: expected 3 crust interfaces, got $(interfaces)")
        moho = interfaces[3]
        push!(models, (moho = moho, layers = layers))
        empty!(cur_rows)
    end
    for ln in readlines(path)
        s = strip(ln)
        isempty(s) && continue
        if startswith(s, "% model")
            flush_model!()
        elseif occursin(r"^ *[0-9]", s) && !occursin(r"[A-Za-z]", s)
            v = parse.(Float64, split(s))
            length(v) == 4 || continue          # skip misfit / vpvs lines
            push!(cur_rows, (v[1], v[2], v[3], v[4]))
        end
    end
    flush_model!()
    models
end

# ---- invB -------------------------------------------------------------------
# one line: d1 vs1 r1 d2 vs2 r2 d3 vs3 r3 vs4 r4 misfit
# README: d_i = layer THICKNESS (km), vs_i = S velocity (km/s), r = vp/vs;
# the model ends with a mantle half-space (vs4, r4). Moho = d1+d2+d3.
# Densities not in the inversion: Birch's law from vp (fit to invA).
function parse_invB(path)
    models = []
    for ln in readlines(path)
        s = strip(ln)
        (isempty(s) || !occursin(r"^ *[0-9]", s)) && continue
        v = try parse.(Float64, split(s)) catch; continue; end
        length(v) == 12 || continue
        t1, vs1, r1, t2, vs2, r2, t3, vs3, r3, vs4, r4 = v[1:11]
        d1, d2, d3 = t1, t1 + t2, t1 + t2 + t3
        birch(vp) = (0.0262 + 0.4389 * vp) * 1000.0    # kg/m3
        layers = Tuple{Float64,Float64,Float64,Float64,Float64}[]
        for (dtop, dbot, vs, r) in ((0.0, d1, vs1, r1), (d1, d2, vs2, r2),
                                    (d2, d3, vs3, r3), (d3, 120.0, vs4, r4))
            dbot <= dtop && continue
            vp = vs * r
            push!(layers, (dtop, min(dbot, 120.0), vp * 1000.0,
                          vs * 1000.0, birch(vp)))
        end
        push!(models, (moho = d3, layers = layers))
    end
    models
end

# ---- quantise + emit ----------------------------------------------------------
function emit(models; tag, fh)
    # 8 members evenly spaced in Moho rank (includes min and max)
    order = sortperm([m.moho for m in models])
    idx = [order[max(1, round(Int, k * (length(order) - 1) / 7))] for k in 0:7]
    println(fh, "# ---- $(tag) members (rank 1,8 of $(length(order)) by Moho) ----")
    ens_models = "ENSEMBLE_MODELS_$tag = Dict{String,Vector{NTuple{5,Float64}}}(\n"
    ens_moho = "ENSEMBLE_MOHO_$tag = Dict{String,Float64}(\n"
    ens_bnds = "ENSEMBLE_BOUNDS_$tag = Dict{String,Vector{Float64}}(\n"
    stats = Float64[]
    for (k, i) in enumerate(idx)
        m = models[i]
        push!(stats, m.moho)
        layers = Tuple{Float64,Float64,Float64,Float64,Float64}[]
        bnds = Float64[]
        moho_q = q02(m.moho)
        for L in m.layers
            dtop = q02(max(L[1], 0.0))
            dbot = q02(min(L[2], 120.0))
            dbot <= dtop && continue
            push!(layers, (dtop, dbot, L[3], L[4], L[5]))
            dtop > 0 && dtop < moho_q && push!(bnds, dtop)
        end
        name = "$(tag)_r$(k)"
        println(fh, "# $name: Moho $(round(m.moho, digits=3)) km -> $(moho_q) km, " *
                    "crust interfaces $(round.(bnds, digits=1)) km")
        ens_models *= "    \"$name\" => [\n"
        for L in layers
            ens_models *= "        ($(join(L, ", "))),\n"
        end
        ens_models *= "    ],\n"
        ens_moho *= "    \"$name\" => $moho_q,\n"
        ens_bnds *= "    \"$name\" => $bnds,\n"
    end
    ens_models *= ")\n"; ens_moho *= ")\n"; ens_bnds *= ")\n"
    write(fh, ens_models); write(fh, ens_moho); write(fh, ens_bnds)
    @printf("%s: Moho min/med/max = %.1f / %.1f / %.1f km (full ensemble)\n",
            tag, minimum([m.moho for m in models]),
            ([sort([m.moho for m in models])[div(end,2)]]...)[1],
            maximum([m.moho for m in models]))
    return [models[i] for i in idx]
end

if abspath(PROGRAM_FILE) == @__FILE__   # script mode: regenerate models_gen.jl
A = parse_invA(INV_A)
B = parse_invB(INV_B)
@printf("parsed: %d invA models, %d invB models\n", length(A), length(B))

open(joinpath(@__DIR__, "..", "config", "models_gen.jl"), "w") do fh
    println(fh, "# ===========================================================================")
    println(fh, "# GENERATED by tools/ensemble_pick.jl — DO NOT HAND-EDIT")
    println(fh, "# Source: InSight joint-Ps-RF crust ensembles, Knapmeyer-Endrun et al.")
    println(fh, "# 2021, Science 373 (doi:10.18715/IPGP.2021.kqwc4yaf); 8 members per")
    println(fh, "# ensemble, evenly spaced in Moho-depth rank of 5000. Layers:")
    println(fh, "# (dtop_km, dbot_km, vp_m_s, vs_m_s, rho_kg_m3), 200 m quantised,")
    println(fh, "# clipped at 120 km. invB densities via Birch fit to invA.")
    println(fh, "# ===========================================================================")
    emit(A; tag="invA", fh=fh)
    emit(B; tag="invB", fh=fh)
    write(fh, """
const ENSEMBLE_MODELS = merge(ENSEMBLE_MODELS_invA, ENSEMBLE_MODELS_invB)
const ENSEMBLE_MOHO   = merge(ENSEMBLE_MOHO_invA, ENSEMBLE_MOHO_invB)
const ENSEMBLE_BOUNDS = merge(ENSEMBLE_BOUNDS_invA, ENSEMBLE_BOUNDS_invB)
""")
end
println("wrote config/models_gen.jl")
end                                     # script mode
