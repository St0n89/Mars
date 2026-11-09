# ===========================================================================
# RUN: Mars FDTD — slanted teleseismic plane P wave through a selectable
# 1-D Mars model (TAYAK or InSight RF-ensemble crusts), surface array.
#
# Usage:
#   MARS_MODEL=invA_r3 MARS_SLANT=10 MARS_RUN=runs/invA_r3 \
#   MARS_TMAX=38 MARS_FRAMES="" julia -t 2 run_mars.jl
# Env vars:
#   MARS_MODEL    model name from config/models.jl   (default TAYAK)
#   MARS_SLANT    plane-wave incidence from vertical, deg (default 10)
#   MARS_ARRAY_X0 first receiver x, m                (default 53000)
#   MARS_TMAX     simulation length, s                (default 60)
#   MARS_RUN      output dir  (default mars_run)
#   MARS_FRAMES   frame dir, "" disables  (default mars_frames)
#   MARS_FRAME_EVERY  steps between frames           (default 100)
# Output: <MARS_RUN>/traces.csv (t + receiver traces) and meta.jl.
# ===========================================================================

const RUN_DIR    = get(ENV, "MARS_RUN", "mars_run")
const FRAME_DIR  = get(ENV, "MARS_FRAMES", "mars_frames")
const FRAME_EVERY = parse(Int, get(ENV, "MARS_FRAME_EVERY", "100"))
const tmax = parse(Float64, get(ENV, "MARS_TMAX", "60.0"))

t0 = time()
include("Fdtd2DElastic.jl")
using .Fdtd2DElastic
println("Solver loaded in $(round(time() - t0; digits=1)) s")

include("config/env_mars.jl")
include("config/sr_mars.jl")

t0 = time()
result = run_fdtd(grid, sources, bc, tmax;
    monopoles      = monopoles,
    receivers      = receivers,
    frame_every    = FRAME_EVERY,
    pml_width      = pml_cells,
    f_ref          = f0,
    check_f_max    = check_f_max,
    frame_dir      = isempty(FRAME_DIR) ? nothing : FRAME_DIR,
    progress_every = 1000,
    preflight      = :warn)

runtime = round(time() - t0; digits=1)
nt = length(result.receiver_traces[1])
println("Run finished in $(runtime) s: $(nt) steps, dt = $(result.dt) s")

mkpath(RUN_DIR)
open(joinpath(RUN_DIR, "traces.csv"), "w") do io
    hdr = ["t_s"; ["recv_$(round(Int, x[r.x]))_m" for r in receivers]]
    println(io, join(hdr, ","))
    for n in 1:nt
        t = result.dt * (n - 1)
        vals = [tr[n] for tr in result.receiver_traces]
        println(io, join(round.([t; vals]; sigdigits=8), ","))
    end
end

open(joinpath(RUN_DIR, "meta.jl"), "w") do io
    println(io, "const model = \"", MODEL, "\"")
    println(io, "const moho_km = ", MOHO_KM)
    println(io, "const slant_deg = ", SLANT_DEG)
    println(io, "const p_slant = ", P_SLANT)
    println(io, "const x_center = ", X_CENTER)
    println(io, "const dt = ", result.dt)
    println(io, "const nt = ", nt)
    println(io, "const nx = ", nx)
    println(io, "const nz = ", nz)
    println(io, "const dx = ", dx)
    println(io, "const dz = ", dz)
    println(io, "const f0 = ", f0)
end
println("Wrote $(RUN_DIR)/traces.csv + meta.jl")
