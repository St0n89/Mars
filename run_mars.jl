# ===========================================================================
# RUN: first Mars FDTD — TAYAK slice + plane P wave + surface array
#
# Usage:  julia -t 2 run_mars.jl
# Output: mars_run/traces.csv (t + 13 receiver traces),
#         mars_run/meta.jl   (dt, grid, f0),
#         mars_frames/       (Float32 snapshots for the movie).
# ===========================================================================

const tmax = 60.0
const frame_every = 100

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
    frame_every    = frame_every,
    pml_width      = pml_cells,
    f_ref          = f0,
    check_f_max    = check_f_max,
    frame_dir      = "mars_frames",
    progress_every = 500,
    preflight      = :warn)

runtime = round(time() - t0; digits=1)
nt = length(result.receiver_traces[1])
println("Run finished in $(runtime) s: $(nt) steps, dt = $(result.dt) s")

mkpath("mars_run")
open("mars_run/traces.csv", "w") do io
    hdr = ["t_s"; ["recv_$(round(Int, x[r.x]))_m" for r in receivers]]
    println(io, join(hdr, ","))
    for n in 1:nt
        t = result.dt * (n - 1)
        vals = [tr[n] for tr in result.receiver_traces]
        println(io, join(round.([t; vals]; sigdigits=8), ","))
    end
end

open("mars_run/meta.jl", "w") do io
    println(io, "const dt = ", result.dt)
    println(io, "const nt = ", nt)
    println(io, "const nx = ", nx)
    println(io, "const nz = ", nz)
    println(io, "const dx = ", dx)
    println(io, "const dz = ", dz)
    println(io, "const f0 = ", f0)
end
println("Wrote mars_run/traces.csv + meta.jl")
