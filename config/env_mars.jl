# ===========================================================================
# ENVIRONMENT: Mars — TAYAK reference model slice (crust + upper mantle)
#
# 2-D column through the Martian crust at the InSight landing site:
# TAYAK 1-D model (Khan et al. 2016, PEPI 258; Smrekar et al. 2019),
# quantised to the 200 m grid. z = 0 is the FREE SURFACE (Mars surface),
# z increases downward.
#
# Layers (depth-km, vp-m/s, vs-m/s, rho-kg/m3, Qs):
#    0.0 -  1.0    3677.71  1739.80  1865.33   500   (megaregolith)
#    1.0 - 10.0    4952.25  2780.97  2271.82   500
#   10.0 - 77.4    5846.66  3281.16  2681.72   500   (lower crust)
#   77.4 - 80.0    7400.90  4246.33  3388.58   118   (mantle top = Moho)
#   80.0 - 100.0   7413.64  4247.35  3393.07   118
#  100.0 - 110.0   7436.12  4180.19  3399.95   118
#  110.0 - 120.0   7446.65  4181.17  3403.31   118
# Qp = 10000 everywhere (TAYAK QKP). Interfaces are SHARP (no smoothing):
# the PmP/Ps conversions at the Moho and the two crustal interfaces are
# the physics under test.
# ===========================================================================

dx = 200.0
dz = 200.0
xmax = 40000.0            # 40 km lateral (plane wave in the centre)
zmax = 120000.0           # 120 km depth (Moho at 77.4 km, source at 100 km)
pml_cells = 40            # 8 km absorbing rim (l/r/bottom)

nx = Int(round(xmax / dx)) + 1        # 201
nz = Int(round(zmax / dz)) + 1        # 601
x = collect(0.0:dx:xmax)
z = collect(0.0:dz:zmax)

# (d_top_km, d_bot_km, vp, vs, rho, Qs) — grid-quantised TAYAK
const MARS_LAYERS = [
    (0.0,    1.0,   3677.71, 1739.80, 1865.33, 500.0),
    (1.0,    10.0,  4952.25, 2780.97, 2271.82, 500.0),
    (10.0,   77.4,  5846.66, 3281.16, 2681.72, 500.0),
    (77.4,   80.0,  7400.90, 4246.33, 3388.58, 118.2),
    (80.0,   100.0, 7413.64, 4247.35, 3393.07, 118.2),
    (100.0,  110.0, 7436.12, 4180.19, 3399.95, 118.4),
    (110.0,  120.0, 7446.65, 4181.17, 3403.31, 118.4),
]

lam = zeros(nx, nz)
mu  = zeros(nx, nz)
rho = zeros(nx, nz)
Qp  = fill(10000.0, nx, nz)          # TAYAK QKP
Qs  = zeros(nx, nz)

for j in 1:nz
    zk = z[j] / 1000.0
    for L in MARS_LAYERS
        if L[1] <= zk < L[2] || (j == nz && zk == L[2])
            vp, vs, r, qs = L[3], L[4], L[5], L[6]
            @inbounds for i in 1:nx
                rho[i,j] = r
                mu[i,j]  = r * vs^2
                lam[i,j] = r * (vp^2 - 2.0 * vs^2)
                Qs[i,j]  = qs
            end
            break
        end
    end
end

grid = ElasticGrid(nx, nz, dx, dz, lam, mu, rho, Qp, Qs)
println("Environment loaded: Mars TAYAK slice, $(nx) x $(nz) cells, dx=dz=$(dx) m, " *
        "Moho at 77.4 km, free surface on top")
