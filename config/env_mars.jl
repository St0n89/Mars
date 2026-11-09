# ===========================================================================
# ENVIRONMENT: Mars crustal slice, model selectable via MARS_MODEL.
# Default "TAYAK" (Khan et al. 2016); ensemble crusts from the InSight
# receiver-function inversion (Knapmeyer-Endrun et al. 2021, Science 373)
# are selected with MARS_MODEL=invA_r1 ... invB_r8 (see config/models.jl).
#
# Domain is 120 km wide (nx=601), 120 km deep (nz=601), dx=dz=200 m.
# A slanted plane wave (MARS_SLANT) drifts laterally as it rises, and the
# Moho multiple of a 10 deg-incidence wave returns ~29 km to the side of
# its surface bounce. The wide domain keeps direct P, the conversion
# package and PPmP co-recorded at the array (see config/sr_mars.jl for
# the footprint budget and the PML-corner instability notes).
#
# Attenuation: Qs = 500 above the Moho, TAYAK mantle values below
# (118.2 to 100 km, 118.4 deeper), Qp = 10000. The RF ensembles are
# elastic, so Q is assigned by depth relative to each model's Moho.
# ===========================================================================

using Printf

const MODEL = get(ENV, "MARS_MODEL", "TAYAK")
include(joinpath(@__DIR__, "models.jl"))

dx = 200.0
dz = 200.0
xmax = 120000.0
zmax = 120000.0
pml_cells = 40

nx = Int(round(xmax / dx)) + 1        # 601
nz = Int(round(zmax / dz)) + 1        # 601
x = collect(0.0:dx:xmax)
z = collect(0.0:dz:zmax)

layers = MARS_MODELS[MODEL]
MOHO_KM = MARS_MOHO[MODEL]
moho_m = MOHO_KM * 1000.0

lam = zeros(nx, nz)
mu  = zeros(nx, nz)
rho = zeros(nx, nz)
Qp  = fill(10000.0, nx, nz)
Qs  = zeros(nx, nz)

# fill column by column: layer rows are (dtop_km, dbot_km, vp, vs, rho)
# with velocities in m/s and density in kg/m^3 (200 m quantised). Layer
# boundaries sit ON grid lines: cell iz spans [(iz-1)dz, iz*dz), so a
# layer fills cells round(z0/dz)+1 .. round(z1/dz).
for L in layers
    z0 = max(L[1] * 1000.0, 0.0)
    z1 = min(L[2] * 1000.0, zmax)
    z1 <= z0 && continue
    vp, vs, r = L[3], L[4], L[5]
    lam_l = r * (vp^2 - 2.0 * vs^2)
    mu_l  = r * vs^2
    qs_l  = z1 <= moho_m ? 500.0 : (z1 <= 100000.0 ? 118.2 : 118.4)
    iz0 = Int(round(z0 / dz)) + 1
    iz1 = Int(round(z1 / dz))
    z1 >= zmax - 1e-6 && (iz1 = nz)   # extend deepest layer to the bottom row
    for iz in iz0:min(iz1, nz)
        for ix in 1:nx
            lam[ix, iz] = lam_l
            mu[ix, iz]  = mu_l
            rho[ix, iz] = r
            Qs[ix, iz]  = qs_l
        end
    end
end
@assert minimum(rho) > 0 "unfilled cells in model $(MODEL)"

grid = ElasticGrid(nx, nz, dx, dz, lam, mu, rho, Qp, Qs)

@printf("Env: %s — Moho %.1f km, %d layers, %d x %d cells, dx=dz=%.0f m, free surface on top\n",
        MODEL, MOHO_KM, length(layers), nx, nz, dx)
