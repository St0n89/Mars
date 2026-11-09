# ===========================================================================
# SOURCE + RECEIVERS: slanted teleseismic plane P wave (RF geometry)
#
# A line of monopoles at 100 km depth, one per grid cell, each delayed by
#   tau(x) = (x - x_center) * p,   p = sin(SLANT_DEG)/vp(100 km)
# which is exactly a plane wavefront with horizontal slowness p. The
# line is Hann-tapered over its outer 12 km (full-amplitude plane zone
# x = 22-72 km): the taper suppresses edge diffraction off the line
# ends, which otherwise feeds the top corners (free surface meets side
# PML) and detonates a grazing-wave corner instability by ~45-50 s
# (verified at 10 and 15 deg, untapered; the bottom PML stays clean).
#
# IMPORTANT, verified on TAYAK: even tapered, a numerically growing
# surface-wave train (phase speed ~2.5 km/s, e-folding ~1 s, marching
# outward along the free surface from the line centre) is present in
# ALL slant runs from ~t=30-35 s; it multiplies by ~40x from the line
# centre to the domain corners. Windows that end before ~35 s (direct
# P, the conversion package, PPmP off the 10-km interface) are clean.
# Run slant experiments with MARS_TMAX <= 50 and check late_env_ratio
# in analyze_mars.jl output; PPmPM (onset ~43.4 s) is only trustworthy
# when the array sits at the line centre (the wave passes it early,
# while still tiny) — i.e. in vertical-control geometry, not at slant
# where the array must sit in the shifted PPmPM footprint.
# ===========================================================================

const SLANT_DEG = parse(Float64, get(ENV, "MARS_SLANT", "10"))
const SRC_DEPTH = 100000.0

# source pulse: 2-cycle Hann-windowed sine at f0; the derivative is
# injected as the pressure rate (4-cycle pulses smear the P coda into
# the conversion window; 2 cycles were validated in the first run)
const f0 = 0.5
const T_pulse = 4.0 / f0
const src_amp = 1.0e3
const check_f_max = 1.0

# horizontal slowness (s/m) from the P velocity at source depth
vp_src = 7413.64
for L in layers
    if L[1] <= SRC_DEPTH / 1000.0 < L[2]
        global vp_src = L[3]
        break
    end
end
const P_SLANT = sind(SLANT_DEG) / vp_src      # s/m
const X_CENTER = 46000.0                      # line centre, m

# Hann-windowed sine pressure pulse and its derivative (injection rate)
function hann_sin(t)
    (t <= 0.0 || t >= T_pulse) && return 0.0
    sin(2.0 * pi * f0 * t) * 0.5 * (1.0 - cos(2.0 * pi * t / T_pulse))
end
function hann_sin_rate(t)
    (t <= 0.0 || t >= T_pulse) && return 0.0
    w  = 0.5 * (1.0 - cos(2.0 * pi * t / T_pulse))
    dw = pi / T_pulse * sin(2.0 * pi * t / T_pulse)
    src_amp * (2.0 * pi * f0 * cos(2.0 * pi * f0 * t) * w +
               sin(2.0 * pi * f0 * t) * dw)
end

# source line at z = 100 km (iz = 501), x = 10-82 km, one monopole per
# cell. Conjugate source points for the default array [53, 65] km at
# 10 deg (drift 100 km -> surface: P 14.3, PsM 9.7, PPmPM 28.6 km):
# P [38.7, 50.7], PsM [43.3, 55.3], PPmPM [24.4, 36.4] km — all inside
# the full-amplitude plane zone [22, 72] for every ensemble model
# (thin invA crusts shift PPmPM conjugates left to ~23 km, still full
# amplitude). At SLANT=0 use MARS_ARRAY_X0=40000 so the array sits at
# the line centre (see the instability note above).
const TAPER_M = 12000.0
taper_w(xs) = let d = min(xs - 10000.0, 82000.0 - xs)
    d >= TAPER_M ? 1.0 : 0.5 * (1.0 + cos(pi * (TAPER_M - d) / TAPER_M))
end
iz_src = Int(round(SRC_DEPTH / dz)) + 1
sources = Source[]
monopoles = MonopoleSource[]
for ix in 51:411
    xs = (ix - 1) * dx
    dly = (xs - X_CENTER) * P_SLANT   # plane-wave moveout delay
    w = taper_w(xs)
    w > 0.0 && push!(monopoles, MonopoleSource(ix, iz_src, t -> w * hann_sin_rate(t - dly)))
end

bc = BC(:pml, :pml, :free, :pml)     # free surface on top = Mars surface

# receivers at z = 100 m (iz = 2), 13 of them, 1 km spacing from
# MARS_ARRAY_X0 (default 53 km: conjugates inside the tapered line's
# full-amplitude zone for all phases and ensemble models)
const ARR_X0 = parse(Float64, get(ENV, "MARS_ARRAY_X0", "53000"))
receivers = Receiver[]
for xr in ARR_X0:1000.0:(ARR_X0 + 12000.0)
    push!(receivers, Receiver(Int(round(xr / dx)) + 1, 2))
end

@printf("Source: %d tapered monopoles at z=%.0f km (x=10-82 km, centre %.0f km), slant=%.1f deg (p=%.3e s/m), f0=%.2f Hz\n",
        length(monopoles), SRC_DEPTH / 1000, X_CENTER / 1000, SLANT_DEG, P_SLANT, f0)
@printf("Array:  %d receivers, x=%.0f-%.0f m, z=100 m\n",
        length(receivers), ARR_X0, ARR_X0 + 12000.0)
