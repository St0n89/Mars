# ===========================================================================
# SOURCE + RECEIVERS: teleseismic-style plane P wave from the mantle
#
# A line of monopole pressure sources at 100 km depth spanning x = 10-30 km
# approximates a vertically incident plane P wave, like the P onset of a
# distant marsquake (e.g. S1222a) arriving under InSight. The upgoing P
# converts to Ps at the Moho (77.4 km), the 10 km and 1 km interfaces,
# and produces PPmP Moho multiples — the receiver-function physics that
# InSight used to constrain the Martian crust (Knapmeyer-Endrun et al.
# 2021, Science).
#
# Pulse: Hann-windowed sine, 2 cycles at f0 = 0.5 Hz (VBB band). Short
# enough that the P coda does not smear into the Ps(Moho) window.
# The injected rate is the time derivative of the pressure pulse.
# ===========================================================================

const f0 = 0.5                      # Hz, dominant frequency
const check_f_max = 1.0             # Hz, upper band edge for the preflight
src_amp = 1.0e3                     # Pa, pressure pulse amplitude (arbitrary)
T_pulse = 2.0 / f0                  # 4 s source duration

# Hann-windowed sine pressure pulse and its derivative (injection rate)
hann_sin(t) = begin
    if t < 0.0 || t > T_pulse
        0.0
    else
        sin(2.0 * pi * f0 * t) * 0.5 * (1.0 - cos(2.0 * pi * t / T_pulse))
    end
end
hann_sin_rate(t) = begin
    if t <= 0.0 || t >= T_pulse
        0.0
    else
        w  = 0.5 * (1.0 - cos(2.0 * pi * t / T_pulse))
        dw = pi / T_pulse * sin(2.0 * pi * t / T_pulse)
        src_amp * (2.0 * pi * f0 * cos(2.0 * pi * f0 * t) * w +
                   sin(2.0 * pi * f0 * t) * dw)
    end
end

# Source line at z = 100 km (iz = 501), x = 10-30 km (ix = 51-151)
const src_depth = 100000.0
iz_src = Int(round(src_depth / dz)) + 1
sources = Source[]
monopoles = MonopoleSource[]
for ix in 51:151
    push!(monopoles, MonopoleSource(ix, iz_src, hann_sin_rate))
end

bc = BC(:pml, :pml, :free, :pml)     # free surface on top = Mars surface

# Receivers at z = 100 m (iz = 2), x = 14-26 km, 1 km spacing
receivers = Receiver[]
for xr in 14000.0:1000.0:26000.0
    push!(receivers, Receiver(Int(round(xr / dx)) + 1, 2))
end

println("Source/receivers loaded: $(length(monopoles)) monopoles at 100 km, " *
        "f0=$(f0) Hz, $(length(receivers)) surface receivers (x=14-26 km)")
