# Mars FDTD — synthetic seismograms with a hand-written elastic solver

First forward-modeling experiment from the "Mars open data" track: a
teleseismic-style plane P wave sent through the Martian crust at the InSight
landing site, computed with a private 2-D elastic FDTD solver
(Fdtd2DElastic.jl, spatial order 4, free surface, CPML, GSLS attenuation)
and validated against 1-D ray theory.

## The idea

InSight's SEIS recorded marsquakes and used P-to-S receiver functions to
constrain the Martian crust. Before touching the real waveforms, this repo
answers the sanity question: does the suite reproduce the kinematics of a
plane P wave in the TAYAK reference model (Khan et al. 2016) to the
precision the InSight receiver-function literature works at?

Setup, all open data or open models:

- Model: TAYAK 1-D Martian reference model, quantised to a 200 m grid.
  Crust: 3 layers (0-1, 1-10, 10-77.4 km), Moho at 77.4 km, sharp
  interfaces, Qs 500 crust / 118 mantle, Qp 10000 (TAYAK QKP).
- Source: line of monopoles at 100 km depth, x = 10-30 km, approximating a
  vertically incident plane P wave. Pulse: Hann-windowed sine, 2 cycles at
  0.5 Hz (InSight VBB band).
- "Array": 13 receivers at 100 m depth, x = 14-26 km, free surface on top.
- Grid: 201 x 601 cells (40 x 120 km), dt = 14.65 ms, 4096 steps, 60 s.

## Results

Arrival onsets measured by cross-correlation against the direct-P wavelet
(all phases share the source pulse, so the estimator is bias-free):

| phase    | predicted (s) | measured (s) | error (s) | note                          |
|----------|---------------|--------------|-----------|-------------------------------|
| P        | 16.64         | 16.63        | -0.01     | direct P from 100 km          |
| Ps77     | 27.34         | 27.11        | -0.24     | Ps conversion at the Moho     |
| PPmP77   | 43.82         | 43.70        | -0.12     | Moho multiple (27.2 s delay)  |
| Ps10     | 18.33         | -            | -         | inside direct-P coda          |
| PPmP10   | 20.76         | -            | -         | inside direct-P coda          |

- Direct P through 5 layers: within one time sample of the ray prediction.
- Moho Ps conversion and PPmP multiple: sub-1% timing error over delays of
  10.7 and 27.2 s. The kinematics hold; the suite is usable for Mars RF
  forward modeling.
- Physics check that came out right by accident: at exactly vertical
  incidence, P-to-S conversion vanishes, and indeed the Ps amplitudes are
  ~3% of P, near the numerical coda. Real receiver functions rely on
  oblique incidence (~15-30 deg), so the next experiment slants the plane
  wave with per-source delays.

`analyze_mars.jl` (pure-stdlib Julia: ray predictions, cross-correlation
picks, and a hand-rolled PNG writer with CRC32/Adler32 and stored DEFLATE
blocks, no zlib) regenerates `picks.csv` and `section.png` (record section
with predicted arrivals overlaid) from `mars_run/traces.csv`.

## Layout

```
config/env_mars.jl   TAYAK slice: layers, grid, ElasticGrid
config/sr_mars.jl    plane-wave line source, pulse, receivers, BCs
run_mars.jl          driver (writes mars_run/traces.csv + meta.jl)
analyze_mars.jl      ray predictions, cross-correlation picks, PNG
results/             picks.csv + meta.jl from the validated run
```

Fdtd2DElastic.jl itself lives in a private repo and is not vendored here;
drop it next to `run_mars.jl` and the configs run unmodified.

## Reproduce

```
julia -t 2 run_mars.jl     # ~25 s on 2 cores, 121k cells, 4096 steps
julia analyze_mars.jl
```

## Next steps

1. Slanted plane wave (per-source delays, 20 deg incidence) to bring Ps
   conversions up to realistic amplitudes.
2. Swap the crust for the InSight receiver-function ensembles (invA/invB,
   Knapmeyer-Endrun et al. 2021, doi:10.18715/IPGP.2021.kqwc4yaf) and watch
   Ps77/PPmP77 move; that is the forward kernel of a crustal inversion.
3. Real data: S1222a (M 5, the largest marsquake, 2190 km from InSight,
   sol 1256) from the NASA PDS GeoCSV archive, plus a hand-written
   miniSEED decoder for the IPGP FDSN service. The MQS catalog v14 and
   featured-event GeoCSV/mSEED pairs give free decoder validation.

## References

- Khan, A. et al. (2016), TAYAK spherical Martian reference model,
  PEPI 258. doi:10.1016/j.pepi.2016.07.001
- Knapmeyer-Endrun, B. et al. (2021), Thickness and structure of the
  Martian crust from InSight seismic data, Science 373.
- InSight Marsquake Service catalog v14 (2023),
  https://www.seis-insight.eu/en/public/marsquake-catalog
- InSight SEIS data: NASA PDS Geosciences Node,
  https://pds-geosciences.wustl.edu/insight/urn-nasa-pds-insight_seis/
