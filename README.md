# Mars FDTD — synthetic seismograms with a hand-written elastic solver

Forward-modeling experiments from the "Mars open data" track: teleseismic
plane P waves sent through the Martian crust at the InSight landing site,
computed with a private 2-D elastic FDTD solver (Fdtd2DElastic.jl,
spatial order 4, free surface, CPML, GSLS attenuation) and validated
against 1-D ray theory, then run across the InSight receiver-function
crust ensembles.

## The idea

InSight's SEIS recorded marsquakes and used P-to-S receiver functions to
constrain the Martian crust. Before touching the real waveforms, this repo
answers the sanity question: does the suite reproduce the kinematics of a
plane P wave in Martian crust models to the precision the InSight
receiver-function literature works at? And a second question that turned
out to be the interesting one: at InSight frequencies (0.5 Hz), which
crust models can a receiver function resolve at all?

All open data or open models:

- TAYAK 1-D reference model (Khan et al. 2016): 3-layer crust, Moho at
  77.4 km, quantised to a 200 m grid.
- invA/invB receiver-function crust ensembles (Knapmeyer-Endrun et al.
  2021, Science 373): 5000 posterior models each, 3-layer crust over a
  half-space, 8 members per ensemble picked evenly in Moho-depth rank.

## Geometry (v3)

- Domain 120 x 120 km, 601 x 601 cells, dx = dz = 200 m, PML 40 cells,
  free surface on top, dt = 14.65 ms.
- Source: line of monopoles at 100 km depth, x = 10-82 km, one per cell,
  Hann-tapered over the outer 12 km (full-amplitude plane zone 22-72 km).
  A slanted plane wave is built from per-monopole delays
  `tau(x) = (x - x_center) * p`, `p = sin(slant)/vp(100 km)`.
  Pulse: Hann-windowed sine, 0.5 Hz, 8 s (InSight VBB band).
- Array: 13 receivers at 100 m depth, 1 km spacing, default x = 53-65 km
  (off the line centre, so slant runs see the plane wave the way a real
  station sees a teleseism: with moveout).
- Qs 500 crust / ~118 mantle, Qp 10000 (TAYAK QKP). The RF ensembles are
  elastic, so Q is assigned by depth relative to each model's own Moho.

## Experiment 1 — vertical control (TAYAK, slant 0, tmax 50 s)

Array at the line centre. At exactly vertical incidence P-to-S conversion
vanishes, so this run isolates the pure P kinematics:

| phase  | predicted (s) | measured (s) | error (s) |
|--------|---------------|--------------|-----------|
| P      | 16.64         | 16.63        | -0.01     |
| PPmPM  | delay 27.18    | delay 27.31  | +0.13     |

Moho multiple timing holds to half a percent. The PPmPM *amplitude* is
not clean in this run: a numerical surface-wave train (below) grows into
the window from ~44 s (late-time envelope ratio 1.12). Timing yes,
amplitude no.

## Experiment 2 — slanted plane wave (TAYAK, 10 deg, tmax 38 s)

Slant brings the receiver-function physics up: p = 2.335e-5 s/m,
measured from the array moveout 2.342e-5 (+0.3%).

The first finding is about the phase everyone calls "PsM". The Moho
P->S conversion does not travel to the surface as S: it re-converts
S->P at every crust interface it crosses on the way up. The observable
is a *package* of PsP@b members (one per crust interface) plus PsM
itself. On TAYAK the package spans 9.06 s (PsP at the 10-km interface,
the onset) to 10.76 s (PsM). A single-phase pick against the naive PsM
prediction is therefore biased early by the package spread; the
analyzer (analyze_mars.jl) tracks the package onset and reports both.
The package is measured with a matched filter: the measured P envelope
placed at every predicted member onset, cross-correlated over +-1.5 s.
A bare P-wavelet correlation fails here because the PPmP@10 coda tail
interferes with the package onset.

Measured package delay 10.49 s — inside the predicted package
[9.06, 10.76] s, +1.43 s past the onset, i.e. consistent with the
energy centroid of the member superposition. Package amplitude ~10% of
P. Run is clean (late-time envelope ratio 0.10).

## Experiment 3 — ensemble sweep (16 models, forward-operator study)

Every ensemble member (invA Moho 14.8-69.4 km, invB 24.4-76.4 km) run
at 10 deg slant, tmax 38 s, same array. Full table in
`results/ensemble_sweep.csv`. Headline:

- All 18 runs: horizontal slowness recovered from moveout within 0.6%,
  P timing at ray-theory level, late-time QA clean everywhere (0.03-0.19).
- The Moho conversion package resolves only when its onset clears the
  8-s direct-P coda (>5.5 s). That happens exactly for the thick-crust
  members:

| model    | Moho (km) | package onset (s) | measured (s) | err (s) |
|----------|-----------|-------------------|--------------|---------|
| TAYAK    | 77.4      | 9.06              | 10.49        | +1.43   |
| invA_r8  | 69.4      | 5.94              | 6.21         | +0.26   |
| invB_r8  | 76.4      | 5.63              | 5.63         | +0.00   |

- The other 13 ensemble members (Moho 14.8-52.8 km) have package onsets
  of 0.9-4.5 s: inside the direct-P coda, honestly flagged unresolved.
  At 0.5 Hz with an 8-s wavelet, a receiver function simply cannot
  separate a thin Martian crust — the conversion package rides the P
  coda. This is the forward-operator statement behind why the InSight
  RF inversions carry thick-crust solutions: the data cannot rule
  thin crust in, only out where delays exceed the pulse span.
- Two marginal members (invA_r7, invB_r7, Moho 48.6/52.8 km) sit at
  the coda edge; their matched-filter hits the search boundary, so they
  are reported unresolved rather than half-picked.

invB gotcha, documented in `tools/ensemble_pick.jl`: the ensemble file
stores layer THICKNESSES, VS, and VP/VS ratios (not depths and VP);
densities are not in the inversion and come from a Birch fit
rho = 0.0262 + 0.4389*vp. The parser emits SI units
(vp/vs m/s, rho kg/m3), 200 m quantised, Moho = d1+d2+d3.

## The numerical free-surface instability

Verified with Float32 frame dumps and region energy probes at 10 and
15 deg: slant runs carry a numerically growing surface-wave train
(phase speed ~2.5 km/s, e-folding ~1 s) that marches outward along the
free surface from the source-line centre, seeded around t = 25-35 s,
amplifying ~40x from line centre to the domain corners, and detonating
the free-surface/side-PML corners by ~50-60 s. Hann-tapering the source
line reduces the seeding but does not stop it. The bottom PML absorbs
cleanly.

Practical rules that fell out of it:

- Slant experiments are trustworthy only to ~35-38 s; all windows used
  here (P, conversion package, PPmP@10) close by then.
- It exists even at vertical incidence (the tmax-50 control picks it
  up from ~44 s, asymmetrically, left of the array first).
- `late_env_ratio` in every summary.csv (late envelope over direct-P
  envelope) is the QA metric: <0.5 clean, ~1 contaminated. All 16
  slant runs sit at 0.03-0.19.
- Discriminator if in doubt: the train sweeps the array at ~0.4 s/km
  apparent slowness, 17x steeper than the plane-wave phases, so a
  delay-vs-x profile separates them immediately.

## Layout

```
config/env_mars.jl       model-selectable slice env, grid, ElasticGrid
config/sr_mars.jl        plane-wave line source, taper, receivers, BCs
config/models.jl         TAYAK layers + merge of the ensemble members
config/models_gen.jl     GENERATED ensemble members (tools/ensemble_pick.jl)
tools/ensemble_pick.jl   invA/invB parser: 5000-model ensembles -> 8+8 members
run_mars.jl              driver (MARS_MODEL / MARS_SLANT / MARS_TMAX env vars)
analyze_mars.jl          slant ray predictions, matched-filter picks, PNG
results/ensemble_sweep.csv   one summary row per run, all 18 runs
runs/<model>/            summary.csv, picks.csv, meta.jl per run
```

Fdtd2DElastic.jl itself lives in a private repo and is not vendored
here; drop it next to `run_mars.jl` and the configs run unmodified.
Section PNGs regenerate with `analyze_mars.jl` (they are binary, so
they are not committed).

## Reproduce

```
julia -t 2 run_mars.jl                      # TAYAK, 10 deg slant, 38 s (~33 s)
MARS_RUN=runs/tayak_slant10 julia analyze_mars.jl

MARS_MODEL=invA_r8 MARS_RUN=runs/invA_r8 MARS_TMAX=38 julia -t 2 run_mars.jl
MARS_RUN=runs/invA_r8 julia analyze_mars.jl

julia tools/ensemble_pick.jl                # regenerate models_gen.jl
```

## Next steps

1. Real data: S1222a (M 5, the largest marsquake, 2190 km from InSight,
   sol 1256) from the NASA PDS GeoCSV archive, plus a hand-written
   miniSEED decoder for the IPGP FDSN service. The MQS catalog v14 and
   featured-event GeoCSV/mSEED pairs give free decoder validation.
2. Shorter pulse / higher f0 to sharpen the package onset and pull the
   thin-crust members out of the P coda (bounded by grid resolution:
   min_ppw=6 at 1 Hz on this grid).
3. Invert: the sweep is the forward operator; feed measured delays and
   package shapes back into the ensemble posterior.

## References

- Khan, A. et al. (2016), TAYAK spherical Martian reference model,
  PEPI 258. doi:10.1016/j.pepi.2016.07.001
- Knapmeyer-Endrun, B. et al. (2021), Thickness and structure of the
  Martian crust from InSight seismic data, Science 373.
  doi:10.18715/IPGP.2021.kqwc4yaf (ensemble data)
- Birch, F. (1961), the velocity-density relation for crustal rocks.
- InSight Marsquake Service catalog v14 (2023),
  https://www.seis-insight.eu/en/public/marsquake-catalog
- InSight SEIS data: NASA PDS Geosciences Node,
  https://pds-geosciences.wustl.edu/insight/urn-nasa-pds-insight_seis/
