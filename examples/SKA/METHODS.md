# SKA drift-scan series — how the machinery actually works

A practitioner's tour of what the four notebooks in this directory really do:
the sky and beam models, the TOD simulation, the map-maker and the exact
recipe it is called with, the metrics, and the pitfalls we hit. It documents
*only* what the notebooks use — no polarisation, no MPI, no patchbeam — and
assumes you know the concepts but have not run this pipeline.

Companion documents: `HANDOFF.md` (current state, results, queued work),
`README.md` (experiment index), and the notebooks themselves.

---

## 1. The shape of the whole thing

Every experiment is a **closed-loop simulation**:

```
known sky ──> TODSim (beam convolution + noise) ──> TOD
                                                     │
truth ◄──────────── compare ◄── HPW_mapmaking ◄──────┘
```

We simulate time-ordered data from a sky we know perfectly, invert it back to
a map, and difference against the input. Because every ingredient is known and
seeded, any residual is *attributable*: we can switch off noise components,
swap beams, or change the scan pattern one at a time and read off the cost of
each. That attribution game — not realism — is the design principle behind
most choices below.

Fixed across all experiments: SKA-Mid-like site (Karoo, lat −30.713°), one
15 m dish, one frequency channel at **350 MHz** (HI at z = 3.06 — the extreme
low end of Band 1, hence the very large beam), 2 s sampling, 3 h of
integration total, nside 64 maps.

## 2. Sky model

`ska_common.gdsm_equatorial_sky_model` wraps pygdsm's `GlobalSkyModel`
(GSM2008): generate at 350 MHz in the model's native *galactic* frame, rotate
to equatorial with `hp.Rotator(coord=["G","C"])` at native resolution, then
`ud_grade` to nside 64.

The rotation is the whole point of the wrapper. limTOD's own
`GDSM_sky_model` returns the galactic-frame map unrotated (still true in
v1.8.0), while the simulator and map-maker both point the beam in RA/Dec —
using it directly puts the galactic plane along Dec ≈ 0. This is on the
upstream bug list; until it lands, always use the wrapper.

The two fields, which differ only in LST window:

| field | patch mean | patch rms | role |
|---|---|---|---|
| b ≈ +50° (exp 001/002) | 20.6 K | 2.3 K | typical IM field — the feasibility numbers |
| galactic plane, l ≈ 43° (exp 003/004) | 161 K | 74 K | stress test — 8× brighter, 32× more structure |

Both error mechanisms in this pipeline scale with sky brightness (§4, §7), so
plane results are pessimistic by construction and the two fields must never be
mixed in one quoted number.

## 3. Beams

Two models, both in `ska_common.py`, both **azimuthally symmetric**, both
chromatic (FWHM ∝ λ/D, D = 15 m), both returned as sum-normalised 1-D HEALPix
maps at nside 64 following limTOD's `beam_func(*, freq, nside)` protocol:

- **`ska_beam_func` ("gauss")** — symmetric Gaussian, FWHM = 1.22 λ/D =
  **3.99°** at 350 MHz. The idealised control.
- **`ska_airy_beam_func` ("airy")** — tapered circular aperture with edge
  illumination 0.2 (−14 dB edge taper), evaluated in closed form from Bessel
  functions. FWHM **3.84°**, first sidelobe ring at ~**−23 dB**, rings spaced
  ~λ/D ≈ 3.3°, back hemisphere zeroed. The "realistic dish".

Both simulation and map-making truncate the rotated beam at **1e-5 of peak
(−50 dB)** — deliberately loose, so the airy sidelobe rings survive. (The DSA
notebooks' 1e-3 would amputate everything past the first ring.)

Three scenarios per experiment:

1. **gauss** — Gaussian sky, Gaussian operator (matched, idealised);
2. **airy** — tapered sky, tapered operator (matched, realistic);
3. **mismatch** — tapered-aperture *sky* solved through the *Gaussian*
   operator: unmodelled sidelobes. Because both TODs share noise seeds and the
   operator/prior/mask are identical, `mismatch − gauss` maps isolate pure
   sidelobe pickup. Result: ~3 mK off-plane, ~8 K on the plane.

One convention note: limTOD's beam coordinate convention (boresight at the
map's north pole, φ measured from "up" toward "right") matters a great deal in
general and is documented at length upstream — but it is **invisible to these
beams**, which depend on θ only. That is not an accident; it removes one
entire class of error from the study.

## 4. TOD simulation (`TODSim.generate_TOD`)

Per time sample, the pointing chain is:

1. UTC (`start_time_utc` + `time_list` offsets) → LST via astropy;
2. (LST, latitude, azimuth, elevation) → Euler angles;
3. rotate the beam's alm to that pointing (`hp.rotate_alm`), synthesise on the
   sky grid, truncate at 1e-5;
4. TOD sample = beam-weighted sum over the sky map.

At ~60–70 samples/s this is the expensive step: ~2 min per 5400-sample
strategy per beam, from cold.

The data model as configured here (no receiver temperature, no additive
`Tsys_others` — the only "signal" is sky):

```
TOD(t) = (1 + g(t)) · T_sky(t) · (1 + η(t))
```

- **`g(t)` — 1/f gain drift.** `flicker_model.sim_noise` with the generator
  defaults `[f0, fc, alpha] = [1.335e-5, 1.099e-3, 2]` (angular frequencies),
  i.e. a pure f⁻² drift across the whole observable band. Crucially it is
  **fractional**: the contamination it injects is `g(t)·T_sky(t)`, so it
  scales with sky brightness. Same realization on a 10× brighter field → 10×
  the contamination. This single fact drives most of the plane-vs-off-plane
  behaviour.
- **`η(t)` — white noise**, also fractional, variance `WHITE_VAR = 2.5e-6`
  (so ~0.16% of signal per sample).

**Seeding is the load-bearing wall of the analysis.** Before each night/pass
`i`, the notebooks call `np.random.seed(SEED + i)` with `SEED = 0`;
`generate_TOD` then draws the 1/f sequence first and the white sequence
second, from the legacy global RNG. Consequences we exploit everywhere:

- gauss and airy TODs of the same strategy share *identical* noise
  realizations → scenario differences are purely the beam;
- experiments 001/002 and 003 share time grids and seeds → identical `g(t)` on
  fields of different brightness → a controlled measurement of the brightness
  scaling;
- any cached TOD can be *decomposed after the fact* by replaying the two
  draws and dividing them out (`sky = TOD / ((1+g)(1+η))`, verified to 4e-14
  relative). Every "no 1/f" / "no white" / "noiseless floor" number in the
  series is produced this way, without re-running the beam convolution.

## 5. Scan strategies

Both strategies produce exactly **5400 samples = 3 h**, so comparisons are at
fixed observing time.

**Drift (`drift_scan_night`).** Parked at azimuth 0°, elevation 52/50/48° on
three consecutive "nights", 1 h each; night *n* is offset by *n* sidereal days
so all three nights see the same LST (hence RA) window, and the 2° elevation
steps (≈ half-beam) stack three Dec strips. Off-plane starts
2024-04-15 19:00 UTC (LST ≈ 151°); the plane run only changes the start to
03:34 UTC (LST ≈ 280°).

The geometric fact that shapes everything: a parked dish is fixed in the
rotating Earth frame, so its boresight traces a **constant-Dec circle swept in
+RA** — track position angle 90° at *every* azimuth (verified numerically at
az 0/45/90/270/315). Azimuth selects *which* Dec; it cannot change the scan
direction. **A drift scan cannot cross-link.**

**Constant-elevation raster (`constant_elevation_scan`).** MeerKLASS-style:
fixed elevation 38.19° (chosen so the throw centre lands on the drift patch's
Dec ≈ +9.3°), triangle-wave azimuth ±7.5° about az 45° (rising pass, east of
meridian) and az 315° (setting pass, west), 150 s per half-sweep = 6 arcmin/s,
1.5 h per pass. The two passes cross the same sky at **~74°** — that crossing
is the one thing the drift cannot do. Cost: the raster spreads the same 3 h
over ~2.16× more pixels (684 vs 317 off-plane), so per-pixel depth is ~1.5×
shallower; depth-matched comparisons rescale its noise by 1/√2.16.

## 6. Map-making (`HPW_mapmaking`)

**Construction (expensive, cached).** The constructor replays the pointing
chain of §4 to build the linear model: for each time sample, rotate the beam
to the pointing and take its response over the map pixels. Pixels enter the
map if the stacked, peak-normalised beam response exceeds `threshold = 0.05`
anywhere in the scan — that defines `pixel_indices` (the patch: 317–684
pixels here). The result is one operator **A** per TOD, shape
(n_time × n_pix); rows sum to ≈ 1 (beam-weighted average of a
sum-normalised beam). Construction cost is the same per-sample rotation loop
(~2.5 min per strategy), so instances are pickled and reloaded; the caches are
regenerated if deleted.

**Solve (cheap).** Calling the instance solves the regularised normal
equations

```
(Aᵀ N⁻¹ A + S⁻¹ + λI) x = Aᵀ N⁻¹ d + S⁻¹ μ
```

jointly over all TODs (multi-night/multi-pass = block-stacked A and d). The
notebooks' exact recipe, identical in every experiment:

| ingredient | value | comment |
|---|---|---|
| prior mean μ | truth smoothed with the *assumed* beam FWHM, on the patch | idealised — stands in for "an external low-resolution model" |
| prior covariance S | diagonal, σ = std(truth on patch), uniform | one number; weakly informative |
| noise variance N | per-sample: `WHITE_VAR·(A@truth)² + floor` | the white-noise level the sample actually has; floor avoids zeros |
| regularization λ | 1e-12 | numerical only |
| high-pass (optional) | Butterworth order 2, cutoff 2 mHz, filtfilt | applied as a matrix to **both** d and A, so the model stays consistent |

Two honest caveats about this recipe:

- The prior mean is derived from the truth. That is fine for *attribution*
  (everything is differenced against the same truth) but it flatters absolute
  performance; a real survey would use an external sky model with its own
  errors. The solve does genuinely beat its prior — prior alone sits 26.6 K
  from the plane truth, the solve 21.4 K — so the machine is not just
  returning μ.
- The noise variance uses the truth too (`A@truth` = the noiseless TOD). In
  practice you would estimate it from the data; here it removes one more
  confounder.

**The high-pass** deserves its own note because its verdict flips with
strategy. It is a temporal filter: it removes TOD modes below 2 mHz *from
both data and operator*, so sky structure that only appears below 2 mHz in
the TOD becomes invisible to the solve and is surrendered to the prior. A
drift scan crosses the sky once per hour — 99.8% of its sky power sits below
2 mHz — so the filter strips signal and *worsens* the plane residual by 36%.
The raster modulates sky at the 3.3 mHz sweep rate, keeps 47% of sky power
above the cutoff, and gains −10 to −13%. Same filter, opposite sign, purely
from scan geometry. (Also: a temporal Butterworth on a raster is our
consistency choice, not field practice — MeerKLASS-type pipelines use
PCA/SVD cleaning — so raster HP numbers are indicative only.)

**Resolution.** nside 64 (0.92° pixels) = **4.4 pixels across the FWHM**,
comfortably above the ~3 px/FWHM sampling rule of thumb, and it is the *only*
correct choice in this codebase: `HPW_mapmaking` hardcodes
`normalize_beam=False`, so any `nside_target ≠ nside_beam` silently scales
every operator row by (nside_target/nside_beam)² — measured row sums 1/16 and
1/4 at nside 16/32 — and produces garbage maps. Even with a natively built
coarse beam, nside 32 undersamples (2.2 px/FWHM, 3× worse floor) and nside 16
cannot represent the beam at all (37% of beam amplitude beyond the band
limit). Keep 64/64.

## 7. Metrics — what the quoted numbers mean

**Residual = `np.std(est − truth)` on the patch.** `std`, not `rms`: the
monopole (a real +0.5–3.7 K offset) is removed, which is both small (0.2–2%
effect) and appropriate, since IM discards the monopole anyway. Always stated
with two qualifiers:

- **Scale.** *Per-pixel* residual on a beam-oversampled grid is dominated
  (75–82%) by sub-beam structure the data never constrained — legitimate but
  pessimistic. *Beam-scale* averages the residual onto nside 16 (~3.7°) cells
  that lie entirely inside the patch, and is the error at scales the survey
  measures. Off-plane: 1.30 K per pixel vs **0.31 K** beam-scale. The choice
  changes strategy comparisons by 4× (cross-linking: −10% per pixel, −38%
  beam-scale), so it is never optional to say which one you mean.
- **Reference.** All headline numbers are against *unsmoothed* truth. Against
  beam-smoothed truth the plane floor drops 21.4 → 17.7 K; we quote the harder
  one.

**Noise budget by ablation.** Using the RNG replay (§4), re-solve the same
operator with the same data minus one component: `full`, `no 1/f`, `no
white`, `clean`. Differences in quadrature give per-component map errors, and
they close to <0.1% — plane: 21.4 K floor ⊕ 7.3 K white ⊕ 4.8 K 1/f = 23.1 K
total. The "clean" solve *is* the **beam + prior floor**: the residual with
noiseless data, i.e. sky the strategy never measured, filled by the prior.

**Mode counting.** The grid-independent summary of a strategy's information
content: eigenvalues of AᵀN⁻¹A (no prior). The drift constrains **~17**
modes above 1% of λ_max (3.9 effective); the crossed raster **~34** (5.5
effective) — identical at nside 32 and 64, over patches of 300–700 pixels.
Everything outside those modes is prior. This is also why the naive
radiometer estimate Tsys/√(samples-per-pixel) is the wrong benchmark here: it
assumes every pixel is an independent measurement, overstating the
information by ~20×.

## 8. Pitfalls log

Things that actually bit us, in decreasing order of danger:

1. **`nside_target ≠ nside_beam` is silently wrong** (row-sum bug, §6). No
   error, no warning, plausible-looking garbage. On the upstream report list.
2. **`GDSM_sky_model` frame** (§2). Use the wrapper.
3. **A "beam-matched coarse grid" is a trap.** The 21.5× pixels-per-beam
   ratio at nside 64 looks like overparameterisation, but coarsening trades a
   prior-filled null space (harmless, prior's job) for beam undersampling
   (real information loss). Diagnose with mode counts, not pixel counts.
4. **Per-pixel residuals on an oversampled grid mislead** in both directions:
   they overstate absolute error (~4×) and understate strategy differences
   (~4×). Quote beam-scale alongside.
5. **VSCode notebooks:** close a notebook before executing it externally
   (`jupyter nbconvert --execute --inplace`), or the editor's cached copy
   overwrites the result.
6. **The editable install means the checked-out branch selects the limTOD
   version.** Old branches silently revert the library.
7. **Caches** (`simulated_TODs_*.npz`, `mapmaker_ops_*.pkl`) are keyed by
   filename only — nothing hashes the config. Change geometry or seeds
   without renaming the cache and you will analyse stale data.

## 9. Where the numbers live

| what | where |
|---|---|
| shared helpers (sky, beams, scans, geometry) | `ska_common.py` |
| off-plane drift, beams, mismatch (exp 001/002) | `ska_drift_scan.ipynb` |
| plane drift + 1/f budget (exp 003) | `ska_drift_galplane.ipynb` |
| plane raster vs drift (exp 004) | `ska_meerklass_scan.ipynb` |
| off-plane raster, temperature scale, summary figures | `ska_results_summary.ipynb` |
| current state, results tables, queued work | `HANDOFF.md` |
| archived per-experiment PDFs | `results/` |
