# SKA drift-scan feasibility — handoff (2026-07-30)

Goal of the series: **evaluate whether drift scanning is viable for SKA-Mid-like
single-dish intensity mapping.** Everything here is 350 MHz, 15 m dish, Karoo
site, GDSM sky, nside 64 maps.

---

## RESOLVED 2026-07-31: the nside 16 re-run was the wrong call — keep nside 64

The 2026-07-30 handoff made "re-run at beam-matched nside 16" the top priority,
on the grounds that nside 64 under a 3.99° beam is 21.5× overparameterised.
**That recommendation was wrong and has been withdrawn.** nside 64 is the
correct grid. Nothing in the results below needs retracting.

Three independent lines all point the same way:

1. **Sampling rule of thumb (~3 pixels across the FWHM).** Pixel size vs the
   3.99° beam: nside 16 → 1.1 px/FWHM, nside 32 → 2.2, nside 64 → **4.4**.
   Coarsening *undersamples* the beam and throws away real information. The
   earlier argument conflated sampling with conditioning; degenerate sub-beam
   modes are the Wiener prior's job, and it was already doing it.
2. **Band-limiting.** A 3.99° Gaussian retains 37% amplitude at nside 16's band
   edge (lmax 47) but only 1.8% at nside 32's (lmax 95). nside 16 cannot
   represent this beam at all.
3. **Measured residuals.** Coarser is monotonically worse — the drift's
   noiseless floor goes 21.4 K (ns64) → 70.0 K (ns32, native beam) →
   84.7 K (ns16). At nside 32 the high-pass verdict even flips sign, which is
   how you can tell the coarse grids are unreliable rather than merely blunt.

### Why the residual is large — the real answer

Not overparameterisation. Eigen-analysis of `AᵀN⁻¹A` (no prior) gives, **and
this is identical at nside 32 and nside 64**:

| | effective modes (Σλ/λmax) | modes > 1% of λmax | pixels solved |
|---|---|---|---|
| drift | 3.9 | **17** | 321 |
| scan | 5.5 | **34** | 681 |

The strategy — not the pixelisation — sets how much sky is measurable. The
drift constrains ~17 independent modes over that patch however finely you
pixelate; the rest of the map is prior. So "Tsys/√(samples per pixel)" was
never the right benchmark: it assumes every pixel is independently measured and
overstates the information by ~19×. Coarsening to nside 16 would still have
measured 17 modes — the grid was never the problem.

This also gives a clean, grid-independent statement of why cross-linking helps:
**it doubles the number of measured modes** (17 → 34).

### How to quote residuals

Absolute temperatures are meaningful, but state the metric. Per-pixel RMS
against unsmoothed truth is legitimate but pessimistic — 75–82% of it lives in
sub-beam modes the data never constrained. Averaged onto beam-sized cells
(nside 16 cells fully inside the patch):

| no HP, gauss | per-pixel | beam-scale |
|---|---|---|
| drift total | 23.1 K | **5.7 K** |
| drift floor | 21.4 K | 5.1 K |
| scan total | 26.8 K | **4.4 K** |
| scan floor | 17.0 K | 3.7 K |

The Slack caveat that absolute temperatures are meaningless was too strong, and
the "deconvolution issue we'll fix tomorrow" framing was wrong — coarsening the
grid does not fix it because there was nothing to fix. The honest correction is
the mode-count argument above.

---

## Where things stand

Branch `ska-drift-analysis`, **3 commits ahead of origin, unpushed**:

```
d605f46 SKA exp 004: add total-residual comparison plot
b444d80 SKA: add MeerKLASS-style constant-elevation scan (experiment 004)
ce10a29 Relabel "reconstruction floor" as "beam + prior floor"
```

Pushing needs no special scope now — `.github/workflows/publish.yml` was
dropped in `a4dd646` (upstream's PyPI workflow, would only ever fail on a fork,
and GitHub demands `workflow` OAuth scope to push it). It will come back as a
conflict on the next upstream merge; resolve with
`git rm .github/workflows/publish.yml`.

The repo was merged up to **upstream limTOD 1.8.0** (`972ce74`). Verified before
merging that this changes nothing: upstream reproduces the SKA sky TODs
bit-identically and matches the cached TODs to 4e-14 relative.

### Notebooks

| notebook | what it is |
|---|---|
| `ska_drift_scan.ipynb` | exp 001/002 — drift at b ≈ +50°, off-plane baseline |
| `ska_drift_galplane.ipynb` | exp 003 — same drift across the galactic plane, plus the 1/f noise budget |
| `ska_meerklass_scan.ipynb` | exp 004 — constant-elevation raster, rising + setting |
| `ska_results_summary.ipynb` | small notebook that parses the others' printed output and plots 3 summary figures |

---

## Results so far (all nside 64 — relative only)

**1/f is not the limiter.** The plane residual decomposes in quadrature (to
<0.1%) into a 21.4 K beam + prior floor, 7.3 K white, 4.8 K 1/f.

**The 2 mHz high-pass backfires on the drift (+36%) and helps the raster
(−10.6%).** Reason, measured: 99.8% of the drift's sky power sits below 2 mHz
(one traverse per 1 h pass, fundamental 0.28 mHz), so the filter removes the
signal, not the drift. The raster modulates sky at the 3.33 mHz sweep rate
(300 s per sweep) where only 53% sits below the cutoff. **Do not quote HP
numbers for the raster** — a temporal Butterworth on a raster is not standard
practice (MeerKLASS-type pipelines use PCA/SVD cleaning plus fast scanning);
2 mHz is not optimised for it.

**A drift scan cannot cross-link, at any azimuth.** A parked dish is fixed in
the rotating Earth frame, so its boresight traces a constant-Dec circle swept in
+RA: track position angle is 90° at az = 0/45/90/270/315. Azimuth only picks
*which* Dec. Verified numerically.

**Cross-linking is the lever (headline result).** Rising + setting raster passes
cross at ~75°. On the 321 shared pixels the beam + prior floor drops
**21.4 → 15.7 K (−26.5%)**.

**Total residual, matched for depth.** As run the raster looks slightly worse
(26.8 vs 23.1 K) purely because it covers 681 px to the drift's 321 in the same
3 h — 2.12× shallower. Correcting for that (noise scaled by 1/√2.12, exact for
white, approximate for 1/f):

| comparison, no HP both sides | drift | scan @ depth | |
|---|---|---|---|
| total residual | 23081 | 21857 | −5.3% |
| beam + prior floor | 21368 | 15702 | −26.5% |

Only 5.3% because at 3 h the scan's budget is half floor (15.7 K), half noise
(15.2 K) — a 1.03:1 ratio, so neither term alone dominates.

**How to improve it: more time, not more speed.**

| total integration (split evenly across the 2 passes) | total | vs drift |
|---|---|---|
| 3 h (1.5 + 1.5, as run) | 21857 | −5.3% |
| 6 h | 19030 | −17.6% |
| 12 h | 17445 | −24.4% |
| ∞ | 15702 | −32.0% |

Faster scanning attacks the wrong term: 1/f is only 27% of the white-noise
term, so eliminating it entirely buys **1.6%**. Scan speed is already solved —
6 arcmin/s puts sky at 3.33 mHz where 1/f power is ~140× below the drift's
fundamental.

When extending integration, **repeat the same 1.5 h window on more sidereal
days** rather than lengthening each pass — a longer pass drifts further in RA,
enlarges the patch and reintroduces the depth dilution.

---

## Queued next

0. **Raise the number of measured modes** — this is now the headline lever,
   replacing the nside re-run. The drift measures 17, the scan 34; everything
   else is prior. More crossing angles and more elevations both add modes,
   whereas more integration time on the same tracks does not.
1. **Depth-matched raster done properly** — narrow the azimuth throw so the scan
   natively selects ~321 px, instead of the analytic noise scaling used above.
2. **More crossing angles** — currently only two (~75° apart). More elevations
   or azimuths should push the floor below 15.7 K. Track the mode count as the
   figure of merit, not just the residual.
3. **Report to the limTOD author (new, 2026-07-31): `HPW_mapmaking` produces a
   silently wrong operator whenever `nside_target != nside_beam`.**
   `HPW_filter.py:449` (and the `num_tods == 1` branch at :470) hardcodes
   `normalize_beam=False`, and `generate_sky2sys_projection` then synthesises
   the beam alm — built at `nside_beam` — straight onto the `nside_target`
   grid. Sum-normalisation is an nside-dependent convention, so every operator
   row is scaled by `(nside_target/nside_beam)²`. Measured row sums with an
   nside-64 beam: **0.0631 at nside_target 16, 0.2525 at 32, 1.0096 at 64** —
   i.e. exactly 1/16, 1/4, 1. The map-maker then scales the sky up to
   compensate, giving nonsense (drift floor 1921 K at nside 16 vs 21 K at 64).
   Multiplying the operator by `(nside_beam/nside_target)²` restores every grid
   to the same 1.0096, confirming the diagnosis.
   Note `TODSim.generate_TOD` also defaults to `normalize_beam=False`, so
   simulation and map-making agree *only* when the nsides match — which is why
   our own 64/64 runs are unaffected. Suggested fix: normalise per pointing, or
   scale by the nside ratio, or refuse `nside_target != nside_beam`.
   Also worth telling the author that a *natively* built coarse beam is still
   not a workaround: at nside 32 (2.2 px/FWHM) per-pointing beam shape error
   inflates the floor 3× and flips the high-pass verdict, and a scalar gain
   correction does not remove it (0.90–1.02 sweep moves the floor <2%).
4. Report to the limTOD author: `GDSM_sky_model` still returns the unrotated
   galactic-frame map in 1.8.0 (so `ska_common.gdsm_equatorial_sky_model`
   remains necessary; note upstream also switched to GSM16 while ours uses
   GSM08, so adopting a fixed upstream version would change the sky model and
   break cache reproduction). Also `limTOD.patchbeam`'s `beam_solid_angle`
   integrates `dl·dm` where the measure is `dl·dm/√(1−l²−m²)` — ~0.01% for a
   main-lobe-dominated beam, but it does *not* cancel against the HEALPix
   map-makers, which have no patchbeam awareness.

---

## Practical notes

- Run everything with `/home/geoff/limTOD/.venv/bin/python` from the repo root.
  The install is editable, so **the checked-out branch decides the limTOD
  version** — checking out a pre-merge branch silently reverts to 1.1.0.
- Notebooks: `.venv/bin/python -m jupyter nbconvert --to notebook --execute
  --inplace <nb>`. **Close the notebook in VSCode first** or it will overwrite
  the executed file with its cached copy.
- Caches (`simulated_TODs_*.npz`, `mapmaker_ops_*.pkl`) are regenerated if
  deleted. Cold cost is ~11 min per beam for TOD generation plus operator
  construction; with them warm a full re-run is a few minutes.
- `CLAUDE.md`, `examples/**/figures/` and `examples/SKA/results/` are gitignored
  (some files under the latter two are already tracked, so `git add` on them
  errors unless forced).
- The noise-budget decomposition works by replaying the seeded RNG draws
  `generate_TOD` makes — 1/f first via `sim_noise`, then white via
  `np.random.normal`, seeded `SEED + index`. This reproduces the exact
  realizations in the cached TODs, verified to 4e-14. It is how every "clean" /
  "no 1/f" / "no white" number is produced.
- Residuals are quoted as `np.std(est − truth)`, i.e. RMS about the residual's
  own mean. There is a real monopole offset of 1.8–3.7 K; including it changes
  the numbers by only 0.2–2.1%, and IM discards the monopole anyway.
