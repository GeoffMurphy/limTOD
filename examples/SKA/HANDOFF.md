# SKA drift-scan feasibility — handoff (2026-07-30)

Goal of the series: **evaluate whether drift scanning is viable for SKA-Mid-like
single-dish intensity mapping.** Everything here is 350 MHz, 15 m dish, Karoo
site, GDSM sky, nside 64 maps.

---

## TOP PRIORITY TOMORROW

**Re-run everything at beam-matched resolution (nside 16).** This is not a
nice-to-have — it is the thing standing between these results and a real
feasibility number.

The map is solved at nside 64 (0.92° pixels) under a 3.99° beam, so the
inversion is **21.5× overparameterised** (= beam solid angle 18.04 deg² / pixel
area 0.839 deg²). It is effectively deconvolving, which inflates *both* error
terms by ~2 orders of magnitude:

| | pixels | samples/px | radiometer limit | measured white term | inflation |
|---|---|---|---|---|---|
| drift | 321 | 16.8 | 81 mK | 7277 mK | **90×** |
| scan | 681 | 7.9 | 117 mK | 15977 mK | **136×** |

At nside 16 the patch holds ~20 pixels with ~269 samples each, a 20 mK
radiometer limit, and only **1.3×** overparameterisation — a well-posed solve.

Two things to watch:

1. The overparameterisation factor is identical (21.5×) for both strategies —
   it depends only on nside and beam — but the *consequence* is not: the scan
   is penalised ~1.5× harder (136× vs 90×). Beam-matching should therefore
   favour the scan, if anything.
2. **The relative conclusions may not survive.** Cross-linking helps largely by
   conditioning an ill-posed inversion. Make the inversion well-posed and its
   benefit may shrink. This is untested and is the main risk to the headline
   result below.

Absolute temperatures from any run so far must **not** be quoted as achievable
map error. Quote percentages. This caveat has already gone out on Slack.

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

## Queued after the nside 16 re-run

1. **Depth-matched raster done properly** — narrow the azimuth throw so the scan
   natively selects ~321 px, instead of the analytic noise scaling used above.
2. **More crossing angles** — currently only two (~75° apart). More elevations
   or azimuths should push the floor below 15.7 K.
3. Report to the limTOD author: `GDSM_sky_model` still returns the unrotated
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
