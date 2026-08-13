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

Branch `ska-drift-analysis`, **pushed and in sync with origin** (2026-08-13).
Everything through experiment 005 is committed: the sweep code and notebook,
the write-ups here and in `METHODS.md`, and `ANALYSIS_LOG.md` (the chronological
record of reasoning and maths — results belong here, mechanisms in `METHODS.md`,
derivations there).

**Deliberately not in git**, so do not be surprised by a thin checkout:

- `mapmaker_ops_*.pkl` (25–100 MB each) and `simulated_TODs_*.npz` — gitignored,
  rebuilt by `run_freq_sweep.sh` or the notebooks if deleted. Four early SKA
  caches are tracked from before that rule and were left alone rather than
  rewritten out of history.
- `figures/` is gitignored repo-wide; the freqsweep PNGs regenerate from
  `plot_freqsweep_maps.py`.
- `ANALYSIS_LOG.pdf` is built from `ANALYSIS_LOG.md`.

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
| `ska_freq_sweep.ipynb` | exp 005 — off-plane drift swept across Band 1 (5 channels, shared nside 128); config + cached heavy steps in `ska_freq_sweep.py`, driver `run_freq_sweep.sh` |
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

## Off-plane raster (added 2026-07-31)

Experiment 004's raster existed only on the plane. `ska_results_summary.ipynb`
now simulates an **off-plane counterpart** — same site, frequency, Gaussian
beam, sky, noise model, seeds, 3 h and recipe as the off-plane drift; only the
pass start times shift, by the −8.574 sidereal hours between the two fields'
RA. Caches: `simulated_TODs_ska_meerklass_offplane_gauss.npz` and
`mapmaker_ops_ska_meerklass_offplane_gauss_ns64.pkl`. Gauss only — the
tapered-aperture beam matters for the sidelobe question, which is a bright-field
result.

Geometry check: both passes centre on RA ≈ 157.9° against the drift patch's
158.30°; track position angles 126.6° and 52.7°, a **73.9° crossing** (the plane
version was ~75°). Raster selects 684 px to the drift's 317, so it is 2.16×
shallower at fixed clock time.

| off-plane, no HP, gauss | per-pixel | beam-scale |
|---|---|---|
| drift floor → raster floor | 1.287 → 1.158 K (**−10.0%**) | 0.300 → 0.187 K (**−37.8%**) |
| drift total → raster total @ matched depth | 1.300 → 1.241 K (−4.5%) | 0.313 → 0.175 K (**−44.0%**) |

**Cross-linking is worth ~40% off-plane at beam scale, but only ~10% per
pixel.** The per-pixel number is diluted by sub-beam structure *neither*
strategy measures, which is why beam-scale is the fair comparison — and why the
plane's headline −26.5% floor gain (per-pixel) understated the effect.

Note the residuals carry a positive monopole (drift +0.754 K, raster +0.542 K
off-plane). All quoted numbers are `np.std`, so it does not enter them.

### Flat-prior audit of the above (added 2026-08-13)

Every number in this file uses the beam-smoothed-truth prior, which turned out
to supply most of the apparent frequency trend in experiment 005. Re-solved
under a flat prior (constant at the patch mean) with everything else held fixed
— `audit_flat_prior.py`, a re-solve off the caches, ~1 min. The smoothed-prior
arm reproduces the published numbers to three decimals, so the comparison is
like for like.

**The cross-linking result survives, and strengthens** — the opposite of what
happened to the frequency trend, where the prior manufactured the effect. Here
it was *masking* it:

| off-plane, no HP, gauss | per-pixel | beam-scale |
|---|---|---|
| smoothed prior, floor | 1.287 → 1.158 K (−10.0%) | 0.300 → 0.187 K (−37.8%) |
| **flat prior, floor** | 1.806 → 1.218 K (**−32.6%**) | 0.392 → 0.210 K (**−46.4%**) |
| **flat prior, total @ depth** | 1.819 → 1.294 K (−28.9%) | 0.402 → 0.196 K (−51.3%) |

**The reason is that the drift leans on the prior far harder than the raster
does.** Withdrawing the prior's structure costs the drift 40% per-pixel
(1.287 → 1.806 K) and the raster only 5% (1.158 → 1.218 K). That is what the
mode-count argument predicts — the drift measures 17 modes to the raster's 34,
so it has twice the null space for a truth-derived prior to fill flatteringly.
Prior-independently the floor is 0.786 × sky structure for the drift against
0.530 for the raster.

So quote the **flat-prior** gains when the claim is about what the survey
measures. The case for cross-linking is not just that the raster's residual is
lower, it is that the drift's number was propped up by a prior no real survey
has.

**Unresolved, and pre-existing.** At beam scale the raster's floor exceeds its
total (0.187 vs 0.181 K smoothed, 0.210 vs 0.199 K flat) — adding noise slightly
*reduces* the beam-scale residual, so the floor/noise quadrature split that holds
on the plane to <0.1% does not hold off-plane at beam scale. Present in the
published numbers too. Understand this before quoting a beam-scale noise term
for this field.

---

## Frequency sweep across Band 1 (experiment 005, added 2026-08-05)

Everything above is 350 MHz, the bottom edge of Band 1. `ska_freq_sweep.ipynb`
repeats the **off-plane** drift geometry of experiment 001 at five equally
spaced channels — **350, 525, 700, 875, 1050 MHz** — changing nothing but the
frequency. Gaussian beam only (sidelobes are a bright-field question). All five
solved on a **shared nside 128 grid** so residuals compare like for like;
config and the cached heavy steps live in `ska_freq_sweep.py`, warmed by
`./run_freq_sweep.sh` (~1.5 h for the sweep on 12 cores).

**Frequency is itself a lever on the mode count — the headline result.** Queued
item 0 below asks how to raise the number of measured modes. Going up in
frequency does it without touching the strategy:

| | 350 MHz | 525 | 700 | 875 | 1050 |
|---|---|---|---|---|---|
| FWHM | 3.99° | 2.66° | 2.00° | 1.60° | 1.33° |
| modes > 1% λmax | 17 | 27 | 37 | 46 | **55** |
| observed area [deg²] | 270 | 184 | 150 | 133 | 119 |
| modes / 100 deg² | 6.3 | 14.7 | 24.6 | 34.5 | **46.2** |
| residual / sky rms | 0.525 | 0.452 | 0.427 | 0.437 | 0.457 |

**3.2× more modes across the band, 7.3× per unit sky** — against cross-linking's
factor of 2 (17 → 34). The mode count is robust: the trend holds at every
eigenvalue threshold (>10%: 8 → 38; >1%: 17 → 55; >0.1%: 23 → 69) and is
unchanged if the truth-dependent noise weighting is dropped for a plain `AᵀA`
(17 → 56). The patch *shrinks* 2.3× at the same time, so the cost of going up in
frequency is **area, not map quality**.

**The residual row above is biased high at the top channels — corrected here.**
On the shared nside 128 grid the ratios (0.53, 0.45, 0.43, 0.44, 0.46) look
flat, but 875 and 1050 MHz sit at 3.5 and 2.9 px/FWHM and are undersampled.
Re-run at nside 256 they improve substantially, and the **flatness turns out to
be entirely an artefact**:

| residual ÷ sky rms | 350 | 525 | 700 | 875 | 1050 |
|---|---|---|---|---|---|
| shared nside 128 | 0.525 | 0.452 | 0.427 | 0.437 | 0.457 |
| **adequately sampled** | 0.525 | 0.452 | 0.427 | **0.350** | **0.336** |
| (grid used) | ns128 | ns128 | ns128 | ns256 | ns256 |

**The corrected trend is monotonic: the drift reconstructs the sky better, in
relative terms, the higher in Band 1 you go** — 0.525 → 0.336, a 36%
improvement. The mode count is unaffected by the regrid at either channel
(875: 46 → 46; 1050: 55 → 56), which is the fourth independent confirmation of
its grid-independence.

Caveat when quoting this: the *absolute* residual falls ~27× across the band
(1.242 → 0.034 K) but the sky structure itself falls ~24× (2.37 → 0.10 K), so
almost all of the absolute drop is the sky dimming as ν^−2.7. The reconstruction
improvement is the 1.56× in the ratio, not the 27×.

**The elevation ladder stops working above 700 MHz.** 52°/50°/48° is 2° Dec
steps — half a beam at 350 MHz, but a *whole* beam at 700 MHz and 1.5 beams at
1050 MHz, so the three strips separate into ribbons with unobserved gaps. This
does not degrade the residual (the mode gain outweighs it) but it is why the
observed area falls faster than the beam does. A ladder spaced in *beams* is the
obvious fix.

**Off-plane the 2 mHz high-pass is a non-issue at every channel** (effect within
±5%, and it does *not* track the crossover). Two separate reasons, and the first
corrects a claim made earlier in this file.

*The sky band does not move with frequency.* The tempting prediction — a beam
crossing takes FWHM / (15 deg/hr), so the sky band should climb from 1.05 to
3.14 mHz and cross the fixed 2 mHz cutoff near 700 MHz — is **wrong at the
premise**, and measuring the cached sky TODs refutes it:

| | 350 | 525 | 700 | 875 | 1050 MHz |
|---|---|---|---|---|---|
| beam edge 1/(crossing time) | 1.05 | 1.57 | 2.09 | 2.62 | 3.14 mHz |
| **median-power frequency** | **0.28** | **0.28** | **0.28** | **0.28** | **0.28 mHz** |
| 90% of power below | 1.67 | 1.67 | 1.67 | 1.67 | 1.67 mHz |
| power below 2 mHz | 93.0% | 92.7% | 92.9% | 93.2% | 93.6% |

The beam edge triples across the band; **nothing in the measured spectrum moves
with it.** Diffuse GDSM emission is red, so the large-scale gradient along the
drift dominates however narrow the beam gets — the sky power sits at the pass
fundamental (0.28 mHz, one 1 h traverse) at every channel. This is consistent
with the original empirical note above ("99.8% of the drift's sky power sits
below 2 mHz ... fundamental 0.28 mHz"); the beam crossing time bounds the
*support* of the sky signal, not where its power lives, and a high-pass cares
about the latter. `beam_crossing_freq_hz` is documented accordingly — do not
compare it against a filter cutoff.

*And the term it would govern is not the limiting one anyway.* The residual is
96–100% beam + prior floor and only 5–8% 1/f, so perfectly removing 1/f could
not move the total by more than a few percent. The cutoff question belongs to
the plane, where 1/f is a much larger share of a much larger residual. Note
also that the **+36% high-pass penalty quoted above is the galactic-plane
drift** — off-plane at nside 64 the same filter already helped slightly
(−9.7%).

**Cross-check, and a third grid for the mode-count claim.** The 350 MHz channel
is experiment 001 re-run at nside 128 (4.06× the pixels):

| 350 MHz, same field/seeds | nside 64 | nside 128 |
|---|---|---|
| pixels solved | 317 | 1288 |
| per-pixel residual | 1.300 K | 1.242 K |
| beam-scale residual | 0.313 K | 0.178 K |
| **modes > 1%** | **17** | **17** |
| effective modes | 3.7 | 3.7 |

The mode count is identical, as the grid-independence argument requires — it now
rests on nside 32, 64 and 128 (and 256, below). The beam-scale numbers differ
only because
`beam_cells` keeps cells whose pixels *all* lie inside the patch, and that
criterion tightens 4× with the grid; beam-scale values are comparable within the
sweep but **not** against the nside-64 numbers quoted earlier in this file.

### Prior dependence: the residual trend was the prior, not the survey

Raised by a reviewer after the maps figure was circulated, and they were right.
The Wiener prior mean is the truth smoothed with **that channel's beam**, so as
the beam narrows the prior sharpens by itself — the prior is not held fixed
across the sweep. Measured (residual ÷ sky rms):

| | 350 | 525 | 700 | 875 | 1050 |
|---|---|---|---|---|---|
| beam-smoothed-truth prior, **no data** | 0.511 | 0.470 | 0.438 | 0.406 | 0.360 |
| full reconstruction, that prior | 0.525 | 0.452 | 0.427 | 0.350 | 0.336 |
| what the data bought | 0.97× | 1.04× | 1.03× | 1.16× | 1.07× |
| **flat prior**, reconstruction | 0.769 | 0.685 | 0.650 | 0.597 | 0.597 |
| what the data bought | 1.30× | 1.46× | 1.54× | 1.67× | 1.67× |

**The 0.525 → 0.336 trend is the prior sharpening, not the map improving** —
the data contributes 0.97–1.16×, and at 350 MHz it makes the map marginally
*worse* by injecting noise into the modes it measures. With a flat prior
(constant at the patch mean) the data's own contribution is visible and does
grow with frequency: residual 0.769 → 0.597, gain 1.30× → 1.67×, saturating
above 875 MHz.

**So the conclusion survives but the evidence changed.** Quote the flat-prior
numbers or the mode count — both prior-independent — when the claim is about
what the survey measures. `plot_freqsweep_maps.py` emits both map versions plus
`figures/freqsweep_prior_dependence.png`, which is the figure to send anyone
who asks about priors.

Note this is not specific to experiment 005: every residual in this file is
quoted with the truth-derived prior. What is specific is that a *frequency*
trend is confounded, because the prior moves with the beam. On the plane the
data does beat its prior comfortably (1.24×) — off-plane it barely does.

### Grid check at 1050 MHz: ~3 px/FWHM is NOT a safe floor

The top channel re-run at nside 256 (5.8 px/FWHM) against the shared nside 128
grid (2.9 px/FWHM). Cost ~1.2 h; run it detached (`setsid nohup`) — background
shell jobs do not survive a session ending, which killed two earlier attempts.

| 1050 MHz | ns128 (2.9 px/FWHM) | ns256 (5.8) | |
|---|---|---|---|
| modes > 1% | 55 | 56 | +1.8% |
| effective modes | 18.95 | 18.99 | +0.2% |
| observed area | 119.0 | 117.0 deg² | −1.7% |
| sky structure rms | 0.0985 | 0.1000 K | +1.5% |
| **per-pixel residual** | **0.0451** | **0.0336 K** | **−25.5%** |
| beam + prior floor | 0.0432 | 0.0328 K | −24.1% |

**Two different answers in one table.** The *mode count* is grid-independent —
that claim now rests on nside 32, 64, 128 and 256. The *residual* is not: it
falls ~25% on the finer grid while the sky structure it is measured against
barely moves (+1.5%), so this is a genuine reconstruction improvement, not a
normalisation artefact. Note it moves *opposite* to the naive expectation, since
a finer grid has more sub-beam structure to get wrong.

So **`HANDOFF.md`'s ~3 px/FWHM rule of thumb is too generous** — 2.9 px/FWHM
cost 25% in residual here. Use ~5 px/FWHM when the residual (rather than the
mode count) is the quantity of interest. The nside-64 conclusions elsewhere in
this file are unaffected: they sit at 4.4 px/FWHM and their headline claims are
mode-count based.

875 MHz was re-run the same way and behaves identically: modes 46 → 46, sky rms
+1.5%, residual 0.0723 → 0.0587 K (−19%). Both re-runs are cached
(`*_f0875_ns256.pkl`, `*_f1050_ns256.pkl`); 700 MHz at 4.4 px/FWHM is the
remaining unchecked channel, and by the trend of the other two is likely biased
high by <10%.

---

## Queued next

0. **Raise the number of measured modes** — this is now the headline lever,
   replacing the nside re-run. The drift measures 17, the scan 34; everything
   else is prior. More crossing angles and more elevations both add modes,
   whereas more integration time on the same tracks does not. **Experiment 005
   adds a third lever: frequency** — 55 modes at 1050 MHz on the same tracks,
   at the price of a 2.3× smaller patch.
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
