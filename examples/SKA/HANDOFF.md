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

## HI signal recovery (experiment 006, added 2026-08-18)

The first experiment in the series containing a 21 cm signal. 350-400 MHz, 32
contiguous channels, nside 64, both strategies, on the 277-px patch common to
every channel. HI from fastbox (Tb 0.518 mK, b_HI 1.559 at z_eff 2.80),
projected onto HEALPix by 3D lightcone interpolation. Code in
`ska_hi_mock.py`, `ska_hi_experiment.py`, `ska_hi_analysis.py`,
`run_hi_experiment.py`; results in
`results/hi_experiment_f350_400_nc32_ns64.npz`. Full reasoning in
`ANALYSIS_LOG.md`.

**Cross-linking improves HI recovery at every scale.** T(k_par), fraction of
map-made HI surviving the clean, 20 mocks:

| k_par | 0.012 | 0.024 | 0.036 | 0.048 | plateau |
|---|---|---|---|---|---|
| drift, 4 modes | 0.228 | 0.469 | 0.732 | 0.854 | ~0.90 |
| raster, 4 modes | 0.449 | 0.822 | 0.965 | 0.961 | ~0.96 |
| drift, 8 modes | 0.095 | 0.138 | 0.256 | 0.451 | ~0.80 |
| raster, 8 modes | 0.034 | 0.142 | 0.437 | 0.769 | ~0.94 |

T(k) is depth-independent (each arm's ratio against its own injected response),
so this is fair despite the raster being 2.3x shallower on the shared patch.

**The drift map-maker makes HI look like foreground.** Splitting the HI by
frequency structure, the drift keeps only 38% of the frequency-*varying* part
while *amplifying* the frequency-coherent part 3.14x — it projects HI onto a
nearly frequency-constant subspace, which PCA then removes. Raster: 0.532 and
1.11x. A plain rms ratio (0.924) hides this completely; do not quote it.

**No HI is recoverable here, and it is the floor, not the noise.** Ablation
(common-resolution, interior pixels, re-run 2026-08-26): noiseless "floor" and
full "total" track each other to a few per cent in amplitude at every mode
count, while noise alone sits **52x** (drift) / **10x** (raster) below the floor
at 4 modes. The blocker is entirely the beam + prior floor, whose frequency
structure comes from the chromatic null space — not spectrally smooth, so PCA
cannot touch it.

| 4 modes, common resolution | floor | total | noise alone |
|---|---|---|---|
| drift | 335.0x | 347.3x | 6.4x |
| raster | 22.7x | 24.0x | 2.2x |

The older "agree to 0.2%, noise ~38x below" was the un-reconvolved full-patch
run, and the 0.2% compared the two *median* summary numbers rather than the
per-k residuals (which differed by ~3% even then). Quote the table above.

**Residual/HI, revised 2026-08-19** after adding the common-resolution
reconvolution the pipeline had been missing (see below). Amplitude, median over
k, at 4 modes removed:

| | as run, 277 px | as run, 119 interior px | **common resolution, 119 px** |
|---|---|---|---|
| drift | 570.6x | 322.1x | 347.3x |
| raster | 65.9x | 48.5x | **24.0x** |
| cross-linking gain | 8.7x | 6.6x | **14.4x** |

Two corrections are folded in here and should not be confused. Restricting to
pixels >3 deg from the patch edge (needed because the reconvolution smooths a
zero-padded patch) removes an edge inflation that was present in the original
full-patch numbers. The reconvolution itself then **halves the raster** and does
essentially nothing for the drift, which has too few measured modes for a
cleaner beam to help. At 10 modes the raster reaches **13.4x**. Net: the
cross-linking case is considerably stronger than first published.

**Do not read `p_corrected_*` from the results file as a recovered spectrum.**
With the residual 13-350x the signal the cross-power estimator has no signal in
it. `residual_over_hi_*` is stored alongside to make that explicit. The transfer
function itself is unaffected — injection is differential and mock-averaged.

### Prior variance sweep (added 2026-09-21)

**PLAN item 2 is narrower than PLAN describes, for experiment 006 only.** The
flat-prior audit it asks for is already done here by construction: 006's prior
*mean* is flat, a constant at the patch mean, set in `run_hi_experiment.py` and
`ska_hi_ablation.py` (both commented `# flat prior`). The beam-smoothed-truth
mean that `audit_flat_prior.py` exists to test belongs to experiments 004 and
005, and the prior-dependence section above is about those. The only
`hp.smoothing` in the 006 chain is `common_resolution`, which reconvolves
*output* cubes and is not a prior.

What is left in 006 is the prior *variance*. The solve uses
`S^-1 = I / var(fg_truth_patch)` — one scalar on every pixel, no correlation, no
angular or spectral structure — plus the truth-weighted `N^-1`. That scalar is
also the crossover in `tr(WA) = sum lambda/(lambda + S^-1)`, so it moves the
measured mode count and the floor *together*. That matters more for the design
figure than for any single residual, because the design figure puts the mode
count on an axis.

`ska_hi_priorvar.py` sweeps it via `A.PRIOR_VAR_SCALE` (new, default 1.0, so
nothing published changes). Results in
`results/hi_priorvar_f350_400_nc32_ns64.npz`, figure `figures/hi_priorvar.png`.
**Scale 1 reproduces the published 335.0/347.3 and 22.7/24.0 exactly**, so this
is the published pipeline with one scalar changed.

| prior var / truth var | drift modes/100 deg² | drift floor | drift total | raster modes/100 deg² | raster floor | raster total |
|---|---|---|---|---|---|---|
| 0.01 | 9.17 | 262.7 | 262.5 | 10.45 | 8.5 | **8.6** |
| 0.1 | 11.45 | 271.6 | 275.5 | 15.00 | 15.0 | 15.3 |
| **1 (published)** | **14.62** | **335.0** | **347.3** | **22.31** | **22.7** | **24.0** |
| 10 | 18.42 | 389.9 | 404.3 | 33.56 | 49.4 | 52.5 |
| 100 | 21.93 | 514.3 | 643.0 | 47.10 | 96.2 | 136.7 |
| 1000 | 27.72 | 972.1 | 2627.0 | 61.89 | 191.3 | 531.6 |

(residual/HI at 4 modes removed, common resolution, interior pixels.)

**1. The trade space is not a picture of the prior.** The prior track has
log-log slope **+1.8 (drift) / +2.2 (raster)** in the (modes, residual) plane.
Every measured survey relation runs the other way: the ladder at **−1.28**
(drift → drift12, adding dishes) and crossing the tracks at **−22.4**
(drift8 → raster). Opposite sign either way — the prior tracks run *across*
the survey relations rather than along them. This was the question the design
figure needed answered before it could be drawn honestly.

> **Superseded 2026-09-22, corrected 2026-09-28.** This paragraph originally
> compared the prior tracks against a single **−6.3** survey curve fitted
> through the two published anchors. That slope was withdrawn below ("adding
> strips" vs "crossing the tracks") as an artefact of fitting one line across
> the discontinuity between two regimes, and `collect` no longer emits the
> `slope`/`ratio` keys it needed. `fig_priorvar` kept drawing it behind a
> `d=None` default, so the reference line silently vanished from
> `hi_priorvar.png` while the panel title still claimed "the tracks cross the
> line". Line dropped and panel retitled 2026-09-28; the conclusion is
> unchanged, since the sign comparison survives against the measured slopes.

**2. But "modes per unit sky" is not a sufficient statistic on its own.** Two
numbers settle it:

- At 100× prior variance the drift reaches **21.9 modes/100 deg²**, essentially
  the raster's published 22.3 — and its residual is **643×** against the
  raster's 24×.
- At 10× the raster reaches **33.6 modes/100 deg²**, right at the 36.9 the
  design figure's two-point extrapolation says would bring the residual to 1.
  Its actual residual is **52.5×**.
  (The 36.9 target came from the retired −6.3 extrapolation, so treat it as
  illustrative, not a design number — the measured point, 33.6 modes at
  residual 52.5×, stands on its own.)

Modes bought by loosening the prior are worth nothing; modes bought by
cross-linking are worth 14×. Physically that is the right statement —
cross-linking adds rows to `A`, the prior only reweights rows already there —
but **the design figure's x axis is a common currency only at fixed prior, and
must say so**, or a reader will take the extrapolation as a design target. Add
this to PLAN's "Figure honesty" list.

**3. The loose end stops being floor-dominated.** At 1000× the drift's total
runs 2.7× its floor (2627 vs 972) and the raster's 2.8× (532 vs 191), so the
residual has crossed from beam+prior-floor-dominated to noise-dominated. That is
the trade the prior is making, and it is the first place in this series where
the floor is not the whole story.

**4. The tight end survives T(k) — checked 2026-09-21, 10 mocks.**
`ska_hi_priorvar.py --mocks N` measures the two terms that could have faked it:
`P_mapmade / P_true` (what the map-maker does to the HI) and `T(k_par)` (what
the clean does to what is left). Results in
`results/hi_priorvar_tk_f350_400_nc32_ns64.npz`.

| | drift | | | raster | | |
|---|---|---|---|---|---|---|
| prior var | HI kept | T(k) | resid ÷ surviving HI | HI kept | T(k) | resid ÷ surviving HI |
| 0.01 | 0.344 | 0.745 | 352 (0.84) | 0.365 | 0.804 | **10.7** (0.39) |
| 0.1 | 0.359 | 0.788 | 350 (0.83) | 0.389 | 0.846 | 18.1 (0.66) |
| **1 (published)** | 0.366 | 0.829 | **419** (1.00) | 0.411 | 0.873 | **27.5** (1.00) |
| 10 | 0.371 | 0.857 | 472 (1.13) | 0.431 | 0.888 | 59.2 (2.15) |
| 100 | 0.373 | 0.871 | 738 (1.76) | 0.450 | 0.901 | 152 (5.50) |
| 1000 | 0.375 | 0.908 | 2893 (6.91) | 0.459 | 0.923 | 576 (20.9) |

At 0.01× the total surviving HI (`kept × T`) falls only **16% (drift) / 18%
(raster)** while the residual falls 24% and 2.8×, so the gain is real, not a
denominator collapse. Corrected, the raster gains **2.6×** (27.5 → 10.7) and the
drift **16%** (419 → 352, saturating below 0.1×). Above scale 1 both `kept` and
`T` are flat, so all the residual growth at the loose end is genuine noise.

Control: the 10-mock `T` at scale 1 is 0.829 (drift) / 0.873 (raster) against
the published 20-mock 0.846 / 0.882 — 2.0% and 1.0% low, consistent with
halving the mock count. The same 10 mocks (seeds 1000-1009, the first half of
the published set) are reused at every scale, so the cross-scale comparison is
clean regardless.

**Three reasons this is not a headline.** (a) It does not change the
conclusion — 10.7× is still an order of magnitude above the HI. (b) It is an
analysis choice, not a survey-design lever: the prior mean is the patch mean, so
tightening is shrinkage toward a constant, trading bias for variance, and it
wins only because it suppresses the foreground residual faster than the HI.
(c) It moves the mode count the *wrong* way, which is finding 2 above. Worth a
sentence in the paper's systematics discussion; not worth a figure.


### The declination ladder, and what cross-linking actually is (added 2026-09-22)

**PLAN item 6 is answered, and the answer overturns the design figure's
premise.** Item 6 asked for a demonstration that deliberately raising the
measured mode count moves the residual. It does -- but only linearly, and the
steep slope between the published drift and raster anchors turns out not to be
a mode-count relationship at all.

**Item 6 splits in two, and only half applies to the drift.** A parked dish
cannot cross-link at any azimuth (see above), so "more crossing angles" is a
raster-only lever. The drift's *only* geometric freedom is the declination
ladder: how many strips, and how far apart. That is also the multi-dish
question, because `drift_scan_night` offsets night *n* by a whole sidereal day,
so every strip covers the same LST window -- **N strips on N nights is
arithmetically the same operator as N dishes parked at N elevations for one
pass.** The published 3-night drift already *is* a 3-dish configuration;
multi-dish buys wall-clock, not modes.

#### Tier 1 -- the mode count alone (`ska_hi_ladder.py`)

Operator only, no TOD, one channel: ~10 min for the whole sweep.
`results/hi_ladder_f350_ns64.npz`, figures `hi_ladder.png` and
`hi_ladder_geometry.png`. Control: the published 52/50/48 geometry rebuilt from
scratch gives tr(WA) 38.77 against the cached operator's 38.78.

| strips (0.5 beams) | tr(WA) | area [deg²] | modes/100 deg² |
|---|---|---|---|
| 1 | 13.44 | 182.1 | 7.38 |
| 2 | 26.12 | 222.4 | 11.74 |
| **3 (published)** | **38.77** | **265.2** | **14.62** |
| 5 | 63.76 | 359.2 | 17.75 |
| 8 | 100.05 | 499.4 | 20.03 |

**The mechanism in one line: the first strip costs 182 deg² to buy 13.4 modes;
every later strip buys ~12 modes for ~44 deg².** Ninety-two per cent of the
information for a quarter of the sky. Density therefore rises 2.72x from 1 to 8
dishes.

**Extended to N = 20 (`results/hi_ladder_long_ns64.npz`, added 2026-09-22),
and the ladder does NOT run into diminishing returns.** The marginal return per
strip is constant across the whole range -- 12.68, 12.65, 12.50, 12.01, 11.69,
11.65, 12.21 modes per added strip from N = 2 to N = 20, against 43-47 deg² of
extra sky each time. Densities: 14.62 (N=3), 17.75 (5), 21.13 (10), 21.81 (12),
22.51 (15), 23.73 (20). The density curve flattens purely because area grows
linearly alongside modes, toward 100 x 12 / 44 ~ **27.3 per 100 deg²** -- an
arithmetic asymptote, not a physical saturation. (An earlier note here said
"saturating toward ~26.6", fitted on N <= 8; the constant-marginal-return
statement is the correct one.)

**This also brackets the raster's 22.31 modes/100 deg²**: N = 12 gives 21.81
(2.2% below) and N = 15 gives 22.51 (0.9% above), so a parked drift reaches the
cross-linked raster's mode density at N ~ 13-14. `drift12` is therefore a
matched-density comparison that needs no extrapolation, which is what the
`matched` Tier 2 run measures.

Spacing, at 3 strips: 0.25 beams gives 15.89, 0.5 gives 14.62, 0.75 gives
12.74, 1.0 gives 11.41, 1.5 gives 9.34 per 100 deg². Beyond ~0.75 beams the
total mode count flattens at exactly N x a single strip (3 x 13.44 = 40.3
predicted, 41.4 measured at 1.5 beams) -- the strips have become independent
surveys. This is the same curve as the "ladder stops working above 700 MHz"
note above, with the fixed 2° step sliding rightward in beam units as the beam
narrows.

#### Tier 2 -- does it move the residual? (ilifu, jobs 925104 / 925290)

Full 32-channel HI experiment plus ablation for two new geometries, each with
its own strategy tag so the published caches are untouched:
`drift8` (8 strips, 0.5 beams, 8 dish-hours) and `drift3t` (3 strips, 0.25
beams, 3 dish-hours -- cost-matched to the published drift). Results in
`results/hi_{experiment,ablation}_{mech,costmatched}_*.npz`.

**Mechanism test** (`mech`: drift, drift8, raster), all on the same 119
interior pixels, one HI realisation, one set of noise seeds. Common
resolution, 4 modes removed:

| strategy | modes/100 deg² | **floor** | total | noise | 1/f |
|---|---|---|---|---|---|
| drift | 14.62 | **405.9** | 409.8 | 9.0 | 5.7 |
| drift8 | 20.03 | **299.3** | 306.4 | 6.5 | 2.4 |
| raster | 22.31 | **26.7** | 27.7 | 2.6 | 0.4 |

The floor is 97-99% of the total throughout, so this is entirely a floor
result. **Two regimes:**

- **adding strips** (drift → drift8, no cross-linking): floor 405.9 → 299.3
  for 1.37x the mode density. Local slope **-0.97**.
- **crossing the tracks** (drift8 → raster): floor 299.3 → 26.7 for 1.11x the
  mode density. Local slope **-22.4**.

**At essentially the same mode density, cross-linking lowers the floor 11x.**
Mode density is therefore *not* a sufficient statistic, and the -6.3 slope in
the design-figure prototype was an artefact of fitting one line across a
discontinuity between two regimes.

**This sharpens the physics rather than weakening it.** Cross-linking is not
"more modes": it changes the *structure* of the null space. Every drift track
has position angle 90°, so its unmeasured modes share an angular structure that
correlates across frequency -- which is exactly why the residual is spectrally
non-smooth and PCA cannot remove it. Eight parallel strips do not break that
correlation; two crossed passes do. For a survey designer: **you cannot buy your
way out of a parked drift with dishes.** Eight dishes at eight elevations buy
26%; two crossed passes buy a factor of 11.

**Cost-matched test** (`costmatched`: drift, drift3t, 92 interior px) is
**inconclusive and should not be quoted as a gain.** Tightening 0.5 → 0.25
beams at fixed 3 dish-hours moves the floor 330.3 → 271.8 on the
common-resolution metric (18% better) but 526.2 → 672.7 as-run (28% worse) and
299.2 → 356.9 on interior pixels (19% worse). A sign that flips with the metric
is not a result; the mode-density gain there is only 8.7%, comparable to
realisation scatter. Note `drift3t`'s patch (222 deg²) is barely larger than
the common patch (~202 deg²), so it has little margin for the reconvolution's
boundary effects -- a plausible but unverified explanation for why only the
reconvolved metric favours it. **Practical reading: the ladder's value is in
adding strips, not in repacking the ones you have.**

#### HI realisation scatter -- PLAN item 7, half-answered by accident

`build_hi_cube` sizes the fastbox box from the *union* of the strategies' pixel
sets, so a run with a different strategy list realises a different HI field at
the same grid resolution (20 Mpc/cell either way). The `mech` run's drift
therefore gives 409.8x where the published run gives 347.3x **on identical
pixels** -- an 18% shift from the HI realisation alone, with the raster showing
15% (27.6x vs 24.0x). Consequences: (a) numbers are comparable *within* a run
and not *across* runs, which is why `drift` is included in both Tier 2 runs as
the internal reference; (b) the headline residuals need an error bar, and ~15-18%
is the scale. The machinery to do this properly is a seed loop over
`HIBandConfig`.


### Frequency repeat at 675-725 MHz (experiment 007, added 2026-09-22)

The series asserted a negative result for Band 1 having tested only its bottom
edge. `f700` repeats the whole HI experiment at **675-725 MHz, nside 128**
(4.2-4.5 px/FWHM, matching the 350-400 run's 3.8-4.4 -- `check_sampling()` now
reports this before every run). Same geometry, same seeds, only the band
changes. Results in `results/hi_{experiment,ablation}_f700_f675_725_nc32_ns128.npz`;
Ilifu job 13858047, 1 h 19 min.

**The negative result survives, which is the headline.** At 4 modes removed,
common resolution, 171 interior pixels (35.9 deg²):

| | modes/100 deg² (patch) | (local) | floor | total | noise |
|---|---|---|---|---|---|
| drift | 41.87 | 57.69 | 16.4 | **16.5** | 1.3 |
| raster | 76.55 | 124.48 | 36.2 | **38.1** | 1.2 |

Both remain entirely floor-limited. HI is not recoverable anywhere in Band 1,
not just at its bottom edge.

**But the drift/raster ordering INVERTS, and the reason is not geometric.**
Decomposing `residual/HI` into its two factors (`median sqrt(p_clean_auto_4)`
and `median sqrt(p_mapmade)`):

| band | strategy | HI kept | post-clean FG residual | resid/HI |
|---|---|---|---|---|
| 350-400 | drift | 0.362 | 6.55e-2 | 409.8 |
| | raster | 0.394 | 6.43e-3 | 27.6 |
| 675-725 | drift | 0.298 | **6.25e-4** | 16.3 |
| | raster | 0.401 | **2.60e-3** | 38.0 |

1. **The drift's foreground residual improves 105x with frequency; the raster's
   only 2.5x.** At 350 MHz the drift is crippled by its null space, so a
   narrower beam helps it enormously; the raster was already well-conditioned
   and has little to gain. Instrumentally the drift overtakes the raster.
2. **The HI signal is 4.3x fainter in the upper band** (`sqrt(p_true)`
   7.42e-4 at z~2.8 vs 1.73e-4 at z~1.03). That is cosmology, and it penalises
   the *ratio* for both strategies equally.

Net: drift 105 / 4.3 = 24x better ratio; raster 2.5 / 4.3 = 1.7x worse. Both
reproduce the measurements.

**Consequence, and it is a warning about this series' headline metric:
`residual/HI` is not comparable across bands.** It conflates instrument
performance with signal amplitude. Within a band it is fine; across bands,
quote the foreground residual and the HI amplitude separately.

**Three explanations were proposed and refuted before the decomposition**, all
worth not re-proposing: (a) the raster's cross-scan sampling degrading with
frequency -- refuted, its mode density is 1.83x the drift's at 700 MHz;
(b) the evaluation region sitting in the raster's beam wings -- refuted, both
fields are co-pointed at dec 9.3° (the raster's declination needs the full
`sin d = sin(phi) sin(h) + cos(phi) cos(h) cos(A)`, not the az=0 form);
(c) patch-average vs local mode density -- refuted by `local_modes.py`, which
computes `sum_{i in region} (WA)_ii` and preserves the ordering (raster 124.48
vs drift 57.69 per 100 deg² locally, with the worse residual).

`local_modes.py` is worth keeping regardless: it is the honest local version of
the mode count, and it shows the patch average understates the evaluation
region's mode density by 1.4-2.4x in every case.


### The 500-550 MHz point, and where the crossover is (added 2026-09-23)

`f500` fills in the middle of Band 1 at **500-550 MHz, nside 128** (272
interior pixels). Ilifu job **13858845**, 2 h 53 min, exit 0, clean log.
Results in `results/hi_{experiment,ablation}_f500_f500_550_nc32_ns128.npz`.

At 4 modes removed, common resolution, relative to HI:

| band | strategy | floor | total | noise | white | gain |
|---|---|---|---|---|---|---|
| 350-400 (ns64) | drift | 399.6 | 413.5 | 8.27 | 6.22 | 5.10 |
| | raster | 25.6 | 27.1 | 2.33 | 2.21 | 0.41 |
| 500-550 (ns128) | drift | 214.3 | 213.2 | 1.74 | 1.54 | 0.64 |
| | raster | 56.2 | 76.7 | 1.17 | 1.12 | 0.16 |
| 675-725 (ns128) | drift | 16.1 | 16.3 | 1.28 | 1.27 | 0.15 |
| | raster | 36.0 | 38.0 | 1.15 | 1.15 | 0.15 |

Floor-limited in all six cases (`floor` ~ `total`), so nothing here is a noise
statement. **The negative result holds at every band tested.**

**The inversion is real, not a resolution artefact.** This was the obvious
worry, since the 350 band is nside 64 and the other two are nside 128. It is
answered by the 500/675 pair, which share a grid, share a geometry and share
seeds: across that pair alone the drift goes from 3.8x *worse* than the raster
to 2.2x *better*. The ordering flips inside the matched-grid comparison.

**Both strategies improve monotonically with frequency; the ratio does not.**
Separating the two factors (`fg` = median `sqrt(p_clean_auto_4)`):

| band | drift fg | raster fg | drift/raster | HI `sqrt(p_true)` | drift resid/HI | raster resid/HI |
|---|---|---|---|---|---|---|
| 350-400 | 6.54e-2 | 6.45e-3 | 10.1 | 7.81e-4 | 413.5 | 27.1 |
| 500-550 | 1.48e-2 | 5.90e-3 | 2.51 | 4.52e-4 | 213.2 | 76.7 |
| 675-725 | 6.25e-4 | 2.60e-3 | 0.24 | 1.73e-4 | 16.3 | 38.0 |

Every instrumental column is monotonic: the drift's foreground residual falls
105x across Band 1, the raster's 2.5x, and the HI falls 4.5x on its own. Yet
the raster's `residual/HI` goes **27 -> 77 -> 38** -- non-monotonic, out of two
monotonic quantities. That is the cleanest available demonstration of the
warning the 700 MHz run raised: **`residual/HI` is a quotient of an instrument
number and a cosmology number and must not be read as a trend across bands.**
The third point turns that from an argument into a picture.

**Where the drift overtakes the raster: ~600 MHz.** Log-interpolating the
foreground residual ratio between the two nside-128 bands (2.51 at 525 MHz,
0.24 at 700 MHz) puts the crossing at **~594 MHz**. Treat it as "somewhere in
the upper half of Band 1", not a measured number -- it is two points and the
350 MHz point does not lie on the same straight line.

**Consequence for the paper, now with a number attached.** Every
"cross-linking is the lever" claim is a *bottom-of-Band-1* statement. Above
~600 MHz the drift's null space has narrowed enough that it wins on
foreground residual, and cross-linking stops being the discriminator. The
qualifier "at 350 MHz" must be attached in `paper/ska_drift/main.tex` and in
the earlier sections of this file (still outstanding).

**`ska_hi_design.py` now draws all eight figures in one run.** `main()` had
only ever been wired to the first four -- `hi_frequency_design`,
`hi_cleaning_depth`, `hi_chromatic` and `hi_realisations` were drawn by hand
on 2026-09-22 and would silently have gone stale. `fig_cleaning_depth` also
gained the 500 MHz band as a third line style. Both fixed 2026-09-23.


### Design figures, and two gaps closed (added 2026-09-22 evening)

`ska_hi_design.py` draws all of these from cached results in seconds. They are
prototypes for the paper, not final art; every number on them is measured.

| figure | what it is for |
|---|---|
| `hi_design.png` | the headline. Two regimes at 350-400 MHz: the ladder is a shallow power law (slope -1.28 over 3 → 12 dishes), cross-linking is a **9.3x step change at 2%-matched mode density**. Plus the lever chart, now on the FLOOR rather than the mode count. |
| `hi_frequency_design.png` | why `residual/HI` must not be compared across bands: instrument, cosmology, and their misleading quotient, in three panels. Three bands as of 2026-09-23 (350, 500, 675 MHz), so the panels show curves rather than line segments. |
| `hi_cleaning_depth.png` | the referee's first question. Residual and surviving HI both fall with PCA depth, the residual faster, so the net never turns over. **There is no optimal cleaning depth: the limit is the floor.** |
| `hi_ladder.png` | the drift's only geometry lever, out to 20 dishes, with the two measured floors overlaid. |
| `hi_ladder_geometry.png` | the ladder as sky geometry, for explaining it to others. |
| `hi_chromatic.png` | **why the floor is immune to PCA** (see below). |
| `hi_priorvar.png` | prior-amplitude robustness. |
| `hi_realisations.png` | the error bar. Nine HI realisations: 13.4% scatter on the drift, 8.7% on the raster, **6.4% on the paired ratio** -- which is why this series quotes ratios between strategies, not absolute residuals. |

#### Why the floor survives foreground cleaning (`hi_chromatic.png`)

Built from `results/hi_rank_f350_400_nc32_ns64.npz`, which already existed --
no new compute. Foreground cleaning assumes the contaminant is spectrally
smooth and therefore low-rank. Three facts, all measured:

1. **Rank.** The foreground is rank 1 (one mode carries 99% of its variance).
   The floor needs 2 (raster) to 6 (drift). The HI needs 23-25 of 32 -- it is
   the one component that is genuinely not compressible.
2. **Shape.** The foreground's first eigenvector is a clean monotonic power
   law; the floor's first three are progressively structured; the HI's is
   noise-like from the start. That progression is what "chromatic" means
   concretely -- `WA` varies across the band, so the floor carries spectral
   structure the sky never had.
3. **Overlap.** Principal angles between the floor's subspace and the HI's have
   cosines **0.42-0.78**. They are nowhere near orthogonal, so projecting out
   the floor necessarily removes HI. This is the reason deeper cleaning cannot
   win, and it is the same statement `hi_cleaning_depth.png` shows empirically.

#### HI realisation scatter -- the error bar (PLAN item 7, DONE)

`run_hi_experiment.py --hi-seed N` varies the true HI seed alone, leaving the
noise seeds and pixels fixed, so the spread isolates HI realisation scatter.
Six such runs (seeds 11-16, tag `seed<N>`) plus the three existing runs give
**nine realisations on the same 119 interior pixels**. Figure
`hi_realisations.png`; results in `results/hi_experiment_seed*_*.npz`.

| | N | mean | sd | range |
|---|---|---|---|---|
| drift | 9 | 391.5 | **52.5 (13.4%)** | 329.6-511.4 |
| raster | 9 | 26.1 | **2.3 (8.7%)** | 22.7-29.9 |
| **drift / raster** | 9 | **14.9** | **1.0 (6.4%)** | -- |

**Quote ratios between strategies, not absolute residuals.** The paired ratio
is tighter than either absolute because the realisation partly cancels when
both strategies see the same HI field. The cross-linking comparison therefore
carries a ~6% error bar, not ~13%.

**The published headline is on the optimistic side.** 347.3 sits 0.85 sd below
the mean of 391.5, and the raster's 24.0 likewise below 26.1. The published
*ratio*, 14.5, is representative of the mean 14.9 -- another argument for the
ratio being the quotable quantity.

**Consistency check:** the median transfer function is identical across all
seeds (0.829 drift, 0.873 raster), because the mocks use fixed seeds 1000+j and
are unaffected by the true HI seed. So the spread is genuinely in the
realisation and is not leaking in from the transfer-function estimate.

**Scale check, and why this matters:** 13% is the same size as several effects
this series has reported. It is why the cost-matched ladder test (0.5 -> 0.25
beams, an 8.7% mode-density change) came out inconclusive, and it is the bar
any future sub-20% claim has to clear.


### PICK UP HERE — state at 2026-10-01

**Checked 2026-10-01, 7 h 45 min in. Two done, one still running.**

| job | state | outcome |
|---|---|---|
| **13920752** density confirmation | COMPLETED 2 h 16 min | **75.72 modes/100 deg^2 at 0.579 beams, 1.08% below the raster's 76.55** -- a tighter match than drift12 managed at 350 MHz (2.2%). Pulled to `results/hi_ladder_f700d_ns128.npz`. |
| **13920776** warm array, drift12m | COMPLETED 2 h 31 min | caches built |
| **13920795** `matched700b` | **RUNNING** | see below |

**The validity check has ALREADY PASSED.** The running job's log shows the
drift arm finishing on **171 interior pixels** -- exactly the count of the
published f700 run, not the 57 that made `matched700` void. The wider
footprint did what it was meant to, so this run is usable. Its drift arm also
gives residual/HI = 16.6x at 4 modes, against the published 16.3x: consistent,
and not identical only because the HI realisation differs.

**TIMING RISK, and it is the thing to check first.** At 7 h 45 min the job had
done drift (1872 s) and was 24/32 channels into drift12m at ~1052 s/channel --
12 strips over a 13.2 deg patch is far more solve than anything run before.
Extrapolating: ~2.3 h more for drift12m, ~0.5 h for raster, then the ablation.
Total plausibly 13-16 h against a **16 h wall limit**.

`run_hi_freq.sbatch` runs the experiment and the ablation as two separate
python calls, so a timeout would leave the experiment npz written and **no
ablation** -- and the ablation is the half that separates floor from total,
i.e. the half the matched-density question needs. So:

```bash
ssh ilifu "bash -lc 'sacct -j 13920795 --format=State,Elapsed'"
ssh ilifu "bash -lc 'ls -la ~/ska/limTOD-src/examples/SKA/results/*matched700b*'"
```

If the experiment landed but the ablation did not, re-run the ablation alone
rather than the whole job -- the caches and the experiment npz are already
there:

```bash
sbatch run_ablation.sbatch matched700b drift drift12m raster   # check its env sources freq_env.sh
```

**Pull when they land:**

```bash
cd ~/limTOD/examples/SKA
for f in hi_ladder_f700d_ns128.npz \
         hi_experiment_matched700b_f675_725_nc32_ns128.npz \
         hi_ablation_matched700b_f675_725_nc32_ns128.npz; do
  ssh ilifu "bash -lc 'cat ~/ska/limTOD-src/examples/SKA/results/$f'" > results/$f
done
```

**CHECK THIS FIRST, before reading any floor number.** The previous attempt
(`matched700`, drift6m) was void because its 2.58 deg footprint became the
limiting one and the interior collapsed to 57 px over 0.60 deg — 0.29 beams.
`drift12m` spans 13.16 deg, so the interior should be back near the **171 px**
of the published f700 run:

```python
import numpy as np
e = np.load("results/hi_experiment_matched700b_f675_725_nc32_ns128.npz")
print(int(np.sum(e["interior"])))     # want ~171, NOT ~57
```

If it is not, the run is void for the same reason and must not be quoted.

**Then the question it answers:** at 350 MHz, matched mode density still leaves
cross-linking worth 9.3x. At 675-725 MHz the drift beats the raster outright
(16.3 against 38.0). Does the drift still lose at matched density where it wins
overall? Read `drift12m_cr_floor_4 / raster_cr_floor_4` from the ablation.
Ratios within the run are exact; do NOT compare against the published f700
numbers, because three strategies in the union means a different HI realisation.

---

**Everything else is committed and pushed.** `ska-drift-analysis` is at
`d12b243` on GitHub, 0 commits ahead.

**Three files are uncommitted**, all from the 2026-10-01 observable-sky fix:
`HANDOFF.md`, `paper/ska_drift/main.tex`, `ska_hi_kmodes.py`. The fix replaced
a hardcoded 41253 deg² (the full sphere) with the sky SKA-Mid can reach —
30 717 deg² at el > 30 deg from the Karoo, derived from `SKA_LAT` with
`EL_MIN_DEG` as the knob. No verdict in the figure changed; the label had
overstated the available sky by a third.

**Paper workflow changed 2026-09-29.** Geoff edits `main.tex` on **Overleaf**
and copy-pastes changes across. So: keep editing the local file, but **hand him
the LaTeX to paste** rather than assuming the local copy is what he is reading.
A bundle was exported to `~/ska_drift_overleaf.zip` (verified to compile
standalone, 11 pages). **Figures are generated and gitignored**, so any figure
that changes must be re-uploaded to Overleaf separately — currently that means
`hi_kmodes.png` (the line moved), `hi_maps.png` and `hi_cleaned.png`.

**A copy of the whole repo sits at `C:\Users\geoff\code\limTOD`** (2.1 GB,
`.venv` excluded), made 2026-09-30 so Claude Science could read it from
Windows. It is a snapshot, not a link — re-run the rsync to refresh it.
`~/research_overview.md` describes this project plus SBI_2026, IM_Gibbs and
SPHEREx_Gibbs for that purpose; it is not in any repository.

**Notebook hazard, seen twice.** `ska_hi_paper_figures.ipynb` is regenerated by
`tools_regen_figures_notebook.py` whenever a figure changes, and a VSCode window
holding an older copy will silently write it back. Git is the recoverable
truth: `git checkout HEAD -- examples/SKA/ska_hi_paper_figures.ipynb`. Do not
confuse it with `ska_hi_figures.ipynb`, which is the superseded 2026-09-23 set.



> **Restored 2026-10-01.** The four sections below were dropped from the
> working copy when the 2026-09-23 "PICK UP HERE" block was rewritten in
> place; they are results records, not status, so they are back verbatim
> from the last commit. Some of the matched-density story is now also
> summarised in PICK UP HERE above.

### Matched density at 675-725 MHz: YES, and the prediction was wrong (2026-09-29)

Tier 1 sweep, job **13889513**, 8 h 04 min, nside 128, spacing 0.25 beams at
675.78 MHz (0.517 deg). Results in `results/hi_ladder_f700_ns128.npz`.

| N strips | tr(WA) | area deg^2 | modes/100 deg^2 | vs raster (76.55) |
|---|---|---|---|---|
| control, published drift | 63.70 | 152.1 | 41.87 | below |
| 3 | 56.38 | 95.5 | 59.05 | below |
| **8** | 126.42 | 153.0 | **82.65** | **ABOVE** |
| 12 | 179.06 | 188.4 | 95.03 | above |
| 20 | 290.55 | 269.4 | 107.84 | above |

**Matched density IS reachable, and comfortably.** Eight strips already exceed
the raster, and the density does not saturate -- it climbs to 107.84 at N = 20.
The crossing is at **N ~ 7** (linear between the two bracketing points, so
worth pinning rather than trusting).

**Why the prediction failed, because the reasoning is worth not repeating.**
This file predicted the drift would top out near 68 and never reach 76.55. That
extrapolated the 350 MHz *density* curve, which saturates near 24 -- but that
curve was measured at **0.5-beam** spacing, where strips are nearly disjoint
and area grows as fast as modes. This sweep is at **0.25 beams**, the
overlapping regime, where area grows far more slowly: N = 3 -> 20 multiplies
modes by 5.2 but area by only 2.8, so density keeps rising. The saturation is a
property of the spacing, not of ladders. Tightening at fixed N is also worth
much more here than at 350 MHz: 1.0 -> 0.25 beams takes 41.87 -> 59.05 (+41%),
against the +8.7% that 0.5 -> 0.25 bought at 350 MHz.

**Next, if the matched-density test is wanted.**

1. **DONE, job 13891716 (2 h 24 min).** N = 6 gives **76.47** modes/100 deg^2
   against the raster's **76.55** -- **0.10% apart**, a tighter match than
   drift12 manages at 350 MHz (2.2%). N = 7 overshoots at 80.12. Results in
   `results/hi_ladder_f700b_ns128.npz`.

2. **Strategy added: `drift6m`** = `ladder(6, 0.12977)`. The spacing is quoted
   in beams at the band bottom because that is what `ladder()` does, so 0.12977
   there IS 0.25 beams at 675.78 MHz. Elevations 48.708 ... 51.292, verified
   against the sweep's own output before submitting.

3. **RUNNING: warm array 13892711** (drift6m, 32 channels, nside 128), then
   experiment + ablation **13892723**, chained `afterok`. Tag `matched700`,
   strategies `drift drift6m raster` in one run so all three share an HI
   realisation.

**What it answers.** At 350 MHz, matching mode density leaves cross-linking
worth 9.3x. At 675-725 MHz the drift beats the raster outright on
`residual/HI` (16.3 against 38.0). This asks whether the drift still loses at
*matched density* even where it wins outright -- i.e. whether cross-linking is
a permanent structural advantage or one that dissolves once the beam narrows.
Either answer is a result; they say different things about what a survey
designer should do above 600 MHz.

**Do not compare against the published f700 numbers.** Three strategies in the
union means a different HI realisation. Ratios within the run are exact.

**RESULT (jobs 13892711 / 13892723, 1 h 37 min) -- AND IT IS NOT USABLE.**
57 interior pixels, common resolution, 4 modes:

| strategy | floor | total | noise |
|---|---|---|---|
| drift (3 strips) | 10.2 | 10.8 | 0.77 |
| drift6m (matched density) | 57.4 | 65.6 | 1.92 |
| raster | 33.7 | 34.2 | 1.12 |

Taken at face value: drift6m/raster = **1.70x** (the drift still loses at
matched density), drift/raster = 0.30x (the known inversion), and drift6m is
**5.6x worse than the 3-strip drift** -- the longer ladder actively hurt.

**Why it must not be quoted.** `drift6m` reaches the raster's density by
packing six strips into **2.58 deg** of declination -- *less than the 3 deg
interior margin*. The common patch across the three strategies is therefore
thin, and after masking the evaluation region is a sliver **0.60 deg** across,
i.e. **0.29 beams** at this band, against 1.17 beams for the published f700
run. Every pixel in it is sub-beam correlated, and it sits where drift6m's own
coverage is falling off. The 5.6x is most likely that, not geometry.

**This is a design error, not a code error.** The sweep matched *density*, and
density is achievable at 700 MHz only by concentrating strips -- which shrinks
the footprint until there is nothing left to evaluate on. Matching density and
keeping an evaluable patch are separate constraints and only one was imposed.

**The fix, if the test is still wanted.** Find a configuration at the same
density but a much wider footprint: more strips at wider spacing. From the
first sweep at 0.25 beams, N = 12 gives 95.03 over a 5.68 deg span, so
widening the spacing at N = 12 should drop the density toward 76.55 while
spanning >10 deg. Needs a spacing sweep at fixed N = 12 (~36 passes,
~660 s each, so ~6.6 h) before any TOD, then warm + run (~4 h). Roughly a day.

**FIXED, and the fix is running (2026-09-30).** Spacing sweep at N = 12,
job 13898503 (4 h 42 min), `results/hi_ladder_f700c_ns128.npz`:

| N = 12, spacing @700 | modes/100 deg^2 | area deg^2 | dec span |
|---|---|---|---|
| 0.50 beams | 82.06 | 297.1 | 11.4 deg |
| 0.70 beams | 68.09 | 377.3 | 16.0 deg |
| 0.90 beams | 56.45 | 462.5 | 20.4 deg |

The raster's 76.55 is crossed at **0.579 beams**, i.e. `ladder(12, 0.3005)`,
added as strategy **`drift12m`**. Elevations 43.418 ... 56.582, **span
13.16 deg**.

**Why this one works where drift6m did not.** The common patch is set by the
*narrowest* footprint among the strategies in the run. drift6m spanned 2.58 deg
-- narrower than drift+raster's own common patch -- so it became the limiting
footprint and the interior collapsed to 0.60 deg. drift12m spans 13.16 deg,
comfortably wider than the published drift's 4 deg, so the common patch stays
set by drift+raster and the interior should return to the ~171 px of the
published f700 run. **Check that first when the results land: if the interior
is not back around 171 px, the run is void for the same reason as before.**

**Queued:** warm array **13920776** (drift12m, 32 channels, nside 128), then
experiment + ablation **13920795**, chained `afterok`, tag `matched700b`,
strategies `drift drift12m raster`. Job **13920752** independently re-measures
the density at exactly 0.579 beams so the matched figure can be quoted rather
than interpolated.

**Superseded:** `matched700` (drift6m) is discarded. Its npz files stay in
`results/` but must not be quoted.

**Diagnostic worth running either way:** `local_modes.py` over the 57 pixels.
If drift6m's *local* mode density there is far below its patch-average 76.47,
that settles it -- the premise fails locally and the number is void. The
operators are on Ilifu only; they are not in the local `hi_cache/`.


### Tying modes to science: n_eff(k_perp) and a sufficiency test (2026-09-28)

`ska_hi_kmodes.py`, figure `hi_kmodes.png`, results in
`results/hi_kmodes_f350_ns64.npz`. Reads cached operators only -- no new runs.

**The decomposition.** Eigen-decompose `A^T N^-1 A`; each mode contributes
`w_i = lambda_i/(lambda_i + S^-1)` to `tr(WA)`. Give each eigenvector an
angular power spectrum by direct non-uniform Fourier transform on the tangent
plane (no regridding), normalised to sum to 1 over k bins. Then
`n_eff(k) = sum_i w_i f_i(k)` and **`sum_k n_eff(k) = tr(WA)` exactly** --
verified at 38.78 and 128.09, so this is a decomposition of the published
number, not a new statistic.

**Result 1: the raster's advantage is nearly scale-independent**, so scale is
not where it lives either.

| | ratio raster/drift |
|---|---|
| large scales, abs(k) < 30 (theta > 12 deg) | 1.44 |
| near the beam, abs(k) > 60 (theta < 6 deg) | 1.81 |
| whole range | 1.30 - 2.11 |

**This is the third axis on which mode count fails to explain the floor.** The
raster beats the drift 1.5x in mode density, 1.8x at best by scale -- and 9.3x
on the floor at matched density. Density, depth (the 2026-09-24 run) and now
angular scale have each been eliminated. What is left is the null-space
*structure*: every drift track has position angle 90 deg.

**Result 2: the drift turns over before the beam.** Its density peaks at
abs(k) ~ 65 (theta ~ 5.5 deg, about 1.4 beams) and falls, while the raster
climbs to the beam at abs(k) = 90.4. A single drift strip cannot sample the
cross-scan direction, so it never reaches beam-scale information.

**Result 3, the sufficiency test -- the useful half for a survey designer.**
For `N_m` independent modes in a k bin, an auto-spectrum gives
`sigma_P/P = sqrt(2/N_m)(1 + P_fl/P_HI)` and a cross-correlation against a
tracer `sigma/P_x = sqrt((P_fl/P_HI)/N_m)/r`. The difference is a **square**,
because in auto the floor is a **bias** (deterministic given the strategy, so
sky area cannot average it away) while in cross it is dominated by unmeasured
*foreground*, uncorrelated with the tracer, hence variance only.

Sky area needed for `sigma/P = 1` at the 350-400 MHz floor:

| | needed | vs now | verdict |
|---|---|---|---|
| raster, cross | **2 937 deg^2** | 5.1x | **buildable -- MeerKLASS scale** |
| drift, cross | 1.1e6 deg^2 | 4 100x | > whole sky |
| raster, auto | 3.9e6 deg^2 | 6 700x | > whole sky |
| drift, auto | 3.5e11 deg^2 | 1.3e9x | not a survey question |

So "raster is better" sharpens to: **the raster is within a factor of 5 of a
useful cross-correlation and the drift is 4 000x away** -- in units a designer
acts on, rather than the 9.3x ratio.

> **SUPERSEDED 2026-10-01 -- do not quote this table.** The two auto rows are
> wrong in kind: the floor is a bias there, so no area suffices and
> `3.9e6`/`3.5e11` are only lower bounds that assume the floor's mean power is
> known and subtracted. The cross rows now carry the cross-power's own sample
> variance and a shot-noise term; see "Fig. 12 right panel rebuilt". Note also
> that "4 000x" above is **vs now** -- 4 100x the 265 deg^2 simulated, NOT 4000x
> the raster. The drift/raster ratio is **370x**. That column was misread into
> the paper's abstract notes and has been corrected there.

**Assumptions to state if this is quoted.** `r = 1` for the cross-correlation
coefficient; Gaussian statistics; the patch-total `N_eff` used as if it were
one k bin, which is optimistic (the per-bin version is the left panel and has
not been folded into the right one); floor/HI taken patch-averaged at 4 modes,
common resolution. The foreground part of the floor is uncorrelated with a
tracer but the signal-loss part is not -- that is what `T(k)` already corrects.


### Depth-matched raster (PLAN item 5, 2026-09-24)

**Tier 1, done.** `ska_hi_raster_depth.py` sweeps the azimuth throw and counts
modes off the operator alone (job 13875729, 20 min).
Results in `results/hi_raster_depth_ns64.npz`.

| throw | px | samp/px vs drift | modes/100 deg² |
|---|---|---|---|
| 7.50° (published) | 684 | 0.46x | 22.31 |
| 5.00° | 550 | 0.57x | 21.22 |
| 4.00° | 485 | 0.65x | 21.20 |
| 3.00° | 420 | 0.75x | 21.13 |
| 2.50° | 411 | 0.77x | 19.96 |
| **2.00°** | **377** | **0.84x** | **19.86** |
| *drift* | *316* | *1.00x* | *14.62* |

**The raster cannot be fully depth-matched by narrowing the throw.** Coverage
is set mostly by the 22.5° of RA that drifts past in 1.5 h and by the beam
(FWHM 3.98°), not by the throw, so 2° -- already sub-beam -- still selects
377 px. It closes the depth gap from 0.46x to 0.84x, which is as far as the
geometry goes. Do not read the remaining 16% as a free parameter.

**The mode-density advantage survives depth-matching**: 19.86 against the
drift's 14.62 at the deepest throw. This is depth-*inclusive*, since more
samples per pixel raise the eigenvalues and hence tr(WA) directly, so the
comparison accounts for depth rather than assuming it away.

**Tier 2, running.** `rasternarrow` (2°, same seeds as the published raster, so
the noise realisation is controlled) added to `ska_hi_experiment.py`. Cache
warm array **13876390** (32 tasks), then experiment + ablation **13876408**,
chained with `--dependency=afterok` so it cannot run against a partial cache.
Tag `depth`, strategies `drift raster rasternarrow` -- all three in ONE run so
they share a single HI realisation and the three-way comparison is internally
consistent.

**Do not compare `depth`'s drift/raster numbers with the published ones.**
`build_hi_cube` sizes its box from the union of a run's strategies, so a
different strategy list is a different HI realisation (~15-18%, see the
realisation-scatter section). Within the run it is exact; across runs use
ratios.

**Tier 2 result (jobs 13876390 / 13876408, warm 96/96 COMPLETED, run 36 min).**
One HI realisation, 96 interior pixels, common resolution, 4 modes removed:

| strategy | depth vs drift | modes/100 deg² | **floor** | total | noise | 1/f |
|---|---|---|---|---|---|---|
| drift | 1.00x | 14.62 | **338.6** | 348.9 | 6.57 | 3.96 |
| raster (7.5°) | 0.46x | 22.31 | **22.6** | 24.6 | 2.17 | 0.40 |
| rasternarrow (2°) | 0.84x | 19.86 | **36.0** | 36.9 | 1.78 | 0.57 |

**The answer is geometry, not depth, and it is unambiguous.** The depth
confound, if it were real, predicts that giving the raster the drift's depth
shrinks its advantage. Instead:

- **At essentially matched depth (0.84x) the raster's floor is still 9.4x
  better than the drift's.** That is the number the referee question asks for.
- **Making the raster deeper made it WORSE**, 22.6 -> 36.0, a factor 1.60 --
  because narrowing the throw costs modes (22.31 -> 19.86 per 100 deg²) and
  the floor follows the modes, not the depth. Depth and floor move in
  *opposite* directions here.

**Two independent controls now agree.** Matching mode density by lengthening
the drift's ladder (`drift12` vs `raster`) gives 9.33x; matching depth by
narrowing the raster (`drift` vs `rasternarrow`) gives 9.40x. They control
different things and are not equivalent -- the depth-matched pair still differs
in mode density -- but both land on ~9.4x from opposite directions, where the
uncontrolled comparison gives 15.0x.

*Sanity check that the depth change took effect:* the noise term falls
2.17 -> 1.78 (x0.82) against x0.74 predicted by 1/sqrt(samples per pixel).
Same direction, right order; the residual difference is the changed mode
content, not an error.

*Do not compare these absolutes with the published run* -- 96 interior pixels
against 119, and a different HI realisation. Ratios within the run are exact.


### Section 5 is blocked on a prior, not on figures (found 2026-09-24)

Filling the paper's empty "Foreground reconstruction" section (experiments
001-005) ran into this: **every `summary_*.png` is dated 2026-07-31 and is
therefore a smoothed-prior solve**, while the flat-prior audit of 2026-08-13
concluded "quote the **flat-prior** gains when the claim is about what the
survey measures". The two disagree by a lot:

| off-plane, no HP, gauss | per-pixel | beam-scale |
|---|---|---|
| smoothed prior (what the figures show) | -10.0% | -37.8% |
| flat prior (what the text should quote) | **-32.6%** | **-46.4%** |

So dropping `summary_offplane_strategy.png` in and writing the flat-prior
numbers around it would reproduce exactly the `hi_ablation` caption failure:
a caption describing numbers its own figure does not show.

`audit_flat_prior.py` re-solves this in ~1 min but **writes no figures** --
it is numbers only, and only for the experiment-004 cross-linking comparison,
not for the noise budget or the residual RMS. The three figures §5 wants
(`summary_residual_rms`, `summary_noise_budget`, `summary_offplane_strategy`)
all come from `ska_results_summary.ipynb`, on the smoothed prior.

Three ways out, none of them free:

1. **Regenerate the three under the flat prior.** Extend `audit_flat_prior.py`
   to emit the budget and RMS, then re-plot. Makes §5 self-consistent and is
   the only option that lets the text quote the numbers HANDOFF says to quote.
2. **Caption them explicitly as smoothed-prior** and quote smoothed-prior
   numbers in §5. Cheap and honest, but it contradicts the 2026-08-13 guidance
   and invites the question of why two sections use two priors.
3. **Write §5 from numbers, no figures.** Defers the choice.

Note the qualitative claim is safe under either prior and in fact strengthens
under the flat one: withdrawing the prior's structure costs the drift 40% and
the raster 5%, so the floor grows and "1/f is not the limiter" gets *stronger*.
It is the quoted percentages that move.

**Added meanwhile:** `meerklass_scan_tracks.png` into Sec. 2.3. It is pure
geometry -- boresight tracks plus the drift patch -- so it is prior-independent
and safe, and it is the paper's only picture of the crossing that the whole
argument rests on.

**Also unresolved, and pre-existing** (from the 2026-08-13 audit, repeated here
because §5 is where it would be written): off-plane at beam scale the raster's
floor *exceeds* its total, 0.187 vs 0.181 K smoothed and 0.210 vs 0.199 K flat.
Adding noise slightly reduces the beam-scale residual, so the floor/noise
quadrature split does not hold there. Understand this before quoting any
beam-scale noise term for this field.


**Still to do in the text:** the earlier sections of *this file*, above
"Frequency repeat at 675-725 MHz", still state the cross-linking result
without the band qualifier. The prose of the paper is still a skeleton --
this pass changed figures, captions and numbers, not the argument.

**Higher bands would need nside 256.** 850-900 MHz gives 3.4 px/FWHM at nside
128 and 1025-1075 gives 2.8, both below this series' own safe floor (see the
grid check above). That is a much larger job than the 700 MHz repeat and is
not queued. Three points across Band 1 is arguably enough to make the
frequency argument without them.

**Ilifu reminders.** The ControlMaster expires after 8 h, so start with
`ssh ilifu` **in a real terminal** (the OTP prompt needs a TTY; an agent
cannot do it). `squeue`/`sacct` are not on the transfer node and the login
shell is not a login shell by default, so wrap: `ssh ilifu "bash -lc '...'"`.
Result npz files are small enough to pull with
`ssh ilifu "bash -lc 'cat <remote>'" > <local>`.

### Fig. 12 right panel rebuilt (added 2026-10-01)

`hi_kmodes`'s right panel used to quote a sky area for the auto-spectrum, which
contradicted its own caption: `sufficiency()` applied the variance-only count
`2(1+R)^2` to a term the same figure called a bias. It now says what is true --
**the auto rows carry no area at all**, only the bias `R = P_floor/P_HI` (655x
for the raster, 1.6e5x for the parked drift). The `2(1+R)^2` count survives in
the returned dict under `auto_lower_bound_if_floor_known`, because that is all
it ever was: a bound that assumes the floor's mean power is known and
subtracted, i.e. that you know the sky the survey failed to measure. **Do not
quote 3.9e6 deg^2 as an auto requirement.**

The cross rows now use the full Gaussian cross-spectrum variance,

    (dP_x/P_x)^2 = [1 + (1 + R)(1 + 1/(nbar P_g)) / r^2] / N_m,

parameterised by `R_CROSS`, `NPG = (inf, 1.0)` and `SNR = (1, 3, 5)`. The extra
leading 1 (the cross-power's own sample variance) raises the count from `R` to
`2 + R`, which is 0.31% for the raster -- `main()` asserts the 2937 deg^2 the
draft quotes survives, and it does, at 2946.

| row | nPg | 1 sigma | 3 sigma | 5 sigma |
|---|---|---|---|---|
| raster | inf | **2 946** | **26 515** | 73 653 |
| raster | 1 | 5 888 | 52 990 | 147 193 |
| 12-strip drift | inf | 2.63e5 | 2.36e6 | 6.56e6 |
| 3-strip drift | inf | 1.09e6 | 9.83e6 | 2.73e7 |

Bold is inside the reachable sky (30 717 deg^2), so the raster cross-correlation
is buildable **to 3 sigma, not 5**. The twelve-strip drift is the new row and
the point of it: at the raster's mode density it still needs 90x the raster's
sky, so matching mode count does not buy the science back.

**One number in the draft was wrong and is now fixed**: the abstract notes said
the drift needs "4000x more sky than the raster" in cross-correlation. It is
**370x** (1.09e6 / 2946) -- 4000 does not come out of any version of the
formula. Both planning comment blocks in `main.tex` were rewritten.

drift12 has no cached operator (the ladder sweep discards its operators), so
`drift12_geometry()` reads n_eff 147.18 / 674.8 deg^2 from
`results/hi_ladder_long_ns64.npz`, falling back to the tab:modes values. It
appears in the right panel only -- there is no k decomposition for it.

Notation: the x label is now a plain area, the dP/P definition moved into the
legend, and the reference line is the reachable sky with 4pi drawn faint behind
it. `MUTED` is now within a hair of `INK`, so reachability is carried by alpha
(`DIM`), not colour.

**Band-parameterised, so 675-725 MHz can be drawn too (2026-10-01).**
`ska_hi_kmodes.py` now has a `Band` dataclass and `BANDS = {"f350", "f700"}`,
and `main()` takes `--band` and `--draw-only`. Everything band-specific is in
the config: frequency, nside, the matched ablation npz the floors come from,
the strategy names (`drift12` at 350, `drift12m` at 700), the ladder fallback
row, the output npz, the figure suffix and the axis limits. `f350` is
bit-for-bit unchanged -- both PNGs reproduce byte-identically, and the 2 946
deg^2 assertion still passes at 0.31%.

Two things worth knowing:

* `FLOOR_AMP` is gone. Floors are now read from the band's ablation npz by
  `floor_amp()` and **rounded to 1 dp**, which reproduces the old literals
  exactly (399.6 / 239.3 / 25.6) while making them impossible to leave stale.
  The rounding is deliberate: the ratio enters squared, so further digits are
  noise against 15-18% realisation scatter.
* `compute()` needs the cached operators, so at 700 MHz it runs on Ilifu.
  `save()`/`load()` split the work: everything the figures need goes into
  `results/hi_kmodes_<band>_*.npz`, and `--draw-only` redraws from that alone.
  Pull the npz, draw locally.

**Queued: job 13985514 (`ska-kmodes`), `--dependency=afterok:13920795`.** It
fires when the matched700b ablation lands and writes
`results/hi_kmodes_f700_ns128.npz`. `scancel 13985514` to drop it. Then:

```bash
# pull, then
/home/geoff/limTOD/.venv/bin/python ska_hi_kmodes.py --band f700 --draw-only
```

All three f700 operators are confirmed present on Ilifu (`op_{drift,drift12m,
raster}_f675.781_ns128.pkl`), and the ladder cross-check for the 12-strip row
resolves to n_eff 251.20 over 331.73 deg^2. Note Ilifu's `ska_hi_plots.py` was
**stale** (no `MONO_STRAT`) and has been updated along with the new
`ska_hi_kmodes.py` and `run_kmodes.sbatch`; `ska_hi_ladder.py`, `ska_common.py`
and `ska_hi_experiment.py` already matched. Anything run there needs
`export FASTBOX_PATH=$HOME/ska/FastBox` or `ska_hi_mock` refuses to import.

**What to expect at 700 MHz.** Required area goes as (2 + R)/density with
R the SQUARED floor. The floors invert between bands (drift 16.1 vs raster 36.0
at 675-725, against 399.6 vs 25.6 at 350-400), so the drift rows should come
out BELOW the raster and the figure should reorder. Do not present the 350 MHz
version as a general statement -- see the band table above.

**There is no mock tracer, and the figure now says so (2026-10-01).** Nothing
in this repo simulates a galaxy survey -- no catalogue, no n(z), no bias, no
footprint. The only mocks are HI: the fastbox `HITracer` cube (b_HI 1.559,
Tb 0.518 mK, with RSD) and the independent HI mock `X_m` that `T(k)` uses. The
cross-correlation is pure algebra, and the tracer enters through exactly two
numbers, `R_CROSS` and `NPG`. So "galaxy survey" was replaced by **"ideal
tracer"** in the row labels and the legend, and the caption now states plainly
that no tracer is simulated.

Assumptions a referee will go for, in order: the tracer is granted the HI map's
own `N_eff` over the same patch (no window of its own); the whole patch is one
band power, with no k binning, which is the most optimistic aggregation
possible; Gaussian errors; floor/HI patch-averaged at 4 modes, common
resolution; the HI map's thermal noise is NOT in the cross variance (10-52x
below the floor, small but absent).

**The redshift is the weak point.** 350-400 MHz is z = 2.55-3.06, D_c = 6155
Mpc. There is no wide-area spectroscopic *galaxy* survey there -- the plausible
tracers at z ~ 2.8 are the Lyman-alpha forest, quasars, or LBGs. `nbar P_g = 1`
over the accessible k (k_perp <= 0.016, k_par 0.012-0.192 Mpc^-1) needs nbar of
roughly 1e-5 to 4e-4 Mpc^-3 depending on k and bias (back-of-envelope with the
pipeline's own cosmology; not in the paper). If the claim is ever to name a
tracer, that is the number to check against a real n(z).

**Split into two figures, same day.** The sufficiency panel was too dense beside
a log-log panel, so `fig_kmodes` is now the angular-scale panel alone
(`hi_kmodes.png`) and `fig_sufficiency` is the sky-area panel alone
(`hi_sufficiency.png`); both call shared `_panel_*` helpers, the same split
pattern as `fig_design`/`fig_levers`. In the paper they are `fig:kmodes`
(single column) and `fig:suff` (`figure*`, 0.92 textwidth). The generator
registry has both, so a notebook regeneration keeps them.

**main.tex was re-synced from Overleaf on 2026-10-01.** Geoff writes there and
codes here; the local copy had drifted badly (different captions, no `\iffalse`
blocks, plain `HI` where the draft uses `\textsc{Hi}`). The local file is now a
mirror of the Overleaf version plus the two new figure blocks, the corrected
abstract sentence and the corrected Discussion comment. **Do not edit prose
here** -- hand Geoff tex to paste instead, or the two diverge again. Note the
draft includes `hi_patch_maps` and `hi_transfer_function` as `.pdf`, and
`_save()` writes only `STYLE["fmt"]`, so a `.png` regeneration leaves those two
stale without saying so.

### The ladder is penalised for its own extent (found 2026-10-01)

**The prior is set per strategy, from that strategy's own patch.**
`ska_hi_analysis._solve_args` builds `S^-1 = I / (PRIOR_VAR_SCALE *
var(fg_truth_patch))`, and `fg_truth_patch` is the pixels THAT strategy
observes. So a longer declination ladder spans more sky, picks up more
structure, gets a larger sigma_p, and is handed a LOOSER prior -- which
`fig_priorvar` shows raises tr(WA) and *worsens* the floor. A ladder is
therefore penalised for its extent, independently of its geometry. This is a
confound in every drift-vs-drift12 comparison in the series.

Measured on the real footprint at 675-725 MHz (RA 150.9-165.9 deg, the 15 deg
of drift, plus a beam margin):

| | dec band | npix | <T> | sd | max |
|---|---|---|---|---|---|
| drift (3 strips) | 7.30-11.28 | 480 | 3.476 K | **0.321 K** | 4.83 |
| drift12m (12) | 2.70-15.87 | 1200 | 3.571 K | **0.401 K** | 4.83 |

So **not a brighter sky** -- 3% in the mean, identical maximum, and over the
full RA ring the wider band is marginally COOLER (8.08 vs 8.14 K). What differs
is structure: sigma_p is 1.25x wider, i.e. the prior variance is **1.56x
looser** for drift12m than for drift.

Direction of the effect is right for the puzzle: at 675-725 drift12m's floor is
20.8x against the 3-strip drift's 16.5x, i.e. WORSE with 12 strips (and with 4x
the integration time -- the ladder gives one strip per sidereal night, so 12
strips is 12 nights against 3). Whether 1.56x of prior accounts for a 26% floor
difference is not yet known: the prior sweep is at 350 MHz over 1e-2 to 1e3, and
the local slope at 1.0 has not been read off.

**ANSWERED 2026-10-01 by job 13990883 (`ska-pin`), 3:36:15, MaxRSS 31.8 G.**
`run_ablation_pin.sbatch pin700 drift drift drift12m raster` ->
`results/hi_ablation_pin700_f675_725_nc32_ns128.npz`. Same 713/171 px as
matched700b, verified identical.

| strategy | prior scale | unpinned floor | pinned floor | change |
|---|---|---|---|---|
| drift (reference) | 1.0000 | 16.48 | **16.48** | **+0.0%** |
| drift12m | 0.6817 | 20.79 | 19.27 | -7.3% |
| raster | 0.5421 | 30.56 | 27.17 | -11.1% |

**The control is exact**: drift reproduces 16.48 to 2 dp and every printed digit
of its log line (16.8 / 16.9 / 1.3 / 0.1 / 1.3). So the comparison is sound.

**The confound is real but explains only about a third of the anomaly.**
drift12m/drift goes 1.261 -> 1.169, i.e. 35% of the excess is the prior and
**two thirds is not**. At a common prior the twelve-strip ladder is still 17%
worse than three strips at this band, with 4x the integration time. The
ordering is unchanged -- drift < drift12m < raster -- so the 675-725 conclusion
survives the correction.

**The raster is penalised MOST (-11.1%), because its patch is the widest.**
That matters more for 350-400 than anything here: pinning there would improve
the raster more than the drift, so **the published cross-linking gains are
conservative** -- 9.3x and 14.4x are lower bounds, not best cases. NOT yet
measured at 350; a pinned run at nside 64 is cheap and would settle it.

**This also bears on 350-400**, where drift12 spans the same 13.16 deg. There
drift12 came out BETTER (239.3 vs 399.6), so the looser prior did not stop it
winning -- but the measured cross-linking gain is understated by whatever this
costs, and Table 2's mode densities are quoted at each strategy's own prior.

### One fixed prior per band (added 2026-10-02)

**Definition.** For each band and channel, `sigma_p^2` is the variance of the
GDSM truth over the **three-strip drift's own patch**, stored in
`results/prior_ref_f{lo}_{hi}_ns{nside}.npz` and read by every run in that
band. Not the common patch: that is the intersection over whichever strategies
are in the run, so it moves when the strategy list does and is not a property
of the band. At 675-725 the two agree to **0.6%** (0.3397 vs 0.3416 K at
channel 0), because the drift is the narrowest strategy, so `pin700` stands as
measured. `PRIOR_VAR_SCALE` stays 1.0.

**New module `ska_hi_prior.py`**, with a single switch `prior.PINNED`, OFF by
default so every published number reproduces untouched -- merely building a
reference file cannot move Table 2. Runs opt in with `--pin-prior`.

**Every place S^-1 is set, and what happened to it:**

| site | what it feeds | routed? |
|---|---|---|
| `ska_hi_analysis._solve_args` | every map-making solve | yes |
| `ska_hi_modes.py:90` | Table `tab:modes` | yes |
| `ska_hi_ladder.py:145` | the ladder sweeps | yes |
| `ska_hi_kmodes.eig_decomposition` | hi_kmodes, hi_sufficiency | yes |
| `ska_hi_design.freq_mode_density` | hi_design, hi_levers | yes |
| `local_modes.py:54` | prototype, not a paper figure | **no** |
| `plot_freqsweep_maps.py:106` | experiment 005, different geometry | **no** |
| `audit_flat_prior.py:103` | experiments 004/005 audit | **no** |

The last three are outside the HI comparison chain and run at bands with no
reference; they are listed so nobody has to re-derive this. `ska_hi_priorvar.py`
and `ska_hi_plots.py:847` only *consume* `s_inv`, so they are untouched and the
prior sweep still works.

`solve`, `response` and `validate_linearity` now take `freq_mhz`; a pinned solve
without one **raises** rather than silently using the patch, because that is the
failure mode worth catching. Verified: unpinned, `hi_kmodes.png` and
`hi_sufficiency.png` are byte-identical to before the routing.

**Step 0 finding.** Job 14001763 (850-900 MHz) was submitted before any of this
and would have used the per-strategy prior. Cancelled and replaced.

**Runs queued 2026-10-02** (all via the new `run_hi_band.sbatch
<bandenv|default> <pin|nopin> <tag> <strategy>...`, experiment + ablation):

| job | tag | band | nside | strategies | unpinned reference |
|---|---|---|---|---|---|
| 14002898 | `pin350` | 350-400 | 64 | drift, drift12, raster | `matched` (same set) |
| 14002899 | `pin500` | 500-550 | 128 | drift, raster | `f500` (same set) |
| 14002900 | `f875` | 850-900 | 128 | drift, raster | -- it IS the reference |
| 14002901 | `f875pin` | 850-900 | 128 | drift, raster | 14002900, `afterok` |
| 14002942 | mode check | both | -- | -- | prints own vs pinned |

`pin700` is done (job 13990883) and used the three-strip drift's patch, which is
the definition above. **Control to check on every pinned run: the three-strip
drift must reproduce its unpinned numbers exactly**, as it did at 700 MHz.

**Step 3 DONE (job 14002942, 2:59).** Mode counts at channel 0, own prior vs
pinned. Threshold counts are prior-independent and reproduce Table 2 exactly
(17/37/52/33), and the three-strip drift is +0.0% at both bands -- control.

| band | strategy | area | n_eff own | n_eff pin | change | /100 deg^2 own -> pin |
|---|---|---|---|---|---|---|
| 350-400 | drift | 265.2 | 38.78 | 38.78 | **+0.0%** | 14.62 -> 14.62 |
| | drift8 | 499.4 | 100.05 | 98.94 | -1.1% | 20.03 -> 19.81 |
| | drift12 | 674.8 | 147.18 | 146.09 | -0.7% | 21.81 -> 21.65 |
| | raster | 574.1 | 128.09 | 119.39 | **-6.8%** | 22.31 -> 20.80 |
| 675-725 | drift | 152.1 | 63.70 | 63.70 | **+0.0%** | 41.87 -> 41.87 |
| | drift12m | 331.7 | 251.19 | 240.84 | -4.1% | 75.72 -> 72.60 |
| | raster | 383.3 | 293.45 | 267.09 | **-9.0%** | 76.55 -> 69.67 |

**The raster loses the most mode count**, at both bands, for the same reason it
gained the most floor: its patch is the widest, so its prior was the loosest and
it was counting the most prior-filled directions as data-dominated.

**The density match holds, but crosses over.** drift12 / raster density goes
0.9775 (-2.2%) -> **1.0410 (+4.1%)** at 350-400, and drift12m / raster goes
0.9892 (-1.1%) -> **1.0420 (+4.2%)** at 675-725. Both inside the ~5% tolerance,
so **no ladder re-tuning**. But the sign flips: the twelve-strip drift no longer
*approaches* the raster's density from below, it slightly **exceeds** it. If the
pinned floors hold up, the matched-density argument gets stronger, not weaker --
more modes per unit sky AND a far worse floor.

`prior_ref_f675_725_ns128.npz` reports sigma_p 0.2818-0.3416 K, identical to
what job 13990883 printed, which **confirms pin700 used this definition**.

**Failures 2026-10-02, all diagnosed:**

* `pin500` (14002899) **FAILED** -- and the guard is why. `ska_hi_ablation.py`
  line 158 called `A.response(mm, hi, truth)` with no `freq_mhz`, so the pinned
  solve raised instead of silently using the patch. One site, missed because it
  is an `A.response` and I had patched only the `A.solve` calls. Found every
  remaining one by AST rather than grep (multi-line calls hide from grep).
  Its **experiment stage completed** and wrote `hi_experiment_pin500_*`, so only
  the ablation is re-run, via the new `STAGE=ablation` env on
  `run_hi_band.sbatch`.
* `pin350` (14002898) **TIMEOUT** at 4 h: drift finished in 609 s, drift12 did
  not. My 4 h was far too tight -- the unpinned `matched` run (13857817) took
  3:33:41 for the experiment alone. Resubmitted at 12 h.
* `f875` (14002900) **TIMEOUT** at 6 h having never reached the first strategy:
  it was still building operators. **Cache warming at nside 128 is ~2.5 h per
  channel**, which is why 700 MHz was warmed by a 32-task ARRAY (`ska-warm7`,
  13920776), not in the run. So 850-900 needs
  `sbatch --array=0-31 warm_caches_f875.sbatch drift` and the same for raster --
  64 tasks -- before any run. Script is in place; NOT submitted pending a
  decision on the cost.

**Re-queued:** 14024831 `pin350` (12 h), 14024832 `pin500` ablation-only (4 h).

**Expectation, to be checked not assumed:** at 350-400 the twelve-strip drift
has the widest patch (674.8 deg^2 against the raster's 574.1 and the three-strip
drift's 265.2), so pinning should help it most and the **9.3x should shrink**
while the **14.4x grows**.

### Pinned-prior results, all bands (added 2026-10-05)

Jobs 14024831 (`pin350`, 7:12:31) and 14024832 (`pin500` ablation, 42:44)
COMPLETED; `pin700` was job 13990883. All at 4 modes removed, common
resolution, interior pixels, against an unpinned run with the SAME strategy set
and the same pixel set (verified). **All ratios below are AMPLITUDE.**

| band | strategy | floor/HI own | pinned | change | total/HI own | pinned | change | N_eff/100deg2 own | pinned | change |
|---|---|---|---|---|---|---|---|---|---|---|
| 350-400 | drift | 399.61 | 399.61 | +0.0% | 413.49 | 413.49 | -0.0% | 14.62 | 14.62 | +0.0% |
| 350-400 | drift12 | 239.29 | 236.20 | -1.3% | 268.49 | 262.97 | -2.1% | 21.81 | 21.65 | -0.7% |
| 350-400 | raster | 25.64 | 26.79 | +4.5% | 27.08 | 37.94 | +40.1% | 22.31 | 20.80 | -6.8% |
| 500-550 | drift | 214.30 | 214.30 | +0.0% | 213.17 | 213.17 | +0.0% | -- | -- | -- |
| 500-550 | raster | 56.24 | 47.87 | -14.9% | 76.72 | 69.47 | -9.5% | -- | -- | -- |
| 675-725 | drift | 16.48 | 16.48 | +0.0% | 16.63 | 16.63 | +0.0% | 41.87 | 41.87 | +0.0% |
| 675-725 | drift12m | 20.79 | 19.27 | -7.3% | 22.17 | 20.74 | -6.4% | 75.72 | 72.60 | -4.1% |
| 675-725 | raster | 30.56 | 27.17 | -11.1% | 32.21 | 29.14 | -9.5% | 76.55 | 69.67 | -9.0% |

| ratio (AMPLITUDE, not power) | band | unpinned | pinned | change |
|---|---|---|---|---|
| 12-strip drift / raster, floor (the **9.3x**) | 350-400 | 9.33x | 8.82x | -5.5% |
| 3-strip drift / raster, total, common res (the **14.4x**) | 350-400 | 15.27x | 10.90x | -28.6% |
| 12-strip / 3-strip drift, floor | 675-725 | 1.26x | 1.17x | -7.3% |
| 3-strip drift / raster, total, common res | 500-550 | 2.78x | 3.07x | +10.4% |

**The control holds at every band**: the three-strip drift is +0.0% on floor,
total and N_eff, to every digit. So the mechanism is sound and every change
below is the prior and nothing else.

**What moved, and the direction.** The reference is the drift's patch, which is
the NARROWEST, so pinning TIGHTENS the prior for everything else (at 350-400 the
raster's own sigma_p is 2.7693 K against the reference 2.2867, i.e. prior
variance x0.682). Tighter means fewer data-dominated modes -- hence N_eff falls
for every non-reference strategy, most for the raster (-6.8%, -9.0%), which has
the widest patch.

**The 9.3x shrinks to 8.82x (-5.5%), as predicted.** The twelve-strip drift's
own floor barely moves (-1.3%) while the raster's rises, so the matched-density
gap narrows. Combined with the mode counts above -- drift12 now at 1.0410x the
raster's density rather than 0.9775x -- the matched-density claim is intact and
if anything cleaner: MORE modes per unit sky and still ~9x the floor.

**The 14.4x does NOT grow. It falls to 10.90x (-28.6%).** That is the opposite
of the expectation, and it rests on a number that is **not yet trustworthy**.

#### RESOLVED 2026-10-05: pinning redistributes noise in k_par, at 350 MHz only

The +334% was never a bug. These arrays are per **k_par bin** (16 of them), and
the series' summary statistic is the median over bins -- so a change in the
SHAPE of the spectrum moves the median even when the total power falls.

Per-bin ratio of `cr_noiseonly_4` POWER, pinned / unpinned:

| band | strategy | low k | mid k | high k | bins with total > 1.2x floor (pin / unpin) |
|---|---|---|---|---|---|
| 350-400 | drift12 | 0.83 | 3.59 | 4.97 | 1/16 / 2/16 |
| 350-400 | **raster** | **0.56** | **22.9** | **52.7** | **9/16 / 0/16** |
| 500-550 | raster | 0.61 | 0.79 | 1.82 | 11/16 / 7/16 |
| 675-725 | drift12m | 0.51 | 0.52 | 0.48 | 0/16 / 0/16 |
| 675-725 | raster | 0.43 | 0.41 | 0.37 | 1/16 / 2/16 |

The drift control is **exactly 1.00 in every bin at every band**, so the
mechanism is sound.

**At 675-725 pinning is clean**: the noise falls uniformly across k (0.37-0.52),
which is the expected damping from a tighter prior, and both strategies stay
floor-limited. The direct one-channel test agrees -- map-space |noisy - quiet|
rms 0.531 -> 0.424.

**At 350-400 it redistributes**: low k falls as expected (0.56) but mid and high
k rise 20-50x from a very low base, and the raster stops being floor-limited
above k ~ 0.06 -- 9 of 16 bins, against 0 of 16 unpinned. drift12 shows the same
signature ~10x weaker, so it scales with how far the prior moved, and it is a
property of the BAND, not of the raster.

Not the explanation: the spectral slope of the reference. The raster/drift
sigma_p ratio runs 1.2031-1.2760 across 350-400, a smooth monotonic 6% -- which
would show up at LOW k, not high. More likely the prior change moves modes
across the data/prior crossover (the raster loses 8.7 of 128 modes), and
marginal modes are the ones most sensitive to the per-channel operator changing
as the beam narrows 12.5% across this band -- the largest chromaticity in the
series. **That is a hypothesis; the measurement above is not.**

**Operational rule:**

* Mode counts and FLOORS: use the pinned numbers, every band. The 350 floor
  redistributes mildly too (about +/-15% across k, median +4.5% for the raster),
  so quote 8.82x rather than 9.33x but do not treat the third digit as safe.
* TOTAL residual at 350-400: **use the unpinned run.** The pinned total is
  contaminated by noise the prior moved into the measured k range, and 10.90x is
  not a number to put in the paper.
* 500-550 and 675-725: pinned is fine throughout.


### Amplitude-only prior, and the shared white-noise bug (added 2026-10-05)

**Two independent corrections, both now implemented.**

**(1) The pinned prior shared the wrong thing.** Imposing the reference's
per-channel sigma_p(nu) also imposes its spectral SHAPE. Measured, noise power
before any cleaning, pinned/own:

| | low k | mid k | high k | map rms |
|---|---|---|---|---|
| shape-pinned | 0.68 | **13.54** | **132.70** | -21.3% |
| amplitude-only | 0.62 | 0.61 | 0.58 | -21.7% |

Same total noise removed; only the shape-pin redistributes it. Cause: the
Wiener gain lambda/(lambda + sigma_p^-2) depends on lambda, which varies
channel to channel as the beam narrows 12.5% across 350-400, so a
frequency-varying scale (0.6142-0.6908, a 12.5% swing) changes which modes sit
near the data/prior crossover per channel. The per-mode bound
lambda/(lambda+sigma_p^-2)^2 constrains TOTAL noise power, not its distribution
in k_par, so it forbids none of this.

New definition: `sigma_p,s^2(nu) = sigma_p,s^2,own(nu) * c_s`, `c_s` the
geometric mean over channels of ref/own. Controls pass: drift `c = 1.000000`
and bit-identical (0.000e+00); raster `c = 0.654902`, residual shape spread
**6.1%** -- what we deliberately leave in, and what the paper should quote.

**(2) The white noise is identical in every channel.** `simulate_foreground`
re-seeds with `p["seed"]` (= SEED + night, NO frequency dependence) and calls
`generate_TOD` one channel at a time, so every channel draws the same white
realisation. The CACHED TODs carry it too -- `replay_noise` is faithful, not
out of step. Empirical: ch0/ch1 TOD ratio scatters 1.21e-4 where independent
draws would give 2.24e-3, 18x larger.

A common 1/f gain across the band is right for one receiver and is KEPT.
Identical thermal noise is not: following the sky's smooth spectrum it is
nearly rank-1, so PCA removes most of it and the published post-clean noise is
optimistic. New arms `noiseonly_indep` and `total_indep` redraw the white part
per (pass, channel) with `SeedSequence([seed, channel])` -- not `hash()`, which
is signed and out of range for `np.random.seed`. Published arms untouched.

**What each correction touches** (they are NOT the same set):

| | amplitude-only prior | shared white noise |
|---|---|---|
| floors, N_eff, sufficiency/sigma levels | **redo** | safe |
| total residual, noise arms, hi_frequency_ratio | **redo** | **redo** |
| dead levers (time, Tsys, dishes) | -- | **at risk** |

The sigma levels are floor-derived (`sufficiency()` takes only floor_amp, n_eff,
area; both the floor arm and `p_hi` are noiseless), so the noise bug does not
reach them. The dead-levers claim is what the noise bug most threatens: it
holds because the residual is 96-99% floor, and noise/floor is already 7.7% at
675-725 where the drift/raster inversion lives.

**Ablation-only re-runs are sufficient.** The ablation's `cr_total_4/cr_p_hi`
reproduces the experiment's `p_clean_auto_4/p_mapmade` to **0.0%** at both
bands checked, and every quantity in the report comes from the ablation
(N_eff from the mode check, operators only). That drops pin350 from ~7 h to
~1-2 h. NOTE `ska_hi_design.load_bands` reads the EXPERIMENT npz, so
`hi_frequency_ratio` will keep plotting shared-noise numbers until it is
repointed -- left alone deliberately, no figure edits until the table is
reviewed.

**Queued 2026-10-05:** 14053756 `amp350` (drift, drift12, raster), 14053757
`amp500` (drift, raster), 14053758 `amp700` (drift, drift12m, raster), all
ablation-only with `--pin-prior`; 14053380 `f875pin` released; 14053379 `f875`
waits on the last warm tasks. f875/f875pin pick up the new code, so both bands
get the independent-noise arms.

### Section A/B of the minor TODOs, and two controls (added 2026-10-07)

Claude Science sent a list of minor checks. All of Section A (numbers) and
Section B (figures) are done. The write-ups live in `~/handover/` --
`paper1_minor_todos_A.md` and the MANIFEST inside
`~/handover/ska_figures_sectionB.zip` -- because they are addressed to the
draft, not to this repo. Only what changes the CODE or would be re-derived
from scratch is recorded here.

#### The ablation noise model (the one real decision)

`fig_ablation`'s old caption said noise sits 10-52x BELOW the floor, so
integration time is a dead lever. That reading came from the `whiteonly` /
`noiseonly` arms, and those are **channel-identical**: `simulate_foreground`
reseeds with `p["seed"]` (= SEED + night, no frequency dependence) once per
pass and calls `generate_TOD` one channel at a time, so every channel draws the
SAME white realisation. Channel-identical noise is spectrally smooth, lives in
the same few modes as the foreground, and PCA removes most of it. The arm was
measuring a cleaning artefact.

`noiseonly_indep` / `total_indep` (added 2026-10-05) draw per channel from
`SeedSequence([seed, channel])`. At 350-400, 4 modes removed, cr_ interior:

| | floor | noise, white shared | noise, white independent |
|---|---|---|---|
| 3-strip drift | 399.6 | 8.3 | 1025.7 |
| 1-dish raster | 22.3 | 1.8 | 2050.1 |

A factor of **124x** (drift) and **1161x** (raster) purely from the correlation
assumption, and the honest arm sits ABOVE the floor, not below it.

**Decision: draw both.** `fig_ablation` now plots `floor`, `total`,
`gainonly`, `noiseonly` (dashed) and `noiseonly_indep` (solid, HIGHLIGHT red),
and the gap between the last two is the figure's point. The floor argument is
untouched -- the floor is *defined* as what does not integrate down -- but
"noise is a dead lever" is gone: at 3 h the maps are noise-dominated, and
integration stays live until noise crosses the floor.

Pairing trap, hit and fixed before rendering: `noiseonly_indep` applies
`(1 + g)(1 + w_ind)`, so it carries the 1/f term and its counterpart is
`noiseonly`, NOT `whiteonly`. 1/f is negligible against independent white
(5.1 vs 1025.7 for the drift) so the curve hardly moves, but the label would
have been false.

#### Matched-depth control (PLAN item 5) -- 9.4x becomes 9.93x

New ablation `rnarrow`, strategies `drift raster rasternarrow`, amplitude-only
prior, 350-400 MHz nside 64. Job 14102350, **45 min** (the 2 h 17 estimate was
pessimistic). Pure replay -- all 32 operators and the TOD caches already
existed, no new TOD.

| | floor | drift / x |
|---|---|---|
| drift | 338.58 | -- |
| raster (full throw) | 19.50 | **17.37x** (geometry + depth) |
| rasternarrow | 34.09 | **9.93x** (geometry alone) |

So ~9.9x of the drift's disadvantage is cross-linking geometry and the rest
(~1.75) is depth. Two things that must travel with the number: "matched depth"
is **0.84x, not 1.00x** (the raster cannot be brought further -- coverage is
set by the 22.5 deg of RA drifting past in 1.5 h and by the beam, not the
throw, so even a sub-beam 2 deg half-width takes 377 px against the drift's
316), and **this run has its own common patch and HI realisation** because its
strategy list differs from `amp350`'s. Its drift floor is 338.58 against
`amp350`'s 399.61 and its drift/raster is 17.37x against 17.92x. Quote the
matched-depth comparison from inside this run only.

#### 8-strip row of tab:modes

Job 14102348, 2 min 41 s, `ska_hi_modecheck.py --strategies drift drift8`.
`drift8`'s amplitude constant had never been computed: **c = 0.848088**.

| strategy | area | npix | N_eff own | N_eff pin | /100 deg2 own | /100 deg2 pin |
|---|---|---|---|---|---|---|
| drift (3) | 265.2 | 316 | 38.78 | 38.78 | 14.62 | 14.62 |
| drift8 | 499.4 | 595 | 100.05 | 98.32 | 20.03 | **19.69** |

Use 19.69. Only -1.7% because `drift8` shares the reference's 0.5-beam spacing
and differs only in extent; spacing is what moves the constant. The 3-strip row
reproduces exactly, which is why it was in the list.

Results: `results/hi_{ablation,experiment}_rnarrow_f350_400_nc32_ns64.npz` and
`results/hi_modecheck_f350_400_ns64.npz`, all pulled from ilifu.

#### Figure changes (details in the zip's MANIFEST)

`hi_frequency_ratio` on the floor arm with four bands; `hi_cleaning_depth`
rewired and its duplicate legend entry deduped on label; **new**
`hi_noise_terms.png`; adopted-prior axvline on `hi_priorvar`; `hi_ladder` cut
to two panels; `hi_ladder_geometry` cut to three columns; `_f350` aliases for
`hi_kmodes` and `hi_sufficiency`.

Two numbers from those worth not re-deriving:

- **HI dims 16.08x from 350-400 to 850-900**, against only 4.52x to 675-725.
  Any residual/HI comparison ACROSS bands is partly a statement about the
  denominator.
- **The ladder's marginal return is flat**, 11.65-12.50 modes for 43.1-47.0
  deg2 per added strip from n=3 to n=20. The right panel of `hi_ladder` looks
  like diminishing returns only because of the per-100-deg2 normalisation.

`hi_ladder_geometry` was asked for at 0.5/1/2 beams but the sweep only measured
out to 1.5, so a 2-beam column would carry no measured density. It uses **1.5
beams**; a true 2-beam point is a short `ska_hi_ladder.py` job if wanted.

---

### Line-of-sight modes in the forecast, and the matched white level (2026-10-08)

Two tasks from `claude_code_nm_whitevar.md`, both pure post-processing of the
npz in `results/`. Full write-up and figures in `~/handover/ska_nm_kpar.zip`;
the numbers regenerate with `python ska_nm_kpar.py`.

#### The forecast was throwing away the line of sight

`sufficiency()` set `N_m = N_eff = tr(WA)`, the angular dof of ONE channel,
against a floor/Hi ratio taken as the MEDIAN over k_par bins. The cube has 32
channels and the forecast already assumes an ideal spectroscopic tracer, so it
was entitled to the radial modes too. Replaced by an inverse-variance sum,

    R_b = P_fl,b / (T_b^2 P_HI,b),
    Omega_req = SNR^2 / [rho sum_b n_b/(2 + R_b)],   rho = N_eff/Omega_sim

with `T` squared because `transfer_function` returns a cross-power over the
mock's AUTO-power, so T is an amplitude ratio. That retires the standing
"areas optimistic by 1/T^2" caveat -- the correction is now carried, not noted.

| 1 sigma, ideal tracer | published | per-k_par |
|---|---|---|
| 350-400, 3-strip drift | 1.09e6 | 36,648 |
| 350-400, 12-strip drift | 2.54e5 | 6,729 |
| 350-400, raster | 2,419 | **112** |
| 675-725, 3-strip drift | 655 | **12.9** |
| 675-725, 12-strip drift | 516 | **7.8** |
| 675-725, raster | 1,058 | **35.5** |

**Ordering unchanged at both bands**; the gaps narrow (drift/raster 451x ->
329x at 350-400). Everything at 675-725 MHz now sits BELOW its own simulated
patch, which predicts a 3.3-6.5 sigma cross-detection in data already on disk.
Not checked -- it was not in the brief -- but it is directly checkable and
should be before the number goes in the draft.

**Everything rests on bins 3-11.** They carry 87% of `sum_b n_b/(2 + R_b)`;
bins 1-2 carry 0.4% between them, because that is where the floor is worst
(`R_1 = 86,080` vs a median 621 for the raster) AND where `T_b` is lowest
(0.51, 0.81 vs ~0.88), and `T` enters squared. Convenient: the forecast is
insensitive to exactly the two bins whose `T` is least well measured.

**dof convention, the thing to get right.** `pk_par` subtracts the mean along
frequency and drops DC, so **31 real dof survive, not 32**: `n_b = 2` for bins
1-15, `n_16 = 1`. Counting one dof per bin instead halves every answer -- the
~2x ambiguity is real. `n_b` is right: it counts real dof, which is what
`tr(WA)` counts in the angular direction. **The earlier 12.8/15.9 for
`N_par,eff` could not be reproduced from any of four natural conventions, and
its drift/raster ordering is the reverse of all four.** If it matters, the
recipe has to come from `paper1_frozen_numbers.md`, which is not in this repo.

**Why summing over bins is allowed** (the referee question): the floor is
deterministic, so its bins are not independent draws -- but in CROSS the
variance is the chance correlation between a fixed floor and the tracer, and
the tracer supplies fresh randomness in every (angular mode, k_par) pair. The
objection bites only in the auto-spectrum, where the floor is already treated
as a bias.

**Auto-spectrum bound moved much further than the cross** -- 356,633 -> 1,895
deg^2 for the drift at 675-725, a factor 188 against `N_par,eff` = 70, because
the auto weight carries `(1 + R_b)^-2` and rewards good bins quadratically.
All three strategies now sit inside the reachable sky on that bound. Do not
lean on it: it still requires knowing the sky the survey failed to measure.

**T had to be borrowed across prior conventions.** The amplitude-only runs were
ablation-only, so no amplitude-only `T` exists. Measured cost at 350-400:
`T_b` agrees to 0.6% (12-strip) and 0.4% (raster) between shape-pinned and own
prior, against 5-9% mock scatter, and the drift (the reference) is identical to
1e-9. T is a property of the PCA clean, which never sees the prior.

**500-550 and 850-900 could not be forecast.** `rho` needs `N_eff` and the
patch area; neither band has an `hi_kmodes` npz, a modecheck or cached
operators (`hi_cache/` is 350 MHz only). Stopped rather than submit. Their
`T_b`/`R_b`/`N_par,eff` are in the report. 850-900 is worse off still: no
amplitude-only ablation, and the drift's `T` is unmeasurable above bin 10
(scatter 1.2-12.9 on 20 mocks, `T > 1` in bins 11-13).

#### Matched white noise: the draft's column was right by accident

The draft rescales the whole noise arm by `sqrt(3.2e-7/2.5e-6) = 0.35777`. Only
the WHITE part may be rescaled -- 1/f is a gain fluctuation and does not obey
the radiometer equation -- and `noiseonly_indep` is `(1 + g)(1 + w_ind)`, so it
carries both. Split (`white^2 = noise^2 - (1/f)^2`), rescale, recombine:
**all ten `tab:timetofloor` values reproduce to better than 0.02%**, because
1/f is 0.03-0.5% of the noise. Keep the split anyway: the column is correct by
accident, not by construction, and stops being correct if 1/f is ever raised.

`t_x = 3 h (noise/floor)^2`, the integration time to reach the floor, at 4
modes: 350-400 drift **2.5 h** (so the 3-strip drift is ALREADY floor-limited at
350 MHz in the 3 h simulated), 12-strip 13.4 h, raster 3,246 h; 675-725 drift
1,559 h, 12-strip 648 h, raster 851 h; 850-900 raster 5,668 h. The inversion is
the point: the better the floor, the longer it takes to integrate down to it.

`fig_noise_terms` relabelled: the dashed arm said "white noise (independent)"
but carries 1/f -- the same mis-pairing `fig_ablation` was corrected for on
2026-10-07. Now "all noise (white independent)".

#### Code

`ska_nm_kpar.py` is new and self-contained. `ska_hi_kmodes.sufficiency()` takes
an optional `kpar=` term dict and `load()` a `kpar=` flag, behind a `--kpar`
CLI switch; `ska_hi_design` gains `WHITE_VAR_MATCHED`/`noise_matched()`.
**The single-bin path is untouched** -- `--bands` without `--kpar` still passes
all ten `BANDS_EXPECT` guards and renders byte-identically, re-checked after
every patch.

Figures regenerated per-k_par, each with its single-bin twin kept beside it:
`hi_sufficiency_bands` / `_1bin`, `hi_sufficiency` (= `_f350`) / `_1bin`,
`hi_sufficiency_f700` / `_1bin`, `hi_noise_terms` / `_white2.5e-6`. The brief
named only the merged figure, but **`main.tex:751` includes
`hi_sufficiency.png`**, the single-band f350 panel -- so that one had to move
too or the figure would have contradicted its own caption. `hi_kmodes*.png` are
untouched: they decompose `tr(WA)` by angular scale and carry no `N_m`.

### The per-k_par forecast, checked against the simulated patch (2026-10-08)

`claude_code_forecast_check.md`. Post-processing only. Full write-up in
`forecast_check.md`; regenerate with `python ska_forecast_check.py`.

**The forecast holds, and tr(WA) is now measured, not assumed.**

The brief wanted the 20 Hi mocks re-cleaned against the noiseless floor arm.
**Neither is stored** -- `run_hi_experiment.py` keeps only `T_b` from the
mocks, and the ablation keeps power spectra, not cubes. Two substitutes on
disk do the same job:

1. `p_clean_cross_4` IS the wanted estimator. It is `P(clean(data), W A s_HI)`,
   the cleaned map crossed with the map-made true Hi, which for r = 1 with no
   shot noise is exactly the ideal tracer. So
   `A_b = p_clean_cross_b / (T_b * p_mapmade_b)` should be 1 in every bin and
   its scatter about 1 is the forecast's error bar. Stored for every run and
   every depth. **Remember this exists** -- it was built in 2026-08 and never
   used for anything.
2. `seed11`..`seed16` vary ONLY the Hi field (checked: the cleaned auto-power,
   99% floor, is identical across them to 0.2-0.8% per bin). With the base run
   that is **seven** independent draws of the floor-times-tracer chance
   correlation -- the quantity the variance is made of.

| 350-400 MHz | predicted SNR | mean A | pull | empirical SNR | ratio |
|---|---|---|---|---|---|
| 3-strip drift | 0.08 | -6.44 +/- 4.22 | -0.61 | -- (null) | -- |
| raster | 1.95 | 1.08 +/- 0.25 | +0.16 | 1.65 +/- 0.79 | **0.85** |

Pull = (A - 1) x predicted SNR. **All six band-strategy rows are within 0.7**,
so the forecast predicts its own error bar correctly six times out of six. At
675-725 the amplitude comes back 1.05 / 1.04 / 0.99 against predicted 1-sigma
errors of 0.29 / 0.17 / 0.37.

**tr(WA) measured.** Inverting the per-bin variance,
`N_ang = R_b / (n_b var(A_b))`, gives a mode count that owes nothing to
mode-counting theory:

| | tr(WA) | N_ang measured | ratio |
|---|---|---|---|
| 3-strip drift | 38.8 | 44.3 | **1.14** |
| raster | 118.5 | 135.6 | **1.14** |

1.14 for both, although the counts differ by 3x; rms(z) = 0.95 for each over
7 x 16 samples. The whole chain checks out at once: angular count, radial
count n_b, Blackman binning, and the (2 + R_b) variance. The taper's 1.78-bin
ENBW, flagged as a risk in `nm_kpar_results.md`, does not bite.

**Trap: use tr(WA), NOT rho * Omega_meas, for anything measured on the
interior pixels.** Spectra live on the interior of the common patch (119 px,
100 deg2 at 350-400; 171 px, 36 deg2 at 675-725), not on each strategy's own
patch, so the naive count is low by 3.0x (drift) and 6.6x (raster). The
interior cut throws away the patch EDGE, whose pixels have small lambda_i and
so carry almost none of tr(WA) = sum lambda/(lambda + sigma_p^-2). It costs
area but not modes. **No area in `nm_kpar_results.md` moves**: Omega_req uses
rho = N_eff/Omega_sim, which is right for tiling whole patches, and the
in-patch SNR is sqrt(N_eff sum_b n_b/(2+R_b)) either way.

**Two things not to misread.** The predicted SNRs in `forecast_check.md`
(1.95; 3.41/5.96/2.67) are NOT those in `nm_kpar_results.md` (2.27;
3.43/6.53/3.29): `R_b` there is measured from the same cleaned map as the
cross-power, so it is the own-prior TOTAL arm, not the amplitude-only FLOOR
arm. 10-20% shift, deliberate, so no prediction can fail for a bookkeeping
reason. And the single-realisation rows get **no** error-bar test: with one Hi
field the 16 bins are not 16 independent draws, because the floor is
spectrally low rank (participation ratio 2-5), so rms(z) over bins scatters
0.6-3.2 and a single-draw implied mode count is biased high by ~3.6x.

#### Figure and draft

`hi_sufficiency_bands.png` (per-k_par) re-rendered without the open
simulated-patch markers or their legend entry; the bar runs from 1 sigma
outwards. `fig_sufficiency_bands` gained `patch_markers=`, passed False by
`--kpar`; the single-bin path still renders byte-identically.

**The repo copy of `main.tex` is a 2026-10-01 snapshot** (last commit
`d12b243` 2026-09-30, mtime 2026-10-01 17:41, uncommitted since). It has none
of the 2026-10-08 open-items box, no `hi_sufficiency_bands`, no per-k_par
equation, and still includes `hi_sufficiency.png` at line 751. It is seven
days behind the delivered draft, not a divergent branch -- overwrite it from
the delivered copy. Because it still named them, the single-band
`hi_sufficiency.png` / `_f700.png` were regenerated per-k_par on 2026-10-08;
the delivered draft uses the merged figure, so those two are spare but
consistent.

### Why the raster "loses" at 675-725 MHz: it does not (2026-10-08)

`claude_code_raster_sampling.md`. Diagnostic only -- stored pointings, cached
operators, result npz. Write-up in `raster_sampling.md`, regenerate with
`python ska_raster_sampling.py` (~5 min).

**The premise was wrong.** The raster's post-clean foreground residual -- the
INSTRUMENT half, with the Hi's own dimming divided out -- improves **2.5x**
from 350-400 to 675-725. The three-strip drift's improves **105x**. The
inversion is the drift catching up, not the raster degrading. (That table was
already in HANDOFF under "Both strategies improve monotonically with
frequency"; this is the first time the question was asked of it.)

All four sampling candidates ruled out, measured on the interior pixels:

| | 350-400 | 675-725 |
|---|---|---|
| raster sweep spacing | 0.769 deg = **0.19 beams** | 0.769 deg = **0.37 beams** |
| largest gap (any interior px to nearest sample) | 0.08 beams | 0.10 beams |
| raster two-direction coverage | **100%** | **100%** |
| raster samples/beam (vs 3-strip drift) | 503 (0.68x) | 161 (0.56x) |
| raster patch sky / interior sky | 1.048 | 1.068 |

- **(i) gappy sweeps: ruled out.** The sweep separation is geometry, so it is
  0.769 deg at BOTH bands; only the beam changes. Tighter than the 0.5-0.6
  beams the brief's envelope expected. **"6'/s" in the code is an AZIMUTH
  rate** -- `sweep_s = 150` s is one half-sweep of the full 15 deg throw, so
  0.1 deg/s in az, 4.7'/s on sky after cos(el = 38.19).
- **(ii) single-direction interior pixels: ruled out.** 100% two-direction in
  both bands. Raster passes sit at PA 126 deg (rising) and 54 deg (setting),
  72 deg apart; all drift strips at PA 89.95 with 0.01 deg scatter, i.e. 100%
  single-direction by construction.
- **(iii) sampling density: right sign, 20x too small.** raster/drift
  samples/beam moves 0.68 -> 0.56 (factor 1.2) while floor/Hi moves 0.056x ->
  1.64x (factor 29). And at 350-400 the raster is ALREADY the sparsest
  strategy and wins by 17.9x. **Prior dominance runs the wrong way**: mean
  diag(WA) = tr(WA)/npix exactly, so no operator is needed -- raster 0.173 vs
  drift 0.123 at 350-400, 0.146 vs 0.088 at 675-725. The raster is the LESS
  prior-dominated strategy in both bands. With the cached 350-400 operators
  the interior can be done properly and says it harder: 0.328 (raster) vs
  0.204 (drift), both ~1.7x their patch averages.
- **(iv) brighter footprint: ruled out as the cause.** 7% excess at 675-725,
  and `drift12m` carries 5.9% while beating the raster by 1.4x.

**What actually changed: the drift's ladder spacing.** Its fixed 2 deg step is
0.50 beams at 350.78 MHz and 1.00 beams at the 675-725 centre -- heavily
redundant to contiguous and non-redundant. Same mechanism the draft already
uses for drift vs drift12m at 700 (1.00 vs 0.60) and for the drift's rise at
850-900 (1.25 beams, gapped). The drift-raster inversion near 600 MHz is the
third instance, not a separate effect. Three points, spacing and density not
varied independently, so still a consistency argument and not a measurement of
an optimum.

**Two reusable things.**

- **Per-pixel floor maps can be rebuilt without map-making**: PCA the stored
  `data_cube` at common resolution, take the per-pixel rms over channels.
  Validated -- reproduces `cr_p_clean_auto_4` to 1 part in 1e6. Medians
  (real-space, all k_par, NOT the headline numbers, but they reproduce the
  ordering and the inversion): 350-400 drift 437 / drift12 250 / raster 66;
  675-725 drift 10.7 / drift12m 15.6 / raster 21.1. The raster's correlation
  with exposure flips sign between bands (-0.24 -> +0.62) -- do NOT read that
  as a cause: samples/beam spans only 156-166 across its 675-725 interior
  while the per-pixel floor spans 13.4-34.6.
- **Footprints can be reconstructed from the boresight track alone** (stack a
  Gaussian of the band FWHM, keep the brightest npix pixels). Checked against
  the two cached 350-400 operators: **94%** of the true pixel set. Enough for
  brightness statistics at 675-725, where no operator is cached.

**Still missing: interior diag(WA) at 675-725.** `hi_cache/` holds 350-400
nside 64 only. One channel per strategy at nside 128 is a ~30-60 min local
operator build; not run, since the brief excluded operator construction.

## Queued next

**DONE 2026-08-18: experiment 006 is plotted.** `ska_hi_plots.py` draws four
figures into `figures/` (gitignored, regenerate with
`/home/geoff/gibbs_venv_312/bin/python ska_hi_plots.py`, seconds — it only reads
cached results):

| figure | what it shows |
|---|---|
| `hi_transfer_function.png` | T(k_par), drift vs raster, one panel per mode count, ±1σ over 20 mocks. The cross-linking result. |
| `hi_residual_vs_modes.png` | residual/HI vs modes removed, log scale, k-range bars — never approaches 1. |
| `hi_ablation.png` | floor vs total vs **two** noise arms. Full-data rings sit *on* the floor line. Common resolution + interior pixels since 2026-08-26 (`fig_ablation(..., tag="cr_")`). **The "noise is 10-50x below" reading was wrong and was fixed 2026-10-07 — see "The ablation noise model" in the 2026-10-07 section above.** |
| `hi_frequency_structure.png` | radial power surviving the map-maker, plus the frequency-structure split that explains why the drift fails. |
| `hi_patch_maps.png` | the GDSM field and the nested drift/raster/interior footprints. |
| `hi_maps.png` | map domain: true HI, drift map-made, raster map-made, and the cleaned data on a 50x wider scale. |

Two notes for whoever revises them. A single-line-of-sight panel was tried for
the fourth figure and **abandoned** — the traces are too noisy to read "flat in
frequency" off by eye; the two rms ratios state it far better. The map-domain
panel was skipped at first on the grounds that the cleaned map is pure residual;
**that judgement was wrong and it was added 2026-08-19** — putting the cleaned
panel beside the others on its own 50x colour scale is the clearest single
statement of the result in the set.

Map colour rules, since healpy's defaults break them: signed fields (HI,
residuals) use the diverging blue-grey-red pair so zero recedes and sign reads;
positive-definite temperature uses the one-hue sequential ramp; footprint
membership uses categorical slots. Never a rainbow — it would encode magnitude
as hue and invent structure.

Colours are the validated categorical slots (blue drift, orange raster, aqua for
the third ablation arm) — those three clear all-pairs CVD separation, and every
series is legended or directly labelled so identity never rests on hue.

The rest of the queue is the *improvement* work.

0. **Raise the number of measured modes** — this is now the headline lever,
   replacing the nside re-run. The drift measures 17, the scan 34; everything
   else is prior. More crossing angles and more elevations both add modes,
   whereas more integration time on the same tracks does not. **Experiment 005
   adds a third lever: frequency** — 55 modes at 1050 MHz on the same tracks,
   at the price of a 2.3× smaller patch.
0b. **Follow-ups from experiment 006 (HI).** The negative result is bounded by
   mode count, so the levers are the same ones: (a) repeat at 675-725 MHz,
   where k_perp reaches 0.053 rather than 0.016 so T(k) becomes a genuine 2D
   measurement and the null space is less pathological — needs nside 128 and
   the elevation ladder is marginal there (2 deg steps = 1.0 FWHM); (b) push
   the floor down with more crossing angles and check whether the residual/HI
   ratio follows; (c) the 006 raster is not depth-matched, which does not
   affect T(k) but does affect the residual comparison — fold in item 1 below.

0c. **Would a different cleaning basis help? Measured 2026-08-25: no.**
   `ska_hi_rank.py`. Modes of the channel-channel covariance needed to carry
   90/99/99.9% of the variance:

   | | 90% | 99% | 99.9% | lambda_10/lambda_1 |
   |---|---|---|---|---|
   | GDSM foreground | 1 | 1 | 1 | 2.6e-17 |
   | floor, drift | 2 | 6 | 10 | 6.3e-4 |
   | floor, raster | 1 | 2 | 5 | 5.7e-5 |

   Two things follow. **The foreground is rank 1** over this 50 MHz block — one
   PCA mode removes essentially all of it, and the ablation agrees (after 1 mode
   the drift residual is already 896x the HI). So from mode 2 onward the filter
   is not cleaning foregrounds at all, it is fighting the instrumental floor.
   And **the floor is also low-rank**, so a subspace method *can* capture it.

   ICA, GMCA, NMF and kernel PCA all remove a rank-N subspace and differ only in
   how they pick it — those differences matter for a mixture of astrophysical
   components with different statistics, which is not what this is. Polynomial
   modes should be *worse*: they impose smoothness the floor does not have,
   where PCA at least adapts. Cheap to falsify — one line in
   `ska_hi_analysis.pca_clean`, since fastbox ships `ica_filter`, `nmf_filter`,
   `kernel_pca_filter` and `gpr_filter` — so worth an afternoon as a null
   result, not a research direction.

   **Correction.** An earlier version of this item put GPR first, on the grounds
   that the blocker was a *non-smooth* floor PCA could not represent. The rank
   measurement falsifies that: the floor is low-rank enough for PCA. The blocker
   is that the floor's subspace **overlaps the HI's** — both are the smooth
   low-k_par modes — so removing 6-10 modes costs 83-94% of the HI at low k
   (T = 0.056-0.17 at 10 modes). No basis escapes an overlap.

   What follows instead, in order:

   (i) **Model the floor rather than filter it.** The floor is `(I - WA)s`, and
   in simulation `W` and `A` are known exactly. Compute its covariance and fold
   it into the noise model instead of blindly projecting. Blind separation
   throws away the one real advantage a simulation study has — that the
   instrument is known.

   (ii) **Shrink the floor with geometry.** The raster's floor is already rank 2
   against the drift's 6: more measured modes give a more compressible floor.
   The same lever as item 0, now visible in the rank.

   (iii) GPR is still worth trying, but it does not deserve top billing.

0d. **DONE 2026-08-26: `hi_ablation` now runs at common resolution.**
   `ska_hi_ablation.py` computes the interior mask and reconvolves each arm,
   and stores all three conventions the way `run_hi_experiment.py` does --
   `""` (as run, 277 px), `int_` (as run, 119 interior px) and `cr_` (common
   resolution, 119 px). `fig_ablation` defaults to `tag="cr_"`, so figure 3 and
   figure 2 are now on the same footing; it falls back to the as-run keys with
   a warning if handed an older `.npz`.

   Two checks worth keeping. The as-run column reproduces the previous run
   exactly (drift 572.4x, raster 66.0x at 4 modes), so the solves are
   untouched and only the post-processing changed. And the `cr_` "total" arm
   agrees with the residual/HI table it should match: drift 347.3x vs 347.3x,
   raster 24.0x vs 24.0x.

   What the re-run changed in the *wording*, not the conclusion: the noise
   floor separation is 52x (drift) / 10x (raster) in amplitude rather than
   "~38x", and floor-vs-total agreement is "a few per cent", not 0.2%. The
   ANALYSIS_LOG entry and `paper/main.tex` still carry the old phrasing and
   were deliberately left alone -- the paper is being written by hand from
   here.

0e. **DONE 2026-08-26: the 1/f results are k-resolved and bounded.** Three
   pieces, all on the common-resolution interior pixels at 4 modes removed, so
   they sit on the same footing as Figures 1-3.

   *(i) The ablation splits the noise.* `ska_hi_ablation.py` now carries
   `gainonly` and `whiteonly` alongside `noiseonly`. The noise is
   multiplicative, so the split is exact bar the `g*w` cross term; the run
   prints the check and it closes to 1.9e-4 (drift) / 2.3e-4 (raster) of the
   noise amplitude.

   | 4 modes, common resolution | floor | total | 1/f | white |
   |---|---|---|---|---|
   | drift | 336.3x | 349.4x | **4.0x** | 5.2x |
   | raster | 22.7x | 24.0x | **0.4x** | 2.1x |

   For the raster 1/f alone sits *below the HI*. That is much stronger than the
   old "5-8% of the residual" and it is in the paper's own units.

   *(ii) How much worse would 1/f have to be?* `ska_hi_1f_scan.py` (default
   `--mode knee`) scans the knee -- where 1/f gain power crosses white -- which
   is at **0.95 mHz** for `GAIN_PARAMS`. Note where that sits: the drift's
   fundamental (0.28 mHz) is *below* the knee, the raster's sweep (3.33 mHz)
   *above* it. Knee at which 1/f reaches the floor:

   | | alpha 1.5 | alpha 2.0 | alpha 2.5 |
   |---|---|---|---|
   | drift | 316 mHz (333x) | 80 mHz (84x) | 33 mHz (35x) |
   | raster | 184 mHz (194x) | 57 mHz (60x) | 21 mHz (22x) |

   **Quote the raster at alpha 2.5: 22x, the hardest case to argue with.**

   *Do not present the knee axis as six measurements.* Scaling one realisation
   scales the whole pipeline linearly -- `solve` is affine and differences out,
   and `pca_clean` is scale-invariant since scaling a cube scales its
   covariance and leaves the eigenvectors alone -- so
   `ratio ~ knee**(alpha/2)` exactly, verified against the grid to 2e-4. The
   figure draws lines for the scaling and a filled marker only at the measured
   anchor. The genuinely measured content is three points per panel.

   *(iii) Does the cut-off matter?* `--mode cutoff` scans `w_c`, which changes
   the correlation *shape*, so nothing factors out and every point is a fresh
   draw. Fiducial `w_c` = 0.175 mHz, a 1.6 h timescale -- about one pass.
   Extending the spectrum down to 16 h grows the gain rms 2.6x and moves the
   residual 12% (drift 3.98 -> 4.44x). **Slow gain drift is degenerate with the
   foreground**: it multiplies a rank-1 component, so the cleaning removes it
   too. That is the answer to "you assumed a convenient cut-off". Above the
   fiducial the residual falls, but only because raising `w_c` deletes 1/f
   power (rms drops 6x) -- do not read that as a result.

   *Known soft spot.* The two lowest-`w_c` raster points are non-monotonic in
   rms (3.28e-4 then 3.80e-4) where the drift's are clean: a 3 h observation
   barely samples a 16 h correlation time, so a single draw of the lowest modes
   is not representative. Averaging a few seeds there (~5 min) would give error
   bars; the flat trend itself is not in doubt.

   *Still open from the original note:* scatter over seeds generally, and a
   sentence on what frequency structure the gain model does and does not have
   (the draw is common to all channels by construction of the replay).

   *Two traps when writing it up.* The **+36% high-pass penalty is the
   galactic-plane drift**; off-plane at nside 64 the same 2 mHz filter helped
   slightly (-9.7%), and HP numbers should not be quoted for the raster at all.
   And `beam_crossing_freq_hz` bounds the *support* of the sky signal, not
   where its power lives -- do not compare it against a filter cutoff.

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

5. **Retire the word "published" from anything a reader sees.** Nothing in
   this series is published; the term entered on 2026-08-13 (`6828e1b`) meaning
   "the experiment-006 baseline the flat-prior audit had to reproduce", and
   spread to mean "the canonical run, as opposed to a sweep variant". Fine in
   code, wrong in a figure. Only one live instance reaches the paper:
   `ska_hi_design.py:1251`, the annotation `"$\star$ = the published run"` in
   `fig_realisations`, explaining the stars at lines 1250/1271 (realisation
   index 0 — the seed used in experiment 006). Replace with "fiducial run" and
   regenerate `hi_realisations.png` + its notebook cell. Two related points:
   `ska_hi_kmodes.py:175` also plots a star but it marks the *required* sky
   area, not a fiducial run, so the same glyph means two different things
   across `hi_kmodes` and `hi_realisations` — worth separating. And the
   commented-out annotations (`ska_hi_design.py` 517/571/655/784, and the
   `_title` blocks at design 593 and kmodes 199) all carry the word, so fix
   them before any of them is uncommented. Docstring/comment uses
   (`ska_hi_plots.py` 336/613, design 78/449/484/494, kmodes 19) can stay.

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
