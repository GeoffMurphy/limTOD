# SKA drift-scan series — analysis log

A running record of **why each step was taken and the maths behind it**, in the
order it happened. Written to be re-readable months later without reconstructing
the reasoning from code.

## How this relates to the other three documents

| file | answers | shape |
|---|---|---|
| `README.md` | what is in this directory | index |
| `METHODS.md` | how the machinery works | atemporal reference |
| `HANDOFF.md` | where things stand, what the numbers are | current state |
| **`ANALYSIS_LOG.md`** | **why we did it and how the maths goes** | **chronological** |

Rule of thumb: a *result* goes in `HANDOFF.md`, a *mechanism* goes in
`METHODS.md`, and the *derivation and the argument* go here. Numbers are
repeated here only where the derivation needs them — `HANDOFF.md` stays the
single source of truth for values.

## Entry format

Each entry: **the question**, **the method** (with the maths written out),
**what was computed**, **the result**, and **what it changed**. Dated, newest
last.

---

## Notation

Fixed throughout, matching the code in `ska_freq_sweep.py` and
`plot_freqsweep_maps.py`.

| symbol | meaning | in code |
|---|---|---|
| $s$ | true sky on the solved pixels | `truth` |
| $\hat{x}$ | reconstructed map | `est[kind]` |
| $\mu$ | Wiener prior mean | `Tsky_prior_mean` |
| $A$ | sky → TOD projection operator | built in `HPW_mapmaking.__init__` |
| $N$ | TOD noise covariance | diagonal |
| $S$ | prior covariance on the sky | `Tsky_prior_inv_cov_diag` (diagonal) |
| $\sigma_s$ | sky *structure* rms, `np.std(s)` | `sky` |

The map-maker solves

$$(A^{\top}N^{-1}A + S^{-1} + \lambda I)\,\hat{x} = A^{\top}N^{-1}d + S^{-1}\mu$$

so $\hat{x}$ is a precision-weighted blend of what the data measures and what
the prior asserts. Every claim about "what the survey measures" has to survive
the second term — which is the subject of the first entry below.

The headline metric everywhere is the **residual ratio**

$$r \;=\; \frac{\operatorname{std}(\hat{x} - s)}{\sigma_s}
      \;=\; \frac{\text{residual rms}}{\text{sky structure rms}}$$

normalised because the sky itself dims $\sim\nu^{-2.7}$ (a factor ~24 across
Band 1), so absolute residuals are not comparable between channels.

---

## 2026-08-07 — Is the frequency trend the survey, or the prior?

**Question.** Experiment 005 reported $r$ improving monotonically across Band 1
(0.525 → 0.336, adequately sampled). Raised by a reviewer: is that the *survey*
reconstructing better at high frequency, or just the prior getting sharper?

**Why it was a fair challenge.** The series' standard recipe sets the prior mean
to the beam-smoothed true sky,

$$\mu = B_\nu \, s$$

where $B_\nu$ is that channel's beam. The beam narrows as $\nu^{-1}$, so
$B_\nu \to I$ with increasing frequency and $\mu \to s$. **The prior is not held
fixed across the sweep — it sharpens on its own, with no data involved.** Any
frequency trend in $r$ is therefore confounded by construction.

**Method.** Separate the two contributions by evaluating the prior with the data
switched off, and by re-solving under a prior that cannot improve with
frequency.

*Prior alone.* Setting the data term to zero gives $\hat{x} = \mu$, so

$$r_{\text{prior}} = \frac{\operatorname{std}(\mu - s)}{\sigma_s}$$

- beam-smoothed prior: $r_{\text{prior}} = \operatorname{std}((B_\nu - I)s)/\sigma_s$, which falls as the beam narrows.
- flat prior ($\mu = \langle s \rangle$, a constant): $\mu - s = -(s - \langle s\rangle)$, so $r_{\text{prior}} \equiv 1$ **exactly, at every frequency** — which is precisely why it is the right control.

*What the data bought.* The gain of the full reconstruction over its own prior:

$$g = \frac{r_{\text{prior}}}{r_{\text{recon}}}$$

$g > 1$ means the data improved on the prior; $g < 1$ means it made the map
worse. Under the flat prior $r_{\text{prior}} = 1$, so $g = 1/r_{\text{recon}}$
and the metric reads directly off the reconstruction.

**Computed.** `plot_freqsweep_maps.py::solve_channel` re-solves all five
channels for both prior choices, holding geometry, seeds, noise and grid fixed —
only $\mu$ changes. Each channel on a grid that samples its beam adequately
(nside 128 up to 700 MHz, nside 256 above; see the sampling entry in
`METHODS.md` §9).

**Result.** See `HANDOFF.md` for the full table. The shape of it:

- beam-smoothed prior alone already scores 0.511 → 0.360 across the band
- full reconstruction with that prior scores 0.525 → 0.336
- so $g = 0.97$–$1.16$, and **at 350 MHz $g < 1$ — the data makes the map
  marginally worse**, injecting noise into the modes it measures
- under the flat prior, $r$ = 0.769 → 0.597 and $g$ = 1.30 → 1.67, saturating
  above 875 MHz

**What it changed.** The conclusion survived but the evidence did not: the
0.525 → 0.336 trend is mostly the prior sharpening. Quote **flat-prior numbers
or the mode count** — both prior-independent — whenever the claim is about what
the survey measures. Two map figures are emitted rather than one, plus
`figures/freqsweep_prior_dependence.png`.

**Caveat carried forward.** This is not specific to experiment 005: *every*
residual in `HANDOFF.md` is quoted with the truth-derived prior. What is
specific to 005 is that a **frequency** trend is confounded, because the prior
moves with the beam. A comparison at fixed frequency (e.g. cross-linking) is not
confounded in the same way — but that has not been demonstrated, only argued.
See the next entry.

---

## 2026-08-12 — Position, and what limits a drift scan

**Where the analysis actually stands.** Everything measured so far is
*reconstruction fidelity of the diffuse sky*: a known GDSM sky in, a map out,
residual against truth. The power spectra in the notebooks are `hp.anafast` of
the **residual**, i.e. a diagnostic of reconstruction error — not an HI
measurement. There is no 21 cm signal anywhere in the simulation and no
foreground separation step.

So "a drift scan looks reasonable" currently means *the map-maker recovers the
bright foreground sky reasonably*. Necessary, but a long way from sufficient.

**The limiter already measured — and it is not sensitivity.** Eigen-analysis of
$A^{\top}N^{-1}A$ gives 17 modes above 1% of $\lambda_{\max}$ for the drift
against 34 for the raster (grid-independent: verified at nside 32/64/128/256).
Everything outside that span is prior. The integration table follows from it:
3 h → ∞ moves the total only 21.9 → 15.7 K, because the asymptote *is* the beam
+ prior floor. Eliminating 1/f entirely buys 1.6%.

**A drift scan is therefore information-limited, not noise-limited, and the
information is set by geometry.** This is the cleanest feasibility statement the
series currently supports.

### Queued, in the order I would do them

**1. Audit the cross-linking headline against a flat prior.** — **done, see the
next entry.** The previous entry
established that the truth-derived prior did the work in the frequency trend,
and `HANDOFF.md` notes every residual in the file uses that prior — including
the −26.5% floor and −44% beam-scale cross-linking gains. The argument that
cross-linking is *not* confounded (both arms sit at the same frequency, so the
prior does not move between them) is sound but untested. Operators are cached
(`mapmaker_ops_ska_meerklass_offplane_gauss_ns64.pkl`), so this is a re-solve —
minutes, and `plot_freqsweep_maps.py` already has the flat-prior machinery to
copy. Cheap, and it is the obvious question once the prior-dependence figure is
circulating.

**2. HI plus foreground cleaning.** The real feasibility question, and
experiment 005 is what makes it urgent rather than generic: the measured mode
set changes drastically across Band 1 (17 → 56 modes, patch shrinking 2.3×, the
elevation ladder breaking into ribbons above 700 MHz). **The map's null space is
strongly chromatic.** Foreground cleaning assumes foregrounds are spectrally
smooth *after* the instrument and pipeline; a frequency-dependent null space
imprints spectral structure on the foreground residual that PCA/SVD cannot
remove, so it eats HI instead. The sweep that produced the good news is exactly
the setup for measuring the bad news.

Design: inject an HI cube, solve every channel on a **common grid and common
patch** (per-channel grids are themselves a chromatic effect and would confound
the test), PCA-clean, and measure the HI transfer function against the number of
modes removed. Drift vs raster.

**3. Sidelobe pickup of the galactic plane.** Drift-specific and unavoidable —
a parked dish cannot schedule around the plane transiting its sidelobes.
`ska_tapered_beam` exists and was used on the plane, but the off-plane sweep is
Gaussian only.

**4. Polarisation leakage.** Absent entirely. limTOD accepts `(3, npix)` I,Q,U
beams. Polarised synchrotron plus Faraday rotation is the classic IM killer
precisely because it manufactures spectral structure out of a smooth foreground.

**5. Calibration and bandpass** — worth doing as a *positive*: fixed elevation
gives a very stable ground/spillover signature, a genuine drift-scan advantage
not yet claimed.

Also outstanding from `HANDOFF.md`: the 700 MHz grid check, the properly
depth-matched raster, more crossing angles, and the two limTOD author reports
(the `nside_target != nside_beam` operator bug is real and will bite others).

---

## 2026-08-13 — Does cross-linking survive a flat prior?

**Question.** The 2026-08-07 entry showed the beam-smoothed-truth prior supplied
most of the apparent frequency trend, and `HANDOFF.md` records that *every*
residual in the series uses that prior — including the cross-linking headline
(floor −26.5% on the plane, −37.8%/−44.0% off-plane at beam scale). The argument
that a fixed-frequency comparison is not confounded the same way, because the
prior does not move between the two arms, was made but never demonstrated.

**Method.** `audit_flat_prior.py`. Re-solve the off-plane drift-vs-raster
comparison under both prior means, changing nothing else — same operators, same
cached TODs, same seeded noise replay, same prior covariance
$S^{-1} = I/\sigma_s^2$, same regularisation. Only $\mu$ differs:

- smoothed: $\mu = B_\nu s$, the series' recipe
- flat: $\mu = \langle s \rangle$, constant at the patch mean

Both strategies restricted to the 317 shared pixels; beam-scale numbers use
cells lying entirely inside the shared region, so the two are compared on
identical cells.

**Validation first.** The smoothed-prior arm reproduces the published numbers to
three decimals (floor 1.287 → 1.158 K, −10.0% per-pixel; 0.300 → 0.187 K,
−37.8% beam-scale; total at matched depth −4.5% and −44.0%). The re-solve is
therefore sound and any difference below is the prior, not the pipeline.

**Result — the headline survives, and strengthens.**

| | per-pixel | beam-scale |
|---|---|---|
| smoothed, floor | 1.287 → 1.158 K (**−10.0%**) | 0.300 → 0.187 K (**−37.8%**) |
| smoothed, total @ depth | 1.300 → 1.241 K (−4.5%) | 0.313 → 0.175 K (−44.0%) |
| **flat**, floor | 1.806 → 1.218 K (**−32.6%**) | 0.392 → 0.210 K (**−46.4%**) |
| **flat**, total @ depth | 1.819 → 1.294 K (−28.9%) | 0.402 → 0.196 K (−51.3%) |

This is the **opposite** of what happened to the frequency trend. There the
prior manufactured the effect; here it was *masking* it. Removing the prior's
structure roughly triples the per-pixel cross-linking gain (−10.0% → −32.6%).

**Why — the drift leans on the prior far harder than the raster does.** Taking
the ratio flat/smoothed on the same quantity (1.00 would mean the prior was
doing nothing):

| | drift | raster |
|---|---|---|
| per-pixel floor | 1.287 → 1.806 K (**1.40×**) | 1.158 → 1.218 K (**1.05×**) |
| beam-scale floor | 0.300 → 0.392 K (1.30×) | 0.187 → 0.210 K (1.12×) |

Withdrawing the prior costs the drift 40% and the raster 5%. That is exactly
what the mode-count argument predicts: the drift measures 17 modes to the
raster's 34, so it has twice the null space for the prior to fill, and a
truth-derived prior fills it flatteringly. Prior-independently, on the shared
pixels the floor is 0.786 × sky structure for the drift against 0.530 for the
raster.

**What it changed.** The cross-linking conclusion is now demonstrated rather
than argued, and it is stronger than published. More importantly the *reason*
to prefer cross-linking is sharper: it is not merely that the raster's residual
is lower, it is that the drift's number was being propped up by a prior no real
survey has. Quote the flat-prior gains when the claim is about what the survey
measures.

**Loose end, pre-existing.** At beam scale the raster's "floor" (0.187 K) exceeds
its "total" (0.181 K) under the smoothed prior, and likewise under the flat one
(0.210 vs 0.199 K) — adding noise slightly *reduces* the beam-scale residual, so
the quadrature decomposition into floor and noise that holds on the plane to
<0.1% does not hold here. This is in the published numbers too, not introduced
by this audit. Most likely the beam-scale residual is bias-dominated and the
noise partially decorrelates it. Worth understanding before quoting a beam-scale
noise term off this field.
