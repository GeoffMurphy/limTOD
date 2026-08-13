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

**1. Audit the cross-linking headline against a flat prior.** The previous entry
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
