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

---

## 2026-08-18 — Experiment 006: does any HI survive the pipeline?

**Question.** Every residual in this series so far measures reconstruction
fidelity of the diffuse *foreground* sky. There is no 21 cm signal anywhere in
experiments 001-005 and no foreground separation step, so "a drift scan looks
reasonable" has only ever meant "the map-maker recovers the bright foreground
sky reasonably". This is the queued test that asks the real question.

The specific worry, raised by experiment 005: the map's null space is strongly
chromatic (17 measured modes at 350 MHz against 56 at 1050 MHz, patch shrinking
2.3x). Foreground cleaning assumes foregrounds are spectrally smooth *after*
the instrument and pipeline. A frequency-dependent null space imprints spectral
structure on the foreground residual that PCA cannot remove, so it eats HI
instead.

**Setup.** 350-400 MHz, 32 contiguous channels, nside 64 — the series anchor, so
the established off-plane drift-vs-raster geometry carries over. New simulation
work rather than a re-solve: experiment 005's five channels span the whole of
Band 1, which is useless for PCA. One operator and one foreground TOD per
channel per strategy, 64 of each.

HI from fastbox (`ska_hi_mock.py`): Gaussian density -> HI bias -> log-normal ->
linear RSD -> Tb(z), giving Tb = 0.518 mK and b_HI = 1.559 at z_eff = 2.80. The
box-to-HEALPix bridge is a 3D lightcone interpolation, each (pixel, channel)
sample placed at its true comoving position, not a tangent-plane projection
(which would stretch scales 4.3% at the raster patch edge).

**What this band can measure.** At z ~ 2.8 a 15 m dish resolves 401 Mpc
transverse, so k_perp reaches only 0.0157 Mpc^-1 while k_par runs 0.012 to
0.192. The accessible 3D k-space is a sliver near the k_par axis, and the
transfer function is quoted against k_par for that reason. This is a property
of the survey, not of the simulation.

**Two shortcuts, both validated.** Foregrounds are simulated properly by
`TODSim` (that is what produces the beam + prior floor), but HI is injected as
`A s` through the map-maker's own forward operator `mm.Tsys_operators`, which
reproduces a simulated sky TOD to 4-7% rms (corr 0.998) — the gap being the
sub-beam structure it cannot represent. This is what makes a mock-averaged
transfer function affordable: a solve costs 0.57 s against ~5 min for a TOD
simulation. The Wiener filter is affine in the data, so
`solve(d + As) - solve(d) = response(s)` exactly; measured at 5e-6. Prior is
flat, following the 2026-08-13 audit.

Common patch is the intersection across all channels and both strategies: 277
px, set by the drift's narrowest channel and entirely inside the raster's. The
chromatic shrinkage is visible directly — drift 316 px at 350 MHz falling to
277 px at 399 MHz, raster 684 -> 636.

### Result 1: cross-linking substantially improves HI recovery

T(k_par), the fraction of *map-made* HI surviving the clean, 20 mocks:

| k_par | 0.012 | 0.024 | 0.036 | 0.048 | 0.060 | plateau |
|---|---|---|---|---|---|---|
| drift, 4 modes | 0.228 | 0.469 | 0.732 | 0.854 | 0.880 | ~0.90 |
| raster, 4 modes | 0.449 | 0.822 | 0.965 | 0.961 | 0.956 | ~0.96 |
| drift, 8 modes | 0.095 | 0.138 | 0.256 | 0.451 | 0.554 | ~0.80 |
| raster, 8 modes | 0.034 | 0.142 | 0.437 | 0.769 | 0.899 | ~0.94 |

Loss concentrates at low k_par, as expected — those modes look most like
foregrounds. The raster is better nearly everywhere and the gap widens with
aggressive cleaning. **T(k) is depth-independent**, being each arm's ratio
against its own injected response, so this comparison is fair even though the
raster is ~2.3x shallower on the shared patch.

### Result 2: the drift map-maker turns HI into something that looks like foreground

The headline rms ratio is a trap. Split by frequency structure:

| drift | true -> map-made | ratio |
|---|---|---|
| full rms | 0.2257 -> 0.2085 mK | 0.924 |
| frequency-mean removed | 0.2173 -> 0.0825 mK | **0.380** |
| the frequency-mean itself | 0.0609 -> 0.1915 mK | **3.14** |

The drift *amplifies* the frequency-coherent part of the HI 3.1x while keeping
38% of the frequency-varying part. With ~17 measured modes it projects HI onto
a nearly frequency-constant subspace — so the surviving signal looks like
foreground and PCA removes it. That is the chromatic null space biting the
signal directly, and it explains the collapse of drift T(k) at low k. The
raster gives 1.11x and 0.532, far healthier.

### Result 3: no HI is recoverable, and the blocker is the floor, not noise

After cleaning, the residual sits far above the HI it contains (amplitude
ratio, median over k):

| | 4 modes | 8 modes |
|---|---|---|
| drift | **572x** | 152x |
| raster | **65x** | 46x |

Ablation settles what the residual is. Solving with noiseless data ("floor")
against the full data ("total") gives 326900 vs 327600 (drift, 4 modes) and
4262 vs 4357 (raster) in power — **identical to 0.2%**. The noise alone
contributes 15x HI (drift) and 5.8x (raster), some 38x below the floor. The
blocker is therefore entirely the **beam + prior floor**: sky the strategy
never measured, filled by the prior, whose frequency structure is set by the
chromatic null space. PCA cannot remove it because it is not spectrally smooth.
This is the concern of the 2026-08-12 entry, confirmed and quantified.

**Consequence for the deliverable.** The transfer function is well measured —
injection is differential and mock-averaged, so the floor does not enter it —
but there is nothing in the data to correct. The cross-power estimator meant to
isolate surviving HI has no signal in it and wanders between -44 and +0.6 times
the map-made HI power, which is why `p_corrected_*` in the results file is
negative and unstable. **It must not be read as a recovered HI spectrum**;
`residual_over_hi_*` is recorded alongside it to make the reason explicit.

**What it changed.** The series can now say something about HI rather than only
about foreground reconstruction, and it is a negative result at this band:
a 15 m dish at z ~ 2.8, on 3 h of data, leaves a foreground-cleaning residual
50-500x the HI, dominated by the beam + prior floor. Cross-linking improves
every part of this — 9x lower residual, better T(k) at every k, and a
map-making response that is not degenerate with foregrounds — which is now the
strongest argument in the series for slewing in azimuth. Whether the residual
can be pushed below the HI is a question about mode count, and mode count is
what experiment 005 showed frequency buys.

### The figures, and the maths behind them

Drawn by `ska_hi_plots.py` from cached results; nothing below re-solves.
Every display equation is written on one source line because `md2pdf.py`'s
mathtext backend cannot span two.

**The solve.** The map-maker is a Wiener filter. With $A$ the pointing-and-beam
operator, $N$ the noise covariance, $S$ the prior covariance and $\mu$ the prior
mean, one channel's map is

$$m = (A^\top N^{-1} A + S^{-1})^{-1}(A^\top N^{-1} d + S^{-1}\mu)$$

which is **affine** in the data $d$: a fixed linear map $W = (A^\top N^{-1}A + S^{-1})^{-1}A^\top N^{-1}$ plus a constant prior term. Two consequences run through the whole experiment. Injecting a sky component adds linearly, $m(d + As) - m(d) = WAs$; and solving that component's own TOD with $\mu = 0$ kills the constant term and returns $WAs$ directly, at one solve instead of two. `validate_linearity` measures the identity at $5\times10^{-6}$.

The injected data is built as

$$d = \left[A(s_\text{fg} + s_\text{HI})\right](1 + g)(1 + \eta)$$

with the multiplicative 1/f gain $g$ and white $\eta$ replayed from their seeds, so the foreground TOD is the properly simulated one and only the HI rides in through $A$.

**The radial power spectrum.** Only $k_\parallel$ is well sampled here, so `pk_par` estimates a line-of-sight spectrum, averaged over the $N_\text{pix}$ lines of sight of the common patch. With $w_i$ a Blackman taper over $N_c$ channels of comoving depth $\Delta r$, $L = N_c\Delta r$ and $\bar{w^2}$ the mean squared taper,

$$P_{ab}(k) = \frac{\Delta r^2}{L\,\overline{w^2}}\left\langle \operatorname{Re}\left[\tilde a(k)\,\tilde b^{*}(k)\right]\right\rangle_\text{pix}, \qquad k = 2\pi f_\text{rfft}$$

The taper is not cosmetic: the foregrounds are $\sim\!10^4$ times the HI, so band-edge leakage would otherwise swamp every high-$k_\parallel$ bin.

**The clean.** PCA removes the highest-variance modes of the channel-channel covariance $C = \operatorname{cov}(x - \bar x)$, with $\bar x$ the per-channel mean over pixels. Keeping the leading $n$ eigenvectors as the columns of $U$,

$$x_\text{clean} = x - \left[U U^\top (x - \bar x) + \bar x\right]$$

This is the one step in the chain that is **not** linear in the data — $U$ is estimated from the data itself — which is exactly why the signal loss has to be measured by injection rather than computed.

**The transfer function.** Following Cunnington et al. (2023), with $X$ the data cube and $X_m$ an independent mock HI response,

$$T(k) = \frac{P\left(\mathcal{C}_n[X + X_m] - \mathcal{C}_n[X],\ X_m\right)}{P(X_m,\,X_m)}$$

where $\mathcal{C}_n$ is the PCA clean at $n$ modes. The numerator is a *cross*-power with the known mock; an auto-power of the difference would carry a positive noise bias. Since $X_m$ is each strategy's own map-made response $WAs_m$, $T$ is a ratio within one arm and is therefore **depth-independent** — which is what makes drift and raster comparable here even though the raster is $2.3\times$ shallower on the shared patch.

<div class="figblock">
<img class="figure" src="figures/hi_transfer_function.png" alt="Transfer function, drift against raster">
<p class="caption"><em>Figure 1 — $T(k_\parallel)$ per mode count, band $=\pm1\sigma$ over 20 mocks (Result 1). The raster is above the drift almost everywhere and the gap widens with aggressive cleaning. Loss piles up at low $k_\parallel$: those are the smoothest modes along the line of sight, so PCA cannot tell them from foreground.</em></p>
</div>

<div class="figblock">
<img class="figure" src="figures/hi_residual_vs_modes.png" alt="Residual against HI, versus modes removed">
<p class="caption"><em>Figure 2 — the post-clean residual in units of the HI it contains, $\sqrt{P_\text{clean}/P_\text{map made}}$, median over $k_\parallel$ (Result 3). Dashed is the pipeline as first run, solid after reconvolving to a common beam; both on the same 119 interior pixels, so the gap is the reconvolution alone. It halves the raster and does nothing for the drift. Removing more modes lowers the residual and costs signal, but never brings it near 1.</em></p>
</div>

**The field itself.** Before the power spectra, what is actually being
observed and solved. The drift's patch is a stack of three constant-Dec strips;
the raster's azimuth throw widens it 2.2x and wholly contains it, which is why
the two can be compared on a shared patch at all.

<div class="figblock">
<img class="figure" src="figures/hi_patch_maps.png" alt="The field and the two footprints">
<p class="caption"><em>Figure 5 — left: the GDSM foreground the survey sees, ~20 K against 0.5 mK of HI. Right: the nested footprints. The 119-pixel interior is what the common-resolution numbers are evaluated on, since reconvolving a zero-padded patch cannot be trusted near its boundary. Sequential ramp for temperature, categorical colours for membership — a diverging or rainbow scale would imply an ordering neither field has.</em></p>
</div>

**And the HI itself.** A map-domain panel was skipped when these figures were
first drawn, on the grounds that the cleaned map is pure residual and makes a
poor picture. That was half right: the cleaned panel *is* residual, and putting
it beside the others on its own colour scale turns out to be the clearest
statement of the result in the whole set.

<div class="figblock">
<img class="figure" src="figures/hi_maps.png" alt="HI in the map domain">
<p class="caption"><em>Figure 6 — one channel, frequency mean removed, which is the part a survey can use. Common resolution on the 119 interior pixels, the same footing as the quantitative figures. True HI, the drift's map-made version, the raster's, and the cleaned data on a scale 82x wider — at 4 modes removed nothing in that fourth image is signal. Diverging blue-grey-red, so zero recedes and sign is readable. The colour ratio is a 99th-percentile display scale at one channel, not the residual/HI ratio of Figures 1-3.</em></p>
</div>

**The ablation.** With the noise replayable, the same operator is re-solved three ways at fixed HI injection: `floor` from noiseless data, so the residual is only sky the strategy never measured and the prior filled in; `total` from the full data; and `noiseonly` as the difference of a noisy and a noiseless foreground solve. If `floor` and `total` agree, the noise plays no part.

<div class="figblock">
<img class="figure" src="figures/hi_ablation.png" alt="Ablation: floor against total against noise">
<p class="caption"><em>Figure 3 — they agree to 0.2% in power: the open rings (full data) sit on the floor line at every mode count, for both strategies, while noise alone runs $\sim\!38\times$ lower. The blocker is structural — the beam + prior floor — not sensitivity.</em></p>
</div>

**The frequency-structure split.** Write each cube as a per-pixel mean over channels plus a fluctuation, $x(\nu, p) = \bar x(p) + \delta x(\nu, p)$ with $\bar x(p) = N_c^{-1}\sum_\nu x(\nu, p)$. The two ratios that matter are

$$r_\text{coh} = \frac{\operatorname{rms}_p\,\bar x_\text{map made}}{\operatorname{rms}_p\,\bar x_\text{true}}, \qquad r_\text{var} = \frac{\operatorname{rms}\,\delta x_\text{map made}}{\operatorname{rms}\,\delta x_\text{true}}$$

A plain rms over the whole cube mixes the two and reported $0.924$ for the drift, which reads as near-perfect recovery and is badly misleading. Split apart it is $r_\text{var} = 0.380$ and $r_\text{coh} = 3.14$ — the drift discards most of the frequency-varying HI and *inflates* the frequency-coherent part threefold. Against the raster's $0.532$ and $1.11$ this is the mechanism behind Figure 1: with $\sim\!17$ measured modes the drift projects HI onto a nearly frequency-constant subspace, and anything frequency-constant is precisely what a foreground filter is built to delete.

<div class="figblock">
<img class="figure" src="figures/hi_frequency_structure.png" alt="Why the drift fails">
<p class="caption"><em>Figure 4 — left: the fraction of radial HI power surviving the map-maker, before any cleaning. Right: the split above, on a shared axis because both are map-made/true rms ratios. Result 2.</em></p>
</div>

### Why does even the cross-linked raster not reach the HI?

Asked 2026-08-19, and it splits into a configuration part and a structural one.

**A real gap in our pipeline, now closed.** Experiment 006 as first run never
reconvolved the channels to a common angular resolution. Every intensity-mapping
pipeline does this before cleaning, precisely so the beam cannot imprint
spectral structure: our beam narrows **12.5% across 350-400 MHz**, and leaving
that in attributes to the survey a chromatic term any real analysis would have
removed. `ska_hi_analysis.common_resolution` now brings each channel to the
widest beam in the band with a Gaussian kernel of
$\theta_k = \sqrt{\theta_\text{max}^2 - \theta_\nu^2}$ — exact for the Gaussian
beam used here, approximate only for the tapered-aperture beam, which these runs
do not use. Because it smooths a patch embedded in a zero-filled sphere, the
result is only trusted on pixels more than 3° from the boundary, which costs
well over half the patch; the un-smoothed control is therefore evaluated on that
same subset rather than on all 277 pixels.

**Re-run 2026-08-19 with all three variants stored.** Two corrections, which
must not be confused with each other — residual/HI in amplitude, median over k,
4 modes removed:

| | as run, 277 px | as run, 119 interior px | common resolution, 119 px |
|---|---|---|---|
| drift | 570.6 | 322.1 | 347.3 |
| raster | 65.9 | 48.5 | **24.0** |
| cross-linking gain | 8.7$\times$ | 6.6$\times$ | **14.4$\times$** |

The interior restriction alone accounts for a large part of the original
numbers: the full-patch values carried an edge inflation that had nothing to do
with the physics. The reconvolution then **halves the raster** (48.5 → 24.0, and
33.6 → 17.2 at 8 modes, 29.4 → 13.4 at 10) and does **nothing for the drift**,
which has too few measured modes for a cleaner beam to help. The transfer
function drops slightly under reconvolution (raster median $T$ 0.920 → 0.882 at
4 modes) because smoothing correlates the channels a little and PCA takes more
with it — but the residual falls twice as fast, so the net is a clear gain.

**The cross-linking case is therefore stronger than first published**, 14.4$\times$
rather than 8.7$\times$ at 4 modes. Everything quoted before 2026-08-19 was the
as-run full-patch column; prefer the last one.

**What is left is not a sensitivity problem.** From the ablation, the raster at
8 modes is floor $45.9\times$ and noise $4.9\times$ the HI, which combine in
quadrature to $46.2$ against a measured $47.0$. Drive the noise to zero with
infinite integration and the residual moves to $45.9\times$ — a **2%**
improvement. More dishes behave identically: 64 dishes on the same track give
64 times the samples of the *same* modes, buying noise and not information. The
usual "integrate longer, add dishes" argument does not engage with this blocker
at all. What does move it is mode count, which is geometry: 17 modes to 34 took
the floor from $572\times$ to $65\times$.

**And our band is close to a worst case.** Measured on this patch:

| | 350 MHz | 400 | 700 | 1000 |
|---|---|---|---|---|
| foreground rms | 2.152 K | 1.515 | 0.341 | 0.129 |
| HI $T_b$ | 0.541 mK | 0.492 | 0.270 | 0.149 |
| **foreground / HI** | **3978** | 3082 | 1266 | **862** |
| beam change over 50 MHz | **12.5%** | — | 6.7% | **4.8%** |

Foreground contrast is $4.6\times$ worse and the beam $2.6\times$ more chromatic
than at 1 GHz. Both stack against us, and both were chosen deliberately — 350 MHz
is the series anchor — but it means these numbers should not be read as typical
of Band 1.

**So should we expect to reach the HI, as MeerKLASS does in autocorrelation?**
Not in *this* configuration, and this experiment does not bound the method. Two
differences matter more than the ones usually cited, and neither is sensitivity:

1. **Their strongest results are cross-correlations.** Our residual is a bias
   uncorrelated with the HI, so against an external tracer it largely nulls,
   inflating the variance rather than biasing the answer. That is a far weaker
   requirement than residual < signal, which is what the numbers above measure.
2. **Survey volume does the averaging.** This patch holds roughly 13 beam areas
   times 16 $k_\parallel$ bins, of order 200 independent modes — which is exactly
   why our own cross-power estimator was too noisy to use (see the deliverable
   note above). A survey over thousands of square degrees averages a residual
   down in a way 277 pixels never can.

The honest scope of the result is therefore: **drift scanning at the bottom of
Band 1, on 3 h, is information-limited and cannot reach the HI, and no amount of
integration or extra dishes changes that.** Extending it to a MeerKLASS-like
configuration is a separate experiment that has not been run. Specific MeerKLASS
survey parameters are deliberately not quoted here — check them before putting
any in the write-up.

### Would a different cleaning basis help? (2026-08-25)

**Question.** If PCA is not removing the floor, would ICA, GMCA, NMF, a
polynomial expansion or GPR do better? This is answerable rather than a matter
of taste, because every one of those except GPR removes a **rank-N subspace** of
the channel-channel covariance and they differ only in how they choose it. So
the measurement to make is how many modes each component actually occupies.

**Method.** `ska_hi_rank.py`. Eigen-decompose the channel-channel covariance of
three cubes on the common patch: the GDSM truth, the beam + prior floor (the
residual of a *noiseless* solve, i.e. sky the strategy never measured), and the
map-made HI response. Modes needed to carry a given fraction of the variance:

| | 90% | 99% | 99.9% | $\lambda_{10}/\lambda_1$ |
|---|---|---|---|---|
| GDSM foreground | **1** | **1** | **1** | 2.6e-17 |
| floor, drift | 2 | 6 | 10 | 6.3e-4 |
| floor, raster | 1 | 2 | 5 | 5.7e-5 |
| **HI, drift** | **14** | **23** | **28** | — |
| **HI, raster** | **17** | **25** | **30** | — |

<div class="figblock">
<img class="figure" src="figures/hi_rank.png" alt="Eigenspectra of foreground, floor and HI">
<p class="caption"><em>Figure 7 — left: eigenspectrum of the channel-channel covariance, normalised to the first mode. Right: cumulative variance. The foreground falls off a cliff after mode 1; the floors decay fast; the HI decays slowly and is still climbing at mode 20.</em></p>
</div>

<div class="figblock">
<img class="figure" src="figures/hi_eigenvectors.png" alt="The spectral shapes, and how much each carries">
<p class="caption"><em>Figure 8 — top: the leading eigenvectors, i.e. the spectral shapes each component is built from. The foreground's mode 1 is a single smooth power law; its modes 2--4 look like structure but carry $\lambda/\lambda_1\sim10^{-17}$ and are floating-point noise. Bottom: how much each mode carries, on a <em>linear</em> axis. The foreground's modes 2--6 total $2\times10^{-8}$ of mode 1 — a blank column rather than a short one, which is what rank 1 looks like. The floor is down to 0.08 by mode 2; the HI is still at 0.51 by mode 6.</em></p>
</div>

**Result 1: the foreground is rank 1.** Over a 50 MHz block the GDSM sky is a
single spectral shape — $\lambda_{10}/\lambda_1 = 3\times10^{-17}$. One PCA mode
removes essentially all of it, and the ablation agrees: after **1** mode the
drift residual is already 896x the HI. **So from mode 2 onward the filter is not
cleaning foregrounds at all — it is fighting the instrumental floor.**
"Foreground cleaning" is a misnomer for what this pipeline does.

**Result 2: the HI is the incompressible component.** It needs 23-25 of 32 modes
for 99%, against 2-6 for the floor. That is expected — it is a stochastic field
with a short line-of-sight correlation length — but it is the fact that makes
the whole thing work at all: floor and signal *are* distinguishable by rank.

**Result 3: and yet removing the floor still costs the signal.** Mean cosine of
the principal angles between the leading rank-6 subspaces of floor and HI is
0.597 (drift) and 0.737 (raster) — angles of 53 and 42 degrees, far from
orthogonal. The floor's few modes are the smooth ones, and while the HI's power
is spread over ~25 modes, its **low-$k_\parallel$** power sits in the same smooth
few. That is exactly the measured $T(k)$ shape: 0.06 at the lowest
$k_\parallel$ rising to 0.95 at the highest. (Read the rank-6 overlap with
care for the raster, whose floor is rank 2 — modes 3-6 there carry almost no
variance, so the comparison is noisy.)

**What it changed.** The answer to the basis question is **no**. ICA, GMCA, NMF
and kernel PCA would find much the same subspace; the differences between them
matter for a mixture of astrophysical components with different statistics,
which is not what this is — one rank-1 sky plus a deterministic instrumental
floor. Polynomial modes should be *worse*, since they impose a smoothness the
floor does not have where PCA at least adapts. All are cheap to falsify (one
line in `ska_hi_analysis.pca_clean`; fastbox ships `ica_filter`, `nmf_filter`,
`kernel_pca_filter`, `gpr_filter`), so worth an afternoon as a null result.

**Correction to the queued work.** `HANDOFF.md` previously put GPR first, on the
reasoning that the blocker was a *non-smooth* floor PCA could not represent.
This measurement falsifies that: the floor is low-rank and PCA can represent it
fine. The blocker is the **overlap** with the signal, which no choice of basis
escapes. What follows instead: (i) **model** the floor rather than filter it —
it is $(I - WA)s$ and both $W$ and $A$ are known exactly in simulation, so its
covariance can go into the noise model rather than being blindly projected out;
(ii) **shrink** it with geometry — the raster's floor is already rank 2 against
the drift's 6, so more measured modes give a more compressible floor. GPR stays
on the list but not at the top.

### Position after 006 — consolidate before improving

Stated 2026-08-18, at the end of the session that produced the above. The
result is a negative one, so **the next step was figures and a presentable
write-up of the non-recoverability, not another attempt to improve recovery.**
**Done the same day:** `ska_hi_plots.py` draws four figures from cached results
(listed in `HANDOFF.md`), and `ska_hi_ablation.py` now saves the floor/total/
noise split k-resolved rather than printing medians to a console, which is what
made the third figure drawable. The ablation re-run reproduced the original
numbers exactly (drift 571.8x / 572.4x at 4 modes, raster 65.3x / 66.0x).

The reasoning is worth recording because it cuts against the obvious instinct.
A negative result that is *understood* — the blocker is the beam + prior floor,
it is not noise (floor and total agree to 0.2%), and PCA cannot touch it
because the chromatic null space makes it spectrally non-smooth — is a
publishable statement about drift scanning on its own. Rushing to improve
recovery before that is drawn up risks losing the clean version of the claim,
and risks tuning against a result nobody has yet inspected visually. The
machinery is documented in `METHODS.md` §11 so the presentation work does not
have to re-derive it.

One consequence for the improvement work when it does start: since the blocker
is a *non-smooth* residual, a filter that does not assume spectral smoothness
(GPR, or a foreground model built from the null space itself) is a more
promising lever than tuning the PCA mode count, which the residual/HI table
above shows saturating.


### 2026-08-26 — the map-domain figure, and how much the convention moves it

Figure 6 was still built from the as-run 277-pixel cube at native per-channel
resolution, while Figures 1-3 had moved to common resolution on the 119
interior pixels. Switching it over, and what that cost.

**Three of the four panels did not need it.** Each panel is a single channel,
so there is one beam in it and nothing chromatic to correct. The exception is
the cleaned panel, where PCA ran *across* frequency: on as-run data it shows a
residual the current pipeline no longer produces, which is the real reason to
switch.

**Reconvolving costs a factor ~3 in the HI itself.** At channel 16, 99th
percentile of $|\delta T_b|$:

| | true HI | raster map-made | cleaned (4 modes) | colour ratio |
|---|---|---|---|---|
| as run, 277 px | 0.699 µK | 0.341 | 35.8 | 51x |
| as run, 119 interior | 0.632 | 0.243 | 25.6 | 40x |
| **common resolution, 119 px** | **0.223** | 0.164 | 19.7 | **82x** |

The signal is small-scale, so smoothing to the widest beam in the band eats
most of it while the residual only falls ~1.8x. The displayed ratio therefore
*rises*. That is the honest cost of a standard procedure, not a defect.

**That ratio is not the residual/HI number.** 82x is a 99th-percentile display
scale at one channel; the 24x in Figures 1-3 is a median over $k$ of a power
ratio. They measure different things and will never agree, so the label reads
"the shared scale" rather than anything that looks like a physical ratio.

**One claim did not survive the change, and it is worth flagging.** The
frequency-structure split was quoted as drift $r_\text{var} = 0.380$,
$r_\text{coh} = 3.14$ against raster $0.532$, $1.11$ — the basis for "the drift
makes HI look like foreground" and for the Figure 6 caption's "visibly washed
out". At common resolution on the interior those become:

| | $r_\text{var}$ | $r_\text{coh}$ |
|---|---|---|
| drift | 0.380 → **0.721** | 3.14 → **2.61** |
| raster | 0.532 → **0.776** | 1.11 → **1.04** |

The drift's apparent loss of frequency-varying structure was largely the
small-scale structure the reconvolution removes from *both* cubes. On matched
resolution the two strategies retain nearly the same varying fraction (0.72 vs
0.78), and the mechanism rests on the coherent term instead: the drift still
inflates the frequency-*constant* part 2.6x, and a foreground filter is built
to delete exactly that. The headline is untouched — residual/HI is 347x drift
against 24x raster, and both of those are already common-resolution numbers —
but Figure 4 and its subtitle still quote 0.38/3.14 from the as-run cubes and
need the same treatment before the paper uses them.


### 2026-08-26 — the 1/f term, in the paper's own units and with a bound

The series has said "1/f is not the limiter" since experiment 003, on the
strength of a quadrature split of the galactic-plane residual: 21.4 K beam +
prior floor, 7.3 K white, 4.8 K 1/f, closing to better than 0.1%. That claim
was true but unusable in the paper, because it is a map-domain rms in kelvin,
on the full patch, at native per-channel resolution, with no HI in the
simulation. Every other number in the paper is a $k$-resolved ratio against
the HI at common resolution on the 119 interior pixels. Three steps put 1/f on
that footing and turn the claim into a bound.

**Splitting the noise costs two solves per channel.** The ablation already
re-solved the same operators with components removed; its `noiseonly` arm was
the difference between a noisy and a noiseless solve of the foreground. Since
the noise enters multiplicatively,

$$d = (\mathbf{A}s)\,(1 + g)\,(1 + w),$$

the gain and white contributions can be isolated by applying one factor at a
time, and the two pieces sum to the joint term up to the cross product $g w$.
Both draws are of order $10^{-3}$, so that term is $\sim\!10^{-6}$ relative;
measured, `gainonly + whiteonly` reproduces `noiseonly` to $1.9\times10^{-4}$
(drift) and $2.3\times10^{-4}$ (raster) of the noise amplitude. The script
prints the check on every run, since it is the assumption the split rests on.

At 4 modes removed, common resolution, as amplitude ratios to the HI:

| | floor | total | 1/f | white |
|---|---|---|---|---|
| drift | 336.3x | 349.4x | **4.0x** | 5.2x |
| raster | 22.7x | 24.0x | **0.4x** | 2.1x |

For the raster the 1/f term alone sits *below the HI it is trying to hide*.
That is a considerably sharper statement than "5-8% of the residual", and it
is in the same units as the transfer function and the residual/HI table.

**What the model actually is.** `limTOD.flicker_model` draws the gain from a
Gaussian process with a closed-form autocorrelation rather than by sampling a
spectrum, so the realisation is exact and has no periodicity artefacts. The
correlation is the cosine transform of a cut-off power law,

$$C(\tau) = \frac{1}{\pi}\int_{\omega_c}^{\infty}
           \left(\frac{\omega_0}{\omega}\right)^{\alpha}\cos(\omega\tau)\,
           \mathrm{d}\omega ,$$

evaluated through the upper incomplete gamma function, with
$C(0) = (\omega_c/\pi)(\omega_0/\omega_c)^{\alpha}/(\alpha - 1) + \sigma_w^2$;
the $1/(\alpha-1)$ is why $\alpha = 1$, the pure $1/f$ exponent, is rejected
rather than returning infinity. Verified against a numerical cosine transform
to $\sim\!10^{-4}$ out to $\tau = 300$ s. **$\omega_0$ and $\omega_c$ are
angular frequencies** — `GAIN_PARAMS` = $[1.335\times10^{-5},\,
1.099\times10^{-3},\, 2]$ is 2.13e-6 Hz and 1.75e-4 Hz, and quoting the raw
numbers as Hz is wrong by $2\pi$.

**The knee is the useful parameterisation.** The white gain term has flat
spectrum $\sigma_w^2 \Delta t$ in the same convention, so the two cross at

$$\omega_{\rm knee} = \omega_0 \left(\sigma_w^2 \Delta t\right)^{-1/\alpha}
  = 0.95\ \text{mHz} .$$

That number locates the assumption against the scan: the drift's fundamental
(0.28 mHz, one traverse per 1 h pass) sits *below* the knee, in the
1/f-dominated regime, while the raster's sweep rate (3.33 mHz, 300 s per
sweep) sits 3.5x *above* it. This is the "scan speed is already solved"
argument stated as a crossing frequency rather than a power ratio.

**How much worse would it have to be?** Scanning the knee and reading off
where the 1/f arm meets the floor:

| | $\alpha = 1.5$ | $\alpha = 2.0$ | $\alpha = 2.5$ |
|---|---|---|---|
| drift | 316 mHz (333x) | 80 mHz (84x) | 33 mHz (35x) |
| raster | 184 mHz (194x) | 57 mHz (60x) | 21 mHz (22x) |

Quote the raster at $\alpha = 2.5$: even in the least favourable corner, the
knee has to move 22x before 1/f reaches the floor.

*That axis is analytic, and the write-up must say so.* Scaling one realisation
scales the entire pipeline linearly: `solve` is affine, so the constant term
cancels in the noisy-minus-quiet difference, and `pca_clean` is
scale-invariant, because scaling a cube scales its covariance and leaves the
eigenvectors untouched. Hence
$\mathrm{ratio}(\omega_{\rm knee}) \propto \omega_{\rm knee}^{\alpha/2}$
exactly — confirmed against the computed grid to $2\times10^{-4}$. Only the
three points at the fiducial knee are independent measurements, and the figure
encodes that with lines for the scaling and a filled marker at the anchor.
Drawing a marker at every knee would claim eighteen runs where there were
three.

**The cut-off, which is the honest test.** $\omega_c$ changes the *shape* of
the correlation rather than its amplitude, so nothing factors out and every
point is a fresh draw. The fiducial 0.175 mHz is a 1.6 h timescale — about one
pass — which is a suspiciously convenient place for the model to stop adding
correlated power. Extending it downward:

| $\omega_c$ | timescale | gain rms | drift | raster |
|---|---|---|---|---|
| 0.0175 mHz | 16 h | 7.12e-4 | 4.44x | 0.40x |
| 0.055 | 5 h | 4.13e-4 | 4.30x | 0.45x |
| **0.175 (assumed)** | 1.6 h | 2.69e-4 | 3.98x | 0.38x |
| 0.55 | 30 min | 1.51e-4 | 3.20x | 0.37x |
| 1.75 | 10 min | 7.96e-5 | 2.31x | 0.15x |
| 5.5 | 3 min | 4.21e-5 | 0.22x | 0.09x |

Left of the assumed value the curve is **flat while the gain rms grows 2.6x**.
The reason is a degeneracy, and it is the same one that runs through the whole
experiment: gain drift slower than a pass acts as a near-constant
multiplicative error on the foreground, and the foreground is rank 1, so the
product is rank 1 too and the cleaning removes it along with the foreground.
The result therefore does not rest on where the spectrum is truncated. Right
of the assumed value the residual falls, but only because raising $\omega_c$
deletes 1/f power outright (rms drops 6x) — that side is not a finding.

**One soft spot, recorded rather than hidden.** The two lowest-$\omega_c$
raster points are non-monotonic in rms (3.28e-4 then 3.80e-4) where the
drift's are clean. A 3 h observation barely samples a 16 h correlation time,
so a single draw of the lowest-frequency modes is not representative. The flat
trend is not in doubt, but those individual points carry realisation scatter
that one draw cannot quantify; a few seeds would give error bars. Related, and
worth a sentence in the paper: the gain draw is common to all channels by
construction of the RNG replay, so the model has no frequency structure in the
gain — and frequency structure is exactly what would decide whether PCA could
remove it.
