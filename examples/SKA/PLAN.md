# Plan to first draft

**Target: complete draft of the drift-scan paper by end of November 2026.**
Written 2026-09-17. This is the forward-looking doc — scope, cut list and
schedule. It deliberately holds almost no numbers; those live in `HANDOFF.md`
and are quoted from there, because this series has already been bitten once by
a duplicated table going stale (`HANDOFF.md`'s old "Reference numbers"
section). If a number matters here, follow the pointer.

Sibling docs and their remits: `README.md` (index), `METHODS.md` (how the
machinery works, atemporal), `HANDOFF.md` (state and numbers),
`ANALYSIS_LOG.md` (chronological record and maths). `paper/main.tex` is being
written by hand from 2026-08-26 onward.

## Where we are

Six experiments are complete and plotted: off-plane drift (001/002),
galactic-plane drift with a 1/f budget (003), MeerKLASS-style constant-elevation
raster (004), the Band 1 frequency sweep (005) and the HI recovery experiment
(006). Six figures regenerate from cache in seconds via `ska_hi_plots.py`. The
science content is, in one line each:

- **No HI is recoverable at 350–400 MHz.** Post-clean residual is 347× (drift)
  and 24× (raster) the signal at 4 modes removed.
- **It is not a noise problem.** Ablation attributes essentially all of it to
  the beam + prior floor; 1/f alone sits *below the HI* for the raster.
- **The binding constraint is the measured mode count.** The drift measures
  ~17 sky modes over the patch, the raster ~34; everything else is prior-filled.
  Frequency, cross-linking and the elevation ladder are all levers on that one
  number, and integration time is not.

That is a coherent, publishable negative result. The remaining work is about
making it useful to somebody designing a survey, and about defending it.

## The decision that sets the scope

**Which paper is this?** The answer changes the cut list below and it should be
settled before any more runs are queued.

- **(A) The negative result, with a design figure as its synthesis.** "Drift
  scanning is not viable for SKA-Mid Band 1 single-dish IM, here is why, and
  here is the chart that tells you which levers move it." Mostly written.
  Makes November.
- **(B) A survey-design framework**, with the Band 1 case as the worked
  example. Wants the full lever set — real beams, RFI, spillover — and each of
  those is new runs plus a new robustness section. Does not make November.

**Recommendation: (A)**, with the design figure doing the work that makes it
more than a null result. (B) is the natural follow-up and inherits everything
deferred below.

## Needed for the draft

| # | item | why it is in | effort |
|---|---|---|---|
| 1 | **The design / trade-space figure** — residual/HI against measured modes per unit sky, with the 1/f family and the beam+prior floor asymptote. See "Figure honesty" below. | This is the primary result. Nothing currently puts the levers on one canvas. | days |
| 2 | **Flat-prior audit of every headline number** | The prior has already manufactured one result (the 005 frequency trend) and masked another (cross-linking). Non-negotiable before anything becomes a headline. | ~1 day, `audit_flat_prior.py` exists |
| 3 | **Scan period sweep** | `f_scan / f_knee` is the parameter a designer actually controls, and it gives the cleanest statement of the cross-linking result: the drift's fundamental sits *below* the knee, the raster's sweep above it. Partly quantified already (HANDOFF 0e(ii)). `sweep_s` is already an argument to `constant_elevation_scan`. | cheap |
| 4 | **Dish count, stated analytically** | Co-pointed dishes add samples, not rows of `A`, so they cannot move a floor-dominated residual — it is the same statement as "integration time does not add modes". First question any referee or survey designer asks. Worth a paragraph and no new runs. Note the one real exception: dishes at *different elevations* buy the ladder for free. | hours |
| 5 | **Depth-matched raster** (queued item 1) | The 006 raster is not depth-matched; does not affect T(k) but does affect the residual comparison that the paper leads with. | queued, scoped |
| 6 | **More crossing angles** (queued item 2 / item 0) | The paper claims mode count is the binding constraint. It needs at least one demonstration that deliberately raising it moves the residual. | days |
| 7 | **Seed scatter / error bars** | The known soft spot at low `w_c`, and generally: single-realisation numbers in a headline table will be challenged. ~5 min of compute for the worst case. | cheap |

## Explicitly deferred

Not in the November draft. Each is a genuine question, and each is a robustness
section that could swallow the schedule.

- **Measured (holographic) beam.** Right idea — the spectral ripple from the
  feed/subreflector standing wave is non-smooth in frequency and so immune to
  PCA, which is exactly this paper's failure mode, and a real beam should make
  it worse. But `katbeam`/MeerKAT holography covers UHF 580–1015 and L-band
  856–1712 MHz; the HI experiment is at 350–400, where no measured MeerKAT beam
  exists. Doing it properly means re-running at 700–1050 MHz. That is the first
  item of the follow-up, not a check to bolt on here.
- **RFI flagging fraction.** Missing data changes the null space, so it hits
  the floor rather than the noise — a first-class lever, not a nuisance. Needs
  a flagging model this series does not have.
- **Spillover / ground pickup versus elevation.** The hidden cost of the
  elevation ladder, which currently looks free.
- **Is 1/f common-mode across dishes?** Determines whether multi-dish averaging
  helps at all. One parameter, but only meaningful once there is a dish count.
- **Model the floor instead of filtering it** (queued 0c(i)). The strongest
  idea in the queue and a research direction in its own right — `(I - WA)` is
  known exactly in simulation. This is the next paper.
- **675–725 MHz repeat** (queued 0b(a)), GPR and alternative cleaning bases
  (0c(iii)) — the latter is an afternoon and a null result if anyone wants it.

## Figure honesty

Constraints the design figure has to respect, all established:

- **Most of the knee axis is a scaling law, not data.** `ratio ~ knee^(α/2)`
  holds exactly (verified to 2e-4) because the pipeline is affine in the gain
  draw. Draw lines for the scaling, filled markers only at the measured
  anchors — three per panel. The existing `fig_1f_scan` already does this; the
  new figure must not quietly undo it.
- **Do not put the high-pass cutoff on the chart** without `T(k)` beside it.
  It is the only lever whose ratio crosses 1, and it almost certainly crosses
  because the filter removes the HI along with the 1/f. Unqualified, a reader
  reads off "use a 5.5 mHz high-pass and the survey works".
- **Normalise modes per unit sky.** Frequency buys modes and shrinks the patch
  ~2.3× across the band; unnormalised, the chart recommends 1050 MHz
  unconditionally.
- **Carry a cost axis.** Fixed telescope-hours, or time as an explicit axis. A
  trade-space plot without a cost cannot trade.
- **Put the dead levers on it.** Tsys, integration time and co-pointed dish
  count are all confirmed non-levers because the residual is floor-dominated.
  Publishing which knobs do *not* work is rarer and more useful than another
  sensitivity curve.
- `ratio = 1` is necessary, not sufficient.

## Schedule

Ten weeks. Two are held back deliberately — this series has revised its
headline numbers twice (2026-08-19 reconvolution, 2026-08-26 common
resolution), and both times it was the check, not the run, that found it.

| when | what |
|---|---|
| **w/c 22 Sep** | Settle (A) vs (B). Items 2 and 4 — both cheap, both risk-reducing, and item 2 gates everything downstream. |
| **Oct, first half** | Items 3, 5, 7. Prototype the design figure against cached results; iterate on form while the underlying runs are still moving. |
| **Oct, second half** | Item 6 (more crossing angles) — the only item needing real compute. Design figure to final form. |
| **Nov, first half** | Write. `paper/main.tex` is hand-written from here; `ANALYSIS_LOG.md` and `main.tex` still carry pre-2026-08-26 phrasing that needs reconciling with the common-resolution numbers. |
| **Nov, second half** | Buffer, internal read, figure polish. |

## Risks

1. **The 2026-08-26 wording drift.** `ANALYSIS_LOG.md` and `paper/main.tex`
   were deliberately left carrying the old phrasing ("~38x", "0.2%") after the
   common-resolution re-run. Reconciling them is a writing task, not a physics
   one, but it is easy to forget and it is wrong in the current text.
2. **Item 6 is the only schedule risk with compute in it.** If more crossing
   angles do not move the residual, that is a *result* and the paper absorbs
   it — but it needs to be known by end of October, not in November.
3. **Scope creep toward (B).** Every deferred item above is genuinely
   interesting. The list exists so they can be declined quickly.
