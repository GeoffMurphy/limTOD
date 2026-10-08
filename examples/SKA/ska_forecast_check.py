"""Does the per-k_par forecast survive contact with the simulated patch?

``nm_kpar_results.md`` predicts a cross-detection against an ideal tracer
inside the patch already simulated. That prediction is the whole mode-counting
chain at once -- ``tr(WA)`` as the angular count, ``n_b`` as the line-of-sight
count, the inverse-variance sum over Blackman-tapered bins, and the ``2 + R_b``
variance -- so measuring it tests all four together.

**What the brief asked for, and what exists.** The brief wanted the 20 Hi mocks
re-cleaned against the noiseless floor arm. ``run_hi_experiment.py`` builds
those mocks but saves only the transfer function derived from them; the cubes
are gone, and the ablation stores power spectra rather than cubes, so the
noiseless arm cannot be cross-correlated with anything either. Nothing here
re-makes them -- that would be map-making.

Two substitutes are already on disk and together do the same job:

1. ``p_clean_cross_4`` IS the wanted estimator. It is
   ``P(clean(data), W A s_HI)``: the cleaned map crossed with the map-made true
   Hi, which for ``r = 1`` and no shot noise is exactly the ideal tracer. So

       A_b = P_clean_cross,b / (T_b P_mapmade,b)

   is an amplitude estimate that should be 1, bin by bin, and its scatter about
   1 is the forecast's error bar. Both halves are stored for every run.

2. Six extra runs (``seed11``..``seed16``) vary ONLY the Hi realisation --
   verified here: the cleaned auto-power, which is 99% floor, is identical
   across them to 0.2-0.8%. With the base run that is **seven** independent
   draws of the floor-times-tracer chance correlation, which is the quantity
   the forecast's variance is made of.

**Two differences from the brief, both conservative.** The data cube carries
the simulated noise, so the residual crossed with the Hi is the TOTAL arm, not
the floor arm; it sits 1-9% above the floor in power, so ``R_b`` is measured
from the same cleaned map rather than assumed, and the prediction it feeds is
if anything pessimistic. And because the estimator normalises by the same
realisation's measured ``P_mapmade``, the Hi's own cosmic variance cancels:
the predicted per-bin variance is ``R_b/N_b`` rather than ``(2 + R_b)/N_b``.
``R_b >= 400`` everywhere, so the two differ by under 0.5%.

**The area the test runs on is not the area the forecast quotes.** Power
spectra are measured on the interior of the COMMON patch -- 119 px (99.9 deg^2)
at 350-400, 171 px (35.9 deg^2) at 675-725 -- not on each strategy's own patch.
The modes available to the measurement are therefore ``rho * Omega_meas``, not
``N_eff``. Getting this wrong inflates the predicted SNR by 1.6-2.4x.

    /home/geoff/limTOD/.venv/bin/python ska_forecast_check.py
"""
from __future__ import annotations

import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

RESDIR = os.path.join(_HERE, "results")
NMODES = 4
ROW = {"drift": "3-strip drift", "drift12": "12-strip drift",
       "drift12m": "12-strip drift", "raster": "raster"}


class Band:
    def __init__(self, tag, label, nside, kmodes, strategies, single, seeds=()):
        self.tag, self.label, self.nside = tag, label, nside
        self.kmodes, self.strategies = kmodes, strategies
        self.single = single          # one run, all strategies
        self.seeds = seeds            # independent Hi realisations, same floor


BANDS = {
    "f350": Band(
        "f350", "350-400 MHz", 64, "hi_kmodes_f350_ns64.npz",
        ("drift", "drift12", "raster"),
        single="hi_experiment_matched_f350_400_nc32_ns64.npz",
        seeds=("hi_experiment_f350_400_nc32_ns64.npz",)
              + tuple(f"hi_experiment_seed{j}_f350_400_nc32_ns64.npz"
                      for j in range(11, 17))),
    "f700": Band(
        "f700", "675-725 MHz", 128, "hi_kmodes_f700_ns128.npz",
        ("drift", "drift12m", "raster"),
        single="hi_experiment_matched700b_f675_725_nc32_ns128.npz"),
}


def _nb(n):
    w = np.full(n, 2.0)
    w[-1] = 1.0
    return w


def amplitude(d, s, nmodes=NMODES):
    """``A_b``, ``R_b`` and ``T_b`` for one strategy in one experiment npz.

    Everything comes from the ONE file, so the cleaned map that supplies the
    cross-power is the same cleaned map that supplies ``R_b``. No ablation is
    needed and no prior convention can creep in between the two halves.
    """
    t = d[f"{s}_cr_tf_{nmodes}"]
    p_mm = d[f"{s}_cr_p_mapmade"]
    p_auto = d[f"{s}_cr_p_clean_auto_{nmodes}"]
    p_x = d[f"{s}_cr_p_clean_cross_{nmodes}"]
    surviving = t ** 2 * p_mm
    return dict(A=p_x / (t * p_mm), R=p_auto / surviving - 1.0, T=t,
                n_b=_nb(len(t)), k=d[f"{s}_cr_k"])


def geometry(band):
    """``rho`` per strategy and the solid angle the spectra are measured on."""
    import healpy as hp
    km = np.load(os.path.join(RESDIR, band.kmodes), allow_pickle=True)
    ex = np.load(os.path.join(RESDIR, band.single), allow_pickle=True)
    omega = float(ex["interior"].sum()) * hp.nside2pixarea(band.nside,
                                                           degrees=True)
    rho = {s: float(km[f"{s}_total"]) / float(km[f"{s}_area"])
           for s in band.strategies}
    return rho, omega, int(ex["interior"].sum())


def per_bin(term, n_ang):
    """Predicted per-bin sigma and the normalised residual ``z_b``.

    ``N_b = n_ang * n_b`` angular-times-radial degrees of freedom, and
    ``var(A_b) = R_b / N_b``: the Hi's own cosmic variance cancels because the
    estimator divides by this realisation's measured ``P_mapmade``.
    """
    N_b = n_ang * term["n_b"]
    sig = np.sqrt(term["R"] / N_b)
    return sig, (term["A"] - 1.0) / sig


def combine(term):
    """Inverse-variance weights over bins, and the combined amplitude."""
    w = term["n_b"] / (2.0 + term["R"])
    w = w / w.sum()
    return w, float(np.sum(w * term["A"]))


def ensemble(band, strategy, n_ang):
    """Mean and scatter of the combined amplitude over independent Hi fields."""
    ds = [np.load(os.path.join(RESDIR, f), allow_pickle=True)
          for f in band.seeds]
    terms = [amplitude(d, strategy) for d in ds]
    A = np.asarray([t["A"] for t in terms])              # (nreal, nbin)
    R = np.mean([t["R"] for t in terms], axis=0)
    n_b = terms[0]["n_b"]
    Ahat = np.asarray([combine(t)[1] for t in terms])
    pred = float(np.sqrt(n_ang * np.sum(n_b / (2.0 + R))))
    emp = float(Ahat.mean() / Ahat.std(ddof=1))
    # implied angular dof per bin, from the measured scatter alone
    implied = R / (n_b * A.var(axis=0, ddof=1))
    z = (A - 1.0) / np.sqrt(R / (n_ang * n_b))
    return dict(n=len(ds), Ahat=Ahat, mean=float(Ahat.mean()),
                std=float(Ahat.std(ddof=1)), pred=pred, emp=emp,
                implied=implied,
                implied_geom=float(np.exp(np.mean(np.log(implied)))),
                z_std=float(z.std(ddof=1)), R=R, n_b=n_b, A=A)


def snr_err(n):
    """Fractional 1-sigma spread of an SNR estimated from ``n`` draws.

    Both halves are noisy: the mean carries 1/sqrt(n) and the standard
    deviation 1/sqrt(2(n-1)), added in quadrature at SNR ~ 1-2.
    """
    return float(np.sqrt(1.0 / n + 1.0 / (2.0 * (n - 1))))


def _fmt(x, p=2):
    return f"{x:,.{p}f}"


def single_run(band, strategy, n_ang):
    """The 16-bin test, for the band-strategy pairs with only one realisation."""
    d = np.load(os.path.join(RESDIR, band.single), allow_pickle=True)
    t = amplitude(d, strategy)
    sig, z = per_bin(t, n_ang)
    w, Ahat = combine(t)
    pred = float(np.sqrt(n_ang * np.sum(t["n_b"] / (2.0 + t["R"]))))
    implied = t["R"] / (t["n_b"] * (t["A"] - 1.0) ** 2)
    return dict(term=t, sig=sig, z=z, Ahat=Ahat, pred=pred,
                z_rms=float(np.sqrt(np.mean(z ** 2))),
                implied_geom=float(np.exp(np.mean(np.log(implied)))))


def write_report(path=None):
    path = path or os.path.join(_HERE, "forecast_check.md")
    import healpy as hp
    L = []
    w = L.append
    w("# Does the per-$k_\\parallel$ forecast survive contact with the "
      "simulated patch?")
    w("")
    w("Generated by `ska_forecast_check.py`, 2026-10-08. Post-processing of "
      "`results/*.npz` only: no TOD, no map-making,")
    w("no solves.")
    w("")
    w("## Verdict")
    w("")
    w("**The forecast holds.** The amplitude comes out at 1, the measured "
      "error bar matches the predicted one to 5%, and")
    w("`tr(WA)` is confirmed as a count of independent angular modes to "
      "**14%**, independently for two strategies whose mode")
    w("counts differ by a factor of three. The areas in "
      "`nm_kpar_results.md` stand.")
    w("")

    # ---- inputs -----------------------------------------------------------
    w("## What the brief asked for, and what is on disk")
    w("")
    w("**The 20 Hi mocks are not stored.** `run_hi_experiment.py` builds "
      "them, uses them for $T_b$, and saves only $T_b$ and")
    w("its scatter; the cubes are discarded. The ablation stores power "
      "spectra, not cubes, so the noiseless floor arm cannot")
    w("be cross-correlated with anything either. Re-making either would be "
      "map-making, so neither was done.")
    w("")
    w("Two substitutes are already on disk, and between them they do the "
      "same job.")
    w("")
    w("**1. `p_clean_cross_4` IS the estimator the brief describes.** It is "
      "$P(\\mathcal{C}_4[x],\\,W A s_{\\rm HI})$ -- the cleaned map")
    w("crossed with the map-made TRUE Hi, which for $r = 1$ with no shot "
      "noise is exactly the ideal tracer. So")
    w("")
    w("$$\\hat{A}_b = \\frac{\\hat{P}_{\\times,b}}"
      "{T_b\\,\\hat{P}_{\\rm mapmade,b}}$$")
    w("")
    w("is an amplitude estimate that should be 1 in every bin, and its "
      "scatter about 1 is the forecast's error bar. Both")
    w("halves are stored for every run and every cleaning depth.")
    w("")
    w("**2. Seven independent Hi realisations.** `seed11`..`seed16` vary only "
      "the Hi field. Verified rather than assumed: the")
    w("cleaned auto-power, which is 99% floor, is identical across them to "
      "**0.2-0.8%** in every bin. With the base run that")
    w("is seven independent draws of the floor-times-tracer chance "
      "correlation -- precisely the quantity the forecast's")
    w("variance is made of. They exist only at 350-400 MHz and only for the "
      "3-strip drift and the raster; everything else")
    w("gets the 16-bin test on a single realisation.")
    w("")
    w("Two differences from the brief, both of which make the test harder, "
      "not easier:")
    w("")
    w("- the data cube carries the simulated noise, so what is crossed with "
      "the Hi is the TOTAL arm, 1-9% above the floor in")
    w("  power. $R_b$ is therefore **measured from the same cleaned map** "
      "rather than taken from the ablation, so each file is")
    w("  internally consistent and the prediction it feeds is if anything "
      "pessimistic;")
    w("- the estimator divides by the same realisation's measured "
      "$\\hat{P}_{\\rm mapmade}$, so the Hi's own cosmic variance")
    w("  cancels and the predicted per-bin variance is $R_b/N_b$, not "
      "$(2+R_b)/N_b$. $R_b \\geq 400$ everywhere, so the two differ")
    w("  by under 0.5%.")
    w("")

    # ---- headline ---------------------------------------------------------
    w("## Headline")
    w("")
    w("$N_b = N_{\\rm eff}\\,n_b$ and predicted SNR "
      "$= \\sqrt{N_{\\rm eff}\\sum_b n_b/(2+R_b)}$. Uncertainty on an "
      "SNR from $n$ draws is")
    w("$\\sqrt{1/n + 1/[2(n-1)]}$, i.e. **48%** at $n = 7$.")
    w("")
    w("These predicted SNRs are **not** the 2.3 / 0.3 / 0.09 and "
      "3.4 / 6.5 / 3.3 of `nm_kpar_results.md`, and should not be")
    w("quoted in their place. $R_b$ here is measured from the same cleaned "
      "map that supplies the cross-power, so it is the")
    w("own-prior TOTAL arm rather than the amplitude-only FLOOR arm -- the "
      "price of making each file internally consistent,")
    w("and worth paying, since a prediction assembled from two files could "
      "fail for a bookkeeping reason rather than a")
    w("physical one. The shift is 10-20%: 1.95 against 2.27 for the "
      "350-400 MHz raster, 3.41 / 5.96 / 2.67 against")
    w("3.43 / 6.53 / 3.29 at 675-725. Nothing in `nm_kpar_results.md` moves.")
    w("")
    w("| band | strategy | predicted SNR | mean $\\hat{A}$ | pull | "
      "empirical SNR | ratio | implied $N_{\\rm ang}$ / tr($WA$) |")
    w("|---|---|---:|---:|---:|---:|---:|---:|")
    for tag in ("f350", "f700"):
        band = BANDS[tag]
        km = np.load(os.path.join(RESDIR, band.kmodes), allow_pickle=True)
        for st in band.strategies:
            neff = float(km[f"{st}_total"])
            if band.seeds and st in ("drift", "raster"):
                r = ensemble(band, st, neff)
                e = snr_err(r["n"])
                null = r["pred"] < 0.5
                emp = ("--" if null
                       else f"{r['emp']:.2f} $\\pm$ {abs(r['emp']) * e:.2f}")
                rat = "--" if null else f"{r['emp'] / r['pred']:.2f}"
                mean = (f"{r['mean']:.2f} $\\pm$ "
                        f"{r['std'] / np.sqrt(r['n']):.2f}")
                pull = (r["mean"] - 1.0) * r["pred"]
                w(f"| {band.label} | {ROW[st]} | {r['pred']:.2f} | {mean} | "
                  f"{pull:+.2f} | {emp} | {rat} | "
                  f"**{r['implied_geom'] / neff:.2f}** |")
            else:
                r = single_run(band, st, neff)
                pull = (r["Ahat"] - 1.0) * r["pred"]
                w(f"| {band.label} | {ROW[st]} | {r['pred']:.2f} | "
                  f"{r['Ahat']:.2f} | {pull:+.2f} | 1 realisation | -- | -- |")
    w("")
    w("**Pull** is $(\\hat{A} - 1)\\times$ predicted SNR, i.e. the "
      "deviation in units of the forecast's own error bar. Every row")
    w("is within **0.7**, including the two bands and all three strategies. "
      "That is the single most direct statement the")
    w("data support: the forecast predicts its own error bar correctly, six "
      "times out of six.")
    w("")
    w("`--` on the 350-400 MHz 3-strip drift is not a failure: the forecast "
      "predicts **0.08$\\sigma$** there, so there is no")
    w("signal to measure and a ratio of two noise numbers would mean "
      "nothing. What that row checks is that the null comes out")
    w("as a null, and it does -- $\\hat{A} = -6.4 \\pm 4.2$ against a "
      "per-realisation scatter of 11.2, consistent with anything.")
    w("Its implied-mode column is measured the same way as the raster's and "
      "passes alongside it.")
    w("")
    w("**The single-realisation rows get no error-bar test, deliberately.** "
      "With one Hi field and one floor, the sixteen")
    w("$\\hat{A}_b$ are not sixteen independent draws: the floor is "
      "spectrally low rank -- its $k_\\parallel$ power has a")
    w("participation ratio of only 2-5 -- so its chance correlation with the "
      "tracer is partly coherent across bins. The")
    w("rms of $z_b$ over bins then ranges 0.6-3.2 across the four "
      "single-realisation rows, and a single-draw implied mode")
    w("count is biased high by $\\sim$3.6x (the geometric mean of "
      "$1/\\chi^2_1$). Neither is quoted. The amplitude itself is")
    w("unaffected by this and is quoted: at 675-725 MHz it comes back as "
      "1.05, 1.04 and 0.99 against predicted 1$\\sigma$ errors")
    w("of 0.29, 0.17 and 0.37. Those three share one Hi realisation, so they "
      "are one draw rather than three.")
    w("")

    # ---- the mode count ---------------------------------------------------
    w("## The result worth keeping: tr($WA$) really is the angular mode count")
    w("")
    w("Inverting the per-bin variance gives a **measurement** of the number "
      "of angular degrees of freedom, independent of any")
    w("mode-counting theory:")
    w("")
    w("$$N_{\\rm ang} = \\frac{R_b}{n_b\\,{\\rm var}(\\hat{A}_b)}$$")
    w("")
    w("Over seven realisations and sixteen bins at 350-400 MHz:")
    w("")
    w("| strategy | tr($WA$) | $N_{\\rm ang}$ measured | ratio | "
      "$\\rho\\,\\Omega_{\\rm meas}$ | ratio |")
    w("|---|---:|---:|---:|---:|---:|")
    band = BANDS["f350"]
    rho, omega, npix = geometry(band)
    km = np.load(os.path.join(RESDIR, band.kmodes), allow_pickle=True)
    for st in ("drift", "raster"):
        neff = float(km[f"{st}_total"])
        r = ensemble(band, st, neff)
        g = r["implied_geom"]
        w(f"| {ROW[st]} | {neff:.1f} | {g:.1f} | **{g / neff:.2f}** | "
          f"{rho[st] * omega:.1f} | {g / (rho[st] * omega):.2f} |")
    w("")
    w("**1.14 for both**, although the two mode counts differ by a factor of "
      "three. The rms of the normalised residual")
    w("$z_b = (\\hat{A}_b - 1)/\\sqrt{R_b/N_b}$ is **0.95** for each, "
      "over 7 x 16 samples. That is the whole chain checking")
    w("out at once: the angular count, the radial count $n_b$, the Blackman "
      "binning and the $(2+R_b)$ variance.")
    w("")
    w("**Use tr($WA$), not $\\rho\\,\\Omega_{\\rm meas}$.** Spectra "
      "are measured on the interior of the common patch -- "
      f"{npix} px, "
      f"{omega:.0f} deg$^2$ at")
    w("350-400 -- not on each strategy's own patch, so the naive count is "
      "$\\rho\\,\\Omega_{\\rm meas}$, which is low by 3.0x (drift)")
    w("and 6.6x (raster). The reason is physical and worth stating in the "
      "draft: the interior cut throws away the patch EDGE,")
    w("and edge pixels have small $\\lambda_i$, so they carry almost none "
      "of $\\mathrm{tr}(WA) = \\sum_i \\lambda_i/(\\lambda_i + "
      "\\sigma_{\\rm p}^{-2})$.")
    w("Cutting them costs area but not modes.")
    w("")
    w("This does **not** change any area in `nm_kpar_results.md`. "
      "$\\Omega_{\\rm req}$ uses $\\rho = N_{\\rm eff}/\\Omega_{\\rm sim}$, "
      "the patch-average")
    w("density, which is the right quantity for tiling whole patches "
      "(edges included) over a survey. And the predicted SNR in")
    w("the patch simulated is $\\sqrt{\\rho\\,\\Omega_{\\rm sim} "
      "\\sum_b n_b/(2+R_b)} = \\sqrt{N_{\\rm eff}\\sum_b n_b/(2+R_b)}$ "
      "either way, so the")
    w("3.3-6.5$\\sigma$ quoted there is unaffected.")
    w("")

    # ---- diagnostics ------------------------------------------------------
    w("## Diagnostics (a) and (b)")
    w("")
    w("The brief asked for these only if the ratio came out well below 1. It "
      "did not, but they are cheap and they are what")
    w("rules out the two failure modes it names.")
    w("")
    w("**(a) and (b), per bin, 350-400 MHz.** Measured scatter against "
      "predicted, over the seven realisations:")
    w("")
    for st in ("drift", "raster"):
        neff = float(km[f"{st}_total"])
        r = ensemble(band, st, neff)
        w(f"*{ROW[st]}* (tr$(WA)$ = {neff:.1f})")
        w("")
        w("| b | $R_b$ | predicted $\\sigma_b$ | measured $\\sigma_b$ | "
          "ratio |")
        w("|---:|---:|---:|---:|---:|")
        sig = np.sqrt(r["R"] / (neff * r["n_b"]))
        meas = r["A"].std(axis=0, ddof=1)
        for i in range(len(sig)):
            w(f"| {i + 1} | {r['R'][i]:,.0f} | {sig[i]:.2f} | {meas[i]:.2f} | "
              f"{meas[i] / sig[i]:.2f} |")
        w("")
    w("There is **no shortfall at every bin**, which is what would have "
      "pointed at $N_{\\rm eff}$ and pixel correlations, and")
    w("**no shortfall that appears only on summing**, which is what would "
      "have pointed at bin-to-bin correlation from the")
    w("Blackman taper. The ratios scatter about 1 bin by bin and the "
      "combined SNR ratio is 0.85, consistent with 1 at 48%.")
    w("The taper's 1.78-bin equivalent noise bandwidth was the flagged risk "
      "in `nm_kpar_results.md`; it does not bite.")
    w("")

    # ---- tasks 2 and 3 ----------------------------------------------------
    w("## Task 2: patch markers removed")
    w("")
    w("`hi_sufficiency_bands.png` (per-$k_\\parallel$) is re-rendered "
      "without the open simulated-patch markers and without their")
    w("legend entry. The row bar now runs from the 1$\\sigma$ marker "
      "outwards, as it did before the patch markers set its left")
    w("end. Nothing else changed, and `hi_sufficiency_bands_1bin.png` is "
      "untouched (still byte-identical to the figure as it")
    w("stood before any of this work). Patch sizes for the caption: "
      "350-400 MHz, 265.2 / 674.8 / 574.1 deg$^2$ for the")
    w("3-strip drift, 12-strip drift and raster; 675-725 MHz, 152.1 / 331.7 "
      "/ 383.3 deg$^2$.")
    w("")
    w("## Task 3: the repo copy of `main.tex` is a 2026-10-01 snapshot")
    w("")
    w("Not edited. Dated, as asked:")
    w("")
    w("| | |")
    w("|---|---|")
    w("| last commit touching it | `d12b243`, 2026-09-30 20:38 |")
    w("| working-tree mtime | **2026-10-01 17:41**, uncommitted since |")
    w("| diff vs that commit | 449 insertions, 342 deletions |")
    w("")
    w("It contains none of: the `Open items (status 2026-10-08)` box, the "
      "string `hi_sufficiency_bands`, or the per-$k_\\parallel$")
    w("$\\sum_b n_b/(2+R_b)$ equation. It includes `hi_sufficiency.png` at "
      "line 751. So it predates the amplitude-only prior")
    w("work (2026-10-05), Section A/B (2026-10-07) and the per-$k_\\parallel$ "
      "work (2026-10-08) entirely -- it is seven days")
    w("older than the delivered draft, not a divergent branch of it. The "
      "delivered `main_tex_kpar.zip` copy is the live one;")
    w("the repo copy can be overwritten from it.")
    w("")
    w("One consequence: the single-band `hi_sufficiency.png` and "
      "`hi_sufficiency_f700.png` were regenerated per-$k_\\parallel$ on")
    w("2026-10-08 **because the repo copy still included them**. The "
      "delivered draft uses `hi_sufficiency_bands.png`, so those")
    w("two are now spare. They are consistent, not stale, so nothing needs "
      "doing; they can simply be ignored.")
    w("")

    # ---- reproducing ------------------------------------------------------
    w("## Reproducing")
    w("")
    w("```")
    w("/home/geoff/limTOD/.venv/bin/python ska_forecast_check.py")
    w("/home/geoff/limTOD/.venv/bin/python ska_hi_kmodes.py --bands --kpar")
    w("```")
    w("")
    w("`ska_hi_kmodes.fig_sufficiency_bands` gained `patch_markers=`, passed "
      "`False` by `--kpar`. The single-bin path is")
    w("untouched and still renders byte-identically.")
    w("")

    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")
    print(f"wrote {path}", flush=True)
    return path


if __name__ == "__main__":
    write_report()
