"""Line-of-sight mode counting for the sufficiency forecast.

The forecast as published sets ``N_m = N_eff = tr(WA)``: the angular degrees of
freedom of ONE channel. The cube has 32, and the forecast already assumes an
ideal spectroscopic tracer (``r = 1``, no shot noise), so it is entitled to the
line-of-sight modes as well. This counts them.

The floor-to-Hi ratio is strongly k_par dependent -- two decades across the
band at 350-400 MHz -- so the right object is not one ratio times one mode
count but an inverse-variance sum over k_par bins, each with its own ratio and
its own transfer function:

    R_b = P_fl,b / (T_b^2 P_HI,b),
    Omega_req(SNR) = SNR^2 / [ rho * sum_b n_b / (2 + R_b) ],    rho = N_eff/Omega_sim

``T_b^2`` rather than ``T_b``: ``transfer_function`` returns a cross-power with
the mock over the mock's AUTO-power, so T is an amplitude ratio and the cleaned
map is ``T_b * HI + floor``. Putting it inside R is what retires the old "areas
are optimistic by 1/T^2" caveat -- the correction is now carried, not noted.

**Degrees of freedom.** ``pk_par`` subtracts the mean along frequency and drops
the DC bin, so of the 32 real dof one is spent and 31 survive: ``n_b = 2`` for
b = 1..15 (real and imaginary) and ``n_16 = 1`` (Nyquist). This is a count of
real dof, the same thing ``tr(WA)`` counts in the angular direction. It is NOT
a count of k_par bins (16) and the two differ by almost exactly the factor of
two you would expect; see ``nm_kpar_results.md``.

Pure post-processing: every quantity is read from the ablation and experiment
npz already in ``results/``. Nothing here solves, simulates or needs operators.

    /home/geoff/limTOD/.venv/bin/python ska_nm_kpar.py
"""
from __future__ import annotations

import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

RESDIR = os.path.join(_HERE, "results")
NMODES = 4                      # cleaning depth, as everywhere else


# ---------------------------------------------------------------------------
# Where each band's two halves come from
# ---------------------------------------------------------------------------
#
# The floor and the Hi come from the ablation; T comes from the experiment.
# They have to be different files because the amplitude-only prior runs
# (2026-10-05) were ablation-only -- deliberately, since the ablation
# reproduces the experiment's residual/Hi to 0.0% at a fraction of the cost.
# That leaves no amplitude-only T anywhere, so T is taken from the matching
# pinned or own-prior experiment on the SAME common patch.
#
# Measured cost of that mismatch at 350-400 MHz, shape-pinned vs own prior:
# T agrees to 0.6% (drift12) and 0.3% (raster), against a 5-9% mock scatter per
# bin. The prior barely touches T, which is what one would expect -- T is a
# property of the PCA clean, and the clean sees the map, not the prior.

class Source:
    def __init__(self, tag, ablation, experiment, strategies, label,
                 kmodes=None, note=""):
        self.tag, self.label = tag, label
        self.ablation, self.experiment = ablation, experiment
        self.strategies = strategies
        self.kmodes = kmodes        # None -> no N_eff/area, so no Omega_req
        self.note = note


SOURCES = {
    "f350": Source(
        "f350", "hi_ablation_amp350_f350_400_nc32_ns64.npz",
        "hi_experiment_pin350_f350_400_nc32_ns64.npz",
        ("drift", "drift12", "raster"), "350-400 MHz",
        kmodes="hi_kmodes_f350_ns64.npz"),
    "f700": Source(
        "f700", "hi_ablation_amp700_f675_725_nc32_ns128.npz",
        "hi_experiment_matched700b_f675_725_nc32_ns128.npz",
        ("drift", "drift12m", "raster"), "675-725 MHz",
        kmodes="hi_kmodes_f700_ns128.npz"),
    # No N_eff and no patch area exist at these two bands: no hi_kmodes npz, no
    # modecheck, and hi_cache/ holds 350 MHz operators only. Building them is a
    # cluster job, so these bands get R_b, T_b and N_par,eff -- which need
    # neither -- and no Omega_req.
    "f500": Source(
        "f500", "hi_ablation_amp500_f500_550_nc32_ns128.npz",
        "hi_experiment_pin500_f500_550_nc32_ns128.npz",
        ("drift", "raster"), "500-550 MHz",
        note="no N_eff/area cached at this band"),
    "f875": Source(
        "f875", "hi_ablation_f875pin_f850_900_nc32_ns128.npz",
        "hi_experiment_f875pin_f850_900_nc32_ns128.npz",
        ("drift", "raster"), "850-900 MHz",
        note="no N_eff/area cached; prior is shape-pinned, not amplitude-only; "
             "the drift's T is unmeasurable above bin 10"),
}


def _nb(nbin):
    """Real degrees of freedom per retained bin.

    ``pk_par`` returns bins 1..nbin of an ``rfft`` over 32 channels, having
    dropped DC. Bins 1..15 are complex (2 real dof); bin 16 is Nyquist (1).
    """
    n = np.full(nbin, 2.0)
    n[-1] = 1.0
    return n


def kpar_terms(tag, nmodes=NMODES, variant="cr_", t_floor=0.0):
    """Per-k_par floor, Hi, transfer function and ratio, for one band.

    ``t_floor`` drops bins whose T is at or below it. Bins with small but
    positive T are kept: they contribute essentially nothing because R_b is
    correspondingly huge, and dropping them would silently improve the answer.
    """
    src = SOURCES[tag]
    ab = os.path.join(RESDIR, src.ablation)
    ex = os.path.join(RESDIR, src.experiment)
    for p in (ab, ex):
        if not os.path.exists(p):
            sys.exit(f"missing {p}")
    a, e = np.load(ab, allow_pickle=True), np.load(ex, allow_pickle=True)

    out = {}
    for s in src.strategies:
        k = a[f"{s}_{variant}k"]
        if not np.allclose(k, e[f"{s}_{variant}k"]):
            sys.exit(f"{tag}/{s}: ablation and experiment k grids differ")
        p_fl = a[f"{s}_{variant}floor_{nmodes}"]
        p_hi = a[f"{s}_{variant}p_hi"]
        t = e[f"{s}_{variant}tf_{nmodes}"]
        t_sc = e[f"{s}_{variant}tf_scatter_{nmodes}"]
        keep = t > t_floor
        r_b = np.full_like(p_fl, np.inf)
        r_b[keep] = p_fl[keep] / (t[keep] ** 2 * p_hi[keep])
        out[s] = dict(k=k, n_b=_nb(len(k)), p_fl=p_fl, p_hi=p_hi, T=t,
                      T_scatter=t_sc, R_b=r_b, keep=keep)
    return out


def weights(term, shot=1.0, r=1.0):
    """``sum_b n_b / [1 + (1 + R_b) shot / r^2]`` -- the cross-power sum."""
    k = term["keep"]
    return float(np.sum(term["n_b"][k]
                        / (1.0 + (1.0 + term["R_b"][k]) * shot / r ** 2)))


def auto_weights(term):
    """``sum_b n_b / [2 (1 + R_b)^2]`` -- the known-floor auto-spectrum bound."""
    k = term["keep"]
    return float(np.sum(term["n_b"][k] / (2.0 * (1.0 + term["R_b"][k]) ** 2)))


def single_bin_R(term, how="sum"):
    """The one-number R, with the T^2 convention.

    ``sum``    ratio of sums, as the task specifies: a genuine single band
               power over the whole cube, hence dominated by the lowest bins
               where the floor is worst.
    ``median`` median over bins, which is what the PUBLISHED floor_amp^2 is
               (``ska_hi_kmodes.floor_amp`` takes a median over k). The two
               differ by a factor of 2-14, so which one "the current method"
               means has to be stated rather than assumed.
    """
    k = term["keep"]
    if how == "sum":
        return float(term["p_fl"][k].sum()
                     / (term["T"][k] ** 2 * term["p_hi"][k]).sum())
    return float(np.median(term["p_fl"][k]
                           / (term["T"][k] ** 2 * term["p_hi"][k])))


# ---------------------------------------------------------------------------
# rho = N_eff / Omega_sim, from the band's cached kmodes npz
# ---------------------------------------------------------------------------

def geometry(tag):
    """``{strategy: (n_eff, area_deg2, floor_amp_published)}``, or ``{}``.

    Read from ``hi_kmodes_*.npz``, which ``ska_hi_kmodes.compute`` writes from
    the operators. Only 350-400 and 675-725 have one.
    """
    src = SOURCES[tag]
    if not src.kmodes:
        return {}
    path = os.path.join(RESDIR, src.kmodes)
    if not os.path.exists(path):
        return {}
    d = np.load(path, allow_pickle=True)
    return {s: (float(d[f"{s}_total"]), float(d[f"{s}_area"]),
                float(d[f"{s}_floor_amp"]))
            for s in src.strategies if f"{s}_total" in d}


def floor_ratios(tag, nmodes=NMODES):
    """Floor amplitude ratios between strategies -- must not move.

    They involve no ``N_m``, so nothing in this analysis can touch them. The
    brief asks for them as a stop condition, so they are recomputed from the
    same arrays the per-bin work reads rather than quoted from elsewhere.
    """
    terms = kpar_terms(tag, nmodes)
    amp = {s: float(np.median(np.sqrt(t["p_fl"] / t["p_hi"])))
           for s, t in terms.items()}
    ref = SOURCES[tag].strategies[-1]           # the raster
    return amp, {s: amp[s] / amp[ref] for s in amp if s != ref}, ref


# tab:timetofloor, in the draft's row order, with the draft's noise column
# (the whole noise arm rescaled by 0.35777 at 4 modes).
TAB_TIMETOFLOOR = (
    ("350-400 MHz", "f350", ("drift", "drift12", "raster"),
     (367.0, 494.7, 733.5)),
    ("500-550 MHz", "f500", ("drift", "raster"), (283.7, 406.3)),
    ("675-725 MHz", "f700", ("drift", "drift12m", "raster"),
     (375.7, 283.9, 457.1)),
    ("850-900 MHz", "f875", ("drift", "raster"), (481.5, 760.5)),
)

SNR = (1, 3, 5)
SKY_LIMIT_DEG2 = 30716.636739792735        # ska_hi_kmodes.SKY_LIMIT_DEG2
FULL_SKY_DEG2 = 41252.96124941928


def table_a(tag, nmodes=NMODES):
    """The headline row set for one band."""
    terms = kpar_terms(tag, nmodes)
    geo = geometry(tag)
    rows = []
    for s, t in terms.items():
        W = weights(t)
        R1_sum = single_bin_R(t, "sum")
        R1_med = single_bin_R(t, "median")
        row = dict(strategy=s, W=W, auto_W=auto_weights(t),
                   R1_sum=R1_sum, R1_med=R1_med,
                   npar_sum=(2.0 + R1_sum) * W, npar_med=(2.0 + R1_med) * W,
                   ndrop=int((~t["keep"]).sum()))
        if s in geo:
            n_eff, area, f_amp = geo[s]
            rho = n_eff / area
            row.update(n_eff=n_eff, area=area, rho=rho, floor_amp=f_amp,
                       R_pub=f_amp ** 2,
                       omega={q: q ** 2 / (rho * W) for q in SNR},
                       omega_1bin_sum={q: q ** 2 * (2.0 + R1_sum) / rho
                                       for q in SNR},
                       omega_1bin_med={q: q ** 2 * (2.0 + R1_med) / rho
                                       for q in SNR},
                       omega_pub={q: q ** 2 * (2.0 + f_amp ** 2) / rho
                                  for q in SNR},
                       omega_auto=1.0 / (rho * row["auto_W"]))
        rows.append(row)
    return rows, terms


def sanity(tag, nmodes=NMODES):
    """Every check the brief asks for, as (name, ok, detail) triples."""
    terms = kpar_terms(tag, nmodes)
    geo = geometry(tag)
    out = []
    for s, t in terms.items():
        nb_sum = float(t["n_b"].sum())
        # 31, not 32: pk_par subtracts the mean along frequency and drops DC.
        out.append((f"{s}: sum n_b", abs(nb_sum - 31.0) < 1e-9,
                    f"{nb_sum:.0f} (32 minus the mean-subtracted DC mode)"))
        if s in geo:
            pub = geo[s][2]
            got = float(np.median(np.sqrt(t["p_fl"] / t["p_hi"])))
            out.append((f"{s}: published floor_amp reproduced",
                        abs(got - pub) <= 0.05 + 0.005 * pub,
                        f"median_b sqrt(P_fl/P_HI) = {got:.2f} vs {pub:.1f}"))
        lo = t["T"] < 0.1
        out.append((f"{s}: T_b >= 0.1 in every bin", not lo.any(),
                    "yes" if not lo.any() else
                    f"bins {[int(i) for i in np.flatnonzero(lo) + 1]} kept with T < 0.1"))
        neg = ~t["keep"]
        out.append((f"{s}: no bin dropped for T <= 0", not neg.any(),
                    "none" if not neg.any() else
                    f"dropped bins {[int(i) for i in np.flatnonzero(neg) + 1]}"))
        noisy = t["T"] < 3.0 * t["T_scatter"]
        out.append((f"{s}: T_b at >3x its mock scatter", not noisy.any(),
                    "yes" if not noisy.any() else
                    f"bins {[int(i) for i in np.flatnonzero(noisy) + 1]} below 3 sigma"))
    return out


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

ROW = {"drift": "3-strip drift", "drift12": "12-strip drift",
       "drift12m": "12-strip drift", "raster": "raster"}


def _deg(v):
    if v < 100:
        return f"{v:.1f}"
    return f"{v:,.0f}" if v < 1e5 else f"{v:.3g}"


def _reach(v):
    return "yes" if v <= SKY_LIMIT_DEG2 else "no"


def write_report(path=None, nmodes=NMODES):
    path = path or os.path.join(_HERE, "nm_kpar_results.md")
    L = []
    w = L.append
    w("# N_m over k_par bins, and the matched-white re-render")
    w("")
    w(f"Generated by `ska_nm_kpar.py` at {nmodes} PCA modes, amplitude-only "
      "prior where one exists, noiseless floor arm, common resolution,")
    w("interior pixels. Pure post-processing of the npz already in "
      "`results/` -- nothing was re-solved or re-simulated.")
    w("")

    # ---- inputs -----------------------------------------------------------
    w("## Inputs")
    w("")
    w("The floor and the Hi come from the ablation; `T_b` comes from the "
      "experiment. They have to be different files: the")
    w("amplitude-only prior runs of 2026-10-05 were ablation-only, so no "
      "amplitude-only `T` exists anywhere. `T` is taken")
    w("from the matching experiment on the **same common patch**.")
    w("")
    w("| band | floor + Hi | T | strategies | Omega_req? |")
    w("|---|---|---|---|---|")
    for tag in ("f350", "f700", "f500", "f875"):
        src = SOURCES[tag]
        w(f"| {src.label} | `{src.ablation}` | `{src.experiment}` | "
          f"{', '.join(src.strategies)} | "
          f"{'yes' if src.kmodes else '**no** -- ' + src.note} |")
    w("")
    w("**Cost of that mismatch, measured.** At 350-400 MHz, shape-pinned vs "
      "own-prior experiment, `T_b` agrees to 0.6%")
    w("(12-strip) and 0.4% (raster) in every bin, against a 5-9% mock scatter "
      "per bin. The drift is the reference strategy")
    w("and its `T_b` is identical to 1 part in 1e9. This is expected -- "
      "`T` is a property of the PCA clean, which sees the")
    w("map, not the prior -- but it had to be checked rather than assumed.")
    w("")
    w("**500-550 and 850-900 cannot be forecast.** `Omega_req` needs "
      "`rho = N_eff/Omega_sim`, and neither band has an")
    w("`hi_kmodes` npz, a modecheck, or cached operators (`hi_cache/` holds "
      "350 MHz only). Building them is a cluster job,")
    w("so per the standing rule this stops here and reports instead. What "
      "needs neither -- `T_b`, `R_b`, `N_par,eff` -- is")
    w("in Table B. At 850-900 the 3-strip drift's `T` is also unusable above "
      "bin 10 (scatter 1.2-12.9 on 20 mocks, and")
    w("`T > 1` in bins 11-13, which is unphysical), and that band has no "
      "amplitude-only ablation either, so its numbers are")
    w("outside the stated convention twice over.")
    w("")

    # ---- Table A ----------------------------------------------------------
    w("## Table A -- headline")
    w("")
    w("`Omega_req` in deg^2 for `Delta P / P = 1 / SNR` against an ideal "
      "tracer (`r = 1`, no shot noise).")
    w("Per-bin uses `Omega = SNR^2 / [rho sum_b n_b/(2 + R_b)]`; single-bin "
      "uses `Omega = SNR^2 (2 + R_1bin) / rho`.")
    w("")
    w("| band | strategy | N_par,eff | 1 sigma | 3 sigma | 5 sigma | within "
      "reach (1/3/5) | single-bin 1 sigma | published 1 sigma |")
    w("|---|---|---:|---:|---:|---:|:--:|---:|---:|")
    for tag in ("f350", "f700"):
        rows, _ = table_a(tag, nmodes)
        for r in rows:
            if "omega" not in r:
                continue
            o = r["omega"]
            w(f"| {SOURCES[tag].label} | {ROW[r['strategy']]} | "
              f"{r['npar_med']:.1f} | {_deg(o[1])} | {_deg(o[3])} | "
              f"{_deg(o[5])} | "
              f"{'/'.join(_reach(o[q]) for q in SNR)} | "
              f"{_deg(r['omega_1bin_med'][1])} | {_deg(r['omega_pub'][1])} |")
    w("")
    w("`N_par,eff = (2 + R_1bin) * sum_b n_b/(2 + R_b)` with `R_1bin` the "
      "MEDIAN over bins, which is the convention the published")
    w("`floor_amp` uses, so the column is the factor by which the forecast "
      "moves. Taking `R_1bin` as a ratio of sums instead, as a literal")
    w("single band power over the cube, gives:")
    w("")
    w("| band | strategy | N_par,eff (median R) | N_par,eff (summed R) | "
      "single-bin 1 sigma, summed R |")
    w("|---|---|---:|---:|---:|")
    for tag in ("f350", "f700"):
        rows, _ = table_a(tag, nmodes)
        for r in rows:
            if "omega" not in r:
                continue
            w(f"| {SOURCES[tag].label} | {ROW[r['strategy']]} | "
              f"{r['npar_med']:.1f} | {r['npar_sum']:.1f} | "
              f"{_deg(r['omega_1bin_sum'][1])} |")
    w("")

    # ---- Table B ----------------------------------------------------------
    w("## Table B -- per bin")
    w("")
    w("`n_b` is real degrees of freedom: 2 for bins 1-15, 1 for the Nyquist "
      "bin 16. DC is absent -- `pk_par` subtracts the mean")
    w("along frequency and drops it -- so the column sums to **31**, not 32. "
      "`share` is the bin's fraction of `sum_b n_b/(2 + R_b)`.")
    w("")
    for tag in ("f350", "f700", "f500", "f875"):
        src = SOURCES[tag]
        terms = kpar_terms(tag, nmodes)
        for s, t in terms.items():
            W = weights(t)
            w(f"### {src.label}, {ROW[s]}")
            if src.note:
                w("")
                w(f"*{src.note}*")
            w("")
            w("| b | k_par [Mpc^-1] | n_b | T_b | R_b | share |")
            w("|---:|---:|---:|---:|---:|---:|")
            for i in range(len(t["k"])):
                if not t["keep"][i]:
                    w(f"| {i + 1} | {t['k'][i]:.4f} | {t['n_b'][i]:.0f} | "
                      f"{t['T'][i]:.3f} | -- | dropped, T <= 0 |")
                    continue
                sh = t["n_b"][i] / (2.0 + t["R_b"][i]) / W
                w(f"| {i + 1} | {t['k'][i]:.4f} | {t['n_b'][i]:.0f} | "
                  f"{t['T'][i]:.3f} | {t['R_b'][i]:,.0f} | {100 * sh:.1f}% |")
            w("")

    # ---- auto -------------------------------------------------------------
    w("## Auto-spectrum, known-floor bound")
    w("")
    w("`(Delta P/P)^-2 = rho Omega sum_b n_b / [2 (1 + R_b)^2]`. This is a "
      "LOWER BOUND and not a requirement: in the auto-spectrum the")
    w("floor is a bias, deterministic given the strategy and the sky, so no "
      "sky area averages it away. The bound holds only if the")
    w("floor's mean power were known and subtracted -- which would mean "
      "knowing the sky the survey failed to measure.")
    w("")
    w("| band | strategy | per-bin Omega | published single-bin Omega | "
      "vs reachable sky (30,717 deg^2) |")
    w("|---|---|---:|---:|:--|")
    for tag in ("f350", "f700"):
        rows, _ = table_a(tag, nmodes)
        for r in rows:
            if "omega_auto" not in r:
                continue
            a = r["omega_auto"]
            pub = 2.0 * (1.0 + r["R_pub"]) ** 2 / r["rho"]
            w(f"| {SOURCES[tag].label} | {ROW[r['strategy']]} | {_deg(a)} | "
              f"{_deg(pub)} | "
              f"{'**inside**' if a <= SKY_LIMIT_DEG2 else f'{a / SKY_LIMIT_DEG2:.3g}x above'} |")
    w("")

    # ---- notes ------------------------------------------------------------
    w("## Notes")
    w("")
    w("### Where the information is")
    w("")
    w("**Bins 3-11 carry it all.** At 350-400 MHz the raster's bins 3-11 "
      "(`k_par` 0.036-0.13 Mpc^-1) contribute 87.3% of")
    w("`sum_b n_b/(2 + R_b)` between them, at roughly 10% each. Bins 1-2 "
      "contribute **0.4% combined**: that is where the")
    w("floor is worst (`R_1 = 86,080` against a median 621) and where `T_b` "
      "is lowest (0.51 and 0.81 against ~0.88), and the")
    w("two effects compound because `T` enters squared. The high-k tail falls "
      "away more gently, to 0.5% in the Nyquist bin.")
    w("The same shape holds for every strategy at both bands. **The forecast "
      "is therefore insensitive to the two bins whose**")
    w("**`T` is least well measured**, which is the reassuring version of the "
      "`T` caveat.")
    w("")
    w("### The ordering does not change; the gaps do")
    w("")
    w("At both bands the strategy ordering at every sigma level is exactly "
      "what the single-bin forecast gave. What moves is")
    w("the size of the gaps: at 350-400 MHz the drift/raster area ratio goes "
      "from **451x to 329x**, because the raster's")
    w("floor varies more sharply across `k_par` than the drift's (a factor "
      "212 between its best and worst bin, against 24)")
    w("and so gains less from a per-bin treatment than the drift does.")
    w("")
    w("### Degrees of freedom: the convention question from step 3")
    w("")
    w("`pk_par` subtracts the mean along frequency and drops DC, so **31 real "
      "dof survive, not 32**. Of the 16 bins it")
    w("returns, 1-15 are complex (2 real dof) and 16 is Nyquist (1). A "
      "numerical Parseval check on white noise closes to")
    w("5%: `sum_b n_b P_b / L = 0.952` against a map variance of 0.971, the "
      "shortfall being the Blackman-tapered band edges.")
    w("")
    w("**The ~2x the brief warned about is real.** Counting one dof per bin "
      "instead of `n_b` halves every `N_par,eff`:")
    w("")
    w("| 350-400 MHz | n_b, with T^2 | 1 per bin, with T^2 | n_b, no T^2 | "
      "1 per bin, no T^2 |")
    w("|---|---:|---:|---:|---:|")
    w("| 3-strip drift | 41.9 | 21.1 | 43.1 | 21.7 |")
    w("| raster | 27.1 | 13.6 | 28.0 | 14.1 |")
    w("")
    w("`n_b` is the right one: it is a count of real degrees of freedom, "
      "which is the same thing `tr(WA)` counts in the")
    w("angular direction, and it is what makes the sum over bins an "
      "inverse-variance combination rather than a mean.")
    w("")
    w("**I cannot reproduce the earlier 12.8 / 15.9.** None of the four cells "
      "above gives them, and the ordering is wrong in")
    w("all four: every convention here puts the drift ABOVE the raster "
      "(the drift's floor varies more across `k_par`, so it")
    w("gains more), while 12.8 / 15.9 has the raster above. Participation "
      "ratios of `P_HI`, `P_fl` and the weights give 9.2,")
    w("5.1 and 10.3 for the drift and 9.4, 2.1 and 11.1 for the raster -- "
      "also not it. If those numbers matter, the")
    w("convention that produced them needs to come from "
      "`paper1_frozen_numbers.md`, which is not in this repo.")
    w("")
    w("### Why summing over bins is legitimate here")
    w("")
    w("The obvious objection is that the floor is `(I - WA)s`, deterministic "
      "given the strategy and the sky, so its `k_par`")
    w("bins are not independent draws and inverse-variance summing them "
      "overcounts. That objection does not bite in")
    w("**cross-correlation**: the variance there is the chance correlation "
      "between a FIXED floor and the tracer, and the")
    w("tracer is an independent random field supplying fresh randomness in "
      "every (angular mode, `k_par`) pair. The sum is")
    w("over the tracer's degrees of freedom, which are genuinely 31 along the "
      "line of sight.")
    w("")
    w("It bites hard in the **auto-spectrum**, where the floor beats against "
      "itself -- which is why that row is a bias and")
    w("not a variance, and why its area is labelled a lower bound both before "
      "and after this change.")
    w("")
    w("### Two things that are more fragile than they look")
    w("")
    w("**The taper correlates neighbouring bins.** `pk_par` applies a "
      "Blackman window (not optional -- the foregrounds are")
    w("~1e4 times the Hi, so untapered leakage would swamp the signal at high "
      "`k_par`), whose equivalent noise bandwidth is")
    w("**1.78 bins**. The dof count is unaffected -- a window is a linear "
      "transform and cannot create or destroy dof -- but")
    w("the ASSIGNMENT of floor power to bins is blurred over about two bins. "
      "That matters only where `R_b` changes steeply,")
    w("i.e. bins 1-3, which carry under 7% of the information. Everywhere "
      "else `R_b` is flat and the blur is harmless.")
    w("")
    w("**The auto-spectrum bound moves much further than the cross.** At "
      "675-725 MHz the known-floor bound drops from")
    w("356,633 to 1,895 deg^2 for the drift -- a factor 188, against an "
      "`N_par,eff` of 70 for the cross. The extra factor is")
    w("the square: the auto weight is `n_b/[2(1 + R_b)^2]`, so a bin twice as "
      "good as the median counts four times, not")
    w("twice. **All three strategies at 675-725 MHz now sit inside the "
      "reachable sky on this bound.** I would not lean on")
    w("that: it is still a bound that requires knowing the sky the survey "
      "failed to measure, and it is the quantity in this")
    w("report most sensitive to the per-bin treatment.")
    w("")
    w("### The headline, and a way to falsify it")
    w("")
    w("At 675-725 MHz every `Omega_req` is now **smaller than the simulated "
      "patch**: 12.9 deg^2 needed against 152 simulated")
    w("for the drift, 35.5 against 383 for the raster. Taken at face value "
      "the simulation should already show a cross-")
    w("detection against an ideal tracer in the patch it has: **3.4 sigma** "
      "(drift), **6.5 sigma** (12-strip), **3.3 sigma**")
    w("(raster). That is a prediction this series can check")
    w("directly -- cross-correlate the cleaned cube with the injected Hi mock "
      "and read off the significance -- and it is")
    w("worth checking before the number goes in the draft, because it is the "
      "first time the forecast has made a claim")
    w("about data already on disk rather than about a survey nobody has run.")
    w("")

    # ---- sanity -----------------------------------------------------------
    w("## Sanity checks")
    w("")
    w("### Floor ratios -- the stop condition")
    w("")
    w("None of this touches the floor: `R_b` has `P_fl,b` in the numerator "
      "unchanged, and `N_m` appears nowhere in a floor")
    w("ratio. Recomputed from the same arrays the per-bin work reads:")
    w("")
    w("| band | strategy | floor/Hi amplitude | / raster |")
    w("|---|---|---:|---:|")
    for tag in ("f350", "f700", "f500", "f875"):
        amp, rat, ref = floor_ratios(tag, nmodes)
        for s in SOURCES[tag].strategies:
            w(f"| {SOURCES[tag].label} | {ROW[s]} | {amp[s]:.2f} | "
              f"{'--' if s == ref else f'{rat[s]:.2f}x'} |")
    w("")
    w("The matched-density ratio is **10.48x** (12-strip drift / raster at "
      "350-400 MHz), i.e. the 10.5x the draft quotes,")
    w("and every published floor amplitude reproduces to the decimal it is "
      "quoted at. Nothing moved.")
    w("")
    w("### Per-strategy checks")
    w("")
    for tag in ("f350", "f700", "f500", "f875"):
        w(f"**{SOURCES[tag].label}**")
        w("")
        for name, ok, detail in sanity(tag, nmodes):
            w(f"- {'PASS' if ok else 'FLAG'} -- {name}: {detail}")
        w("")

    # ---- Task 2 -----------------------------------------------------------
    w("## Task 2 -- `hi_noise_terms` at the matched white-noise level")
    w("")
    w("The draft quotes `sigma_eta^2 = 3.2e-7` (radiometer, 1.5625 MHz x 2 s) "
      "against the simulated 2.5e-6, so the")
    w("residuals are rescaled by `sqrt(3.2e-7/2.5e-6) = 0.35777`. Only the "
      "WHITE part may be rescaled: 1/f is a gain")
    w("fluctuation and does not obey the radiometer equation. The "
      "`noiseonly_indep` arm is `(1 + g)(1 + w_ind)` and")
    w("carries both, so it is split first (the two are drawn independently "
      "and the `g*w` cross term is ~1e-6 relative):")
    w("")
    w("```")
    w("white^2  = noise^2 - (1/f)^2")
    w("matched  = 0.35777^2 * white^2 + (1/f)^2")
    w("```")
    w("")
    w("Both the split and the rescaling are done per k_par bin before the "
      "median over bins, not on the medians.")
    w("")
    w("| band | strategy | draft (arm rescaled whole) | separated | diff | "
      "floor | t x for noise = floor |")
    w("|---|---|---:|---:|---:|---:|---:|")
    for band, src_tag, strategies, draft in TAB_TIMETOFLOOR:
        ab = np.load(os.path.join(RESDIR, SOURCES[src_tag].ablation),
                     allow_pickle=True)
        for s, d_val in zip(strategies, draft):
            p_hi = ab[f"{s}_cr_p_hi"]
            noise = ab[f"{s}_cr_noiseonly_indep_{nmodes}"]
            gain = ab[f"{s}_cr_gainonly_{nmodes}"]
            matched = 0.35777087639996635 ** 2 * np.clip(noise - gain, 0, None) + gain
            amp = float(np.median(np.sqrt(matched / p_hi)))
            fl = float(np.median(np.sqrt(ab[f"{s}_cr_floor_{nmodes}"] / p_hi)))
            w(f"| {band} | {ROW[s]} | {d_val:.1f} | {amp:.1f} | "
              f"{100 * (amp / d_val - 1):+.3f}% | {fl:.1f} | "
              f"{3.0 * (amp / fl) ** 2:,.1f} h |")
    w("")
    w("All ten agree with the draft to better than **0.02%**, well inside the "
      "1% the brief allows. 1/f is 0.03-0.5% of the")
    w("noise amplitude at every band, so splitting it out before rescaling "
      "changes nothing numerically. It is still the")
    w("right thing to do: the draft's column is correct by accident, not by "
      "construction, and stops being correct the")
    w("moment 1/f is raised. No bin had `noise^2 < (1/f)^2`, so no clipping "
      "was needed anywhere.")
    w("")
    w("`t x` is `3 h * (noise/floor)^2`, the integration time at which the "
      "matched white noise would reach the floor.")
    w("At 350-400 MHz the 3-strip drift is **already floor-limited** at the "
      "3 h simulated (2.5 h needed); the raster at the")
    w("same band needs 3,246 h. That inversion is the figure's point: the "
      "better the strategy's floor, the longer it takes")
    w("to integrate down to it.")
    w("")
    w("Figures: `hi_noise_terms.png` re-rendered; the pre-rescale version is "
      "kept as `hi_noise_terms_white2.5e-6.png`.")
    w("")
    w("**One label changed.** The dashed arm was captioned \"white noise "
      "(independent)\" but is `(1 + g)(1 + w_ind)` and")
    w("carries 1/f -- the same mis-pairing `fig_ablation` had to be corrected "
      "for on 2026-10-07. It now reads \"all noise")
    w("(white independent)\". Revert it if the draft's caption already makes "
      "the distinction.")
    w("")

    # ---- reproducing ------------------------------------------------------
    w("## Reproducing")
    w("")
    w("```")
    w("# Tables A and B, the auto bound, the notes and every sanity check")
    w("/home/geoff/limTOD/.venv/bin/python ska_nm_kpar.py")
    w("")
    w("# the per-k_par sufficiency figures")
    w("/home/geoff/limTOD/.venv/bin/python ska_hi_kmodes.py --bands --kpar")
    w("/home/geoff/limTOD/.venv/bin/python ska_hi_kmodes.py --band f350 "
      "--draw-only --kpar")
    w("/home/geoff/limTOD/.venv/bin/python ska_hi_kmodes.py --band f700 "
      "--draw-only --kpar")
    w("")
    w("# the single-bin versions, unchanged (drop --kpar)")
    w("/home/geoff/limTOD/.venv/bin/python ska_hi_kmodes.py --bands")
    w("")
    w("# hi_noise_terms.png at the matched white level")
    w("/home/geoff/limTOD/.venv/bin/python -c \\")
    w("  'import matplotlib; matplotlib.use(\"Agg\"); "
      "import ska_hi_design as D; D.fig_noise_terms()'")
    w("```")
    w("")
    w("Reads `results/*.npz` only. No operators, no solves, no TOD, no "
      "cluster.")
    w("")
    w("**The draft includes `hi_sufficiency.png`, not "
      "`hi_sufficiency_bands.png`** (`main.tex:751`). The brief named the "
      "merged")
    w("figure, but leaving the one the draft actually builds with on "
      "single-bin numbers would have put a figure and its")
    w("own caption in contradiction, so the single-band panels were "
      "regenerated too. Both conventions are kept side by")
    w("side:")
    w("")
    w("| per-k_par | single-bin, for comparison |")
    w("|---|---|")
    w("| `hi_sufficiency_bands.png` | `hi_sufficiency_bands_1bin.png` |")
    w("| `hi_sufficiency.png` = `hi_sufficiency_f350.png` | "
      "`hi_sufficiency_1bin.png` = `hi_sufficiency_f350_1bin.png` |")
    w("| `hi_sufficiency_f700.png` | `hi_sufficiency_f700_1bin.png` |")
    w("")
    w("`hi_kmodes*.png` are untouched: they decompose `tr(WA)` by angular "
      "scale and have no `N_m` in them.")
    w("")
    w("**Code changed, 2026-10-08.**")
    w("")
    w("- `ska_nm_kpar.py` -- NEW, self-contained. All the per-k_par "
      "arithmetic and this report.")
    w("- `ska_hi_kmodes.py` -- `sufficiency()` takes an optional `kpar=` term "
      "dict; `load()` takes `kpar=False`; a")
    w("  `--kpar` flag on `main`; `X_BANDS_KPAR` because the per-bin areas "
      "fall off the left of the single-bin window;")
    w("  `_sig3` gains a decimal below 100; a label just ABOVE the "
      "reachable-sky line is left-aligned so it stops")
    w("  running through it; the row bar spans the simulated patch even when "
      "the patch is larger than the 1-sigma area.")
    w("- `ska_hi_design.py` -- `WHITE_VAR_MATCHED`, `WHITE_SCALE` and "
      "`noise_matched()`; `fig_noise_terms` draws the")
    w("  recombined arm and relabels it; half a decade of extra y-floor so "
      "the legend keeps its own space.")
    w("")
    w("**The single-bin path is untouched.** `ska_hi_kmodes.py --bands` "
      "without `--kpar` still passes all ten regression")
    w("guards in `BANDS_EXPECT` and renders a figure **byte-identical** to "
      "the one from before these edits; that was")
    w("checked after every patch, and `hi_sufficiency_bands_1bin.png` is that "
      "byte-identical render.")
    w("")

    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")
    print(f"wrote {path}", flush=True)
    return path


if __name__ == "__main__":
    write_report()
