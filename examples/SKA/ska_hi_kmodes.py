"""Where the measured sky modes live in angular scale -- and what they buy.

``tr(WA)`` says *how many* sky modes a strategy measures. It does not say
*which*, and the series has already shown that the count alone is not a
sufficient statistic: at 22 modes/100 deg^2 the cross-linked raster's floor is
9.3x below a twelve-strip parked drift's. This decomposes the count by angular
scale, which is the axis a science requirement is actually written on.

The decomposition is exact rather than a projection onto a non-orthogonal
basis. Eigen-decompose ``A^T N^-1 A`` into modes ``u_i`` with eigenvalues
``lambda_i``; each contributes ``w_i = lambda_i / (lambda_i + S^-1)`` to
``tr(WA)``, between 0 (prior-filled) and 1 (data-dominated). Every ``u_i`` is
then given an angular power spectrum ``f_i(k)`` normalised to sum to one over
the k bins, by direct (non-uniform) Fourier transform on the tangent plane --
no regridding, so no interpolation artefacts. Then

    n_eff(k) = sum_i w_i f_i(k),      sum_k n_eff(k) = tr(WA) exactly,

so the curve is a decomposition of the published mode count, not a new
statistic that has to be reconciled with it.

The second panel turns that into a sufficiency statement. For a k bin holding
``N_m`` independent modes, an auto-spectrum measures ``sigma_P/P =
sqrt(2/N_m)(1 + P_fl/P_HI)`` and a cross-correlation against a tracer
``sigma/P_x = sqrt((P_fl/P_HI)/N_m)/r``. The difference is a square, and it is
the whole story: in auto the floor is a **bias** -- deterministic given the
strategy, so no amount of sky averages it away -- while in cross it is
dominated by unmeasured *foreground*, which is uncorrelated with the tracer and
therefore contributes variance only.

    /home/geoff/limTOD/.venv/bin/python ska_hi_kmodes.py

Reads the cached operators in ``hi_cache/`` only; no TOD, no new runs.
"""
from __future__ import annotations

import os
import pickle
import sys

import healpy as hp
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import ska_hi_experiment as X
import ska_hi_plots as P
from ska_hi_plots import INK, INK2, MUTED, SURFACE, BASELINE, GRID
from ska_common import gdsm_equatorial_sky_model, ska_beam_fwhm_deg
from ska_hi_ladder import FREQ_MHZ, operator_list

# Comoving distance implied by the series' own anchor: a 3.98 deg beam at
# 350 MHz subtends 401 Mpc at z_eff = 2.80, i.e. k_perp = 2 pi / 401 = 0.0157
# (Sec. "Injected signal and geometry"). Deriving D_c this way rather than
# recomputing it guarantees this figure agrees with the number the paper quotes.
FWHM_RAD = np.radians(ska_beam_fwhm_deg(FREQ_MHZ))
D_COMOVING_MPC = (2 * np.pi / FWHM_RAD) / 0.0157

# floor/HI in AMPLITUDE, 4 modes removed, common resolution, matched run
FLOOR_AMP = {"drift": 399.6, "raster": 25.6}
R_CROSS = 1.0          # cross-correlation coefficient assumed for the estimate


def eig_decomposition(mm):
    """``lambda_i``, ``u_i`` and ``S^-1`` -- the recipe ``mode_count`` uses."""
    pix = np.asarray(mm.pixel_indices)
    truth = gdsm_equatorial_sky_model(freq=FREQ_MHZ, nside=X.NSIDE)[pix]
    nv_floor = X.WHITE_VAR * (1e-3 * float(np.mean(truth))) ** 2
    M = None
    for o in operator_list(mm):
        A = np.asarray(o, float)
        var = X.WHITE_VAR * (A @ truth) ** 2 + nv_floor
        term = A.T @ (A / var[:, None])
        M = term if M is None else M + term
    lam, U = np.linalg.eigh(M)
    lam = np.maximum(lam[::-1], 0.0)
    U = U[:, ::-1]
    s_inv = 1.0 / max(float(np.std(truth)), 1e-3) ** 2
    return lam, U, s_inv, pix


def tangent_xy(pix):
    """Tangent-plane offsets in radians about the patch centre."""
    th, ph = hp.pix2ang(X.NSIDE, pix)
    ra, dec = np.unwrap(ph), np.pi / 2 - th
    ra0, dec0 = ra.mean(), dec.mean()
    return (ra - ra0) * np.cos(dec0), dec - dec0


def n_eff_of_k(mm, kedges, n_ang=24):
    """``n_eff`` per |k| bin. Sums to ``tr(WA)`` by construction."""
    lam, U, s_inv, pix = eig_decomposition(mm)
    w = lam / (lam + s_inv)                       # (nmode,), sums to tr(WA)
    x, y = tangent_xy(pix)
    kmid = np.sqrt(kedges[:-1] * kedges[1:])
    angles = np.linspace(0.0, np.pi, n_ang, endpoint=False)

    Pk = np.zeros((len(kmid), U.shape[1]))
    for b, kk in enumerate(kmid):
        acc = np.zeros(U.shape[1])
        for a in angles:
            phase = np.exp(-1j * (kk * np.cos(a) * x + kk * np.sin(a) * y))
            acc += np.abs(phase @ U) ** 2
        Pk[b] = acc / n_ang
    # normalise each mode's spectrum to unit total, so the k sum is exact
    Pk /= np.maximum(Pk.sum(axis=0, keepdims=True), 1e-300)
    return kmid, Pk @ w, float(w.sum())


def sufficiency(floor_amp, n_eff, area):
    """Modes and sky area needed for sigma/P = 1, auto and cross."""
    ratio = floor_amp ** 2                        # P_floor / P_HI
    need = {"auto": 2.0 * (1.0 + ratio) ** 2,
            "cross": ratio / R_CROSS ** 2}
    return {k: dict(modes=v, factor=v / n_eff, area=area * v / n_eff)
            for k, v in need.items()}


def fig_kmodes(res, text=None):
    """Where the modes are, and what they are enough for."""
    fig, (axA, axB) = plt.subplots(
        1, 2, figsize=P._figsize("hi_kmodes", (13.2, 5.4)),
        gridspec_kw=dict(width_ratios=[1.0, 1.05], wspace=0.34))

    # --- left: the mode count, decomposed by angular scale -----------------
    P._tidy(axA); axA.set_xscale("log"); axA.set_yscale("log")
    for s_ in ("drift", "raster"):
        st = P.MONO_STRAT[s_]
        r = res[s_]
        dens = 100.0 * r["n"] / r["area"]         # per 100 deg^2, the currency
        axA.plot(r["k"], dens, color=st["color"], ls=st["ls"],
                 marker=st["marker"], ms=6, mfc=SURFACE, mew=1.6, zorder=4,
                 label=f"{s_} \u2014 {100 * r['total'] / r['area']:.1f} total")
    axA.legend(loc="upper left", fontsize=10.2)
    k_beam = 2 * np.pi / FWHM_RAD
    axA.axvline(k_beam, color=INK, lw=1.4, ls=(0, (4, 3)), zorder=3)
    axA.annotate("beam", (k_beam, axA.get_ylim()[0]),
                 textcoords="offset points", xytext=(-6, 8), ha="right",
                 va="bottom", color=INK, fontsize=9.9, rotation=90)
    axA.set_xlabel("$|k_\\perp|$  [rad$^{-1}$]   $\\simeq \\ell$")
    axA.set_ylabel("measured modes per 100 deg$^2$, per bin")
    sec = axA.secondary_xaxis(
        "top", functions=(lambda k: k / D_COMOVING_MPC,
                          lambda kp: kp * D_COMOVING_MPC))
    sec.set_xlabel("$k_\\perp$  [Mpc$^{-1}$]  at $z = 2.80$", fontsize=10.6)

    # --- right: what that is enough for ------------------------------------
    P._tidy(axB); axB.set_xscale("log"); axB.grid(axis="y", visible=False)
    rows = [("raster", "cross"), ("drift", "cross"),
            ("raster", "auto"), ("drift", "auto")]
    ys = np.arange(len(rows))[::-1]
    for y, (s_, kind) in zip(ys, rows):
        a = res[s_]["suff"][kind]["area"]
        ok = a <= 41253.0
        axB.plot([res[s_]["area"], a], [y, y], color=INK if ok else MUTED,
                 lw=2.2, zorder=3, solid_capstyle="round")
        axB.plot([res[s_]["area"]], [y], marker=P.MONO_STRAT[s_]["marker"],
                 ms=8, mfc=SURFACE, mec=INK, mew=1.8, zorder=5)
        axB.plot([a], [y], marker="*", ms=16,
                 color=INK if ok else MUTED, zorder=5)
        axB.annotate(f"{a:,.0f} deg$^2$" if a < 1e5 else f"{a:.1g} deg$^2$",
                     (a, y), textcoords="offset points", xytext=(10, 0),
                     va="center", ha="left", color=INK if ok else MUTED,
                     fontsize=10.2, fontweight="semibold" if ok else "normal")
    axB.axvline(41253.0, color=INK, lw=1.6, zorder=4)
    axB.annotate("whole sky", (41253.0, len(rows) - 0.45),
                 textcoords="offset points", xytext=(-6, 0), ha="right",
                 va="center", color=INK, fontsize=9.9)
    axB.set_yticks(ys)
    axB.set_yticklabels([f"{s_}\n{kind}" for s_, kind in rows], fontsize=10.2)
    axB.set_ylim(-0.7, len(rows) - 0.25)
    axB.set_xlim(80.0, 1e9)
    axB.set_xlabel("sky area needed for $\\sigma_P/P = 1$   "
                   "(marker = observed now)")
    '''axB.set_title("Only one of the four is a survey you could build",
                  color=INK2, fontsize=11.0, loc="left", pad=6)'''

    '''P._title(fig, "How many modes is enough? Only the cross-correlation asks a "
             "reachable number",
             "Left: tr($WA$) split by angular scale \u2014 the k sum reproduces "
             "the published 14.6 and 22.3 per 100 deg$^2$ exactly.\n"
             "Right: in an auto-spectrum the floor is a BIAS that sky area "
             "cannot average away; against a tracer it is uncorrelated "
             "foreground, so the requirement drops by a square.",
             x=0.012, y_title=1.075, y_sub=0.965, override=text)'''
    fig.tight_layout(rect=(0, 0, 1, 0.80))
    plt.show()
    return P._save(fig, "hi_kmodes")


def compute(kedges=None):
    """The per-strategy decomposition. ~1 min: two eigendecompositions."""
    if kedges is None:
        kedges = np.geomspace(8.0, 1.25 * 2 * np.pi / FWHM_RAD, 13)
    res = {"_kedges": kedges}
    for s in ("drift", "raster"):
        path = os.path.join(_HERE, "hi_cache",
                            f"op_{s}_f{FREQ_MHZ:07.3f}_ns{X.NSIDE}.pkl")
        with open(path, "rb") as f:
            mm = pickle.load(f)
        k, n_of_k, tot = n_eff_of_k(mm, kedges)
        area = len(np.asarray(mm.pixel_indices)) * hp.nside2pixarea(
            X.NSIDE, degrees=True)
        res[s] = dict(k=k, n=n_of_k, total=tot, area=area,
                      suff=sufficiency(FLOOR_AMP[s], tot, area))
        print(f"{s:7s} tr(WA)={tot:7.2f}  sum over k={n_of_k.sum():7.2f}  "
              f"area={area:6.1f} deg^2", flush=True)
    return res


def main():
    res = compute()
    strat = [k for k in res if not k.startswith("_")]
    np.savez(os.path.join(_HERE, "results", "hi_kmodes_f350_ns64.npz"),
             kedges=res["_kedges"],
             **{f"{s}_{a}": res[s][a] for s in strat for a in ("k", "n")},
             **{f"{s}_total": res[s]["total"] for s in strat},
             **{f"{s}_area": res[s]["area"] for s in strat},
             d_comoving_mpc=D_COMOVING_MPC)
    print("wrote results/hi_kmodes_f350_ns64.npz", flush=True)
    for s_ in strat:
        for kind, v in res[s_]["suff"].items():
            print(f"  {s_:7s} {kind:5s}: {v['modes']:10.3g} modes "
                  f"({v['factor']:9.3g}x current) -> {v['area']:11.4g} deg^2"
                  + ("" if v["area"] <= 41253 else "   > whole sky"), flush=True)
    print("wrote", fig_kmodes(res), flush=True)
    return res


if __name__ == "__main__":
    main()
