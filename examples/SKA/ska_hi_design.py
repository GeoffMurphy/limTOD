"""PROTOTYPE of the design / trade-space figure (PLAN.md item 1).

**This is a form prototype, not a result.** It exists so the shape of the
primary figure can be argued about while the runs behind it are still moving
(PLAN.md, "Oct, first half"). Items 2, 3, 5, 6 and 7 of the needed-for-draft
table are all still open, and item 6 in particular is the one that would put a
third point on the left panel. Nothing here should be quoted.

Two panels, because the question splits in two:

  left   what does a measured mode per unit sky buy you?  Residual/HI against
         mode density. Only TWO points on it are measured -- the drift and the
         raster at 350 MHz, both at 3 telescope-hours -- so the curve through
         them is a two-point power law and is drawn as such.
  right  which levers move that axis, and by how much. Including the ones that
         do not move it at all, which is the more useful half (PLAN.md,
         "Figure honesty": *put the dead levers on it*).

Honesty constraints this respects, and where it is still short:

* **Modes are counted threshold-free**, as ``tr(WA) = sum lambda/(lambda + S^-1)``
  over each strategy's own patch, then divided by that patch's area. The
  1%-of-lambda_max count is NOT used: it makes the raster look *worse* per unit
  sky than the drift (5.75 vs 6.41 per 100 deg^2) while its residual is 14x
  better, because the raster covers 2.2x the sky. See ``ska_hi_modes.py``.
* **Normalised per unit sky**, so the frequency lever cannot win simply by
  shrinking the patch.
* **Cost is carried** -- but trivially, because the two anchors already cost the
  same 3 telescope-hours (drift 3 x 3600 s, raster 2 x 5400 s). That is stated
  on the panel rather than drawn as an axis. It stops being trivial the moment
  item 5 (depth-matched raster) or a time axis lands.
* **The high-pass cutoff is absent**, deliberately -- it is the one lever whose
  ratio crosses 1, and it must not appear without T(k) beside it.
* **The extrapolation is marked as one.** The slope between the two anchors is
  -6.3, which is steep enough that the extrapolation to ratio = 1 (~37 modes
  per 100 deg^2) is the least defensible number on the chart. It is drawn
  because it is exactly the claim item 6 exists to test, not because it is
  believed.

Known gaps, all of which would change the figure:

* The left panel's y is measured at 350-400 MHz only. The frequency rug along
  its foot shows where frequency lands on the x axis -- 12.2 to 75.1 modes per
  100 deg^2 -- but there is no HI run at those channels, so it carries no y.
  Drawing it as a rug rather than as points is the whole of the honesty here.
* The mode densities for the two anchors come from each strategy's own patch,
  while the residual is measured on the 119 common interior pixels. Those are
  not the same region.
* Single realisation, prior-dependent (item 2 is the flat-prior audit).

Reads ``results/hi_experiment_*``, ``results/hi_ablation_*``,
``results/hi_1f_scan_*``, ``results/hi_modes_*`` and -- computing it on first
run from the experiment-005 operator cache, ~25 s -- ``results/hi_freqmodes.npz``.

    /home/geoff/limTOD/.venv/bin/python ska_hi_design.py

(This one only reads caches, so it runs under either venv. The rest of
experiment 006 needs ``gibbs_venv_312``, which is the one with pyccl.)
"""

from __future__ import annotations

import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

# Palette, rcParams and the _save/_tidy/_title helpers all come from the
# figure module so this prototype cannot drift from the published figures.
import ska_hi_plots as P
from ska_hi_plots import (DRIFT, RASTER, THIRD, SURFACE, INK, INK2, MUTED,
                          GRID, BASELINE, STRAT, RESDIR)

NMODES = 4          # PCA modes removed, the mode count the series leads with
ALPHA_1F = 2.0      # 1/f slope for the whisker family


# ---------------------------------------------------------------------------
# The frequency lever, in the same currency as the anchors
# ---------------------------------------------------------------------------

def freq_mode_density(nside=128, refresh=False):
    """tr(WA) per 100 deg^2 for the five experiment-005 channels.

    Same recipe as ``ska_hi_modes.py`` -- the map-maker's own truth-dependent
    noise weighting and the flat prior the solves use -- applied to the cached
    experiment-005 operators, so the frequency lever is measured in the same
    currency as the two anchors rather than in threshold counts.

    Sanity check on first run: the 1%-of-lambda_max counts this also returns
    reproduce HANDOFF.md's frequency table exactly (6.3, 14.7, 24.6, 34.5,
    46.2 per 100 deg^2).
    """
    path = os.path.join(RESDIR, f"hi_freqmodes_ns{nside}.npz")
    if os.path.exists(path) and not refresh:
        return np.load(path)

    import healpy as hp
    import ska_freq_sweep as S
    from ska_common import gdsm_equatorial_sky_model

    freqs, areas, n_eff, n_thr = [], [], [], []
    for f in S.CHANNELS_MHZ:
        mm = S.build_operator(f, nside)
        pix = np.asarray(mm.pixel_indices)
        truth = gdsm_equatorial_sky_model(freq=f, nside=nside)[pix]
        nv_floor = S.WHITE_VAR * (1e-3 * float(np.mean(truth))) ** 2
        M = None
        for o in mm.Tsys_operators:
            Aop = np.asarray(o, float)
            var = S.WHITE_VAR * (Aop @ truth) ** 2 + nv_floor
            term = Aop.T @ (Aop / var[:, None])
            M = term if M is None else M + term
        w = np.maximum(np.linalg.eigvalsh(M)[::-1], 0.0)
        s_inv = 1.0 / max(float(np.std(truth)), 1e-3) ** 2
        freqs.append(f)
        areas.append(len(pix) * hp.nside2pixarea(nside, degrees=True))
        n_eff.append(float(np.sum(w / (w + s_inv))))
        n_thr.append(int((w / w[0] > 0.01).sum()))
        print(f"  {f:5.0f} MHz  {areas[-1]:6.1f} deg^2  tr(WA) {n_eff[-1]:7.2f}"
              f"  ({100 * n_eff[-1] / areas[-1]:6.2f}/100deg^2)"
              f"  n>1%: {n_thr[-1]:3d}", flush=True)

    out = dict(freq_mhz=np.asarray(freqs, float),
               area_deg2=np.asarray(areas, float),
               n_eff=np.asarray(n_eff, float),
               n_above_1pct=np.asarray(n_thr, float), nside=nside)
    np.savez(path, **out)
    print(f"wrote {path}", flush=True)
    return np.load(path)


# ---------------------------------------------------------------------------
# The numbers the figure is made of, gathered in one place
# ---------------------------------------------------------------------------

def collect(mabl, ladder, modes, fm, nmodes=NMODES, extra_abl=None):
    """Every scalar the design figure draws.

    Primary input is the **matched** run (`hi_ablation_matched_*`), whose three
    strategies -- drift, drift12, raster -- share one HI realisation, one pixel
    set and one set of noise seeds. That matters: ``build_hi_cube`` sizes its box
    from the union of a run's strategies, so numbers are comparable WITHIN a run
    and not across runs (~15-18% realisation scatter; see HANDOFF).

    ``extra_abl`` optionally supplies the **mech** run so drift8 can be shown as
    a third ladder point. It is a *different* realisation, but the two runs'
    shared `drift` reference differs by only 0.9% (409.8 vs 413.5), so the point
    is plotted open and labelled as coming from another run.

    The quantity is the beam + prior FLOOR, which is 96-99% of the total for
    every strategy here and is the term the geometry acts on.

    Everything on this figure is **350-400 MHz**. The 675-725 MHz repeat
    (HANDOFF, experiment 007) inverts the drift/raster ordering, so none of it
    generalises across the band.
    """
    d = {}
    def _floor(abl, s_):
        p_hi = abl[f"{s_}_cr_p_hi"]
        return float(np.median(np.sqrt(abl[f"{s_}_cr_floor_{nmodes}"] / p_hi)))
    for s_ in ("drift", "drift12", "raster"):
        d[s_] = dict(floor=_floor(mabl, s_))
    if extra_abl is not None and "drift8_cr_p_hi" in extra_abl:
        d["drift8"] = dict(floor=_floor(extra_abl, "drift8"), other_run=True)

    lad = {int(n): float(dd) for n, sp, dd
           in zip(ladder["n_strips"], ladder["spacing_beams"], ladder["dens"])
           if abs(float(sp) - 0.5) < 1e-9}
    d["drift"]["dens"] = lad[3]
    if "drift8" in d:
        d["drift8"]["dens"] = lad[8]
    d["drift12"]["dens"] = 21.81       # from hi_ladder_long_ns64.npz
    d["raster"]["dens"] = (100.0 * float(modes["raster_n_eff"])
                           / float(modes["raster_area_deg2"]))
    d["ladder"] = sorted(lad.items())

    # The ladder slope, now from two points in ONE run (drift -> drift12).
    d["slope_ladder"] = (np.log(d["drift12"]["floor"] / d["drift"]["floor"])
                         / np.log(d["drift12"]["dens"] / d["drift"]["dens"]))
    # No extrapolation any more: drift12 sits at the raster's density (21.81 vs
    # 22.31, 2.2% apart), so the cross-linking gain is measured directly.
    d["crosslink_gain"] = d["drift12"]["floor"] / d["raster"]["floor"]
    d["freq"] = dict(mhz=np.asarray(fm["freq_mhz"], float),
                     dens=100.0 * np.asarray(fm["n_eff"], float)
                     / np.asarray(fm["area_deg2"], float))
    return d


# ---------------------------------------------------------------------------
# The design figure
# ---------------------------------------------------------------------------

def fig_design(d, text=None):
    fig, (axA, axB) = plt.subplots(
        1, 2, figsize=P._figsize("hi_design", (14.8, 6.1)),
        gridspec_kw=dict(width_ratios=[1.25, 1.0], wspace=0.52))
    _panel_tradespace(axA, d)
    _panel_levers(axB, d)
    P._title(fig,
             "Cross-linking is not more modes \u2014 it is a different regime",
             "Beam + prior floor at 4 modes removed, common resolution, 119 "
             "interior pixels. Filled points share one HI realisation, one pixel "
             "set and one set of noise seeds; the open point is a different run "
             "(~15% realisation scatter). Mode density is tr($WA$) per unit sky.",
             x=0.012, y_title=1.04, y_sub=0.975, override=text)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    return P._save(fig, "hi_design")


def _panel_tradespace(ax, d):
    P._tidy(ax)
    ax.set_xscale("log"); ax.set_yscale("log")
    xlo, xhi, ylo, yhi = 11.0, 30.0, 0.55, 1200.0
    ax.set_xlim(xlo, xhi); ax.set_ylim(ylo, yhi)

    dr, d12, ra = d["drift"], d["drift12"], d["raster"]

    # --- the ladder family: measured, one run, slope from two of its points --
    xs = np.geomspace(xlo, xhi, 200)
    ys = dr["floor"] * (xs / dr["dens"]) ** d["slope_ladder"]
    seg = (xs >= dr["dens"]) & (xs <= d12["dens"])
    ax.plot(xs, ys, color=DRIFT, lw=1.2, ls=(0, (5, 3)), alpha=0.6, zorder=2)
    ax.plot(xs[seg], ys[seg], color=DRIFT, lw=2.8, zorder=3)

    # drift8 comes from the OTHER run, so it is drawn open and said so.
    if "drift8" in d:
        d8 = d["drift8"]
        ax.plot([d8["dens"]], [d8["floor"]], marker="o", ms=10, mfc=SURFACE,
                mec=DRIFT, mew=2.2, zorder=5)
        ax.annotate("8 dishes\n(other run)", (d8["dens"], d8["floor"]),
                    textcoords="offset points", xytext=(-9, -6), ha="right",
                    va="top", color=MUTED, fontsize=8.2, linespacing=1.4)

    for k, lab, off in (("drift", "3 dishes", (-10, 4)),
                        ("drift12", "12 dishes", (12, 12))):
        ax.plot([d[k]["dens"]], [d[k]["floor"]], marker="o", ms=13,
                color=DRIFT, mec=SURFACE, mew=2.0, zorder=6)
        ax.annotate(lab, (d[k]["dens"], d[k]["floor"]),
                    textcoords="offset points", xytext=off,
                    ha="right" if off[0] < -8 else "left",
                    color=DRIFT, fontsize=9.4, fontweight="semibold")
    ax.annotate(f"parked drift, adding dishes\nslope ${d['slope_ladder']:.2f}$",
                (xlo * 1.04, 110.0), color=DRIFT, fontsize=8.8, ha="left",
                va="center", linespacing=1.45)

    # --- the raster, essentially the same mode density, a decade lower -------
    ax.plot([ra["dens"]], [ra["floor"]], marker="o", ms=13, color=RASTER,
            mec=SURFACE, mew=2.0, zorder=6)
    ax.annotate("raster\ncross-linked", (ra["dens"], ra["floor"]),
                textcoords="offset points", xytext=(10, -2), ha="left",
                va="top", color=RASTER, fontsize=9.4, fontweight="semibold",
                linespacing=1.4)
    # The measured gap. Both endpoints are data; nothing here is extrapolated.
    xbar = np.sqrt(d12["dens"] * ra["dens"])
    ax.plot([xbar, xbar], [ra["floor"], d12["floor"]], color=INK, lw=1.8,
            zorder=4)
    for y in (ra["floor"], d12["floor"]):
        ax.plot([xbar * 0.985, xbar * 1.015], [y, y], color=INK, lw=1.8,
                zorder=4)
    ax.annotate(f"$\\times${d['crosslink_gain']:.0f}\nat the same\nmode density",
                (xbar * 0.97, np.sqrt(ra["floor"] * d12["floor"])), ha="right",
                va="center", color=INK, fontsize=9.6, fontweight="semibold",
                linespacing=1.45)
    ax.annotate(f"{d12['dens']:.1f} vs {ra['dens']:.1f} modes/100 deg$^2$ "
                "\u2014 2% apart",
                (xbar, ylo * 2.2), ha="center", va="center", color=MUTED,
                fontsize=8.2)

    ax.axhline(1.0, color=INK, lw=1.6, zorder=3)
    ax.annotate("HI level", (xlo * 1.04, 1.0), textcoords="offset points",
                xytext=(0, 6), color=INK, fontsize=8.8)

    ax.set_xticks([12, 15, 20, 25, 30])
    ax.set_xticklabels(["12", "15", "20", "25", "30"])
    ax.minorticks_off()
    ax.set_xlabel("measured sky modes per 100 deg$^2$   "
                  "[$\\mathrm{tr}(WA)$ / patch area]")
    ax.set_ylabel(f"beam + prior floor / HI   ({NMODES} modes removed)")
    ax.set_title("Two regimes, measured \u2014 350$-$400 MHz", color=INK2,
                 fontsize=9.6, loc="left", pad=6)


# --- right panel: what each lever does to the FLOOR -------------------------

# (label, factor on the beam+prior floor, kind, note). `None` = not measured.
# The axis is the floor, not the mode count: Tier 2 showed the two are not
# interchangeable in either direction -- the ladder moves modes 1.49x and the
# floor 1.54x, while cross-linking moves modes 1.02x and the floor 10x.
def _levers(d):
    dr, d12, ra = d["drift"], d["drift12"], d["raster"]
    return [
        ("Cross-linking (at matched density)",
         ra["floor"] / d12["floor"], "measured",
         f"drift12 vs raster, {d12['dens']:.1f} vs {ra['dens']:.1f} "
         "modes/100 deg$^2$ \u2014 a step change, not a slope"),
        ("Declination ladder, 3 $\\rightarrow$ 12 dishes",
         d12["floor"] / dr["floor"], "measured",
         f"12 dish-hours against 3; floor slope ${d['slope_ladder']:.2f}$, "
         "linear over 3$-$12"),
        ("Frequency, 350 $\\rightarrow$ 700 MHz (drift)", 1.0 / 105.0,
         "measured",
         "foreground residual only; the HI is 4.3$\\times$ fainter there, so "
         "the RATIO gains 24$\\times$"),
        ("Frequency, 350 $\\rightarrow$ 700 MHz (raster)", 1.0 / 2.5,
         "measured",
         "the raster barely gains \u2014 it was already well conditioned"),
        ("Ladder spacing, 0.5 $\\rightarrow$ 0.25 beams", None, "ambiguous",
         "cost-matched; sign flips with the metric, so not a result"),
        ("Integration time", 1.0, "dead", "adds samples, not rows of $A$"),
        ("System temperature", 1.0, "dead",
         "the floor is 96$-$99% of the residual"),
        ("Co-pointed dishes", 1.0, "dead",
         "same statement as integration time"),
    ]


def _panel_levers(ax, d):
    P._tidy(ax)
    rows = _levers(d)
    y = np.arange(len(rows))[::-1]
    xlo, xhi = 0.005, 3.0
    ax.set_xscale("log")
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(-0.8, len(rows) - 0.3)
    ax.axvline(1.0, color=INK, lw=1.6, zorder=4)
    ax.grid(axis="y", visible=False)

    for yi, (label, factor, kind, note) in zip(y, rows):
        if kind == "measured":
            ax.barh(yi, factor - 1.0, left=1.0, height=0.32, color=THIRD,
                    zorder=3)
            txt = f"$\\times${factor:.2f}" if factor >= 0.1 else \
                  f"$\\times${factor:.3f}"
            ax.annotate(txt, (factor, yi), textcoords="offset points",
                        xytext=(6, 0), ha="left", va="center", color=INK,
                        fontsize=9.4, fontweight="semibold",
                        bbox=dict(facecolor=SURFACE, edgecolor="none", pad=1.0))
        elif kind in ("open", "ambiguous"):
            ax.add_patch(Rectangle((xlo, yi - 0.16), 1.0 - xlo, 0.32,
                                   facecolor="none", edgecolor=BASELINE,
                                   lw=1.1, hatch="////", zorder=3))
            ax.annotate("inconclusive" if kind == "ambiguous" else "not measured",
                        (0.8, yi), va="center", ha="right", color=MUTED,
                        fontsize=8.4, style="italic", zorder=5,
                        bbox=dict(facecolor=SURFACE, edgecolor="none", pad=1.5))
        else:
            ax.plot([1.0], [yi], marker="o", ms=9, color=MUTED, mec=SURFACE,
                    mew=1.6, zorder=5)
            ax.annotate("1.00$\\times$ exactly", (1.12, yi), va="center",
                        color=MUTED, fontsize=8.6)
        ax.annotate(note, (xlo * 1.15, yi - 0.30), color=MUTED, fontsize=7.3,
                    va="center", ha="left")

    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=8.8)
    for tick, (_l, _f, kind, _n) in zip(ax.get_yticklabels(), rows):
        tick.set_color(MUTED if kind == "dead" else INK)
    ax.set_xticks([0.01, 0.05, 0.2, 1.0])
    ax.set_xticklabels(["0.01", "0.05", "0.2", "1$\\times$"])
    ax.minorticks_off()
    ax.set_xlabel("factor change in the beam + prior floor   (left is better)")
    ax.set_title("Frequency is the biggest lever, and it is band-dependent",
                 color=INK2, fontsize=9.6, loc="left", pad=6)


def main():
    matched = os.path.join(RESDIR, "hi_ablation_matched_f350_400_nc32_ns64.npz")
    mech = os.path.join(RESDIR, "hi_ablation_mech_f350_400_nc32_ns64.npz")
    if not os.path.exists(matched):
        sys.exit(f"need {os.path.basename(matched)} (Ilifu job 13857817); "
                 "pull it before drawing the design figure")
    d = collect(np.load(matched),
                np.load(os.path.join(RESDIR, "hi_ladder_f350_ns64.npz")),
                np.load(os.path.join(RESDIR, "hi_modes_f350_400_nc32_ns64.npz")),
                freq_mode_density(),
                extra_abl=np.load(mech) if os.path.exists(mech) else None)
    for s_ in ("drift", "drift8", "drift12", "raster"):
        if s_ in d:
            tag = "  (other run)" if d[s_].get("other_run") else ""
            print(f"{s_:8s} {d[s_]['dens']:6.2f} modes/100deg^2   "
                  f"floor {d[s_]['floor']:7.1f}x{tag}")
    print(f"ladder slope {d['slope_ladder']:.2f} (drift -> drift12); "
          f"cross-linking at matched density is worth "
          f"{1 / d['crosslink_gain']:.1f}x on the floor "
          f"({d['crosslink_gain']:.1f}x lower)")
    print("wrote", fig_design(d))

    lad_long = os.path.join(RESDIR, "hi_ladder_long_ns64.npz")
    print("wrote", fig_ladder(np.load(os.path.join(
        RESDIR, "hi_ladder_f350_ns64.npz"))))
    print("wrote", fig_ladder_geometry(np.load(os.path.join(
        RESDIR, "hi_ladder_f350_ns64.npz"))))
    pv = os.path.join(RESDIR, "hi_priorvar_tk_f350_400_nc32_ns64.npz")
    if os.path.exists(pv):
        print("wrote", fig_priorvar(np.load(pv)))


def fig_priorvar(pv, d=None, nmodes=NMODES, text=None):
    """Does the trade space survive the prior, or is it a picture of the prior?

    Left: the residual against the prior variance, floor and total, both
    strategies. Right: the same sweep drawn as a *track* across the trade
    space, against the two-point power law the design figure draws through the
    published anchors. The question the right panel asks is the only one that
    matters for the design figure:

      * a track that runs ALONG the power law means the curve is a statement
        about the prior. Loosening the prior would then buy you residual for
        free, and the chart is recommending a prior, not a survey.
      * a track that runs ACROSS it means the mode count and the residual are
        separately real, and the curve is a statement about the survey.

    ``d`` is the design-figure dict from ``collect``; passing it draws the
    published anchors and the power law for reference.
    """
    scales = np.asarray(pv["scales"], float)
    fig, (axA, axB) = plt.subplots(
        1, 2, figsize=P._figsize("hi_priorvar", (13.4, 5.6)),
        gridspec_kw=dict(width_ratios=[1.0, 1.0], wspace=0.26))

    # --- left: residual against the prior scale ----------------------------
    P._tidy(axA)
    for s, style in STRAT.items():
        axA.plot(scales, pv[f"{s}_total_{nmodes}"], color=style["color"],
                 marker="o", ms=7, mfc=SURFACE, mew=2.0, zorder=4)
        axA.plot(scales, pv[f"{s}_floor_{nmodes}"], color=style["color"],
                 ls=(0, (4, 2)), lw=1.6, alpha=0.75, zorder=3)
        axA.annotate(s, (scales[-1], pv[f"{s}_total_{nmodes}"][-1]),
                     textcoords="offset points", xytext=(9, 0), va="center",
                     color=style["color"], fontsize=9.5, fontweight="semibold")
    axA.axhline(1.0, color=INK, lw=1.4, zorder=2)
    axA.annotate("HI level", (scales[0], 1.0), textcoords="offset points",
                 xytext=(2, 6), color=INK, fontsize=8.8)
    axA.axvline(1.0, color=MUTED, lw=1.0, ls=(0, (1, 2)), zorder=1)
    axA.annotate("as published", (1.0, axA.get_ylim()[1]),
                 textcoords="offset points", xytext=(4, -10), color=MUTED,
                 fontsize=8.4, rotation=90, va="top")
    axA.set_xscale("log"); axA.set_yscale("log")
    axA.set_xlim(scales.min() * 0.5, scales.max() * 3.0)
    axA.set_xlabel("prior variance / truth patch variance   "
                   "(right = looser prior)")
    axA.set_ylabel(f"post-clean residual / HI   ({nmodes} modes removed)")
    axA.legend(handles=[
        Line2D([], [], color=MUTED, marker="o", ms=7, mfc=SURFACE, mew=2.0,
               label="total"),
        Line2D([], [], color=MUTED, ls=(0, (4, 2)), lw=1.6,
               label="beam + prior floor alone")],
        loc="lower right", bbox_to_anchor=(1.0, 0.07), fontsize=8.8)
    axA.set_title("What the prior's amplitude is worth", color=INK2,
                  fontsize=9.6, loc="left", pad=6)

    # --- right: the same sweep as a track across the trade space -----------
    P._tidy(axB)
    axB.set_xscale("log"); axB.set_yscale("log")
    bx0, bx1, by0, by1 = 7.0, 95.0, 3.0, 6e3
    axB.set_xlim(bx0, bx1); axB.set_ylim(by0, by1)
    if d is not None:
        # Clipped to the panel: unclipped it spans five decades and flattens
        # everything else. Both published anchors lie ON it by construction --
        # it is the line through them -- so what the eye is being asked to
        # judge is the ANGLE the prior tracks cross it at.
        xs = np.geomspace(bx0, bx1, 400)
        ys = d["raster"]["ratio"] * (xs / d["raster"]["dens"]) ** d["slope"]
        keep = (ys > by0) & (ys < by1)
        axB.plot(xs[keep], ys[keep], color=MUTED, lw=1.6, ls=(0, (5, 3)),
                 zorder=2)
        axB.annotate("two-point power law from the\ndesign figure "
                     f"(slope ${d['slope']:.1f}$)",
                     (bx0 * 1.05, 2600), color=MUTED, fontsize=8.4, ha="left",
                     va="center", linespacing=1.45)
    lab_off = {"drift": (0, -30, "center", "top"),
               "raster": (12, -12, "left", "top")}
    for s_, style in STRAT.items():
        x = np.asarray(pv[f"{s_}_dens"], float)
        y = np.asarray(pv[f"{s_}_total_{nmodes}"], float)
        axB.plot(x, y, color=style["color"], lw=2.2, zorder=4)
        axB.plot(x, y, color=style["color"], ls="none", marker="o", ms=6,
                 mfc=SURFACE, mew=1.8, zorder=5)
        j = int(np.argmin(np.abs(scales - 1.0)))
        axB.plot([x[j]], [y[j]], color=style["color"], marker="o", ms=13,
                 mec=SURFACE, mew=2.0, zorder=6)
        dx, dy, ha, va = lab_off[s_]
        axB.annotate(f"{s_}\nas published", (x[j], y[j]),
                     textcoords="offset points", xytext=(dx, dy), ha=ha, va=va,
                     color=style["color"], fontsize=9.2,
                     fontweight="semibold", linespacing=1.4)
        axB.annotate("tight", (x[0], y[0]), textcoords="offset points",
                     xytext=(-7, 0), ha="right", va="center", color=MUTED,
                     fontsize=8.0)
        axB.annotate("loose", (x[-1], y[-1]), textcoords="offset points",
                     xytext=(7, 0), ha="left", va="center", color=MUTED,
                     fontsize=8.0)
    axB.set_xticks([10, 15, 20, 30, 50, 70])
    axB.set_xticklabels(["10", "15", "20", "30", "50", "70"])
    axB.minorticks_off()
    axB.set_xlabel("measured sky modes per 100 deg$^2$   "
                   "[$\\mathrm{tr}(WA)$ / patch area]")
    axB.set_ylabel(f"post-clean residual / HI   ({nmodes} modes removed)")
    axB.set_title("The tracks cross the line, they do not run along it",
                  color=INK2, fontsize=9.6, loc="left", pad=6)

    P._title(fig, "Loosening the prior buys modes that are worth less than nothing",
             f"Prior variance swept over {scales.min():g}$-${scales.max():g}$\\times$ "
             "the truth patch variance, everything else fixed (the prior MEAN was "
             "already flat; this is its amplitude). Scale 1 reproduces the published "
             "335.0/347.3 and 22.7/24.0 exactly. PROTOTYPE.",
             x=0.012, y_title=1.05, y_sub=0.982, override=text)
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    return P._save(fig, "hi_priorvar")

# ---------------------------------------------------------------------------
# The drift's Dec ladder (ska_hi_ladder.py)
# ---------------------------------------------------------------------------

def fig_ladder(ld, ld_long=None, floors=None, text=None):
    """What the drift's only geometry lever is worth, out to 20 dishes.

    Three panels:

      left    density against strip SPACING. Below one beam the strips overlap
              and densify the cross-scan direction a single drift strip cannot
              sample; at or above one beam they are disjoint and you buy area.
      middle  density against strip COUNT, now to N = 20, with N = 1 (a single
              parked dish) as the baseline.
      right   the decomposition, normalised to N = 1. The marginal return per
              strip is CONSTANT -- 11.7-12.7 modes for 43-47 deg² every time --
              so the density curve flattens purely because area grows linearly
              alongside modes, not because the ladder saturates.

    ``floors`` optionally overlays the two measured HI floors (3 and 12 dishes)
    on the middle panel, which is the only place mode count and residual appear
    together for the same geometry family.
    """
    n = np.asarray(ld["n_strips"], int)
    sp = np.asarray(ld["spacing_beams"], float)
    dens = np.asarray(ld["dens"], float)
    n_eff = np.asarray(ld["n_eff"], float)
    area = np.asarray(ld["area"], float)
    ref_n, ref_sp = int(ld["ref_count"]), float(ld["ref_spacing_beams"])

    # Merge the long sweep, so the count arm runs 1..20 rather than 1..8.
    if ld_long is not None:
        m = np.abs(np.asarray(ld_long["spacing_beams"], float) - ref_sp) < 1e-9
        for k, v, a, d in zip(np.asarray(ld_long["n_strips"], int)[m],
                              np.asarray(ld_long["n_eff"], float)[m],
                              np.asarray(ld_long["area"], float)[m],
                              np.asarray(ld_long["dens"], float)[m]):
            if k not in set(n[np.abs(sp - ref_sp) < 1e-9]):
                n = np.r_[n, k]; sp = np.r_[sp, ref_sp]
                n_eff = np.r_[n_eff, v]; area = np.r_[area, a]; dens = np.r_[dens, d]

    A = np.argsort(sp[np.abs(sp - sp) < 1])  # placeholder, replaced below
    iA = np.where(n == ref_n)[0][np.argsort(sp[n == ref_n])]
    iB = np.where(np.abs(sp - ref_sp) < 1e-9)[0]
    iB = iB[np.argsort(n[iB])]

    fig, axes = plt.subplots(1, 3, figsize=P._figsize("hi_ladder", (14.6, 4.9)),
                             gridspec_kw=dict(wspace=0.30))

    ax = axes[0]; P._tidy(ax)
    ax.plot(sp[iA], dens[iA], color=DRIFT, marker="o", ms=8, mfc=SURFACE,
            mew=2.0, zorder=4)
    j0 = int(np.where((n[iA] == ref_n) & (np.abs(sp[iA] - ref_sp) < 1e-9))[0][0])
    ax.plot([sp[iA][j0]], [dens[iA][j0]], color=DRIFT, marker="o", ms=13,
            mec=SURFACE, mew=2.0, zorder=5)
    ax.annotate("published\n52/50/48", (sp[iA][j0], dens[iA][j0]),
                textcoords="offset points", xytext=(0, 16), ha="center",
                color=DRIFT, fontsize=9.0, fontweight="semibold",
                linespacing=1.4)
    ax.axvline(1.0, color=INK, lw=1.4, ls=(0, (4, 3)), zorder=2)
    ax.annotate("one beam \u2014 strips stop overlapping", (1.0, dens[iA].min()),
                textcoords="offset points", xytext=(-6, 4), rotation=90,
                ha="right", va="bottom", color=INK, fontsize=8.4)
    ax.set_xlabel("strip spacing  [beam FWHM]")
    ax.set_ylabel("measured sky modes per 100 deg$^2$")
    ax.set_title(f"{ref_n} strips, varying spacing", color=INK2, fontsize=9.6,
                 loc="left", pad=6)

    ax = axes[1]; P._tidy(ax)
    ax.plot(n[iB], dens[iB], color=DRIFT, marker="o", ms=7, mfc=SURFACE,
            mew=2.0, zorder=4)
    j1 = int(np.where(n[iB] == 1)[0][0])
    ax.plot([1], [dens[iB][j1]], color=MUTED, marker="o", ms=11, mec=SURFACE,
            mew=1.8, zorder=5)
    ax.annotate("one parked dish\n(no ladder)", (1, dens[iB][j1]),
                textcoords="offset points", xytext=(10, -2), ha="left",
                va="top", color=MUTED, fontsize=8.6, linespacing=1.4)
    if floors:
        for nn, fl in floors.items():
            k = int(np.where(n[iB] == nn)[0][0])
            ax.plot([nn], [dens[iB][k]], color=DRIFT, marker="o", ms=13,
                    mec=SURFACE, mew=2.0, zorder=6)
            ax.annotate(f"{nn} dishes\nfloor {fl:.0f}$\\times$",
                        (nn, dens[iB][k]), textcoords="offset points",
                        xytext=(0, 15), ha="center", color=DRIFT,
                        fontsize=8.8, fontweight="semibold", linespacing=1.4)
    ax.set_xlabel("number of strips  (= dishes, observing simultaneously)")
    ax.set_ylabel("measured sky modes per 100 deg$^2$")
    ax.set_title(f"spacing fixed at {ref_sp:g} beams", color=INK2,
                 fontsize=9.6, loc="left", pad=6)

    ax = axes[2]; P._tidy(ax)
    ax.plot(n[iB], n_eff[iB] / n_eff[iB][j1], color=DRIFT, marker="o", ms=7,
            mfc=SURFACE, mew=2.0, zorder=4, label="measured modes, tr($WA$)")
    ax.plot(n[iB], area[iB] / area[iB][j1], color=THIRD, marker="s", ms=6,
            mfc=SURFACE, mew=2.0, zorder=4, label="sky area")
    ax.axhline(1.0, color=BASELINE, lw=1.0, ls=(0, (4, 3)), zorder=2)
    ax.set_xlabel("number of strips")
    ax.set_ylabel("relative to a single parked dish")
    ax.legend(loc="upper left", fontsize=8.8)
    ax.set_title("Constant marginal return: ~12 modes per 44 deg$^2$, every strip",
                 color=INK2, fontsize=9.6, loc="left", pad=6)

    P._title(fig, "The drift cannot cross-link, so the Dec ladder is its only "
             "geometry lever",
             "Mode count only, one channel at 350.8 MHz, no TOD. N strips on N "
             "sidereal days is arithmetically the same as N dishes parked at N "
             "elevations for one hour. The ladder does not saturate \u2014 but a "
             "mode count is not a residual, and cross-linking beats it 9$\\times$ "
             "at matched density.",
             x=0.012, y_title=1.05, y_sub=0.975, override=text)
    fig.tight_layout(rect=(0, 0, 1, 0.89))
    return P._save(fig, "hi_ladder")


def fig_ladder_geometry(ld=None, text=None):
    """The ladder as geometry, for explaining the mode-count result to others.

    Parked due north, a dish at elevation ``el`` points at declination
    ``dec = 90 + phi - el`` (phi = -30.713 at the Karoo), so elevation steps are
    declination steps one-for-one, and each strip is a band of constant Dec
    swept through 15 deg of RA by one hour of Earth rotation.

    Top row: the footprints. Bands are drawn one FWHM thick and semi-transparent,
    so overlap darkens. Bottom row: the same thing as a cross-scan cut -- the
    beam response against Dec, individual strips light, their sum bold. The sum
    is the honest picture of what the ladder does: tight spacing gives a narrow,
    smoothly covered band; wide spacing gives separated ribbons with troughs
    between them, which is the regime where the strips have become independent
    surveys and the ladder has stopped buying information.

    The bands show FWHM only. The map-maker keeps every pixel above 5 per cent
    of beam peak, which reaches considerably further, so the areas quoted from
    ``ska_hi_ladder.py`` are larger than these rectangles -- they are annotated
    from the measurement, never derived from the drawing.
    """
    from matplotlib.patches import Rectangle

    lat, fwhm = -30.713, 3.983
    ra_span = 15.0                       # one hour of drift
    sigma = fwhm / 2.3548

    # (n_strips, spacing_beams) -- the spacing arm, plus one point of the
    # count arm so both levers appear in the same picture.
    panels = [(3, 0.25), (3, 0.5), (3, 1.0), (3, 1.5), (8, 0.5)]
    meas = {}
    if ld is not None:
        for n, sp, d, a in zip(ld["n_strips"], ld["spacing_beams"],
                               ld["dens"], ld["area"]):
            meas[(int(n), round(float(sp), 2))] = (float(d), float(a))

    fig, axes = plt.subplots(2, len(panels), figsize=P._figsize(
        "hi_ladder_geometry", (15.0, 6.4)),
        gridspec_kw=dict(height_ratios=[1.25, 1.0], hspace=0.38, wspace=0.22))

    decs_all = []
    for n, sp in panels:
        els = EL_CENTRE_ = 50.0 + (np.arange(n) - (n - 1) / 2.0) * sp * fwhm
        decs_all.append(90.0 + lat - els)
    dlo = min(d.min() for d in decs_all) - fwhm
    dhi = max(d.max() for d in decs_all) + fwhm

    grid = np.linspace(dlo, dhi, 900)
    ymax = max(np.sum([np.exp(-0.5 * ((grid - d) / sigma) ** 2) for d in decs],
                      axis=0).max() for decs in decs_all)

    for col, ((n, sp), decs) in enumerate(zip(panels, decs_all)):
        # --- top: footprints on the sky ---------------------------------
        ax = axes[0, col]; P._tidy(ax)
        for d in decs:
            ax.add_patch(Rectangle((0, d - fwhm / 2), ra_span, fwhm,
                                   facecolor=DRIFT, edgecolor="none",
                                   alpha=0.42, zorder=3))
            ax.plot([0, ra_span], [d, d], color=SURFACE, lw=0.8, zorder=4)
        ax.set_xlim(-1.0, ra_span + 1.0); ax.set_ylim(dlo, dhi)
        ax.set_xlabel("RA drift  [deg]", fontsize=8.6)
        if col == 0:
            ax.set_ylabel("declination  [deg]")
        else:
            ax.set_yticklabels([])
        lab = f"{n} strips, {sp:g} beam" + ("" if sp == 1 else "s")
        if (n, sp) == (3, 0.5):
            lab += "\n(published 52/50/48)"
        elif (n, sp) == (8, 0.5):
            lab += "\n(8 dishes)"
        ax.set_title(lab, color=INK, fontsize=9.4, fontweight="semibold",
                     loc="left", pad=6, linespacing=1.35)

        # --- bottom: the cross-scan cut ---------------------------------
        ax = axes[1, col]; P._tidy(ax)
        tot = np.zeros_like(grid)
        for d in decs:
            g = np.exp(-0.5 * ((grid - d) / sigma) ** 2)
            tot += g
            ax.plot(grid, g, color=DRIFT, lw=1.0, alpha=0.40, zorder=3)
        ax.plot(grid, tot, color=DRIFT, lw=2.4, zorder=5)
        ax.set_xlim(dlo, dhi); ax.set_ylim(0, ymax * 1.10)
        ax.set_xlabel("declination  [deg]", fontsize=8.6)
        if col == 0:
            ax.set_ylabel("summed beam response")
        else:
            ax.set_yticklabels([])
        # Only report a trough where one exists. Below ~0.5 beams the sum is
        # a single peak and "trough depth" measured between the outer strip
        # centres is meaningless -- it reads as if tighter packing were worse.
        lo = np.r_[False, (tot[1:-1] < tot[:-2]) & (tot[1:-1] < tot[2:]), False]
        inner = lo & (grid > decs.min()) & (grid < decs.max())
        if n > 1 and inner.any():
            note = f"dips to {tot[inner].min() / tot.max():.2f} of peak"
        elif n > 1:
            note = "no dip — strips fully merged"
        else:
            note = "single strip"
        ax.annotate(note, (0.5, 0.94), xycoords="axes fraction", ha="center",
                    color=MUTED, fontsize=8.0)
        if (n, round(sp, 2)) in meas:
            d_, a_ = meas[(n, round(sp, 2))]
            ax.annotate(f"{d_:.1f} modes/100 deg$^2$\n{a_:.0f} deg$^2$ measured",
                        (0.5, 0.74), xycoords="axes fraction", ha="center",
                        color=INK, fontsize=8.6, fontweight="semibold",
                        linespacing=1.4)

    P._title(fig, "The declination ladder: overlapping strips concentrate "
             "information, separated strips just add sky",
             "Parked due north, $\\delta = 90\\degree + \\phi - \\mathrm{el}$, so "
             "elevation steps are declination steps one-for-one, and one hour of "
             "Earth rotation sweeps each strip through 15$\\degree$ of RA. Bands "
             "are one FWHM thick; the map-maker keeps everything above 5% of beam "
             "peak, so measured areas (annotated) exceed the rectangles.",
             x=0.012, y_title=1.035, y_sub=0.985, override=text)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    return P._save(fig, "hi_ladder_geometry")


# ---------------------------------------------------------------------------
# Why residual/HI must not be compared across bands
# ---------------------------------------------------------------------------

BANDS = [
    ("350$-$400", "hi_experiment_matched_f350_400_nc32_ns64.npz", 375.0),
    ("500$-$550", "hi_experiment_f500_f500_550_nc32_ns128.npz", 525.0),
    ("675$-$725", "hi_experiment_f700_f675_725_nc32_ns128.npz", 700.0),
]


def load_bands(nmodes=NMODES, strategies=("drift", "raster")):
    """Foreground residual, HI amplitude and their ratio, per band.

    Only bands whose npz is present are returned, so this works with two bands
    and improves with three.
    """
    out = []
    for label, fname, fc in BANDS:
        path = os.path.join(RESDIR, fname)
        if not os.path.exists(path):
            continue
        e = np.load(path)
        rec = dict(label=label, fc=fc)
        rec["hi"] = float(np.median(np.sqrt(e["drift_cr_p_true"])))
        for s_ in strategies:
            if f"{s_}_cr_p_true" not in e:
                continue
            rec[s_] = dict(
                fg=float(np.median(np.sqrt(e[f"{s_}_cr_p_clean_auto_{nmodes}"]))),
                ratio=float(np.median(np.sqrt(
                    e[f"{s_}_cr_p_clean_auto_{nmodes}"]
                    / e[f"{s_}_cr_p_mapmade"]))))
        out.append(rec)
    return out


def fig_frequency(bands, text=None):
    """Three panels: the numerator, the denominator, and why their ratio lies.

    ``residual/HI`` is the quantity this series quotes throughout, and within a
    band it is the right one. Across bands it is not: the HI's own brightness
    falls with redshift, so the denominator moves for reasons that have nothing
    to do with the instrument. Plotting the two factors separately is the only
    honest way to show a frequency dependence.
    """
    if len(bands) < 2:
        return None
    x = [b["fc"] for b in bands]
    fig, axes = plt.subplots(1, 3, figsize=P._figsize("hi_frequency_design",
                                                      (14.2, 4.6)),
                             gridspec_kw=dict(wspace=0.32))

    panels = [
        ("fg", "post-clean foreground residual", "what the INSTRUMENT does"),
        (None, "HI amplitude  $\\sqrt{P_{\\mathrm{true}}}$",
         "what COSMOLOGY does"),
        ("ratio", f"residual / HI   ({NMODES} modes removed)",
         "their quotient — and it inverts"),
    ]
    for ax, (key, ylab, title) in zip(axes, panels):
        P._tidy(ax); ax.set_yscale("log")
        if key is None:
            y = [b["hi"] for b in bands]
            ax.plot(x, y, color=INK, marker="o", ms=8, mfc=SURFACE, mew=2.0,
                    zorder=4)
            ax.annotate(f"{y[0] / y[-1]:.1f}$\\times$ fainter\nacross this range",
                        (x[0], y[-1]), textcoords="offset points",
                        xytext=(8, 0), ha="left", va="center", color=INK,
                        fontsize=8.8, linespacing=1.45)
        else:
            for s_, style in STRAT.items():
                if s_ not in bands[0]:
                    continue
                y = [b[s_][key] for b in bands]
                ax.plot(x, y, color=style["color"], marker="o", ms=8,
                        mfc=SURFACE, mew=2.0, zorder=4)
                ax.annotate(s_, (x[-1], y[-1]), textcoords="offset points",
                            xytext=(8, 0), va="center", color=style["color"],
                            fontsize=9.4, fontweight="semibold")
                fac = y[0] / y[-1]
                ax.annotate(f"{fac:.0f}$\\times$" if fac >= 10
                            else f"{fac:.1f}$\\times$",
                            (x[-1], y[-1]), textcoords="offset points",
                            xytext=(8, 13), va="center", color=MUTED,
                            fontsize=8.4)
            if key == "ratio":
                ax.axhline(1.0, color=INK, lw=1.4, zorder=3)
                ax.annotate("HI level", (x[0], 1.0),
                            textcoords="offset points", xytext=(2, 6),
                            color=INK, fontsize=8.4)
        ax.set_xlabel("band centre  [MHz]")
        ax.set_ylabel(ylab, fontsize=9.2)
        ax.set_title(title, color=INK2, fontsize=9.6, loc="left", pad=6)
        ax.set_xticks(x); ax.set_xticklabels([b["label"] for b in bands],
                                             fontsize=8.4)
        ax.set_xlim(min(x) - 60, max(x) + 110)
        ax.margins(y=0.18)

    P._title(fig, "residual / HI is not comparable across bands",
             "Left: the drift's foreground residual improves far faster with "
             "frequency than the raster's. Middle: the HI itself dims with "
             "redshift, which is cosmology, not instrument. Right: dividing one "
             "by the other inverts the ordering — so quote the two factors "
             "separately whenever bands are compared.",
             x=0.012, y_title=1.06, y_sub=0.975, override=text)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    return P._save(fig, "hi_frequency_design")


# ---------------------------------------------------------------------------
# Would cleaning harder help?
# ---------------------------------------------------------------------------

def fig_cleaning_depth(text=None):
    """The first question asked of any floor-limited result: clean deeper?

    Removing more PCA modes lowers the residual and destroys more HI, and the
    two are not separable -- the floor is spectrally non-smooth precisely
    because the null space is chromatic, so it overlaps the signal. The useful
    quantity is therefore not the residual but the residual measured against
    the HI that SURVIVES the clean, ``ratio / T``. If that has a minimum, there
    is an optimal cleaning depth; if it falls monotonically, the limit is not
    the cleaning.

    Both bands, so the answer can be seen to be band-dependent or not.
    """
    runs = [("350$-$400", "hi_experiment_matched_f350_400_nc32_ns64.npz",
             ("drift", "drift12", "raster"), "-"),
            ("675$-$725", "hi_experiment_f700_f675_725_nc32_ns128.npz",
             ("drift", "raster"), (0, (4, 2)))]
    colour = {"drift": DRIFT, "drift12": THIRD, "raster": RASTER}

    fig, (axA, axB) = plt.subplots(
        1, 2, figsize=P._figsize("hi_cleaning_depth", (12.6, 5.0)),
        gridspec_kw=dict(wspace=0.28))
    seen, strat_seen = [], []
    for band, fname, strategies, ls in runs:
        path = os.path.join(RESDIR, fname)
        if not os.path.exists(path):
            continue
        e = np.load(path)
        nm = [int(v) for v in e["nmodes_grid"]]
        for s_ in strategies:
            if f"{s_}_cr_tf_{nm[0]}" not in e:
                continue
            ratio = [float(np.median(np.sqrt(
                e[f"{s_}_cr_p_clean_auto_{k}"] / e[f"{s_}_cr_p_mapmade"])))
                for k in nm]
            tf = [float(np.median(e[f"{s_}_cr_tf_{k}"])) for k in nm]
            net = [r / t for r, t in zip(ratio, tf)]
            c = colour[s_]
            axA.plot(nm, ratio, color=c, ls=ls, marker="o", ms=6, mfc=SURFACE,
                     mew=1.8, zorder=4)
            axB.plot(nm, net, color=c, ls=ls, marker="o", ms=6, mfc=SURFACE,
                     mew=1.8, zorder=4)
            seen.append((band, ls))
            if s_ not in [t[0] for t in strat_seen]:
                strat_seen.append((s_, c))
    for ax, ylab, title in (
            (axA, f"post-clean residual / HI",
             "What deeper cleaning buys"),
            (axB, "residual / HI that SURVIVES the clean",
             "What it costs — and whether it is worth it")):
        P._tidy(ax); ax.set_yscale("log")
        ax.axhline(1.0, color=INK, lw=1.5, zorder=3)
        ax.annotate("HI level", (12.0, 1.0), textcoords="offset points",
                    xytext=(0, 6), ha="right", color=INK, fontsize=8.6)
        ax.set_xlabel("PCA modes removed")
        ax.set_ylabel(ylab, fontsize=9.4)
        ax.set_title(title, color=INK2, fontsize=9.6, loc="left", pad=6)
        ax.set_xticks([1, 2, 3, 4, 6, 8, 10])
        ax.set_xlim(0.5, 12.2)
        ax.set_ylim(bottom=0.7)
    bands = list(dict.fromkeys(seen))
    leg1 = axA.legend(handles=[Line2D([], [], color=c, lw=2.4, label=s_)
                               for s_, c in strat_seen],
                      loc="upper right", bbox_to_anchor=(1.0, 1.0),
                      fontsize=8.8, title="strategy", title_fontsize=8.4,
                      labelspacing=0.35)
    axA.add_artist(leg1)
    axA.legend(handles=[Line2D([], [], color=MUTED, ls=ls, lw=2.0,
                               label=f"{b} MHz") for b, ls in bands],
               loc="upper right", bbox_to_anchor=(1.0, 0.62),
               fontsize=8.8, title="band", title_fontsize=8.4,
               labelspacing=0.35)

    P._title(fig, "Cleaning harder always helps — which is the point",
             "The residual keeps falling with mode count and so does the HI, "
             "but the residual falls faster, so the net never turns over within "
             "the range tested. There is no optimal depth to find: the limit is "
             "the beam + prior floor, not the cleaning. Common resolution, "
             "interior pixels.",
             x=0.012, y_title=1.05, y_sub=0.975, override=text)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    return P._save(fig, "hi_cleaning_depth")


# ---------------------------------------------------------------------------
# Why the floor is immune to foreground cleaning
# ---------------------------------------------------------------------------

def fig_chromatic(rank, text=None):
    """The physical heart of the result, in three panels.

    Foreground cleaning works by assuming the contaminant is spectrally smooth
    and therefore low-rank in frequency. That is true of the foreground itself,
    and it is why one PCA mode removes it. It is NOT true of the beam + prior
    floor: the floor is ``(I - WA)s``, and ``WA`` is chromatic because the beam
    (and hence the measured subspace) changes across the band. The floor
    therefore inherits spectral structure the sky never had.

      left    rank. The foreground is rank 1. The floor needs 2-6 modes. The HI
              needs 23-25 of 32 -- it is the one component that is genuinely
              not compressible.
      middle  the eigenvectors themselves, which is what "chromatic" means
              concretely: the foreground's first mode is a clean power law, the
              floor's get progressively more structured, and the HI's are
              noise-like from the start.
      right   principal angles between the floor's subspace and the HI's. They
              are not orthogonal, so projecting out the floor necessarily takes
              HI with it. This is why deeper cleaning cannot win.
    """
    f = np.asarray(rank["freqs"], float)
    arms = [("fg", "foreground", INK), ("floor", "beam + prior floor", DRIFT),
            ("hi", "HI", RASTER)]
    fig, axes = plt.subplots(1, 3, figsize=P._figsize("hi_chromatic", (14.4, 4.7)),
                             gridspec_kw=dict(wspace=0.30))

    # --- left: how many modes each component occupies ----------------------
    ax = axes[0]; P._tidy(ax)
    for s_, ls in (("drift", "-"), ("raster", (0, (4, 2)))):
        for arm, lab, c in arms:
            w = np.asarray(rank[f"{s_}_w_{arm}"], float)
            cum = np.cumsum(w) / w.sum()
            ax.plot(np.arange(1, len(cum) + 1), cum, color=c, ls=ls, lw=2.0,
                    zorder=4, label=lab if s_ == "drift" else None)
    ax.axhline(0.99, color=BASELINE, lw=1.0, ls=(0, (4, 3)), zorder=2)
    ax.annotate("99% of variance", (len(f), 0.99), textcoords="offset points",
                xytext=(-2, -12), ha="right", color=MUTED, fontsize=8.4)
    ax.set_xlim(1, 30); ax.set_ylim(0.3, 1.02)
    ax.set_xlabel("eigenmodes retained")
    ax.set_ylabel("cumulative fraction of variance")
    leg = ax.legend(loc="lower right", fontsize=8.8)
    ax.add_artist(leg)
    ax.legend(handles=[Line2D([], [], color=MUTED, ls="-", lw=2.0,
                              label="drift"),
                       Line2D([], [], color=MUTED, ls=(0, (4, 2)), lw=2.0,
                              label="raster")],
              loc="center right", fontsize=8.4, labelspacing=0.3)
    ax.set_title("The HI is the only thing that is not low-rank",
                 color=INK2, fontsize=9.6, loc="left", pad=6)

    # --- middle: the spectral shapes ---------------------------------------
    ax = axes[1]; P._tidy(ax)
    V = {a: np.asarray(rank[f"drift_V_{a}"], float) for a, _, _ in arms}
    def norm(v):
        v = v - v.mean()
        return v / (np.abs(v).max() or 1.0)
    ax.plot(f, norm(V["fg"][:, 0]), color=INK, lw=2.4, zorder=5,
            label="foreground, mode 1")
    for k, alpha in zip(range(3), (1.0, 0.65, 0.42)):
        ax.plot(f, norm(V["floor"][:, k]) + 2.4 * (k + 1), color=DRIFT, lw=2.0,
                alpha=alpha, zorder=4,
                label="floor, modes 1$-$3" if k == 0 else None)
    ax.plot(f, norm(V["hi"][:, 0]) + 2.4 * 4, color=RASTER, lw=1.6, zorder=4,
            label="HI, mode 1")
    ax.set_yticks([]); ax.set_xlabel("frequency  [MHz]")
    ax.set_ylabel("eigenvector  (offset, normalised)")
    ax.set_xlim(f[0] - 1, f[-1] + 17)
    for y0, lab, c in ((0.0, "foreground\nmode 1", INK),
                       (2.4, "floor\nmodes 1$-$3", DRIFT),
                       (2.4 * 4, "HI\nmode 1", RASTER)):
        ax.annotate(lab, (f[-1] + 1.5, y0), ha="left", va="center", color=c,
                    fontsize=8.4, fontweight="semibold", linespacing=1.35)
    ax.set_title("Smooth, then structured, then noise-like",
                 color=INK2, fontsize=9.6, loc="left", pad=6)

    # --- right: the floor and the HI are not orthogonal ---------------------
    ax = axes[2]; P._tidy(ax)
    for s_, style in STRAT.items():
        c = np.asarray(rank[f"{s_}_overlap"], float)
        ax.plot(np.arange(1, len(c) + 1), c, color=style["color"], marker="o",
                ms=6, mfc=SURFACE, mew=1.8, zorder=4, label=s_)
    ax.axhline(0.0, color=INK, lw=1.4, zorder=3)
    ax.annotate("0 would mean the floor could be removed\nwithout touching the HI",
                (1, 0.06), color=MUTED, fontsize=8.4, ha="left", va="bottom",
                linespacing=1.45)
    ax.set_ylim(-0.05, 1.0)
    ax.set_xlabel("principal angle index")
    ax.set_ylabel("cos(principal angle), floor vs HI")
    ax.legend(loc="upper left", fontsize=8.8)
    ax.set_title("They overlap, so cleaning one costs the other",
                 color=INK2, fontsize=9.6, loc="left", pad=6)

    P._title(fig, "The floor is low-rank but not spectrally smooth — which "
             "is why PCA cannot remove it",
             "Foreground cleaning assumes the contaminant is smooth in "
             "frequency. The foreground is (rank 1). The beam + prior floor is "
             "not: $WA$ is chromatic, so the floor carries structure the sky "
             "never had, and its subspace overlaps the HI's at "
             "cos $=0.42-0.78$. 350$-$400 MHz.",
             x=0.012, y_title=1.05, y_sub=0.975, override=text)
    fig.tight_layout(rect=(0, 0, 1, 0.89))
    return P._save(fig, "hi_chromatic")


# ---------------------------------------------------------------------------
# The error bar nobody had measured
# ---------------------------------------------------------------------------

SEED_RUNS = ["hi_experiment_f350_400_nc32_ns64.npz",
             "hi_experiment_mech_f350_400_nc32_ns64.npz",
             "hi_experiment_matched_f350_400_nc32_ns64.npz"] + \
            [f"hi_experiment_seed{sd}_f350_400_nc32_ns64.npz"
             for sd in (11, 12, 13, 14, 15, 16)]


def load_realisations(nmodes=NMODES):
    """residual/HI per HI realisation, for the strategies present in each run.

    Nine runs on the same 119 interior pixels. Six vary only the true HI seed
    (``--hi-seed``, noise seeds fixed); the other three differ because
    ``build_hi_cube`` sizes its box from the union of a run's strategies, so a
    different strategy list is also a different realisation.
    """
    out = {"drift": [], "raster": []}
    for fn in SEED_RUNS:
        p = os.path.join(RESDIR, fn)
        if not os.path.exists(p):
            continue
        e = np.load(p)
        for s_ in out:
            k = f"{s_}_cr_residual_over_hi_{nmodes}"
            if k in e:
                out[s_].append(float(np.median(e[k])))
    return {k: np.asarray(v) for k, v in out.items()}


def fig_realisations(vals, text=None):
    """How much of a single-realisation number is the realisation?

    Left: every run, per strategy, with mean and +/-1 sd. Right: the PAIRED
    drift/raster ratio, which is tighter than either absolute because the
    realisation partly cancels when both strategies see the same HI field.
    That is the practical conclusion -- quote ratios between strategies, not
    absolute residuals.
    """
    dr, ra = vals["drift"], vals["raster"]
    if len(dr) < 3:
        return None
    fig, (axA, axB) = plt.subplots(
        1, 2, figsize=P._figsize("hi_realisations", (11.8, 4.8)),
        gridspec_kw=dict(width_ratios=[1.35, 1.0], wspace=0.30))

    P._tidy(axA); axA.set_yscale("log")
    rng = np.random.default_rng(0)
    for i, (s_, v) in enumerate((("drift", dr), ("raster", ra))):
        c = STRAT[s_]["color"]
        x = i + 1 + rng.uniform(-0.09, 0.09, len(v))
        axA.plot(x, v, ls="none", marker="o", ms=7, mfc=SURFACE, mec=c,
                 mew=1.8, zorder=4)
        m, sd = v.mean(), v.std(ddof=1)
        axA.plot([i + 1 - 0.28, i + 1 + 0.28], [m, m], color=c, lw=2.6,
                 zorder=5)
        axA.add_patch(Rectangle((i + 1 - 0.28, m - sd), 0.56, 2 * sd,
                                facecolor=c, alpha=0.13, edgecolor="none",
                                zorder=2))
        axA.annotate(f"{m:.0f} $\\pm$ {sd:.0f}\n({100 * sd / m:.1f}%)",
                     (i + 1 + 0.33, m), va="center", ha="left", color=c,
                     fontsize=9.2, fontweight="semibold", linespacing=1.4)
        axA.plot([i + 1], [v[0]], marker="*", ms=14, color=INK, zorder=6)
    axA.annotate("$\\star$ = the published run", (0.03, 0.05),
                 xycoords="axes fraction", color=INK, fontsize=8.6)
    axA.set_xlim(0.55, 2.75); axA.set_xticks([1, 2])
    axA.set_xticklabels(["drift", "raster"])
    axA.set_ylabel(f"post-clean residual / HI   ({NMODES} modes removed)")
    axA.set_title(f"{len(dr)} HI realisations, same pixels, same noise seeds",
                  color=INK2, fontsize=9.6, loc="left", pad=6)

    P._tidy(axB)
    r = dr / ra
    x = 1 + rng.uniform(-0.09, 0.09, len(r))
    axB.plot(x, r, ls="none", marker="o", ms=7, mfc=SURFACE, mec=THIRD,
             mew=1.8, zorder=4)
    m, sd = r.mean(), r.std(ddof=1)
    axB.plot([0.72, 1.28], [m, m], color=THIRD, lw=2.6, zorder=5)
    axB.add_patch(Rectangle((0.72, m - sd), 0.56, 2 * sd, facecolor=THIRD,
                            alpha=0.13, edgecolor="none", zorder=2))
    axB.annotate(f"{m:.1f} $\\pm$ {sd:.1f}\n({100 * sd / m:.1f}%)",
                 (1.35, m), va="center", ha="left", color=THIRD, fontsize=9.6,
                 fontweight="semibold", linespacing=1.4)
    axB.plot([1], [r[0]], marker="*", ms=14, color=INK, zorder=6)
    axB.set_xlim(0.5, 2.1); axB.set_xticks([1])
    axB.set_xticklabels(["drift / raster"])
    axB.set_ylabel("ratio between strategies")
    axB.set_title("The paired ratio is tighter than either number",
                  color=INK2, fontsize=9.6, loc="left", pad=6)

    P._title(fig, "A single-realisation residual carries a 9$-$13% error bar",
             "Only the true HI seed varies; noise seeds and pixels are fixed, so "
             "this is realisation scatter alone. It is the same size as several "
             "effects this series has reported, which is why the cost-matched "
             "ladder test came out inconclusive. Quote ratios, not absolutes.",
             x=0.012, y_title=1.05, y_sub=0.975, override=text)
    fig.tight_layout(rect=(0, 0, 1, 0.89))
    return P._save(fig, "hi_realisations")


if __name__ == "__main__":
    main()
