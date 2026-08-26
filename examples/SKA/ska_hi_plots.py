"""Experiment 006 figures: the non-recoverability result, made presentable.

Four figures, all drawn from cached results so nothing is re-solved:

    figures/hi_transfer_function.png   T(k_par), drift vs raster, per mode count
    figures/hi_residual_vs_modes.png   how far the residual sits above the HI
    figures/hi_ablation.png            floor vs total vs noise-only
    figures/hi_frequency_structure.png why the drift fails

Reads ``results/hi_experiment_*.npz`` (figures 1, 2, 4) and
``results/hi_ablation_*.npz`` (figure 3).

    /home/geoff/gibbs_venv_312/bin/python ska_hi_plots.py

Palette is the validated categorical default: slot 1 blue for the drift, slot 2
orange for the raster, slot 3 aqua for the third ablation arm. Those three
validate all-pairs for colourblind separation (worst min(protan, deutan) OKLab
dE 9.2, worst normal-vision dE 24.0, light surface), so identity never rests on
hue alone -- every series is also directly labelled or legended.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import ska_hi_analysis as A
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

_HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(_HERE, "figures")
RESDIR = os.path.join(_HERE, "results")

# --- palette (validated categorical slots 1-3 + chart chrome) ---------------
DRIFT, RASTER, THIRD = "#2a78d6", "#eb6834", "#1baf7a"
SURFACE = "#fcfcfb"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, BASELINE = "#e1e0d9", "#c3c2b7"

STRAT = {"drift": dict(color=DRIFT, label="Drift (parked, no cross-linking)"),
         "raster": dict(color=RASTER, label="Raster (cross-linked, ~74$\\degree$)")}

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "font.size": 9.5, "axes.labelsize": 10, "axes.titlesize": 10.5,
    "axes.edgecolor": BASELINE, "axes.labelcolor": INK2,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "legend.fontsize": 9,
    "lines.linewidth": 2.0, "figure.dpi": 110, "savefig.dpi": 200,
})


def _tidy(ax):
    ax.set_axisbelow(True)
    ax.tick_params(length=3, width=0.6)
    return ax


def _load():
    exp = np.load(os.path.join(
        RESDIR, "hi_experiment_f350_400_nc32_ns64.npz"))
    abl_path = os.path.join(RESDIR, "hi_ablation_f350_400_nc32_ns64.npz")
    abl = np.load(abl_path) if os.path.exists(abl_path) else None
    rk_path = os.path.join(RESDIR, "hi_rank_f350_400_nc32_ns64.npz")
    rk = np.load(rk_path) if os.path.exists(rk_path) else None
    return exp, abl, rk


# ---------------------------------------------------------------------------
# Figure 1 -- the transfer function
# ---------------------------------------------------------------------------

def fig_transfer_function(exp, var="cr_"):
    """T(k_par). Defaults to the common-resolution, interior-pixel variant --
    the pipeline a real analysis would run."""
    nmodes = [int(n) for n in exp["nmodes_grid"]]
    k = exp[f"drift_{var}k"]
    fig, axes = plt.subplots(2, 4, figsize=(13.4, 6.4), sharex=True, sharey=True)

    for ax, nm in zip(axes.ravel(), nmodes):
        _tidy(ax)
        ax.axhline(1.0, color=BASELINE, lw=1.0, ls=(0, (4, 3)), zorder=1)
        for s, style in STRAT.items():
            t = exp[f"{s}_{var}tf_{nm}"]
            e = exp[f"{s}_{var}tf_scatter_{nm}"]
            ax.fill_between(k, t - e, t + e, color=style["color"], alpha=0.16,
                            lw=0, zorder=2)
            ax.plot(k, t, color=style["color"], zorder=3,
                    solid_capstyle="round")
        ax.set_xscale("log")
        ax.set_title(f"{nm} mode{'s' if nm > 1 else ''} removed",
                     color=INK, pad=6)
        ax.set_ylim(-0.05, 1.15)
        ax.set_xlim(k.min() * 0.9, k.max() * 1.1)
        # The cell below the last top-row panel holds the legend, so that panel
        # has to carry its own x axis or the column loses it entirely.
        if ax is axes[0, -1]:
            ax.tick_params(labelbottom=True)
            ax.set_xlabel("$k_\\parallel$  [Mpc$^{-1}$]")

    # Direct labels on the first panel, so identity is not colour-alone.
    a0 = axes[0, 0]
    a0.annotate("drift", (k[7], exp[f"drift_{var}tf_{nmodes[0]}"][7]),
                textcoords="offset points", xytext=(0, -16),
                color=DRIFT, fontsize=9, fontweight="semibold", ha="center")
    a0.annotate("raster", (k[7], exp[f"raster_{var}tf_{nmodes[0]}"][7]),
                textcoords="offset points", xytext=(0, 8),
                color=RASTER, fontsize=9, fontweight="semibold", ha="center")

    # Spare cell carries the legend and the reading instruction.
    spare = axes.ravel()[len(nmodes)]
    spare.set_axis_off()
    spare.legend(handles=[Line2D([], [], color=v["color"], lw=2.4,
                                 label=v["label"]) for v in STRAT.values()]
                 + [Line2D([], [], color=BASELINE, lw=1.0, ls=(0, (4, 3)),
                           label="$T=1$, no signal lost")],
                 loc="upper left", bbox_to_anchor=(-0.02, 0.95),
                 handlelength=1.8, labelspacing=0.9)
    spare.text(-0.02, 0.30,
               "Band = $\\pm1\\sigma$ over 20\nindependent HI mocks.\n\n"
               "Loss concentrates at low\n$k_\\parallel$: those modes are the\n"
               "smoothest along the line of\nsight, so PCA mistakes them\n"
               "for foreground.",
               transform=spare.transAxes, va="top", color=INK2, fontsize=8.6,
               linespacing=1.5)

    for ax in axes[1, :len(nmodes) - 4]:
        ax.set_xlabel("$k_\\parallel$  [Mpc$^{-1}$]")
    # One shared y label — per-axes labels on a shared axis collide between rows.
    fig.supylabel("$T(k_\\parallel)$   —   fraction of HI power surviving the clean",
                  color=INK2, fontsize=10, x=0.005)

    fig.suptitle("Cross-linking preserves HI through foreground cleaning at every scale",
                 color=INK, fontsize=13, fontweight="semibold", x=0.055,
                 ha="left", y=0.985)
    fig.text(0.055, 0.938,
             "350$-$400 MHz, 32 channels, nside 64. Channels reconvolved to a "
             "common 3.99$\\degree$ beam, 119 interior pixels. $T$ is each "
             "strategy's ratio against its own injected response, so it is "
             "depth-independent.",
             color=INK2, fontsize=9.4, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.925))
    out = os.path.join(FIGDIR, "hi_transfer_function.png")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# Figure 2 -- residual against the signal
# ---------------------------------------------------------------------------

def fig_residual_vs_modes(exp):
    """Residual/HI against modes removed, as run and after reconvolution.

    Both variants are on the same 119 interior pixels, so the improvement is
    the reconvolution alone and not the pixel-set change it forces.
    """
    nmodes = [int(n) for n in exp["nmodes_grid"]]
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    _tidy(ax)

    VAR = [("int_", (0, (5, 2)), "none", "as run"),
           ("cr_", "-", SURFACE, "common resolution")]
    for s, style in STRAT.items():
        for var, ls, mfc, _lab in VAR:
            med = [np.median(np.sqrt(exp[f"{s}_{var}p_clean_auto_{nm}"]
                                     / exp[f"{s}_{var}p_mapmade"]))
                   for nm in nmodes]
            ax.plot(nmodes, med, color=style["color"], ls=ls, marker="o", ms=6,
                    mfc=mfc, mec=style["color"], mew=2.0,
                    alpha=1.0 if var == "cr_" else 0.55, zorder=3)
        ax.annotate(s, (nmodes[-1], med[-1]), textcoords="offset points",
                    xytext=(12, 0), color=style["color"], fontsize=9.5,
                    fontweight="semibold", va="center")

    ax.axhline(1.0, color=INK, lw=1.4)
    ax.annotate("HI level — residual would have to reach here",
                (nmodes[0], 1.0), textcoords="offset points", xytext=(2, 7),
                color=INK, fontsize=8.8)
    ax.set_yscale("log")
    ax.set_xlabel("PCA modes removed")
    ax.set_ylabel("post-clean residual / HI   (amplitude)")
    ax.set_xticks(nmodes); ax.set_xticklabels([str(n) for n in nmodes])
    ax.set_xlim(0.4, nmodes[-1] + 2.6)
    ax.legend(handles=[Line2D([], [], color=MUTED, ls=ls, marker="o", ms=6,
                              mfc=(SURFACE if var == "cr_" else "none"),
                              mec=MUTED, mew=2.0, label=lab)
                       for var, ls, _m, lab in VAR],
              loc="center left", bbox_to_anchor=(0.02, 0.30),
              handlelength=2.6, labelspacing=0.8)
    ax.set_title("Reconvolving to a common beam halves the raster's residual, "
                 "and does nothing for the drift",
                 color=INK, fontsize=12, fontweight="semibold", loc="left",
                 pad=26)
    ax.text(0, 1.045,
            "Both variants on the same 119 interior pixels. Even corrected, the "
            "cross-linked raster is still 13$-$30$\\times$ above the HI.",
            transform=ax.transAxes, color=INK2, fontsize=9.2)
    fig.tight_layout()
    out = os.path.join(FIGDIR, "hi_residual_vs_modes.png")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# Figure 3 -- the ablation
# ---------------------------------------------------------------------------

def fig_ablation(abl):
    nmodes = [int(n) for n in abl["nmodes_grid"]]
    arms = [("floor", DRIFT, "Beam + prior floor (noiseless data)", "-"),
            ("total", INK, "Full data (floor + 1/f + white)", "none"),
            ("noiseonly", THIRD, "Noise alone (1/f + white)", "-")]
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.9), sharey=True)

    for ax, s in zip(axes, ("drift", "raster")):
        _tidy(ax)
        p_hi = abl[f"{s}_p_hi"]
        for arm, colour, label, ls in arms:
            med = [np.median(np.sqrt(abl[f"{s}_{arm}_{nm}"] / p_hi))
                   for nm in nmodes]
            if arm == "total":
                # Open rings on top of `floor`: the point is that they coincide.
                ax.plot(nmodes, med, ls="none", marker="o", ms=9,
                        mfc="none", mec=colour, mew=1.6, zorder=4)
            else:
                ax.plot(nmodes, med, color=colour, ls=ls, marker="o", ms=5,
                        mfc=SURFACE, mew=1.8, zorder=3)
        ax.axhline(1.0, color=INK, lw=1.4, zorder=2)
        ax.set_yscale("log")
        ax.set_xticks(nmodes); ax.set_xticklabels([str(n) for n in nmodes])
        ax.set_xlabel("PCA modes removed")
        ax.set_title(STRAT[s]["label"], color=STRAT[s]["color"], pad=6,
                     fontweight="semibold")
    axes[0].set_ylabel("cleaned residual / HI   (amplitude)")
    axes[0].annotate("HI level", (nmodes[0], 1.0), textcoords="offset points",
                     xytext=(2, 7), color=INK, fontsize=8.8)

    # Figure-level and horizontal: inside either panel it collides with the
    # noise curve or the HI-level line.
    fig.legend(handles=[Line2D([], [], color=c, ls=(ls if ls != "none" else "none"),
                               marker="o", ms=(9 if a == "total" else 5),
                               mfc=("none" if a == "total" else SURFACE),
                               mec=c, mew=1.6, lw=(0 if ls == "none" else 2),
                               label=lab)
                        for a, c, lab, ls in arms],
               loc="upper left", bbox_to_anchor=(0.045, 0.90), ncols=3,
               handlelength=2.2, columnspacing=2.4)

    fig.suptitle("The blocker is the floor, not the noise",
                 color=INK, fontsize=13, fontweight="semibold", x=0.045,
                 ha="left", y=1.03)
    fig.text(0.045, 0.965,
             "Noiseless and full data coincide to 0.2% in power — the rings sit "
             "on the line. Noise alone is ~38$\\times$ below the floor, so it is "
             "irrelevant to HI recovery.",
             color=INK2, fontsize=9.4, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.84))
    out = os.path.join(FIGDIR, "hi_ablation.png")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# Figure 4 -- why the drift fails
# ---------------------------------------------------------------------------

def fig_frequency_structure(exp):
    freqs = exp["freqs"]
    k = exp["drift_k"]
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.0),
                             gridspec_kw=dict(width_ratios=[1, 1.15]))

    # (a) radial power retained by the map-maker
    ax = _tidy(axes[0])
    for s, style in STRAT.items():
        ratio = exp[f"{s}_p_mapmade"] / exp[f"{s}_p_true"]
        ax.plot(k, ratio, color=style["color"], marker="o", ms=4,
                mfc=SURFACE, mew=1.4, label=s)
        # Label at the left end, where the two curves are furthest apart and
        # there is clear space; at the right end they collide with the markers.
        ax.annotate(s, (k[0], ratio[0]), textcoords="offset points",
                    xytext=(6, 10), color=style["color"], fontsize=9.5,
                    fontweight="semibold", ha="left")
    ax.set_xscale("log")
    ax.set_ylim(0, 0.42)
    ax.set_xlabel("$k_\\parallel$  [Mpc$^{-1}$]")
    ax.set_ylabel("$P_\\mathrm{map\\,made}\\,/\\,P_\\mathrm{true}$")
    ax.set_title("Radial HI power surviving the map-maker", color=INK, pad=6)

    # (b) the frequency-structure split, which is the actual mechanism.
    # A single line of sight was tried here and did not show it -- the traces
    # are too noisy to read "flat in frequency" off by eye. The two ratios
    # state it directly, and share an axis because both are map-made/true.
    ax = _tidy(axes[1])
    metrics, groups = [], ["frequency-varying\npart retained",
                           "frequency-coherent\npart (the mean)"]
    for s in STRAT:
        t, m = exp[f"{s}_hi_true"], exp[f"{s}_hi_mapmade"]
        varying = ((m - m.mean(axis=0)).std() / (t - t.mean(axis=0)).std())
        coherent = m.mean(axis=0).std() / t.mean(axis=0).std()
        metrics.append((varying, coherent))
    x = np.arange(2)
    w = 0.34
    for i, (s, style) in enumerate(STRAT.items()):
        off = (i - 0.5) * (w + 0.02)
        bars = ax.bar(x + off, metrics[i], width=w, color=style["color"],
                      edgecolor=SURFACE, linewidth=2.0, zorder=3,
                      label=style["label"])
        for b, v in zip(bars, metrics[i]):
            ax.annotate(f"{v:.2f}$\\times$",
                        (b.get_x() + b.get_width() / 2, v),
                        textcoords="offset points", xytext=(0, 4),
                        ha="center", color=INK, fontsize=9.5,
                        fontweight="semibold")
    ax.axhline(1.0, color=INK, lw=1.4, zorder=4)
    # Left of the bars: at the right it collides with the raster bar's label.
    ax.annotate("1.0 — faithful", (-0.47, 1.0), textcoords="offset points",
                xytext=(0, 6), color=INK, fontsize=8.8, ha="left")
    ax.set_xticks(x); ax.set_xticklabels(groups, color=INK2)
    ax.set_ylabel("map-made / true   (rms ratio)")
    ax.set_ylim(0, 3.7)
    ax.set_xlim(-0.55, 1.55)   # explicit, else the reference-line label clips
    ax.set_title("The drift suppresses structure and inflates the mean",
                 color=INK, pad=6)
    ax.grid(axis="x", visible=False)
    ax.legend(loc="upper left", bbox_to_anchor=(0.015, 0.99))

    fig.suptitle("Why the drift fails: it turns HI into something that looks like foreground",
                 color=INK, fontsize=13, fontweight="semibold", x=0.045,
                 ha="left", y=1.01)
    fig.text(0.045, 0.945,
             "The drift keeps only 38% of the frequency-varying HI while "
             "amplifying the frequency-coherent part 3.14$\\times$ "
             "(raster: 0.53 and 1.11$\\times$). With ~17 measured modes it "
             "projects HI onto a nearly frequency-constant subspace.",
             color=INK2, fontsize=9.4, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    out = os.path.join(FIGDIR, "hi_frequency_structure.png")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# Figures 5 and 6 -- the map domain
# ---------------------------------------------------------------------------

# Colour maps come from the same validated palette as the line charts.
# Signed fields (HI fluctuations, residuals) get the DIVERGING pair -- two
# opposite-reading poles with a NEUTRAL GREY midpoint, so "zero" recedes.
# Positive-definite fields (sky temperature) get the one-hue SEQUENTIAL ramp,
# light to dark. Never a rainbow: healpy's default would encode magnitude as
# hue and invent structure that is not in the data.
def _cmaps():
    from matplotlib.colors import LinearSegmentedColormap
    div = LinearSegmentedColormap.from_list(
        "blue_grey_red", ["#184f95", "#2a78d6", "#f0efec", "#e34948", "#a02b2b"])
    seq = LinearSegmentedColormap.from_list(
        "blues", ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf",
                  "#184f95", "#0d366b"])
    for c in (div, seq):
        c.set_bad(SURFACE)          # unobserved sky reads as page, not as data
    return div, seq


def _project(values, pix, nside, reso, xsize, rot):
    """One HEALPix patch -> a 2D array, with everything unobserved masked."""
    import healpy as hp
    full = np.full(hp.nside2npix(nside), hp.UNSEEN)
    full[np.asarray(pix)] = values
    return hp.gnomview(full, rot=rot, reso=reso, xsize=xsize, no_plot=True,
                       return_projected_map=True)


def _mapshow(ax, img, cmap, vmin, vmax, title, extent):
    im = ax.imshow(img, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax,
                   extent=extent, interpolation="nearest")
    ax.set_title(title, color=INK, fontsize=10, pad=5)
    ax.set_xlabel("offset  [deg]", fontsize=8.5)
    ax.grid(False)
    ax.tick_params(labelsize=8)
    return im


def fig_patch_maps(exp, nside=64, channel=None):
    """Figure 5 -- the field, and which pixels each strategy actually solves."""
    import healpy as hp
    from ska_common import gdsm_equatorial_sky_model

    div, seq = _cmaps()
    common = exp["common"]
    interior = exp["interior"]
    ops = {s: A.channel_operators(s, verbose=False)[0]
           for s in ("drift", "raster")}
    rot = (158.30, 9.375)
    reso, xsize = 8.0, 380
    half = reso * xsize / 60.0 / 2.0
    extent = (-half, half, -half, half)

    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.6), layout="constrained")

    # (a) the sky this experiment actually observes
    allpix = np.union1d(np.asarray(ops["drift"].pixel_indices),
                        np.asarray(ops["raster"].pixel_indices))
    sky = gdsm_equatorial_sky_model(freq=350.0, nside=nside)[allpix]
    img = _project(sky, allpix, nside, reso, xsize, rot)
    im = _mapshow(axes[0], img, seq, np.nanpercentile(sky, 1),
                  np.nanpercentile(sky, 99),
                  "GDSM foreground at 350 MHz", extent)
    cb = fig.colorbar(im, ax=axes[0], shrink=0.62, aspect=15, pad=0.02)
    cb.set_label("$T_b$  [K]", fontsize=8.5); cb.ax.tick_params(labelsize=8)
    axes[0].set_ylabel("offset  [deg]", fontsize=8.5)

    # (b) coverage. Categorical, so it uses the categorical slots, not a ramp.
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.patches import Patch
    code = np.zeros(hp.nside2npix(nside))
    code[np.asarray(ops["raster"].pixel_indices)] = 1
    code[common] = 2
    code[np.asarray(common)[interior]] = 3
    img = _project(code[allpix], allpix, nside, reso, xsize, rot)
    cmap = ListedColormap(["#d8d7d0", RASTER, DRIFT, "#0d366b"])
    axes[1].imshow(img, origin="lower", cmap=cmap,
                   norm=BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], 4),
                   extent=extent, interpolation="nearest")
    axes[1].set_title("What each strategy solves", color=INK, fontsize=10, pad=5)
    axes[1].set_xlabel("offset  [deg]", fontsize=8.5)
    axes[1].grid(False); axes[1].tick_params(labelsize=8)
    n_rast_only = len(ops["raster"].pixel_indices) - len(common)
    axes[1].legend(handles=[
        Patch(facecolor=RASTER, label=f"raster only ({n_rast_only} px)"),
        Patch(facecolor=DRIFT, label=f"common to both ({len(common)} px)"),
        Patch(facecolor="#0d366b",
              label=f"interior, used for common-res ({interior.sum()} px)")],
        loc="upper left", bbox_to_anchor=(0.0, -0.12), fontsize=8.6)
    for ax in axes:                 # the patch is a wide strip; square axes
        ax.set_ylim(-11, 11)        # would be mostly empty page

    fig.suptitle("The field and the two footprints", color=INK, fontsize=13,
                 fontweight="semibold", x=0.01, ha="left", y=1.13)
    fig.text(0.01, 1.04,
             "Gnomonic projection about RA 158.3$\\degree$, Dec +9.4$\\degree$. The drift's "
             "patch is a stack of three constant-Dec strips; the raster's azimuth "
             "throw widens it 2.2$\\times$ and wholly contains it "
             f"({len(ops['raster'].pixel_indices)} px against "
             f"{len(ops['drift'].pixel_indices)}).",
             color=INK2, fontsize=9.4, ha="left")
    out = os.path.join(FIGDIR, "hi_patch_maps.png")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def fig_hi_maps(exp, nside=64, nmodes=4):
    """Figure 6 -- the HI itself, in the map domain.

    The frequency mean is removed from every panel. That is deliberate: the
    frequency-*varying* part is the HI a survey can actually use, and it is the
    part the drift loses (0.38x, while inflating the coherent part 3.14x). Raw
    channel maps would show the drift looking healthy for the wrong reason.
    """
    div, _ = _cmaps()
    common = exp["common"]
    freqs = exp["freqs"]
    ch = len(freqs) // 2
    rot = (158.30, 9.375)
    reso, xsize = 5.0, 330
    half = reso * xsize / 60.0 / 2.0
    extent = (-half, half, -half, half)

    def demean(c):
        return c - c.mean(axis=0, keepdims=True)

    true = demean(exp["drift_hi_true"])[ch] * 1e3          # uK
    dmap = demean(exp["drift_hi_mapmade"])[ch] * 1e3
    rmap = demean(exp["raster_hi_mapmade"])[ch] * 1e3
    clean = demean(A.pca_clean(exp["raster_data_cube"], nmodes))[ch] * 1e3

    v = float(np.nanpercentile(np.abs(true), 99))
    fig, axes = plt.subplots(1, 4, figsize=(16.2, 3.9), layout="constrained")
    panels = [(true, "True HI", v), (dmap, "Drift, map-made", v),
              (rmap, "Raster, map-made", v),
              (clean, f"Raster, cleaned data ({nmodes} modes)",
               float(np.nanpercentile(np.abs(clean), 99)))]
    ims = []
    for ax, (vals, title, vv) in zip(axes, panels):
        img = _project(vals, common, nside, reso, xsize, rot)
        ims.append(_mapshow(ax, img, div, -vv, vv, title, extent))
        ax.set_ylim(-9, 9)
    # One bar for the three panels that share a scale, one for the outlier --
    # four identical bars would imply four different scales.
    cb = fig.colorbar(ims[0], ax=list(axes[:3]), shrink=0.72, aspect=16,
                      pad=0.015)
    cb.set_label("$\\delta T_b$  [$\\mu$K]", fontsize=8.5)
    cb.ax.tick_params(labelsize=8)
    cb2 = fig.colorbar(ims[3], ax=axes[3], shrink=0.72, aspect=16, pad=0.03)
    cb2.set_label("$\\delta T_b$  [$\\mu$K] — 50$\\times$ the scale at left",
                  fontsize=8.5)
    cb2.ax.tick_params(labelsize=8)
    axes[0].set_ylabel("offset  [deg]", fontsize=8.5)

    fig.suptitle("The HI in the map domain, and why you cannot see it",
                 color=INK, fontsize=13, fontweight="semibold", x=0.01,
                 ha="left", y=1.10)
    fig.text(0.01, 1.02,
             f"Channel {ch} ({freqs[ch]:.1f} MHz), frequency mean removed. First "
             "three share a colour scale; the fourth needs its own, and that is "
             "the result — the cleaned map is residual, not signal.",
             color=INK2, fontsize=9.4, ha="left")
    out = os.path.join(FIGDIR, "hi_maps.png")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# Figure 7 -- how compressible is each component?
# ---------------------------------------------------------------------------

def fig_rank(rank, nshow=20):
    """Eigenspectra of the channel-channel covariance.

    The question this answers is whether swapping PCA for ICA/GMCA/NMF could
    help. All of them remove a rank-N subspace, so what matters is how many
    modes each component actually occupies.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.8), layout="constrained")
    series = [("drift_w_fg", INK, "GDSM foreground", "-"),
              ("drift_w_floor", DRIFT, "beam + prior floor, drift", "-"),
              ("raster_w_floor", RASTER, "beam + prior floor, raster", "-"),
              ("drift_w_hi", THIRD, "HI (map-made)", "-")]

    ax = _tidy(axes[0])
    for key, colour, label, ls in series:
        w = rank[key][:nshow]
        ax.plot(np.arange(1, len(w) + 1), np.maximum(w, 1e-18), color=colour,
                ls=ls, marker="o", ms=3.5, mfc=SURFACE, mew=1.2, label=label)
    ax.set_yscale("log")
    ax.set_ylim(1e-18, 3)
    ax.set_xlabel("eigenmode of the channel-channel covariance")
    ax.set_ylabel("$\\lambda_i\\,/\\,\\lambda_1$")
    ax.set_title("Eigenspectrum", color=INK, pad=6)
    ax.legend(loc="center right", fontsize=8.4)  # lower left sits on the GDSM line

    ax = _tidy(axes[1])
    for key, colour, label, ls in series:
        w = rank[key]
        c = np.cumsum(w) / w.sum()
        ax.plot(np.arange(1, min(nshow, len(c)) + 1), c[:nshow], color=colour,
                ls=ls, marker="o", ms=3.5, mfc=SURFACE, mew=1.2)
    for frac, lbl in ((0.90, "90%"), (0.99, "99%")):
        ax.axhline(frac, color=BASELINE, lw=1.0, ls=(0, (4, 3)))
        ax.annotate(lbl, (nshow, frac), textcoords="offset points",
                    xytext=(-2, 4), ha="right", fontsize=8.2, color=MUTED)
    ax.set_ylim(0, 1.04)
    ax.set_xlabel("eigenmodes retained")
    ax.set_ylabel("cumulative fraction of variance")
    ax.set_title("How many modes each component occupies", color=INK, pad=6)

    fig.suptitle("A different cleaning basis cannot help: the floor is already low-rank",
                 color=INK, fontsize=13, fontweight="semibold", x=0.01,
                 ha="left", y=1.10)
    fig.text(0.01, 1.02,
             "The foreground is rank 1 — one mode removes it. The floor needs 6 "
             "(drift) or 2 (raster). The HI needs 23$-$25 of 32: it is the one "
             "component that is not compressible.",
             color=INK2, fontsize=9.4, ha="left")
    out = os.path.join(FIGDIR, "hi_rank.png")
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    exp, abl, rk = _load()
    made = [fig_transfer_function(exp), fig_residual_vs_modes(exp),
            fig_frequency_structure(exp), fig_patch_maps(exp), fig_hi_maps(exp)]
    if abl is not None:
        made.insert(2, fig_ablation(abl))
    else:
        print("!! ablation results missing -- run ska_hi_ablation.py", flush=True)
    if rk is not None:
        made.append(fig_rank(rk))
    else:
        print("!! rank results missing -- run ska_hi_rank.py", flush=True)
    for m in made:
        print(f"wrote {os.path.relpath(m, _HERE)}", flush=True)


if __name__ == "__main__":
    main()
