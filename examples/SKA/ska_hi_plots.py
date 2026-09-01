"""Experiment 006 figures: the non-recoverability result, made presentable.

Eight figures, all drawn from cached results so nothing is re-solved:

    figures/hi_transfer_function.png   T(k_par), drift vs raster, per mode count
    figures/hi_residual_vs_modes.png   how far the residual sits above the HI
    figures/hi_ablation.png            floor vs total vs noise-only
    figures/hi_frequency_structure.png why the drift fails
    figures/hi_patch_maps.png          the field and the two footprints
    figures/hi_maps.png                HI in the map domain
    figures/hi_rank.png                eigenspectra: is the floor compressible?
    figures/hi_eigenvectors.png        the spectral shapes themselves

Reads ``results/hi_experiment_*.npz`` (figures 1, 2, 4, 6),
``results/hi_ablation_*.npz`` (figure 3) and ``results/hi_rank_*.npz``
(figures 7, 8). Figure 5 also loads the cached operators for the footprints.

**These are report figures, not publication figures.** See METHODS.md section 12
for what would have to change -- chiefly that every title and subtitle is baked
into the PNG and would need to move to a LaTeX caption.

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
    # Axes furniture in ink rather than the recessive grey the
    # data-viz default suggests -- override via STYLE["rc"] to go back.
    "axes.edgecolor": INK, "axes.labelcolor": INK,
    "xtick.color": INK, "ytick.color": INK,
    "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "legend.fontsize": 9,
    "lines.linewidth": 2.0, "figure.dpi": 110, "savefig.dpi": 200,
})


# ---------------------------------------------------------------------------
# Everything a caller might want to change, in one place.
# ---------------------------------------------------------------------------
# `ska_hi_figures.ipynb` edits this dict and calls apply_style(), so the
# notebook never duplicates plotting code. For publication the two that matter
# are draw_titles=False (the headline moves into the LaTeX caption) and
# fmt="pdf". See METHODS.md section 12.
STYLE = dict(
    draw_titles=True,      # False strips suptitle + subtitle from every figure
    fmt="png",             # "pdf" for vector text and lines
    dpi=200,
    figsize={},            # per-figure override, e.g. {"hi_rank": (7.0, 3.2)}
    rc={},                 # extra rcParams, merged over the defaults above
)

_BASE_RC = dict(plt.rcParams)


def apply_style():
    """Re-apply rcParams after STYLE['rc'] has been edited."""
    plt.rcParams.update(_BASE_RC)
    plt.rcParams.update(STYLE.get("rc", {}))


def _figsize(name, default):
    return STYLE.get("figsize", {}).get(name, default)


def _title(fig, title, subtitle, x=0.01, y_title=1.05, y_sub=0.99,
           override=None):
    """Headline + subtitle.

    Suppressed entirely when STYLE['draw_titles'] is False. `override` is the
    `text=` argument every fig_* accepts, so a caller can retitle a figure
    without editing this module:

        P.fig_rank(rk, text={'title': '...', 'subtitle': '...'})
    """
    if not STYLE.get("draw_titles", True):
        return
    o = override or {}
    title = o.get("title", title)
    subtitle = o.get("subtitle", subtitle)
    fig.suptitle(title, color=INK, fontsize=13, fontweight="semibold",
                 x=x, ha="left", y=y_title)
    fig.text(x, y_sub, subtitle, color=INK2, fontsize=9.4, ha="left")


def _save(fig, name):
    """Write figures/<name>.<fmt>, honouring STYLE."""
    out = os.path.join(FIGDIR, f"{name}.{STYLE.get('fmt','png')}")
    fig.savefig(out, bbox_inches="tight", dpi=STYLE.get("dpi", 200))
    plt.close(fig)
    return out


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

def fig_transfer_function(exp, var="cr_", text=None):
    """T(k_par). Defaults to the common-resolution, interior-pixel variant --
    the pipeline a real analysis would run."""
    nmodes = [int(n) for n in exp["nmodes_grid"]]
    k = exp[f"drift_{var}k"]
    fig, axes = plt.subplots(2, 4, figsize=_figsize("hi_transfer_function", (13.4, 6.4)), sharex=True, sharey=True)

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

    _title(fig, "Cross-linking preserves HI through foreground cleaning at every scale",
             "350$-$400 MHz, 32 channels, nside 64. Channels reconvolved to a "
             "common 3.99$\\degree$ beam, 119 interior pixels. $T$ is each "
             "strategy's ratio against its own injected response, so it is "
             "depth-independent.",
           x=0.055, y_title=0.985, y_sub=0.938, override=text)
    fig.tight_layout(rect=(0, 0, 1, 0.925))
    return _save(fig, "hi_transfer_function")


# ---------------------------------------------------------------------------
# Figure 2 -- residual against the signal
# ---------------------------------------------------------------------------

def fig_residual_vs_modes(exp, text=None):
    """Residual/HI against modes removed, as run and after reconvolution.

    Both variants are on the same 119 interior pixels, so the improvement is
    the reconvolution alone and not the pixel-set change it forces.
    """
    nmodes = [int(n) for n in exp["nmodes_grid"]]
    fig, ax = plt.subplots(figsize=_figsize("hi_residual_vs_modes", (8.6, 5.6)))
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
    if STYLE.get("draw_titles", True):   # axes-level here, not fig-level
        _o = text or {}
        ax.set_title(_o.get("title",
                     "Reconvolving to a common beam halves the raster's "
                     "residual, and does nothing for the drift"),
                     color=INK, fontsize=12, fontweight="semibold", loc="left",
                     pad=26)
        ax.text(0, 1.045,
                _o.get("subtitle",
                       "Both variants on the same 119 interior pixels. Even "
                       "corrected, the cross-linked raster is still "
                       "13$-$30$\\times$ above the HI."),
                transform=ax.transAxes, color=INK2, fontsize=9.2)
    fig.tight_layout()
    return _save(fig, "hi_residual_vs_modes")


# ---------------------------------------------------------------------------
# Figure 3 -- the ablation
# ---------------------------------------------------------------------------

def fig_ablation(abl, text=None, tag="cr_"):
    """`tag` selects the convention: "cr_" (common resolution, interior pixels,
    matching figure 2), "int_" (as run, interior pixels) or "" (as run, full
    patch -- the originally published numbers). Older result files hold only
    the "" keys, so fall back to those rather than raising."""
    if f"drift_{tag}p_hi" not in abl:
        print(f"!! ablation results predate the '{tag}' variant -- "
              f"using as-run full-patch keys", flush=True)
        tag = ""
    nmodes = [int(n) for n in abl["nmodes_grid"]]
    # key, colour, label, linestyle, marker, markersize -- one marker per line
    # so the arms stay separable in greyscale and for CVD readers.
    #
    # FOURTH is categorical slot 4 of the validated palette (the next free slot:
    # slot 2, orange, is reserved for "raster" across this figure set, so using
    # it for a noise term inside a per-strategy panel would clash). Checked with
    # the skill's validator ported to Python -- blue/aqua/yellow on all pairs:
    # CVD dE 9.1 (target 8), normal-vision dE 22.9 (floor 15). Yellow's contrast
    # against the surface is 2.11, a WARN, which is why every line is legended
    # AND carries its own marker rather than resting on hue.
    FOURTH = "#eda100"
    arms = [("floor", DRIFT, "Beam + prior floor (noiseless data)", "-", "o", 5),
            ("total", INK, "Full data (floor + 1/f + white)", "none", "x", 9),
            ("gainonly", THIRD, "1/f gain alone", "-", "^", 6),
            ("whiteonly", FOURTH, "White noise alone", "-", "s", 5)]
    fig, axes = plt.subplots(1, 2, figsize=_figsize("hi_ablation", (11.6, 4.9)), sharey=True)

    for ax, s in zip(axes, ("drift", "raster")):
        _tidy(ax)
        p_hi = abl[f"{s}_{tag}p_hi"]
        for arm, colour, label, ls, mk, ms in arms:
            med = [np.median(np.sqrt(abl[f"{s}_{tag}{arm}_{nm}"] / p_hi))
                   for nm in nmodes]
            if arm == "total":
                # Markers only, over `floor`: the point is that they coincide.
                ax.plot(nmodes, med, ls="none", marker=mk, ms=ms,
                        mfc="none", mec=colour, mew=1.6, zorder=4)
            else:
                ax.plot(nmodes, med, color=colour, ls=ls, marker=mk, ms=ms,
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
                               marker=mk, ms=ms,
                               mfc=("none" if a == "total" else SURFACE),
                               mec=c, mew=1.6, lw=(0 if ls == "none" else 2),
                               label=lab)
                        for a, c, lab, ls, mk, ms in arms],
               loc="upper left", bbox_to_anchor=(0.045, 0.90), ncols=4,
               handlelength=2.0, columnspacing=1.8)

    _title(fig, "The blocker is the floor, and 1/f is the smallest term in it",
             "Noiseless and full data track each other to a few per cent. Split "
             "apart, 1/f sits below the white noise, which is itself far below "
             "the floor — so neither is what limits HI recovery.",
           x=0.045, y_title=1.03, y_sub=0.965, override=text)
    fig.tight_layout(rect=(0, 0, 1, 0.84))
    return _save(fig, "hi_ablation")


# ---------------------------------------------------------------------------
# Figure 4 -- why the drift fails
# ---------------------------------------------------------------------------

def fig_frequency_structure(exp, text=None):
    freqs = exp["freqs"]
    k = exp["drift_k"]
    fig, axes = plt.subplots(1, 2, figsize=_figsize("hi_frequency_structure", (12.4, 5.0)),
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

    _title(fig, "Why the drift fails: it turns HI into something that looks like foreground",
             "The drift keeps only 38% of the frequency-varying HI while "
             "amplifying the frequency-coherent part 3.14$\\times$ "
             "(raster: 0.53 and 1.11$\\times$). With ~17 measured modes it "
             "projects HI onto a nearly frequency-constant subspace.",
           x=0.045, y_title=1.01, y_sub=0.945, override=text)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    return _save(fig, "hi_frequency_structure")


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


def _data_window(img, extent, pad=1.0):
    """(xlim, ylim) enclosing the observed pixels, with a small margin."""
    good = np.isfinite(img) & (img != hp_unseen())
    rows = np.flatnonzero(good.any(axis=1))
    cols = np.flatnonzero(good.any(axis=0))
    if not len(rows) or not len(cols):
        return (extent[0], extent[1]), (extent[2], extent[3])
    x0, x1, y0, y1 = extent
    dx = (x1 - x0) / img.shape[1]
    dy = (y1 - y0) / img.shape[0]
    return ((x0 + dx * cols[0] - pad, x0 + dx * (cols[-1] + 1) + pad),
            (y0 + dy * rows[0] - pad, y0 + dy * (rows[-1] + 1) + pad))


def hp_unseen():
    import healpy as hp
    return hp.UNSEEN


def _mapshow(ax, img, cmap, vmin, vmax, title, extent):
    im = ax.imshow(img, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax,
                   extent=extent, interpolation="nearest")
    ax.set_title(title, color=INK, fontsize=10, pad=5)
    ax.set_xlabel("offset  [deg]", fontsize=8.5)
    ax.grid(False)
    ax.tick_params(labelsize=8)
    return im


def fig_patch_maps(exp, nside=64, channel=None, text=None):
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

    fig, axes = plt.subplots(1, 2, figsize=_figsize("hi_patch_maps", (11.6, 4.6)), layout="constrained")

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

    _title(fig, "The field and the two footprints",
             "Gnomonic projection about RA 158.3$\\degree$, Dec +9.4$\\degree$. The drift's "
             "patch is a stack of three constant-Dec strips; the raster's azimuth "
             "throw widens it 2.2$\\times$ and wholly contains it "
             f"({len(ops['raster'].pixel_indices)} px against "
             f"{len(ops['drift'].pixel_indices)}).",
           x=0.01, y_title=1.13, y_sub=1.04, override=text)
    return _save(fig, "hi_patch_maps")


def fig_hi_maps(exp, nside=64, nmodes=4, text=None, variant="cr_"):
    """Figure 6 -- the HI itself, in the map domain.

    The frequency mean is removed from every panel. That is deliberate: the
    frequency-*varying* part is the HI a survey can actually use, and it is the
    part the drift loses (0.38x, while inflating the coherent part 3.14x). Raw
    channel maps would show the drift looking healthy for the wrong reason.

    `variant="cr_"` reconvolves to common resolution and keeps the 119 interior
    pixels, matching every quantitative figure; `"as_run"` is the full 277-pixel
    patch at native per-channel resolution, as first published. Each panel is a
    single channel, so three of the four are internally consistent either way --
    but the cleaned panel is not, because PCA ran *across* frequency, and on
    as-run data it shows a residual the current pipeline no longer produces.

    Reconvolving costs a factor ~3 in true HI amplitude (0.699 -> 0.223 uK at
    the 99th percentile): the signal is small-scale, so smoothing to the widest
    beam in the band eats most of it while the residual only falls ~1.8x. The
    colour-scale ratio therefore *rises*, 51x -> 82x. That is the real cost of
    a standard procedure, not a defect in the figure.

    That ratio is a 99th-percentile *display* scale at one channel, and is not
    the residual/HI ratio of Figures 1-3 (a median over k of a power ratio,
    24x for the raster). They measure different things and will never agree --
    hence the label says "the shared scale", not a physical ratio.
    """
    div, _ = _cmaps()
    common = np.asarray(exp["common"])
    freqs = exp["freqs"]
    ch = len(freqs) // 2
    rot = (158.30, 9.375)
    reso, xsize = 5.0, 330
    half = reso * xsize / 60.0 / 2.0
    extent = (-half, half, -half, half)

    if variant in ("cr_", "cr"):
        interior = np.asarray(exp["interior"], bool)
        pix = common[interior]
        def prep(cube):
            return A.common_resolution(np.asarray(cube), freqs, common,
                                       nside)[:, interior]
    else:
        pix = common
        def prep(cube):
            return np.asarray(cube)

    def demean(c):
        return c - c.mean(axis=0, keepdims=True)

    true = demean(prep(exp["drift_hi_true"]))[ch] * 1e3          # uK
    dmap = demean(prep(exp["drift_hi_mapmade"]))[ch] * 1e3
    rmap = demean(prep(exp["raster_hi_mapmade"]))[ch] * 1e3
    # Reconvolve first, then clean -- the order run_hi_experiment.py uses.
    clean = demean(A.pca_clean(prep(exp["raster_data_cube"]), nmodes))[ch] * 1e3

    v = float(np.nanpercentile(np.abs(true), 99))
    v_clean = float(np.nanpercentile(np.abs(clean), 99))
    # Explicit gridspec rather than plt.subplots: a colorbar attached to one
    # panel steals space from that panel alone, leaving the 2x2 visibly
    # unequal. Dedicated cells for both bars keep the four panels identical.
    fig = plt.figure(figsize=_figsize("hi_maps", (9.6, 5.9)),
                     layout="constrained")
    gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 0.055],
                          width_ratios=[1, 1, 0.035])
    flat = [fig.add_subplot(gs[r, c]) for r in (0, 1) for c in (0, 1)]
    cax_shared = fig.add_subplot(gs[0:2, 2])
    cax_clean = fig.add_subplot(gs[2, 1])
    panels = [(true, "True HI", v), (dmap, "Drift, map-made", v),
              (rmap, "Raster, map-made", v),
              (clean, f"Raster, cleaned data ({nmodes} modes)", v_clean)]
    ims = []
    win = None
    for ax, (vals, title, vv) in zip(flat, panels):
        img = _project(vals, pix, nside, reso, xsize, rot)
        ims.append(_mapshow(ax, img, div, -vv, vv, title, extent))
        # Size the window to the patch rather than fixing it: the interior is
        # much shorter than the full common patch, and dead paper around it
        # doubles in a 2x2 grid.
        win = win or _data_window(img, extent)
        ax.set_xlim(*win[0]); ax.set_ylim(*win[1])
    # Axis furniture only on the outer edges of the grid.
    for ax in flat[:2]:
        ax.set_xlabel("")
        ax.tick_params(labelbottom=False)   # bottom row carries the axis
    for ax in (flat[1], flat[3]):
        ax.tick_params(labelleft=False)
    for ax in (flat[0], flat[2]):
        ax.set_ylabel("offset  [deg]", fontsize=8.5)
    # One bar for the three panels that share a scale, one for the outlier --
    # four identical bars would imply four different scales. The odd one is
    # horizontal under its own panel so the split reads before the numbers do.
    cb = fig.colorbar(ims[0], cax=cax_shared)
    cb.set_label("$\\delta T_b$  [$\\mu$K]", fontsize=8.5)
    cb.ax.tick_params(labelsize=8)
    cb2 = fig.colorbar(ims[3], cax=cax_clean, orientation="horizontal")
    # Derived, not hardcoded: the ratio moves with `nmodes` (51x at 4 modes,
    # 186x at 2, 18x at 8), and a fixed "50x" would quietly go wrong.
    cb2.set_label(f"$\\delta T_b$  [$\\mu$K] — {v_clean / v:.0f}$\\times$ "
                  "the shared scale", fontsize=8.5)
    cb2.ax.tick_params(labelsize=8)

    _title(fig, "The HI in the map domain, and why you cannot see it",
             f"Channel {ch} ({freqs[ch]:.1f} MHz), frequency mean removed. First "
             "three share a colour scale; the fourth needs its own, and that is "
             "the result — the cleaned map is residual, not signal.",
           x=0.01, y_title=1.10, y_sub=1.02, override=text)
    return _save(fig, "hi_maps")


# ---------------------------------------------------------------------------
# Figure 7 -- how compressible is each component?
# ---------------------------------------------------------------------------

def fig_rank(rank, nshow=20, text=None):
    """Eigenspectra of the channel-channel covariance.

    The question this answers is whether swapping PCA for ICA/GMCA/NMF could
    help. All of them remove a rank-N subspace, so what matters is how many
    modes each component actually occupies.
    """
    fig, axes = plt.subplots(1, 2, figsize=_figsize("hi_rank", (12.4, 4.8)), layout="constrained")
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

    _title(fig, "A different cleaning basis cannot help: the floor is already low-rank",
             "The foreground is rank 1 — one mode removes it. The floor needs 6 "
             "(drift) or 2 (raster). The HI needs 23$-$25 of 32: it is the one "
             "component that is not compressible.",
           x=0.01, y_title=1.10, y_sub=1.02, override=text)
    return _save(fig, "hi_rank")


def fig_1f_scan(scan, abl=None, text=None, cutoff=None):
    """Figure 9 -- how bad would 1/f have to be before it mattered?

    The ablation fixes 1/f at one parameter set; this asks what it would take
    for that term to reach the beam + prior floor. x is the knee frequency --
    where the 1/f gain power crosses the white gain power -- so the fiducial
    model is one point on it and everything to the right is "worse than
    assumed". The floor is drawn as the thing 1/f would have to reach.
    """
    knees = np.asarray(scan["knees_mhz"], float)
    alphas = [float(a) for a in scan["alphas"]]
    nm = int(scan["nmodes"])
    # Slope is an ordinal family, not an identity: one hue, light to dark,
    # rather than three categorical slots that would imply unrelated series.
    ramp = ["#86b6ef", "#2a78d6", "#184f95"]
    nrow = 2 if cutoff is not None else 1
    fig, axes = plt.subplots(nrow, 2, squeeze=False,
                             figsize=_figsize("hi_1f_scan",
                                              (11.6, 4.9 * nrow)))

    for ax, s in zip(axes[0], ("drift", "raster")):
        _tidy(ax)
        # Honest encoding: the knee axis is an exact scaling, not a set of
        # independent runs. Scaling one realisation scales the whole pipeline
        # linearly -- `solve` is affine and differences out, and `pca_clean` is
        # scale-invariant (scaling a cube scales its covariance, leaving the
        # eigenvectors alone) -- so ratio ~ knee**(alpha/2) exactly, verified
        # against the measured grid to 2e-4. So: LINE for the scaling, and a
        # filled MARKER only at the fiducial knee, which is the measured point.
        kfid = float(scan["fiducial_knee_mhz"])
        for a, colour in zip(alphas, ramp):
            y = np.array([scan[f"{s}_ratio_a{a}_k{k}"] for k in knees], float)
            ax.plot(knees, y, color=colour, lw=2.0, zorder=3,
                    label=f"$\\alpha = {a}$")
            j = int(np.argmin(np.abs(knees - kfid)))
            ax.plot([knees[j]], [y[j]], color=colour, marker="o", ms=7,
                    mfc=colour, mec=SURFACE, mew=1.6, zorder=5)
        # What 1/f is being compared against: the floor at the same mode count.
        if abl is not None:
            p_hi = abl[f"{s}_cr_p_hi"]
            floor = np.median(np.sqrt(abl[f"{s}_cr_floor_{nm}"] / p_hi))
            ax.axhline(floor, color=INK, lw=1.4, zorder=4)
            ax.annotate(f"beam + prior floor ({floor:.0f}$\\times$)",
                        (knees[0], floor), textcoords="offset points",
                        xytext=(2, -13), color=INK, fontsize=8.8)
        ax.axhline(1.0, color=BASELINE, lw=1.0, ls=(0, (4, 3)), zorder=2)
        ax.annotate("HI level", (knees[0], 1.0), textcoords="offset points",
                    xytext=(2, 4), color=MUTED, fontsize=8.4)
        # The fiducial model, and the frequency the raster actually scans at.
        ax.axvline(float(scan["fiducial_knee_mhz"]), color=MUTED, lw=1.0,
                   ls=(0, (1, 2)), zorder=1)
        ax.annotate("assumed", (float(scan["fiducial_knee_mhz"]), ax.get_ylim()[1]),
                    textcoords="offset points", xytext=(3, -11),
                    color=MUTED, fontsize=8.4, rotation=90, va="top")
        ax.set_xscale("log"); ax.set_yscale("log")
        # Margin either side: the measured points sit at the left edge of the
        # grid and would otherwise be clipped by the spine.
        ax.set_xlim(knees.min() * 0.7, knees.max() * 1.4)
        ax.set_xlabel("1/f knee frequency  [mHz]")
        ax.set_title(STRAT[s]["label"], color=STRAT[s]["color"], pad=6,
                     fontweight="semibold")
    axes[0, 0].set_ylabel(f"cleaned 1/f residual / HI   ({nm} modes removed)")

    if cutoff is not None:
        # w_c changes the SHAPE of the correlation, so unlike the knee axis
        # nothing factors out and every point here is an independent draw.
        wcs = np.asarray(cutoff["cutoffs_mhz"], float)
        wc_fid = float(cutoff["fiducial_wc_mhz"])
        for ax, s in zip(axes[1], ("drift", "raster")):
            _tidy(ax)
            y = [cutoff[f"{s}_ratio_wc{w}"] for w in wcs]
            ax.plot(wcs, y, color=DRIFT if s == "drift" else RASTER,
                    marker="o", ms=6, mfc=SURFACE, mew=1.8, zorder=3)
            if abl is not None:
                p_hi = abl[f"{s}_cr_p_hi"]
                floor = np.median(np.sqrt(abl[f"{s}_cr_floor_{nm}"] / p_hi))
                ax.axhline(floor, color=INK, lw=1.4, zorder=4)
                ax.annotate(f"floor ({floor:.0f}$\\times$)", (wcs[0], floor),
                            textcoords="offset points", xytext=(2, -13),
                            color=INK, fontsize=8.8)
            ax.axhline(1.0, color=BASELINE, lw=1.0, ls=(0, (4, 3)), zorder=2)
            ax.axvline(wc_fid, color=MUTED, lw=1.0, ls=(0, (1, 2)), zorder=1)
            ax.annotate("assumed", (wc_fid, ax.get_ylim()[1]),
                        textcoords="offset points", xytext=(3, -11),
                        color=MUTED, fontsize=8.4, rotation=90, va="top")
            ax.set_xscale("log"); ax.set_yscale("log")
            ax.set_xlim(wcs.min() * 0.7, wcs.max() * 1.4)
            ax.set_xlabel("1/f low-frequency cut-off $\\omega_c$  [mHz]")
        axes[1, 0].set_ylabel("cleaned 1/f residual / HI")
    # Lower right: the curves run bottom-left to top-right, so upper left is
    # where the "assumed" marker lives and lower right is the free corner.
    axes[0, 0].legend(loc="lower right", fontsize=8.6)

    sub = ("Top: knee = where 1/f gain power crosses white. Filled circles are "
           "measured, one per slope; lines are the exact "
           "$\\propto\\,$knee$^{\\alpha/2}$ scaling, since the pipeline is linear "
           "in the noise amplitude.")
    if cutoff is not None:
        sub += ("\nBottom: the low-frequency cut-off, every point an "
                "independent draw. Left of the assumed value the curve is flat "
                "while the gain rms grows 2.6$\\times$ — power added on "
                "timescales longer than a pass multiplies a rank-1 foreground, "
                "so the cleaning takes it out too.")
    _title(fig, "1/f would have to be orders of magnitude worse to matter", sub,
           x=0.045, y_title=1.03, y_sub=0.965, override=text)
    fig.tight_layout(rect=(0, 0, 1, 0.86 if nrow == 1 else 0.93))
    return _save(fig, "hi_1f_scan")


def fig_eigenvectors(rank, nmodes=4, nbars=6, text=None):
    """Figure 8 -- the spectral shapes, and how much each one carries.

    "Rank" is abstract until you look at the eigenvectors: these ARE the shapes
    each component is built from. The foreground's first mode is a clean power
    law; the floor's get progressively wigglier; the HI's are noise-like from
    the start, which is why it needs ~25 of them.

    The lower row is the point of the figure. On a LINEAR axis the foreground
    has one bar and then nothing -- mode 2 is 2e-8 of mode 1, which is a blank
    column, not a short one. That is what rank 1 looks like, and it is the
    thing the log-scale eigenspectrum in fig_rank cannot show.
    """
    from matplotlib.colors import LinearSegmentedColormap
    ramp = LinearSegmentedColormap.from_list(
        "blues", ["#9ec5f4", "#5598e7", "#2a78d6", "#184f95", "#0d366b"])
    freqs = rank["freqs"]
    panels = [("drift_V_fg", "drift_w_fg", INK,
               "GDSM foreground\n(rank 1)"),
              ("drift_V_floor", "drift_w_floor", DRIFT,
               "Beam + prior floor, drift\n(rank 6)"),
              ("drift_V_hi", "drift_w_hi", THIRD,
               "HI, map-made\n(rank 23)")]
    fig, axes = plt.subplots(2, 3, figsize=_figsize("hi_eigenvectors", (13.2, 7.0)),
                             sharex="row", layout="constrained")

    for col, (vkey, wkey, colour, title) in enumerate(panels):
        # --- top: the shapes -------------------------------------------
        ax = _tidy(axes[0, col])
        V = rank[vkey]
        for m in range(min(nmodes, V.shape[1])):
            v = V[:, m]
            # Sign of an eigenvector is arbitrary; fix it so the panels are
            # comparable rather than flipping at random between components.
            if v[np.argmax(np.abs(v))] < 0:
                v = -v
            ax.plot(freqs, v, color=ramp(m / max(nmodes - 1, 1)), lw=1.8,
                    label=f"mode {m + 1}")
        ax.axhline(0.0, color=BASELINE, lw=1.0)
        ax.set_xlabel("frequency  [MHz]")
        ax.set_title(title, color=INK, pad=6, fontsize=10)

        # --- bottom: how much each shape carries -----------------------
        ax = _tidy(axes[1, col])
        w = rank[wkey][:nbars]
        ax.bar(np.arange(1, len(w) + 1), w, color=colour, width=0.62,
               edgecolor=SURFACE, linewidth=1.4, zorder=3)
        ax.set_ylim(0, 1.12)
        ax.set_xticks(np.arange(1, len(w) + 1))
        ax.set_xlabel("eigenmode")
        ax.grid(axis="x", visible=False)
        # Spell out the invisible bars rather than leaving them a mystery.
        if w[1] < 1e-3:
            ax.annotate(f"modes 2$-${len(w)} total\n"
                        f"{w[1:].sum():.0e} of mode 1",
                        (2.6, 0.5), ha="left", va="center",
                        fontsize=8.6, color=INK2)
        else:
            for i, v in enumerate(w[1:], start=2):
                ax.annotate(f"{v:.2f}", (i, v), textcoords="offset points",
                            xytext=(0, 3), ha="center", fontsize=7.6,
                            color=INK2)
    axes[0, 0].set_ylabel("eigenvector amplitude  (arb.)")
    axes[0, 0].legend(loc="upper right", fontsize=8.4, ncols=2)
    axes[1, 0].set_ylabel("$\\lambda_i\\,/\\,\\lambda_1$   (linear)")

    _title(fig, "What the components are actually made of",
             "Top: leading eigenvectors of the channel-channel covariance. The "
             "foreground's mode 1 is a single smooth power law; its modes 2$-$4 look "
             "like structure but are floating-point noise.\n"
             "Bottom: how much each mode carries, on a LINEAR axis — the "
             "foreground's modes 2+ are a blank column, not a short one. That is "
             "what rank 1 looks like.",
           x=0.01, y_title=1.11, y_sub=1.02, override=text)
    return _save(fig, "hi_eigenvectors")


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
        made.append(fig_eigenvectors(rk))
    else:
        print("!! rank results missing -- run ska_hi_rank.py", flush=True)
    for m in made:
        print(f"wrote {os.path.relpath(m, _HERE)}", flush=True)


if __name__ == "__main__":
    main()
