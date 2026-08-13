"""Send-ready figures for the experiment-005 frequency sweep.

Emits three PNGs into `figures/`:

* `freqsweep_maps_sampled.png`   — maps, **beam-smoothed-truth prior** (the
  series' standard recipe; this is the version first circulated)
* `freqsweep_maps_flatprior.png` — maps, **flat prior**, so the panels show
  what the *data* recovers rather than what the prior supplies
* `freqsweep_prior_dependence.png` — the quantitative prior answer

Two things differ from the notebook's `figures/freqsweep_maps.png`, both aimed
at a reader who has not read `METHODS.md`:

1. **Each channel is shown on a grid where its beam is adequately sampled.**
   The notebook figure puts all five channels on one nside 128 grid so the
   residuals compare like for like, which is right for the sweep — but that
   grid undersamples the top two channels (3.5 and 2.9 px/FWHM) and inflates
   their residuals by 19% and 25%. Here 875 and 1050 MHz use nside 256.
2. **The numbers are printed on the panels.** Colour scales are necessarily
   per-row (the sky dims ~24x across the band, so a shared scale would black
   out the top rows), so the images alone cannot show the residual falling.

**Why the flat-prior version exists.** The standard recipe sets the Wiener
prior mean to the *beam-smoothed true sky*, and the beam narrows with
frequency — so the prior sharpens across the band on its own, with no data
involved. Measured, it supplies essentially all of the apparent improvement:
the beam-smoothed prior alone scores 0.511 -> 0.360 across the band while the
full reconstruction scores 0.525 -> 0.336, i.e. the data buys 0.97-1.16x. With
a flat prior the data's own contribution is visible and does grow with
frequency (1.30x -> 1.67x). Quote the flat-prior numbers, or the mode count,
when the claim is about what the survey measures.

Run: `.venv/bin/python examples/SKA/plot_freqsweep_maps.py` (needs the caches
warmed by `run_freq_sweep.sh`, plus the nside-256 runs for 875 and 1050).
"""

from __future__ import annotations

import os
import sys

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))

import ska_freq_sweep as S                                    # noqa: E402
from ska_common import gdsm_equatorial_sky_model              # noqa: E402

# Grid per channel: the shared nside 128 everywhere it is adequately sampled,
# nside 256 for the two channels where it is not (see pitfall 4 in METHODS.md).
GRID = {350: 128, 525: 128, 700: 128, 875: 256, 1050: 256}

FOV_X_DEG, FOV_Y_DEG = 22.0, 16.0   # roomy enough that no row's patch clips
FIGDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")

# Okabe-Ito. Two hues for the two prior choices; line style separates
# "prior alone" from "reconstruction" within each, so identity never rests on
# a pair of near-neighbour colours.
C_INFO, C_FLAT, C_GREY = "#0072B2", "#D55E00", "0.45"

PRIORS = {
    "smoothed": dict(
        label="beam-smoothed-truth prior",
        short="beam-smoothed truth",
        file="freqsweep_maps_sampled.png"),
    "flat": dict(
        label="flat prior (patch mean, no structure)",
        short="flat",
        file="freqsweep_maps_flatprior.png"),
}


def prior_mean(kind, truth_full, truth, pix, freq):
    if kind == "smoothed":
        return hp.smoothing(truth_full,
                            fwhm=np.radians(S.ska_beam_fwhm_deg(freq)))[pix]
    if kind == "flat":
        # Constant at the patch mean: carries the level but no structure, so
        # every spatial mode in the answer has to come from the data.
        return np.full_like(truth, float(np.mean(truth)))
    raise ValueError(kind)


def solve_channel(freq, nside, kinds=("smoothed", "flat")):
    """Full-noise, no-high-pass solve for each prior choice."""
    tod = S.simulate_channel(freq, nside, verbose=False)
    mm = S.build_operator(freq, nside, tod, verbose=False)
    pix = mm.pixel_indices
    truth_full = gdsm_equatorial_sky_model(freq=freq, nside=nside)
    truth = truth_full[pix]
    sky = float(np.std(truth))
    nv_floor = S.WHITE_VAR * (1e-3 * float(np.mean(truth)))**2
    noise = [S.WHITE_VAR * (np.asarray(o) @ truth)**2 + nv_floor
             for o in mm.Tsys_operators]
    tg = [np.asarray(t, np.float64) for t in tod["TOD_group"]]

    out = dict(pix=pix, truth=truth, nside=nside, sky=sky, est={}, prior={})
    for kind in kinds:
        mu = prior_mean(kind, truth_full, truth, pix, freq)
        out["est"][kind] = mm(
            TOD_group=tg, dtime=S.DT, Tsky_prior_mean=mu,
            Tsky_prior_inv_cov_diag=np.ones_like(truth) / max(sky, 1e-3)**2,
            noise_variance=noise, regularization=1e-12,
            return_full_cov=False)[0]
        out["prior"][kind] = mu

    A = np.vstack([np.asarray(o) for o in mm.Tsys_operators])
    Ninv = np.concatenate([1.0 / n for n in noise])
    ev = np.linalg.eigvalsh(A.T @ (Ninv[:, None] * A))[::-1]
    out["modes"] = int((ev > 0.01 * ev.max()).sum())
    return out


def project(nside, pix, vals, centre):
    """Gnomonic patch as a plain 2-D array, so the panels can be laid out and
    annotated normally instead of through gnomview's own subplot handling."""
    m = np.full(hp.nside2npix(nside), hp.UNSEEN)
    m[pix] = vals
    reso = 3.0 * (128 / nside)                       # arcmin, constant FOV
    img = hp.gnomview(m, rot=centre, reso=reso,
                      xsize=int(FOV_X_DEG * 60 / reso),
                      ysize=int(FOV_Y_DEG * 60 / reso),
                      return_projected_map=True, no_plot=True)
    # Keep UNSEEN masked. Dropping the mask would clip those pixels to the
    # bottom of the colour scale, so unobserved sky would render as though it
    # were real cold sky.
    return np.ma.masked_less(np.ma.filled(img, hp.UNSEEN), -1e29)


def maps_figure(res, chans, centre, kind):
    cmap_sky = plt.get_cmap("cividis").copy()
    cmap_res = plt.get_cmap("coolwarm").copy()
    for cm in (cmap_sky, cmap_res):
        cm.set_bad("0.90")

    fig, axes = plt.subplots(len(chans), 3, figsize=(14.5, 3.0 * len(chans)),
                             constrained_layout=True)
    for r, f in enumerate(chans):
        d = res[f]
        est = d["est"][kind]
        resid = est - d["truth"]
        lo, hi = np.percentile(d["truth"], [1, 99])
        rlim = float(np.percentile(np.abs(resid), 99))
        panels = [(d["truth"], cmap_sky, lo, hi),
                  (est, cmap_sky, lo, hi),
                  (resid, cmap_res, -rlim, rlim)]
        ims = []
        for c, (vals, cmap, vmin, vmax) in enumerate(panels):
            img = project(d["nside"], d["pix"], vals, centre)
            ax = axes[r, c]
            ims.append(ax.imshow(img, origin="lower", cmap=cmap, vmin=vmin,
                                 vmax=vmax, interpolation="nearest",
                                 extent=[-FOV_X_DEG / 2, FOV_X_DEG / 2,
                                         -FOV_Y_DEG / 2, FOV_Y_DEG / 2]))
            ax.set_xticks([]), ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_edgecolor("0.7")
            if r == 0:
                ax.set_title(["input sky", "reconstruction (same scale)",
                              "residual (rec − input)"][c], fontsize=12)

        # One colourbar per axes: attaching to a list of axes places the bar
        # between columns and misaligns it across rows.
        for im, ax in ((ims[1], axes[r, 1]), (ims[2], axes[r, 2])):
            cb = fig.colorbar(im, ax=ax, fraction=0.038, pad=0.015)
            cb.set_label("K", fontsize=9)
            cb.ax.tick_params(labelsize=8)

        box = dict(boxstyle="round,pad=0.32", fc="white", ec="0.6", alpha=0.88)
        axes[r, 2].annotate(
            f"residual rms {float(np.std(resid)):.3f} K\n"
            f"= {float(np.std(resid)) / d['sky']:.2f} × sky structure",
            xy=(0.03, 0.03), xycoords="axes fraction", fontsize=9.5,
            va="bottom", ha="left", bbox=box)
        axes[r, 0].annotate(
            f"{f} MHz\nFWHM {S.ska_beam_fwhm_deg(f):.2f}°\n"
            f"{d['modes']} sky modes",
            xy=(0.03, 0.97), xycoords="axes fraction", fontsize=9.5,
            va="top", ha="left", bbox=box)

    if kind == "flat":
        head = ("SKA-Mid-like drift scan across Band 1, reconstructed with a "
                "FLAT prior:\nwhat the data alone recovers, and it recovers "
                "more at higher frequency")
        prior_note = (
            "Wiener prior mean is a CONSTANT (the patch mean), so every "
            "spatial structure shown had to come from the data. Residual/sky "
            "falls 0.77 → 0.60 across the band;\nthe data's improvement over "
            "its own prior grows 1.30× → 1.67×. The companion figure uses the "
            "series' standard beam-smoothed-truth prior, which is sharper and "
            "flatters the maps.")
    else:
        head = ("SKA-Mid-like drift scan across Band 1: the beam narrows, the "
                "survey measures more sky modes,\nand the reconstruction "
                "improves relative to the sky it is measuring")
        prior_note = (
            "Wiener prior mean is the BEAM-SMOOTHED INPUT SKY, so this is not "
            "a blind recovery — and because the beam narrows with frequency "
            "the prior sharpens across the band by itself.\nThat prior alone, "
            "with no data, scores 0.51 → 0.36, so it supplies most of the "
            "trend below; see the flat-prior companion figure for what the "
            "data contributes.")

    fig.suptitle(head, fontsize=13.5)
    fig.supxlabel(
        "Off-plane field (b ≈ +50°), 3 h of drift on three sidereal nights, "
        "identical geometry and noise seeds at every channel — only the "
        "frequency changes.\n"
        "Colour scales are PER ROW: the sky itself dims ~24× across the band "
        "(ν$^{-2.7}$), so most of the fall in absolute residual is the fainter "
        "sky, not better map-making.\n"
        + prior_note + "\n"
        "Each channel is shown on a grid that samples its beam adequately "
        "(nside 128 up to 700 MHz, nside 256 above). The shrinking patch is "
        "the fixed 2° elevation ladder becoming wider than a beam.",
        fontsize=8.6, color="0.25")

    out = os.path.join(FIGDIR, PRIORS[kind]["file"])
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def prior_figure(res, chans):
    """The direct answer to 'how much of this is the prior?'."""
    def ratio(vals, f):
        return float(np.std(vals - res[f]["truth"])) / res[f]["sky"]

    p_info = [ratio(res[f]["prior"]["smoothed"], f) for f in chans]
    r_info = [ratio(res[f]["est"]["smoothed"], f) for f in chans]
    p_flat = [ratio(res[f]["prior"]["flat"], f) for f in chans]
    r_flat = [ratio(res[f]["est"]["flat"], f) for f in chans]

    fig, (axa, axb) = plt.subplots(1, 2, figsize=(11.5, 4.4))

    axa.plot(chans, p_info, "o--", color=C_INFO, lw=2, ms=7,
             label="beam-smoothed-truth prior, NO data")
    axa.plot(chans, r_info, "o-", color=C_INFO, lw=2.4, ms=8,
             label="reconstruction with that prior")
    axa.plot(chans, p_flat, "s--", color=C_FLAT, lw=2, ms=7,
             label="flat prior, NO data")
    axa.plot(chans, r_flat, "s-", color=C_FLAT, lw=2.4, ms=8,
             label="reconstruction with flat prior")
    axa.set_xlabel("Frequency [MHz]")
    axa.set_ylabel("residual rms ÷ sky structure rms")
    axa.set_title("With the standard prior, the prior alone already\n"
                  "explains the trend (blue dashed ≈ blue solid)", fontsize=10)
    axa.set_ylim(0, 1.12)
    axa.legend(frameon=False, fontsize=8.5, loc="lower left")
    axa.grid(alpha=0.3)

    gain_info = [p / r for p, r in zip(p_info, r_info)]
    gain_flat = [p / r for p, r in zip(p_flat, r_flat)]
    axb.plot(chans, gain_info, "o-", color=C_INFO, lw=2.4, ms=8,
             label="beam-smoothed-truth prior")
    axb.plot(chans, gain_flat, "s-", color=C_FLAT, lw=2.4, ms=8,
             label="flat prior")
    axb.axhline(1.0, color=C_GREY, ls=":", lw=1.5)
    axb.annotate("data adds nothing", (chans[-1], 1.0),
                 textcoords="offset points", xytext=(-6, -15), fontsize=8.5,
                 color=C_GREY, ha="right")
    axb.set_xlabel("Frequency [MHz]")
    axb.set_ylabel("prior-only residual ÷ reconstruction residual")
    axb.set_title("What the DATA buys over its own prior —\n"
                  "visible only once the prior stops carrying the answer",
                  fontsize=10)
    axb.legend(frameon=False, fontsize=8.5, loc="upper left")
    axb.grid(alpha=0.3)

    fig.suptitle("How much of the frequency trend is the prior?", fontsize=13)
    fig.tight_layout()
    out = os.path.join(FIGDIR, "freqsweep_prior_dependence.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")

    print(f"\n{'freq':>6s} | {'beam-smoothed-truth prior':>26s} | "
          f"{'flat prior':>26s}")
    print(f"{'':6s} | {'prior':>8s} {'recon':>8s} {'gain':>8s} | "
          f"{'prior':>8s} {'recon':>8s} {'gain':>8s}")
    for i, f in enumerate(chans):
        print(f"{f:6d} | {p_info[i]:8.3f} {r_info[i]:8.3f} "
              f"{gain_info[i]:7.2f}x | {p_flat[i]:8.3f} {r_flat[i]:8.3f} "
              f"{gain_flat[i]:7.2f}x")


def main():
    chans = S.CHANNELS_MHZ
    res = {f: solve_channel(f, GRID[f]) for f in chans}

    # One field of view for every row, centred on the widest patch, so the
    # shrinking Dec coverage is visible rather than normalised away.
    p0 = res[chans[0]]
    v = np.array(hp.pix2vec(p0["nside"], p0["pix"])).mean(axis=1)
    lon, lat = hp.vec2dir(v, lonlat=True)
    centre = (float(lon), float(lat))

    os.makedirs(FIGDIR, exist_ok=True)
    for kind in ("smoothed", "flat"):
        maps_figure(res, chans, centre, kind)
    prior_figure(res, chans)


if __name__ == "__main__":
    main()
