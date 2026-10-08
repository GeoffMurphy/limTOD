"""Why does the one-dish raster lose at 675-725 MHz?

At 350-400 MHz the raster's floor is 10.5x BELOW a mode-density-matched drift
ladder. At 675-725 MHz it is above both drifts: 27.1 against 16.5 (three-strip)
and 19.3 (``drift12m``). The drift-drift difference there is already explained
by ladder spacing -- 1.00 beams against 0.60. This asks whether the raster's
loss is also a sampling effect of its fixed geometry (el 38.19 deg, +/-7.5 deg
azimuth throw about 45/315, one half-sweep per 150 s) now that the beam is half
as wide: FWHM 3.98 deg at 350.78 MHz against 2.07 deg at 675.78 MHz.

Diagnostic only. Everything here is read from the stored pointings, the cached
operators and the result npz; nothing is simulated, solved or map-made.

**What is reconstructed, and how well.** The operators exist only at
350-400 MHz, nside 64 (``hi_cache/``), so each strategy's footprint at
675-725 MHz is rebuilt from its boresight track: stack a Gaussian of the band's
FWHM at every sample and keep the brightest ``npix`` pixels, with ``npix``
taken from the run that measured it. Checked against the two cached operators
at 350-400 MHz, that reproduces **94%** of the true pixel set. It is used only
for patch-brightness statistics (item 4); everything else runs on the interior
pixel set, which is stored exactly.

    /home/geoff/limTOD/.venv/bin/python ska_raster_sampling.py
"""
from __future__ import annotations

import os
import pickle
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

RESDIR = os.path.join(_HERE, "results")
SIGMA_PER_FWHM = 1.0 / 2.3548200450309493
ROW = {"drift": "3-strip drift", "drift12": "12-strip drift",
       "drift12m": "12-strip drift", "raster": "1-dish raster"}


class Band:
    def __init__(self, tag, label, freq, nside, experiment, strategies, npix):
        self.tag, self.label, self.freq, self.nside = tag, label, freq, nside
        self.experiment, self.strategies, self.npix = experiment, strategies, npix


# ``npix`` is channel 0's operator pixel count, from the runs that measured it:
# hi_modecheck_* for the drifts at both bands and the raster at 675-725, and
# the cached operators for 350-400.
BANDS = {
    "f350": Band("f350", "350-400 MHz", 350.78125, 64,
                 "hi_experiment_matched_f350_400_nc32_ns64.npz",
                 ("drift", "drift12", "raster"),
                 dict(drift=316, drift12=804, raster=684)),
    "f700": Band("f700", "675-725 MHz", 675.78125, 128,
                 "hi_experiment_matched700b_f675_725_nc32_ns128.npz",
                 ("drift", "drift12m", "raster"),
                 dict(drift=725, drift12m=1581, raster=1827)),
}


# ---------------------------------------------------------------------------
# Boresight geometry, from the stored pointings alone
# ---------------------------------------------------------------------------

def _location():
    from astropy.coordinates import EarthLocation
    from astropy import units as u
    from ska_common import SKA_LAT, SKA_LON, SKA_HGT
    return EarthLocation(lat=SKA_LAT * u.deg, lon=SKA_LON * u.deg,
                         height=SKA_HGT * u.m)


def tracks(strategy):
    """ICRS boresight track per pass: ``[(ra_deg, dec_deg, az_deg), ...]``."""
    import ska_hi_experiment as X
    from ska_common import azel_to_radec

    loc = _location()
    out = []
    for p in X.STRATEGIES[strategy]["pointings"]():
        az = np.asarray(p["azlist"], float)
        ra, dec = azel_to_radec(az, np.full(len(az), float(p["el"])),
                                np.asarray(p["tlist"], float),
                                X.STRATEGIES[strategy]["utc"], loc)
        out.append((np.asarray(ra), np.asarray(dec), az))
    return out


def _vecs(ra, dec):
    import healpy as hp
    return hp.ang2vec(ra, dec, lonlat=True)


def sweeps(ra, dec, az):
    """Split one raster pass into half-sweeps by the sign of the azimuth step."""
    d = np.sign(np.diff(az))
    d[d == 0] = 1.0
    edges = np.flatnonzero(np.diff(d) != 0) + 1
    idx = np.split(np.arange(len(az)), edges + 1)
    return [i for i in idx if len(i) > 5]


def sweep_spacing(strategy, fwhm):
    """Perpendicular separation between consecutive SAME-direction sweeps.

    For every point of sweep j, the nearest point of sweep j+2; the median over
    the sweep is that pair's separation. Reported over all pairs of all passes.
    The tracks are geometry, so the separation in degrees is identical at both
    bands -- only the beam it is measured against changes.
    """
    seps = []
    for ra, dec, az in tracks(strategy):
        sw = sweeps(ra, dec, az)
        V = _vecs(ra, dec)
        for j in range(len(sw) - 2):
            a, b = V[sw[j]], V[sw[j + 2]]
            cos = np.clip(a @ b.T, -1.0, 1.0)
            seps.append(float(np.median(np.degrees(np.arccos(cos)).min(axis=1))))
    s = np.asarray(seps)
    return dict(n=len(s), median_deg=float(np.median(s)),
                max_deg=float(s.max()), min_deg=float(s.min()),
                median_beams=float(np.median(s) / fwhm),
                max_beams=float(s.max() / fwhm))


def beam_exposure(strategy, pix, nside, fwhm, per_pass=True):
    """``sum_t exp(-theta^2/2 sigma^2)`` at each pixel, per pass.

    With the Gaussian normalised to unit peak this sum IS samples per beam
    solid angle: ``int w dOmega = 2 pi sigma^2``, the Gaussian beam area.
    """
    import healpy as hp
    sig = fwhm * SIGMA_PER_FWHM
    pv = np.asarray(hp.pix2vec(nside, np.asarray(pix))).T       # (npix, 3)
    out = []
    for ra, dec, _az in tracks(strategy):
        V = _vecs(ra, dec)                                       # (nsamp, 3)
        acc = np.zeros(len(pv))
        peak = np.zeros(len(pv))
        for i0 in range(0, len(V), 2000):                        # chunked
            cos = np.clip(pv @ V[i0:i0 + 2000].T, -1.0, 1.0)
            w = np.exp(-0.5 * (np.degrees(np.arccos(cos)) / sig) ** 2)
            acc += w.sum(axis=1)
            peak = np.maximum(peak, w.max(axis=1))
        out.append((acc, peak))
    return out if per_pass else (sum(a for a, _ in out),
                                 np.max([p for _, p in out], axis=0))


def nearest_track_deg(strategy, pix, nside):
    """Angular distance from each pixel to the nearest boresight sample."""
    import healpy as hp
    pv = np.asarray(hp.pix2vec(nside, np.asarray(pix))).T
    best = np.full(len(pv), np.inf)
    for ra, dec, _az in tracks(strategy):
        V = _vecs(ra, dec)
        for i0 in range(0, len(V), 2000):
            cos = np.clip(pv @ V[i0:i0 + 2000].T, -1.0, 1.0)
            best = np.minimum(best, np.degrees(np.arccos(cos)).min(axis=1))
    return best


def patch_pixels(strategy, band):
    """The strategy's footprint, as the ``npix`` most-exposed pixels."""
    import healpy as hp
    sig = band.freq and ska_fwhm(band) * SIGMA_PER_FWHM
    m = np.zeros(hp.nside2npix(band.nside))
    rad = np.radians(3.5 * ska_fwhm(band))
    for ra, dec, _az in tracks(strategy):
        for v in _vecs(ra, dec):
            ip = hp.query_disc(band.nside, v, rad)
            cos = np.clip(np.asarray(hp.pix2vec(band.nside, ip)).T @ v, -1, 1)
            m[ip] += np.exp(-0.5 * (np.degrees(np.arccos(cos)) / sig) ** 2)
    return np.argsort(m)[::-1][:band.npix[strategy]]


def ska_fwhm(band):
    from ska_common import ska_beam_fwhm_deg
    return float(ska_beam_fwhm_deg(band.freq))


def interior_pixels(band):
    d = np.load(os.path.join(RESDIR, band.experiment), allow_pickle=True)
    return np.asarray(d["common"])[d["interior"]], d


def position_angles(strategy):
    """Median track position angle (deg E of N) per pass, and per sweep
    direction for the raster. The drift's is 90 deg at every azimuth."""
    out = []
    for ra, dec, az in tracks(strategy):
        dra = np.diff(np.unwrap(np.radians(ra))) * np.cos(np.radians(dec[:-1]))
        pa = np.degrees(np.arctan2(dra, np.diff(np.radians(dec)))) % 180.0
        if strategy.startswith("raster"):
            sw = sweeps(ra, dec, az)
            f = np.concatenate([sw[j] for j in range(0, len(sw), 2)])
            b = np.concatenate([sw[j] for j in range(1, len(sw), 2)])
            f, b = f[f < len(pa)], b[b < len(pa)]
            out.append((float(np.median(pa[f])), float(np.median(pa[b]))))
        else:
            out.append((float(np.median(pa)), float(np.median(pa))))
    return out


def collect(band):
    """Every number the report needs, for one band."""
    import healpy as hp
    from ska_common import gdsm_equatorial_sky_model

    fwhm = ska_fwhm(band)
    inter, d = interior_pixels(band)
    sky = gdsm_equatorial_sky_model(freq=band.freq, nside=band.nside)
    out = dict(fwhm=fwhm, n_interior=len(inter),
               area_interior=len(inter) * hp.nside2pixarea(band.nside,
                                                           degrees=True),
               sky_interior_mean=float(sky[inter].mean()),
               sky_interior_rms=float(sky[inter].std()), strategies={})
    for s in band.strategies:
        pa = position_angles(s)
        # two DIRECTIONS, not two passes: a pass counts separately only if its
        # track position angle differs from another's by more than 20 deg.
        pp = beam_exposure(s, inter, band.nside, fwhm)
        seen = np.asarray([p > 0.5 for _, p in pp])
        ndir = np.zeros(len(inter), int)
        for i in range(len(pp)):
            fresh = seen[i] & ~np.any(
                [seen[j] & (abs(pa[i][0] - pa[j][0]) < 20.0)
                 for j in range(i)] or [np.zeros(len(inter), bool)], axis=0)
            ndir += fresh
        rec = dict(
            n_pass=len(pp),
            n_samples=int(sum(len(t[0]) for t in tracks(s))),
            pa=[p[0] for p in pa],
            samp_per_beam=float(sum(a for a, _ in pp).mean()),
            two_dir=float(np.mean(ndir >= 2)),
            one_dir=float(np.mean(ndir == 1)),
            no_cover=float(np.mean(ndir == 0)),
            max_gap_beams=float(nearest_track_deg(s, inter,
                                                  band.nside).max() / fwhm),
            npix=band.npix[s])
        if s.startswith("raster"):
            rec["sweep"] = sweep_spacing(s, fwhm)
        patch = patch_pixels(s, band)
        rec["sky_patch_mean"] = float(sky[patch].mean())
        rec["sky_patch_rms"] = float(sky[patch].std())
        rec["sky_patch_over_interior"] = (rec["sky_patch_mean"]
                                          / out["sky_interior_mean"])
        out["strategies"][s] = rec
    return out


# Measured elsewhere in the series, quoted so the conclusion is self-contained.
# Median sqrt(p_clean_auto_4), cr_, interior -- the INSTRUMENT half of
# residual/Hi, with the cosmological half (the Hi's own dimming) divided out.
# Source: HANDOFF.md, "Both strategies improve monotonically with frequency".
FG_RESIDUAL = {"350-400": dict(drift=6.54e-2, raster=6.45e-3),
               "500-550": dict(drift=1.48e-2, raster=5.90e-3),
               "675-725": dict(drift=6.25e-4, raster=2.60e-3)}
# Floor/Hi amplitude at 4 modes, amplitude-only prior (nm_kpar_results.md).
FLOOR = {"f350": dict(drift=399.6, drift12=233.7, raster=22.3),
         "f700": dict(drift=16.5, drift12m=19.3, raster=27.1)}
# tr(WA) and patch pixel count, from hi_kmodes_* / hi_modecheck_*. The mean of
# diag(WA) over a patch IS tr(WA)/npix, so this is patch-average prior
# dominance without touching an operator.
TRWA = {"f350": dict(drift=38.78, drift12=145.35, raster=118.52),
        "f700": dict(drift=63.70, drift12m=240.65, raster=266.92)}
# Measured from the cached 350-400 operators (the only ones that exist here).
DIAG_INTERIOR_F350 = dict(drift=0.2038, raster=0.3283)


def write_report(path=None):
    path = path or os.path.join(_HERE, "raster_sampling.md")
    res = {t: collect(BANDS[t]) for t in ("f350", "f700")}
    L = []
    w = L.append
    w("# Why does the one-dish raster lose at 675-725 MHz?")
    w("")
    w("Generated by `ska_raster_sampling.py`, 2026-10-08. Diagnostic only: "
      "stored pointings, cached operators and result")
    w("npz. Nothing simulated, solved or map-made.")
    w("")
    w("## Answer")
    w("")
    w("**It does not.** The raster's post-clean foreground residual -- the "
      "instrument half of floor/Hi, with the Hi's own")
    w("dimming divided out -- *improves* by 2.5x from 350-400 to "
      "675-725 MHz. The three-strip drift's improves by **105x**")
    w("over the same interval. The inversion is the drift catching up and "
      "overtaking, not the raster degrading, and all")
    w("four sampling candidates in the brief are ruled out below.")
    w("")
    w("| | 350-400 | 500-550 | 675-725 | 350 -> 725 |")
    w("|---|---:|---:|---:|---:|")
    for s in ("drift", "raster"):
        v = [FG_RESIDUAL[b][s] for b in ("350-400", "500-550", "675-725")]
        w(f"| {ROW[s]}, foreground residual | {v[0]:.2e} | {v[1]:.2e} | "
          f"{v[2]:.2e} | **{v[0] / v[2]:,.0f}x better** |"
          if v[0] / v[2] >= 10 else
          f"| {ROW[s]}, foreground residual | {v[0]:.2e} | {v[1]:.2e} | "
          f"{v[2]:.2e} | **{v[0] / v[2]:.1f}x better** |")
    w("")
    w("That table is not new -- it is in HANDOFF under \"Both strategies "
      "improve monotonically with frequency\" -- but it is")
    w("the quantity the question should be asked of. `floor/Hi` is a quotient "
      "of an instrument number and a cosmology")
    w("number, and the Hi dims 4.5x over the same interval on its own.")
    w("")

    # ---- the table --------------------------------------------------------
    w("## Per band")
    w("")
    w("Measured on the interior of the common patch, the pixel set the floor "
      "is measured on. `samples/beam` is")
    w("$\\sum_t \\exp(-\\theta_t^2/2\\sigma^2)$ at each pixel, which is "
      "literally samples per Gaussian beam solid angle")
    w("($\\int w\\,d\\Omega = 2\\pi\\sigma^2$), averaged over interior "
      "pixels. `max gap` is the largest distance from any interior pixel")
    w("to the nearest boresight sample.")
    w("")
    for tag in ("f350", "f700"):
        r = res[tag]
        b = BANDS[tag]
        w(f"**{b.label}** -- FWHM {r['fwhm']:.3f} deg, interior "
          f"{r['n_interior']} px = {r['area_interior']:.1f} deg$^2$, "
          f"mean sky {r['sky_interior_mean']:.2f} K")
        w("")
        w("| strategy | floor/Hi | samples | patch px | samples/beam | "
          "track PA | 2-direction | max gap | patch sky / interior | "
          "tr($WA$)/px |")
        w("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
        for s in b.strategies:
            v = r["strategies"][s]
            w(f"| {ROW[s]} | {FLOOR[tag][s]:.1f} | {v['n_samples']:,} | "
              f"{v['npix']:,} | {v['samp_per_beam']:.0f} | "
              f"{v['pa'][0]:.0f} deg | {100 * v['two_dir']:.0f}% | "
              f"{v['max_gap_beams']:.2f} beams | "
              f"{v['sky_patch_over_interior']:.3f} | "
              f"{TRWA[tag][s] / v['npix']:.4f} |")
        w("")

    # ---- item 1 -----------------------------------------------------------
    rs = res["f350"]["strategies"]["raster"]["sweep"]
    rs7 = res["f700"]["strategies"]["raster"]["sweep"]
    w("## 1. Sweep spacing: the sweeps are not gappy, by a wide margin")
    w("")
    w(f"**\"6 arcmin/s\" is an AZIMUTH rate, not an on-sky rate.** "
      f"`constant_elevation_scan`'s `sweep_s = 150` s is one")
    w("half-sweep across the full 15 deg throw, so the azimuth rate is "
      "$2 \\times 7.5/150 = 0.1$ deg/s = 6'/s, and the on-sky")
    w("rate is that times $\\cos(\\mathrm{el}) = \\cos 38.19^\\circ$, i.e. "
      "4.7'/s. The triangle period is 300 s.")
    w("")
    w(f"Perpendicular separation between consecutive same-direction sweeps, "
      f"measured over {rs['n']} pairs across both passes:")
    w("")
    w("| | degrees | beams at 350-400 | beams at 675-725 |")
    w("|---|---:|---:|---:|")
    w(f"| median | {rs['median_deg']:.3f} | {rs['median_beams']:.3f} | "
      f"{rs7['median_beams']:.3f} |")
    w(f"| maximum | {rs['max_deg']:.3f} | {rs['max_beams']:.3f} | "
      f"{rs7['max_beams']:.3f} |")
    w("")
    w("The separation is geometry, so it is the same in degrees at both "
      "bands; only the beam it is measured against changes.")
    w("**0.37 beams at 675-725**, against the 0.5-0.6 the brief's envelope "
      "expected -- tighter, not gappier. The direct")
    w("coverage measure agrees: the largest distance from any interior pixel "
      "to the nearest boresight sample is")
    w(f"**{res['f700']['strategies']['raster']['max_gap_beams']:.2f} beams** "
      f"at 675-725 and "
      f"{res['f350']['strategies']['raster']['max_gap_beams']:.2f} at "
      "350-400. There are no holes. **Candidate (i) is ruled out.**")
    w("")

    # ---- item 2 -----------------------------------------------------------
    w("## 2. Two-direction coverage: 100%, in both bands")
    w("")
    w("A pixel counts as seen by a direction if some sample of a pass with a "
      "distinct track position angle (within 20 deg)")
    w("puts it above half the beam peak. The raster's two passes sit at "
      "PA 126 deg (rising, az 45) and 54 deg (setting,")
    w("az 315), 72 deg apart. **Every interior pixel is seen by both, in "
      "both bands.** For contrast, all three drift strips")
    w("sit at PA 89.9 deg with a scatter of 0.01 deg, so the drift is 100% "
      "single-direction by construction -- which is")
    w("Section \"Cross-linking and drift scans\" of the draft, measured. "
      "**Candidate (ii) is ruled out.**")
    w("")

    # ---- item 3 -----------------------------------------------------------
    d3, d7 = res["f350"]["strategies"], res["f700"]["strategies"]
    w("## 3. Sampling density: the right sign, the wrong size")
    w("")
    w("The raster is the sparsest strategy in both bands, and it does get "
      "relatively sparser:")
    w("")
    w("| | 350-400 | 675-725 |")
    w("|---|---:|---:|")
    w(f"| raster samples/beam | {d3['raster']['samp_per_beam']:.0f} | "
      f"{d7['raster']['samp_per_beam']:.0f} |")
    w(f"| 3-strip drift samples/beam | {d3['drift']['samp_per_beam']:.0f} | "
      f"{d7['drift']['samp_per_beam']:.0f} |")
    w(f"| raster / drift | {d3['raster']['samp_per_beam'] / d3['drift']['samp_per_beam']:.2f} | "
      f"{d7['raster']['samp_per_beam'] / d7['drift']['samp_per_beam']:.2f} |")
    w("")
    w("The cause is geometric: the raster's footprint is set by a throw fixed "
      "in degrees, so it barely shrinks as the beam")
    w("narrows, while the drift's cross-scan extent is set by the beam. From "
      "350-400 to 675-725 the raster's patch grows")
    w("from 2.2x to 2.5x the drift's, so the same 5400 samples are spread "
      "thinner.")
    w("")
    w("**But the size is wrong by an order of magnitude.** The raster/drift "
      "sampling ratio moves 0.68 -> 0.56, a factor")
    w("1.2, while floor/Hi moves from 0.056x the drift's to 1.64x -- a factor "
      "29. And at 350-400 the raster is *already*")
    w("the sparsest strategy and still wins by 17.9x. Sampling density "
      "plainly does not set the floor.")
    w("")
    w("**Prior dominance runs the wrong way too.** The mean of "
      "$\\mathrm{diag}(WA)$ over a patch is exactly "
      "tr$(WA)/N_{\\rm pix}$,")
    w("so it needs no operator:")
    w("")
    w("| | 350-400 | 675-725 |")
    w("|---|---:|---:|")
    for s3, s7 in (("drift", "drift"), ("raster", "raster")):
        w(f"| {ROW[s3]} | {TRWA['f350'][s3] / d3[s3]['npix']:.4f} | "
          f"{TRWA['f700'][s7] / d7[s7]['npix']:.4f} |")
    w("")
    w("The raster is **less** prior-dominated than the drift in both bands, "
      "not more. Where the operators exist (350-400,")
    w("nside 64) the interior can be done properly, and it says the same "
      "thing more strongly: mean diag$(WA)$ over the 119")
    w(f"interior pixels is **{DIAG_INTERIOR_F350['raster']:.3f}** for the "
      f"raster against **{DIAG_INTERIOR_F350['drift']:.3f}** for the drift, "
      "both about 1.7x their own")
    w("patch averages (the interior cut keeps the best-measured core -- the "
      "same effect as in `forecast_check.md`).")
    w("**Candidate (iii) is consistent in direction but cannot carry the "
      "effect, and its prior-dominance half is ruled out.**")
    w("")

    # ---- item 4 -----------------------------------------------------------
    w("## 4. Footprint brightness: a 5% effect against a 1.6x difference")
    w("")
    w("| | 350-400 | 675-725 |")
    w("|---|---:|---:|")
    for s3, s7 in (("drift", "drift"), ("drift12", "drift12m"),
                   ("raster", "raster")):
        w(f"| {ROW[s3]} patch sky / interior | "
          f"{d3[s3]['sky_patch_over_interior']:.3f} | "
          f"{d7[s7]['sky_patch_over_interior']:.3f} |")
    w("")
    w("The raster does reach brighter sky, and slightly more so at the higher "
      "band (1.048 -> 1.068), which is the right")
    w("direction. It is still only a 7% excess, and the twelve-strip drift's "
      "footprint is 5.9% brighter at 675-725 while")
    w("beating the raster by 1.4x. A 5-7% brightness difference cannot "
      "produce a 1.6x floor difference. **Candidate (iv) is")
    w("ruled out as the cause**, though it is a real and correctly-signed "
      "sub-effect.")
    w("")

    # ---- item 5 -----------------------------------------------------------
    w("## 5. Per-pixel floor map")
    w("")
    w("No per-pixel floor map is stored, but one can be rebuilt without "
      "map-making: PCA the stored `data_cube` at common")
    w("resolution and take the per-pixel rms over channels. The chain was "
      "validated first -- it reproduces the stored")
    w("`cr_p_clean_auto_4` to 1 part in $10^6$ in every bin.")
    w("")
    w("| band | strategy | median floor/Hi per px | 10-90% | "
      "corr(log samples/beam, log floor/Hi) |")
    w("|---|---|---:|---:|---:|")
    w("| 350-400 | 3-strip drift | 437 | 196-789 | -0.38 |")
    w("| 350-400 | 12-strip drift | 250 | 96-633 | -0.16 |")
    w("| 350-400 | 1-dish raster | 66 | 26-155 | -0.24 |")
    w("| 675-725 | 3-strip drift | 10.7 | 3.7-29.4 | -0.16 |")
    w("| 675-725 | 12-strip drift | 15.6 | 7.8-38.9 | -0.10 |")
    w("| 675-725 | 1-dish raster | 21.1 | 10.0-44.4 | **+0.62** |")
    w("")
    w("These are real-space rms ratios over all $k_\\parallel$, so they are "
      "not the headline numbers and should not be")
    w("quoted as such; what matters is that they reproduce the ordering and "
      "the inversion (66 < 250 < 437 at 350-400;")
    w("10.7 < 15.6 < 21.1 at 675-725), which makes them a valid relative "
      "probe.")
    w("")
    w("The raster's correlation with exposure **flips sign** between bands, "
      "-0.24 to +0.62. I would not read that as a")
    w("cause. Across the raster's interior at 675-725 the samples/beam spans "
      "only 156 to 166 -- a 6% range -- while the")
    w("per-pixel floor spans 13.4 to 34.6, a factor 2.6. A 6% lever cannot "
      "move a 2.6x spread, so the correlation is with")
    w("something that varies across the patch alongside exposure, not with "
      "exposure. What that something is, these")
    w("numbers do not isolate.")
    w("")

    # ---- verdict ----------------------------------------------------------
    w("## Which candidate survives")
    w("")
    w("| candidate | verdict |")
    w("|---|---|")
    w("| (i) gappy sweeps | **ruled out.** 0.37-beam sweep spacing, "
      "0.10-beam largest gap at 675-725. |")
    w("| (ii) single-pass interior pixels | **ruled out.** 100% "
      "two-direction coverage in both bands. |")
    w("| (iii) sampling density / prior dominance | **right sign, far too "
      "small.** 1.2x against a 29x swing; and the raster is the LESS "
      "prior-dominated strategy in both bands. |")
    w("| (iv) brighter footprint | **ruled out as the cause.** 7% excess, "
      "and the strategy that beats it carries 6%. |")
    w("")
    w("None of the four explains it, because the premise is wrong: the "
      "raster's sampling of the sky is essentially")
    w("unchanged between the bands -- same tracks, same 100% two-direction "
      "coverage, no gaps at either beam width -- and its")
    w("instrumental residual improves 2.5x. **What changed is the drift.** "
      "Its fixed 2 deg ladder step is 0.50 beams at")
    w("350.78 MHz and 1.00 beams at the 675-725 band centre: from heavily "
      "redundant to contiguous and non-redundant. That")
    w("is the same mechanism the draft already uses for drift vs "
      "`drift12m` at this band (1.00 against 0.60 beams) and for")
    w("the drift's rise at 850-900 (1.25 beams, gapped). The drift-raster "
      "inversion near 600 MHz is the third instance of")
    w("it, not a separate effect.")
    w("")
    w("This is a consistency argument from three points, not a measurement "
      "of an optimum. Spacing and mode density are not")
    w("varied independently anywhere in the series, so \"near one beam\" "
      "remains the hypothesis the draft already states.")
    w("")

    # ---- limits -----------------------------------------------------------
    w("## What could not be done")
    w("")
    w("- **Interior diag$(WA)$ at 675-725.** Operators are cached only at "
      "350-400, nside 64 (`hi_cache/`). Building one at")
    w("  nside 128 is operator construction, which the brief excludes, so "
      "the patch average (tr$(WA)/N_{\\rm pix}$, exact) is")
    w("  used there instead. It is a ~30-60 min local job for one channel "
      "per strategy if the interior value is wanted.")
    w("- **Footprints at 675-725** are reconstructed from the boresight "
      "track by stacking a Gaussian of the band FWHM and")
    w("  keeping the brightest $N_{\\rm pix}$ pixels, with $N_{\\rm pix}$ "
      "from the run that measured it. Checked against the two")
    w("  cached 350-400 operators, that recovers **94%** of the true pixel "
      "set. Used only for the brightness column.")
    w("")
    w("## Reproducing")
    w("")
    w("```")
    w("/home/geoff/limTOD/.venv/bin/python ska_raster_sampling.py")
    w("```")
    w("")
    w("The item-5 table needs `pyccl` for `ska_hi_analysis`, so it was run "
      "under `/home/geoff/gibbs_venv_312/bin/python`")
    w("and is transcribed here rather than recomputed.")

    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")
    print(f"wrote {path}", flush=True)
    return path


if __name__ == "__main__":
    write_report()
