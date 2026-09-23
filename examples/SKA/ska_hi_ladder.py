"""Tier 1: what does the drift's Dec ladder buy? (PLAN item 6, drift half.)

**A parked drift cannot cross-link, at any azimuth** (HANDOFF.md: the boresight
traces a constant-Dec circle swept in +RA, so every track has position angle 90
deg and azimuth only chooses *which* Dec). So PLAN item 6 -- "more crossing
angles" -- applies to the raster only. The drift's *one* geometry lever is the
Dec ladder: how many strips, and how far apart.

That is also the multi-dish question. ``drift_scan_night`` offsets night n by a
full sidereal day, so every strip already covers the same LST window --
**the published 3-night drift is arithmetically identical to three dishes
parked at three elevations for one hour.** Multi-dish buys the factor of 3 in
wall-clock, not in modes. Dishes at different elevations are the one exception
to "co-pointed dishes are a dead lever" (PLAN item 4), and this is what that
exception is worth.

Tier 1 because it needs the map-making OPERATOR but not the TOD: the mode count
is the eigenspectrum of ``A^T N^-1 A``, and the only thing the TOD is needed for
is the residual. That drops the cost from ~660 s per pass per channel to ~125 s,
and one channel suffices because the mode count is what is being measured.
Budget ~65 min for the whole sweep on 12 cores.

Two one-dimensional sweeps, crossing at the published geometry:

  A  three strips, spacing 0.25 - 1.5 beams. The published 52/50/48 deg ladder
     is 2 deg steps = 0.50 beams at 350 MHz, so it appears as one point here.
  B  spacing fixed at 0.50 beams, 1 to 8 strips. N = 1 is the honest baseline:
     it is what a single parked dish measures, and the ladder's whole value is
     the gap between it and N = 3.

**The quantity that matters is modes per unit sky, not modes.** Adding strips
adds area as well as modes, and the design figure normalises. The fork this
sweep is here to settle:

  * spacing >= 1 beam -> strips are disjoint, you buy area only, modes/deg^2
    flat. HANDOFF already records the symptom at the top of Band 1, where the
    fixed 2 deg ladder becomes a whole beam and "the three strips separate into
    ribbons with unobserved gaps".
  * spacing < 1 beam -> strips overlap and densify the *cross-scan* sampling,
    which is the direction a single drift strip cannot sample at all. That
    should raise modes per unit sky.

Caveats, both real:

* A mode count is not a residual. The prior-variance sweep
  (``ska_hi_priorvar.py``, HANDOFF 2026-09-21) showed modes and residual can
  move in opposite directions, so a ladder that adds modes per unit sky still
  needs a Tier 2 run before anyone claims it moves the HI residual.
* Spillover and ground pickup versus elevation are not modelled anywhere in
  this series (PLAN, deferred). A wide ladder reaches lower elevations and
  therefore looks free here when it is not.

Writes ``results/hi_ladder_f350_ns64.npz``.

    /home/geoff/gibbs_venv_312/bin/python ska_hi_ladder.py
"""

from __future__ import annotations

import argparse
import os
import pickle
import sys
import time

import healpy as hp
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import ska_hi_experiment as X
from ska_common import (SKA_LAT, SKA_LON, SKA_HGT, ska_beam_func,
                        ska_beam_fwhm_deg, gdsm_equatorial_sky_model,
                        drift_scan_night)

FREQ_MHZ = 350.78125            # channel 0, the frequency hi_modes.py quotes
EL_CENTRE = 50.0                # the published ladder is 52/50/48
SPACINGS_BEAMS = (0.25, 0.5, 0.75, 1.0, 1.5)
COUNTS = (1, 2, 3, 5, 8)
REF_SPACING, REF_COUNT = 0.5, 3   # the published geometry, in these units


def elevations(n_strips, spacing_beams, fwhm_deg):
    """Symmetric ladder about EL_CENTRE. n=3, 0.5 beams -> 48/50/52."""
    step = spacing_beams * fwhm_deg
    return EL_CENTRE + (np.arange(n_strips) - (n_strips - 1) / 2.0) * step


def build_operator(els):
    """One drift pass per elevation, each on its own sidereal day.

    Night n is a whole sidereal day later, so every strip sees the same LST
    window -- which is exactly what N dishes observing simultaneously would
    give. No TOD is simulated; only the LSTs the operator needs.
    """
    from limTOD import HPW_mapmaking
    from limTOD.simulator import generate_LSTs_deg

    lst_g, az_g, el_g = [], [], []
    for night, el in enumerate(els):
        tlist, azlist = drift_scan_night(X.DRIFT_NIGHT_S, dt=X.DT, night=night,
                                         azimuth_deg=X.DRIFT_AZ)
        lst_g.append(np.asarray(generate_LSTs_deg(
            SKA_LAT, SKA_LON, SKA_HGT, tlist, X.DRIFT_UTC), float))
        az_g.append(np.asarray(azlist, float))
        el_g.append(np.full(len(tlist), float(el)))
    return HPW_mapmaking(
        beam_map=ska_beam_func(freq=FREQ_MHZ, nside=X.NSIDE),
        LST_deg_list_group=lst_g, lat_deg=SKA_LAT,
        azimuth_deg_list_group=az_g, elevation_deg_list_group=el_g,
        threshold=X.THRESHOLD, nside_target=X.NSIDE,
        beam_truncate_frac_thres=X.TRUNCATE_FRAC_THRES)


def operator_list(mm):
    """``mm.Tsys_operators`` as a list, whatever the group count.

    ``HPW_mapmaking`` stores a *list* of per-TOD operators when ``num_tods > 1``
    but a single stacked ndarray when ``num_tods == 1`` (HPW_filter.py, and its
    own comment says so). Iterating it blindly then yields TOD *rows* instead of
    operators. Nothing else in this series has tripped over it because every
    other geometry here has at least two scans; the N = 1 rung of this ladder is
    the first single-group operator the series has ever built.
    """
    ops = mm.Tsys_operators
    if isinstance(ops, np.ndarray) and ops.ndim == 2:
        return [ops]
    return list(ops)


def mode_count(mm):
    """tr(WA) and friends, the same recipe as ska_hi_modes.py."""
    pix = np.asarray(mm.pixel_indices)
    truth = gdsm_equatorial_sky_model(freq=FREQ_MHZ, nside=X.NSIDE)[pix]
    nv_floor = X.WHITE_VAR * (1e-3 * float(np.mean(truth))) ** 2
    M = None
    for o in operator_list(mm):
        A = np.asarray(o, float)
        var = X.WHITE_VAR * (A @ truth) ** 2 + nv_floor
        term = A.T @ (A / var[:, None])
        M = term if M is None else M + term
    w = np.maximum(np.linalg.eigvalsh(M)[::-1], 0.0)
    s_inv = 1.0 / max(float(np.std(truth)), 1e-3) ** 2
    area = len(pix) * hp.nside2pixarea(X.NSIDE, degrees=True)
    return dict(n_eff=float(np.sum(w / (w + s_inv))),
                n_1pct=int((w / w[0] > 0.01).sum()),
                npix=len(pix), area=area)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--control-only", action="store_true",
                    help="just re-count the cached published operator")
    ap.add_argument("--counts", type=int, nargs="*", default=list(COUNTS),
                    help="strip counts to sweep at the reference spacing")
    ap.add_argument("--spacings", type=float, nargs="*",
                    default=list(SPACINGS_BEAMS),
                    help="spacings (in beams) to sweep at the reference count")
    ap.add_argument("--out", default=None,
                    help="output npz basename; defaults to hi_ladder_f350")
    args = ap.parse_args()

    fwhm = ska_beam_fwhm_deg(FREQ_MHZ)
    print(f"{FREQ_MHZ:.2f} MHz, FWHM {fwhm:.3f} deg, nside {X.NSIDE}, "
          f"ladder centred on {EL_CENTRE:.0f} deg", flush=True)

    # --- control: the cached published operator, free -----------------------
    cache = os.path.join(_HERE, "hi_cache",
                         f"op_drift_f{FREQ_MHZ:07.3f}_ns{X.NSIDE}.pkl")
    with open(cache, "rb") as f:
        r = mode_count(pickle.load(f))
    print(f"CONTROL published 52/50/48 (cached op): tr(WA) {r['n_eff']:.2f}, "
          f"{r['area']:.1f} deg^2, {100 * r['n_eff'] / r['area']:.2f}/100deg^2, "
          f"n>1% {r['n_1pct']}   [hi_modes.py: 38.78, 265.2, 14.62, 17]",
          flush=True)
    if args.control_only:
        return

    # --- the two sweeps, sharing the point where they cross -----------------
    jobs = {(REF_COUNT, s) for s in args.spacings}
    jobs |= {(n, REF_SPACING) for n in args.counts}
    jobs = sorted(jobs, key=lambda j: (j[0], j[1]))
    npass = sum(n for n, _ in jobs)
    print(f"\n{len(jobs)} configurations, {npass} passes total "
          f"(~{npass * 125 / 60:.0f} min)\n", flush=True)

    out = {k: [] for k in ("n_strips", "spacing_beams", "spacing_deg",
                           "n_eff", "n_1pct", "npix", "area", "dens",
                           "el_min", "el_max")}
    for n, sp in jobs:
        t0 = time.time()
        els = elevations(n, sp, fwhm)
        r = mode_count(build_operator(els))
        dens = 100.0 * r["n_eff"] / r["area"]
        for k, v in (("n_strips", n), ("spacing_beams", sp),
                     ("spacing_deg", sp * fwhm), ("dens", dens),
                     ("el_min", float(els.min())), ("el_max", float(els.max()))):
            out[k].append(v)
        for k in ("n_eff", "n_1pct", "npix", "area"):
            out[k].append(r[k])
        tag = "  <- published" if (n, sp) == (REF_COUNT, REF_SPACING) else ""
        print(f"  N={n} spacing={sp:.2f} beams ({sp * fwhm:4.2f} deg, "
              f"el {els.min():.1f}-{els.max():.1f})  "
              f"tr(WA) {r['n_eff']:7.2f}  {r['area']:6.1f} deg^2  "
              f"{dens:6.2f}/100deg^2  n>1% {r['n_1pct']:3d}"
              f"  ({time.time() - t0:.0f} s){tag}", flush=True)

    out = {k: np.asarray(v) for k, v in out.items()}
    out["freq_mhz"] = FREQ_MHZ
    out["fwhm_deg"] = fwhm
    out["ref_spacing_beams"], out["ref_count"] = REF_SPACING, REF_COUNT
    base = args.out or "hi_ladder_f350"
    path = os.path.join(_HERE, "results", f"{base}_ns{X.NSIDE}.npz")
    np.savez(path, **out)
    print(f"\nwrote {path}", flush=True)


if __name__ == "__main__":
    main()
