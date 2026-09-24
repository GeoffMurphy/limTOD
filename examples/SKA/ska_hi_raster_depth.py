"""How narrow must the raster's azimuth throw be to match the drift's depth?

PLAN item 5 / HANDOFF queued item 1. The published comparison is not
depth-matched: both strategies produce 5400 samples, but the drift spreads them
over 316 pixels and the raster over 684, so the raster is ~2.2x shallower per
pixel. Every residual comparison in the paper therefore confounds *geometry*
with *depth*, and the obvious referee question is whether cross-linking is
doing the work or whether the raster is simply being read at a different noise
level.

The earlier answer was an analytic noise rescaling. This measures it instead:
narrow ``RASTER_AZ_HALFWIDTH`` until the scan *natively* selects the drift's
pixel count, then hand that geometry to the full HI experiment.

Tier 1, like ``ska_hi_ladder.py``: the operator alone, one channel, no TOD.
The pixel count and tr(WA) are properties of ``A``, so the TOD is not needed to
choose the throw -- only to measure the residual afterwards.

Two things to watch in the output, because narrowing the throw is not a pure
depth knob:

  * it removes sky area, so modes/deg^2 can *rise* while total modes fall;
  * the two passes still cross at 73.9 deg at any throw -- azimuth centre sets
    the crossing angle, not the width -- but a narrow throw shrinks the region
    where both passes overlap, which is where cross-linking actually acts.

Reported per throw so the choice is made on measurement, not on the assumption
that only depth changed.

    python ska_hi_raster_depth.py --halfwidths 7.5 5 4 3 2.5 2
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import healpy as hp
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import ska_hi_experiment as X
from ska_common import (SKA_LAT, SKA_LON, SKA_HGT, ska_beam_func,
                        ska_beam_fwhm_deg, constant_elevation_scan)
# the operator/eigen machinery is identical to the drift ladder's; importing it
# keeps one definition of "how this series counts modes"
from ska_hi_ladder import FREQ_MHZ, mode_count

# what we are matching: the published drift, from hi_modes_f350_400_nc32_ns64
DRIFT_NPIX, DRIFT_AREA, DRIFT_NEFF = 316, 265.2, 38.78
HALFWIDTHS = (7.5, 5.0, 4.0, 3.0, 2.5, 2.0)


def build_operator(halfwidth_deg):
    """The two published raster passes, with the azimuth throw overridden."""
    from limTOD import HPW_mapmaking
    from limTOD.simulator import generate_LSTs_deg

    lst_g, az_g, el_g = [], [], []
    for az_centre, start_h in X.RASTER_PASSES:
        tlist, azlist, _ = constant_elevation_scan(
            X.RASTER_PASS_S, az_centre, halfwidth_deg, X.RASTER_SWEEP_S,
            dt=X.DT, t0_s=start_h * 3600.0)
        lst_g.append(np.asarray(generate_LSTs_deg(
            SKA_LAT, SKA_LON, SKA_HGT, tlist, X.RASTER_UTC), float))
        az_g.append(np.asarray(azlist, float))
        el_g.append(np.full(len(tlist), float(X.RASTER_EL)))
    return HPW_mapmaking(
        beam_map=ska_beam_func(freq=FREQ_MHZ, nside=X.NSIDE),
        LST_deg_list_group=lst_g, lat_deg=SKA_LAT,
        azimuth_deg_list_group=az_g, elevation_deg_list_group=el_g,
        threshold=X.THRESHOLD, nside_target=X.NSIDE,
        beam_truncate_frac_thres=X.TRUNCATE_FRAC_THRES)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--halfwidths", type=float, nargs="*",
                    default=list(HALFWIDTHS))
    ap.add_argument("--out", default="hi_raster_depth")
    args = ap.parse_args()

    print(f"{FREQ_MHZ:.2f} MHz, FWHM {ska_beam_fwhm_deg(FREQ_MHZ):.3f} deg, "
          f"nside {X.NSIDE}, el {X.RASTER_EL}, "
          f"{len(X.RASTER_PASSES)} passes of {X.RASTER_PASS_S:.0f} s",
          flush=True)
    print(f"TARGET (published drift): {DRIFT_NPIX} px, {DRIFT_AREA:.1f} deg^2, "
          f"tr(WA) {DRIFT_NEFF:.2f}\n", flush=True)

    keys = ("halfwidth_deg", "n_eff", "n_1pct", "npix", "area", "dens",
            "samples_per_px")
    out = {k: [] for k in keys}
    for hw in args.halfwidths:
        t0 = time.time()
        r = mode_count(build_operator(hw))
        dens = 100.0 * r["n_eff"] / r["area"]
        # both passes together, which is what sets the depth per pixel
        nsamp = len(X.RASTER_PASSES) * X.RASTER_PASS_S / X.DT
        for k, v in (("halfwidth_deg", hw), ("dens", dens),
                     ("samples_per_px", nsamp / r["npix"])):
            out[k].append(v)
        for k in ("n_eff", "n_1pct", "npix", "area"):
            out[k].append(r[k])
        tag = "  <- published" if abs(hw - X.RASTER_AZ_HALFWIDTH) < 1e-9 else ""
        if abs(r["npix"] - DRIFT_NPIX) <= 0.05 * DRIFT_NPIX:
            tag += "  <== DEPTH-MATCHED (within 5% of the drift)"
        print(f"  halfwidth {hw:4.2f} deg  tr(WA) {r['n_eff']:7.2f}  "
              f"{r['npix']:4d} px  {r['area']:6.1f} deg^2  "
              f"{dens:6.2f}/100deg^2  n>1% {r['n_1pct']:3d}  "
              f"{nsamp / r['npix']:6.1f} samp/px  "
              f"({time.time() - t0:.0f} s){tag}", flush=True)

    out = {k: np.asarray(v) for k, v in out.items()}
    out["freq_mhz"], out["nside"] = FREQ_MHZ, X.NSIDE
    out["drift_npix"], out["drift_area"] = DRIFT_NPIX, DRIFT_AREA
    out["published_halfwidth_deg"] = X.RASTER_AZ_HALFWIDTH
    path = os.path.join(_HERE, "results", f"{args.out}_ns{X.NSIDE}.npz")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    np.savez(path, **out)
    print(f"\nwrote {path}", flush=True)


if __name__ == "__main__":
    main()
