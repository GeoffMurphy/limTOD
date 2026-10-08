"""Mode counts with and without the pinned prior, from cached operators only.

Table ``tab:modes`` quotes tr(WA) at channel 0 for each strategy, and the paper
says in Sec. "Accessible sky modes" that modes per unit sky is a fair
comparison only at *fixed* sigma_p. It was not: each strategy used its own
patch variance (see ska_hi_prior). This reports both, side by side, so the
cost of the correction is visible rather than inferred.

No TOD, no solves -- one eigendecomposition per strategy.

    python ska_hi_modecheck.py --strategies drift drift8 drift12 raster
    source freq_env.sh && python ska_hi_modecheck.py \
        --strategies drift drift12m raster
"""
from __future__ import annotations

import os
import sys

import healpy as hp
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import ska_hi_experiment as X
import ska_hi_analysis as A
import ska_hi_modes as MD
import ska_hi_prior as prior
from ska_common import gdsm_equatorial_sky_model

THRESH = 0.01       # the ">1% of lambda_max" count Table 2's second column uses


def main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--strategies", nargs="*",
                    default=["drift", "drift8", "drift12", "raster"])
    ap.add_argument("--channel", type=int, default=0)
    args = ap.parse_args()

    freqs = X.channel_freqs()
    f = freqs[args.channel]
    print(f"band {X.F_LO_MHZ:.0f}-{X.F_HI_MHZ:.0f} MHz, nside {X.NSIDE}, "
          f"channel {args.channel} = {f:.3f} MHz", flush=True)

    if prior.load() is None:
        if prior.REF_STRATEGY not in args.strategies:
            sys.exit(f"no band reference yet, and '{prior.REF_STRATEGY}' is "
                     "not in --strategies, so it cannot be built here")
        prior.build(A.channel_operators(prior.REF_STRATEGY))
    ref = prior.load(required=True)
    print(f"reference sigma_p = {np.sqrt(ref['sigma_p2'][args.channel]):.5f} K "
          f"({ref['ref_strategy']}, channel {args.channel})\n", flush=True)

    hdr = (f"{'strategy':10s} {'area':>8s} {'thr':>4s} "
           f"{'n_eff own':>10s} {'n_eff pin':>10s} {'change':>8s} "
           f"{'/100deg2 own':>13s} {'/100deg2 pin':>13s} {'change':>8s}")
    print(hdr)
    print("-" * len(hdr))
    rows = {}
    for s in args.strategies:
        s_ops = A.channel_operators(s)
        mm = s_ops[args.channel]
        pix = np.asarray(mm.pixel_indices)
        truth = gdsm_equatorial_sky_model(freq=f, nside=X.NSIDE)[pix]
        w = np.maximum(np.linalg.eigvalsh(MD.information_matrix(mm, truth))[::-1],
                       0.0)
        area = len(pix) * hp.nside2pixarea(X.NSIDE, degrees=True)
        n_thr = int((w / w[0] > THRESH).sum())     # prior-independent
        out = {}
        for key, pinned in (("own", False), ("pin", True)):
            prior.PINNED = pinned
            prior.set_strategy(s, s_ops, verbose=pinned)
            si = prior.s_inv(f, truth, label=s)
            out[key] = float(np.sum(w / (w + si)))
        prior.PINNED = False
        d_own, d_pin = 100 * out["own"] / area, 100 * out["pin"] / area
        rows[s] = dict(area=area, n_thr=n_thr, npix=len(pix),
                       n_eff_own=out["own"], n_eff_pin=out["pin"],
                       dens_own=d_own, dens_pin=d_pin)
        print(f"{s:10s} {area:8.1f} {n_thr:4d} {out['own']:10.2f} "
              f"{out['pin']:10.2f} {100 * (out['pin'] / out['own'] - 1):+7.1f}% "
              f"{d_own:13.2f} {d_pin:13.2f} "
              f"{100 * (d_pin / d_own - 1):+7.1f}%")

    # the density match that the matched-density argument rests on
    if "raster" in rows:
        print()
        for s in args.strategies:
            if s == "raster":
                continue
            for key in ("own", "pin"):
                r = rows[s][f"dens_{key}"] / rows["raster"][f"dens_{key}"]
                print(f"  {s:9s} / raster density, {key}: {r:6.4f} "
                      f"({100 * (r - 1):+.1f}%)")
    np.savez(os.path.join(_HERE, "results",
                          f"hi_modecheck_f{X.F_LO_MHZ:.0f}_{X.F_HI_MHZ:.0f}"
                          f"_ns{X.NSIDE}.npz"),
             channel=args.channel, freq_mhz=f,
             strategies=np.array(args.strategies, dtype=object),
             **{f"{s}_{k}": v for s, d in rows.items() for k, v in d.items()})
    print("\nwrote results/hi_modecheck_"
          f"f{X.F_LO_MHZ:.0f}_{X.F_HI_MHZ:.0f}_ns{X.NSIDE}.npz", flush=True)


if __name__ == "__main__":
    main()
