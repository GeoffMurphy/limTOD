"""Modes measured *in the region the residual is evaluated on*.

``tr(WA)`` as reported everywhere in this series is a PATCH AVERAGE over each
strategy's own footprint. The residual, however, is measured on the small
common interior region. If a strategy's measured degrees of freedom are not
uniform across its footprint, the patch average says nothing about how well
that particular region is constrained -- and the two strategies have very
different footprint sizes (at 700 MHz the interior is 24% of the drift's patch
but 9% of the raster's).

The local quantity is exact and needs no new simulation:

    N_eff(region) = sum_{i in region} (WA)_ii,   WA = (M + S^-1)^-1 M,

the degrees of freedom the data constrains *there*. Summed over the whole patch
it reduces to the usual tr(WA).

    python local_modes.py <experiment-npz> <strategy> [<strategy> ...]
"""
import os, sys, pickle
import numpy as np, healpy as hp
import ska_hi_experiment as X, ska_hi_ladder as L
from ska_common import gdsm_equatorial_sky_model, ska_beam_fwhm_deg

def info_matrix(mm, truth):
    nv = X.WHITE_VAR * (1e-3 * float(np.mean(truth))) ** 2
    M = None
    for o in L.operator_list(mm):
        A = np.asarray(o, float)
        var = X.WHITE_VAR * (A @ truth) ** 2 + nv
        t = A.T @ (A / var[:, None])
        M = t if M is None else M + t
    return M

def main():
    exp_path, strategies = sys.argv[1], sys.argv[2:]
    e = np.load(exp_path)
    common = np.asarray(e["common"]); interior = np.asarray(e["interior"], bool)
    region = common[interior]
    pxa = hp.nside2pixarea(X.NSIDE, degrees=True)
    f = X.channel_freqs()[0]
    print(f"{exp_path}")
    print(f"channel 0 = {f:.2f} MHz, nside {X.NSIDE}, FWHM {ska_beam_fwhm_deg(f):.3f} deg")
    print(f"evaluation region: {len(region)} px = {len(region)*pxa:.1f} deg^2\n")
    print(f"{'strategy':8s}{'patch tr(WA)':>14}{'patch deg2':>12}{'patch/100':>11}"
          f"{'LOCAL N_eff':>13}{'LOCAL/100':>11}{'region/patch':>14}")
    for s in strategies:
        p = os.path.join(X.CACHE_DIR, f"op_{s}_f{f:07.3f}_ns{X.NSIDE}.pkl")
        with open(p, "rb") as fh:
            mm = pickle.load(fh)
        pix = np.asarray(mm.pixel_indices)
        truth = gdsm_equatorial_sky_model(freq=f, nside=X.NSIDE)[pix]
        M = info_matrix(mm, truth)
        s_inv = 1.0 / max(float(np.std(truth)), 1e-3) ** 2
        WA = np.linalg.solve(M + s_inv * np.eye(len(pix)), M)
        d = np.diag(WA)
        sel = np.isin(pix, region)
        patch_area = len(pix) * pxa
        loc, loc_area = float(d[sel].sum()), int(sel.sum()) * pxa
        print(f"{s:8s}{d.sum():14.2f}{patch_area:12.1f}{100*d.sum()/patch_area:11.2f}"
              f"{loc:13.2f}{100*loc/loc_area:11.2f}{loc_area/patch_area:14.3f}")

if __name__ == "__main__":
    main()
