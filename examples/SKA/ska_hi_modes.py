"""How many independent sky modes does each strategy actually measure?

The mode count is the mechanism behind everything else in the series: it is why
the beam + prior floor is what it is, why the drift's floor needs more PCA
modes than the raster's, and why the improvement lever is more crossing angles
rather than more time. It has been quoted throughout (drift ~17, raster ~34 at
350 MHz) without a figure.

The quantity is the eigenspectrum of the map-maker's own information matrix,

    M = A^T N^-1 A = sum_scans A_s^T diag(1/sigma_s^2) A_s,

with the *same* truth-dependent noise weighting the solves use, so the count
describes the estimator actually run rather than an idealised one. Modes are
counted above a fraction of lambda_max; the series quotes the 1% threshold, and
10% / 0.1% are recorded to show the contrast is not an artefact of where the
line is drawn.

The primary number is threshold-free. The map is m = W A s + prior, so the
resolution matrix W A has trace

    N_eff = tr(W A) = sum_i lambda_i / (lambda_i + S^-1),

the degrees of freedom the DATA constrains: each mode counts between 0 and 1
according to whether data or prior dominates it. It is exactly complementary to
the floor, r_floor = (I - W A) s, so N_pix - N_eff counts the prior-filled
directions the floor is built from. It is prior-dependent by construction --
S^-1 is the flat prior actually used, 1/var(s_truth) -- which is honest rather
than a defect, since the floor depends on the prior too.

A fixed fraction of lambda_max is NOT used for the headline: it is arbitrary,
and the drift-to-raster factor moves with it (1.75x at 10%, 1.94x at 1%, 2.39x
at 0.1%) where tr(WA) gives 3.30x. The threshold counts are kept as a
robustness line, with the participation ratio as a second threshold-free
measure.

Writes ``results/hi_modes_<band>.npz``.

    /home/geoff/gibbs_venv_312/bin/python ska_hi_modes.py
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
from ska_common import gdsm_equatorial_sky_model

THRESHOLDS = (0.10, 0.01, 0.001)


def information_matrix(mm, truth_patch):
    """A^T N^-1 A with the map-maker's own diagonal noise weights."""
    nv_floor = X.WHITE_VAR * (1e-3 * float(np.mean(truth_patch))) ** 2
    M = None
    for o in mm.Tsys_operators:
        Aop = np.asarray(o, float)
        var = X.WHITE_VAR * (Aop @ truth_patch) ** 2 + nv_floor
        term = Aop.T @ (Aop / var[:, None])
        M = term if M is None else M + term
    return M


def main():
    freqs = X.channel_freqs()
    ch = 0                      # 350.78 MHz, the band edge the series quotes
    out = dict(freq_mhz=freqs[ch], thresholds=np.asarray(THRESHOLDS))
    print(f"channel {ch}: {freqs[ch]:.2f} MHz", flush=True)

    for strategy in ("drift", "raster"):
        mm = A.channel_operators(strategy)[ch]
        pix = np.asarray(mm.pixel_indices)
        truth = gdsm_equatorial_sky_model(freq=freqs[ch], nside=X.NSIDE)[pix]
        M = information_matrix(mm, truth)
        w = np.linalg.eigvalsh(M)[::-1]          # descending
        w = np.maximum(w, 0.0)
        rel = w / w[0]
        # The prior the solves actually use: flat and diagonal at the patch
        # variance. Being a multiple of the identity it shifts every eigenvalue
        # equally, so tr(WA) follows directly in the eigenbasis.
        s_inv = 1.0 / max(float(np.std(truth)), 1e-3) ** 2
        n_eff = float(np.sum(w / (w + s_inv)))
        area = len(pix) * hp.nside2pixarea(X.NSIDE, degrees=True)
        part = float(w.sum() ** 2 / (w ** 2).sum())
        counts = {t: int((rel > t).sum()) for t in THRESHOLDS}
        out[f"{strategy}_eig"] = rel
        out[f"{strategy}_eig_abs"] = w
        out[f"{strategy}_s_inv"] = s_inv
        out[f"{strategy}_n_eff"] = n_eff
        out[f"{strategy}_n_above_prior"] = int((w > s_inv).sum())
        out[f"{strategy}_area_deg2"] = area
        out[f"{strategy}_npix"] = len(pix)
        out[f"{strategy}_participation"] = part
        for t in THRESHOLDS:
            out[f"{strategy}_n{t}"] = counts[t]
        print(f"  {strategy:7s} {len(pix):3d} px, {area:5.1f} deg^2  "
              + " ".join(f">{t*100:g}%: {counts[t]:3d}" for t in THRESHOLDS)
              + f"  |  tr(WA) {n_eff:6.1f} ({n_eff / area:.3f}/deg^2)"
              + f"  participation {part:.1f}", flush=True)

    os.makedirs(os.path.join(_HERE, "results"), exist_ok=True)
    path = os.path.join(_HERE, "results",
                        f"hi_modes_f{X.F_LO_MHZ:.0f}_{X.F_HI_MHZ:.0f}"
                        f"_nc{X.NCHAN}_ns{X.NSIDE}.npz")
    np.savez(path, **out)
    print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    main()
