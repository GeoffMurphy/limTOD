"""How compressible is the beam + prior floor, and does it overlap the HI?

Motivated by a straightforward question: if PCA is not removing the floor,
would a different basis -- ICA, GMCA, NMF, a polynomial expansion -- do better?

That is answerable rather than a matter of taste. Every one of those methods
removes a rank-N subspace of the channel-channel covariance; they differ only in
how they choose it. So there are two things to measure:

1. **How low-rank is each component?** If the floor needs many modes, no
   subspace method can capture it and the answer is a different *kind* of
   filter. If it needs few, they can, and the basis is a detail.

2. **How much does the floor's subspace overlap the HI's?** If they are nearly
   the same modes, removing the floor necessarily removes the signal, and no
   choice of basis escapes that -- the fix has to be a smaller floor, or a floor
   that is modelled rather than filtered.

Principal angles answer (2): given orthonormal bases U (floor) and V (HI) for
the leading subspaces, the singular values of U^T V are the cosines of the
principal angles. cos = 1 means a shared direction; cos = 0 means orthogonal.

Writes ``results/hi_rank_<band>.npz``. Run under ``gibbs_venv_312``.
"""

from __future__ import annotations

import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import ska_hi_experiment as X
import ska_hi_analysis as A
from ska_common import gdsm_equatorial_sky_model
from ska_hi_mock import build_hi_cube, project_to_healpix


def spectrum(cube):
    """Eigenvalues (descending, normalised to the first) of the channel-channel
    covariance, plus the orthonormal eigenvectors."""
    d = np.asarray(cube, float)
    d = d - d.mean(axis=0, keepdims=True)
    C = np.cov(d)
    w, V = np.linalg.eigh(C)
    order = np.argsort(w)[::-1]
    w, V = w[order], V[:, order]
    return w / w[0], V


def n_modes_for(w, frac):
    """Modes needed to carry `frac` of the total variance."""
    return int(np.searchsorted(np.cumsum(w) / w.sum(), frac)) + 1


def principal_angles(U, V, n):
    """cos of the principal angles between the leading n-dim subspaces."""
    s = np.linalg.svd(U[:, :n].T @ V[:, :n], compute_uv=False)
    return np.clip(s, 0.0, 1.0)


def main():
    freqs = X.channel_freqs()
    cfg = X.band_config()
    ops = {s: A.channel_operators(s, verbose=False)
           for s in ("drift", "raster")}
    common = A.common_patch([ops[s] for s in ("drift", "raster")])
    union = np.unique(np.concatenate(
        [np.asarray(m.pixel_indices) for s in ops for m in ops[s]]))
    print(f"common patch {len(common)} px", flush=True)

    cube, box = build_hi_cube(cfg, union, verbose=False)
    hi_u = project_to_healpix(cube, box, cfg, union)
    look = {p: i for i, p in enumerate(union)}
    fg = {f: gdsm_equatorial_sky_model(freq=f, nside=X.NSIDE) for f in freqs}

    out = dict(freqs=freqs, common=common)
    for strategy in ("drift", "raster"):
        gain, white = X.replay_noise(strategy)
        truth_c, floor_c, hi_c = [], [], []
        for i, f in enumerate(freqs):
            mm = ops[strategy][i]
            pix = np.asarray(mm.pixel_indices)
            truth = fg[f][pix]
            sel = np.isin(pix, common)
            hi = hi_u[i, [look[p] for p in pix]]
            tods = X.simulate_foreground(strategy, f, verbose=False)
            clean = [np.asarray(t, float) / ((1 + gain[j]) * (1 + white[j]))
                     for j, t in enumerate(tods["TOD_group"])]
            mu = np.full_like(truth, float(np.mean(truth)))
            # Noiseless solve: the residual is sky the strategy never measured.
            m = A.solve(mm, clean, truth, mu)
            truth_c.append(truth[sel])
            floor_c.append((m - truth)[sel])
            hi_c.append(A.response(mm, hi, truth)[sel])
            if (i + 1) % 8 == 0:
                print(f"  [{strategy}] {i + 1}/{len(freqs)}", flush=True)

        w_fg, _ = spectrum(np.asarray(truth_c))
        w_fl, V_fl = spectrum(np.asarray(floor_c))
        w_hi, V_hi = spectrum(np.asarray(hi_c))
        # Overlap of the leading subspaces, at each rank up to 12.
        cosines = np.array([principal_angles(V_fl, V_hi, n).mean()
                            for n in range(1, 13)])

        out[f"{strategy}_w_fg"] = w_fg
        out[f"{strategy}_w_floor"] = w_fl
        out[f"{strategy}_w_hi"] = w_hi
        out[f"{strategy}_overlap"] = cosines
        # Eigenvectors too: the spectral SHAPES each component is built from.
        # Reading them is what makes "rank" concrete -- the foreground's first
        # mode should be a clean power law, the floor's should get wigglier.
        out[f"{strategy}_V_fg"] = spectrum(np.asarray(truth_c))[1][:, :8]
        out[f"{strategy}_V_floor"] = V_fl[:, :8]
        out[f"{strategy}_V_hi"] = V_hi[:, :8]
        for lbl, w in (("foreground", w_fg), ("floor", w_fl), ("HI", w_hi)):
            print(f"  {strategy:7s} {lbl:11s} modes for 90/99/99.9%: "
                  f"{n_modes_for(w,0.90):3d} {n_modes_for(w,0.99):3d} "
                  f"{n_modes_for(w,0.999):3d}", flush=True)
        print(f"  {strategy:7s} mean cos(principal angle), floor vs HI, "
              f"rank 6: {cosines[5]:.3f}", flush=True)

    os.makedirs(os.path.join(_HERE, "results"), exist_ok=True)
    path = os.path.join(_HERE, "results",
                        f"hi_rank_f{X.F_LO_MHZ:.0f}_{X.F_HI_MHZ:.0f}"
                        f"_nc{X.NCHAN}_ns{X.NSIDE}.npz")
    np.savez(path, **out)
    print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    main()
