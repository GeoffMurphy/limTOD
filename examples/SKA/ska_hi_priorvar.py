"""How much of experiment 006 is the prior variance? (PLAN item 2, narrowed.)

PLAN.md item 2 asks for a flat-prior audit of every headline number. For
experiment 006 that audit is already done by construction: the prior *mean* is
flat, a constant at the patch mean, set in ``run_hi_experiment.py`` and
``ska_hi_ablation.py``. The beam-smoothed-truth mean that ``audit_flat_prior.py``
exists to test belongs to experiments 004 and 005.

What is left is the prior *variance*. The solve uses

    S^-1 = I / var(fg_truth_patch),

one scalar on every pixel, and that scalar is the only remaining place the
prior knows anything about the truth. It is also the crossover in

    tr(WA) = sum_i lambda_i / (lambda_i + S^-1),

so it moves the measured mode count and the beam + prior floor *together* --
which matters more for the design figure than for any single residual, because
the design figure puts the mode count on an axis. If the two anchors slide
along the trade-space curve as the prior loosens, the curve is a statement
about the prior; if they move off it, the curve is a statement about the
survey. That is the question this script answers.

Method. Sweep ``A.PRIOR_VAR_SCALE`` (prior variance in units of the truth patch
variance; 1.0 is every published 006 number) and at each value measure:

  y  post-clean residual / HI at 4 modes removed, common resolution, interior
     pixels -- the design figure's y axis. Re-solved.
  y_floor  the same with the noise realisation switched off, so the beam +
     prior floor can be separated from the total. Re-solved.
  x  tr(WA) per 100 deg^2 at channel 0 -- the design figure's x axis. NOT
     re-solved: the eigenvalues of A^T N^-1 A do not depend on the prior, so
     tr(WA) at any scale follows analytically from the spectrum
     ``ska_hi_modes.py`` already cached. Only the crossover moves.

Cost is 3 solves per channel per strategy per scale (floor, total, HI
response), ~2 min per scale for both strategies. ``--mocks N`` adds N more
solves per channel and measures two further things, which is what it takes to
decide whether a residual ratio that *falls* at tight prior is a real gain:

  T(k_par)  what the PCA clean does to the map-made HI, by mock injection.
  P_mapmade / P_true  what the MAP-MAKER alone does to the HI, before any
     cleaning. ``residual/HI`` divides by P_mapmade, so a prior that suppresses
     the HI in the map flatters the ratio; this is the term that would show it.

Without ``--mocks`` neither is computed and residual/HI is
``sqrt(P_auto / P_mapmade)``, which needs no mocks.

Caveats, both real:

* Loosening the prior is not free. With ``regularization=1e-12`` the solve has
  essentially nothing holding down the directions the data does not constrain,
  so at large scale the unmeasured modes are limited only by round-off. A
  residual that grows at the loose end is the expected behaviour, not a bug --
  it is the other half of the trade the prior is making.
* This sweeps the prior's *amplitude*, not its *shape*. The prior is diagonal
  and isotropic at every scale, so nothing here tests what a correlated or
  angular-structured prior would do.

Writes ``results/hi_priorvar_<band>.npz``.

    /home/geoff/gibbs_venv_312/bin/python ska_hi_priorvar.py [--scales ...]
"""

from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import ska_hi_experiment as X
import ska_hi_analysis as A
from ska_common import gdsm_equatorial_sky_model
from ska_hi_mock import HIBandConfig, build_hi_cube, project_to_healpix

SCALES = (0.01, 0.1, 1.0, 10.0, 100.0, 1000.0)
NMODES_GRID = (1, 2, 4, 10)
MARGIN_DEG = 3.0
STRATEGIES = ("drift", "raster")


def mode_density(modes, strategy, scales):
    """tr(WA) per 100 deg^2 against prior scale, analytically.

    ``A^T N^-1 A`` does not contain the prior, so its eigenvalues are fixed and
    the whole scale dependence is the crossover S^-1 -> S^-1 / scale. No
    re-solve is needed, which is the only reason the x axis is cheap.
    """
    w = np.asarray(modes[f"{strategy}_eig_abs"], float)
    s_inv = float(modes[f"{strategy}_s_inv"])
    area = float(modes[f"{strategy}_area_deg2"])
    n_eff = np.array([float(np.sum(w / (w + s_inv / c))) for c in scales])
    return n_eff, 100.0 * n_eff / area


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scales", type=float, nargs="*", default=list(SCALES),
                    help="prior variance in units of the truth patch variance")
    ap.add_argument("--strategies", nargs="*", default=list(STRATEGIES))
    ap.add_argument("--mocks", type=int, default=0,
                    help="HI realisations for T(k_par); 0 skips it. "
                         "run_hi_experiment.py uses 20")
    args = ap.parse_args()
    scales = np.asarray(sorted(args.scales), float)

    freqs = X.channel_freqs()
    cfg = X.band_config()
    dr = A.channel_dr_mpc(cfg)
    ops = {s: A.channel_operators(s) for s in args.strategies}
    common = A.common_patch([ops[s] for s in args.strategies])
    union = np.unique(np.concatenate(
        [np.asarray(m.pixel_indices) for s in args.strategies for m in ops[s]]))
    interior = A.interior_mask(common, X.NSIDE, margin_deg=MARGIN_DEG)
    print(f"common {len(common)} px, interior {interior.sum()} px, "
          f"dr {dr:.2f} Mpc, scales {list(scales)}", flush=True)

    cube, box = build_hi_cube(cfg, union, verbose=False)
    hi_u = project_to_healpix(cube, box, cfg, union)
    look = {p: i for i, p in enumerate(union)}

    # Independent HI realisations for the transfer function. Same seeds as
    # run_hi_experiment.py (1000 + j), so a shared mock is the same field.
    mocks_u = []
    for j in range(args.mocks):
        mcfg = HIBandConfig(f_lo_mhz=X.F_LO_MHZ, f_hi_mhz=X.F_HI_MHZ,
                            nchan=X.NCHAN, nside=X.NSIDE, seed=1000 + j)
        mc, mb = build_hi_cube(mcfg, union, verbose=False)
        mocks_u.append(project_to_healpix(mc, mb, mcfg, union))
        print(f"  mock {j + 1}/{args.mocks}", flush=True)
    fg = {f: gdsm_equatorial_sky_model(freq=f, nside=X.NSIDE) for f in freqs}

    modes_path = os.path.join(
        _HERE, "results", f"hi_modes_f{X.F_LO_MHZ:.0f}_{X.F_HI_MHZ:.0f}"
                          f"_nc{X.NCHAN}_ns{X.NSIDE}.npz")
    modes = np.load(modes_path)

    out = dict(scales=scales, common=common, interior=interior,
               margin_deg=MARGIN_DEG, nmodes_grid=np.asarray(NMODES_GRID),
               freq_mhz=modes["freq_mhz"])

    for strategy in args.strategies:
        n_eff, dens = mode_density(modes, strategy, scales)
        out[f"{strategy}_n_eff"] = n_eff
        out[f"{strategy}_dens"] = dens
        out[f"{strategy}_area_deg2"] = float(modes[f"{strategy}_area_deg2"])
        print(f"\n=== {strategy} ===  tr(WA)/100deg^2 across the sweep: "
              + "  ".join(f"{d:.1f}" for d in dens), flush=True)

        g, w = X.replay_noise(strategy)
        for si, scale in enumerate(scales):
            A.PRIOR_VAR_SCALE = float(scale)
            t0 = time.time()
            floor_c, total_c, hi_c = [], [], []
            mock_c = [[] for _ in range(args.mocks)]
            for i, f in enumerate(freqs):
                mm = ops[strategy][i]
                pix = np.asarray(mm.pixel_indices)
                truth = fg[f][pix]
                sel = np.isin(pix, common)
                hi = hi_u[i, [look[p] for p in pix]]
                tods = X.simulate_foreground(strategy, f, verbose=False)
                clean = [np.asarray(t, float) / ((1 + g[j]) * (1 + w[j]))
                         for j, t in enumerate(tods["TOD_group"])]
                hitod = [np.asarray(o) @ hi for o in mm.Tsys_operators]
                mu = np.full_like(truth, float(np.mean(truth)))
                floor_c.append(A.solve(
                    mm, [c + h for c, h in zip(clean, hitod)], truth, mu)[sel])
                total_c.append(A.solve(
                    mm, [(c + h) * (1 + g[j]) * (1 + w[j])
                         for j, (c, h) in enumerate(zip(clean, hitod))],
                    truth, mu)[sel])
                hi_c.append(A.response(mm, hi, truth)[sel])
                for j in range(args.mocks):
                    mock_c[j].append(A.response(
                        mm, mocks_u[j][i, [look[p] for p in pix]], truth)[sel])

            def cr(c):
                return A.common_resolution(np.asarray(c), freqs, common,
                                           X.NSIDE)[:, interior]
            floor_cube, total_cube, hi_cube = cr(floor_c), cr(total_c), cr(hi_c)
            kk, p_hi = A.pk_par(hi_cube, dr_mpc=dr)
            if args.mocks:
                # What the MAP-MAKER does to the HI, before any cleaning. This
                # is the denominator of residual/HI, so it has to be reported
                # beside it or a prior that suppresses the HI reads as a win.
                hi_true_int = hi_u[:, [look[p] for p in common]][:, interior]
                _, p_true = A.pk_par(hi_true_int, dr_mpc=dr)
                out.setdefault(f"{strategy}_p_true", p_true)
                out.setdefault("k_par", kk)
                out.setdefault(f"{strategy}_p_mapmade",
                               np.zeros((len(scales), len(kk))))
                out[f"{strategy}_p_mapmade"][si] = p_hi
                out.setdefault(f"{strategy}_hi_kept", np.zeros(len(scales)))
                out[f"{strategy}_hi_kept"][si] = float(
                    np.median(np.sqrt(p_hi / p_true)))
                mock_cubes = [cr(m) for m in mock_c]
            for nm in NMODES_GRID:
                if args.mocks:
                    _, tf, _tfr = A.transfer_function(total_cube, mock_cubes,
                                                      nm, dr, verbose=False)
                    out.setdefault(f"{strategy}_tf_{nm}",
                                   np.zeros((len(scales), len(kk))))
                    out[f"{strategy}_tf_{nm}"][si] = tf
                    out.setdefault(f"{strategy}_tfmed_{nm}",
                                   np.zeros(len(scales)))
                    out[f"{strategy}_tfmed_{nm}"][si] = float(np.median(tf))
                _, p_floor = A.pk_par(A.pca_clean(floor_cube, nm), dr_mpc=dr)
                _, p_total = A.pk_par(A.pca_clean(total_cube, nm), dr_mpc=dr)
                for tag, p in (("floor", p_floor), ("total", p_total)):
                    key = f"{strategy}_{tag}_{nm}"
                    out.setdefault(key, np.zeros(len(scales)))
                    out[key][si] = float(np.median(np.sqrt(p / p_hi)))
            out.setdefault(f"{strategy}_hi_rms", np.zeros(len(scales)))
            out[f"{strategy}_hi_rms"][si] = float(hi_cube.std())
            extra = ""
            if args.mocks:
                extra = (f"  HI kept by mapmaker {out[f'{strategy}_hi_kept'][si]:.3f}"
                         f"  median T(k) at 4 modes "
                         f"{out[f'{strategy}_tfmed_4'][si]:.3f}")
            print(f"  scale {scale:8.3g}  tr(WA)/100deg^2 {dens[si]:6.2f}  "
                  f"residual/HI at 4 modes: floor {out[f'{strategy}_floor_4'][si]:8.1f}x"
                  f"  total {out[f'{strategy}_total_4'][si]:8.1f}x{extra}"
                  f"   ({time.time() - t0:.0f} s)", flush=True)
        A.PRIOR_VAR_SCALE = 1.0

    out["n_mocks"] = args.mocks
    suffix = "_tk" if args.mocks else ""
    path = os.path.join(_HERE, "results",
                        f"hi_priorvar{suffix}_f{X.F_LO_MHZ:.0f}_{X.F_HI_MHZ:.0f}"
                        f"_nc{X.NCHAN}_ns{X.NSIDE}.npz")
    np.savez(path, **out)
    print(f"\nwrote {path}", flush=True)


if __name__ == "__main__":
    main()
