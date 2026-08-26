"""Experiment 006 driver: inject HI, clean, measure and correct the loss.

Assumes the operator and TOD caches are warm (``python ska_hi_experiment.py``).
Everything after that is solves, which cost ~0.6 s each.

    /home/geoff/gibbs_venv_312/bin/python run_hi_experiment.py [--mocks N]

Writes ``results/hi_experiment_<band>.npz`` with the cubes, the transfer
functions and every power spectrum, so the plotting and write-up can be redone
without re-solving.
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

NMODES_GRID = (1, 2, 3, 4, 6, 8, 10)
STRATEGIES = ("drift", "raster")


def hi_on_union(cfg, union_pix, verbose=True):
    """Project one HI realisation onto the union patch: ``(nchan, npix)`` in K."""
    cube, box = build_hi_cube(cfg, union_pix, verbose=verbose)
    return project_to_healpix(cube, box, cfg, union_pix)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mocks", type=int, default=20,
                    help="independent HI realisations for the transfer function")
    ap.add_argument("--strategies", nargs="*", default=list(STRATEGIES))
    ap.add_argument("--margin", type=float, default=3.0,
                    help="degrees from the patch edge to discard when "
                         "evaluating common-resolution cubes (zero-padded "
                         "smoothing pulls the map down near the boundary)")
    args = ap.parse_args()

    freqs = X.channel_freqs()
    cfg = X.band_config()
    dr = A.channel_dr_mpc(cfg)
    print(f"band {X.F_LO_MHZ:.0f}-{X.F_HI_MHZ:.0f} MHz, {X.NCHAN} channels, "
          f"nside {X.NSIDE}, dr = {dr:.2f} Mpc/channel", flush=True)

    # ---- operators and the common patch --------------------------------
    ops = {s: A.channel_operators(s) for s in args.strategies}
    for s in args.strategies:
        n = [len(m.pixel_indices) for m in ops[s]]
        print(f"  {s}: {min(n)}-{max(n)} px per channel", flush=True)
    common = A.common_patch([ops[s] for s in args.strategies])
    union = np.unique(np.concatenate(
        [np.asarray(m.pixel_indices) for s in args.strategies for m in ops[s]]))
    interior = A.interior_mask(common, X.NSIDE, margin_deg=args.margin)
    print(f"  common patch: {len(common)} px; union {len(union)} px; "
          f"interior (>{args.margin:.1f} deg from edge) {interior.sum()} px",
          flush=True)

    # Index of every union pixel, so per-channel patches can be sliced out.
    lookup = {p: i for i, p in enumerate(union)}
    def on_patch(cube_union, pix):
        return cube_union[:, [lookup[p] for p in np.asarray(pix)]]

    # ---- the true HI realisation, plus independent mocks ---------------
    print("\n[HI] building true realisation", flush=True)
    hi_true_union = hi_on_union(cfg, union)
    mocks_union = []
    for j in range(args.mocks):
        mcfg = HIBandConfig(f_lo_mhz=X.F_LO_MHZ, f_hi_mhz=X.F_HI_MHZ,
                            nchan=X.NCHAN, nside=X.NSIDE, seed=1000 + j)
        mocks_union.append(hi_on_union(mcfg, union, verbose=(j == 0)))
        if (j + 1) % 5 == 0:
            print(f"[HI] mock {j + 1}/{args.mocks}", flush=True)

    fg_truth = {f: gdsm_equatorial_sky_model(freq=f, nside=X.NSIDE)
                for f in freqs}

    out = dict(freqs=freqs, common=common, dr_mpc=dr, interior=interior,
               margin_deg=args.margin, nmodes_grid=np.asarray(NMODES_GRID))

    for strategy in args.strategies:
        print(f"\n=== {strategy} ===", flush=True)
        gain, white = X.replay_noise(strategy)
        t0 = time.time()

        data_cube, hi_map_cube = [], []
        mock_cubes = [[] for _ in range(args.mocks)]
        for i, f in enumerate(freqs):
            mm = ops[strategy][i]
            pix = np.asarray(mm.pixel_indices)
            truth = fg_truth[f][pix]
            tods = X.simulate_foreground(strategy, f, verbose=False)

            # Strip the multiplicative noise to recover the clean sky TOD,
            # add the injected HI, then re-apply the same noise realisation.
            hi_patch = on_patch(hi_true_union, pix)[i]
            clean_tod = [np.asarray(t, float) / ((1 + gain[j]) * (1 + white[j]))
                         for j, t in enumerate(tods["TOD_group"])]
            total = [(c + np.asarray(o) @ hi_patch)
                     * (1 + gain[j]) * (1 + white[j])
                     for j, (c, o) in enumerate(zip(clean_tod,
                                                    mm.Tsys_operators))]

            mu = np.full_like(truth, float(np.mean(truth)))   # flat prior
            sel = np.isin(pix, common)
            data_cube.append(A.solve(mm, total, truth, mu)[sel])
            hi_map_cube.append(A.response(mm, hi_patch, truth)[sel])
            for j in range(args.mocks):
                mock_patch = on_patch(mocks_union[j], pix)[i]
                mock_cubes[j].append(A.response(mm, mock_patch, truth)[sel])
            if (i + 1) % 8 == 0:
                print(f"  solved {i + 1}/{len(freqs)} channels "
                      f"({time.time() - t0:.0f} s)", flush=True)

        data_cube = np.asarray(data_cube)
        hi_map_cube = np.asarray(hi_map_cube)
        mock_cubes = [np.asarray(m) for m in mock_cubes]
        hi_true_common = on_patch(hi_true_union, common)

        print(f"  HI: true rms {hi_true_common.std() * 1e3:.4f} mK, "
              f"map-made rms {hi_map_cube.std() * 1e3:.4f} mK "
              f"({hi_map_cube.std() / hi_true_common.std():.3f}x)", flush=True)

        # Three variants, so the reconvolution's effect is never confounded
        # with the pixel-set restriction it forces:
        #   ""     as run, full common patch  (the originally published numbers)
        #   "int_" as run, interior pixels    (the like-for-like control)
        #   "cr_"  common resolution, interior pixels
        print("  reconvolving to common resolution...", flush=True)
        def cr(c):
            return A.common_resolution(c, freqs, common, X.NSIDE)
        variants = {
            "": (data_cube, hi_map_cube, mock_cubes, slice(None)),
            "int_": (data_cube, hi_map_cube, mock_cubes, interior),
            "cr_": (cr(data_cube), cr(hi_map_cube), [cr(m) for m in mock_cubes],
                    interior),
        }

        res = {}
        for tag, (dc, hm, mocks, sel_pix) in variants.items():
            dc, hm = dc[:, sel_pix], hm[:, sel_pix]
            mocks = [m[:, sel_pix] for m in mocks]
            k, p_true = A.pk_par(hi_true_common[:, sel_pix], dr_mpc=dr)
            _, p_mapmade = A.pk_par(hm, dr_mpc=dr)
            res[f"{tag}k"] = k
            res[f"{tag}p_true"] = p_true
            res[f"{tag}p_mapmade"] = p_mapmade
            for nm in NMODES_GRID:
                _, tf, tf_real = A.transfer_function(dc, mocks, nm, dr,
                                                     verbose=False)
                clean = A.pca_clean(dc, nm)
                _, p_auto = A.pk_par(clean, dr_mpc=dr)
                _, p_cross = A.pk_par(clean, hm, dr_mpc=dr)
                res[f"{tag}tf_{nm}"] = tf
                res[f"{tag}tf_scatter_{nm}"] = tf_real.std(axis=0)
                res[f"{tag}p_clean_auto_{nm}"] = p_auto
                res[f"{tag}p_clean_cross_{nm}"] = p_cross
                res[f"{tag}p_corrected_{nm}"] = p_cross / tf

                # How far the post-cleaning residual sits above the HI it
                # contains. This governs whether p_corrected means anything:
                # the cross-power estimator only isolates surviving HI when the
                # residual is not overwhelmingly larger than the signal.
                # Measured here it is tens to hundreds, so p_corrected is
                # noise-dominated and must NOT be read as a recovered HI
                # spectrum. The transfer function itself is unaffected --
                # injection is differential and mock-averaged.
                res[f"{tag}residual_over_hi_{nm}"] = np.sqrt(p_auto / p_mapmade)
            print(f"    [{tag or 'as-run':7s}] {dc.shape[1]:3d} px  "
                  f"residual/HI at 4 modes "
                  f"{np.median(res[f'{tag}residual_over_hi_4']):7.1f}x  "
                  f"median T {np.median(res[f'{tag}tf_4']):.3f}", flush=True)

        for key, val in res.items():
            out[f"{strategy}_{key}"] = val
        out[f"{strategy}_data_cube"] = data_cube
        out[f"{strategy}_hi_mapmade"] = hi_map_cube
        out[f"{strategy}_hi_true"] = hi_true_common
        print(f"  {strategy} done in {time.time() - t0:.0f} s", flush=True)

    os.makedirs(os.path.join(_HERE, "results"), exist_ok=True)
    path = os.path.join(_HERE, "results",
                        f"hi_experiment_f{X.F_LO_MHZ:.0f}_{X.F_HI_MHZ:.0f}"
                        f"_nc{X.NCHAN}_ns{X.NSIDE}.npz")
    np.savez(path, **out)
    print(f"\nwrote {path}", flush=True)


if __name__ == "__main__":
    main()
