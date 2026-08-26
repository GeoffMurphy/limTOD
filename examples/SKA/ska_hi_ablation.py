"""Experiment 006 ablation: what is actually blocking HI recovery?

The post-cleaning residual sits 50-500x above the HI. This asks whether that is
the beam + prior floor or the noise, by re-solving the same operators with
components removed -- the same RNG-replay trick the rest of the series uses to
split residuals, applied here to the cleaned cube instead.

Three arms, all with the HI injected:

* ``floor``     -- noiseless data. The residual is then purely sky the strategy
                   never measured, filled by the prior: the beam + prior floor.
* ``total``     -- the full data, with 1/f gain and white noise.
* ``noiseonly`` -- the noise contribution alone, as the difference between a
                   noisy and a noiseless solve of the foreground.

If ``floor`` and ``total`` agree, the noise is irrelevant to HI recovery and
the blocker is structural. Results are saved k-resolved for every mode count so
the figures can be drawn without re-solving.

Each arm is stored under the three conventions ``run_hi_experiment.py`` uses --
as-run full patch (no prefix), as-run interior pixels (``int_``), and
common-resolution interior pixels (``cr_``). The ablation's *conclusion* holds
under any of them, since the arms are compared within one convention, but the
absolute residual/HI ratios are only comparable with figure 2 under ``cr_``.

    /home/geoff/gibbs_venv_312/bin/python ska_hi_ablation.py
"""

from __future__ import annotations

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
from ska_hi_mock import build_hi_cube, project_to_healpix

ARMS = ("floor", "total", "noiseonly")
MARGIN_DEG = 3.0


def main():
    freqs = X.channel_freqs()
    cfg = X.band_config()
    dr = A.channel_dr_mpc(cfg)

    ops = {s: A.channel_operators(s) for s in ("drift", "raster")}
    common = A.common_patch([ops[s] for s in ("drift", "raster")])
    union = np.unique(np.concatenate(
        [np.asarray(m.pixel_indices) for s in ops for m in ops[s]]))
    interior = A.interior_mask(common, X.NSIDE, margin_deg=MARGIN_DEG)
    print(f"common patch {len(common)} px, interior {interior.sum()} px, "
          f"dr {dr:.2f} Mpc", flush=True)

    cube, box = build_hi_cube(cfg, union, verbose=False)
    hi_u = project_to_healpix(cube, box, cfg, union)
    look = {p: i for i, p in enumerate(union)}
    fg = {f: gdsm_equatorial_sky_model(freq=f, nside=X.NSIDE) for f in freqs}

    out = dict(nmodes_grid=np.asarray(X.__dict__.get("NMODES_GRID",
                                                     (1, 2, 3, 4, 6, 8, 10))),
               common=common, interior=interior, margin_deg=MARGIN_DEG)
    for strategy in ("drift", "raster"):
        g, w = X.replay_noise(strategy)
        t0 = time.time()
        cubes = {a: [] for a in ARMS}
        hi_mm = []
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

            cubes["floor"].append(
                A.solve(mm, [c + h for c, h in zip(clean, hitod)],
                        truth, mu)[sel])
            cubes["total"].append(
                A.solve(mm, [(c + h) * (1 + g[j]) * (1 + w[j])
                             for j, (c, h) in enumerate(zip(clean, hitod))],
                        truth, mu)[sel])
            noisy = A.solve(mm, [c * (1 + g[j]) * (1 + w[j])
                                 for j, c in enumerate(clean)], truth, mu)[sel]
            quiet = A.solve(mm, clean, truth, mu)[sel]
            cubes["noiseonly"].append(noisy - quiet)
            hi_mm.append(A.response(mm, hi, truth)[sel])
            if (i + 1) % 8 == 0:
                print(f"  [{strategy}] {i + 1}/{len(freqs)} "
                      f"({time.time() - t0:.0f} s)", flush=True)

        hi_mm = np.asarray(hi_mm)
        cubes = {a: np.asarray(c) for a, c in cubes.items()}

        # The same three conventions run_hi_experiment.py stores, so figure 3
        # can be read alongside figure 2 instead of against it:
        #   ""     as run, full common patch  (the originally published numbers)
        #   "int_" as run, interior pixels    (the like-for-like control)
        #   "cr_"  common resolution, interior pixels
        def cr(c):
            return A.common_resolution(c, freqs, common, X.NSIDE)
        print("  reconvolving to common resolution...", flush=True)
        cr_cubes = {a: cr(c) for a, c in cubes.items()}
        cr_hi = cr(hi_mm)
        variants = {
            "": (cubes, hi_mm, slice(None)),
            "int_": (cubes, hi_mm, interior),
            "cr_": (cr_cubes, cr_hi, interior),
        }

        for tag, (arm_cubes, hi_ref, sel_pix) in variants.items():
            k, p_hi = A.pk_par(hi_ref[:, sel_pix], dr_mpc=dr)
            out[f"{strategy}_{tag}k"] = k
            out[f"{strategy}_{tag}p_hi"] = p_hi
            for nm in out["nmodes_grid"]:
                for arm in ARMS:
                    cleaned = A.pca_clean(arm_cubes[arm][:, sel_pix], int(nm))
                    _, p = A.pk_par(cleaned, dr_mpc=dr)
                    out[f"{strategy}_{tag}{arm}_{nm}"] = p
                print(f"  [{strategy} {tag or 'as-run':6s}] nmodes={nm}: "
                      + ", ".join(
                          f"{arm} {np.sqrt(np.median(out[f'{strategy}_{tag}{arm}_{nm}'] / p_hi)):.1f}x"
                          for arm in ARMS), flush=True)

    os.makedirs(os.path.join(_HERE, "results"), exist_ok=True)
    path = os.path.join(_HERE, "results",
                        f"hi_ablation_f{X.F_LO_MHZ:.0f}_{X.F_HI_MHZ:.0f}"
                        f"_nc{X.NCHAN}_ns{X.NSIDE}.npz")
    np.savez(path, **out)
    print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    main()
