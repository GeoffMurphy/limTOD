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
* ``gainonly``  -- the 1/f gain contribution alone, same difference with only
                   ``(1 + g)`` applied.
* ``whiteonly`` -- the white contribution alone, only ``(1 + w)`` applied.

The last two split ``noiseonly`` so the 1/f term can be quoted on the same
footing as everything else. The noise is multiplicative, ``(1 + g)(1 + w)``, so
the split is exact only up to the ``g*w`` cross term; both draws are ~1e-3, so
that term is ~1e-6 relative and ``gainonly + whiteonly`` reproduces
``noiseonly`` to well below plotting precision. The run prints the check.

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
import ska_hi_prior as prior
from ska_common import gdsm_equatorial_sky_model
from ska_hi_mock import build_hi_cube, project_to_healpix

ARMS = ("floor", "total", "noiseonly", "gainonly", "whiteonly")
# Independent-white arms, added 2026-10-05. simulate_foreground re-seeds with
# p["seed"] -- which depends on the PASS, not the frequency -- and calls
# generate_TOD one channel at a time, so every channel draws the IDENTICAL
# white realisation. Verified empirically: the ch0/ch1 TOD ratio scatters
# 1.2e-4 where independent draws would give 2.2e-3. A common 1/f gain across
# the band is physically right for one receiver and is kept; identical thermal
# noise is not, and because it follows the sky's smooth spectrum PCA removes
# most of it, so the published post-clean noise is optimistic. These arms
# redraw the white part per channel and change nothing else.
INDEP_ARMS = ("noiseonly_indep", "total_indep")
ALL_ARMS = ARMS + INDEP_ARMS
MARGIN_DEG = 3.0


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--strategies", nargs="*", default=["drift", "raster"])
    ap.add_argument("--tag", default="",
                    help="suffix for the output filename; required for any "
                         "strategy set other than the published drift+raster")
    ap.add_argument("--pin-prior", action="store_true",
                    help="solve every strategy with the BAND's reference "
                         "sigma_p (ska_hi_prior: the three-strip drift's own "
                         "patch, per channel) instead of each strategy's own "
                         "patch variance. Without it a wider footprint picks "
                         "up more sky structure, is handed a looser prior and "
                         "is penalised for its extent -- which also breaks the "
                         "paper's stated condition that modes per unit sky is "
                         "comparable only at fixed sigma_p. The reference "
                         "strategy is solved at its own variance either way, "
                         "so its numbers must reproduce the unpinned run; that "
                         "is the control.")
    args = ap.parse_args()
    if not args.tag and sorted(args.strategies) != ["drift", "raster"]:
        sys.exit("refusing to overwrite the published ablation npz: pass --tag")
    if args.pin_prior:
        if not args.tag:
            sys.exit("--pin-prior changes every number: pass --tag")
        if prior.REF_STRATEGY not in args.strategies:
            sys.exit(f"--pin-prior needs '{prior.REF_STRATEGY}' in "
                     f"--strategies ({' '.join(args.strategies)}): it defines "
                     "the band reference, and including it also keeps the "
                     "common patch and the HI realisation comparable")
        prior.PINNED = True

    freqs = X.channel_freqs()
    cfg = X.band_config()
    dr = A.channel_dr_mpc(cfg)

    ops = {s: A.channel_operators(s) for s in args.strategies}
    common = A.common_patch([ops[s] for s in args.strategies])
    union = np.unique(np.concatenate(
        [np.asarray(m.pixel_indices) for s in ops for m in ops[s]]))
    interior = A.interior_mask(common, X.NSIDE, margin_deg=MARGIN_DEG)
    print(f"common patch {len(common)} px, interior {interior.sum()} px, "
          f"dr {dr:.2f} Mpc", flush=True)

    cube, box = build_hi_cube(cfg, union, verbose=False)
    hi_u = project_to_healpix(cube, box, cfg, union)
    look = {p: i for i, p in enumerate(union)}
    fg = {f: gdsm_equatorial_sky_model(freq=f, nside=X.NSIDE) for f in freqs}

    # Build the band reference now if it does not exist: it needs exactly the
    # operators already loaded here, so this costs nothing extra.
    if args.pin_prior:
        if prior.load() is None:
            prior.build(ops[prior.REF_STRATEGY])
        _ref = prior.load(required=True)
        _sp = np.sqrt(_ref["sigma_p2"])
        print(f"prior pinned to the band reference ({_ref['ref_strategy']}): "
              f"sigma_p {_sp.min():.4f}-{_sp.max():.4f} K", flush=True)

    out = dict(nmodes_grid=np.asarray(X.__dict__.get("NMODES_GRID",
                                                     (1, 2, 3, 4, 6, 8, 10))),
               common=common, interior=interior, margin_deg=MARGIN_DEG,
               arms=np.array(ALL_ARMS, dtype=object),
               pin_prior=bool(args.pin_prior),
               pin_prior_to=prior.REF_STRATEGY if args.pin_prior else "")
    for strategy in args.strategies:
        g, w = X.replay_noise(strategy)
        pass_seeds = [int(p["seed"]) for p in X.STRATEGIES[strategy]["pointings"]()]
        prior.set_strategy(strategy, ops[strategy])
        t0 = time.time()
        cubes = {a: [] for a in ALL_ARMS}
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
            # independent white, one draw per (pass, channel). SeedSequence
            # rather than hash(): reproducible, in range, and no collisions.
            w_ind = [np.random.default_rng(
                        np.random.SeedSequence([sd, i])
                     ).normal(0.0, np.sqrt(X.WHITE_VAR), size=len(w[j]))
                     for j, sd in enumerate(pass_seeds)]

            cubes["floor"].append(
                A.solve(mm, [c + h for c, h in zip(clean, hitod)],
                        truth, mu, freq_mhz=f)[sel])
            cubes["total"].append(
                A.solve(mm, [(c + h) * (1 + g[j]) * (1 + w[j])
                             for j, (c, h) in enumerate(zip(clean, hitod))],
                        truth, mu, freq_mhz=f)[sel])
            quiet = A.solve(mm, clean, truth, mu, freq_mhz=f)[sel]
            noisy = A.solve(mm, [c * (1 + g[j]) * (1 + w[j])
                                 for j, c in enumerate(clean)], truth, mu, freq_mhz=f)[sel]
            gainy = A.solve(mm, [c * (1 + g[j])
                                 for j, c in enumerate(clean)], truth, mu, freq_mhz=f)[sel]
            whitey = A.solve(mm, [c * (1 + w[j])
                                  for j, c in enumerate(clean)], truth, mu, freq_mhz=f)[sel]
            cubes["noiseonly"].append(noisy - quiet)
            cubes["gainonly"].append(gainy - quiet)
            cubes["whiteonly"].append(whitey - quiet)
            noisy_i = A.solve(mm, [c * (1 + g[j]) * (1 + w_ind[j])
                                   for j, c in enumerate(clean)],
                              truth, mu, freq_mhz=f)[sel]
            cubes["noiseonly_indep"].append(noisy_i - quiet)
            cubes["total_indep"].append(
                A.solve(mm, [(c + h) * (1 + g[j]) * (1 + w_ind[j])
                             for j, (c, h) in enumerate(zip(clean, hitod))],
                        truth, mu, freq_mhz=f)[sel])
            hi_mm.append(A.response(mm, hi, truth, freq_mhz=f)[sel])
            if (i + 1) % 8 == 0:
                print(f"  [{strategy}] {i + 1}/{len(freqs)} "
                      f"({time.time() - t0:.0f} s)", flush=True)

        hi_mm = np.asarray(hi_mm)
        cubes = {a: np.asarray(c) for a, c in cubes.items()}
        # The 1/f + white split is exact bar the g*w cross term; check it.
        resid = cubes["gainonly"] + cubes["whiteonly"] - cubes["noiseonly"]
        print(f"  [{strategy}] gain+white vs noiseonly: max residual "
              f"{np.abs(resid).max() / np.abs(cubes['noiseonly']).max():.2e} "
              f"of the noise amplitude", flush=True)

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
                for arm in ALL_ARMS:
                    cleaned = A.pca_clean(arm_cubes[arm][:, sel_pix], int(nm))
                    _, p = A.pk_par(cleaned, dr_mpc=dr)
                    out[f"{strategy}_{tag}{arm}_{nm}"] = p
                print(f"  [{strategy} {tag or 'as-run':6s}] nmodes={nm}: "
                      + ", ".join(
                          f"{arm} {np.sqrt(np.median(out[f'{strategy}_{tag}{arm}_{nm}'] / p_hi)):.1f}x"
                          for arm in ALL_ARMS), flush=True)

    os.makedirs(os.path.join(_HERE, "results"), exist_ok=True)
    _tag = f"_{args.tag}" if args.tag else ""
    path = os.path.join(_HERE, "results",
                        f"hi_ablation{_tag}_f{X.F_LO_MHZ:.0f}_{X.F_HI_MHZ:.0f}"
                        f"_nc{X.NCHAN}_ns{X.NSIDE}.npz")
    np.savez(path, **out)
    print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    main()
