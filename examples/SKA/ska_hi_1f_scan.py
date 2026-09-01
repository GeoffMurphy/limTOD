"""How bad would 1/f have to be before it mattered?

The ablation says the 1/f term sits far below the beam + prior floor -- but at
one parameter set, `GAIN_PARAMS = [w0, wc, alpha] = [1.335e-5, 1.099e-3, 2]`.
alpha = 2 is steeper than the alpha ~ 1 usually quoted for MeerKAT, and a steep
slope pushes power below the scan frequency, which is *favourable* to the
conclusion. So the honest version of "1/f is not the limiter" is a bound: how
much worse would the noise have to be before it reached the floor?

This scans the **knee frequency** -- where the 1/f gain power crosses the white
gain power -- at several slopes, and reports the cleaned 1/f residual against
the HI on the same axis as figure 3.

Two things make it cheap:

1. **The realisation is rescaled, not redrawn.** With `white_n_variance = 0` the
   correlation is exactly proportional to w0**alpha, so a draw made at w0 is a
   valid draw at w0' after multiplying by (w0'/w0)**(alpha/2). One draw per
   (alpha, pass) covers every knee, and holding the random numbers fixed along
   the x axis makes the curves smooth instead of realisation-noisy.
2. **The HI reference is reused** from `hi_experiment_*.npz`, so no HI cube has
   to be rebuilt -- this script never imports fastbox.

Writes ``results/hi_1f_scan_<band>.npz``.

    /home/geoff/gibbs_venv_312/bin/python ska_hi_1f_scan.py
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

# Ordinary frequency, mHz. The fiducial knee works out at 0.95 mHz; the grid
# runs to 100 mHz, ~100x worse, which is well beyond anything reported.
KNEES_MHZ = (0.95, 2.0, 5.0, 10.0, 30.0, 100.0)
ALPHAS = (1.5, 2.0, 2.5)
# Low-frequency cut-off, ordinary frequency in mHz. Fiducial is 0.175 mHz --
# a 1.6 h timescale, about one pass. Below that the model is adding correlated
# drift on timescales LONGER than a scan; above it, faster fluctuation.
# w_c does not enter the knee (that is set by w0 and alpha alone), so this
# axis extends the spectrum downward without moving the crossing with white.
CUTOFFS_MHZ = (0.0175, 0.055, 0.175, 0.55, 1.75, 5.5)
NMODES = 4
MARGIN_DEG = 3.0

W0_FID, WC_FID, ALPHA_FID = X.GAIN_PARAMS


def white_psd_level():
    """Flat PSD of the white gain term, in the convention of `flicker_corr`.

    `flicker_corr` returns C(tau) = (1/pi) \\int P(w) cos(w tau) dw, so a flat
    P = p integrated to the Nyquist angular frequency pi/dt gives C(0) = p/dt.
    The white draw has variance `WHITE_VAR` per sample, hence p = WHITE_VAR*dt.
    """
    return X.WHITE_VAR * X.DT


def omega0_for_knee(knee_mhz, alpha):
    """w0 placing the 1/f-white crossing at `knee_mhz` (ordinary frequency)."""
    w_knee = 2 * np.pi * knee_mhz * 1e-3
    return w_knee * white_psd_level() ** (1.0 / alpha)


def fiducial_knee_mhz(alpha=ALPHA_FID, w0=W0_FID):
    return w0 * white_psd_level() ** (-1.0 / alpha) / (2 * np.pi) * 1e3


def base_realisations(strategy, alpha):
    """One 1/f draw per pass at w0 = W0_FID, on the strategy's own seeds."""
    from limTOD.flicker_model import sim_noise
    out = []
    for p in X.STRATEGIES[strategy]["pointings"]():
        np.random.seed(p["seed"])
        out.append(sim_noise(W0_FID, WC_FID, alpha, p["tlist"],
                             n_samples=1, white_n_variance=0.0)[0])
    return out


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=("knee", "cutoff"), default="knee",
                    help="knee: scan the 1/f-white crossing (amplitude, an "
                         "exact rescaling). cutoff: scan w_c, which changes "
                         "the correlation SHAPE, so every point is a fresh "
                         "draw and nothing factors out.")
    args = ap.parse_args()

    freqs = X.channel_freqs()
    cfg = X.band_config()
    dr = A.channel_dr_mpc(cfg)
    ops = {s: A.channel_operators(s) for s in ("drift", "raster")}
    common = A.common_patch([ops[s] for s in ("drift", "raster")])
    interior = A.interior_mask(common, X.NSIDE, margin_deg=MARGIN_DEG)
    fg = {f: gdsm_equatorial_sky_model(freq=f, nside=X.NSIDE) for f in freqs}
    print(f"common {len(common)} px, interior {interior.sum()} px; "
          f"fiducial knee {fiducial_knee_mhz():.2f} mHz", flush=True)

    exp = np.load(os.path.join(_HERE, "results",
                               f"hi_experiment_f{X.F_LO_MHZ:.0f}_{X.F_HI_MHZ:.0f}"
                               f"_nc{X.NCHAN}_ns{X.NSIDE}.npz"))

    def cr_int(cube):
        return A.common_resolution(np.asarray(cube), freqs, common,
                                   X.NSIDE)[:, interior]

    out = dict(knees_mhz=np.asarray(KNEES_MHZ), alphas=np.asarray(ALPHAS),
               cutoffs_mhz=np.asarray(CUTOFFS_MHZ),
               fiducial_wc_mhz=WC_FID / (2 * np.pi) * 1e3,
               nmodes=NMODES, fiducial_knee_mhz=fiducial_knee_mhz(),
               fiducial_alpha=ALPHA_FID)

    for strategy in ("drift", "raster"):
        t0 = time.time()
        gain_fid, white_fid = X.replay_noise(strategy)

        # Per-channel clean TODs, the noiseless solve, and the truth/prior.
        clean_tods, quiet, truths, mus = [], [], [], []
        for i, f in enumerate(freqs):
            mm = ops[strategy][i]
            truth = fg[f][np.asarray(mm.pixel_indices)]
            mu = np.full_like(truth, float(np.mean(truth)))
            tods = X.simulate_foreground(strategy, f, verbose=False)
            clean = [np.asarray(t, float) / ((1 + gain_fid[j]) * (1 + white_fid[j]))
                     for j, t in enumerate(tods["TOD_group"])]
            clean_tods.append(clean)
            truths.append(truth)
            mus.append(mu)
            quiet.append(A.solve(mm, clean, truth, mu))
        print(f"  [{strategy}] baseline ready ({time.time() - t0:.0f} s)",
              flush=True)

        sel = np.isin(np.asarray(ops[strategy][0].pixel_indices), common)
        hi_ref = cr_int(exp[f"{strategy}_hi_mapmade"])
        k, p_hi = A.pk_par(hi_ref, dr_mpc=dr)
        out[f"{strategy}_k"] = k

        if args.mode == "cutoff":
            from limTOD.flicker_model import sim_noise
            for wc_mhz in CUTOFFS_MHZ:
                wc = 2 * np.pi * wc_mhz * 1e-3
                g = []
                for pt in X.STRATEGIES[strategy]["pointings"]():
                    np.random.seed(pt["seed"])
                    g.append(sim_noise(W0_FID, wc, ALPHA_FID, pt["tlist"],
                                       n_samples=1, white_n_variance=0.0)[0])
                rms = float(np.std(np.concatenate(g)))
                cube = []
                for i in range(len(freqs)):
                    mm = ops[strategy][i]
                    pix = np.asarray(mm.pixel_indices)
                    noisy = A.solve(mm, [c * (1 + g[j]) for j, c
                                         in enumerate(clean_tods[i])],
                                    truths[i], mus[i])
                    cube.append((noisy - quiet[i])[np.isin(pix, common)])
                cleaned = A.pca_clean(cr_int(np.asarray(cube)), NMODES)
                _, p = A.pk_par(cleaned, dr_mpc=dr)
                ratio = float(np.median(np.sqrt(p / p_hi)))
                out[f"{strategy}_wc{wc_mhz}"] = p
                out[f"{strategy}_ratio_wc{wc_mhz}"] = ratio
                out[f"{strategy}_rms_wc{wc_mhz}"] = rms
                print(f"  [{strategy}] w_c={wc_mhz:7.4f} mHz  gain rms "
                      f"{rms:.3e} -> {ratio:8.2f}x HI   "
                      f"({time.time() - t0:.0f} s)", flush=True)
            continue

        for alpha in ALPHAS:
            base = base_realisations(strategy, alpha)
            for knee in KNEES_MHZ:
                # Exact rescaling: C ~ w0**alpha, so the sample scales by
                # (w0'/w0)**(alpha/2) for the same random draw.
                fac = (omega0_for_knee(knee, alpha) / W0_FID) ** (alpha / 2.0)
                cube = []
                for i in range(len(freqs)):
                    mm = ops[strategy][i]
                    pix = np.asarray(mm.pixel_indices)
                    g = [fac * b for b in base]
                    noisy = A.solve(mm, [c * (1 + g[j]) for j, c
                                         in enumerate(clean_tods[i])],
                                    truths[i], mus[i])
                    cube.append((noisy - quiet[i])[np.isin(pix, common)])
                cleaned = A.pca_clean(cr_int(np.asarray(cube)), NMODES)
                _, p = A.pk_par(cleaned, dr_mpc=dr)
                ratio = float(np.median(np.sqrt(p / p_hi)))
                out[f"{strategy}_a{alpha}_k{knee}"] = p
                out[f"{strategy}_ratio_a{alpha}_k{knee}"] = ratio
                print(f"  [{strategy}] alpha={alpha} knee={knee:6.2f} mHz -> "
                      f"{ratio:8.2f}x HI   ({time.time() - t0:.0f} s)", flush=True)

    os.makedirs(os.path.join(_HERE, "results"), exist_ok=True)
    tag = "scan" if args.mode == "knee" else "cutoff"
    path = os.path.join(_HERE, "results",
                        f"hi_1f_{tag}_f{X.F_LO_MHZ:.0f}_{X.F_HI_MHZ:.0f}"
                        f"_nc{X.NCHAN}_ns{X.NSIDE}.npz")
    np.savez(path, **out)
    print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    main()
