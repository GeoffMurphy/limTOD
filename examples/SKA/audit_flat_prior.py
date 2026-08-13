"""Does the cross-linking result survive a flat prior?

Re-solves the off-plane drift-vs-raster comparison (experiment 004) under two
prior means, changing nothing else: the series' beam-smoothed-truth prior, and a
flat prior constant at the patch mean. Operators and TODs come from cache, so
this is a re-solve, not a re-simulation (~1 min).

Motivation is in ANALYSIS_LOG.md: the beam-smoothed-truth prior turned out to
supply most of the apparent frequency trend in experiment 005, and every
residual in HANDOFF.md uses that prior — including the cross-linking headline.
A comparison at fixed frequency should not be confounded the same way, since
the prior does not move between the two arms, but that was argued rather than
demonstrated. This demonstrates it.

Run from examples/SKA with the limTOD venv:

    ../../.venv/bin/python audit_flat_prior.py

Needs the cached operators and TODs for the off-plane drift (ska_drift_gauss)
and the off-plane raster (ska_meerklass_offplane_gauss).
"""
import pickle

import numpy as np
import healpy as hp

from limTOD.flicker_model import sim_noise
from ska_common import (constant_elevation_scan, drift_scan_night,
                        gdsm_equatorial_sky_model, ska_beam_fwhm_deg)

# ---- config, identical to ska_results_summary.ipynb ----
FREQ, NS, DT, WHITE_VAR = 350, 64, 2.0, 2.5e-6
NS_CELL = 16                       # ~3.7 deg cells, about one beam
GAIN_PARAMS = [1.335e-5, 1.099e-3, 2]
SEED = 0
FWHM = ska_beam_fwhm_deg(FREQ)

RASTER_STEM = "ska_meerklass_offplane_gauss"
AZ_HALFWIDTH, SWEEP_S, PASS_S = 7.5, 150.0, 5400.0
RASTER_PASSES = [(45.0, 16.401), (315.0, 20.961)]   # az centre, start hour

sky_truth_full = gdsm_equatorial_sky_model(freq=FREQ, nside=NS)
prior_full = hp.smoothing(sky_truth_full, fwhm=np.radians(FWHM))

# ---- noise replay (seeded; reproduces the realisations in the cached TODs) ----
d_tlists = [drift_scan_night(3600.0, dt=DT, night=n)[0] for n in range(3)]
gain_noise, white_noise = [], []
for n, tl in enumerate(d_tlists):
    np.random.seed(SEED + n)
    gain_noise.append(sim_noise(*GAIN_PARAMS, tl, n_samples=1,
                                white_n_variance=0.0)[0])
    white_noise.append(np.random.normal(0, np.sqrt(WHITE_VAR),
                                        size=(1, len(tl)))[0])

r_tlists = [constant_elevation_scan(PASS_S, az, AZ_HALFWIDTH, SWEEP_S,
                                    dt=DT, t0_s=h * 3600.0)[0]
            for az, h in RASTER_PASSES]
r_gain, r_white = [], []
for i, tl in enumerate(r_tlists):
    np.random.seed(SEED + i)
    r_gain.append(sim_noise(*GAIN_PARAMS, tl, n_samples=1,
                            white_n_variance=0.0)[0])
    r_white.append(np.random.normal(0, np.sqrt(WHITE_VAR),
                                    size=(1, len(tl)))[0])

# ---- caches ----
with open(f"mapmaker_ops_{RASTER_STEM}_ns{NS}.pkl", "rb") as f:
    r_mm = pickle.load(f)
_rd = np.load(f"simulated_TODs_{RASTER_STEM}.npz", allow_pickle=True)
r_sky = [np.asarray(_rd["TOD_group"][i], np.float64)
         / ((1 + r_gain[i]) * (1 + r_white[i])) for i in range(len(r_tlists))]

d_mm = pickle.load(open(f"mapmaker_ops_ska_drift_gauss_ns{NS}.pkl", "rb"))
d_sky = [np.asarray(np.load("simulated_TODs_ska_drift_gauss.npz",
                            allow_pickle=True)["TOD_group"][i], np.float64)
         / ((1 + gain_noise[i]) * (1 + white_noise[i])) for i in range(3)]

common = np.intersect1d(d_mm.pixel_indices, r_mm.pixel_indices)
depth_ratio = len(r_mm.pixel_indices) / len(d_mm.pixel_indices)
print(f"drift {len(d_mm.pixel_indices)} px, raster {len(r_mm.pixel_indices)} px, "
      f"{len(common)} shared -> raster {depth_ratio:.2f}x shallower")


def solve(mm, TOD_group, prior_kind):
    """The experiment notebooks' recipe, no high-pass.

    Only the prior MEAN differs between the two variants. The prior covariance,
    the noise model, the regularisation and the operator are all held fixed, so
    any difference is attributable to the prior mean alone.
    """
    pix = mm.pixel_indices
    truth = sky_truth_full[pix]
    if prior_kind == "smoothed":
        mu = prior_full[pix]
    elif prior_kind == "flat":
        # Carries the level but no structure, so every spatial mode in the
        # answer has to come from the data.
        mu = np.full_like(truth, float(np.mean(truth)))
    else:
        raise ValueError(prior_kind)
    nv_floor = WHITE_VAR * (1e-3 * float(np.mean(truth)))**2
    return mm(TOD_group=TOD_group, dtime=DT, Tsky_prior_mean=mu,
              Tsky_prior_inv_cov_diag=np.ones_like(truth)
              / max(float(np.std(truth)), 1e-3)**2,
              noise_variance=[WHITE_VAR * (np.asarray(o) @ truth)**2 + nv_floor
                              for o in mm.Tsys_operators],
              regularization=1e-12, return_full_cov=False)[0]


def cell_rms(pix, resid, hit_pix):
    """Beam-scale RMS, restricted to cells lying entirely inside `hit_pix` so
    both strategies are compared on identical cells."""
    full = np.zeros(hp.nside2npix(NS))
    hit = np.zeros(hp.nside2npix(NS), bool)
    full[pix], hit[hit_pix] = resid, True
    order = hp.ring2nest(NS, np.arange(hp.nside2npix(NS)))
    fn, hn = np.zeros_like(full), np.zeros_like(hit)
    fn[order], hn[order] = full, hit
    k = (NS // NS_CELL)**2
    fc, hc = fn.reshape(-1, k), hn.reshape(-1, k)
    keep = hc.all(axis=1)
    return float(np.std(fc[keep].mean(axis=1)))


sky_struct = float(np.std(sky_truth_full[common]))

results = {}
for prior_kind in ("smoothed", "flat"):
    for name, mm, sky, g, w in (("drift", d_mm, d_sky, gain_noise, white_noise),
                                ("raster", r_mm, r_sky, r_gain, r_white)):
        sel = np.isin(mm.pixel_indices, common)
        rng = range(len(sky))
        sc = 1 / np.sqrt(depth_ratio)
        variants = {
            "total": [sky[i] * (1 + g[i]) * (1 + w[i]) for i in rng],
            "floor": [sky[i] for i in rng],
        }
        if name == "raster":                   # same depth, not same clock time
            variants["total @ depth"] = [
                sky[i] * (1 + g[i] * sc) * (1 + w[i] * sc) for i in rng]
        for key, tg in variants.items():
            resid = solve(mm, tg, prior_kind) - sky_truth_full[mm.pixel_indices]
            results[(prior_kind, name, key)] = dict(
                pix=float(np.std(resid[sel])),
                beam=cell_rms(mm.pixel_indices, resid, common))

# ---- report ----
for prior_kind in ("smoothed", "flat"):
    tag = ("beam-smoothed-truth prior (the series' recipe)"
           if prior_kind == "smoothed" else "FLAT prior (constant at patch mean)")
    print(f"\n=== {tag} ===")
    print(f"{'strategy':10s} {'quantity':15s} {'per-pixel':>10s} {'beam-scale':>11s}")
    for name in ("drift", "raster"):
        for key in ("total", "floor", "total @ depth"):
            if (prior_kind, name, key) not in results:
                continue
            v = results[(prior_kind, name, key)]
            print(f"{name:10s} {key:15s} {v['pix']:9.3f}K {v['beam']:10.3f}K")

print("\n\n=== does cross-linking survive the prior change? ===")
print(f"{'':31s} {'per-pixel':>22s} {'beam-scale':>22s}")
for prior_kind in ("smoothed", "flat"):
    for label, dkey, rkey in (("floor", "floor", "floor"),
                              ("total @ matched depth", "total", "total @ depth")):
        out = []
        for scname in ("pix", "beam"):
            d = results[(prior_kind, "drift", dkey)][scname]
            r = results[(prior_kind, "raster", rkey)][scname]
            out.append(f"{d:.3f}->{r:.3f} ({(r - d) / d:+6.1%})")
        print(f"{prior_kind:8s} {label:22s} {out[0]:>22s} {out[1]:>22s}")

print("\n=== how much does each strategy lean on the prior? ===")
print("(flat / smoothed on the same quantity: 1.00 = the prior was doing nothing)")
for scname, lbl in (("pix", "per-pixel"), ("beam", "beam-scale")):
    for name in ("drift", "raster"):
        s = results[("smoothed", name, "floor")][scname]
        f = results[("flat", name, "floor")][scname]
        print(f"  {lbl:10s} {name:7s} floor {s:.3f} -> {f:.3f} K  ({f / s:.2f}x)")

print(f"\nsky structure on shared pixels: {sky_struct:.3f} K")
print("residual / sky structure, floor (prior-independent framing):")
for prior_kind in ("smoothed", "flat"):
    d = results[(prior_kind, "drift", "floor")]["pix"] / sky_struct
    r = results[(prior_kind, "raster", "floor")]["pix"] / sky_struct
    print(f"  {prior_kind:8s} drift {d:.3f}, raster {r:.3f}")
