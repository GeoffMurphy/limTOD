"""Experiment 006: HI signal recovery through foreground cleaning.

The queued test from ``ANALYSIS_LOG.md``. Everything measured in experiments
001-005 is reconstruction fidelity of the diffuse *foreground* sky -- there is
no 21 cm signal anywhere in those simulations and no foreground separation
step. This experiment injects one, cleans it, and measures how much HI the
pipeline destroys.

The concern that makes it urgent rather than generic: experiment 005 showed the
map's null space is strongly chromatic (17 modes at 350 MHz against 56 at
1050 MHz, patch shrinking 2.3x). Foreground cleaning assumes foregrounds are
spectrally smooth *after* the instrument and pipeline; a frequency-dependent
null space imprints spectral structure on the foreground residual that PCA
cannot remove, so it eats HI instead.

Configuration
-------------
350-400 MHz, 32 contiguous channels, nside 64 -- the series anchor, so the
established off-plane drift (317 px) vs raster (684 px) comparison carries
over unchanged. Note what this band can and cannot measure: at z ~ 2.8 a 15 m
dish resolves 401 Mpc transverse, so k_perp reaches only 0.016 Mpc^-1 while
k_par runs from 0.012 to 0.19. The accessible 3D k-space is a thin sliver
near the k_par axis, and the transfer function is reported against k_par for
that reason. This is a real property of the survey, not a limitation of the
simulation.

Design decisions worth knowing
------------------------------
* **Common patch.** Operator pixel selection is beam-dependent, so each
  channel selects a slightly different pixel set. Per-channel grids are
  themselves a chromatic effect and would confound the test, so everything is
  analysed on the intersection across all channels and both strategies.

* **Foregrounds are simulated; HI is injected through the operator.** The GDSM
  foreground TOD is generated properly by ``TODSim`` (full-sky, through the
  beam), which is what produces the beam + prior floor. The HI is injected as
  ``A s`` using the map-maker's own forward operator ``mm.Tsys_operators``.
  Measured against a real simulation the forward operator reproduces the sky
  TOD to 4-7% rms (corr 0.998), the gap being exactly the sub-beam and
  outside-patch structure it cannot represent. Injecting through ``A`` is what
  makes a mock-averaged transfer function affordable at all -- a solve costs
  0.57 s against minutes for a TOD simulation -- and the residual approximation
  largely cancels in a ratio. ``validate_injection`` quantifies it.

* **Linearity.** The map-maker is a Wiener filter, linear in the sky at fixed
  operator and noise weights, so an injected component adds linearly to the
  solved map. That is what lets mock injection be done as extra solves rather
  than extra simulations. The PCA step is *not* linear, which is precisely why
  the transfer function has to be measured by injection.

Run under ``gibbs_venv_312`` (has pyccl for fastbox, and reproduces the cached
TODs to 4e-14):

    /home/geoff/gibbs_venv_312/bin/python ska_hi_experiment.py
"""

from __future__ import annotations

import os
import pickle
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from ska_common import (SKA_LAT, SKA_LON, SKA_HGT, ska_beam_func,
                        ska_beam_fwhm_deg, gdsm_equatorial_sky_model,
                        drift_scan_night, constant_elevation_scan)
from ska_hi_mock import HIBandConfig

# ---------------------------------------------------------------------------
# Band and grid
# ---------------------------------------------------------------------------

F_LO_MHZ, F_HI_MHZ, NCHAN = 350.0, 400.0, 32
NSIDE = 64
NS_CELL = 16                        # ~3.7 deg cells, about one beam at 350 MHz

# Held identical to experiments 001-004.
DT = 2.0
SEED = 0
WHITE_VAR = 2.5e-6
TRUNCATE_FRAC_THRES = 1e-5
THRESHOLD = 0.05
GAIN_PARAMS = [1.335e-5, 1.099e-3, 2]

# --- drift (experiments 001/002, off-plane) --------------------------------
DRIFT_AZ = 0.0
DRIFT_ELEVATIONS = [52.0, 50.0, 48.0]
DRIFT_NIGHT_S = 3600.0
DRIFT_UTC = "2024-04-15 19:00:00"

# --- raster (experiment 004, off-plane) ------------------------------------
# BASE_UTC recovered by matching the cached LST of
# simulated_TODs_ska_meerklass_offplane_gauss.npz exactly (0.0000 deg).
RASTER_EL = 38.19
RASTER_AZ_HALFWIDTH = 7.5
RASTER_SWEEP_S = 150.0
RASTER_PASS_S = 5400.0
RASTER_PASSES = [(45.0, 16.401), (315.0, 20.961)]   # az centre, start hour
RASTER_UTC = "2024-04-16 00:00:00"

CACHE_DIR = os.path.join(_HERE, "hi_cache")


def channel_freqs():
    edges = np.linspace(F_LO_MHZ, F_HI_MHZ, NCHAN + 1)
    return 0.5 * (edges[:-1] + edges[1:])


def band_config(seed=SEED):
    return HIBandConfig(f_lo_mhz=F_LO_MHZ, f_hi_mhz=F_HI_MHZ, nchan=NCHAN,
                        nside=NSIDE, seed=seed)


# ---------------------------------------------------------------------------
# Scan geometry
# ---------------------------------------------------------------------------

def drift_pointings():
    """(time, azimuth, elevation) per night for the off-plane drift."""
    out = []
    for night, el in enumerate(DRIFT_ELEVATIONS):
        tlist, azlist = drift_scan_night(DRIFT_NIGHT_S, dt=DT, night=night,
                                         azimuth_deg=DRIFT_AZ)
        out.append(dict(tlist=tlist, azlist=azlist, el=el, seed=SEED + night))
    return out


def raster_pointings():
    """(time, azimuth, elevation) per pass for the off-plane raster."""
    out = []
    for i, (az_centre, start_h) in enumerate(RASTER_PASSES):
        tlist, azlist, _ = constant_elevation_scan(
            RASTER_PASS_S, az_centre, RASTER_AZ_HALFWIDTH, RASTER_SWEEP_S,
            dt=DT, t0_s=start_h * 3600.0)
        out.append(dict(tlist=tlist, azlist=azlist, el=RASTER_EL,
                        seed=SEED + i))
    return out


STRATEGIES = {
    "drift": dict(pointings=drift_pointings, utc=DRIFT_UTC),
    "raster": dict(pointings=raster_pointings, utc=RASTER_UTC),
}


# ---------------------------------------------------------------------------
# Noise replay
# ---------------------------------------------------------------------------

def replay_noise(strategy):
    """Regenerate the seeded 1/f gain and white draws for one strategy.

    ``generate_TOD`` draws 1/f first via ``sim_noise`` and then white via
    ``np.random.normal`` from a stream seeded per night/pass, so replaying in
    the same order reproduces the realisations inside the simulated TODs. The
    noise is multiplicative: TOD = (A s) (1 + g) (1 + w).
    """
    from limTOD.flicker_model import sim_noise

    gain, white = [], []
    for p in STRATEGIES[strategy]["pointings"]():
        np.random.seed(p["seed"])
        gain.append(sim_noise(*GAIN_PARAMS, p["tlist"], n_samples=1,
                              white_n_variance=0.0)[0])
        white.append(np.random.normal(0, np.sqrt(WHITE_VAR),
                                      size=(1, len(p["tlist"])))[0])
    return gain, white


# ---------------------------------------------------------------------------
# The expensive, cached steps
# ---------------------------------------------------------------------------

def _cache(name):
    os.makedirs(CACHE_DIR, exist_ok=True)
    return os.path.join(CACHE_DIR, name)


def simulate_foreground(strategy, freq_mhz, verbose=True):
    """Simulate (or load) the GDSM foreground TOD for one channel."""
    from limTOD import TODSim

    path = _cache(f"tod_{strategy}_f{freq_mhz:07.3f}_ns{NSIDE}.npz")
    try:
        d = np.load(path, allow_pickle=True)
        return {k: [np.asarray(x, np.float64) for x in d[k]]
                for k in ("TOD_group", "LST_deg_list_group",
                          "azimuth_deg_list_group", "elevation_deg_list_group")}
    except (FileNotFoundError, KeyError):
        pass

    if verbose:
        print(f"[{strategy} {freq_mhz:.3f} MHz] simulating foreground TOD...",
              flush=True)
    sim = TODSim(ant_latitude_deg=SKA_LAT, ant_longitude_deg=SKA_LON,
                 ant_height_m=SKA_HGT, beam_func=ska_beam_func,
                 sky_func=gdsm_equatorial_sky_model,
                 beam_nside=NSIDE, sky_nside=NSIDE)
    groups = dict(TOD_group=[], LST_deg_list_group=[],
                  azimuth_deg_list_group=[], elevation_deg_list_group=[])
    for p in STRATEGIES[strategy]["pointings"]():
        np.random.seed(p["seed"])
        tod, _, _, lst = sim.generate_TOD(
            freq_list=[freq_mhz], time_list=p["tlist"],
            azimuth_deg_list=p["azlist"], elevation_deg=p["el"],
            start_time_utc=STRATEGIES[strategy]["utc"],
            white_noise_var=WHITE_VAR, return_LSTs=True,
            normalize_beam=False, truncate_frac_thres=TRUNCATE_FRAC_THRES)
        groups["TOD_group"].append(np.asarray(tod[0], np.float64))
        groups["LST_deg_list_group"].append(np.asarray(lst, np.float64))
        groups["azimuth_deg_list_group"].append(
            np.asarray(p["azlist"], np.float64))
        groups["elevation_deg_list_group"].append(
            p["el"] * np.ones(len(p["tlist"]), np.float64))
    np.savez(path, **{k: np.asarray(v) for k, v in groups.items()})
    return groups


def build_operator(strategy, freq_mhz, tod_groups=None, verbose=True):
    """Build (or load) the HPW map-making operator for one channel.

    ``nside_target`` is deliberately equal to the beam nside: ``HPW_mapmaking``
    hardcodes ``normalize_beam=False`` and mis-scales every operator row by
    ``(nside_target/nside_beam)**2`` when they differ (HANDOFF.md item 3).
    """
    from limTOD import HPW_mapmaking

    path = _cache(f"op_{strategy}_f{freq_mhz:07.3f}_ns{NSIDE}.pkl")
    try:
        with open(path, "rb") as f:
            return pickle.load(f)
    except (FileNotFoundError, EOFError):
        pass

    if tod_groups is None:
        tod_groups = simulate_foreground(strategy, freq_mhz, verbose=verbose)
    if verbose:
        print(f"[{strategy} {freq_mhz:.3f} MHz] building operator "
              f"(FWHM {ska_beam_fwhm_deg(freq_mhz):.2f} deg)...", flush=True)
    mm = HPW_mapmaking(
        beam_map=ska_beam_func(freq=freq_mhz, nside=NSIDE),
        LST_deg_list_group=tod_groups["LST_deg_list_group"],
        lat_deg=SKA_LAT,
        azimuth_deg_list_group=tod_groups["azimuth_deg_list_group"],
        elevation_deg_list_group=tod_groups["elevation_deg_list_group"],
        threshold=THRESHOLD, nside_target=NSIDE,
        beam_truncate_frac_thres=TRUNCATE_FRAC_THRES)
    with open(path, "wb") as f:
        pickle.dump(mm, f)
    if verbose:
        print(f"[{strategy} {freq_mhz:.3f} MHz] operator built, "
              f"{len(mm.pixel_indices)} pixels", flush=True)
    return mm


def prepare_all(strategies=("drift", "raster"), verbose=True):
    """Warm every operator and TOD cache. This is the multi-hour step."""
    for strategy in strategies:
        for f in channel_freqs():
            tod = simulate_foreground(strategy, f, verbose=verbose)
            build_operator(strategy, f, tod, verbose=verbose)
    print("all caches warm", flush=True)


if __name__ == "__main__":
    which = sys.argv[1:] or ["drift", "raster"]
    prepare_all(tuple(which))
