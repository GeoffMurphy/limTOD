"""Configuration and cached heavy steps for the off-plane frequency sweep
(experiment 005).

Experiments 001-004 all sat at 350 MHz, the bottom edge of SKA-Mid Band 1.
This module repeats the *off-plane* drift geometry of experiment 001
(``ska_drift_scan.ipynb``) at five equally spaced channels spanning the whole
of Band 1, holding everything except the frequency fixed: same site, same
parked azimuth, same three elevations on three sidereal nights, same 3 h of
2 s samples, same noise model and seeds, same map-making recipe.

Two things change on their own as the frequency rises, and separating them is
the point of the experiment:

1. **The beam narrows** as 1.22 lambda/D, from 3.99 deg at 350 MHz to 1.33 deg
   at 1050 MHz. The elevation ladder (52/50/48 deg, i.e. 2 deg Dec steps) was
   chosen as *half-beam* spacing at 350 MHz, so above ~600 MHz the three Dec
   strips stop overlapping and the patch breaks into disjoint ribbons.
2. **The sky signal's temporal band does *not* move — a prediction this sweep
   tested and refuted.** The tempting argument: a drift crosses the beam in
   FWHM / (15 deg/hr), 957 s at 350 MHz but 319 s at 1050 MHz, so the sky band
   should climb from ~1.0 to ~3.1 mHz and cross the 2 mHz high-pass cutoff,
   flipping the filter's verdict. **Measured, it does not.** The sky TOD's
   power sits at the *pass fundamental* (0.28 mHz, one 1 h traverse) at every
   channel: median-power frequency 0.28 mHz and ~93% of power below 2 mHz,
   identical from 350 to 1050 MHz. Diffuse GDSM emission is red, so the
   large-scale gradient along the drift dominates the spectrum no matter how
   narrow the beam is; the beam crossing time sets the *upper edge* of the band
   the beam can transfer, not where the power lives.

The whole sweep is solved on a **single nside 128 grid** so map-space residuals
are directly comparable across channels (per-pixel RMS is not comparable across
different grids — see HANDOFF.md). nside 128 keeps at least 2.9 pixels across
the FWHM everywhere, and the top channel is re-run at nside 256 as a check that
the grid is not what limits it.

Gaussian beam only: the tapered-aperture sidelobe question is a bright-field
result, and this sweep is off-plane (b ~ +50 deg), following the same choice
made for the off-plane raster in ``ska_results_summary.ipynb``.
"""

from __future__ import annotations

import os
import pickle

import numpy as np

from ska_common import (SKA_LAT, SKA_LON, SKA_HGT, ska_beam_func,
                        ska_beam_fwhm_deg, gdsm_equatorial_sky_model,
                        drift_scan_night)

# ---------------------------------------------------------------------------
# Sweep configuration
# ---------------------------------------------------------------------------

# Five equally spaced channels across SKA-Mid Band 1 (350-1050 MHz). The first
# is the frequency every earlier experiment in this series used.
CHANNELS_MHZ = [350, 525, 700, 875, 1050]

# One grid for the whole sweep, so residuals compare like for like.
NSIDE_SWEEP = 128

# Grid check on the narrowest beam, where nside 128 gives only 2.9 px/FWHM.
NSIDE_CHECK = 256
CHECK_CHANNEL_MHZ = 1050

# --- held fixed, identical to experiment 001 (off-plane drift) -------------
DRIFT_AZ = 0.0                     # parked on the meridian, facing north
ELEVATIONS = [52.0, 50.0, 48.0]    # one per night; 2 deg Dec steps
NIGHT_DURATION_S = 3600.0
DT = 2.0
START_UTC = "2024-04-15 19:00:00"  # LST ~ 151 deg at Karoo -> b ~ +50 deg
SEED = 0
WHITE_VAR = 2.5e-6                 # fractional white-noise variance
TRUNCATE_FRAC_THRES = 1e-5         # -50 dB beam truncation
THRESHOLD = 0.05                   # beam-response pixel selection
HP_CUTOFF = 2e-3                   # Hz, the series' standard high-pass

# Flicker-noise parameters of limTOD's default gain model, needed to replay
# the seeded draws when decomposing the residual (see ``replay_noise``).
GAIN_PARAMS = [1.335e-5, 1.099e-3, 2]

SIDEREAL_RATE_DEG_PER_HR = 360.0 / 23.9344696


def nights():
    """The three drift nights: time and azimuth lists, one dict per night."""
    out = []
    for night, el in enumerate(ELEVATIONS):
        tlist, azlist = drift_scan_night(NIGHT_DURATION_S, dt=DT, night=night,
                                         azimuth_deg=DRIFT_AZ)
        out.append(dict(night=night, el=el, tlist=tlist, azlist=azlist))
    return out


def beam_crossing_freq_hz(freq_mhz):
    """Upper edge of the temporal band a parked dish can transfer.

    The sky sweeps past at the sidereal rate, so a point source takes
    ``FWHM / 15 deg/hr`` to cross the beam; the reciprocal is the highest
    temporal frequency at which the beam still passes sky structure, since
    anything finer is smoothed away.

    **This is not where the sky power lies, and it is the wrong number to
    compare against ``HP_CUTOFF``.** Measured on the cached TODs, the sky
    spectrum is essentially frequency-independent: power concentrates at the
    pass fundamental (0.28 mHz for a 1 h traverse) and ~93% of it sits below
    2 mHz at every channel from 350 to 1050 MHz, because diffuse emission is
    red and the large-scale gradient along the drift dominates. Use
    ``sky_band_percentiles`` for the question the high-pass actually cares
    about.
    """
    crossing_s = (ska_beam_fwhm_deg(freq_mhz) / SIDEREAL_RATE_DEG_PER_HR
                  * 3600.0)
    return 1.0 / crossing_s


def sky_band_percentiles(sky_tod_group, percentiles=(0.5, 0.9, 0.99)):
    """Where the sky TOD's power actually sits, in Hz.

    Returns ``(freqs, cumulative_power_fraction, {p: freq_at_p})`` for the
    mean periodogram over the group. This is the measured answer to "would a
    high-pass at f_c remove the sky signal?", and it is what
    ``beam_crossing_freq_hz`` is *not*.
    """
    from scipy import signal

    fr, P = signal.periodogram(np.asarray(sky_tod_group), fs=1.0 / DT,
                               detrend="constant", axis=-1)
    P = P.mean(0)
    cum = np.cumsum(P) / P.sum()
    at = {p: float(fr[np.searchsorted(cum, p)]) for p in percentiles}
    return fr, cum, at


def px_per_fwhm(freq_mhz, nside):
    """Pixels across the beam FWHM — the sampling figure of merit.

    HANDOFF.md records ~3 px/FWHM as the floor: at 2.2 the reconstruction
    degrades badly and even the high-pass verdict flips sign.
    """
    import healpy as hp
    return ska_beam_fwhm_deg(freq_mhz) / np.degrees(
        hp.nside2resol(nside))


# ---------------------------------------------------------------------------
# Cache paths
# ---------------------------------------------------------------------------

_HERE = os.path.dirname(os.path.abspath(__file__))


def tod_cache_path(freq_mhz, nside):
    """TOD cache. The nside is in the name because the simulation itself is
    done at that resolution (``sky_nside`` = ``beam_nside`` = nside), so the
    nside-256 check is a genuinely different simulation, not a re-solve."""
    return os.path.join(
        _HERE, f"simulated_TODs_ska_freqsweep_gauss_f{freq_mhz:04d}"
               f"_ns{nside}.npz")


def op_cache_path(freq_mhz, nside):
    return os.path.join(
        _HERE, f"mapmaker_ops_ska_freqsweep_gauss_f{freq_mhz:04d}"
               f"_ns{nside}.pkl")


# ---------------------------------------------------------------------------
# The two expensive steps, both cached
# ---------------------------------------------------------------------------

def simulate_channel(freq_mhz, nside, verbose=True):
    """Simulate (or load) the three nights of TOD for one channel.

    Seeded ``SEED + night``, exactly as experiments 001-004, so the 1/f and
    white realizations can be replayed later by ``replay_noise``.
    """
    from limTOD import TODSim

    cache = tod_cache_path(freq_mhz, nside)
    try:
        data = np.load(cache, allow_pickle=True)
        groups = {k: [np.asarray(x, np.float64) for x in data[k]]
                  for k in ("TOD_group", "LST_deg_list_group",
                            "azimuth_deg_list_group",
                            "elevation_deg_list_group")}
        if verbose:
            print(f"[{freq_mhz} MHz ns{nside}] loaded cached TODs")
        return groups
    except (FileNotFoundError, KeyError):
        pass

    if verbose:
        print(f"[{freq_mhz} MHz ns{nside}] simulating TODs "
              f"(FWHM {ska_beam_fwhm_deg(freq_mhz):.2f} deg, "
              f"{px_per_fwhm(freq_mhz, nside):.1f} px/FWHM)...", flush=True)

    tod_sim = TODSim(
        ant_latitude_deg=SKA_LAT,
        ant_longitude_deg=SKA_LON,
        ant_height_m=SKA_HGT,
        beam_func=ska_beam_func,
        sky_func=gdsm_equatorial_sky_model,
        beam_nside=nside,
        sky_nside=nside,
    )
    groups = dict(TOD_group=[], LST_deg_list_group=[],
                  azimuth_deg_list_group=[], elevation_deg_list_group=[])
    for n in nights():
        np.random.seed(SEED + n["night"])
        tod_array, _, _, lst_deg = tod_sim.generate_TOD(
            freq_list=[freq_mhz],
            time_list=n["tlist"],
            azimuth_deg_list=n["azlist"],
            elevation_deg=n["el"],
            start_time_utc=START_UTC,
            white_noise_var=WHITE_VAR,
            return_LSTs=True,
            normalize_beam=False,
            truncate_frac_thres=TRUNCATE_FRAC_THRES,
        )
        groups["TOD_group"].append(np.asarray(tod_array[0], np.float64))
        groups["LST_deg_list_group"].append(np.asarray(lst_deg, np.float64))
        groups["azimuth_deg_list_group"].append(
            np.asarray(n["azlist"], np.float64))
        groups["elevation_deg_list_group"].append(
            n["el"] * np.ones_like(groups["TOD_group"][-1]))

    np.savez(cache, freq_list=np.asarray([freq_mhz]), **groups)
    if verbose:
        print(f"[{freq_mhz} MHz ns{nside}] TODs saved", flush=True)
    return groups


def build_operator(freq_mhz, nside, tod_groups=None, verbose=True):
    """Build (or load) the HPW map-making operator for one channel.

    ``nside_target`` is deliberately equal to the beam nside: HPW_mapmaking
    hardcodes ``normalize_beam=False`` and mis-scales every operator row by
    ``(nside_target/nside_beam)**2`` when they differ (HANDOFF.md item 3).
    """
    from limTOD import HPW_mapmaking

    cache = op_cache_path(freq_mhz, nside)
    try:
        with open(cache, "rb") as f:
            mm = pickle.load(f)
        if verbose:
            print(f"[{freq_mhz} MHz ns{nside}] loaded cached operator "
                  f"({len(mm.pixel_indices)} pixels)")
        return mm
    except (FileNotFoundError, EOFError):
        pass

    if tod_groups is None:
        tod_groups = simulate_channel(freq_mhz, nside, verbose=verbose)
    if verbose:
        print(f"[{freq_mhz} MHz ns{nside}] building operator...", flush=True)

    mm = HPW_mapmaking(
        beam_map=ska_beam_func(freq=freq_mhz, nside=nside),
        LST_deg_list_group=tod_groups["LST_deg_list_group"],
        lat_deg=SKA_LAT,
        azimuth_deg_list_group=tod_groups["azimuth_deg_list_group"],
        elevation_deg_list_group=tod_groups["elevation_deg_list_group"],
        threshold=THRESHOLD,
        nside_target=nside,
        beam_truncate_frac_thres=TRUNCATE_FRAC_THRES,
    )
    with open(cache, "wb") as f:
        pickle.dump(mm, f)
    if verbose:
        print(f"[{freq_mhz} MHz ns{nside}] operator built, "
              f"{len(mm.pixel_indices)} pixels", flush=True)
    return mm


def replay_noise(n_nights=3):
    """Regenerate the exact 1/f and white draws ``generate_TOD`` made.

    ``generate_TOD`` draws 1/f first via ``sim_noise`` and then white via
    ``np.random.normal``, from a stream seeded ``SEED + night``. Replaying in
    the same order reproduces the realizations inside the cached TODs (verified
    to 4e-14 in earlier experiments), which is what lets the residual be split
    into beam+prior floor, white and 1/f terms. The draws depend only on the
    time grid and seed, so one replay serves every channel.
    """
    from limTOD.flicker_model import sim_noise

    gain, white = [], []
    for n in nights()[:n_nights]:
        np.random.seed(SEED + n["night"])
        gain.append(sim_noise(*GAIN_PARAMS, n["tlist"], n_samples=1,
                              white_n_variance=0.0)[0])
        white.append(np.random.normal(0, np.sqrt(WHITE_VAR),
                                      size=(1, len(n["tlist"])))[0])
    return gain, white


def prepare_all(channels=None, nside=NSIDE_SWEEP, with_check=True):
    """Warm every cache the notebook needs. Safe to re-run; skips what exists.

    This is the multi-hour step — run it from the command line
    (``python ska_freq_sweep.py``) rather than waiting on it in the notebook.
    """
    jobs = [(f, nside) for f in (channels or CHANNELS_MHZ)]
    if with_check:
        jobs.append((CHECK_CHANNEL_MHZ, NSIDE_CHECK))
    for freq, ns in jobs:
        tod = simulate_channel(freq, ns)
        build_operator(freq, ns, tod)
    print("all caches warm")


if __name__ == "__main__":
    prepare_all()
