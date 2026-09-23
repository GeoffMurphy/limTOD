"""Experiment 006 analysis: PCA cleaning and the HI transfer function.

Consumes the operators and foreground TODs warmed by ``ska_hi_experiment.py``
and produces the signal-recovery numbers: how much HI survives the map-maker,
how much more the foreground clean destroys, and what the transfer function
correction puts back.

Four power spectra tell the story, all along the line of sight (see the module
docstring of ``ska_hi_experiment`` for why k_par is the only well-sampled
direction in this band):

1. ``P_true``      -- the injected HI cube, before anything touches it.
2. ``P_mapmade``   -- HI after the map-maker, ``W A s``. The gap to ``P_true``
                      is the beam plus null-space loss this series has been
                      measuring all along, now expressed on the signal rather
                      than on the foreground residual.
3. ``P_clean``     -- HI-attributable power surviving PCA foreground removal.
4. ``P_corrected`` -- ``P_clean / T(k)``. Should recover ``P_mapmade``, not
                      ``P_true``: the transfer function corrects for the
                      cleaning, not for the map-maker.

The transfer function follows Cunnington et al. (2023) [arXiv:2302.07034],
which is also what ``fastbox.filters.pca_transfer_function`` implements:

    T(k) = P(X_m_clean - X_clean, X_m) / P(X_m, X_m)

with the cross-power in the numerator to avoid the positive noise bias an
auto-power of the difference would carry. It is estimated here rather than
called from fastbox because fastbox's version bins in 3D |k| on a Cartesian
box, and our data lives on an irregular HEALPix patch where only the radial
direction is well sampled.

Prior choice: **flat**, constant at the patch mean. The 2026-08-13 audit
concluded that beam-smoothed-truth prior numbers flatter the drift (it leans on
the prior 1.40x against the raster's 1.05x) and that flat-prior numbers are the
ones to quote when the claim is about what the survey measures. A signal
recovery claim is exactly that.
"""

from __future__ import annotations

import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import ska_hi_experiment as X
from ska_common import gdsm_equatorial_sky_model
from ska_hi_mock import (HIBandConfig, build_hi_cube, project_to_healpix,
                         comoving_distances)


# ---------------------------------------------------------------------------
# Patch bookkeeping
# ---------------------------------------------------------------------------

def channel_operators(strategy, verbose=False):
    """All per-channel operators for one strategy, in channel order."""
    return [X.build_operator(strategy, f, verbose=verbose)
            for f in X.channel_freqs()]


def interior_mask(pixel_indices, nside, margin_deg=3.0):
    """Pixels whose `margin_deg` neighbourhood lies entirely inside the patch.

    ``common_resolution`` smooths a patch embedded in a zero-filled sphere, so
    the map is pulled toward zero within roughly a kernel width of the
    boundary. Restricting to the interior is what keeps that out of the
    numbers. It is not free -- at the default margin it costs well over half
    the 277-pixel patch -- so any comparison against un-smoothed cubes must be
    made on this same subset, not on the full patch.
    """
    import healpy as hp

    pix = np.asarray(pixel_indices)
    inside = np.zeros(hp.nside2npix(nside), bool)
    inside[pix] = True
    keep = np.empty(len(pix), bool)
    for i, p in enumerate(pix):
        nb = hp.query_disc(nside, hp.pix2vec(nside, p), np.radians(margin_deg))
        keep[i] = inside[nb].all()
    return keep


def common_resolution(cube, freqs, pixel_indices, nside, target_fwhm_deg=None):
    """Reconvolve every channel to a single angular resolution.

    Standard intensity-mapping practice, and absent from experiment 006 as
    originally run. The beam narrows 12.5% across 350-400 MHz, and an
    unmatched resolution imprints exactly the kind of spectral structure a
    foreground filter cannot remove -- so leaving it out attributes to the
    survey a chromatic term that any real pipeline would have taken out.

    Each channel is convolved with a Gaussian of FWHM
    ``sqrt(target^2 - fwhm(nu)^2)`` to bring it to the widest beam in the
    band. That composition rule is *exact* for the Gaussian beam these
    experiments use (`ska_beam_func`); it would only be approximate for the
    tapered-aperture beam, which is not used here.

    Measured effect (interior pixels, 4 modes removed): the raster's residual
    halves, 48.9x -> 24.0x the HI, while the drift is unchanged at ~300x. The
    drift has too few measured modes for a cleaner beam to matter.
    """
    import healpy as hp
    from ska_common import ska_beam_fwhm_deg

    pix = np.asarray(pixel_indices)
    fwhm = np.array([ska_beam_fwhm_deg(f) for f in freqs])
    target = float(target_fwhm_deg if target_fwhm_deg is not None else fwhm.max())

    out = np.zeros_like(np.asarray(cube, float))
    npix = hp.nside2npix(nside)
    for i in range(out.shape[0]):
        kernel = np.sqrt(max(target**2 - fwhm[i]**2, 0.0))
        if kernel <= 0.01:
            out[i] = cube[i]
            continue
        full = np.zeros(npix)
        full[pix] = cube[i]
        out[i] = hp.smoothing(full, fwhm=np.radians(kernel), verbose=False)[pix]
    return out


def common_patch(operator_sets):
    """Pixels selected by *every* channel of *every* strategy.

    Operator pixel selection is beam-dependent, so each channel picks a
    slightly different set. Per-channel grids are themselves a chromatic
    effect and would confound the test, so the analysis runs on the
    intersection.
    """
    pix = None
    for ops in operator_sets:
        for mm in ops:
            p = np.asarray(mm.pixel_indices)
            pix = p if pix is None else np.intersect1d(pix, p)
    return pix


# ---------------------------------------------------------------------------
# Solving
# ---------------------------------------------------------------------------

# Prior variance in units of the truth patch variance. 1.0 is what every
# published number in experiment 006 used; ``ska_hi_priorvar.py`` sweeps it.
#
# The prior mean is already flat (a constant at the patch mean, set by each
# caller), so this scalar is the ONLY place the prior still carries information
# about the truth. It is also the data-vs-prior crossover in
# ``tr(WA) = sum lambda/(lambda + S^-1)``, which means it moves the measured
# mode count and the beam + prior floor together -- the two axes of the design
# figure. Raising it loosens the prior.
PRIOR_VAR_SCALE = 1.0


def _solve_args(mm, fg_truth_patch):
    """Noise weights and prior, held identical across every solve.

    ``noise_variance`` is the series' truth-dependent weighting (radiometer
    noise scaling with Tsys, which is dominated by the foreground). The HI is
    ~1e-5 of the foreground, so it does not perturb the weights.
    """
    nv_floor = X.WHITE_VAR * (1e-3 * float(np.mean(fg_truth_patch)))**2
    return dict(
        Tsky_prior_inv_cov_diag=np.ones_like(fg_truth_patch)
        / (PRIOR_VAR_SCALE
           * max(float(np.std(fg_truth_patch)), 1e-3)**2),
        noise_variance=[X.WHITE_VAR * (np.asarray(o) @ fg_truth_patch)**2
                        + nv_floor for o in mm.Tsys_operators],
        regularization=1e-12, return_full_cov=False)


def solve(mm, tod_group, fg_truth_patch, prior_mean):
    return mm(TOD_group=tod_group, dtime=X.DT, Tsky_prior_mean=prior_mean,
              **_solve_args(mm, fg_truth_patch))[0]


def response(mm, sky_patch, fg_truth_patch):
    """Map-domain response ``W A s`` to a sky component.

    The map-maker is a Wiener filter, affine in the data: ``solve(d) = W d +
    (A^T N^-1 A + S^-1)^-1 S^-1 mu``. Solving the component's own TOD with a
    zero prior mean kills the second term and returns ``W A s`` exactly, which
    is what ``solve(d + A s) - solve(d)`` would give at twice the cost.
    ``validate_linearity`` checks this identity numerically.
    """
    tod = [np.asarray(o) @ sky_patch for o in mm.Tsys_operators]
    return solve(mm, tod, fg_truth_patch, np.zeros_like(sky_patch))


# ---------------------------------------------------------------------------
# Radial power spectrum
# ---------------------------------------------------------------------------

def _taper(nchan):
    """Blackman-Harris along frequency.

    Not optional in practice: the foregrounds are ~1e4 times the HI here, so
    spectral leakage from the band edges would swamp the signal at high k_par.
    """
    return np.blackman(nchan)


def pk_par(cube_a, cube_b=None, dr_mpc=1.0, taper=True, subtract_mean=True):
    """Line-of-sight (cross-)power spectrum, averaged over pixels.

    ``cube`` is ``(nchan, npix)``. Returns ``(k, P)`` with k in Mpc^-1 and P in
    the square of the map units times Mpc.
    """
    a = np.asarray(cube_a, float)
    b = a if cube_b is None else np.asarray(cube_b, float)
    nchan = a.shape[0]
    if subtract_mean:
        a = a - a.mean(axis=0, keepdims=True)
        b = b - b.mean(axis=0, keepdims=True)
    if taper:
        w = _taper(nchan)[:, None]
        a, b = a * w, b * w
        norm = np.mean(_taper(nchan)**2)
    else:
        norm = 1.0

    L = nchan * dr_mpc
    fa = np.fft.rfft(a, axis=0)
    fb = np.fft.rfft(b, axis=0)
    p = np.real(fa * np.conj(fb)).mean(axis=1) * dr_mpc**2 / (L * norm)
    k = 2.0 * np.pi * np.fft.rfftfreq(nchan, d=dr_mpc)
    return k[1:], p[1:]        # drop the DC mode


def channel_dr_mpc(cfg: HIBandConfig):
    """Mean comoving depth of one channel, in Mpc."""
    r = comoving_distances(cfg)
    return float(abs(r[0] - r[-1]) / (len(r) - 1))


# ---------------------------------------------------------------------------
# PCA and the transfer function
# ---------------------------------------------------------------------------

def pca_clean(cube, nmodes):
    """Remove the ``nmodes`` highest-variance frequency-frequency eigenmodes.

    Thin wrapper over ``fastbox.filters.pca_filter``, which wants frequency on
    the last axis while our cubes are ``(nchan, npix)``.
    """
    from fastbox.filters import pca_filter
    return pca_filter(np.asarray(cube).T, nmodes).T


def transfer_function(data_cube, hi_mocks, nmodes, dr_mpc, verbose=True):
    """Estimate T(k_par) by mock injection.

    ``hi_mocks`` is a list of map-domain HI responses ``W A s_m`` for
    independent realisations -- already through the map-maker, so the transfer
    function isolates the cleaning step alone.
    """
    clean = pca_clean(data_cube, nmodes)
    num, den = [], []
    for i, mock in enumerate(hi_mocks):
        clean_m = pca_clean(data_cube + mock, nmodes)
        k, p_num = pk_par(clean_m - clean, mock, dr_mpc=dr_mpc)
        _, p_den = pk_par(mock, dr_mpc=dr_mpc)
        num.append(p_num)
        den.append(p_den)
        if verbose and (i + 1) % 5 == 0:
            print(f"    mock {i + 1}/{len(hi_mocks)}", flush=True)
    num, den = np.asarray(num), np.asarray(den)
    t_realisations = num / den
    return k, t_realisations.mean(axis=0), t_realisations


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def validate_linearity(mm, fg_truth_patch, tod_group, sky_patch, prior_mean):
    """Check ``solve(d + A s) - solve(d) == response(s)`` to round-off."""
    base = solve(mm, tod_group, fg_truth_patch, prior_mean)
    tod_inj = [np.asarray(t) + np.asarray(o) @ sky_patch
               for t, o in zip(tod_group, mm.Tsys_operators)]
    both = solve(mm, tod_inj, fg_truth_patch, prior_mean)
    direct = response(mm, sky_patch, fg_truth_patch)
    diff = np.max(np.abs((both - base) - direct))
    scale = np.max(np.abs(direct))
    return diff / scale


def validate_injection(strategy, freq_mhz, hi_patch, verbose=True):
    """Compare operator-injected HI against a full TOD simulation of the same sky.

    Quantifies the one real approximation in the design: the HI is injected as
    ``A s`` through the map-maker's own forward model rather than simulated, so
    it misses the sub-beam and outside-patch structure a real simulation would
    carry. Returns the fractional TOD rms difference.
    """
    import healpy as hp
    from limTOD import TODSim

    mm = X.build_operator(strategy, freq_mhz, verbose=False)
    pix = np.asarray(mm.pixel_indices)
    full = np.zeros(hp.nside2npix(X.NSIDE))
    full[pix] = hi_patch

    def hi_sky(*, freq, nside):
        return full if nside == X.NSIDE else hp.ud_grade(full, nside)

    sim = TODSim(ant_latitude_deg=X.SKA_LAT, ant_longitude_deg=X.SKA_LON,
                 ant_height_m=X.SKA_HGT, beam_func=X.ska_beam_func,
                 sky_func=hi_sky, beam_nside=X.NSIDE, sky_nside=X.NSIDE)
    out = []
    for i, p in enumerate(X.STRATEGIES[strategy]["pointings"]()):
        tod, _, _, _ = sim.generate_TOD(
            freq_list=[freq_mhz], time_list=p["tlist"],
            azimuth_deg_list=p["azlist"], elevation_deg=p["el"],
            start_time_utc=X.STRATEGIES[strategy]["utc"],
            white_noise_var=0.0, return_LSTs=True, normalize_beam=False,
            truncate_frac_thres=X.TRUNCATE_FRAC_THRES)
        sim_tod = np.asarray(tod[0], np.float64)
        fwd_tod = np.asarray(mm.Tsys_operators[i]) @ hi_patch
        out.append(float(np.std(sim_tod - fwd_tod) / np.std(sim_tod)))
        if verbose:
            print(f"  [{strategy} pass {i}] injected-vs-simulated HI TOD: "
                  f"{out[-1]:.2%}", flush=True)
    return out
