"""Shared helpers for the SKA-like drift-scan notebooks.

Site constants, a chromatic Gaussian beam for an SKA-Mid-style 15 m dish,
drift-scan timing helpers and the AltAz -> ICRS boresight-track transform.
Plot helpers (``plot_patch``, ``plot_map_compare``) are shared with the DSA
notebooks — import them from ``examples/DSA/dsa_vis.py``.
"""

from __future__ import annotations

import os

import numpy as np
import healpy as hp

from limTOD import example_symm_beam_map

# ---------------------------------------------------------------------------
# SKA-Mid site constants (Karoo, South Africa — same values as the TODSim
# defaults, kept explicit here so the notebooks read unambiguously)
# ---------------------------------------------------------------------------
SKA_LAT = -30.7130  # deg
SKA_LON = 21.4430   # deg
SKA_HGT = 1054.0    # m

SKA_DISH_DIAMETER_M = 15.0  # SKA-Mid dish

# Sidereal day in seconds. Offsetting consecutive nights by this value keeps
# the LST (hence RA) window identical from night to night.
SIDEREAL_DAY_SEC = 86164.0905

C_M_PER_S = 299792458.0


def ska_beam_fwhm_deg(freq_mhz, dish_diameter_m=SKA_DISH_DIAMETER_M):
    """Airy-disc FWHM ~ 1.22 lambda/D, in degrees, for ``freq_mhz`` in MHz."""
    lam = C_M_PER_S / (freq_mhz * 1e6)
    return np.degrees(1.22 * lam / dish_diameter_m)


def ska_beam_func(*, freq, nside):
    """Chromatic symmetric-Gaussian beam for an SKA-Mid-style dish.

    Signature matches ``limTOD.TODSim``'s ``beam_func`` protocol
    (``beam_func(*, freq, nside) -> 1-D HEALPix map``, sum-normalised).
    Unlike the DSA zenith beam this model is chromatic: the FWHM follows
    1.22 lambda/D with D = 15 m, so ``freq`` (MHz) matters.
    """
    fwhm = ska_beam_fwhm_deg(freq)
    return example_symm_beam_map(freq=freq, nside=nside, FWHM=fwhm)


# ---------------------------------------------------------------------------
# Sky model in the correct frame
# ---------------------------------------------------------------------------

def gdsm_equatorial_sky_model(*, freq, nside):
    """GDSM sky map rotated to equatorial (celestial) coordinates.

    ``limTOD.GDSM_sky_model`` returns pygdsm's native *galactic*-frame map
    unrotated, while the TOD simulator and map-maker point the beam in
    *equatorial* coordinates — so using it directly mislabels the sky (the
    galactic plane ends up along Dec ≈ 0 instead of its true celestial
    track). This wrapper rotates galactic → celestial at the native
    resolution before degrading, so RA/Dec in the SKA notebooks mean what
    they say.
    """
    from pygdsm import GlobalSkyModel

    gsm = GlobalSkyModel()
    skymap = gsm.generate(freq)  # native galactic frame
    rot = hp.Rotator(coord=["G", "C"])
    skymap = rot.rotate_map_pixel(skymap)
    return hp.ud_grade(skymap, nside_out=nside)


# ---------------------------------------------------------------------------
# Tapered-aperture ("realistic Airy") beam
# ---------------------------------------------------------------------------

def _tapered_aperture_field(x, taper=0.2):
    """Far-field voltage pattern of a circular aperture with illumination
    ``g(r) = taper + (1 - taper) * (1 - r^2)`` (r = normalised aperture
    radius), evaluated at ``x = pi * D * sin(theta) / lambda``.

    Closed form via Lambda functions: the uniform-disc term is 2 J1(x)/x
    and the (1 - r^2) term is 8 J2(x)/x^2; each tends to 1 as x -> 0 and
    they are weighted by their aperture integrals (1 and 1/2). ``taper``
    is the edge-illumination amplitude: taper = 1 recovers the uniform
    Airy pattern (first sidelobe -17.6 dB); the default 0.2 (-14 dB edge
    taper, typical of real dishes) widens the main lobe to ~1.15 lambda/D
    and lowers the first sidelobe to ~ -23 dB.
    """
    from scipy.special import j1, jn

    x = np.asarray(x, dtype=float)
    small = np.abs(x) < 1e-8
    xs = np.where(small, 1.0, x)
    E_uniform = np.where(small, 1.0, 2.0 * j1(xs) / xs)
    E_tapered = np.where(small, 1.0, 8.0 * jn(2, xs) / xs**2)
    w_uniform = taper
    w_tapered = (1.0 - taper) / 2.0
    return (w_uniform * E_uniform + w_tapered * E_tapered) / (
        w_uniform + w_tapered)


def tapered_aperture_profile(theta_deg, freq_mhz,
                             dish_diameter_m=SKA_DISH_DIAMETER_M,
                             taper=0.2):
    """Power pattern P(theta)/P(0) of the tapered aperture (for plotting)."""
    lam = C_M_PER_S / (freq_mhz * 1e6)
    theta = np.radians(np.asarray(theta_deg, dtype=float))
    x = np.pi * dish_diameter_m / lam * np.sin(np.clip(theta, 0, np.pi / 2))
    power = _tapered_aperture_field(x, taper=taper) ** 2
    return np.where(theta > np.pi / 2, 0.0, power)


def tapered_aperture_fwhm_deg(freq_mhz, dish_diameter_m=SKA_DISH_DIAMETER_M,
                              taper=0.2):
    """Numerical FWHM (degrees) of the tapered-aperture power pattern."""
    lam = C_M_PER_S / (freq_mhz * 1e6)
    first_null_deg = np.degrees(1.6 * lam / dish_diameter_m)
    theta = np.linspace(0.0, first_null_deg, 20001)
    power = tapered_aperture_profile(theta, freq_mhz,
                                     dish_diameter_m=dish_diameter_m,
                                     taper=taper)
    return 2.0 * float(np.interp(0.5, power[::-1], theta[::-1]))


def ska_airy_beam_func(*, freq, nside):
    """Chromatic tapered-aperture beam with realistic sidelobe rings.

    Same ``beam_func`` protocol as ``ska_beam_func`` (sum-normalised 1-D
    HEALPix map), but with the full circular-aperture diffraction pattern
    instead of a Gaussian: sidelobe rings spaced ~lambda/D (≈ 3.3° at
    350 MHz), first ring at ~ -23 dB for the default 0.2 edge taper.
    The back hemisphere (theta > 90°) is set to zero.
    """
    lam = C_M_PER_S / (freq * 1e6)
    theta, _ = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
    x = np.pi * SKA_DISH_DIAMETER_M / lam * np.sin(
        np.clip(theta, 0, np.pi / 2))
    beam_map = _tapered_aperture_field(x) ** 2
    beam_map[theta > np.pi / 2] = 0.0
    beam_map /= beam_map.sum()
    return beam_map


def drift_scan_night(duration_s, dt=2.0, night=0, azimuth_deg=0.0):
    """Time/azimuth lists for one night of a drift scan.

    The telescope is parked (constant azimuth; pass the constant elevation
    separately to ``TODSim.generate_TOD``). Night ``n`` is offset by ``n``
    sidereal days so every night covers the same LST window and the map-maker
    sees the same RA strip at each elevation setting.

    Returns
    -------
    time_list : (ntime,) array
        Offsets in seconds from the survey's single ``start_time_utc``.
    azimuth_deg_list : (ntime,) array
        Constant array at ``azimuth_deg``.
    """
    ntime = int(round(duration_s / dt))
    t_list = night * SIDEREAL_DAY_SEC + np.arange(ntime) * dt
    az_list = np.full(ntime, azimuth_deg)
    return t_list, az_list


def save_results_pdf(number, slug, description, image_paths,
                     rms_table=None, out_dir="results"):
    """Archive one experiment's results as a numbered PDF.

    Writes ``<out_dir>/<NNN>_<slug>.pdf`` with a title page (short setup
    description + optional results table) followed by one page per image.

    Parameters
    ----------
    number : int
        Experiment number; becomes the zero-padded filename prefix so the
        archive sorts chronologically.
    slug : str
        Short kebab-case name, e.g. ``"short-drift-tapered-beam"``.
    description : str
        One-or-two-sentence setup description shown on the title page.
    image_paths : list of str
        Figure PNGs (as saved by the notebook) to embed, one per page.
    rms_table : list of (label, value) pairs, optional
        Rendered monospaced under the description (e.g. residual RMS rows).

    Returns
    -------
    path : str
        The written PDF path.
    """
    import textwrap
    from datetime import date

    import matplotlib.image as mpimg
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f"{number:03d}_{slug}.pdf")

    with PdfPages(path) as pdf:
        # Title page (landscape A4, matching the wide map figures)
        fig = plt.figure(figsize=(11.69, 8.27))
        fig.text(0.06, 0.90, f"Experiment {number:03d} — {slug}",
                 fontsize=18, weight="bold")
        fig.text(0.06, 0.86, date.today().isoformat(), fontsize=10,
                 color="0.4")
        body = "\n".join(textwrap.fill(par, width=95)
                         for par in description.split("\n"))
        fig.text(0.06, 0.80, body, fontsize=12, va="top")
        if rms_table:
            width = max(len(str(k)) for k, _ in rms_table)
            lines = "\n".join(f"{str(k):{width}s}  {v}"
                              for k, v in rms_table)
            fig.text(0.06, 0.55, lines, fontsize=11, va="top",
                     family="monospace")
        pdf.savefig(fig)
        plt.close(fig)

        for img_path in image_paths:
            img = mpimg.imread(img_path)
            fig = plt.figure(figsize=(11.69, 8.27))
            ax = fig.add_axes([0.02, 0.05, 0.96, 0.90])
            ax.imshow(img)
            ax.axis("off")
            ax.set_title(os.path.basename(img_path), fontsize=9,
                         color="0.4")
            pdf.savefig(fig)
            plt.close(fig)

    print(f"Archived results to {path}")
    return path


def azel_to_radec(az_deg, el_deg, time_list_sec, start_time_utc, location):
    """Boresight ICRS track of a parked telescope (inverse of
    ``dsa_vis.radec_to_azel``).

    Parameters
    ----------
    az_deg, el_deg : float
        Fixed horizontal pointing in degrees.
    time_list_sec : array-like
        Time offsets in seconds from ``start_time_utc``.
    start_time_utc : str
        UTC start time parseable by ``astropy.time.Time``.
    location : astropy.coordinates.EarthLocation
        Observatory location.

    Returns
    -------
    ra_deg, dec_deg : np.ndarray
        ICRS coordinates of the boresight at each time sample.
    """
    from astropy.coordinates import AltAz, ICRS, SkyCoord
    from astropy.time import Time, TimeDelta
    from astropy import units as u

    start = Time(start_time_utc)
    times = start + TimeDelta(np.asarray(time_list_sec), format="sec")
    pointing = SkyCoord(
        az=az_deg * u.deg, alt=el_deg * u.deg,
        frame=AltAz(obstime=times, location=location),
    )
    icrs = pointing.transform_to(ICRS())
    return icrs.ra.deg, icrs.dec.deg
