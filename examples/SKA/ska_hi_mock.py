"""Mock HI brightness-temperature signal for the drift-scan feasibility series.

Experiments 001-005 contain no 21 cm signal at all: every residual quoted in
``HANDOFF.md`` measures reconstruction fidelity of the diffuse *foreground*
sky. This module supplies the missing ingredient for the queued signal-recovery
test (``ANALYSIS_LOG.md``, "HI plus foreground cleaning"): a cosmological HI
cube, projected onto the HEALPix patch the map-maker actually solves.

limTOD has no 21 cm model -- ``sky_model.py`` offers ``GDSM_sky_model`` and
``generate_gaussian_field``, and the latter is an Alonso et al. (2014)
*foreground* covariance (zero mean, spectrally smooth by construction), the
opposite of what has to be injected. The signal therefore comes from fastbox
(``/home/geoff/FastBox``), whose ``generate_hi_mock`` runs the standard chain
Gaussian density -> HI bias -> log-normal -> linear RSD -> Tb(z).

Two coordinate facts drive the design:

* fastbox works in a Cartesian comoving box; limTOD works in HEALPix. The
  bridge here is a genuine 3D lightcone interpolation -- each (pixel, channel)
  pair is placed at its true comoving position ``r(z) * n_pix`` and the box is
  sampled there -- rather than a flat-sky tangent-plane projection. Over the
  raster's ~20 deg patch radius a gnomonic projection would stretch scales by
  4.3% at the edge; the 3D route has no such error.

* Two approximations remain, both documented rather than hidden. The field is
  realised at a single effective redshift (no growth across the band, which
  spans dz ~ 0.5 at the band bottom), and fastbox applies RSD along the box
  z-axis while the true radial direction tilts away from it by up to the patch
  radius, so the RSD component is cos(20 deg) = 0.94 of correct at the extreme
  edge. Neither matters for a transfer-function measurement, which is a ratio.

Runs under ``gibbs_venv_312`` (pyccl 3.3.0); the limTOD venv has no pyccl.
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass, field

import numpy as np

FASTBOX_PATH = "/home/geoff/FastBox"
if FASTBOX_PATH not in sys.path:
    sys.path.insert(0, FASTBOX_PATH)

NU21_MHZ = 1420.405752

# Planck-like, matching the defaults used elsewhere in the fastbox examples.
#
# transfer_function: pyccl defaults to 'boltzmann_camb', which needs camb.
# camb 2.0.3 and pyccl 3.3.0 disagree over the NonLinear flag (pyccl asserts
# camb left it at NonLinear_none and camb 2.x no longer does), so the venv is
# left without camb and the analytic Eisenstein & Hu (1998) transfer function
# is used instead. Halofit still supplies the non-linear correction on top.
# Adequate here by a wide margin: this mock exists to be injected and
# recovered, and the transfer function it feeds is a ratio, so percent-level
# differences in the input P(k) cancel.
COSMO_PARAMS = dict(Omega_c=0.25, Omega_b=0.05, h=0.7, n_s=0.96, sigma8=0.8,
                    transfer_function="eisenstein_hu")


@dataclass
class HIBandConfig:
    """A contiguous sub-band to run the signal-recovery test in.

    The sweep of experiment 005 used five channels spanning the whole of
    Band 1, which is useless for PCA foreground cleaning -- that needs many
    *contiguous* channels within one sub-band, so this is new simulation work
    rather than a re-solve of the cached operators.
    """

    f_lo_mhz: float
    f_hi_mhz: float
    nchan: int
    nside: int
    seed: int = 0
    # Grid spacing of the comoving box, in Mpc. Must comfortably resolve both
    # the channel depth and the HEALPix pixel; the beam smooths far coarser.
    cell_mpc: float = 20.0
    # Padding added to every box dimension, in Mpc, so interpolation never
    # reaches the grid edge.
    pad_mpc: float = 120.0
    sigma_nl: float = 120.0
    cosmo_params: dict = field(default_factory=lambda: dict(COSMO_PARAMS))

    @property
    def channel_freqs_mhz(self) -> np.ndarray:
        """Channel centres, low to high frequency."""
        edges = np.linspace(self.f_lo_mhz, self.f_hi_mhz, self.nchan + 1)
        return 0.5 * (edges[:-1] + edges[1:])

    @property
    def redshifts(self) -> np.ndarray:
        return NU21_MHZ / self.channel_freqs_mhz - 1.0

    @property
    def z_centre(self) -> float:
        z = self.redshifts
        return float(0.5 * (z.min() + z.max()))

    @property
    def slug(self) -> str:
        return (f"f{self.f_lo_mhz:04.0f}_{self.f_hi_mhz:04.0f}"
                f"_nc{self.nchan}_ns{self.nside}_s{self.seed}")


def _cosmology(cfg: HIBandConfig):
    import pyccl as ccl
    return ccl.Cosmology(**cfg.cosmo_params)


def comoving_distances(cfg: HIBandConfig) -> np.ndarray:
    """Comoving radial distance to each channel, in Mpc."""
    import pyccl as ccl
    cosmo = _cosmology(cfg)
    a = 1.0 / (1.0 + cfg.redshifts)
    return ccl.comoving_radial_distance(cosmo, a)


def patch_frame(pixel_indices: np.ndarray, nside: int):
    """Orthonormal box frame for a HEALPix patch.

    Returns ``(e_x, e_y, e_z)`` with ``e_z`` pointing at the patch centroid,
    so the box z-axis is the line of sight (the axis fastbox applies RSD
    along) and x/y are transverse.
    """
    import healpy as hp

    vec = np.asarray(hp.pix2vec(nside, np.asarray(pixel_indices)))
    e_z = vec.mean(axis=1)
    e_z /= np.linalg.norm(e_z)

    # Any vector not parallel to e_z; the pole is safe because these patches
    # sit at Dec ~ +9 deg.
    seed_vec = np.array([0.0, 0.0, 1.0])
    e_x = seed_vec - np.dot(seed_vec, e_z) * e_z
    e_x /= np.linalg.norm(e_x)
    e_y = np.cross(e_z, e_x)
    return e_x, e_y, e_z


def lightcone_coordinates(cfg: HIBandConfig, pixel_indices: np.ndarray):
    """Comoving box coordinates of every (channel, pixel) sample.

    Returns ``(u, v, w, frame)`` each of shape ``(nchan, npix)``, in Mpc,
    relative to the patch centroid direction at the band's central distance.
    """
    import healpy as hp

    pixel_indices = np.asarray(pixel_indices)
    e_x, e_y, e_z = patch_frame(pixel_indices, cfg.nside)
    n = np.asarray(hp.pix2vec(cfg.nside, pixel_indices))          # (3, npix)
    r = comoving_distances(cfg)                                    # (nchan,)
    r_c = float(0.5 * (r.min() + r.max()))

    # (nchan, npix) comoving positions projected onto the box axes.
    u = r[:, None] * (e_x @ n)[None, :]
    v = r[:, None] * (e_y @ n)[None, :]
    w = r[:, None] * (e_z @ n)[None, :] - r_c
    return u, v, w, (e_x, e_y, e_z, r_c)


def box_geometry(cfg: HIBandConfig, pixel_indices: np.ndarray):
    """Box side lengths and grid shape that cover the patch with padding."""
    u, v, w, _ = lightcone_coordinates(cfg, pixel_indices)
    spans = []
    for arr in (u, v, w):
        half = max(abs(float(arr.min())), abs(float(arr.max())))
        spans.append(2.0 * half + 2.0 * cfg.pad_mpc)
    nsamp = tuple(int(np.ceil(s / cfg.cell_mpc)) for s in spans)
    return tuple(spans), nsamp


def build_hi_cube(cfg: HIBandConfig, pixel_indices: np.ndarray, verbose=True):
    """Realise one mock HI brightness-temperature box, in mK.

    Thin wrapper over ``fastbox.tracers.generate_hi_mock`` that sizes the box
    to the patch and seeds the global RNG that ``realise_density`` draws from,
    so a given ``cfg.seed`` is reproducible.
    """
    from fastbox.box import CosmoBox
    from fastbox.tracers import HITracer
    from numpy import fft

    box_scale, nsamp = box_geometry(cfg, pixel_indices)
    if verbose:
        print(f"[hi_mock] box {box_scale[0]:.0f} x {box_scale[1]:.0f} x "
              f"{box_scale[2]:.0f} Mpc on {nsamp} grid "
              f"(z_eff = {cfg.z_centre:.3f})", flush=True)

    np.random.seed(cfg.seed)
    box = CosmoBox(cosmo=_cosmology(cfg), box_scale=box_scale, nsamp=nsamp,
                   redshift=cfg.z_centre, realise_now=False)

    # Same chain as fastbox's generate_hi_mock, inlined so the box (needed for
    # the coordinate grids) stays available and the seed is explicit.
    box.realise_density()
    tracer = HITracer(box)
    delta_hi = box.delta_x * tracer.bias_HI()
    delta_ln = box.lognormal(delta_hi)
    vel_k = box.realise_velocity(delta_x=box.delta_x, inplace=True)
    vel_z = fft.ifftn(vel_k[2]).real
    delta_s = box.redshift_space_density(delta_x=delta_ln.real,
                                         velocity_z=vel_z,
                                         sigma_nl=cfg.sigma_nl,
                                         method='linear')
    cube_mk = tracer.signal_amplitude() * (1.0 + delta_s)

    # fastbox divides by k^2 in realise_velocity, so the DC mode warns and can
    # in principle poison the RSD field. It does not in practice, but a silent
    # NaN here would propagate into every power spectrum downstream.
    if not np.all(np.isfinite(cube_mk)):
        raise ValueError("HI cube contains non-finite values "
                         f"({np.count_nonzero(~np.isfinite(cube_mk))} cells)")

    if verbose:
        print(f"[hi_mock] Tb = {tracer.signal_amplitude():.4f} mK, "
              f"b_HI = {tracer.bias_HI():.3f}, "
              f"cube rms = {cube_mk.std():.4f} mK", flush=True)
    return cube_mk, box


def project_to_healpix(cube_mk, box, cfg: HIBandConfig, pixel_indices,
                       order=3):
    """Sample the comoving cube along the lightcone onto HEALPix channels.

    Returns an ``(nchan, npix)`` array in **K**, matching limTOD's sky units
    (``gdsm_equatorial_sky_model`` returns K).
    """
    from scipy.ndimage import map_coordinates

    u, v, w, _ = lightcone_coordinates(cfg, pixel_indices)

    # Box grids are linspace over [-L/2, L/2] with N points, so index = the
    # fractional position on that grid.
    def to_index(coord, axis_vals):
        lo, hi = axis_vals[0], axis_vals[-1]
        n = len(axis_vals)
        return (coord - lo) / (hi - lo) * (n - 1)

    idx = np.stack([to_index(u, box.x), to_index(v, box.y),
                    to_index(w, box.z)])

    for name, arr, n in (("x", idx[0], box.Nx), ("y", idx[1], box.Ny),
                         ("z", idx[2], box.Nz)):
        if arr.min() < 0 or arr.max() > n - 1:
            raise ValueError(
                f"lightcone leaves the box on the {name} axis "
                f"({arr.min():.1f} .. {arr.max():.1f} vs 0 .. {n - 1}); "
                "increase HIBandConfig.pad_mpc")

    maps_mk = map_coordinates(cube_mk, idx, order=order, mode="nearest")
    return maps_mk * 1e-3  # mK -> K


def hi_maps(cfg: HIBandConfig, pixel_indices, verbose=True):
    """Convenience: build a cube and project it. Returns ``(nchan, npix)`` K."""
    cube, box = build_hi_cube(cfg, pixel_indices, verbose=verbose)
    maps_k = project_to_healpix(cube, box, cfg, pixel_indices)
    if verbose:
        print(f"[hi_mock] projected to {maps_k.shape} HEALPix samples, "
              f"rms = {maps_k.std() * 1e3:.4f} mK", flush=True)
    return maps_k
