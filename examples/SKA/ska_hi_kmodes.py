"""Where the measured sky modes live in angular scale -- and what they buy.

``tr(WA)`` says *how many* sky modes a strategy measures. It does not say
*which*, and the series has already shown that the count alone is not a
sufficient statistic: at 22 modes/100 deg^2 the cross-linked raster's floor is
9.3x below a twelve-strip parked drift's. This decomposes the count by angular
scale, which is the axis a science requirement is actually written on.

The decomposition is exact rather than a projection onto a non-orthogonal
basis. Eigen-decompose ``A^T N^-1 A`` into modes ``u_i`` with eigenvalues
``lambda_i``; each contributes ``w_i = lambda_i / (lambda_i + S^-1)`` to
``tr(WA)``, between 0 (prior-filled) and 1 (data-dominated). Every ``u_i`` is
then given an angular power spectrum ``f_i(k)`` normalised to sum to one over
the k bins, by direct (non-uniform) Fourier transform on the tangent plane --
no regridding, so no interpolation artefacts. Then

    n_eff(k) = sum_i w_i f_i(k),      sum_k n_eff(k) = tr(WA) exactly,

so the curve is a decomposition of the published mode count, not a new
statistic that has to be reconciled with it.

The second panel turns that into a sufficiency statement, and the two cases it
shows are not symmetric.

Against an **external tracer** the floor is unmeasured foreground, which
the tracer knows nothing about: it is uncorrelated, so it costs variance only.
Treating the whole patch as a single band power of ``N_m`` modes,

    (dP_x / P_x)^2 = [1 + (1 + R) (1 + 1/(nbar P_g)) / r^2] / N_m,
    R = P_floor / P_HI,

so a detection at significance ``snr`` needs ``N_m = snr^2 [1 + (1 + R)(1 +
1/nPg)/r^2]`` modes, and -- modes scaling with area at fixed strategy -- a sky
area ``area * N_m / n_eff``.

In the **Hi auto-spectrum** there is no such area, and the figure quotes none.
The floor is ``(I - WA)s``: deterministic given the strategy and the sky, so it
enters the auto-spectrum as a *bias*, not as variance, and no amount of sky
averages it down. What the auto rows report is the bias itself, ``R``. The
variance-only count ``2(1 + R)^2`` is retained in the returned dict under
``auto_lower_bound_if_floor_known`` precisely because that is all it is: a lower
bound, valid only if the floor's mean power were known and subtracted -- which
would require knowing the sky the survey failed to measure.

    /home/geoff/limTOD/.venv/bin/python ska_hi_kmodes.py

Reads the cached operators in ``hi_cache/`` only; no TOD, no new runs.
"""
from __future__ import annotations

import os
import pickle
import sys
from dataclasses import dataclass

import healpy as hp
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import ska_hi_experiment as X
import ska_hi_prior as prior
import ska_hi_plots as P
from ska_hi_plots import INK, INK2, MUTED, SURFACE, BASELINE, GRID
from ska_common import gdsm_equatorial_sky_model, ska_beam_fwhm_deg, SKA_LAT
from ska_hi_ladder import operator_list

# --------------------------------------------------------------------------
# Bands. Everything that differs between 350-400 and 675-725 MHz lives here, so
# the figure is drawn by the same code for both and cannot drift between them.
# --------------------------------------------------------------------------

EL_MIN_DEG = 30.0
# The ceiling is the sky SKA-Mid can actually reach, not 4*pi. A source at
# declination d reaches elevation 90 - |lat - d|, so from the Karoo it clears
# EL_MIN only for d <= lat + (90 - EL_MIN). Using the full sphere (41253 deg^2)
# overstates the available sky by a third and is the wrong line to draw.
_DEC_MAX = min(SKA_LAT + (90.0 - EL_MIN_DEG), 90.0)
SKY_LIMIT_DEG2 = float(2 * np.pi * (np.sin(np.radians(_DEC_MAX)) + 1.0)
                       * (180.0 / np.pi) ** 2)
FULL_SKY_DEG2 = 4 * np.pi * (180.0 / np.pi) ** 2        # 41253, drawn faint

R_CROSS = 1.0           # HI-tracer correlation coefficient
NPG = (np.inf, 1.0)     # tracer shot noise: perfect, and nbar P_g = 1
SNR = (1, 3, 5)         # detection significance; dP/P = 1/snr
Q_LABEL = 3             # the second area printed, on the raster row only
DIM = 0.35              # alpha for anything past the reachable sky. MUTED is
                        # within a hair of INK in the current palette, so
                        # transparency, not colour, has to carry the distinction.
NMODES = 4              # cleaning depth the floor/HI ratio is read at

# Shared x range for the merged two-band figure. Both groups must share it or
# the eye compares areas across groups that are not on the same scale.
X_BANDS = (1e2, 3e8)
# Per-k_par (--kpar) the areas drop by one to two orders of magnitude, so the
# single-bin window would clip the 675-725 MHz rows off the left edge. Both
# sky lines still have to fit, hence the top end.
X_BANDS_KPAR = (3.0, 4e6)
# the single-band panels keep their auto rows, whose floor/Hi text sits at the
# left edge, so they need their own window rather than X_BANDS_KPAR
X_SINGLE_KPAR = (3.0, 4e6)

# "1-dish" on the raster is not decoration: this is an autocorrelation paper
# about a vaguely MeerKLASS-like scan, so "raster x ..." invites the reading
# that dishes are being correlated with each other. Nothing here is. The drift
# strip counts are the declination ladder of Sec. "Dish count and the
# declination ladder", where N strips on N sidereal days and N dishes at N
# elevations give the identical operator A.
ROW_LABEL = {"drift": "3-strip drift", "drift12": "12-strip drift",
             "drift12m": "12-strip drift", "raster": "1-dish raster"}


@dataclass(frozen=True)
class Band:
    """One observing band, and where every number for it comes from.

    ``ablation`` must be a SINGLE run containing all three strategies:
    ``build_hi_cube`` sizes its box from the union of a run's strategies, so
    floors are comparable within a run and not across runs (15-18% realisation
    scatter; see HANDOFF). That is why the matched runs are used rather than
    the published ones, at both bands.
    """
    tag: str
    freq_mhz: float                 # channel 0, the frequency the ops are cached at
    nside: int
    z_centre: float
    d_comoving_mpc: float
    ablation: str                   # results/<...>.npz, one run, all strategies
    strategies: tuple               # row order: 3-strip, 12-strip, raster
    ladder: str                     # fallback geometry for the 12-strip row
    ladder_row: tuple               # (n_strips, spacing_beams) in that sweep
    out: str                        # results/<...>.npz this writes
    suffix: str                     # figure-name suffix ("" keeps hi_kmodes.png)
    x_lim: tuple
    auto_text_x: float

    @property
    def fwhm_rad(self):
        return np.radians(ska_beam_fwhm_deg(self.freq_mhz))


BANDS = {
    # d_comoving is NOT recomputed here: it is the series' own anchor, a
    # 3.98 deg beam at 350 MHz subtending 401 Mpc at z = 2.80, i.e.
    # k_perp = 2 pi / 401 = 0.0157 (Sec. "Injected signal and geometry").
    # Deriving it this way guarantees the figure agrees with the paper. It is
    # 6% off the exact value, which is the rounding in the quoted 0.0157.
    "f350": Band(
        tag="f350", freq_mhz=350.78125, nside=64, z_centre=2.80,
        d_comoving_mpc=(2 * np.pi
                        / np.radians(ska_beam_fwhm_deg(350.78125))) / 0.0157,
        ablation="hi_ablation_amp350_f350_400_nc32_ns64.npz",
        strategies=("drift", "drift12", "raster"),
        ladder="hi_ladder_long_ns64.npz", ladder_row=(12, 0.5),
        out="hi_kmodes_f350_ns64.npz", suffix="",
        x_lim=(80.0, 3e8), auto_text_x=100.0),
    # No published anchor exists at this band, so d_comoving is the exact
    # comoving distance to the band centre (z = 1.028) under the pipeline's own
    # cosmology, ska_hi_mock.COSMO_PARAMS. Stated rather than imported so this
    # module keeps running under limTOD/.venv, which has no pyccl.
    "f700": Band(
        tag="f700", freq_mhz=675.78125, nside=128, z_centre=1.028,
        d_comoving_mpc=3371.1,
        ablation="hi_ablation_amp700_f675_725_nc32_ns128.npz",
        strategies=("drift", "drift12m", "raster"),
        ladder="hi_ladder_f700d_ns128.npz", ladder_row=(12, 0.579),
        out="hi_kmodes_f700_ns128.npz", suffix="_f700",
        # everything lands inside the reachable sky at this band, so the axis
        # only has to carry the 4pi line and the legend -- 3e7 was dead space
        x_lim=(30.0, 2e6), auto_text_x=38.0),
}
DEFAULT_BAND = "f350"

# Regression guard on the f350 raster's 1-sigma cross area. NOT a constant of
# nature: it moves whenever the floor or the mode count moves, and the guard
# exists so a refactor cannot shift it silently. Its history:
#
#   2937  own-patch prior -- the value the draft quotes
#   2946  + the cross-power's own sample-variance term (N_m: R -> 2 + R), <0.5%
#   2419  + the amplitude-only prior (2026-10-05). The raster has the widest
#         patch, so correcting the prior helps it most: floor 25.6 -> 22.3.
#
# At 2419 the raster's 3 sigma lands at 2.18e4 deg^2, comfortably inside the
# reachable 30717 rather than just outside it, so "3 sigma at the edge of the
# reachable sky" becomes "3 sigma reachable". The draft needs that change.
DRAFT_RASTER_1SIG_DEG2 = 2419.0


def eig_decomposition(mm, band):
    """``lambda_i``, ``u_i`` and ``S^-1`` -- the recipe ``mode_count`` uses."""
    pix = np.asarray(mm.pixel_indices)
    truth = gdsm_equatorial_sky_model(freq=band.freq_mhz, nside=band.nside)[pix]
    nv_floor = X.WHITE_VAR * (1e-3 * float(np.mean(truth))) ** 2
    M = None
    for o in operator_list(mm):
        A = np.asarray(o, float)
        var = X.WHITE_VAR * (A @ truth) ** 2 + nv_floor
        term = A.T @ (A / var[:, None])
        M = term if M is None else M + term
    lam, U = np.linalg.eigh(M)
    lam = np.maximum(lam[::-1], 0.0)
    U = U[:, ::-1]
    s_inv = prior.s_inv(band.freq_mhz, truth, label=band.tag)
    return lam, U, s_inv, pix


def tangent_xy(pix, nside):
    """Tangent-plane offsets in radians about the patch centre."""
    th, ph = hp.pix2ang(nside, pix)
    ra, dec = np.unwrap(ph), np.pi / 2 - th
    ra0, dec0 = ra.mean(), dec.mean()
    return (ra - ra0) * np.cos(dec0), dec - dec0


def n_eff_of_k(mm, kedges, band, n_ang=24):
    """``n_eff`` per |k| bin. Sums to ``tr(WA)`` by construction."""
    lam, U, s_inv, pix = eig_decomposition(mm, band)
    w = lam / (lam + s_inv)                       # (nmode,), sums to tr(WA)
    x, y = tangent_xy(pix, band.nside)
    kmid = np.sqrt(kedges[:-1] * kedges[1:])
    angles = np.linspace(0.0, np.pi, n_ang, endpoint=False)

    Pk = np.zeros((len(kmid), U.shape[1]))
    for b, kk in enumerate(kmid):
        acc = np.zeros(U.shape[1])
        for a in angles:
            phase = np.exp(-1j * (kk * np.cos(a) * x + kk * np.sin(a) * y))
            acc += np.abs(phase @ U) ** 2
        Pk[b] = acc / n_ang
    # normalise each mode's spectrum to unit total, so the k sum is exact
    Pk /= np.maximum(Pk.sum(axis=0, keepdims=True), 1e-300)
    return kmid, Pk @ w, float(w.sum())


def _need(modes, n_eff, area):
    """Modes, the factor on what is measured now, and the sky area that buys."""
    return dict(modes=float(modes), factor=float(modes / n_eff),
                area=float(area * modes / n_eff))


def sufficiency(floor_amp, n_eff, area, r=R_CROSS, npgs=NPG, snrs=SNR,
                kpar=None):
    """What the measured modes are enough for -- cross-correlation, and auto.

    ``cross[nPg][snr]`` is the sky area needed to detect the Hi x tracer
    cross-power amplitude at ``snr``; see the module docstring for the
    variance. ``auto`` carries no area on purpose: the floor is a bias there,
    so sky area cannot buy a detection, and the only number that means
    anything is the bias itself.

    ``kpar`` is the per-k_par term dict from ``ska_nm_kpar.kpar_terms``. Given
    one, ``N_m`` stops being the angular mode count of a single channel and
    becomes an inverse-variance sum over the cube's line-of-sight degrees of
    freedom, each bin with its OWN floor/Hi ratio and its own transfer
    function. Left at ``None`` every number is bit-for-bit what the single-bin
    forecast gave, which is what the regression guards below check.
    """
    import ska_nm_kpar as NK

    R = float(floor_amp) ** 2                      # P_floor / P_HI
    cross = {}
    for npg in npgs:
        shot = 1.0 if np.isinf(npg) else 1.0 + 1.0 / float(npg)
        if kpar is None:
            at_1sig = 1.0 + (1.0 + R) * shot / float(r) ** 2
        else:
            # Omega = snr^2 / [rho sum_b n_b/(1 + (1+R_b) shot/r^2)], and
            # _need turns "modes" into area as area * modes / n_eff = modes/rho
            at_1sig = 1.0 / NK.weights(kpar, shot=shot, r=float(r))
        cross[float(npg)] = {int(q): _need(q ** 2 * at_1sig, n_eff, area)
                             for q in snrs}
    auto_modes = (2.0 * (1.0 + R) ** 2 if kpar is None
                  else 1.0 / NK.auto_weights(kpar))
    return dict(
        R=R, cross=cross,
        auto=dict(
            bias_ratio=R,
            # NOT a requirement: the variance-only count, which holds only if
            # the floor's mean power is known and subtracted.
            auto_lower_bound_if_floor_known=_need(auto_modes, n_eff, area)))


def floor_amp(band, nmodes=NMODES):
    """floor/HI in AMPLITUDE per strategy, from the band's own matched run.

    The median over channels of sqrt(floor_n / P_HI) at common resolution --
    keys ``{s}_cr_floor_{n}`` and ``{s}_cr_p_hi``, the recipe
    ``ska_hi_design.collect._floor`` uses, so these are the same numbers
    ``fig_design`` and ``fig_ladder`` draw.

    Rounded to one decimal on purpose. The ratio enters squared, so the extra
    digits are noise against the 15-18% realisation scatter, and rounding keeps
    the printed values identical to the band table in HANDOFF and to what the
    paper quotes (399.6 / 239.3 / 25.6 at 350-400, 16.1 / 36.0 at 675-725).
    """
    path = os.path.join(_HERE, "results", band.ablation)
    if not os.path.exists(path):
        sys.exit(f"need results/{band.ablation} for band {band.tag} -- "
                 "pull it before drawing the sufficiency figure")
    d = np.load(path)
    out = {}
    for s in band.strategies:
        if f"{s}_cr_p_hi" not in d:
            sys.exit(f"{band.ablation} has no '{s}' arm "
                     f"(has: {sorted({k.split('_cr_')[0] for k in d})})")
        out[s] = round(float(np.median(np.sqrt(d[f"{s}_cr_floor_{nmodes}"]
                                               / d[f"{s}_cr_p_hi"]))), 1)
    return out


def ladder_geometry(band):
    """``n_eff`` and patch area for the band's twelve-strip row, off the sweep.

    Only used when no operator is cached for that strategy -- true at 350 MHz,
    where the ladder sweep builds its operators and discards them. At 675-725
    an ``op_drift12m_*`` pickle does exist, so the operator wins and this is
    just a cross-check.
    """
    path = os.path.join(_HERE, "results", band.ladder)
    n_want, sp_want = band.ladder_row
    if os.path.exists(path):
        d = np.load(path)
        m = ((d["n_strips"] == n_want)
             & (np.abs(d["spacing_beams"] - sp_want) < 1e-3))
        if m.any():
            i = int(np.argmax(m))
            return float(d["n_eff"][i]), float(d["area"][i])
    sys.exit(f"no {n_want}-strip / {sp_want}-beam row in results/{band.ladder}; "
             "the twelve-strip geometry has to come from somewhere")


def operator_path(band, strategy):
    return os.path.join(_HERE, "hi_cache",
                        f"op_{strategy}_f{band.freq_mhz:07.3f}_ns{band.nside}.pkl")


def _sig3(v):
    """A readable area or ratio: grouped integer below 1e5, else three
    significant figures as a mathtext power of ten. The cut is at 1e5 rather
    than 1e4 so the 3-sigma raster area (26,515) is grouped like the reachable
    sky it is being compared with (30,717), instead of reading 2.65e4."""
    if v < 100:
        return f"{v:,.1f}"          # per-k_par areas reach single digits
    if v < 1e5:
        return f"{v:,.0f}"
    mant, exp = f"{v:.2e}".split("e")
    return f"${mant}\\times 10^{{{int(exp)}}}$"


def _reach_segment(ax, x0, x1, y, **kw):
    """One horizontal run, solid inside the reachable sky and faded past it."""
    lo, hi = min(x0, x1), max(x0, x1)
    if hi <= lo:
        return
    cut = min(max(SKY_LIMIT_DEG2, lo), hi)
    if cut > lo:
        ax.plot([lo, cut], [y, y], color=INK, solid_capstyle="round", **kw)
    if hi > cut:
        ax.plot([cut, hi], [y, y], color=INK, alpha=DIM,
                solid_capstyle="round", **kw)


def _alias_f350(path, stem, band):
    """Also write `<stem>_f350.png` for the band whose suffix is empty.

    f350 is the default band, so its files are `hi_kmodes.png` and
    `hi_sufficiency.png`, but the draft's \\includegraphics asks for
    `hi_kmodes_f350.png`. Writing both keeps either spelling working rather
    than silently leaving a missing figure at build time.
    """
    import shutil
    if band.suffix:
        return path
    alias = path.replace(f"{stem}.", f"{stem}_f350.")
    shutil.copyfile(path, alias)
    return f"{path} (+ {os.path.basename(alias)})"


def _panel_kmodes(ax, res):
    """Where the measured modes sit in angular scale."""
    band = res["_band"]
    P._tidy(ax); ax.set_xscale("log"); ax.set_yscale("log")
    for s_ in [k for k in band.strategies if "k" in res[k]]:
        # same colour and marker per strategy as every other figure, so a
        # colour means one thing throughout the paper
        st = P.PAPER_STRAT[s_]
        r = res[s_]
        dens = 100.0 * r["n"] / r["area"]         # per 100 deg^2, the currency
        ax.plot(r["k"], dens, color=st["color"], ls=st["ls"],
                marker=st["marker"], ms=6, mfc=SURFACE, mew=1.6, zorder=4,
                label=f"{st['label']} — "
                      f"{100 * r['total'] / r['area']:.1f} total")
    ax.legend(loc="upper left", fontsize=10.2)
    k_beam = 2 * np.pi / band.fwhm_rad
    ax.axvline(k_beam, color=INK, lw=1.4, ls=(0, (4, 3)), zorder=3)
    ax.annotate("beam", (k_beam, ax.get_ylim()[0]),
                textcoords="offset points", xytext=(-6, 8), ha="right",
                va="bottom", color=INK, fontsize=9.9, rotation=90)
    ax.set_xlabel("$|k_\\perp|$  [rad$^{-1}$]   $\\simeq \\ell$")
    ax.set_ylabel("measured modes per 100 deg$^2$, per bin")
    sec = ax.secondary_xaxis(
        "top", functions=(lambda k: k / band.d_comoving_mpc,
                          lambda kp: kp * band.d_comoving_mpc))
    sec.set_xlabel(f"$k_\\perp$  [Mpc$^{{-1}}$]  at $z = {band.z_centre:.2f}$",
                   fontsize=10.6)


def _panel_sufficiency(ax, res, xlim=None):
    """What those modes are enough for -- cross-correlation, and auto."""
    band = res["_band"]
    P._tidy(ax); ax.set_xscale("log"); ax.grid(axis="y", visible=False)
    # top to bottom. Only the cross rows get an area; see sufficiency().
    lo, mid, hi = band.strategies[0], band.strategies[1], band.strategies[-1]
    rows = [(hi, "cross"), (mid, "cross"), (lo, "cross"),
            (hi, "auto"), (lo, "auto")]
    ys = np.arange(len(rows))[::-1]
    npg_best, npg_shot = NPG[0], NPG[-1]
    q_lo, q_hi = min(SNR), max(SNR)
    for y, (s_, kind) in zip(ys, rows):
        st = P.MONO_STRAT[s_]
        if kind == "auto":
            R = res[s_]["suff"]["auto"]["bias_ratio"]
            # no area, by construction: the floor is a bias here. The text
            # stays short so it cannot reach either vertical sky line --
            # the caption carries "a bias, so no sky area suffices".
            ax.annotate(f"floor/HI power = {_sig3(R)}",
                    ((xlim or band.x_lim)[0] * band.auto_text_x
                     / band.x_lim[0], y),
                        va="center", ha="left", color=INK, fontsize=10.6,
                        zorder=6)
            continue

        cr = res[s_]["suff"]["cross"]
        area_of = lambda npg, q: cr[float(npg)][int(q)]["area"]
        # ideal tracer out to 5 sigma, then the shot-noise case as a dotted
        # extension -- it costs a factor ~2 in area, not an order of magnitude
        # min/max so the bar still reaches every marker when the simulated
        # patch is LARGER than the 1-sigma area, which per-k_par it can be
        _reach_segment(ax, min(res[s_]["area"], area_of(npg_best, q_lo)),
                       max(res[s_]["area"], area_of(npg_best, q_hi)), y,
                       lw=2.2, zorder=3)
        _reach_segment(ax, area_of(npg_best, q_hi), area_of(npg_shot, q_hi), y,
                       lw=1.4, ls=(0, (1.4, 1.6)), zorder=2)
        ax.plot([res[s_]["area"]], [y], marker=st["marker"], ms=8.5,
                mfc=SURFACE, mec=INK, mew=1.8, zorder=5)
        for q in SNR:
            a = area_of(npg_best, q)
            ax.plot([a], [y], marker=st["marker"], ms=7.5, color=INK,
                    alpha=1.0 if a <= SKY_LIMIT_DEG2 else DIM, zorder=6)
            if s_ == hi:                 # label the scale once, on the top row
                ax.annotate(f"{q}$\\sigma$", (a, y),
                            textcoords="offset points", xytext=(0, 11),
                            ha="center", va="bottom", color=INK,
                            alpha=1.0 if a <= SKY_LIMIT_DEG2 else DIM,
                            fontsize=10.2)
        a1 = area_of(npg_best, q_lo)
        # a label centred on a marker beside the reachable-sky line runs
        # straight through it; grow it away from the line instead
        _near = 0.3 < a1 / SKY_LIMIT_DEG2 <= 1.0
        _past = 1.0 < a1 / SKY_LIMIT_DEG2 < 3.0
        ax.annotate(f"{_sig3(a1)} deg$^2$", (a1, y),
                    textcoords="offset points",
                    xytext=((-4, -13) if _near
                            else (8, -13) if _past else (0, -13)),
                    ha="right" if _near else "left" if _past else "center",
                    va="top", color=INK, fontsize=10.6,
                    alpha=1.0 if a1 <= SKY_LIMIT_DEG2 else DIM,
                    fontweight="semibold" if a1 <= SKY_LIMIT_DEG2 else "normal")
        # The raster is the only row whose 3-sigma point also lands inside the
        # reachable sky, so it is the only one where that number is worth
        # printing. It sits 0.95 decades from the 1-sigma label, which is not
        # enough clearance at this font size, hence the leader.
        if s_ == hi and Q_LABEL in SNR:
            a3 = area_of(npg_best, Q_LABEL)
            near = a3 <= SKY_LIMIT_DEG2
            ax.annotate(f"{_sig3(a3)} deg$^2$", (a3, y),
                        textcoords="offset points",
                        xytext=(34, -34) if not res.get("_kpar") else (26, -20),
                        ha="left", va="top", color=INK, fontsize=10.6,
                        alpha=1.0 if near else DIM,
                        fontweight="semibold" if near else "normal",
                        arrowprops=dict(arrowstyle="-", lw=0.9, color=INK,
                                        alpha=0.75, shrinkA=1, shrinkB=5))

    # both sky lines are labelled in a band above the rows, clear of the sigma
    # labels that sit just above the raster markers
    y_sky = len(rows) + 0.46
    ax.axvline(SKY_LIMIT_DEG2, color=INK, lw=1.6, zorder=4)
    ax.annotate(f"reachable from Karoo\n($\\delta \\leq {_DEC_MAX:.1f}"
                f"\\degree$, {SKY_LIMIT_DEG2:,.0f} deg$^2$)",
                (SKY_LIMIT_DEG2, y_sky), textcoords="offset points",
                xytext=(-9, 0), ha="right", va="center", color=INK,
                fontsize=10.0, linespacing=1.4)
    ax.axvline(FULL_SKY_DEG2, color=INK, lw=1.2, ls=(0, (4, 3)), alpha=DIM,
               zorder=3)
    ax.annotate("$4\\pi$", (FULL_SKY_DEG2, y_sky), textcoords="offset points",
                xytext=(8, 0), ha="left", va="center", color=INK, alpha=0.6,
                fontsize=10.0)

    ax.set_yticks(ys)
    ax.set_yticklabels(
        [ROW_LABEL[s_] + (" × ideal tracer" if kind == "cross"
                          else ", HI auto-spectrum") for s_, kind in rows],
        fontsize=10.6)
    ax.set_ylim(-1.15, len(rows) + 0.85)
    ax.set_xlim(xlim or band.x_lim)
    ax.set_xlabel("sky area needed  [deg$^2$]")
    ax.legend(handles=[
        Line2D([], [], ls="none", marker="o", ms=8.5, mfc=SURFACE, mec=INK,
               mew=1.8, label="open marker = simulated patch"),
        Line2D([], [], color=INK, lw=2.2,
               label="ideal tracer ($r = 1$, no shot noise)"),
        Line2D([], [], color=INK, lw=1.4, ls=(0, (1.4, 1.6)),
               label="$\\bar{n}P_\\mathrm{g} = 1$")],
        loc="lower right", fontsize=9.8, framealpha=0.95,
        title="$\\Delta P/P = 1/\\mathrm{SNR}$", title_fontsize=10.0)


def fig_kmodes(res, text=None):
    """Where the measured modes live in angular scale.

    Split from the sufficiency panel on 2026-10-01: the two ask different
    questions, and the sufficiency half needed the room -- five rows, three
    detection thresholds and two tracer assumptions do not fit beside a
    log-log panel. ``fig_sufficiency`` is the other half.
    """
    band = res["_band"]
    fig, ax = plt.subplots(
        figsize=P._figsize("hi_kmodes" + band.suffix, (7.6, 5.4)))
    _panel_kmodes(ax, res)
    '''P._title(fig, "Where the measured modes are, by angular scale",
             "tr($WA$) split by angular scale — the k sum reproduces the "
             "14.6 and 22.3 per 100 deg$^2$ of Table 2 exactly. The raster's "
             "advantage is close to scale-independent.",
             x=0.012, y_title=1.05, y_sub=0.975, override=text)'''
    fig.tight_layout()      # rect=(0, 0, 1, 0.87) if the title is restored
    plt.show()
    out = P._save(fig, "hi_kmodes" + band.suffix)
    return _alias_f350(out, "hi_kmodes", band)


def fig_sufficiency(res, text=None):
    """How much sky each strategy would need, and what that buys.

    The other half of the old two-panel hi_kmodes (split 2026-10-01).
    """
    band = res["_band"]
    fig, ax = plt.subplots(
        figsize=P._figsize("hi_sufficiency" + band.suffix, (11.0, 5.8)))
    _panel_sufficiency(ax, res, xlim=X_SINGLE_KPAR if res.get("_kpar") else None)
    '''P._title(fig, "How many modes is enough? Only the cross-correlation asks a "
             "reachable number",
             "Against an external tracer the floor is uncorrelated foreground and "
             "costs variance only; in the Hi auto-spectrum it is a BIAS, which "
             "no sky area averages away.",
             x=0.012, y_title=1.05, y_sub=0.975, override=text)'''
    fig.tight_layout()      # rect=(0, 0, 1, 0.87) if the title is restored
    plt.show()
    out = P._save(fig, "hi_sufficiency" + band.suffix)
    return _alias_f350(out, "hi_sufficiency", band)


# ---------------------------------------------------------------------------
# Both bands in one figure
# ---------------------------------------------------------------------------

BAND_GROUPS = (("f350", "350$-$400 MHz", "$z \\simeq 2.8$"),
               ("f700", "675$-$725 MHz", "$z \\simeq 1.0$"))
ROW_ORDER = ("raster", "drift12", "drift")     # top to bottom, within a group
_ALIAS = {"drift12": ("drift12", "drift12m")}  # the 12-strip row is named per band


def _row_strategy(band, row):
    """The strategy in ``band`` that fills ``row`` -- drift12 vs drift12m."""
    for name in _ALIAS.get(row, (row,)):
        if name in band.strategies:
            return name
    return None


def _reach_colour(ax, x0, x1, y, colour, **kw):
    """One horizontal run, full colour inside the reachable sky, faded beyond.

    Faded rather than greyed, so the colour still identifies the strategy past
    the limit -- which is where most of the rows live.
    """
    lo, hi = min(x0, x1), max(x0, x1)
    if hi <= lo:
        return
    cut = min(max(SKY_LIMIT_DEG2, lo), hi)
    if cut > lo:
        ax.plot([lo, cut], [y, y], color=colour, solid_capstyle="round", **kw)
    if hi > cut:
        ax.plot([cut, hi], [y, y], color=colour, alpha=DIM,
                solid_capstyle="round", **kw)


def fig_sufficiency_bands(results, text=None, xlim=None, patch_markers=True):
    """Both bands' sufficiency rows in one panel, six rows in two groups.

    ``results`` maps band tag -> the dict ``compute``/``load`` returns. The
    auto-spectrum rows are gone: they carry no area by construction, so they
    were two rows saying "not applicable". Their floor-to-Hi power ratios go in
    the caption instead, and ``auto_ratios`` returns them.
    """
    groups = [(BANDS[t], lab, z) for t, lab, z in BAND_GROUPS if t in results]
    if len(groups) < 2:
        return None
    nrow = len(ROW_ORDER)
    fig, ax = plt.subplots(
        figsize=P._figsize("hi_sufficiency_bands", (11.4, 6.6)))
    P._tidy(ax); ax.set_xscale("log"); ax.grid(axis="y", visible=False)

    # upper group occupies y = 6,5,4; lower y = 2,1,0; the rule sits at y = 3
    top_of = lambda gi: (len(groups) - 1 - gi) * (nrow + 1) + nrow - 1
    npg_best, npg_shot = NPG[0], NPG[-1]
    q_lo, q_hi = min(SNR), max(SNR)
    yticks, ylabels = [], []

    for gi, (band, glabel, zlabel) in enumerate(groups):
        res = results[band.tag]
        y_top = top_of(gi)
        for ri, row in enumerate(ROW_ORDER):
            s_ = _row_strategy(band, row)
            if s_ is None:
                continue
            y = y_top - ri
            st = P.PAPER_STRAT[row]
            col = st["color"]
            yticks.append(y); ylabels.append(st["label"])
            cr = res[s_]["suff"]["cross"]
            area_of = lambda npg, q: cr[float(npg)][int(q)]["area"]

            # With the patch markers drawn the bar has to reach them, which
            # per-k_par means reaching LEFT of 1 sigma. Without them it runs
            # from 1 sigma outwards, as it did before they set its left end.
            lo_end = (min(res[s_]["area"], area_of(npg_best, q_lo))
                      if patch_markers else area_of(npg_best, q_lo))
            hi_end = (max(res[s_]["area"], area_of(npg_best, q_hi))
                      if patch_markers else area_of(npg_best, q_hi))
            _reach_colour(ax, lo_end, hi_end, y, col, lw=2.4, zorder=3)
            _reach_colour(ax, area_of(npg_best, q_hi),
                          area_of(npg_shot, q_hi), y, col,
                          lw=1.5, ls=(0, (1.4, 1.6)), zorder=2)
            if patch_markers:
                ax.plot([res[s_]["area"]], [y], marker=st["marker"], ms=9,
                        mfc=SURFACE, mec=col, mew=2.0, zorder=5,
                        linestyle="none")
            for q in SNR:
                a = area_of(npg_best, q)
                ax.plot([a], [y], marker=st["marker"], ms=8, color=col,
                        alpha=1.0 if a <= SKY_LIMIT_DEG2 else DIM, zorder=6,
                        linestyle="none")
                # sigma levels labelled once, over the top group's raster row
                if gi == 0 and ri == 0:
                    ax.annotate(f"{q}$\\sigma$", (a, y),
                                textcoords="offset points", xytext=(0, 11),
                                ha="center", va="bottom", color=INK,
                                alpha=1.0 if a <= SKY_LIMIT_DEG2 else DIM,
                                fontsize=10.2)
            # numbers: the 1 sigma area on every row, plus the headline 3 sigma
            labels = [(q_lo, area_of(npg_best, q_lo))]
            if gi == 0 and ri == 0:
                labels.append((3, area_of(npg_best, 3)))
            for q, a in labels:
                near = a <= SKY_LIMIT_DEG2
                # A label centred on a marker that sits just inside the
                # reachable-sky line would run straight through it, so anything
                # within half a decade of the line is right-aligned instead and
                # grows away from it.
                close = 0.3 < a / SKY_LIMIT_DEG2 <= 1.0
                # ...and one just ABOVE it runs through it just as badly, so
                # it is left-aligned and grows away to the right instead
                past = 1.0 < a / SKY_LIMIT_DEG2 < 3.0
                ax.annotate(f"{_sig3(a)}", (a, y), textcoords="offset points",
                            xytext=((-4, -13) if close
                                    else (8, -13) if past else (0, -13)),
                            ha="right" if close
                            else "left" if past else "center", va="top",
                            color=col, fontsize=10.2,
                            alpha=1.0 if near else DIM,
                            fontweight="semibold" if near else "normal")

        ax.annotate(f"{glabel}    {zlabel}", ((xlim or X_BANDS)[0], y_top + 0.62),
                    textcoords="offset points", xytext=(8, 0), ha="left",
                    va="center", color=INK, fontsize=11.6,
                    fontweight="semibold")

    if len(groups) > 1:
        # midway between the bottom row of the upper group and the top of the
        # lower one -- NOT through a row
        ax.axhline(0.5 * ((top_of(0) - nrow + 1) + top_of(1)),
                   color=INK, lw=0.8, alpha=0.45, zorder=1)

    ax.axvline(SKY_LIMIT_DEG2, color=INK, lw=1.6, zorder=4)
    ax.annotate(f"reachable from Karoo\n($\\delta \\leq {_DEC_MAX:.1f}"
                f"\\degree$, {SKY_LIMIT_DEG2:,.0f} deg$^2$)",
                (SKY_LIMIT_DEG2, top_of(0) + 1.35),
                textcoords="offset points", xytext=(-9, 0), ha="right",
                va="center", color=INK, fontsize=10.0, linespacing=1.4)
    ax.axvline(FULL_SKY_DEG2, color=INK, lw=1.2, ls=(0, (4, 3)), alpha=DIM,
               zorder=3)
    ax.annotate("$4\\pi$", (FULL_SKY_DEG2, top_of(0) + 1.35),
                textcoords="offset points", xytext=(8, 0), ha="left",
                va="center", color=INK, alpha=0.6, fontsize=10.0)

    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=10.6)
    ax.set_ylim(-1.35, top_of(0) + 1.95)
    ax.set_xlim(xlim or X_BANDS)
    ax.set_xlabel("sky area needed  [deg$^2$]")
    handles = [Line2D([], [], ls="none", marker="o", ms=9, mfc=SURFACE,
                      mec=INK, mew=2.0, label="open marker = simulated patch")]
    if not patch_markers:
        handles = []                   # the caption carries the patch sizes
    ax.legend(handles=handles + [
        Line2D([], [], color=INK, lw=2.4,
               label="ideal tracer ($r = 1$, no shot noise)"),
        Line2D([], [], color=INK, lw=1.5, ls=(0, (1.4, 1.6)),
               label="$\\bar{n}P_\\mathrm{g} = 1$")],
        loc="lower right", fontsize=9.8, frameon=False,
        title="$\\Delta P/P = 1/\\mathrm{SNR}$", title_fontsize=10.0)
    fig.tight_layout()
    plt.show()
    return P._save(fig, "hi_sufficiency_bands")


# Expected areas for the merged figure, deg^2. Exact values are asserted to
# +/-1 deg^2; the ones written approximately in the spec are checked to 1%.
# The merge must not change any number -- that is the whole point of the check.
BANDS_EXPECT = {
    ("f350", "raster"):   ((1, 2419, 1.0), (3, 21767, 1.0), (5, 60470, "1%")),
    ("f350", "drift12"):  ((1, 2.54e5, "1%"),),
    ("f350", "drift"):    ((1, 1.09e6, "1%"),),
    ("f700", "raster"):   ((1, 1058, 1.0), (3, 9519, 1.0), (5, 26450, "1%")),
    ("f700", "drift12m"): ((1, 516, 1.0),),
    ("f700", "drift"):    ((1, 655, 1.0),),
}


def check_bands(results, strict=True):
    """Assert the merged figure reproduces the two single-band panels."""
    print("\nmerged-figure check (areas in deg^2, ideal tracer):", flush=True)
    bad = []
    for (tag, s_), checks in sorted(BANDS_EXPECT.items()):
        if tag not in results or s_ not in results[tag]:
            continue
        cr = results[tag][s_]["suff"]["cross"][float(NPG[0])]
        for q, want, tol in checks:
            got = cr[int(q)]["area"]
            ok = (abs(got - want) <= 0.01 * want if tol == "1%"
                  else abs(got - want) <= tol)
            mark = "ok " if ok else "** MISMATCH **"
            print(f"  {tag} {s_:9s} {q}sigma  got {got:11.1f}  want {want:11.1f}"
                  f"   {mark}", flush=True)
            if not ok:
                bad.append((tag, s_, q, got, want))
    if bad and strict:
        raise AssertionError(
            "merged figure changed numbers that must not move: "
            + "; ".join(f"{t} {s} {q}sigma {g:.1f} != {w:.1f}"
                        for t, s, q, g, w in bad))
    print("  all values reproduce the single-band panels", flush=True)
    return not bad


def auto_ratios(results):
    """floor/HI POWER per strategy per band -- the auto rows, for the caption."""
    out = {}
    for tag, _, _ in BAND_GROUPS:
        if tag not in results:
            continue
        band = BANDS[tag]
        out[tag] = {s: results[tag][s]["suff"]["auto"]["bias_ratio"]
                    for s in band.strategies}
    return out


def compute(band=None, kedges=None):
    """The per-strategy decomposition for one band.

    Needs the cached operators in ``hi_cache/``, so at 675-725 MHz this runs on
    Ilifu, not here. Everything the figures need is written to the band's npz,
    and ``load()`` reads it back, so the plotting half runs anywhere.
    """
    band = BANDS[band] if isinstance(band, str) else (band or BANDS[DEFAULT_BAND])
    if kedges is None:
        kedges = np.geomspace(8.0, 1.25 * 2 * np.pi / band.fwhm_rad, 13)
    floors = floor_amp(band)
    res = {"_kedges": kedges, "_band": band, "_floor_amp": floors}
    for s in band.strategies:
        path = operator_path(band, s)
        if os.path.exists(path):
            with open(path, "rb") as f:
                mm = pickle.load(f)
            if prior.PINNED:
                import ska_hi_analysis as A
                prior.set_strategy(s, A.channel_operators(s))
            k, n_of_k, tot = n_eff_of_k(mm, kedges, band)
            area = len(np.asarray(mm.pixel_indices)) * hp.nside2pixarea(
                band.nside, degrees=True)
            res[s] = dict(k=k, n=n_of_k, total=tot, area=area)
            print(f"{s:8s} tr(WA)={tot:7.2f}  sum over k={n_of_k.sum():7.2f}  "
                  f"area={area:6.1f} deg^2", flush=True)
        else:
            # the matched-density row: twelve strips reach the raster's mode
            # density and still carry a much larger floor, which is the point
            # of it. No operator cached, so no k decomposition -- right panel
            # only (true at 350 MHz; at 700 the operator exists).
            tot, area = ladder_geometry(band)
            res[s] = dict(total=tot, area=area)
            print(f"{s:8s} tr(WA)={tot:7.2f}  (ladder sweep, no k "
                  f"decomposition)   area={area:6.1f} deg^2", flush=True)
        res[s]["suff"] = sufficiency(floors[s], res[s]["total"], res[s]["area"])
    return res


def _npz_path(band):
    return os.path.join(_HERE, "results", band.out)


def save(res):
    """Everything the two figures need, so they can be drawn without operators."""
    band = res["_band"]
    strat = [k for k in res if not k.startswith("_")]
    out = dict(kedges=res["_kedges"], d_comoving_mpc=band.d_comoving_mpc,
               band_tag=band.tag, freq_mhz=band.freq_mhz, nside=band.nside,
               z_centre=band.z_centre, strategies=np.array(strat, dtype=object),
               suff_snr=np.asarray(SNR, float),
               suff_npg=np.asarray([float(n) for n in NPG]),
               r_cross=float(R_CROSS), sky_limit_deg2=SKY_LIMIT_DEG2,
               full_sky_deg2=FULL_SKY_DEG2)
    for s in strat:
        su = res[s]["suff"]
        out[f"{s}_total"] = res[s]["total"]
        out[f"{s}_area"] = res[s]["area"]
        out[f"{s}_floor_amp"] = res["_floor_amp"][s]
        out[f"{s}_bias_ratio"] = su["auto"]["bias_ratio"]
        out[f"{s}_auto_lower_bound_area"] = (
            su["auto"]["auto_lower_bound_if_floor_known"]["area"])
        # (len(NPG), len(SNR)), rows ordered as suff_npg, columns as suff_snr
        out[f"{s}_cross_area"] = np.asarray(
            [[su["cross"][float(n)][int(q)]["area"] for q in SNR] for n in NPG])
        if "k" in res[s]:
            out[f"{s}_k"], out[f"{s}_n"] = res[s]["k"], res[s]["n"]
    np.savez(_npz_path(band), **out)
    print(f"wrote results/{band.out}", flush=True)
    return _npz_path(band)


def load(band=None, kpar=False):
    """Rebuild what ``compute`` returned, from the npz. No operators needed.

    ``kpar=True`` swaps the single-bin ``N_m = N_eff`` for the per-k_par sum
    (``ska_nm_kpar``), which reads only the band's ablation and experiment
    npz -- still no operators, still nothing re-solved.
    """
    band = BANDS[band] if isinstance(band, str) else (band or BANDS[DEFAULT_BAND])
    path = _npz_path(band)
    if not os.path.exists(path):
        sys.exit(f"need results/{band.out} -- run `compute` for band "
                 f"{band.tag} where the operators are cached, then pull it")
    d = np.load(path, allow_pickle=True)
    floors = {s: float(d[f"{s}_floor_amp"]) for s in band.strategies}
    res = {"_kedges": d["kedges"], "_band": band, "_floor_amp": floors,
           "_kpar": bool(kpar)}
    terms = {}
    if kpar:
        import ska_nm_kpar as NK
        terms = NK.kpar_terms(band.tag, nmodes=NMODES)
    for s in band.strategies:
        res[s] = dict(total=float(d[f"{s}_total"]), area=float(d[f"{s}_area"]))
        if f"{s}_k" in d:
            res[s]["k"], res[s]["n"] = d[f"{s}_k"], d[f"{s}_n"]
        res[s]["suff"] = sufficiency(floors[s], res[s]["total"],
                                     res[s]["area"], kpar=terms.get(s))
    return res


def report(res):
    """The printed checks. Numbers here are what the figure draws."""
    band = res["_band"]
    print(f"\n[{band.tag}] floor/HI amplitude at {NMODES} modes, common "
          f"resolution, from {band.ablation}:", flush=True)
    for s in band.strategies:
        print(f"  {ROW_LABEL[s]:16s} {res['_floor_amp'][s]:8.1f}x", flush=True)

    print(f"\ncross-power detection, r = {R_CROSS:g}: sky area [deg^2] needed "
          f"for dP/P = 1/snr", flush=True)
    print("  " + "row".ljust(16) + "nPg".rjust(5)
          + "".join(f"{q} sigma".rjust(13) for q in SNR), flush=True)
    for s in reversed(band.strategies):
        for n in NPG:
            cells = "".join(
                f"{res[s]['suff']['cross'][float(n)][int(q)]['area']:13.4g}"
                for q in SNR)
            tag = "inf" if np.isinf(n) else f"{n:g}"
            print(f"  {ROW_LABEL[s]:16s}{tag:>5s}{cells}", flush=True)

    print("\nHi auto-spectrum: the floor is a bias, so no area is quoted",
          flush=True)
    for s in (band.strategies[-1], band.strategies[0]):
        a = res[s]["suff"]["auto"]
        print(f"  {ROW_LABEL[s]:16s} floor/HI power = {a['bias_ratio']:10.4g}"
              f"   (lower bound IF the floor were known and subtracted: "
              f"{a['auto_lower_bound_if_floor_known']['area']:.3g} deg^2)",
              flush=True)

    reach = [s for s in band.strategies
             if res[s]["suff"]["cross"][float(NPG[0])][int(min(SNR))]["area"]
             <= SKY_LIMIT_DEG2]
    print(f"\ninside the reachable sky at 1 sigma: "
          f"{', '.join(ROW_LABEL[s] for s in reach) or 'none'}", flush=True)

    if band.tag == "f350" and not res.get("_kpar"):
        got = res["raster"]["suff"]["cross"][float(NPG[0])][int(min(SNR))]["area"]
        rel = abs(got / DRAFT_RASTER_1SIG_DEG2 - 1.0)
        assert rel < 0.01, (f"raster 1-sigma cross area moved from "
                            f"{DRAFT_RASTER_1SIG_DEG2:.0f} to {got:.1f} deg^2 "
                            f"({100 * rel:.1f}%) -- the draft quotes the old value")
        print(f"check: raster 1 sigma, ideal tracer = {got:.1f} deg^2, "
              f"{100 * rel:.2f}% from the {DRAFT_RASTER_1SIG_DEG2:.0f} "
              "in the draft", flush=True)


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--band", default=DEFAULT_BAND, choices=sorted(BANDS),
                    help="which observing band to build")
    ap.add_argument("--draw-only", action="store_true",
                    help="skip the operators and draw from the cached npz")
    ap.add_argument("--bands", action="store_true",
                    help="draw the MERGED two-band sufficiency figure "
                         "(hi_sufficiency_bands) from both bands' cached npz, "
                         "and assert it reproduces the single-band panels. "
                         "Reads results only -- no operators, no cluster.")
    ap.add_argument("--pin-prior", action="store_true",
                    help="use the band's amplitude-only prior (ska_hi_prior), "
                         "matching the ablations the floors are read from. "
                         "Without it n_eff is computed at each strategy's own "
                         "patch variance while the floors are pinned, which is "
                         "inconsistent.")
    ap.add_argument("--kpar", action="store_true",
                    help="count the cube's LINE-OF-SIGHT degrees of freedom "
                         "too, instead of N_m = N_eff (one channel's angular "
                         "modes). Each k_par bin carries its own floor/Hi "
                         "ratio and its own transfer function; see "
                         "ska_nm_kpar. Changes every area, so the single-bin "
                         "regression guards are skipped.")
    ap.add_argument("--no-figures", action="store_true",
                    help="compute and save the npz, draw nothing. For Ilifu, "
                         "which has no figures/ directory and no need for one; "
                         "pull the npz and redraw here with --draw-only")
    args = ap.parse_args(argv)
    band = BANDS[args.band]

    if args.bands:
        results = {t: load(BANDS[t], kpar=args.kpar) for t, _, _ in BAND_GROUPS}
        if args.kpar:
            print("per-k_par N_m: the single-bin regression guards do not "
                  "apply and are skipped", flush=True)
        else:
            check_bands(results)
        print("\nauto floor/HI POWER ratios (for the caption):", flush=True)
        for tag, d in auto_ratios(results).items():
            print(f"  {tag}: "
                  + ", ".join(f"{k} {v:,.0f}" for k, v in d.items()), flush=True)
        print("wrote", fig_sufficiency_bands(
            results, xlim=X_BANDS_KPAR if args.kpar else None,
            patch_markers=not args.kpar), flush=True)
        return results

    if args.pin_prior:
        prior.PINNED = True
    res = load(band, kpar=args.kpar) if args.draw_only else compute(band)
    if not args.draw_only:
        save(res)
    report(res)
    if not args.no_figures:
        print("wrote", fig_kmodes(res), flush=True)
        print("wrote", fig_sufficiency(res), flush=True)
    return res


if __name__ == "__main__":
    main()
