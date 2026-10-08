"""One prior width per band, shared by every strategy.

The solve uses a flat isotropic prior ``S^-1 = I / sigma_p^2``. Until
2026-10-02 every site set ``sigma_p`` from the sky variance over *its own
strategy's patch*, which makes the prior a function of the instrument: a wider
footprint picks up more sky structure, gets a looser prior, and is penalised
for it. Job 13990883 measured the cost at 675-725 MHz -- pinning every strategy
to a common ``sigma_p`` improved the twelve-strip drift by 7.3% and the raster
by 11.1%, while the reference strategy reproduced its unpinned numbers exactly.

That also broke the paper's own stated condition: Sec. "Accessible sky modes"
says modes per unit sky is only a fair comparison at *fixed* ``sigma_p``, and
``tr(WA) = sum lambda_i / (lambda_i + sigma_p^-2)`` moves with it by
construction.

**Definition (amplitude-only, 2026-10-05).** Each strategy keeps its OWN
per-channel ``sigma_p(nu)``, rescaled by a single constant:

    sigma_p,s^2(nu) = sigma_p,s^2,own(nu) * c_s,
    c_s = geometric mean over channels of sigma_p,ref^2(nu) / sigma_p,s^2,own(nu)

with the reference the three-strip drift's patch. ``c_ref = 1`` exactly, so the
reference strategy is untouched and remains the control.

**Why amplitude only.** The confound -- wider patch, looser prior -- is purely
an amplitude effect. The spectral SHAPE of ``sigma_p(nu)`` is a property of the
sky over that patch, and overwriting it with another patch's shape imposes a
scale that varies across the band (0.6142-0.6908 at 350-400, a 12.5% swing).
Because the Wiener gain ``lambda_i / (lambda_i + sigma_p^-2)`` depends on
``lambda_i``, which varies channel to channel as the beam narrows 12.5% across
that band, a frequency-varying scale changes WHICH modes sit near the
data/prior crossover per channel. That injects channel-to-channel structure
into the solution. Measured, noise power before any cleaning, pinned/own:

    shape-pinned   low k 0.68   mid k 13.54   high k 132.70
    amplitude-only low k 0.62   mid k  0.61   high k   0.58

Both remove the same total noise (rms -21.3% vs -21.7%); only the shape-pin
redistributes it. The per-mode bound lambda/(lambda + sigma_p^-2)^2 constrains
total noise power, not its distribution in k_par, so there is nothing in it to
forbid the redistribution.

Why that region rather than the common patch: the common patch is the
intersection over whichever strategies are in the run, so it changes when the
strategy list changes and is not a property of the band. The three-strip drift
is a fixed, named geometry. At 675-725 MHz the two agree to 0.6% anyway,
because the drift is the narrowest strategy, so ``pin700`` stands as measured.

The reference lives in ``results/prior_ref_f{lo}_{hi}_ns{nside}.npz`` and every
run in that band reads it. Nothing recomputes sigma_p from its own patch once
the file exists.

    python ska_hi_prior.py            # build it for the current band
    python ska_hi_prior.py --show     # print it
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

REF_STRATEGY = "drift"          # the three-strip ladder, 52/50/48
_FLOOR = 1e-3                   # the clip every site already applied to std

# OFF by default, so every published number reproduces untouched. A run opts in
# with --pin-prior, which sets this once; every S^-1 site then applies the
# amplitude constant. One switch, honoured everywhere -- so merely building the
# reference file cannot silently move Table 2.
PINNED = False

# c_s for the strategy currently being solved. Callers set it with
# set_strategy() at the top of each strategy's loop; 1.0 means "unchanged",
# which is what the reference strategy gets.
AMP_SCALE = 1.0
_amp_cache = {}
_cache = {}
_warned = set()


def ref_path(f_lo=None, f_hi=None, nside=None):
    f_lo = X.F_LO_MHZ if f_lo is None else f_lo
    f_hi = X.F_HI_MHZ if f_hi is None else f_hi
    nside = X.NSIDE if nside is None else nside
    return os.path.join(_HERE, "results",
                        f"prior_ref_f{f_lo:.0f}_{f_hi:.0f}_ns{nside}.npz")


def build(ops=None, verbose=True):
    """Compute and save the band's reference sigma_p, from the drift's patch.

    ``ops`` is the reference strategy's per-channel operator list if the caller
    already has it (the HI runs do, so this costs nothing extra). Otherwise it
    is loaded, which builds and caches the operators if they do not exist.
    """
    import ska_hi_analysis as A
    if ops is None:
        ops = A.channel_operators(REF_STRATEGY)
    freqs = X.channel_freqs()
    if len(ops) != len(freqs):
        sys.exit(f"{REF_STRATEGY} has {len(ops)} operators for {len(freqs)} "
                 "channels -- band mismatch")
    sig2, npix = [], []
    for f, mm in zip(freqs, ops):
        pix = np.asarray(mm.pixel_indices)
        truth = gdsm_equatorial_sky_model(freq=f, nside=X.NSIDE)[pix]
        sig2.append(max(float(np.std(truth)), _FLOOR) ** 2)
        npix.append(len(pix))
    out = dict(
        freq_mhz=np.asarray(freqs, float), sigma_p2=np.asarray(sig2, float),
        npix=np.asarray(npix, int), ref_strategy=REF_STRATEGY,
        nside=X.NSIDE, f_lo_mhz=X.F_LO_MHZ, f_hi_mhz=X.F_HI_MHZ,
        definition=("sigma_p^2 = var(GDSM truth) over the three-strip drift's "
                    "own patch, per channel. One value per channel, shared by "
                    "every strategy in this band."))
    path = ref_path()
    np.savez(path, **out)
    if verbose:
        s = np.sqrt(out["sigma_p2"])
        print(f"[prior] wrote {os.path.basename(path)}: sigma_p "
              f"{s.min():.4f}-{s.max():.4f} K over {npix[0]}-{npix[-1]} px "
              f"({REF_STRATEGY})", flush=True)
    return out


def load(f_lo=None, f_hi=None, nside=None, required=False):
    """The band's reference, or None if it has not been built."""
    path = ref_path(f_lo, f_hi, nside)
    if path in _cache:
        return _cache[path]
    if not os.path.exists(path):
        if required:
            sys.exit(f"need {os.path.relpath(path, _HERE)} -- build it with "
                     "`python ska_hi_prior.py` where the operators are cached")
        return None
    d = np.load(path, allow_pickle=True)
    _cache[path] = dict(freq_mhz=d["freq_mhz"], sigma_p2=d["sigma_p2"],
                        npix=d["npix"], ref_strategy=str(d["ref_strategy"]))
    return _cache[path]


def sigma2(freq_mhz, patch=None, label=""):
    """The prior VARIANCE to use at this frequency.

    With ``PINNED`` off this IS the caller's own patch variance, i.e. exactly
    what every site did before this module existed. With it on, the band
    reference wins; if that file is missing the patch is used anyway, but
    loudly and once per band, since a silent fallback is how the confound
    survived this long.
    """
    if patch is None:
        sys.exit("ska_hi_prior.sigma2 needs the patch: the amplitude-only "
                 "prior rescales the patch's OWN variance")
    own = max(float(np.std(patch)), _FLOOR) ** 2
    if not PINNED:
        return own
    if AMP_SCALE == 1.0 and label and label not in _warned:
        _warned.add(label)
        print(f"   [prior] {label}: c = 1.000000 (reference, or set_strategy "
              "not called)", flush=True)
    return own * AMP_SCALE


def s_inv(freq_mhz, patch=None, label=""):
    """``S^-1`` -- what the mode counts actually divide by."""
    return 1.0 / sigma2(freq_mhz, patch, label)


def set_strategy(strategy, ops=None, verbose=True):
    """Set ``AMP_SCALE`` to this strategy's ``c_s`` for the current band.

    Call once at the top of each strategy's loop. The reference strategy gets
    exactly 1.0, so it is bit-identical to an unpinned run -- the control.
    """
    global AMP_SCALE
    if not PINNED:
        AMP_SCALE = 1.0
        return AMP_SCALE
    key = (ref_path(), strategy)
    if key not in _amp_cache:
        _amp_cache[key] = 1.0 if strategy == REF_STRATEGY else amp_scale(ops)
    AMP_SCALE = _amp_cache[key]
    if verbose:
        print(f"   [prior] {strategy}: c = {AMP_SCALE:.6f} "
              f"(sigma scaled by {np.sqrt(AMP_SCALE):.6f})", flush=True)
    return AMP_SCALE


def amp_scale(ops):
    """``c_s``: geometric mean over channels of ref variance / own variance."""
    import ska_hi_analysis  # noqa: F401  (import cycle guard)
    ref = load(required=True)
    freqs = X.channel_freqs()
    if ops is None:
        sys.exit("amp_scale needs this strategy's per-channel operators")
    own = np.array([
        max(float(np.std(gdsm_equatorial_sky_model(freq=f, nside=X.NSIDE)
                         [np.asarray(mm.pixel_indices)])), _FLOOR) ** 2
        for f, mm in zip(freqs, ops)])
    return float(np.exp(np.mean(np.log(ref["sigma_p2"] / own))))


def shape_residual(ops):
    """max/min over channels of sigma_p,s(nu)/sigma_p,ref(nu) AFTER rescaling.

    The residual shape difference the amplitude-only prior leaves behind, which
    is the honest thing to quote: it is what we are deliberately NOT removing.
    """
    ref = load(required=True)
    freqs = X.channel_freqs()
    own = np.array([
        max(float(np.std(gdsm_equatorial_sky_model(freq=f, nside=X.NSIDE)
                         [np.asarray(mm.pixel_indices)])), _FLOOR) ** 2
        for f, mm in zip(freqs, ops)])
    r = np.sqrt(own * amp_scale(ops) / ref["sigma_p2"])
    return float(r.min()), float(r.max())


def main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--show", action="store_true",
                    help="print the band's reference instead of building it")
    args = ap.parse_args()
    globals()["PINNED"] = True
    print(f"band {X.F_LO_MHZ:.0f}-{X.F_HI_MHZ:.0f} MHz, nside {X.NSIDE}, "
          f"{X.NCHAN} channels", flush=True)
    ref = load() if args.show else build()
    if ref is None:
        sys.exit(f"{os.path.basename(ref_path())} does not exist")
    s = np.sqrt(ref["sigma_p2"])
    f = ref["freq_mhz"]
    print(f"  reference strategy: {ref.get('ref_strategy', REF_STRATEGY)}")
    for i in (0, len(f) // 2, len(f) - 1):
        print(f"  ch {i:2d}  {f[i]:8.3f} MHz  sigma_p {s[i]:.5f} K  "
              f"S^-1 {1 / ref['sigma_p2'][i]:10.3f}")


if __name__ == "__main__":
    main()
