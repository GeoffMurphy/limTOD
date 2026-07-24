# SKA-like drift-scan strategy

Starting point for simulating an SKA-Mid-like single dish (15 m, Karoo site)
observing in **drift-scan** mode: the telescope is parked at a fixed (Az, El)
and the sky drifts through the beam, with the elevation stepped once per
sidereal night to build up a declination strip.

This complements the strategies in [`examples/DSA/`](../DSA/) (azimuth scan,
stop-and-stare, elevation cascade) — a drift scan is the degenerate case of
the azimuth scan with zero sweep amplitude.

Two beam models run through the identical pipeline so their effect is
isolated: an idealised **Gaussian** (FWHM = 1.22 λ/D ≈ 4.0° at 350 MHz) and a
**tapered-aperture** ("realistic Airy") pattern — circular-aperture
diffraction with a −14 dB edge taper, giving FWHM ≈ 3.8°, sidelobe rings
every ~3.3° (first at −23 dB) and ~1% of beam power beyond 10° off-axis.
Drift scans are the worst case for sidelobes: with the dish parked, off-axis
structure sweeps through the rings coherently and nothing cross-links it away.

| File | Purpose |
|------|---------|
| [ska_drift_scan.ipynb](./ska_drift_scan.ipynb) | End-to-end drift-scan TOD simulation + HPW map-making, Gaussian vs tapered-aperture beam (experiments 001/002) |
| [ska_drift_galplane.ipynb](./ska_drift_galplane.ipynb) | Same design with the LST window shifted onto the galactic plane — sidelobe stress test (experiment 003) |
| [ska_common.py](./ska_common.py) | Site constants, both beam models, equatorial-frame GDSM wrapper, drift-night timing, AltAz→ICRS track |
| `simulated_TODs_ska_drift_{gauss,airy}.npz` | Cached TODs (regenerated if deleted) |
| `mapmaker_ops_ska_drift_{gauss,airy}_ns64.pkl` | Cached map-maker operators (regenerated if deleted) |
| [results/](./results/) | Numbered experiment-record PDFs (`NNN_slug.pdf`), written by `ska_common.save_results_pdf` |

The notebook also runs a **beam-mismatch** scenario — TOD simulated with the
tapered-aperture beam but map made with the Gaussian operator — so the cost
of *unmodelled* sidelobes (the realistic survey case) is separated from the
cost of modelled ones. Noise is seeded per night with the same seed across
beams, so scenario comparisons share identical 1/f and white-noise
realizations.

## Results archive

Each experiment is recorded as a numbered PDF in `results/` (title page with
a short setup description + the residual-RMS table, then one figure per
page), with a companion `.txt` of the same name holding the full
interpretation — findings, caveats, and how to reproduce. Current entries:

- `001_short-drift-matched-beams.{pdf,txt}` — 3-night drift at b ≈ +50°,
  matched forward/inverse beams, Gaussian vs tapered aperture, no-HP vs
  2 mHz HP. HP helps ~10%; matched sidelobes are benign.
- `002_beam-mismatch-unmodelled-sidelobes.{pdf,txt}` — tapered-aperture sky
  through a Gaussian operator, against matched-gauss with identical noise.
  At b ≈ +50° the unmodelled-sidelobe map error is 0.09–0.23 K RMS —
  6–13× below the 1/f-dominated residual — see the txt's caveats before
  generalising.
- `003_galplane-drift-sidelobe-stress.{pdf,txt}` — same drift design swept
  across the galactic plane (crossing at l ≈ 43°). The unmodelled-sidelobe
  map error grows ~35–50× to 4.5–8 K RMS, and the 2 mHz high-pass flips
  from helping to badly hurting (it removes the plane crossing itself).

Plot helpers are imported from [`../DSA/dsa_vis.py`](../DSA/dsa_vis.py).

**Coordinate-frame note.** `limTOD.GDSM_sky_model` returns pygdsm's native
*galactic*-frame map unrotated, while the simulator points the beam in
*equatorial* coordinates (see the repo CHANGELOG's known-issues entry).
These notebooks use `ska_common.gdsm_equatorial_sky_model`, which rotates
G → C at native resolution, so quoted (RA, Dec) patches contain the sky
they claim to.

Key drift-scan-specific consideration (see the notebook's opening cell): the
sky signal occupies temporal frequencies ≲ 1 mHz (one beam crossing ≈ 16 min),
overlapping the 1/f knee — so the aggressive high-pass filtering used by the
DSA azimuth-scan analysis is not available, and the notebook compares
map-making with no high-pass vs. a 2 mHz cutoff.

Survey to-dos are collected in "SKA Scanning Strategy.pdf" (repo root, not yet
folded into these notebooks).
