# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Known issues

- `limTOD.sky_model.GDSM_sky_model` returns pygdsm's native **galactic**-frame
  map without rotating to equatorial, while `TODSim` and `HPW_mapmaking` point
  the beam in equatorial coordinates. Any analysis using it therefore observes
  the GDSM with its galactic coordinates silently reinterpreted as (RA, Dec) —
  self-consistent between forward model and map-maker, but the sky content at
  a quoted (RA, Dec) is wrong (e.g. the galactic plane lies along Dec ≈ 0
  instead of its true celestial track). This affects the patch descriptions in
  `examples/DSA/`. The SKA notebooks now use
  `examples/SKA/ska_common.gdsm_equatorial_sky_model`, which rotates G → C at
  native resolution before degrading; consider applying the same fix (or an
  explicit `coord` argument) to `GDSM_sky_model` itself in a future release.

### Added

- `examples/SKA/` — drift-scan strategy notebook (`ska_drift_scan.ipynb`) for an
  SKA-Mid-like 15 m dish at the Karoo site: telescope parked at fixed (Az, El),
  elevation stepped per sidereal night, map-making with and without the
  high-pass filter. Shared helpers in `examples/SKA/ska_common.py`
  (chromatic 1.22 λ/D Gaussian beam, drift-night timing, AltAz→ICRS track).
- Tapered-aperture ("realistic Airy") beam model in `examples/SKA/ska_common.py`
  (closed-form circular-aperture diffraction with configurable edge taper;
  first sidelobe ≈ −23 dB at the default −14 dB taper). The drift-scan
  notebook now runs the Gaussian and tapered-aperture beams through the
  identical pipeline to isolate the cost of sidelobes in drift-scan mode.
- Beam-mismatch scenario in the drift-scan notebook (tapered-aperture sky
  solved with the Gaussian operator — the cost of *unmodelled* sidelobes),
  per-night noise seeding shared across beams so scenario comparisons use
  identical noise realizations, and a numbered experiment-record PDF archive
  (`examples/SKA/results/`, written by `ska_common.save_results_pdf`) with
  companion `.txt` notes holding findings and caveats.
- `examples/SKA/ska_drift_galplane.ipynb` — galactic-plane sidelobe stress
  test (experiment 003), plus `ska_common.gdsm_equatorial_sky_model`, the
  galactic→celestial-rotated GDSM wrapper that the frame issue above made
  necessary (all SKA experiments re-run with it).
- Power-spectrum comparison in `ska_drift_scan.ipynb`: 1D pseudo-Cℓ per
  scenario (transfer function T_ℓ and truth-correlation r_ℓ, recipe from
  `examples/DSA/scripts/compare_power_spectra.py`) and a 2D flat-sky power
  spectrum of the map residual, which resolves the drift-scan 1/f signature
  Cℓ integrates away — a band at low ℓ_Dec extended in ℓ_RA (the drift is
  coherent across the sidereal-offset Dec rows), suppressed by the 2 mHz
  high-pass.

## [1.2.0] - 2025-10-06

### Added

- Add `CHANGELOG.md`

### Changed

- **BREAKING**: Renamed `TODsim` class to `TODSim` for better Python naming conventions
- Renamed old `limTODsim` class references to `TODSim` throughout codebase
- Updated all import statements and class instantiations to use `TODSim`
- Updated `__init__.py` and `__all__` exports to reflect new class name

### Improved

- 📝 **Documentation**:
  - **Structure**: Moved "Latest Updates" section from README.md to dedicated CHANGELOG.md file
  - 📋 **Table of Contents**: Updated README.md Table of Contents to accurately reflect document structure
  - **Examples**: Removed in flavour of example notebooks to simplify maintainace
- 🔧 **Code Organization**:
  - Move example notebooks to `examples/`
  - Improved consistency in class naming across all files including:
    - Source code (`simulator.py`)
    - Package exports (`__init__.py`)
    - Documentation (`README.md`)
    - Example notebooks (`examples.ipynb`, `mm_example.ipynb`)
    - Change Log (`CHANGELOG.md`)

### Fixed

- Corrected all references to use consistent `TODSim` class name
- Fixed import statements in example notebooks and documentation

## [1.1.0] - 2025-10-05

### Added

- 🎯 **Full Stokes Support**: Added complete polarization handling (I, Q, U, V) for both TOD simulation and map-making
- 🗺️ **Map-Making Pipeline**: Implemented `HPW_mapmaking` class combining high-pass filtering and Wiener filtering for sky reconstruction from TOD
- 🎲 **Gaussian Random Field Generator**: Added generator for correlated sky realizations from frequency-frequency angular power spectra C_ℓ(ν,ν'), enabling realistic simulation of line intensity mapping signals with spectral correlations (credit: Katrine Alice Glasscock, Philip Bull)
- 📓 **Example Notebooks**: Added comprehensive Jupyter notebook demonstrating the full map-making workflow ([examples/mm_example.ipynb](examples/mm_example.ipynb))

### Changed

- **BREAKING**: `beam_func` and `sky_func` now require keyword-only arguments, two of which must be `freq` and `nside`:

  ```python
  # Old: beam_func(freq, nside) and sky_func(freq, nside)
  # New: beam_func(freq=xx, nside=xx) and sky_func(freq=xx, nside=xx)
  ```

- **BREAKING**: Function outputs must be HEALPix maps with specific shapes:
  - 1D array of length npix for unpolarized (**I**) beam/sky
  - 2D array of shape (3, npix) for polarized (**I, Q, U**) beam/sky
  - 2D array of shape (4, npix) for polarized (**I, Q, U, V**) beam/sky

### Fixed

- 🐛 **Bug Fix**: Corrected a critical sign error in coordinate rotation transformations

## [1.0.0] - 2025-09-01

### Initial Release

- Initial release of limTOD: Time-Ordered Data simulation for single-dish radio telescopes

### Key Features

- TOD simulation with realistic noise models (1/f noise, white noise, gain variations)
- Support for asymmetric beam patterns
- Direct beam convolution using HEALPix spherical harmonics rotation and sum to calculate Tsky
- Flexible beam and sky model functions
- Global Sky Model (GDSM) integration
- MPI parallelization support
- Example scanning patterns and beam models
- Example notebooks for getting started
- Documentation and examples
