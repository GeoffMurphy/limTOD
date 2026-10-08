# 850-900 MHz block: the fourth point in Band 1, above the ~600 MHz crossover.
#
# nside 128 gives 3.6-3.4 px/FWHM, which is the COARSEST block in the series
# (350-400 runs 4.4 at nside 64; 500-550 runs 6.1-5.5; 675-725 runs 4.5-4.2)
# and sits below the ~5 px/FWHM that HANDOFF's grid check recommends when the
# RESIDUAL, rather than the mode count, is the quantity of interest. At 875 MHz
# that check measured the drift's residual falling 19% from nside 128 to 256,
# so these numbers are expected to read HIGH by of order that much.
#
# nside 256 is not an option for a full block: the single-channel re-run cost
# ~1.2 h, so 32 channels x 2 strategies is ~77 h. The bias is far too small to
# threaten the ordering (the drift leads the raster 2.2x at 675-725), which is
# what this block is for -- but do not quote its residual/HI at the same
# precision as the other three without a raster grid check at 875.
export SKA_HI_F_LO=850 SKA_HI_F_HI=900 SKA_HI_NCHAN=32 SKA_HI_NSIDE=128 SKA_HI_NS_CELL=32
