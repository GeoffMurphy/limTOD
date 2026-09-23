# Band/grid for the 700 MHz HI repeat. Source before any ska_hi_* script.
# nside 128 keeps 4.2-4.5 px/FWHM, matching the 350-400 MHz run's 3.8-4.4;
# nside 64 would give 2.1-2.3 and undersample the beam (HANDOFF grid check).
export SKA_HI_F_LO=675 SKA_HI_F_HI=725 SKA_HI_NCHAN=32 SKA_HI_NSIDE=128 SKA_HI_NS_CELL=32
