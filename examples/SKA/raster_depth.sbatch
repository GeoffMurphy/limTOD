#!/bin/bash
# PLAN item 5: find the azimuth throw at which the raster natively matches the
# drift's pixel count, so the residual comparison stops confounding geometry
# with depth. Tier 1 -- operator only, one channel, no TOD.
#SBATCH --job-name=ska-rdepth
#SBATCH --partition=Main
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=%x-%j.out
set -euo pipefail
export OMP_NUM_THREADS=1
export FASTBOX_PATH="$HOME/ska/FastBox"
source "$HOME/ska/.venv/bin/activate"
cd "$HOME/ska/limTOD-src/examples/SKA"
echo "host $(hostname)"
python -u ska_hi_raster_depth.py --halfwidths 7.5 5 4 3 2.5 2
echo "ALL DONE raster-depth"
