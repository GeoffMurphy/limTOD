#!/bin/bash
# The fourth Band 1 block: 850-900 MHz, above the ~600 MHz crossover. Does the
# drift/raster ordering stay flipped? nside 128 gives 3.6-3.4 px/FWHM, the
# coarsest in the series -- see freq_env_875.sh for the grid caveat.
#   sbatch run_hi_f875.sbatch <tag> <strategy>...
#SBATCH --job-name=ska-hi875
#SBATCH --partition=Main
#SBATCH --cpus-per-task=8
# Sized from measurement, not inherited: the 500-550 MHz block (same pair,
# same nside, same nchan) took 2:53, and at 850-900 the beam is narrower so
# the patches are SMALLER and the per-channel solve cheaper. 64G/16h came
# from run_hi_freq.sbatch, which needs it for drift12m -- not run here. The
# pinned 3-strategy run including drift12m peaked at 31.8G (job 13990883),
# so 32G is ample for drift+raster. A smaller box backfills far sooner on a
# full partition, which is what left this job with a 2-day start estimate.
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=%x-%j.out
set -euo pipefail
TAG="${1:?usage: sbatch run_hi_freq.sbatch <tag> <strategy>...}"; shift
STRATEGIES=("$@")
export OMP_NUM_THREADS=1
export FASTBOX_PATH="$HOME/ska/FastBox"
source "$HOME/ska/.venv/bin/activate"
cd "$HOME/ska/limTOD-src/examples/SKA"
source ./freq_env_875.sh
echo "host $(hostname)  tag ${TAG}  strategies ${STRATEGIES[*]}"
python -u -c "import ska_hi_experiment as X; X.check_sampling()"
echo "=== experiment ==="
python -u run_hi_experiment.py --tag "$TAG" --strategies "${STRATEGIES[@]}"
echo "=== ablation ==="
python -u ska_hi_ablation.py --tag "$TAG" --strategies "${STRATEGIES[@]}"
echo "ALL DONE ${TAG}"
