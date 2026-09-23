#!/bin/bash
# The 700 MHz HI repeat: does the Band 1 negative result hold above its bottom
# edge? Same pipeline, same geometry, nside 128 to keep 4.2-4.5 px/FWHM.
#   sbatch run_hi_freq.sbatch <tag> <strategy>...
#SBATCH --job-name=ska-hi7
#SBATCH --partition=Main
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=16:00:00
#SBATCH --output=%x-%j.out
set -euo pipefail
TAG="${1:?usage: sbatch run_hi_freq.sbatch <tag> <strategy>...}"; shift
STRATEGIES=("$@")
export OMP_NUM_THREADS=1
export FASTBOX_PATH="$HOME/ska/FastBox"
source "$HOME/ska/.venv/bin/activate"
cd "$HOME/ska/limTOD-src/examples/SKA"
source ./freq_env.sh
echo "host $(hostname)  tag ${TAG}  strategies ${STRATEGIES[*]}"
python -u -c "import ska_hi_experiment as X; X.check_sampling()"
echo "=== experiment ==="
python -u run_hi_experiment.py --tag "$TAG" --strategies "${STRATEGIES[@]}"
echo "=== ablation ==="
python -u ska_hi_ablation.py --tag "$TAG" --strategies "${STRATEGIES[@]}"
echo "ALL DONE ${TAG}"
