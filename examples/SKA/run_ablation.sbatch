#!/bin/bash
# Ablation only, for when the experiment npz already exists.
#   sbatch run_ablation.sbatch <tag> <strategy>...
#SBATCH --job-name=ska-abl
#SBATCH --partition=Main
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=%x-%j.out
set -euo pipefail

TAG="${1:?usage: sbatch run_ablation.sbatch <tag> <strategy>...}"; shift
STRATEGIES=("$@")

export OMP_NUM_THREADS=1
export FASTBOX_PATH="$HOME/ska/FastBox"
source "$HOME/ska/.venv/bin/activate"
cd "$HOME/ska/limTOD-src/examples/SKA"

echo "host $(hostname)  tag ${TAG}  strategies ${STRATEGIES[*]}"
python -u ska_hi_ablation.py --tag "$TAG" --strategies "${STRATEGIES[@]}"
echo "ALL DONE ${TAG}"
