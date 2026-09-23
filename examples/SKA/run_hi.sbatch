#!/bin/bash
# Tier 2: the HI experiment for one strategy set, plus its ablation.
#
#   sbatch run_hi.sbatch <tag> <strategy> [<strategy> ...]
#   sbatch run_hi.sbatch ladder8 drift drift8 raster
#
# The tag keeps the published results npz safe -- the filename otherwise
# encodes only the band and grid. Both drivers refuse to run untagged on a
# non-default strategy set.
#
# The ablation matters more than the experiment here: the drift residual is
# 96% beam+prior floor, so "did the ladder help?" is really "did the floor
# move?", and only the ablation separates floor from total.
#SBATCH --job-name=ska-hi
#SBATCH --partition=Main
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=%x-%j.out
set -euo pipefail

TAG="${1:?usage: sbatch run_hi.sbatch <tag> <strategy>...}"; shift
STRATEGIES=("$@")

export OMP_NUM_THREADS=1
export FASTBOX_PATH="$HOME/ska/FastBox"
source "$HOME/ska/.venv/bin/activate"
cd "$HOME/ska/limTOD-src/examples/SKA"

echo "host $(hostname)  tag ${TAG}  strategies ${STRATEGIES[*]}"
echo "=== experiment ==="
python -u run_hi_experiment.py --tag "$TAG" --strategies "${STRATEGIES[@]}"
echo "=== ablation ==="
python -u ska_hi_ablation.py --tag "$TAG" --strategies "${STRATEGIES[@]}"
echo "ALL DONE ${TAG}"
