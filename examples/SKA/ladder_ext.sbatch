#!/bin/bash
# Extend the Tier 1 ladder sweep to long ladders, to pin the density curve
# where the design figure currently extrapolates it.
#SBATCH --job-name=ska-ladder-ext
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
python -u ska_hi_ladder.py --counts 5 10 12 15 20 --spacings 0.5 --out hi_ladder_long
echo "ALL DONE ladder-ext"
