#!/bin/bash
# Reproduce the published drift mode count on ilifu before trusting anything
# else here. Laptop: tr(WA) = 38.77 over 265.2 deg^2.
#SBATCH --job-name=ska-validate
#SBATCH --partition=Main
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:45:00
#SBATCH --output=%x-%j.out
set -euo pipefail
export OMP_NUM_THREADS=1
export FASTBOX_PATH="$HOME/ska/FastBox"
source "$HOME/ska/.venv/bin/activate"
cd "$HOME/ska/limTOD-src/examples/SKA"
echo "host $(hostname)"
python -u validate_ilifu.py
