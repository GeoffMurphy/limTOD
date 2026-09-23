#!/bin/bash
# Mode density at 700 MHz for both strategies -- the diagnostic for why the
# raster/drift ordering inverts between 350 and 700 MHz.
#SBATCH --job-name=modes700
#SBATCH --partition=Main
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --output=%x-%j.out
set -euo pipefail
export OMP_NUM_THREADS=1
export FASTBOX_PATH="$HOME/ska/FastBox"
source "$HOME/ska/.venv/bin/activate"
cd "$HOME/ska/limTOD-src/examples/SKA"
source ./freq_env.sh
python -u modes700.py
