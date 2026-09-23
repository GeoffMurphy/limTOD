#!/bin/bash
# One HI realisation, drift + raster, for the error bar on residual/HI.
# Only the TRUE HI seed varies; noise seeds are fixed, so the spread measures
# HI realisation scatter alone.
#   sbatch run_hi_seed.sbatch <seed>
#SBATCH --job-name=ska-seed
#SBATCH --partition=Main
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=%x-%j.out
set -euo pipefail
SEED="${1:?usage: sbatch run_hi_seed.sbatch <seed>}"
export OMP_NUM_THREADS=1
export FASTBOX_PATH="$HOME/ska/FastBox"
source "$HOME/ska/.venv/bin/activate"
cd "$HOME/ska/limTOD-src/examples/SKA"
echo "host $(hostname)  HI seed ${SEED}"
python -u run_hi_experiment.py --tag "seed${SEED}" --hi-seed "${SEED}" \
       --strategies drift raster --mocks 10
echo "ALL DONE seed${SEED}"
