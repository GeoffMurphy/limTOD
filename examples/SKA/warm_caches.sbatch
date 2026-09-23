#!/bin/bash
# Warm the experiment-006 TOD + operator caches, one array task per channel.
#
#   sbatch --array=0-31 warm_caches.sbatch drift
#   sbatch --array=0-31 warm_caches.sbatch raster
#
# Tasks are independent and idempotent, so a partial failure is re-submitted
# with --array=<the failed ids> and costs nothing for what already succeeded.
#SBATCH --job-name=ska-warm
#SBATCH --partition=Main
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --output=%x-%A_%a.out
set -euo pipefail

STRATEGY="${1:?usage: sbatch --array=0-31 warm_caches.sbatch <strategy>}"

export OMP_NUM_THREADS=1          # healpy's OpenMP would oversubscribe the 2 cpus
export FASTBOX_PATH="$HOME/ska/FastBox"
source "$HOME/ska/.venv/bin/activate"
cd "$HOME/ska/limTOD-src/examples/SKA"

echo "host $(hostname)  task ${SLURM_ARRAY_TASK_ID}  strategy ${STRATEGY}"
python -u ska_hi_warm.py "$STRATEGY" "${SLURM_ARRAY_TASK_ID}"
