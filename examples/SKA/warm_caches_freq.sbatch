#!/bin/bash
# Warm caches for the 700 MHz HI repeat. One array task per channel.
#   sbatch --array=0-31 warm_caches_freq.sbatch raster
# nside 128 makes every step several times more expensive than the 350 MHz run,
# hence the longer wall limit and more memory.
#SBATCH --job-name=ska-warm7
#SBATCH --partition=Main
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=%x-%A_%a.out
set -euo pipefail
STRATEGY="${1:?usage: sbatch --array=0-31 warm_caches_freq.sbatch <strategy>}"
export OMP_NUM_THREADS=1
export FASTBOX_PATH="$HOME/ska/FastBox"
source "$HOME/ska/.venv/bin/activate"
cd "$HOME/ska/limTOD-src/examples/SKA"
source ./freq_env.sh
echo "host $(hostname)  task ${SLURM_ARRAY_TASK_ID}  strategy ${STRATEGY}"
python -u -c "import ska_hi_experiment as X; X.check_sampling()"
python -u ska_hi_warm.py "$STRATEGY" "${SLURM_ARRAY_TASK_ID}"
