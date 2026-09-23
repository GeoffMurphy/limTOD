#!/bin/bash
#SBATCH --job-name=localmodes
#SBATCH --partition=Main
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.out
set -euo pipefail
export OMP_NUM_THREADS=1
export FASTBOX_PATH="$HOME/ska/FastBox"
source "$HOME/ska/.venv/bin/activate"
cd "$HOME/ska/limTOD-src/examples/SKA"
echo "=== 350-400 MHz (nside 64), mech run ==="
python -u local_modes.py results/hi_experiment_mech_f350_400_nc32_ns64.npz drift drift8 raster
echo
echo "=== 675-725 MHz (nside 128), f700 run ==="
source ./freq_env.sh
python -u local_modes.py results/hi_experiment_f700_f675_725_nc32_ns128.npz drift raster
echo "ALL DONE localmodes"
