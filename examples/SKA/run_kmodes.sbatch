#!/bin/bash
# tr(WA) by angular scale, and the sufficiency areas, for one band.
#   sbatch run_kmodes.sbatch [band]          # band defaults to f700
#
# Reads the cached operators in hi_cache/ and the band's matched ablation npz;
# runs no TOD and no simulation. Minutes, not hours -- the cost is one
# eigendecomposition per strategy (npix ~ 1600 at nside 128).
#
# Writes results/hi_kmodes_<band>_*.npz, which is all the figures need: pull
# that back and draw locally with
#   /home/geoff/limTOD/.venv/bin/python ska_hi_kmodes.py --band f700 --draw-only
#SBATCH --job-name=ska-kmodes
#SBATCH --partition=Main
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=%x-%j.out
set -euo pipefail
BAND="${1:-f700}"
export OMP_NUM_THREADS=1
export FASTBOX_PATH="$HOME/ska/FastBox"
source "$HOME/ska/.venv/bin/activate"
cd "$HOME/ska/limTOD-src/examples/SKA"
# freq_env.sh pins the 675-725 grid for the OTHER modules. The kmodes band
# config is self-contained -- it takes frequency, nside and every input path
# from BANDS[band] -- and the only thing it reads from ska_hi_experiment is the
# band-independent WHITE_VAR. So source it only where it is the right band, and
# never for f350, where it would advertise the wrong grid for no reason.
if [ "$BAND" = "f700" ]; then source ./freq_env.sh; fi
echo "host $(hostname)  band ${BAND}"
# --no-figures: there is no figures/ here, and the npz is the deliverable.
python -u ska_hi_kmodes.py --band "$BAND" --no-figures
echo "ALL DONE ${BAND}"
