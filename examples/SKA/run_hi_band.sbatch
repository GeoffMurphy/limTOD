#!/bin/bash
# Experiment + ablation for one band, pinned or not.
#   sbatch run_hi_band.sbatch <bandenv|default> <pin|nopin> <tag> <strategy>...
#
# bandenv  'default' for 350-400 MHz / nside 64, else a freq_env*.sh to source
# pin      'pin' uses the band's reference sigma_p (ska_hi_prior: the 3-strip
#          drift's own patch, per channel) for EVERY strategy; 'nopin' leaves
#          each strategy on its own patch variance, as every run before
#          2026-10-02 did. A pinned run needs an unpinned one with the SAME
#          strategy set to compare against.
#
# Keep the strategy list identical between the two: common_patch and the HI
# cube's box are both built from it, so a different list means a different
# pixel set AND a different HI realisation.
#
# STAGE=both|experiment|ablation (env, default both) restarts half a run when
# the other half already wrote its npz.
#SBATCH --job-name=ska-band
#SBATCH --partition=Main
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=%x-%j.out
set -euo pipefail
BANDENV="${1:?usage: sbatch run_hi_band.sbatch <bandenv|default> <pin|nopin> <tag> <strategy>...}"
PIN="${2:?need 'pin' or 'nopin'}"
TAG="${3:?need a tag}"
shift 3
STRATEGIES=("$@")
case "$PIN" in
  pin)   PINARG="--pin-prior" ;;
  nopin) PINARG="" ;;
  *) echo "second argument must be 'pin' or 'nopin', got '$PIN'" >&2; exit 2 ;;
esac
export OMP_NUM_THREADS=1
export FASTBOX_PATH="$HOME/ska/FastBox"
source "$HOME/ska/.venv/bin/activate"
cd "$HOME/ska/limTOD-src/examples/SKA"
[ "$BANDENV" = default ] || source "./$BANDENV"
echo "host $(hostname)  band ${BANDENV}  ${PIN}  tag ${TAG}  strategies ${STRATEGIES[*]}"
python -u -c "import ska_hi_experiment as X; X.check_sampling()"
STAGE="${STAGE:-both}"
if [ "$STAGE" = both ] || [ "$STAGE" = experiment ]; then
  echo "=== experiment ==="
  python -u run_hi_experiment.py --tag "$TAG" $PINARG --strategies "${STRATEGIES[@]}"
fi
if [ "$STAGE" = both ] || [ "$STAGE" = ablation ]; then
  echo "=== ablation ==="
  python -u ska_hi_ablation.py --tag "$TAG" $PINARG --strategies "${STRATEGIES[@]}"
fi
echo "ALL DONE ${TAG}"
