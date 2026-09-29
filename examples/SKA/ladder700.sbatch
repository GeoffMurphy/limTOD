#!/bin/bash
# Can a PARKED drift reach the raster's mode density at 675-725 MHz?
#
# At 350 MHz it can: drift12 gets to 21.81 against the raster's 22.31, and the
# matched-density comparison then leaves cross-linking worth 9.3x. At 700 MHz
# the raster is at 76.55 per 100 deg^2 and the published drift at 41.87, a
# bigger gap -- and the 350 MHz sweep showed the drift's density SATURATES
# around 24 (a factor 1.63 over 3 strips), so matching may be impossible.
# Settling that is the point: if the drift cannot reach it, "matched density"
# is not a test that exists at this band, which is itself the answer.
#
# Tier 1 -- operator only, one channel, no TOD. Spacing 0.25 beams AT 700 MHz
# (the tightest sensible, where strips overlap and densify cross-scan);
# ska_hi_ladder quotes spacings in beams at --freq.
#SBATCH --job-name=ska-ladder700
#SBATCH --partition=Main
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=16:00:00
#SBATCH --output=%x-%j.out
set -euo pipefail
export OMP_NUM_THREADS=1
export FASTBOX_PATH="$HOME/ska/FastBox"
source "$HOME/ska/.venv/bin/activate"
cd "$HOME/ska/limTOD-src/examples/SKA"
source ./freq_env.sh
echo "host $(hostname)"
# Counts and output name are arguments so a follow-up sweep does not need a
# second script. --ref-count is pinned to the FIRST count so the spacing sweep
# lands on a point the count sweep already covers, instead of silently adding
# the default N = 3 again.
COUNTS="${@:-3 8 12 20}"
OUT="${OUT:-hi_ladder_f700}"
FIRST="${COUNTS%% *}"
python -u ska_hi_ladder.py --freq 675.78125 --counts $COUNTS \
    --spacings 0.25 --ref-spacing 0.25 --ref-count "$FIRST" --out "$OUT"
echo "ALL DONE ladder700"
