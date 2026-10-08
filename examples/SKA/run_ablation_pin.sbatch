#!/bin/bash
# SUPERSEDED 2026-10-02 by run_hi_band.sbatch, which takes the band as an
# argument and runs the experiment as well. Kept because job 13990883 (pin700)
# was produced by this script and the npz records that provenance.
# Ablation with the prior variance PINNED across strategies.
#   sbatch run_ablation_pin.sbatch <tag> <reference> <strategy>...
#
# Why. The solve sets S^-1 = I / var(fg over THAT strategy's own patch), so a
# ladder spanning more declination picks up more sky structure and is handed a
# looser prior. fig_priorvar shows a looser prior raises tr(WA) and WORSENS the
# floor, so a long ladder is penalised for its extent independently of its
# geometry. At 675-725 MHz the 12-strip drift spans 13.16 deg against the
# 3-strip drift's 3.98 and its patch sigma is ~25% wider, i.e. a prior ~1.6x
# looser -- this run asks how much of its worse floor that accounts for.
#
# The reference strategy is solved at scale 1.0 exactly, so its floor must
# reproduce the unpinned run. That is the control: if it moves, something else
# changed.
#
# Keep the strategy list identical to the run being compared against, because
# common_patch and the HI cube's box are both built from it -- a different list
# means a different pixel set AND a different HI realisation.
#SBATCH --job-name=ska-pin
#SBATCH --partition=Main
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=10:00:00
#SBATCH --output=%x-%j.out
set -euo pipefail
TAG="${1:?usage: sbatch run_ablation_pin.sbatch <tag> <reference> <strategy>...}"; shift
REF="${1:?need a reference strategy to pin the prior to}"; shift
STRATEGIES=("$@")
export OMP_NUM_THREADS=1
export FASTBOX_PATH="$HOME/ska/FastBox"
source "$HOME/ska/.venv/bin/activate"
cd "$HOME/ska/limTOD-src/examples/SKA"
source ./freq_env.sh          # 675-725 MHz, nside 128
echo "host $(hostname)  tag ${TAG}  pin ${REF}  strategies ${STRATEGIES[*]}"
python -u ska_hi_ablation.py --tag "$TAG" --pin-prior-to "$REF" \
       --strategies "${STRATEGIES[@]}"
echo "ALL DONE ${TAG}"
