#!/usr/bin/env bash
# Warm every cache the frequency-sweep notebook needs (experiment 005).
# Channels are independent, so they run concurrently; 2 OpenMP threads each
# keeps 6 jobs inside 12 cores. The nside-256 grid check is the long pole
# (~1 h of TOD generation on its own), so it starts first.
set -u
cd "$(dirname "$0")/../.."
PY=.venv/bin/python
LOG=examples/SKA/logs
mkdir -p "$LOG"

run () {  # freq nside
  OMP_NUM_THREADS=2 $PY -c "
import sys; sys.path.insert(0, 'examples/SKA')
import ska_freq_sweep as S
tod = S.simulate_channel($1, $2)
S.build_operator($1, $2, tod)
" > "$LOG/sweep_f$1_ns$2.log" 2>&1 &
}

run 1050 256
for f in 350 525 700 875 1050; do run $f 128; done
wait
echo "sweep caches complete"
