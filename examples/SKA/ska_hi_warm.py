"""Warm one (strategy, channel) TOD + operator cache. One Slurm array task.

``ska_hi_experiment.prepare_all`` walks every channel of every strategy in one
serial process -- roughly 17 h, which is what makes it a cluster job rather
than a laptop job. The channels are completely independent (each writes its own
``tod_*``/``op_*`` file and reads nothing the others produce), so the whole
thing is embarrassingly parallel over a Slurm array with no coordination.

    python ska_hi_warm.py <strategy> <channel-index>
    python ska_hi_warm.py drift 0

Idempotent: both underlying steps return the cached object if it already
exists, so re-running a failed array task costs nothing for the tasks that
already finished.
"""
import os
import sys
import time

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import ska_hi_experiment as X


def main():
    if len(sys.argv) != 3:
        sys.exit(f"usage: {sys.argv[0]} <strategy> <channel-index>")
    strategy, idx = sys.argv[1], int(sys.argv[2])
    freqs = X.channel_freqs()
    if not 0 <= idx < len(freqs):
        sys.exit(f"channel index {idx} outside 0..{len(freqs) - 1}")
    f = freqs[idx]

    t0 = time.time()
    print(f"[{strategy} ch{idx:02d} {f:.3f} MHz] start", flush=True)
    tod = X.simulate_foreground(strategy, f)
    t1 = time.time()
    print(f"[{strategy} ch{idx:02d}] TOD {t1 - t0:.0f} s", flush=True)
    mm = X.build_operator(strategy, f, tod)
    print(f"[{strategy} ch{idx:02d}] operator {time.time() - t1:.0f} s, "
          f"{len(mm.pixel_indices)} px", flush=True)
    print(f"[{strategy} ch{idx:02d}] done in {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
