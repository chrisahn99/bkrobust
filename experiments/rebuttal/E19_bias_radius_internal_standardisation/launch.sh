#!/usr/bin/env bash
# Full E19 run on the workstation (16 cores, 15 GB): nohup, log, pid file, done-file.
# Usage, from the repository root:  bash experiments/rebuttal/E19_bias_radius_internal_standardisation/launch.sh [WORKERS]
set -euo pipefail
cd "$(dirname "$0")/../../.."
WORKERS="${1:-6}"
if [ "$WORKERS" -gt 8 ]; then echo "refusing: at most 8 workers" >&2; exit 2; fi
OUT=results/rebuttal/E19
LOG="$OUT/run.log"; PID="$OUT/run.pid"; DONE="$OUT/DONE"
mkdir -p "$OUT"
if [ -e "$DONE" ]; then echo "refusing: $DONE exists (the run finished; remove it to rerun)" >&2; exit 1; fi
if [ -e "$PID" ] && kill -0 "$(cat "$PID")" 2>/dev/null; then
  echo "refusing: run already live, pid $(cat "$PID")" >&2; exit 1
fi
[ -x .venv/bin/python ] || { echo "no .venv/bin/python in $(pwd)" >&2; exit 1; }
E=experiments/rebuttal/E19_bias_radius_internal_standardisation
# each worker is one process with ~200 MB; 8 x 200 MB stays far below 15 GB
nohup bash -c "
  .venv/bin/python $E/run.py --out $OUT --workers $WORKERS; rc=\$?
  .venv/bin/python $E/analyse.py --out $OUT > $OUT/analyse.txt 2>&1 || true
  echo \"exit=\$rc \$(date -Iseconds)\" > $DONE
  rm -f $PID
" >> "$LOG" 2>&1 &
echo $! > "$PID"
echo "started pid $(cat "$PID"), workers $WORKERS, log $LOG, done-file $DONE"
