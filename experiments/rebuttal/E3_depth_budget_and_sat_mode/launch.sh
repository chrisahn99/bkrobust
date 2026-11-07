#!/usr/bin/env bash
# Full E3 run on the workstation: serial (one worker), as registered, so timings are uncontended.
# Run from anywhere; it changes to the repository root.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
cd "$ROOT"

OUT="results/rebuttal/E3"
WORKERS=1   # PREREG: serial. Never above 8 on this machine.
LOGDIR="$OUT/full"
LOG="$LOGDIR/run.log"
DONE="$LOGDIR/DONE"
PID="$LOGDIR/run.pid"
mkdir -p "$LOGDIR"

if [ -e "$DONE" ]; then
  echo "refusing: $DONE exists (the run finished; analyse it, do not rerun)"; exit 1
fi
if [ -e "$PID" ] && kill -0 "$(cat "$PID")" 2>/dev/null; then
  echo "refusing: run already live with pid $(cat "$PID")"; exit 1
fi
if [ ! -x .venv/bin/python ]; then
  echo "refusing: .venv/bin/python not found in $ROOT"; exit 1
fi

nohup bash -c "
  .venv/bin/python -u '$HERE/run.py' --out '$OUT' --workers $WORKERS
  rc=\$?
  echo \"exit=\$rc \$(date -Is)\" > '$DONE'
  rm -f '$PID'
" >> "$LOG" 2>&1 &
echo $! > "$PID"
echo "started pid $(cat "$PID"); log $ROOT/$LOG; done-file $ROOT/$DONE"
