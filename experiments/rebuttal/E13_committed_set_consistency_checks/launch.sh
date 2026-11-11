#!/usr/bin/env bash
# E13 full run on the workstation: nohup, from the repository root, with the project venv.
# Refuses to start if the run already finished or is still alive.
set -euo pipefail
cd "$(dirname "$0")/../../.."
OUT=results/rebuttal/E13
WORKERS="${WORKERS:-6}"            # 16 cores, 15 GB: a 4,096-extension G0 can take a few hundred MB per worker
if [ "$WORKERS" -gt 8 ]; then echo "at most 8 workers" >&2; exit 1; fi
mkdir -p "$OUT"
DONE="$OUT/E13.done"; PID="$OUT/E13.pid"; LOG="$OUT/E13.log"
if [ -e "$DONE" ]; then echo "already done: $DONE" >&2; exit 1; fi
if [ -e "$PID" ] && kill -0 "$(cat "$PID")" 2>/dev/null; then echo "still running: pid $(cat "$PID")" >&2; exit 1; fi
test -x .venv/bin/python || { echo "no .venv/bin/python" >&2; exit 1; }
test -f experiments/rebuttal/_data/per_query_records.csv || { echo "missing experiments/rebuttal/_data/per_query_records.csv" >&2; exit 1; }
nohup bash -c "
  PYTHONPATH=src .venv/bin/python experiments/rebuttal/E13_committed_set_consistency_checks/run.py \
    --workers $WORKERS --out $OUT; rc=\$?
  PYTHONPATH=src .venv/bin/python experiments/rebuttal/E13_committed_set_consistency_checks/analyse.py --out $OUT || true
  echo \"exit \$rc \$(date -Is)\" > $DONE
" > "$LOG" 2>&1 &
echo $! > "$PID"
echo "started pid $(cat "$PID"), log $LOG, done-file $DONE"
