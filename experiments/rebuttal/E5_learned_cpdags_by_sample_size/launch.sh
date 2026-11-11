#!/usr/bin/env bash
# E5 full run on the workstation: nohup, the repo's venv, a log, a done-file. Run from anywhere.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
cd "$REPO"

WORKERS="${WORKERS:-6}"            # spec: about 1 h on 6 workers
OUT="${OUT:-results/rebuttal/E5}"
if [ "$WORKERS" -gt 8 ]; then echo "refusing: at most 8 workers (15 GB machine)" >&2; exit 1; fi

mkdir -p "$OUT"
LOG="$OUT/run.log"; PID="$OUT/run.pid"; DONE="$OUT/run.done"
if [ -e "$DONE" ]; then echo "refusing: $DONE exists (finished run); move it away to rerun" >&2; exit 1; fi
if [ -e "$PID" ] && kill -0 "$(cat "$PID")" 2>/dev/null; then
  echo "refusing: run already live, pid $(cat "$PID")" >&2; exit 1
fi
if [ ! -x .venv/bin/python ]; then echo "refusing: .venv/bin/python missing in $REPO" >&2; exit 1; fi
if [ ! -f experiments/rebuttal/_data/pest/graphs.jsonl ]; then
  echo "refusing: data missing under experiments/rebuttal/_data/pest/" >&2; exit 1
fi

nohup bash -c "
  .venv/bin/python '$HERE/run.py' --workers $WORKERS --out '$OUT'; rc=\$?
  .venv/bin/python '$HERE/analyse.py' --out '$OUT' > '$OUT/analyse.txt' 2>&1; ra=\$?
  echo \"run_exit=\$rc analyse_exit=\$ra \$(date -Is)\" > '$DONE'
  rm -f '$PID'
" > "$LOG" 2>&1 &
echo $! > "$PID"
echo "E5 started, pid $(cat "$PID"), log $LOG, done-file $DONE"
