#!/usr/bin/env bash
# Full E2 run on the workstation, serial (the specification says serial; never more than 8 workers).
# Run from anywhere; it cd's to the repository root. Refuses to start if a done-file or a live
# pid file exists. Resumable: rerunning after a crash skips units already written.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
cd "$ROOT"

OUT="${OUT:-results/rebuttal/E2}"
WORKERS="${WORKERS:-1}"
PY="$ROOT/.venv/bin/python"
LOG="$OUT/run.log"
PID="$OUT/run.pid"
DONE="$OUT/DONE"

if [ "$WORKERS" -gt 8 ]; then echo "refusing: WORKERS=$WORKERS > 8" >&2; exit 2; fi
if [ ! -x "$PY" ]; then echo "refusing: no virtual environment at $PY" >&2; exit 2; fi
mkdir -p "$OUT"
if [ -e "$DONE" ]; then echo "refusing: $DONE exists (run finished; move it to rerun)" >&2; exit 2; fi
if [ -e "$PID" ] && kill -0 "$(cat "$PID")" 2>/dev/null; then
  echo "refusing: run already live, pid $(cat "$PID")" >&2; exit 2
fi

nohup bash -c "
  '$PY' '$HERE/run.py' --out '$OUT' --workers $WORKERS
  rc=\$?
  '$PY' '$HERE/analyse.py' --out '$OUT' > '$OUT/analyse.txt' 2>&1 || true
  echo \"exit \$rc \$(date -Is)\" > '$DONE'
  rm -f '$PID'
" > "$LOG" 2>&1 &
echo $! > "$PID"
echo "E2 started: pid $(cat "$PID"), log $LOG, done-file $DONE"
