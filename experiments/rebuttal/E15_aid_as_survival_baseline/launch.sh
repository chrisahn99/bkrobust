#!/usr/bin/env bash
# Full E15 run on the workstation, detached. Run from anywhere; it cd's to the repository root.
set -euo pipefail
cd "$(dirname "$0")/../../.."
WORKERS="${WORKERS:-8}"
if [ "$WORKERS" -gt 8 ]; then echo "at most 8 workers on the workstation" >&2; exit 1; fi
OUT=results/rebuttal/E15
mkdir -p "$OUT"
DONE="$OUT/E15.done"; PID="$OUT/E15.pid"; LOG="$OUT/E15.log"
if [ -e "$DONE" ]; then echo "already finished: $DONE exists" >&2; exit 1; fi
if [ -e "$PID" ] && kill -0 "$(cat "$PID")" 2>/dev/null; then
  echo "already running: pid $(cat "$PID")" >&2; exit 1
fi
nohup bash -c "
  PYTHONPATH=src .venv/bin/python experiments/rebuttal/E15_aid_as_survival_baseline/run.py --out $OUT --workers $WORKERS
  rc=\$?
  # 0 = all rows consistent, 3 = all rows written but some flagged: both are a finished run.
  # Anything else is a crash: no done-file, so relaunching resumes from the rows written.
  if [ \$rc -eq 0 ] || [ \$rc -eq 3 ]; then
    PYTHONPATH=src .venv/bin/python experiments/rebuttal/E15_aid_as_survival_baseline/analyse.py --out $OUT || rc=5
    echo \"exit=\$rc \$(date -u +%Y-%m-%dT%H:%M:%SZ)\" > $DONE
  else
    echo \"run.py exit=\$rc \$(date -u +%Y-%m-%dT%H:%M:%SZ)\" > $OUT/E15.crashed
  fi
  rm -f $PID
" > "$LOG" 2>&1 &
echo $! > "$PID"
echo "started pid $(cat "$PID"), log $LOG, done-file $DONE"
