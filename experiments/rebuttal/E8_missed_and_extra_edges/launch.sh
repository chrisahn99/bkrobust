#!/usr/bin/env bash
# E8 full run on the workstation (16 cores, 15 GB): nohup, log, pid file, done-file.
# Run from the repository root:  bash experiments/rebuttal/E8_missed_and_extra_edges/launch.sh [workers]
set -euo pipefail
WORKERS="${1:-8}"
if [ "$WORKERS" -gt 8 ]; then echo "refusing more than 8 workers" >&2; exit 1; fi
[ -f pyproject.toml ] && [ -d experiments/rebuttal/E8_missed_and_extra_edges ] || {
  echo "run this from the repository root" >&2; exit 1; }
[ -x .venv/bin/python ] || { echo "no .venv/bin/python" >&2; exit 1; }

OUT=results/rebuttal/E8
mkdir -p "$OUT"
DONE="$OUT/E8.done"
PID="$OUT/E8.pid"
LOG="$OUT/E8_$(date +%Y%m%d_%H%M%S).log"
if [ -e "$DONE" ]; then echo "already finished: $DONE exists" >&2; exit 1; fi
if [ -e "$PID" ] && kill -0 "$(cat "$PID")" 2>/dev/null; then
  echo "already running: pid $(cat "$PID")" >&2; exit 1
fi

# Done-file only on exit 0; any other exit writes E8.failed instead, so relaunching resumes the run
# (run.py skips units already written) and a failure is never mistaken for a result.
rm -f "$OUT/E8.failed"
nohup bash -c "
  .venv/bin/python experiments/rebuttal/E8_missed_and_extra_edges/run.py --out $OUT --workers $WORKERS
  rc=\$?
  if [ \$rc -eq 0 ]; then echo \"exit=0 \$(date -Iseconds)\" > $DONE;
  else echo \"exit=\$rc \$(date -Iseconds)\" > $OUT/E8.failed; fi
  rm -f $PID
" > "$LOG" 2>&1 &
echo $! > "$PID"
echo "started pid $(cat "$PID"), log $LOG, done-file $DONE"
