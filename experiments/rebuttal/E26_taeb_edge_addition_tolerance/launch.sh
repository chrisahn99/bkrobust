#!/usr/bin/env bash
# Full E26 run on the workstation, detached (alternative to the Slurm array). Run from the repository root:
#   bash experiments/rebuttal/E26_taeb_edge_addition_tolerance/launch.sh [workers]
# Writes results/rebuttal/E26/{rows.jsonl,config.json,run.log,run.pid,DONE}, then analyses.
set -euo pipefail
MAX_WORKERS=10
WORKERS="${1:-10}"
if [ "$WORKERS" -gt "$MAX_WORKERS" ]; then echo "refusing more than $MAX_WORKERS workers" >&2; exit 1; fi
HERE="experiments/rebuttal/E26_taeb_edge_addition_tolerance"
OUT="results/rebuttal/E26"
PY=".venv/bin/python"
[ -f "$HERE/run.py" ] || { echo "run from the repository root" >&2; exit 1; }
[ -x "$PY" ] || { echo "no $PY" >&2; exit 1; }
mkdir -p "$OUT"
if [ -e "$OUT/DONE" ]; then echo "$OUT/DONE exists: already finished" >&2; exit 1; fi
if [ -f "$OUT/run.pid" ] && kill -0 "$(cat "$OUT/run.pid")" 2>/dev/null; then
  echo "a run is live (pid $(cat "$OUT/run.pid"))" >&2; exit 1
fi
nohup bash -c "
  $PY -u $HERE/run.py --out $OUT --workers $WORKERS; rc=\$?
  $PY -u $HERE/analyse.py --out $OUT; arc=\$?
  echo \"run_exit=\$rc analyse_exit=\$arc finished=\$(date -Iseconds)\" > $OUT/DONE
" >> "$OUT/run.log" 2>&1 &
echo $! > "$OUT/run.pid"
echo "started pid $(cat "$OUT/run.pid"), log $OUT/run.log, done-file $OUT/DONE"
