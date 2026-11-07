#!/usr/bin/env bash
# Full E17 run on the workstation, detached. Run from the repository root:
#   bash experiments/rebuttal/E17_worked_case_and_validity_by_radius/launch.sh [workers]
# Writes results/rebuttal/E17/{rows.*.jsonl,config.*.json,run.log,run.pid,DONE}, then analyses.
set -euo pipefail
WORKERS="${1:-8}"
if [ "$WORKERS" -gt 8 ]; then echo "refusing more than 8 workers" >&2; exit 1; fi
HERE="experiments/rebuttal/E17_worked_case_and_validity_by_radius"
OUT="results/rebuttal/E17"
PY=".venv/bin/python"
DATA="experiments/rebuttal/_data/per_query_records.csv"
[ -f "$HERE/run.py" ] || { echo "run from the repository root" >&2; exit 1; }
[ -x "$PY" ] || { echo "no $PY" >&2; exit 1; }
[ -f "$DATA" ] || { echo "missing data file $DATA" >&2; exit 1; }
mkdir -p "$OUT"
if [ -e "$OUT/DONE" ]; then echo "$OUT/DONE exists: already finished" >&2; exit 1; fi
if [ -f "$OUT/run.pid" ] && kill -0 "$(cat "$OUT/run.pid")" 2>/dev/null; then
  echo "a run is live (pid $(cat "$OUT/run.pid"))" >&2; exit 1
fi
nohup bash -c "
  $PY -u $HERE/run.py --out $OUT --workers $WORKERS; rc=\$?
  $PY -u $HERE/analyse.py --out $OUT > $OUT/analyse.txt 2>&1; arc=\$?
  echo \"run_exit=\$rc analyse_exit=\$arc finished=\$(date -Iseconds)\" > $OUT/DONE
" >> "$OUT/run.log" 2>&1 &
echo $! > "$OUT/run.pid"
echo "started pid $(cat "$OUT/run.pid"), log $OUT/run.log, done-file $OUT/DONE"
