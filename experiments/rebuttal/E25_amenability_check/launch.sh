#!/usr/bin/env bash
# Full E25 run on the workstation, detached. Run from the repository root:
#   bash experiments/rebuttal/E25_amenability_check/launch.sh [workers]
# Writes results/rebuttal/E25/{rows.jsonl,config.json,run.log,run.pid,results.json,DONE}.
set -euo pipefail
MAX_WORKERS=8
WORKERS="${1:-6}"
if [ "$WORKERS" -gt "$MAX_WORKERS" ]; then echo "refusing more than $MAX_WORKERS workers" >&2; exit 1; fi
HERE="experiments/rebuttal/E25_amenability_check"
OUT="results/rebuttal/E25"
PY=".venv/bin/python"
[ -f "$HERE/run.py" ] || { echo "run from the repository root" >&2; exit 1; }
[ -x "$PY" ] || { echo "no $PY" >&2; exit 1; }
for f in experiments/rebuttal/_data/per_query_records.csv experiments/rebuttal/_data/knowledge_states.json \
         experiments/rebuttal/_data/chris_e10/amenability_radius.csv; do
  [ -f "$f" ] || { echo "missing data file $f" >&2; exit 1; }
done
mkdir -p "$OUT"
if [ -e "$OUT/DONE" ]; then echo "$OUT/DONE exists: already finished" >&2; exit 1; fi
if [ -f "$OUT/run.pid" ] && kill -0 "$(cat "$OUT/run.pid")" 2>/dev/null; then
  echo "a run is live (pid $(cat "$OUT/run.pid"))" >&2; exit 1
fi
nohup bash -c "
  $PY -u $HERE/run.py --out $OUT --workers $WORKERS; rc=\$?
  $PY -u $HERE/analyse.py --out $OUT > $OUT/analyse.log 2>&1; arc=\$?
  echo \"run_exit=\$rc analyse_exit=\$arc finished=\$(date -Iseconds)\" > $OUT/DONE
" >> "$OUT/run.log" 2>&1 &
echo $! > "$OUT/run.pid"
echo "started pid $(cat "$OUT/run.pid"), log $OUT/run.log, done-file $OUT/DONE"
