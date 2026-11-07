#!/usr/bin/env bash
# Full E10 run on the workstation, detached. Run from the repository root:
#   bash experiments/rebuttal/E10_radius_across_valid_sets/launch.sh [workers]
# Writes results/rebuttal/E10/{rows_networks.jsonl,rows_census.jsonl,config.json,run.log,run.pid,DONE},
# then analyses. The specification says one process, so the default and the cap are 1.
set -euo pipefail
MAX_WORKERS=1
WORKERS="${1:-1}"
if [ "$WORKERS" -gt "$MAX_WORKERS" ]; then echo "refusing more than $MAX_WORKERS worker(s)" >&2; exit 1; fi
HERE="experiments/rebuttal/E10_radius_across_valid_sets"
OUT="results/rebuttal/E10"
PY=".venv/bin/python"
[ -f "$HERE/run.py" ] || { echo "run from the repository root" >&2; exit 1; }
[ -x "$PY" ] || { echo "no $PY" >&2; exit 1; }
for f in experiments/rebuttal/_data/per_query_records.csv experiments/rebuttal/_data/knowledge_states.json; do
  [ -f "$f" ] || { echo "missing data file $f" >&2; exit 1; }
done
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
