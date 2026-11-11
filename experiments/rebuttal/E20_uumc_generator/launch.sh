#!/usr/bin/env bash
# Full E20 run on the workstation, detached. Run from the repository root:
#   bash experiments/rebuttal/E20_uumc_generator/launch.sh [workers]
# Writes results/rebuttal/E20/{units.shard0of1.jsonl,config.shard0of1.json,results.json,run.log,run.pid,DONE}.
set -euo pipefail
WORKERS="${1:-6}"
# The cost is per query (child CO2->DuctFlow dominates), so more than 8 workers buys nothing.
if [ "$WORKERS" -gt 8 ]; then echo "refusing more than 8 workers" >&2; exit 1; fi
HERE="experiments/rebuttal/E20_uumc_generator"
OUT="results/rebuttal/E20"
PY=".venv/bin/python"
[ -f "$HERE/run.py" ] || { echo "run from the repository root" >&2; exit 1; }
[ -x "$PY" ] || { echo "no $PY" >&2; exit 1; }
for f in results/elicit/knowledge.json results/mean_table_fullspace/instances.jsonl \
         results/mean_table_relevant_v2/x_component/instances.jsonl \
         results/mean_table_relevant_v2/extended/instances.jsonl; do
  [ -f "$f" ] || { echo "missing input $f" >&2; exit 1; }
done
mkdir -p "$OUT"
if [ -e "$OUT/DONE" ]; then echo "$OUT/DONE exists: already finished" >&2; exit 1; fi
if [ -f "$OUT/run.pid" ] && kill -0 "$(cat "$OUT/run.pid")" 2>/dev/null; then
  echo "a run is live (pid $(cat "$OUT/run.pid"))" >&2; exit 1
fi
export PYTHONPATH="src${PYTHONPATH:+:$PYTHONPATH}"
nohup bash -c "
  $PY -u $HERE/run.py --out $OUT --workers $WORKERS; rc=\$?
  $PY -u $HERE/analyse.py --out $OUT; arc=\$?
  echo \"run_exit=\$rc analyse_exit=\$arc finished=\$(date -Iseconds)\" > $OUT/DONE
" >> "$OUT/run.log" 2>&1 &
echo $! > "$OUT/run.pid"
echo "started pid $(cat "$OUT/run.pid"), log $OUT/run.log, done-file $OUT/DONE"
