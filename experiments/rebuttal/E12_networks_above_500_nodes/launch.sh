#!/usr/bin/env bash
# Full E12 run on the workstation, detached. Run from the repository root:
#   bash experiments/rebuttal/E12_networks_above_500_nodes/launch.sh [workers]
# Writes results/rebuttal/E12/{rows.jsonl,config.json,run.log,run.pid,DONE}, then analyses.
# Six workers by default: each has a 2 GB data cap, so six stay under the 15 GB of the machine.
set -euo pipefail
WORKERS="${1:-6}"
if [ "$WORKERS" -gt 8 ]; then echo "refusing more than 8 workers" >&2; exit 1; fi
HERE="experiments/rebuttal/E12_networks_above_500_nodes"
OUT="results/rebuttal/E12"
PY=".venv/bin/python"
[ -f "$HERE/run.py" ] || { echo "run from the repository root" >&2; exit 1; }
[ -x "$PY" ] || { echo "no $PY" >&2; exit 1; }
for f in experiments/rebuttal/_data/knowledge_states.json results/frame/frame.jsonl; do
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
