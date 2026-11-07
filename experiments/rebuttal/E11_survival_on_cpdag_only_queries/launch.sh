#!/usr/bin/env bash
# Full E11 run on the workstation, detached. Run from the repository root:
#   bash experiments/rebuttal/E11_survival_on_cpdag_only_queries/launch.sh [workers]
# Writes results/rebuttal/E11/{units.*.jsonl,config.*.json,run.log,run.pid,DONE}, then analyses.
set -euo pipefail
WORKERS="${1:-8}"
if [ "$WORKERS" -gt 8 ]; then echo "refusing more than 8 workers" >&2; exit 1; fi
HERE="experiments/rebuttal/E11_survival_on_cpdag_only_queries"
OUT="results/rebuttal/E11"
PY=".venv/bin/python"
[ -f "$HERE/run.py" ] || { echo "run from the repository root" >&2; exit 1; }
[ -x "$PY" ] || { echo "no $PY" >&2; exit 1; }
for f in results/frame/frame.jsonl results/axis_robustness_real/frame.jsonl results/elicit/knowledge.json \
         results/axisa3/networks/acquisition_manifest.json; do
  [ -f "$f" ] || { echo "missing input $f" >&2; exit 1; }
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
