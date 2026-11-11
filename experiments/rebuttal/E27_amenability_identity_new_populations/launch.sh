#!/usr/bin/env bash
# Full E27 run on the workstation, detached. Run from the repository root:
#   bash experiments/rebuttal/E27_amenability_identity_new_populations/launch.sh [workers]
# Writes results/rebuttal/E27/{units.*.jsonl,census_population_A.json,config.*.json,run.log,run.pid,DONE}, then analyses.
set -euo pipefail
WORKERS="${1:-6}"
if [ "$WORKERS" -gt 10 ]; then echo "refusing more than 10 workers" >&2; exit 1; fi
HERE="experiments/rebuttal/E27_amenability_identity_new_populations"
OUT="results/rebuttal/E27"
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
  $PY -u $HERE/run.py --out $OUT --workers $WORKERS --mem-gb 1.5; rc=\$?
  $PY -u $HERE/analyse.py --out $OUT > $OUT/analyse.txt 2>&1; arc=\$?
  echo \"run_exit=\$rc analyse_exit=\$arc finished=\$(date -Iseconds)\" > $OUT/DONE
" >> "$OUT/run.log" 2>&1 &
echo $! > "$OUT/run.pid"
echo "started pid $(cat "$OUT/run.pid"), log $OUT/run.log, done-file $OUT/DONE"
