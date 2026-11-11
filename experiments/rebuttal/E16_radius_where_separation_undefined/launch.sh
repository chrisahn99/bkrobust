#!/usr/bin/env bash
# Full E16 run on the workstation, detached. Run from the repository root:
#   bash experiments/rebuttal/E16_radius_where_separation_undefined/launch.sh [workers]
# Writes results/rebuttal/E16/{rows.jsonl,config.json,results.json,run.log,run.pid,DONE}.
set -euo pipefail
WORKERS="${1:-4}"
if [ "$WORKERS" -gt 8 ]; then echo "refusing more than 8 workers" >&2; exit 1; fi
HERE="experiments/rebuttal/E16_radius_where_separation_undefined"
OUT="results/rebuttal/E16"
PY=".venv/bin/python"
[ -f "$HERE/run.py" ] || { echo "run from the repository root" >&2; exit 1; }
[ -x "$PY" ] || { echo "no $PY" >&2; exit 1; }
for f in results/axis_robustness_llm/analysis_units.csv results/axis_robustness_llm_v2/empty_valid_units.csv \
         results/radius_gac_recompute/axis_robustness_llm.csv results/radius_gac_recompute/axis_robustness_llm_large.csv \
         results/radius_gac_recompute/axis_robustness_llm_v2_empty_valid.csv; do
  [ -f "$f" ] || { echo "missing data file $f" >&2; exit 1; }
done
mkdir -p "$OUT"
if [ -e "$OUT/DONE" ]; then echo "$OUT/DONE exists: already finished" >&2; exit 1; fi
if [ -f "$OUT/run.pid" ] && kill -0 "$(cat "$OUT/run.pid")" 2>/dev/null; then
  echo "a run is live (pid $(cat "$OUT/run.pid"))" >&2; exit 1
fi
nohup bash -c "
  PYTHONPATH=src $PY -u $HERE/run.py --out $OUT --workers $WORKERS; rc=\$?
  $PY -u $HERE/analyse.py --out $OUT; arc=\$?
  echo \"run_exit=\$rc analyse_exit=\$arc finished=\$(date -Iseconds)\" > $OUT/DONE
" >> "$OUT/run.log" 2>&1 &
echo $! > "$OUT/run.pid"
echo "started pid $(cat "$OUT/run.pid"), log $OUT/run.log, done-file $OUT/DONE"
