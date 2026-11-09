#!/usr/bin/env bash
# Full E29 run on the workstation, detached, niced. Run from the repository root:
#   bash experiments/rebuttal/E29_radius_vs_degree_sep_undefined/launch.sh [workers] [cpu-list]
# Writes results/rebuttal/E29/{units.*.jsonl,config.*.json,consistency_failures.*,errors.*,run.log,run.pid,
# results.json,analyse.txt,DONE}. Never deletes anything; refuses to start over a live or finished run.
set -euo pipefail
WORKERS="${1:-2}"
CPUS="${2:-}"
if [ "$WORKERS" -gt 6 ]; then echo "refusing more than 6 workers" >&2; exit 1; fi
HERE="experiments/rebuttal/E29_radius_vs_degree_sep_undefined"
OUT="results/rebuttal/E29"
PY=".venv/bin/python"
[ -f "$HERE/run.py" ] || { echo "run from the repository root" >&2; exit 1; }
[ -x "$PY" ] || { echo "no $PY" >&2; exit 1; }
for f in results/axis_robustness_llm/analysis_units.csv results/radius_gac_recompute/axis_robustness_llm.csv \
         results/radius_gac_recompute/axis_robustness_llm_large.csv results/elicit/knowledge.json \
         results/axisa3/networks/acquisition_manifest.json; do
  [ -f "$f" ] || { echo "missing input $f" >&2; exit 1; }
done
mkdir -p "$OUT"
if [ -e "$OUT/DONE" ]; then echo "$OUT/DONE exists: already finished" >&2; exit 1; fi
if [ -f "$OUT/run.pid" ] && kill -0 "$(cat "$OUT/run.pid")" 2>/dev/null; then
  echo "a run is live (pid $(cat "$OUT/run.pid"))" >&2; exit 1
fi
PIN=""
[ -n "$CPUS" ] && PIN="taskset -c $CPUS"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH=src
nohup $PIN nice -n 10 bash -c "
  $PY -u $HERE/run.py --out $OUT --workers $WORKERS --mem-gb 1.5; rc=\$?
  $PY -u $HERE/analyse.py --out $OUT > $OUT/analyse.log 2>&1; arc=\$?
  echo \"run_exit=\$rc analyse_exit=\$arc finished=\$(date -Iseconds)\" > $OUT/DONE
" >> "$OUT/run.log" 2>&1 &
echo $! > "$OUT/run.pid"
echo "started pid $(cat "$OUT/run.pid"), log $OUT/run.log, done-file $OUT/DONE"
