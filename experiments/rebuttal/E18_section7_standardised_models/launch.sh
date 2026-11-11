#!/usr/bin/env bash
# Full E18 run on the workstation, detached. Run from the repository root:
#   bash experiments/rebuttal/E18_section7_standardised_models/launch.sh [workers]
# Writes results/rebuttal/E18/{units.s0of1.jsonl,config.s0of1.json,failures.s0of1.jsonl,
# run.log,run.pid,DONE,results.json}. CPU only; 16 cores / 15 GB, so at most 8 workers.
set -euo pipefail
MAX_WORKERS=8
WORKERS="${1:-6}"
if [ "$WORKERS" -gt "$MAX_WORKERS" ]; then echo "refusing more than $MAX_WORKERS workers" >&2; exit 1; fi
HERE="experiments/rebuttal/E18_section7_standardised_models"
OUT="results/rebuttal/E18"
PY=".venv/bin/python"
[ -f "$HERE/run.py" ] || { echo "run from the repository root" >&2; exit 1; }
[ -x "$PY" ] || { echo "no $PY" >&2; exit 1; }
for f in results/mean_table_relevant_v2/x_component/instances.jsonl \
         results/mean_table_relevant_v2/extended/instances.jsonl \
         results/elicit/knowledge.json; do
  [ -f "$f" ] || { echo "missing input $f" >&2; exit 1; }
done
mkdir -p "$OUT"
if [ -e "$OUT/DONE" ]; then echo "$OUT/DONE exists: already finished" >&2; exit 1; fi
if [ -f "$OUT/run.pid" ] && kill -0 "$(cat "$OUT/run.pid")" 2>/dev/null; then
  echo "a run is live (pid $(cat "$OUT/run.pid"))" >&2; exit 1
fi
nohup bash -c "
  export PYTHONPATH=src OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
  $PY -u $HERE/run.py --out $OUT --workers $WORKERS; rc=\$?
  $PY -u $HERE/analyse.py --out $OUT; arc=\$?
  echo \"run_exit=\$rc analyse_exit=\$arc finished=\$(date -Iseconds)\" > $OUT/DONE
" >> "$OUT/run.log" 2>&1 &
echo $! > "$OUT/run.pid"
echo "started pid $(cat "$OUT/run.pid"), log $OUT/run.log, done-file $OUT/DONE"
