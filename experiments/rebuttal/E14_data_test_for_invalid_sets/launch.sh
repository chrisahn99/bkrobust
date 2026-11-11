#!/usr/bin/env bash
# Full E14 run on the workstation, detached. Run from the repository root:
#   bash experiments/rebuttal/E14_data_test_for_invalid_sets/launch.sh [workers]
# Writes results/rebuttal/E14/{units.jsonl,config.json,run.log,run.pid,DONE}, then analyses.
set -euo pipefail
WORKERS="${1:-8}"
if [ "$WORKERS" -gt 8 ]; then echo "refusing more than 8 workers" >&2; exit 1; fi
HERE="experiments/rebuttal/E14_data_test_for_invalid_sets"
OUT="results/rebuttal/E14"
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
# One BLAS thread per worker: the workers are the parallelism, and 8 x 16 threads would thrash.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
nohup bash -c "
  $PY -u $HERE/run.py --out $OUT --workers $WORKERS; rc=\$?
  $PY -u $HERE/analyse.py --out $OUT; arc=\$?
  echo \"run_exit=\$rc analyse_exit=\$arc finished=\$(date -Iseconds)\" > $OUT/DONE
" >> "$OUT/run.log" 2>&1 &
echo $! > "$OUT/run.pid"
echo "started pid $(cat "$OUT/run.pid"), log $OUT/run.log, done-file $OUT/DONE"
