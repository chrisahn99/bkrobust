#!/usr/bin/env bash
# Second pass of E27 (DEVIATIONS.md item 1). Run from the repository root:
#   nohup bash experiments/rebuttal/E27_amenability_identity_new_populations/resume_after_done.sh >> results/rebuttal/E27/run.log 2>&1 &
set -uo pipefail
HERE="experiments/rebuttal/E27_amenability_identity_new_populations"
OUT="results/rebuttal/E27"
PY=".venv/bin/python"
until [ -e "$OUT/DONE" ]; do sleep 30; done
cp "$OUT/analyse.txt" "$OUT/analyse_pass1.txt" 2>/dev/null
cp "$OUT/results.json" "$OUT/results_pass1.json" 2>/dev/null
echo "pass 2 start $(date -Iseconds)"
$PY -u $HERE/run.py --out $OUT --workers 2 --mem-gb 4; rc=$?
$PY -u $HERE/analyse.py --out $OUT > $OUT/analyse.txt 2>&1; arc=$?
echo "run_exit=$rc analyse_exit=$arc finished=$(date -Iseconds)" > $OUT/DONE2
