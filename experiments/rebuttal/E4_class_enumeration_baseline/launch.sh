#!/usr/bin/env bash
# E4 full run on the workstation (16 cores, 15 GB). From the repository root:
#   bash experiments/rebuttal/E4_class_enumeration_baseline/launch.sh [WORKERS]
set -euo pipefail
cd "$(dirname "$0")/../../.."
WORKERS="${1:-6}"
if [ "$WORKERS" -gt 8 ]; then echo "refusing: at most 8 workers" >&2; exit 1; fi
OUT=results/rebuttal/E4
LOGDIR=$OUT/logs
mkdir -p "$LOGDIR"
DONE=$OUT/E4.done
PID=$OUT/E4.pid
if [ -e "$DONE" ]; then echo "refusing: $DONE exists (run finished)" >&2; exit 1; fi
if [ -e "$PID" ] && kill -0 "$(cat "$PID")" 2>/dev/null; then echo "refusing: live run, pid $(cat "$PID")" >&2; exit 1; fi
[ -x .venv/bin/python ] || { echo "no .venv/bin/python in $(pwd)" >&2; exit 1; }
LOG=$LOGDIR/run_$(date +%Y%m%d_%H%M%S).log
# one BLAS thread per worker: the work is pure Python, threads would only contend
nohup env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 bash -c "
  .venv/bin/python experiments/rebuttal/E4_class_enumeration_baseline/run.py --out $OUT --workers $WORKERS
  rc=\$?
  .venv/bin/python experiments/rebuttal/E4_class_enumeration_baseline/analyse.py $OUT > $OUT/analyse.txt 2>&1 || true
  echo \"exit \$rc \$(date -u +%FT%TZ)\" > $DONE
  rm -f $PID
" > "$LOG" 2>&1 &
echo $! > "$PID"
echo "started pid $(cat "$PID"), log $LOG, done-file $DONE"
