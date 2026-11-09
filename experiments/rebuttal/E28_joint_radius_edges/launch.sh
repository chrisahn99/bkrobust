#!/usr/bin/env bash
# E28 full run on the workstation: six shards, one core each (taskset), nice 10, a 2 GB address-space cap
# per shard, a log per shard, a done-file per shard, then analyse.py and run.done. Run from anywhere.
#   bash experiments/rebuttal/E28_joint_radius_edges/launch.sh [CORES]      # default CORES=2,3,4,5,6,7
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
cd "$REPO"

CORES="${1:-${CORES:-2,3,4,5,6,7}}"
IFS=',' read -r -a CORE_LIST <<< "$CORES"
N="${#CORE_LIST[@]}"
if [ "$N" -gt 6 ]; then echo "refusing: at most 6 cores" >&2; exit 1; fi
OUT="${OUT:-results/rebuttal/E28}"
mkdir -p "$OUT"
LOG="$OUT/run.log"; PID="$OUT/run.pid"; DONE="$OUT/run.done"
if [ -e "$DONE" ]; then echo "refusing: $DONE exists (finished run); move it away to rerun" >&2; exit 1; fi
if [ -e "$PID" ] && kill -0 "$(cat "$PID")" 2>/dev/null; then
  echo "refusing: run already live, pid $(cat "$PID")" >&2; exit 1
fi
if [ ! -x .venv/bin/python ]; then echo "refusing: .venv/bin/python missing in $REPO" >&2; exit 1; fi
if [ ! -f results/rebuttal/E22/rows.shard0of1.jsonl ] || [ ! -f results/rebuttal/E22/run.done ]; then
  echo "refusing: E22 rows missing or E22 not finished (results/rebuttal/E22/)" >&2; exit 1
fi
grep -q "E28_joint_radius_edges/PREREG.md" experiments/rebuttal/_queue/PREREG_SHA256.txt || {
  echo "refusing: PREREG.md hash not in experiments/rebuttal/_queue/PREREG_SHA256.txt" >&2; exit 1; }
want="$(grep "E28_joint_radius_edges/PREREG.md" experiments/rebuttal/_queue/PREREG_SHA256.txt | tail -1 | cut -d' ' -f1)"
have="$(sha256sum "$HERE/PREREG.md" | cut -d' ' -f1)"
if [ "$want" != "$have" ]; then echo "refusing: PREREG.md changed since its hash was recorded" >&2; exit 1; fi

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
nohup bash -c "
  for i in \$(seq 0 $((N - 1))); do
    core=\$(echo '$CORES' | cut -d, -f\$((i + 1)))
    ( nice -n 10 taskset -c \$core .venv/bin/python '$HERE/run.py' --shard \$i --n-shards $N --out '$OUT' --mem-gb 2 \
        > '$OUT/run.shard'\$i'.log' 2>&1; echo \"exit=\$? \$(date -Is)\" > '$OUT/shard'\$i'.done' ) &
  done
  wait
  nice -n 10 .venv/bin/python '$HERE/analyse.py' --out '$OUT' > '$OUT/analyse.txt' 2>&1; ra=\$?
  echo \"shards: \$(cat '$OUT'/shard*.done | tr '\n' ' ') analyse_exit=\$ra \$(date -Is)\" > '$DONE'
  rm -f '$PID'
" > "$LOG" 2>&1 &
echo $! > "$PID"
echo "E28 started, pid $(cat "$PID"), $N shards on cores $CORES, logs $OUT/run.shard*.log, done-file $DONE"
