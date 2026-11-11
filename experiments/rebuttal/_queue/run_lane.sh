#!/bin/bash
# One lane of the rebuttal queue on the workstation: runs the given experiments one after the other through their own
# launchers, pinned to a set of cores, and waits for each done-file. Usage, from the repository root:
#   bash experiments/rebuttal/_queue/run_lane.sh <lane> <cpu-list> <workers> E5 E6 ...
# Writes experiments/rebuttal/_queue/<lane>.log and one line per experiment in STATUS.txt. Never deletes anything.
cd "$(dirname "$0")/../../.." || exit 1
export PATH=$HOME/.local/bin:$PATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH=src
LANE="$1"; CPUS="$2"; W="$3"; shift 3
Q=experiments/rebuttal/_queue
for E in "$@"; do
  DIR=$(ls -d experiments/rebuttal/${E}_* 2>/dev/null | head -1)
  OUT=results/rebuttal/$E
  if [ -z "$DIR" ]; then echo "$(date -Is) $LANE $E no folder" >> $Q/STATUS.txt; continue; fi
  if [ -n "$(find $OUT -maxdepth 2 \( -name DONE -o -name '*.done' \) 2>/dev/null | head -1)" ]; then echo "$(date -Is) $LANE $E already done, skipped" >> $Q/STATUS.txt; continue; fi
  # wait for memory: at least 4 GB available before starting another experiment
  while [ "$(awk '/MemAvailable/{print int($2/1048576)}' /proc/meminfo)" -lt 4 ]; do sleep 60; done
  start=$(date +%s); echo "$(date -Is) $LANE $E start workers=$W cpus=$CPUS" >> $Q/STATUS.txt
  WORKERS=$W taskset -c "$CPUS" bash "$DIR/launch.sh" "$W" >> $Q/$LANE.log 2>&1
  sleep 5
  while [ -z "$(find $OUT -maxdepth 2 \( -name DONE -o -name '*.done' \) 2>/dev/null | head -1)" ]; do
    # stop waiting if the launcher never started anything (no pid file and no log after two minutes)
    if [ $(( $(date +%s) - start )) -gt 120 ] && [ -z "$(find $OUT -maxdepth 2 -name '*.pid' 2>/dev/null | head -1)" ]; then echo "$(date -Is) $LANE $E launcher did not start (see $Q/$LANE.log)" >> $Q/STATUS.txt; break; fi
    sleep 30
  done
  DF=$(find $OUT -maxdepth 2 \( -name DONE -o -name '*.done' \) 2>/dev/null | head -1)
  echo "$(date -Is) $LANE $E end seconds=$(( $(date +%s) - start )) done-file: $(cat "$DF" 2>/dev/null | tr '\n' ' ')" >> $Q/STATUS.txt
done
echo "$(date -Is) $LANE lane finished" >> $Q/STATUS.txt
