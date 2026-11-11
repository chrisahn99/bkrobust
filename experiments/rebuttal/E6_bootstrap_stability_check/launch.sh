#!/usr/bin/env bash
# Full E6 run on the workstation, detached. Run from anywhere; it works from the repository root.
#   bash experiments/rebuttal/E6_bootstrap_stability_check/launch.sh [WORKERS]
# Refuses to start when the done-file exists or a live pid file is present. The run is
# resumable, so after a crash remove the stale pid file and launch again.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
cd "$ROOT"

WORKERS="${1:-6}"
if [ "$WORKERS" -gt 8 ]; then
  echo "at most 8 workers on the 16-core workstation (15 GB)" >&2
  exit 2
fi
OUT="results/rebuttal/E6"
LOG="$OUT/run.log"
PID="$OUT/run.pid"
DONE="$OUT/DONE"
PY="$ROOT/.venv/bin/python"

[ -x "$PY" ] || { echo "no virtual environment at $PY" >&2; exit 2; }
[ -d experiments/rebuttal/_data/pest ] || { echo "missing experiments/rebuttal/_data/pest (chat_e.json, knowledge.json, rows.csv)" >&2; exit 2; }
mkdir -p "$OUT"
if [ -e "$DONE" ]; then
  echo "$DONE exists: the run already finished" >&2
  exit 1
fi
if [ -e "$PID" ] && kill -0 "$(cat "$PID")" 2>/dev/null; then
  echo "a run is live (pid $(cat "$PID"))" >&2
  exit 1
fi

nohup bash -c "
  '$PY' experiments/rebuttal/E6_bootstrap_stability_check/run.py --out '$OUT' --workers '$WORKERS' \
  && '$PY' experiments/rebuttal/E6_bootstrap_stability_check/analyse.py --out '$OUT' \
  && date -u +%FT%TZ > '$DONE'
  echo \"exit \$?\" >> '$LOG'
  rm -f '$PID'
" >> "$LOG" 2>&1 &
echo $! > "$PID"
echo "started pid $(cat "$PID"), log $ROOT/$LOG, done-file $ROOT/$DONE"
