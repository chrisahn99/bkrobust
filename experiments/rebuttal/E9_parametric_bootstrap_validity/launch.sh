#!/usr/bin/env bash
# Full E9 run on the workstation, detached. Run from anywhere; it works from the repository root.
#   bash experiments/rebuttal/E9_parametric_bootstrap_validity/launch.sh [WORKERS]
# Refuses to start when the done-file exists, a live pid file is present, or E6 has no full
# result (PREREG: E9 runs only after E6; set E9_SKIP_E6_GATE=1 to override, and say so).
# The run is resumable, so after a crash remove the stale pid file and launch again.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
cd "$ROOT"

WORKERS="${1:-6}"
if [ "$WORKERS" -gt 8 ] || [ "$WORKERS" -lt 1 ]; then
  echo "between 1 and 8 workers on the 16-core workstation (15 GB)" >&2
  exit 2
fi
OUT="results/rebuttal/E9"
E6="results/rebuttal/E6"
LOG="$OUT/run.log"
PID="$OUT/run.pid"
DONE="$OUT/DONE"
PY="$ROOT/.venv/bin/python"

[ -x "$PY" ] || { echo "no virtual environment at $PY" >&2; exit 2; }
[ -d experiments/rebuttal/_data/pest ] || { echo "missing experiments/rebuttal/_data/pest (chat_e.json, knowledge.json, rows.csv)" >&2; exit 2; }
if [ "${E9_SKIP_E6_GATE:-0}" != "1" ]; then
  "$PY" -c "import json,sys; r=json.load(open('$E6/results.json')); sys.exit(0 if r.get('auc') is not None and not r.get('smoke') else 1)" 2>/dev/null \
    || { echo "E6 has no full results.json with an AUC in $E6; E9 waits for it (E9_SKIP_E6_GATE=1 to override)" >&2; exit 1; }
fi
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
  '$PY' experiments/rebuttal/E9_parametric_bootstrap_validity/run.py --out '$OUT' --workers '$WORKERS' \
  && '$PY' experiments/rebuttal/E9_parametric_bootstrap_validity/analyse.py --out '$OUT' --e6 '$E6' \
  && date -u +%FT%TZ > '$DONE'
  echo \"exit \$?\" >> '$LOG'
  rm -f '$PID'
" >> "$LOG" 2>&1 &
echo $! > "$PID"
echo "started pid $(cat "$PID"), log $ROOT/$LOG, done-file $ROOT/$DONE"
