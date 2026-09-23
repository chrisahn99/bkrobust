#!/usr/bin/env bash
# Survival/ranking sweep on ELICITED background knowledge, across the model panel.
#
# K comes from results/elicit/knowledge.json -- what each model actually
# asserted when asked about the CPDAG's undirected pairs. The programmatic
# benchmarks.measure.select_knowledge path is not used anywhere in this sweep.
#
# The variation axis is the model panel: eight real-naming conditions from a
# 7B up to two 70B/72B models, plus two scrambled-naming controls that are run
# and reported but never pooled with the panel. Each (condition, network) cell
# is one shard and carries its own |K| and its own measured accuracy.
#
# Usage:
#   scripts/run_llm_survival_panel.sh                  # full panel + control
#   scripts/run_llm_survival_panel.sh --workers 4
#   NAMING=real scripts/run_llm_survival_panel.sh      # panel only
#
# Re-invoking is safe and skips completed shards: a shard writes its
# completion marker only when it finishes, and a shard with no marker is
# redone from scratch rather than resumed.

set -euo pipefail

cd "$(dirname "$0")/.."

PY="${PY:-.venv/bin/python}"
OUT_DIR="${OUT_DIR:-results/axis_robustness_llm}"
ELICIT="${ELICIT:-results/elicit/knowledge.json}"
N_DRAWS="${N_DRAWS:-1000}"
WORKERS="${WORKERS:-8}"
N_BOOT="${N_BOOT:-10000}"
LOG_DIR="$OUT_DIR/_logs"

# Optional naming filter: "real" for the panel alone, "scrambled" for the
# control alone, unset for both. macOS ships bash 3.2, where expanding an empty
# array under `set -u` is an error, so every array expansion below uses the
# `${a[@]+"${a[@]}"}` form.
NAMING_ARG=()
if [[ -n "${NAMING:-}" ]]; then
  NAMING_ARG=(--naming "$NAMING")
fi

# Pass-through flags (e.g. --workers 4) override the environment defaults.
EXTRA=("$@")

export PYTHONPATH=src

mkdir -p "$LOG_DIR"
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
MAIN_LOG="$LOG_DIR/panel_${STAMP}.log"

log() { printf '%s | %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a "$MAIN_LOG"; }

log "=== elicited-knowledge survival panel ==="
log "out_dir=$OUT_DIR elicit=$ELICIT n_draws=$N_DRAWS workers=$WORKERS"
log "git=$(git rev-parse HEAD) branch=$(git rev-parse --abbrev-ref HEAD)"
log "python=$($PY --version 2>&1)"

log "--- stage 1/4: manifest ---"
$PY -m bkrobust.robustness.run_llm_survival manifest \
    --out-dir "$OUT_DIR" --elicit "$ELICIT" --n-draws "$N_DRAWS" 2>&1 | tee -a "$MAIN_LOG"

log "--- stage 2/4: pending shards ---"
$PY -m bkrobust.robustness.run_llm_survival status \
    --out-dir "$OUT_DIR" --elicit "$ELICIT" ${NAMING_ARG[@]+"${NAMING_ARG[@]}"} \
    2>&1 | tee -a "$MAIN_LOG"

log "--- stage 3/4: sweep (each shard is its own process; a failure is logged, never fatal) ---"
SWEEP_RC=0
$PY -m bkrobust.robustness.run_llm_survival pool \
    --out-dir "$OUT_DIR" --elicit "$ELICIT" --n-draws "$N_DRAWS" \
    --workers "$WORKERS" ${NAMING_ARG[@]+"${NAMING_ARG[@]}"} ${EXTRA[@]+"${EXTRA[@]}"} \
    2>&1 | tee -a "$MAIN_LOG" || SWEEP_RC=$?
log "sweep exit status: $SWEEP_RC (nonzero means at least one shard failed; see $LOG_DIR)"

log "--- stage 4/4: analysis ---"
$PY -m bkrobust.robustness.llm_analyse \
    --out-dir "$OUT_DIR" --n-boot "$N_BOOT" 2>&1 | tee -a "$MAIN_LOG"

log "--- final status ---"
$PY -m bkrobust.robustness.run_llm_survival status \
    --out-dir "$OUT_DIR" --elicit "$ELICIT" 2>&1 | tee -a "$MAIN_LOG"
log "done. log: $MAIN_LOG"
