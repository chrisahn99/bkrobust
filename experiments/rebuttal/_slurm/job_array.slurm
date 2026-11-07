#!/bin/bash
# Array job for one rebuttal experiment on a Slurm cluster. Submit from the repository root with:
#   mkdir -p logs
#   sbatch --job-name=rb_E11 --array=0-23 --export=ALL,EXP=E11,NSHARDS=24,EXTRA="" experiments/rebuttal/_slurm/job_array.slurm
# BKR_ROOT is the repository checkout on the cluster (default: the submit directory).
# The Python environment is either BKR_VENV (a directory with bin/python) or, if unset,
# $BKR_ROOT/venv.tar.gz unpacked to the node's temporary directory (useful under a small file quota).
# Partition, time and memory below are examples: adjust them to your cluster.
#SBATCH --partition=cpu_med
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=28G
#SBATCH --output=logs/%x_%A_%a.out
set -euo pipefail
R="${BKR_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
if [ -n "${BKR_VENV:-}" ]; then
  PY="$BKR_VENV/bin/python"
else
  T=${TMPDIR:-/tmp}/rb_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-0}; mkdir -p "$T"; tar xzf "$R/venv.tar.gz" -C "$T"
  trap 'rm -rf "$T"' EXIT
  PY="$T/venv/bin/python"
fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH="$R/src"
cd "$R"
DIR=$(ls -d experiments/rebuttal/${EXP}_* | head -1)
echo "exp $EXP shard ${SLURM_ARRAY_TASK_ID:-0}/${NSHARDS:-1} commit $(git rev-parse --short HEAD) start $(date -Is)"
"$PY" -u "$DIR/run.py" --out "results/rebuttal/$EXP" --workers "${SLURM_CPUS_PER_TASK:-8}" --shard "${SLURM_ARRAY_TASK_ID:-0}" --n-shards "${NSHARDS:-1}" ${EXTRA:-}
echo "exit $? end $(date -Is)"
