#!/bin/bash
# Smoke test of one experiment on a Slurm compute node (never on the login node). From the repository root:
#   mkdir -p logs
#   sbatch --job-name=smoke_E11 --export=ALL,EXP=E11 experiments/rebuttal/_slurm/job_smoke.slurm
# BKR_ROOT and BKR_VENV as in job_array.slurm. Partition, time and memory are examples.
#SBATCH --partition=cpu_short
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --output=logs/%x_%j.out
set -euo pipefail
R="${BKR_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
if [ -n "${BKR_VENV:-}" ]; then
  PY="$BKR_VENV/bin/python"
else
  T=${TMPDIR:-/tmp}/rb_${SLURM_JOB_ID}; mkdir -p "$T"; tar xzf "$R/venv.tar.gz" -C "$T"
  trap 'rm -rf "$T"' EXIT
  PY="$T/venv/bin/python"
fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH="$R/src"
cd "$R"
DIR=$(ls -d experiments/rebuttal/${EXP}_* | head -1)
echo "smoke of $EXP, python $("$PY" --version), commit $(git rev-parse --short HEAD)"
"$PY" -u "$DIR/run.py" --smoke --out "results/rebuttal/${EXP}_smoke_slurm" --workers 4
echo "smoke exit $?"
