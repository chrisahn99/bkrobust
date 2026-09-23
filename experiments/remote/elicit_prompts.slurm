#!/bin/bash
#SBATCH --job-name=bk-elicit
#SBATCH --account=<your-account>
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-task=24
#SBATCH --hint=nomultithread
#SBATCH --time=00:45:00
#SBATCH --output=%x-%j.out
# usage: sbatch elicit_prompts.slurm <checkpoint dir> <model id> [--no-think]
set -euo pipefail
# Load a Python environment with vLLM and PyTorch here (site-specific).
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 VLLM_LOGGING_LEVEL=WARNING
# keep vLLM / Triton / Inductor caches off a small home quota
export XDG_CACHE_HOME="${SCRATCH:-$TMPDIR}/cache" VLLM_CACHE_ROOT="${SCRATCH:-$TMPDIR}/cache/vllm"
export TRITON_CACHE_DIR="${SCRATCH:-$TMPDIR}/cache/triton" TORCHINDUCTOR_CACHE_DIR="${SCRATCH:-$TMPDIR}/cache/inductor"
mkdir -p "${SCRATCH:-$TMPDIR}/cache"
cd "$SLURM_SUBMIT_DIR"
mkdir -p "${OUTDIR:-out}"
python elicit_vllm.py --model "$1" --model-id "$2" --prompts "${PROMPTS:-prompts.jsonl}" \
  --out "${OUTDIR:-out}/$2.jsonl" --tp 4 --eager ${3:-}
