#!/bin/bash
#SBATCH --job-name=bk-elicit
#SBATCH --account=zkh@h100
#SBATCH --constraint=h100
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-task=24
#SBATCH --hint=nomultithread
#SBATCH --time=00:45:00
#SBATCH --qos=qos_gpu_h100-dev
#SBATCH --output=%x-%j.out
# usage: sbatch elicit_jz.slurm <checkpoint dir> <model id> [--no-think]
set -euo pipefail
module purge
module load arch/h100
module load pytorch-gpu/py3/2.8.0
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 VLLM_LOGGING_LEVEL=WARNING
# the home directory has a 3 GB quota and vLLM's compile cache overflowed it once
export XDG_CACHE_HOME="$SCRATCH/cache" VLLM_CACHE_ROOT="$SCRATCH/cache/vllm"
export TRITON_CACHE_DIR="$SCRATCH/cache/triton" TORCHINDUCTOR_CACHE_DIR="$SCRATCH/cache/inductor"
mkdir -p "$SCRATCH/cache"
cd "$SLURM_SUBMIT_DIR"
mkdir -p "${OUTDIR:-out}"
python elicit_vllm.py --model "$1" --model-id "$2" --prompts "${PROMPTS:-prompts.jsonl}" \
  --out "${OUTDIR:-out}/$2.jsonl" --tp 4 --eager ${3:-}
