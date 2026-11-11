#!/bin/bash
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 PYTHONPATH=src:rebuttal
mkdir -p results/rebuttal/e3s
caffeinate -dimsu -w $$ &
.venv/bin/python rebuttal/e3s_solvers.py --workers 8 --pops P-SATR,P-E3,P-SCALE > results/rebuttal/e3s/run_main.log 2>&1
date > results/rebuttal/e3s/DONE
