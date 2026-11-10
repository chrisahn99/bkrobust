#!/bin/bash
# Overnight E4L: main pass, then the long-cap validation subsample. Keeps the Mac awake.
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 PYTHONPATH=src:rebuttal
mkdir -p results/rebuttal/e4l
caffeinate -dimsu -w $$ &
.venv/bin/python rebuttal/e4l_listing.py --workers 8 > results/rebuttal/e4l/run_main.log 2>&1
.venv/bin/python rebuttal/e4l_listing.py --long --workers 8 > results/rebuttal/e4l/run_long.log 2>&1
date > results/rebuttal/e4l/DONE
