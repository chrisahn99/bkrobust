#!/bin/bash
cd /home/costaj/latent-causal/x2-locality/code
PY=/home/costaj/e1venv/bin/python
while pgrep -f "run_scan.py [0-9]" >/dev/null 2>&1; do sleep 10; done
mkdir -p ../results/pass1
cp ../results/scan_*.json ../results/pass1/ 2>/dev/null
echo "PASS1 archived $(date)"
# identical seed and identical n_scm as pass 1 -> identical SCMs, plus the
# add/rep operator dimension.  Reproduction of pass 1 is therefore exact.
$PY run_scan.py 28000 licensed ../results/scan_licensed.json 12
$PY run_scan.py 10000 large    ../results/scan_large.json    12
echo "OPSPLIT_DONE $(date)"
$PY run_lemma.py 6 sample 2500 ../results/lemma_p6.json 12
echo "LEMMA6_DONE $(date)"
$PY analyse_x2.py ../results ../results/x2_analysis.json > ../logs/analyse.log 2>&1
echo "ANALYSE_DONE $(date)"
