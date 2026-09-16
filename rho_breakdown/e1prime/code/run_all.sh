#!/usr/bin/env bash
# E1' full run. betelgeuse: 16 cores, no GPU needed, no network.
set -euo pipefail
PY="${PY:-$HOME/e1venv/bin/python}"
NW="${NW:-12}"
cd "$(dirname "$0")"
mkdir -p ../results ../logs

echo "=== ARM 1 : population sweep ==="
$PY run_arm1.py 4000 licensed ../results/arm1_licensed.json  "$NW"
$PY run_arm1.py 4000 original ../results/arm1_original.json  "$NW"
$PY run_arm1.py 3000 large    ../results/arm1_large.json     "$NW"
$PY run_arm1.py 1200 k8       ../results/arm1_k8.json        "$NW"

echo "=== ARM 1 : analysis ==="
$PY analyse_arm1.py ../results/arm1_licensed.json licensed K4 > ../logs/arm1_licensed.txt
$PY analyse_arm1.py ../results/arm1_original.json original K4 > ../logs/arm1_original.txt
$PY analyse_arm1.py ../results/arm1_large.json    large    K4 > ../logs/arm1_large.txt
for K in K4 K6 K8; do
  $PY analyse_arm1.py ../results/arm1_k8.json k8 "$K" > "../logs/arm1_k8_$K.txt"
done

echo "=== ARM 2 : finite-sample analyst simulation ==="
$PY run_arm2.py ../results/arm1_licensed.json ../results/arm2_licensed.json 2500 20 "$NW"
$PY analyse_arm2.py ../results/arm2_licensed.json > ../logs/arm2_licensed.txt

echo "=== DONE ==="
ls -la ../results
