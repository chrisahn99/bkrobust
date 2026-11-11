#!/bin/bash
# Smoke test of every rebuttal experiment, one after the other, each with a time limit. Run from the repository root.
cd "$(dirname "$0")/../../.." || exit 1
export PATH=$HOME/.local/bin:/usr/bin:/bin:$PATH
OUT=experiments/rebuttal/_smoke; : > $OUT/SUMMARY.txt
ONLY="$1"
if [ -z "$ONLY" ] || [[ " $ONLY " == *" E1 "* ]]; then
  start=$(date +%s); timeout 600 bash -c '.venv/bin/python experiments/rebuttal/E1_cpdag_check_on_lower_bound_rows/run.py --smoke && .venv/bin/python experiments/rebuttal/E1_cpdag_check_on_lower_bound_rows/analyse.py --out results/rebuttal/E1_smoke' > $OUT/E1.log 2>&1; rc=$?; end=$(date +%s)
  echo "E1 exit=$rc seconds=$((end-start))" | tee -a $OUT/SUMMARY.txt
fi
if [ -z "$ONLY" ] || [[ " $ONLY " == *" E2 "* ]]; then
  start=$(date +%s); timeout 600 bash -c '.venv/bin/python experiments/rebuttal/E2_running_time_by_graph_size/run.py --smoke && .venv/bin/python experiments/rebuttal/E2_running_time_by_graph_size/analyse.py --out results/rebuttal/E2_smoke' > $OUT/E2.log 2>&1; rc=$?; end=$(date +%s)
  echo "E2 exit=$rc seconds=$((end-start))" | tee -a $OUT/SUMMARY.txt
fi
if [ -z "$ONLY" ] || [[ " $ONLY " == *" E3 "* ]]; then
  start=$(date +%s); timeout 600 bash -c '.venv/bin/python experiments/rebuttal/E3_depth_budget_and_sat_mode/run.py --smoke && .venv/bin/python experiments/rebuttal/E3_depth_budget_and_sat_mode/analyse.py --smoke' > $OUT/E3.log 2>&1; rc=$?; end=$(date +%s)
  echo "E3 exit=$rc seconds=$((end-start))" | tee -a $OUT/SUMMARY.txt
fi
if [ -z "$ONLY" ] || [[ " $ONLY " == *" E4 "* ]]; then
  start=$(date +%s); timeout 600 bash -c '.venv/bin/python experiments/rebuttal/E4_class_enumeration_baseline/run.py --smoke && .venv/bin/python experiments/rebuttal/E4_class_enumeration_baseline/analyse.py results/rebuttal/E4_smoke' > $OUT/E4.log 2>&1; rc=$?; end=$(date +%s)
  echo "E4 exit=$rc seconds=$((end-start))" | tee -a $OUT/SUMMARY.txt
fi
if [ -z "$ONLY" ] || [[ " $ONLY " == *" E5 "* ]]; then
  start=$(date +%s); timeout 600 bash -c '.venv/bin/python experiments/rebuttal/E5_learned_cpdags_by_sample_size/run.py --smoke && .venv/bin/python experiments/rebuttal/E5_learned_cpdags_by_sample_size/analyse.py --out results/rebuttal/E5_smoke' > $OUT/E5.log 2>&1; rc=$?; end=$(date +%s)
  echo "E5 exit=$rc seconds=$((end-start))" | tee -a $OUT/SUMMARY.txt
fi
if [ -z "$ONLY" ] || [[ " $ONLY " == *" E6 "* ]]; then
  start=$(date +%s); timeout 600 bash -c '.venv/bin/python experiments/rebuttal/E6_bootstrap_stability_check/run.py --smoke && .venv/bin/python experiments/rebuttal/E6_bootstrap_stability_check/analyse.py --out results/rebuttal/E6_smoke' > $OUT/E6.log 2>&1; rc=$?; end=$(date +%s)
  echo "E6 exit=$rc seconds=$((end-start))" | tee -a $OUT/SUMMARY.txt
fi
if [ -z "$ONLY" ] || [[ " $ONLY " == *" E7 "* ]]; then
  start=$(date +%s); timeout 600 bash -c '.venv/bin/python experiments/rebuttal/E7_one_hidden_variable/run.py --smoke && .venv/bin/python experiments/rebuttal/E7_one_hidden_variable/analyse.py --out results/rebuttal/E7_smoke' > $OUT/E7.log 2>&1; rc=$?; end=$(date +%s)
  echo "E7 exit=$rc seconds=$((end-start))" | tee -a $OUT/SUMMARY.txt
fi
if [ -z "$ONLY" ] || [[ " $ONLY " == *" E8 "* ]]; then
  start=$(date +%s); timeout 600 bash -c '.venv/bin/python experiments/rebuttal/E8_missed_and_extra_edges/run.py --smoke && .venv/bin/python experiments/rebuttal/E8_missed_and_extra_edges/analyse.py --smoke' > $OUT/E8.log 2>&1; rc=$?; end=$(date +%s)
  echo "E8 exit=$rc seconds=$((end-start))" | tee -a $OUT/SUMMARY.txt
fi
if [ -z "$ONLY" ] || [[ " $ONLY " == *" E9 "* ]]; then
  start=$(date +%s); timeout 600 bash -c '.venv/bin/python experiments/rebuttal/E9_parametric_bootstrap_validity/run.py --smoke && .venv/bin/python experiments/rebuttal/E9_parametric_bootstrap_validity/analyse.py --out results/rebuttal/E9_smoke --e6 results/rebuttal/E6_smoke' > $OUT/E9.log 2>&1; rc=$?; end=$(date +%s)
  echo "E9 exit=$rc seconds=$((end-start))" | tee -a $OUT/SUMMARY.txt
fi
if [ -z "$ONLY" ] || [[ " $ONLY " == *" E10 "* ]]; then
  start=$(date +%s); timeout 600 bash -c '.venv/bin/python experiments/rebuttal/E10_radius_across_valid_sets/run.py --smoke && .venv/bin/python experiments/rebuttal/E10_radius_across_valid_sets/analyse.py --out results/rebuttal/E10_smoke' > $OUT/E10.log 2>&1; rc=$?; end=$(date +%s)
  echo "E10 exit=$rc seconds=$((end-start))" | tee -a $OUT/SUMMARY.txt
fi
if [ -z "$ONLY" ] || [[ " $ONLY " == *" E11 "* ]]; then
  start=$(date +%s); timeout 600 bash -c '.venv/bin/python experiments/rebuttal/E11_survival_on_cpdag_only_queries/run.py --smoke && .venv/bin/python experiments/rebuttal/E11_survival_on_cpdag_only_queries/analyse.py --out results/rebuttal/E11_smoke --n-boot 200 --n-perm 50' > $OUT/E11.log 2>&1; rc=$?; end=$(date +%s)
  echo "E11 exit=$rc seconds=$((end-start))" | tee -a $OUT/SUMMARY.txt
fi
if [ -z "$ONLY" ] || [[ " $ONLY " == *" E12 "* ]]; then
  start=$(date +%s); timeout 600 bash -c '.venv/bin/python experiments/rebuttal/E12_networks_above_500_nodes/run.py --smoke && .venv/bin/python experiments/rebuttal/E12_networks_above_500_nodes/analyse.py --out results/rebuttal/E12/smoke' > $OUT/E12.log 2>&1; rc=$?; end=$(date +%s)
  echo "E12 exit=$rc seconds=$((end-start))" | tee -a $OUT/SUMMARY.txt
fi
if [ -z "$ONLY" ] || [[ " $ONLY " == *" E13 "* ]]; then
  start=$(date +%s); timeout 600 bash -c 'PYTHONPATH=src .venv/bin/python experiments/rebuttal/E13_committed_set_consistency_checks/run.py --smoke && PYTHONPATH=src .venv/bin/python experiments/rebuttal/E13_committed_set_consistency_checks/analyse.py --out results/rebuttal/E13_smoke' > $OUT/E13.log 2>&1; rc=$?; end=$(date +%s)
  echo "E13 exit=$rc seconds=$((end-start))" | tee -a $OUT/SUMMARY.txt
fi
if [ -z "$ONLY" ] || [[ " $ONLY " == *" E14 "* ]]; then
  start=$(date +%s); timeout 600 bash -c '.venv/bin/python experiments/rebuttal/E14_data_test_for_invalid_sets/run.py --smoke && .venv/bin/python experiments/rebuttal/E14_data_test_for_invalid_sets/analyse.py --out results/rebuttal/E14_smoke' > $OUT/E14.log 2>&1; rc=$?; end=$(date +%s)
  echo "E14 exit=$rc seconds=$((end-start))" | tee -a $OUT/SUMMARY.txt
fi
if [ -z "$ONLY" ] || [[ " $ONLY " == *" E15 "* ]]; then
  start=$(date +%s); timeout 600 bash -c 'PYTHONPATH=src .venv/bin/python experiments/rebuttal/E15_aid_as_survival_baseline/run.py --smoke && PYTHONPATH=src .venv/bin/python experiments/rebuttal/E15_aid_as_survival_baseline/analyse.py --out results/rebuttal/E15_smoke --n-boot 200 --n-boot-secondary 50' > $OUT/E15.log 2>&1; rc=$?; end=$(date +%s)
  echo "E15 exit=$rc seconds=$((end-start))" | tee -a $OUT/SUMMARY.txt
fi
if [ -z "$ONLY" ] || [[ " $ONLY " == *" E16 "* ]]; then
  start=$(date +%s); timeout 600 bash -c 'PYTHONPATH=src .venv/bin/python experiments/rebuttal/E16_radius_where_separation_undefined/run.py --smoke && .venv/bin/python experiments/rebuttal/E16_radius_where_separation_undefined/analyse.py --out results/rebuttal/E16_smoke' > $OUT/E16.log 2>&1; rc=$?; end=$(date +%s)
  echo "E16 exit=$rc seconds=$((end-start))" | tee -a $OUT/SUMMARY.txt
fi
if [ -z "$ONLY" ] || [[ " $ONLY " == *" E17 "* ]]; then
  start=$(date +%s); timeout 600 bash -c '.venv/bin/python experiments/rebuttal/E17_worked_case_and_validity_by_radius/run.py --smoke --out results/rebuttal/E17_smoke && .venv/bin/python experiments/rebuttal/E17_worked_case_and_validity_by_radius/analyse.py --out results/rebuttal/E17_smoke' > $OUT/E17.log 2>&1; rc=$?; end=$(date +%s)
  echo "E17 exit=$rc seconds=$((end-start))" | tee -a $OUT/SUMMARY.txt
fi
if [ -z "$ONLY" ] || [[ " $ONLY " == *" E18 "* ]]; then
  start=$(date +%s); timeout 600 bash -c 'PYTHONPATH=src .venv/bin/python experiments/rebuttal/E18_section7_standardised_models/run.py --smoke && PYTHONPATH=src .venv/bin/python experiments/rebuttal/E18_section7_standardised_models/analyse.py --smoke' > $OUT/E18.log 2>&1; rc=$?; end=$(date +%s)
  echo "E18 exit=$rc seconds=$((end-start))" | tee -a $OUT/SUMMARY.txt
fi
if [ -z "$ONLY" ] || [[ " $ONLY " == *" E19 "* ]]; then
  start=$(date +%s); timeout 600 bash -c '.venv/bin/python experiments/rebuttal/E19_bias_radius_internal_standardisation/run.py --smoke && .venv/bin/python experiments/rebuttal/E19_bias_radius_internal_standardisation/analyse.py --smoke' > $OUT/E19.log 2>&1; rc=$?; end=$(date +%s)
  echo "E19 exit=$rc seconds=$((end-start))" | tee -a $OUT/SUMMARY.txt
fi
if [ -z "$ONLY" ] || [[ " $ONLY " == *" E20 "* ]]; then
  start=$(date +%s); timeout 600 bash -c 'PYTHONPATH=src .venv/bin/python experiments/rebuttal/E20_uumc_generator/run.py --smoke --workers 2 && .venv/bin/python experiments/rebuttal/E20_uumc_generator/analyse.py --out results/rebuttal/E20_smoke' > $OUT/E20.log 2>&1; rc=$?; end=$(date +%s)
  echo "E20 exit=$rc seconds=$((end-start))" | tee -a $OUT/SUMMARY.txt
fi
if [ -z "$ONLY" ] || [[ " $ONLY " == *" E21 "* ]]; then
  start=$(date +%s); timeout 600 bash -c '.venv/bin/python experiments/rebuttal/E21_walkthrough_toy/run.py --smoke --out results/rebuttal/E21_smoke && .venv/bin/python experiments/rebuttal/E21_walkthrough_toy/analyse.py --out results/rebuttal/E21_smoke' > $OUT/E21.log 2>&1; rc=$?; end=$(date +%s)
  echo "E21 exit=$rc seconds=$((end-start))" | tee -a $OUT/SUMMARY.txt
fi
touch $OUT/DONE
