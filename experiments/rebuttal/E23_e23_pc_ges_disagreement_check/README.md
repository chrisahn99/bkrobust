# E23: PC and GES agreement as a check for discovery error

Checks whether disagreement between the PC and GES states learned on the same sample flags committed sets that are invalid in the true DAG (E6's bootstrap share gave AUC 0.245). For each certifiable A_ORACLE_R row of E5 on a learned graph, Z is tested with the GAC in the other method's state, same knowledge rule. Registered statistic: AUC of that flag for invalidity at truth, pooled over n and both directions, network-then-query cluster interval, bar 0.65 (PREREG.md).

Needs `experiments/rebuttal/_data/pest/` and the finished E5 run in `results/rebuttal/E5/` (rebuild check). No discovery is rerun.

    .venv/bin/python experiments/rebuttal/E23_e23_pc_ges_disagreement_check/run.py --smoke && .venv/bin/python experiments/rebuttal/E23_e23_pc_ges_disagreement_check/analyse.py --out results/rebuttal/E23_smoke
    bash experiments/rebuttal/E23_e23_pc_ges_disagreement_check/launch.sh 2

Smoke: ecoli70 and magic-niab, seed 0, n 200 and 5000, 27 rows, 2 s. Full: 987 rows (arth150 PC seeds without a GES graph are counted, not paired), about 1 to 2 minutes on 1 core. No Slurm array is needed.
Outputs in `results/rebuttal/E23/`: `rows.*.jsonl` (one line per row), `config.*.json`, `e5_check.*.json`, `run.log`, `analyse.txt`, `results.json`, `DONE`.
Exit 2 from either program means a consistency check failed: read `consistency_failures.*.json`.
Fidelity: `python ../fidelity_gate.py FIDELITY.yaml --no-judge` after the run.
