# E22: which discovery error breaks each committed set

For the 987 certifiable A_ORACLE_R rows of E5 on learned graphs (534 invalid at the truth), rebuilds G0 and Z with E5's functions (exit 2 on any mismatch with the stored row), lists every difference between the learned CPDAG and the true CPDAG with its hop distance from {x, y}, names the failing GAC condition at the truth, and counts invalid rows that satisfy the S3 hypothesis inside An({x, y} u Z). Plan and registered statistic: `PREREG.md`.

Needs the finished E5 run (`results/rebuttal/E5/rows.shard0of1.jsonl`, `run.done`) and `experiments/rebuttal/_data/pest/`. No discovery is run.

    .venv/bin/python experiments/rebuttal/E22_e22_discovery_error_census/run.py --smoke && \
    .venv/bin/python experiments/rebuttal/E22_e22_discovery_error_census/analyse.py --out results/rebuttal/E22_smoke --n-boot 200
    bash experiments/rebuttal/E22_e22_discovery_error_census/launch.sh 2

Duration: smoke under a minute; full run well under 1 h on 2 workers (see the smoke timings), far below 30 core-hours, so no Slurm array.
Outputs in `results/rebuttal/E22/`: `rows.shard*.jsonl` (one row per query, with the full error list `errors`: type, nodes, hop, in_A, for valid and invalid rows; reusable by E24 and the local-correctness tabulation), `config.*.json`, `run.log`, `analyse.txt`, `results.json`, `run.done`.
Fidelity: `python ../fidelity_gate.py FIDELITY.yaml --no-judge` after the run.
