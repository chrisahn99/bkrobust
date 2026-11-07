# E24 — real far-away discovery errors transplanted into the true graph

For each learned CPDAG unit certifiable under E5's repaired oracle arm, hybrid DAGs keep the truth within d hops of {x, y} (d = 1, 2, 3; controls d = none and d = all) and the learned graph's errors elsewhere; the committed set of each hybrid is tested in the true DAG. Registered statistic and prediction: `PREREG.md`.

    .venv/bin/python experiments/rebuttal/E24_e24_far_error_transplant/run.py --smoke       # < 2 min, results/rebuttal/E24_smoke
    bash experiments/rebuttal/E24_e24_far_error_transplant/launch.sh 2           # full run, nice, 2 workers

Inputs: `experiments/rebuttal/_data/pest/` (same files as E5) and E5's finished rows in `results/rebuttal/E5/` (consistency check 1).
Outputs in `results/rebuttal/E24/`: `rows.shard*.jsonl` (one line per unit: `g0`, `anchor`, or `h|...|<d>`), `config.*.json`, `run.log`, `analyse.txt`, `results.json`, `run.done`.
Exit 2 from either program means a consistency check failed: read `consistency_failures.*.json`.
Fidelity: `python ../fidelity_gate.py FIDELITY.yaml --no-judge` after the run.
No discovery runs; closed-form sets and Meek closures only. Expected duration: about 5 min on 2 workers (smoke 2 s; all 35 arth150 graphs took 76 s on 3 cores, probe kept in `results/rebuttal/E24_timing_probe_arth150/`). Far below 30 core-hours, so no Slurm array.
