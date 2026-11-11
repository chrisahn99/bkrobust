# E5 — validity at truth and radius on learned CPDAGs, by sample size

Measures, for committed adjustment sets read off learned CPDAGs (PC/GES, n = 200/1000/5000, every stored seed) and the true CPDAG, the share that is invalid in the true DAG, plus the GAC radius, per knowledge arm. Phase A first re-scores all 560 stored P-EST rows with the submitted code and must reproduce them. Registered statistic and prediction: `PREREG.md`.

Data: the panel files from `results/pest/` of bkrobust (main branch, commit ba7ddc12) are committed in `experiments/rebuttal/_data/pest/`.

    .venv/bin/python experiments/rebuttal/E5_learned_cpdags_by_sample_size/run.py --smoke      # < 1 min, results/rebuttal/E5_smoke
    bash experiments/rebuttal/E5_learned_cpdags_by_sample_size/launch.sh            # full run, 6 workers

Expected duration of the full run: about 1 h on 6 workers. The estimate comes from the stored sweep (560 rows in 320 s) scaled to about 35,000 rows; the per-row wall bounds the worst case.
Outputs in `results/rebuttal/E5/`: `rows.shard*.jsonl` (one row per graph x arm x query), `config.*.json`, `run.log`, `analyse.txt`, `results.json`, `run.done`.
Then `analyse.py --out results/rebuttal/E5`. Exit 2 from either program means a consistency check failed: read `consistency_failures.*.json`.
Fidelity: `python ../fidelity_gate.py FIDELITY.yaml --no-judge` after the run.
