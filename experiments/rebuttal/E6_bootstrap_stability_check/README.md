# E6: bootstrap stability as a companion check for discovery error

Measures whether the share of bootstrap resamples (100 per network, PC rerun on each, same knowledge re-imposed) in which the committed set Z stays GAC-valid flags the 244 certifiable rows of the learned-CPDAG panel whose Z is invalid in the true DAG. Registered statistic: AUC of instability for invalidity at truth, with a network-then-query cluster bootstrap interval (PREREG.md).

Needs `experiments/rebuttal/_data/pest/{chat_e.json,knowledge.json,rows.csv}` (from the second author's `results/pest/`) and causal-learn 0.1.4.8.

    .venv/bin/python experiments/rebuttal/E6_bootstrap_stability_check/run.py --smoke && .venv/bin/python experiments/rebuttal/E6_bootstrap_stability_check/analyse.py --out results/rebuttal/E6_smoke
    bash experiments/rebuttal/E6_bootstrap_stability_check/launch.sh 6

Smoke: ecoli70 and magic-niab, 2 replicates plus identity, under a minute. Full: 404 PC runs, about 2 h on 6 workers (arth150 dominates, 80-230 s per run).
Outputs in `results/rebuttal/E6/`: `config.json`, `reference.jsonl` (rebuilt rows), `units_*.jsonl` (one line per network x replicate), `results.json`, `per_row.json`, `run.log`, `DONE`.
Exit 3 = rebuilt rows disagree with the frozen panel; exit 4 = identity replicate does not reproduce the reference CPDAG. Split with `--shard I --n-shards K`.
