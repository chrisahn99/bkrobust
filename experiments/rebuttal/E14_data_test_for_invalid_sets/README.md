# E14: can a data test across sets valid in G0 flag an invalid committed set?

Runs the coefficient-agreement test of Su and Henckel (2022) across up to 20 sets GAC-valid in the knowledge-induced graph G0 of each distinct released state, on 50 non-Gaussian linear SEMs of the true DAG at n = 1,000 and 10,000. Registered statistic (PREREG.md): detection rate on committed sets invalid at truth minus false-alarm rate on valid ones, at n = 10,000, by released radius. Exits 1 on a consistency failure.

Data: copy `per_query_records.csv` and `knowledge_states.json` from the paper's `release/` into `experiments/rebuttal/_data/` (or pass `--records`, `--knowledge`). Networks come from `results/axisa3/networks/example_models` in this repository.

    .venv/bin/python experiments/rebuttal/E14_data_test_for_invalid_sets/run.py --smoke && \
      .venv/bin/python experiments/rebuttal/E14_data_test_for_invalid_sets/analyse.py --out results/rebuttal/E14_smoke
    bash experiments/rebuttal/E14_data_test_for_invalid_sets/launch.sh 8

Smoke: 8 states of asia and sachs, 3 replicates, under 60 s. Full run: about 2,500 states, estimated 2 to 6 h on 8 cores (dominated by the large networks: candidate-set search and n = 10,000 regressions).
Outputs in `results/rebuttal/E14/`: `units.jsonl` (one line per state, resumable), `config.json`, `run.log`, `DONE`, `results.json` (registered statistic, cluster-bootstrap interval, strata, size, population ceiling). Array jobs: `--shard I --n-shards K`.
