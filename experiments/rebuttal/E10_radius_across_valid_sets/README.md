# E10: the breakdown radius across all valid adjustment sets of a query (criterion gac)

Measures, for each certifiable released query (3,014 rows: ten elicited arms and D_RAND) and for every 4- and 5-node CPDAG state, how much the radius varies across the sets valid in G0 and how often a valid set is certified strictly more robust than the committed optimal set O (registered statistic, PREREG.md). Every set also gets its validity at the true DAG. Exit 1 on a consistency failure (hybrid vs exact BFS, ladder vs search, the two truth oracles, or O against the released row).

Data: copy `per_query_records.csv` and `knowledge_states.json` from the paper's `release/` into `experiments/rebuttal/_data/` (or pass `--records`, `--knowledge`). Networks come from `results/axisa3/networks/example_models` in this repository.

    .venv/bin/python experiments/rebuttal/E10_radius_across_valid_sets/run.py --smoke && \
      .venv/bin/python experiments/rebuttal/E10_radius_across_valid_sets/analyse.py --out results/rebuttal/E10_smoke
    bash experiments/rebuttal/E10_radius_across_valid_sets/launch.sh

Smoke: four released rows (exhaustive, sampled, ladder, large-network paths) and the 3-node census plus one 4-node CPDAG, tiny budgets, about a minute. Full run: estimated 10 to 30 h in one process (the per-set caps bound the worst case); `--part networks|census` and `--shard I --n-shards K` split it.
Outputs in `results/rebuttal/E10/`: `rows_networks.jsonl` and `rows_census.jsonl` (one line per unit, resumable), `config.json`, `run.log`, `DONE`, and `results.json` with the registered fraction, its Clopper-Pearson and network-clustered bootstrap intervals, spread, censoring and the census counts.
