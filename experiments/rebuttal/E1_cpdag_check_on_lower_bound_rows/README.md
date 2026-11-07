# E1: the validity check at the CPDAG on the 372 released lower-bound radii

Measures how many of the 372 rows of `release/per_query_records.csv` reported only as lower bounds (`hop_method = local_up_fast_budget`) are settled as infinite by checking the committed set at the CPDAG (Algorithm 1, line 2), under the GAC (registered, PREREG.md) and the back-door criterion; checks consistency (exit 1 on a contradiction); replays the released search, and reruns unsettled rows with ten times the closure budget.

Data: copy `per_query_records.csv` and `knowledge_states.json` from the paper's `release/` into `experiments/rebuttal/_data/` (or pass `--records`, `--knowledge`). Networks come from `results/axisa3/networks/example_models` in this repository.

    .venv/bin/python experiments/rebuttal/E1_cpdag_check_on_lower_bound_rows/run.py --smoke && \
      .venv/bin/python experiments/rebuttal/E1_cpdag_check_on_lower_bound_rows/analyse.py --out results/rebuttal/E1_smoke
    bash experiments/rebuttal/E1_cpdag_check_on_lower_bound_rows/launch.sh 8

Smoke: 4 rows of magic-irri, tiny budgets, under 60 s. Full run: 1 to 3 h on 8 cores if most rows settle; at most about 8 h if none do (each unsettled search capped at 300 s).
Outputs in `results/rebuttal/E1/`: `rows.jsonl` (one line per row, resumable), `config.json`, `run.log`, `DONE`, and `results.json` with the registered statistic, its Clopper-Pearson interval and the per-stratum table. Array jobs: `--shard I --n-shards K`.
