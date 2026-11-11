# E9: parametric bootstrap of validity under a fitted model

Measures whether simulating from a linear-Gaussian model fitted to one DAG of the analyst's state (200 datasets per network, PC rerun, same knowledge re-imposed, committed set recomputed and checked with the GAC in the simulated truth) flags the 244 certifiable rows of the learned-CPDAG panel whose Z is invalid in the true DAG better than E6's row bootstrap. Registered statistic: AUC, and Delta = AUC(E9) - AUC(E6) with a paired network-then-pair cluster bootstrap interval (PREREG.md).

Needs `experiments/rebuttal/_data/pest/{chat_e.json,knowledge.json,rows.csv}` (same files as E6) and causal-learn 0.1.4.8. Runs only after E6's full run (the launcher checks).

    .venv/bin/python experiments/rebuttal/E9_parametric_bootstrap_validity/run.py --smoke && .venv/bin/python experiments/rebuttal/E9_parametric_bootstrap_validity/analyse.py --out results/rebuttal/E9_smoke --e6 results/rebuttal/E6_smoke
    bash experiments/rebuttal/E9_parametric_bootstrap_validity/launch.sh 6

Smoke: ecoli70 and magic-niab, 2 simulations plus identity, 2 extra extensions, about a minute. Full: 804 PC runs, about twice E6 (roughly 4 to 6 h on 6 workers; arth150 dominates).
Outputs in `results/rebuttal/E9/`: `config.json`, `reference.jsonl`, `models.json` (fitted truths, covariance gaps), `units_*.jsonl` (one line per network x simulation), `results.json`, `per_row.json`, `run.log`, `DONE`.
Since 2026-11-06 (DEVIATIONS.md) the 38 arth150 rows are rebuilt and checked but excluded from the simulation, because the arth150 reference CPDAG has a chordless undirected 4-cycle and so no DAG extension; `exclusions.json` lists them, and `run.py --check-reference --out <dir>` runs the rebuild alone (about 2 s).
Exit 3 = rebuilt rows disagree with the frozen panel; exit 4 = a model invariant broke (no extension, equivalent fits differ, Z invalid in a fitted truth, identity unit not the reference). Split with `--shard I --n-shards K`.
