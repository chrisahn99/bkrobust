# E12: radius and audit on the five networks above 500 nodes

Measures whether the breakdown radius can be computed on link, munin, munin2, munin3 and munin4 (excluded by the released 500-node cap), how often it is censored, how long each call takes, and what the audit at the true DAG shows. 100 frame rows x 14 released conditions (about 1,600 units). Registered statistics in PREREG.md: N_net and the censored fraction C.

    python experiments/rebuttal/E12_networks_above_500_nodes/run.py --smoke && python experiments/rebuttal/E12_networks_above_500_nodes/analyse.py --out results/rebuttal/E12/smoke
    bash experiments/rebuttal/E12_networks_above_500_nodes/launch.sh 6

Needs `experiments/rebuttal/_data/knowledge_states.json` (the released file). Smoke: 2 rows of alarm and link, tiny budgets, one core, under a minute (not verified). Full run: 4 to 8 h on the workstation with 6 workers, bounded by the per-unit wall clocks (search 120 s, hybrid 60 s, claims 120 s).

Outputs in `results/rebuttal/E12/`: `rows.jsonl` (one line per unit, resumable), `config.json`, `oracle_knowledge.<net>.json`, `run.log`, `DONE`, and `results.json` from analyse.py. Exit code 3 means a consistency failure.
