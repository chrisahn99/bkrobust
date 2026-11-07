# E2: running time of the submitted hybrid under the GAC, n = 8 to 400

Times each stage of `bkrobust.hybrid.breakdown_radius(criterion="gac", search_budget=3)` in a fresh process: CPDAG build, top-state check, one validity call, depth-3 search, SAT build, SAT ladder alone, end-to-end hybrid (plus the unbounded search as a cross-check when the SAT stage answers). Two families: the Appendix D grid (n 8-24, 256 units) and sparse Erdos-Renyi with expected degree 2 or 4 (n 16-400, 160 units). Cap 120 s and 4 GB per leg; censored legs are kept. Protocol and registered statistic: PREREG.md. Fidelity: FIDELITY.yaml.

Run from the repository root:

    .venv/bin/python experiments/rebuttal/E2_running_time_by_graph_size/run.py --smoke      # < 60 s, results/rebuttal/E2_smoke
    bash experiments/rebuttal/E2_running_time_by_graph_size/launch.sh            # full run, serial, nohup
    .venv/bin/python experiments/rebuttal/E2_running_time_by_graph_size/analyse.py --out results/rebuttal/E2

Expected duration: 6 to 10 h serial (most units take seconds; the large sparse units hit the 120 s caps until the per-series skip rule stops them). Split with `--shard I --n-shards K` if needed.

Outputs in `results/rebuttal/E2/`: `rows_shard*.jsonl` (one line per unit, with criterion and commit), `run_config_shard*.json`, `run.log`, `analyse.txt`, `results.json`, `DONE`. No input data: instances are generated from seeds.
