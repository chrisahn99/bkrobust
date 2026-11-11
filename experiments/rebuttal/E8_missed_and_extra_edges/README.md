# E8: missed and extra edges, how local is the damage

Measures, on X3's synthetic population with fresh seeds (32,800 SCM draws, population covariance), (1) the
silent-bias rate of the committed adjustment set O* when discovery misses one edge two or more hops from the
query (registered bar 0.002), and (2) whether one extra edge, oriented as the true order allows, ever makes O*
invalid in the true DAG (registered: zero failures). GAC validity, two implementations. See `PREREG.md`, `FIDELITY.yaml`.

    .venv/bin/python experiments/rebuttal/E8_missed_and_extra_edges/run.py --smoke && \
    .venv/bin/python experiments/rebuttal/E8_missed_and_extra_edges/analyse.py --smoke
    bash experiments/rebuttal/E8_missed_and_extra_edges/launch.sh 8   # then analyse.py

Expected duration: smoke under a minute; full run a few hours on 8 cores (not measured; the additions,
exhaustive on p up to 25, dominate). Split with `--shard I --n-shards K` for array jobs.
Outputs in `results/rebuttal/E8/`: `E8_<shard>.jsonl` (one line per SCM), `config_<shard>.json`, the log,
`E8.done` or `E8.failed`, and `results.json` from `analyse.py`.
