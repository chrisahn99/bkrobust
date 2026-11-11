# E26: an edge-addition tolerance in the sense of Taeb et al., Section 6 (criterion gac)

For each of the 3,585 released queries with a committed set, the least number of directed edges that, added to one DAG of the analyst's class under a fixed causal order, make the committed set fail the GAC (exhaustive to depth 3), compared with the released breakdown radius: Kendall tau_b with a network-clustered bootstrap interval, the share of radius-ordered pairs the tolerance does not order alike, and the share of queries at tolerance 1 (registered statistic and prediction in PREREG.md). It answers h8pC-W3: the novelty sentence against Taeb, Guo and Henckel (arXiv:2511.10625v3) needs a number. Exit 1 on a consistency failure (rebuild, D0 validity, witness not confirmed by the submitted and released GAC oracles).

    .venv/bin/python experiments/rebuttal/E26_taeb_edge_addition_tolerance/run.py --smoke && \
      .venv/bin/python experiments/rebuttal/E26_taeb_edge_addition_tolerance/analyse.py --out results/rebuttal/E26_smoke
    # Slurm cluster (registered launch), from the repository root on the cluster:
    sbatch --job-name=rb_E26 --array=0-7 --export=ALL,EXP=E26,NSHARDS=8,EXTRA="" experiments/rebuttal/_slurm/job_array.slurm
    # or the workstation, 10 workers:
    bash experiments/rebuttal/E26_taeb_edge_addition_tolerance/launch.sh 10

Smoke: 12 rows from six networks (4 to 413 nodes), short budgets, about 5 s. Full run: 3,585 rows; most queries are decided at depth 1 in well under a second, the cost is the queries above 1 on large networks (120 s + 60 s cap each).
Outputs in `results/rebuttal/E26/`: `rows[.shardIofK].jsonl` (one line per query, resumable), `config*.json`, and from `analyse.py` `results.json` and `analyse.txt` with the registered statistic, the verdict against the dies-if bar, the cross-tabulation, extension sensitivity and the AUC secondary.
`vendor_ledger.py` is E10's copy, unchanged (rebuild of the released knowledge and committed sets, released truth oracle).
