# E11: survival ranking on queries chosen from the CPDAG alone

Measures the within-state Kendall tau_b between the breakdown radius (GAC) and the survival AUC on the
frozen 648-row CPDAG-only frame (results/frame/frame.jsonl), eight real-name elicited conditions, with the
submitted survival protocol; compares it with the corrected 0.43 [0.07, 0.69] of Section 6.3. The Table 1
frame is rescored with the same code and both criteria as a reference. Registered plan: PREREG.md.

    .venv/bin/python experiments/rebuttal/E11_survival_on_cpdag_only_queries/run.py --smoke && \
    .venv/bin/python experiments/rebuttal/E11_survival_on_cpdag_only_queries/analyse.py --out results/rebuttal/E11_smoke --n-boot 200 --n-perm 50
    bash experiments/rebuttal/E11_survival_on_cpdag_only_queries/launch.sh 8

Duration: smoke under a minute; full run estimated 2 to 12 h on 8 workers (the committed Table 1 panel
took 2.1 CPU-hours; this run scores 648 + 543 pairs with two radii each). Split with --shard I --n-shards K.
Outputs in results/rebuttal/E11/: units.*.jsonl (one line per unit), config.*.json, consistency_failures.*.jsonl,
results.json (registered statistic, secondary analyses, reproduction of the committed panel), run.log, DONE.
