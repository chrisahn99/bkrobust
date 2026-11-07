# E19: the App. E.4 bias-radius study under internal standardisation

Measures whether the App. E.4 shares (threshold never reached on 70.1% of 900 instances at epsilon = 100%; step staircase on 84.1%) hold when the same coefficient draw is internally standardised (iSCM). Each stored instance runs three arms: naive (must reproduce the stored row), post hoc standardised (must equal naive), internal. Criterion gac. Prereg: PREREG.md; rubric: FIDELITY.yaml (checks in fidelity_checks.py).

Smoke (6 fast instances, under a minute, from the repository root):
    .venv/bin/python experiments/rebuttal/E19_bias_radius_internal_standardisation/run.py --smoke && .venv/bin/python experiments/rebuttal/E19_bias_radius_internal_standardisation/analyse.py --smoke

Full run (the workstation, nohup, log, done-file; split with --shard I --n-shards K if wanted):
    bash experiments/rebuttal/E19_bias_radius_internal_standardisation/launch.sh 6

Expected duration: not measured. The stored rows record 17 s in total for the 900 instances (max 1.4 s each), so 3 arms should take minutes, unless the internal arm traverses much further (cap: --timeout 60 s per arm).

Outputs, in results/rebuttal/E19/: units.shard*of*.jsonl (one line per instance, all arms), run_config.shard*.json, run.log, analyse.txt, results.json, DONE.
Exit 3 from run.py means a consistency check failed (naive differs from the stored study, post hoc differs from naive, strategies disagree, variance not 1, r_val differs across arms).
