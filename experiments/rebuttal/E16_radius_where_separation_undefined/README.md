# E16: radius ranking where separation is undefined

Measures Kendall tau_b of the GAC breakdown radius with survival AUC on the 1,033 pooled real-name Table 1
units where the separation proxy is undefined (network-cluster bootstrap, 10,000 resamples), on the submitted
and corrected cohorts, plus the within-state paired lead over separation where it is defined. Pure re-analysis
of committed per-unit tables; no solver runs. It also re-derives published Table 1 points and exits 2 on any mismatch.

    PYTHONPATH=src .venv/bin/python experiments/rebuttal/E16_radius_where_separation_undefined/run.py --smoke
    bash experiments/rebuttal/E16_radius_where_separation_undefined/launch.sh 4

Smoke: under a minute, writes results/rebuttal/E16_smoke/ (analyse with `analyse.py --out results/rebuttal/E16_smoke`).
Full: a few minutes on 4 cores. Outputs in results/rebuttal/E16/: rows.jsonl (one line per cohort x statistic),
config.json, results.json (registered statistic, verdict, all cells), run.log, DONE. Plan: PREREG.md; gate: FIDELITY.yaml.
