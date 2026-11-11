# E7 — one hidden variable, near or far from the query

Measures how often hiding one true variable before PC (n = 5000, seeds 0-4, four Gaussian networks of the learned-CPDAG panel) makes the committed adjustment set invalid in the true DAG, paired with the run where nothing is hidden, split by the hidden variable's distance to the query (hop 1 vs hop 2+). Registered statistic and prediction: PREREG.md.

    .venv/bin/python experiments/rebuttal/E7_one_hidden_variable/run.py --smoke && .venv/bin/python experiments/rebuttal/E7_one_hidden_variable/analyse.py --out results/rebuttal/E7_smoke
    bash experiments/rebuttal/E7_one_hidden_variable/launch.sh 6

Inputs: experiments/rebuttal/_data/pest/{frame.jsonl,structure.csv,graphs.jsonl} (from the panel's results/pest), and panel_samples/, the panel's own samples for the two cells that do not regenerate bit for bit on the workstation (see DEVIATIONS.md; `run.py --check-samples` checks every cell in a few seconds). Smoke: under a minute on one core. Full: about 3 h on 6 workers (arth150's PC runs dominate, about 2.5 min each).
Outputs in results/rebuttal/E7: config.json, units_*.jsonl (one line per PC run with all its query rows), results.json, run.log, DONE. Exit 3 means a consistency check of PREREG.md failed.
