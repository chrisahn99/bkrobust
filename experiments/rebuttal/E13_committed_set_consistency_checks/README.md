# E13: consistency checks on committed sets and the error budget

Measures (a) whether the 898 real-name units with "no optimal set" (status_split.json) are amenable
at G0; (b) on all 8,308 released rows, whether the closed-form committed set equals the optimal set
read off every DAG extension of G0; (c) whether every row invalid at the truth has the truth reversing
at least r_val closure orientations of G0 (Corollary 1(c)), plus r_val recomputed by Corollary 1(c)
where Ext(CPDAG) has at most 256 DAGs. GAC throughout. Registered statistics: PREREG.md.

    PYTHONPATH=src .venv/bin/python experiments/rebuttal/E13_committed_set_consistency_checks/run.py --smoke
    bash experiments/rebuttal/E13_committed_set_consistency_checks/launch.sh

Needs experiments/rebuttal/_data/per_query_records.csv (the release file). Expected duration: minutes
to a few hours on 6 workers, dominated by G0s with thousands of extensions on the largest networks.
Outputs in results/rebuttal/E13/: rows.jsonl (one line per unit), config.json, results.json
(from analyse.py), E13.log, E13.done. Exit 1 = a consistency failure (see HARD lines in the log).
