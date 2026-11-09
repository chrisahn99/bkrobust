# E29: the radius against the undirected degree of X where separation is undefined

On E16's 1,033 units where the separation proxy is undefined (submitted Table 1 cohort), compares the GAC
radius with the number of undirected edges at X in the CPDAG and with the size of X's chain component, and
records for each radius-1 unit which of the analyst's orientations its nearest failing states retract and
why the committed set fails there. Registered statistics: Kendall tau_b(radius, undirected degree of X), and
the share of radius-1 units whose only orientation at X is the retracted one. Secondary: the tie ceiling of
E16's tau_b, the failure mechanism, and the within-state paired lead of the radius over the degree as rankers
of survival. Registered plan: PREREG.md; gate: FIDELITY.yaml.

    PYTHONPATH=src .venv/bin/python experiments/rebuttal/E29_radius_vs_degree_sep_undefined/run.py --smoke && \
    PYTHONPATH=src .venv/bin/python experiments/rebuttal/E29_radius_vs_degree_sep_undefined/analyse.py --out results/rebuttal/E29_smoke
    bash experiments/rebuttal/E29_radius_vs_degree_sep_undefined/launch.sh 2 <cpu-list>

Duration: smoke under a minute (4 of 94 states); full run expected a few minutes on 2 workers (the smoke scored 28 units in 0.3 s; the radius
is computed twice per unit, once with the witness and once by the submitted search, under a 600 s wall clock
and a 1.5 GB data cap per worker). Outputs in results/rebuttal/E29/: units.*.jsonl (one line per unit),
config.*.json, consistency_failures.*.jsonl, errors.*.jsonl, results.json, analyse.txt, run.log, DONE.
