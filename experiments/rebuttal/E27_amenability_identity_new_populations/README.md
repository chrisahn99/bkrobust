# E27: the amenability identity on populations the first author did not score

Measures, per query, whether the radius of the committed set (r_val, the submitted `breakdown_radius`,
criterion gac) equals the amenability radius (r_amen, the fewest retractions of the analyst's orientations
after which (X, Y) is no longer amenable; the first author's `e10b_amenability.py` definition), and the time
each takes. Population A: loadable networks outside the paper's 33 under the CPDAG-only rule (census; all
seven such files have a fully oriented CPDAG or are not DAGs, so it is expected to be empty). Population B:
the frozen 648-pair CPDAG-only frame x the ten elicited knowledge states, states built as in E11. Both radii
are checked by plain enumeration when G0 commits at most 12 orientations. Registered plan: PREREG.md.

    .venv/bin/python experiments/rebuttal/E27_amenability_identity_new_populations/run.py --smoke && \
    .venv/bin/python experiments/rebuttal/E27_amenability_identity_new_populations/analyse.py --out results/rebuttal/E27_smoke
    bash experiments/rebuttal/E27_amenability_identity_new_populations/launch.sh 6

Duration: smoke under a minute; full run expected under 3 h on 6 workers (E11's median r_val call was 2 ms;
r_amen and the enumeration are the unknown part; every call has a 600 s wall clock and a 1.5 GB data cap per worker).
Outputs in results/rebuttal/E27/: units.*.jsonl (one line per unit), census_population_A.json, config.*.json,
consistency_failures.*.jsonl, errors.*.jsonl, results.json, analyse.txt, run.log, DONE.

Not covered here (card V2): the 388 E12 units on the five networks above 500 nodes with the synthetic arms,
and the E5 learned graphs. The E11 frame already contains link, munin and munin1-4 rows with elicited knowledge.
