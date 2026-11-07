# E17: worked applied case, and validity at truth by radius

Serves qQYf-S2. (1) On Acid_1996 / Polzer_2012, one D_LLM knowledge state with a radius-1 query (named claim to confirm) and an unreachable or large-radius query, recomputed with the submitted code and matched to the released rows. (2) Share of committed sets GAC-invalid at the true DAG by GAC radius bin (1, 2, >=3, unreached), network-clustered, D_LLM against the matched-|K| random states of results/matched_control; registered statistic and prediction in PREREG.md.

Smoke (two small networks, about 30 s):
    .venv/bin/python experiments/rebuttal/E17_worked_case_and_validity_by_radius/run.py --smoke --out results/rebuttal/E17_smoke && .venv/bin/python experiments/rebuttal/E17_worked_case_and_validity_by_radius/analyse.py --out results/rebuttal/E17_smoke

Full (the workstation, 8 workers, expected 20-50 min, dominated by the six largest networks timing out at 120 s per unit):
    bash experiments/rebuttal/E17_worked_case_and_validity_by_radius/launch.sh 8

Needs experiments/rebuttal/_data/per_query_records.csv (the released file). Outputs in results/rebuttal/E17/: rows.shard*.jsonl (one row per unit), config.shard*.json, results.json, analyse.txt, run.log, DONE. run.py exits 2 on any failed consistency check.
