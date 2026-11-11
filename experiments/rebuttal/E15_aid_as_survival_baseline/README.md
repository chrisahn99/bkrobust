# E15: adjustment identification distance (AID) as a survival baseline

Measures whether the AID of Henckel et al. (2024), from the analyst's graph G0 to the CPDAG (analyst-side)
and to the true DAG (oracle), ranks survival (stored AUC_frac) as well as the GAC breakdown radius, on the
480 untargeted and 4,104 targeted stored instances, each rebuilt from its seed and checked against the
stored x, Z*, radius and separation. Registered statistic and prediction: PREREG.md. Our AID uses the
submitted GAC predicate (aid.py); gadjid's oset/parent/ancestor AID is computed beside it.

    PYTHONPATH=src .venv/bin/python experiments/rebuttal/E15_aid_as_survival_baseline/run.py --smoke && \
      PYTHONPATH=src .venv/bin/python experiments/rebuttal/E15_aid_as_survival_baseline/analyse.py --out results/rebuttal/E15_smoke --n-boot 200
    WORKERS=8 experiments/rebuttal/E15_aid_as_survival_baseline/launch.sh   # full run + analysis

Expected: under one hour on 8 cores (the earlier GAC recompute of the same instances took ~4 min serially; AID adds ~182 pairs per graph).
Outputs in results/rebuttal/E15/: units_shard*.jsonl (one row per instance), config_shard*.json, results.json, E15.log, E15.done.
Install gadjid first (pip install gadjid); a full run refuses to start without it unless --allow-no-gadjid.
