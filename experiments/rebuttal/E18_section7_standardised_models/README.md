# E18: Section 7 profiles under standardised models, ten coefficient draws

Measures whether the mean-bias profiles r_mu and ceilings B(Chat)/|theta_Z| of Tables 2 and 6 change
when the linear-Gaussian model is standardised post hoc or internally (iSCM), over ten coefficient
draws, on the 49 queries of Table 6 (reviewer C2Bd-Q1). Plan and decision rules: PREREG.md.

    PYTHONPATH=src .venv/bin/python experiments/rebuttal/E18_section7_standardised_models/run.py --smoke && PYTHONPATH=src .venv/bin/python experiments/rebuttal/E18_section7_standardised_models/analyse.py --smoke
    bash experiments/rebuttal/E18_section7_standardised_models/launch.sh 6

Smoke: 5 queries x 2 draws x 3 arms, under a minute. Full: 1,470 units; expected well under an
hour on 6 workers (dominated by the one-off slice enumeration of child CO2->DuctFlow, ~5 min, and
the four Kampen_2014 queries, ~1 min each).
Outputs in results/rebuttal/E18/: units.s*.jsonl (one row per query, draw, arm), config.s*.json,
failures.s*.jsonl (must be empty), results.json and the table printed by analyse.py. run.py exits 2
on any consistency failure (posthoc differs from naive, draw 0 does not reproduce the stored rows).
