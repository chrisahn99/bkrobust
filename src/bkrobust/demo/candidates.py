"""Hand-designed candidate scenarios, screened against the design gates.

Six candidates, all variations on a cardiovascular-prevention theme: does statin
therapy reduce cardiovascular events, given observational data on a handful of
clinical covariates? They differ in where the undirected component sits relative
to the treatment, which is what determines whether a perturbation can do damage.

The candidates are hand-written rather than sampled because the witness graphs
have to be *contestable*: the report asks a domain expert to look at a failure
and judge whether the orientation claim behind it is plausible, and that is only
possible if the variables mean something. Random DAGs pass every structural gate
and fail that one.

The screening results (which passed, which were rejected and why) are reported
in ``report.md``.
"""

from __future__ import annotations

CANDIDATES: list[dict[str, object]] = [
    {
        "name": "C1_block_far",
        "story": (
            "Undirected component sits among upstream lifestyle variables, far "
            "from the treatment. Included as a deliberate negative control: if "
            "the component cannot touch a back-door path, nothing can break."
        ),
        "x": "Statin",
        "y": "CVD",
        "edges": [
            ("SES", "Smoke"),
            ("SES", "Diet"),
            ("Smoke", "Diet"),
            ("Smoke", "BMI"),
            ("Diet", "BMI"),
            ("BMI", "Chol"),
            ("Chol", "Statin"),
            ("Statin", "BP"),
            ("Age", "BP"),
            ("Statin", "CVD"),
            ("BP", "CVD"),
            ("Age", "CVD"),
        ],
    },
    {
        "name": "C2_treatment_in_block",
        "story": (
            "Statin therapy and cardiovascular events, with the treatment itself "
            "inside the undirected component. A single edge orientation decides "
            "whether inflammation (CRP) is a pre-treatment confounder or a "
            "post-treatment mediator -- the canonical applied error."
        ),
        "x": "Statin",
        "y": "CVD",
        "edges": [
            ("Smoke", "BMI"),
            ("BMI", "Chol"),
            ("BMI", "CRP"),
            ("Chol", "CRP"),
            ("CRP", "Statin"),
            ("Chol", "Statin"),
            ("Statin", "BP"),
            ("Age", "BP"),
            ("Statin", "CVD"),
            ("BP", "CVD"),
            ("Age", "CVD"),
            ("Smoke", "CVD"),
        ],
    },
    {
        "name": "C3_clique4",
        "story": (
            "As C2 but the block is a 4-clique on {Smoke, BMI, Chol, Statin}, "
            "giving the largest interesting space (a clique's Meek-closed "
            "orientations are exactly its linear orders)."
        ),
        "x": "Statin",
        "y": "CVD",
        "edges": [
            ("Smoke", "BMI"),
            ("Smoke", "Chol"),
            ("BMI", "Chol"),
            ("Smoke", "Statin"),
            ("BMI", "Statin"),
            ("Chol", "Statin"),
            ("Statin", "BP"),
            ("Age", "BP"),
            ("Statin", "CVD"),
            ("BP", "CVD"),
            ("Age", "CVD"),
        ],
    },
    {
        "name": "C4_crp_pendant",
        "story": (
            "As C2 but CRP hangs off the block as a pendant rather than sitting "
            "in a triangle, shrinking the component to test the lower end of the "
            "4-6 edge window."
        ),
        "x": "Statin",
        "y": "CVD",
        "edges": [
            ("Smoke", "BMI"),
            ("BMI", "Chol"),
            ("Chol", "CRP"),
            ("CRP", "Statin"),
            ("Chol", "Statin"),
            ("Statin", "BP"),
            ("Age", "BP"),
            ("Statin", "CVD"),
            ("BP", "CVD"),
            ("Age", "CVD"),
            ("Smoke", "CVD"),
        ],
    },
    {
        "name": "C5_two_confounder_routes",
        "story": (
            "Two distinct back-door routes from Statin to CVD (via Chol and via "
            "Smoke), so several genuinely different valid adjustment sets exist "
            "and the robustness-vs-efficiency comparison has something to say."
        ),
        "x": "Statin",
        "y": "CVD",
        "edges": [
            ("Smoke", "BMI"),
            ("BMI", "Chol"),
            ("Smoke", "Chol"),
            ("BMI", "CRP"),
            ("Chol", "CRP"),
            ("CRP", "Statin"),
            ("Chol", "Statin"),
            ("Statin", "BP"),
            ("Age", "BP"),
            ("Statin", "CVD"),
            ("BP", "CVD"),
            ("Age", "CVD"),
            ("Smoke", "CVD"),
            ("Chol", "CVD"),
        ],
    },
    {
        "name": "C6_mediator_and_confounder",
        "story": (
            "As C2 with an explicit measured mediator (BP) and a second "
            "pre-treatment marker, so the adjustment-set lattice contains both "
            "over-adjustment and under-adjustment failures."
        ),
        "x": "Statin",
        "y": "CVD",
        "edges": [
            ("Smoke", "BMI"),
            ("BMI", "Chol"),
            ("BMI", "CRP"),
            ("Chol", "CRP"),
            ("CRP", "Statin"),
            ("Chol", "Statin"),
            ("Statin", "BP"),
            ("Age", "BP"),
            ("Age", "Chol"),
            ("Statin", "CVD"),
            ("BP", "CVD"),
            ("Age", "CVD"),
            ("Smoke", "CVD"),
        ],
    },
]
