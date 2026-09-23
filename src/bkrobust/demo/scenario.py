"""The chosen running example and the three analyst knowledge states.

Selected as candidate ``C4_crp_pendant`` after screening six hand-designed
candidates against the five design gates in :mod:`bkrobust.demo.example`; four
were rejected. See ``report.md`` for the screening table.

The scenario
------------
An analyst asks whether statin therapy (``Statin``) reduces cardiovascular
events (``CVD``), from observational data on eight variables. The estimated
CPDAG leaves a five-edge undirected component on
``{Smoke, BMI, Chol, CRP, Statin}`` -- the data cannot orient those edges -- and
six edges are compelled by v-structures and therefore frozen.

Why this example
----------------
The treatment sits *inside* the undirected component, so a single edge
orientation decides whether C-reactive protein is a pre-treatment confounder or
a post-treatment mediator. That is a live clinical question -- statins do lower
CRP -- and getting it backwards converts a legitimate covariate into one that
must never be adjusted for. The witness graphs are therefore contestable by a
domain expert, which is what the demonstration needs.
"""

from __future__ import annotations

from bkrobust.demo.graph import MPDAG

Edge = tuple[str, str]

#: Ground-truth DAG. Directed edges only.
TRUE_EDGES: list[Edge] = [
    ("Smoke", "BMI"),  # smoking raises BMI
    ("BMI", "Chol"),  # BMI raises cholesterol
    ("Chol", "CRP"),  # cholesterol drives inflammation
    ("CRP", "Statin"),  # inflammation drives prescribing
    ("Chol", "Statin"),  # cholesterol drives prescribing
    ("Statin", "BP"),  # statin lowers blood pressure
    ("Age", "BP"),  # age raises blood pressure
    ("Statin", "CVD"),  # the effect of interest
    ("BP", "CVD"),  # blood pressure raises event risk
    ("Age", "CVD"),  # age raises event risk
    ("Smoke", "CVD"),  # smoking raises event risk
]

TREATMENT = "Statin"
OUTCOME = "CVD"

#: The four orientation claims a clinician would volunteer. All true. Only two
#: (Smoke->BMI and CRP->Statin) are strictly needed -- Meek forces the rest,
#: which is the cascade the report points at.
K_TRUE: list[Edge] = [
    ("Smoke", "BMI"),
    ("BMI", "Chol"),
    ("Chol", "CRP"),
    ("CRP", "Statin"),
]

#: The load-bearing claim. Scenarios A and B both target it -- A withholds it,
#: B asserts its opposite -- so the two are directly comparable.
PIVOT: Edge = ("CRP", "Statin")


def true_dag() -> MPDAG:
    """The ground-truth DAG."""
    nodes = sorted({n for e in TRUE_EDGES for n in e})
    return MPDAG(nodes, directed=TRUE_EDGES)


def scenarios() -> dict[str, dict[str, object]]:
    """The three analyst knowledge states, all imposed on the same CPDAG.

    Returns:
        Mapping from scenario label to a dict with ``label``, ``description``,
        ``knowledge`` (the asserted orientations) and ``relation`` (how it
        differs from ``K_TRUE``, in words).
    """
    a = [e for e in K_TRUE if e != PIVOT]
    b = [e for e in K_TRUE if e != PIVOT] + [(PIVOT[1], PIVOT[0])]
    # C omits Smoke->BMI, not BMI->Chol: omitting BMI->Chol is a no-op because
    # Meek re-derives it, which would make C identical to B. Verified in
    # tests/demo/test_scenario.py::test_scenarios_are_distinct.
    c = [e for e in K_TRUE if e not in (PIVOT, ("Smoke", "BMI"))] + [(PIVOT[1], PIVOT[0])]
    return {
        "A": {
            "label": "A - mild (one assertion omitted)",
            "description": (
                "The analyst declines to say whether inflammation precedes or "
                "follows statin therapy, and asserts the other three claims. "
                "Knows less than the truth, but says nothing false."
            ),
            "knowledge": a,
            "relation": "K_true minus {CRP->Statin}",
        },
        "B": {
            "label": "B - wrong (one assertion flipped)",
            "description": (
                "The analyst believes statins lower CRP, so asserts "
                "Statin->CRP. Consistent with the CPDAG, and false. This "
                "converts CRP from a confounder into a mediator."
            ),
            "knowledge": b,
            "relation": "K_true with CRP->Statin flipped to Statin->CRP",
        },
        "C": {
            "label": "C - compound (one flip plus one omission)",
            "description": (
                "As B, and the analyst also declines to commit on whether "
                "smoking raises BMI. The realistic messy case: one false belief "
                "and one gap. Strictly B plus ambiguity."
            ),
            "knowledge": c,
            "relation": "K_true with CRP->Statin flipped and Smoke->BMI omitted",
        },
    }
