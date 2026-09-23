r"""Baseline comparison (naive Hamming distance) and calibration for the certificate.

Two independent pieces, per the brief:

Part 1 -- naive assertion counting.
    The model-oriented distance used throughout the rest of the demonstration is
    BFS hop count on the covering-relation neighbour graph
    (:mod:`bkrobust.demo.space`). The obvious alternative an analyst without that
    machinery would reach for is much simpler: Hamming distance between
    knowledge *sets* -- how many orientation claims would need to be added,
    dropped or flipped to get from one state of background knowledge to
    another. :func:`naive_distance` is that count; :func:`naive_vs_model_frame`
    lines the two metrics up side by side over the enumerated space so they can
    be compared directly, and :func:`inconsistent_fraction` answers a related
    question the brief asks for explicitly: of all the single-assertion edits an
    analyst might make to their knowledge, what fraction produce something
    outright inconsistent with the data (Meek closure FAILs) rather than a
    merely-different MPDAG?

    IMPORTANT FRAMING, stated plainly and not editorialised in either
    direction: on the three scenarios this module was actually run against
    (48-element space each), the model radius (``r_val``) and the naive radius
    are *close but never identical* -- ``model_radius - naive_radius == 1`` in
    every one of scenarios A, B and C (3 vs 2, 3 vs 2, 2 vs 1). Rank
    correlation (Spearman) between the two distance columns is 0.73-0.77 and
    Pearson is 0.75-0.79 across the three -- a real, strong association, but
    with substantial scatter: about 14 of 48 elements (~29%) in every scenario
    differ by 2 or more model-hops from their naive-Hamming count, some by as
    much as 4 (see :func:`disagreement_examples`, and the report this module
    feeds, for concrete cases -- in every disagreement observed here, the
    naive count *underestimates* the model distance, because Meek propagation
    chains a single re-assertion into several forced covering-relation steps
    that a raw Hamming count on K cannot see). So a naive analyst using
    Hamming distance on this example would not have been catastrophically
    misled -- the two rankings broadly track each other, and the naive radius
    happens to be off by exactly one step in every scenario tried -- but they
    also would not have gotten the model radius exactly right, and the gap is
    systematic in one direction (naive is optimistic: it reports the nearest
    failure as closer than it model-theoretically is), not noise. Whether that
    holds on other examples is untested; the value of the model-oriented
    distance is that it is *derived* from the actual model-inclusion structure
    rather than an incidental property of the assertion count, so it does not
    depend on this near-agreement continuing to hold elsewhere.

Part 2 -- calibration.
    The running example has a known ground truth, so the certificate
    (``r_val`` covering a given adjustment set ``Z``) can actually be checked
    rather than merely trusted. :func:`calibration_table` treats every element
    of the enumerated space in turn as a hypothetical truth and records whether
    it falls inside the certified radius and whether ``Z`` is genuinely valid
    there; :func:`calibration_summary` reduces that to the two numbers the brief
    asks for -- coverage (should be 1.0 by construction; a bug if not) and
    conservativeness (how often the certificate is stricter than it needed to
    be); :func:`calibration_for_actual_truth` reports the single real data
    point, the actual ground-truth DAG, handling explicitly the case where that
    DAG is not itself an element of the enumerated space.

Written for Python 3.9: ``from __future__ import annotations`` defers every
annotation to a string, so builtin generics (``dict[str, int]``) and ``X | Y``
unions are used freely *in annotations only*; no such syntax appears in a
runtime position (default values, isinstance checks, etc).
"""

from __future__ import annotations

import warnings
from collections.abc import Iterable
from itertools import product

import pandas as pd

from bkrobust.demo.evaluate import is_valid_adjustment_set_dag, is_valid_adjustment_set_mpdag
from bkrobust.demo.graph import MPDAG, Edge, canon
from bkrobust.demo.meek import apply_orientations
from bkrobust.demo.pipeline import UNREACHED
from bkrobust.demo.scenario import OUTCOME, TREATMENT
from bkrobust.demo.space import atomic_moves

# --------------------------------------------------------------------------------
# Part 1: naive assertion counting
# --------------------------------------------------------------------------------


def naive_distance(k1: Iterable[Edge], k2: Iterable[Edge]) -> int:
    """Hamming distance between two knowledge sets: ``|K1 symmetric-difference K2|``.

    Each knowledge set is a collection of ``(tail, head)`` orientation
    assertions. This is a genuine metric on the space of such sets (symmetric,
    zero iff equal, satisfies the triangle inequality -- all three checked
    directly in the test suite): it is exactly the Hamming distance between the
    sets' indicator vectors over the universe of possible assertions.

    Args:
        k1: A knowledge set.
        k2: Another knowledge set.

    Returns:
        ``len(set(k1) ^ set(k2))``.
    """
    return len(set(k1) ^ set(k2))


def k_sets_for_space(cpdag: MPDAG, k_assumed: Iterable[Edge], space: list) -> dict:
    """For each element of ``space``, the cheapest knowledge-set edit that reaches it.

    For every element ``G`` of ``space``, this is the minimum
    :func:`naive_distance` from ``k_assumed`` to any knowledge set ``K'`` whose
    Meek closure (:func:`~bkrobust.demo.meek.apply_orientations` on ``cpdag``)
    equals ``G``.

    ``K'`` is drawn from the candidate assertion pool implied by ``cpdag``'s own
    undirected edges: for each of ``cpdag``'s ``k`` undirected edges, either
    orient it one way, orient it the other way, or omit it entirely from
    ``K'``. That is ``3**k`` combinations; ``k`` is at most 5 in this
    demonstration, so brute force is cheap and exact -- no candidate ``K'``
    outside this pool can reach an element of ``space`` at all, since every
    element of ``space`` is itself one of ``cpdag``'s ``3**k`` orientation
    assignments (see :func:`bkrobust.demo.space.enumerate_space`).

    Args:
        cpdag: The background-knowledge CPDAG the space is anchored at.
        k_assumed: The analyst's actual asserted knowledge (the scenario's
            ``"knowledge"``).
        space: The enumerated space, as returned by
            :func:`bkrobust.demo.space.enumerate_space`.

    Returns:
        A mapping from each element of ``space`` to the minimum Hamming
        distance found, or :data:`~bkrobust.demo.pipeline.UNREACHED` if no
        candidate ``K'`` in the pool closes to that element (should not happen
        in practice, since every space element is itself reachable by asserting
        exactly its own orientations, but checked rather than assumed).
    """
    undirected = sorted(cpdag.undirected_edges)
    k = len(undirected)
    k_assumed_set = set(k_assumed)

    best: dict = {}
    for state in product((0, 1, 2), repeat=k):
        assertion: list = []
        for (a, b), s in zip(undirected, state):  # noqa: B905
            if s == 1:
                assertion.append((a, b))
            elif s == 2:
                assertion.append((b, a))
        closure = apply_orientations(cpdag, assertion)
        if closure is None:
            continue
        dist = naive_distance(k_assumed_set, assertion)
        if closure not in best or dist < best[closure]:
            best[closure] = dist

    return {g: best.get(g, UNREACHED) for g in space}


def naive_vs_model_frame(result: dict) -> pd.DataFrame:
    """One row per space element, comparing model distance and naive distance.

    Columns (exact names, consumed by subagent S4's fig5 -- do not rename):
    ``edge_string``, ``model_distance``, ``naive_distance``, ``z_valid``,
    ``z_is_optimal``, ``mean_abs_bias``.

    Args:
        result: Output of :func:`bkrobust.demo.pipeline.run_scenario`.

    Returns:
        A DataFrame with one row per element of ``result["space"]``, in the
        same order.
    """
    cpdag = result["cpdag"]
    space = result["space"]
    rows = result["rows"]
    k_assumed = result["spec"]["knowledge"]
    k_sets = k_sets_for_space(cpdag, k_assumed, space)

    records = []
    for g, row in zip(space, rows):  # noqa: B905 - equal length by construction
        records.append(
            {
                "edge_string": g.edge_string(),
                "model_distance": row["shell"],
                "naive_distance": k_sets[g],
                "z_valid": row["z_valid"],
                "z_is_optimal": row["z_is_optimal"],
                "mean_abs_bias": row["mean_abs_bias"],
            }
        )
    return pd.DataFrame.from_records(
        records,
        columns=[
            "edge_string",
            "model_distance",
            "naive_distance",
            "z_valid",
            "z_is_optimal",
            "mean_abs_bias",
        ],
    )


def naive_radius(result: dict) -> int:
    """The naive-metric analogue of ``r_val``: min naive distance to a ``z_valid=False`` element.

    Args:
        result: Output of :func:`bkrobust.demo.pipeline.run_scenario`.

    Returns:
        The minimum ``naive_distance`` over rows with ``z_valid`` False and a
        defined (not :data:`~bkrobust.demo.pipeline.UNREACHED`) naive distance;
        :data:`~bkrobust.demo.pipeline.UNREACHED` if there is no such row.
    """
    frame = naive_vs_model_frame(result)
    invalid = frame[~frame["z_valid"]]
    reachable = invalid[invalid["naive_distance"] != UNREACHED]
    if reachable.empty:
        return UNREACHED
    return int(reachable["naive_distance"].min())


def _defined_pairs(frame: pd.DataFrame) -> pd.DataFrame:
    """Rows of ``frame`` where both distance columns are defined (not UNREACHED)."""
    return frame[(frame["model_distance"] != UNREACHED) & (frame["naive_distance"] != UNREACHED)]


def _agreement_from_frame(frame: pd.DataFrame) -> dict:
    """Correlation and disagreement stats between ``model_distance`` and ``naive_distance``.

    Split out from :func:`compare_radii` so it is directly unit-testable on a
    hand-built frame, in particular to verify that rows carrying
    :data:`~bkrobust.demo.pipeline.UNREACHED` in either distance column are
    excluded from every aggregate rather than silently treated as a real (and
    very small, since ``UNREACHED == -1``) distance.

    Args:
        frame: A frame shaped like :func:`naive_vs_model_frame`'s output (at
            minimum, ``model_distance`` and ``naive_distance`` columns).

    Returns:
        A dict with ``spearman``, ``pearson`` (``nan`` if fewer than two
        comparable rows, or if either column is constant on the comparable
        rows), ``n_compared`` (rows with both distances defined) and
        ``n_disagree_by_2_or_more``.
    """
    defined = _defined_pairs(frame)
    n_compared = len(defined)

    model_col = defined["model_distance"]
    naive_col = defined["naive_distance"]
    if n_compared >= 2 and model_col.nunique() > 1 and naive_col.nunique() > 1:
        pearson = float(model_col.corr(naive_col, method="pearson"))
        # Spearman rho without a scipy dependency: Pearson correlation of the
        # (average, for ties) ranks is exactly Spearman's rho.
        spearman = float(model_col.rank().corr(naive_col.rank(), method="pearson"))
    else:
        pearson = float("nan")
        spearman = float("nan")

    if n_compared:
        diff = (defined["model_distance"] - defined["naive_distance"]).abs()
        n_disagree = int((diff >= 2).sum())
    else:
        n_disagree = 0

    return {
        "spearman": spearman,
        "pearson": pearson,
        "n_compared": n_compared,
        "n_disagree_by_2_or_more": n_disagree,
    }


def compare_radii(result: dict) -> dict:
    """Compare the model radius and the naive radius, plus overall agreement stats.

    Args:
        result: Output of :func:`bkrobust.demo.pipeline.run_scenario`.

    Returns:
        A dict with ``model_radius`` (``result["radii"]["r_val"]``),
        ``naive_radius`` (:func:`naive_radius`), ``radii_agree`` (whether the
        two are equal, including both being
        :data:`~bkrobust.demo.pipeline.UNREACHED`), ``spearman``, ``pearson``
        and ``n_compared`` (see :func:`_agreement_from_frame`), and
        ``n_disagree_by_2_or_more`` -- the count of elements where the two
        distances differ by 2 or more, among elements where both are defined.
    """
    frame = naive_vs_model_frame(result)
    model_r = result["radii"]["r_val"]
    naive_r = naive_radius(result)
    agreement = _agreement_from_frame(frame)
    return {
        "model_radius": model_r,
        "naive_radius": naive_r,
        "radii_agree": model_r == naive_r,
        **agreement,
    }


def disagreement_examples(result: dict, n: int = 3) -> list:
    """The elements where the model distance and naive distance differ most.

    Args:
        result: Output of :func:`bkrobust.demo.pipeline.run_scenario`.
        n: How many examples to return (fewer if there are not that many
            elements with both distances defined).

    Returns:
        Up to ``n`` dicts, most-disagreeing first, each with ``edge_string``,
        ``model_distance``, ``naive_distance``, ``difference``,
        ``atomic_moves`` (the move sequence from ``G0`` to that element, from
        :func:`bkrobust.demo.space.atomic_moves`) and a plain-language
        ``description`` of the case.
    """
    cpdag = result["cpdag"]
    space = result["space"]
    rows = result["rows"]
    g0 = result["g0"]
    k_assumed = result["spec"]["knowledge"]
    k_sets = k_sets_for_space(cpdag, k_assumed, space)

    candidates = []
    for g, row in zip(space, rows):  # noqa: B905 - equal length by construction
        model_d = row["shell"]
        naive_d = k_sets[g]
        if model_d == UNREACHED or naive_d == UNREACHED:
            continue
        candidates.append((abs(model_d - naive_d), g, model_d, naive_d))

    candidates.sort(key=lambda t: (-t[0], t[1].edge_string()))

    examples = []
    for diff, g, model_d, naive_d in candidates[:n]:
        moves = atomic_moves(g0, g)
        move_text = ", ".join(moves) if moves else "no moves (identical to G0)"
        if model_d > naive_d:
            description = (
                f"{g.edge_string()!r} sits {model_d} model-hops from G0 "
                f"({move_text}), but the cheapest knowledge-set edit that reaches "
                f"it costs only naive distance {naive_d}. A single re-assertion "
                f"there looks cheap by the naive count, yet Meek propagation "
                f"chains it into {model_d} covering-relation steps -- the naive "
                f"metric understates how far this graph actually is in the "
                f"model-inclusion order."
            )
        elif naive_d > model_d:
            description = (
                f"{g.edge_string()!r} is only {model_d} model-hops from G0 via "
                f"the covering relation ({move_text}), yet no consistent "
                f"knowledge-set edit reaches it in fewer than naive distance "
                f"{naive_d} assertions. The model-oriented distance says this "
                f"graph is close; a naive Hamming count on K would say it is far."
            )
        else:
            description = f"{g.edge_string()!r}: model and naive distances agree at {model_d}."
        examples.append(
            {
                "edge_string": g.edge_string(),
                "model_distance": model_d,
                "naive_distance": naive_d,
                "difference": diff,
                "atomic_moves": moves,
                "description": description,
            }
        )
    return examples


def inconsistent_fraction(cpdag: MPDAG, k_assumed: Iterable[Edge]) -> dict:
    """Fraction of single-assertion perturbations of ``k_assumed`` that are inconsistent.

    Enumerates every single-assertion perturbation of ``k_assumed``:

    * ``remove``: drop one existing assertion.
    * ``add``: assert one currently-unasserted undirected edge of ``cpdag``, in
      either direction (two perturbations per unasserted edge).
    * ``flip``: reverse one existing assertion.

    For each, :func:`~bkrobust.demo.meek.apply_orientations` is run against
    ``cpdag``; a ``None`` result means that perturbation is inconsistent with
    ``cpdag`` (no valid MPDAG extends it) rather than merely landing on a
    different (but valid) MPDAG.

    Args:
        cpdag: The background-knowledge CPDAG.
        k_assumed: The knowledge set to perturb.

    Returns:
        A dict with ``counts`` (per perturbation type: ``total`` and
        ``inconsistent``), ``total_perturbations``, ``total_inconsistent`` and
        ``fraction_inconsistent`` (``nan`` if there are no perturbations to
        make at all).
    """
    k_list = list(dict.fromkeys(k_assumed))
    k_set = set(k_list)
    asserted_canon = {canon(a, b) for a, b in k_list}
    unasserted = sorted(e for e in cpdag.undirected_edges if e not in asserted_canon)

    perturbations: dict = {"remove": [], "add": [], "flip": []}
    for e in sorted(k_list):
        perturbations["remove"].append(sorted(k_set - {e}))
    for a, b in unasserted:
        perturbations["add"].append(sorted(k_set | {(a, b)}))
        perturbations["add"].append(sorted(k_set | {(b, a)}))
    for tail, head in sorted(k_list):
        flipped = (k_set - {(tail, head)}) | {(head, tail)}
        perturbations["flip"].append(sorted(flipped))

    counts: dict = {}
    total = 0
    total_inconsistent = 0
    for kind, plist in perturbations.items():
        n_total = len(plist)
        n_bad = sum(1 for p in plist if apply_orientations(cpdag, p) is None)
        counts[kind] = {"total": n_total, "inconsistent": n_bad}
        total += n_total
        total_inconsistent += n_bad

    return {
        "counts": counts,
        "total_perturbations": total,
        "total_inconsistent": total_inconsistent,
        "fraction_inconsistent": (total_inconsistent / total) if total else float("nan"),
    }


# --------------------------------------------------------------------------------
# Part 2: calibration
# --------------------------------------------------------------------------------


def _inside_radius(d: int, r_val: int) -> bool:
    """Whether a shell distance ``d`` counts as "inside" a certified radius ``r_val``.

    Two edge cases handled explicitly:

    * ``d == UNREACHED``: the element is not reachable from G0 at all in the
      neighbour graph, so it is never "inside" -- there is no finite distance
      to compare against the radius.
    * ``r_val == UNREACHED``: no failure of the certified property was found
      anywhere in the (finite, fully enumerated) space, i.e. the radius is
      effectively infinite. Every reachable element then counts as inside.
    """
    if d == UNREACHED:
        return False
    if r_val == UNREACHED:
        return True
    return d < r_val


def calibration_table(result: dict) -> pd.DataFrame:
    """Treat every element of the space as a hypothetical truth and check the certificate.

    For each element ``G_true`` of ``result["space"]``: ``distance`` is its
    model (shell) distance from ``G0``, ``inside_radius`` is whether that
    distance is inside the certified ``r_val`` (see :func:`_inside_radius`),
    and ``z_valid`` is whether ``result["z"]`` is actually a valid adjustment
    set in ``G_true``, checked directly via
    :func:`~bkrobust.demo.evaluate.is_valid_adjustment_set_mpdag` (i.e. valid
    in every DAG extension of ``G_true``) -- not read off the memoised
    per-element rows already in ``result``, so this is an independent check.

    Args:
        result: Output of :func:`bkrobust.demo.pipeline.run_scenario`.

    Returns:
        A DataFrame with columns ``edge_string``, ``distance``,
        ``inside_radius``, ``z_valid``, one row per element of
        ``result["space"]``.
    """
    space = result["space"]
    shells = result["shells"]
    r_val = result["radii"]["r_val"]
    z = result["z"]

    records = []
    for g in space:
        d = shells.get(g, UNREACHED)
        records.append(
            {
                "edge_string": g.edge_string(),
                "distance": d,
                "inside_radius": _inside_radius(d, r_val),
                "z_valid": is_valid_adjustment_set_mpdag(g, TREATMENT, OUTCOME, z),
            }
        )
    return pd.DataFrame.from_records(
        records, columns=["edge_string", "distance", "inside_radius", "z_valid"]
    )


def calibration_summary(result: dict) -> dict:
    """Coverage and conservativeness of the certified radius, from :func:`calibration_table`.

    * COVERAGE: among elements inside the radius, the fraction where ``Z`` is
      actually valid. This should be 1.0 *by construction* -- ``r_val`` is
      defined as the smallest shell containing a validity failure, so nothing
      strictly closer than ``r_val`` can fail. If it is not 1.0, that is a bug
      in the radius computation or in this calibration check, not a finding
      about the certificate's tightness; a loud warning is raised in that case
      and ``coverage_is_bug`` is set True.
    * CONSERVATIVENESS: among elements outside the radius *that are actually
      reachable* (excluding elements with no defined distance at all, which
      are a separate "unreached" category reported alongside), the fraction
      where ``Z`` nonetheless remained valid. This measures how pessimistic
      the certificate is -- a high value means the radius is a conservative
      (safely small) bound rather than a tight one.

    Args:
        result: Output of :func:`bkrobust.demo.pipeline.run_scenario`.

    Returns:
        A dict with ``coverage``, ``coverage_is_bug``, ``conservativeness``,
        and ``counts`` (raw 2x2-style counts: ``inside_valid``,
        ``inside_invalid``, ``outside_reachable_valid``,
        ``outside_reachable_invalid``, ``n_unreached``, ``n_total``).
    """
    table = calibration_table(result)
    inside = table[table["inside_radius"]]
    reachable = table[table["distance"] != UNREACHED]
    outside_reachable = reachable[~reachable["inside_radius"]]

    coverage = float(inside["z_valid"].mean()) if len(inside) else float("nan")
    coverage_is_bug = len(inside) > 0 and coverage < 1.0
    if coverage_is_bug:
        warnings.warn(
            f"calibration BUG: coverage is {coverage:.4f}, not 1.0 -- some element inside the "
            "certified radius has an invalid adjustment set. r_val is defined as the "
            "smallest shell containing a validity failure, so this should be "
            "impossible; treat this as a correctness bug in the radius computation "
            "or in is_valid_adjustment_set_mpdag, not as a robustness finding.",
            stacklevel=2,
        )

    conservativeness = (
        float(outside_reachable["z_valid"].mean()) if len(outside_reachable) else float("nan")
    )

    counts = {
        "inside_valid": int(inside["z_valid"].sum()),
        "inside_invalid": int((~inside["z_valid"]).sum()),
        "outside_reachable_valid": int(outside_reachable["z_valid"].sum()),
        "outside_reachable_invalid": int((~outside_reachable["z_valid"]).sum()),
        "n_unreached": int((table["distance"] == UNREACHED).sum()),
        "n_total": len(table),
    }

    return {
        "coverage": coverage,
        "coverage_is_bug": coverage_is_bug,
        "conservativeness": conservativeness,
        "counts": counts,
    }


def calibration_for_actual_truth(result: dict) -> dict:
    """The single real calibration data point: the actual ground-truth DAG.

    Reports the model distance from ``G0`` to ``result["truth"]``
    (:func:`~bkrobust.demo.scenario.true_dag`), whether that distance falls
    inside the certified ``r_val``, and whether ``Z`` is genuinely a valid
    adjustment set in the truth (checked directly at the DAG level via
    :func:`~bkrobust.demo.evaluate.is_valid_adjustment_set_dag`, since the
    truth is a fully oriented DAG).

    The true DAG need not be an element of the enumerated space at all (it
    always should be, in this demonstration, since it is one particular
    consistent DAG extension of the CPDAG the space is built from -- but that
    is verified here, not assumed): ``truth_in_space`` and ``truth_in_shells``
    report this explicitly rather than silently treating an absent truth as
    "distance 0" or crashing.

    Args:
        result: Output of :func:`bkrobust.demo.pipeline.run_scenario`.

    Returns:
        A dict with ``truth_in_space``, ``truth_in_shells``,
        ``distance_to_truth``, ``inside_radius``, ``z_valid_in_truth`` and
        ``r_val``.
    """
    truth = result["truth"]
    space = result["space"]
    shells = result["shells"]
    r_val = result["radii"]["r_val"]
    z = result["z"]

    truth_in_space = truth in space
    truth_in_shells = truth in shells
    d = shells.get(truth, UNREACHED)

    z_valid_in_truth = is_valid_adjustment_set_dag(truth, TREATMENT, OUTCOME, z)

    return {
        "truth_in_space": truth_in_space,
        "truth_in_shells": truth_in_shells,
        "distance_to_truth": d,
        "inside_radius": _inside_radius(d, r_val),
        "z_valid_in_truth": z_valid_in_truth,
        "r_val": r_val,
    }
