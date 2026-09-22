"""Survival on **elicited** background knowledge, across the model panel.

What this module is for, and why it is a new file
-------------------------------------------------
:mod:`bkrobust.robustness.real_survival` builds the analyst's claim set with
``benchmarks.measure.select_knowledge(dag, cpdag, coverage)``. That function
reads the **ground-truth DAG**: at every coverage level it returns claims that
are correct by construction, and the only wrongness in that design is the
synthetic ``base_wrongness`` dial applied on top. The elicited claims in
``results/elicit/knowledge.json`` -- what a model actually answered when asked
about each undirected pair -- were never on that path.

This module puts them on it. ``K`` is read from ``knowledge.json`` for a
``(condition, network)`` pair and ``select_knowledge`` is never called, so
``G0`` is the Meek closure of claims an analyst really supplied, wrong ones
included, and ``shd_truth`` is a measured error rather than a dial setting.

``real_survival`` is **not modified**. Its docstring records the repository's
standing rule -- every number in ``report_fragility_and_pareto.md`` and
``results/axis_robustness_real/`` must stay reproducible from the code that
produced it -- and the committed 831-row corpus depends on it. This module is
its elicited-knowledge sibling: it imports ``real_survival``'s grid, seeding,
closure-cache, radius and cell machinery wholesale so the two arms cannot
drift apart on any definition, and it owns a separate results tree.

The model panel is the variation axis
--------------------------------------
On the ``select_knowledge`` corpus, ``|K|`` is a deterministic function of
``(network, coverage)`` and ``shd_truth`` is identically zero at
``base_wrongness = 0``, so within a stratum there is nothing to rank. The
panel supplies both missing dimensions empirically:

* ``|K|`` varies **within a network** across conditions -- ``sachs`` ranges
  over eight distinct sizes from 1 to 14 claims, ``ecoli70`` over ten from 3
  to 23 -- because a bigger model asserts on more of the pairs it is asked
  about instead of declining.
* per-condition accuracy over the frame's networks runs from about 0.61 to
  about 0.74, so ``G0`` is genuinely wrong at a rate that differs by model.

Neither dimension is a knob anyone set. Both are recorded per row.

Accuracy is a measured property of ``K``, never an input
---------------------------------------------------------
:func:`k_accuracy` compares each claim to the true DAG orientation. It is a
**diagnostic column**, written next to the survival endpoint so the analysis
can condition on it. Nothing in the corruption process, the grid, the seeding
or the endpoint reads it, and no ``G0`` is ever built from the truth.

Degenerate cases are recorded, never skipped
---------------------------------------------
A ``(condition, network)`` cell can fail to produce a usable ``G0`` three
ways, and they are three different facts:

``n_k_zero``
    The model declined every pair it was asked about. A real elicitation
    outcome, not an error -- ``D_LLM_32B`` did this on ``Acid_1996``.
``k_contradictory``
    The asserted claims admit no consistent MPDAG. None were observed on this
    corpus; the branch exists because absence is a finding that must be
    measured rather than assumed.
``o_g0_extensions_intractable``
    ``G0`` has too many undirected edges left for the extension enumeration.
    A **measurement limit** of this machine, kept distinct from the two
    structural reasons above.

``pathfinder`` carries three frame rows and is absent from ``knowledge.json``
entirely -- it was never elicited. Its rows produce no shard and the driver
stamps the exclusion in its manifest rather than letting the network vanish
silently from a network count.

Sentinel, RNG and gate discipline are inherited unchanged from
``real_survival``: ``UNREACHED`` is a status and never averaged, every draw is
seeded by a SHA-256 derivation so no global RNG is touched, and every row
carries the literal ``"fast_gate"``.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from bkrobust.demo.graph import MPDAG
from bkrobust.robustness import real_survival as rs

Edge = tuple[str, str]

#: The committed elicitation bundle every condition is read from.
DEFAULT_ELICIT_PATH = Path("results/elicit/knowledge.json")

#: Conditions whose node names are the real domain names. These are the panel
#: proper: the ranking question is about knowledge an analyst would actually
#: supply, and a scrambled-name arm is a semantic control, not an analyst.
REAL_NAMING = "real"

#: The scrambled-name control arm. Run and reported **beside** the panel, never
#: pooled into it: its claims are answers about relabelled variables, so its
#: accuracy measures what survives when the domain semantics are destroyed.
SCRAMBLED_NAMING = "scrambled"


# --- the bundle ---------------------------------------------------------------


def load_bundle(path: Path = DEFAULT_ELICIT_PATH) -> tuple[dict[str, Any], str]:
    """Read ``knowledge.json`` and its digest, together.

    The digest is stamped on every shard marker so a results tree can be tied
    back to the exact elicitation bundle it was computed from.

    Args:
        path: The bundle file.

    Returns:
        ``(bundle, sha256_hex)``.

    Raises:
        FileNotFoundError: If the bundle is absent.
    """
    if not path.is_file():
        raise FileNotFoundError(f"elicitation bundle missing: {path}")
    raw = path.read_bytes()
    return json.loads(raw), hashlib.sha256(raw).hexdigest()


def panel_conditions(bundle: dict[str, Any], *, naming: str | None = REAL_NAMING) -> list[str]:
    """The condition ids of the panel, in a deterministic order.

    Args:
        bundle: The loaded bundle.
        naming: ``"real"`` for the panel proper, ``"scrambled"`` for the
            control arm, or ``None`` for every condition in the bundle.

    Returns:
        Sorted condition ids.
    """
    if naming is None:
        return sorted(bundle)
    return sorted(c for c in bundle if bundle[c].get("naming") == naming)


def condition_meta(bundle: dict[str, Any], condition: str) -> dict[str, Any]:
    """The descriptive columns every row of a condition carries.

    Args:
        bundle: The loaded bundle.
        condition: A condition id.

    Returns:
        ``condition``, ``model``, ``family``, ``naming``, ``elicit_arm``,
        ``order_seed``, ``questionnaire_sha256`` and ``true_dag_on_path``.

    Raises:
        KeyError: If the condition is not in the bundle.
    """
    c = bundle[condition]
    return {
        "condition": condition,
        "model": c.get("model"),
        "family": c.get("family"),
        "naming": c.get("naming"),
        "elicit_arm": c.get("arm"),
        "order_seed": c.get("order_seed"),
        "questionnaire_sha256": c.get("questionnaire_sha256"),
        "true_dag_on_path": c.get("true_dag_on_path"),
    }


def has_network(bundle: dict[str, Any], condition: str, network: str) -> bool:
    """Whether a condition was ever asked about a network.

    Args:
        bundle: The loaded bundle.
        condition: A condition id.
        network: A network name.

    Returns:
        ``True`` if the network appears under that condition.
    """
    return network in bundle.get(condition, {}).get("networks", {})


def elicited_k(
    bundle: dict[str, Any], condition: str, network: str
) -> tuple[tuple[Edge, ...], dict[str, Any]]:
    """The claims a condition actually asserted for a network, and its audit.

    This is the whole point of the module: ``K`` comes from the questionnaire
    answers, so it may be wrong, partial or empty. ``select_knowledge`` is not
    called and the ground-truth DAG is not read here.

    Args:
        bundle: The loaded bundle.
        condition: A condition id.
        network: A network name.

    Returns:
        ``(k, elicitation_columns)``. ``k`` is sorted so the seed derivation
        and the ``k_b_sha256`` column are order-independent.

    Raises:
        KeyError: If the condition never covered that network.
    """
    nets = bundle[condition]["networks"]
    if network not in nets:
        raise KeyError(f"{condition} has no network {network!r}")
    rec = nets[network]
    k = tuple(sorted((str(a), str(b)) for a, b in rec["k"]))
    n_asked = rec.get("n_asked")
    cols = {
        "k_source": "knowledge.json",
        "n_asked": n_asked,
        "n_asserted": rec.get("n_asserted"),
        "n_asserted_pre": rec.get("n_asserted_pre"),
        "n_declined": rec.get("declined"),
        "n_parse_fail": rec.get("parse_fail"),
        "n_not_reached": rec.get("not_reached"),
        "n_invented": rec.get("invented"),
        "n_contradicts_order": rec.get("contradicts_order"),
        "n_dropped_for_consistency": rec.get("dropped_for_consistency"),
        "meek_consistent_pre": rec.get("meek_consistent_pre"),
        "assert_rate": (len(k) / n_asked) if n_asked else None,
        "k_sha256_elicit": rec.get("k_sha256"),
    }
    return k, cols


# --- measured properties of an elicited K -------------------------------------


def k_accuracy(k: Sequence[Edge], dag: MPDAG) -> dict[str, Any]:
    """How many of the asserted claims point the way the true DAG does.

    A **diagnostic**, computed after ``K`` is fixed and never fed back into the
    corruption process, the grid, the seeding or the endpoint.

    Args:
        k: The elicited claims.
        dag: The ground-truth DAG.

    Returns:
        ``n_k_correct``, ``n_k_wrong`` and ``k_accuracy`` -- the last is
        ``None``, never zero, when ``K`` is empty.
    """
    truth = set(dag.directed_edges)
    n_correct = sum(1 for e in k if tuple(e) in truth)
    return {
        "n_k_correct": n_correct,
        "n_k_wrong": len(k) - n_correct,
        "k_accuracy": (n_correct / len(k)) if k else None,
    }


def off_skeleton(k: Sequence[Edge], cpdag: MPDAG) -> list[Edge]:
    """Claims whose unordered pair is not an undirected edge of the CPDAG.

    The questionnaire only ever asked about the CPDAG's undirected pairs, so
    this should be empty for every condition and network. It is checked on
    every shard rather than trusted, because a non-empty result would mean the
    answers were parsed against a different graph than the one being scored.

    Args:
        k: The elicited claims.
        cpdag: The oracle CPDAG.

    Returns:
        The offending claims, in input order.
    """
    und = {tuple(sorted(e)) for e in cpdag.undirected_edges}
    return [tuple(e) for e in k if tuple(sorted(e)) not in und]


def k_diagnostics(k: Sequence[Edge], dag: MPDAG, cpdag: MPDAG) -> dict[str, Any]:
    """Every measured column an elicited ``K`` carries, in one call.

    Args:
        k: The elicited claims.
        dag: The ground-truth DAG.
        cpdag: The oracle CPDAG.

    Returns:
        :func:`k_accuracy`'s columns plus ``n_off_skeleton``,
        ``off_skeleton_claims`` and ``n_cpdag_undirected``.
    """
    off = off_skeleton(k, cpdag)
    return {
        **k_accuracy(k, dag),
        "n_off_skeleton": len(off),
        "off_skeleton_claims": [list(e) for e in off] or None,
        "n_cpdag_undirected": len(cpdag.undirected_edges),
    }


# --- the instance population --------------------------------------------------


def distinct_pairs(frame: Sequence[dict[str, Any]], network: str) -> list[dict[str, Any]]:
    """One frame row per distinct ``(X, Y)`` pair of a network.

    Elicited ``K`` is a property of ``(condition, network)`` and carries no
    coverage parameter, so a network's frame rows -- which repeat a pair once
    per coverage level -- are deduplicated exactly the way
    ``run_real_survival``'s tiered arm deduplicates them. The first row in
    frame order wins, which is deterministic because the frame is frozen.

    Args:
        frame: The frozen frame rows.
        network: The network to select.

    Returns:
        The deduplicated rows, in frame order.
    """
    seen: set[tuple[str, str]] = set()
    out: list[dict[str, Any]] = []
    for r in frame:
        if r["network"] != network:
            continue
        key = (r["x"], r["y"])
        if key in seen:
            continue
        seen.add(key)
        out.append(r)
    return out


def panel_family(condition: str, network: str) -> str:
    """The seed-derivation family for a ``(condition, network)`` shard.

    Distinct from every family ``run_real_survival`` uses -- those are built
    from ``(network, coverage, base_wrongness, replicate)`` and never contain a
    condition id -- so no draw of this sweep can collide with a committed one.

    Args:
        condition: The condition id.
        network: The network name.

    Returns:
        The family string.
    """
    return f"llm|{condition}|{network}"


def build_state(
    bundle: dict[str, Any], condition: str, network: str, dag: MPDAG, cpdag: MPDAG
) -> tuple[tuple[Edge, ...], MPDAG | None, str, dict[str, Any]]:
    """Read ``K`` for one ``(condition, network)`` and Meek-close it.

    Args:
        bundle: The loaded bundle.
        condition: The condition id.
        network: The network name.
        dag: The ground-truth DAG, read only for the accuracy diagnostic.
        cpdag: The oracle CPDAG.

    Returns:
        ``(k, g0_or_None, g0_status, shared_columns)``. ``g0_status`` is one of
        ``"ok"``, ``"n_k_zero"``, ``"k_contradictory"`` or
        ``"o_g0_extensions_intractable"``.
    """
    k, elicit_cols = elicited_k(bundle, condition, network)
    g0, reason = rs.build_g0(cpdag, k)
    shared: dict[str, Any] = {
        **condition_meta(bundle, condition),
        **elicit_cols,
        **k_diagnostics(k, dag, cpdag),
        "n_k": len(k),
        "k_b_sha256": rs.sha_of(sorted(k)),
        "g0_status": reason,
        "g0_sha256": rs.sha_of(sorted(g0.directed_edges)) if g0 is not None else None,
        "g0_undirected_edges": len(g0.undirected_edges) if g0 is not None else None,
        "k_g0": rs.commitment_size(g0, cpdag) if g0 is not None else None,
        "shd_truth": rs.directed_symdiff(g0, dag) if g0 is not None else None,
        # Held for schema compatibility with the committed real arm, where they
        # are the synthetic corruption dials. On this arm there is no dial: the
        # wrongness is whatever the model said, counted in `n_k_wrong`.
        "base_wrongness": None,
        "bw_abs": None,
        "bw_is_inert": False,
        "n_claims_actually_wrong": len(k) - k_accuracy(k, dag)["n_k_correct"],
        "depth_grid": "frac10",
    }
    return k, g0, reason, shared


__all__ = [
    "DEFAULT_ELICIT_PATH", "REAL_NAMING", "SCRAMBLED_NAMING",
    "build_state", "condition_meta", "distinct_pairs", "elicited_k",
    "has_network", "k_accuracy", "k_diagnostics", "load_bundle",
    "off_skeleton", "panel_conditions", "panel_family",
]
