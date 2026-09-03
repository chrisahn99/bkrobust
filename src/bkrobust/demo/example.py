"""The running example: candidate scenarios and the design gates they must pass.

The brief requires a single worked example, chosen deliberately rather than
sampled, with interpretable variable names. Random DAGs would satisfy the
structural criteria and fail the interpretability one -- the whole point of the
witness graphs is that a domain expert can read the failure and disagree with
it -- so candidates are hand-designed from plausible clinical stories and then
screened computationally against the gates below.

Design gates (all must pass):

G1  8-10 vertices.
G2  The CPDAG has an undirected (chordal) component with 4-6 edges.
G3  X and Y admit at least two distinct valid adjustment sets in the truth.
G4  At least one undirected edge lies on, or adjacent to, a back-door path
    from X to Y -- otherwise no perturbation can change anything and the
    demonstration is vacuous.
G5  Sanity gate: at least one SINGLE atomic perturbation of ``K_true`` changes
    the validity status of at least one adjustment set.

``check_candidate`` reports pass/fail per gate so the rejected candidates can be
counted honestly in the report.
"""

from __future__ import annotations

from bkrobust.demo.graph import MPDAG, canon, undirected_components

Edge = tuple[str, str]


def dag_to_cpdag(dag: MPDAG) -> MPDAG:
    """Return the CPDAG of a DAG: keep the v-structures, then close under Meek.

    The essential graph of the Markov equivalence class. Every edge that is not
    compelled -- by a v-structure or by Meek propagation from one -- comes back
    undirected, and those undirected edges are exactly what background knowledge
    is in a position to orient.

    Args:
        dag: A fully oriented, acyclic MPDAG.

    Returns:
        The CPDAG.

    Raises:
        ValueError: If ``dag`` is not a DAG, or if the closure unexpectedly
            FAILs (which would indicate a bug in the closure, not bad input).
    """
    from bkrobust.demo.meek import meek_closure

    if not dag.is_dag():
        raise ValueError("dag_to_cpdag expects a fully oriented acyclic graph")

    v_structures: set[Edge] = set()
    for b in dag.nodes:
        parents = sorted(dag.parents(b))
        for i, a in enumerate(parents):
            for c in parents[i + 1 :]:
                if not dag.has_edge(a, c):
                    v_structures.add((a, b))
                    v_structures.add((c, b))

    skeleton = dag.skeleton()
    start = MPDAG(
        dag.nodes,
        directed=v_structures,
        undirected=[e for e in skeleton if e not in {canon(a, b) for a, b in v_structures}],
    )
    closed = meek_closure(start)
    if closed is None:
        raise ValueError("Meek closure FAILed on a CPDAG construction -- closure bug")
    return closed


def knowledge_to_recover(dag: MPDAG, cpdag: MPDAG) -> set[Edge]:
    """A minimal set of orientations that turns ``cpdag`` back into ``dag``.

    Greedy and deterministic: repeatedly assert the true orientation of an
    undirected edge, take the Meek closure, and stop when nothing is left
    undirected. Greedy minimality is enough here -- ``K_true`` only has to be a
    small honest set of claims an analyst might supply, not a certified minimum.

    Args:
        dag: The ground-truth DAG.
        cpdag: Its CPDAG.

    Returns:
        The asserted orientations, as ``(tail, head)`` pairs.
    """
    from bkrobust.demo.meek import meek_closure

    current = cpdag
    asserted: set[Edge] = set()
    while current.undirected_edges:
        best: Edge | None = None
        best_gain = -1
        for a, b in sorted(current.undirected_edges):
            tail, head = (a, b) if dag.is_directed_edge(a, b) else (b, a)
            trial = meek_closure(current.oriented(tail, head))
            if trial is None:
                continue
            gain = len(current.undirected_edges) - len(trial.undirected_edges)
            if gain > best_gain:
                best_gain, best = gain, (tail, head)
        if best is None:
            raise ValueError("could not orient remaining undirected edges toward the truth")
        asserted.add(best)
        closed = meek_closure(current.oriented(*best))
        assert closed is not None
        current = closed
    return asserted


def backdoor_adjacent_undirected(dag: MPDAG, cpdag: MPDAG, x: str, y: str) -> set[Edge]:
    """Undirected edges of ``cpdag`` touching a node on some back-door path in ``dag``.

    Gate G4. An undirected edge far from every back-door path cannot change the
    validity of any adjustment set, so an example whose undirected component is
    disjoint from the relevant paths would be vacuous.
    """
    relevant: set[str] = set()
    for node in dag.nodes:
        if node in (x, y):
            continue
        # A node matters if it can sit on a back-door path: it is an ancestor of
        # x or of y (or is adjacent to one), and is not a descendant of x.
        if node in dag.descendants(x):
            continue
        if node in dag.ancestors(x) or node in dag.ancestors(y):
            relevant.add(node)
    relevant |= {x, y}
    return {e for e in cpdag.undirected_edges if e[0] in relevant or e[1] in relevant}


def check_candidate(
    name: str,
    story: str,
    edges: list[Edge],
    x: str,
    y: str,
) -> dict[str, object]:
    """Screen one candidate scenario against the five design gates.

    Args:
        name: Short candidate identifier.
        story: One-line description of the clinical scenario.
        edges: The ground-truth DAG as ``(tail, head)`` pairs.
        x: Treatment.
        y: Outcome.

    Returns:
        A report dict with a boolean per gate, the derived objects, and a
        ``passed`` flag that is True only if every gate passed.
    """
    from bkrobust.demo.evaluate import all_valid_adjustment_sets_mpdag
    from bkrobust.demo.meek import apply_orientations

    nodes = sorted({n for e in edges for n in e})
    dag = MPDAG(nodes, directed=edges)
    cpdag = dag_to_cpdag(dag)
    comps = undirected_components(cpdag)
    comp_edge_counts = [
        len([e for e in cpdag.undirected_edges if e[0] in c and e[1] in c]) for c in comps
    ]

    g1 = 8 <= len(nodes) <= 10
    g2 = any(4 <= c <= 6 for c in comp_edge_counts)

    valid_true = all_valid_adjustment_sets_mpdag(dag, x, y)
    g3 = len(valid_true) >= 2

    g4 = len(backdoor_adjacent_undirected(dag, cpdag, x, y)) > 0

    # G5: does any single atomic perturbation of K_true flip a validity status?
    k_true = knowledge_to_recover(dag, cpdag)
    g5 = False
    g5_witness: str | None = None
    if g3:
        for drop in sorted(k_true):
            remaining = sorted(k_true - {drop})
            g0 = apply_orientations(cpdag, remaining)
            if g0 is None:
                continue
            for z in valid_true:
                from bkrobust.demo.evaluate import is_valid_adjustment_set_mpdag

                if not is_valid_adjustment_set_mpdag(g0, x, y, z):
                    g5 = True
                    g5_witness = "omitting {}->{} invalidates Z={}".format(
                        drop[0], drop[1], sorted(z) or "{}"
                    )
                    break
            if g5:
                break

    gates = {
        "G1_size": g1,
        "G2_component": g2,
        "G3_two_sets": g3,
        "G4_adjacent": g4,
        "G5_sanity": g5,
    }
    return {
        "name": name,
        "story": story,
        "gates": gates,
        "passed": all(gates.values()),
        "n_nodes": len(nodes),
        "component_edge_counts": comp_edge_counts,
        "n_valid_sets_true": len(valid_true),
        "valid_sets_true": valid_true,
        "k_true": sorted(k_true),
        "dag": dag,
        "cpdag": cpdag,
        "g5_witness": g5_witness,
    }
