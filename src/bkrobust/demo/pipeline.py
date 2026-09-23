"""Per-scenario driver: build the space, walk it, compute the radii.

One function, :func:`run_scenario`, does the whole pipeline for one analyst
knowledge state and returns everything the report and the figures need.

Design decisions made here (documented in ``report.md``):

*The adjustment set under test.* ``Z = O*(G0)`` -- the optimal adjustment set of
the graph the analyst actually ends up with. That is what a sophisticated
analyst would adopt, and it is the set the theory recommends. Radii for every
*other* valid adjustment set are computed too, in :func:`robustness_frontier`,
because the comparison between them is the point of the robust-selection
section.

*Bias needs a fully specified model.* A space element ``G`` is generally not a
DAG, so it does not determine a linear-Gaussian SEM on its own. For each ``G``
we therefore instantiate over **every DAG in [G] crossed with several
coefficient draws**, and report the mean and the maximum over that product. The
maximum is a *sampled proxy* for a worst case, not a supremum -- coefficients
are drawn, not optimised over. The report says so plainly.

*Radii.* ``r_val`` is the smallest shell index containing an element where ``Z``
is not valid; ``r_opt`` the smallest containing an element where ``Z`` is not
the optimal set; ``r_eps`` the smallest containing an element whose **mean**
absolute bias reaches ``eps``.

The mean, not the max, keys ``r_eps``. The max over coefficient draws is
unstable -- it moved from 4.36 to 2.68 between two runs differing only in the
number of draws -- because it is an extreme order statistic of a random sample,
not an optimisation over coefficients. The max is still reported per element as
a diagnostic, with that caveat attached. If
no element in the whole (finite, fully enumerated) space breaks the property,
the radius is reported as :data:`UNREACHED` rather than as a large number --
"never breaks in this space" and "breaks far away" are different statements.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from bkrobust.demo.evaluate import (
    LinearSEM,
    adjusted_estimand,
    all_valid_adjustment_sets_mpdag,
    asymptotic_variance,
    is_valid_adjustment_set_dag,
    optimal_adjustment_set_dag,
    optimal_adjustment_set_mpdag,
    random_sem,
)
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.scenario import OUTCOME, TREATMENT, scenarios, true_dag

#: Reported when no element of the enumerated space breaks the property.
UNREACHED = -1

#: Coefficient draws per DAG. Fixed seed; see ``report.md`` for reproduction.
N_DRAWS = 100
SEED = 20260919


def _sem_stats(
    dag: MPDAG,
    z: frozenset[str],
    rng: np.random.Generator,
    n_draws: int = N_DRAWS,
) -> tuple[float, float, float]:
    """Mean |bias|, max |bias| and mean |true effect| over coefficient draws on one DAG."""
    biases: list[float] = []
    effects: list[float] = []
    for _ in range(n_draws):
        sem: LinearSEM = random_sem(dag, rng)
        tau = sem.true_total_effect(TREATMENT, OUTCOME)
        est = adjusted_estimand(sem, TREATMENT, OUTCOME, z)
        biases.append(abs(est - tau))
        effects.append(abs(tau))
    return float(np.mean(biases)), float(np.max(biases)), float(np.mean(effects))


def evaluate_element(
    g: MPDAG,
    z: frozenset[str],
    rng: np.random.Generator,
    n_draws: int = N_DRAWS,
) -> dict[str, Any]:
    """Evaluate ``Z`` at one element of the space.

    Args:
        g: The space element.
        z: The adjustment set under test.
        rng: Seeded generator; the sole source of randomness.
        n_draws: Coefficient draws per DAG in ``[g]``.

    Returns:
        Validity, optimality and bias summaries for ``g``.
    """
    exts = enumerate_dag_extensions(g)
    valid = bool(exts) and all(is_valid_adjustment_set_dag(d, TREATMENT, OUTCOME, z) for d in exts)
    opts = [frozenset(optimal_adjustment_set_dag(d, TREATMENT, OUTCOME)) for d in exts]
    is_opt = bool(opts) and all(o == z for o in opts)

    means, maxes, effs = [], [], []
    for d in exts:
        m, mx, ef = _sem_stats(d, z, rng, n_draws)
        means.append(m)
        maxes.append(mx)
        effs.append(ef)
    return {
        "edge_string": g.edge_string(),
        "n_directed": len(g.directed_edges),
        "n_undirected": len(g.undirected_edges),
        "n_extensions": len(exts),
        "z_valid": valid,
        "z_is_optimal": is_opt,
        "mean_abs_bias": float(np.mean(means)) if means else float("nan"),
        "max_abs_bias": float(np.max(maxes)) if maxes else float("nan"),
        "mean_abs_effect": float(np.mean(effs)) if effs else float("nan"),
    }


def _first_failure(
    rows: list[dict[str, Any]],
    shells: dict[MPDAG, int],
    space: list[MPDAG],
    predicate: str,
    eps: float | None = None,
) -> tuple[int, MPDAG | None]:
    """Smallest shell index at which ``predicate`` fails, and a witness there."""
    best_shell, witness = UNREACHED, None
    for g, row in zip(space, rows):  # noqa: B905 - equal by construction; strict= is 3.10+
        if g not in shells:
            continue
        if predicate == "valid":
            broken = not row["z_valid"]
        elif predicate == "optimal":
            broken = not row["z_is_optimal"]
        else:
            broken = row["mean_abs_bias"] >= float(eps)
        if not broken:
            continue
        s = shells[g]
        if s == 0:
            continue  # G0 itself is the analyst's own graph, not a failure
        if best_shell == UNREACHED or s < best_shell:
            best_shell, witness = s, g
    return best_shell, witness


def run_scenario(
    label: str,
    n_draws: int = N_DRAWS,
    seed: int = SEED,
) -> dict[str, Any]:
    """Run the full pipeline for one scenario.

    Args:
        label: ``"A"``, ``"B"`` or ``"C"``.
        n_draws: Coefficient draws per DAG.
        seed: Root seed. All randomness descends from it.

    Returns:
        A dict with the space, shells, per-element rows, radii, witnesses, the
        metric check, and the robustness frontier.

    Raises:
        ValueError: If the scenario's knowledge is inconsistent with the CPDAG.
    """
    from bkrobust.demo.space import (
        all_pairs_distances,
        atomic_moves,
        bfs_distances,
        check_metric_axioms,
        covering_pairs,
        enumerate_space,
        neighbour_graph,
        represented_dags,
    )

    truth = true_dag()
    cpdag = dag_to_cpdag(truth)
    spec = scenarios()[label]
    g0 = apply_orientations(cpdag, spec["knowledge"])
    if g0 is None:
        raise ValueError(f"scenario {label} is inconsistent with the CPDAG")

    space = enumerate_space(cpdag)
    reps = represented_dags(space)
    covers = covering_pairs(space, reps)
    nbrs = neighbour_graph(space, covers)
    shells = bfs_distances(nbrs, g0)
    metric = check_metric_axioms(all_pairs_distances(nbrs), space)

    z_opt = optimal_adjustment_set_mpdag(g0, TREATMENT, OUTCOME)
    z = frozenset(z_opt) if z_opt is not None else frozenset()

    rng = np.random.default_rng(seed)
    rows = [evaluate_element(g, z, rng, n_draws) for g in space]
    for g, row in zip(space, rows):  # noqa: B905 - equal by construction; strict= is 3.10+
        row["shell"] = shells.get(g, UNREACHED)
        row["in_component"] = g in shells

    mean_eff = float(np.nanmean([r["mean_abs_effect"] for r in rows]))
    # One epsilon tied to a plausible standard error, so the radius can be read
    # against what a 95% interval at n=1000 would actually notice; two tied to
    # the effect magnitude itself.
    se_rng = np.random.default_rng(seed + 1)  # one generator, drawn from 200 times
    se = float(
        np.sqrt(
            np.mean(
                [
                    asymptotic_variance(random_sem(truth, se_rng), TREATMENT, OUTCOME, z)
                    for _ in range(200)
                ]
            )
            / 1000.0
        )
    )
    eps_grid = {
        "eps_se": 2.0 * se,
        "eps_small": 0.05 * mean_eff,
        "eps_medium": 0.25 * mean_eff,
    }

    radii: dict[str, Any] = {}
    witnesses: dict[str, Any] = {}
    for name, pred, ev in (
        ("r_val", "valid", None),
        ("r_opt", "optimal", None),
        ("r_eps_se", "eps", eps_grid["eps_se"]),
        ("r_eps_small", "eps", eps_grid["eps_small"]),
        ("r_eps_medium", "eps", eps_grid["eps_medium"]),
    ):
        r, w = _first_failure(rows, shells, space, pred, ev)
        radii[name] = r
        witnesses[name] = (
            None
            if w is None
            else {
                "edge_string": w.edge_string(),
                "shell": shells[w],
                "moves": atomic_moves(g0, w),
            }
        )

    return {
        "label": label,
        "spec": spec,
        "cpdag": cpdag,
        "truth": truth,
        "g0": g0,
        "z": z,
        "space": space,
        "shells": shells,
        "covers": covers,
        "neighbours": nbrs,
        "rows": rows,
        "radii": radii,
        "witnesses": witnesses,
        "eps_grid": eps_grid,
        "metric": metric,
        "mean_abs_effect": mean_eff,
        "standard_error_n1000": se,
    }


def robustness_frontier(
    result: dict[str, Any],
    n_draws: int = 25,
    seed: int = SEED,
) -> list[dict[str, Any]]:
    """``r_val`` and asymptotic variance for EVERY valid adjustment set of ``G0``.

    The efficiency-robustness tradeoff. If the optimal adjustment set turns out
    to be the least robust, that is the headline of the demonstration.

    Args:
        result: Output of :func:`run_scenario`.
        n_draws: Coefficient draws for the variance estimate.
        seed: Root seed.

    Returns:
        One row per valid adjustment set: the set, its size, its ``r_val``, its
        mean asymptotic variance, and whether it is ``O*(G0)``.
    """
    g0, space, shells = result["g0"], result["space"], result["shells"]
    truth = result["truth"]
    rng = np.random.default_rng(seed)
    out: list[dict[str, Any]] = []
    for zs in all_valid_adjustment_sets_mpdag(g0, TREATMENT, OUTCOME):
        z = frozenset(zs)
        rows = []
        for g in space:
            exts = enumerate_dag_extensions(g)
            valid = bool(exts) and all(
                is_valid_adjustment_set_dag(d, TREATMENT, OUTCOME, z) for d in exts
            )
            rows.append({"z_valid": valid, "z_is_optimal": False, "max_abs_bias": 0.0})
        r_val, _ = _first_failure(rows, shells, space, "valid")
        variances = [
            asymptotic_variance(random_sem(truth, rng), TREATMENT, OUTCOME, z)
            for _ in range(n_draws)
        ]
        out.append(
            {
                "z": sorted(z),
                "size": len(z),
                "r_val": r_val,
                "asymptotic_variance": float(np.mean(variances)),
                "is_optimal_at_g0": z == result["z"],
            }
        )
    return out
