"""The end-to-end certificate, and the units it reports in.

The unit test at the bottom of this file exists because the bug it pins was real:
the staircase is computed in effect units while the certificate reports in
relative units, and returning one where the other was expected made a sound
certificate look violated by a factor of ``|theta_Z|``. Anything that compares
``worst_case_at``, ``realised_error`` and ``beta_top`` has to have them in the
same units.
"""

from __future__ import annotations

import numpy as np
import pytest

from bkrobust.core.conventions import UNREACHED
from bkrobust.demo.evaluate import LinearSEM, optimal_adjustment_set_mpdag, random_sem
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations
from bkrobust.epsilon.bias import knowledge_of
from bkrobust.epsilon.certify import certify, describe


def _setup(seed: int = 5) -> tuple[MPDAG, MPDAG, MPDAG, frozenset[str], LinearSEM]:
    rng = np.random.default_rng(seed)
    dag = MPDAG(
        ["V0", "V1", "V2", "V3", "V4"],
        directed=[
            ("V0", "V1"),
            ("V0", "V2"),
            ("V1", "V3"),
            ("V2", "V3"),
            ("V3", "V4"),
            ("V2", "V4"),
        ],
    )
    cpdag = dag_to_cpdag(dag)
    k = sorted(set(dag.directed_edges) - set(cpdag.directed_edges))
    g0 = apply_orientations(cpdag, k)
    assert g0 is not None
    z = optimal_adjustment_set_mpdag(g0, "V2", "V4")
    assert z is not None
    return dag, cpdag, g0, frozenset(z), random_sem(dag, rng)


def test_certificate_runs_end_to_end_from_a_sem() -> None:
    _dag, cpdag, g0, z, sem = _setup()
    cert = certify(cpdag, None, "V2", "V4", sem=sem, g0=g0, z=z, units="relative")
    assert cert.status == "ok"
    assert cert.r_val >= 1
    assert cert.n_knowledge == len(knowledge_of(cpdag, g0))
    assert set(cert.r_eps) == set(cert.epsilons)
    assert isinstance(describe(cert), str) and describe(cert)


def test_certificate_runs_from_a_covariance_with_no_true_graph() -> None:
    """The analyst path: no SEM, no truth, and the audit field is honestly ``nan``."""
    dag, cpdag, g0, z, sem = _setup()
    cert = certify(
        cpdag,
        None,
        "V2",
        "V4",
        sigma=sem.covariance(),
        nodes=dag.nodes,
        g0=g0,
        z=z,
        units="absolute",
    )
    assert cert.status == "ok"
    assert np.isnan(cert.realised_error), "no truth was supplied, so none may be reported"


def test_worst_case_at_is_in_the_same_units_as_realised_error() -> None:
    """The regression test for the unit bug described in the module docstring."""
    _dag, cpdag, g0, z, sem = _setup()
    rel = certify(cpdag, None, "V2", "V4", sem=sem, g0=g0, z=z, units="relative")
    absolute = certify(cpdag, None, "V2", "V4", sem=sem, g0=g0, z=z, units="absolute")
    k = rel.n_knowledge
    # The same staircase, read in two units, must differ by exactly |theta_Z|.
    assert rel.worst_case_at(k) == pytest.approx(absolute.worst_case_at(k) / abs(rel.theta_z))
    assert rel.beta_top == pytest.approx(absolute.beta_top / abs(rel.theta_z))
    # And the ceiling must dominate every point of the staircase, in either unit.
    for cert in (rel, absolute):
        assert cert.worst_case_at(cert.n_knowledge) <= cert.beta_top + 1e-9


def test_the_certificate_bounds_the_realised_error_at_the_right_budget() -> None:
    """Theorem D, audited -- and audited in the unit the budget is actually denominated in.

    The budget is the number of orientations of the CLOSURE that are false of the
    truth, not the number of asserted claims: Meek's rules cascade, so a smaller
    number of wrong assertions can produce more wrong closure orientations.
    Getting this wrong is what made the certificate look violated once.
    """
    rng = np.random.default_rng(0)
    checked = 0
    for seed in range(40):
        dag, cpdag, g0, _z, sem = _setup(seed)
        # Give the analyst some genuinely false knowledge by flipping one claim.
        k0 = knowledge_of(cpdag, g0)
        if not k0:
            continue
        flipped = [(b, a) if (a, b) == k0[0] else (a, b) for a, b in k0]
        g0b = apply_orientations(cpdag, flipped)
        if g0b is None:
            continue
        zb = optimal_adjustment_set_mpdag(g0b, "V2", "V4")
        if zb is None:
            continue
        cert = certify(cpdag, None, "V2", "V4", sem=sem, g0=g0b, z=frozenset(zb), units="absolute")
        if cert.r_val == UNREACHED:
            continue
        k_false = len([e for e in knowledge_of(cpdag, g0b) if e not in set(dag.directed_edges)])
        assert cert.realised_error <= cert.worst_case_at(k_false) + 1e-9, (
            f"certificate violated at seed {seed}: realised {cert.realised_error} "
            f"> beta_up({k_false}) = {cert.worst_case_at(k_false)}"
        )
        checked += 1
        _ = rng  # determinism is via _setup's own seed; nothing else draws here
    assert checked > 0, "the audit did not actually exercise any instance"


def test_relative_units_are_refused_when_the_estimate_is_zero() -> None:
    """A percentage of zero is not a number, and must not silently become one."""
    dag = MPDAG(["A", "B", "X", "Y"], directed=[("A", "X"), ("B", "X"), ("A", "Y")])
    cpdag = dag_to_cpdag(dag)
    sem = random_sem(dag, np.random.default_rng(1))
    g0 = apply_orientations(cpdag, [])
    assert g0 is not None
    # X has no causal path to Y here, so the reported estimate is zero.
    with pytest.raises(ValueError, match="relative units are undefined"):
        certify(cpdag, None, "X", "Y", sem=sem, g0=g0, z=frozenset({"A"}), units="relative")


def test_describe_never_invents_a_guarantee_when_r_val_did_not_complete() -> None:
    _dag, cpdag, g0, z, sem = _setup()
    cert = certify(cpdag, None, "V2", "V4", sem=sem, g0=g0, z=z)
    cert.status = "r_val_not_exact"
    assert "must not be read as a guarantee" in describe(cert)
