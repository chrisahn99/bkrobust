"""The sampler contract: consistent, false, at the requested radius, of the requested kind.

Everything downstream assumes these four properties. If a sampled knowledge set
is inconsistent, it is caught by the check practitioners already run and is not
what the paper is about; if it is accidentally true, it belongs in the control
arm; if it is off-radius, every curve plotted against delta is mislabelled.

These are the tests to run first when a result looks surprising.
"""

from __future__ import annotations

import pytest

from bkrobust.knowledge import perturb
from bkrobust.knowledge.taxonomy import MisspecificationKind


def test_sample_is_consistent(mpdag, chain, rng):
    """Every sample passes Meek's Algorithm 1 -- part 1 of the contract."""
    from bkrobust.graphs import consistency

    cpdag, true_graph = mpdag(chain, as_cpdag=True), mpdag(chain)
    with pytest.raises(NotImplementedError):
        bk = perturb.sample_consistent_but_false(
            true_graph, cpdag, delta=1, kind=MisspecificationKind.ORIENTATION, rng=rng
        )
        assert consistency.is_consistent(cpdag, bk)


def test_sample_is_false(mpdag, chain, rng):
    """Every sample contradicts the true DAG -- part 2 of the contract."""
    from bkrobust.graphs import consistency

    cpdag, true_graph = mpdag(chain, as_cpdag=True), mpdag(chain)
    with pytest.raises(NotImplementedError):
        bk = perturb.sample_consistent_but_false(
            true_graph, cpdag, delta=1, kind=MisspecificationKind.ORIENTATION, rng=rng
        )
        assert not consistency.is_true(true_graph, bk)


def test_sample_is_at_the_requested_radius(mpdag, chain, rng):
    """The realised radius equals the requested one -- part 3 of the contract."""
    cpdag, true_graph = mpdag(chain, as_cpdag=True), mpdag(chain)
    with pytest.raises(NotImplementedError):
        bk = perturb.sample_consistent_but_false(
            true_graph,
            cpdag,
            delta=1,
            kind=MisspecificationKind.ORIENTATION,
            rng=rng,
            exact_radius=True,
        )
        assert bk.metadata["realised_delta"] == 1


def test_sample_is_of_the_requested_kind(mpdag, chain, rng):
    """Every false constraint is of the requested kind -- part 4 of the contract."""
    cpdag, true_graph = mpdag(chain, as_cpdag=True), mpdag(chain)
    with pytest.raises(NotImplementedError):
        bk = perturb.sample_consistent_but_false(
            true_graph, cpdag, delta=1, kind=MisspecificationKind.ORIENTATION, rng=rng
        )
        assert bk.metadata["kind"] == MisspecificationKind.ORIENTATION


def test_delta_zero_returns_true_knowledge(mpdag, chain, rng):
    """delta=0 is the control arm: true knowledge, and the falsity requirement is waived."""
    from bkrobust.graphs import consistency

    cpdag, true_graph = mpdag(chain, as_cpdag=True), mpdag(chain)
    with pytest.raises(NotImplementedError):
        bk = perturb.sample_consistent_but_false(
            true_graph, cpdag, delta=0, kind=MisspecificationKind.ORIENTATION, rng=rng
        )
        assert consistency.is_true(true_graph, bk)


def test_sampling_is_reproducible(mpdag, chain):
    """Two generators from the same seed give the same sample.

    Rejection sampling consumes a variable number of draws, so a sampler that
    reaches for any RNG other than the one passed in will fail this.
    """
    import numpy as np

    cpdag, true_graph = mpdag(chain, as_cpdag=True), mpdag(chain)
    with pytest.raises(NotImplementedError):
        a = perturb.sample_consistent_but_false(
            true_graph, cpdag, 1, MisspecificationKind.ORIENTATION, np.random.default_rng(0)
        )
        b = perturb.sample_consistent_but_false(
            true_graph, cpdag, 1, MisspecificationKind.ORIENTATION, np.random.default_rng(0)
        )
        assert a == b


def test_exhaustion_raises_rather_than_returning_off_radius(mpdag, collider, rng):
    """When the target set is empty, the sampler raises rather than approximating.

    The collider's CPDAG is fully oriented, so there is nothing consistent left
    to get wrong. Returning a near-miss here would put mislabelled rows into
    every downstream table.
    """
    cpdag, true_graph = mpdag(collider, as_cpdag=True), mpdag(collider)
    with pytest.raises((NotImplementedError, perturb.PerturbationExhaustedError)):
        perturb.sample_consistent_but_false(
            true_graph,
            cpdag,
            delta=1,
            kind=MisspecificationKind.ORIENTATION,
            rng=rng,
            max_rejections=50,
        )


def test_rejection_count_is_recorded(mpdag, chain, rng):
    """Rejection counts are data about how determined the CPDAG is, and are kept."""
    cpdag, true_graph = mpdag(chain, as_cpdag=True), mpdag(chain)
    with pytest.raises(NotImplementedError):
        bk = perturb.sample_consistent_but_false(
            true_graph, cpdag, delta=1, kind=MisspecificationKind.ORIENTATION, rng=rng
        )
        assert "n_rejections" in bk.metadata


def test_true_arm_matches_false_arm_in_size(mpdag, chain, rng):
    """The control arm asserts as many constraints as the treatment arm.

    Knowledge shrinks the equivalence class whether or not it is correct. Without
    size matching, a difference between arms could be informativeness rather
    than falsity.
    """
    cpdag, true_graph = mpdag(chain, as_cpdag=True), mpdag(chain)
    with pytest.raises(NotImplementedError):
        false_bk = perturb.sample_consistent_but_false(
            true_graph, cpdag, delta=1, kind=MisspecificationKind.ORIENTATION, rng=rng
        )
        true_bk = perturb.sample_consistent_and_true(
            true_graph, cpdag, len(false_bk), MisspecificationKind.ORIENTATION, rng
        )
        assert len(true_bk) == len(false_bk)


def test_enumeration_agrees_with_sampling(mpdag, chain, rng):
    """Sampled knowledge sets are members of the enumerated ball."""
    cpdag, true_graph = mpdag(chain, as_cpdag=True), mpdag(chain)
    with pytest.raises(NotImplementedError):
        ball = list(
            perturb.enumerate_consistent_but_false(
                true_graph, cpdag, delta=1, kind=MisspecificationKind.ORIENTATION
            )
        )
        bk = perturb.sample_consistent_but_false(
            true_graph, cpdag, delta=1, kind=MisspecificationKind.ORIENTATION, rng=rng
        )
        assert bk in ball
