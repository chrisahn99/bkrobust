"""Guards for the knowledge-selection modes and the three-valued verdict.

The property that matters and that no other test covers: a coverage sweep is
only a retraction sequence if the smaller knowledge set is a subset of the
larger one. The original even stride does not guarantee that, and on this
corpus it fails on four networks, so a sweep over it measures a change of
knowledge rather than a loss of it. The nested mode guarantees it and agrees
with the stride at full coverage, which is what lets the committed rows stand.

These run offline against the cached corpus and skip cleanly without it.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from bkrobust.benchmarks.describe import parse_file, to_mpdag
from bkrobust.benchmarks.measure import select_knowledge
from bkrobust.demo.example import dag_to_cpdag

CACHE = Path("results/axisa3/networks")

#: Small enough to build a CPDAG for inside the test budget.
SMALL = ("asia", "child", "sachs", "insurance")
#: Measured on the committed corpus: the stride keeps indices at 0.25 that are
#: not among those it keeps at 0.5 on these four networks.
STRIDE_NOT_NESTED = ("diabetes", "ecoli70", "magic-irri", "munin1")
COVERAGES = (1.0, 0.5, 0.25)


@pytest.fixture(scope="module")
def graphs():
    """The small networks as ``(dag, cpdag)`` pairs, from the cached corpus."""
    if not CACHE.exists():
        pytest.skip("benchmark cache unavailable")
    out = {}
    for name in SMALL:
        match = [p for p in sorted(CACHE.glob("example_models/*")) if p.name.startswith(name)]
        if not match:
            continue
        digest = hashlib.sha256(match[0].read_bytes()).hexdigest()
        parsed = parse_file(match[0], digest)
        dag = to_mpdag(parsed)
        out[name] = (dag, dag_to_cpdag(dag))
    if not out:
        pytest.skip("no small network found in the cache")
    return out


def test_modes_agree_at_full_coverage(graphs) -> None:
    """Nothing is selected away at coverage 1.0, so the modes cannot differ there."""
    for dag, cpdag in graphs.values():
        stride = select_knowledge(dag, cpdag, 1.0, mode="stride")
        nested = select_knowledge(dag, cpdag, 1.0, mode="nested")
        assert stride == nested


def test_nested_mode_is_nested(graphs) -> None:
    """Lower coverage is a subset of higher coverage, which is what makes a sweep a retraction."""
    for dag, cpdag in graphs.values():
        sets = {c: set(select_knowledge(dag, cpdag, c, mode="nested")) for c in COVERAGES}
        assert sets[0.25] <= sets[0.5] <= sets[1.0]


def test_nested_mode_is_deterministic(graphs) -> None:
    """Same seed, same selection: no global RNG anywhere on this path."""
    for dag, cpdag in graphs.values():
        first = select_knowledge(dag, cpdag, 0.5, mode="nested", seed=7)
        second = select_knowledge(dag, cpdag, 0.5, mode="nested", seed=7)
        assert first == second


def test_sizes_match_between_modes(graphs) -> None:
    """The modes select the same number of claims; only which ones differs."""
    for dag, cpdag in graphs.values():
        for c in COVERAGES:
            stride = select_knowledge(dag, cpdag, c, mode="stride")
            nested = select_knowledge(dag, cpdag, c, mode="nested")
            assert len(stride) == len(nested)


def test_unknown_mode_is_refused(graphs) -> None:
    """A typo in the mode must not silently fall back to the committed behaviour."""
    dag, cpdag = next(iter(graphs.values()))
    with pytest.raises(ValueError):
        select_knowledge(dag, cpdag, 0.5, mode="evenly-spaced")
