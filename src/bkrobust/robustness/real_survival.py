"""Survival under corrupted background knowledge, on real network structure ([RE-11]).

The synthetic counterpart is :mod:`bkrobust.robustness.survival`, which is left
**unmodified** so every number in ``report_fragility_and_pareto.md`` stays
reproducible from the code that produced it (the precedent is commit
``748f2e1``). This module is its real-structure sibling and shares none of its
files.

Design is fixed by ``results/axis_robustness_real/PREREGISTRATION.md``; nothing
here may deviate from it without a dated appendix in that file.

The one structural difference from synthetic, and it drives everything
---------------------------------------------------------------------
On the real corpus the analyst's knowledge is a property of the **network and
the coverage level**, not of the query: ``benchmarks.measure.select_knowledge``
returns one claim set per ``(network, coverage)`` and every admissible
``(X, Y)`` pair of that network reads its adjustment set off the *same* ``G0``.
So a corrupted knowledge state is shared by all of a network's pairs, and the
Monte-Carlo draw is taken **once per network per grid point** and scored against
every pair.

Three consequences, all recorded rather than smoothed over:

* ``|K|`` cannot vary within a ``(network, coverage)`` cell -- ``flip``
  preserves list length -- so the ``n_k`` baseline can rank networks but not
  queries on this corpus.
* The contradiction rate at a grid point is a property of the network, not of
  the pair, and is therefore identical across a network's rows.
* Survival outcomes for two pairs of one network at one ``(grid point, rep)``
  are perfectly dependent. That is exactly why the pre-registration bootstraps
  over **networks** and prints no per-row interval.

Gate discipline
---------------
``benchmarks.measure.fast_gate`` is the only gate, and it is inherited: every
row of the frame was admitted by it in ``results/axisa3/instances.jsonl``. Every
row this module writes carries the literal string ``"fast_gate"`` in a ``gate``
column. ``synth.runner.gate`` is never imported -- it calls the exponential
``all_valid_adjustment_sets_mpdag`` -- and neither is
``all_valid_adjustment_sets_mpdag`` itself, on any path.

Sentinel discipline
-------------------
``UNREACHED = -1`` is a status, never a radius: it is never averaged, never
plotted numerically, never fed to a correlation. A wall-capped computation
carries ``wall_until_timeout_s`` and **never** a key a real measurement uses --
the ``{"timed_out_at_s": CAP, "seconds": <wall clock>}`` shape that let a
timeout masquerade as a datum in session 4 cannot occur here because no row
this module writes has a key named ``seconds``. A contradictory corrupted ``K``
is status ``corrupted_k_contradictory`` and is **not** a survival-0 outcome;
``S`` is computed over non-contradictory draws and the contradiction rate is its
own series.

RNG discipline
--------------
No global RNG is ever touched. Every draw goes through an
``np.random.Generator`` built by :func:`derived_seed` -- ``hashlib.sha256`` over
the parts, never Python's salted built-in ``hash()`` -- so results are identical
across ``PYTHONHASHSEED`` values. That is checked by the driver, not assumed.

The closure cache
-----------------
``apply_orientations`` is a pure function of ``(cpdag, k_cor)`` and the GAC
predicate is a pure function of ``(g, x, y, Z*)``. Within one shard, draws that
happen to produce the *same* corrupted claim set are memoised. This changes no
statistic -- every draw is still drawn and still counted, only the graph
computation is reused -- and it matters because at small ``|K|`` the number of
distinct corruptions at a given depth is tiny (``|K| = 1`` at depth 1 admits
exactly one). Cache hits are recorded per cell so the saving is auditable.
"""

from __future__ import annotations

import hashlib
from typing import Any, Iterable, Sequence

import numpy as np

from bkrobust.benchmarks import describe as _describe
from bkrobust.benchmarks.acquire import load_cached
from bkrobust.benchmarks.measure import (
    MAX_G0_UNDIRECTED_FOR_EXTENSIONS,
    O_INTRACTABLE,
    select_knowledge,
)
from bkrobust.core.conventions import UNREACHED
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations
from bkrobust.gac import is_gac_valid_mpdag
from bkrobust.hybrid import breakdown_radius
from bkrobust.synth.knowledge import flip, tiered

Edge = tuple[str, str]

#: Literal gate name stamped on every output row.
GATE_NAME = "fast_gate"

#: The fixed fractional grid the ``AUC_frac`` endpoint averages over
#: (``results/axis_robustness/PREREGISTRATION.md`` §4).
FRAC_GRID: tuple[float, ...] = tuple(round(i / 10, 1) for i in range(1, 11))

#: The tiered arm's corruption-rate grid, identical to the synthetic sweep's.
RATE_GRID: tuple[float, ...] = tuple(round(0.05 * i, 2) for i in range(0, 11))

#: The cross-arm's tiered relocation grid (rate 0 carries no comparison).
XARM_TIERED_RATE_GRID: tuple[float, ...] = tuple(round(0.05 * i, 2) for i in range(1, 11))

#: Minimum non-contradictory draws a grid point needs to enter a ``_usable``
#: endpoint or a paired cross-arm bin.
MIN_USABLE_N = 30

#: Shared-intensity bin width for the cross-arm.
BIN_WIDTH = 0.05

#: Wall cap on a single radius computation, matching the committed corpus.
RADIUS_TIME_LIMIT_S = 300.0

#: Pre-registered draw count.
N_DRAWS_DEFAULT = 1000

_OK = "ok"


# --- seeding -----------------------------------------------------------------


def derived_seed(*parts: Any) -> int:
    """A stable, ``PYTHONHASHSEED``-independent seed derived from ``parts``.

    Uses ``hashlib.sha256`` over the string form of ``parts`` -- never Python's
    built-in ``hash()``, which is salted per process.

    Args:
        *parts: Anything with a stable ``str()``. Order matters.

    Returns:
        A nonnegative integer suitable for ``np.random.default_rng``.
    """
    payload = "|".join(str(p) for p in parts).encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big")


def sha_of(obj: Any) -> str:
    """A short, stable SHA-256 digest of ``obj``'s sorted string form.

    Args:
        obj: Anything with a stable ``str()`` after sorting, e.g. a sorted list
            of edges.

    Returns:
        The first 16 hex characters of the digest.
    """
    return hashlib.sha256(repr(obj).encode("utf-8")).hexdigest()[:16]


# --- networks ----------------------------------------------------------------

_NET_CACHE: dict[str, tuple[MPDAG, MPDAG]] = {}


def network_filename(name: str, available: Iterable[str]) -> str:
    """The corpus filename for a network, by the repository's naming rule.

    Args:
        name: Network name as it appears in ``results/axisa3/instances.jsonl``.
        available: The filenames present in the acquisition cache.

    Returns:
        The matching filename.

    Raises:
        KeyError: If no candidate filename is present.
    """
    have = set(available)
    for candidate in (f"{name}.bif.gz", f"{name}.json", f"{name}.txt"):
        if candidate in have:
            return candidate
    raise KeyError(f"no corpus file for network {name!r}")


def load_network(name: str) -> tuple[MPDAG, MPDAG]:
    """Parse one network and build its oracle CPDAG, memoised per process.

    Args:
        name: Network name.

    Returns:
        ``(dag, cpdag)``.
    """
    if name not in _NET_CACHE:
        acq = load_cached()
        recs = {r.name: r for r in acq.files}
        fname = network_filename(name, recs)
        parsed = _describe.parse_file(recs[fname].path, recs[fname].sha256)
        dag = _describe.to_mpdag(parsed)
        _NET_CACHE[name] = (dag, dag_to_cpdag(dag))
    return _NET_CACHE[name]


# --- predictors ---------------------------------------------------------------


def directed_symdiff(g1: MPDAG, g2: MPDAG) -> int:
    """``|dir(g1) Δ dir(g2)|``, the directed-edge symmetric difference.

    Recorded on sampled states as ``symdiff_proxy_not_distance``. It is an
    upper-bound surrogate (``THEOREMS.md`` §10-§11), never a distance.
    """
    return len(g1.directed_edges ^ g2.directed_edges)


def commitment_size(g0: MPDAG, cpdag: MPDAG) -> int:
    """``k_g0`` -- orientations ``g0`` commits to that ``cpdag`` left undirected.

    Args:
        g0: The analyst's Meek-closed graph.
        cpdag: The oracle CPDAG.

    Returns:
        The commitment size. This is **not** "the knowledge size"; the stated
        size is ``|K|`` (``docs/PAPER_NARRATIVE.md`` §11).
    """
    und = cpdag.undirected_edges
    return sum(1 for (a, b) in g0.directed_edges if tuple(sorted((a, b))) in und)


# --- radius -------------------------------------------------------------------


def radius_columns(
    cpdag: MPDAG, g0: MPDAG, x: str, y: str, z_star: Iterable[str],
    *, time_limit_s: float = RADIUS_TIME_LIMIT_S,
) -> dict[str, Any]:
    """Compute the breakdown radius and flatten it, with sentinel discipline.

    A timed-out ladder leaves ``radius`` ``None`` -- never a number -- and
    stamps ``wall_until_timeout_s``. ``UNREACHED`` is stored as the sentinel
    integer only when the search *completed* and genuinely found no failure, and
    is flagged by ``r_status`` so it can never be silently averaged.

    Args:
        cpdag: The oracle CPDAG.
        g0: The analyst's graph.
        x: Treatment.
        y: Outcome.
        z_star: The committed adjustment set.
        time_limit_s: Wall cap on the search.

    Returns:
        A flat dict of columns.
    """
    out = breakdown_radius(cpdag, None, x, y, frozenset(z_star), g0=g0, time_limit_s=time_limit_s)
    row: dict[str, Any] = {
        "dispatch_leg": out.method,
        "oracle": out.oracle,
        "exact": out.exact,
        "assumes": out.assumes,
        "r_search_seconds": out.search_seconds,
        "r_total_seconds": out.total_seconds,
    }
    if not out.exact:
        row["radius"] = None
        row["r_status"] = "censored_wall_cap"
        row["wall_until_timeout_s"] = out.total_seconds
    elif out.radius == UNREACHED:
        row["radius"] = UNREACHED
        row["r_status"] = "unreached"
        row["wall_until_timeout_s"] = None
    else:
        row["radius"] = out.radius
        row["r_status"] = _OK
        row["wall_until_timeout_s"] = None
    return row


# --- grids --------------------------------------------------------------------


def depth_grid(n_k: int) -> tuple[int, ...]:
    """The claim depths the flip arm sweeps, for an instance with ``n_k`` claims.

    Exactly the depths the pre-registered ``AUC_frac`` endpoint averages over:
    ``clamp(round(f · n_k), 1, n_k)`` for ``f`` on :data:`FRAC_GRID`,
    deduplicated. For ``n_k <= 10`` this is the dense sweep ``1 … n_k``; above
    it, the ten fractional grid points. Nothing is measured that the endpoint
    discards and nothing the endpoint needs is missing.

    Args:
        n_k: The stated knowledge size.

    Returns:
        The ascending, deduplicated depth grid.

    Raises:
        ValueError: If ``n_k < 1``.
    """
    if n_k < 1:
        raise ValueError(f"n_k must be >= 1, got {n_k}")
    return tuple(sorted({min(n_k, max(1, round(f * n_k))) for f in FRAC_GRID}))


# --- the committed state ------------------------------------------------------


def build_g0(
    cpdag: MPDAG, k: Sequence[Edge]
) -> tuple[MPDAG | None, str]:
    """Meek-close ``k`` onto ``cpdag`` and screen the result for tractability.

    Args:
        cpdag: The oracle CPDAG.
        k: The analyst's orientation claims.

    Returns:
        ``(g0, "ok")``, or ``(None, reason)`` with reason
        ``"k_contradictory"``, ``"n_k_zero"`` or
        ``"o_g0_extensions_intractable"``. The last is a **measurement limit**,
        deliberately kept distinct from every structural reason.
    """
    if len(k) == 0:
        return None, "n_k_zero"
    g0 = apply_orientations(cpdag, list(k))
    if g0 is None:
        return None, "k_contradictory"
    if len(g0.undirected_edges) > MAX_G0_UNDIRECTED_FOR_EXTENSIONS:
        return None, O_INTRACTABLE
    return g0, _OK


def commit_z_star(g0: MPDAG, x: str, y: str) -> tuple[frozenset[str] | None, str]:
    """Read the adjustment set off ``g0`` once, and check it is valid there.

    Args:
        g0: The analyst's graph.
        x: Treatment.
        y: Outcome.

    Returns:
        ``(z_star, "ok")``, or ``(None, reason)``.
    """
    z = optimal_adjustment_set_mpdag(g0, x, y)
    if not z:
        return None, "optimal_set_undefined"
    z = frozenset(z)
    if not is_gac_valid_mpdag(g0, x, y, z):
        return None, "z_invalid_at_g0"
    return z, _OK


# --- corruption draws ---------------------------------------------------------


class ClosureCache:
    """Memoises ``apply_orientations`` and the GAC verdict within one shard.

    Both memoised functions are pure, so this changes no statistic: every draw
    is still drawn and still counted, only the graph computation is reused.

    Args:
        cpdag: The oracle CPDAG every corrupted state is closed onto.
    """

    def __init__(self, cpdag: MPDAG) -> None:
        self._cpdag = cpdag
        self._graphs: dict[tuple[Edge, ...], MPDAG | None] = {}
        self._verdicts: dict[tuple[tuple[Edge, ...], str, str, frozenset], bool] = {}
        self.graph_hits = 0
        self.graph_misses = 0
        self.verdict_hits = 0
        self.verdict_misses = 0

    def graph(self, k_cor: Sequence[Edge]) -> tuple[tuple[Edge, ...], MPDAG | None]:
        """The Meek closure of ``k_cor``, and the key it was cached under.

        Args:
            k_cor: A corrupted claim set.

        Returns:
            ``(key, graph_or_None)``. ``None`` means ``k_cor`` admits no
            consistent MPDAG -- a contradiction, not a failure.
        """
        key = tuple(sorted(k_cor))
        if key in self._graphs:
            self.graph_hits += 1
        else:
            self.graph_misses += 1
            self._graphs[key] = apply_orientations(self._cpdag, list(key))
        return key, self._graphs[key]

    def survived(
        self, key: tuple[Edge, ...], g: MPDAG, x: str, y: str, z_star: frozenset
    ) -> bool:
        """Whether ``z_star`` is still GAC-valid for ``(x, y)`` at ``g``.

        Args:
            key: The cache key :meth:`graph` returned for ``g``.
            g: The corrupted state.
            x: Treatment.
            y: Outcome.
            z_star: The committed adjustment set, held fixed since depth 0.

        Returns:
            The verdict of the polynomial GAC predicate.
        """
        vkey = (key, x, y, z_star)
        if vkey in self._verdicts:
            self.verdict_hits += 1
        else:
            self.verdict_misses += 1
            self._verdicts[vkey] = is_gac_valid_mpdag(g, x, y, z_star)
        return self._verdicts[vkey]


def flip_draw(
    k_b: Sequence[Edge], d: int, rep: int, *, family: str
) -> tuple[list[Edge], int]:
    """One flip-arm draw: reverse exactly ``d`` of ``k_b``'s claims.

    Args:
        k_b: The analyst's (possibly base-wrongness-corrupted) claim set.
        d: Target depth, ``1 <= d <= len(k_b)``.
        rep: Repetition index.
        family: The instance-family id the seed is derived from.

    Returns:
        ``(k_cor, seed)``.

    Raises:
        AssertionError: If the realised number of reversals is not ``d``. The
            realised count is asserted, never assumed.
    """
    n = len(k_b)
    assert 1 <= d <= n, f"d={d} out of range for n_k={n}"
    seed = derived_seed(family, "flip", d, rep)
    rng = np.random.default_rng(seed)
    k_cor = flip(list(k_b), rng, d / n)
    realised = sum(1 for e in k_cor if e not in set(k_b))
    assert realised == d, f"flip realised {realised} reversals, wanted {d}"
    return k_cor, seed


def tiered_draw(
    dag: MPDAG, cpdag: MPDAG, n_tiers: int, rate: float, rep: int, *, family: str
) -> tuple[list[Edge], int]:
    """One tiered-arm draw: relocate ``rate`` of the nodes, re-derive the claims.

    Args:
        dag: The ground-truth DAG, read only for the depth-based tier order.
        cpdag: The oracle CPDAG.
        n_tiers: Number of tiers.
        rate: Fraction of nodes misfiled into a random other tier.
        rep: Repetition index.
        family: The instance-family id the seed is derived from.

    Returns:
        ``(k_cor, seed)``.
    """
    seed = derived_seed(family, "tiered", rate, rep)
    rng = np.random.default_rng(seed)
    return tiered(dag, cpdag, rng, n_tiers, rate), seed


# --- cells --------------------------------------------------------------------


def cell_statistics(
    n_draws: int, n_contradictory: int, n_survived: int
) -> dict[str, Any]:
    """Aggregate one ``(instance, grid point)`` cell.

    ``S`` is undefined -- ``None``, never a number and never zero -- when every
    draw at the grid point was contradictory.

    Args:
        n_draws: Draws taken at this grid point.
        n_contradictory: Draws whose corrupted ``K`` admitted no MPDAG.
        n_survived: Non-contradictory draws at which ``Z*`` stayed valid.

    Returns:
        ``n_eval``, ``contradiction_rate``, ``S`` and ``S_contra_as_fail``.
    """
    n_eval = n_draws - n_contradictory
    return {
        "n_draws": n_draws,
        "n_eval": n_eval,
        "n_contradictory": n_contradictory,
        "n_survived": n_survived,
        "contradiction_rate": (n_contradictory / n_draws) if n_draws else None,
        "S": (n_survived / n_eval) if n_eval else None,
        "S_contra_as_fail": (n_survived / n_draws) if n_draws else None,
    }


def auc_over_grid(cells: dict[Any, dict[str, Any]], targets: Sequence[Any]) -> float | None:
    """Mean ``S`` over ``targets``, taking the nearest grid point with a defined ``S``.

    Ties break toward the smaller grid point, matching
    :func:`bkrobust.robustness.survival.auc_frac` exactly.

    Args:
        cells: ``{grid_point: cell_statistics(...)}``.
        targets: The grid points the endpoint averages over.

    Returns:
        The mean, or ``None`` if no grid point has a defined ``S``.
    """
    defined = sorted(g for g, c in cells.items() if c["S"] is not None)
    if not defined:
        return None
    vals = [cells[min(defined, key=lambda g: (abs(g - t), g))]["S"] for t in targets]
    return sum(vals) / len(vals)


def auc_usable(
    cells: dict[Any, dict[str, Any]],
    targets: Sequence[Any] | None = None,
    *,
    min_n: int = MIN_USABLE_N,
) -> float | None:
    """The conservative ``_usable`` endpoint: the raw endpoint, recomputed after
    dropping grid points backed by fewer than ``min_n`` evaluable draws.

    With ``targets`` this is :func:`auc_over_grid` on the filtered cells, which
    is exactly how the synthetic analysis defines ``AUC_frac_usable``
    (``robustness/analyse.py`` recomputes ``auc_frac`` on the restricted curve,
    it does not take a plain mean). Without ``targets`` -- the tiered arm, whose
    grid is fixed and identical across instances, so no nearest-point mapping is
    needed -- it is the plain mean over the surviving grid points, which is how
    ``AUC_rate_usable`` is defined.

    It is reported **beside** the raw endpoint, never instead of it: Appendix J
    of the synthetic pre-registration established that it conditions on
    ``n_eval``, which is correlated with the radius, so it is not a conservative
    control at low draw counts. The gap between the two is the diagnostic.

    Args:
        cells: ``{grid_point: cell_statistics(...)}``.
        targets: The grid points the endpoint averages over, for the flip arm.
        min_n: The evaluable-draw threshold.

    Returns:
        The endpoint, or ``None`` if no grid point qualifies.
    """
    kept = {g: c for g, c in cells.items() if c["n_eval"] >= min_n and c["S"] is not None}
    if not kept:
        return None
    if targets is not None:
        return auc_over_grid(kept, targets)
    return sum(c["S"] for c in kept.values()) / len(kept)


def intensity_bin(intensity: float) -> float:
    """Left edge of the :data:`BIN_WIDTH`-wide bin ``intensity`` falls into.

    No ceiling: a heavily corrupted state can disagree with ``G0`` on more
    orientations than ``G0`` had, and that is reported rather than clipped.

    Args:
        intensity: ``|dir(G0) Δ dir(G)| / |dir(G0)|``.

    Returns:
        The bin's left edge.
    """
    import math

    return round(math.floor(round(intensity / BIN_WIDTH, 9)) * BIN_WIDTH, 2)


__all__ = [
    "GATE_NAME", "FRAC_GRID", "RATE_GRID", "XARM_TIERED_RATE_GRID",
    "MIN_USABLE_N", "BIN_WIDTH", "RADIUS_TIME_LIMIT_S", "N_DRAWS_DEFAULT",
    "derived_seed", "sha_of", "network_filename", "load_network",
    "directed_symdiff", "commitment_size", "radius_columns", "depth_grid",
    "build_g0", "commit_z_star", "ClosureCache", "flip_draw", "tiered_draw",
    "cell_statistics", "auc_over_grid", "auc_usable", "intensity_bin",
]
