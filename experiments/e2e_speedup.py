"""One end-to-end speedup figure: space -> search -> hybrid, on ONE sweep.

Section 5.3 needs a single composable speedup number. Two prior sweeps measured
the two legs separately and cannot be composed:

* ``results/axisb2/scaling/scaling.csv`` -- full-space BFS vs
  :func:`bkrobust.search.exact.radius_local_up`, both under the enumeration
  oracle, at n=5..7.
* ``results/axisb4/hybrid_envelope.jsonl`` -- hybrid
  (:func:`bkrobust.hybrid.breakdown_radius`) vs the frozen search, at n=8..24,
  with the enumeration oracle. Its generating script was never committed.

This module measures every leg -- full-space BFS, the frozen search, the
criterion-accelerated search, and the hybrid dispatcher -- on the SAME 256
instances and the SAME ``(x, y, Z)`` query, in one sweep, so the two ratios
compose into one exact identity:

    speedup_e2e = speedup_space_vs_search * speedup_search_vs_hybrid

Query-rule reconstruction
--------------------------

The committed ``hybrid_envelope.jsonl`` records ``x``, ``y`` and ``z_size`` per
``(n, edge_prob, seed)`` but not the rule that produced them, and its
generating script was never committed either. This module first tries to
reproduce it: build ``G0`` via :func:`bkrobust.gac.sweep.dense_instance` with
the literal ``(n, seed, edge_prob)`` fields from the committed file, and take
``Z = optimal_adjustment_set_mpdag(g0, x, y)`` using the committed ``x, y``.

That reconstruction was checked against all 256 committed rows before this
module was written (see the smoke-test log referenced in ``manifest.json``)
and it reproduces ``z_size`` on only 6/256 rows -- consistent with noise, not
a matching rule. A stronger, graph-topology-only check confirms the instance
itself differs: for row 0 (n=8, edge_prob=0.3, seed=0), the committed query is
``x=V0, y=V3, z_size=3``, but ``dense_instance(8, 0, 0.3)`` places node ``V0``
with zero edges in the underlying DAG (it is an isolated node under
``erdos_renyi_dag``'s topological-order construction), which makes any
adjustment set of size 3 for a query rooted at ``V0`` impossible under any
valid definition of validity -- there are no backdoor paths through an
isolated node. So the committed envelope was not built by literally calling
``dense_instance(n, seed, edge_prob)`` with the fields it stored; the mapping
from those fields to the actual generator seed (or generator itself) used in
session 4 could not be recovered from the repository.

Given that, this driver falls back to the conservative option the revision
brief anticipates: it builds its own 256-instance grid from
``dense_instance(n, seed, edge_prob)`` with the literal grid values below, and
picks its own deterministic, documented query per instance -- the first
``(x, y)`` pair (sorted-node permutation order, matching the iteration order
in ``bkrobust.search.scaling.measure_instance``) whose
``optimal_adjustment_set_mpdag`` is defined, non-empty, and valid under the
enumeration oracle. This keeps every leg on the same instances and the same
query within THIS sweep, which is what the composed speedup number needs; it
does not attempt to match the old envelope's specific rows.

Resample gate (``--gate``)
---------------------------

Default behaviour (``--gate`` unset) is unchanged from the committed
``results/e2e_speedup`` run: attempt 0, i.e. the literal
``dense_instance(n, seed, edge_prob)`` graph, is used for every cell, and a
cell with no valid query is skipped (not written), exactly as before -- so the
committed run stays reproducible with its recorded command line.

``--gate`` turns on a deterministic resample gate whose OUTPUT IS NOT the
session-4 envelope. It is a k>=1-gated draw of the SAME grid: for each grid
cell ``(n, edge_prob, seed)`` it tries attempt ``a = 0, 1, 2, ..., 100``.
Attempt 0 uses ``effective_seed = seed`` (so any cell that already passes
keeps its attempt-0 instance, bit-identical to the ungated run). Attempt
``a >= 1`` uses ``effective_seed = seed + 10_000 * a`` -- no other source of
randomness is touched, and no global RNG is mutated. The first attempt whose
CPDAG has ``k_undirected >= 1`` AND whose ``choose_query`` finds a query is
accepted; ``attempt`` and ``effective_seed`` are recorded on every row (jsonl
and csv) regardless of whether the gate is enabled. If no attempt in
``0..100`` satisfies the gate, the cell is recorded as a ``gate_exhausted``
error row -- it is never silently dropped. Worker subprocesses rebuild the
instance from ``(n, effective_seed, edge_prob)`` so they are bit-identical to
what the driver selected.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import os
import platform
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"

# --------------------------------------------------------------------------
# Grid and constants
# --------------------------------------------------------------------------

N_GRID: tuple[int, ...] = (8, 10, 12, 14, 16, 18, 20, 24)
EDGE_PROB_GRID: tuple[float, ...] = (0.3, 0.5, 0.7, 0.85)
SEED_GRID: tuple[int, ...] = tuple(range(8))

CAP_S = 120.0
LEGS = ("space", "search_frozen", "search_fast", "hybrid")
GATE_MAX_ATTEMPTS = 100  # attempts 0..GATE_MAX_ATTEMPTS inclusive


def instance_grid() -> list[tuple[int, float, int]]:
    """The full 256-cell grid, in deterministic (n, edge_prob, seed) order."""
    return [
        (n, ep, seed)
        for n in N_GRID
        for ep in EDGE_PROB_GRID
        for seed in SEED_GRID
    ]


def instance_key(n: int, edge_prob: float, seed: int) -> str:
    return f"n{n}_ep{edge_prob:g}_seed{seed}"


# --------------------------------------------------------------------------
# Instance construction (shared by parent process and worker subprocess so
# both build bit-identical graphs from (n, seed, edge_prob) alone).
# --------------------------------------------------------------------------


def build_cpdag_and_g0(n: int, seed: int, edge_prob: float):
    """Re-derive the CPDAG and G0 with the identical calls ``dense_instance`` makes.

    ``dense_instance`` itself returns only ``G0``; the space-building legs need
    the CPDAG too, so this reproduces the same sequence of calls rather than
    importing a private helper, to keep both bit-identical to what
    :func:`bkrobust.gac.sweep.dense_instance` returns.
    """
    from bkrobust.demo.example import dag_to_cpdag, knowledge_to_recover
    from bkrobust.demo.meek import apply_orientations
    from bkrobust.synth.generators import erdos_renyi_dag

    rng = np.random.default_rng(seed)
    dag = erdos_renyi_dag(n, rng, edge_prob)
    cpdag = dag_to_cpdag(dag)
    knowledge = sorted(knowledge_to_recover(dag, cpdag))
    g0 = apply_orientations(cpdag, knowledge)
    if g0 is None:
        raise RuntimeError(f"orientations FAILed on n={n}, seed={seed}, p={edge_prob}")
    return cpdag, g0


def choose_query(cpdag, g0):
    """First valid, non-empty, sorted-permutation-order query on ``g0``.

    Mirrors the iteration order in
    :func:`bkrobust.search.scaling.measure_instance`: sorted nodes,
    ``itertools.permutations(nodes, 2)``, first pair with a defined, non-empty,
    valid optimal adjustment set under the enumeration oracle.

    Returns:
        ``(x, y, z)`` or ``None`` if no instance in this CPDAG admits one.
    """
    from bkrobust.core.oracle import is_valid
    from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag

    nodes = sorted(g0.nodes)
    for x, y in itertools.permutations(nodes, 2):
        o = optimal_adjustment_set_mpdag(g0, x, y)
        if o is None:
            continue
        z = frozenset(o)
        if not z or not is_valid(z, g0, x, y):
            continue
        return x, y, z
    return None


def graph_stats(cpdag, g0) -> dict[str, Any]:
    """Topology fields independent of the query: k_undirected, k_g0, largest_component."""
    import networkx as nx

    from bkrobust.search.exact import retractable_edges

    k_undirected = len(cpdag.undirected_edges)
    k_g0 = len(retractable_edges(cpdag, g0))

    graph = nx.Graph()
    graph.add_nodes_from(cpdag.nodes)
    graph.add_edges_from(cpdag.undirected_edges)
    comps = [len(c) for c in nx.connected_components(graph) if len(c) > 1]
    largest_component = max(comps) if comps else 0

    return {"k_undirected": k_undirected, "k_g0": k_g0, "largest_component": largest_component}


# --------------------------------------------------------------------------
# Worker: runs exactly one leg for one instance, in a fresh subprocess.
# --------------------------------------------------------------------------


def _worker_space(n: int, seed: int, edge_prob: float, x: str, y: str, z: list[str]) -> dict[str, Any]:
    from bkrobust.core.oracle import clear_cache as oracle_clear_cache
    from bkrobust.core.oracle import is_valid
    from bkrobust.core.spacelib import distances_from, radius
    from bkrobust.mpdag_criterion.criterion import clear_cache as crit_clear_cache
    from bkrobust.mpdag_criterion.criterion import is_valid_mpdag
    from bkrobust.search.space_fixed import build_space_correct

    cpdag, g0 = build_cpdag_and_g0(n, seed, edge_prob)
    zf = frozenset(z)

    def fails_enum(g):
        return not is_valid(zf, g, x, y)

    def fails_crit(g):
        return not is_valid_mpdag(g, x, y, zf)

    oracle_clear_cache()
    crit_clear_cache()
    t0 = time.perf_counter()
    space = build_space_correct(cpdag)
    build_s = time.perf_counter() - t0

    dists = distances_from(space, g0)

    t0 = time.perf_counter()
    r_enum, _w1 = radius(space, dists, fails_enum)
    bfs_query_enum_s = time.perf_counter() - t0

    t0 = time.perf_counter()
    r_crit, _w2 = radius(space, dists, fails_crit)
    bfs_query_crit_s = time.perf_counter() - t0

    return {
        "build_s": build_s,
        "bfs_query_enum_s": bfs_query_enum_s,
        "bfs_query_crit_s": bfs_query_crit_s,
        "space_size": len(space),
        "r_space_enum": r_enum,
        "r_space_crit": r_crit,
    }


def _worker_search_frozen(n: int, seed: int, edge_prob: float, x: str, y: str, z: list[str]) -> dict[str, Any]:
    from bkrobust.core.oracle import clear_cache as oracle_clear_cache
    from bkrobust.core.oracle import is_valid
    from bkrobust.search.exact import SearchStats, radius_local_up

    cpdag, g0 = build_cpdag_and_g0(n, seed, edge_prob)
    zf = frozenset(z)

    def fails_enum(g):
        return not is_valid(zf, g, x, y)

    oracle_clear_cache()
    st = SearchStats()
    t0 = time.perf_counter()
    res = radius_local_up(cpdag, g0, fails_enum, stats=st)
    search_s = time.perf_counter() - t0

    return {
        "search_s": search_s,
        "radius": res.radius,
        "exact": res.exact,
        "method": res.method,
        **st.as_dict(),
    }


def _worker_search_fast(n: int, seed: int, edge_prob: float, x: str, y: str, z: list[str]) -> dict[str, Any]:
    from bkrobust.mpdag_criterion.criterion import clear_cache as crit_clear_cache
    from bkrobust.mpdag_criterion.criterion import is_valid_mpdag
    from bkrobust.search.exact import SearchStats
    from bkrobust.search.exact_fast import radius_local_up_fast

    cpdag, g0 = build_cpdag_and_g0(n, seed, edge_prob)
    zf = frozenset(z)

    def fails_crit(g):
        return not is_valid_mpdag(g, x, y, zf)

    crit_clear_cache()
    st = SearchStats()
    t0 = time.perf_counter()
    res = radius_local_up_fast(cpdag, g0, fails_crit, stats=st)
    search_s = time.perf_counter() - t0

    return {
        "search_s": search_s,
        "radius": res.radius,
        "exact": res.exact,
        "method": res.method,
        **st.as_dict(),
    }


def _worker_hybrid(n: int, seed: int, edge_prob: float, x: str, y: str, z: list[str]) -> dict[str, Any]:
    from bkrobust.mpdag_criterion.criterion import clear_cache as crit_clear_cache
    from bkrobust.hybrid import breakdown_radius

    cpdag, g0 = build_cpdag_and_g0(n, seed, edge_prob)
    zf = frozenset(z)

    crit_clear_cache()
    res = breakdown_radius(cpdag, None, x, y, zf, g0=g0)

    return {
        "radius": res.radius,
        "method": res.method,
        "oracle": res.oracle,
        "exact": res.exact,
        "search_seconds": res.search_seconds,
        "ladder_seconds": res.ladder_seconds,
        "total_seconds": res.total_seconds,
    }


_WORKERS = {
    "space": _worker_space,
    "search_frozen": _worker_search_frozen,
    "search_fast": _worker_search_fast,
    "hybrid": _worker_hybrid,
}


def run_worker(leg: str, payload: dict[str, Any]) -> None:
    """Entry point invoked in the fresh child subprocess. Prints one JSON line."""
    fn = _WORKERS[leg]
    out = fn(payload["n"], payload["seed"], payload["edge_prob"], payload["x"], payload["y"], payload["z"])
    print(json.dumps(out))


# --------------------------------------------------------------------------
# Parent-side subprocess dispatch, with the 120s cap and censoring discipline.
# --------------------------------------------------------------------------


def run_leg(leg: str, n: int, seed: int, edge_prob: float, x: str, y: str, z: list[str], cap_s: float) -> dict[str, Any]:
    payload = {"n": n, "seed": seed, "edge_prob": edge_prob, "x": x, "y": y, "z": z}
    cmd = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--worker",
        leg,
        "--worker-payload",
        json.dumps(payload),
    ]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(SRC)

    wall0 = time.perf_counter()
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=cap_s,
            env=env,
            cwd=str(ROOT),
        )
    except subprocess.TimeoutExpired:
        wall = time.perf_counter() - wall0
        return {"censored": True, "wall_until_timeout_s": wall}

    if proc.returncode != 0:
        return {"error": (proc.stderr or "").strip()[-4000:], "returncode": proc.returncode}

    lines = [ln for ln in proc.stdout.strip().splitlines() if ln.strip()]
    if not lines:
        return {"error": "worker produced no stdout", "stderr": (proc.stderr or "").strip()[-4000:]}
    try:
        return json.loads(lines[-1])
    except json.JSONDecodeError as exc:
        return {"error": f"worker stdout not JSON: {exc}", "stdout": lines[-1][:2000]}


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------


UNREACHED = -1


def build_instance_row(
    n: int,
    edge_prob: float,
    seed: int,
    gate: bool = False,
    max_attempts: int = GATE_MAX_ATTEMPTS,
) -> dict[str, Any] | None:
    """Everything computable without timing: query, keys, topology.

    ``gate=False`` (default): the original, ungated behaviour -- attempt 0
    only (``effective_seed = seed``), returns ``None`` if degenerate (no
    valid query). This keeps the committed ungated run reproducible.

    ``gate=True``: tries attempt ``a = 0, 1, ..., max_attempts``. Attempt 0
    uses ``effective_seed = seed``; attempt ``a >= 1`` uses
    ``effective_seed = seed + 10_000 * a``. Accepts the first attempt whose
    CPDAG has ``k_undirected >= 1`` AND has a valid query. If none of the
    ``max_attempts + 1`` attempts qualifies, returns a ``gate_exhausted``
    error-row dict (never ``None`` -- the cell is never silently dropped).
    Every returned row carries ``attempt`` and ``effective_seed``.
    """
    if not gate:
        cpdag, g0 = build_cpdag_and_g0(n, seed, edge_prob)
        q = choose_query(cpdag, g0)
        if q is None:
            return None
        x, y, z = q
        row = {
            "n": n,
            "edge_prob": edge_prob,
            "seed": seed,
            "attempt": 0,
            "effective_seed": seed,
            "x": x,
            "y": y,
            "z": sorted(z),
        }
        row.update(graph_stats(cpdag, g0))
        return row

    for a in range(0, max_attempts + 1):
        effective_seed = seed if a == 0 else seed + 10_000 * a
        cpdag, g0 = build_cpdag_and_g0(n, effective_seed, edge_prob)
        gstats = graph_stats(cpdag, g0)
        if gstats["k_undirected"] < 1:
            continue
        q = choose_query(cpdag, g0)
        if q is None:
            continue
        x, y, z = q
        row = {
            "n": n,
            "edge_prob": edge_prob,
            "seed": seed,
            "attempt": a,
            "effective_seed": effective_seed,
            "x": x,
            "y": y,
            "z": sorted(z),
        }
        row.update(gstats)
        return row

    return {
        "n": n,
        "edge_prob": edge_prob,
        "seed": seed,
        "attempt": max_attempts,
        "effective_seed": None,
        "gate_exhausted": True,
        "error": (
            f"gate exhausted: no attempt in 0..{max_attempts} had "
            "k_undirected >= 1 and a valid query"
        ),
    }


def load_completed(jsonl_path: Path) -> set[str]:
    if not jsonl_path.exists():
        return set()
    done = set()
    with jsonl_path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            done.add(instance_key(r["n"], r["edge_prob"], r["seed"]))
    return done


def radius_agreement(row: dict[str, Any]) -> dict[str, Any]:
    """Compare radii across completed legs. Never averages UNREACHED."""
    candidates: dict[str, Any] = {}
    space = row["legs"].get("space", {})
    if "r_space_enum" in space:
        candidates["space_enum"] = space["r_space_enum"]
    if "r_space_crit" in space:
        candidates["space_crit"] = space["r_space_crit"]
    sf = row["legs"].get("search_frozen", {})
    if "radius" in sf:
        candidates["search_frozen"] = sf["radius"]
    sfast = row["legs"].get("search_fast", {})
    if "radius" in sfast:
        candidates["search_fast"] = sfast["radius"]
    hyb = row["legs"].get("hybrid", {})
    if "radius" in hyb:
        candidates["hybrid"] = hyb["radius"]

    values = set(candidates.values())
    return {
        "values": candidates,
        "agree": len(values) <= 1,
        "n_legs_completed": len(candidates),
    }


def run_instance(row: dict[str, Any], cap_s: float) -> dict[str, Any]:
    # Legs are dispatched on the SAME (n, effective_seed, edge_prob) the driver
    # used to pick this row's query, so worker subprocesses rebuild a
    # bit-identical instance (effective_seed == seed when the gate is off, or
    # for attempt 0 under the gate).
    n, edge_prob, effective_seed, x, y, z = (
        row["n"],
        row["edge_prob"],
        row.get("effective_seed", row["seed"]),
        row["x"],
        row["y"],
        row["z"],
    )
    legs: dict[str, Any] = {}
    for leg in LEGS:
        legs[leg] = run_leg(leg, n, effective_seed, edge_prob, x, y, z, cap_s)
    out = dict(row)
    out["legs"] = legs
    out["radius_agreement"] = radius_agreement(out)
    return out


_CSV_FIELDS = (
    "n", "edge_prob", "seed", "attempt", "effective_seed",
    "x", "y", "z_size", "k_undirected", "k_g0", "largest_component",
    "space_build_s", "space_bfs_query_enum_s", "space_bfs_query_crit_s",
    "space_censored", "space_error", "space_total_crit_s", "space_total_enum_s",
    "search_frozen_s", "search_frozen_censored", "search_frozen_error",
    "search_fast_s", "search_fast_censored", "search_fast_error",
    "hybrid_total_s", "hybrid_censored", "hybrid_error", "hybrid_method",
    "radius_agree",
    "speedup_space_vs_search", "speedup_space_vs_search_lb",
    "speedup_search_vs_hybrid", "speedup_search_vs_hybrid_lb",
    "speedup_e2e", "speedup_e2e_lb",
    "speedup_space_vs_search_enum", "speedup_search_vs_hybrid_enum", "speedup_e2e_enum",
    "gate_exhausted", "gate_error",
)


def flatten_csv_row(row: dict[str, Any]) -> dict[str, Any]:
    if row.get("gate_exhausted"):
        # Never dropped -- recorded as an all-None error row with the same
        # column set as a normal row, so csv.DictWriter's fieldnames (taken
        # from row 0) are stable regardless of grid position.
        out = {k: None for k in _CSV_FIELDS}
        out.update(
            {
                "n": row["n"],
                "edge_prob": row["edge_prob"],
                "seed": row["seed"],
                "attempt": row["attempt"],
                "effective_seed": row.get("effective_seed"),
                "gate_exhausted": True,
                "gate_error": row.get("error"),
            }
        )
        return out

    legs = row["legs"]
    space = legs.get("space", {})
    sf = legs.get("search_frozen", {})
    sfast = legs.get("search_fast", {})
    hyb = legs.get("hybrid", {})

    def censored(leg: dict[str, Any]) -> bool:
        return bool(leg.get("censored"))

    def errored(leg: dict[str, Any]) -> bool:
        return "error" in leg

    space_total_crit = None
    space_total_enum = None
    if not censored(space) and not errored(space):
        space_total_crit = space["build_s"] + space["bfs_query_crit_s"]
        space_total_enum = space["build_s"] + space["bfs_query_enum_s"]

    search_frozen_s = sf.get("search_s") if not censored(sf) and not errored(sf) else None
    search_fast_s = sfast.get("search_s") if not censored(sfast) and not errored(sfast) else None
    hybrid_total_s = hyb.get("total_seconds") if not censored(hyb) and not errored(hyb) else None

    out: dict[str, Any] = {
        "n": row["n"],
        "edge_prob": row["edge_prob"],
        "seed": row["seed"],
        "attempt": row.get("attempt", 0),
        "effective_seed": row.get("effective_seed", row["seed"]),
        "x": row["x"],
        "y": row["y"],
        "z_size": len(row["z"]),
        "k_undirected": row["k_undirected"],
        "k_g0": row["k_g0"],
        "largest_component": row["largest_component"],
        "space_build_s": space.get("build_s"),
        "space_bfs_query_enum_s": space.get("bfs_query_enum_s"),
        "space_bfs_query_crit_s": space.get("bfs_query_crit_s"),
        "space_censored": censored(space),
        "space_error": errored(space),
        "space_total_crit_s": space_total_crit,
        "space_total_enum_s": space_total_enum,
        "search_frozen_s": search_frozen_s,
        "search_frozen_censored": censored(sf),
        "search_frozen_error": errored(sf),
        "search_fast_s": search_fast_s,
        "search_fast_censored": censored(sfast),
        "search_fast_error": errored(sfast),
        "hybrid_total_s": hybrid_total_s,
        "hybrid_censored": censored(hyb),
        "hybrid_error": errored(hyb),
        "hybrid_method": hyb.get("method") if not censored(hyb) and not errored(hyb) else None,
        "radius_agree": row["radius_agreement"]["agree"],
    }

    def ratio(num, den):
        if num is None or den is None or den == 0:
            return None
        return num / den

    def lb_ratio(num_leg, den):
        # numerator censored -> ratio is only a LOWER bound, using the cap as numerator
        if den is None or den == 0:
            return None
        if not censored(num_leg):
            return None
        return num_leg["wall_until_timeout_s"] / den

    out["speedup_space_vs_search"] = ratio(space_total_crit, search_fast_s)
    out["speedup_space_vs_search_lb"] = lb_ratio(space, search_fast_s)
    out["speedup_search_vs_hybrid"] = ratio(search_fast_s, hybrid_total_s)
    out["speedup_search_vs_hybrid_lb"] = lb_ratio(sfast, hybrid_total_s)
    out["speedup_e2e"] = ratio(space_total_crit, hybrid_total_s)
    out["speedup_e2e_lb"] = lb_ratio(space, hybrid_total_s)

    out["speedup_space_vs_search_enum"] = ratio(space_total_enum, search_frozen_s)
    out["speedup_search_vs_hybrid_enum"] = ratio(search_frozen_s, hybrid_total_s)
    out["speedup_e2e_enum"] = ratio(space_total_enum, hybrid_total_s)

    out["gate_exhausted"] = False
    out["gate_error"] = None

    return out


def write_manifest(out_dir: Path, args: argparse.Namespace, total_wall_s: float, n_rows: int) -> None:
    try:
        sha = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=str(ROOT), capture_output=True, text=True, check=True
        ).stdout.strip()
    except Exception:
        sha = None
    try:
        dirty = bool(
            subprocess.run(
                ["git", "status", "--porcelain"], cwd=str(ROOT), capture_output=True, text=True, check=True
            ).stdout.strip()
        )
    except Exception:
        dirty = None

    import importlib.metadata as im

    def pkg_version(name: str) -> str | None:
        try:
            return im.version(name)
        except Exception:
            return None

    manifest = {
        "created_utc": __import__("datetime").datetime.now(__import__("datetime").timezone.utc).isoformat(),
        "git_sha": sha,
        "git_dirty": dirty,
        "python": sys.version,
        "platform": platform.platform(),
        "packages": {
            "numpy": pkg_version("numpy"),
            "networkx": pkg_version("networkx"),
            "ortools": pkg_version("ortools"),
        },
        "grid": {
            "n": list(N_GRID),
            "edge_prob": list(EDGE_PROB_GRID),
            "seed": list(SEED_GRID),
            "n_instances_grid": len(instance_grid()),
            "n_rows_written": n_rows,
        },
        "cap_s": CAP_S,
        "legs": list(LEGS),
        "oracles": {
            "space_leg": "both enumeration (bkrobust.core.oracle.is_valid) and criterion (bkrobust.mpdag_criterion.criterion.is_valid_mpdag), same BFS distances reused for both",
            "search_frozen": "enumeration oracle, bkrobust.search.exact.radius_local_up, unbounded",
            "search_fast": "criterion oracle, bkrobust.search.exact_fast.radius_local_up_fast, unbounded (no max_depth)",
            "hybrid": "bkrobust.hybrid.breakdown_radius defaults: criterion oracle, search_budget=3",
        },
        "subprocess_timing_discipline": (
            "Every leg for every instance runs in its own fresh Python subprocess "
            "(PYTHONPATH=src), timed inside the child with time.perf_counter() "
            "around the computation only (imports and instance construction "
            "excluded from the timed interval). Legs are run strictly serially, "
            "one at a time, never concurrently, to avoid wall-clock contamination "
            "between timed processes. A leg exceeding the 120s cap is killed by "
            "subprocess.run(timeout=...) and recorded as "
            "{censored: true, wall_until_timeout_s: <wall>} with no seconds key."
        ),
        "query_rule": (
            "The committed results/axisb4/hybrid_envelope.jsonl query rule could "
            "not be reconstructed from dense_instance(n, seed, edge_prob) with its "
            "own (n, seed, edge_prob) fields -- see the e2e_speedup.py module "
            "docstring for the topology-level proof (row 0: x=V0 is isolated in "
            "dense_instance(8, 0, 0.3), incompatible with its committed z_size=3). "
            "This run instead builds its own 256-instance grid and its own "
            "deterministic query per instance: first (x, y) in sorted-node "
            "permutation order with a non-empty, valid optimal_adjustment_set_mpdag "
            "under the enumeration oracle, matching the iteration order in "
            "bkrobust.search.scaling.measure_instance. Verification against the "
            "committed envelope (same n/edge_prob/seed, committed x/y, freshly "
            "computed z_size) matched on only 6/256 rows, consistent with chance."
        ),
        "command_line": " ".join(sys.argv),
        "total_wall_s": total_wall_s,
        "gate": {
            "enabled": bool(args.gate),
            "description": (
                "k>=1-gated draw of the same grid -- NOT a reproduction of "
                "results/axisb4/hybrid_envelope.jsonl (see module docstring "
                "and query_rule_verification.json)."
            ),
            "max_attempts": getattr(args, "gate_max_attempts", GATE_MAX_ATTEMPTS),
            "rule": (
                "For each grid cell (n, edge_prob, seed), attempt a = 0, 1, 2, ..., "
                "max_attempts. Attempt 0 uses effective_seed = seed (dense_instance(n, "
                "seed, edge_prob), so any cell that already passes keeps its prior "
                "instance). Attempt a >= 1 uses effective_seed = seed + 10_000 * a "
                "(no other source of randomness; no global RNG is touched). Accept the "
                "first attempt where the CPDAG has k_undirected >= 1 AND choose_query() "
                "finds a query (first sorted-permutation (x, y) with a non-empty, valid "
                "optimal adjustment set). If no attempt in 0..max_attempts qualifies, the "
                "cell is recorded as a gate_exhausted error row -- never silently dropped. "
                "Worker subprocesses rebuild the instance from (n, effective_seed, "
                "edge_prob), so they are bit-identical to the driver's choice."
            ),
        },
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, default=str) + "\n")


def verify_query_rule_against_envelope() -> dict[str, Any]:
    """Check reconstruction against results/axisb4/hybrid_envelope.jsonl. Read-only."""
    envelope_path = ROOT / "results" / "axisb4" / "hybrid_envelope.jsonl"
    if not envelope_path.exists():
        return {"checked": False, "reason": "envelope file not found"}

    from bkrobust.core.oracle import is_valid
    from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag

    rows = [json.loads(line) for line in envelope_path.read_text().splitlines() if line.strip()]
    n_total = len(rows)
    n_match_own_rule = 0
    n_match_fallback_zsize = 0
    n_topology_impossible_checked = 0
    for r in rows:
        cpdag, g0 = build_cpdag_and_g0(r["n"], r["seed"], r["edge_prob"])
        q = choose_query(cpdag, g0)
        if q is not None:
            x, y, z = q
            if x == r["x"] and y == r["y"] and len(z) == r["hybrid"]["z_size"]:
                n_match_own_rule += 1
        o = optimal_adjustment_set_mpdag(g0, r["x"], r["y"])
        if o is not None and len(o) == r["hybrid"]["z_size"]:
            n_match_fallback_zsize += 1
        n_topology_impossible_checked += 1

    return {
        "checked": True,
        "n_committed_rows": n_total,
        "n_match_own_query_rule": n_match_own_rule,
        "n_match_fallback_zsize_only": n_match_fallback_zsize,
        "n_rows_checked_for_zsize_fallback": n_topology_impossible_checked,
        "conclusion": (
            "dense_instance(n, seed, edge_prob) with the committed fields does not "
            "reproduce the committed envelope's instances or queries; this driver "
            "uses its own self-consistent grid and query rule instead (see manifest.json)."
        ),
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", default="results/e2e_speedup")
    p.add_argument("--limit", type=int, default=None, help="only process the first N grid cells (smoke test)")
    p.add_argument("--resume", action="store_true", help="skip instances already present in the output jsonl")
    p.add_argument("--cap-s", type=float, default=CAP_S)
    p.add_argument("--worker", default=None, choices=list(_WORKERS.keys()), help=argparse.SUPPRESS)
    p.add_argument("--worker-payload", default=None, help=argparse.SUPPRESS)
    p.add_argument("--skip-verify", action="store_true", help="skip the (slow) envelope verification pass")
    p.add_argument(
        "--gate",
        action="store_true",
        help=(
            "enable the deterministic k>=1 resample gate (see module docstring). "
            "Default OFF, so the committed ungated run in results/e2e_speedup "
            "stays reproducible with its recorded command line. This is a "
            "k>=1-gated draw of the SAME grid, NOT a reproduction of "
            "results/axisb4/hybrid_envelope.jsonl."
        ),
    )
    p.add_argument(
        "--gate-max-attempts",
        type=int,
        default=GATE_MAX_ATTEMPTS,
        help="max resample attempts per cell under --gate (attempts 0..N inclusive)",
    )
    args = p.parse_args()

    if args.worker is not None:
        run_worker(args.worker, json.loads(args.worker_payload))
        return

    out_dir = ROOT / args.out
    out_dir.mkdir(parents=True, exist_ok=True)
    jsonl_path = out_dir / "e2e_speedup.jsonl"
    csv_path = out_dir / "e2e_speedup.csv"

    already = load_completed(jsonl_path) if args.resume else set()

    grid = instance_grid()
    if args.limit is not None:
        grid = grid[: args.limit]

    wall0 = time.perf_counter()

    verification = None
    if not args.skip_verify:
        print("[e2e_speedup] verifying query-rule reconstruction against committed envelope...", flush=True)
        verification = verify_query_rule_against_envelope()
        (out_dir / "query_rule_verification.json").write_text(json.dumps(verification, indent=2) + "\n")
        print(f"[e2e_speedup] verification: {verification}", flush=True)

    mode = "a" if (args.resume and jsonl_path.exists()) else "w"
    n_skipped_degenerate = 0
    n_gate_exhausted = 0
    n_written = 0
    with jsonl_path.open(mode) as jf:
        for i, (n, ep, seed) in enumerate(grid):
            key = instance_key(n, ep, seed)
            if key in already:
                continue
            row = build_instance_row(n, ep, seed, gate=args.gate, max_attempts=args.gate_max_attempts)
            if row is None:
                # Only reachable with --gate unset (old, ungated skip behaviour).
                n_skipped_degenerate += 1
                print(f"[e2e_speedup] ({i+1}/{len(grid)}) {key}: no valid query found, skipping", flush=True)
                continue
            if row.get("gate_exhausted"):
                # Never silently dropped -- recorded as an error row with no legs run.
                jf.write(json.dumps({**row, "legs": {}, "radius_agreement": {"values": {}, "agree": None, "n_legs_completed": 0}}, default=str) + "\n")
                jf.flush()
                n_written += 1
                n_gate_exhausted += 1
                print(
                    f"[e2e_speedup] ({i+1}/{len(grid)}) {key}: GATE EXHAUSTED after "
                    f"{args.gate_max_attempts + 1} attempts, recorded as error row",
                    flush=True,
                )
                continue
            t_inst0 = time.perf_counter()
            result = run_instance(row, args.cap_s)
            t_inst = time.perf_counter() - t_inst0
            jf.write(json.dumps(result, default=str) + "\n")
            jf.flush()
            n_written += 1
            agree = result["radius_agreement"]["agree"]
            print(
                f"[e2e_speedup] ({i+1}/{len(grid)}) {key} attempt={row.get('attempt', 0)} "
                f"eff_seed={row.get('effective_seed', seed)} x={row['x']} y={row['y']} "
                f"|z|={len(row['z'])} k_g0={row['k_g0']} k_und={row['k_undirected']} "
                f"agree={agree} took={t_inst:.1f}s",
                flush=True,
            )

    # Rebuild the CSV from the full jsonl (so --resume runs regenerate it fully).
    all_rows = []
    with jsonl_path.open() as jf:
        for line in jf:
            line = line.strip()
            if line:
                all_rows.append(json.loads(line))
    csv_rows = [flatten_csv_row(r) for r in all_rows]
    if csv_rows:
        fieldnames = list(csv_rows[0].keys())
        with csv_path.open("w", newline="") as cf:
            writer = csv.DictWriter(cf, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(csv_rows)

    total_wall_s = time.perf_counter() - wall0
    write_manifest(out_dir, args, total_wall_s, len(all_rows))
    write_summary(out_dir, all_rows, csv_rows, verification, n_skipped_degenerate, args.gate, args.gate_max_attempts)

    print(
        f"[e2e_speedup] done: {n_written} instances written this run, "
        f"{len(all_rows)} total rows, {n_skipped_degenerate} degenerate skips, "
        f"{n_gate_exhausted} gate-exhausted error rows, {total_wall_s:.1f}s wall",
        flush=True,
    )


def _median(values: list[float]) -> float | None:
    vals = [v for v in values if v is not None]
    if not vals:
        return None
    return float(np.median(vals))


def write_summary(
    out_dir: Path,
    all_rows: list[dict[str, Any]],
    csv_rows: list[dict[str, Any]],
    verification: dict[str, Any] | None,
    n_skipped_degenerate: int,
    gate_enabled: bool = False,
    gate_max_attempts: int = GATE_MAX_ATTEMPTS,
) -> None:
    gate_exhausted_rows = [r for r in csv_rows if r.get("gate_exhausted")]
    ok_rows = [r for r in csv_rows if not r.get("gate_exhausted")]

    attempts_histogram: dict[str, int] = {}
    for r in ok_rows:
        a = r.get("attempt")
        if a is None:
            continue
        attempts_histogram[str(a)] = attempts_histogram.get(str(a), 0) + 1

    by_k: dict[int, list[dict[str, Any]]] = {}
    for r in ok_rows:
        by_k.setdefault(r["k_undirected"], []).append(r)

    per_k = {}
    for k, rs in sorted(by_k.items()):
        space_completed = [r for r in rs if not r["space_censored"] and not r["space_error"]]
        space_censored = [r for r in rs if r["space_censored"]]
        sf_completed = [r for r in rs if not r["search_frozen_censored"] and not r["search_frozen_error"]]
        sfast_completed = [r for r in rs if not r["search_fast_censored"] and not r["search_fast_error"]]
        hyb_completed = [r for r in rs if not r["hybrid_censored"] and not r["hybrid_error"]]

        per_k[str(k)] = {
            "n_rows": len(rs),
            "space": {
                "n_completed": len(space_completed),
                "n_censored": len(space_censored),
                "median_space_total_crit_s": _median([r["space_total_crit_s"] for r in space_completed]),
            },
            "search_frozen": {
                "n_completed": len(sf_completed),
                "n_censored": sum(1 for r in rs if r["search_frozen_censored"]),
                "median_s": _median([r["search_frozen_s"] for r in sf_completed]),
            },
            "search_fast": {
                "n_completed": len(sfast_completed),
                "n_censored": sum(1 for r in rs if r["search_fast_censored"]),
                "median_s": _median([r["search_fast_s"] for r in sfast_completed]),
            },
            "hybrid": {
                "n_completed": len(hyb_completed),
                "n_censored": sum(1 for r in rs if r["hybrid_censored"]),
                "median_s": _median([r["hybrid_total_s"] for r in hyb_completed]),
            },
            "median_speedup_space_vs_search": _median([r["speedup_space_vs_search"] for r in rs]),
            "median_speedup_search_vs_hybrid": _median([r["speedup_search_vs_hybrid"] for r in rs]),
            "median_speedup_e2e": _median([r["speedup_e2e"] for r in rs]),
            "n_lower_bound_rows": sum(
                1
                for r in rs
                if r["speedup_space_vs_search_lb"] is not None
                or r["speedup_search_vs_hybrid_lb"] is not None
                or r["speedup_e2e_lb"] is not None
            ),
        }

    # Radius agreement is only meaningful for rows that ran legs; gate_exhausted
    # rows carry agree=None (vacuous) and are excluded here, not counted as
    # disagreements.
    agreement_rows = [r for r in all_rows if r["radius_agreement"]["agree"] is not None]
    n_agree = sum(1 for r in agreement_rows if r["radius_agreement"]["agree"])
    n_disagree = len(agreement_rows) - n_agree
    disagreements = [
        {"n": r["n"], "edge_prob": r["edge_prob"], "seed": r["seed"], "values": r["radius_agreement"]["values"]}
        for r in agreement_rows
        if not r["radius_agreement"]["agree"]
    ]

    summary = {
        "n_rows": len(all_rows),
        "n_skipped_degenerate": n_skipped_degenerate,
        "gate": {
            "enabled": bool(gate_enabled),
            "max_attempts": gate_max_attempts,
            "description": "k>=1-gated draw of the same grid -- NOT the session-4 envelope",
            "n_cells_gate_exhausted": len(gate_exhausted_rows),
            "n_cells_resolved": len(ok_rows),
            "attempts_histogram": dict(sorted(attempts_histogram.items(), key=lambda kv: int(kv[0]))),
            "n_cells_needing_resample": sum(v for a, v in attempts_histogram.items() if a != "0"),
        },
        "per_k_undirected": per_k,
        "overall": {
            "median_speedup_space_vs_search": _median([r["speedup_space_vs_search"] for r in ok_rows]),
            "median_speedup_search_vs_hybrid": _median([r["speedup_search_vs_hybrid"] for r in ok_rows]),
            "median_speedup_e2e": _median([r["speedup_e2e"] for r in ok_rows]),
        },
        "radius_agreement": {
            "n_agree": n_agree,
            "n_disagree": n_disagree,
            "disagreements": disagreements[:50],
        },
        "query_rule_verification": verification,
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")


if __name__ == "__main__":
    main()
