"""E4: listing the Markov equivalence class as a baseline for the breakdown radius.

Corollary 1(c): r = min |K_G0 \\ dir(D)| over DAGs D in [C^] where Z fails. This
driver times that formula, computed by listing [C^], against
bkrobust.hybrid.breakdown_radius on three fixed instance families, and checks
that both give the same radius. See PREREG.md for the registered statistic.

Run from the repository root:
    .venv/bin/python experiments/rebuttal/E4_class_enumeration_baseline/run.py --smoke
"""

from __future__ import annotations

import argparse
import itertools
import json
import multiprocessing as mp
import os
import platform
import subprocess
import sys
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "src"))

import numpy as np  # noqa: E402

from bkrobust.core.conventions import UNREACHED  # noqa: E402
from bkrobust.core.oracle import clear_cache as oracle_clear_cache  # noqa: E402
from bkrobust.core.oracle import is_valid as backdoor_valid_mpdag  # noqa: E402
from bkrobust.demo.evaluate import is_valid_adjustment_set_dag, optimal_adjustment_set_mpdag  # noqa: E402
from bkrobust.demo.example import dag_to_cpdag, knowledge_to_recover  # noqa: E402
from bkrobust.demo.graph import MPDAG, undirected_components  # noqa: E402
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions, is_consistent_extension  # noqa: E402
from bkrobust.gac.dag_level import is_gac_valid_dag  # noqa: E402
from bkrobust.gac.mpdag_level import clear_cache as gac_clear_cache  # noqa: E402
from bkrobust.gac.mpdag_level import is_gac_valid_mpdag  # noqa: E402
from bkrobust.hybrid import DEFAULT_SEARCH_BUDGET, breakdown_radius  # noqa: E402
from bkrobust.search.conjecture_study import all_cpdags  # noqa: E402
from bkrobust.search.space_fixed import enumerate_space_correct  # noqa: E402
from bkrobust.synth.generators import erdos_renyi_dag  # noqa: E402

FOUR_NODE_EXPECTED = 1044  # App. D; reported against, never enforced


class DeadBranch(RuntimeError):
    """Orienting an undirected edge of a closed MPDAG failed: contradicts Meek (1995), so a bug."""


class _Timeout(Exception):
    pass


def git_commit() -> str:
    try:
        out = subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"], capture_output=True, text=True, timeout=10)
        sha = out.stdout.strip() or "unknown"
        dirty = subprocess.run(
            ["git", "-C", str(REPO), "status", "--porcelain", "--untracked-files=no"],
            capture_output=True, text=True, timeout=10,
        ).stdout.strip()
        return sha + ("-dirty" if dirty else "")
    except Exception:
        return "unknown"


# --------------------------------------------------------------------------- instances


def _edges(g: MPDAG) -> dict:
    return {"nodes": list(g.nodes), "directed": sorted(g.directed_edges), "undirected": sorted(g.undirected_edges)}


def _graph(d: dict) -> MPDAG:
    return MPDAG(d["nodes"], [tuple(e) for e in d["directed"]], [tuple(e) for e in d["undirected"]])


def four_node_units(filter_criterion: str):
    """Every (C^, G0, X, Y) on 4 nodes with a non-empty optimal Z valid in G0.

    Same loop as tests/hybrid/test_hybrid_gac.py, which is the closest written
    trace of how App. D's 1,044 instances were drawn.
    """
    i = 0
    for cpdag in all_cpdags(4):
        if not cpdag.undirected_edges:
            continue
        for g0 in enumerate_space_correct(cpdag):
            for x, y in itertools.permutations(sorted(cpdag.nodes), 2):
                o = optimal_adjustment_set_mpdag(g0, x, y)
                if o is None or not o:
                    continue
                z = frozenset(o)
                ok = backdoor_valid_mpdag(z, g0, x, y) if filter_criterion == "backdoor" else is_gac_valid_mpdag(g0, x, y, z)
                if not ok:
                    continue
                yield {
                    "unit_id": f"four:{i:05d}", "family": "four", "cpdag": _edges(cpdag), "g0": _edges(g0),
                    "x": x, "y": y, "z": sorted(z), "filter_criterion": filter_criterion,
                }
                i += 1


def highk_units(path: Path) -> list[dict]:
    out = []
    for i, line in enumerate(path.read_text().splitlines()):
        if not line.strip():
            continue
        r = json.loads(line)
        out.append({
            "unit_id": f"highk:{i:04d}", "family": "highk", "n": r["n"], "seed": r["seed"], "edge_prob": r["edge_prob"],
            "x": r["x"], "y": r["y"], "z": None, "stored": {"z_size": r["z_size"], "m_undirected": r["m_undirected"], "k_g0": r["k_g0"]},
        })
    return out


def timing_units(path: Path) -> list[dict]:
    out = []
    for i, line in enumerate(path.read_text().splitlines()):
        if not line.strip():
            continue
        r = json.loads(line)
        if "x" not in r:  # gate_exhausted error rows carry no query
            continue
        out.append({
            "unit_id": f"timing:{i:04d}", "family": "timing", "n": r["n"], "seed": r["effective_seed"], "edge_prob": r["edge_prob"],
            "x": r["x"], "y": r["y"], "z": sorted(r["z"]), "stored": {"m_undirected": r["k_undirected"], "k_g0": r["k_g0"]},
        })
    return out


def build_instance(u: dict):
    """(cpdag, g0, x, y, z, mismatch) with mismatch a list of fields that differ from the stored row."""
    if u["family"] == "four":
        return _graph(u["cpdag"]), _graph(u["g0"]), u["x"], u["y"], frozenset(u["z"]), []
    # Same calls, in the same order, as bkrobust.gac.sweep.dense_instance and experiments/e2e_speedup.py.
    rng = np.random.default_rng(u["seed"])
    dag = erdos_renyi_dag(u["n"], rng, u["edge_prob"])
    cpdag = dag_to_cpdag(dag)
    g0 = apply_orientations(cpdag, sorted(knowledge_to_recover(dag, cpdag)))
    if g0 is None:
        return cpdag, None, u["x"], u["y"], frozenset(), ["g0_failed"]
    x, y = u["x"], u["y"]
    mismatch = []
    if u["z"] is None:
        o = optimal_adjustment_set_mpdag(g0, x, y)
        z = frozenset(o) if o is not None else frozenset()
        if o is None or len(z) != u["stored"]["z_size"]:
            mismatch.append("z_size")
    else:
        z = frozenset(u["z"])
    if len(cpdag.undirected_edges) != u["stored"]["m_undirected"]:
        mismatch.append("m_undirected")
    if len(g0.directed_edges - cpdag.directed_edges) != u["stored"]["k_g0"]:
        mismatch.append("k_g0")
    return cpdag, g0, x, y, z, mismatch


# --------------------------------------------------------------------------- class size (independent of the listing)


def _count_uccg(nodes: frozenset, edges: frozenset, memo: dict, deadline: float) -> int:
    """Number of DAGs of a connected chordal undirected graph, by root-picking.

    He, Jia and Yu (2015): each such DAG has one source; fixing v as the source
    and closing leaves chordal chain components whose counts multiply.
    """
    if not edges:
        return 1
    key = (nodes, edges)
    if key in memo:
        return memo[key]
    if time.perf_counter() > deadline:
        raise _Timeout
    g = MPDAG(nodes, [], edges)
    total = 0
    for v in sorted(nodes):
        h = apply_orientations(g, [(v, u) for u in sorted(g.adjacent(v))])
        if h is None:
            raise DeadBranch(f"root {v} not consistent in a chordal component")
        prod = 1
        for comp in undirected_components(h):
            ce = frozenset(e for e in h.undirected_edges if e[0] in comp and e[1] in comp)
            prod *= _count_uccg(frozenset(comp), ce, memo, deadline)
        total += prod
    memo[key] = total
    return total


def class_size(cpdag: MPDAG, limit_s: float) -> int | None:
    deadline = time.perf_counter() + limit_s
    memo: dict = {}
    size = 1
    try:
        for comp in undirected_components(cpdag):
            ce = frozenset(e for e in cpdag.undirected_edges if e[0] in comp and e[1] in comp)
            size *= _count_uccg(frozenset(comp), ce, memo, deadline)
    except _Timeout:
        return None
    return size


# --------------------------------------------------------------------------- listing [C^]


def iter_extensions(g: MPDAG, prefer: frozenset = frozenset()):
    """Depth-first listing of [g]: orient the first undirected edge both ways, close, recurse.

    In a closed MPDAG every undirected edge can be oriented both ways (Meek 1995),
    so no branch is dead and every leaf is a distinct member of [g].
    """
    if not g.undirected_edges:
        yield g
        return
    a, b = min(g.undirected_edges)
    for o in (((b, a), (a, b)) if (b, a) in prefer else ((a, b), (b, a))):
        h = apply_orientations(g, [o])
        if h is None:
            raise DeadBranch(f"orientation {o} failed")
        yield from iter_extensions(h, prefer)


def leg_enum(cpdag, x, y, z, k0, valid_dag, cap, limit_s) -> dict:
    """Corollary 1(c) by plain listing; no pruning."""
    t0 = time.perf_counter()
    best, n, n_fail, stop = None, 0, 0, None
    it = iter_extensions(cpdag)
    for d in it:
        n += 1
        if not valid_dag(d, x, y, z):
            n_fail += 1
            c = len(k0 - d.directed_edges)
            best = c if best is None else min(best, c)
        if n >= cap:
            if next(it, None) is not None:  # one more leaf exists: the listing is incomplete
                stop = "cap"
            break
        if time.perf_counter() - t0 > limit_s:
            stop = "time"
            break
    sec = time.perf_counter() - t0
    return {"radius": UNREACHED if best is None else best, "found_failing": best is not None, "n_dags": n,
            "n_failing": n_fail, "complete": stop is None, "stop": stop, "seconds": sec}


def leg_bnb(cpdag, x, y, z, k0, valid_dag, cap, limit_s) -> dict:
    """Same listing, pruned: reversed K_G0 orientations only grow down a branch."""
    k_rev = frozenset((b, a) for a, b in k0)
    st = {"best": None, "leaves": 0, "nodes": 0, "stop": None}
    t0 = time.perf_counter()

    def rec(g):
        if st["stop"]:
            return
        st["nodes"] += 1
        if time.perf_counter() - t0 > limit_s:
            st["stop"] = "time"
            return
        lb = len(k_rev & g.directed_edges)
        if st["best"] is not None and lb >= st["best"]:
            return
        if not g.undirected_edges:
            st["leaves"] += 1
            if not valid_dag(g, x, y, z):
                st["best"] = lb  # lb < best here, and lb is exact at a leaf
            if st["leaves"] >= cap:
                st["stop"] = "cap"  # conservative: may be complete, never claimed so
            return
        a, b = min(g.undirected_edges)
        for o in (((b, a), (a, b)) if (b, a) in k0 else ((a, b), (b, a))):
            h = apply_orientations(g, [o])
            if h is None:
                raise DeadBranch(f"orientation {o} failed")
            rec(h)

    rec(cpdag)
    sec = time.perf_counter() - t0
    best = st["best"]
    return {"radius": UNREACHED if best is None else best, "found_failing": best is not None, "leaves": st["leaves"],
            "tree_nodes": st["nodes"], "complete": st["stop"] is None, "stop": st["stop"], "seconds": sec}


def leg_lib(cpdag, x, y, z, k0, valid_dag) -> dict:
    """The library's brute force over 2^k orientations, as an independent listing."""
    t0 = time.perf_counter()
    exts = enumerate_dag_extensions(cpdag)
    best = None
    for d in exts:
        if not valid_dag(d, x, y, z):
            c = len(k0 - d.directed_edges)
            best = c if best is None else min(best, c)
    sec = time.perf_counter() - t0
    return {"radius": UNREACHED if best is None else best, "n_dags": len(exts), "seconds": sec,
            "dag_set": sorted(d.edge_string() for d in exts)}


def leg_hybrid(cpdag, g0, x, y, z, criterion, rung_limit_s) -> dict:
    gac_clear_cache()
    oracle_clear_cache()
    t0 = time.perf_counter()
    res = breakdown_radius(cpdag, None, x, y, z, g0=g0, criterion=criterion, time_limit_s=rung_limit_s)
    wall = time.perf_counter() - t0
    return {"radius": res.radius, "method": res.method, "oracle": res.oracle, "exact": res.exact,
            "total_seconds": res.total_seconds, "wall_seconds": wall, "search_budget": DEFAULT_SEARCH_BUDGET}


# --------------------------------------------------------------------------- one unit


def run_unit(job: tuple) -> dict:
    u, cfg = job
    row = {"unit_id": u["unit_id"], "family": u["family"], "criterion": cfg["criterion"], "commit": cfg["commit"],
           "workers": cfg["workers"], "smoke": cfg["smoke"], "cap": cfg["cap"], "enum_limit_s": cfg["enum_limit_s"]}
    for k in ("n", "seed", "edge_prob", "x", "y", "filter_criterion"):
        if k in u:
            row[k] = u[k]
    failures: list[str] = []
    try:
        cpdag, g0, x, y, z, mismatch = build_instance(u)
        if mismatch:
            row.update(status="excluded", exclusion=mismatch, failures=[], consistent=True)
            return row
        valid_dag = is_gac_valid_dag if cfg["criterion"] == "gac" else is_valid_adjustment_set_dag
        k0 = frozenset(g0.directed_edges - cpdag.directed_edges)
        row.update(z=sorted(z), n_nodes=len(cpdag.nodes), m_undirected=len(cpdag.undirected_edges), k_g0=len(k0),
                   largest_component=max((len(c) for c in undirected_components(cpdag)), default=1))

        t0 = time.perf_counter()
        size = class_size(cpdag, cfg["count_limit_s"])
        row.update(class_size=size, class_count_seconds=time.perf_counter() - t0)

        hyb = leg_hybrid(cpdag, g0, x, y, z, cfg["criterion"], cfg["rung_limit_s"])
        enum = leg_enum(cpdag, x, y, z, k0, valid_dag, cfg["cap"], cfg["enum_limit_s"])
        row.update(hybrid=hyb, enum=enum)
        if "bnb" in cfg["legs"]:
            row["bnb"] = leg_bnb(cpdag, x, y, z, k0, valid_dag, cfg["cap"], cfg["enum_limit_s"])
        lib = None
        if "lib" in cfg["legs"] and row["m_undirected"] <= cfg["lib_max_k"]:
            lib = leg_lib(cpdag, x, y, z, k0, valid_dag)

        # ---- consistency, outside every timed region
        hr, h_ok = hyb["radius"], hyb["exact"]
        if not h_ok:
            row["flag_hybrid_inexact"] = True
        if enum["complete"]:
            if h_ok and enum["radius"] != hr:
                failures.append(f"enum_vs_hybrid: enum {enum['radius']} hybrid {hr}")
            if size is not None and enum["n_dags"] != size:
                failures.append(f"enum_count_vs_class_size: {enum['n_dags']} vs {size}")
        elif enum["found_failing"] and h_ok and (hr == UNREACHED or hr > enum["radius"]):
            failures.append(f"enum_upper_bound: capped enum found {enum['radius']} below hybrid {hr}")
        b = row.get("bnb")
        if b is not None:
            if b["complete"] and h_ok and b["radius"] != hr:
                failures.append(f"bnb_vs_hybrid: bnb {b['radius']} hybrid {hr}")
            if not b["complete"] and b["found_failing"] and h_ok and (hr == UNREACHED or hr > b["radius"]):
                failures.append(f"bnb_upper_bound: {b['radius']} below hybrid {hr}")
        if lib is not None:
            mine = sorted(d.edge_string() for d in iter_extensions(cpdag))
            dag_set = lib.pop("dag_set")
            lib["same_dag_set"] = mine == dag_set
            row["lib"] = lib
            if not lib["same_dag_set"]:
                failures.append(f"lib_dag_set: lib {len(dag_set)} vs listing {len(mine)}")
            if h_ok and lib["radius"] != hr:
                failures.append(f"lib_vs_hybrid: lib {lib['radius']} hybrid {hr}")
        n_ver = 0
        for d in itertools.islice(iter_extensions(cpdag), cfg["verify_first"]):
            n_ver += 1
            if not is_consistent_extension(d, cpdag):
                failures.append(f"not_an_extension: {d.edge_string()}")
                break
        row["verified_extensions"] = n_ver
        row["status"] = "ok"
    except DeadBranch as e:
        failures.append(f"dead_branch: {e}")
        row["status"] = "error"
    except Exception:
        failures.append("exception")
        row.update(status="error", traceback=traceback.format_exc()[-2000:])
    row["failures"] = failures
    row["consistent"] = not failures
    return row


# --------------------------------------------------------------------------- driver


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--smoke", action="store_true", help="tiny deterministic subset, under 60 s on one core")
    p.add_argument("--out", default=None, help="output dir (default results/rebuttal/E4, or results/rebuttal/E4_smoke with --smoke)")
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--shard", type=int, default=0)
    p.add_argument("--n-shards", type=int, default=1)
    p.add_argument("--families", default="four,highk,timing")
    p.add_argument("--legs", default="hybrid,enum,bnb,lib", help="hybrid and enum always run")
    p.add_argument("--criterion", choices=("gac", "backdoor"), default="gac")
    p.add_argument("--four-filter", choices=("backdoor", "gac"), default="backdoor",
                   help="validity filter selecting four-node instances (backdoor reproduces App. D's selection)")
    p.add_argument("--highk", default="results/axisb3/highk.jsonl")
    p.add_argument("--timing", default="results/e2e_speedup_gated_v2/e2e_speedup.jsonl")
    p.add_argument("--cap", type=int, default=100_000, help="DAGs listed per class before stopping")
    p.add_argument("--enum-limit-s", type=float, default=120.0, help="wall limit per listing leg")
    p.add_argument("--count-limit-s", type=float, default=120.0, help="wall limit for the class-size count")
    p.add_argument("--rung-limit-s", type=float, default=300.0, help="hybrid SAT rung limit")
    p.add_argument("--lib-max-k", type=int, default=14, help="run the 2^k library listing only up to this k")
    p.add_argument("--verify-first", type=int, default=200, help="listed DAGs re-checked as consistent extensions")
    a = p.parse_args()

    os.chdir(REPO)  # every default path is relative to the repository root
    if a.smoke:
        # cap 20 is below the largest four-node class (24) and most highk classes, so the
        # capped branches and the upper-bound check run too
        a.cap, a.enum_limit_s, a.count_limit_s, a.lib_max_k, a.verify_first = 20, 10.0, 10.0, 12, 50
    out = Path(a.out or ("results/rebuttal/E4_smoke" if a.smoke else "results/rebuttal/E4"))
    out.mkdir(parents=True, exist_ok=True)
    fams = set(a.families.split(","))
    commit = git_commit()

    t_build = time.perf_counter()
    units: list[dict] = []
    n_four = None
    if "four" in fams:
        gen = four_node_units(a.four_filter)
        four = list(itertools.islice(gen, 12)) if a.smoke else list(gen)
        n_four = None if a.smoke else len(four)
        units += four
        if n_four is not None and n_four != FOUR_NODE_EXPECTED:
            print(f"note: four-node selection gives {n_four} instances, App. D prints {FOUR_NODE_EXPECTED}", flush=True)
    if "highk" in fams:
        hk = highk_units(Path(a.highk))
        units += sorted(hk, key=lambda u: (u["stored"]["m_undirected"], u["unit_id"]))[:2] if a.smoke else hk
    if "timing" in fams:
        tm = timing_units(Path(a.timing))
        units += sorted(tm, key=lambda u: (u["stored"]["m_undirected"], u["unit_id"]))[:3] if a.smoke else tm
    mine = [u for i, u in enumerate(units) if i % a.n_shards == a.shard]

    tag = f"shard{a.shard:03d}of{a.n_shards:03d}"
    rows_path = out / f"rows_{tag}.jsonl"
    done = set()
    if rows_path.exists():
        for line in rows_path.read_text().splitlines():
            try:
                done.add(json.loads(line)["unit_id"])
            except (json.JSONDecodeError, KeyError):
                pass  # a line cut by a kill; its unit is redone
        # drop a truncated tail so appended rows start on a fresh line
        txt = rows_path.read_text()
        if txt and not txt.endswith("\n"):
            rows_path.write_text(txt[: txt.rfind("\n") + 1])
    todo = [u for u in mine if u["unit_id"] not in done]

    cfg = {"criterion": a.criterion, "commit": commit, "workers": a.workers, "smoke": a.smoke, "cap": a.cap,
           "enum_limit_s": a.enum_limit_s, "count_limit_s": a.count_limit_s, "rung_limit_s": a.rung_limit_s,
           "lib_max_k": a.lib_max_k, "verify_first": a.verify_first, "legs": a.legs.split(",")}
    manifest = {"experiment": "E4", "config": cfg, "four_filter": a.four_filter, "families": sorted(fams),
                "inputs": {"highk": a.highk, "timing": a.timing}, "n_units_total": len(units), "n_units_shard": len(mine),
                "n_four_node": n_four, "four_node_expected": FOUR_NODE_EXPECTED, "shard": a.shard, "n_shards": a.n_shards,
                "python": sys.version.split()[0], "platform": platform.platform(), "argv": sys.argv,
                "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "unit_build_seconds": time.perf_counter() - t_build}
    (out / f"manifest_{tag}.json").write_text(json.dumps(manifest, indent=1))
    print(f"E4 {tag}: {len(mine)} units, {len(done)} already written, {len(todo)} to run, commit {commit}", flush=True)

    jobs = [(u, cfg) for u in todo]
    t0 = last = time.perf_counter()
    n_done = n_bad = 0
    with rows_path.open("a") as fh, mp.get_context("spawn").Pool(max(1, a.workers), maxtasksperchild=50) as pool:
        it = pool.imap_unordered(run_unit, jobs, chunksize=1)
        while n_done < len(jobs):
            try:
                row = it.next(timeout=30)
            except mp.TimeoutError:
                row = None
            if row is not None:
                fh.write(json.dumps(row, default=str) + "\n")
                fh.flush()
                n_done += 1
                if not row["consistent"]:
                    n_bad += 1
                    print(f"CONSISTENCY FAILURE {row['unit_id']}: {row['failures']}", flush=True)
            now = time.perf_counter()
            if now - last >= 60 or n_done == len(jobs):
                print(f"[{now - t0:7.0f} s] {n_done}/{len(jobs)} done, {n_bad} inconsistent", flush=True)
                last = now

    # resumed runs must also fail on rows written before the restart
    bad = 0
    for line in rows_path.read_text().splitlines():
        try:
            bad += not json.loads(line)["consistent"]
        except (json.JSONDecodeError, KeyError):
            pass
    print(f"E4 {tag}: finished, {bad} inconsistent rows in {rows_path}", flush=True)
    return 2 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
