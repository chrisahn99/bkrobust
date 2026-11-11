"""E13: consistency checks on committed sets and the orientation-error budget.

Part a: the 898 real-name units recorded with no optimal set (status_split.json, classification
"none"): amenability of G0, the closed-form set, and a recount of the extensions' optimal sets.

Parts b and c: every released row (per_query_records.csv), rebuilt and checked against its hashes.
(b) closed-form committed set against the optimal sets O_D of every DAG extension of G0.
(c) on rows invalid at the truth, Corollary 1(c): the truth reverses at least r_val closure
orientations of G0; and, where Ext(CPDAG) is small, r_val recomputed by Corollary 1(c) itself.

One JSON line per unit. See PREREG.md for the registered statistics.

    python experiments/rebuttal/E13_committed_set_consistency_checks/run.py --smoke
    python experiments/rebuttal/E13_committed_set_consistency_checks/run.py --workers 8
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import multiprocessing as mp
import os
import subprocess
import sys
import time
import traceback
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
for p in (ROOT / "src", HERE):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from bkrobust.benchmarks.measure import component_of  # noqa: E402
from bkrobust.core.conventions import UNREACHED  # noqa: E402
from bkrobust.demo.evaluate import optimal_adjustment_set_dag  # noqa: E402
from bkrobust.demo.graph import MPDAG  # noqa: E402
from bkrobust.demo.meek import (  # noqa: E402
    apply_orientations,
    enumerate_dag_extensions,
    is_consistent_extension,
)
from bkrobust.gac.dag_level import is_gac_valid_dag  # noqa: E402
from bkrobust.gac.mpdag_level import clear_cache, is_gac_valid_mpdag  # noqa: E402
from bkrobust.mpdag_criterion.criterion import is_amenable  # noqa: E402
from bkrobust.mpdag_criterion.paths import possible_descendants  # noqa: E402
from vendor_audit import is_valid_adjustment_set_gac_dag  # noqa: E402
from vendor_ledger import load, rebuild_knowledge, sha  # noqa: E402
from vendor_optimal import committed_adjustment_set  # noqa: E402

CRITERION = "gac"
DATA = ROOT / "experiments" / "rebuttal" / "_data"
RELEASE_SEED = 20260913
SMALL_NETWORKS = ("asia", "sachs", "mediator", "paths", "Acid_1996", "Didelez_2010", "Shrier_2008")


def repo_commit() -> tuple[str, bool | None]:
    try:
        head = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "--", "src"],
                               capture_output=True, text=True, check=True).stdout.strip()
        return head, bool(dirty)
    except (OSError, subprocess.CalledProcessError):
        return "unknown", None


def file_sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha16(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()[:16]


# ------------------------------------------------------------------ extensions


def extensions(g: MPDAG, cap: int, time_limit: float) -> tuple[list[MPDAG], str, int]:
    """Every DAG extension of ``g``, by orienting one undirected edge at a time and Meek-closing.

    Brute force (the submitted enumerate_dag_extensions) visits 2**k orientations; this visits
    about two closures per extension, which is what makes G0s with 20-45 undirected edges and few
    extensions reachable. Leaves are checked with is_consistent_extension, and on small k the
    result is compared with brute force by the caller.
    Returns (extensions sorted by edge string, "complete" | "capped" | "timeout", bad leaves).
    """
    deadline = time.perf_counter() + time_limit
    out: list[MPDAG] = []
    bad = 0
    stack = [g]
    while stack:
        if time.perf_counter() > deadline:
            return out, "timeout", bad
        h = stack.pop()
        und = sorted(h.undirected_edges)
        if not und:
            if is_consistent_extension(h, g):
                out.append(h)
                if len(out) > cap:
                    return out, "capped", bad
            else:
                bad += 1
            continue
        a, b = und[0]
        for o in ((b, a), (a, b)):
            h2 = apply_orientations(h, [o])
            if h2 is not None:
                stack.append(h2)
    out.sort(key=lambda d: d.edge_string())
    return out, "complete", bad


_EXT: dict = {}
_KNOW: dict = {}


def knowledge(path: str) -> dict:
    """The elicitation bundle, read once per worker rather than pickled into every task."""
    if path not in _KNOW:
        _KNOW[path] = json.loads(Path(path).read_text())
    return _KNOW[path]


def cached_extensions(key, g: MPDAG, cap: int, time_limit: float, brute_max: int) -> dict:
    """Extensions of ``g`` once per worker, with the brute-force cross-check when ``g`` is small.

    Only the latest G0 and the latest CPDAG are kept: tasks are whole G0 groups, and 4,096
    extensions of a 400-node graph are a few hundred MB.
    """
    if key not in _EXT:
        for old in [k for k in _EXT if k[0] == key[0]]:
            del _EXT[old]
        t = time.perf_counter()
        exts, stop, bad = extensions(g, cap, time_limit)
        rec = {"exts": exts, "stop": stop, "bad_leaves": bad, "n": len(exts),
               "seconds": round(time.perf_counter() - t, 3), "brute_agrees": None}
        if stop == "complete" and len(g.undirected_edges) <= brute_max:
            brute = enumerate_dag_extensions(g)
            rec["brute_agrees"] = [d.edge_string() for d in brute] == [d.edge_string() for d in exts]
        _EXT[key] = rec
    return _EXT[key]


def optimal_sets(exts: list[MPDAG], x: str, y: str) -> list[frozenset]:
    return [frozenset(optimal_adjustment_set_dag(d, x, y)) for d in exts]


def summarise_sets(sets: list[frozenset], z: frozenset | None) -> dict:
    distinct = sorted({tuple(sorted(s)) for s in sets})
    agreed = list(distinct[0]) if len(distinct) == 1 else None
    return {
        "n_distinct_O_D": len(distinct),
        "O_D_agreed": agreed,
        "O_D_examples": [list(s) for s in distinct[:4]],
        "z_equals_agreed": None if z is None or agreed is None else sorted(z) == agreed,
        "z_in_some_O_D": None if z is None else tuple(sorted(z)) in distinct,
    }


def closure_orientations(g0: MPDAG, cpdag: MPDAG) -> list[tuple[str, str]]:
    """K_G0: the orientations of G0 on edges the CPDAG leaves undirected."""
    und = cpdag.undirected_edges
    return sorted((a, b) for (a, b) in g0.directed_edges if tuple(sorted((a, b))) in und)


def reversed_count(kg0: list, d: MPDAG) -> int:
    """|K_G0 \\ dir(D)| for a fully directed D: G0 orientations that D points the other way."""
    return sum(1 for (a, b) in kg0 if not d.is_directed_edge(a, b))


# ------------------------------------------------------------------ part a


def run_group_a(task: tuple) -> list[dict]:
    (cond, net), units, cfg = task
    from bkrobust.robustness import llm_survival as ls
    from bkrobust.robustness import real_survival as rs

    k, _ = ls.elicited_k(knowledge(cfg["knowledge"]), cond, net)
    _dag, cpdag = rs.load_network(net)
    g0, reason = rs.build_g0(cpdag, k)
    out = []
    ext = None
    if g0 is not None:
        # status_split enumerated by brute force; so does this recount, then the fast path is compared.
        brute = enumerate_dag_extensions(g0)
        ext = {"exts": brute, "n": len(brute)}
        # Cap at the brute-force count: anything beyond it is itself a disagreement.
        fast, stop, bad = extensions(g0, len(brute), cfg["enum_time_limit"])
        ext["fast_agrees"] = None if stop == "timeout" else (
            stop == "complete" and bad == 0
            and [d.edge_string() for d in fast] == [d.edge_string() for d in brute])
    for u in units:
        t0 = time.perf_counter()
        x, y = u["x"], u["y"]
        rec = {"part": "a", "key": f"a|{cond}|{net}|{x}|{y}", "condition": cond, "network": net,
               "x": x, "y": y, "naming": u.get("naming"), "real_condition": cond in cfg["real"],
               "released_classification": u["classification"], "criterion": CRITERION,
               "commit": cfg["commit"], "src_dirty": cfg["src_dirty"], "smoke": cfg["smoke"]}
        try:
            if g0 is None:
                rec.update(status="rebuild_mismatch", mismatch=[f"g0_none:{reason}"], hard=[])
            else:
                rec.update(unit_a(g0, x, y, ext))
        except Exception:
            rec.update(status="error", error=traceback.format_exc()[-2000:])
        rec["seconds"] = round(time.perf_counter() - t0, 3)
        out.append(rec)
    clear_cache()
    return out


def unit_a(g0: MPDAG, x: str, y: str, ext: dict) -> dict:
    hard: list[str] = []
    res = {"g0_undirected": len(g0.undirected_edges), "n_extensions": ext["n"],
           "fast_enumeration_agrees": ext["fast_agrees"]}
    if ext["fast_agrees"] is False:
        hard.append("fast enumeration differs from brute force")
    sets = optimal_sets(ext["exts"], x, y)
    res["y_possible_descendant"] = y in possible_descendants(g0, x)
    res["amenable"] = is_amenable(g0, x, y)
    cs = committed_adjustment_set(g0, x, y)
    z = cs.z
    res.update({"closed_form_verdict": cs.verdict, "z": None if z is None else sorted(z),
                "z_valid_at_g0": None if z is None else is_gac_valid_mpdag(g0, x, y, z)})
    res.update(summarise_sets(sets, z))
    res["reproduces_none"] = res["n_distinct_O_D"] > 1 or ext["n"] == 0
    if z is not None:
        per_ext = [is_gac_valid_dag(d, x, y, z) for d in ext["exts"]]
        res["z_valid_in_every_extension"] = all(per_ext)
        if res["z_valid_at_g0"] and not all(per_ext):
            hard.append("GAC-valid at G0 but invalid in an extension")
    # Amenability is necessary for any valid adjustment set in an MPDAG.
    if z is not None and not res["amenable"]:
        hard.append("a committed set on a non-amenable G0")
    res["status"] = "ok" if res["reproduces_none"] else "rebuild_mismatch"
    if not res["reproduces_none"]:
        res["mismatch"] = ["extensions agree on O_D; status_split said none"]
    res["hard"] = hard
    return res


# ------------------------------------------------------------------ parts b, c

_NETS: dict = {}


def network(models: Path, name: str):
    if name not in _NETS:
        _NETS[name] = load(models, name)
    return _NETS[name]


def run_group_bc(task: tuple) -> list[dict]:
    (net, arm, rep), rows, cfg = task
    dag, cpdag, _names = network(Path(cfg["models"]), net)
    k = rebuild_knowledge(net, arm, rep, cpdag, knowledge(cfg["knowledge"]), cfg["seed"])
    g0 = None if k is None else apply_orientations(cpdag, k)
    out = []
    for row in rows:
        t0 = time.perf_counter()
        rec = {"part": "bc", "key": bc_key(row), **{f: row[f] for f in ("network", "x", "y", "arm", "rep")},
               "criterion": CRITERION, "commit": cfg["commit"], "src_dirty": cfg["src_dirty"],
               "smoke": cfg["smoke"],
               "released": {f: row[f] for f in ("status", "committed_verdict", "z_size", "g0_undirected",
                                                "k_g0", "r_hop", "hop_method", "hop_exact",
                                                "hop_lower_bound", "z_valid_at_truth")}}
        try:
            rec.update(unit_bc(row, dag, cpdag, k, g0, cfg))
        except Exception:
            rec.update(status="error", error=traceback.format_exc()[-2000:])
        rec["seconds"] = round(time.perf_counter() - t0, 3)
        out.append(rec)
    clear_cache()
    return out


def unit_bc(row: dict, dag: MPDAG, cpdag: MPDAG, k, g0, cfg: dict) -> dict:
    mism: list[str] = []
    hard: list[str] = []
    res: dict = {"rebuild": {}}
    x, y = row["x"], row["y"]

    def expect(name, got, want):
        res["rebuild"][name] = {"got": got, "released": want}
        if str(got) != str(want):
            mism.append(name)

    expect("chat_source_sha256", sha(cpdag), row["chat_source_sha256"])
    expect("generating_dag_sha256", sha(dag), row["generating_dag_sha256"])
    if k is None or g0 is None:
        return {**res, "status": "rebuild_mismatch", "mismatch": mism + ["knowledge_or_g0"], "hard": hard}
    expect("k_sha256", sha16(json.dumps(sorted(k))), row["k_sha256"])
    expect("n_k", len(k), row["n_k"])
    expect("g0_undirected", len(g0.undirected_edges), row["g0_undirected"])
    kg0 = closure_orientations(g0, cpdag)
    expect("k_g0", len(kg0), row["k_g0"])
    cs = committed_adjustment_set(g0, x, y)
    expect("committed_verdict", cs.verdict, row["committed_verdict"])
    z = cs.z
    if z is not None:
        expect("z_size", len(z), row["z_size"])
        comp = component_of(cpdag, x)
        expect("component_size", len(comp) if comp else 0, row["component_size"])
    if mism:
        return {**res, "status": "rebuild_mismatch", "mismatch": mism, "hard": hard}
    res.update({"status": "ok", "verdict": cs.verdict, "z": None if z is None else sorted(z),
                "n_kg0": len(kg0)})

    # (b) closed form against the extensions of G0
    ext = cached_extensions(("g0", row["network"], row["arm"], row["rep"]), g0,
                            cfg["max_extensions"], cfg["enum_time_limit"], cfg["brute_max"])
    b = {"enumeration": ext["stop"], "n_extensions": ext["n"], "enum_seconds": ext["seconds"],
         "bad_leaves": ext["bad_leaves"], "brute_agrees": ext["brute_agrees"]}
    if ext["brute_agrees"] is False:
        hard.append("fast enumeration differs from brute force on G0")
    if ext["stop"] == "complete":
        if ext["n"] == 0:
            hard.append("G0 has no DAG extension")
        else:
            b.update(summarise_sets(optimal_sets(ext["exts"], x, y), z))
            if cs.verdict == "closed_form":
                b["agree"] = b["z_equals_agreed"] is True
                # The per-extension GAC is the costly leg; capped so big G0s stay in minutes.
                if ext["n"] <= cfg["max_validity_extensions"]:
                    per_ext = [is_gac_valid_dag(d, x, y, z) for d in ext["exts"]]
                    b["z_valid_in_every_extension"] = all(per_ext)
                    if not all(per_ext):
                        hard.append("GAC-valid at G0 but invalid in an extension of G0")
                else:
                    b["z_valid_in_every_extension"] = None
            if cs.verdict == "y_not_possible_descendant" and b["O_D_agreed"] != []:
                hard.append("y not a possible descendant, yet some O_D is non-empty")
    res["b"] = b
    if z is None:
        res["hard"] = hard
        return res

    # (c) the truth against the orientation-error budget
    if not is_gac_valid_mpdag(g0, x, y, z):
        hard.append("committed set invalid at G0, though read off G0")
    valid_nx = is_valid_adjustment_set_gac_dag(dag, x, y, z)
    valid_sub = is_gac_valid_dag(dag, x, y, z)
    if valid_nx != valid_sub:
        hard.append("the two DAG-level GAC implementations disagree at the truth")
    if str(int(valid_nx)) != row["z_valid_at_truth"]:
        hard.append("z_valid_at_truth not reproduced")
    c = reversed_count(kg0, dag)
    r_hop = int(row["r_hop"]) if row["r_hop"] not in ("", None) else None
    exact = row["hop_exact"] == "1"
    if r_hop is not None and r_hop != UNREACHED and exact:
        bound, bound_kind = r_hop, "exact"
    elif r_hop == UNREACHED and exact:
        bound, bound_kind = None, "unreached"  # infinity: Z valid in every extension of the CPDAG
    elif row["hop_lower_bound"] not in ("", None):
        bound, bound_kind = int(float(row["hop_lower_bound"])), "lower_bound"
    else:
        bound, bound_kind = None, "upper_only"  # an inexact finite r_hop bounds r_val from above only
    cc = {"valid_at_truth": valid_nx, "reversed_by_truth": c, "bound": bound, "bound_kind": bound_kind}
    if not valid_nx:
        if bound_kind == "upper_only":
            cc["violation"] = None
        else:
            cc["violation"] = bound_kind == "unreached" or c < bound
            cc["slack"] = None if bound is None else c - bound
            if cc["violation"]:
                hard.append(f"Corollary 1(c): truth reverses {c} orientations, bound {bound_kind} {bound}")
    res["c"] = cc

    # Secondary: r_val by Corollary 1(c) itself, over every extension of the CPDAG.
    ce = cached_extensions(("chat", row["network"]), cpdag, cfg["max_chat_extensions"],
                           cfg["enum_time_limit"], cfg["brute_max"])
    cor = {"enumeration": ce["stop"], "n_extensions": ce["n"], "brute_agrees": ce["brute_agrees"]}
    if ce["brute_agrees"] is False:
        hard.append("fast enumeration differs from brute force on the CPDAG")
    if ce["stop"] == "complete" and ce["n"]:
        failing = [reversed_count(kg0, d) for d in ce["exts"] if not is_gac_valid_dag(d, x, y, z)]
        r_enum = min(failing) if failing else None
        cor["r_enum"] = r_enum
        if bound_kind == "exact":
            cor["match"] = r_enum == bound
        elif bound_kind == "unreached":
            cor["match"] = r_enum is None
        elif bound_kind == "lower_bound":
            cor["match"] = r_enum is None or r_enum >= bound
        else:
            cor["match"] = None
        if cor["match"] is False:
            hard.append(f"r_val by enumeration {r_enum} vs released {bound_kind} {bound}")
    res["corollary"] = cor
    res["hard"] = hard
    return res


# ------------------------------------------------------------------ driver


def bc_key(r: dict) -> str:
    return "bc|" + "|".join((r["network"], r["x"], r["y"], r["arm"], str(r["rep"])))


def done_keys(out: Path) -> set[str]:
    keys = set()
    for f in sorted(out.glob("rows*.jsonl")):
        for line in f.open():
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue  # a line cut by a kill; the unit is rerun
            if rec.get("status") != "error":
                keys.add(rec["key"])
    return keys


def smoke_rows(rows: list[dict]) -> list[dict]:
    """Two rows per category on small networks: every verdict, and invalid-at-truth rows."""
    cats = {
        "ok_valid_enum": lambda r: r["status"] == "ok" and r["z_valid_at_truth"] == "1" and r["g0_undirected"] != "0",
        "ok_invalid": lambda r: r["status"] == "ok" and r["z_valid_at_truth"] == "0",
        "not_amenable": lambda r: r["status"] == "not_amenable",
        "y_not_pd": lambda r: r["status"] == "y_not_possible_descendant",
        "degen": lambda r: r["arm"] == "D_DEGEN" and r["status"] == "ok",
        "rand": lambda r: r["arm"] == "D_RAND" and r["status"] == "ok",
    }
    picked, seen = [], defaultdict(int)
    for r in rows:
        if r["network"] not in SMALL_NETWORKS:
            continue
        for name, test in cats.items():
            if seen[name] < 2 and test(r):
                seen[name] += 1
                picked.append(r)
                break
    return picked


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--records", type=Path, default=DATA / "per_query_records.csv")
    ap.add_argument("--knowledge", type=Path, default=ROOT / "results/elicit/knowledge.json",
                    help="byte-identical to release/knowledge_states.json")
    ap.add_argument("--status-split", type=Path,
                    default=ROOT / "results/axis_robustness_llm_v2/status_split.json")
    ap.add_argument("--models", type=Path, default=ROOT / "results/axisa3/networks/example_models")
    ap.add_argument("--out", type=Path, default=ROOT / "results/rebuttal/E13")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--parts", default="a,bc", help="comma list of a, bc")
    ap.add_argument("--smoke", action="store_true", help="a dozen rows on small networks, < 60 s")
    ap.add_argument("--seed", type=int, default=RELEASE_SEED, help="D_RAND seed of the released sweep")
    ap.add_argument("--max-extensions", type=int, default=4096, help="cap on Ext(G0) for part b")
    ap.add_argument("--max-validity-extensions", type=int, default=512,
                    help="check Z in every extension of G0 only up to this many extensions")
    ap.add_argument("--max-chat-extensions", type=int, default=256, help="cap on Ext(CPDAG) for the r_val recount")
    ap.add_argument("--enum-time-limit", type=float, default=600.0, help="seconds per enumeration")
    ap.add_argument("--brute-max", type=int, default=10, help="brute-force cross-check up to this many undirected edges")
    args = ap.parse_args(argv)
    if not 0 <= args.shard < args.n_shards:
        ap.error("--shard must be in [0, --n-shards)")
    parts = {p for p in args.parts.split(",") if p}
    if args.smoke:
        args.enum_time_limit = min(args.enum_time_limit, 5.0)
        args.max_extensions = min(args.max_extensions, 256)
        if args.out == ROOT / "results/rebuttal/E13":
            args.out = ROOT / "results/rebuttal/E13_smoke"

    commit, dirty = repo_commit()
    for name in ("records", "knowledge", "status_split", "models", "out"):
        setattr(args, name, getattr(args, name).resolve())
    # rs.load_network reads the network cache by a relative path (results/axisa3/networks),
    # so the run works from the repository root; spawned workers inherit this directory.
    os.chdir(ROOT)

    from bkrobust.robustness import llm_survival as ls
    bundle, bundle_sha = ls.load_bundle(args.knowledge)
    real = set(ls.panel_conditions(bundle, naming=ls.REAL_NAMING))
    base = {"commit": commit, "src_dirty": dirty, "smoke": args.smoke,
            "max_extensions": args.max_extensions, "max_chat_extensions": args.max_chat_extensions,
            "enum_time_limit": args.enum_time_limit, "brute_max": args.brute_max,
            "max_validity_extensions": args.max_validity_extensions}

    groups: list[tuple] = []  # (sort key, function, task)
    n_sel = {"a": 0, "bc": 0}
    if "a" in parts:
        split = json.loads(args.status_split.read_text())
        units = sorted((u for u in split["units"] if u["classification"] == "none"),
                       key=lambda u: (u["condition"] not in real, u["condition"], u["network"], u["x"], u["y"]))
        if args.smoke:
            first = (units[0]["condition"], units[0]["network"])
            units = [u for u in units if (u["condition"], u["network"]) == first][:3]
        n_sel["a"] = len(units)
        by: dict = defaultdict(list)
        for u in units:
            by[(u["condition"], u["network"])].append(u)
        cfg_a = {**base, "knowledge": str(args.knowledge), "real": real}
        for g, us in sorted(by.items()):
            groups.append((("a",) + g, run_group_a, (g, us, cfg_a), [f"a|{g[0]}|{g[1]}|{u['x']}|{u['y']}" for u in us]))
    if "bc" in parts:
        rows = list(csv.DictReader(args.records.open()))
        rows.sort(key=bc_key)
        if args.smoke:
            rows = smoke_rows(rows)
        n_sel["bc"] = len(rows)
        by = defaultdict(list)
        for r in rows:
            by[(r["network"], r["arm"], int(r["rep"]))].append(r)
        cfg_bc = {**base, "knowledge": str(args.knowledge), "models": str(args.models), "seed": args.seed}
        for g, rs_ in sorted(by.items()):
            groups.append((("bc",) + tuple(map(str, g)), run_group_bc, (g, rs_, cfg_bc), [bc_key(r) for r in rs_]))

    # Shard by group, so that a G0 and its extensions are built once.
    groups.sort(key=lambda t: t[0])
    groups = [t for i, t in enumerate(groups) if i % args.n_shards == args.shard]
    args.out.mkdir(parents=True, exist_ok=True)
    tag = "" if args.n_shards == 1 else f".shard{args.shard}of{args.n_shards}"
    config = {"experiment": "E13", "criterion": CRITERION, "commit": commit, "src_dirty": dirty,
              "smoke": args.smoke, "parts": sorted(parts), "n_selected": n_sel,
              "shard": args.shard, "n_shards": args.n_shards, "seed": args.seed,
              "max_extensions": args.max_extensions, "max_chat_extensions": args.max_chat_extensions,
              "enum_time_limit": args.enum_time_limit, "brute_max": args.brute_max,
              "max_validity_extensions": args.max_validity_extensions,
              "workers": args.workers,
              "records_sha256": file_sha(args.records) if "bc" in parts else None,
              "knowledge_sha256": bundle_sha,
              "status_split_sha256": file_sha(args.status_split) if "a" in parts else None}
    (args.out / f"config{tag}.json").write_text(json.dumps(config, indent=1))

    skip = done_keys(args.out)
    tasks = []
    for _k, fn, (g, items, cfg), keys in groups:
        keep = [it for it, key in zip(items, keys) if key not in skip]
        if keep:
            tasks.append((fn, (g, keep, cfg)))
    n_units = sum(len(t[1][1]) for t in tasks)
    print(f"[E13] selected {n_sel}; shard {args.shard}/{args.n_shards}: {len(tasks)} groups, "
          f"{n_units} units to run, {len(skip)} already done; commit {commit[:12]}"
          f"{' (src dirty)' if dirty else ''}", flush=True)

    n_hard = n_err = n_mis = done = 0
    t0 = last = time.perf_counter()
    with (args.out / f"rows{tag}.jsonl").open("a") as fh:
        pool = mp.get_context("spawn").Pool(max(1, args.workers))
        pending = [pool.apply_async(fn, (task,)) for fn, task in tasks]
        open_ = list(range(len(pending)))
        while open_:
            still = []
            for i in open_:
                if not pending[i].ready():
                    still.append(i)
                    continue
                try:
                    recs = pending[i].get()
                except Exception:  # a whole group failed; its units are retried on resume
                    n_err += len(tasks[i][1][1])
                    print(f"[E13] GROUP ERROR {tasks[i][1][0]}: {traceback.format_exc().strip().splitlines()[-1]}",
                          flush=True)
                    continue
                for rec in recs:
                    fh.write(json.dumps(rec, default=str) + "\n")
                    done += 1
                    n_err += rec.get("status") == "error"
                    n_mis += rec.get("status") == "rebuild_mismatch"
                    if rec.get("hard"):
                        n_hard += 1
                        print(f"[E13] HARD {rec['key']}: {rec['hard']}", flush=True)
                    if rec.get("status") == "error":
                        print(f"[E13] ERROR {rec['key']}: {rec['error'].strip().splitlines()[-1]}", flush=True)
                fh.flush()
            open_ = still
            now = time.perf_counter()
            if now - last >= 30 or not open_:
                last = now
                print(f"[E13] {done}/{n_units} units, {len(open_)} groups open, rebuild mismatches {n_mis}, "
                      f"hard {n_hard}, errors {n_err}, {now - t0:.0f} s", flush=True)
            if open_:
                time.sleep(0.5)
        pool.close()
        pool.join()
    if n_hard:
        print(f"[E13] CONSISTENCY FAILURE on {n_hard} units", flush=True)
        return 1
    if n_err:
        print(f"[E13] {n_err} units raised; rerun to retry them", flush=True)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
