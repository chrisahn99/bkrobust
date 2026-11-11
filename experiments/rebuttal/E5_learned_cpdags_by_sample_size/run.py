#!/usr/bin/env python3
"""E5: validity at truth and the radius on learned CPDAGs, by sample size.

Phase A recomputes every row of the stored P-EST panel (results/pest/rows.csv) on the reference
learned CPDAG with the submitted code and compares field by field. Phase B scores every stored
learned graph (PC and GES, n in 200/1000/5000, every seed) and the true CPDAG of each network, on
the same 80 frame queries and the same knowledge arms. See PREREG.md.

    .venv/bin/python experiments/rebuttal/E5_learned_cpdags_by_sample_size/run.py --smoke
    .venv/bin/python experiments/rebuttal/E5_learned_cpdags_by_sample_size/run.py --workers 6

Exit 0 when every consistency check holds, 2 when one fails (rows are still written), 1 on error.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import multiprocessing as mp
import os
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path

# one BLAS thread per worker, or six workers oversubscribe sixteen cores
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
for _p in (REPO / "src", HERE):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import numpy as np  # noqa: E402

from bkrobust.benchmarks.describe import parse_file, to_mpdag  # noqa: E402
from bkrobust.benchmarks.measure import component_of, separation  # noqa: E402
from bkrobust.core.conventions import UNREACHED  # noqa: E402
from bkrobust.demo.example import dag_to_cpdag  # noqa: E402
from bkrobust.demo.graph import MPDAG  # noqa: E402
from bkrobust.demo.meek import apply_orientations, is_consistent_extension  # noqa: E402
from bkrobust.gac.mpdag_level import is_gac_valid_mpdag  # noqa: E402
from bkrobust.hybrid import DEFAULT_SEARCH_BUDGET, breakdown_radius  # noqa: E402
from bkrobust.mpdag_criterion import is_valid_mpdag  # noqa: E402
from bkrobust.search.exact import retractable_edges  # noqa: E402
from bkrobust.search.exact_fast import radius_local_up_fast  # noqa: E402
from vendor_bkrobust import (  # noqa: E402
    EXHAUSTIVE_MAX_RETRACTABLE,
    LADDER_MAX_NODES,
    NO_RETRACTABLE_EDGES,
    RUNGS,
    TIER,
    audit_row,
    breaking_depths,
    chance_level,
    claim_lower_bound,
    committed_adjustment_set,
    graph_from_json,
    sha,
)

NETWORKS = ("ecoli70", "arth150", "magic-niab", "magic-irri")
N_GRID = (200, 1000, 5000)
CRITERION = "gac"
DRAND_SEED = 20260913  # the stored sweep's --seed, so D_RAND reproduces
STORED_RAND_REPS = 3
SMOKE_NETWORKS = ("magic-niab", "magic-irri")
SMOKE_GRAPHS = {
    "magic-niab": [("pc", 200, 0), ("pc", 5000, 0), ("ges", 1000, 0), ("true_cpdag", None, None)],
    "magic-irri": [("pc", 5000, 0), ("true_cpdag", None, None)],
}
INF = float("inf")

# Phase A: fields the stored sweep computed deterministically, compared exactly
COMPARE_FIELDS = (
    "tier", "true_dag_on_path", "chat_source_sha256", "generating_dag_sha256", "k_sha256",
    "n_k", "status", "committed_verdict", "z_size", "g0_undirected", "k_g0", "separation",
    "separation_status", "component_size", "r_claim", "r_claim_status", "r_id", "phi_1",
    "witness", "n_wrong_claims", "n_omitted", "z_valid_at_truth", "z_backdoor_at_truth",
    "severity", "bucket_B2_free_rule", "bucket_B4_kg0", "bucket_B5_k", "bucket_B7_claim",
    "truth_status", "truth_extends_g0", "n_claims_on_absent_edges", "y_true_descendant",
    "r_cert_claim", "cpdag_source",
)


def stable_seed(name: str) -> int:
    """experiments/questionnaire_build.py::stable_seed, copied to avoid importing that module."""
    return int(hashlib.sha256(name.encode()).hexdigest()[:8], 16)


def canon(e) -> tuple[str, str]:
    return tuple(sorted(e))


def k_hash(k: list) -> str:
    return hashlib.sha256(json.dumps(sorted(k)).encode()).hexdigest()[:16]


def file_sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def repo_commit() -> str:
    try:
        out = subprocess.run(
            ["git", "-C", str(REPO), "rev-parse", "HEAD"], capture_output=True, text=True, timeout=20
        )
        return out.stdout.strip() or "unknown"
    except (OSError, subprocess.TimeoutExpired):
        return "unknown"


# --- inputs ------------------------------------------------------------------------


def load_inputs(data_dir: Path, models_dir: Path, networks: tuple[str, ...]) -> dict:
    """Everything a worker needs, read once. Fails loudly on a missing or mismatched file."""
    chat = json.loads((data_dir / "chat_e.json").read_text())
    know = json.loads((data_dir / "knowledge.json").read_text())
    frame = [json.loads(line) for line in (data_dir / "frame.jsonl").open() if line.strip()]
    learned = []
    for name in ("graphs.jsonl", "graphs_retry.jsonl"):
        p = data_dir / name
        if p.exists():
            learned += [json.loads(line) for line in p.open() if line.strip()]
    stored: dict[tuple, dict] = {}
    with (data_dir / "rows.csv").open() as fh:
        for r in csv.DictReader(fh):
            stored[(r["network"], r["arm"], int(r["rep"]), r["x"], r["y"])] = r
    dags, names = {}, {}
    for net in networks:
        path = models_dir / f"{net}.json"
        parsed = parse_file(path, file_sha(path))
        dags[net] = to_mpdag(parsed)
        names[net] = dict(parsed.name_map)
    graphs: dict[str, dict] = {}
    for net in networks:
        e = chat["networks"][net]
        g = graph_from_json(e)
        if sha(g) != e["sha256"]:
            raise SystemExit(f"{net}: chat_e.json reference CPDAG hash mismatch")
        graphs[f"{net}|reference"] = {"net": net, "method": "pc", "n": chat["meta"]["n"],
                                      "seed": e["seed"], "kind": "reference", "g": g}
        cp = dag_to_cpdag(dags[net])
        graphs[f"{net}|true_cpdag|true|-"] = {"net": net, "method": "true_cpdag", "n": None,
                                              "seed": None, "kind": "true_cpdag", "g": cp}
    for e in learned:
        if e["network"] not in networks:
            continue
        gid = f"{e['network']}|{e['method']}|{e['n']}|{e['seed']}"
        if gid in graphs:
            raise SystemExit(f"duplicate learned graph {gid}")
        graphs[gid] = {"net": e["network"], "method": e["method"], "n": e["n"], "seed": e["seed"],
                       "kind": "learned", "g": graph_from_json(e)}
    return {"chat": chat, "know": know, "frame": frame, "stored": stored, "dags": dags,
            "names": names, "graphs": graphs}


# --- knowledge arms ----------------------------------------------------------------


def greedy_consistent(g: MPDAG, claims: list) -> tuple[list, int]:
    """Keep claims in order, skipping any whose addition makes the Meek closure fail."""
    kept, dropped = [], 0
    for e in claims:
        if apply_orientations(g, [*kept, e]) is None:
            dropped += 1
        else:
            kept.append(e)
    return kept, dropped


def truthful(dag: MPDAG, edges) -> list:
    return [(a, b) if dag.is_directed_edge(a, b) else (b, a) for (a, b) in edges if dag.has_edge(a, b)]


def plans_stored(inp: dict, net: str, rand_reps: int) -> list[tuple]:
    """The stored sweep's arms on the reference CPDAG, exactly as pest_panel.stage_sweep built them."""
    dag, g = inp["dags"][net], inp["graphs"][f"{net}|reference"]["g"]
    base = [tuple(e) for e in inp["know"]["D_LLM"]["networks"].get(net, {}).get("k", [])]
    plans = [("D_LLM", 0, base, "false", {})]
    if base:
        for rep in range(rand_reps):
            kk = chance_level(g, base, np.random.default_rng([DRAND_SEED, stable_seed(net), rep]))
            if kk is not None:
                plans.append(("D_RAND", rep, kk, "false", {}))
    plans.append(("D_DEGEN", 0, [], "false", {}))
    plans.append(("A_TRUE", 0, truthful(dag, base), "via_K", {}))
    plans.append(("A_ORACLE", 0, truthful(dag, sorted(g.undirected_edges)), "via_K", {}))
    return plans


def plans_transfer(inp: dict, gid: str, rand_reps: int) -> list[tuple]:
    """The same arms on another graph: the supplier's claims restricted to the graph's open edges."""
    info = inp["graphs"][gid]
    net, g, dag = info["net"], info["g"], inp["dags"][info["net"]]
    stored = [tuple(e) for e in inp["know"]["D_LLM"]["networks"].get(net, {}).get("k", [])]
    on_graph = [e for e in stored if canon(e) in g.undirected_edges]
    base, dropped = greedy_consistent(g, on_graph)
    note = {"n_llm_stored": len(stored), "n_llm_off_graph": len(stored) - len(on_graph),
            "n_llm_dropped": dropped}
    plans = [("D_LLM", 0, base, "false", note)]
    if base:
        for rep in range(rand_reps):
            kk = chance_level(g, base, np.random.default_rng([DRAND_SEED, stable_seed(net), rep]))
            if kk is not None:
                plans.append(("D_RAND", rep, kk, "false", note))
    plans.append(("D_DEGEN", 0, [], "false", {}))
    plans.append(("A_TRUE", 0, truthful(dag, base), "via_K", note))
    oracle = truthful(dag, sorted(g.undirected_edges))
    plans.append(("A_ORACLE", 0, oracle, "via_K", {}))
    # the stored A_ORACLE is unusable whenever one truthful claim contradicts a wrongly learned
    # arrow; the repaired arm keeps the truthful claims that close, so every graph contributes
    kept, dropped_o = greedy_consistent(g, oracle)
    plans.append(("A_ORACLE_R", 0, kept, "via_K", {"n_oracle_dropped": dropped_o}))
    return plans


# --- the radius --------------------------------------------------------------------


class WallError(Exception):
    pass


def _fire(signum, frame):
    raise WallError


def with_wall(seconds: float, fn, *a, **kw):
    """(value, timed_out). SIGALRM, so it only interrupts Python code; CP-SAT has its own limit."""
    signal.signal(signal.SIGALRM, _fire)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        return fn(*a, **kw), False
    except WallError:
        return None, True
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


def hop_radius(cpdag, g0, x, y, z, criterion, n_retractable, radius_limit, ladder_max_nodes,
               use_criterion=True) -> dict:
    """The submitted breakdown_radius; above ladder_max_nodes the same dispatch without the ladder.

    The ladder's model is built over every node triple and does not return on graphs of a hundred
    nodes, which is why the stored sweep turned it off above 60 nodes; the submitted function has
    no such switch, so the no-ladder branch repeats its first three legs with the same calls.
    """
    budget = n_retractable if n_retractable <= EXHAUSTIVE_MAX_RETRACTABLE else DEFAULT_SEARCH_BUDGET
    if len(cpdag.nodes) <= ladder_max_nodes:
        res = breakdown_radius(cpdag, None, x, y, z, g0=g0, search_budget=budget,
                               use_criterion=use_criterion, criterion=criterion,
                               time_limit_s=radius_limit)
        r, method, exact = res.radius, res.method, bool(res.exact)
    else:
        if criterion == "gac":
            def fails(g):
                return not is_gac_valid_mpdag(g, x, y, z)
        else:
            def fails(g):
                return not is_valid_mpdag(g, x, y, z)
        if not fails(cpdag):
            r, method, exact = UNREACHED, "top_state", True
        else:
            found = radius_local_up_fast(cpdag, g0, fails, max_depth=budget)
            if found.radius != UNREACHED:
                r, method, exact = found.radius, ("degenerate" if found.radius == 0 else "local_up_fast"), True
            elif found.exact:
                r, method, exact = UNREACHED, "local_up_fast", True
            else:
                r, method, exact = UNREACHED, "local_up_fast_budget", False
    # an inexact answer still certifies that every retraction up to the search budget was clean
    lower = None if exact else budget + 1
    return {"r": r, "method": method, "exact": int(exact), "lower": lower, "budget": budget}


def interval(r, exact, lower) -> tuple[float, float]:
    """A radius answer as the set of values it allows; UNREACHED is +inf."""
    if r is None:
        return (1.0, INF)
    r = int(r)
    if int(exact):
        return (INF, INF) if r == UNREACHED else (float(r), float(r))
    if r == UNREACHED:
        return (float(lower) if lower not in (None, "") else 1.0, INF)
    return (1.0, INF)  # an inexact finite ladder value bounds nothing we rely on


def overlap(a: tuple, b: tuple) -> bool:
    return max(a[0], b[0]) <= min(a[1], b[1])


# --- one row -----------------------------------------------------------------------


def score_row(cfg: dict, inp: dict, block: dict, k: list, g0, x: str, y: str) -> dict:
    gid, info = block["gid"], inp["graphs"][block["gid"]]
    net, cpdag = info["net"], info["g"]
    dag, names = inp["dags"][net], inp["names"][net]
    rec: dict = {
        "phase": block["phase"], "gid": gid, "network": net, "method": info["method"],
        "n": info["n"], "seed": info["seed"], "graph_kind": info["kind"], "tier": TIER.get(net, "T1"),
        "x": x, "y": y, "arm": block["arm"], "rep": block["rep"],
        "true_dag_on_path": block["provenance"], "chat_source_sha256": sha(cpdag),
        "generating_dag_sha256": sha(dag), "k_sha256": k_hash(k), "n_k": len(k),
        "cpdag_source": "pc_fisherz_n5000" if info["kind"] == "reference" else gid,
        "criterion": CRITERION, "repo_commit": cfg["repo_commit"], **block["note"], **block["graph_stats"],
    }
    t0 = time.perf_counter()
    fails: list[str] = []
    # truth fields that do not need a committed set, as the stored sweep wrote them on every row
    absent = sum(1 for (a, b) in k if not dag.has_edge(a, b))
    y_desc = y in dag.descendants(x)
    flags = (["claim_on_absent_edge"] if absent else []) + ([] if y_desc else ["y_not_true_descendant"])
    rec["truth_status"] = ";".join(flags) or "ok"
    rec["n_claims_on_absent_edges"] = absent
    rec["y_true_descendant"] = int(y_desc)
    rec["truth_extends_g0"] = int(g0 is not None and is_consistent_extension(dag, g0))
    if g0 is None:
        rec["status"] = "knowledge_inconsistent"
        return finish(rec, fails, t0)

    # ---- observables only, from here to the audit
    rec["g0_undirected"] = len(g0.undirected_edges)
    rec["k_g0"] = sum(1 for (a, b) in g0.directed_edges if canon((a, b)) in cpdag.undirected_edges)
    cs = committed_adjustment_set(g0, x, y)
    rec["committed_verdict"] = cs.verdict
    if cs.z is None:
        rec["status"] = cs.verdict
        return finish(rec, fails, t0)
    z = cs.z
    rec["status"] = "ok"
    rec["z_size"] = len(z)
    sep, sep_status = separation(cpdag, x, z)
    comp = component_of(cpdag, x)
    rec.update(separation=sep, separation_status=sep_status, component_size=len(comp) if comp else 0)
    n_ret = len(retractable_edges(cpdag, g0))
    rec["n_retractable"] = n_ret

    def hop(criterion, use_criterion=True):
        if not n_ret:
            return {"r": UNREACHED, "method": NO_RETRACTABLE_EDGES, "exact": 1, "lower": None, "budget": 0}, False
        return with_wall(cfg["wall_s"], hop_radius, cpdag, g0, x, y, z, criterion, n_ret,
                         cfg["radius_limit"], cfg["ladder_max_nodes"], use_criterion)

    h, h_to = hop("gac")
    rec["hop_wall_timeout"] = int(h_to)
    if h is not None:
        rec.update(r_hop=h["r"], hop_method=h["method"], hop_exact=h["exact"],
                   hop_lower_bound=h["lower"], hop_search_budget=h["budget"])
    if block["phase"] == "A":
        b, b_to = hop("backdoor")
        rec["bd_wall_timeout"] = int(b_to)
        if b is not None:
            rec.update(bd_r_hop=b["r"], bd_hop_method=b["method"], bd_hop_exact=b["exact"],
                       bd_hop_lower_bound=b["lower"])
    # second method: the same radius by enumerating DAG extensions, where that is affordable
    if (h is not None and len(cpdag.undirected_edges) <= cfg["crosscheck_max_undirected"]
            and len(cpdag.nodes) <= cfg["ladder_max_nodes"]):
        e, e_to = hop("gac", use_criterion=False)
        rec["xc_wall_timeout"] = int(e_to)
        if e is not None:
            rec.update(xc_r_hop=e["r"], xc_hop_exact=e["exact"])
            if h["exact"] and e["exact"] and h["r"] != e["r"]:
                fails.append(f"criterion_vs_enumeration:{h['r']}!={e['r']}")

    def claims():
        r_claim, status, r_id, witness = breaking_depths(cpdag, k, x, y, z)
        n_break = sum(
            1 for i in range(len(k))
            if (gi := apply_orientations(cpdag, [e for j, e in enumerate(k) if j != i])) is not None
            and not is_gac_valid_mpdag(gi, x, y, z)
        )
        return r_claim, status, r_id, witness, n_break

    c, c_to = with_wall(cfg["wall_s"], claims)
    rec["claim_wall_timeout"] = int(c_to)
    if c is not None:
        r_claim, r_status, r_id, witness, n_break = c
        rec.update(
            r_claim=r_claim, r_claim_status=r_status, r_id=r_id,
            phi_1=round(n_break / len(k), 4) if k else None,
            witness=" and ".join(f"{names.get(a, a)}->{names.get(b, b)}" for a, b in witness),
        )
    floor = claim_lower_bound(c[0], c[1]) if c is not None else None
    cert = {
        "B2_free_rule": sep if sep_status == "measured" else 1,
        "B4_kg0": rec["k_g0"],
        "B5_k": len(k),
        "B6_hop": None if h is None else (h["lower"] if not h["exact"] else (None if h["r"] == UNREACHED else h["r"])),
        "B7_claim": floor,
    }

    # ---- the truth
    rec.update(audit_row(dag, cpdag, k, x, y, z, cert))
    if h is None:
        rec["bucket_B6_hop"] = "wall_timeout"
    if c is None:
        rec["bucket_B7_claim"] = "wall_timeout"
    else:
        rec["r_cert_claim"] = UNREACHED if floor is None else floor
    # if the truth is an extension of G0, a set valid in every extension of G0 is valid in it
    if rec["truth_extends_g0"] and not rec["z_valid_at_truth"]:
        fails.append("truth_extends_g0_but_set_invalid_at_truth")
    return finish(rec, fails, t0)


def finish(rec: dict, fails: list, t0: float) -> dict:
    rec["seconds"] = round(time.perf_counter() - t0, 3)
    rec["consistency_fail"] = fails
    return rec


def norm(v):
    if v is None or v == "":
        return None
    if isinstance(v, bool):
        return float(v)
    try:
        return float(v)
    except (TypeError, ValueError):
        return str(v)


def compare_stored(rec: dict, stored: dict | None) -> list[str]:
    """Phase A: every deterministic field equal, the back-door radius compatible."""
    if stored is None:
        return ["stored_row_missing"]
    out = []
    timed = rec.get("claim_wall_timeout") or rec.get("hop_wall_timeout") or rec.get("bd_wall_timeout")
    for f in COMPARE_FIELDS:
        if timed and f in ("r_claim", "r_claim_status", "r_id", "phi_1", "witness", "bucket_B7_claim", "r_cert_claim"):
            continue
        if norm(rec.get(f)) != norm(stored.get(f)):
            out.append(f"{f}:{rec.get(f)!r}!={stored.get(f)!r}")
    if rec.get("status") == "ok" and not rec.get("bd_wall_timeout"):
        new = interval(rec.get("bd_r_hop"), rec.get("bd_hop_exact"), rec.get("bd_hop_lower_bound"))
        old = interval(stored.get("r_hop") or None, stored.get("hop_exact") or 0, stored.get("hop_lower_bound"))
        if not overlap(new, old):
            out.append(f"r_hop_backdoor:{new}!~{old}")
    return out


# --- workers -----------------------------------------------------------------------

_INP: dict = {}
_CFG: dict = {}


def init_worker(cfg: dict) -> None:
    _CFG.update(cfg)
    _INP.update(load_inputs(Path(cfg["data_dir"]), Path(cfg["models_dir"]), tuple(cfg["networks"])))


def run_block(block: dict) -> list[dict]:
    cpdag = _INP["graphs"][block["gid"]]["g"]
    k = [tuple(e) for e in block["k"]]
    g0 = apply_orientations(cpdag, k)
    rows = []
    for x, y in block["pairs"]:
        rec = score_row(_CFG, _INP, block, k, g0, x, y)
        if block["phase"] == "A":
            st = _INP["stored"].get((rec["network"], rec["arm"], rec["rep"], x, y))
            rec["consistency_fail"] = rec["consistency_fail"] + compare_stored(rec, st)
        rec["unit"] = unit_key(block, x, y)
        rows.append(rec)
    return rows


def unit_key(block: dict, x: str, y: str) -> str:
    return f"{block['phase']}|{block['gid']}|{block['arm']}|{block['rep']}|{x}|{y}"


# --- planning ----------------------------------------------------------------------


def graph_stats(dag: MPDAG, g: MPDAG) -> dict:
    ts, gs = dag.skeleton(), g.skeleton()
    return {"g_undirected": len(g.undirected_edges), "g_directed": len(g.directed_edges),
            "g_skeleton_missing": len(ts - gs), "g_skeleton_extra": len(gs - ts),
            "g_class_correct": int(is_consistent_extension(dag, g))}


def make_blocks(inp: dict, args, smoke_pairs: dict | None) -> list[dict]:
    pairs_by_net: dict[str, list] = {}
    for r in inp["frame"]:
        pairs_by_net.setdefault(r["network"], []).append((r["x"], r["y"]))
    if smoke_pairs is not None:
        pairs_by_net = smoke_pairs
    blocks = []
    nets = [n for n in args.networks if n in pairs_by_net]
    if args.phase in ("A", "all"):
        for net in nets:
            gid = f"{net}|reference"
            gs = graph_stats(inp["dags"][net], inp["graphs"][gid]["g"])
            reps = 1 if args.smoke else STORED_RAND_REPS
            for arm, rep, k, prov, note in plans_stored(inp, net, reps):
                blocks.append({"phase": "A", "gid": gid, "arm": arm, "rep": rep, "k": k,
                               "provenance": prov, "note": note, "graph_stats": gs,
                               "pairs": pairs_by_net[net]})
    if args.phase in ("B", "all"):
        for gid, info in sorted(inp["graphs"].items()):
            if info["kind"] == "reference" or info["net"] not in nets:
                continue
            if args.smoke and (info["method"], info["n"], info["seed"]) not in SMOKE_GRAPHS.get(info["net"], []):
                continue
            if info["kind"] == "learned" and (info["n"] not in args.n_grid or info["method"] not in args.methods):
                continue
            gs = graph_stats(inp["dags"][info["net"]], info["g"])
            for arm, rep, k, prov, note in plans_transfer(inp, gid, args.rand_reps):
                blocks.append({"phase": "B", "gid": gid, "arm": arm, "rep": rep, "k": k,
                               "provenance": prov, "note": note, "graph_stats": gs,
                               "pairs": pairs_by_net[info["net"]]})
    size = {n: len(inp["dags"][n].nodes) for n in inp["dags"]}
    blocks.sort(key=lambda b: (b["phase"], b["gid"], b["arm"], b["rep"]))
    mine = [b for i, b in enumerate(blocks) if i % args.n_shards == args.shard]
    mine.sort(key=lambda b: (-size[inp["graphs"][b["gid"]]["net"]], b["phase"], b["gid"], b["arm"], b["rep"]))
    return mine


def smoke_selection(inp: dict) -> dict:
    """Per smoke network, two frame pairs: one whose stored rows reach the slowest radius leg
    (the SAT ladder on magic-niab, the budgeted search on magic-irri), and one not certifiable."""
    out = {}
    st = inp["stored"]
    arms = ("D_LLM", "D_RAND", "D_DEGEN", "A_TRUE", "A_ORACLE")
    for net in SMOKE_NETWORKS:
        pairs = [(r["x"], r["y"]) for r in inp["frame"] if r["network"] == net]

        def legs(p):
            return {st.get((net, a, 0, *p), {}).get("hop_method") for a in arms}

        ok = next((p for leg in ("e1_ladder", "local_up_fast_budget", "local_up_fast")
                   for p in pairs if leg in legs(p)), pairs[0])
        bad = next((p for p in pairs if st.get((net, "D_LLM", 0, *p), {}).get("status") != "ok"), pairs[-1])
        out[net] = [ok, bad] if bad != ok else [ok]
    return out


def done_units(out: Path) -> tuple[set, list]:
    """Units already written by any shard, and every consistency failure among them."""
    done, fails = set(), []
    for p in sorted(out.glob("rows.*.jsonl")):
        for line in p.open():
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue  # a line cut by a crash; its unit is recomputed
            done.add(r["unit"])
            if r.get("consistency_fail"):
                fails.append((r["unit"], r["consistency_fail"]))
    return done, fails


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--smoke", action="store_true", help="two networks, two pairs, a few graphs")
    ap.add_argument("--out", default="results/rebuttal/E5", help="output directory (relative to the repo root)")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--phase", choices=["A", "B", "all"], default="all")
    ap.add_argument("--data-dir", default="experiments/rebuttal/_data/pest",
                    help="graphs.jsonl, graphs_retry.jsonl, knowledge.json, frame.jsonl, chat_e.json, rows.csv")
    ap.add_argument("--models-dir", default="results/axisa3/networks/example_models")
    ap.add_argument("--networks", default=",".join(NETWORKS))
    ap.add_argument("--n-grid", default="200,1000,5000")
    ap.add_argument("--methods", default="pc,ges")
    ap.add_argument("--rand-reps", type=int, default=3)
    ap.add_argument("--radius-limit", type=float, default=5.0, help="ladder per-rung limit, s (stored: 5)")
    ap.add_argument("--ladder-max-nodes", type=int, default=LADDER_MAX_NODES)
    ap.add_argument("--wall-s", type=float, default=120.0, help="per-row wall for each radius computation")
    ap.add_argument("--crosscheck-max-undirected", type=int, default=8)
    args = ap.parse_args(argv)
    if args.workers > 8:
        ap.error("at most 8 workers (the workstation has 15 GB)")
    if not 0 <= args.shard < args.n_shards:
        ap.error("--shard must be in [0, --n-shards)")

    def rel(p):
        p = Path(p)
        return p if p.is_absolute() else REPO / p

    data_dir, models_dir, out = rel(args.data_dir), rel(args.models_dir), rel(args.out)
    if args.smoke and args.out == "results/rebuttal/E5":
        out = REPO / "results/rebuttal/E5_smoke"
    args.networks = list(SMOKE_NETWORKS) if args.smoke else [n for n in args.networks.split(",") if n]
    args.n_grid = [int(v) for v in args.n_grid.split(",")]
    args.methods = [m for m in args.methods.split(",") if m]
    if args.smoke:
        args.rand_reps = 1
    out.mkdir(parents=True, exist_ok=True)

    # dots replaced so FIDELITY.yaml can address each hash by a dotted key
    inputs = {p.name.replace('.', '_'): file_sha(p) for p in sorted(data_dir.iterdir()) if p.is_file()}
    cfg = {
        "experiment": "E5_learned_cpdags_by_sample_size", "criterion": CRITERION,
        "repo_commit": repo_commit(), "vendor_sha256": file_sha(HERE / "vendor_bkrobust.py"),
        "run_sha256": file_sha(Path(__file__)), "smoke": bool(args.smoke), "phase": args.phase,
        "networks": args.networks, "n_grid": args.n_grid, "methods": args.methods,
        "rand_reps": args.rand_reps, "drand_seed": DRAND_SEED, "radius_limit": args.radius_limit,
        "ladder_max_nodes": args.ladder_max_nodes, "wall_s": args.wall_s,
        "crosscheck_max_undirected": args.crosscheck_max_undirected,
        "default_search_budget": DEFAULT_SEARCH_BUDGET,
        "exhaustive_max_retractable": EXHAUSTIVE_MAX_RETRACTABLE, "rungs": list(RUNGS),
        "data_dir": str(data_dir), "models_dir": str(models_dir), "input_sha256": inputs,
        "shard": args.shard, "n_shards": args.n_shards, "workers": args.workers,
    }
    tag = f"shard{args.shard}of{args.n_shards}"
    (out / f"config.{tag}.json").write_text(json.dumps(cfg, indent=1))

    inp = load_inputs(data_dir, models_dir, tuple(args.networks))
    blocks = make_blocks(inp, args, smoke_selection(inp) if args.smoke else None)
    done, _ = done_units(out)
    todo = []
    for b in blocks:
        left = [p for p in b["pairs"] if unit_key(b, *p) not in done]
        if left:
            todo.append({**b, "pairs": left})
    n_rows = sum(len(b["pairs"]) for b in todo)
    print(f"E5 {tag}: {len(blocks)} blocks, {len(todo)} to run ({n_rows} rows), "
          f"{args.workers} workers, criterion {CRITERION}, commit {cfg['repo_commit'][:12]}", flush=True)

    t0 = last = time.perf_counter()
    written = 0
    stop = threading.Event()

    def heartbeat():
        # a block of twenty rows can take minutes, so the minute line cannot wait for a block
        while not stop.wait(60):
            print(f"  [alive] {written}/{n_rows} rows, {time.perf_counter() - t0:.0f} s", flush=True)

    threading.Thread(target=heartbeat, daemon=True).start()
    fh = (out / f"rows.{tag}.jsonl").open("a")

    def write(rows):
        nonlocal written
        for r in rows:
            fh.write(json.dumps(r, default=str) + "\n")
        fh.flush()
        written += len(rows)

    if args.workers <= 1:
        init_worker(cfg)
        for i, b in enumerate(todo, 1):
            write(run_block(b))
            if time.perf_counter() - last >= 30 or i == len(todo):
                last = time.perf_counter()
                print(f"  {i}/{len(todo)} blocks, {written}/{n_rows} rows, {last - t0:.0f} s", flush=True)
    else:
        with mp.get_context("spawn").Pool(args.workers, initializer=init_worker, initargs=(cfg,)) as pool:
            it = pool.imap_unordered(run_block, todo)
            i = 0
            while i < len(todo):
                try:
                    rows = it.next(timeout=60)
                except mp.TimeoutError:
                    print(f"  ... {i}/{len(todo)} blocks, {written}/{n_rows} rows, "
                          f"{time.perf_counter() - t0:.0f} s (waiting)", flush=True)
                    continue
                i += 1
                write(rows)
                if time.perf_counter() - last >= 30 or i == len(todo):
                    last = time.perf_counter()
                    print(f"  {i}/{len(todo)} blocks, {written}/{n_rows} rows, {last - t0:.0f} s", flush=True)
    fh.close()
    stop.set()

    _, fails = done_units(out)
    print(f"done: {written} rows written in {time.perf_counter() - t0:.0f} s; "
          f"{len(fails)} rows with a consistency failure", flush=True)
    if fails:
        (out / f"consistency_failures.{tag}.json").write_text(json.dumps(fails, indent=1))
        for u, f in fails[:20]:
            print(f"  FAIL {u}: {f}", flush=True)
        return 2
    if args.n_shards == 1:
        import analyse  # the cross-phase check needs both phases, so only on an unsharded run
        bad = analyse.cross_phase_failures(analyse.load_rows(out))
        if bad:
            print(f"cross-phase check failed on {len(bad)} rows; see analyse.py", flush=True)
            return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
