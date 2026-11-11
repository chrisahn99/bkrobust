"""E25: are the nearest failing states of the committed set exactly the non-amenable states?

For each query (see PREREG.md for the populations) the knowledge states above G0 are enumerated
level by level (level j = the states reached by j upward covers from G0, i.e. j retracted
orientations), and on every state up to the level where both events have happened we evaluate

* fails: the committed set Z (valid at G0) is not GAC-valid there;
* not amenable: (x, y) is not amenable there (submitted ``is_amenable``).

r_val is the first level holding a failing state, r_amen the first level holding a non-amenable
state. Registered: (a) every failing state at level r_val is non-amenable; (b) r_val = r_amen.
Beside them: the submitted hybrid radius and an upward search for the amenability radius, each
timed, a one-call amenability test at G0, timed, the literal b-amenability of Perković et al.
(2017, Def. 3.1 and 4.3) at every nearest failing state, and the validity of Pa(x, G) at every
amenable state visited (Perković 2020, Prop. 4.2).

One JSON line per query, flushed as it goes, skipped on resume. Exit 1 on a consistency failure,
2 if some units raised.

    .venv/bin/python experiments/rebuttal/E25_amenability_check/run.py --smoke
    .venv/bin/python experiments/rebuttal/E25_amenability_check/run.py --workers 6
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import multiprocessing as mp
import subprocess
import sys
import threading
import time
import traceback
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
for p in (ROOT / "src", HERE):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from bkrobust.benchmarks.describe import parse_file, to_mpdag  # noqa: E402
from bkrobust.core import oracle as core_oracle  # noqa: E402
from bkrobust.core.conventions import UNREACHED  # noqa: E402
from bkrobust.demo.example import dag_to_cpdag  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.gac.mpdag_level import clear_cache, is_gac_valid_mpdag, why_invalid_gac  # noqa: E402
from bkrobust.hybrid import breakdown_radius  # noqa: E402
from bkrobust.mpdag_criterion.criterion import is_amenable  # noqa: E402
from bkrobust.search.exact import retractable_edges  # noqa: E402
from bkrobust.search.exact_fast import local_up_covers_fast  # noqa: E402

import vendor_states as vs  # noqa: E402
from b_amenability import b_amenable  # noqa: E402
from vendor_ledger import RELEASE_SEED, committed_adjustment_set, k_sha, rebuild_knowledge, sha  # noqa: E402
from vendor_search import bfs_radius  # noqa: E402

CRITERION = "gac"
DATA = ROOT / "experiments" / "rebuttal" / "_data"
ELICITED_ARMS = (
    "D_LLM", "D_LLM_INSTR", "D_LLM_SRC", "D_SCRAMBLED", "D_LLM_72B", "D_LLM_72B_INSTR",
    "D_SCRAMBLED_72B", "D_LLM_32B", "D_LLM_SRC_70B", "D_LLM_GEMMA_27B",
)
RELEASED_ARMS = (*ELICITED_ARMS, "D_RAND")
LADDER_MAX_NODES = 60          # as E10: above it the hybrid's SAT ladder does not return
EXHAUSTIVE_MAX_RETRACTABLE = 12
DEFAULT_DEPTH = 3
SOURCES = ("released", "e11", "synthetic")


def file_sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def repo_commit() -> tuple[str, bool | None]:
    try:
        head = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "--", "src"],
                               capture_output=True, text=True, check=True).stdout.strip()
        return head, bool(dirty)
    except (OSError, subprocess.CalledProcessError):
        return "unknown", None


def clear_all() -> None:
    clear_cache()
    core_oracle.clear_cache()


# ------------------------------------------------------------------ the measurements


def kset(cpdag, g) -> frozenset:
    return frozenset(retractable_edges(cpdag, g))


def level_scan(cpdag, g0, x, y, z, cfg) -> dict:
    """Enumerate the up-set of G0 level by level until both r_val and r_amen are found.

    A level is complete when every state on it was generated and evaluated within the budget.
    r_val (r_amen) is decided when its first level is found with every earlier level complete,
    or when the up-set is exhausted. (a) needs the whole r_val level.
    """
    t0 = time.perf_counter()
    # DEVIATIONS.md (1): when Z is valid and (x, y) amenable at the top state (the CPDAG), both
    # radii are infinite by the submitted top-state rule; the scan then only checks that rule on
    # up-sets small enough to exhaust (top_scan_max_states, top_scan_time_limit).
    top_ok = is_gac_valid_mpdag(cpdag, x, y, z) and is_amenable(cpdag, x, y)
    max_states = cfg["top_scan_max_states"] if top_ok else cfg["scan_max_states"]
    deadline = t0 + (cfg["top_scan_time_limit"] if top_ok else cfg["scan_time_limit"])
    k0 = kset(cpdag, g0)
    level = [g0]
    seen = {g0.edge_string()}
    n_states = 0
    hard: list[str] = []
    out = {"r_val": None, "r_amen": None, "r_val_decided": False, "r_amen_decided": False,
           "fail_level_complete": False, "amen_level_complete": False, "stop": None}
    j = 0
    prop42 = {"amenable_states": 0, "pa_invalid": 0, "y_in_pa": 0, "examples": []}
    while level:
        complete = True
        fail_states, nonamen_states = [], []
        for g in level:
            if n_states >= max_states or time.perf_counter() > deadline:
                complete = False
                out["stop"] = "states" if n_states >= max_states else "time"
                break
            n_states += 1
            kg = kset(cpdag, g)
            if not kg <= k0 or len(k0 - kg) != j:
                hard.append(f"level {j}: state with |K0 minus K| = {len(k0 - kg)}, subset {kg <= k0}")
            fo = not is_gac_valid_mpdag(g, x, y, z)
            na = not is_amenable(g, x, y)
            if fo:
                fail_states.append((g, na))
            if na:
                nonamen_states.append((g, fo))
            if not na and out["r_val"] is None:
                # Perković (2020) Prop. 4.2: amenable and y not in Pa(x, G) => Pa(x, G) is valid.
                pa = frozenset(g.parents(x))
                prop42["amenable_states"] += 1
                if y in pa:
                    prop42["y_in_pa"] += 1
                elif not is_gac_valid_mpdag(g, x, y, pa):
                    prop42["pa_invalid"] += 1
                    if len(prop42["examples"]) < 3:
                        prop42["examples"].append(g.edge_string())
        if j == 0 and (fail_states or nonamen_states):
            hard.append("committed set fails, or (x, y) is not amenable, at G0")
        if out["r_val"] is None and fail_states:
            out["r_val"], out["r_val_decided"], out["fail_level_complete"] = j, True, complete
            out["fail_level"] = summarise_fail_level(g0, fail_states, x, y, z, cfg, complete)
        if out["r_amen"] is None and nonamen_states:
            out["r_amen"], out["r_amen_decided"], out["amen_level_complete"] = j, True, complete
            out["amen_level"] = {"n_nonamenable": len(nonamen_states),
                                 "n_nonamenable_and_failing": sum(fo for _, fo in nonamen_states),
                                 "complete": complete}
        if not complete:
            break
        if out["r_val"] is not None and out["r_amen"] is not None:
            out["stop"] = "both_found"
            break
        # next level
        nxt = []
        for g in level:
            if time.perf_counter() > deadline:
                complete = False
                out["stop"] = "time"
                break
            for h in local_up_covers_fast(cpdag, g):
                key = h.edge_string()
                if key not in seen:
                    seen.add(key)
                    nxt.append(h)
        if not complete:
            break
        level = nxt
        j += 1
        if len(seen) > 4 * max_states:
            out["stop"] = "states"
            break
    if not level and out["stop"] is None:
        out["stop"] = "exhausted"
        for k in ("r_val", "r_amen"):
            if out[k] is None:
                out[k], out[f"{k}_decided"] = UNREACHED, True
    if top_ok and (out["r_val"] not in (None, UNREACHED) or out["r_amen"] not in (None, UNREACHED)):
        hard.append("top state valid and amenable, yet a failing or non-amenable state was found")
    if top_ok and out["stop"] != "exhausted":
        out["decided_by_top_state"] = True
        for k in ("r_val", "r_amen"):
            if out[k] is None:
                out[k], out[f"{k}_decided"] = UNREACHED, True
    out["top_state_ok"] = top_ok
    out.update({"levels_scanned": j, "n_states": n_states, "prop42": prop42,
                "seconds": round(time.perf_counter() - t0, 4), "hard_failures": hard})
    return out


def summarise_fail_level(g0, fail_states, x, y, z, cfg, complete) -> dict:
    reasons = Counter()
    n_lit_agree = n_lit_undecided = 0
    exceptions = []
    pa_has_y = 0
    for g, na in fail_states:
        reasons[why_invalid_gac(g, x, y, z) or "valid?"] += 1
        lit = b_amenable(g, x, y, max_steps=cfg["bpc_max_steps"])
        if lit is None:
            n_lit_undecided += 1
        elif lit == (not na):
            n_lit_agree += 1
        if y in g.parents(x):
            pa_has_y += 1
        if not na and len(exceptions) < 5:
            exceptions.append({"state": g.edge_string(), "reason": why_invalid_gac(g, x, y, z),
                               "b_amenable_literal": lit})
    return {"n_failing": len(fail_states), "n_failing_nonamenable": sum(na for _, na in fail_states),
            "n_failing_amenable": sum(not na for _, na in fail_states), "complete": complete,
            "reasons": dict(reasons), "literal_agree": n_lit_agree, "literal_undecided": n_lit_undecided,
            "y_parent_of_x": pa_has_y, "exceptions": exceptions,
            "g0_if_exception": g0.edge_string() if exceptions else None}


def hybrid_radius(cpdag, g0, x, y, z, cfg) -> dict:
    """The submitted radius for Z, as E10's set_radius (ladder only up to 60 nodes)."""
    n_ret = len(retractable_edges(cpdag, g0))
    depth = n_ret if n_ret <= EXHAUSTIVE_MAX_RETRACTABLE else DEFAULT_DEPTH
    t0 = time.perf_counter()
    if len(cpdag.nodes) <= LADDER_MAX_NODES:
        hb = breakdown_radius(cpdag, None, x, y, z, g0=g0, search_budget=depth, criterion=CRITERION,
                              time_limit_s=cfg["ladder_time_limit"])
        rec = {"r": hb.radius if hb.exact else None, "exact": hb.exact, "method": hb.method}
    else:
        if is_gac_valid_mpdag(cpdag, x, y, z):
            rec = {"r": UNREACHED, "exact": True, "method": "top_state"}
        else:
            b = bfs_radius(cpdag, g0, lambda g: not is_gac_valid_mpdag(g, x, y, z), max_depth=None,
                           max_closures=cfg["closure_budget"], time_limit_s=cfg["search_time_limit"])
            rec = {"r": b["radius"] if b["exact"] else None, "exact": b["exact"], "method": "bfs_budget",
                   "lb": b["lb"]}
    rec["seconds"] = round(time.perf_counter() - t0, 5)
    return rec


def amen_radius(cpdag, g0, x, y, cfg) -> dict:
    """Amenability radius as the first author's r_amen: top check, then the same upward search."""
    t0 = time.perf_counter()
    if is_amenable(cpdag, x, y):
        rec = {"r": UNREACHED, "exact": True, "method": "top_state"}
    else:
        b = bfs_radius(cpdag, g0, lambda g: not is_amenable(g, x, y), max_depth=None,
                       max_closures=cfg["closure_budget"], time_limit_s=cfg["search_time_limit"])
        rec = {"r": b["radius"] if b["exact"] else None, "exact": b["exact"], "method": "bfs_budget",
               "lb": b["lb"]}
    rec["seconds"] = round(time.perf_counter() - t0, 5)
    return rec


def measure(cpdag, g0, x, y, z, cfg) -> dict:
    t0 = time.perf_counter()
    a0 = is_amenable(g0, x, y)
    t_call = time.perf_counter() - t0
    hb = hybrid_radius(cpdag, g0, x, y, z, cfg)
    clear_all()
    am = amen_radius(cpdag, g0, x, y, cfg)
    clear_all()
    sc = level_scan(cpdag, g0, x, y, z, cfg)
    hard = list(sc.pop("hard_failures"))
    if not a0:
        hard.append("not amenable at G0 although Z is valid there")
    # two radius computations must agree where both are exact
    if hb["exact"] and sc["r_val_decided"] and hb["r"] != sc["r_val"]:
        hard.append(f"hybrid radius {hb['r']} vs level scan {sc['r_val']}")
    if am["exact"] and sc["r_amen_decided"] and am["r"] != sc["r_amen"]:
        hard.append(f"amenability search {am['r']} vs level scan {sc['r_amen']}")
    rec = {"amenable_g0": a0, "t_amen_call": round(t_call, 6), "hybrid": hb, "amen_search": am,
           "scan": sc, "n_nodes": len(cpdag.nodes), "n_retractable": len(retractable_edges(cpdag, g0)),
           "g0_undirected": len(g0.undirected_edges), "hard_failures": hard}
    exc = sc.get("fail_level", {}).get("n_failing_amenable", 0) > 0
    r_v = sc["r_val"] if sc["r_val_decided"] else (hb["r"] if hb["exact"] else None)
    r_a = sc["r_amen"] if sc["r_amen_decided"] else (am["r"] if am["exact"] else None)
    if exc or hard or (r_v is not None and r_a is not None and r_v != r_a):
        # the graph of every exception is kept, so it can be listed
        rec["exception_graphs"] = {"cpdag": cpdag.edge_string(), "g0": g0.edge_string(),
                                   "x": x, "y": y, "z": sorted(z)}
    return rec


# ------------------------------------------------------------------ the three sources

_NETS: dict[str, tuple] = {}
_KNOW: dict[str, dict] = {}


def released_network(models: Path, name: str):
    if name not in _NETS:
        path = next(p for p in sorted(models.iterdir()) if p.name.startswith(name))
        dag = to_mpdag(parse_file(path, file_sha(path)))
        _NETS.clear()
        _NETS[name] = (dag, dag_to_cpdag(dag))
    return _NETS[name]


def compute_released(u: dict, cfg: dict) -> dict:
    row = u["row"]
    dag, cpdag = released_network(Path(cfg["models"]), row["network"])
    if row["knowledge"] not in _KNOW:
        _KNOW[row["knowledge"]] = json.loads(Path(row["knowledge"]).read_text())
    k = rebuild_knowledge(row["network"], row["arm"], int(row["rep"]), dag, cpdag,
                          _KNOW[row["knowledge"]], cfg["seed"])
    g0 = apply_orientations(cpdag, k) if k is not None else None
    mism = []
    if sha(cpdag) != row["chat_source_sha256"]:
        mism.append("chat_source_sha256")
    if g0 is None or k_sha(k) != row["k_sha256"]:
        mism.append("k_sha256")
    if mism:
        return {"status": "rebuild_mismatch", "hard_failures": mism}
    cs = committed_adjustment_set(g0, row["x"], row["y"])
    if cs.z is None:
        return {"status": "rebuild_mismatch", "hard_failures": ["no committed set"]}
    rec = {"status": "ok", "z": sorted(cs.z), **measure(cpdag, g0, row["x"], row["y"], cs.z, cfg)}
    if row["hop_exact"] == "1" and row["r_hop"] != "" and rec["hybrid"]["exact"] \
            and int(row["r_hop"]) != rec["hybrid"]["r"]:
        rec["hard_failures"].append(f"radius {rec['hybrid']['r']} vs released r_hop {row['r_hop']}")
    rec["released_r_hop"] = row["r_hop"]
    return rec


def compute_e11(u: dict, cfg: dict) -> dict:
    dag, cpdag, g0, status = vs.e11_state(cfg["elicit"], cfg["networks_dir"], u["condition"], u["network"])
    if g0 is None:
        return {"status": status, "hard_failures": []}
    z, reason = vs.e11_commit(g0, u["x"], u["y"])
    if z is None:
        return {"status": reason, "hard_failures": []}
    return {"status": "ok", "z": sorted(z), **measure(cpdag, g0, u["x"], u["y"], z, cfg)}


def compute_synthetic(u: dict, cfg: dict) -> dict:
    row = u["row"]
    inst, why = vs.rebuild(row)
    if inst is None:
        return {"status": "rebuild_failed", "reason": str(why), "hard_failures": ["rebuild_failed"]}
    z_stored = frozenset(s for s in row["z_star"].split(",") if s)
    hard = []
    if inst["x"] != row["x"] or frozenset(inst["z_star"]) != z_stored:
        hard.append("x or z_star differs from the stored instance")
    rec = {"status": "ok", "x": inst["x"], "y": inst["y"], "z": sorted(inst["z_star"]),
           **measure(inst["cpdag"], inst["g0"], inst["x"], inst["y"], frozenset(inst["z_star"]), cfg)}
    rec["hard_failures"] += hard
    rec["stored_r_val"] = row.get("r_val")
    return rec


def run_unit(task: tuple) -> dict:
    u, cfg = task
    t0 = time.perf_counter()
    rec = {k: v for k, v in u.items() if k != "row"}
    rec.update({"criterion": CRITERION, "commit": cfg["commit"], "smoke": cfg["smoke"]})
    try:
        fn = {"released": compute_released, "e11": compute_e11, "synthetic": compute_synthetic}[u["source"]]
        rec.update(fn(u, cfg))
    except Exception:
        rec["status"] = "error"
        rec["error"] = traceback.format_exc()[-2000:]
    clear_all()
    rec["seconds"] = round(time.perf_counter() - t0, 3)
    return rec


# ------------------------------------------------------------------ unit lists


def released_units(records: Path, knowledge: Path, smoke: bool) -> list[dict]:
    rows = [r for r in csv.DictReader(records.open()) if r["status"] == "ok" and r["arm"] in RELEASED_ARMS]
    units = []
    for r in rows:
        r = dict(r, knowledge=str(knowledge))
        units.append({"key": "released|" + "|".join((r["network"], r["arm"], r["rep"], r["x"], r["y"])),
                      "source": "released", "network": r["network"], "arm": r["arm"], "rep": r["rep"],
                      "x": r["x"], "y": r["y"], "row": r})
    units.sort(key=lambda u: u["key"])
    if smoke:
        pick, seen = [], set()
        for u in units:
            tag = (u["network"] in ("asia", "sachs", "magic-irri", "Sebastiani_2005"), u["arm"] == "D_RAND")
            if u["network"] in ("asia", "sachs", "magic-irri", "Sebastiani_2005") and (u["network"], tag) not in seen:
                seen.add((u["network"], tag))
                pick.append(u)
        units = pick[:6]
    return units


def e11_units(elicit: str, frames: dict, smoke: bool) -> list[dict]:
    conds = vs.e11_conditions(elicit)
    seen, units = set(), []
    for fname, path in frames.items():
        for fu in vs.frame_units(fname, vs.read_jsonl(Path(path))):
            for c in conds:
                key = "e11|" + "|".join((c, fu["network"], fu["x"], fu["y"]))
                if key in seen:  # pairs in both frames are scored once, under the first frame
                    continue
                seen.add(key)
                units.append({"key": key, "source": "e11", "frame": fname, "condition": c,
                              "network": fu["network"], "x": fu["x"], "y": fu["y"]})
    units.sort(key=lambda u: u["key"])
    if smoke:
        units = [u for u in units if u["network"] in ("asia", "sachs") and u["condition"] == "D_LLM"][:6]
    return units


def synthetic_units(targeted: Path, untargeted: Path, smoke: bool) -> list[dict]:
    units = []
    for which, path in (("targeted", targeted), ("untargeted", untargeted)):
        rows = list(csv.DictReader(path.open()))
        if smoke:
            got, pick = Counter(), []
            for r in rows:
                if r["component_size"] == "6" and got[r["arm"]] < 2:
                    got[r["arm"]] += 1
                    pick.append(r)
            rows = pick
        for r in rows:
            units.append({"key": f"synthetic|{which}|{r['instance_id']}", "source": "synthetic",
                          "which": which, "instance_id": r["instance_id"], "arm": r["arm"], "row": r})
    return units


def done_keys(out: Path) -> set[str]:
    keys = set()
    for f in sorted(out.glob("rows*.jsonl")):
        for line in f.open():
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            if rec.get("status") != "error":
                keys.add(rec["key"])
    return keys


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--sources", default=",".join(SOURCES))
    ap.add_argument("--records", type=Path, default=DATA / "per_query_records.csv")
    ap.add_argument("--knowledge", type=Path, default=DATA / "knowledge_states.json")
    ap.add_argument("--models", type=Path, default=ROOT / "results/axisa3/networks/example_models")
    ap.add_argument("--elicit", default=str(ROOT / "results/elicit/knowledge.json"))
    ap.add_argument("--networks-dir", default=str(ROOT / "results/axisa3/networks"))
    ap.add_argument("--targeted", type=Path, default=ROOT / vs.TARGETED_CSV)
    ap.add_argument("--untargeted", type=Path, default=ROOT / vs.UNTARGETED_CSV)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--seed", type=int, default=RELEASE_SEED)
    ap.add_argument("--scan-max-states", type=int, default=200_000)
    ap.add_argument("--scan-time-limit", type=float, default=300.0)
    ap.add_argument("--closure-budget", type=int, default=20_000)
    ap.add_argument("--search-time-limit", type=float, default=120.0)
    ap.add_argument("--ladder-time-limit", type=float, default=60.0)
    ap.add_argument("--bpc-max-steps", type=int, default=200_000)
    ap.add_argument("--top-scan-max-states", type=int, default=2_000)
    ap.add_argument("--top-scan-time-limit", type=float, default=20.0)
    args = ap.parse_args(argv)
    if args.workers > 8:
        ap.error("at most 8 workers (the workstation is shared; 12 cores at most)")
    if args.smoke:
        args.scan_max_states, args.scan_time_limit = 5_000, 20.0
        args.closure_budget, args.search_time_limit, args.ladder_time_limit = 2_000, 10.0, 5.0
        args.bpc_max_steps = 20_000
        args.top_scan_max_states, args.top_scan_time_limit = 500, 5.0
    out = args.out or ROOT / ("results/rebuttal/E25_smoke" if args.smoke else "results/rebuttal/E25")
    out.mkdir(parents=True, exist_ok=True)
    commit, src_dirty = repo_commit()
    sources = [s for s in args.sources.split(",") if s]
    frames = {n: str(ROOT / p) for n, (p, _) in vs.FRAMES.items()}
    cfg = {
        "experiment": "E25", "criterion": CRITERION, "commit": commit, "src_dirty": src_dirty,
        "smoke": args.smoke, "sources": sources, "seed": args.seed, "models": str(args.models),
        "elicit": args.elicit, "networks_dir": args.networks_dir,
        "scan_max_states": args.scan_max_states, "scan_time_limit": args.scan_time_limit,
        "closure_budget": args.closure_budget, "search_time_limit": args.search_time_limit,
        "ladder_time_limit": args.ladder_time_limit, "bpc_max_steps": args.bpc_max_steps,
        "top_scan_max_states": args.top_scan_max_states, "top_scan_time_limit": args.top_scan_time_limit,
        "ladder_max_nodes": LADDER_MAX_NODES, "workers": args.workers,
    }
    units: list[dict] = []
    if "released" in sources:
        cfg["records_sha256"] = file_sha(args.records)
        cfg["knowledge_sha256"] = file_sha(args.knowledge)
        ru = released_units(args.records, args.knowledge, args.smoke)
        cfg["n_released"] = len(ru)
        units += ru
    if "e11" in sources:
        cfg["frame_sha256"] = {n: file_sha(Path(p)) for n, p in frames.items()}
        cfg["frame_sha_ok"] = all(cfg["frame_sha256"][n] == vs.FRAMES[n][1] for n in frames)
        cfg["elicit_sha256"] = file_sha(Path(args.elicit))
        eu = e11_units(args.elicit, frames, args.smoke)
        cfg["n_e11"] = len(eu)
        units += eu
    if "synthetic" in sources:
        cfg["targeted_sha256"] = file_sha(args.targeted)
        cfg["untargeted_sha256"] = file_sha(args.untargeted)
        su = synthetic_units(args.targeted, args.untargeted, args.smoke)
        cfg["n_synthetic"] = len(su)
        units += su
    cfg["n_units"] = len(units)
    (out / "config.json").write_text(json.dumps(cfg, indent=1, sort_keys=True))

    have = done_keys(out)
    todo = [(u, cfg) for u in units if u["key"] not in have]
    print(f"[E25] {len(units)} units, {len(units) - len(todo)} done, {len(todo)} to run; "
          f"commit {commit[:12]}, src dirty {src_dirty}", flush=True)
    fh = (out / "rows.jsonl").open("a")
    state = {"done": 0, "hard": 0, "err": 0, "t0": time.perf_counter(), "stop": False}

    def heartbeat():
        while not state["stop"]:
            time.sleep(60)
            if not state["stop"]:
                print(f"[E25] {state['done']}/{len(todo)} done, hard {state['hard']}, errors {state['err']}, "
                      f"{time.perf_counter() - state['t0']:.0f} s", flush=True)

    threading.Thread(target=heartbeat, daemon=True).start()

    def consume(rec: dict):
        fh.write(json.dumps(rec, sort_keys=True) + "\n")
        fh.flush()
        state["done"] += 1
        if rec.get("status") == "error":
            state["err"] += 1
            print(f"[E25] ERROR {rec['key']}: {rec['error'].strip().splitlines()[-1]}", flush=True)
        if rec.get("hard_failures"):
            state["hard"] += 1
            print(f"[E25] CONSISTENCY {rec['key']}: {rec['hard_failures'][:3]}", flush=True)
        if args.smoke:
            sc = rec.get("scan", {})
            print(f"[E25] {rec['key']} {rec.get('status')} r_val {sc.get('r_val')} r_amen {sc.get('r_amen')} "
                  f"hyb {rec.get('hybrid', {}).get('r')} {rec.get('seconds')} s", flush=True)

    try:
        if args.workers <= 1:
            for t in todo:
                consume(run_unit(t))
        else:
            with mp.get_context("spawn").Pool(args.workers, maxtasksperchild=200) as pool:
                for rec in pool.imap_unordered(run_unit, todo, chunksize=1):
                    consume(rec)
    finally:
        state["stop"] = True
        fh.close()
    print(f"[E25] finished {state['done']} units in {time.perf_counter() - state['t0']:.0f} s; "
          f"hard {state['hard']}, errors {state['err']}", flush=True)
    if state["hard"]:
        print("[E25] CONSISTENCY FAILURE (see rows with hard_failures)", flush=True)
        return 1
    return 2 if state["err"] else 0


if __name__ == "__main__":
    sys.exit(main())
