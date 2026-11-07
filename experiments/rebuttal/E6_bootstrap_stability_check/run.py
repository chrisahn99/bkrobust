"""E6: does bootstrap stability of the committed set flag sets invalid at the truth?

For each network of the learned-CPDAG panel, resample the reference n = 5000
sample B times, rerun PC exactly as the panel did, impose each arm's knowledge
on the learned graph, and record per certifiable row whether its committed set Z
is still GAC-valid there and whether the graph around the query moved. One JSON
line per unit (network, replicate); unit -1 is the identity replicate, which
must give back the panel's reference CPDAG.

    python experiments/rebuttal/E6_bootstrap_stability_check/run.py --smoke
    python experiments/rebuttal/E6_bootstrap_stability_check/run.py --workers 6

Exit codes: 0 done, 2 bad arguments or a config clash with an existing output,
3 the rebuilt reference rows disagree with the frozen panel rows, 4 a replicate
broke an invariant (identity replicate not the reference, or Z invalid on it).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import multiprocessing as mp
import os
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path

# PC is single-threaded by design here; BLAS threads would oversubscribe the workers.
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
for _p in (ROOT / "src", HERE):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import numpy as np  # noqa: E402

from bkrobust.demo.graph import MPDAG  # noqa: E402
from bkrobust.demo.meek import apply_orientations, meek_closure  # noqa: E402
from bkrobust.gac.dag_level import is_gac_valid_dag  # noqa: E402
from bkrobust.gac.mpdag_level import clear_cache, is_gac_valid_mpdag  # noqa: E402
from vendor_pest import (  # noqa: E402
    ALPHA,
    WALL_S,
    chance_level,
    committed_adjustment_set,
    data_sha,
    from_endpoint_matrix,
    graph_from_json,
    learn_graph,
    load_gaussian,
    sample_for,
    sha,
    stable_seed,
)

CRITERION = "gac"
NETWORKS = ("ecoli70", "arth150", "magic-niab", "magic-irri")
SMOKE_NETWORKS = ("ecoli70", "magic-niab")
PANEL_N = 5000
BOOT_SEED = 20261106
# The panel's D_RAND draws, so the chance arm's knowledge is regenerated rather than stored.
RAND_SEED = 20260913
RAND_REPS = 3
IDENTITY = -1


class _Wall(Exception):
    pass


def _fire(signum, frame):
    raise _Wall


def repo_commit() -> str:
    try:
        out = subprocess.run(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True, timeout=10
        )
        return out.stdout.strip() or "unknown"
    except (OSError, subprocess.TimeoutExpired):
        return "unknown"


def file_sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def incident(g: MPDAG, nodes: set) -> list[str]:
    """Every edge touching ``nodes``, with its mark; the query's neighbourhood in a state."""
    out = [f"{a}->{b}" for a, b in g.directed_edges if a in nodes or b in nodes]
    out += [f"{a}-{b}" for a, b in g.undirected_edges if a in nodes or b in nodes]
    return sorted(out)


# --- the reference: every certifiable row rebuilt and checked against the frozen panel


def build_reference(pest_dir: Path, models: Path, networks: tuple[str, ...]) -> dict:
    """Rebuild K, G0 and Z for every certifiable row; stop on any disagreement.

    The frozen rows carry Z's size but not Z, so Z is recomputed with the
    panel's own functions, and the rebuild must reproduce the K hash, the size,
    the verdict and the validity at truth (decided here by the submitted
    ``gac.dag_level``, there by an independent networkx audit).
    """
    chat = json.loads((pest_dir / "chat_e.json").read_text())
    know = json.loads((pest_dir / "knowledge.json").read_text())
    frozen = [
        r
        for r in csv.DictReader((pest_dir / "rows.csv").open())
        if r["status"] == "ok" and r["network"] in networks
    ]
    problems: list[str] = []
    out: dict = {}
    seen = set()
    for net in networks:
        dag, weights, noise = load_gaussian(models, net)
        entry = chat["networks"][net]
        cpdag = graph_from_json(entry)
        if sha(cpdag) != entry["sha256"]:
            problems.append(f"{net}: reference CPDAG hash {sha(cpdag)} != {entry['sha256']}")
            continue
        data = sample_for(dag, weights, noise, net, PANEL_N, entry["seed"])
        if data_sha(data) != entry["sample_sha256"]:
            problems.append(f"{net}: reference sample did not regenerate (numpy version?)")
            continue
        base_llm = [tuple(e) for e in know["D_LLM"]["networks"].get(net, {}).get("k", [])]
        plans: list[tuple[str, int, list]] = [("D_LLM", 0, base_llm)]
        if base_llm:
            for rep in range(RAND_REPS):
                rng = np.random.default_rng([RAND_SEED, stable_seed(net), rep])
                kk = chance_level(cpdag, base_llm, rng)
                if kk is not None:
                    plans.append(("D_RAND", rep, kk))
        plans.append(("D_DEGEN", 0, []))
        plans.append(
            (
                "A_TRUE",
                0,
                [(a, b) if dag.is_directed_edge(a, b) else (b, a) for (a, b) in base_llm if dag.has_edge(a, b)],
            )
        )
        plans.append(
            (
                "A_ORACLE",
                0,
                [
                    (a, b) if dag.is_directed_edge(a, b) else (b, a)
                    for (a, b) in sorted(cpdag.undirected_edges)
                    if dag.has_edge(a, b)
                ],
            )
        )
        groups = []
        for arm, rep, k in plans:
            k_sha = hashlib.sha256(json.dumps(sorted(k)).encode()).hexdigest()[:16]
            mine = [r for r in frozen if r["network"] == net and r["arm"] == arm and int(r["rep"]) == rep]
            if not mine:
                continue
            g0 = apply_orientations(cpdag, k)
            rows = []
            for r in mine:
                rid = f"{net}|{arm}|{rep}|{r['x']}|{r['y']}"
                seen.add(rid)
                if r["k_sha256"] != k_sha:
                    problems.append(f"{rid}: K hash {k_sha} != frozen {r['k_sha256']}")
                    continue
                if g0 is None:
                    problems.append(f"{rid}: knowledge inconsistent on rebuild, frozen row says ok")
                    continue
                cs = committed_adjustment_set(g0, r["x"], r["y"])
                if cs.z is None or len(cs.z) != int(r["z_size"]) or cs.verdict != r["committed_verdict"]:
                    problems.append(f"{rid}: committed set rebuilt as {cs.verdict}/{cs.z}, frozen {r['committed_verdict']}/{r['z_size']}")
                    continue
                truth = is_gac_valid_dag(dag, r["x"], r["y"], cs.z)
                if int(truth) != int(r["z_valid_at_truth"]):
                    problems.append(f"{rid}: GAC at truth {int(truth)} (gac.dag_level) != frozen {r['z_valid_at_truth']} (networkx audit)")
                    continue
                rows.append(
                    {
                        "rid": rid,
                        "network": net,
                        "arm": arm,
                        "rep": rep,
                        "x": r["x"],
                        "y": r["y"],
                        "z": sorted(cs.z),
                        "verdict": cs.verdict,
                        "z_valid_at_truth": int(truth),
                        "r_cert_claim": r["r_cert_claim"],
                        "adj_x": sorted(cpdag.adjacent(r["x"])),
                        "adj_y": sorted(cpdag.adjacent(r["y"])),
                        "nb_state": incident(g0, {r["x"], r["y"]}),
                    }
                )
            groups.append({"arm": arm, "rep": rep, "k": [list(e) for e in k], "k_sha256": k_sha, "rows": rows})
        out[net] = {"ref_seed": entry["seed"], "cpdag": entry, "groups": groups}
    missing = {f"{r['network']}|{r['arm']}|{int(r['rep'])}|{r['x']}|{r['y']}" for r in frozen} - seen
    problems += [f"{rid}: frozen certifiable row has no rebuilt plan" for rid in sorted(missing)]
    return {"networks": out, "problems": problems, "n_frozen": len(frozen)}


# --- one unit: one replicate of one network -----------------------------------------

_CACHE: dict = {}


def unit_task(task: dict) -> dict:
    """Resample, learn, impose every arm's K, and score every row of the network."""
    net, b = task["network"], task["b"]
    if net not in _CACHE:
        dag, weights, noise = load_gaussian(Path(task["models"]), net)
        _CACHE[net] = (dag, sample_for(dag, weights, noise, net, PANEL_N, task["ref_seed"]))
    dag, data = _CACHE[net]
    if b == IDENTITY:
        data_b = data
    else:
        rng = np.random.default_rng([BOOT_SEED, stable_seed(net), b])
        data_b = data[rng.integers(0, data.shape[0], size=data.shape[0])]
    rec: dict = {
        "unit": f"{net}|{b}",
        "network": net,
        "b": b,
        "criterion": CRITERION,
        "repo_commit": task["commit"],
        "causal_learn_version": task["cl_version"],
        "alpha": ALPHA,
        "n": int(data_b.shape[0]),
        "resample_sha256": data_sha(data_b),
        "wall_s": task["wall"],
    }
    signal.signal(signal.SIGALRM, _fire)
    signal.setitimer(signal.ITIMER_REAL, task["wall"])
    t0 = time.perf_counter()
    try:
        matrix = learn_graph(data_b, "pc", ALPHA)
    except _Wall:
        rec.update(status="timeout", seconds=round(time.perf_counter() - t0, 1), rows=[])
        return rec
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    rec["seconds"] = round(time.perf_counter() - t0, 1)
    g, counts = from_endpoint_matrix(matrix, dag.nodes)
    closure = meek_closure(g)
    # the panel used the raw graph where the closure failed; kept, and flagged
    used = closure if closure is not None else g
    ref = graph_from_json(task["cpdag"])
    rec.update(
        status="ok",
        closure_ok=int(closure is not None),
        n_conflict_edges=counts["conflict"],
        cpdag_sha256=sha(used),
        is_reference_cpdag=int(sha(used) == task["cpdag"]["sha256"]),
        skeleton_missing_vs_truth=len(dag.skeleton() - used.skeleton()),
        skeleton_extra_vs_truth=len(used.skeleton() - dag.skeleton()),
        skeleton_changed_vs_reference=len(used.skeleton() ^ ref.skeleton()),
    )
    rows = []
    for grp in task["groups"]:
        k = [tuple(e) for e in grp["k"]]
        try:
            state = apply_orientations(used, k)
            state_err = ""
        except Exception as e:  # a non-closed raw graph can break the closure code
            state, state_err = None, type(e).__name__
        for r in grp["rows"]:
            x, y, z = r["x"], r["y"], frozenset(r["z"])
            out = {"rid": r["rid"], "state_ok": int(state is not None), "err": state_err}
            try:
                out["z_valid"] = int(state is not None and is_gac_valid_mpdag(state, x, y, z))
            except Exception as e:  # counted as not valid, never dropped
                out["z_valid"], out["err"] = 0, type(e).__name__
            try:
                cs = committed_adjustment_set(state, x, y) if state is not None else None
                out["z_same"] = int(cs is not None and cs.z is not None and cs.z == z)
            except Exception as e:
                out["z_same"], out["err"] = 0, out["err"] or type(e).__name__
            out["skel_nb_stable"] = int(
                sorted(used.adjacent(x)) == r["adj_x"] and sorted(used.adjacent(y)) == r["adj_y"]
            )
            out["state_nb_stable"] = int(state is not None and incident(state, {x, y}) == r["nb_state"])
            rows.append(out)
    rec["rows"] = rows
    # the closure caches are keyed by graph and never shrink; each replicate's graphs are new
    clear_cache()
    return rec


# --- driver -----------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--smoke", action="store_true", help="two small networks, two replicates")
    ap.add_argument("--out", default=None, help="output folder (default results/rebuttal/E6)")
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--B", type=int, default=100, help="bootstrap replicates per network")
    ap.add_argument("--networks", default=",".join(NETWORKS))
    ap.add_argument("--wall", type=int, default=WALL_S, help="seconds per PC run before a timeout")
    ap.add_argument(
        "--pest-dir",
        default=str(ROOT / "experiments" / "rebuttal" / "_data" / "pest"),
        help="folder with chat_e.json, knowledge.json, rows.csv of the learned-CPDAG panel",
    )
    ap.add_argument(
        "--models",
        default=str(ROOT / "results" / "axisa3" / "networks" / "example_models"),
        help="folder with the four Gaussian network JSON files",
    )
    args = ap.parse_args(argv)

    networks = tuple(s for s in args.networks.split(",") if s)
    B = args.B
    if args.smoke:
        networks, B = SMOKE_NETWORKS, 2
    workers = args.workers or (1 if args.smoke else 6)
    out = Path(args.out) if args.out else ROOT / "results" / "rebuttal" / ("E6_smoke" if args.smoke else "E6")
    if not (0 <= args.shard < args.n_shards) or B < 1 or workers < 1 or set(networks) - set(NETWORKS):
        print("bad --shard/--n-shards/--B/--workers/--networks", file=sys.stderr)
        return 2
    pest_dir, models = Path(args.pest_dir), Path(args.models)
    out.mkdir(parents=True, exist_ok=True)

    try:
        cl_version = importlib.metadata.version("causal-learn")
    except importlib.metadata.PackageNotFoundError:
        print("causal-learn is not installed (see requirements.txt)", file=sys.stderr)
        return 2
    commit = repo_commit()
    config = {
        "experiment": "E6",
        "criterion": CRITERION,
        "method": "pc",
        "indep_test": "fisherz",
        "alpha": ALPHA,
        "panel_n": PANEL_N,
        "B": B,
        "boot_seed": BOOT_SEED,
        "networks": list(networks),
        "wall_s": args.wall,
        "smoke": bool(args.smoke),
        "resample": "rows_with_replacement",
        "identity_replicate": True,
        "causal_learn_version": cl_version,
        "numpy_version": np.__version__,
        "inputs_sha256": {f: file_sha(pest_dir / f) for f in ("chat_e.json", "knowledge.json", "rows.csv")},
        "vendor_sha256": file_sha(HERE / "vendor_pest.py"),
        # a later edit of the preregistration then shows up as a config clash
        "prereg_sha256": file_sha(HERE / "PREREG.md"),
    }
    cfg_path = out / "config.json"
    if cfg_path.exists():
        old = json.loads(cfg_path.read_text())
        if {k: v for k, v in old.items() if k != "repo_commit"} != config:
            print(f"{cfg_path} holds a different configuration; use another --out", file=sys.stderr)
            return 2
    cfg_path.write_text(json.dumps({**config, "repo_commit": commit}, indent=1))

    t0 = time.perf_counter()
    ref = build_reference(pest_dir, models, networks)
    if ref["problems"]:
        (out / "reference_problems.txt").write_text("\n".join(ref["problems"]) + "\n")
        print(f"CONSISTENCY FAILURE: {len(ref['problems'])} rebuilt rows disagree with the frozen panel", file=sys.stderr)
        for p in ref["problems"][:20]:
            print("  " + p, file=sys.stderr)
        return 3
    n_rows = sum(len(g["rows"]) for e in ref["networks"].values() for g in e["groups"])
    tmp = out / "reference.jsonl.tmp"
    with tmp.open("w") as fh:
        for net, e in ref["networks"].items():
            for g in e["groups"]:
                for r in g["rows"]:
                    fh.write(json.dumps({**r, "k_sha256": g["k_sha256"], "criterion": CRITERION, "repo_commit": commit}) + "\n")
    tmp.replace(out / "reference.jsonl")
    print(f"reference rebuilt: {n_rows} certifiable rows on {len(networks)} networks, all match the frozen panel ({time.perf_counter() - t0:.0f} s)", flush=True)

    # identity replicates first, so a pipeline that does not reproduce the panel stops early
    units = [(net, b) for b in [IDENTITY, *range(B)] for net in networks]
    units = [u for i, u in enumerate(units) if i % args.n_shards == args.shard]
    units_path = out / f"units_{args.shard:03d}_of_{args.n_shards:03d}.jsonl"
    done = set()
    if units_path.exists():
        for line in units_path.read_text().splitlines():
            try:
                done.add(json.loads(line)["unit"])
            except (json.JSONDecodeError, KeyError):
                pass  # a line cut by a kill; the unit is simply redone
    todo = [u for u in units if f"{u[0]}|{u[1]}" not in done]
    print(f"shard {args.shard}/{args.n_shards}: {len(units)} units, {len(done)} already written, {len(todo)} to run on {workers} workers", flush=True)
    tasks = [
        {
            "network": net,
            "b": b,
            "models": str(models),
            "ref_seed": ref["networks"][net]["ref_seed"],
            "cpdag": ref["networks"][net]["cpdag"],
            "groups": ref["networks"][net]["groups"],
            "commit": commit,
            "cl_version": cl_version,
            "wall": args.wall,
        }
        for net, b in todo
    ]

    state = {"done": 0, "last": "-"}
    stop = threading.Event()

    def heartbeat():
        while not stop.wait(30):
            print(f"  .. {state['done']}/{len(tasks)} units, last {state['last']}, {time.perf_counter() - t0:.0f} s", flush=True)

    threading.Thread(target=heartbeat, daemon=True).start()
    failures: list[str] = []
    with units_path.open("a") as fh, mp.get_context("spawn").Pool(workers) as pool:
        for rec in pool.imap_unordered(unit_task, tasks):
            # the line carries every row, so the unit is written whole or not at all
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
            state["done"] += 1
            state["last"] = f"{rec['unit']} {rec['status']} {rec.get('seconds')} s"
            print(f"  [{state['done']}/{len(tasks)}] {state['last']}", flush=True)
            if rec["b"] == IDENTITY and rec["status"] == "ok":
                bad = [r["rid"] for r in rec["rows"] if not r["z_valid"]]
                if not rec["is_reference_cpdag"]:
                    failures.append(f"{rec['unit']}: identity replicate learned {rec['cpdag_sha256']}, reference is {ref['networks'][rec['network']]['cpdag']['sha256']}")
                if bad:
                    failures.append(f"{rec['unit']}: Z not valid on the reference state for {len(bad)} rows, e.g. {bad[0]}")
                if failures:
                    pool.terminate()
                    break
    stop.set()
    if failures:
        print("CONSISTENCY FAILURE:", file=sys.stderr)
        for f in failures:
            print("  " + f, file=sys.stderr)
        return 4
    print(f"shard {args.shard}: finished in {time.perf_counter() - t0:.0f} s -> {units_path}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
