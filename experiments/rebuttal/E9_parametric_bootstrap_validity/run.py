"""E9: does a parametric bootstrap under a fitted model flag sets invalid at the truth?

For each network of the learned-CPDAG panel, fit a linear-Gaussian SEM to one DAG
of the analyst's state on the reference n = 5000 sample, simulate S datasets from
it, rerun PC exactly as the panel did, impose each arm's knowledge, recompute the
committed set Z_b, and check Z_b with the GAC in the simulated truth. One JSON line
per unit (network, simulation); unit -1 runs PC on the reference sample itself and
must give back the panel's reference CPDAG.

    python experiments/rebuttal/E9_parametric_bootstrap_validity/run.py --smoke
    python experiments/rebuttal/E9_parametric_bootstrap_validity/run.py --workers 6

Exit codes: 0 done, 2 bad arguments or a config clash with an existing output,
3 the rebuilt reference rows disagree with the frozen panel, 4 a model or replicate
invariant broke (no DAG extension, Markov-equivalent fits that differ, Z invalid in
a simulated truth, identity unit not the reference).
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

# PC is single-threaded here; BLAS threads would oversubscribe the workers.
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
from bkrobust.demo.meek import (  # noqa: E402
    apply_orientations,
    is_chordal_components,
    is_consistent_extension,
    meek_closure,
)
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
    sample_ancestral,
    sample_for,
    sha,
    stable_seed,
)

CRITERION = "gac"
NETWORKS = ("ecoli70", "arth150", "magic-niab", "magic-irri")
SMOKE_NETWORKS = ("ecoli70", "magic-niab")
PANEL_N = 5000
SIM_SEED = 20261109
EXT_SEED = 20261110
# The panel's D_RAND draws, so the chance arm's knowledge is regenerated rather than stored.
RAND_SEED = 20260913
RAND_REPS = 3
IDENTITY = -1
COV_RTOL = 1e-6


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


# --- the simulated truth: one DAG of the state, and its least-squares SEM ---------------


def dor_tarsi(g: MPDAG, rng: np.random.Generator | None = None) -> MPDAG | None:
    """A consistent DAG extension of a PDAG (Dor and Tarsi 1992), or None if none exists.

    Repeatedly removes a node with no outgoing directed edge whose undirected
    neighbours are adjacent to all its other neighbours, orienting its undirected
    edges into it. With ``rng`` None the first such node in sorted order is taken,
    so the extension is a deterministic function of the graph; with an ``rng`` the
    choice is random, which gives other members of the class (not uniformly).
    """
    left = set(g.nodes)
    par = {v: set(g.parents(v)) for v in g.nodes}
    ch = {v: set(g.children(v)) for v in g.nodes}
    und = {v: set(g.neighbors(v)) for v in g.nodes}
    arcs = set(g.directed_edges)
    while left:
        cands = []
        for x in sorted(left):
            if ch[x]:
                continue
            adj_x = und[x] | par[x]
            if all((adj_x - {y}) <= (und[y] | par[y] | ch[y]) for y in und[x]):
                cands.append(x)
        if not cands:
            return None
        x = cands[0] if rng is None else cands[int(rng.integers(len(cands)))]
        for y in und[x]:
            arcs.add((y, x))
            und[y].discard(x)
        for p in par[x]:
            ch[p].discard(x)
        left.discard(x)
        und[x], par[x] = set(), set()
    return MPDAG(g.nodes, directed=arcs, undirected=[])


def fit_sem(dag: MPDAG, data: np.ndarray) -> tuple[dict, dict]:
    """Least-squares coefficients and ML residual variances (RSS / n), no intercept.

    The panel's networks have zero means and the simulator draws zero-mean data, so
    an intercept would only add a parameter the simulation ignores.
    """
    idx = {v: i for i, v in enumerate(dag.nodes)}
    weights: dict = {}
    noise: dict = {}
    n = data.shape[0]
    for v in dag.nodes:
        ps = sorted(dag.parents(v))
        y = data[:, idx[v]]
        if ps:
            X = data[:, [idx[p] for p in ps]]
            beta, *_ = np.linalg.lstsq(X, y, rcond=None)
            resid = y - X @ beta
            for p, b in zip(ps, beta):
                weights[(p, v)] = float(b)
        else:
            resid = y
        noise[v] = float(resid @ resid / n)
    return weights, noise


def implied_cov(dag: MPDAG, weights: dict, noise: dict) -> np.ndarray:
    """Covariance of the SEM: X = X B + E gives Sigma = A^T Omega A with A = (I - B)^-1."""
    idx = {v: i for i, v in enumerate(dag.nodes)}
    p = len(idx)
    B = np.zeros((p, p))
    for (a, b), w in weights.items():
        B[idx[a], idx[b]] = w
    A = np.linalg.inv(np.eye(p) - B)
    return A.T @ np.diag([noise[v] for v in dag.nodes]) @ A


def incident(g: MPDAG, nodes: set) -> list[str]:
    out = [f"{a}->{b}" for a, b in g.directed_edges if a in nodes or b in nodes]
    out += [f"{a}-{b}" for a, b in g.undirected_edges if a in nodes or b in nodes]
    return sorted(out)


# --- the reference: every certifiable row rebuilt and checked against the frozen panel
# Same rebuild as E6 (experiments/rebuttal/E6_bootstrap_stability_check/run.py), so the
# two experiments score the same rows; plus the fitted model and its checks.


def no_extension_reason(g: MPDAG, what: str) -> str:
    """Why a state has no consistent DAG extension, in words, for the exclusion record."""
    parts = []
    if not is_chordal_components(g):
        parts.append("an undirected component is not chordal")
    if not g.is_acyclic():
        parts.append("the directed part has a cycle")
    return f"{what} has no consistent DAG extension ({'; '.join(parts) or 'Dor-Tarsi finds no removable node'}), so there is no simulated truth"


def build_reference(pest_dir: Path, models: Path, networks: tuple[str, ...], n_ext: int) -> dict:
    """Rebuild every certifiable row as E6 does, fit the simulated truths, and sort the rows.

    Every frozen row is rebuilt and checked against the panel (K hash, committed set,
    GAC at the truth) whether or not it can be simulated. A row whose state has no
    consistent DAG extension cannot have a simulated truth: it is kept in ``excluded``
    with the reason instead of in a simulation group (DEVIATIONS.md, 2026-11-06).
    """
    chat = json.loads((pest_dir / "chat_e.json").read_text())
    know = json.loads((pest_dir / "knowledge.json").read_text())
    frozen = [
        r
        for r in csv.DictReader((pest_dir / "rows.csv").open())
        if r["status"] == "ok" and r["network"] in networks
    ]
    problems: list[str] = []
    model_problems: list[str] = []
    out: dict = {}
    excluded: list[dict] = []
    excluded_networks: dict[str, str] = {}
    seen = set()
    for net in networks:
        dag, weights, noise = load_gaussian(models, net)
        entry = chat["networks"][net]
        cpdag = graph_from_json(entry)
        if sha(cpdag) != entry["sha256"]:
            problems.append(f"{net}: reference CPDAG hash {sha(cpdag)} != {entry['sha256']}")
            continue
        if cpdag.nodes != dag.nodes:
            problems.append(f"{net}: reference CPDAG and true DAG have different node sets")
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
            ("A_TRUE", 0, [(a, b) if dag.is_directed_edge(a, b) else (b, a) for (a, b) in base_llm if dag.has_edge(a, b)])
        )
        plans.append(
            (
                "A_ORACLE",
                0,
                [(a, b) if dag.is_directed_edge(a, b) else (b, a) for (a, b) in sorted(cpdag.undirected_edges) if dag.has_edge(a, b)],
            )
        )

        # The simulated world: one DAG of the D_DEGEN state (the CPDAG itself), fitted.
        # Without it the network cannot be simulated at all; its rows are still rebuilt
        # and checked below, then excluded with the reason.
        d0_state = apply_orientations(cpdag, []) or cpdag
        d0 = dor_tarsi(d0_state)
        net_reason = None
        if d0 is None or not is_consistent_extension(d0, cpdag):
            net_reason = no_extension_reason(d0_state, f"{net}: the reference CPDAG (D_DEGEN state)")
            excluded_networks[net] = net_reason
            d0 = None
        else:
            w0, n0 = fit_sem(d0, data)
            cov0 = implied_cov(d0, w0, n0)
            scale = float(np.abs(cov0).max())

        groups = []
        for arm, rep, k in plans:
            k_sha = hashlib.sha256(json.dumps(sorted(k)).encode()).hexdigest()[:16]
            mine = [r for r in frozen if r["network"] == net and r["arm"] == arm and int(r["rep"]) == rep]
            if not mine:
                continue
            g0 = apply_orientations(cpdag, k)
            truths: list[MPDAG] = []
            cov_rel = None
            reason = net_reason
            if g0 is not None and reason is None:
                dg = dor_tarsi(g0)
                if dg is None or not is_consistent_extension(dg, g0):
                    reason = no_extension_reason(g0, f"{net}|{arm}|{rep}: the state")
                else:
                    truths.append(dg)
                    ext_rng = np.random.default_rng([EXT_SEED, stable_seed(f"{net}|{arm}|{rep}")])
                    for _ in range(n_ext):
                        de = dor_tarsi(g0, ext_rng)
                        if de is None or not is_consistent_extension(de, g0):
                            model_problems.append(f"{net}|{arm}|{rep}: a random extension is not consistent")
                            break
                        truths.append(de)
                    for t in truths:
                        wt, nt = fit_sem(t, data)
                        rel = float(np.abs(implied_cov(t, wt, nt) - cov0).max() / scale)
                        cov_rel = rel if cov_rel is None else max(cov_rel, rel)
                    if cov_rel is not None and cov_rel > COV_RTOL:
                        # simulating once per network is only right if every truth has the same law
                        model_problems.append(f"{net}|{arm}|{rep}: fitted covariance differs from D_0 by {cov_rel:.2e} (relative)")
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
                    problems.append(f"{rid}: GAC at truth {int(truth)} != frozen {r['z_valid_at_truth']}")
                    continue
                row = {
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
                }
                if reason is not None:
                    excluded.append({**row, "k_sha256": k_sha, "truth_sha256": None, "excluded": reason})
                    continue
                bad = [i for i, t in enumerate(truths) if not is_gac_valid_dag(t, r["x"], r["y"], cs.z)]
                if bad:
                    model_problems.append(f"{rid}: Z is GAC-valid on the state but not in its extension(s) {bad}")
                rows.append(row)
            if reason is not None:
                continue
            groups.append(
                {
                    "arm": arm,
                    "rep": rep,
                    "k": [list(e) for e in k],
                    "k_sha256": k_sha,
                    "truths": [sorted(map(list, t.directed_edges)) for t in truths],
                    "truth_sha256": sha(truths[0]) if truths else None,
                    "cov_rel_vs_d0": cov_rel,
                    "rows": rows,
                }
            )
        if d0 is None or not groups:
            continue
        out[net] = {
            "ref_seed": entry["seed"],
            "cpdag": entry,
            "d0": sorted(map(list, d0.directed_edges)),
            "d0_sha256": sha(d0),
            "weights": [[a, b, w] for (a, b), w in sorted(w0.items())],
            "noise": n0,
            "groups": groups,
        }
    missing = {f"{r['network']}|{r['arm']}|{int(r['rep'])}|{r['x']}|{r['y']}" for r in frozen} - seen
    problems += [f"{rid}: frozen certifiable row has no rebuilt plan" for rid in sorted(missing)]
    return {
        "networks": out,
        "problems": problems,
        "model_problems": model_problems,
        "excluded": excluded,
        "excluded_networks": excluded_networks,
        "n_frozen": len(frozen),
    }


# --- one unit: one simulated dataset of one network ---------------------------------------

_CACHE: dict = {}


def _world(task: dict) -> dict:
    """Per-worker cache of the fitted model and the simulated truths of one network."""
    net = task["network"]
    if net not in _CACHE:
        nodes = task["cpdag"]["nodes"]
        d0 = MPDAG(nodes, directed=[tuple(e) for e in task["d0"]], undirected=[])
        weights = {(a, b): float(w) for a, b, w in task["weights"]}
        truths = {
            (g["arm"], g["rep"]): [MPDAG(nodes, directed=[tuple(e) for e in t], undirected=[]) for t in g["truths"]]
            for g in task["groups"]
        }
        _CACHE[net] = {"d0": d0, "weights": weights, "noise": task["noise"], "truths": truths, "ref": None}
    return _CACHE[net]


def unit_task(task: dict) -> dict:
    net, b = task["network"], task["b"]
    w = _world(task)
    if b == IDENTITY:
        dag, weights, noise = load_gaussian(Path(task["models"]), net)
        data_b = sample_for(dag, weights, noise, net, PANEL_N, task["ref_seed"])
    else:
        rng = np.random.default_rng([SIM_SEED, stable_seed(net), b])
        data_b = sample_ancestral(w["d0"], w["weights"], w["noise"], PANEL_N, rng)
    rec: dict = {
        "unit": f"{net}|{b}",
        "network": net,
        "b": b,
        "criterion": CRITERION,
        "repo_commit": task["commit"],
        "causal_learn_version": task["cl_version"],
        "alpha": ALPHA,
        "n": int(data_b.shape[0]),
        "data_sha256": data_sha(data_b),
        "simulated_from": "reference_sample" if b == IDENTITY else task["d0_sha256"],
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
    g, counts = from_endpoint_matrix(matrix, w["d0"].nodes)
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
        skeleton_missing_vs_model=len(w["d0"].skeleton() - used.skeleton()),
        skeleton_extra_vs_model=len(used.skeleton() - w["d0"].skeleton()),
        skeleton_changed_vs_reference=len(used.skeleton() ^ ref.skeleton()),
    )
    rows = []
    for grp in task["groups"]:
        k = [tuple(e) for e in grp["k"]]
        truths = w["truths"][(grp["arm"], grp["rep"])]
        try:
            state = apply_orientations(used, k)
            state_err = ""
        except Exception as e:  # a non-closed raw graph can break the closure code
            state, state_err = None, type(e).__name__
        for r in grp["rows"]:
            x, y, z = r["x"], r["y"], frozenset(r["z"])
            o = {"rid": r["rid"], "state_ok": int(state is not None), "err": state_err}
            zb = None
            try:
                cs = committed_adjustment_set(state, x, y) if state is not None else None
                zb = cs.z if cs is not None else None
                o["zb_verdict"] = cs.verdict if cs is not None else "no_state"
            except Exception as e:  # counted as no set, never dropped
                o["zb_verdict"], o["err"] = "error", o["err"] or type(e).__name__
            o["zb_exists"] = int(zb is not None)
            o["zb_same"] = int(zb is not None and zb == z)
            try:
                o["zb_valid"] = int(zb is not None and bool(truths) and is_gac_valid_dag(truths[0], x, y, zb))
                o["zb_valid_ext"] = sum(int(zb is not None and is_gac_valid_dag(t, x, y, zb)) for t in truths[1:])
            except Exception as e:
                o["zb_valid"], o["zb_valid_ext"], o["err"] = 0, 0, o["err"] or type(e).__name__
            o["n_ext"] = max(len(truths) - 1, 0)
            try:
                o["zfix_valid_state"] = int(state is not None and is_gac_valid_mpdag(state, x, y, z))
            except Exception as e:
                o["zfix_valid_state"], o["err"] = 0, o["err"] or type(e).__name__
            rows.append(o)
    rec["rows"] = rows
    # the closure caches are keyed by graph and never shrink; each replicate's graphs are new
    clear_cache()
    return rec


# --- driver -----------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--smoke", action="store_true", help="two small networks, two simulations, two extensions")
    ap.add_argument("--out", default=None, help="output folder (default results/rebuttal/E9)")
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--S", type=int, default=200, help="simulated datasets per network")
    ap.add_argument("--n-ext", type=int, default=5, help="random extra DAG extensions per state (secondary)")
    ap.add_argument("--networks", default=",".join(NETWORKS))
    ap.add_argument("--wall", type=int, default=WALL_S, help="seconds per PC run before a timeout")
    ap.add_argument(
        "--check-reference",
        action="store_true",
        help="only rebuild the reference rows and fit the models, write reference.jsonl and exclusions.json, run no PC",
    )
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
    S, n_ext = args.S, args.n_ext
    if args.smoke:
        networks, S, n_ext = SMOKE_NETWORKS, 2, 2
    workers = args.workers or (1 if args.smoke else 6)
    out = Path(args.out) if args.out else ROOT / "results" / "rebuttal" / ("E9_smoke" if args.smoke else "E9")
    if not (0 <= args.shard < args.n_shards) or S < 1 or n_ext < 0 or workers < 1 or set(networks) - set(NETWORKS):
        print("bad --shard/--n-shards/--S/--n-ext/--workers/--networks", file=sys.stderr)
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
        "experiment": "E9",
        "criterion": CRITERION,
        "method": "pc",
        "indep_test": "fisherz",
        "alpha": ALPHA,
        "panel_n": PANEL_N,
        "S": S,
        "n_ext": n_ext,
        "sim_seed": SIM_SEED,
        "ext_seed": EXT_SEED,
        "networks": list(networks),
        "wall_s": args.wall,
        "smoke": bool(args.smoke),
        "resample": "parametric_linear_gaussian_ols_ml_variance",
        "simulated_truth": "dor_tarsi_sorted_extension_of_state",
        "cov_rtol": COV_RTOL,
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
    print(f"rebuilding the reference rows and fitting the models on {len(networks)} networks", flush=True)
    ref = build_reference(pest_dir, models, networks, n_ext)
    # both lists are written before either stops the run, so one cannot hide the other
    if ref["model_problems"]:
        (out / "model_problems.txt").write_text("\n".join(ref["model_problems"]) + "\n")
    if ref["problems"]:
        (out / "reference_problems.txt").write_text("\n".join(ref["problems"]) + "\n")
        print(f"CONSISTENCY FAILURE: {len(ref['problems'])} rebuilt rows disagree with the frozen panel", file=sys.stderr)
        for p in ref["problems"][:20]:
            print("  " + p, file=sys.stderr)
        if ref["model_problems"]:
            print(f"  (and {len(ref['model_problems'])} fitted-model problems, in model_problems.txt)", file=sys.stderr)
        return 3
    if ref["model_problems"]:
        print(f"CONSISTENCY FAILURE: {len(ref['model_problems'])} fitted-model invariants broke", file=sys.stderr)
        for p in ref["model_problems"][:20]:
            print("  " + p, file=sys.stderr)
        return 4
    excluded = ref["excluded"]
    reasons: dict[str, int] = {}
    for r in excluded:
        reasons[r["excluded"]] = reasons.get(r["excluded"], 0) + 1
    (out / "exclusions.json").write_text(
        json.dumps(
            {
                "n_frozen": ref["n_frozen"],
                "n_excluded": len(excluded),
                "excluded_networks": ref["excluded_networks"],
                "reasons": reasons,
                "rows": [r["rid"] for r in excluded],
            },
            indent=1,
        )
    )
    n_rows = sum(len(g["rows"]) for e in ref["networks"].values() for g in e["groups"])
    tmp = out / "reference.jsonl.tmp"
    with tmp.open("w") as fh:
        for net, e in ref["networks"].items():
            for g in e["groups"]:
                for r in g["rows"]:
                    fh.write(
                        json.dumps(
                            {
                                **r,
                                "k_sha256": g["k_sha256"],
                                "truth_sha256": g["truth_sha256"],
                                "excluded": None,
                                "criterion": CRITERION,
                                "repo_commit": commit,
                            }
                        )
                        + "\n"
                    )
        # every frozen row has a line; an excluded row carries its reason and is never simulated
        for r in excluded:
            fh.write(json.dumps({**r, "criterion": CRITERION, "repo_commit": commit}) + "\n")
    tmp.replace(out / "reference.jsonl")
    models_out = {
        net: {
            "d0_sha256": e["d0_sha256"],
            "n_arcs": len(e["d0"]),
            "weights_sha256": hashlib.sha256(json.dumps(e["weights"]).encode()).hexdigest()[:16],
            "groups": {f"{g['arm']}|{g['rep']}": {"truth_sha256": g["truth_sha256"], "n_truths": len(g["truths"]), "cov_rel_vs_d0": g["cov_rel_vs_d0"]} for g in e["groups"]},
        }
        for net, e in ref["networks"].items()
    }
    (out / "models.json").write_text(json.dumps(models_out, indent=1))
    worst = max(((g["cov_rel_vs_d0"] or 0.0) for e in ref["networks"].values() for g in e["groups"]), default=0.0)
    print(
        f"reference rebuilt: {n_rows + len(excluded)} of {ref['n_frozen']} certifiable rows, all match the frozen panel; "
        f"{n_rows} to simulate, {len(excluded)} excluded; worst fitted-covariance gap {worst:.1e} ({time.perf_counter() - t0:.0f} s)",
        flush=True,
    )
    for why, c in sorted(reasons.items()):
        print(f"  excluded {c} rows: {why}", flush=True)
    if args.check_reference:
        return 0

    # identity units first, so a pipeline that does not reproduce the panel stops early
    # a network with no simulated truth (all its rows excluded above) has no units
    units = [(net, b) for b in [IDENTITY, *range(S)] for net in networks if net in ref["networks"]]
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
    tasks = []
    for net, b in todo:
        e = ref["networks"][net]
        tasks.append(
            {
                "network": net,
                "b": b,
                "models": str(models),
                "ref_seed": e["ref_seed"],
                "cpdag": e["cpdag"],
                "d0": e["d0"],
                "d0_sha256": e["d0_sha256"],
                "weights": e["weights"],
                "noise": e["noise"],
                "groups": e["groups"],
                "commit": commit,
                "cl_version": cl_version,
                "wall": args.wall,
            }
        )

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
                if not rec["is_reference_cpdag"]:
                    failures.append(f"{rec['unit']}: PC on the reference sample learned {rec['cpdag_sha256']}, reference is {ref['networks'][rec['network']]['cpdag']['sha256']}")
                else:
                    # on the reference state the rule must recommit Z, and Z must hold in its own extension
                    bad = [r["rid"] for r in rec["rows"] if not (r["zb_same"] and r["zb_valid"] and r["zb_valid_ext"] == r["n_ext"])]
                    if bad:
                        failures.append(f"{rec['unit']}: Z not recommitted or not valid in its extensions for {len(bad)} rows, e.g. {bad[0]}")
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
