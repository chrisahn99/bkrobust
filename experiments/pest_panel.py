"""P-EST: the estimated-CPDAG panel, Table 7.

Every ledger row so far was computed on the CPDAG of the true DAG and stamped
``via_Chat_only``. This panel replaces that CPDAG with one estimated from data on
the four networks that carry real Gaussian parameters, so the frame, the
questionnaire, the knowledge and every instrument see only what an analyst
would see: samples, a discovery algorithm's output and a supplier's answers.
The true DAG draws the samples and audits the frozen rows, and is on no other
path, so the deployment rows carry ``true_dag_on_path = "false"`` and their
``chat_source_sha256`` differs from ``generating_dag_sha256``.

Stages, each writing under ``results/pest/``:

``structure``     ancestral samples at n in {200, 1000, 5000}, seeds 0..9; PC with
                  Fisher-z at alpha 0.01 and GES with BIC; every learned graph
                  converted to the repository's MPDAG, checked Meek-closed and
                  scored on undirected fraction, chain-component partition,
                  skeleton precision and recall against the truth, and the
                  adjusted Rand index of partitions across seeds. One reference
                  run per network, PC at n = 5000 seed 0, is C-hat_E.
``frame``         the candidate frame from C-hat_E under ``frame_build``'s rules.
``questionnaire`` ``questionnaire_build``'s renderer on C-hat_E's components.
``elicit``        the primary supplier through ``elicit_run``, cached; the prompts
                  are exported instead when the host is unreachable.
``sweep``         ``ledger_sweep``'s instruments on every row and arm, frozen and
                  hashed, then the audit, with two statuses the oracle panel
                  cannot have and the realised-bias calibration A-CAL-2.

The pre-registration in ``results/pest/PREREGISTRATION.md`` fixes the settings
and the compute bound; ``pest_table.py`` prints Table 7.

    python experiments/pest_panel.py structure --workers 6
    python experiments/pest_panel.py frame
    python experiments/pest_panel.py questionnaire
    python experiments/pest_panel.py elicit
    python experiments/pest_panel.py sweep
"""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import importlib.metadata
import itertools
import json
import multiprocessing as mp
import os
import signal
import subprocess
import sys
import time
from pathlib import Path
from types import FrameType

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT / "src", ROOT / "experiments"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import numpy as np  # noqa: E402
from sklearn.metrics import adjusted_rand_score  # noqa: E402

import elicit_export  # noqa: E402
import elicit_run  # noqa: E402
import questionnaire_build  # noqa: E402
from bkrobust.benchmarks.describe import parse_file, to_mpdag  # noqa: E402
from bkrobust.benchmarks.discovery import (  # noqa: E402
    METHODS,
    from_endpoint_matrix,
    gaussian_parameters,
    learn_graph,
    sample_ancestral,
    structure_scores,
)
from bkrobust.core.conventions import UNREACHED  # noqa: E402
from bkrobust.demo.evaluate import LinearSEM  # noqa: E402
from bkrobust.demo.example import dag_to_cpdag  # noqa: E402
from bkrobust.demo.graph import MPDAG, undirected_components  # noqa: E402
from bkrobust.demo.meek import apply_orientations, is_consistent_extension  # noqa: E402
from bkrobust.mpdag_criterion.paths import possible_descendants  # noqa: E402
from frame_build import MIN_COMPONENT, treatments  # noqa: E402
from ledger_sweep import RUNGS, TIER, audit_row, chance_level, instrument_row, sha  # noqa: E402
from questionnaire_build import stable_seed  # noqa: E402

MODELS = ROOT / "results" / "axisa3" / "networks" / "example_models"
OUT = ROOT / "results" / "pest"

NETWORKS = ("ecoli70", "arth150", "magic-niab", "magic-irri")
N_GRID = (200, 1000, 5000)
SEEDS = tuple(range(10))
ALPHA = 0.01
#: Wall per discovery run, and the one cell run on fewer seeds; both are in the
#: pre-registration and a run that hits the wall is a counted ``timeout`` row.
WALL_S = 900
GES_SEED_CAP = {"arth150": 3}
REFERENCE = {"method": "pc", "n": 5000, "seed": 0}
PANEL_SEED = 20260912
Z975 = 1.959963984540054


# --- the substrate -----------------------------------------------------------


def load_gaussian(name: str) -> tuple[MPDAG, dict, dict, dict[str, str]]:
    """The true DAG, its fitted coefficients and residual variances, and the raw names."""
    path = MODELS / f"{name}.json"
    parsed = parse_file(path, hashlib.sha256(path.read_bytes()).hexdigest())
    dag = to_mpdag(parsed)
    weights, noise = gaussian_parameters(path)
    if set(weights) != set(dag.directed_edges):
        raise ValueError(f"{name}: coefficients and arcs disagree")
    return dag, weights, noise, dict(parsed.name_map)


def sample_for(dag: MPDAG, weights: dict, noise: dict, name: str, n: int, seed: int) -> np.ndarray:
    """The panel's sample for one cell, reproducible from its four seeds."""
    rng = np.random.default_rng([PANEL_SEED, stable_seed(name), n, seed])
    return sample_ancestral(dag, weights, noise, n, rng)


def data_sha(data: np.ndarray) -> str:
    """A short hash of a sample, so a later stage can prove it regenerated the same one."""
    return hashlib.sha256(np.ascontiguousarray(data).tobytes()).hexdigest()[:16]


def graph_from_json(e: dict) -> MPDAG:
    """An MPDAG from the ``nodes``/``directed``/``undirected`` fields of a record."""
    return MPDAG(
        e["nodes"],
        directed=[tuple(x) for x in e["directed"]],
        undirected=[tuple(x) for x in e["undirected"]],
    )


def graph_to_json(g: MPDAG) -> dict:
    """The inverse of :func:`graph_from_json`."""
    return {
        "nodes": list(g.nodes),
        "directed": sorted(g.directed_edges),
        "undirected": sorted(g.undirected_edges),
    }


# --- stage: structure ----------------------------------------------------------


class _WallError(Exception):
    """Raised inside a discovery run that exceeded its wall."""


def _fire(signum: int, frame: FrameType | None) -> None:
    raise _WallError


def structure_task(task: tuple[str, int, int, str]) -> dict:
    """One discovery run: sample, learn, convert, score. Runs in a worker process."""
    name, n, seed, method = task
    dag, weights, noise, _ = load_gaussian(name)
    data = sample_for(dag, weights, noise, name, n, seed)
    rec: dict = {
        "network": name,
        "n": n,
        "seed": seed,
        "method": method,
        "wall_s": WALL_S,
        "n_nodes": len(dag.nodes),
        "n_true_edges": len(dag.directed_edges),
        "sample_sha256": data_sha(data),
    }
    signal.signal(signal.SIGALRM, _fire)
    signal.setitimer(signal.ITIMER_REAL, WALL_S)
    t0 = time.perf_counter()
    try:
        matrix = learn_graph(data, method, ALPHA)
    except _WallError:
        rec.update(status="timeout", seconds=round(time.perf_counter() - t0, 1))
        return rec
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    rec["seconds"] = round(time.perf_counter() - t0, 1)
    g, counts = from_endpoint_matrix(matrix, dag.nodes)
    scores, used = structure_scores(g, dag, dag_to_cpdag(dag))
    rec.update(scores)
    rec["status"] = "ok"
    rec["n_conflict_edges"] = counts["conflict"]
    rec["n_other_endpoints"] = counts["other"]
    rec["graph"] = graph_to_json(used)
    return rec


STRUCTURE_FIELDS = (
    "network",
    "n",
    "seed",
    "method",
    "status",
    "seconds",
    "wall_s",
    "n_nodes",
    "n_true_edges",
    "n_learned_edges",
    "n_directed",
    "n_undirected",
    "undirected_fraction",
    "n_conflict_edges",
    "n_other_endpoints",
    "meek_closed_on_arrival",
    "closure_ok",
    "n_forced_by_closure",
    "acyclic",
    "skeleton_tp",
    "skeleton_precision",
    "skeleton_recall",
    "n_directed_in_truth",
    "oriented_correct",
    "n_components_ge2",
    "max_component",
    "ari_vs_oracle",
    "is_reference",
    "sample_sha256",
    "components",
)


def stage_structure(workers: int) -> None:
    """Run every discovery cell, score it, and designate the reference CPDAGs."""
    tasks = [
        (name, n, seed, method)
        for name in NETWORKS
        for method in METHODS
        for n in N_GRID
        for seed in SEEDS
        if method != "ges" or seed < GES_SEED_CAP.get(name, len(SEEDS))
    ]
    size = {name: len(load_gaussian(name)[0].nodes) for name in NETWORKS}
    tasks.sort(key=lambda t: (-size[t[0]], -t[1], t[3]))  # largest first, for balance
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"{len(tasks)} discovery runs on {workers} workers, wall {WALL_S} s each", flush=True)
    t0 = time.perf_counter()
    records: list[dict] = []
    with mp.get_context("spawn").Pool(workers) as pool:
        for i, rec in enumerate(pool.imap_unordered(structure_task, tasks), 1):
            records.append(rec)
            print(
                f"  [{i}/{len(tasks)}] {rec['network']} {rec['method']} n={rec['n']} "
                f"seed={rec['seed']} {rec['status']} {rec['seconds']} s "
                f"({time.perf_counter() - t0:.0f} s elapsed)",
                flush=True,
            )
    records.sort(key=lambda r: (r["network"], r["method"], r["n"], r["seed"]))

    # the oracle CPDAG per network, for the comparison every prediction needs
    oracle: dict[str, dict] = {}
    dags: dict[str, MPDAG] = {}
    for name in NETWORKS:
        dag = load_gaussian(name)[0]
        cp = dag_to_cpdag(dag)
        dags[name] = dag
        comps = undirected_components(cp)
        oracle[name] = {
            "n_nodes": len(dag.nodes),
            "n_edges": len(dag.directed_edges),
            "n_undirected": len(cp.undirected_edges),
            "undirected_fraction": round(len(cp.undirected_edges) / len(dag.directed_edges), 6),
            "n_components_ge2": len(comps),
            "max_component": max((len(c) for c in comps), default=0),
            "sha256": sha(cp),
        }

    # the reference run: PC at n = 5000, the first seed whose closure succeeded
    chat: dict = {
        "meta": {
            **REFERENCE,
            "indep_test": "fisherz",
            "alpha": ALPHA,
            "causal_learn_version": importlib.metadata.version("causal-learn"),
            "panel_seed": PANEL_SEED,
            "true_dag_on_path": "false",
        },
        "networks": {},
    }
    for name in NETWORKS:
        cands = sorted(
            (
                r
                for r in records
                if r["network"] == name
                and r["method"] == REFERENCE["method"]
                and r["n"] == REFERENCE["n"]
                and r["status"] == "ok"
                and r["closure_ok"]
            ),
            key=lambda r: r["seed"],
        )
        if not cands:
            raise SystemExit(f"{name}: no usable reference run")
        ref = cands[0]
        ref["is_reference"] = 1
        g = graph_from_json(ref["graph"])
        chat["networks"][name] = {
            **ref["graph"],
            "seed": ref["seed"],
            "reference_fallback": int(ref["seed"] != REFERENCE["seed"]),
            "sha256": sha(g),
            "generating_dag_sha256": sha(dags[name]),
            "sample_sha256": ref["sample_sha256"],
            **{
                k: ref[k]
                for k in (
                    "meek_closed_on_arrival",
                    "closure_ok",
                    "n_conflict_edges",
                    "undirected_fraction",
                    "skeleton_precision",
                    "skeleton_recall",
                    "oriented_correct",
                    "n_components_ge2",
                    "max_component",
                    "ari_vs_oracle",
                )
            },
        }
    (OUT / "chat_e.json").write_text(json.dumps(chat, indent=1))

    with (OUT / "graphs.jsonl").open("w") as fh:
        for r in records:
            if r["status"] == "ok":
                fh.write(
                    json.dumps({k: r[k] for k in ("network", "n", "seed", "method")} | r["graph"])
                    + "\n"
                )
    with (OUT / "structure.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(STRUCTURE_FIELDS))
        w.writeheader()
        for r in records:
            row = {k: r.get(k) for k in STRUCTURE_FIELDS}
            row["is_reference"] = r.get("is_reference", 0)
            row["components"] = json.dumps(r.get("components", []))
            w.writerow(row)

    # per (network, method, n): the seed-level summary and the ARI across seeds
    groups: dict[tuple, list[dict]] = collections.defaultdict(list)
    for r in records:
        groups[(r["network"], r["method"], r["n"])].append(r)
    cells: dict[str, dict] = {}
    for (name, method, n), rs in sorted(groups.items()):
        ok = [r for r in rs if r["status"] == "ok"]
        pairs = [
            float(adjusted_rand_score(a["partition"], b["partition"]))
            for a, b in itertools.combinations(ok, 2)
        ]

        def mean(key: str, rows: list[dict] = ok) -> float | None:
            vals = [r[key] for r in rows if r.get(key) is not None]
            return round(float(np.mean(vals)), 6) if vals else None

        cells[f"{name}|{method}|{n}"] = {
            "network": name,
            "method": method,
            "n": n,
            "runs": len(rs),
            "runs_ok": len(ok),
            "runs_timeout": sum(1 for r in rs if r["status"] == "timeout"),
            "seconds_mean": mean("seconds", rs),
            "undirected_fraction_mean": mean("undirected_fraction"),
            "undirected_fraction_min": min(
                (r["undirected_fraction"] for r in ok if r["undirected_fraction"] is not None),
                default=None,
            ),
            "undirected_fraction_max": max(
                (r["undirected_fraction"] for r in ok if r["undirected_fraction"] is not None),
                default=None,
            ),
            "skeleton_precision_mean": mean("skeleton_precision"),
            "skeleton_recall_mean": mean("skeleton_recall"),
            "oriented_correct_mean": mean("oriented_correct"),
            "meek_closed_on_arrival_share": mean("meek_closed_on_arrival"),
            "closure_ok_share": mean("closure_ok"),
            "conflict_edges_mean": mean("n_conflict_edges"),
            "components_ge2_mean": mean("n_components_ge2"),
            "ari_across_seeds_mean": round(float(np.mean(pairs)), 6) if pairs else None,
            "ari_across_seeds_min": round(float(min(pairs)), 6) if pairs else None,
            "ari_vs_oracle_mean": mean("ari_vs_oracle"),
        }
    summary = {
        "seconds": round(time.perf_counter() - t0, 1),
        "runs": len(records),
        "timeouts": sum(1 for r in records if r["status"] == "timeout"),
        "settings": {
            "networks": NETWORKS,
            "n_grid": N_GRID,
            "seeds": SEEDS,
            "methods": METHODS,
            "alpha": ALPHA,
            "wall_s": WALL_S,
            "ges_seed_cap": GES_SEED_CAP,
            "reference": REFERENCE,
            "panel_seed": PANEL_SEED,
            "causal_learn_version": importlib.metadata.version("causal-learn"),
        },
        "oracle": oracle,
        "cells": cells,
        "reference": {
            name: {k: v for k, v in e.items() if k not in ("nodes", "directed", "undirected")}
            for name, e in chat["networks"].items()
        },
    }
    (OUT / "structure_summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary["reference"], indent=1))


def stage_structure_retry(workers: int) -> None:
    """Rerun the cells that hit the wall, under the same wall, and record them apart.

    The first pass ran six workers at once and GES on the largest network hit
    the wall under that contention. This reruns exactly those cells with the
    same wall and writes them to ``structure_retry.csv`` and
    ``graphs_retry.jsonl``. The first pass's rows are left as they are, and the
    reference CPDAG is unaffected, since it is a PC run.
    """
    rows = list(csv.DictReader((OUT / "structure.csv").open()))
    tasks = [
        (r["network"], int(r["n"]), int(r["seed"]), r["method"])
        for r in rows
        if r["status"] == "timeout"
    ]
    print(f"{len(tasks)} timed-out cells rerun on {workers} workers, wall {WALL_S} s", flush=True)
    t0 = time.perf_counter()
    records: list[dict] = []
    with mp.get_context("spawn").Pool(workers) as pool:
        for i, rec in enumerate(pool.imap_unordered(structure_task, tasks), 1):
            records.append(rec)
            print(
                f"  [{i}/{len(tasks)}] {rec['network']} {rec['method']} n={rec['n']} "
                f"seed={rec['seed']} {rec['status']} {rec['seconds']} s",
                flush=True,
            )
    records.sort(key=lambda r: (r["network"], r["method"], r["n"], r["seed"]))
    with (OUT / "graphs_retry.jsonl").open("w") as fh:
        for r in records:
            if r["status"] == "ok":
                fh.write(
                    json.dumps({k: r[k] for k in ("network", "n", "seed", "method")} | r["graph"])
                    + "\n"
                )
    with (OUT / "structure_retry.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[*STRUCTURE_FIELDS, "attempt"])
        w.writeheader()
        for r in records:
            row = {k: r.get(k) for k in STRUCTURE_FIELDS}
            row["is_reference"] = 0
            row["components"] = json.dumps(r.get("components", []))
            row["attempt"] = 2
            w.writerow(row)
    print(f"done in {time.perf_counter() - t0:.0f} s", flush=True)


# --- stage: frame ----------------------------------------------------------------


def stage_frame(cap: int) -> None:
    """The candidate frame from C-hat_E, under the same rules as ``frame_build``."""
    chat = json.loads((OUT / "chat_e.json").read_text())
    meta = chat["meta"]
    source = f"estimated_cpdag:{meta['method']}_{meta['indep_test']}_a{meta['alpha']}_n{meta['n']}"
    rows_out = (OUT / "frame.jsonl").open("w")
    per_net: list[dict] = []
    totals: collections.Counter = collections.Counter()
    for name in NETWORKS:
        e = chat["networks"][name]
        cpdag = graph_from_json(e)
        cands, comp_size = treatments(cpdag)
        frame: list[tuple[str, str]] = []
        for x in sorted(cands):
            for y in sorted(possible_descendants(cpdag, x) - {x}):
                frame.append((x, y))
        frame.sort()
        step = max(1, len(frame) // cap) if len(frame) > cap else 1
        sampled = frame[::step][:cap]
        inclusion = len(sampled) / len(frame) if frame else 0.0
        per_net.append(
            {
                "network": name,
                "nodes": len(cpdag.nodes),
                "undirected_edges": len(cpdag.undirected_edges),
                "components_ge2": sum(
                    1 for c in undirected_components(cpdag) if len(c) >= MIN_COMPONENT
                ),
                "candidate_treatments": len(cands),
                "frame_pairs": len(frame),
                "sampled_pairs": len(sampled),
                "inclusion_probability": round(inclusion, 6),
            }
        )
        totals["frame_pairs"] += len(frame)
        totals["sampled_pairs"] += len(sampled)
        for x, y in sampled:
            rows_out.write(
                json.dumps(
                    {
                        "network": name,
                        "x": x,
                        "y": y,
                        "component_size": comp_size.get(x, 0),
                        "inclusion_probability": round(inclusion, 6),
                        "frame_source": f"{source}_seed{e['seed']}",
                        "chat_source_sha256": e["sha256"],
                        "true_dag_on_path": "false",
                    }
                )
                + "\n"
            )
    rows_out.close()
    with (OUT / "per_network.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(per_net[0]))
        w.writeheader()
        w.writerows(per_net)
    summary = {"cap": cap, "min_component": MIN_COMPONENT, "source": source, **totals}
    (OUT / "frame_summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


# --- stage: questionnaire, elicit ------------------------------------------------


def stage_questionnaire() -> None:
    """Render the questionnaire on C-hat_E's components with the shared renderer."""
    questionnaire_build.main(
        [
            "--frame",
            str(OUT / "frame.jsonl"),
            "--cpdag-source",
            str(OUT / "chat_e.json"),
            "--out",
            str(OUT),
        ]
    )


def host_reachable() -> bool:
    """Whether the local supplier's host answers, with no model call made."""
    try:
        proc = subprocess.run(
            [
                "ssh",
                "-o",
                "BatchMode=yes",
                "-o",
                "ConnectTimeout=8",
                elicit_run.HOST,
                "curl",
                "-s",
                "-m",
                "5",
                "http://localhost:11434/api/tags",
            ],
            capture_output=True,
            text=True,
            timeout=40,
        )
    except (subprocess.TimeoutExpired, OSError):
        return False
    return proc.returncode == 0 and '"models"' in proc.stdout


def stage_elicit(export_only: bool) -> bool:
    """Ask the primary supplier, or export the prompts when it cannot be reached."""
    if export_only or not host_reachable():
        n = elicit_export.export(OUT / "questionnaire.json", OUT / "prompts.jsonl")
        print(
            f"supplier host unreachable or export requested: {n} prompts exported to "
            f"{OUT / 'prompts.jsonl'}; the panel stops before elicitation",
            flush=True,
        )
        return False
    elicit_run.main(
        [
            "--conditions",
            "D_LLM",
            "--questionnaire",
            str(OUT / "questionnaire.json"),
            "--cpdag-source",
            str(OUT / "chat_e.json"),
            "--out",
            str(OUT),
        ]
    )
    return True


# --- stage: sweep ----------------------------------------------------------------


def ols_fwl(
    sigma: np.ndarray, idx: dict[str, int], x: str, y: str, z: frozenset, n: int
) -> tuple[float, float]:
    """The OLS coefficient of ``x`` in the regression of ``y`` on ``{x} + z``, with its SE.

    Both from a sample covariance matrix. The standard error is the
    Frisch-Waugh-Lovell form ``sqrt(s2_{y|x,z} / (n * s2_{x|z}))``, the same
    partialling-out formula the repository's asymptotic variance uses, with the
    sample moments in place of the population ones.
    """
    reg = [x, *sorted(set(z) - {x, y})]
    ri = [idx[a] for a in reg]
    yi = idx[y]
    a = sigma[np.ix_(ri, ri)]
    b = sigma[np.ix_(ri, [yi])].ravel()
    try:
        coef = np.linalg.solve(a, b)
    except np.linalg.LinAlgError:
        coef = np.linalg.lstsq(a, b, rcond=None)[0]
    beta = float(coef[0])
    s_y = float(sigma[yi, yi] - b @ coef)
    zi = ri[1:]
    if zi:
        az = sigma[np.ix_(zi, zi)]
        bz = sigma[np.ix_(zi, [idx[x]])].ravel()
        try:
            cz = np.linalg.solve(az, bz)
        except np.linalg.LinAlgError:
            cz = np.linalg.lstsq(az, bz, rcond=None)[0]
        s_x = float(sigma[idx[x], idx[x]] - bz @ cz)
    else:
        s_x = float(sigma[idx[x], idx[x]])
    if s_x <= 0 or s_y <= 0:
        return beta, float("inf")
    return beta, float(np.sqrt(s_y / (n * s_x)))


LEDGER_FIELDS = (
    "network",
    "tier",
    "x",
    "y",
    "arm",
    "rep",
    "true_dag_on_path",
    "chat_source_sha256",
    "generating_dag_sha256",
    "k_sha256",
    "n_k",
    "status",
    "committed_verdict",
    "z_size",
    "g0_undirected",
    "k_g0",
    "separation",
    "separation_status",
    "component_size",
    "r_claim",
    "r_claim_status",
    "r_id",
    "r_hop",
    "hop_method",
    "hop_exact",
    "hop_lower_bound",
    "phi_1",
    "witness",
    "n_wrong_claims",
    "n_omitted",
    "z_valid_at_truth",
    "z_backdoor_at_truth",
    "severity",
    *[f"bucket_{r}" for r in RUNGS],
)
PANEL_FIELDS = (
    "cpdag_source",
    "truth_status",
    "truth_extends_g0",
    "n_claims_on_absent_edges",
    "y_true_descendant",
    "r_cert_claim",
    "panel_n",
    "beta_hat",
    "true_effect",
    "bias_abs",
    "se_fwl",
    "band_lo",
    "band_hi",
    "band_covers_truth",
)


def stage_sweep(rand_reps: int, seed: int, radius_limit: float) -> None:
    """Every arm on every frame row: instruments, freeze and hash, then the audit."""
    chat = json.loads((OUT / "chat_e.json").read_text())
    know = json.loads((OUT / "knowledge.json").read_text())
    frame = [json.loads(line) for line in (OUT / "frame.jsonl").open()]
    by_net: dict[str, list] = collections.defaultdict(list)
    for r in frame:
        by_net[r["network"]].append(r)
    arms = ["D_LLM", "D_RAND", "D_DEGEN", "A_TRUE", "A_ORACLE"]
    source = f"{chat['meta']['method']}_{chat['meta']['indep_test']}_n{chat['meta']['n']}"
    t0 = time.perf_counter()

    # ---- pass 1: the instruments, from observables; the truth enters only to
    # build the two audit arms' knowledge, which is what `via_K` marks
    pending: list[tuple[dict, list, frozenset | None, dict | None, MPDAG | None]] = []
    substrate: dict[str, tuple] = {}
    arm_notes: dict[str, dict] = {}
    for net in sorted(by_net):
        dag, weights, noise, names = load_gaussian(net)
        cpdag = graph_from_json(chat["networks"][net])
        if sha(cpdag) != chat["networks"][net]["sha256"]:
            raise SystemExit(f"{net}: reference CPDAG hash mismatch")
        substrate[net] = (dag, weights, noise, cpdag)
        dag_sha, chat_sha = sha(dag), sha(cpdag)
        entry = know["D_LLM"]["networks"].get(net, {})
        base_llm = [tuple(e) for e in entry.get("k", [])]
        plans: list[tuple[str, int, list, str]] = [("D_LLM", 0, base_llm, "false")]
        if base_llm:
            for rep in range(rand_reps):
                rng = np.random.default_rng([seed, stable_seed(net), rep])
                kk = chance_level(cpdag, base_llm, rng)
                if kk is not None:
                    plans.append(("D_RAND", rep, kk, "false"))
        plans.append(("D_DEGEN", 0, [], "false"))
        # A_TRUE: the truthful orientation of the model's edges that exist in the
        # truth; a claim on an absent edge has no truthful orientation and is dropped
        truthful = [
            (a, b) if dag.is_directed_edge(a, b) else (b, a)
            for (a, b) in base_llm
            if dag.has_edge(a, b)
        ]
        plans.append(("A_TRUE", 0, truthful, "via_K"))
        # A_ORACLE: the truthful orientation of every undirected edge of C-hat_E
        # that exists in the truth
        oracle_k = [
            (a, b) if dag.is_directed_edge(a, b) else (b, a)
            for (a, b) in sorted(cpdag.undirected_edges)
            if dag.has_edge(a, b)
        ]
        plans.append(("A_ORACLE", 0, oracle_k, "via_K"))
        true_skel, est_skel = dag.skeleton(), cpdag.skeleton()
        arm_notes[net] = {
            "llm_claims": len(base_llm),
            "llm_claims_on_absent_edges": len(base_llm) - len(truthful),
            "chat_undirected": len(cpdag.undirected_edges),
            "chat_undirected_in_truth": len(oracle_k),
            "chat_directed": len(cpdag.directed_edges),
            "chat_directed_wrong_way": sum(
                1 for (a, b) in cpdag.directed_edges if dag.is_directed_edge(b, a)
            ),
            "chat_directed_absent": sum(
                1 for (a, b) in cpdag.directed_edges if not dag.has_edge(a, b)
            ),
            "chat_skeleton_missing": len(true_skel - est_skel),
            "chat_skeleton_extra": len(est_skel - true_skel),
        }
        for arm, rep, k, provenance in plans:
            k_sha = hashlib.sha256(json.dumps(sorted(k)).encode()).hexdigest()[:16]
            g0 = apply_orientations(cpdag, k)
            for row in by_net[net]:
                x, y = row["x"], row["y"]
                rec: dict = {
                    "network": net,
                    "tier": TIER.get(net, "T1"),
                    "x": x,
                    "y": y,
                    "arm": arm,
                    "rep": rep,
                    "true_dag_on_path": provenance,
                    "chat_source_sha256": chat_sha,
                    "generating_dag_sha256": dag_sha,
                    "k_sha256": k_sha,
                    "n_k": len(k),
                    "cpdag_source": source,
                }
                if g0 is None:
                    rec["status"] = "knowledge_inconsistent"
                    pending.append((rec, k, None, None, None))
                    continue
                irec, z, cert = instrument_row(cpdag, g0, k, x, y, names, radius_limit)
                rec.update(irec)
                pending.append((rec, k, z, cert, g0))
        print(
            f"  {net}: {len(pending)} instrument rows, {time.perf_counter() - t0:.0f} s",
            flush=True,
        )

    # ---- the freeze: every instrument row written and hashed before the audit
    frozen_text = json.dumps([rec for rec, *_ in pending], indent=1, sort_keys=True)
    (OUT / "frozen_rows.json").write_text(frozen_text)
    frozen_sha = hashlib.sha256(frozen_text.encode()).hexdigest()
    (OUT / "frozen_rows.sha256").write_text(frozen_sha + "\n")
    print(f"frozen {len(pending)} rows, sha256 {frozen_sha}", flush=True)

    # ---- pass 2: the truth revealed
    fields = [*LEDGER_FIELDS, *PANEL_FIELDS]
    fh = (OUT / "rows.csv").open("w", newline="")
    w = csv.DictWriter(fh, fieldnames=fields)
    w.writeheader()
    counts: dict[str, collections.Counter] = collections.defaultdict(collections.Counter)
    truth_counts: dict[str, collections.Counter] = collections.defaultdict(collections.Counter)
    sems: dict[str, tuple] = {}
    for rec, k, z, cert, g0 in pending:
        net, x, y = rec["network"], rec["x"], rec["y"]
        dag, weights, noise, cpdag = substrate[net]
        # the premise of the ledger's structural argument, tested exactly: is the
        # generating DAG a consistent extension of the analyst's graph?
        rec["truth_extends_g0"] = int(g0 is not None and is_consistent_extension(dag, g0))
        absent = sum(1 for (a, b) in k if not dag.has_edge(a, b))
        y_desc = y in dag.descendants(x)
        flags = []
        if absent:
            flags.append("claim_on_absent_edge")
        if not y_desc:
            flags.append("y_not_true_descendant")
        rec["truth_status"] = ";".join(flags) or "ok"
        rec["n_claims_on_absent_edges"] = absent
        rec["y_true_descendant"] = int(y_desc)
        truth_counts[rec["arm"]][rec["truth_status"]] += 1
        counts[rec["arm"]][rec["status"]] += 1
        if z is None or cert is None:
            w.writerow(rec)
            continue
        rec.update(audit_row(dag, cpdag, k, x, y, z, cert))
        floor = cert["B7_claim"]
        rec["r_cert_claim"] = UNREACHED if floor is None else floor
        # A-CAL-2: the adjustment estimate on the reference sample against the true effect
        if net not in sems:
            ref_seed = chat["networks"][net]["seed"]
            data = sample_for(dag, weights, noise, net, REFERENCE["n"], ref_seed)
            if data_sha(data) != chat["networks"][net]["sample_sha256"]:
                raise SystemExit(f"{net}: reference sample did not regenerate")
            sems[net] = (
                LinearSEM(dag=dag, weights=weights, noise_var=noise),
                np.cov(data, rowvar=False),
                {v: i for i, v in enumerate(dag.nodes)},
            )
        sem, sigma, idx = sems[net]
        beta, se = ols_fwl(sigma, idx, x, y, z, REFERENCE["n"])
        truth = sem.true_total_effect(x, y)
        rec.update(
            {
                "panel_n": REFERENCE["n"],
                "beta_hat": round(beta, 6),
                "true_effect": round(truth, 6),
                "bias_abs": round(abs(beta - truth), 6),
                "se_fwl": round(se, 6),
                "band_lo": round(beta - Z975 * se, 6),
                "band_hi": round(beta + Z975 * se, 6),
                "band_covers_truth": int(beta - Z975 * se <= truth <= beta + Z975 * se),
            }
        )
        w.writerow(rec)
    fh.close()

    summary = {
        "frame_rows": len(frame),
        "networks": len(by_net),
        "arms": arms,
        "evaluations": len(pending),
        "seconds": round(time.perf_counter() - t0, 1),
        "frozen_rows_sha256": frozen_sha,
        "status_by_arm": {a: dict(c) for a, c in sorted(counts.items())},
        "truth_status_by_arm": {a: dict(c) for a, c in sorted(truth_counts.items())},
        "knowledge_per_network": arm_notes,
        "settings": {
            "rand_reps": rand_reps,
            "seed": seed,
            "radius_limit_s": radius_limit,
            "panel_n": REFERENCE["n"],
            "cpdag_source": source,
            "instruments": "ledger_sweep.instrument_row, unchanged",
        },
    }
    (OUT / "sweep_summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps({k: summary[k] for k in ("status_by_arm", "truth_status_by_arm")}, indent=1))


# --- entry point -------------------------------------------------------------------


def main(argv: list[str] | None = None) -> None:
    """Run one stage of the panel, or all of them in order."""
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "stage",
        choices=[
            "structure",
            "structure-retry",
            "frame",
            "questionnaire",
            "elicit",
            "sweep",
            "all",
        ],
    )
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--cap", type=int, default=20)
    ap.add_argument("--rand-reps", type=int, default=3)
    ap.add_argument("--seed", type=int, default=20260913)
    ap.add_argument("--radius-limit", type=float, default=5.0)
    ap.add_argument("--export-only", action="store_true", help="export prompts, never call")
    args = ap.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    stages = (
        ["structure", "frame", "questionnaire", "elicit", "sweep"]
        if args.stage == "all"
        else [args.stage]
    )
    for stage in stages:
        print(f"\n== {stage}", flush=True)
        if stage == "structure":
            stage_structure(args.workers)
        elif stage == "structure-retry":
            stage_structure_retry(args.workers)
        elif stage == "frame":
            stage_frame(args.cap)
        elif stage == "questionnaire":
            stage_questionnaire()
        elif stage == "elicit":
            if not stage_elicit(args.export_only):
                return
        elif stage == "sweep":
            stage_sweep(args.rand_reps, args.seed, args.radius_limit)


if __name__ == "__main__":
    main()
