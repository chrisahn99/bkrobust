"""The matched control for reviewer Point 7.

The paper's final table (``experiments/final_table.py`` -> ``results/final_table/``)
observes that the validity radius r_val concentrates at 1 on real networks at
LLM-elicited background knowledge (46 of 53 *informative* instances -- see
``results/final_table/instances.jsonl``, reproduced and re-verified by this
module). The reviewer's objection: concentration at 1 might be a property of
the real-network CPDAGs, of the *number* of claims asserted (|K|), or of the
elicitation format -- not of the model's knowledge *content*. This module
builds the two controls needed to tell those apart, on the SAME (network, x, y)
queries D_LLM used, and compares all arms.

Arms
----
D_LLM             -- treatment. Reused verbatim from results/final_table/
                      (committed; never re-run or modified by this module).
D_SCRAMBLED,
D_SCRAMBLED_72B   -- semantic control. Same questionnaire, same graphs,
                      variable names scrambled so domain knowledge cannot be
                      exploited. Run via experiments/final_table.py itself, as
                      an unmodified subprocess (see run_table_arm below), with
                      identical --seed/--queries-per-network/--skip.
D_RANDOM_MATCHED  -- quantity control (implemented here). For each network,
                      draw |K| orientations uniformly at random from the
                      CPDAG's undirected edges, |K| set EXACTLY to that
                      network's D_LLM count, redraw on Meek-inconsistency, and
                      run the SAME certify() pipeline final_table.py uses.
                      >= 20 independent seeds per network; every seed and every
                      rejected draw is recorded (see draws.jsonl per network).

Everything here IMPORTS final_table.py's functions/CLI; it never edits that
file, src/bkrobust/hybrid.py, experiments/elicit_*.py, results/final_table/ or
results/elicit/ -- see the module docstring in NOTES.md for the standing
precedent (commit 748f2e1) this respects.

Interpreter note (IMPORTANT -- read before touching this file)
----------------------------------------------------------------
This module, and every arm it runs, uses ``.venv/bin/python`` (3.14.7,
numpy 2.5.3, scipy 1.18.1), NOT the repo's usual pinned ``/usr/bin/python3``
(3.9.6). Why: src/bkrobust/epsilon/certify.py uses zip(..., strict=True), a
Python 3.10+ keyword, which raises under 3.9.6. Git history settles which
interpreter actually produced the committed D_LLM table: certify.py's
strict=True line landed in commit 4276d15c (2026-09-17 11:38:10+02:00) and
results/final_table/instances.jsonl was committed afterwards, the same day,
in commit b546cedd (2026-09-17 22:35:36+02:00) -- so that run necessarily
executed the current certify.py and therefore ran on a 3.10+ interpreter, not
3.9.6. Confirmed directly: replaying the committed D_LLM instance
(Acid_1996, x1->x10, seed 20260917) under ``PYTHONPATH=src .venv/bin/python``,
with NO shim of any kind, reproduces r_val=1, n_knowledge=1,
theta_z=0.11769815101763414 and the full r_eps grid exactly.

An earlier version of this module instead patched Python 3.9.6's zip() at the
builtins level (in-process) and via a sitecustomize.py shim (for the
final_table.py subprocess), reasoning -- wrongly -- that a byte-for-byte match
on one instance under the shimmed 3.9.6 was sufficient. It was abandoned
because it put the treatment arm (committed D_LLM, run on a 3.10+ interpreter
with numpy 2.5.3) and the control arms (3.9.6 + a builtins patch, numpy 2.0.2)
on two different interpreters and numpy majors -- exactly the confound a
matched-control study must not have. There is no shim in this file any more;
every arm, including D_LLM's re-verification, runs on .venv/bin/python, and
the interpreter/numpy/scipy versions are recorded per arm in manifest.json.

Usage::

    PYTHONPATH=src:experiments .venv/bin/python experiments/matched_control.py verify
    PYTHONPATH=src:experiments .venv/bin/python experiments/matched_control.py \\
        table-arm --condition D_SCRAMBLED
    PYTHONPATH=src:experiments .venv/bin/python experiments/matched_control.py \\
        random --networks net1,net2 --seeds 20 --shard-tag a
    PYTHONPATH=src:experiments .venv/bin/python experiments/matched_control.py merge-random
    PYTHONPATH=src:experiments .venv/bin/python experiments/matched_control.py compare
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
EXPDIR = ROOT / "experiments"
for p in (SRC, EXPDIR):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

import numpy as np  # noqa: E402
import scipy  # noqa: E402

import final_table as ft  # noqa: E402  (experiments/final_table.py: imported, never modified)
from bkrobust.core.conventions import UNREACHED  # noqa: E402
from bkrobust.demo.evaluate import random_sem  # noqa: E402
from bkrobust.demo.graph import MPDAG  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402

OUT_ROOT = ROOT / "results" / "matched_control"
ARMS_DIR = OUT_ROOT / "arms"
LOGS_DIR = OUT_ROOT / "_logs"
COMPARISON_DIR = OUT_ROOT / "comparison"
KNOWLEDGE_PATH = ROOT / "results" / "elicit" / "knowledge.json"
FINAL_TABLE_D_LLM = ROOT / "results" / "final_table"
PYTHON_BIN = str(ROOT / ".venv" / "bin" / "python")


def interpreter_versions() -> dict[str, str]:
    """Python/numpy/scipy versions of THIS process, for manifest.json provenance."""
    return {
        "python": platform.python_version(),
        "python_executable": sys.executable,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
    }


DEFAULT_SEED = 20260917
DEFAULT_QPN = 5
DEFAULT_SKIP = "pathfinder"
DEFAULT_TIME_LIMIT = 120.0
DEFAULT_SEARCH_BUDGET = 3
RANDOM_CONDITION = "D_RANDOM_MATCHED"


def _log(msg: str, log_path: Path | None = None) -> None:
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, file=sys.stderr, flush=True)
    if log_path is not None:
        with log_path.open("a") as f:
            f.write(line + "\n")


# --------------------------------------------------------------------------
# Outcome taxonomy, matching final_table.py / the paper's own convention
# exactly (verified against results/final_table/instances.jsonl: of 160
# committed D_LLM instances, 66 are "ok & r_val != 0"; of those, 13 are
# UNREACHED (r_val == -1, robust throughout the reachable space) and the
# remaining 53 are genuinely finite -- 46 of which sit at r_val == 1. That
# 46/53 is exactly the paper's headline figure, so "informative" below means
# finite, non-degenerate r_val, with UNREACHED broken out as its own bucket.
# --------------------------------------------------------------------------
def classify(record: dict[str, Any]) -> str:
    """One of: informative, unreached, degenerate, timeout, error."""
    if record["status"] == "timeout":
        return "timeout"
    if record["status"] == "error":
        return "error"
    if record["status"] != "ok":
        return "error"
    r_val = record["r_val"]
    if r_val == 0:
        return "degenerate"
    if r_val == UNREACHED:
        return "unreached"
    return "informative"


# --------------------------------------------------------------------------
# verify: confirm the arms are matched by construction before spending any
# compute on them.
# --------------------------------------------------------------------------
def cmd_verify(args: argparse.Namespace) -> None:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOGS_DIR / "verify.log"

    knowledge_raw = json.loads(KNOWLEDGE_PATH.read_text())
    conditions = ["D_LLM", "D_SCRAMBLED", "D_SCRAMBLED_72B"]
    qshas = {c: knowledge_raw[c]["questionnaire_sha256"] for c in conditions}
    _log(f"questionnaire_sha256 by condition: {qshas}", log_path)
    assert len(set(qshas.values())) == 1, f"questionnaire hashes differ: {qshas}"
    _log("PASS: D_LLM, D_SCRAMBLED, D_SCRAMBLED_72B share one questionnaire hash.", log_path)

    # Queries are loaded from frame.jsonl only -- no condition argument at all
    # in load_queries's signature -- so they are identical across conditions
    # by construction. Confirm the file itself is untouched (sha256) and that
    # load_queries is deterministic.
    frame_sha = hashlib.sha256(ft.FRAME_PATH.read_bytes()).hexdigest()
    q1 = ft.load_queries(DEFAULT_QPN)
    q2 = ft.load_queries(DEFAULT_QPN)
    assert q1 == q2, "load_queries is not deterministic"
    _log(f"PASS: frame.jsonl sha256={frame_sha}; load_queries({DEFAULT_QPN}) deterministic, "
         f"{sum(len(v) for v in q1.values())} queries over {len(q1)} networks.", log_path)

    # Networks: same CPDAGs regardless of condition (parsing does not read
    # knowledge.json at all).
    knowledge = {c: ft.load_knowledge(c) for c in conditions}
    net_sets = {c: set(knowledge[c]) for c in conditions}
    _log(f"network counts per condition: {[(c, len(net_sets[c])) for c in conditions]}", log_path)
    assert net_sets["D_LLM"] == net_sets["D_SCRAMBLED"] == net_sets["D_SCRAMBLED_72B"], (
        "conditions do not cover the same networks"
    )
    _log("PASS: all three conditions cover the identical 32-network universe.", log_path)

    # |K| per network per condition, for the record (this is exactly why the
    # quantity control is needed -- D_LLM and D_SCRAMBLED are NOT matched on
    # |K|).
    k_counts = {
        c: {net: len(v) for net, v in knowledge[c].items()} for c in conditions
    }
    totals = {c: sum(k_counts[c].values()) for c in conditions}
    _log(f"total |K| per condition: {totals}", log_path)

    skip = {s.strip() for s in DEFAULT_SKIP.split(",") if s.strip()}
    parsed = ft.load_networks(names=None, skip=skip, max_nodes=0)
    universe = sorted(set(knowledge["D_LLM"]) & set(q1) & set(parsed) - skip)
    _log(f"final universe (knowledge & queries & parsed - skip): {len(universe)} networks", log_path)

    undirected_counts = {}
    insufficient = []
    for net in universe:
        cpdag = parsed[net]["cpdag"]
        n_undirected = len(cpdag.undirected_edges)
        target_k = k_counts["D_LLM"][net]
        undirected_counts[net] = {
            "n_nodes": len(cpdag.nodes),
            "n_undirected": n_undirected,
            "target_k_D_LLM": target_k,
            "sufficient": target_k <= n_undirected,
        }
        if target_k > n_undirected:
            insufficient.append(net)
    if insufficient:
        _log(f"WARNING: networks where |K|_D_LLM exceeds |undirected edges|: {insufficient}", log_path)
    else:
        _log("PASS: every network has enough undirected CPDAG edges to draw |K|_D_LLM from "
             "without replacement.", log_path)

    # SEM seeding: derive_seed does not depend on condition either, so the SEM
    # (and hence theta_z, covariance) is identical across arms for the same
    # (network, x, y) -- confirmed by inspecting final_table.derive_seed's
    # signature (base_seed, network, x, y only).
    s1 = ft.derive_seed(DEFAULT_SEED, "asia", "x", "y")
    s2 = ft.derive_seed(DEFAULT_SEED, "asia", "x", "y")
    assert s1 == s2
    _log(f"PASS: derive_seed(base_seed, network, x, y) has no condition argument -- the SEM is "
         f"identical across arms for a matched query by construction (sample: {s1}).", log_path)

    # Interpreter provenance: re-verify, on THIS process's interpreter, that the
    # committed D_LLM instance reproduces exactly. See the module docstring's
    # "Interpreter note" for why this matters (a prior version of this check
    # passed under a shimmed 3.9.6 while the committed run actually executed on
    # a 3.10+ interpreter -- a confound between the treatment and control arms
    # that this repeated check exists to catch).
    from bkrobust.benchmarks.describe import parse_file, to_mpdag
    from bkrobust.demo.example import dag_to_cpdag
    from bkrobust.epsilon.certify import certify

    acid_path = next(p for p in sorted(ft.NETWORKS_DIR.iterdir()) if "Acid" in p.name)
    acid_parsed = parse_file(acid_path, hashlib.sha256(acid_path.read_bytes()).hexdigest())
    acid_dag = to_mpdag(acid_parsed)
    acid_cpdag = dag_to_cpdag(acid_dag)
    acid_seed = ft.derive_seed(DEFAULT_SEED, "Acid_1996", "x1", "x10")
    acid_sem = random_sem(acid_dag, np.random.default_rng(acid_seed))
    acid_cert = certify(
        acid_cpdag, [("x4", "x1")], "x1", "x10", sem=acid_sem,
        epsilons=(0.01, 0.05, 0.1, 0.25, 0.5, 1.0), units="relative",
        method="semilocal", search_budget=DEFAULT_SEARCH_BUDGET,
    )
    expected = {"r_val": 1, "n_knowledge": 1, "theta_z": 0.11769815101763414}
    reproduces = (
        acid_cert.r_val == expected["r_val"]
        and acid_cert.n_knowledge == expected["n_knowledge"]
        and abs(acid_cert.theta_z - expected["theta_z"]) < 1e-9
        and all(v == 1 for v in acid_cert.r_eps.values())
    )
    versions = interpreter_versions()
    _log(f"interpreter provenance check on {versions}: reproduces committed "
         f"Acid_1996 x1->x10 instance = {reproduces} "
         f"(got r_val={acid_cert.r_val} n_knowledge={acid_cert.n_knowledge} "
         f"theta_z={acid_cert.theta_z} r_eps={acid_cert.r_eps})", log_path)
    if not reproduces:
        raise SystemExit(
            "FATAL: this interpreter does NOT reproduce the committed D_LLM "
            "instance -- stop and escalate before running any arm. See "
            f"{log_path}."
        )
    _log("PASS: committed D_LLM instance reproduces exactly on this interpreter, no shim.", log_path)

    manifest = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "questionnaire_sha256": qshas["D_LLM"],
        "frame_jsonl_sha256": frame_sha,
        "base_seed": DEFAULT_SEED,
        "queries_per_network": DEFAULT_QPN,
        "skip": DEFAULT_SKIP,
        "time_limit_s": DEFAULT_TIME_LIMIT,
        "search_budget": DEFAULT_SEARCH_BUDGET,
        "universe": universe,
        "n_networks": len(universe),
        "k_counts": k_counts,
        "k_totals": totals,
        "undirected_edge_audit": undirected_counts,
        "insufficient_edges_for_quantity_control": insufficient,
        "arms": {
            "D_LLM": {
                "role": "treatment",
                "source": "results/final_table/ (committed; reused verbatim, not re-run)",
                "interpreter_versions": (
                    "not re-run; provenance established by git history (see "
                    "environment_note) and by the interpreter_provenance_check "
                    "above, run on THIS process's interpreter"
                ),
            },
            "D_SCRAMBLED": {
                "role": "semantic control",
                "source": "results/matched_control/arms/D_SCRAMBLED/ (this module, via final_table.py subprocess)",
                "interpreter_versions": versions,
            },
            "D_SCRAMBLED_72B": {
                "role": "semantic control (larger model)",
                "source": "results/matched_control/arms/D_SCRAMBLED_72B/ (this module, via final_table.py subprocess)",
                "interpreter_versions": versions,
            },
            RANDOM_CONDITION: {
                "role": "quantity control (|K| matched to D_LLM per network)",
                "source": "results/matched_control/arms/D_RANDOM_MATCHED/ (this module, custom loop)",
                "interpreter_versions": versions,
            },
        },
        "environment_note": (
            "The repo's usual pinned interpreter, /usr/bin/python3 3.9.6, cannot "
            "run src/bkrobust/epsilon/certify.py as committed to this tree "
            "(zip(..., strict=True) is Python 3.10+). Git history shows "
            "certify.py's strict=True line (commit 4276d15c, "
            "2026-09-17 11:38:10+02:00) predates the commit of "
            "results/final_table/instances.jsonl (b546cedd, "
            "2026-09-17 22:35:36+02:00), so the committed D_LLM run necessarily "
            "executed on a 3.10+ interpreter, not 3.9.6. All arms in this "
            "study -- D_LLM's re-verification, D_SCRAMBLED, D_SCRAMBLED_72B and "
            "D_RANDOM_MATCHED -- therefore run on PYTHONPATH=src .venv/bin/python "
            "(see interpreter_versions below), with no shim of any kind. "
            "Re-verified on this process before writing this manifest: the "
            "committed D_LLM instance (Acid_1996, x1, x10; seed 20260917) "
            "reproduces r_val=1, n_knowledge=1, theta_z=0.11769815101763414 and "
            "the full r_eps grid exactly under this interpreter."
        ),
        "interpreter_versions": versions,
    }
    (OUT_ROOT / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True))
    _log(f"wrote {OUT_ROOT / 'manifest.json'}", log_path)


# --------------------------------------------------------------------------
# table-arm: run D_SCRAMBLED / D_SCRAMBLED_72B (or any other knowledge.json
# condition) through final_table.py itself, unmodified, as a subprocess.
# --------------------------------------------------------------------------
def cmd_table_arm(args: argparse.Namespace) -> None:
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    out_dir = ARMS_DIR / args.condition
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = LOGS_DIR / f"table_arm_{args.condition}.log"

    env = dict(os.environ)
    env["PYTHONPATH"] = str(SRC)
    cmd = [
        PYTHON_BIN,
        str(EXPDIR / "final_table.py"),
        "--condition", args.condition,
        "--seed", str(args.seed),
        "--queries-per-network", str(args.queries_per_network),
        "--skip", args.skip,
        "--time-limit", str(args.time_limit),
        "--out", str(out_dir.relative_to(ROOT)),
    ]
    _log(f"launching (unmodified final_table.py) on {PYTHON_BIN}: {' '.join(cmd)}", log_path)
    t0 = time.perf_counter()
    with log_path.open("a") as logf:
        proc = subprocess.run(cmd, cwd=ROOT, env=env, stdout=logf, stderr=subprocess.STDOUT)
    elapsed = time.perf_counter() - t0
    _log(f"final_table.py --condition {args.condition} exited {proc.returncode} in {elapsed:.1f}s", log_path)
    if proc.returncode != 0:
        raise SystemExit(f"table-arm {args.condition} failed; see {log_path}")


# --------------------------------------------------------------------------
# random: the quantity control this module implements directly.
# --------------------------------------------------------------------------
def _draw_seed(base_seed: int, net: str, seed_idx: int) -> int:
    payload = f"{base_seed}:randomK:{net}:{seed_idx}".encode()
    digest = hashlib.sha256(payload).digest()
    return int.from_bytes(digest[:8], "big") % (2**32)


def draw_random_k(
    rng: np.random.Generator,
    undirected_edges: list[tuple[str, str]],
    target_k: int,
    cpdag: MPDAG,
    max_attempts: int,
) -> tuple[list[tuple[str, str]] | None, int, int]:
    """Draw ``target_k`` random orientations from ``undirected_edges``, Meek-checked.

    Args:
        rng: Seeded generator for this (network, seed_idx).
        undirected_edges: The CPDAG's undirected edges, canonically ordered and
            sorted (so the draw is reproducible from ``rng`` alone).
        target_k: Exactly how many orientations to draw (= |K|_D_LLM for this
            network).
        cpdag: The estimated CPDAG, for the Meek-consistency check.
        max_attempts: Redraws before giving up.

    Returns:
        ``(k, n_rejected, n_attempts)``. ``k is None`` iff every attempt up to
        ``max_attempts`` was Meek-inconsistent -- ``n_rejected == max_attempts``
        in that case.
    """
    n = len(undirected_edges)
    for attempt in range(1, max_attempts + 1):
        idx = rng.choice(n, size=target_k, replace=False)
        chosen_edges = [undirected_edges[i] for i in idx]
        k: list[tuple[str, str]] = []
        for a, b in chosen_edges:
            if rng.random() < 0.5:
                k.append((a, b))
            else:
                k.append((b, a))
        g0 = apply_orientations(cpdag, k)
        if g0 is not None:
            return k, attempt - 1, attempt
    return None, max_attempts, max_attempts


def cmd_random(args: argparse.Namespace) -> None:
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    out_dir = ARMS_DIR / RANDOM_CONDITION
    out_dir.mkdir(parents=True, exist_ok=True)
    tag = args.shard_tag
    log_path = LOGS_DIR / f"random_{tag}.log"
    instances_path = out_dir / f"shard_{tag}.instances.jsonl"
    draws_path = out_dir / f"shard_{tag}.draws.jsonl"

    skip = {s.strip() for s in args.skip.split(",") if s.strip()}
    requested = {n.strip() for n in args.networks.split(",") if n.strip()} if args.networks else None

    knowledge = ft.load_knowledge("D_LLM")
    queries = ft.load_queries(args.queries_per_network)
    parsed = ft.load_networks(names=requested, skip=skip, max_nodes=0)
    universe = sorted(set(knowledge) & set(queries) & set(parsed) - skip)
    if requested is not None:
        universe = [n for n in universe if n in requested]

    _log(
        f"shard={tag} networks={universe} seeds={args.seeds} "
        f"max_attempts={args.max_attempts} time_limit={args.time_limit} "
        f"interpreter={interpreter_versions()}",
        log_path,
    )
    (out_dir / f"shard_{tag}.meta.json").write_text(json.dumps({
        "shard_tag": tag, "networks": universe, "seeds": args.seeds,
        "max_attempts": args.max_attempts, "time_limit": args.time_limit,
        "interpreter_versions": interpreter_versions(),
    }, indent=2))

    t_start = time.perf_counter()
    with instances_path.open("w") as inst_f, draws_path.open("w") as draw_f:
        for net in universe:
            dag = parsed[net]["dag"]
            cpdag = parsed[net]["cpdag"]
            target_k = len(knowledge[net])
            undirected_edges = sorted(cpdag.undirected_edges)
            net_queries = queries[net]
            _log(
                f"  {net}: |V|={len(dag.nodes)} target_k={target_k} "
                f"n_undirected={len(undirected_edges)} n_queries={len(net_queries)} "
                f"elapsed={time.perf_counter() - t_start:.1f}s",
                log_path,
            )
            for seed_idx in range(args.seeds):
                draw_seed = _draw_seed(args.seed, net, seed_idx)
                rng = np.random.default_rng(draw_seed)
                k, n_rejected, n_attempts = draw_random_k(
                    rng, undirected_edges, target_k, cpdag, args.max_attempts
                )
                draw_record = {
                    "network": net,
                    "seed_idx": seed_idx,
                    "draw_seed": draw_seed,
                    "target_k": target_k,
                    "n_undirected": len(undirected_edges),
                    "n_rejected": n_rejected,
                    "n_attempts": n_attempts,
                    "exhausted": k is None,
                    "k": k,
                }
                draw_f.write(json.dumps(draw_record) + "\n")
                draw_f.flush()
                if k is None:
                    _log(
                        f"    {net} seed_idx={seed_idx}: EXHAUSTED after {n_attempts} attempts "
                        f"(all Meek-inconsistent) -- no instance run for this draw.",
                        log_path,
                    )
                    continue
                _log(
                    f"    {net} seed_idx={seed_idx}: accepted draw after "
                    f"{n_attempts} attempt(s) ({n_rejected} rejected)",
                    log_path,
                )
                for x, y in net_queries:
                    # Same SEM seed final_table.py would use for this (net, x, y)
                    # -- condition-independent by construction (see cmd_verify) --
                    # so the covariance is identical to D_LLM/D_SCRAMBLED for the
                    # matched query; only K differs.
                    instance_seed = ft.derive_seed(args.seed, net, x, y)
                    record = ft.run_instance(
                        RANDOM_CONDITION, net, x, y, dag, cpdag, k,
                        tuple(float(e) for e in args.epsilons.split(",")),
                        instance_seed, args.time_limit, DEFAULT_SEARCH_BUDGET,
                    )
                    record["seed_idx"] = seed_idx
                    record["draw_seed"] = draw_seed
                    record["n_rejected_draws"] = n_rejected
                    record["k"] = k
                    inst_f.write(json.dumps(record) + "\n")
                    inst_f.flush()
                    _log(
                        f"      {net} seed_idx={seed_idx} {x}->{y} status={record['status']} "
                        f"r_val={record['r_val']} seconds={record['seconds']:.2f}"
                        + (f" error={record['error']}" if record["error"] else ""),
                        log_path,
                    )
    _log(f"shard {tag} done in {time.perf_counter() - t_start:.1f}s -> {instances_path}", log_path)


def cmd_merge_random(args: argparse.Namespace) -> None:
    out_dir = ARMS_DIR / RANDOM_CONDITION
    shard_instance_files = sorted(out_dir.glob("shard_*.instances.jsonl"))
    shard_draw_files = sorted(out_dir.glob("shard_*.draws.jsonl"))
    instances: list[dict[str, Any]] = []
    for path in shard_instance_files:
        instances.extend(json.loads(line) for line in path.read_text().splitlines() if line)
    draws: list[dict[str, Any]] = []
    for path in shard_draw_files:
        draws.extend(json.loads(line) for line in path.read_text().splitlines() if line)
    # Deterministic merge order: by (network, seed_idx, x, y).
    instances.sort(key=lambda r: (r["network"], r["seed_idx"], r["x"], r["y"]))
    draws.sort(key=lambda r: (r["network"], r["seed_idx"]))

    (out_dir / "instances.jsonl").write_text(
        "\n".join(json.dumps(r) for r in instances) + ("\n" if instances else "")
    )
    (out_dir / "draws.jsonl").write_text(
        "\n".join(json.dumps(r) for r in draws) + ("\n" if draws else "")
    )

    n_networks = len({d["network"] for d in draws})
    n_exhausted = sum(1 for d in draws if d["exhausted"])
    n_rejected_total = sum(d["n_rejected"] for d in draws)
    n_attempts_total = sum(d["n_attempts"] for d in draws)
    summary = {
        "n_shards": len(shard_instance_files),
        "n_networks": n_networks,
        "n_draws": len(draws),
        "n_draws_exhausted": n_exhausted,
        "n_instances": len(instances),
        "rejection_rate": (n_rejected_total / n_attempts_total) if n_attempts_total else None,
        "n_rejected_draws_total": n_rejected_total,
        "n_attempts_total": n_attempts_total,
        "status_counts": dict_count([classify(r) for r in instances]),
        "interpreter_versions": interpreter_versions(),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


def dict_count(values: list[str]) -> dict[str, int]:
    out: dict[str, int] = defaultdict(int)
    for v in values:
        out[v] += 1
    return dict(out)


# --------------------------------------------------------------------------
# compare: cross-arm comparison on the matched (network, x, y) queries.
# --------------------------------------------------------------------------
def _load_jsonl(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    return [json.loads(line) for line in path.read_text().splitlines() if line]


def _paired_bootstrap(
    diffs_by_network: dict[str, list[float]],
    n_resamples: int,
    rng: np.random.Generator,
) -> dict[str, float]:
    """Cluster (by network) bootstrap CI for a mean paired difference.

    Resample networks with replacement; within a resampled network, resample
    its instance-level differences with replacement too, then average
    everything pooled. This respects the paper's own precedent (PAPER_NARRATIVE
    in final_table.py) that networks, not instance pairs, are the clustering
    unit -- an instance-level-only bootstrap would treat a 20-query network as
    20 independent votes when it is one graph's worth of evidence.
    """
    nets = sorted(diffs_by_network)
    all_diffs = [d for net in nets for d in diffs_by_network[net]]
    point = float(np.mean(all_diffs)) if all_diffs else float("nan")
    boot = []
    for _ in range(n_resamples):
        sampled_nets = rng.choice(nets, size=len(nets), replace=True)
        pooled = []
        for net in sampled_nets:
            vals = diffs_by_network[net]
            if not vals:
                continue
            idx = rng.integers(0, len(vals), size=len(vals))
            pooled.extend(vals[i] for i in idx)
        if pooled:
            boot.append(float(np.mean(pooled)))
    boot.sort()
    if not boot:
        return {"point": point, "ci_lo": float("nan"), "ci_hi": float("nan"), "n_resamples": 0}
    lo = boot[int(0.025 * len(boot))]
    hi = boot[min(len(boot) - 1, int(0.975 * len(boot)))]
    return {"point": point, "ci_lo": lo, "ci_hi": hi, "n_resamples": len(boot)}


def cmd_compare(args: argparse.Namespace) -> None:
    COMPARISON_DIR.mkdir(parents=True, exist_ok=True)
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOGS_DIR / "compare.log"

    d_llm = _load_jsonl(FINAL_TABLE_D_LLM / "instances.jsonl")
    d_scr = _load_jsonl(ARMS_DIR / "D_SCRAMBLED" / "instances.jsonl")
    d_scr72 = _load_jsonl(ARMS_DIR / "D_SCRAMBLED_72B" / "instances.jsonl")
    d_rand = _load_jsonl(ARMS_DIR / RANDOM_CONDITION / "instances.jsonl")

    _log(f"loaded: D_LLM={len(d_llm)} D_SCRAMBLED={len(d_scr)} "
         f"D_SCRAMBLED_72B={len(d_scr72)} D_RANDOM_MATCHED={len(d_rand)}", log_path)

    # --- per-arm outcome denominators, exactly as they are (never silently dropped) ---
    def denom_table(records: list[dict[str, Any]]) -> dict[str, Any]:
        classes = [classify(r) for r in records]
        counts = dict_count(classes)
        informative = [r for r in records if classify(r) == "informative"]
        at1 = sum(1 for r in informative if r["r_val"] == 1)
        return {
            "n_total": len(records),
            "counts": counts,
            "n_informative": len(informative),
            "n_at_r_val_1": at1,
            "share_at_1_of_informative": (at1 / len(informative)) if informative else None,
        }

    arm_denoms = {
        "D_LLM": denom_table(d_llm),
        "D_SCRAMBLED": denom_table(d_scr),
        "D_SCRAMBLED_72B": denom_table(d_scr72),
        RANDOM_CONDITION: denom_table(d_rand),
    }
    _log("per-arm denominators: " + json.dumps(arm_denoms, indent=2), log_path)

    # --- matched-query comparison: D_LLM vs D_SCRAMBLED, informative-in-D_LLM queries ---
    d_llm_by_q = {(r["network"], r["x"], r["y"]): r for r in d_llm}
    d_scr_by_q = {(r["network"], r["x"], r["y"]): r for r in d_scr}
    matched_keys = sorted(set(d_llm_by_q) & set(d_scr_by_q))
    _log(f"matched (network,x,y) keys present in both D_LLM and D_SCRAMBLED: {len(matched_keys)}", log_path)

    rows = []
    diffs_at1_by_net: dict[str, list[float]] = defaultdict(list)
    for key in matched_keys:
        net, x, y = key
        rl, rs = d_llm_by_q[key], d_scr_by_q[key]
        cl, cs = classify(rl), classify(rs)
        row = {
            "network": net, "x": x, "y": y,
            "D_LLM": {"status": rl["status"], "class": cl, "r_val": rl["r_val"], "len_k": rl["len_k"]},
            "D_SCRAMBLED": {"status": rs["status"], "class": cs, "r_val": rs["r_val"], "len_k": rs["len_k"]},
        }
        rows.append(row)
        if cl == "informative" and cs == "informative":
            # Both arms informative on this query, so the pair is defined. An
            # arm that is informative but not at radius 1 contributes 0.0, not
            # NaN: NaN here would propagate through the mean and void the whole
            # paired statistic, which is exactly what it did before this fix.
            at1_llm = 1.0 if rl["r_val"] == 1 else 0.0
            at1_scr = 1.0 if rs["r_val"] == 1 else 0.0
            diffs_at1_by_net[net].append(at1_llm - at1_scr)
    (COMPARISON_DIR / "matched_queries_llm_vs_scrambled.json").write_text(
        json.dumps(rows, indent=2)
    )

    rng = np.random.default_rng(args.bootstrap_seed)
    paired_at1 = _paired_bootstrap(diffs_at1_by_net, args.n_resamples, rng)
    n_paired_informative_both = sum(len(v) for v in diffs_at1_by_net.values())
    _log(f"paired share-at-1 (D_LLM - D_SCRAMBLED), informative-in-both, n={n_paired_informative_both} "
         f"pairs over {len(diffs_at1_by_net)} networks: {paired_at1}", log_path)

    # --- quantity control comparison: D_LLM vs D_RANDOM_MATCHED, per (network,x,y),
    #     averaged over the >=20 random draws for that network ---
    d_rand_by_net_q: dict[tuple[str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for r in d_rand:
        d_rand_by_net_q[(r["network"], r["x"], r["y"])].append(r)

    rand_rows = []
    diffs_at1_rand_by_net: dict[str, list[float]] = defaultdict(list)
    for key in sorted(set(d_llm_by_q) & set(d_rand_by_net_q)):
        net, x, y = key
        rl = d_llm_by_q[key]
        draws = d_rand_by_net_q[key]
        cl = classify(rl)
        draw_classes = [classify(r) for r in draws]
        n_draws = len(draws)
        n_at1 = sum(1 for r, c in zip(draws, draw_classes) if c == "informative" and r["r_val"] == 1)
        n_informative_draws = sum(1 for c in draw_classes if c == "informative")
        share_at1_draws = (n_at1 / n_informative_draws) if n_informative_draws else None
        rand_rows.append({
            "network": net, "x": x, "y": y,
            "D_LLM": {"status": rl["status"], "class": cl, "r_val": rl["r_val"], "len_k": rl["len_k"]},
            "D_RANDOM_MATCHED": {
                "n_draws": n_draws,
                "class_counts": dict_count(draw_classes),
                "n_informative_draws": n_informative_draws,
                "n_at_r_val_1": n_at1,
                "share_at_1_of_informative_draws": share_at1_draws,
                "len_k_target": draws[0]["len_k"] if draws else None,
            },
        })
        if cl == "informative" and n_informative_draws > 0:
            at1_llm = 1.0 if rl["r_val"] == 1 else 0.0
            diffs_at1_rand_by_net[net].append(at1_llm - share_at1_draws)
    (COMPARISON_DIR / "matched_queries_llm_vs_random.json").write_text(
        json.dumps(rand_rows, indent=2)
    )
    paired_at1_rand = _paired_bootstrap(diffs_at1_rand_by_net, args.n_resamples, rng)
    n_paired_rand = sum(len(v) for v in diffs_at1_rand_by_net.values())
    _log(f"paired share-at-1 (D_LLM - mean(D_RANDOM_MATCHED draws)), n={n_paired_rand} pairs over "
         f"{len(diffs_at1_rand_by_net)} networks: {paired_at1_rand}", log_path)

    # --- elicitation accuracy: compelled_control right/wrong, D_LLM vs D_SCRAMBLED(_72B) ---
    knowledge_raw = json.loads(KNOWLEDGE_PATH.read_text())

    def control_accuracy(cond: str) -> dict[str, Any]:
        nets = knowledge_raw[cond]["networks"]
        # compelled_control is absent (None) for networks where no compelled
        # question was reached, e.g. the large graphs the elicitation never got
        # through. Those contribute nothing to the tally and are counted
        # separately rather than treated as zeros, which would silently deflate
        # the accuracy denominator.
        blocks = [v.get("compelled_control") for v in nets.values()]
        present = [b for b in blocks if isinstance(b, dict)]
        n_networks_without_block = len(blocks) - len(present)

        def total(field: str) -> int:
            return sum(int(b.get(field) or 0) for b in present)

        right, wrong = total("right"), total("wrong")
        n = right + wrong
        return {
            "right": right, "wrong": wrong,
            "declined": total("declined"), "not_reached": total("not_reached"),
            "n_networks_with_compelled_control": len(present),
            "n_networks_without_compelled_control": n_networks_without_block,
            "accuracy_of_answered": (right / n) if n else None,
            "true_dag_on_path": knowledge_raw[cond].get("true_dag_on_path"),
        }

    accuracy = {c: control_accuracy(c) for c in ["D_LLM", "D_SCRAMBLED", "D_SCRAMBLED_72B"]}
    _log("compelled_control accuracy by condition: " + json.dumps(accuracy, indent=2), log_path)

    # Two-proportion (right vs wrong) bootstrap CI on the accuracy difference,
    # D_LLM minus D_SCRAMBLED, clustered by network (paired: same compelled
    # questions are asked under both conditions -- see questionnaire_build.py).
    def per_network_accuracy_diffs(cond_a: str, cond_b: str) -> dict[str, list[float]]:
        nets_a = knowledge_raw[cond_a]["networks"]
        nets_b = knowledge_raw[cond_b]["networks"]
        out: dict[str, list[float]] = {}
        for net in sorted(set(nets_a) & set(nets_b)):
            ca, cb = nets_a[net].get("compelled_control"), nets_b[net].get("compelled_control")
            # A network whose elicitation never reached a compelled question
            # carries no block under that condition; the pair is undefined there
            # and the network is skipped, never imputed as agreement.
            if not isinstance(ca, dict) or not isinstance(cb, dict):
                continue
            na = int(ca.get("right") or 0) + int(ca.get("wrong") or 0)
            nb = int(cb.get("right") or 0) + int(cb.get("wrong") or 0)
            if na == 0 or nb == 0:
                continue
            out[net] = [(int(ca["right"]) / na) - (int(cb["right"]) / nb)]
        return out

    acc_diffs = per_network_accuracy_diffs("D_LLM", "D_SCRAMBLED")
    paired_accuracy = _paired_bootstrap(acc_diffs, args.n_resamples, rng)
    _log(f"paired per-network accuracy diff (D_LLM - D_SCRAMBLED), n_networks_with_data="
         f"{len(acc_diffs)}: {paired_accuracy}", log_path)

    comparison = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "bootstrap_seed": args.bootstrap_seed,
        "n_resamples": args.n_resamples,
        "arm_denominators": arm_denoms,
        "matched_queries": {
            "D_LLM_vs_D_SCRAMBLED": {
                "n_matched_keys": len(matched_keys),
                "n_paired_both_informative": n_paired_informative_both,
                "paired_diff_share_at_1": paired_at1,
            },
            "D_LLM_vs_D_RANDOM_MATCHED": {
                "n_matched_keys": len(rand_rows),
                "n_paired_both_informative": n_paired_rand,
                "paired_diff_share_at_1": paired_at1_rand,
            },
        },
        "compelled_control_accuracy": accuracy,
        "paired_accuracy_diff_D_LLM_minus_D_SCRAMBLED": paired_accuracy,
    }
    (COMPARISON_DIR / "summary.json").write_text(json.dumps(comparison, indent=2, sort_keys=True))
    print(json.dumps(comparison, indent=2, sort_keys=True))
    _log(f"wrote {COMPARISON_DIR / 'summary.json'}", log_path)


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_verify = sub.add_parser("verify", help="confirm arms are matched by construction; write manifest.json")
    p_verify.set_defaults(func=cmd_verify)

    p_table = sub.add_parser("table-arm", help="run a knowledge.json condition via final_table.py (subprocess)")
    p_table.add_argument("--condition", required=True)
    p_table.add_argument("--seed", type=int, default=DEFAULT_SEED)
    p_table.add_argument("--queries-per-network", type=int, default=DEFAULT_QPN)
    p_table.add_argument("--skip", default=DEFAULT_SKIP)
    p_table.add_argument("--time-limit", type=float, default=DEFAULT_TIME_LIMIT)
    p_table.set_defaults(func=cmd_table_arm)

    p_random = sub.add_parser("random", help="quantity control: random-at-matched-|K|")
    p_random.add_argument("--networks", default="", help="comma list; default all")
    p_random.add_argument("--seeds", type=int, default=20)
    p_random.add_argument("--max-attempts", type=int, default=200)
    p_random.add_argument("--shard-tag", required=True)
    p_random.add_argument("--seed", type=int, default=DEFAULT_SEED)
    p_random.add_argument("--queries-per-network", type=int, default=DEFAULT_QPN)
    p_random.add_argument("--skip", default=DEFAULT_SKIP)
    p_random.add_argument("--time-limit", type=float, default=DEFAULT_TIME_LIMIT)
    p_random.add_argument("--epsilons", default=ft.DEFAULT_EPSILONS)
    p_random.set_defaults(func=cmd_random)

    p_merge = sub.add_parser("merge-random", help="merge random-arm shards deterministically")
    p_merge.set_defaults(func=cmd_merge_random)

    p_compare = sub.add_parser("compare", help="cross-arm comparison on matched queries")
    p_compare.add_argument("--bootstrap-seed", type=int, default=DEFAULT_SEED + 1)
    p_compare.add_argument("--n-resamples", type=int, default=10000)
    p_compare.set_defaults(func=cmd_compare)

    return parser


def main() -> None:
    args = build_arg_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
