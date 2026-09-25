"""Reviewer 4.7 control: matched random-knowledge accuracy comparison.

Reviewer 4.7 (ICLR 2027 v2 revision) objects that the scrambled-name control
does not separate from real-name conditions on claim accuracy, which weakens
the claim that LLM elicitation reproduces meaningful domain knowledge, and
asks for "a matched random-knowledge control".

That control already exists on disk as ``results/matched_control`` (built for
a *different* reviewer point -- whether radius-1 concentration is explained
by knowledge *content*). Its ``D_RANDOM_MATCHED`` arm draws, per network, the
same number of orientations |K| that D_LLM asserted, uniformly at random from
the CPDAG's undirected edges, Meek-checked, over 20 seeds (2 on six large/slow
networks) -- exactly the random-orientation control this task asks for. This
script does not re-run elicitation or the survival panel; it computes claim
accuracy (k_acc: share of oriented claims agreeing with the data-generating
DAG) for that existing random arm and compares it, per network, to k_acc of
the primary LLM condition (D_LLM) and to k_acc pooled over all real-name LLM
conditions, both already recorded in
``results/axis_robustness_llm/analysis_units.csv``.

Outputs (this task's own files, written fresh, nothing under
results/matched_control or results/axis_robustness_llm is modified):
    results/axis_robustness_llm_v2/controls_k_accuracy.json

Usage:
    PYTHONPATH=src:experiments .venv/bin/python experiments/llm_controls_v2.py
"""

from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "experiments"))

import final_table as ft  # noqa: E402  (experiments/final_table.py: read-only use)
from bkrobust.robustness.llm_survival import k_accuracy  # noqa: E402

MATCHED_CONTROL_DIR = REPO / "results" / "matched_control"
DRAWS_PATH = MATCHED_CONTROL_DIR / "arms" / "D_RANDOM_MATCHED" / "draws.jsonl"
ANALYSIS_UNITS_PATH = REPO / "results" / "axis_robustness_llm" / "analysis_units.csv"
OUT_DIR = REPO / "results" / "axis_robustness_llm_v2"
OUT_PATH = OUT_DIR / "controls_k_accuracy.json"

PRIMARY_CONDITION = "D_LLM"
N_BOOT = 10000
SEED = 20260925


def load_jsonl(path: Path) -> list[dict[str, Any]]:
    with path.open() as f:
        return [json.loads(line) for line in f if line.strip()]


def random_k_accuracy_by_network(dags: dict[str, dict[str, Any]]) -> dict[str, list[float]]:
    """Per-network list of per-seed k_acc for the D_RANDOM_MATCHED draws."""
    draws = load_jsonl(DRAWS_PATH)
    out: dict[str, list[float]] = defaultdict(list)
    for d in draws:
        if d.get("exhausted"):
            continue
        net = d["network"]
        if net not in dags:
            continue
        dag = dags[net]["dag"]
        k = [tuple(e) for e in d["k"]]
        acc = k_accuracy(k, dag)["k_accuracy"]
        if acc is not None:
            out[net].append(acc)
    return dict(out)


def real_condition_k_accuracy_by_network() -> dict[str, dict[str, float]]:
    """{condition: {network: k_accuracy}} read off the committed real-network panel.

    k_accuracy is a per-(condition, network) property (same claims answer
    every query in that network), so this collapses analysis_units.csv's
    per-query rows to one value per (condition, network) and asserts they
    agree, rather than silently averaging over duplicates.
    """
    with ANALYSIS_UNITS_PATH.open() as f:
        rows = list(csv.DictReader(f))
    out: dict[str, dict[str, float]] = defaultdict(dict)
    for r in rows:
        cond, net = r["condition"], r["network"]
        val = r["k_accuracy"]
        if val in ("", None):
            continue
        val = float(val)
        prev = out[cond].get(net)
        if prev is not None and abs(prev - val) > 1e-9:
            raise ValueError(f"k_accuracy disagreement within {cond}/{net}: {prev} vs {val}")
        out[cond][net] = val
    return dict(out)


def paired_network_bootstrap(
    pairs: list[tuple[str, float, float]], n_boot: int = N_BOOT, seed: int = SEED
) -> dict[str, Any]:
    """Bootstrap CI for the mean paired difference a-b, resampling networks.

    Args:
        pairs: (network, value_a, value_b) triples, one row per network.
    """
    nets = [p[0] for p in pairs]
    a = np.array([p[1] for p in pairs])
    b = np.array([p[2] for p in pairs])
    diff = a - b
    point = float(diff.mean())
    rng = np.random.default_rng(seed)
    n = len(pairs)
    boots = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        boots.append(diff[idx].mean())
    boots = np.array(boots)
    return {
        "n_networks": n,
        "networks": nets,
        "point_estimate_mean_diff": point,
        "ci_lo_2p5": float(np.percentile(boots, 2.5)),
        "ci_hi_97p5": float(np.percentile(boots, 97.5)),
        "frac_diff_gt_0": float((boots > 0).mean()),
    }


def cluster_bootstrap_mean(values_by_network: dict[str, list[float]], n_boot: int = N_BOOT, seed: int = SEED) -> dict[str, Any]:
    """Network-cluster bootstrap CI for a pooled mean (e.g. pooled k_acc)."""
    nets = sorted(values_by_network)
    rng = np.random.default_rng(seed)
    n = len(nets)
    all_vals = [v for net in nets for v in values_by_network[net]]
    point = float(np.mean(all_vals)) if all_vals else None
    boots = []
    for _ in range(n_boot):
        draw_nets = rng.choice(nets, size=n, replace=True)
        vals = [v for net in draw_nets for v in values_by_network[net]]
        if vals:
            boots.append(np.mean(vals))
    boots = np.array(boots)
    return {
        "n_networks": n,
        "n_items": len(all_vals),
        "point_estimate": point,
        "ci_lo_2p5": float(np.percentile(boots, 2.5)) if len(boots) else None,
        "ci_hi_97p5": float(np.percentile(boots, 97.5)) if len(boots) else None,
    }


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    draws = load_jsonl(DRAWS_PATH)
    networks_needed = sorted(set(d["network"] for d in draws))
    dags = ft.load_networks(names=set(networks_needed), skip=set(), max_nodes=0)
    missing = set(networks_needed) - set(dags)
    if missing:
        print(f"WARNING: could not load DAGs for {missing}", file=sys.stderr)

    random_acc = random_k_accuracy_by_network(dags)
    real_acc = real_condition_k_accuracy_by_network()

    primary_acc = real_acc.get(PRIMARY_CONDITION, {})
    real_name_conditions = [
        c for c in real_acc
        if c not in ("D_SCRAMBLED", "D_SCRAMBLED_72B")
    ]

    # --- (1) primary LLM (D_LLM) vs random, paired by network ---
    common_primary = sorted(set(primary_acc) & set(random_acc))
    pairs_primary = [
        (net, primary_acc[net], float(np.mean(random_acc[net])))
        for net in common_primary
    ]
    paired_primary = paired_network_bootstrap(pairs_primary)

    # --- (2) pooled real-name conditions vs random, paired by network ---
    # per network, pool k_acc across every real-name condition that has a
    # value there (each condition contributes one number per network, same
    # weight regardless of |K|).
    pooled_real_by_network: dict[str, list[float]] = defaultdict(list)
    for cond in real_name_conditions:
        for net, v in real_acc[cond].items():
            pooled_real_by_network[net].append(v)
    common_pooled = sorted(set(pooled_real_by_network) & set(random_acc))
    pairs_pooled = [
        (net, float(np.mean(pooled_real_by_network[net])), float(np.mean(random_acc[net])))
        for net in common_pooled
    ]
    paired_pooled = paired_network_bootstrap(pairs_pooled)

    # --- descriptive pooled means with network-cluster CIs ---
    primary_acc_by_net_list = {net: [v] for net, v in primary_acc.items()}
    desc_primary = cluster_bootstrap_mean(primary_acc_by_net_list)
    desc_pooled_real = cluster_bootstrap_mean(pooled_real_by_network)
    desc_random = cluster_bootstrap_mean(random_acc)

    out = {
        "control": "D_RANDOM_MATCHED (results/matched_control), reused verbatim -- not re-run",
        "primary_condition": PRIMARY_CONDITION,
        "real_name_conditions_pooled": real_name_conditions,
        "n_boot": N_BOOT,
        "seed": SEED,
        "descriptive": {
            "D_LLM_pooled_k_acc": desc_primary,
            "pooled_real_name_k_acc": desc_pooled_real,
            "D_RANDOM_MATCHED_k_acc": desc_random,
        },
        "paired_LLM_minus_random": paired_primary,
        "paired_pooled_real_minus_random": paired_pooled,
        "per_network_random_k_acc_n_seeds": {
            net: len(vals) for net, vals in random_acc.items()
        },
    }
    OUT_PATH.write_text(json.dumps(out, indent=2))
    print(f"wrote {OUT_PATH}")
    print(json.dumps(out["descriptive"], indent=2))
    print("paired D_LLM - random:", json.dumps(paired_primary, indent=2)[:500])
    print("paired pooled_real - random:", json.dumps(paired_pooled, indent=2)[:500])


if __name__ == "__main__":
    main()
