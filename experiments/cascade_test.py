"""The cascade test: how many edges each asserted claim compels, per supplier.

The protocol predicted, before any supplier was scored, that the greedy
recovering set used by the committed pipeline compels more edges per assertion
than an elicited set, because it is chosen to maximise what propagates. The
ledger carries both numbers on every row: the number of asserted claims and the
number of edges the closure orients that the CPDAG left open. Their ratio is
constant across the rows of one network under one arm, so the test is paired by
network. If the paired interval of the oracle minus the model contains zero,
the claim that the committed generator cascades unrealistically is withdrawn.

Writes ``results/ledger/cascade_test.json`` and prints the table.

    python experiments/cascade_test.py
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ROWS = ROOT / "results" / "ledger" / "rows.csv"
OUT = ROOT / "results" / "ledger" / "cascade_test.json"
B = 4000


def boot(diff: dict[str, float], reps: int = B, seed: int = 20260913) -> tuple[float, float, float]:
    """Mean of per-network differences with a percentile interval over networks."""
    names = sorted(diff)
    vals = np.array([diff[n] for n in names])
    rng = np.random.default_rng(seed)
    draws = np.sort(
        [vals[rng.integers(0, len(names), size=len(names))].mean() for _ in range(reps)]
    )
    return float(vals.mean()), float(draws[int(0.025 * reps)]), float(draws[int(0.975 * reps) - 1])


def main() -> None:
    """Closure per assertion, per arm and network, and the paired contrasts."""
    ratio: dict[str, dict[str, float]] = {}
    for r in csv.DictReader(ROWS.open()):
        if r["status"] == "knowledge_inconsistent" or not r["k_g0"] or int(r["n_k"]) == 0:
            continue
        ratio.setdefault(r["arm"], {})[r["network"]] = int(r["k_g0"]) / int(r["n_k"])
    out: dict = {"per_arm": {}, "paired": {}}
    for arm, per in sorted(ratio.items()):
        vals = np.array(list(per.values()))
        out["per_arm"][arm] = {
            "networks": len(per),
            "mean": float(vals.mean()),
            "median": float(np.median(vals)),
            "max": float(vals.max()),
        }
    for a, b in (
        ("A_ORACLE", "D_LLM"),
        ("A_ORACLE", "D_LLM_72B"),
        ("A_ORACLE", "D_LLM_SRC_70B"),
        ("A_TRUE", "D_LLM"),
        ("D_LLM_72B", "D_LLM"),
    ):
        if a in ratio and b in ratio:
            shared = sorted(set(ratio[a]) & set(ratio[b]))
            if shared:
                m, lo, hi = boot({n: ratio[a][n] - ratio[b][n] for n in shared})
                larger = sum(1 for n in shared if ratio[a][n] > ratio[b][n])
                out["paired"][f"{a} - {b}"] = {
                    "mean_diff": m,
                    "ci": [lo, hi],
                    "networks": len(shared),
                    "a_larger_on": larger,
                    "verdict": "separated"
                    if lo > 0
                    else ("reversed" if hi < 0 else "not separated"),
                }
    OUT.write_text(json.dumps(out, indent=1))
    print(f"{'arm':16s} {'nets':>4} {'mean':>6} {'median':>6} {'max':>6}")
    for arm, s in out["per_arm"].items():
        print(f"{arm:16s} {s['networks']:4d} {s['mean']:6.2f} {s['median']:6.2f} {s['max']:6.2f}")
    print("\npaired, closure per assertion, difference of per-network values:")
    for k, v in out["paired"].items():
        ci = f"[{v['ci'][0]:+.3f}, {v['ci'][1]:+.3f}]"
        print(
            f"  {k:28s} {v['mean_diff']:+.3f} {ci}  larger on {v['a_larger_on']}/{v['networks']}"
            f"  -> {v['verdict']}"
        )
    if math.isnan(out["paired"].get("A_ORACLE - D_LLM", {}).get("mean_diff", math.nan)):
        print("A_ORACLE or D_LLM absent")


if __name__ == "__main__":
    main()
