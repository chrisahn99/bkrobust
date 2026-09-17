"""The monotonicity audit: failure at a retraction set implies failure at every superset.

The claim certificate is a minimum minimal correction set only if failure is
upward-closed under retraction: once retracting a set of the analyst's claims
makes the committed set invalid, retracting more claims cannot make it valid
again, because re-closure from fewer claims represents more DAGs. That is
Theorem 1 of the theory as implemented, and the protocol makes it a halting
check: a single counterexample stops everything. This samples failing subsets
on the ledger's rows and random supersets of each, and counts.

Writes ``results/ledger/monotonicity_audit.json`` and prints the count.

    python experiments/monotonicity_audit.py --pairs 5000
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments"))

from bkrobust.demo.example import knowledge_to_recover  # noqa: E402
from bkrobust.demo.graph import MPDAG  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.gac.mpdag_level import is_gac_valid_mpdag  # noqa: E402
from bkrobust.mpdag_criterion.optimal import committed_adjustment_set  # noqa: E402
from ledger_sweep import KNOWLEDGE, load  # noqa: E402

ROWS = ROOT / "results" / "ledger" / "rows.csv"
OUT = ROOT / "results" / "ledger" / "monotonicity_audit.json"


def failing_subsets(
    cpdag: MPDAG,
    k: list[tuple[str, str]],
    x: str,
    y: str,
    z: frozenset,
    max_depth: int = 3,
    budget: int = 200,
) -> list[frozenset]:
    """Retraction subsets of size at most ``max_depth`` that break validity, up to a budget."""
    found: list[frozenset] = []
    used = 0
    for depth in range(1, min(max_depth, len(k)) + 1):
        for subset in itertools.combinations(range(len(k)), depth):
            used += 1
            if used > budget:
                return found
            g = apply_orientations(cpdag, [e for i, e in enumerate(k) if i not in subset])
            if g is not None and not is_gac_valid_mpdag(g, x, y, z):
                found.append(frozenset(subset))
    return found


def main(argv: list[str] | None = None) -> None:
    """Sample failing subsets and their supersets, count violations."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=20260913)
    args = ap.parse_args(argv)
    rng = np.random.default_rng(args.seed)
    know = json.loads(KNOWLEDGE.read_text())
    rows = [
        r
        for r in csv.DictReader(ROWS.open())
        if r["status"] == "ok" and r["arm"] != "D_RAND" and int(r["n_k"]) >= 2
    ]
    rng.shuffle(rows)
    graphs: dict[str, tuple] = {}
    checked = violations = rows_used = 0
    examples: list[dict] = []
    t0 = time.perf_counter()
    for r in rows:
        if checked >= args.pairs:
            break
        net, arm = r["network"], r["arm"]
        if net not in graphs:
            graphs[net] = load(net)
        dag, cpdag, _ = graphs[net]
        if arm in know:
            k = [tuple(e) for e in know[arm]["networks"][net]["k"]]
        elif arm == "A_ORACLE":
            k = sorted(knowledge_to_recover(dag, cpdag))
        elif arm == "A_TRUE":
            base = [tuple(e) for e in know["D_LLM"]["networks"][net]["k"]]
            k = [(a, b) if dag.is_directed_edge(a, b) else (b, a) for (a, b) in base]
        else:
            continue
        g0 = apply_orientations(cpdag, k)
        if g0 is None:
            continue
        z = committed_adjustment_set(g0, r["x"], r["y"]).z
        if z is None:
            continue
        fails = failing_subsets(cpdag, k, r["x"], r["y"], z)
        if not fails:
            continue
        rows_used += 1
        for s in fails[:4]:
            rest = [i for i in range(len(k)) if i not in s]
            for _ in range(3):
                if checked >= args.pairs or not rest:
                    break
                extra = rng.choice(
                    rest, size=int(rng.integers(1, min(3, len(rest)) + 1)), replace=False
                )
                sup = s | set(int(i) for i in extra)
                g = apply_orientations(cpdag, [e for i, e in enumerate(k) if i not in sup])
                checked += 1
                if g is not None and is_gac_valid_mpdag(g, r["x"], r["y"], z):
                    violations += 1
                    examples.append(
                        {
                            "network": net,
                            "arm": arm,
                            "x": r["x"],
                            "y": r["y"],
                            "S": sorted(s),
                            "S_sup": sorted(sup),
                        }
                    )
    summary = {
        "pairs_checked": checked,
        "rows_used": rows_used,
        "violations": violations,
        "one_sided_95_upper_bound": round(1 - 0.05 ** (1 / checked), 6)
        if checked and not violations
        else None,
        "seconds": round(time.perf_counter() - t0, 1),
        "examples": examples[:20],
    }
    OUT.write_text(json.dumps(summary, indent=1))
    print(json.dumps({k: v for k, v in summary.items() if k != "examples"}))
    if violations:
        print("HALT: monotonicity violated; see examples in the json")


if __name__ == "__main__":
    main()
