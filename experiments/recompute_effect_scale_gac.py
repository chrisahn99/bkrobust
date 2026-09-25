"""Recompute r_val under the GAC for the effect-scale queries of results/final_table.

``results/final_table/instances.jsonl`` stores r_val computed by
``bkrobust.hybrid.breakdown_radius`` when it decided validity with the back-door
criterion. ``breakdown_radius`` now defaults to the generalised adjustment
criterion (GAC); this script rebuilds each solved instance exactly as
``experiments/final_table.py`` does (same networks, same LLM knowledge, Z the
optimal set of G0) and compares the stored radius with the GAC radius.

Run with::

    PYTHONPATH=src:experiments .venv/bin/python experiments/recompute_effect_scale_gac.py
"""

from __future__ import annotations

import json
import signal
import sys
import time
from pathlib import Path

sys.argv = sys.argv[:1]
import final_table as ft  # noqa: E402

from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.hybrid import breakdown_radius  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "results" / "final_table" / "instances.jsonl"
OUT = ROOT / "results" / "radius_gac_recompute" / "effect_scale_final_table.json"
TIME_LIMIT_S = 30


class _Timeout(Exception):
    pass


def _alarm(*_: object) -> None:
    raise _Timeout


def _same(stored: object, new: object) -> bool:
    try:
        return float(str(stored)) == float(new)
    except (TypeError, ValueError):
        return str(stored) == str(new)


def main() -> None:
    records = [json.loads(line) for line in SRC.open()]
    solved = [r for r in records if r["status"] == "ok" and r["r_val"] not in (None, "None")]
    knowledge = {c: ft.load_knowledge(c) for c in {r["condition"] for r in solved}}
    nets = ft.load_networks(names={r["network"] for r in solved}, skip=set(), max_nodes=0)
    signal.signal(signal.SIGALRM, _alarm)

    n_equal = n_differ = n_skipped = 0
    differences: list[dict[str, object]] = []
    started = time.time()
    for r in solved:
        net = nets.get(r["network"])
        if net is None:
            n_skipped += 1
            continue
        cpdag = net["cpdag"]
        k = knowledge[r["condition"]].get(r["network"], [])
        signal.alarm(TIME_LIMIT_S)
        try:
            g0 = apply_orientations(cpdag, k)
            z = frozenset(optimal_adjustment_set_mpdag(g0, r["x"], r["y"]))
            res = breakdown_radius(cpdag, None, r["x"], r["y"], z, g0=g0, criterion="gac")
            signal.alarm(0)
        except _Timeout:
            n_skipped += 1
            continue
        if _same(r["r_val"], res.radius):
            n_equal += 1
        else:
            n_differ += 1
            differences.append(
                {"network": r["network"], "x": r["x"], "y": r["y"],
                 "stored_backdoor": r["r_val"], "gac": res.radius, "method": res.method}
            )

    summary = {
        "source": str(SRC.relative_to(ROOT)),
        "n_units": len(solved),
        "n_equal": n_equal,
        "n_differ": n_differ,
        "n_skipped": n_skipped,
        "differences": differences,
        "seconds": round(time.time() - started, 1),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
