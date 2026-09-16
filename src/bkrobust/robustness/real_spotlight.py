"""The discordance spotlight, and the draw-count verification it must survive.

Session 8 published a spotlight whose endpoint rested on depths with ``n_eval``
of 35, 10, 1 and 1, and withdrew it on verification. The pre-registration for
this campaign therefore requires that **every** grid point contributing to a
spotlight carry ``n_eval >= 30``, checked before the case is written up rather
than after. This module performs that check and writes the result, so the
verification travels with the claim.

The case: ``paths``, the flip arm, base wrongness 0.00.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import os
from pathlib import Path
from typing import Any

FIELDS = [
    "network", "x", "y", "n_k", "k_g0", "leverage", "shd_truth", "radius",
    "separation", "separation_status", "z_size", "grid_point", "n_draws",
    "n_eval", "n_contradictory", "S", "assumes",
]


def load_case(out_dir: Path, network: str) -> list[dict[str, Any]]:
    """Every flip-arm, base-wrongness-zero row of one network, with its cells.

    Args:
        out_dir: The results subtree.
        network: The network to spotlight.

    Returns:
        One row per (instance, grid point).
    """
    done = {os.path.basename(f)[:-5] for f in glob.glob(str(out_dir / "_done" / "*.json"))}
    insts: dict[tuple[str, str], dict[str, Any]] = {}
    for f in glob.glob(str(out_dir / "shards" / f"flip__{network}__*.instances.jsonl")):
        if os.path.basename(f)[: -len(".instances.jsonl")] not in done:
            continue
        for line in open(f):
            r = json.loads(line)
            if r["base_wrongness"] == 0.0 and r["status"] == "ok":
                insts[(r["x"], r["y"])] = r
    rows: list[dict[str, Any]] = []
    for f in glob.glob(str(out_dir / "shards" / f"flip__{network}__*.cells.jsonl")):
        if os.path.basename(f)[: -len(".cells.jsonl")] not in done:
            continue
        for line in open(f):
            c = json.loads(line)
            if c["base_wrongness"] != 0.0 or c["status"] != "ok":
                continue
            i = insts.get((c["x"], c["y"]))
            if i is None:
                continue
            rows.append({
                "network": network, "x": c["x"], "y": c["y"], "n_k": i["n_k"],
                "k_g0": i["k_g0"], "leverage": i["k_g0"] / i["n_k"] if i["n_k"] else None,
                "shd_truth": i["shd_truth"], "radius": i["radius"],
                "separation": i["separation"], "separation_status": i["separation_status"],
                "z_size": i["z_size"], "grid_point": c["grid_point"],
                "n_draws": c["n_draws"], "n_eval": c["n_eval"],
                "n_contradictory": c["n_contradictory"], "S": c["S"],
                "assumes": i["assumes"],
            })
    return sorted(rows, key=lambda r: (-(r["radius"] or 0), r["x"], r["y"], r["grid_point"]))


def main(argv: list[str] | None = None) -> int:
    """Entry point.

    Args:
        argv: Argument vector, or ``None`` for ``sys.argv``.

    Returns:
        Process exit status.
    """
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-dir", default="results/axis_robustness_real")
    p.add_argument("--network", default="paths")
    p.add_argument("--min-n-eval", type=int, default=30)
    args = p.parse_args(argv)

    out_dir = Path(args.out_dir)
    rows = load_case(out_dir, args.network)
    path = out_dir / f"spotlight_{args.network}.csv"
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    min_eval = min((r["n_eval"] for r in rows), default=0)
    radii = sorted({r["radius"] for r in rows})
    survivals = sorted({r["S"] for r in rows})
    summary = {
        "network": args.network,
        "n_instances": len({(r["x"], r["y"]) for r in rows}),
        "n_grid_point_rows": len(rows),
        "min_n_eval_over_all_contributing_grid_points": min_eval,
        "min_n_eval_required": args.min_n_eval,
        "spotlight_survives_draw_count_check": min_eval >= args.min_n_eval,
        "n_k": rows[0]["n_k"] if rows else None,
        "k_g0": rows[0]["k_g0"] if rows else None,
        "leverage": rows[0]["leverage"] if rows else None,
        "shd_truth": rows[0]["shd_truth"] if rows else None,
        "radius_values": radii,
        "survival_values": survivals,
        "total_contradictory_draws": sum(r["n_contradictory"] for r in rows),
        "assumes": rows[0]["assumes"] if rows else None,
    }
    (out_dir / f"spotlight_{args.network}.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
