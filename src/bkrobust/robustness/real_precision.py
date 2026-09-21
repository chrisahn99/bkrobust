"""How much of each cell's corruption space the sweep actually covered.

A draw count is not the same thing as an amount of evidence. In the flip arm the
corruption space at depth ``d`` from ``|K|`` claims has exactly ``C(|K|, d)``
elements, and on this corpus ``|K|`` is small on most networks: at ``|K| = 1``
there is **one** possible corruption, so a thousand draws are a thousand copies
of the same state. That is not a weakness -- it makes the cell **exhaustive**,
and an exhaustive cell carries *no Monte-Carlo error at all* rather than a small
one -- but it must be reported as what it is, because "1000 draws" reads as a
thousand independent pieces of evidence and at ``|K| = 1`` it is not.

The count of distinct states actually computed is recorded for free by the
closure cache in :class:`bkrobust.robustness.real_survival.ClosureCache`, whose
``graph_misses`` is exactly the number of distinct corrupted claim sets a shard
saw. This module compares that against the combinatorial size of the space.

Nothing here changes a statistic. It describes the evidence behind one.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import os
from pathlib import Path
from typing import Any

FIELDS = [
    "shard_id", "arm", "network", "coverage", "base_wrongness", "bw_abs",
    "n_k", "n_draws_per_grid_point", "n_grid_points", "space_size",
    "n_distinct_states_computed", "n_draws_total", "coverage_of_space",
    "is_exhaustive",
]


def flip_space_size(n_k: int, depths: list[int]) -> int | None:
    """The number of distinct corrupted claim sets the flip arm can produce.

    Args:
        n_k: The stated knowledge size.
        depths: The depth grid actually swept.

    Returns:
        ``sum(C(n_k, d))`` over the swept depths, or ``None`` if ``n_k`` is 0.
    """
    if n_k <= 0:
        return None
    return sum(math.comb(n_k, d) for d in depths if 0 <= d <= n_k)


def main(argv: list[str] | None = None) -> int:
    """Entry point.

    Args:
        argv: Argument vector, or ``None`` for ``sys.argv``.

    Returns:
        Process exit status.
    """
    from bkrobust.robustness.real_survival import depth_grid

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-dir", default="results/axis_robustness_real")
    args = p.parse_args(argv)
    out_dir = Path(args.out_dir)

    rows: list[dict[str, Any]] = []
    for f in sorted(glob.glob(str(out_dir / "_done" / "flip__*.json"))):
        m = json.loads(Path(f).read_text())
        if m.get("g0_status") != "ok" or not m.get("n_cells"):
            continue
        n_k = int(m["n_k"])
        depths = list(depth_grid(n_k)) if n_k >= 1 else []
        space = flip_space_size(n_k, depths)
        distinct = int(m["cache"]["graph_misses"])
        n_draws = int(m["n_draws"])
        total = n_draws * len(depths)
        spec = m["spec"]
        rows.append({
            "shard_id": m["shard_id"], "arm": "flip", "network": spec["network"],
            "coverage": spec["coverage"], "base_wrongness": spec.get("base_wrongness"),
            "bw_abs": spec.get("bw_abs"), "n_k": n_k,
            "n_draws_per_grid_point": n_draws, "n_grid_points": len(depths),
            "space_size": space, "n_distinct_states_computed": distinct,
            "n_draws_total": total,
            "coverage_of_space": (distinct / space) if space else None,
            "is_exhaustive": bool(space is not None and distinct >= space),
        })

    path = out_dir / "precision_space_coverage.csv"
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    exhaustive = [r for r in rows if r["is_exhaustive"]]
    sampled = [r for r in rows if not r["is_exhaustive"]]
    summary = {
        "n_flip_shards_scored": len(rows),
        "n_exhaustive_shards": len(exhaustive),
        "n_sampled_shards": len(sampled),
        "exhaustive_networks": sorted({r["network"] for r in exhaustive}),
        "sampled_networks": sorted({r["network"] for r in sampled}),
        "max_n_k_that_is_exhaustive": max((r["n_k"] for r in exhaustive), default=None),
        "min_n_k_that_is_sampled": min((r["n_k"] for r in sampled), default=None),
        "min_space_coverage_among_sampled": min(
            (r["coverage_of_space"] for r in sampled if r["coverage_of_space"] is not None),
            default=None,
        ),
        "note": (
            "An exhaustive shard enumerated its entire corruption space, so its "
            "survival numbers carry no Monte-Carlo error at all. A sampled shard "
            "drew n_draws per grid point from a larger space."
        ),
    }
    (out_dir / "precision_space_coverage.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
