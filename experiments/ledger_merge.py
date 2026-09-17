"""Join the per-process parts of a split ledger sweep back into one ledger.

``ledger_sweep.py --networks ... --out results/ledger/parts/<name>`` writes one
``rows.csv`` and one ``sweep_summary.json`` per part. This concatenates the
rows in network order, sums the status counts, takes the longest part as the
wall time, unions the size exclusions, and checks that no network appears in
two parts and that every part ran under the same settings.

    python experiments/ledger_merge.py results/ledger/parts/*
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "ledger"


def main(parts: list[str]) -> None:
    """Merge the parts into ``results/ledger/{rows.csv,sweep_summary.json}``."""
    rows: list[dict] = []
    fields: list[str] | None = None
    seen: dict[str, str] = {}
    status: dict[str, dict[str, int]] = {}
    settings = None
    arms = None
    excluded: dict[str, dict] = {}
    seconds = 0.0
    for part in sorted(parts):
        d = Path(part)
        with (d / "rows.csv").open(newline="") as fh:
            reader = csv.DictReader(fh)
            fields = fields or list(reader.fieldnames or [])
            if list(reader.fieldnames or []) != fields:
                raise SystemExit(f"{d}: columns differ from the first part")
            for r in reader:
                if seen.setdefault(r["network"], d.name) != d.name:
                    raise SystemExit(f"{r['network']} appears in {seen[r['network']]} and {d.name}")
                rows.append(r)
        s = json.loads((d / "sweep_summary.json").read_text())
        if settings is None:
            settings, arms = s["settings"], s["arms"]
        elif s["settings"] != settings or s["arms"] != arms:
            raise SystemExit(f"{d}: settings or arms differ from the first part")
        for arm, c in s["status_by_arm"].items():
            for k, v in c.items():
                status.setdefault(arm, {})[k] = status.get(arm, {}).get(k, 0) + v
        for e in s.get("excluded_for_size", []):
            excluded[e["network"]] = e
        seconds = max(seconds, s["seconds"])
    rows.sort(key=lambda r: r["network"])  # stable: keeps arm and row order within a network
    assert fields is not None
    with (OUT / "rows.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    networks = sorted(seen)
    frame_rows = len({(r["network"], r["x"], r["y"]) for r in rows})
    summary = {
        "frame_rows": frame_rows,
        "networks": len(networks),
        "arms": arms,
        "evaluations": len(rows),
        "seconds": round(seconds, 1),
        "parts": [Path(p).name for p in sorted(parts)],
        "status_by_arm": {a: status[a] for a in sorted(status)},
        "settings": settings,
        "excluded_for_size": sorted(excluded.values(), key=lambda e: e["network"]),
    }
    (OUT / "sweep_summary.json").write_text(json.dumps(summary, indent=1))
    print(f"{len(rows)} rows over {len(networks)} networks from {len(parts)} parts")
    print(json.dumps(summary["status_by_arm"], indent=1))


if __name__ == "__main__":
    main(sys.argv[1:])
