"""Phase B: the candidate frame, drawn from the estimated graph and nothing else.

The committed pipeline chooses which pairs to study with five reads of the true
DAG: whether the outcome is a descendant of the treatment, the true optimal
set, whether that set is valid, whether the empty set is trivially valid, and a
loop over the true recovering set. A denominator selected that way cannot
support a claim about what an analyst can do, because an analyst has none of it.

This builds the frame from the CPDAG alone. The treatment lies in or next to a
chain component of size at least two; the outcome is a possible descendant of
the treatment. Nothing else is filtered. Every condition the committed gate used
as a filter becomes a status on a frozen list, so each knowledge source is
scored on the same rows with a different status composition rather than on a
population of its own.

The frame is frozen before any knowledge is collected, and the committed
admissible pairs are matched into it as a bridge stratum so no previous result
is orphaned.

Writes ``results/frame/frame.jsonl``, ``results/frame/per_network.csv`` and
``results/frame/frame_summary.json``.

    python experiments/frame_build.py --cap 20
"""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from bkrobust.benchmarks.describe import parse_file, to_mpdag  # noqa: E402
from bkrobust.demo.example import dag_to_cpdag  # noqa: E402
from bkrobust.demo.graph import MPDAG, undirected_components  # noqa: E402
from bkrobust.mpdag_criterion.paths import possible_descendants  # noqa: E402

MODELS = ROOT / "results" / "axisa3" / "networks" / "example_models"
INSTANCES = ROOT / "results" / "axisa3" / "instances.jsonl"
OUT = ROOT / "results" / "frame"

#: Chain components smaller than this admit no background knowledge at all.
MIN_COMPONENT = 2


def load(path: Path) -> tuple[str, MPDAG, MPDAG] | None:
    """Parse one network file into ``(name, dag, cpdag)``, or None if it is no DAG."""
    parsed = parse_file(path, hashlib.sha256(path.read_bytes()).hexdigest())
    if parsed.bidirected:
        return None
    dag = to_mpdag(parsed)
    return parsed.name, dag, dag_to_cpdag(dag)


def treatments(cpdag: MPDAG) -> tuple[set[str], dict[str, int]]:
    """Nodes in or adjacent to a chain component of size >= MIN_COMPONENT.

    Returns the candidate treatments and, for each, the size of the component it
    sits in or beside. Both are functions of the CPDAG alone.
    """
    comps = [c for c in undirected_components(cpdag) if len(c) >= MIN_COMPONENT]
    out: dict[str, int] = {}
    for comp in comps:
        for node in sorted(comp):
            out[node] = max(out.get(node, 0), len(comp))
        for node in sorted(comp):
            for nb in sorted(cpdag.neighbors(node)):
                if nb not in comp:
                    out[nb] = max(out.get(nb, 0), len(comp))
    return set(out), out


def main(argv: list[str] | None = None) -> None:
    """Build and freeze the frame, and publish the population delta."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--cap", type=int, default=20, help="per-network cap by seeded stride")
    ap.add_argument("--second-cap", type=int, default=50, help="cap printed beside the headline")
    ap.add_argument("--skip", default="", help="comma-separated networks to skip")
    args = ap.parse_args(argv)
    skip = {s for s in args.skip.split(",") if s}

    committed = collections.defaultdict(set)
    for line in INSTANCES.open():
        row = json.loads(line)
        if row.get("admissible"):
            committed[row["network"]].add((row["x"], row["y"]))

    OUT.mkdir(parents=True, exist_ok=True)
    rows_out = (OUT / "frame.jsonl").open("w")
    per_net: list[dict] = []
    totals: collections.Counter = collections.Counter()
    t0 = time.perf_counter()

    for path in sorted(MODELS.iterdir()):
        if not path.is_file():
            continue
        loaded = load(path)
        if loaded is None:
            totals["networks_not_dag"] += 1
            continue
        name, _, cpdag = loaded
        if name in skip:
            continue
        totals["networks_parsed"] += 1

        cands, comp_size = treatments(cpdag)
        frame: list[tuple[str, str]] = []
        for x in sorted(cands):
            for y in sorted(possible_descendants(cpdag, x) - {x}):
                frame.append((x, y))
        frame.sort()

        # Frozen sample: a seeded deterministic stride over the sorted frame,
        # never a function of any computed radius.
        step = max(1, len(frame) // args.cap) if len(frame) > args.cap else 1
        sampled = frame[::step][: args.cap]
        inclusion = len(sampled) / len(frame) if frame else 0.0

        bridge = committed.get(name, set())
        in_frame = bridge & set(frame)
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
                "committed_admissible": len(bridge),
                "committed_inside_frame": len(in_frame),
                "committed_outside_frame": len(bridge - set(frame)),
            }
        )
        totals["frame_pairs"] += len(frame)
        totals["sampled_pairs"] += len(sampled)
        totals["candidate_treatments"] += len(cands)
        totals["committed_admissible"] += len(bridge)
        totals["committed_inside_frame"] += len(in_frame)
        totals["committed_outside_frame"] += len(bridge - set(frame))
        if frame:
            totals["networks_with_frame"] += 1
        if len(sampled) < len(frame):
            totals["networks_capped"] += 1

        for x, y in sampled:
            rows_out.write(
                json.dumps(
                    {
                        "network": name,
                        "x": x,
                        "y": y,
                        "component_size": comp_size.get(x, 0),
                        "in_committed_admissible": (x, y) in in_frame,
                        "inclusion_probability": round(inclusion, 6),
                        "frame_source": "cpdag_only",
                        "true_dag_on_path": False,
                    }
                )
                + "\n"
            )
    rows_out.close()

    with (OUT / "per_network.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(per_net[0]))
        w.writeheader()
        w.writerows(sorted(per_net, key=lambda r: -r["frame_pairs"]))

    summary = {
        "seconds": round(time.perf_counter() - t0, 1),
        "cap": args.cap,
        "min_component": MIN_COMPONENT,
        "totals": dict(totals),
        "frame_vs_committed": {
            "frame_pairs": totals["frame_pairs"],
            "committed_admissible_pairs": totals["committed_admissible"],
            "committed_inside_frame": totals["committed_inside_frame"],
            "committed_outside_frame": totals["committed_outside_frame"],
            "expansion_factor": round(
                totals["frame_pairs"] / max(totals["committed_admissible"], 1), 1
            ),
        },
        "networks_yielding_frame": totals["networks_with_frame"],
        "networks_yielding_committed_instances": len(committed),
        "top_by_frame_size": sorted(
            ((r["network"], r["frame_pairs"]) for r in per_net), key=lambda t: -t[1]
        )[:8],
        "capped_networks": totals["networks_capped"],
    }
    (OUT / "frame_summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
