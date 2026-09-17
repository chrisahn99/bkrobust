"""Phase B: the applied papers, each on the effect its own authors set out to estimate.

Every other query in this project is one we sampled. The dagitty files carry
their authors' own declaration of the exposure and the outcome, and that is the
only query in the corpus a domain expert chose. Until now the parser captured
those annotations in a regular expression group and dropped them, so the
pre-registration's domain-documented stratum had no rows.

For each declared pair this reports whether it is in the observables-only frame,
whether it is certifiable, what the radius is in claims and in hops, and the
first-failing retraction set in the authors' own variable names.

Writes ``results/frame/declared_rows.csv`` and ``results/frame/declared_summary.json``.

    python experiments/declared_queries.py
"""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import itertools
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from bkrobust.benchmarks.describe import parse_file, to_mpdag  # noqa: E402
from bkrobust.benchmarks.measure import separation  # noqa: E402
from bkrobust.demo.evaluate import (  # noqa: E402
    is_valid_adjustment_set_dag,
    optimal_adjustment_set_mpdag,
)
from bkrobust.demo.example import dag_to_cpdag, knowledge_to_recover  # noqa: E402
from bkrobust.demo.graph import MPDAG, undirected_components  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.gac.mpdag_level import is_gac_valid_mpdag  # noqa: E402
from bkrobust.hybrid import breakdown_radius  # noqa: E402
from bkrobust.mpdag_criterion.criterion import is_amenable  # noqa: E402
from bkrobust.mpdag_criterion.paths import possible_descendants  # noqa: E402

MODELS = ROOT / "results" / "axisa3" / "networks" / "example_models"
OUT = ROOT / "results" / "frame"
MAX_DEPTH = 3


def first_breaking(cpdag: MPDAG, k: list, x: str, y: str, z: frozenset) -> tuple[int | None, tuple]:
    """Smallest retraction that invalidates ``z``, and the first such subset."""
    for depth in range(1, MAX_DEPTH + 1):
        if depth > len(k):
            break
        for subset in itertools.combinations(range(len(k)), depth):
            drop = set(subset)
            g = apply_orientations(cpdag, [e for i, e in enumerate(k) if i not in drop])
            if g is None:
                continue
            if not is_gac_valid_mpdag(g, x, y, z):
                return depth, tuple(k[i] for i in subset)
    return None, ()


def main(argv: list[str] | None = None) -> None:
    """Evaluate every author-declared exposure-outcome pair in the corpus."""
    ap = argparse.ArgumentParser()
    ap.parse_args(argv)

    rows: list[dict] = []
    for path in sorted(MODELS.glob("*.txt")):
        parsed = parse_file(path, hashlib.sha256(path.read_bytes()).hexdigest())
        if not parsed.roles:
            continue
        exposures = [n for n, r in parsed.roles.items() if "exposure" in r]
        outcomes = [n for n, r in parsed.roles.items() if "outcome" in r]
        adjusted = sorted(n for n, r in parsed.roles.items() if "adjusted" in r)
        latent = sorted(n for n, r in parsed.roles.items() if "latent" in r)
        selected = sorted(n for n, r in parsed.roles.items() if "selected" in r)
        names = dict(parsed.name_map)

        record: dict = {
            "network": parsed.name,
            "is_dag": not parsed.bidirected,
            "exposure": names.get(exposures[0], exposures[0]) if exposures else "",
            "outcome": names.get(outcomes[0], outcomes[0]) if outcomes else "",
            "author_adjusted": " ".join(names.get(a, a) for a in adjusted),
            "author_latent": len(latent),
            "author_selected": " ".join(names.get(a, a) for a in selected),
            "status": "",
            "in_frame": "",
            "n_claims": "",
            "r_claim": "",
            "r_hop": "",
            "separation": "",
            "retraction_set": "",
            "z_size": "",
            "z_valid_at_truth": "",
        }
        if parsed.bidirected:
            record["status"] = "not_a_dag"
            rows.append(record)
            continue
        if not exposures or not outcomes:
            record["status"] = "no_declared_pair"
            rows.append(record)
            continue

        dag = to_mpdag(parsed)
        cpdag = dag_to_cpdag(dag)
        x, y = exposures[0], outcomes[0]

        comps = [c for c in undirected_components(cpdag) if len(c) >= 2]
        near = any(x in c or any(cpdag.has_edge(x, v) for v in sorted(c)) for c in comps)
        record["in_frame"] = int(near and y in possible_descendants(cpdag, x) - {x})
        if not record["in_frame"]:
            record["status"] = (
                "treatment_has_no_free_edge" if not near else "outcome_not_possibly_descendant"
            )
            rows.append(record)
            continue

        k = sorted(knowledge_to_recover(dag, cpdag))
        record["n_claims"] = len(k)
        g0 = apply_orientations(cpdag, k)
        if g0 is None:
            record["status"] = "knowledge_inconsistent"
            rows.append(record)
            continue
        if not is_amenable(g0, x, y):
            record["status"] = "not_amenable"
            rows.append(record)
            continue
        o = optimal_adjustment_set_mpdag(g0, x, y)
        if o is None:
            record["status"] = "optimal_not_identified"
            rows.append(record)
            continue
        if not o:
            record["status"] = "optimal_empty"
            rows.append(record)
            continue
        z = frozenset(o)
        if not is_gac_valid_mpdag(g0, x, y, z):
            record["status"] = "z_invalid_at_g0"
            rows.append(record)
            continue

        sep, sep_status = separation(cpdag, x, z)
        hop = breakdown_radius(cpdag, None, x, y, z, g0=g0, time_limit_s=30.0)
        depth, subset = first_breaking(cpdag, k, x, y, z)
        record.update(
            {
                "status": "ok",
                "r_claim": depth if depth is not None else f">{MAX_DEPTH}",
                "r_hop": "unreached" if hop.radius == -1 else hop.radius,
                "separation": sep if sep_status == "measured" else "undefined",
                "retraction_set": (
                    " and ".join(f"{names.get(a, a)} -> {names.get(b, b)}" for a, b in subset)
                    if subset
                    else "nothing inside the depth budget"
                ),
                "z_size": len(z),
                "z_valid_at_truth": int(is_valid_adjustment_set_dag(dag, x, y, z)),
            }
        )
        rows.append(record)

    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "declared_rows.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    summary = {
        "files_with_declared_roles": len(rows),
        "status": dict(collections.Counter(r["status"] for r in rows).most_common()),
        "in_frame": sum(1 for r in rows if r["in_frame"] == 1),
        "certifiable": sum(1 for r in rows if r["status"] == "ok"),
        "with_author_adjustment_set": sum(1 for r in rows if r["author_adjusted"]),
        "with_selection_node": sum(1 for r in rows if r["author_selected"]),
        "with_latent_nodes": sum(1 for r in rows if r["author_latent"]),
        "certifiable_and_nothing_breaks": sum(
            1 for r in rows if r["status"] == "ok" and r["r_claim"] == f">{MAX_DEPTH}"
        ),
        "certifiable_and_breaks_inside_budget": sum(
            1 for r in rows if r["status"] == "ok" and r["r_claim"] != f">{MAX_DEPTH}"
        ),
        "selection_note": (
            "The committed pipeline rejects a pair when no single retraction of the recovering "
            "set changes validity, which is the status "
            "`no_atomic_perturbation_changes_validity`. The author-declared queries fall "
            "overwhelmingly on that side, so the population every radius in this project was "
            "measured on is close to the complement of the queries these papers were written to "
            "answer."
        ),
        "rows": rows,
    }
    (OUT / "declared_summary.json").write_text(json.dumps(summary, indent=1))

    print(f"{len(rows)} files declare a query; {summary['certifiable']} are certifiable\n")
    header = f"{'network':18s} {'exposure':22s} {'outcome':14s} {'status':22s}"
    print(f"{header} {'r_claim':>7} {'r_hop':>5}")
    for r in rows:
        print(
            f"{r['network']:18s} {str(r['exposure'])[:22]:22s} {str(r['outcome'])[:14]:14s} "
            f"{r['status']:22s} {r['r_claim']!s:>7} {r['r_hop']!s:>5}"
        )
    breaks = [r for r in rows if r["status"] == "ok" and r["r_claim"] != f">{MAX_DEPTH}"]
    robust = [r for r in rows if r["status"] == "ok" and r["r_claim"] == f">{MAX_DEPTH}"]
    print(f"\nnothing inside a depth-{MAX_DEPTH} retraction breaks the set on {len(robust)} of")
    print(f"the {summary['certifiable']} certifiable declared queries:")
    for r in robust:
        print(f"  {r['network']}: {r['exposure']} on {r['outcome']}, {r['n_claims']} claims")
    print("\nand where something does break it:")
    for r in breaks:
        print(
            f"  {r['network']}: the effect of {r['exposure']} on {r['outcome']} holds "
            f"unless {r['retraction_set']} is retracted"
        )


if __name__ == "__main__":
    main()
