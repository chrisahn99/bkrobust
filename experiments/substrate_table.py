"""Table 2 substrate: every network of the corpus with its structure, tier and sortability.

One row per file under ``example_models``. Structure comes from the frozen
descriptive sweep, ``results/axisa3/descriptive_structure.json``, and is
cross-checked against ``results/frame/per_network.csv``; the two sortabilities
come from ``results/substrate/sortability.json``; the ledger row count from
``results/ledger/rows.csv``; the tier from ``ledger_sweep.TIER``. A fixture the
descriptive sweep does not know, meaning one added after the acquisition, is
described at run time and flagged as outside the frozen frame.

Mean degree is ``2 * edges / nodes``, the skeleton's. Undirected fraction is the
share of edges the CPDAG leaves undirected, the budget background knowledge has
to work with. The chain-component column lists component sizes (at least two
nodes) in descending order with multiplicities.

Writes ``results/substrate/TABLE2_SUBSTRATE.md`` and ``table2_substrate.csv``,
with the illustration tier T0 in its own table above a rule.

    python experiments/substrate_table.py
"""

from __future__ import annotations

import collections
import csv
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
for entry in (ROOT, ROOT / "src"):
    if str(entry) not in sys.path:
        sys.path.insert(0, str(entry))

from bkrobust.benchmarks.describe import describe_network, parse_file  # noqa: E402
from experiments.ledger_sweep import TIER  # noqa: E402

MODELS = ROOT / "results" / "axisa3" / "networks" / "example_models"
DESCRIPTIVE = ROOT / "results" / "axisa3" / "descriptive_structure.json"
PER_NETWORK = ROOT / "results" / "frame" / "per_network.csv"
SORTABILITY = ROOT / "results" / "substrate" / "sortability.json"
ROWS = ROOT / "results" / "ledger" / "rows.csv"
OUT_DIR = ROOT / "results" / "substrate"

#: Networks the sweep excluded for cost (PROTOCOL_RUN.md, Stage 5).
SWEEP_NODE_CAP = 500
#: The questionnaire's per-component cap (PROTOCOL_RUN.md, departures).
QUESTIONNAIRE_COMPONENT_CAP = 20

COLUMNS = (
    "network",
    "tier",
    "nodes",
    "edges",
    "mean_degree",
    "undirected_fraction",
    "chain_components",
    "var_sortability",
    "r2_sortability",
    "ledger_rows",
    "note",
)


def component_string(sizes: list[int]) -> str:
    """Sizes in descending order with multiplicities, e.g. ``19, 14, 5, 2x3``."""
    if not sizes:
        return "none"
    counts = collections.Counter(sizes)
    return ", ".join(
        f"{size}x{count}" if count > 1 else str(size)
        for size, count in sorted(counts.items(), reverse=True)
    )


def sortability_cell(entry: dict[str, object] | None, key: str) -> str:
    """``mean +- sd`` over seeds, or the reason it is undefined."""
    if entry is None:
        return "UNDEFINED (no probe entry)"
    if entry.get("status") != "defined":
        return "UNDEFINED"
    stats = entry[key]
    assert isinstance(stats, dict)
    return f"{stats['mean']:.3f} +- {stats['sd']:.3f}"


def network_key(file_name: str) -> str:
    """The network name a file is known by, ``alarm`` for ``alarm.bif.gz``."""
    return file_name[:-7] if file_name.endswith(".bif.gz") else Path(file_name).stem


def describe_unknown(path: Path) -> dict[str, object]:
    """Describe a fixture the frozen sweep does not list, in the sweep's own vocabulary."""
    parsed = parse_file(path, hashlib.sha256(path.read_bytes()).hexdigest())
    row = describe_network(parsed)
    return {
        "network": row.network,
        "n_nodes": row.n_nodes,
        "n_edges": row.n_edges,
        "n_undirected_edges": row.n_undirected_edges,
        "undirected_fraction": row.undirected_fraction,
        "component_sizes": list(row.component_sizes),
        "status": row.status,
        "note": row.note,
    }


def note_for(desc: dict[str, object], ledger_rows: int, frozen: bool) -> str:
    """The one-line remark that says why a network is where it is."""
    notes: list[str] = []
    if not frozen:
        notes.append("added 2026-09-12; outside the frozen frame and the ledger")
    if desc["status"] == "not_a_dag":
        notes.append(str(desc["note"]))
    elif int(desc["n_nodes"]) > SWEEP_NODE_CAP:
        notes.append(f"above {SWEEP_NODE_CAP} nodes; excluded from the sweep")
    elif not desc["component_sizes"]:
        notes.append("no chain component; no frame pair")
    elif frozen and ledger_rows == 0:
        sizes = [int(s) for s in desc["component_sizes"]]
        if min(sizes) > QUESTIONNAIRE_COMPONENT_CAP:
            notes.append(
                f"every component above the questionnaire cap of {QUESTIONNAIRE_COMPONENT_CAP}"
            )
        else:
            notes.append("no ledger row")
    name = str(desc["network"])
    if name == "Didelez_2010":
        notes.append("S is a selection node (dagitty 'selected'); see DECISIONS.md")
    if name == "sachs":
        notes.append(
            "17 arcs identical to the Sachs et al. 2005 inferred network; see DECISIONS.md"
        )
    if name == "M-bias":
        notes.append("expanded as M-structure; see DECISIONS.md")
    return "; ".join(notes)


def main() -> None:
    """Assemble the table and write both renderings."""
    descriptive = {
        r["network"]: r for r in json.loads(DESCRIPTIVE.read_text(encoding="utf-8"))["networks"]
    }
    sortability = json.loads(SORTABILITY.read_text(encoding="utf-8"))["networks"]
    ledger_rows: collections.Counter[str] = collections.Counter()
    with ROWS.open(encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            ledger_rows[row["network"]] += 1

    mismatches: list[str] = []
    with PER_NETWORK.open(encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            desc = descriptive.get(row["network"])
            if desc is None:
                mismatches.append(f"{row['network']}: in per_network.csv, not in the sweep")
                continue
            if (
                int(row["nodes"]) != desc["n_nodes"]
                or int(row["undirected_edges"]) != (desc["n_undirected_edges"])
            ):
                mismatches.append(f"{row['network']}: per_network.csv disagrees with the sweep")

    table: list[dict[str, object]] = []
    for path in sorted(MODELS.iterdir()):
        if not path.is_file() or path.name.startswith("_"):
            continue
        name = network_key(path.name)
        frozen = name in descriptive
        desc = descriptive[name] if frozen else describe_unknown(path)
        probe = sortability.get(name)
        if probe is not None and probe.get("status") == "defined":
            if probe["n_nodes"] != desc["n_nodes"] or probe["n_edges"] != desc["n_edges"]:
                mismatches.append(f"{name}: sortability.json disagrees with the sweep")
        nodes = int(desc["n_nodes"])
        edges = int(desc["n_edges"])
        fraction = desc["undirected_fraction"]
        table.append(
            {
                "network": name,
                "tier": TIER.get(name, "T1"),
                "nodes": nodes,
                "edges": edges,
                "mean_degree": round(2 * edges / nodes, 2) if nodes else 0.0,
                "undirected_fraction": "n/a" if fraction is None else f"{fraction:.3f}",
                "chain_components": (
                    "n/a"
                    if desc["status"] == "not_a_dag"
                    else component_string([int(s) for s in desc["component_sizes"]])
                ),
                "var_sortability": sortability_cell(probe, "var_sortability"),
                "r2_sortability": sortability_cell(probe, "r2_sortability"),
                "ledger_rows": ledger_rows.get(name, 0),
                "note": note_for(desc, ledger_rows.get(name, 0), frozen),
            }
        )

    order = {"T0": 0, "T1": 1, "T2": 2, "T3": 3}
    table.sort(key=lambda r: (order[str(r["tier"])], str(r["network"]).lower()))
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    with (OUT_DIR / "table2_substrate.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(COLUMNS))
        writer.writeheader()
        writer.writerows(table)

    header = (
        "| network | tier | nodes | edges | mean deg. | undirected | chain components | "
        "var-sort. | R2-sort. | ledger rows | note |\n"
        "|---|---|---:|---:|---:|---:|---|---|---|---:|---|\n"
    )

    def render(rows: list[dict[str, object]]) -> str:
        return header + "".join(
            "| " + " | ".join(str(r[c]) for c in COLUMNS) + " |\n" for r in rows
        )

    t0 = [r for r in table if r["tier"] == "T0"]
    rest = [r for r in table if r["tier"] != "T0"]
    probe_meta = json.loads(SORTABILITY.read_text(encoding="utf-8"))
    lines = [
        "# Table 2 substrate: the corpus, per network\n",
        "Generated by `experiments/substrate_table.py`. Structure from the frozen descriptive",
        "sweep (`results/axisa3/descriptive_structure.json`), cross-checked against",
        "`results/frame/per_network.csv`; sortability from `results/substrate/sortability.json`",
        f"(n = {probe_meta['n']}, {len(probe_meta['seeds'])} seeds, mean +- sd over seeds);",
        "ledger rows from `results/ledger/rows.csv`; tier from `ledger_sweep.TIER`.",
        "Mean degree is `2 * edges / nodes`. Undirected is the fraction of edges the CPDAG",
        "leaves undirected. Chain components are the sizes of the undirected components of at",
        "least two nodes, descending, with multiplicities. Sortability is UNDEFINED wherever",
        "the file carries no continuous parameters. Decisions on the fixtures are in",
        "`DECISIONS.md`.\n",
        "## T0, illustration\n",
        render(t0),
        "\n---\n",
        "## T1 (BIF networks), T2 (applied dagitty DAGs), T3 (Gaussian networks)\n",
        render(rest),
    ]
    if mismatches:
        lines += ["\n## Cross-check\n"] + [f"- {m}" for m in mismatches] + [""]
    else:
        lines += [
            "\n## Cross-check\n",
            "`per_network.csv` and `sortability.json` agree with the descriptive sweep on every",
            "node and edge count they share.\n",
        ]
    (OUT_DIR / "TABLE2_SUBSTRATE.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"{len(table)} rows ({len(t0)} in T0); mismatches: {len(mismatches)}")
    for m in mismatches:
        print("  " + m)


if __name__ == "__main__":
    main()
