"""Table 8, the decision table, on the elicited rows.

For every certifiable query under the primary supplier whose certificate is
exact, the first-failing retraction set named in the analyst's own variables,
the number of claims it takes, what the free hop-count rule said, and, in the
audit block that reads the truth, whether that claim was in fact false and how
the committed set fares at the truth. This is the column no count baseline can
produce: a sentence of the form "the estimate holds unless this claim is
retracted", in domain terms.

Writes ``results/ledger/TABLE8_DECISION.md`` and ``decision_table.csv``.

    python experiments/decision_table.py --arm D_LLM
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments"))

from ledger_sweep import KNOWLEDGE, TIER, load  # noqa: E402

LEDGER = ROOT / "results" / "ledger"


def main(argv: list[str] | None = None) -> None:
    """Assemble the decision table for one arm."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="D_LLM")
    ap.add_argument("--max-rows", type=int, default=60, help="rows printed in the markdown")
    args = ap.parse_args(argv)
    know = json.loads(KNOWLEDGE.read_text())
    rows = [
        r
        for r in csv.DictReader((LEDGER / "rows.csv").open())
        if r["arm"] == args.arm and r["status"] == "ok"
    ]
    names_of: dict[str, dict] = {}
    dags: dict[str, object] = {}
    out_rows = []
    for r in rows:
        net = r["network"]
        if net not in names_of:
            dag, _, names = load(net)
            names_of[net], dags[net] = names, dag
        names = names_of[net]
        dag = dags[net]
        k = {tuple(e) for e in know[args.arm]["networks"][net]["k"]}
        exact = r["r_claim_status"] == "exact"
        witness = r["witness"]
        # was the named claim false at the truth? (audit column)
        witness_false = None
        if exact and witness:
            inv = {f"{names.get(a, a)}->{names.get(b, b)}": (a, b) for a, b in k}
            parts = [p.strip() for p in witness.split(" and ")]
            claims = [inv[p] for p in parts if p in inv]
            if claims:
                witness_false = any(not dag.is_directed_edge(a, b) for a, b in claims)
        out_rows.append(
            {
                "network": net,
                "tier": TIER.get(net, "T1"),
                "query": f"{names.get(r['x'], r['x'])} -> {names.get(r['y'], r['y'])}",
                "r_claim": r["r_claim"] if exact else r["r_claim_status"],
                "safe_moves": (int(r["r_claim"]) - 1) if exact else "",
                "first_failing_retraction": witness if exact else "",
                "free_rule": r["separation"]
                if r["separation_status"] == "measured"
                else "undefined",
                "r_hop": r["r_hop"] if r["hop_exact"] == "1" else f">= {r['hop_lower_bound']}",
                "claims_asserted": r["n_k"],
                "audit_named_claim_false": witness_false,
                "audit_set_at_truth": "valid" if r["z_valid_at_truth"] == "1" else "invalid",
                "audit_severity": r["severity"],
            }
        )
    with (LEDGER / "decision_table.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)
    exact_rows = [o for o in out_rows if o["first_failing_retraction"]]
    unique = sum(1 for o in exact_rows if " and " not in o["first_failing_retraction"])
    named_false = sum(1 for o in exact_rows if o["audit_named_claim_false"])
    lines = [
        f"# Table 8 — the decision table, primary supplier (`{args.arm}`)",
        "",
        f"{len(out_rows)} certifiable queries, {len(exact_rows)} with an exact certificate. The"
        " first-failing retraction set is named in the analyst's variables; it is a single claim on"
        f" {unique} of the {len(exact_rows)} exact rows. Columns left of the double rule are"
        " computed from observables; the audit columns read the truth and say whether the named"
        f" claim was false, which it was on {named_false} of {len(exact_rows)}."
        f" First {min(args.max_rows, len(exact_rows))} exact rows shown; the full table is"
        " `decision_table.csv`.",
        "",
    ]
    head = [
        "network",
        "query",
        "r_claim",
        "holds unless retracted",
        "free rule",
        "r_hop",
        "‖ named claim false",
        "set at truth",
    ]
    lines.append("| " + " | ".join(head) + " |")
    lines.append("|" + "---|" * len(head))
    for o in exact_rows[: args.max_rows]:
        cells = [
            o["network"],
            o["query"],
            o["r_claim"],
            o["first_failing_retraction"],
            o["free_rule"],
            o["r_hop"],
            "yes" if o["audit_named_claim_false"] else "no",
            o["audit_set_at_truth"],
        ]
        lines.append("| " + " | ".join(str(c) for c in cells) + " |")
    (LEDGER / "TABLE8_DECISION.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[:14]))


if __name__ == "__main__":
    main()
