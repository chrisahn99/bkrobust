"""Candidate instances for the paper's first figure, read off the ledger.

The figure wants one real query under elicited knowledge where the certificate
in claims and a count in hops disagree, with the first-failing retraction set
named in the analyst's own variables. The ledger has those rows: the hop
instrument certified a safe move that the truth says does not exist, while the
claim certificate said the committed set is one claim from failure. This
script lists them, with the supplier's claims marked against the truth, the
committed set, and the witness, so that a human can choose one. It reads the
audit columns and is therefore a post-hoc tool, not an instrument.

Writes ``results/ledger/figure1_candidates.md``.

    python experiments/figure1_candidates.py
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments"))

from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.mpdag_criterion.optimal import committed_adjustment_set  # noqa: E402
from ledger_sweep import KNOWLEDGE, TIER, load  # noqa: E402

ROWS = ROOT / "results" / "ledger" / "rows.csv"
OUT = ROOT / "results" / "ledger" / "figure1_candidates.md"
ELICITED = ("D_LLM", "D_LLM_INSTR", "D_LLM_SRC", "D_SCRAMBLED")


def main() -> None:
    """List the unit-gap rows and the free-rule over-certifications under elicited knowledge."""
    know = json.loads(KNOWLEDGE.read_text())
    rows = [r for r in csv.DictReader(ROWS.open()) if r["status"] == "ok" and r["arm"] in ELICITED]
    graphs: dict[str, tuple] = {}
    lines = [
        "# Figure 1 candidates, from the ledger",
        "",
        "Rows under elicited knowledge where an instrument certified a move the truth does not",
        "allow. Claims are the supplier's after the consistency repair; a claim marked `false`",
        "is wrong at the generating graph. `Z` is the committed set on the analyst's graph; the",
        "witness is the first-failing retraction set the claim certificate names. Post hoc: this",
        "file reads the audit columns and exists to choose an illustration, not to measure.",
        "",
    ]
    for title, pick in (
        (
            "## The unit gap: hop certificate dangerous, claim certificate safe",
            lambda r: r["bucket_B6_hop"] == "dangerous" and r["bucket_B7_claim"] == "safe",
        ),
        (
            "## The free rule over-certifies, the claim certificate does not",
            lambda r: r["bucket_B2_free_rule"] == "dangerous" and r["bucket_B7_claim"] == "safe",
        ),
    ):
        lines += [title, ""]
        chosen = [r for r in rows if pick(r)]
        chosen.sort(key=lambda r: (TIER.get(r["network"], "T1") != "T2", r["network"], r["arm"]))
        for r in chosen:
            net, arm = r["network"], r["arm"]
            if net not in graphs:
                graphs[net] = load(net)
            dag, cpdag, names = graphs[net]
            k = [tuple(e) for e in know[arm]["networks"][net]["k"]]
            g0 = apply_orientations(cpdag, k)
            z = committed_adjustment_set(g0, r["x"], r["y"]).z or frozenset()

            claims = ", ".join(
                f"{names.get(a, a)}->{names.get(b, b)}"
                + ("" if dag.is_directed_edge(a, b) else " (false)")
                for a, b in k
            )
            query = f"{names.get(r['x'], r['x'])} -> {names.get(r['y'], r['y'])}"
            zs = ", ".join(sorted(names.get(v, v) for v in z))
            truth = "valid" if r["z_valid_at_truth"] == "1" else "invalid"
            lines += [
                f"- **{net}** ({TIER.get(net, 'T1')}), `{arm}`: query **{query}**;"
                f" r_claim = {r['r_claim'] or r['r_claim_status']}, r_hop = {r['r_hop']},"
                f" free rule = {r['separation']}, |K| = {r['n_k']}, wrong = {r['n_wrong_claims']},"
                f" set at truth: {truth}, severity {r['severity']}.",
                f"  Z = {{{zs}}}; witness: {r['witness'] or '(none)'}; claims: {claims}.",
            ]
        lines += ["", f"{len(chosen)} rows.", ""]
    OUT.write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
