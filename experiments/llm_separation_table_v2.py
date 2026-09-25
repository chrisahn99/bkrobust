"""Separation as a Table 1 predictor on the LLM panel: pooled, balanced 8/8, within state.

Reuses experiments/llm_paired_v2.py so the panel definitions and the
network-cluster bootstrap match the committed table exactly.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import llm_paired_v2 as L  # noqa: E402

OUT = L.OUT_DIR / "separation_table.json"


def main() -> int:
    units = L.load_units()
    pooled = L.real_naming(units)
    balanced = L.subset_at_least(units, L.ok_conditions_by_triple(units), 8)
    out = {}
    for label, rows in (("pooled", pooled), ("balanced_8of8", balanced)):
        out[label] = {
            "separation": L.tau_for_stratum(rows, "separation", L.ENDPOINT,
                                            n_boot=L.N_BOOT, seed=L.SEED),
            "paired_radius_minus_separation": L.paired_comparison(
                rows, "radius", "separation", L.ENDPOINT, "network",
                n_boot=L.N_BOOT, seed=L.SEED),
        }
    OUT.write_text(json.dumps(out, indent=2, default=L._json_default))
    for label, d in out.items():
        s, p = d["separation"], d["paired_radius_minus_separation"]
        print(label, "sep tau", s["tau_b"], s["ci_lo_2p5"], s["ci_hi_97p5"], "n", s["n"],
              "nets", s.get("n_networks"), "| paired", p["delta_point"], p["ci_lo_2p5"],
              p["ci_hi_97p5"], "n", p["n"], "radius tau on matched", p["tau_a_point"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
