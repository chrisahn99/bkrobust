"""Changed-draws-only survival ranking, per stratum, beside the whole-sample value.

Reads the per-(instance, arm, intensity) aggregate written by
``experiments/paired_synth_v2.py`` (``_changed_state_agg.tsv``) and rebuilds
each instance's changed-only survival AUC with the same ``auc_frac_from_S``.
Reports Kendall tau_b(r_val, AUC) per paired-table stratum and pooled per arm,
for the whole sample (committed ``AUC_frac``) and for changed draws only, with
instance-bootstrap 95% intervals.
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import paired_synth_v2 as P  # noqa: E402

N_BOOT = 2000
OUT = P.OUT_DIR / "changed_state_strata.json"


def tau(x: np.ndarray, y: np.ndarray) -> float:
    return P._kendall_tau_b(x, y)


def main() -> int:
    rows = P.load_instances()
    by_inst: dict[str, dict[int, float | None]] = defaultdict(dict)
    for line in (P.OUT_DIR / "_changed_state_agg.tsv").open():
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 5:
            continue
        k = parts[0].split("\x1c")
        if len(k) != 3:
            continue
        n_changed, n_surv = int(parts[3]), int(parts[4])
        by_inst[k[0]][int(k[2])] = (n_surv / n_changed) if n_changed else None
    recs = []
    for r in rows:
        if r["r_val"] in ("", None) or r["AUC_frac"] in ("", None):
            continue
        auc_c = P.auc_frac_from_S(by_inst.get(r["instance_id"], {}), int(r["n_k"]))
        if auc_c is None:
            continue
        recs.append((P.stratum_of(r), r["arm"], float(r["r_val"]), float(r["AUC_frac"]), auc_c))
    groups = defaultdict(list)
    for s, arm, rv, a, ac in recs:
        groups[s].append((rv, a, ac))
        groups[f"pooled_{arm}"].append((rv, a, ac))
    rng = np.random.default_rng(0)
    out = {}
    for g, v in sorted(groups.items()):
        arr = np.array(v)
        pt = [tau(arr[:, 0], arr[:, 1]), tau(arr[:, 0], arr[:, 2])]
        boots = []
        for _ in range(N_BOOT):
            i = rng.integers(0, len(arr), len(arr))
            boots.append([tau(arr[i, 0], arr[i, 1]), tau(arr[i, 0], arr[i, 2])])
        b = np.array(boots)
        out[g] = {"n": len(arr), "tau_whole": pt[0], "tau_changed": pt[1],
                  "ci_whole": list(np.nanpercentile(b[:, 0], [2.5, 97.5])),
                  "ci_changed": list(np.nanpercentile(b[:, 1], [2.5, 97.5]))}
        print(f"{g:24s} n={len(arr):5d} whole={pt[0]:.3f} changed={pt[1]:.3f} "
              f"ci_changed=[{out[g]['ci_changed'][0]:.3f},{out[g]['ci_changed'][1]:.3f}]")
    OUT.write_text(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
