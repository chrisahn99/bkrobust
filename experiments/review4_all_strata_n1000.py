"""[RE-12] in full: re-run every flip stratum at N = 1000, not just the two marked cells.

Both cells the published table flagged as unstable moved materially when re-run at
N = 1000, and one changed sign. That is 2 for 2, so the remaining flip strata are
still reported at a draw count the project has itself shown distorts the restricted
endpoint. This closes the gap by rebinding `nullcell`'s stratum constants and reusing
the committed driver unchanged, exactly as the single-cell re-run did.

Each cell carries its own reproduction check: the driver recomputes tau on the first
200 of its 1000 draws, which must land on the published N = 200 value.
"""
from __future__ import annotations

import json
import pathlib
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from bkrobust.robustness import nullcell as nc  # noqa: E402
from bkrobust.robustness import run_nullcell as rn  # noqa: E402

OUT = ROOT / "results" / "axis_robustness" / "review4_strata_n1000"
# The two already re-run are skipped; the other four flip cells are the gap.
CELLS = [(0.5, 0.0), (0.5, 0.1), (1.0, 0.1), (1.0, 0.25)]


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    summary = []
    for cov, bw in CELLS:
        tag = f"cov{cov}_bw{bw}"
        d = OUT / tag
        if (d / "null_tau_comparison.csv").exists():
            print(f"{tag}: already present, skipping")
            continue
        d.mkdir(parents=True, exist_ok=True)
        nc.COVERAGE, nc.BASE_WRONGNESS = cov, bw
        print(f"\n=== {tag} at N=1000 ===", flush=True)
        t0 = time.perf_counter()
        rc = rn.run_full_mode(d)
        summary.append({"cell": tag, "coverage": cov, "base_wrongness": bw,
                        "returncode": rc, "seconds": round(time.perf_counter() - t0, 1)})
        (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
        print(f"{tag}: rc={rc} in {summary[-1]['seconds']}s", flush=True)
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
