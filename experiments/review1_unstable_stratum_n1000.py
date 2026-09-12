"""Review round 1, weakness 3 / [RE-12]: retire the second unstable cell.

`table_tau_comparisons.md` carries an "unstable at N = 200" marker on two cells.
The null cell (flip, coverage 0.5, base wrongness 0.25) was re-run at N = 1000 and
is now a confirmed null. The other, flip coverage 1.0 base wrongness 0.00, has
never been re-run, and the paper therefore cannot state a verdict for it.

`nullcell` hardcodes the stratum as module constants. This rebinds them and
reuses the committed driver unchanged, so the design, the seeds, the endpoint
and the bootstrap are identical to the published null-cell run; only which
stratum is measured differs.
"""
from __future__ import annotations

import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from bkrobust.robustness import nullcell as nc  # noqa: E402
from bkrobust.robustness import run_nullcell as rn  # noqa: E402

COVERAGE = 1.0
BASE_WRONGNESS = 0.0
OUT = ROOT / "results" / "axis_robustness" / "review1_unstable_cell"


def main() -> int:
    nc.COVERAGE = COVERAGE
    nc.BASE_WRONGNESS = BASE_WRONGNESS
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"re-running stratum flip coverage={COVERAGE} base_wrongness={BASE_WRONGNESS} at N=1000")
    return rn.run_full_mode(OUT)


if __name__ == "__main__":
    raise SystemExit(main())
