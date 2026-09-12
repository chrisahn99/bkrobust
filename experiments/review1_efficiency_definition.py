"""Review round 1, weakness 2: is efficiency really invariant under truthful
background knowledge, or is that an artefact of how the optimal set is defined?

Phase 2 measured utility through `optimal_adjustment_set_mpdag`, which returns a
set only when every DAG extension agrees and `None` otherwise. Under that rule
any G0 admitting the true DAG either returns the true DAG's own optimal set or
returns nothing, so achievable variance cannot move and "0 of 24 strata vary" is
forced. This script re-runs the same bases, the same proposals and the same SEM
draws under the standard graphical definition

    O(x, y, G) = pa(cn(x, y, G)) \\ forb(x, y, G)

(`bkrobust.demo.optimal_mpdag`), which returns a set whenever G is amenable, and
asks the same question again.

Outputs `results/axis_robustness/review1_efficiency_definition.csv`, one row per
(base, proposal), with the variance under both definitions.
"""
from __future__ import annotations

import csv
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))

from bkrobust.demo.evaluate import (  # noqa: E402
    asymptotic_variance,
    optimal_adjustment_set_mpdag,
)
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.demo.optimal_mpdag import optimal_adjustment_set_mpdag_hpm  # noqa: E402
from bkrobust.gac import is_gac_valid_mpdag  # noqa: E402
from bkrobust.robustness import pareto as P  # noqa: E402

OUT = pathlib.Path(__file__).resolve().parent.parent / "results" / "axis_robustness"
N_SEM = 8  # variance is a population quantity here; draws only vary the weights


def avar(sem, x, y, z):
    try:
        v = asymptotic_variance(sem, x, y, z)
    except Exception:
        return None
    return v if v is not None and v == v else None


def main() -> None:
    bases = P.select_base_cpdags()
    print(f"{len(bases)} base CPDAGs")
    rows = []
    for bi, base in enumerate(bases):
        x, y = base.treatment, base.outcome
        menu = P.truthful_orientation_menu(base)
        if not menu:
            continue
        desc_x = base.dag_obj.descendants(x)
        anc_y = base.dag_obj.ancestors(y) | {y}
        proposals = P.sample_proposals(base, menu, desc_x, anc_y)
        sems = P.draw_sems(base, 0, n=N_SEM)
        for prop in proposals:
            g0 = apply_orientations(base.cpdag_obj, list(prop.orientations))
            if g0 is None:
                continue
            z_agree = optimal_adjustment_set_mpdag(g0, x, y)
            z_hpm = optimal_adjustment_set_mpdag_hpm(g0, x, y)
            hpm_valid = (
                is_gac_valid_mpdag(g0, x, y, z_hpm) if z_hpm is not None else False
            )
            for s in sems:
                rows.append(
                    {
                        "base_id": base.base_id,
                        "component_size": base.component_size,
                        "separation": base.separation,
                        "proposal_id": prop.proposal_id,
                        "proposal_size": prop.size,
                        "sem_k": s.k,
                        "z_agree": "" if z_agree is None else "|".join(sorted(z_agree)),
                        "z_agree_defined": int(z_agree is not None),
                        "avar_agree": avar(s.sem, x, y, z_agree) if z_agree is not None else "",
                        "z_hpm": "" if z_hpm is None else "|".join(sorted(z_hpm)),
                        "z_hpm_defined": int(z_hpm is not None),
                        "z_hpm_gac_valid": int(hpm_valid),
                        "avar_hpm": avar(s.sem, x, y, z_hpm) if z_hpm is not None else "",
                    }
                )
        print(f"  [{bi + 1}/{len(bases)}] {base.base_id}: {len(proposals)} proposals")

    path = OUT / "review1_efficiency_definition.csv"
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {path} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
