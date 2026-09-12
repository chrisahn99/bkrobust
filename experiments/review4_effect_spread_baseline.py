"""Round-3 review, W9: benchmark against the diagnostic a reader would reach for first.

The paper tests the radius against SHD to the truth (uncomputable in practice), the
count of asserted claims (pre-registered as inert) and the separation. It never tests
it against the obvious competitor from the literature it cites: the spread of total
effects still possible given the analyst's knowledge, i.e. the multiset
`{effect of X on Y in D : D in [G0]}` (Guo & Perkovic; IDA with background knowledge).
That quantity is already a structural-risk diagnostic, is computable from G0 and the
sample, and needs no breakdown-radius machinery.

For each survival instance we rebuild G0 deterministically, enumerate its DAG
extensions, compute the optimal-set-adjusted total effect in each under a shared SEM,
and take the spread. We then score that predictor against the same survival AUC the
radius is scored against, with the same paired bootstrap.

A wide spread means the knowledge leaves the effect badly undetermined; the radius
instead measures how far the knowledge can be perturbed before the *chosen set* breaks.
Whether the second buys anything over the first is the question.
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd
from scipy.stats import kendalltau

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from bkrobust.core.conventions import UNREACHED  # noqa: E402
from bkrobust.demo.evaluate import (  # noqa: E402
    optimal_adjustment_set_dag,
    random_sem,
)
from bkrobust.demo.meek import enumerate_dag_extensions  # noqa: E402
from bkrobust.robustness import analyse as an  # noqa: E402
from bkrobust.robustness import build_tau_comparisons as btc  # noqa: E402
from bkrobust.robustness import survival as sv  # noqa: E402

R = ROOT / "results" / "axis_robustness"
OUT = R / "review4_effect_spread.csv"
MAX_EXT = 4000      # extension enumeration cap; instances above it are reported, not dropped
N_SEM = 4           # SEM draws per instance; the spread is a population quantity


def adjusted_effect(sem, x, y, z):
    """Population estimand for the effect of ``x`` on ``y`` adjusting for ``z``."""
    from bkrobust.demo.evaluate import adjusted_estimand

    try:
        return float(adjusted_estimand(sem, x, y, frozenset(z)))
    except Exception:
        return float("nan")


def spread_for(inst) -> tuple[float, int]:
    """Range of achievable adjusted effects over the DAG extensions of ``G0``."""
    g0, dag, x, y = inst["g0"], inst["dag"], inst["x"], inst["y"]
    exts = enumerate_dag_extensions(g0)
    if not exts or len(exts) > MAX_EXT:
        return float("nan"), len(exts)
    vals = []
    for k in range(N_SEM):
        rng = np.random.default_rng(abs(hash((inst["instance_id"], "sem", k))) % (2**32))
        sem = random_sem(dag, rng)
        eff = [adjusted_effect(sem, x, y, optimal_adjustment_set_dag(d, x, y)) for d in exts]
        eff = [e for e in eff if e == e]
        if len(eff) >= 2:
            vals.append(max(eff) - min(eff))
    return (float(np.median(vals)) if vals else float("nan")), len(exts)


def main() -> None:
    inst_df = pd.read_csv(R / "survival_instances.csv")
    flip = inst_df[inst_df.arm == "flip"]
    print(f"{len(flip)} flip instances to rebuild")

    rows = []
    for i, r in enumerate(flip.itertuples(), 1):
        built, reason = sv.build_flip_instance(
            int(r.component_size), int(r.separation), float(r.coverage),
            float(r.base_wrongness), int(r.seed))
        if built is None or built["instance_id"] != r.instance_id:
            rows.append({"instance_id": r.instance_id, "effect_spread": float("nan"),
                         "n_extensions": -1, "status": "rebuild_mismatch"})
            continue
        sp, ne = spread_for(built)
        rows.append({"instance_id": r.instance_id, "effect_spread": sp,
                     "n_extensions": ne, "status": "ok" if sp == sp else "too_many_extensions"})
        if i % 200 == 0:
            print(f"  {i}/{len(flip)}", flush=True)

    sp_df = pd.DataFrame(rows)
    sp_df.to_csv(OUT, index=False)
    ok = sp_df[sp_df.status == "ok"]
    print(f"\nrebuilt {len(sp_df)}, usable spread on {len(ok)} "
          f"({(sp_df.status == 'rebuild_mismatch').sum()} rebuild mismatches)")

    # --- score it against the same endpoint the radius is scored against ---------
    curves = pd.read_csv(R / "survival_curves.csv")
    d = an.attach_auc_frac_usable(inst_df, curves).merge(sp_df, on="instance_id", how="left")
    print(f"\n{'stratum':30s} {'tau(r_hop)':>11s} {'tau(spread)':>12s} {'n':>5s}")
    out = []
    for key, g in d[d.arm == "flip"].groupby(["coverage", "base_wrongness"]):
        s = g[["r_val", "effect_spread", "AUC_frac_usable"]].copy()
        s = s[s.r_val != UNREACHED].dropna()
        if len(s) < 30:
            continue
        t_r = kendalltau(s.r_val, s.AUC_frac_usable, variant="b")[0]
        t_s = kendalltau(s.effect_spread, s.AUC_frac_usable, variant="b")[0]
        out.append({"stratum": str(key), "n": len(s), "tau_r_hop": t_r,
                    "tau_effect_spread": t_s})
        print(f"{str(key):30s} {t_r:+11.3f} {t_s:+12.3f} {len(s):5d}")
    pd.DataFrame(out).to_csv(R / "review4_effect_spread_tau.csv", index=False)
    print(f"\nwrote {OUT} and review4_effect_spread_tau.csv")


if __name__ == "__main__":
    main()
