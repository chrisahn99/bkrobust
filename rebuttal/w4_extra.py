"""C2Bd-W4: in-text numbers of App. H (within-state paragraph, cohort paragraph)
on the original and the full (corrected) cohort, beyond results/rebuttal/w4/tables.json.

LOO-network range of the within-state tau_b, networks with positive value,
states with C>D / C<D / C==D, permutation test within states, share of
within-state units with separation undefined, r_val==s agreement share,
balanced 8/8 paired leads, pooled paired lead over separation. Same cohort
builder and helpers as w4_corrected_tables.py (10,000 resamples, seed 0).
On the original cohort it reproduces every published value exactly.

Output: results/rebuttal/w4/extra.json. Laptop timing: about 18 min with
NPERM=20 (dominated by the 18 leave-one-network-out bootstraps per cohort),
plus about 2 min for the 10,000 permutations.

    PYTHONPATH=src:rebuttal .venv/bin/python rebuttal/w4_extra.py [NPERM]
"""
import json
import sys
from collections import defaultdict

import numpy as np

from common import OUT, ROOT  # noqa: E402

sys.path.insert(0, str(ROOT / "experiments"))
import w4_corrected_tables as w4  # noqa: E402
import llm_paired_v2 as base  # noqa: E402
from bkrobust.robustness.paired_resample import paired_comparison  # noqa: E402

E, NB, SEED = "AUC_frac", 10_000, 0
NPERM = int(sys.argv[1]) if len(sys.argv) > 1 else 10_000  # paper: 10,000


def stau(st):
    return base._stratified_tau_b(st, "radius", E)


def extra(units):
    st = w4.states(units)
    out = {"n_states_total_any": len(w4.states(units, require_var=False)),
           "n_states_var": len(st), "n_units": sum(map(len, st.values())),
           "tau": stau(st)["tau_b"]}
    nets = sorted({k[0] for k in st})
    loo = {}
    for n in nets:
        sub = {k: v for k, v in st.items() if k[0] != n}
        loo[n] = w4.strat(sub, "radius")
    lo_n = min(loo, key=lambda n: loo[n]["tau_b"])
    hi_n = max(loo, key=lambda n: loo[n]["tau_b"])
    out["loo_min"] = {"net": lo_n, **{k: loo[lo_n][k] for k in ("tau_b", "lo", "hi")}}
    out["loo_max"] = {"net": hi_n, **{k: loo[hi_n][k] for k in ("tau_b", "lo", "hi")}}
    out["loo_any_ci_includes_zero"] = any(v["lo"] <= 0 for v in loo.values())
    pernet = {n: stau({k: v for k, v in st.items() if k[0] == n})["tau_b"] for n in nets}
    out["n_networks"] = len(nets)
    out["n_networks_positive"] = sum(1 for v in pernet.values() if v is not None and v > 0)
    cd = [(stau({k: v})["C"], stau({k: v})["D"]) for k, v in st.items()]
    out["states_C_gt_D"] = sum(c > d for c, d in cd)
    out["states_C_lt_D"] = sum(c < d for c, d in cd)
    out["states_C_eq_D"] = sum(c == d for c, d in cd)
    # permutation: shuffle radius within states
    rng = np.random.default_rng(SEED)
    obs = out["tau"]
    ge = 0
    for _ in range(NPERM):
        perm = {}
        for k, v in st.items():
            r = [u["radius"] for u in v]
            rng.shuffle(r)
            perm[k] = [{**u, "radius": ri} for u, ri in zip(v, r)]
        if stau(perm)["tau_b"] >= obs:
            ge += 1
    out["perm_n"], out["perm_ge_obs"] = NPERM, ge
    sep_units = [u for v in st.values() for u in v if u.get("separation") is not None]
    out["share_within_sep_undefined"] = 1 - len(sep_units) / out["n_units"]
    out["agree_r_eq_s_within_sepdef"] = float(np.mean(
        [float(u["radius"]) == float(u["separation"]) for u in sep_units]))
    okt = base.ok_conditions_by_triple(units)
    bal = base.subset_at_least(units, okt, 8)
    out["balanced_paired"] = {}
    for other in ("k_g0", "n_k"):
        r = paired_comparison(bal, "radius", other, E, "network", n_boot=NB, seed=SEED)
        out["balanced_paired"][other] = {k: r.get(k) for k in ("delta_point", "ci_lo_2p5", "ci_hi_97p5", "n")}
    pooled = base.real_naming(units)
    r = paired_comparison(pooled, "radius", "separation", E, "network", n_boot=NB, seed=SEED)
    out["pooled_paired_sep"] = {k: r.get(k) for k in ("delta_point", "ci_lo_2p5", "ci_hi_97p5", "n")}
    return out


if __name__ == "__main__":
    old, comb, _ = w4.cohort()
    res = {"original": extra(old), "corrected": extra(comb)}
    out = OUT / "w4" / "extra.json"
    out.write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps(res, indent=1, default=float))
