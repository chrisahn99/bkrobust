"""
Paired K4 -> K8 test (the marginal CIs printed by analyse_ksweep.py overlap by
construction and are the wrong comparison; K4 and K8 are measured on the SAME
SCM with K4 subset K8, so the difference is paired).

Reports, per ensemble and per criterion:
  - censored fraction at K4, K6, K8 on the paired subset
  - the PAIRED difference K8 - K4 with a bootstrap CI over SCMs
  - McNemar exact test on the discordant pairs
  - the same on the UNPAIRED per-arm subsets (what a deployed method would see)

Usage: python paired_test.py
"""
import json
import sys

import numpy as np
from scipy.stats import binomtest

sys.path.insert(0, ".")
from analyse import breakdown_radius  # noqa: E402

MAX_RHO = 3
RNG = np.random.default_rng(20260814)


def rstar(r, key, crit):
    a = r["arms"][key]
    return breakdown_radius(a, r["tau"], a["est0"], MAX_RHO, crit)


def boot_diff(d, B=20000):
    d = np.asarray(d, float)
    idx = RNG.integers(0, len(d), size=(B, len(d)))
    bs = d[idx].mean(axis=1)
    return float(d.mean()), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))


def run(label):
    res = json.load(open(f"../results/ksweep_{label}.json"))
    keys = ["K4", "K6", "K8"]
    paired = [r for r in res
              if all(r["arms"].get(k, {}).get("mpdag_amenable") for k in keys)]
    out = {"label": label, "n_raw": len(res), "n_paired": len(paired)}

    # composition of the rejection-sampled ensemble
    from collections import Counter
    out["paired_deg_dist"] = dict(sorted(Counter(r["deg"] for r in paired).items()))
    out["paired_p_dist"] = dict(sorted(Counter(r["p"] for r in paired).items()))

    for crit in ("sign", "rel50", "ident"):
        c = {}
        cens = {k: np.array([rstar(r, k, crit) > 3 for r in paired], float) for k in keys}
        for k in keys:
            c[f"censored_{k}_paired"] = round(float(cens[k].mean()), 4)
        d = cens["K8"] - cens["K4"]
        m, lo, hi = boot_diff(d)
        c["paired_diff_K8_minus_K4"] = f"{m:+.4f} [{lo:+.4f}, {hi:+.4f}]"
        b = int(((cens["K4"] == 1) & (cens["K8"] == 0)).sum())   # K4 cens, K8 not
        cc = int(((cens["K4"] == 0) & (cens["K8"] == 1)).sum())  # K8 cens, K4 not
        c["mcnemar_b_K4cens_K8free"] = b
        c["mcnemar_c_K8cens_K4free"] = cc
        c["mcnemar_p"] = (round(float(binomtest(b, b + cc, 0.5).pvalue), 5)
                          if b + cc else None)
        # UNPAIRED: each arm on every SCM where THAT arm is amenable
        for k in keys:
            sub = [r for r in res if r["arms"].get(k, {}).get("mpdag_amenable")]
            v = np.array([rstar(r, k, crit) > 3 for r in sub], float)
            c[f"censored_{k}_unpaired"] = f"{v.mean():.4f} (n={len(sub)})"
        out[crit] = c
    return out


def main():
    allout = {}
    for label in ("original", "licensed", "large"):
        allout[label] = run(label)
    print(json.dumps(allout, indent=1))
    with open("../results/paired_test.json", "w") as f:
        json.dump(allout, f, indent=1)


if __name__ == "__main__":
    main()
