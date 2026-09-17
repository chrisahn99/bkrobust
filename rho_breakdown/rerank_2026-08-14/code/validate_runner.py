"""
Faithfulness check: does run_ksweep.py's re-implementation reproduce the pilot's
own numbers when configured to the pilot's regime?

Config: pilot grid (p in 5..8, deg in {1.5,2,2.5}), |K| = 4 only, no
n_undirected >= 8 rejection (require >= 4 so |K|=4 is drawable). Compare the
rho=1 four-way decomposition against the pilot's published generic-arm figures
(RESULTS.md / summary.json): 33.0% Meek-inconsistent, 14.1% silent bias,
changed_still_valid = 0.

The two ensembles are NOT identical -- the pilot also admits |K|=3 SCMs, which
this arm excludes -- so agreement is expected to be close, not exact.
"""
import json
from multiprocessing import Pool

import numpy as np

import run_ksweep as R

R.K_SIZES = (4,)
R.NEED_U = 4


def main():
    rng = np.random.default_rng(20260814)
    jobs = [(int(s), int(rng.choice(list(range(5, 9)))),
             float(rng.choice([1.5, 2.0, 2.5]))) for s in range(20000)]
    with Pool(8) as pool:
        res = []
        for r in pool.imap(R.analyse_one, jobs, chunksize=16):
            if r and r["arms"].get("K4", {}).get("mpdag_amenable"):
                res.append(r)
                if len(res) >= 600:
                    break
        pool.terminate()

    ms = [m for r in res for m in r["arms"]["K4"]["members"] if m["rho"] == 1]
    n = len(ms)
    inc = sum(1 for m in ms if not m["consistent"])
    cons = [m for m in ms if m["consistent"]]
    na = sum(1 for m in cons if not m["amenable"])
    am = [m for m in cons if m["amenable"]]
    unch = sum(1 for m in am if not m["ostar_changed"])
    cv = sum(1 for m in am if m["ostar_changed"] and m["ostar_valid"])
    inv = sum(1 for m in am if m["ostar_changed"] and not m["ostar_valid"])

    out = {
        "n_scm": len(res), "n_rho1_members": n,
        "meek_inconsistent": round(inc / n, 4),
        "not_amenable": round(na / n, 4),
        "ostar_unchanged": round(unch / n, 4),
        "changed_still_valid": cv,
        "silent_bias": round(inv / n, 4),
        "PILOT_published": {"meek_inconsistent": 0.330, "silent_bias": 0.141,
                            "changed_still_valid": 0, "ostar_untouched": 0.63,
                            "loud_non_amenable": 0.23},
    }
    print(json.dumps(out, indent=1))
    with open("../results/validate_runner.json", "w") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
