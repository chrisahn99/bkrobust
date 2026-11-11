"""E10 analysis: radius spread across valid sets, and how often the committed set O is beaten.

Reads rows_networks*.jsonl and rows_census*.jsonl in --out, computes the registered statistic
of PREREG.md with its Clopper-Pearson and network-clustered bootstrap intervals, prints a small
table and writes results.json beside the rows. Works on the smoke output too.

    .venv/bin/python experiments/rebuttal/E10_radius_across_valid_sets/analyse.py --out results/rebuttal/E10
"""

from __future__ import annotations

import argparse
import collections
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import beta

ROOT = Path(__file__).resolve().parents[3]
UNREACHED = -1  # bkrobust.core.conventions.UNREACHED; the sentinel, never a number
THRESHOLD = 0.05  # PREREG.md: O beaten on more than 5% of queries goes against us
N_BOOT = 2000


def more_robust(a: int, b: int) -> bool:
    if a == UNREACHED:
        return b != UNREACHED
    if b == UNREACHED:
        return False
    return a > b


def beats(a: dict, o: dict) -> bool | None:
    """Whether set a is certified strictly more robust than O; None when censoring hides it."""
    if o["exact"]:
        if o["r"] == UNREACHED:
            return False
        if a["exact"]:
            return more_robust(a["r"], o["r"])
        return True if a["lb"] is not None and a["lb"] > o["r"] else None
    if a["exact"] and a["r"] != UNREACHED and o["lb"] is not None and a["r"] <= o["lb"]:
        return False
    return None


def differ(a: dict, b: dict) -> bool | None:
    if a["exact"] and b["exact"]:
        return a["r"] != b["r"]
    for u, v in ((a, b), (b, a)):
        if u["exact"] and not v["exact"] and u["r"] != UNREACHED and v["lb"] is not None and v["lb"] > u["r"]:
            return True
    return None


def verdict(values) -> str:
    values = list(values)
    if any(v is True for v in values):
        return "yes"
    if all(v is False for v in values):
        return "no"
    return "undecided"


def clopper_pearson(k: int, n: int, level: float = 0.95) -> list[float | None]:
    if n == 0:
        return [None, None]
    a = (1 - level) / 2
    lo = 0.0 if k == 0 else float(beta.ppf(a, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(1 - a, k + 1, n - k))
    return [lo, hi]


def cluster_bootstrap(units: list[dict], seed: int = 0) -> list[float | None]:
    """Percentile interval of beaten/decided, resampling networks with replacement."""
    by_net: dict[str, list[int]] = collections.defaultdict(lambda: [0, 0])
    for u in units:
        if u["o_beaten"] != "undecided":
            by_net[u["network"]][0] += u["o_beaten"] == "yes"
            by_net[u["network"]][1] += 1
    nets = sorted(by_net)
    if len(nets) < 2:
        return [None, None]
    k = np.array([by_net[n][0] for n in nets])
    m = np.array([by_net[n][1] for n in nets])
    rng = np.random.default_rng(seed)
    draws = rng.integers(len(nets), size=(N_BOOT, len(nets)))
    num, den = k[draws].sum(1), m[draws].sum(1)
    frac = num[den > 0] / den[den > 0]
    return [float(np.quantile(frac, 0.025)), float(np.quantile(frac, 0.975))]


def load(out: Path, part: str) -> list[dict]:
    rows: dict[str, dict] = {}
    for f in sorted(out.glob(f"rows_{part}*.jsonl")):
        for line in f.open():
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            if rec.get("status") == "error" and rec["key"] in rows:
                continue  # a later failed retry never hides an earlier success
            rows[rec["key"]] = rec
    return list(rows.values())


# ------------------------------------------------------------------ networks


def summarise_unit(rec: dict, include_empty: bool) -> dict:
    sets = rec["sets"]
    o = sets[0]
    alts = [s for s in sets[1:] if include_empty or s["size"] > 0]
    exact_r = [s["r"] for s in sets if s["exact"]]
    finite = [r for r in exact_r if r != UNREACHED]
    pairs = [differ(a, b) for i, a in enumerate(sets) for b in sets[i + 1:]]
    return {
        "network": rec["network"], "arm": rec["arm"], "mode": rec["mode"], "n_sets": len(sets),
        "n_alts": len(alts),
        "o_beaten": verdict(beats(a, o) for a in alts) if alts else "no_alternative",
        "n_beating": sum(beats(a, o) is True for a in alts),
        "nonflat": verdict(pairs) if pairs else "no_alternative",
        "spread_finite": (max(finite) - min(finite)) if finite else None,
        "mixed_unreached": bool(finite) and UNREACHED in exact_r,
        "n_censored": sum(not s["exact"] for s in sets),
        "o_exact": o["exact"], "o_valid_at_truth": o["valid_at_truth"],
        "n_alts_invalid_at_truth": sum(not s["valid_at_truth"] for s in alts),
    }


def network_stats(units: list[dict]) -> dict:
    elig = [u for u in units if u["o_beaten"] != "no_alternative"]
    dec = [u for u in elig if u["o_beaten"] != "undecided"]
    k = sum(u["o_beaten"] == "yes" for u in dec)
    n_und = len(elig) - len(dec)
    nonflat_dec = [u for u in elig if u["nonflat"] != "undecided"]
    spreads = [u["spread_finite"] for u in elig if u["spread_finite"] is not None]
    return {
        "n_units": len(units), "n_with_alternative": len(elig), "n_decided": len(dec),
        "n_undecided": n_und, "n_o_beaten": k,
        "frac_o_beaten": k / len(dec) if dec else None,
        "cp95": clopper_pearson(k, len(dec)),
        "cluster_boot95": cluster_bootstrap(elig),
        "worst_case_frac": (k + n_und) / len(elig) if elig else None,
        "frac_undecided": n_und / len(elig) if elig else None,
        "n_nonflat": sum(u["nonflat"] == "yes" for u in nonflat_dec),
        "n_nonflat_decided": len(nonflat_dec),
        "frac_nonflat": (sum(u["nonflat"] == "yes" for u in nonflat_dec) / len(nonflat_dec))
        if nonflat_dec else None,
        "max_finite_spread": max(spreads) if spreads else None,
        "spread_hist": dict(sorted(collections.Counter(spreads).items())),
        "n_mixed_unreached": sum(u["mixed_unreached"] for u in elig),
        "n_networks": len({u["network"] for u in elig}),
    }


def analyse_networks(rows: list[dict]) -> dict:
    ok = [r for r in rows if r.get("status") == "ok"]
    n_sets = sum(r["n_sets"] for r in ok)
    n_cens = sum(not s["exact"] for r in ok for s in r["sets"])
    res: dict = {
        "n_rows": len(rows),
        "status": dict(collections.Counter(r.get("status") for r in rows)),
        "n_hard": sum(bool(r.get("hard_failures")) for r in rows),
        "n_sets": n_sets, "n_sets_censored": n_cens,
        "frac_sets_censored": n_cens / n_sets if n_sets else None,
        "modes": dict(collections.Counter(r["mode"] for r in ok)),
        "methods": dict(collections.Counter(s["method"] for r in ok for s in r["sets"])),
    }
    for label, inc in (("registered_nonempty_alternatives", False), ("secondary_with_empty", True)):
        units = [summarise_unit(r, inc) for r in ok]
        res[label] = network_stats(units)
        if not inc:
            res["by_arm_family"] = {
                fam: network_stats([u for u in units if (u["arm"] == "D_RAND") == (fam == "D_RAND")])
                for fam in ("elicited", "D_RAND")}
            res["by_mode"] = {m: network_stats([u for u in units if u["mode"] == m])
                              for m in sorted({u["mode"] for u in units})}
            res["by_network"] = {n: {k: v for k, v in network_stats([u for u in units if u["network"] == n]).items()
                                     if k in ("n_with_alternative", "n_decided", "n_o_beaten", "n_nonflat",
                                              "max_finite_spread")}
                                 for n in sorted({u["network"] for u in units})}
    vc = [r for r in ok if "variance_cost" in r]
    ratios = [x for r in vc for x in r["variance_cost"]["avar_ratio_best_over_o"] if x is not None]
    res["beaten_variance_cost"] = {
        "n_units": len(vc),
        "median_avar_ratio": float(np.median(ratios)) if ratios else None,
        "frac_best_valid_at_truth": (sum(r["variance_cost"]["best_valid_at_truth"] for r in vc) / len(vc))
        if vc else None,
        "examples": [{"key": r["key"], "z_o": r["z_o"], **r["variance_cost"]} for r in vc[:10]],
    }
    return res


# ------------------------------------------------------------------ census


def analyse_census(rows: list[dict]) -> dict:
    ok = [r for r in rows if r.get("status") == "ok"]
    out: dict = {"n_rows": len(rows), "status": dict(collections.Counter(r.get("status") for r in rows)),
                 "n_hard": sum(bool(r.get("hard_failures")) for r in rows), "by_n": {}}
    for n in sorted({r["n"] for r in ok}):
        c = collections.Counter()
        examples = []
        for r in (r for r in ok if r["n"] == n):
            c["cpdags"] += 1
            c["censored_hybrid"] += r.get("n_censored", 0)
            for inst in r["instances"]:
                # set = [Z, hybrid radius or None, method, exact space-BFS radius]; the BFS radius
                # is exact everywhere, so comparisons use it.
                c["instances_ge2_sets"] += 1
                radii = {s[0]: s[3] for s in inst["sets"]}
                nonempty = {z: r_ for z, r_ in radii.items() if z}
                if len(nonempty) >= 2:
                    c["frontier_instances"] += 1
                    c["frontier_nonflat"] += len(set(nonempty.values())) > 1
                if len(set(radii.values())) > 1:
                    c["allsets_nonflat"] += 1
                for tag, o in (("closed_form", inst["o"] if inst["verdict"] == "closed_form" else None),
                               ("enum", inst["o_enum"])):
                    if o is None or o not in radii:
                        continue
                    r_o = radii[o]
                    # frontier convention: O non-empty, alternatives non-empty, at least two sets
                    if o and len(nonempty) >= 2:
                        c[f"frontier_o_{tag}_defined"] += 1
                        c[f"frontier_o_{tag}_unreached"] += r_o == UNREACHED
                        if any(more_robust(v, r_o) for z, v in nonempty.items() if z != o):
                            c[f"frontier_o_{tag}_beaten"] += 1
                            if tag == "closed_form" and len(examples) < 5:
                                examples.append({"cpdag": r["cpdag"], **inst})
                    c[f"allsets_o_{tag}_defined"] += 1
                    if any(more_robust(v, r_o) for z, v in radii.items() if z != o):
                        c[f"allsets_o_{tag}_beaten"] += 1
        d = dict(sorted(c.items()))
        if d.get("frontier_instances"):
            d["frontier_frac_nonflat"] = d.get("frontier_nonflat", 0) / d["frontier_instances"]
        d["beaten_examples"] = examples
        out["by_n"][str(n)] = d
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", type=Path, default=ROOT / "results/rebuttal/E10")
    args = ap.parse_args(argv)
    net_rows, cen_rows = load(args.out, "networks"), load(args.out, "census")
    if not net_rows and not cen_rows:
        print(f"no rows in {args.out}", file=sys.stderr)
        return 1
    cfg_files = sorted(args.out.glob("config*.json"))
    res: dict = {"experiment": "E10", "criterion": "gac",
                 "config": json.loads(cfg_files[0].read_text()) if cfg_files else None,
                 "networks": analyse_networks(net_rows) if net_rows else None,
                 "census": analyse_census(cen_rows) if cen_rows else None}

    reg = (res["networks"] or {}).get("registered_nonempty_alternatives")
    census_beaten = sum(d.get("frontier_o_closed_form_beaten", 0) for d in (res["census"] or {"by_n": {}})["by_n"].values())
    res["registered"] = {
        "statistic": "fraction of network queries (>= 1 non-empty valid alternative, comparison decided) "
                     "where a listed valid set is certified strictly more robust than O",
        "value": reg["frac_o_beaten"] if reg else None,
        "cp95": reg["cp95"] if reg else None,
        "cluster_boot95": reg["cluster_boot95"] if reg else None,
        "census_o_beaten": census_beaten,
        "against_authors": bool((reg and reg["frac_o_beaten"] is not None and reg["frac_o_beaten"] > THRESHOLD)
                                or census_beaten > 0),
        "more_than_half_undecided": bool(reg and reg["frac_undecided"] is not None and reg["frac_undecided"] > 0.5),
        "consistency_failures": ((res["networks"] or {}).get("n_hard", 0) + (res["census"] or {}).get("n_hard", 0)),
    }
    (args.out / "results.json").write_text(json.dumps(res, indent=1, sort_keys=True))

    def f(x):
        return "-" if x is None else (f"{x:.4f}" if isinstance(x, float) else str(x))

    print("E10, criterion gac")
    if reg:
        print(f"{'stratum':<22}{'queries':>8}{'decided':>8}{'beaten':>8}{'frac':>8}  {'CP95':<18}{'nonflat':>9}{'maxspread':>10}")
        strata = [("all (registered)", reg)]
        strata += [(f"arm {k}", v) for k, v in res["networks"]["by_arm_family"].items()]
        strata += [(f"mode {k}", v) for k, v in res["networks"]["by_mode"].items()]
        strata += [("with empty set", res["networks"]["secondary_with_empty"])]
        for name, s in strata:
            ci = "[" + ", ".join(f(v) for v in s["cp95"]) + "]"
            print(f"{name:<22}{s['n_with_alternative']:>8}{s['n_decided']:>8}{s['n_o_beaten']:>8}"
                  f"{f(s['frac_o_beaten']):>8}  {ci:<18}{f(s['frac_nonflat']):>9}{f(s['max_finite_spread']):>10}")
        n = res["networks"]
        print(f"cluster bootstrap 95%: {reg['cluster_boot95']}; worst case (undecided = beaten): "
              f"{f(reg['worst_case_frac'])}; undecided {reg['n_undecided']}; censored sets "
              f"{n['n_sets_censored']}/{n['n_sets']}; consistency failures {n['n_hard']}")
    if res["census"]:
        print(f"{'census n':<10}{'inst>=2':>9}{'frontier':>9}{'nonflat':>9}{'O def':>7}{'O beaten':>9}{'O enum beaten':>14}")
        for k, d in res["census"]["by_n"].items():
            print(f"{k:<10}{d.get('instances_ge2_sets', 0):>9}{d.get('frontier_instances', 0):>9}"
                  f"{d.get('frontier_nonflat', 0):>9}{d.get('frontier_o_closed_form_defined', 0):>7}"
                  f"{d.get('frontier_o_closed_form_beaten', 0):>9}{d.get('frontier_o_enum_beaten', 0):>14}")
        print(f"census consistency failures (hybrid vs space BFS): {res['census']['n_hard']}")
    r = res["registered"]
    print(f"REGISTERED: O beaten fraction {f(r['value'])} CP95 {r['cp95']}; census O beaten {r['census_o_beaten']}; "
          f"against authors: {r['against_authors']}")
    print(f"wrote {args.out / 'results.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
