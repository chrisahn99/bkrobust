"""E15: adjustment identification distance (AID) as a survival baseline.

AID (Henckel, Wuertzen & Weichwald 2024, optimal-set strategy), implemented
from its definition with the submitted GAC predicate: for every ordered pair
(a, b) of distinct nodes, the guess graph G proposes either "no effect"
(b not a possible descendant of a), "not identified" (G not amenable), or the
optimal set O(a, b, G). It is a mistake in the reference graph T when
  - G says no effect but b is a possible descendant of a in T, or
  - G says not identified but T is amenable, or
  - G proposes O and O is not a valid adjustment set in T (GAC), or
  - G proposes O but T is not amenable.
AID = mistakes / (n (n - 1)).

Per instance, with G = G0 (the analyst's state):
  aid_oracle      T = data-generating DAG   (needs the truth, like SHD)
  aid_analyst     T = CPDAG C_hat            (available to the analyst)
  aidx_oracle / aidx_analyst: the same restricted to pairs with treatment X
  (the query-specific variant).
Instances are rebuilt from their stored seeds with the submitted builders and
must reproduce the stored x, Z and r_val; failures are counted and dropped.
"""

from __future__ import annotations

import argparse
import csv
import multiprocessing as mp

from common import OUT, ROOT, UNREACHED, hybrid, is_gac_valid_mpdag
from oset import optimal_set_closed_form

from bkrobust.mpdag_criterion import is_amenable, possible_descendants
from bkrobust.robustness import survival as sv
from bkrobust.robustness.unconfounded_generator import build_unconfounded_tiered_instance

FILES = {
    "targeted": ROOT / "results/axis_robustness_p6/survival_instances.csv",
    "untargeted": ROOT / "results/axis_robustness_unconfounded_v2/survival_instances.csv",
}


def guess(g, a, b):
    if b not in possible_descendants(g, a):
        return ("zero", None)
    o, st = optimal_set_closed_form(g, a, b)
    if st == "not_amenable":
        return ("nid", None)
    if st == "y_not_descendant":  # possible but not definite descendant, amenable
        o = frozenset()
    return ("adj", o)


def mistake(t, a, b, gs):
    kind, o = gs
    if kind == "zero":
        return b in possible_descendants(t, a) and a != b and (b in possible_descendants(t, a))
    if kind == "nid":
        return is_amenable(t, a, b)
    if not is_amenable(t, a, b):
        return True
    return not is_gac_valid_mpdag(t, a, b, o)


def aid(g, t, only_x=None):
    nodes = list(g.nodes)
    m = n = 0
    for a in nodes:
        if only_x is not None and a != only_x:
            continue
        for b in nodes:
            if a == b:
                continue
            n += 1
            m += mistake(t, a, b, guess(g, a, b))
    return m / n if n else 0.0


def rebuild(r, which):
    if which == "untargeted":
        inst, why = build_unconfounded_tiered_instance(int(r["component_size"]), int(r["seed"]),
                                                       int(r["n_tiers"]))
    elif r["arm"] == "flip":
        inst, why = sv.build_flip_instance(int(r["component_size"]), int(r["separation"]),
                                           float(r["coverage"]), float(r["base_wrongness"]),
                                           int(r["seed"]))
    else:
        inst, why = sv.build_tiered_instance(int(r["component_size"]), int(r["separation"]),
                                             int(r["seed"]), int(r["n_tiers"]))
    return inst, why


def work(args):
    r, which = args
    out = {"which": which, "instance_id": r["instance_id"], "arm": r["arm"],
           "coverage": r["coverage"], "base_wrongness": r["base_wrongness"],
           "n_tiers": r["n_tiers"], "AUC_frac": r["AUC_frac"], "r_val_stored": r["r_val"],
           "separation": r["separation"], "shd_truth": r["shd_truth"], "n_k": r["n_k"],
           "k_g0": r["k_g0"], "phi_1": r["phi_1"]}
    inst, why = rebuild(r, which)
    if inst is None:
        out["reproduced"] = f"rebuild_failed:{why}"
        return out
    zs = ",".join(sorted(inst["z_star"]))
    if inst["instance_id"] != r["instance_id"] or inst["x"] != r["x"] or zs != r["z_star"]:
        out["reproduced"] = "mismatch_query"
        return out
    h = hybrid(inst["cpdag"], inst["g0"], inst["x"], inst["y"], inst["z_star"])
    rv = "-1" if h.radius == UNREACHED else str(h.radius)
    out["r_val"] = rv
    out["reproduced"] = "ok" if rv == r["r_val"] else f"radius_mismatch:{rv}"
    g0, dag, cpdag, x = inst["g0"], inst["dag"], inst["cpdag"], inst["x"]
    out["aid_oracle"] = aid(g0, dag)
    out["aid_analyst"] = aid(g0, cpdag)
    out["aidx_oracle"] = aid(g0, dag, only_x=x)
    out["aidx_analyst"] = aid(g0, cpdag, only_x=x)
    out["n_nodes"] = len(g0.nodes)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=3)
    a = ap.parse_args()
    todo = []
    for which, p in FILES.items():
        with open(p) as f:
            todo += [(r, which) for r in csv.DictReader(f)]
    print(len(todo), "instances", flush=True)
    res = []
    with mp.Pool(a.workers, maxtasksperchild=200) as P:
        for i, o in enumerate(P.imap_unordered(work, todo, chunksize=8)):
            res.append(o)
            if i % 500 == 0:
                print(i, o["which"], o.get("reproduced"), o.get("aid_oracle"), flush=True)
    keys = sorted({k for o in res for k in o})
    path = OUT / "e15" / "aid_instances.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(res)
    print("wrote", path)
