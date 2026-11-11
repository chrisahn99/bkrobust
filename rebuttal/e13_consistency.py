"""E13: consistency checks on committed sets and the error budget.

(a) Survival-panel units classified "no optimal set" (status_split.json, the
    898 real-name units): amenability at G0, closed-form O, its GAC validity.
(b) Closed-form O against the submitted enumeration-based O, and against the
    release row's status and |Z|, on every elicited release row.
(c) Corollary 1(c) on release rows invalid at truth: the data-generating DAG
    must contradict at least r_val closure orientations of G0.
"""

from __future__ import annotations

import collections
import json
import signal

from common import OUT, ROOT, dump, is_gac_valid_mpdag, k_of
from oset import optimal_set_closed_form
from release_rows import knowledge, rows, state

from bkrobust.demo.evaluate import optimal_adjustment_set_dag
from bkrobust.demo.meek import enumerate_dag_extensions
from bkrobust.robustness import llm_survival as ls
from bkrobust.robustness import real_survival as rs

ENUM_LIMIT = 3000


class TO(Exception):
    pass


def _h(*a):
    raise TO


signal.signal(signal.SIGALRM, _h)


def enum_o(g, x, y):
    """Enumeration O: 'disagree', a frozenset, or 'too_many'."""
    exts = enumerate_dag_extensions(g, limit=ENUM_LIMIT + 1)
    if len(exts) > ENUM_LIMIT:
        return "too_many"
    sets = {frozenset(optimal_adjustment_set_dag(d, x, y)) for d in exts}
    return next(iter(sets)) if len(sets) == 1 else "disagree"


def part_b_c():
    rs_ = rows(set(knowledge()))
    cb = collections.Counter()
    cc = collections.Counter()
    viol = []
    examples = []
    for r in rs_:
        dag, cpdag, k, g0, _ = state(r["arm"], r["network"])
        x, y = r["x"], r["y"]
        o, st = optimal_set_closed_form(g0, x, y)
        rel = r["status"]
        relmap = {"ok": "ok", "not_amenable": "not_amenable",
                  "y_not_possible_descendant": "y_not_descendant"}
        cb[f"status_match={relmap[rel] == st}"] += 1
        if rel == "ok" and st == "ok":
            cb[f"zsize_match={len(o) == int(r['z_size'])}"] += 1
            cb[f"gac_valid_at_g0={is_gac_valid_mpdag(g0, x, y, o)}"] += 1
            signal.alarm(10)
            try:
                e = enum_o(g0, x, y)
            except TO:
                e = "timeout"
            signal.alarm(0)
            if isinstance(e, frozenset):
                cb[f"enum_equal={e == o}"] += 1
                if e != o and len(examples) < 5:
                    examples.append({"net": r["network"], "x": x, "y": y,
                                     "closed": sorted(o), "enum": sorted(e)})
            else:
                cb[f"enum_{e}"] += 1
            # (c)
            if r["z_valid_at_truth"] == "0" and r["r_hop"] not in ("", "-1"):
                rr = int(r["r_hop"])
                n_contra = len(k_of(cpdag, g0) - set(dag.directed_edges))
                ok = n_contra >= rr
                cc[f"cor1c_holds={ok}"] += 1
                cc["margin_" + ("0" if n_contra == rr else ">0")] += 1
                if not ok:
                    viol.append({"arm": r["arm"], "net": r["network"], "x": x, "y": y,
                                 "r": rr, "contradicted": n_contra})
                truth_valid = is_gac_valid_mpdag(dag, x, y, o)
                cc[f"recheck_invalid_at_truth={not truth_valid}"] += 1
        elif st != relmap[rel] and len(examples) < 10:
            examples.append({"net": r["network"], "x": x, "y": y,
                             "release_status": rel, "closed_status": st})
    return dict(cb), dict(cc), viol, examples


def part_a():
    ss = json.loads((ROOT / "results/axis_robustness_llm_v2/status_split.json").read_text())
    real = set(ss["meta"]["real_conditions"])
    units = [u for u in ss["units"] if u["condition"] in real]
    bundle, _ = ls.load_bundle()
    cache = {}
    ca = collections.Counter()
    per_class = collections.defaultdict(collections.Counter)
    for u in units:
        net, cond, x, y = u["network"], u["condition"], u["x"], u["y"]
        key = (cond, net)
        if key not in cache:
            dag, cpdag = rs.load_network(net)
            k, g0, why, _ = ls.build_state(bundle, cond, net, dag, cpdag)
            cache[key] = (dag, cpdag, g0, why)
        dag, cpdag, g0, why = cache[key]
        if g0 is None:
            ca["g0_missing"] += 1
            continue
        o, st = optimal_set_closed_form(g0, x, y)
        cls = u["classification"]
        if st == "ok":
            v = is_gac_valid_mpdag(g0, x, y, o)
            lab = f"amenable_O_size{'0' if not o else '>0'}_valid={v}"
        else:
            lab = st
        per_class[cls][lab] += 1
        ca["n"] += 1
    return dict(ca), {k: dict(v) for k, v in per_class.items()}


if __name__ == "__main__":
    ca, pa = part_a()
    print("(a)", ca, json.dumps(pa, indent=1))
    cb, cc, viol, ex = part_b_c()
    print("(b)", cb)
    print("(c)", cc, "violations:", len(viol))
    print("examples:", json.dumps(ex[:8], indent=1))
    dump({"a_counts": ca, "a_by_class": pa, "b": cb, "c": cc, "c_violations": viol,
          "examples": ex}, OUT / "e13" / "consistency.json")
