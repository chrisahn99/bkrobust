"""
Machinery tests for the C-FIREWALL census.  A census that measures the wrong
object is worse than no census, so every assumption the scan rests on is checked
by brute force at p <= 4 (and p = 5 for the cheap ones).

T1  bk_assert on a NON-ADJACENT pair == "add the edge, then Meek-close"
    (the re-stamp is a no-op for a single statement) and it really grows the
    skeleton, so skeleton(H) != skeleton(D) for every D in [G0].
T2  [G0] as computed (consistent_dag_extensions(G0, v_structures(C))) equals the
    brute-force set { D in MEC(C) : D agrees with every directed edge of G0 },
    where MEC(C) is itself brute-forced from all labelled DAGs.
T3  [G0] is non-empty for every reachable MPDAG, and G0 is exactly the common
    orientation graph of [G0]  (i.e. reachable_mpdags really enumerates MPDAGs).
T4  ostar_and_paths == adjust.optimal_adjustment_set on every reachable MPDAG.
T5  BASELINE / oracle power: O*(X,Y,G0) IS valid in every D in [G0]
    (Henckel-Perkovic-Maathuis class invariance) -- the validity checker returns
    True where it must.
T6  ORACLE POWER against the null: the SAME validity predicate, applied to
    O*(X,Y,H) computed from graphs H obtained by asserting a->b and STAMPING IT
    ANYWAY when a coherence check FAILS, does fire.  Without T6 a 0/N result is
    consistent with a checker that never returns False.
"""
import sys
import importlib.util
from itertools import product

import numpy as np

BASE = "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious"
sys.path.insert(0, BASE + "/x2/code")
from graphs import (common_orientation_graph, consistent_dag_extensions,
                    dag_agrees_with, dag_to_cpdag, directed_edges,
                    has_directed_cycle, meek_closure, skeleton, v_structures)
import adjust
import x2lib as X
from run_lemma import all_dags, reachable_mpdags

_spec = importlib.util.spec_from_file_location("x1_ops", BASE + "/x1/code/x1_ops.py")
x1_ops = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(x1_ops)

fails = []


def chk(name, ok, extra=""):
    print(("PASS " if ok else "FAIL ") + name + ("" if ok else "  " + extra))
    if not ok:
        fails.append(name)


def run(p):
    dags = list(all_dags(p))
    cpdags = {}
    for D in dags:
        C = dag_to_cpdag(D)
        cpdags.setdefault(C.tobytes(), C)
    # brute-force MEC index
    mec = {}
    for D in dags:
        mec.setdefault(dag_to_cpdag(D).tobytes(), []).append(D)

    t1 = t2 = t3 = t4 = t5 = True
    t6_fired = 0
    t6_tot = 0
    n_mpdag = 0
    for kb, C in cpdags.items():
        ref = v_structures(C)
        for G0 in reachable_mpdags(C):
            n_mpdag += 1
            ext = consistent_dag_extensions(G0, ref_vstructs=ref)
            brute = [D for D in mec[kb] if dag_agrees_with(D, G0)]
            if sorted(D.tobytes() for D in ext) != sorted(D.tobytes() for D in brute):
                t2 = False
            if not ext:
                t3 = False
            elif not np.array_equal(common_orientation_graph(ext, skeleton(G0)), G0):
                t3 = False
            for (a0, b0) in x1_ops.nonadjacent_pairs(G0):
                for (a, b) in ((a0, b0), (b0, a0)):
                    H, info = x1_ops.bk_assert(G0, [(a, b)])
                    Gm = G0.copy()
                    Gm[a, b] = 1
                    Gm[b, a] = 0
                    if not np.array_equal(H, meek_closure(Gm)):
                        t1 = False
                    if np.array_equal(skeleton(H), skeleton(G0)):
                        t1 = False
                    ok = ((not info["conflict"]) and (not has_directed_cycle(H))
                          and x1_ops.pdag_extendable(H))
                    if ok:
                        continue
                    # ---- T6: the control arm -- stamp anyway, same predicate
                    for x in range(p):
                        for y in range(p):
                            if x == y:
                                continue
                            O1, np1, am1 = X.ostar_and_paths(H, x, y)
                            if not am1:
                                continue
                            for D in ext:
                                t6_tot += 1
                                if not adjust.is_valid_adjustment_set(D, x, y, set(O1)):
                                    t6_fired += 1
            for x in range(p):
                for y in range(p):
                    if x == y:
                        continue
                    O, npa, am = X.ostar_and_paths(G0, x, y)
                    lib = adjust.optimal_adjustment_set(G0, x, y)
                    if (O if am else None) != lib:
                        t4 = False
                    if am:
                        for D in ext:
                            if not adjust.is_valid_adjustment_set(D, x, y, set(O)):
                                t5 = False
    print(f"--- p={p}: {len(cpdags)} cpdags, {n_mpdag} mpdags")
    chk(f"T1 p={p} bk_assert == add+meek_closure, skeleton grows", t1)
    chk(f"T2 p={p} [G0] == brute-force MEC members agreeing with G0", t2)
    chk(f"T3 p={p} [G0] non-empty and G0 == common orientation graph", t3)
    chk(f"T4 p={p} ostar_and_paths == optimal_adjustment_set", t4)
    chk(f"T5 p={p} O*(G0) valid in every D of [G0] (HPM class invariance)", t5)
    chk(f"T6 p={p} validity predicate FIRES in the incoherent control arm "
        f"({t6_fired}/{t6_tot})", t6_fired > 0, "checker may be vacuous")


if __name__ == "__main__":
    for p in (3, 4):
        run(p)
    print("\nFAILURES:", fails if fails else "none")
    sys.exit(1 if fails else 0)
