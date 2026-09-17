"""
Validation of the new machinery against the pilot's already-validated machinery.

T-A: meek_closure_rules(G, ALL_RULES) == meek_closure(G) on 500 random CPDAGs
     with random partial orientations.
T-B: apply_background_knowledge_rules(C, K, ALL_RULES) agrees with
     apply_background_knowledge(C, K) -- same graph, same MeekFail verdict --
     on every knowledge set drawn below.
T-C: R1-only closure is a SUBGRAPH-of-orientations of the full closure
     (it can never orient an edge the full rule set leaves undirected).
T-D: the rule counter is non-degenerate -- REPORTED, and R1/R2/R4 asserted.
     NOTE the counter is FIRST-RULE-WINS (it credits the first rule in the
     R1->R2->R3->R4 order that can orient a given edge, exactly as the pilot's
     meek_closure evaluates them), so it UNDER-attributes R3 and R4: whenever R1
     also applies to the same edge, R1 takes the credit. The counter is therefore
     a lower bound per rule and is used only descriptively. The R1-sufficiency
     result in run_e2.py does NOT use the counter -- it compares the two closed
     GRAPHS for equality, which is attribution-free.
"""
import numpy as np

from graphs import (MeekFail, apply_background_knowledge, dag_to_cpdag,
                    meek_closure, random_dag, undirected_edges)
from meek_rules import (ALL_RULES, apply_background_knowledge_rules,
                        meek_closure_rules, orientation_status)


def main():
    rng = np.random.default_rng(20260814)
    n_a = n_b = n_c = 0
    fails_match = 0
    counter_total = {}
    for t in range(1500):
        p = int(rng.integers(5, 13))
        deg = float(rng.choice([1.5, 2.0, 2.5, 3.0, 4.0, 5.0]))
        D = random_dag(p, deg, rng)
        C = dag_to_cpdag(D)
        U = undirected_edges(C)
        if len(U) < 2:
            continue

        # T-A: random partial orientation, then closure both ways
        G = C.copy()
        for (i, j) in U:
            r = rng.random()
            if r < 0.3:
                G[j, i] = 0
            elif r < 0.6:
                G[i, j] = 0
        assert np.array_equal(meek_closure_rules(G, ALL_RULES), meek_closure(G)), \
            f"T-A mismatch at t={t}"
        n_a += 1

        # T-B / T-C: knowledge sets
        for _ in range(3):
            k = int(rng.integers(1, min(len(U), 6) + 1))
            idx = rng.choice(len(U), size=k, replace=False)
            K = []
            for i2 in idx:
                (i, j) = U[i2]
                K.append((i, j) if rng.random() < 0.5 else (j, i))
            try:
                Gref = apply_background_knowledge(C, K)
                ok_ref = True
            except MeekFail:
                Gref, ok_ref = None, False
            cnt = {}
            try:
                Gnew = apply_background_knowledge_rules(C, K, ALL_RULES, counter=cnt)
                ok_new = True
            except MeekFail:
                Gnew, ok_new = None, False
            assert ok_ref == ok_new, f"T-B verdict mismatch at t={t}"
            if ok_ref:
                assert np.array_equal(Gref, Gnew), f"T-B graph mismatch at t={t}"
                for kk, vv in cnt.items():
                    counter_total[kk] = counter_total.get(kk, 0) + vv
                # T-C: R1-only never orients more than the full rule set
                try:
                    G1 = apply_background_knowledge_rules(C, K, ("R1",))
                except MeekFail:
                    G1 = None
                if G1 is not None:
                    for (i, j) in U:
                        s1 = orientation_status(G1, i, j)
                        sf = orientation_status(Gnew, i, j)
                        assert s1 == 0 or s1 == sf, f"T-C violation at t={t} edge {(i,j)}"
                    n_c += 1
            else:
                fails_match += 1
            n_b += 1

    print(f"T-A closure equality   : {n_a} graphs        OK")
    print(f"T-B alg-1 equality     : {n_b} knowledge sets OK ({fails_match} matched MeekFails)")
    print(f"T-C R1 subset-of-full  : {n_c} MPDAGs         OK")
    print(f"T-D rule firings (first-rule-wins, LOWER BOUND per rule): {counter_total}")
    missing = [r for r in ("R1", "R2", "R4") if counter_total.get(r, 0) == 0]
    assert not missing, f"T-D: rules never fired -> selector untestable: {missing}"
    if counter_total.get("R3", 0) == 0:
        print("T-D NOTE: R3 credited 0 times -- expected under first-rule-wins "
              "attribution; R3's precondition is a superset of a configuration "
              "where R1 usually also applies. Does NOT affect R1-sufficiency, "
              "which is measured by graph equality.")
    print("ALL TESTS PASSED")


if __name__ == "__main__":
    main()
