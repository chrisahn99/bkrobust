#!/usr/bin/env python3
"""inverse-2 audit, part 2: is 'data-refuted' really refuted BY THE DATA?

The definition's second exact statement says a data-refuted claim 'cannot be wrong in
any direction consistent with the true equivalence class'.  But G_{k-bar} is computed
from (K \ {k}) u {reverse(k)} -- it conditions on the REST OF THE ELICITED CLAIMS,
which may themselves be false.  Test: does the CPDAG ALONE forbid the reversal?
"""
import sys
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
from bkrobust.demo.scenario import true_dag, scenarios, TREATMENT, OUTCOME
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag, is_valid_adjustment_set_mpdag

dag = true_dag(); cpdag = dag_to_cpdag(dag); X, Y = TREATMENT, OUTCOME
for lab, sc in scenarios().items():
    K = [tuple(e) for e in sc["knowledge"]]
    g0 = apply_orientations(cpdag, K); z = frozenset(optimal_adjustment_set_mpdag(g0, X, Y))
    print(f"\n=== scenario {lab}")
    for k in K:
        rest = [e for e in K if e != k]
        g_rev_cond = apply_orientations(cpdag, rest + [(k[1], k[0])])   # as defined
        g_rev_alone = apply_orientations(cpdag, [(k[1], k[0])])         # data alone
        if g_rev_cond is None:
            verdict = "flagged DATA-REFUTED"
            if g_rev_alone is None:
                truth = "  and the CPDAG alone really does forbid it   -> honest"
            else:
                bad = not is_valid_adjustment_set_mpdag(g_rev_alone, X, Y, z)
                truth = (f"  BUT the CPDAG ALONE PERMITS the reversal "
                         f"(Z invalid there? {bad})   -> label is FALSE")
            print(f"  {str(k):22s} {verdict}{truth}")
