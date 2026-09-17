import sys
sys.path.insert(0, "/Users/josecosta/bkrobust/src")
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-09-07_radius-redef")
from bkrobust.demo.scenario import true_dag, scenarios, TREATMENT, OUTCOME
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag
from final_inv2_L import partition
dag = true_dag(); cpdag = dag_to_cpdag(dag); x, y = TREATMENT, OUTCOME
for lab, sc in scenarios().items():
    K = [tuple(e) for e in sc["knowledge"]]
    g0 = apply_orientations(cpdag, K)
    z = optimal_adjustment_set_mpdag(g0, x, y)
    L, Fr, I, D = partition(cpdag, K, x, y, frozenset(z), g0)
    P = [k for k in K if k[0] in (x, y) or k[1] in (x, y)]
    print(f"{lab}: K={K}  Z={sorted(z)}")
    print(f"    L={L}  free={Fr}  inert={I}  data-refuted(Chat-alone)={D}  P(hop0)={P}")
