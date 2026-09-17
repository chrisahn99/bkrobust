import sys
sys.path.insert(0,"/Users/josecosta/bkrobust/src"); sys.path.insert(0,"/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from harness import random_dag
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import enumerate_dag_extensions, is_consistent_extension
from bkrobust.demo.evaluate import (optimal_adjustment_set_dag, adjusted_estimand,
                                    is_valid_adjustment_set_dag, random_sem)
rng = np.random.default_rng(7)
bad=0; tot=0; badext=0
for _ in range(60):
    dag = random_dag(rng,6,0.35)
    if not dag.directed_edges: continue
    cp = dag_to_cpdag(dag)
    exts = enumerate_dag_extensions(cp)
    inx = any(set(d.directed_edges)==set(dag.directed_edges) for d in exts)
    if not inx: badext+=1
    sem = random_sem(dag,rng)
    for x in sorted(dag.nodes):
        for y in sorted(dag.nodes):
            if x==y: continue
            tot+=1
            zs = optimal_adjustment_set_dag(dag,x,y)
            tau = sem.true_total_effect(x,y)
            th = adjusted_estimand(sem,x,y,zs)
            if abs(th-tau) > 1e-8:
                bad+=1
                if bad<6:
                    print(f"MISMATCH x={x} y={y} tau={tau:.4f} theta={th:.4f} O*={sorted(zs)} "
                          f"valid={is_valid_adjustment_set_dag(dag,x,y,zs)} "
                          f"desc={y in dag.descendants(x)}")
print(f"true DAG missing from extensions of its own CPDAG: {badext}/60")
print(f"O*(true dag) fails to identify tau on {bad}/{tot} queries")
