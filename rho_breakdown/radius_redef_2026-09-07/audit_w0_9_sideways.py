import sys; sys.path.insert(0,'.')
from audit_w0_rk import *
from bkrobust.demo.space import enumerate_space, represented_dags, covering_pairs, neighbour_graph, bfs_distances
rng=np.random.default_rng(8080); probs,_=gen(rng,150,n=7,p=0.40,k_claims=3,n_false=0)
side=0; side_bad=0; tot=0; pws=0
for pr in probs:
    sp=enumerate_space(pr['cpdag']); reps=represented_dags(sp)
    cov=covering_pairs(sp,reps); nb=neighbour_graph(sp,cov)
    g0=next((g for g in sp if g==pr['g0']),None)
    if g0 is None: continue
    dist=bfs_distances(nb,g0)
    D0=set(enumerate_dag_extensions(g0)); hs=False
    for g in sp:
        if g==g0 or g not in dist: continue
        tot+=1; Dg=set(enumerate_dag_extensions(g))
        if not (Dg<=D0 or D0<=Dg):
            side+=1; hs=True
            if not is_valid_adjustment_set_mpdag(g,pr['x'],pr['y'],pr['z']): side_bad+=1
    pws+=hs
print(f"non-G0 reachable elements examined        : {tot}")
print(f"  SIDEWAYS (incomparable to [G0])         : {side}  ({side/max(1,tot):.3f})")
print(f"    of those, Z INVALID                   : {side_bad}")
print(f"  problems containing >=1 sideways element: {pws}/{len(probs)}")
