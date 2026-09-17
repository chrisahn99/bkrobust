import sys, itertools, numpy as np
sys.path.insert(0,'/private/tmp/claude-501/-Users-josecosta-mugango/eee5422c-2e7e-4a55-8d52-0664260abb06/scratchpad/ref')
sys.path.insert(0,'/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/code')
import indep as I
from graphs import random_dag
from adjust import possibly_causal_paths

CASES=[("licensed",10442,7,2.0,3,1,(0,2)),
       ("large",1847,23,2.0,9,21,(6,22)),
       ("licensed",2650,9,2.0,4,2,(8,3)),
       ("large",3610,20,2.5,15,18,(17,4)),
       ("large",4802,20,2.0,3,10,(11,6)),
       ("licensed",11746,9,2.0,3,5,(2,6))]
for ens,seed,p,deg,x,y,stmt in CASES:
    rng=np.random.default_rng(seed); Dm=random_dag(p,deg,rng)
    dag=[(i,j) for i in range(p) for j in range(p) if Dm[i,j]==1 and Dm[j,i]==0]
    Dg=I.G(p); Dg.d=set(dag)
    C=I.dag_to_cpdag(dag,p)
    pairs=[(a,b) for a in range(p) for b in range(p) if a!=b and len(possibly_causal_paths(Dm,a,b))>0]
    xq,yq=pairs[int(rng.integers(len(pairs)))]
    dd=I.dist(C,[x,y])
    u,v=stmt
    print('=== %s seed=%d p=%d deg=%.1f  query (%d,%d) [regen query %s] stmt %d->%d'%(ens,seed,p,deg,x,y,(xq,yq),u,v))
    print('   delta(u)=%s delta(v)=%s -> dmin=%s dmax=%s ; true in D? %s ; edge undirected in C? %s'%(
        dd[u],dd[v],min(dd[u],dd[v]),max(dd[u],dd[v]),(u,v) in Dg.d, C.isund(u,v)))
    # search over base K sets of true statements that make G0 amenable and leave (u,v) open
    Uedges=[tuple(sorted(e)) for e in C.u]
    Kall=[(a,b) if (a,b) in Dg.d else (b,a) for (a,b) in Uedges]
    found=None
    rng2=np.random.default_rng(0)
    for trial in range(4000):
        k=int(rng2.integers(1,5))
        idx=rng2.permutation(len(Kall))[:k]
        Kb=[Kall[i] for i in idx]
        if any(tuple(sorted(s))==tuple(sorted((u,v))) for s in Kb): continue
        try: G0=I.apply_K(C,Kb)
        except I.Inconsistent: continue
        if not G0.isund(u,v): continue
        if not I.amenable(G0,x,y): continue
        try: G1=I.apply_K(C,Kb+[(u,v)])
        except I.Inconsistent: continue
        npaths=len(I.pcp(G1,x,y))
        if (not I.amenable(G1,x,y)) and npaths>0:
            found=(Kb,G0,G1,npaths); break
    if found:
        Kb,G0,G1,npaths=found
        print('   FOUND base K =',Kb)
        print('   G0 amenable=%s  O*(G0)=%s'%(I.amenable(G0,x,y),sorted(I.Ostar(G0,x,y))))
        print('   after asserting %d->%d : amenable=%s, n possibly-causal paths=%d  => IDENTIFICATION LOST'%(u,v,I.amenable(G1,x,y),npaths))
        pc=I.pcp(G1,x,y)
        bad=[q for q in pc if not G1.isdir(q[0],q[1])]
        print('   offending path(s) not starting with a directed edge out of X:',bad[:3])
        print('   edges changed by the assertion:', sorted(set(G1.d)-set(G0.d)))
    else:
        print('   (no base K found by random search in 4000 tries)')
