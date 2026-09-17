import sys, itertools, time
sys.path.insert(0,"/Users/josecosta/bkrobust/src"); sys.path.insert(0,"/Users/josecosta/mugango/output/2026-09-07_radius-redef")
import numpy as np
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.space import enumerate_space, represented_dags, covering_pairs, neighbour_graph, bfs_distances
from bkrobust.demo.evaluate import optimal_adjustment_set_dag, random_sem, adjusted_estimand
from harness import random_dag

def tau_range(g,sem,x,y,c):
    k=g.edge_string()
    if k in c: return c[k]
    v=[]
    for D in enumerate_dag_extensions(g):
        try: v.append(adjusted_estimand(sem,x,y,optimal_adjustment_set_dag(D,x,y)))
        except Exception: pass
    c[k]=(min(v),max(v)) if v else None; return c[k]

def rhig(cp,K,j,sem,x,y,c):
    m=len(K); lo,hi=np.inf,-np.inf; ok=False
    for S in itertools.combinations(range(m),m-j):
        g=apply_orientations(cp,[K[i] for i in S])
        if g is None: continue
        r=tau_range(g,sem,x,y,c)
        if r is None: continue
        ok=True; lo=min(lo,r[0]); hi=max(hi,r[1])
    return (lo,hi) if ok else None

rng=np.random.default_rng(11); pairs=[]; bigger=0; nbig=0; t0=time.time(); tried=0
while len(pairs)<150 and tried<80000 and time.time()-t0<300:
    tried+=1
    dag=random_dag(rng,7,0.32)
    if not dag.directed_edges: continue
    cp=dag_to_cpdag(dag); und=sorted(tuple(sorted(e)) for e in cp.undirected_edges)
    if not (3<=len(und)<=6): continue
    nodes=sorted(dag.nodes); cand=[(x,y) for x in nodes for y in nodes if x!=y]; rng.shuffle(cand)
    sem=random_sem(dag,rng); c={}
    for x,y in cand[:6]:
        b=tau_range(cp,sem,x,y,c)
        if b is None or b[1]-b[0]<1e-8: continue
        # truthful orientations of all undirected edges, in a fixed random order
        perm=list(rng.permutation(len(und)))
        full=[]
        for t in perm:
            a,bb=und[t]; full.append((a,bb) if (a,bb) in dag.directed_edges else (bb,a))
        w={}
        for m in (1,2,3):
            if m>len(full): continue
            K=full[:m]
            if apply_orientations(cp,K) is None: continue
            r=rhig(cp,K,1,sem,x,y,c)
            if r: w[m]=r[1]-r[0]
        if len(w)>=2:
            pairs.append((w,b[1]-b[0]))
        # BIG_1 vs RHIG_1 on m=3
        if 3 in w:
            K=full[:3]; g0=apply_orientations(cp,K)
            sp=enumerate_space(cp); rp=represented_dags(sp); cov=covering_pairs(sp,rp); nb=neighbour_graph(sp,cov)
            g0s=next((g for g in sp if g==g0),None)
            if g0s is not None:
                d=bfs_distances(nb,g0s); lo,hi=np.inf,-np.inf
                for g,dd in d.items():
                    if dd is not None and dd<=1:
                        r=tau_range(g,sem,x,y,c)
                        if r: lo=min(lo,r[0]); hi=max(hi,r[1])
                if lo<np.inf:
                    nbig+=1
                    if hi-lo > w[3]+1e-9: bigger+=1
        break
print(f"n={len(pairs)} tried={tried} {time.time()-t0:.1f}s")
for a,bb in ((1,2),(2,3),(1,3)):
    s=[(p[0][a],p[0][bb]) for p in pairs if a in p[0] and bb in p[0]]
    if not s: continue
    A=np.array([u for u,_ in s]); B=np.array([v for _,v in s])
    print(f"RHIG_1 with m={a} vs m={bb} claims, SAME graph+query, all claims TRUE (n={len(s)}):")
    print(f"   mean width {A.mean():.4f} -> {B.mean():.4f} | wider with more knowledge {100*np.mean(B>A+1e-9):.1f}% | narrower {100*np.mean(B<A-1e-9):.1f}% | equal {100*np.mean(np.abs(B-A)<1e-9):.1f}%")
print(f"\nBIG_1 strictly wider than RHIG_1 (m=3): {bigger}/{nbig}")
