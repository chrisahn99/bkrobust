"""ADVERSARIAL search for a counterexample to CONJECTURE S1:
  G0 MPDAG, D in [G0], R(G0,x,y) != REFUSE (identified & sound);
  H = bk_assert(G0, [a->b]) with (a,b) NON-adjacent, H acyclic+extendable+no conflict
  ==>  R(H,x,y) sound ?
Adversarial bias: dense graphs, and spurious statements aimed AT cn / near x,y.
Also: rho up to 4 statements, and the same question for state ZERO."""
import sys
from collections import defaultdict
sys.path.insert(0,"/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/falsegain/code")
from fglib import *
import fglib, numpy as np, importlib.util
from graphs import random_dag, is_undirected
from adjust import causal_nodes
spec=importlib.util.spec_from_file_location("x1_ops","/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x1/code/x1_ops.py")
x1_ops=importlib.util.module_from_spec(spec); spec.loader.exec_module(x1_ops)

def true_bk(C,D,rng):
    G=C.copy(); U=undirected_edges(G)
    if not U: return G
    for t in rng.permutation(len(U))[:int(rng.integers(0,len(U)+1))]:
        u,v=U[t]
        if not is_undirected(G,u,v): continue
        if D[u,v]==1: G[v,u]=0
        else: G[u,v]=0
        G=meek_closure(G)
    return G

def main(n,seed0):
    rng=np.random.default_rng(seed0); acc=defaultdict(int); bad=[]
    for it in range(n):
        p=int(rng.choice([5,6,7,8,9,10])); deg=float(rng.choice([1.5,2.0,2.5,3.0,4.0,5.0]))
        D=random_dag(p,deg,rng); C=dag_to_cpdag(D); G0=true_bk(C,D,rng)
        NA=x1_ops.nonadjacent_pairs(C)
        if not NA: continue
        qs=[(x,y) for x in range(p) for y in range(p) if x!=y]
        rng.shuffle(qs)
        for (x,y) in qs[:14]:
            s0,O0=report_state(G0,x,y)
            if s0==REFUSE: continue
            n0=sound(D,x,y,s0,O0)
            acc[("base",s0,n0)]+=1
            cnset=set(causal_nodes(G0,x,y)) | {x,y}
            # adversarial pool: at least one endpoint in cn u {x,y}
            adv=[pr for pr in NA if pr[0] in cnset or pr[1] in cnset]
            pool = adv if adv else NA
            for rho in (1,2,3,4):
                for rep in range(3):
                    if len(pool)<rho: continue
                    idx=rng.choice(len(pool),size=rho,replace=False)
                    K=[]
                    for t in idx:
                        a,b=pool[int(t)]
                        K.append((int(a),int(b)) if rng.random()<0.5 else (int(b),int(a)))
                    H,info=x1_ops.bk_assert(G0,K)
                    ok=(not info["conflict"]) and (not has_directed_cycle(H)) and x1_ops.pdag_extendable(H)
                    s1,O1=report_state(H,x,y); n1=sound(D,x,y,s1,O1)
                    acc[("S1" if ok else "S1rej",rho,s0,s1,n1)]+=1
                    if ok and n1 is False and len(bad)<6:
                        bad.append(dict(p=p,D=D.tolist(),G0=G0.tolist(),H=H.tolist(),x=x,y=y,K=K,s0=s0,s1=s1))
        if len(fglib._valid_cache)>200000: fglib._valid_cache.clear(); fglib._cn_cache.clear()
        if it%100==0: print("  it",it,flush=True)
    return acc,bad

if __name__=="__main__":
    acc,bad=main(int(sys.argv[1]),int(sys.argv[2]))
    for k,v in sorted(acc.items(),key=lambda z:str(z[0])): print("|".join(map(str,k)),v)
    print("COUNTEREXAMPLES", len(bad))
    import json; json.dump(bad,open("/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/falsegain/code/adv_bad.json","w"),indent=1)
