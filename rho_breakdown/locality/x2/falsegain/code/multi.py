"""RHO >= 2: does the immunity of an already-identified query survive MANY false statements?
Base G0 = M(C,K_true).  Then K' adds rho statements of the chosen class.
Class R : Meek semantics (apply_background_knowledge on C with K_true+K_false); MeekFail -> caught.
Class S : bk_assert with rho spurious required edges; S1 = extendable & acyclic & no conflict.
"""
import sys
from collections import defaultdict
sys.path.insert(0,"/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/falsegain/code")
from fglib import *
import fglib, numpy as np, importlib.util
from graphs import random_dag, is_undirected, apply_background_knowledge, MeekFail
spec=importlib.util.spec_from_file_location("x1_ops","/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x1/code/x1_ops.py")
x1_ops=importlib.util.module_from_spec(spec); spec.loader.exec_module(x1_ops)

def main(n, seed0):
    rng=np.random.default_rng(seed0); acc=defaultdict(int); wit=[]
    for it in range(n):
        p=int(rng.choice([5,6,7,8])); deg=float(rng.choice([1.5,2.0,2.5,3.0]))
        D=random_dag(p,deg,rng); C=dag_to_cpdag(D); U=undirected_edges(C)
        if len(U)<3: continue
        Ktrue=[(u,v) if D[u,v]==1 else (v,u) for (u,v) in U]
        for rho in (1,2,3):
            for rep in range(4):
                nt=int(rng.integers(0,max(len(U)-rho,0)+1))
                perm=rng.permutation(len(U))
                base_idx=list(perm[:nt]); flip_idx=list(perm[nt:nt+rho])
                if len(flip_idx)<rho: continue
                Kb=[Ktrue[i] for i in base_idx]
                Kf=[(Ktrue[i][1],Ktrue[i][0]) for i in flip_idx]
                try: G0=apply_background_knowledge(C,Kb)
                except MeekFail: continue
                # ---- class R, Meek semantics
                try:
                    H=apply_background_knowledge(C,Kb+Kf); ok=True
                except MeekFail:
                    ok=False
                acc[("R_meek",rho,"accepted" if ok else "CAUGHT")]+=1
                if ok:
                    for x in range(p):
                        for y in range(p):
                            if x==y: continue
                            s0,O0=report_state(G0,x,y); n0=sound(D,x,y,s0,O0)
                            s1,O1=report_state(H,x,y); n1=sound(D,x,y,s1,O1)
                            acc[("R",rho,s0,n0,s1,n1)]+=1
                            if s0!=REFUSE and n1 is False and len(wit)<5:
                                wit.append(("R",rho,p,D.tolist(),G0.tolist(),H.tolist(),x,y,Kb,Kf))
                # ---- class S
                NA=x1_ops.nonadjacent_pairs(C)
                if len(NA)>=rho:
                    idx=rng.choice(len(NA),size=rho,replace=False)
                    Ks=[]
                    for t in idx:
                        a,b=NA[int(t)]
                        Ks.append((int(a),int(b)) if rng.random()<0.5 else (int(b),int(a)))
                    HS,info=x1_ops.bk_assert(G0,Ks)
                    okS=(not info["conflict"]) and (not has_directed_cycle(HS)) and x1_ops.pdag_extendable(HS)
                    acc[("S_ext",rho,"S1" if okS else "S1rej")]+=1
                    for x in range(p):
                        for y in range(p):
                            if x==y: continue
                            s0,O0=report_state(G0,x,y); n0=sound(D,x,y,s0,O0)
                            s1,O1=report_state(HS,x,y); n1=sound(D,x,y,s1,O1)
                            acc[("S1" if okS else "S1rej",rho,s0,n0,s1,n1)]+=1
                            if okS and s0!=REFUSE and n1 is False and len(wit)<5:
                                wit.append(("S1",rho,p,D.tolist(),G0.tolist(),HS.tolist(),x,y,Kb,Ks))
        if len(fglib._valid_cache)>200000: fglib._valid_cache.clear(); fglib._cn_cache.clear()
        if it%100==0: print("  it",it,flush=True)
    return acc,wit

if __name__=="__main__":
    acc,wit=main(int(sys.argv[1]),int(sys.argv[2]))
    for k,v in sorted(acc.items(),key=lambda z:str(z[0])): print("|".join(map(str,k)),v)
    print("N_WITNESSES_POS_CORRUPTED", len(wit))
    import json; json.dump(wit,open("/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/falsegain/code/multi_wit.json","w"),indent=1,default=str)
