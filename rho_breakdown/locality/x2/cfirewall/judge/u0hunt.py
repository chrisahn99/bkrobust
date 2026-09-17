"""Directed hunt for a counterexample to LEMMA U0 (the named gap for clause 1):
      pa_H(cn_H) n poss_de_{G0}(X)  subset of  cn_H u {X}
   plus the derived clause-(1) statement  O*(H) n de_D(X) = {}."""
import sys, json, time
from collections import defaultdict
from multiprocessing import Pool
import numpy as np
from jlib import *

def one(arg):
    seed,p,deg=arg
    rng=np.random.default_rng(seed); acc=defaultdict(int); wit=[]
    D0=random_dag(p,deg,rng); C=dag_to_cpdag(D0); ref=v_structures(C); G0=C.copy()
    for _ in range(int(rng.integers(0,5))):
        UU=undirected_edges(G0)
        if not UU: break
        u,v=UU[int(rng.integers(len(UU)))]
        a,b=(u,v) if rng.random()<0.5 else (v,u)
        T=G0.copy(); T[b,a]=0; T=meek_closure(T)
        if has_directed_cycle(T) or v_structures(T)!=ref: continue
        G0=T
    NA=nonadjacent_pairs(G0)
    if not NA: return dict(acc),wit
    Q=[]
    for x in range(p):
        for y in range(p):
            if x==y: continue
            O0,n0,am0,cn0,fb0=ostar(G0,x,y)
            if am0: Q.append((x,y,O0,poss_de(G0,{x})))
    if not Q: return dict(acc),wit
    for (a0,b0) in NA:
        for (a,b) in ((a0,b0),(b0,a0)):
            H,info=bk_assert(G0,[(a,b)])
            if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
            for (x,y,O0,pdX0) in Q:
                O1,n1,am1,cn1,fb1=ostar(H,x,y)
                if not am1: continue
                acc["report"]+=1
                moved=(O1!=O0); acc["moved"]+=moved
                bad=(parents_of_set(H,cn1)&pdX0)-(cn1|{x})
                acc["U0_nonvac"]+= bool(parents_of_set(H,cn1)&pdX0)
                if bad:
                    acc["U0_FAIL"]+=1
                    if len(wit)<5: wit.append(dict(kind="U0",G0=G0.tolist(),H=H.tolist(),x=x,y=y,a=a,b=b,bad=sorted(bad)))
                if not (pdX0<=poss_de(H,{x})):
                    acc["Mfail"]+=1; acc["Mfail_moved"]+=moved
                    if moved and len(wit)<8: wit.append(dict(kind="Mmoved",G0=G0.tolist(),H=H.tolist(),x=x,y=y,a=a,b=b))
                # derived clause-(1) statement against the truth
                if O1 & pdX0:
                    acc["C1_STRONG_FAIL"]+=1
                    if len(wit)<10: wit.append(dict(kind="C1",G0=G0.tolist(),H=H.tolist(),x=x,y=y,a=a,b=b,O1=sorted(O1)))
    return dict(acc),wit

def main():
    n=int(sys.argv[1]); nw=int(sys.argv[2]); s0=int(sys.argv[3])
    jobs=[(s0+i,int(5+(i%7)),float([1.5,2.0,2.5,3.0,4.0][i%5])) for i in range(n)]
    acc=defaultdict(int); W=[]; t0=time.time()
    with Pool(nw) as pool:
        for a,w in pool.imap_unordered(one,jobs,chunksize=16):
            for k,v in a.items(): acc[k]+=v
            W.extend(w[:1])
    acc["elapsed"]=round(time.time()-t0,1)
    print(json.dumps(dict(sorted(acc.items())),indent=1))
    if W: json.dump(W[:20],open("U0_WIT.json","w"),indent=1); print("!!! WITNESS",len(W))
if __name__=="__main__": main()
