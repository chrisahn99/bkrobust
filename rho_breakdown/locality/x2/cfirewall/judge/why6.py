"""Why-distribution of the certificate-uncovered set at p=6..9, and whether
   'D+(a->b) acyclic'  =>  covered by certA v certL v certC."""
import sys, json
from collections import defaultdict
from multiprocessing import Pool
import numpy as np
from jlib import *

def one(arg):
    seed,p,deg=arg
    rng=np.random.default_rng(seed); acc=defaultdict(int); wit=[]
    D0=random_dag(p,deg,rng); C=dag_to_cpdag(D0); ref=v_structures(C); G0=C.copy()
    for _ in range(int(rng.integers(0,4))):
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
            if am0: Q.append((x,y,O0))
    if not Q: return dict(acc),wit
    ext=consistent_dag_extensions(G0,ref_vstructs=ref,limit=120)
    if not ext: return dict(acc),wit
    DQ={}
    for di,D in enumerate(ext):
        for (x,y,_) in Q:
            if (di,x,y) in DQ: continue
            cnD=causal_nodes(D,x,y); fbD=(poss_de(D,cnD)|{x}) if cnD else {x}
            R=reach_outside_forb(D,x,y,cnD,fbD) if cnD else set()
            need=((parents_of_set(D,cnD)-fbD)&R) if cnD else None
            DQ[(di,x,y)]=(cnD,fbD,need,{u for u in range(p) if is_directed(D,u,x)},poss_de(D,{x}))
    for (a0,b0) in NA:
        for (a,b) in ((a0,b0),(b0,a0)):
            H,info=bk_assert(G0,[(a,b)])
            if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
            liftc={}
            for (x,y,O0) in Q:
                O1,n1,am1,cn1,fb1=ostar(H,x,y)
                if not am1: continue
                Z=set(O1)
                for di,D in enumerate(ext):
                    if di not in liftc: liftc[di]=lifts(D,a,b,H)[0::2]
                    ok,why=liftc[di]
                    cnD,fbD,need,paX,deX=DQ[(di,x,y)]
                    acc["n"]+=1
                    if ok: continue
                    if (need is not None) and (not(Z&fbD)) and (need<=Z): continue
                    if (paX<=Z) and not (Z&deX): continue
                    acc["OPEN"]+=1; acc["open_"+why]+=1
                    v=is_valid_adjustment_set(D,x,y,Z)
                    acc["OPEN_INVALID"]+=(not v)
                    if why!="cycle":
                        acc["OPEN_ACYCLIC"]+=1
                        if len(wit)<4: wit.append(dict(G0=G0.tolist(),H=H.tolist(),D=D.tolist(),x=x,y=y,a=a,b=b,O0=sorted(O0),O1=sorted(O1),why=why,valid=bool(v)))
    return dict(acc),wit

if __name__=="__main__":
    n=int(sys.argv[1]); nw=int(sys.argv[2]); s0=int(sys.argv[3])
    jobs=[(s0+i,int(6+(i%4)),float([2.0,2.5,3.0][i%3])) for i in range(n)]
    acc=defaultdict(int); W=[]
    with Pool(nw) as pool:
        for a,w in pool.imap_unordered(one,jobs,chunksize=8):
            for k,v in a.items(): acc[k]+=v
            W.extend(w[:2])
    print(json.dumps(dict(sorted(acc.items())),indent=1))
    if W: json.dump(W[:20],open("OPEN_ACYCLIC.json","w"),indent=1); print("acyclic-open witnesses:",len(W))
