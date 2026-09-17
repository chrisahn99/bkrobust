"""Adversarial search CONCENTRATED on the certificate-uncovered residual.
Given certA/certL/certC are sound, ANY counterexample to C-FIREWALL must lie here."""
import sys, json, time
from collections import defaultdict
from multiprocessing import Pool
import numpy as np
from jlib import *

def one(arg):
    seed,p,deg,fam=arg
    rng=np.random.default_rng(seed); acc=defaultdict(int); wit=[]
    if fam=="er": D0=random_dag(p,deg,rng)
    elif fam=="chain":                      # long chordless path + a few chords
        order=rng.permutation(p); D0=np.zeros((p,p),np.int8)
        for i in range(p-1): D0[order[i],order[i+1]]=1
        for _ in range(int(rng.integers(0,3))):
            i,j=sorted(rng.choice(p,2,replace=False))
            if j>i+1: D0[order[i],order[j]]=1
    elif fam=="star":
        order=rng.permutation(p); D0=np.zeros((p,p),np.int8)
        h=order[0]
        for i in range(1,p):
            if rng.random()<0.5: D0[h,order[i]]=1
            else: D0[order[i],h]=1
        for _ in range(int(rng.integers(0,p))):
            i,j=sorted(rng.choice(range(1,p),2,replace=False))
            D0[order[i],order[j]]=1
    else: D0=random_dag(p,deg,rng)
    C=dag_to_cpdag(D0); ref=v_structures(C)
    U=undirected_edges(C); G0=C.copy()
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
    ext=consistent_dag_extensions(G0,ref_vstructs=ref,limit=200)
    if not ext: return dict(acc),wit
    acc["G0"]=1
    DQ={}
    for di,D in enumerate(ext):
        for (x,y,_) in Q:
            if (di,x,y) in DQ: continue
            cnD=causal_nodes(D,x,y); fbD=(poss_de(D,cnD)|{x}) if cnD else {x}
            R=reach_outside_forb(D,x,y,cnD,fbD) if cnD else set()
            need=((parents_of_set(D,cnD)-fbD)&R) if cnD else None
            paX={u for u in range(p) if is_directed(D,u,x)}
            DQ[(di,x,y)]=(fbD,need,paX,poss_de(D,{x}))
    for (a0,b0) in NA:
        for (a,b) in ((a0,b0),(b0,a0)):
            H,info=bk_assert(G0,[(a,b)])
            if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
            liftc={}
            for (x,y,O0) in Q:
                O1,n1,am1,cn1,fb1=ostar(H,x,y)
                if not am1: continue
                acc["report"]+=1; S1=set(O1)
                for di,D in enumerate(ext):
                    if di not in liftc: liftc[di]=lifts(D,a,b,H)[0]
                    fbD,need,paX,deX=DQ[(di,x,y)]
                    acc["n"]+=1
                    if liftc[di]: continue
                    if (need is not None) and (not(S1&fbD)) and (need<=S1): continue
                    if (paX<=S1) and not (S1&deX): continue
                    acc["UNCOV"]+=1
                    acc["UNCOV_moved"]+= (O1!=O0)
                    if not is_valid_adjustment_set(D,x,y,S1):
                        acc["UNCOV_INVALID"]+=1
                        if len(wit)<5: wit.append(dict(G0=G0.tolist(),H=H.tolist(),D=D.tolist(),x=x,y=y,a=a,b=b,O1=sorted(O1),O0=sorted(O0)))
    return dict(acc),wit

def main():
    n=int(sys.argv[1]); nw=int(sys.argv[2]); s0=int(sys.argv[3])
    jobs=[]
    fams=["er","chain","star"]
    for i in range(n):
        p=int(6+(i%5)); deg=[2.0,2.5,3.0,3.5][i%4]; fam=fams[i%3]
        jobs.append((s0+i,p,deg,fam))
    acc=defaultdict(int); W=[]; t0=time.time()
    with Pool(nw) as pool:
        for a,w in pool.imap_unordered(one,jobs,chunksize=8):
            for k,v in a.items(): acc[k]+=v
            W.extend(w[:2])
    acc["elapsed"]=round(time.time()-t0,1)
    print(json.dumps(dict(sorted(acc.items())),indent=1))
    if W: json.dump(W[:20],open("HUNT_WITNESS.json","w"),indent=1); print("!!!! COUNTEREXAMPLE",len(W))
if __name__=="__main__": main()
