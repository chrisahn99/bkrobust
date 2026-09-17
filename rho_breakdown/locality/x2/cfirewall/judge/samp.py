"""Sampled census at p>=6: full sweep of each sampled G0 (all queries x all non-adj pairs x all D)
   + certificate coverage + the verification battery."""
import sys, json, time
from collections import defaultdict
from multiprocessing import Pool
import numpy as np
from jlib import *

def one(arg):
    seed,p,deg=arg
    rng=np.random.default_rng(seed)
    acc=defaultdict(int); wit=[]
    D0=random_dag(p,deg,rng); C=dag_to_cpdag(D0); ref=v_structures(C)
    # random MPDAG: orient a random subset of undirected edges consistently
    U=undirected_edges(C)
    G0=C.copy()
    k=int(rng.integers(0,min(3,len(U))+1)) if U else 0
    for _ in range(k):
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
    ext=consistent_dag_extensions(G0,ref_vstructs=ref,limit=64)
    if not ext: return dict(acc),wit
    acc["G0"]=1; acc["extsum"]=len(ext)
    DQ={}
    for di,D in enumerate(ext):
        for (x,y,_) in Q:
            if (di,x,y) in DQ: continue
            cnD=causal_nodes(D,x,y); fbD=(poss_de(D,cnD)|{x}) if cnD else {x}
            R=reach_outside_forb(D,x,y,cnD,fbD) if cnD else set()
            need=((parents_of_set(D,cnD)-fbD)&R) if cnD else None
            paX={u for u in range(p) if is_directed(D,u,x)}
            deX=poss_de(D,{x})
            DQ[(di,x,y)]=(cnD,fbD,need,paX,deX)
    for (a0,b0) in NA:
        for (a,b) in ((a0,b0),(b0,a0)):
            H,info=bk_assert(G0,[(a,b)])
            c2=not has_directed_cycle(H); c3=pdag_extendable(H)
            acc["stmt"]+=1
            if info["conflict"]: acc["c1_fail"]+=1
            if c2 and not c3: acc["acyc_notext"]+=1
            if c3 and not c2: acc["ext_cyc"]+=1
            if not((not info["conflict"]) and c2 and c3): continue
            liftc={}
            for (x,y,O0) in Q:
                O1,n1,am1,cn1,fb1=ostar(H,x,y)
                acc["trial"]+=1
                if not am1:
                    acc["abort"]+=1; acc["abort_nopath"]+=(n1==0); continue
                acc["report"]+=1; S1=set(O1); acc["moved"]+=(O1!=O0)
                pdX0=poss_de(G0,{x})
                if (parents_of_set(H,cn1)&pdX0)-(cn1|{x}): acc["U0_FAIL"]+=1
                if not (pdX0<=poss_de(H,{x})):
                    acc["Mfail"]+=1; acc["Mfail_moved"]+=(O1!=O0)
                nb=0
                for di,D in enumerate(ext):
                    if di not in liftc: liftc[di]=lifts(D,a,b,H)[0]
                    A=liftc[di]
                    cnD,fbD,need,paX,deX=DQ[(di,x,y)]
                    L=(need is not None) and (not(S1&fbD)) and (need<=S1)
                    Cc=(paX<=S1) and not (S1&deX)
                    acc["n"]+=1; acc["A"]+=A; acc["L"]+=L; acc["C"]+=Cc
                    acc["ALC"]+=(A or L or Cc)
                    if S1&deX: acc["deX_hit"]+=1
                    v=is_valid_adjustment_set(D,x,y,S1)
                    if not v:
                        nb+=1; acc["VIOL"]+=1
                        if len(wit)<5: wit.append(dict(G0=G0.tolist(),H=H.tolist(),D=D.tolist(),x=x,y=y,a=a,b=b,O1=sorted(O1)))
                    if not (A or L or Cc):
                        acc["UNCOV"]+=1; acc["UNCOV_invalid"]+=(not v)
                if nb and nb!=len(ext): acc["split"]+=1
    return dict(acc),wit

def main():
    p=int(sys.argv[1]); n=int(sys.argv[2]); nw=int(sys.argv[3]); deg=float(sys.argv[4]); s0=int(sys.argv[5])
    acc=defaultdict(int); W=[]; t0=time.time()
    with Pool(nw) as pool:
        for a,w in pool.imap_unordered(one,[(s0+i,p,deg) for i in range(n)],chunksize=8):
            for k,v in a.items(): acc[k]+=v
            W.extend(w[:2])
    acc["elapsed"]=round(time.time()-t0,1)
    print("p=%d n=%d deg=%s"%(p,n,deg),json.dumps(dict(sorted(acc.items())),indent=1))
    if W: json.dump(W[:20],open("sampwit_p%d.json"%p,"w"),indent=1); print("!!! VIOLATION WITNESSES",len(W))
if __name__=="__main__": main()
