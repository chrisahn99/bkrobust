"""Sampled coverage of the two proved routes at p=6,7 and for k>=2 spurious edges."""
import sys, json
from collections import defaultdict
from itertools import combinations
import numpy as np
sys.path.insert(0,"/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/cfirewall/proof-pathmono")
from lib import *

def random_mpdag(C,ref,rng):
    G=C.copy(); t=int(rng.integers(0,len(undirected_edges(C))+1))
    for _ in range(t):
        U=undirected_edges(G)
        if not U: break
        u,v=U[int(rng.integers(len(U)))]
        a,b=(u,v) if rng.random()<0.5 else (v,u)
        Hh=G.copy(); Hh[b,a]=0; Hh=meek_closure(Hh)
        if has_directed_cycle(Hh) or v_structures(Hh)!=ref: continue
        G=Hh
    return G

def block(args):
    p,n,seed,maxu,k=args
    rng=np.random.default_rng(seed); acc=defaultdict(int)
    for _ in range(n):
        deg=float(rng.choice([1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0]))
        Dr=random_dag(p,deg,rng); C=dag_to_cpdag(Dr); ref=v_structures(C)
        G0=random_mpdag(C,ref,rng)
        if len(undirected_edges(G0))>maxu: continue
        NA=nonadjacent_pairs(G0)
        if len(NA)<k: continue
        ext=consistent_dag_extensions(G0,ref_vstructs=ref)
        if not ext: continue
        Q=[]
        for x in range(p):
            for y in range(p):
                if x==y: continue
                O0,_,am0=X.ostar_and_paths(G0,x,y)
                if am0: Q.append((x,y,O0))
        if not Q: continue
        stmts=[]
        for pair in combinations(range(len(NA)),k):
            for signs in range(2**k):
                K=[]
                for t,ix in enumerate(pair):
                    u,v=NA[ix]
                    K.append((u,v) if (signs>>t)&1 else (v,u))
                stmts.append(K)
        if len(stmts)>60:
            sel=rng.choice(len(stmts),60,replace=False); stmts=[stmts[i] for i in sel]
        ODc={}
        for K in stmts:
            H,info=bk_assert(G0,K)
            if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
            acc["coherent_stmt"]+=1
            for (x,y,O0) in Q:
                O1,_,am1=X.ostar_and_paths(H,x,y)
                if not am1: acc["abort"]+=1; continue
                for iD,D in enumerate(ext):
                    acc["trial"]+=1
                    if (iD,x,y) not in ODc:
                        ODc[(iD,x,y)]=(adjust.optimal_adjustment_set(D,x,y),adjust.forb(D,x,y))
                    OD,fbD=ODc[(iD,x,y)]
                    Dp=D.copy()
                    for (u,v) in K: Dp[u,v]=1; Dp[v,u]=0
                    A=(not has_directed_cycle(Dp)) and in_class(Dp,H)
                    S=(OD is not None) and set(OD)<=set(O1)
                    acc["A"]+=int(A); acc["S"]+=int(S); acc["AorS"]+=int(A or S)
                    if set(O1)&fbD: acc["OH_meets_fbD"]+=1
                    if not (A or S):
                        acc["resid"]+=1
                        acc["resid_valid"]+=int(adjust.is_valid_adjustment_set(D,x,y,set(O1)))
    return dict(acc)

if __name__=="__main__":
    p=int(sys.argv[1]); ndraw=int(sys.argv[2]); k=int(sys.argv[3]); nw=int(sys.argv[4]); maxu=int(sys.argv[5]) if len(sys.argv)>5 else 10
    from multiprocessing import Pool
    nb=200; per=max(1,ndraw//nb)
    tot=defaultdict(int)
    with Pool(nw) as pool:
        for acc in pool.imap_unordered(block,[(p,per,1000+7919*i,maxu,k) for i in range(nb)]):
            for kk,v in acc.items(): tot[kk]+=v
    print("p=%d k=%d draws=%d"%(p,k,nb*per),json.dumps(dict(sorted(tot.items())),indent=1),flush=True)
