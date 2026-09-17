"""BIG CONJECTURE at p>=5 by sampling monotone MPDAG refinements H of an amenable MPDAG G0.
   Also reports coverage of the two proved certificates (lift / dominance)."""
import sys; sys.path.insert(0,'/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/cfirewall/proof')
from common import *
from collections import defaultdict
import numpy as np, json, time
from multiprocessing import Pool

def ostar_raw(G,x,y):
    paths=X.pcp_capped(G,x,y)
    if not paths: return None,False
    amen=all(is_directed(G,q[0],q[1]) for q in paths)
    cnv=set()
    for q in paths: cnv.update(q[1:])
    fb=(adjust.poss_de(G,cnv)|{x}) if cnv else {x}
    return frozenset(adjust.parents_of_set(G,cnv)-fb),amen

def coherent(H):
    return (not has_directed_cycle(H)) and pdag_extendable(H) and np.array_equal(meek_closure(H),H)

def job(args):
    p,n,seed=args
    rng=np.random.default_rng(seed); acc=defaultdict(int); wit=[]
    for _ in range(n):
        deg=float(rng.choice([1.5,2.0,2.5,3.0,3.5,4.0]))
        D0=random_dag(p,deg,rng); C=dag_to_cpdag(D0); ref=v_structures(C)
        G0=C.copy(); t=int(rng.integers(0,len(undirected_edges(C))+1))
        for _ in range(t):
            U=undirected_edges(G0)
            if not U: break
            uu,vv=U[int(rng.integers(len(U)))]
            aa,bb=(uu,vv) if rng.random()<0.5 else (vv,uu)
            Hh=G0.copy(); Hh[bb,aa]=0; Hh=meek_closure(Hh)
            if has_directed_cycle(Hh) or v_structures(Hh)!=ref: continue
            G0=Hh
        if len(undirected_edges(G0))>13: continue
        QA=[]
        for x in range(p):
            for y in range(p):
                if x==y: continue
                O0,a0=ostar_raw(G0,x,y)
                if O0 is not None and a0: QA.append((x,y,O0))
        if not QA: continue
        ext=consistent_dag_extensions(G0,ref_vstructs=ref)
        if not ext: continue
        acc["draw"]+=1
        # build monotone refinements H
        for _h in range(8):
            H=G0.copy(); nadd=0
            for _step in range(int(rng.integers(1,4))):
                cand=[]
                for i in range(p):
                    for j in range(p):
                        if i==j: continue
                        if H[i,j]==0 and H[j,i]==0: cand.append(("add",i,j))
                        elif is_undirected(H,i,j): cand.append(("ori",i,j))
                if not cand: break
                k,i,j=cand[int(rng.integers(len(cand)))]
                H2=H.copy(); H2[i,j]=1; H2[j,i]=0; H2=meek_closure(H2)
                for (u,v) in directed_edges(G0):
                    if not is_directed(H2,u,v): H2=None; break
                if H2 is None: continue
                if not coherent(H2): continue
                H=H2; nadd+= (k=="add")
            if np.array_equal(H,G0): continue
            skG=skeleton(G0)
            if not all(is_directed(H,u,v) for (u,v) in directed_edges(G0)): continue
            for (x,y,O0) in QA:
                OH,aH=ostar_raw(H,x,y)
                if OH is None or not aH: continue
                acc["trial"]+=1
                if OH!=O0: acc["moved"]+=1
                extra=[(u,v) for (u,v) in directed_edges(H)+undirected_edges(H) if skG[u,v]==0]
                for D in ext:
                    acc["check"]+=1
                    cnD=adjust.causal_nodes(D,x,y); fbD=adjust.forb(D,x,y)
                    OD=adjust.parents_of_set(D,cnD)-fbD
                    certB = bool(cnD) and (OD<=OH) and not (OH & fbD)
                    Dp=D.copy(); okA=True
                    for (u,v) in extra:
                        if is_directed(H,u,v): Dp[u,v]=1
                        else: okA=False
                    certA = okA and (not has_directed_cycle(Dp)) and all(is_directed(Dp,u,v) for (u,v) in directed_edges(H)) and v_structures(Dp)==v_structures(H)
                    if certA: acc["certA"]+=1
                    if certB: acc["certB"]+=1
                    if certA or certB: acc["covered"]+=1
                    else: acc["residual"]+=1
                    if not adjust.is_valid_adjustment_set(D,x,y,set(OH)):
                        acc["VIOL"]+=1
                        if not (certA or certB): acc["VIOL_residual"]+=1
                        if len(wit)<4: wit.append((G0.tolist(),H.tolist(),D.tolist(),int(x),int(y),sorted(int(z) for z in OH)))
    return dict(acc),wit

def main():
    p=int(sys.argv[1]); nw=int(sys.argv[2]); n=int(sys.argv[3]); seed=int(sys.argv[4]) if len(sys.argv)>4 else 3
    nb=100; jobs=[(p,max(1,n//nb),seed+7919*i) for i in range(nb)]
    tot=defaultdict(int); W=[]; t0=time.time()
    with Pool(nw) as pool:
        for acc,w in pool.imap_unordered(job,jobs,chunksize=1):
            for k,v in acc.items(): tot[k]+=v
            W.extend(w[:max(0,4-len(W))])
    print("p=%d %.0fs"%(p,time.time()-t0)); print(json.dumps(dict(sorted(tot.items())),indent=1))
    for it in W: print("WIT",json.dumps(it))
if __name__=="__main__": main()
