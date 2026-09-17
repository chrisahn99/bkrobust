"""Coverage of the two PROVED certificates, and a hunt in the residual.
CERT-A (lift):  D+(a->b) in [H]                      -> HPM(H) + deletion lemma
CERT-B (dom):   cn_D != {} and O*_D subset O*_H      -> superset lemma (+ O*_H cap forb_D = {})
RESIDUAL: neither.  Any violation must live there."""
import sys; sys.path.insert(0,'/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/cfirewall/proof')
from common import *
from collections import defaultdict
import numpy as np, json, time
from multiprocessing import Pool

def ostar_raw(G,x,y):
    paths = X.pcp_capped(G,x,y)
    if not paths: return None,0,False
    amen = all(is_directed(G,pth[0],pth[1]) for pth in paths)
    cnv=set()
    for pth in paths: cnv.update(pth[1:])
    fb = (adjust.poss_de(G,cnv) | {x}) if cnv else {x}
    return frozenset(adjust.parents_of_set(G,cnv)-fb), len(paths), amen

def scan_mpdag(C,G0,ref,acc,wit,maxu=15):
    p=G0.shape[0]
    NA = nonadjacent_pairs(G0)
    if not NA: return
    if len(undirected_edges(G0))>maxu: acc["skip_bigclass"]+=1; return
    QA=[]
    for x in range(p):
        for y in range(p):
            if x==y: continue
            O0,np0,am0 = ostar_raw(G0,x,y)
            if O0 is not None and am0: QA.append((x,y,O0,adjust.forb(G0,x,y)))
    if not QA: return
    ext = consistent_dag_extensions(G0, ref_vstructs=ref)
    if not ext: return
    Dcache={}
    for di,D in enumerate(ext):
        for (x,y,O0,fb0) in QA:
            cnD=adjust.causal_nodes(D,x,y); fbD=adjust.forb(D,x,y)
            Dcache[(di,x,y)]=(cnD,fbD,adjust.parents_of_set(D,cnD)-fbD)
    for (u,v) in NA:
        for (a,b) in ((u,v),(v,u)):
            H,info = bk_assert(G0,[(a,b)])
            if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
            vH=v_structures(H); dirH_S=[(s,t) for (s,t) in directed_edges(H) if not (s==a and t==b)]
            lift=[]
            for D in ext:
                Dp=D.copy(); Dp[a,b]=1
                lift.append((not has_directed_cycle(Dp)) and all(is_directed(Dp,s,t) for (s,t) in dirH_S) and v_structures(Dp)==vH)
            for (x,y,O0,fb0) in QA:
                OH,nph,amH = ostar_raw(H,x,y)
                if OH is None or not amH: continue
                acc["trial"]+=1
                f1 = not (OH & fb0)
                if not f1: acc["F1_viol"]+=1
                for di,D in enumerate(ext):
                    cnD,fbD,OD = Dcache[(di,x,y)]
                    acc["check"]+=1
                    A = lift[di]
                    B = bool(cnD) and (OD<=OH) and not (OH & fbD)
                    if A: acc["certA"]+=1
                    if B: acc["certB"]+=1
                    if A or B: acc["covered"]+=1; continue
                    acc["residual"]+=1
                    ok = adjust.is_valid_adjustment_set(D,x,y,set(OH))
                    if not ok:
                        acc["RESIDUAL_VIOL"]+=1
                        if len(wit)<5: wit.append((C.tolist(),G0.tolist(),D.tolist(),int(x),int(y),int(a),int(b),sorted(int(z) for z in OH)))

def job_census(kb_p):
    kb,p=kb_p
    C=np.frombuffer(kb,dtype=np.int8).reshape(p,p).copy()
    ref=v_structures(C); acc=defaultdict(int); wit=[]
    for G0 in reachable_mpdags(C): scan_mpdag(C,G0,ref,acc,wit)
    return dict(acc),wit

def job_sample(args):
    p,n,seed=args
    rng=np.random.default_rng(seed); acc=defaultdict(int); wit=[]
    for _ in range(n):
        deg=float(rng.choice([1.5,2.0,2.5,3.0,3.5,4.0]))
        D=random_dag(p,deg,rng); C=dag_to_cpdag(D); ref=v_structures(C)
        G=C.copy(); t=int(rng.integers(0,len(undirected_edges(C))+1))
        for _ in range(t):
            U=undirected_edges(G)
            if not U: break
            uu,vv=U[int(rng.integers(len(U)))]
            aa,bb=(uu,vv) if rng.random()<0.5 else (vv,uu)
            Hh=G.copy(); Hh[bb,aa]=0; Hh=meek_closure(Hh)
            if has_directed_cycle(Hh) or v_structures(Hh)!=ref: continue
            G=Hh
        acc["draw"]+=1
        scan_mpdag(C,G,ref,acc,wit)
    return dict(acc),wit

def main():
    mode=sys.argv[1]; p=int(sys.argv[2]); nw=int(sys.argv[3])
    if mode=="census":
        cps={}
        for D in all_dags(p):
            C=dag_to_cpdag(D); cps.setdefault(C.tobytes(),C)
        jobs=[(k,p) for k in cps]; fn=job_census
    else:
        n=int(sys.argv[4]); seed=int(sys.argv[5]) if len(sys.argv)>5 else 1
        nb=200; jobs=[(p,max(1,n//nb),seed+7919*i) for i in range(nb)]; fn=job_sample
    tot=defaultdict(int); W=[]; t0=time.time(); done=0
    with Pool(nw) as pool:
        for acc,w in pool.imap_unordered(fn,jobs,chunksize=1):
            for k,v in acc.items(): tot[k]+=v
            W.extend(w[:max(0,5-len(W))]); done+=1
            if done%max(1,len(jobs)//10)==0:
                print("  %d/%d %.0fs check=%d resid=%d VIOL=%d"%(done,len(jobs),time.time()-t0,tot["check"],tot["residual"],tot["RESIDUAL_VIOL"]),flush=True)
    print(mode,p,json.dumps(dict(sorted(tot.items())),indent=1),flush=True)
    for it in W: print("RESIDUAL_WIT",json.dumps(it),flush=True)
if __name__=="__main__": main()
