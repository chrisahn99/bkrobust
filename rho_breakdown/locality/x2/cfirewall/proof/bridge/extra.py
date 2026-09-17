"""(i) L1 over ALL H (incoherent included): (C1)&(C2) ==> D+ preserves H ==> H acyclic&extendable
   (ii) per-TRIAL coverage of the 'for every D' version
   (iii) also: is G0-amenability needed?  test THEOREM A conclusion on NON-amenable G0 queries."""
from common import *
from collections import defaultdict
from itertools import combinations
import json,sys,time

def uc(W):
    p=W.shape[0]; out=set()
    for v in range(p):
        pars=[u for u in range(p) if is_directed(W,u,v)]
        for u,w in combinations(sorted(pars),2):
            if not adjacent(W,u,w): out.add((u,v,w))
    return out

def scan(C):
    p=C.shape[0]; ref=v_structures(C); acc=defaultdict(int)
    for G0 in reachable_mpdags(C):
        NA=nonadjacent_pairs(G0)
        if not NA: continue
        ext=consistent_dag_extensions(G0,ref_vstructs=ref)
        Qam=[];Qnon=[]
        for x in range(p):
            for y in range(p):
                if x==y: continue
                O0,n0,am0=X.ostar_and_paths(G0,x,y)
                (Qam if am0 else Qnon).append((x,y,O0,am0))
        for (a0,b0) in NA:
            for (a,b) in ((a0,b0),(b0,a0)):
                H,info=bk_assert(G0,[(a,b)])
                coh=(not info["conflict"]) and (not has_directed_cycle(H)) and pdag_extendable(H)
                P=G0.copy(); P[a,b]=1;P[b,a]=0; ucP=uc(P); dH=directed_edges(H)
                good=[]
                for D in ext:
                    Dp=D.copy(); Dp[a,b]=1
                    ok=(not has_directed_cycle(Dp)) and (uc(Dp)<=ucP)
                    if ok:
                        acc["ok_pairs"]+=1
                        if not all(is_directed(Dp,i,j) for (i,j) in dH): acc["L1_FAIL_ALL"]+=1
                        if not coh: acc["INCOH_but_ok"]+=1
                    good.append(ok)
                allok=all(good)
                # (ii) per-trial coverage, on coherent H only
                if coh:
                    for (x,y,O0,am0) in Qam:
                        O1,n1,am1=X.ostar_and_paths(H,x,y)
                        if not am1: continue
                        acc["trial"]+=1
                        acc["trial_allok" if allok else "trial_not"]+=1
                    # (iii) G0 NOT amenable, H amenable: is O*_H still valid in every D?
                    for (x,y,O0,am0) in Qnon:
                        O1,n1,am1=X.ostar_and_paths(H,x,y)
                        if not am1: continue
                        for i,D in enumerate(ext):
                            acc["nonamenG0_trialD"]+=1
                            if good[i]: acc["nonamenG0_thmA"]+=1
                            if not adjust.is_valid_adjustment_set(D,x,y,set(O1)):
                                acc["nonamenG0_VIOLATION"]+=1
                                if good[i]: acc["nonamenG0_VIOL_inThmA"]+=1
    return dict(acc)

def _job(t):
    kb,p=t
    return scan(np.frombuffer(kb,dtype=np.int8).reshape(p,p).copy())

def main():
    p=int(sys.argv[1]); nw=int(sys.argv[2]) if len(sys.argv)>2 else 10
    cp={}
    for D in all_dags(p):
        C=dag_to_cpdag(D); cp.setdefault(C.tobytes(),C)
    tot=defaultdict(int)
    from multiprocessing import Pool
    with Pool(nw) as pool:
        for acc in pool.imap_unordered(_job,[(k,p) for k in cp],chunksize=4):
            for k,v in acc.items(): tot[k]+=v
    print(json.dumps(dict(sorted(tot.items())),indent=1))

if __name__=='__main__': main()
