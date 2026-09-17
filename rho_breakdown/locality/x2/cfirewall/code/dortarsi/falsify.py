"""Falsification check of Theorem A:
   (E1) directed(H) = directed(G0) + {a->b}   (Meek adds nothing)
   (E2) a not in poss_de(b,G0)
  claim: for EVERY D in [G0], D+(a->b) is in [H]; hence (if H amenable) O*(X,Y,H)
  is valid in D -- WITHOUT any amenability assumption on G0.
Also re-checks Lemma M (edge addition) directly."""
import sys
from collections import defaultdict
import numpy as np
BASE="/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious"
sys.path.insert(0,BASE+"/x2/code"); sys.path.insert(0,BASE+"/x2/cfirewall/code")
from graphs import (consistent_dag_extensions, dag_to_cpdag, dag_agrees_with,
                    has_directed_cycle, v_structures, directed_edges)
import adjust, x2lib as X
from run_lemma import all_dags, reachable_mpdags
import cfirewall as CF
x1=CF.x1_ops
def scan(C):
    ref=v_structures(C); acc=defaultdict(int); wit=[]
    for G0 in reachable_mpdags(C):
        p=G0.shape[0]; NA=x1.nonadjacent_pairs(G0)
        if not NA: continue
        ext=consistent_dag_extensions(G0,ref_vstructs=ref)
        if not ext: continue
        for (a0,b0) in NA:
            for (a,b) in ((a0,b0),(b0,a0)):
                H,info=x1.bk_assert(G0,[(a,b)])
                E1 = set(directed_edges(H))==set(directed_edges(G0))|{(a,b)}
                E2 = a not in adjust.poss_de(G0,{b})
                if not (E1 and E2): continue
                acc["stmt_E1E2"]+=1
                acc["coh_auto"] += (not info["conflict"]) and (not has_directed_cycle(H)) and x1.pdag_extendable(H)
                for x in range(p):
                    for y in range(p):
                        if x==y: continue
                        O0,n0,am0=X.ostar_and_paths(G0,x,y)
                        Z,n1,am1=X.ostar_and_paths(H,x,y)
                        if not am1: continue
                        lbl = "G0amen" if am0 else ("G0nopath" if n0==0 else "G0unamen")
                        nbad=sum(1 for D in ext if not adjust.is_valid_adjustment_set(D,x,y,set(Z)))
                        acc[(lbl,"rep")]+=1
                        acc[(lbl,"BAD")] += (nbad>0)
                        if nbad and len(wit)<4:
                            wit.append(dict(G0=G0.tolist(),H=H.tolist(),a=int(a),b=int(b),x=int(x),y=int(y),
                                            Z=sorted(map(int,Z)),lbl=lbl))
                # Lemma M direct: valid in D+ => valid in D, over all queries/subsets
                import itertools
                for D in ext:
                    Dp=D.copy()
                    if has_directed_cycle(np.where(np.eye(p,dtype=bool),0,D)) : pass
                    Dp[a,b]=1; Dp[b,a]=0
                    if has_directed_cycle(Dp): continue
                    for x in range(p):
                        for y in range(p):
                            if x==y: continue
                            oth=[v for v in range(p) if v not in (x,y)]
                            for r in range(len(oth)+1):
                                for S in itertools.combinations(oth,r):
                                    vp=adjust.is_valid_adjustment_set(Dp,x,y,set(S))
                                    if vp:
                                        acc["M_chk"]+=1
                                        if not adjust.is_valid_adjustment_set(D,x,y,set(S)):
                                            acc["M_FAIL"]+=1
    return {str(k):v for k,v in acc.items()},wit
def _job(arg):
    kb,p=arg
    return scan(np.frombuffer(kb,dtype=np.int8).reshape(p,p).copy())
if __name__=="__main__":
    p=int(sys.argv[1]); cp={}
    for D in all_dags(p): C=dag_to_cpdag(D); cp.setdefault(C.tobytes(),C)
    from multiprocessing import Pool
    tot=defaultdict(int); W=[]
    with Pool(4) as pool:
        for acc,w in pool.imap_unordered(_job,[(k,p) for k in cp],chunksize=4):
            for k,v in acc.items(): tot[k]+=v
            if len(W)<4: W.extend(w[:4-len(W)])
    for k in sorted(tot,key=str): print(f"{k:28s} {tot[k]}")
    if W: print("WITNESS",W[0])
