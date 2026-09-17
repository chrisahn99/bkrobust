"""(1) is H the common orientation graph of [H]  (i.e. a bona fide MPDAG)?
   (2) does (E1) no-propagation + (E2) a not in poss_de(b,G0) force D+(a->b) in [H] for all D?
       and what does it cover?
   (3) if G0 is NOT amenable but H is, is the report ALWAYS invalid?  (sharpness of the
       G0-amenability hypothesis)"""
import sys, json
from collections import defaultdict
import numpy as np
BASE="/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious"
sys.path.insert(0,BASE+"/x2/code"); sys.path.insert(0,BASE+"/x2/cfirewall/code")
from graphs import (consistent_dag_extensions, dag_to_cpdag, dag_agrees_with,
                    common_orientation_graph, has_directed_cycle, skeleton,
                    v_structures, directed_edges, is_directed)
import adjust, x2lib as X
from run_lemma import all_dags, reachable_mpdags
import cfirewall as CF
x1=CF.x1_ops

def scan(C):
    ref=v_structures(C); acc=defaultdict(int); wit=[]
    for G0 in reachable_mpdags(C):
        p=G0.shape[0]
        NA=x1.nonadjacent_pairs(G0)
        if not NA: continue
        ext=consistent_dag_extensions(G0,ref_vstructs=ref)
        if not ext: continue
        QA=[]; QN=[]
        for x in range(p):
            for y in range(p):
                if x==y: continue
                O0,n0,am0=X.ostar_and_paths(G0,x,y)
                (QA if am0 else QN).append((x,y,O0))
        for (a0,b0) in NA:
            for (a,b) in ((a0,b0),(b0,a0)):
                H,info=x1.bk_assert(G0,[(a,b)])
                cyc=has_directed_cycle(H); extok=x1.pdag_extendable(H)
                coh=(not info["conflict"]) and (not cyc) and extok
                arm="T" if coh else "K"
                acc[arm+"_stmt"]+=1
                # (1) H a bona fide MPDAG?
                extH=consistent_dag_extensions(H,ref_vstructs=v_structures(H))
                acc[arm+"_extH_empty"]+= (len(extH)==0)
                if extH:
                    cog=common_orientation_graph(extH,skeleton(H))
                    acc[arm+"_H_is_MPDAG"]+= bool(np.array_equal(cog,H))
                    if not np.array_equal(cog,H) and len(wit)<4:
                        wit.append(dict(kind="notMPDAG",G0=G0.tolist(),H=H.tolist(),cog=cog.tolist(),a=int(a),b=int(b)))
                # (2) E1/E2
                E1 = set(directed_edges(H)) == set(directed_edges(G0))|{(a,b)}
                E2 = a not in adjust.poss_de(G0,{b})
                acc[arm+"_E1"]+=E1; acc[arm+"_E2"]+=E2; acc[arm+"_E1E2"]+=(E1 and E2)
                if E1 and E2:
                    allin=True
                    for D in ext:
                        Dp=D.copy(); Dp[a,b]=1; Dp[b,a]=0
                        if has_directed_cycle(Dp) or not dag_agrees_with(Dp,H) or v_structures(Dp)!=v_structures(H):
                            allin=False
                    acc[arm+"_E1E2_allin"]+= allin
                    if not allin and len(wit)<4:
                        wit.append(dict(kind="E12fail",G0=G0.tolist(),H=H.tolist(),a=int(a),b=int(b)))
                # coverage of E1E2 among MOVED reports
                for (x,y,O0) in QA:
                    Z,n1,am1=X.ostar_and_paths(H,x,y)
                    if not am1: continue
                    mv = (Z!=O0)
                    acc[(arm,"cov",int(E1 and E2),int(mv))]+=1
                # (3) G0 not amenable, H amenable
                for (x,y,O0) in QN:
                    Z,n1,am1=X.ostar_and_paths(H,x,y)
                    if not am1: continue
                    acc[arm+"_gain_report"]+=1
                    nbad=sum(1 for D in ext if not adjust.is_valid_adjustment_set(D,x,y,set(Z)))
                    if nbad==0: acc[arm+"_gain_ALLVALID"]+=1
                    elif nbad==len(ext): acc[arm+"_gain_allbad"]+=1
                    else: acc[arm+"_gain_somebad"]+=1
    return {str(k):v for k,v in acc.items()},wit

def _job(arg):
    kb,p=arg
    return scan(np.frombuffer(kb,dtype=np.int8).reshape(p,p).copy())

if __name__=="__main__":
    p=int(sys.argv[1]); cp={}
    for D in all_dags(p): C=dag_to_cpdag(D); cp.setdefault(C.tobytes(),C)
    from multiprocessing import Pool
    tot=defaultdict(int); W=[]
    with Pool(8) as pool:
        for acc,w in pool.imap_unordered(_job,[(k,p) for k in cp],chunksize=4):
            for k,v in acc.items(): tot[k]+=v
            if len(W)<4: W.extend(w[:4-len(W)])
    for k in sorted(tot,key=str): print(f"{k:35s} {tot[k]}")
    json.dump(dict(counts={k:int(v) for k,v in tot.items()},wit=W),open(f"final_p{p}.json","w"))
