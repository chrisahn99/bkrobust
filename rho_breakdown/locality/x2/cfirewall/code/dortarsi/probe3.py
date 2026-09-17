"""Path-level: for each proper non-causal definite-status path in G0, does the
canonical Case-A blocker w_j (j = max index outside cn_H) do the job?"""
import sys, json
from collections import defaultdict
import numpy as np
BASE="/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious"
sys.path.insert(0,BASE+"/x2/code"); sys.path.insert(0,BASE+"/x2/cfirewall/code")
sys.path.insert(0,"/private/tmp/claude-501/-Users-josecosta-mugango/bceac0a4-93a7-4757-b3b3-a3d4e6c42d55/scratchpad/dt")
from graphs import (dag_to_cpdag, has_directed_cycle, v_structures, is_directed,
                    is_undirected, adjacent)
import adjust, x2lib as X
from run_lemma import all_dags, reachable_mpdags
from gac import all_paths, node_status, is_causal
import cfirewall as CF
x1=CF.x1_ops

def scan(C):
    ref=v_structures(C); acc=defaultdict(int); wit=[]
    for G0 in reachable_mpdags(C):
        p=G0.shape[0]
        NA=x1.nonadjacent_pairs(G0)
        if not NA: continue
        Q=[]
        for x in range(p):
            for y in range(p):
                if x==y: continue
                O0,n0,am0=X.ostar_and_paths(G0,x,y)
                if not am0: continue
                F0=adjust.forb(G0,x,y); CN0=adjust.causal_nodes(G0,x,y)
                P=[]
                for path in all_paths(G0,x,y):
                    if len(path)<2: continue
                    st=[node_status(G0,path[i-1],path[i],path[i+1]) for i in range(1,len(path)-1)]
                    if any(s=="undef" for s in st): continue
                    if is_causal(G0,path): continue
                    P.append((path,st))
                Q.append((x,y,O0,F0,CN0,P))
        if not Q: continue
        for (a0,b0) in NA:
            for (a,b) in ((a0,b0),(b0,a0)):
                H,info=x1.bk_assert(G0,[(a,b)])
                coh=(not info["conflict"]) and (not has_directed_cycle(H)) and x1.pdag_extendable(H)
                arm="T" if coh else "K"
                for (x,y,O0,F0,CN0,P) in Q:
                    Z,n1,am1=X.ostar_and_paths(H,x,y)
                    if not am1: continue
                    CNH=adjust.causal_nodes(H,x,y); FH=adjust.forb(H,x,y)
                    acc[arm+"_rep"]+=1
                    acc[arm+"_P1_fail"]+= bool(set(Z)&F0)
                    for (path,st) in P:
                        acc[arm+"_path"]+=1
                        k=len(path)-1
                        idx=[i for i in range(k+1) if path[i] not in CNH]
                        if not idx:      # every node incl X in cn_H : impossible
                            acc[arm+"_noj"]+=1; continue
                        j=max(idx)
                        wj, wj1 = path[j], path[j+1]
                        caseA = is_directed(H,wj,wj1)
                        inZ = wj in Z
                        # is the path blocked at all (truth)?
                        blocked=False; how=None
                        for i,s in enumerate(st,start=1):
                            w=path[i]
                            if s=="dnc" and w in Z: blocked=True; how="dnc"; break
                            if s=="col" and not (adjust.poss_de(G0,{w}) & set(Z)): blocked=True; how="col"; break
                        acc[(arm,"blocked",int(blocked))]+=1
                        acc[(arm,"case",int(caseA),int(inZ),int(blocked))]+=1
                        if blocked and not (caseA and inZ) and len(wit)<400:
                            wit.append(dict(G0=G0.tolist(),H=H.tolist(),a=int(a),b=int(b),x=int(x),y=int(y),
                                Z=sorted(map(int,Z)),O0=sorted(map(int,O0)),path=list(map(int,path)),st=st,
                                j=j,caseA=bool(caseA),how=how,arm=arm,
                                cnH=sorted(map(int,CNH)),cn0=sorted(map(int,CN0)),
                                fH=sorted(map(int,FH)),f0=sorted(map(int,F0))))
    return {str(k):v for k,v in acc.items()}, wit

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
            if len(W)<400: W.extend(w[:400-len(W)])
    for k in sorted(tot,key=str): print(f"{k:45s} {tot[k]}")
    json.dump(dict(counts={k:int(v) for k,v in tot.items()},wit=W),open(f"probe3_p{p}.json","w"))
