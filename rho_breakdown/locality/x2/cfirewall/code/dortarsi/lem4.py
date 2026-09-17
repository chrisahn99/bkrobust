"""L4: does every proper non-causal definite-status path in an amenable MPDAG G0
have (i) a collider in forb(G0), or (ii) a definite non-collider in O*(X,Y,G0)?"""
import sys, json
from collections import defaultdict
import numpy as np
BASE="/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious"
sys.path.insert(0,BASE+"/x2/code")
sys.path.insert(0,"/private/tmp/claude-501/-Users-josecosta-mugango/bceac0a4-93a7-4757-b3b3-a3d4e6c42d55/scratchpad/dt")
from graphs import dag_to_cpdag, v_structures, is_directed
import adjust, x2lib as X
from run_lemma import all_dags, reachable_mpdags
from gac import all_paths, node_status, is_causal

def scan(C):
    acc=defaultdict(int); wit=[]
    for G0 in reachable_mpdags(C):
        p=G0.shape[0]
        for x in range(p):
            for y in range(p):
                if x==y: continue
                O0,n0,am0=X.ostar_and_paths(G0,x,y)
                if not am0: continue
                F0=adjust.forb(G0,x,y); CN0=adjust.causal_nodes(G0,x,y)
                for path in all_paths(G0,x,y):
                    if len(path)<2: continue
                    st=[node_status(G0,path[i-1],path[i],path[i+1]) for i in range(1,len(path)-1)]
                    if any(s=="undef" for s in st): continue
                    if is_causal(G0,path): continue
                    acc["path"]+=1
                    colF = any(s=="col" and path[i] in F0 for i,s in enumerate(st,1))
                    dncO = any(s=="dnc" and path[i] in O0 for i,s in enumerate(st,1))
                    acc[("L4",int(colF),int(dncO))]+=1
                    if not (colF or dncO):
                        acc["L4_FAIL"]+=1
                        if len(wit)<6: wit.append(dict(G0=G0.tolist(),x=int(x),y=int(y),
                            path=list(map(int,path)),st=st,O0=sorted(map(int,O0)),
                            F0=sorted(map(int,F0)),CN0=sorted(map(int,CN0))))
                    # canonical blocker
                    idx=[i for i in range(len(path)) if path[i] not in CN0]
                    j=max(idx); wj=path[j]
                    caseA = j+1<len(path) and is_directed(G0,wj,path[j+1])
                    acc[("canon",int(caseA),int(wj in O0),int(colF))]+=1
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
            if len(W)<6: W.extend(w[:6-len(W)])
    for k in sorted(tot,key=str): print(f"{k:35s} {tot[k]}")
    json.dump(dict(counts={k:int(v) for k,v in tot.items()},wit=W),open(f"lem4_p{p}.json","w"))
