"""Parallel p=5 census of the two remaining sharpness questions:
   (Sb) reorienting an UNDIRECTED edge of G0 against D  -> can O*_H be invalid?
   (Sc) H coherent but NOT amenable, report pa(cn)\\forb anyway -> invalid?"""
from common import *
from collections import defaultdict
import json,sys

def ostar_raw(G,x,y):
    cn=adjust.causal_nodes(G,x,y)
    if not cn: return frozenset()
    return frozenset(adjust.parents_of_set(G,cn)-adjust.forb(G,x,y))

def scan(C):
    p=C.shape[0]; ref=v_structures(C); acc=defaultdict(int); wit=[]
    for G0 in reachable_mpdags(C):
        ext=consistent_dag_extensions(G0,ref_vstructs=ref)
        Q=[]
        for x in range(p):
            for y in range(p):
                if x!=y:
                    O0,n0,am0=X.ostar_and_paths(G0,x,y); Q.append((x,y,O0,am0))
        # Sc : non-adjacent statement, coherent H, H NOT amenable
        for (a0,b0) in nonadjacent_pairs(G0):
            for (a,b) in ((a0,b0),(b0,a0)):
                H,info=bk_assert(G0,[(a,b)])
                if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
                for (x,y,O0,am0) in Q:
                    if not am0: continue
                    O1,n1,am1=X.ostar_and_paths(H,x,y)
                    if am1: continue
                    Or=ostar_raw(H,x,y)
                    for D in ext:
                        acc["Sc_n"]+=1
                        if not adjust.is_valid_adjustment_set(D,x,y,set(Or)):
                            acc["Sc_VIOL"]+=1
                            if len(wit)<3: wit.append(("Sc",C.tolist(),G0.tolist(),H.tolist(),D.tolist(),a,b,x,y,sorted(map(int,Or))))
        # Sb : reorient an UNDIRECTED edge of G0 the wrong way
        for (u,v) in undirected_edges(G0):
            for (a,b) in ((u,v),(v,u)):
                H,info=bk_assert(G0,[(a,b)])
                if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
                for (x,y,O0,am0) in Q:
                    if not am0: continue
                    O1,n1,am1=X.ostar_and_paths(H,x,y)
                    if not am1: continue
                    for D in ext:
                        if not is_directed(D,b,a): continue   # statement must be FALSE of D
                        acc["Sb_n"]+=1
                        if not adjust.is_valid_adjustment_set(D,x,y,set(O1)):
                            acc["Sb_VIOL"]+=1
                            if len(wit)<3: wit.append(("Sb",C.tolist(),G0.tolist(),H.tolist(),D.tolist(),a,b,x,y,sorted(map(int,O1))))
    return dict(acc),wit

def _job(t):
    kb,p=t
    return scan(np.frombuffer(kb,dtype=np.int8).reshape(p,p).copy())

def main():
    p=int(sys.argv[1]); nw=int(sys.argv[2]) if len(sys.argv)>2 else 10
    cp={}
    for D in all_dags(p):
        C=dag_to_cpdag(D); cp.setdefault(C.tobytes(),C)
    tot=defaultdict(int); W=[]
    from multiprocessing import Pool
    with Pool(nw) as pool:
        for acc,w in pool.imap_unordered(_job,[(k,p) for k in cp],chunksize=4):
            for k,v in acc.items(): tot[k]+=v
            if len(W)<3: W.extend(w[:3-len(W)])
    print(json.dumps(dict(sorted(tot.items())),indent=1))
    json.dump(W,open(f"sharpwit_p{p}.json","w"),indent=1)
if __name__=='__main__': main()
