"""L2: for every path pi in D^pbd from X to Y, pivot i=min{j>=1: v_j in cn_D};
   if v_{i-1}->v_i in D  then  v_{i-1} in O*_H ?   (Case A of the proof)"""
import sys; sys.path.insert(0,'/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/cfirewall/proof')
from common import *
from collections import defaultdict
import numpy as np, json
from multiprocessing import Pool

def all_paths(A,x,y):
    """all simple paths x..y in the *undirected* sense over the edge set of amat A (DAG)."""
    p=A.shape[0]; S=((A+A.T)>0)
    out=[]
    def dfs(path,vis):
        v=path[-1]
        if v==y: out.append(list(path)); return
        for w in range(p):
            if S[v,w] and w not in vis:
                vis.add(w); path.append(w); dfs(path,vis); path.pop(); vis.remove(w)
    dfs([x],{x}); return out

def blocked(D,path,Z):
    """is path blocked given Z in DAG D?"""
    n=len(path)
    for k in range(1,n-1):
        u,v,w=path[k-1],path[k],path[k+1]
        coll = D[u,v]==1 and D[w,v]==1
        if coll:
            # descendant in Z?
            de=set([v]); st=[v]
            while st:
                a=st.pop()
                for c in np.flatnonzero(D[a]):
                    c=int(c)
                    if c not in de: de.add(c); st.append(c)
            if not (de & set(Z)): return True
        else:
            if v in Z: return True
    return False

def scan_cpdag(kb_p):
    kb,p = kb_p
    C = np.frombuffer(kb,dtype=np.int8).reshape(p,p).copy()
    ref = v_structures(C); acc=defaultdict(int); wit=[]
    for G0 in reachable_mpdags(C):
        NA = nonadjacent_pairs(G0)
        if not NA: continue
        Q=[]
        for x in range(p):
            for y in range(p):
                if x==y: continue
                O0,np0,am0 = ostar(G0,x,y)
                if am0: Q.append((x,y,O0))
        if not Q: continue
        ext = consistent_dag_extensions(G0, ref_vstructs=ref)
        if not ext: continue
        for (a0,b0) in NA:
            for (a,b) in ((a0,b0),(b0,a0)):
                H,info = bk_assert(G0,[(a,b)])
                if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
                for (x,y,O0) in Q:
                    OH,nph,amH = ostar(H,x,y)
                    if not amH: continue
                    for D in ext:
                        cnD = cn(D,x,y)
                        if not cnD: continue
                        Dp=D.copy()
                        for w in cnD:
                            if is_directed(Dp,x,w): Dp[x,w]=0
                        for path in all_paths(Dp,x,y):
                            idx=[j for j in range(1,len(path)) if path[j] in cnD]
                            if not idx: continue
                            i=idx[0]; u=path[i-1]; v=path[i]
                            if not (Dp[u,v]==1 and Dp[v,u]==0):   # need u->v (Case A)
                                acc["caseB"]+=1; continue
                            acc["caseA"]+=1
                            if u not in OH:
                                acc["caseA_pivot_not_in_OH"]+=1
                                if len(wit)<6: wit.append(("PIVOTMISS",C.tolist(),G0.tolist(),D.tolist(),x,y,a,b,path,sorted(OH),sorted(cnD),sorted(cn(H,x,y)),sorted(forb(H,x,y))))
                                if not blocked(Dp,path,OH):
                                    acc["caseA_pivot_not_in_OH_and_OPEN"]+=1
                                    if len(wit)<3: wit.append((C.tolist(),G0.tolist(),D.tolist(),x,y,a,b,path,sorted(OH)))
    return dict(acc), wit

def main():
    p=int(sys.argv[1]); nw=int(sys.argv[2]) if len(sys.argv)>2 else 10
    cps={}
    for D in all_dags(p):
        C=dag_to_cpdag(D); cps.setdefault(C.tobytes(),C)
    jobs=[(k,p) for k in cps]
    tot=defaultdict(int); W=[]
    with Pool(nw) as pool:
        for acc,w in pool.imap_unordered(scan_cpdag,jobs,chunksize=4):
            for k,v in acc.items(): tot[k]+=v
            W.extend(w[:max(0,3-len(W))])
    print("p",p,json.dumps(dict(sorted(tot.items())),indent=1))
    for it in W[:3]: print("WIT",json.dumps(it))
if __name__=="__main__": main()
