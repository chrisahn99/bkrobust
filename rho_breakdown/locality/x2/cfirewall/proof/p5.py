"""Full census at p: verify L1 (Meek soundness), L2 (D+ in [H]), E7 (O*_H = O(D+)),
   coverage of THEOREM A, forb_D <= forb_H, and the residual size."""
from common import *
from collections import defaultdict
import sys, json, time
from itertools import combinations

def ucolliders(W):
    p=W.shape[0]; out=set()
    for v in range(p):
        pars=[u for u in range(p) if is_directed(W,u,v)]
        for u,w in combinations(sorted(pars),2):
            if not adjacent(W,u,w): out.add((u,v,w))
    return out

def scan_cpdag(C):
    p=C.shape[0]; ref=v_structures(C)
    acc=defaultdict(int); wit=[]
    for G0 in reachable_mpdags(C):
        NA=nonadjacent_pairs(G0)
        if not NA: continue
        Q=[]
        for x in range(p):
            for y in range(p):
                if x==y: continue
                O0,n0,am0=X.ostar_and_paths(G0,x,y)
                if am0: Q.append((x,y,O0))
        if not Q: continue
        ext=consistent_dag_extensions(G0,ref_vstructs=ref)
        for (a0,b0) in NA:
            for (a,b) in ((a0,b0),(b0,a0)):
                H,info=bk_assert(G0,[(a,b)])
                if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
                P=G0.copy(); P[a,b]=1; P[b,a]=0
                ucP=ucolliders(P); ucH=ucolliders(H); dirH=directed_edges(H)
                # precompute per-D structural facts (independent of query)
                info_D=[]
                for D in ext:
                    Dp=D.copy(); Dp[a,b]=1
                    cyc=has_directed_cycle(Dp)
                    c2 = (not cyc) and (ucolliders(Dp)<=ucP)
                    pres = all(is_directed(Dp,i,j) for (i,j) in dirH) if not cyc else False
                    inext = pres and (not cyc) and (ucolliders(Dp)<=ucH)
                    presD = all(is_directed(D,i,j) for (i,j) in dirH if (i,j)!=(a,b))
                    info_D.append((D,Dp,cyc,c2,pres,inext,presD))
                    # L1 : (not cyc) and c2  ==>  pres
                    if (not cyc) and c2 and not pres: acc["L1_FAIL"]+=1
                    if (not cyc) and c2 and pres and not inext: acc["L2_FAIL"]+=1
                for (x,y,O0) in Q:
                    O1,n1,am1=X.ostar_and_paths(H,x,y)
                    if not am1: continue
                    fbH=adjust.forb(H,x,y); cnH=adjust.causal_nodes(H,x,y)
                    moved=(O1!=O0)
                    for (D,Dp,cyc,c2,pres,inext,presD) in info_D:
                        acc["trialD"]+=1
                        thmA = (not cyc) and c2
                        acc["thmA" if thmA else ("res_same" if not moved else "res_moved")]+=1
                        if not thmA:
                            acc["res_cyc" if cyc else "res_newuc"]+=1
                            if presD: acc["res_presD"]+=1
                        if inext:
                            Od=adjust.optimal_adjustment_set(Dp,x,y)
                            if Od is None: acc["E7_none"]+=1
                            elif Od==O1: acc["E7_eq"]+=1
                            else: acc["E7_NE"]+=1
                        fbD=adjust.forb(D,x,y); cnD=adjust.causal_nodes(D,x,y)
                        if not (fbD<=fbH): acc["forbD_NOTSUB_fbH"]+=1
                        if not (cnD<=cnH): acc["cnD_notsub"]+=1
                        if O1 & fbD: acc["CLAUSE1_FAIL"]+=1
                        if not cnD: acc["cnD_EMPTY"]+=1
    return dict(acc),wit

def _job(kb_p):
    kb,p=kb_p
    return scan_cpdag(np.frombuffer(kb,dtype=np.int8).reshape(p,p).copy())

def main():
    p=int(sys.argv[1]); nw=int(sys.argv[2]) if len(sys.argv)>2 else 10
    cp={}
    for D in all_dags(p):
        C=dag_to_cpdag(D); cp.setdefault(C.tobytes(),C)
    print(f"p={p}: {len(cp)} CPDAGs",flush=True)
    tot=defaultdict(int); t0=time.time(); done=0
    from multiprocessing import Pool
    with Pool(nw) as pool:
        for acc,w in pool.imap_unordered(_job,[(k,p) for k in cp],chunksize=4):
            for kk,vv in acc.items(): tot[kk]+=vv
            done+=1
            if done % max(1,len(cp)//20)==0:
                print(f"  {done}/{len(cp)} {time.time()-t0:.0f}s trialD={tot['trialD']}",flush=True)
    print(json.dumps(dict(sorted(tot.items())),indent=1))
    json.dump(dict(tot),open(f"p{p}_summary.json","w"),indent=1)
if __name__ == '__main__':
    main()
