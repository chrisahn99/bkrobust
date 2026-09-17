"""Fixed k, random statements drawn EITHER all from non-adjacent pairs (the C-FIREWALL pool)
   OR all from existing undirected edges (ordinary background knowledge).  Base = CPDAG C."""
import sys; sys.path.insert(0,'/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/cfirewall/proof')
from common import *
from collections import defaultdict
import numpy as np, json, time
from multiprocessing import Pool

def ostar_raw(G,x,y):
    paths=X.pcp_capped(G,x,y)
    if not paths: return None,0,False
    amen=all(is_directed(G,q[0],q[1]) for q in paths)
    cnv=set()
    for q in paths: cnv.update(q[1:])
    fb=(adjust.poss_de(G,cnv)|{x}) if cnv else {x}
    return frozenset(adjust.parents_of_set(G,cnv)-fb),len(paths),amen

def job(args):
    p,n,seed,K=args
    rng=np.random.default_rng(seed); acc=defaultdict(int); wit={}
    for _ in range(n):
        deg=float(rng.choice([2.0,2.5,3.0,3.5,4.0,5.0,6.0]))
        D=random_dag(p,deg,rng); C=dag_to_cpdag(D); ref=v_structures(C)
        if len(undirected_edges(C))>13: continue
        ext=consistent_dag_extensions(C,ref_vstructs=ref)
        if not ext: continue
        QA=[]
        for x in range(p):
            for y in range(p):
                if x==y: continue
                O0,n0,a0=ostar_raw(C,x,y)
                if O0 is not None and a0: QA.append((x,y,O0))
        if not QA: continue
        NA=nonadjacent_pairs(C); UN=undirected_edges(C)
        for tag,pairs in (("NA",NA),("UND",UN)):
            if len(pairs)<K: continue
            for _rep in range(6):
                sel=rng.choice(len(pairs),size=K,replace=False)
                stmts=[]
                for s in sel:
                    u,v=pairs[int(s)]
                    stmts.append((u,v) if rng.random()<0.5 else (v,u))
                H,info=bk_assert(C,stmts)
                if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
                for (x,y,O0) in QA:
                    OH,nh,amH=ostar_raw(H,x,y)
                    if OH is None or not amH: continue
                    acc[tag+"_trial"]+=1
                    if OH!=O0: acc[tag+"_moved"]+=1
                    for Dd in ext:
                        acc[tag+"_check"]+=1
                        if not adjust.is_valid_adjustment_set(Dd,x,y,set(OH)):
                            acc[tag+"_VIOL"]+=1
                            if tag not in wit:
                                wit[tag]=dict(p=p,K=K,C=C.tolist(),H=H.tolist(),D=Dd.tolist(),
                                              x=int(x),y=int(y),stmts=[[int(i),int(j)] for i,j in stmts],
                                              OH=sorted(int(z) for z in OH),O0=sorted(int(z) for z in O0))
    return dict(acc),wit

def main():
    p=int(sys.argv[1]); nw=int(sys.argv[2]); n=int(sys.argv[3]); K=int(sys.argv[4]); seed=int(sys.argv[5]) if len(sys.argv)>5 else 11
    nb=80; jobs=[(p,max(1,n//nb),seed+7919*i,K) for i in range(nb)]
    tot=defaultdict(int); W={}; t0=time.time()
    with Pool(nw) as pool:
        for acc,w in pool.imap_unordered(job,jobs,chunksize=1):
            for k,v in acc.items(): tot[k]+=v
            for k,v in w.items(): W.setdefault(k,v)
    print("p=%d K=%d  %s  (%.0fs)"%(p,K,json.dumps(dict(sorted(tot.items()))),time.time()-t0),flush=True)
    for k in sorted(W): print("WIT",k,json.dumps(W[k]),flush=True)
if __name__=="__main__": main()
