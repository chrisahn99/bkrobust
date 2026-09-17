"""k=2 statements: NA+NA vs UND+UND vs mixed.  Does non-adjacency matter?"""
import sys; sys.path.insert(0,'/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/cfirewall/proof')
from common import *
from collections import defaultdict
from itertools import permutations
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

def scan(C,G0,ref,acc,wit,kmax=2):
    p=G0.shape[0]
    if len(undirected_edges(G0))>13: return
    QA=[]
    for x in range(p):
        for y in range(p):
            if x==y: continue
            O0,n0,a0=ostar_raw(G0,x,y)
            if O0 is not None and a0: QA.append((x,y,O0))
    if not QA: return
    ext=consistent_dag_extensions(G0,ref_vstructs=ref)
    if not ext: return
    NA=[(a,b) for (u,v) in nonadjacent_pairs(G0) for (a,b) in ((u,v),(v,u))]
    UN=[(a,b) for (u,v) in undirected_edges(G0) for (a,b) in ((u,v),(v,u))]
    def pool_pairs(P1,P2,tag):
        for s1 in P1:
            for s2 in P2:
                if {s1[0],s1[1]}=={s2[0],s2[1]}: continue
                H,info=bk_assert(G0,[s1,s2])
                if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
                for (x,y,O0) in QA:
                    OH,nh,amH=ostar_raw(H,x,y)
                    if OH is None or not amH: continue
                    acc[tag+"_trial"]+=1
                    if OH!=O0: acc[tag+"_moved"]+=1
                    for D in ext:
                        acc[tag+"_check"]+=1
                        if not adjust.is_valid_adjustment_set(D,x,y,set(OH)):
                            acc[tag+"_VIOL"]+=1
                            if tag not in wit:
                                wit[tag]=dict(C=C.tolist(),G0=G0.tolist(),D=D.tolist(),H=H.tolist(),
                                              x=int(x),y=int(y),s1=[int(z) for z in s1],s2=[int(z) for z in s2],
                                              OH=sorted(int(z) for z in OH),O0=sorted(int(z) for z in O0))
    pool_pairs(NA,NA,"NAxNA"); pool_pairs(UN,UN,"UNxUN"); pool_pairs(NA,UN,"NAxUN")

def job(args):
    p,n,seed=args
    rng=np.random.default_rng(seed); acc=defaultdict(int); wit={}
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
        scan(C,G,ref,acc,wit)
    return dict(acc),wit

def main():
    p=int(sys.argv[1]); nw=int(sys.argv[2]); n=int(sys.argv[3]); seed=int(sys.argv[4]) if len(sys.argv)>4 else 7
    nb=100; jobs=[(p,max(1,n//nb),seed+7919*i) for i in range(nb)]
    tot=defaultdict(int); W={}; t0=time.time(); done=0
    with Pool(nw) as pool:
        for acc,w in pool.imap_unordered(job,jobs,chunksize=1):
            for k,v in acc.items(): tot[k]+=v
            for k,v in w.items(): W.setdefault(k,v)
            done+=1
            if done%10==0: print("  %d/%d %.0fs %s"%(done,nb,time.time()-t0,{k:v for k,v in tot.items() if 'VIOL' in k or 'check' in k}),flush=True)
    print("p",p,json.dumps(dict(sorted(tot.items())),indent=1))
    for k in sorted(W): print("WIT",k,json.dumps(W[k]))
if __name__=="__main__": main()
