"""Sharpness of the BIG CONJECTURE at p=4, exhaustive over MPDAG pairs.
 arm A : dir(G0) subseteq dir(H), G0 amenable, H amenable          -> the conjecture
 arm B : H CONTRADICTS G0 on some edge (drop orientation-monotonicity)
 arm C : dir(G0) subseteq dir(H) but G0 NOT amenable
 arm D : dir(G0) subseteq dir(H), G0 amenable, H NOT amenable (report anyway)"""
import sys; sys.path.insert(0,'/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/cfirewall/proof')
from common import *
from itertools import product
from collections import defaultdict
import numpy as np, json, time
from multiprocessing import Pool
p=4
def enum_mpdags(p):
    idx=[(i,j) for i in range(p) for j in range(i+1,p)]; out=[]
    for state in product([0,1,2,3],repeat=len(idx)):
        H=np.zeros((p,p),dtype=np.int8)
        for (i,j),s in zip(idx,state):
            if s==1: H[i,j]=1
            elif s==2: H[j,i]=1
            elif s==3: H[i,j]=H[j,i]=1
        if has_directed_cycle(H): continue
        if not np.array_equal(meek_closure(H),H): continue
        if not pdag_extendable(H): continue
        out.append(H)
    return out
MP=enum_mpdags(p)
def ostar_raw(G,x,y):
    paths=X.pcp_capped(G,x,y)
    if not paths: return None,False
    amen=all(is_directed(G,q[0],q[1]) for q in paths)
    cnv=set()
    for q in paths: cnv.update(q[1:])
    fb=(adjust.poss_de(G,cnv)|{x}) if cnv else {x}
    return frozenset(adjust.parents_of_set(G,cnv)-fb),amen
def job(gi):
    G0=MP[gi]; acc=defaultdict(int); wit={}
    ext=consistent_dag_extensions(G0,ref_vstructs=v_structures(G0))
    if not ext: return dict(acc),wit
    QQ=[]
    for x in range(p):
        for y in range(p):
            if x==y: continue
            O0,a0=ostar_raw(G0,x,y)
            if O0 is not None: QQ.append((x,y,O0,a0))
    if not QQ: return dict(acc),wit
    skG=skeleton(G0); dirG=directed_edges(G0)
    for H in MP:
        if not np.all(skeleton(H)>=skG): continue
        mono=all(is_directed(H,u,v) for (u,v) in dirG)
        contra=any(is_directed(H,v,u) for (u,v) in dirG)
        for (x,y,O0,a0) in QQ:
            OH,aH=ostar_raw(H,x,y)
            if OH is None: continue
            arm = ("A" if (mono and a0 and aH) else "B" if (contra and a0 and aH) else
                   "C" if (mono and (not a0) and aH) else "D" if (mono and a0 and (not aH)) else None)
            if arm is None: continue
            acc[arm+"_trial"]+=1
            for D in ext:
                acc[arm+"_check"]+=1
                if not adjust.is_valid_adjustment_set(D,x,y,set(OH)):
                    acc[arm+"_VIOL"]+=1
                    if arm not in wit:
                        wit[arm]=dict(G0=G0.tolist(),H=H.tolist(),D=D.tolist(),x=int(x),y=int(y),
                                      OH=sorted(int(z) for z in OH),O0=sorted(int(z) for z in O0))
    return dict(acc),wit
def main():
    tot=defaultdict(int); W={}; t0=time.time()
    with Pool(7) as pool:
        for acc,w in pool.imap_unordered(job,range(len(MP)),chunksize=4):
            for k,v in acc.items(): tot[k]+=v
            for k,v in w.items(): W.setdefault(k,v)
    print("MPDAGs",len(MP),"%.0fs"%(time.time()-t0))
    print(json.dumps(dict(sorted(tot.items())),indent=1))
    for k in sorted(W): print("WIT",k,json.dumps(W[k]))
if __name__=="__main__": main()
