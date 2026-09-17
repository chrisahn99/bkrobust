"""Hunt sharpness witnesses at p=5 (census) for:
   abl1: {a,b} ADJACENT in G0 (orient an existing undirected edge) -- drop non-adjacency
   abl2: H NOT amenable, report anyway -- drop the visible abort"""
import sys; sys.path.insert(0,'/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/cfirewall/proof')
from common import *
from collections import defaultdict
import numpy as np, json, time
from multiprocessing import Pool

def ostar_raw(G,x,y):
    paths = X.pcp_capped(G,x,y)
    if not paths: return None,0,False
    amen = all(is_directed(G,pth[0],pth[1]) for pth in paths)
    cnv=set()
    for pth in paths: cnv.update(pth[1:])
    fb=(adjust.poss_de(G,cnv)|{x}) if cnv else {x}
    return frozenset(adjust.parents_of_set(G,cnv)-fb), len(paths), amen

def job(kb_p):
    kb,p=kb_p
    C=np.frombuffer(kb,dtype=np.int8).reshape(p,p).copy()
    ref=v_structures(C); acc=defaultdict(int); wit={}
    for G0 in reachable_mpdags(C):
        if len(undirected_edges(G0))>14: continue
        QA=[]
        for x in range(p):
            for y in range(p):
                if x==y: continue
                O0,n0,a0=ostar_raw(G0,x,y)
                if O0 is not None and a0: QA.append((x,y,O0))
        if not QA: continue
        ext=consistent_dag_extensions(G0,ref_vstructs=ref)
        if not ext: continue
        for kind,pairs in (("abl1",undirected_edges(G0)),("abl2",nonadjacent_pairs(G0))):
            for (u,v) in pairs:
                for (a,b) in ((u,v),(v,u)):
                    H,info=bk_assert(G0,[(a,b)])
                    if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
                    for (x,y,O0) in QA:
                        OH,nh,amH=ostar_raw(H,x,y)
                        if OH is None: continue
                        if kind=="abl1" and not amH: continue
                        if kind=="abl2" and amH: continue
                        for D in ext:
                            acc[kind+"_check"]+=1
                            if not adjust.is_valid_adjustment_set(D,x,y,set(OH)):
                                acc[kind+"_viol"]+=1
                                if kind not in wit:
                                    wit[kind]=dict(C=C.tolist(),G0=G0.tolist(),D=D.tolist(),x=int(x),y=int(y),
                                                   a=int(a),b=int(b),Ostar_H=sorted(int(z) for z in OH),
                                                   Ostar_G0=sorted(int(z) for z in O0),H=H.tolist())
    return dict(acc),wit

def main():
    p=int(sys.argv[1]); nw=int(sys.argv[2])
    cps={}
    for D in all_dags(p):
        C=dag_to_cpdag(D); cps.setdefault(C.tobytes(),C)
    jobs=[(k,p) for k in cps]
    tot=defaultdict(int); W={}; t0=time.time(); done=0
    with Pool(nw) as pool:
        for acc,w in pool.imap_unordered(job,jobs,chunksize=2):
            for k,v in acc.items(): tot[k]+=v
            for k,v in w.items(): W.setdefault(k,v)
            done+=1
            if done%max(1,len(jobs)//8)==0: print("  %d/%d %.0fs %s"%(done,len(jobs),time.time()-t0,dict(tot)),flush=True)
    print(json.dumps(dict(sorted(tot.items())),indent=1))
    for k in sorted(W): print("WIT",k,json.dumps(W[k]))
if __name__=="__main__": main()
