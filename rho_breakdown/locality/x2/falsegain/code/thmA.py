"""THEOREM A check, exhaustive: for EVERY reachable MPDAG H at p<=5, every query,
every consistent DAG extension D in [H]:  R(H,x,y) is SOUND against D."""
import sys
from collections import defaultdict
from multiprocessing import Pool
sys.path.insert(0,"/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/falsegain/code")
from fglib import *
import fglib, numpy as np
from graphs import consistent_dag_extensions

def job(arg):
    Cb,p=arg
    C=np.frombuffer(Cb,dtype=np.int8).reshape(p,p).copy()
    acc=defaultdict(int); bad=[]
    for H in reachable_mpdags(C):
        ext=consistent_dag_extensions(H)
        for x in range(p):
            for y in range(p):
                if x==y: continue
                s,O=report_state(H,x,y)
                for D in ext:
                    sd=sound(D,x,y,s,O)
                    acc[(s,sd)]+=1
                    if sd is False and len(bad)<5:
                        bad.append((H.tolist(),D.tolist(),x,y,s))
    fglib._valid_cache.clear(); fglib._cn_cache.clear()
    return dict(acc),bad

if __name__=="__main__":
    p=int(sys.argv[1]); nw=int(sys.argv[2])
    from run_lemma import all_dags
    cp={}
    for D in all_dags(p):
        C=dag_to_cpdag(D); cp.setdefault(C.tobytes(),C)
    keys=list(cp.keys()); print(f"p={p} {len(keys)} CPDAGs",flush=True)
    tot=defaultdict(int); bads=[]
    with Pool(nw) as pool:
        for acc,b in pool.imap_unordered(job,[(k,p) for k in keys],chunksize=4):
            for k,v in acc.items(): tot[k]+=v
            bads+=b
    for k,v in sorted(tot.items(),key=lambda z:str(z[0])): print("|".join(map(str,k)),v)
    print("VIOLATIONS", len(bads))
