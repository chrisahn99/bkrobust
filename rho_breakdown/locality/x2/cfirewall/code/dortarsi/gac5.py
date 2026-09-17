import sys, itertools, json
from collections import defaultdict
import numpy as np
sys.path.insert(0,"/private/tmp/claude-501/-Users-josecosta-mugango/bceac0a4-93a7-4757-b3b3-a3d4e6c42d55/scratchpad/dt")
from gac import gac, truth, all_dags, reachable_mpdags, consistent_dag_extensions, dag_to_cpdag, v_structures
import x2lib as X
def job(arg):
    kb,p=arg
    C=np.frombuffer(kb,dtype=np.int8).reshape(p,p).copy()
    ref=v_structures(C); acc=defaultdict(int); wit=[]
    for G0 in reachable_mpdags(C):
        ext=consistent_dag_extensions(G0,ref_vstructs=ref)
        if not ext: continue
        for x in range(p):
            for y in range(p):
                if x==y: continue
                O0,n0,am0=X.ostar_and_paths(G0,x,y)
                if not am0: continue
                others=[v for v in range(p) if v not in (x,y)]
                for r in range(len(others)+1):
                    for Z in itertools.combinations(others,r):
                        t=truth(G0,ext,x,y,Z)
                        for var in ("poss","dir"):
                            g=gac(G0,x,y,Z,var)
                            acc[(var,int(t),int(g))]+=1
                            if t!=g and len(wit)<4: wit.append((var,G0.tolist(),x,y,list(Z),t,g))
    return {str(k):v for k,v in acc.items()},wit
if __name__=="__main__":
    p=5; cp={}
    for D in all_dags(p): C=dag_to_cpdag(D); cp.setdefault(C.tobytes(),C)
    from multiprocessing import Pool
    tot=defaultdict(int); W=[]
    with Pool(8) as pool:
        for acc,w in pool.imap_unordered(job,[(k,p) for k in cp],chunksize=4):
            for k,v in acc.items(): tot[k]+=v
            W.extend(w)
    print(json.dumps({k:int(v) for k,v in sorted(tot.items())},indent=1)); print("MISMATCH",W[:2])
