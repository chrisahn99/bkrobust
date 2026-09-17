"""BIG CONJECTURE.  G0 an MPDAG amenable rel (X,Y).  H ANY MPDAG with
   skeleton(H) >= skeleton(G0) and dir(H) >= dir(G0)  (an orientation-monotone refinement,
   possibly with extra edges), H amenable rel (X,Y).
   Then is O*(X,Y,H) valid for (X,Y) in EVERY D in [G0]?"""
import sys; sys.path.insert(0,'.')
from common import *
from itertools import product
from collections import defaultdict
import numpy as np, json, time
from multiprocessing import Pool
p=int(sys.argv[1]); nw=int(sys.argv[2]) if len(sys.argv)>2 else 8

def enum_mpdags(p):
    idx=[(i,j) for i in range(p) for j in range(i+1,p)]
    out=[]
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

def job(gi):
    G0=MP[gi]; acc=defaultdict(int); wit=[]
    ext=consistent_dag_extensions(G0,ref_vstructs=v_structures(G0))
    if not ext: return dict(acc),wit
    QA=[]
    for x in range(p):
        for y in range(p):
            if x==y: continue
            O0=adjust.optimal_adjustment_set(G0,x,y)
            if O0 is not None: QA.append((x,y,O0))
    if not QA: return dict(acc),wit
    skG=skeleton(G0); dirG=directed_edges(G0)
    for H in MP:
        if not np.all(skeleton(H)>=skG): continue
        if not all(is_directed(H,u,v) for (u,v) in dirG): continue
        nadd=int((skeleton(H).sum()-skG.sum())//2)
        for (x,y,O0) in QA:
            OH=adjust.optimal_adjustment_set(H,x,y)
            if OH is None: continue
            acc["trial"]+=1; acc["trial_add%d"%min(nadd,3)]+=1
            if OH!=O0: acc["moved"]+=1
            for D in ext:
                acc["check"]+=1
                if not adjust.is_valid_adjustment_set(D,x,y,set(OH)):
                    acc["VIOL"]+=1; acc["VIOL_add%d"%min(nadd,3)]+=1
                    if len(wit)<4: wit.append((G0.tolist(),H.tolist(),D.tolist(),int(x),int(y),sorted(int(z) for z in OH),sorted(int(z) for z in O0)))
    return dict(acc),wit

def main():
  print("MPDAGs on %d nodes: %d"%(p,len(MP)),flush=True)
  tot=defaultdict(int); W=[]; t0=time.time(); done=0
  with Pool(nw) as pool:
      for acc,w in pool.imap_unordered(job,range(len(MP)),chunksize=4):
        for k,v in acc.items(): tot[k]+=v
        W.extend(w[:max(0,4-len(W))]); done+=1
        if done%max(1,len(MP)//10)==0: print("  %d/%d %.0fs check=%d VIOL=%d"%(done,len(MP),time.time()-t0,tot["check"],tot["VIOL"]),flush=True)
  print(json.dumps(dict(sorted(tot.items())),indent=1))
  for it in W: print("WIT",json.dumps(it))

if __name__=="__main__": main()
