"""Verify the whole proof chain, treatment AND control arm.

A  (P1)  Z n forb(G0) = {}
B  (L4)  every proper non-causal definite-status path in G0 has a collider in
         forb(G0)  OR  its canonical node w_j (j = max i : w_i not in cn(G0))
         satisfies w_j -> w_{j+1} in G0 and w_j in O*(G0)
C  (KC)  if such a path has NO collider in forb(G0) then w_j in Z
D  (GAC) Z satisfies the MPDAG adjustment criterion in G0     [= truth]
E  containments cn0/cnH, F0/FH
"""
import sys, json
from collections import defaultdict
import numpy as np
BASE="/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious"
sys.path.insert(0,BASE+"/x2/code"); sys.path.insert(0,BASE+"/x2/cfirewall/code")
sys.path.insert(0,"/private/tmp/claude-501/-Users-josecosta-mugango/bceac0a4-93a7-4757-b3b3-a3d4e6c42d55/scratchpad/dt")
from graphs import dag_to_cpdag, has_directed_cycle, v_structures, is_directed
import adjust, x2lib as X
from run_lemma import all_dags, reachable_mpdags
from gac import all_paths, node_status, is_causal, gac
import cfirewall as CF
x1=CF.x1_ops

def scan(C):
    acc=defaultdict(int); wit=defaultdict(list)
    for G0 in reachable_mpdags(C):
        p=G0.shape[0]
        NA=x1.nonadjacent_pairs(G0)
        if not NA: continue
        Q=[]
        for x in range(p):
            for y in range(p):
                if x==y: continue
                O0,n0,am0=X.ostar_and_paths(G0,x,y)
                if not am0: continue
                F0=adjust.forb(G0,x,y); CN0=adjust.causal_nodes(G0,x,y)
                P=[]
                for path in all_paths(G0,x,y):
                    if len(path)<2: continue
                    st=[node_status(G0,path[i-1],path[i],path[i+1]) for i in range(1,len(path)-1)]
                    if any(s=="undef" for s in st) or is_causal(G0,path): continue
                    colF=any(s=="col" and path[i] in F0 for i,s in enumerate(st,1))
                    j=max(i for i in range(len(path)) if path[i] not in CN0)
                    P.append((path,st,colF,j))
                    acc["L4_fail"] += (not colF) and not (j+1<len(path) and is_directed(G0,path[j],path[j+1]) and path[j] in O0)
                Q.append((x,y,O0,F0,CN0,P))
        if not Q: continue
        for (a0,b0) in NA:
            for (a,b) in ((a0,b0),(b0,a0)):
                H,info=x1.bk_assert(G0,[(a,b)])
                coh=(not info["conflict"]) and (not has_directed_cycle(H)) and x1.pdag_extendable(H)
                arm="T" if coh else "K"
                for (x,y,O0,F0,CN0,P) in Q:
                    Z,n1,am1=X.ostar_and_paths(H,x,y)
                    if not am1: acc[arm+"_abort"]+=1; continue
                    acc[arm+"_rep"]+=1
                    CNH=adjust.causal_nodes(H,x,y); FH=adjust.forb(H,x,y)
                    A = not (set(Z)&F0)
                    acc[arm+"_A_fail"] += (not A)
                    Cok=True
                    for (path,st,colF,j) in P:
                        if colF: continue
                        if path[j] not in Z:
                            Cok=False; acc[arm+"_C_fail_path"]+=1
                            if len(wit[arm+"C"])<5:
                                wit[arm+"C"].append(dict(G0=G0.tolist(),H=H.tolist(),a=int(a),b=int(b),
                                  x=int(x),y=int(y),path=list(map(int,path)),st=st,j=j,
                                  Z=sorted(map(int,Z)),O0=sorted(map(int,O0)),
                                  cn0=sorted(map(int,CN0)),cnH=sorted(map(int,CNH)),
                                  F0=sorted(map(int,F0)),FH=sorted(map(int,FH))))
                    acc[arm+"_C_fail"] += (not Cok)
                    D = gac(G0,x,y,Z,"poss")
                    acc[arm+"_D_fail"] += (not D)
                    acc[(arm,"AC_vs_D",int(A and Cok),int(D))]+=1
                    acc[arm+"_cn0_sub_cnH"] += set(CN0)<=set(CNH)
                    acc[arm+"_cnH_sub_cn0"] += set(CNH)<=set(CN0)
                    acc[arm+"_F0_sub_FH"]   += set(F0)<=set(FH)
                    acc[arm+"_cn0_sub_FH"]  += set(CN0)<=set(FH)
                    if not set(F0)<=set(FH) and len(wit[arm+"F"])<5:
                        wit[arm+"F"].append(dict(G0=G0.tolist(),H=H.tolist(),a=int(a),b=int(b),x=int(x),y=int(y),
                          Z=sorted(map(int,Z)),F0=sorted(map(int,F0)),FH=sorted(map(int,FH)),
                          cn0=sorted(map(int,CN0)),cnH=sorted(map(int,CNH))))
    return {str(k):v for k,v in acc.items()},{k:v for k,v in wit.items()}

def _job(arg):
    kb,p=arg
    return scan(np.frombuffer(kb,dtype=np.int8).reshape(p,p).copy())

if __name__=="__main__":
    p=int(sys.argv[1]); cp={}
    for D in all_dags(p): C=dag_to_cpdag(D); cp.setdefault(C.tobytes(),C)
    from multiprocessing import Pool
    tot=defaultdict(int); W=defaultdict(list)
    with Pool(8) as pool:
        for acc,w in pool.imap_unordered(_job,[(k,p) for k in cp],chunksize=4):
            for k,v in acc.items(): tot[k]+=v
            for k,v in w.items():
                if len(W[k])<5: W[k].extend(v[:5-len(W[k])])
    for k in sorted(tot,key=str): print(f"{k:35s} {tot[k]}")
    json.dump(dict(counts={k:int(v) for k,v in tot.items()},wit={k:v for k,v in W.items()}),open(f"verify_p{p}.json","w"))
