"""What does Dor-Tarsi buy?  Structural probe over the full census."""
import sys, json, importlib.util
from collections import defaultdict
import numpy as np
BASE = "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious"
sys.path.insert(0, BASE + "/x2/code"); sys.path.insert(0, BASE + "/x2/cfirewall/code")
from graphs import (consistent_dag_extensions, dag_to_cpdag, dag_agrees_with,
                    has_directed_cycle, meek_closure, random_dag, skeleton,
                    undirected_edges, v_structures, directed_edges)
import adjust, x2lib as X
from run_lemma import all_dags, reachable_mpdags
import cfirewall as CF
x1 = CF.x1_ops

def ancestors(D, v):
    p = D.shape[0]; seen={v}; st=[v]
    while st:
        u=st.pop()
        for w in range(p):
            if D[w,u]==1 and D[u,w]==0 and w not in seen:
                seen.add(w); st.append(w)
    return seen

def scan_cpdag(C):
    ref = v_structures(C); acc = defaultdict(int); wit=[]
    for G0 in reachable_mpdags(C):
        p = G0.shape[0]
        NA = x1.nonadjacent_pairs(G0)
        if not NA: continue
        Q=[]
        for x in range(p):
            for y in range(p):
                if x==y: continue
                O0,n0,am0 = X.ostar_and_paths(G0,x,y)
                if am0: Q.append((x,y,O0,adjust.forb(G0,x,y)))
        if not Q: continue
        ext = consistent_dag_extensions(G0, ref_vstructs=ref)
        if not ext: continue
        for (a0,b0) in NA:
            for (a,b) in ((a0,b0),(b0,a0)):
                H, info = x1.bk_assert(G0,[(a,b)])
                cyc = has_directed_cycle(H); ext_ok = x1.pdag_extendable(H)
                acc["stmt"]+=1
                if info["conflict"]: acc["conflict"]+=1
                if cyc!=(not ext_ok): acc["cycle_ne_notext"]+=1
                if cyc: acc["cyc"]+=1
                if not ext_ok: acc["notext"]+=1
                coherent = (not info["conflict"]) and (not cyc) and ext_ok
                vH = v_structures(H); dirH = directed_edges(H)
                # how does H relate to G0 orientation-wise
                for (x,y,O0,F0) in Q:
                    Z,n1,am1 = X.ostar_and_paths(H,x,y)
                    tag = "T" if coherent else "K"     # treatment / control
                    acc[tag+"_trial"]+=1
                    if not am1: acc[tag+"_abort"]+=1; continue
                    acc[tag+"_report"]+=1
                    moved = (Z!=O0)
                    acc[tag+"_moved"]+= moved
                    if Z & F0: acc[tag+"_hits_forbG0"]+=1
                    for D in ext:
                        acc[tag+"_chk"]+=1
                        ok = adjust.is_valid_adjustment_set(D,x,y,set(Z))
                        if not ok:
                            acc[tag+"_bad"]+=1
                            if len(wit)<8 and tag=="K":
                                wit.append(dict(G0=G0.tolist(),a=int(a),b=int(b),x=int(x),y=int(y),
                                                Z=sorted(map(int,Z)),D=D.tolist(),cyc=bool(cyc)))
                        # D+ analysis
                        anc_a = ancestors(D,a)
                        acyc = b not in anc_a
                        acc[tag+"_ab_acyclic_inD"]+= acyc
                        if acyc:
                            Dp = D.copy(); Dp[a,b]=1; Dp[b,a]=0
                            agree = dag_agrees_with(Dp,H)
                            vsame = (v_structures(Dp)==vH)
                            inH = agree and vsame
                            acc[tag+"_Dp_agree"]+= agree
                            acc[tag+"_Dp_vsame"]+= vsame
                            acc[tag+"_Dp_inH"]+= inH
                            okp = adjust.is_valid_adjustment_set(Dp,x,y,set(Z))
                            acc[tag+"_Dp_valid"]+= okp
                            if okp and not ok: acc[tag+"_MONO_FAIL"]+=1
                        else:
                            acc[tag+"_ab_cyclic_inD"]+=1
                            if not ok: acc[tag+"_bad_when_cyclic"]+=1
    return dict(acc), wit

def main():
    p=int(sys.argv[1]); nw=int(sys.argv[2]) if len(sys.argv)>2 else 8
    cp={}
    for D in all_dags(p):
        C=dag_to_cpdag(D); cp.setdefault(C.tobytes(),C)
    from multiprocessing import Pool
    tot=defaultdict(int); W=[]
    with Pool(nw) as pool:
        for acc,w in pool.imap_unordered(_job,[(k,p) for k in cp],chunksize=4):
            for k,v in acc.items(): tot[k]+=v
            W.extend(w)
    print(json.dumps({k:int(v) for k,v in sorted(tot.items())},indent=1))
    json.dump(dict(counts={k:int(v) for k,v in tot.items()},wit=W[:8]),open(f"/private/tmp/claude-501/-Users-josecosta-mugango/bceac0a4-93a7-4757-b3b3-a3d4e6c42d55/scratchpad/dt/probe1_p{p}.json","w"))

def _job(arg):
    kb,p=arg
    C=np.frombuffer(kb,dtype=np.int8).reshape(p,p).copy()
    return scan_cpdag(C)

if __name__=="__main__": main()
