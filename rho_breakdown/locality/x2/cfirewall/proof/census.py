"""Consolidated census: Theorem-1 verification + F1 + sharpness ablations."""
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
    fb = (adjust.poss_de(G,cnv) | {x}) if cnv else {x}
    return frozenset(adjust.parents_of_set(G,cnv)-fb), len(paths), amen

def scan_cpdag(kb_p):
    kb,p = kb_p
    C = np.frombuffer(kb,dtype=np.int8).reshape(p,p).copy()
    ref = v_structures(C); acc=defaultdict(int); wit={}
    def rec(tag,payload):
        if tag not in wit: wit[tag]=payload
    for G0 in reachable_mpdags(C):
        ext = consistent_dag_extensions(G0, ref_vstructs=ref)
        if not ext: continue
        QA=[];QN=[]
        for x in range(p):
            for y in range(p):
                if x==y: continue
                O0,np0,am0 = ostar_raw(G0,x,y)
                if O0 is None: continue
                (QA if am0 else QN).append((x,y,O0,adjust.forb(G0,x,y)))
        if not QA and not QN: continue
        NA = nonadjacent_pairs(G0); UND = undirected_edges(G0)
        vcache={}
        def isvalid(D,x,y,Z,di):
            k=(di,x,y,Z)
            v=vcache.get(k)
            if v is None:
                v=adjust.is_valid_adjustment_set(D,x,y,set(Z)); vcache[k]=v
            return v
        for kind,pairs in (("NA",NA),("UND",UND)):
            for (u,v) in pairs:
                for (a,b) in ((u,v),(v,u)):
                    H,info = bk_assert(G0,[(a,b)])
                    c1=not info["conflict"]; c2=not has_directed_cycle(H); c3=pdag_extendable(H)
                    coherent = c1 and c2 and c3
                    dirH_S=[(s,t) for (s,t) in directed_edges(H) if not (s==a and t==b)]
                    lift=[]
                    if coherent:
                        vH=v_structures(H)
                        for D in ext:
                            if kind=="NA":
                                Dp=D.copy(); Dp[a,b]=1
                                ok=(not has_directed_cycle(Dp)) and all(is_directed(Dp,s,t) for (s,t) in dirH_S) and v_structures(Dp)==vH
                            else:
                                ok=all(is_directed(D,s,t) for (s,t) in directed_edges(H))
                            lift.append(ok)
                    for (amenG0,QS) in ((True,QA),(False,QN)):
                        for (x,y,O0,fb0) in QS:
                            OH,nph,amH = ostar_raw(H,x,y)
                            if OH is None: continue
                            if kind=="NA" and coherent and amH and amenG0:
                                acc["main_trial"]+=1
                                if OH!=O0: acc["main_moved"]+=1
                                if OH & fb0:
                                    acc["F1_viol"]+=1; rec("F1",(C.tolist(),G0.tolist(),x,y,a,b,sorted(OH),sorted(fb0)))
                            for di,D in enumerate(ext):
                                ok = isvalid(D,x,y,OH,di)
                                fbD=adjust.forb(D,x,y)
                                if kind=="NA" and coherent and amH and amenG0:
                                    acc["main_check"]+=1
                                    if not (fbD<=fb0): acc["fbD_notsub_fb0"]+=1
                                    if not adjust.causal_nodes(D,x,y): acc["cnD_empty"]+=1
                                    if not ok:
                                        acc["MAIN_VIOL"]+=1; rec("MAIN",(C.tolist(),G0.tolist(),D.tolist(),x,y,a,b,sorted(OH)))
                                # Theorem 1: liftable => valid (any query, amenable G0 or not)
                                if coherent and amH and lift and lift[di]:
                                    acc["T1_%s_check"%kind]+=1
                                    if not amenG0: acc["T1_%s_check_G0nonamen"%kind]+=1
                                    if not ok:
                                        acc["T1_%s_VIOL"%kind]+=1
                                        rec("T1"+kind,(C.tolist(),G0.tolist(),D.tolist(),x,y,a,b,sorted(OH)))
                                if coherent and amH and lift and not lift[di] and kind=="NA" and amenG0:
                                    acc["nolift_check"]+=1
                                    if not ok: acc["nolift_VIOL"]+=1
                                # ablations
                                if kind=="UND" and coherent and amH and amenG0:
                                    acc["abl1_check"]+=1
                                    if not ok:
                                        acc["abl1_viol"]+=1; rec("abl1",(C.tolist(),G0.tolist(),D.tolist(),x,y,a,b,sorted(OH),sorted(O0)))
                                if kind=="NA" and coherent and (not amH) and amenG0:
                                    acc["abl2_check"]+=1
                                    if not ok:
                                        acc["abl2_viol"]+=1; rec("abl2",(C.tolist(),G0.tolist(),D.tolist(),x,y,a,b,sorted(OH)))
                                if kind=="NA" and (not coherent) and amH and amenG0:
                                    acc["abl3_check"]+=1
                                    if not ok:
                                        acc["abl3_viol"]+=1; rec("abl3",(C.tolist(),G0.tolist(),D.tolist(),x,y,a,b,sorted(OH),[c1,c2,c3]))
                                if kind=="NA" and coherent and amH and (not amenG0):
                                    acc["abl4_check"]+=1
                                    if not ok:
                                        acc["abl4_viol"]+=1; rec("abl4",(C.tolist(),G0.tolist(),D.tolist(),x,y,a,b,sorted(OH)))
    return dict(acc), wit

def main():
    p=int(sys.argv[1]); nw=int(sys.argv[2]) if len(sys.argv)>2 else 8
    lim=int(sys.argv[3]) if len(sys.argv)>3 else 0
    cps={}
    for D in all_dags(p):
        C=dag_to_cpdag(D); cps.setdefault(C.tobytes(),C)
    jobs=[(k,p) for k in cps]
    if lim: jobs=jobs[:lim]
    print("cpdags",len(jobs),flush=True)
    tot=defaultdict(int); W={}; t0=time.time(); done=0
    with Pool(nw) as pool:
        for acc,w in pool.imap_unordered(scan_cpdag,jobs,chunksize=2):
            for k,v in acc.items(): tot[k]+=v
            for k,v in w.items(): W.setdefault(k,v)
            done+=1
            if done%max(1,len(jobs)//20)==0: print("  %d/%d %.0fs main=%d viol=%d"%(done,len(jobs),time.time()-t0,tot["main_trial"],tot["MAIN_VIOL"]),flush=True)
    out=dict(p=p,counts={k:int(v) for k,v in sorted(tot.items())},witnesses=W,elapsed=time.time()-t0)
    with open("census_p%d.json"%p,"w") as f: json.dump(out,f,indent=1)
    print(json.dumps(out["counts"],indent=1))
    for k in sorted(W): print("WIT",k,json.dumps(W[k])[:400])
if __name__=="__main__": main()
