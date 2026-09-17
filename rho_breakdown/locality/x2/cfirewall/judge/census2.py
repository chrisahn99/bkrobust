import sys, json, time
from collections import defaultdict
from multiprocessing import Pool
import numpy as np
from jlib import *

MAXW = 40
def scan_cpdag(args):
    C, do_mpdagcheck = args
    p = C.shape[0]
    ref = v_structures(C)
    acc = defaultdict(int); wit=[]; res=[]
    for G0 in reachable_mpdags(C):
        NA = nonadjacent_pairs(G0)
        if not NA: continue
        Q = []
        for x in range(p):
            for y in range(p):
                if x == y: continue
                O0, n0, am0, cn0, fb0 = ostar(G0, x, y)
                if am0: Q.append((x, y, O0, cn0, fb0))
        if not Q: continue
        ext = consistent_dag_extensions(G0, ref_vstructs=ref)
        acc["mpdag"] += 1
        if not ext: acc["ext_empty"] += 1; continue
        Dcn = {}
        for di,D in enumerate(ext):
            for (x,y,_,_,_) in Q:
                if (di,x,y) in Dcn: continue
                cnD = causal_nodes(D,x,y)
                fbD = (poss_de(D,cnD)|{x}) if cnD else {x}
                Dcn[(di,x,y)] = (cnD, fbD, parents_of_set(D,cnD)-fbD if cnD else set())
        vcache = {}
        for (a0,b0) in NA:
            for (a,b) in ((a0,b0),(b0,a0)):
                acc["stmt_all"] += 1
                H, info = bk_assert(G0, [(a,b)])
                c1 = not info["conflict"]; c2 = not has_directed_cycle(H); c3 = pdag_extendable(H)
                acc["c1_fail"] += (not c1); acc["c2_fail"] += (not c2); acc["c3_fail"] += (not c3)
                if c2 and not c3: acc["acyc_notext"] += 1
                if c3 and not c2: acc["ext_cyc"] += 1
                if not (c1 and c2 and c3): acc["stmt_rej"] += 1; continue
                acc["stmt_ok"] += 1
                if do_mpdagcheck:
                    extH = consistent_dag_extensions(H, ref_vstructs=v_structures(H))
                    if not extH: acc["H_ext_empty_BUG"] += 1
                    elif not np.array_equal(common_orientation_graph(extH, skeleton(H)), H):
                        acc["H_not_maximal"] += 1
                    acc["mpdagchecked"] += 1
                liftcache = {}
                pdX_H = poss_de(H, {a for a in []} or set())  # placeholder
                for (x,y,O0,cn0,fb0) in Q:
                    O1, n1, am1, cn1, fb1 = ostar(H, x, y)
                    acc["trial"] += 1
                    if not am1:
                        acc["abort"] += 1; acc["abort_nopath"] += (n1==0); continue
                    acc["report"] += 1
                    moved = (O1 != O0); acc["moved"] += moved
                    if O1 & fb0: acc["F1_fail"] += 1
                    pdX0 = poss_de(G0,{x}); pdX1 = poss_de(H,{x})
                    if not (pdX0 <= pdX1): acc["M_fail"] += 1
                    if cn1 and ((parents_of_set(H,cn1) & pdX1) - fb1): acc["PD_fail_H"] += 1
                    nbad=0
                    for di,D in enumerate(ext):
                        kk=(di,x,y,O1)
                        v = vcache.get(kk)
                        if v is None:
                            v = is_valid_adjustment_set(D,x,y,set(O1)); vcache[kk]=v
                        lk = di
                        lo = liftcache.get(lk)
                        if lo is None:
                            lo = lifts(D,a,b,H)[0::2]; liftcache[lk]=lo
                        ok_lift, why = lo
                        cnD, fbD, paD = Dcn[(di,x,y)]
                        certB = (not (O1 & fbD)) and bool(cnD) and (paD <= set(O1))
                        acc["checks"] += 1
                        acc["lift_ok"] += ok_lift
                        if not ok_lift: acc["lift_fail_"+why]+=1
                        acc["certB_ok"] += certB
                        cov = ok_lift or certB
                        acc["covered"] += cov
                        if not cov:
                            acc["residual"] += 1
                            if not v: acc["residual_bad"] += 1
                            if len(res)<MAXW:
                                res.append(dict(G0=G0.tolist(),H=H.tolist(),D=D.tolist(),x=x,y=y,a=a,b=b,
                                                O1=sorted(O1),valid=bool(v),why=why))
                        if O1 & fbD: acc["clause1_fail"] += 1
                        if not v:
                            nbad+=1; acc["violation"]+=1
                            if len(wit)<MAXW: wit.append(dict(G0=G0.tolist(),H=H.tolist(),D=D.tolist(),
                                                              x=x,y=y,a=a,b=b,O1=sorted(O1)))
                    if nbad and nbad!=len(ext): acc["split_validity"]+=1
    return dict(acc), wit, res

def main():
    p=int(sys.argv[1]); nw=int(sys.argv[2]) if len(sys.argv)>2 else 10
    mp = (len(sys.argv)>3 and sys.argv[3]=="mpdagcheck")
    Cs=cpdags(p); t0=time.time()
    acc=defaultdict(int); WIT=[]; RES=[]
    with Pool(nw) as pool:
        for a,w,r in pool.imap_unordered(scan_cpdag, [(C,mp) for C in Cs], chunksize=2):
            for k,v in a.items(): acc[k]+=v
            WIT.extend(w[:3]); RES.extend(r[:3])
    acc["ncpdag"]=len(Cs); acc["elapsed"]=round(time.time()-t0,1)
    print(json.dumps(dict(sorted(acc.items())),indent=1))
    json.dump(RES[:400], open("res_p%d.json"%p,"w"))
    if WIT: json.dump(WIT[:40], open("wit_p%d.json"%p,"w"),indent=1); print("!!! WITNESSES", len(WIT))

if __name__=="__main__": main()
