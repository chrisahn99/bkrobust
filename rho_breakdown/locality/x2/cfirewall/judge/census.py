import sys, json, time
from collections import defaultdict
from multiprocessing import Pool
import numpy as np
from jlib import *

def scan_cpdag(C):
    p = C.shape[0]
    ref = v_structures(C)
    acc = defaultdict(int)
    wit = []
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
        if not ext:
            acc["ext_empty"] += 1; continue
        # per-D cached cn/forb
        Dinfo = []
        for D in ext:
            Dinfo.append(D)
        for (a0,b0) in NA:
            for (a,b) in ((a0,b0),(b0,a0)):
                acc["stmt_all"] += 1
                H, info = bk_assert(G0, [(a,b)])
                c1 = not info["conflict"]; c2 = not has_directed_cycle(H); c3 = pdag_extendable(H)
                acc["c1_fail"] += (not c1); acc["c2_fail"] += (not c2); acc["c3_fail"] += (not c3)
                if c2 and not c3: acc["acyc_notext"] += 1
                if not (c1 and c2 and c3):
                    acc["stmt_rej"] += 1
                    continue
                acc["stmt_ok"] += 1
                # H is Meek-closed?
                if not np.array_equal(H, meek_closure(H)): acc["H_not_meekclosed"] += 1
                # H == mc(G0 + a->b)?
                Gp = G0.copy(); Gp[a,b]=1; Gp[b,a]=0
                if not np.array_equal(H, meek_closure(Gp)): acc["restamp_not_noop"] += 1
                # H is an MPDAG (= common orientation graph of its own extensions)?
                extH = consistent_dag_extensions(H, ref_vstructs=v_structures(H))
                if not extH: acc["H_ext_empty_BUG"] += 1
                else:
                    cog = common_orientation_graph(extH, skeleton(H))
                    if not np.array_equal(cog, H): acc["H_not_maximal"] += 1
                for (x,y,O0,cn0,fb0) in Q:
                    O1, n1, am1, cn1, fb1 = ostar(H, x, y)
                    acc["trial"] += 1
                    if not am1:
                        acc["abort"] += 1
                        acc["abort_nopath"] += (n1 == 0)
                        continue
                    acc["report"] += 1
                    moved = (O1 != O0)
                    acc["moved"] += moved
                    # unconditional structural facts
                    if O1 & fb0: acc["F1_fail"] += 1           # O*(H) n forb(G0) != {}
                    pdX_G0 = poss_de(G0, {x}); pdX_H = poss_de(H, {x})
                    if not (pdX_G0 <= pdX_H): acc["M_fail"] += 1
                    # Lemma PD on H itself
                    if cn1 and (parents_of_set(H, cn1) & pdX_H) - fb1: acc["PD_fail_H"] += 1
                    nbad = 0; nlift=0; ncertB=0; nres=0; nres_bad=0
                    for D in ext:
                        v = is_valid_adjustment_set(D, x, y, set(O1))
                        ok_lift, Dp, why = lifts(D, a, b, H)
                        cnD = causal_nodes(D, x, y)
                        fbD = (poss_de(D, cnD) | {x}) if cnD else {x}
                        certB = (not (O1 & fbD)) and cnD and ((parents_of_set(D, cnD) - fbD) <= set(O1))
                        acc["checks"] += 1
                        acc["lift_ok"] += ok_lift
                        acc["lift_fail_"+why] += (not ok_lift)
                        acc["certB_ok"] += bool(certB)
                        cov = ok_lift or bool(certB)
                        acc["covered"] += cov
                        if not cov:
                            acc["residual"] += 1
                            nres += 1
                            if not v: nres_bad += 1
                        if O1 & fbD: acc["clause1_fail"] += 1
                        if not v:
                            nbad += 1
                            acc["violation"] += 1
                            if len(wit) < 20:
                                wit.append(dict(G0=G0.tolist(), H=H.tolist(), D=D.tolist(),
                                                x=x,y=y,a=a,b=b,O1=sorted(O1),O0=sorted(O0) if O0 is not None else None,
                                                lift=ok_lift, certB=bool(certB)))
                    if nbad and nbad != len(ext): acc["split_validity"] += 1
    return dict(acc), wit

def main():
    p = int(sys.argv[1]); nw = int(sys.argv[2]) if len(sys.argv)>2 else 8
    Cs = cpdags(p)
    t0=time.time()
    acc = defaultdict(int); WIT=[]
    with Pool(nw) as pool:
        for a, w in pool.imap_unordered(scan_cpdag, Cs, chunksize=4):
            for k,v in a.items(): acc[k]+=v
            WIT.extend(w[:5])
    acc["ncpdag"]=len(Cs); acc["elapsed"]=round(time.time()-t0,1)
    print(json.dumps(dict(sorted(acc.items())), indent=1))
    if WIT:
        json.dump(WIT[:20], open("wit_p%d.json"%p,"w"), indent=1)
        print("WITNESSES WRITTEN", len(WIT))

if __name__=="__main__": main()
