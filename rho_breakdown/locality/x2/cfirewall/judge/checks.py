"""Judge's verification battery."""
import sys, json, time
from collections import defaultdict
from multiprocessing import Pool
import numpy as np
from jlib import *

def scan(C):
    p=C.shape[0]; ref=v_structures(C); acc=defaultdict(int); wit=[]
    for G0 in reachable_mpdags(C):
        NA=nonadjacent_pairs(G0)
        if not NA: continue
        QA=[];QN=[]
        for x in range(p):
            for y in range(p):
                if x==y: continue
                O0,n0,am0,cn0,fb0=ostar(G0,x,y)
                (QA if am0 else QN).append((x,y,O0,cn0,fb0,n0))
        if not (QA or QN): continue
        ext=consistent_dag_extensions(G0,ref_vstructs=ref)
        if not ext: continue
        # (A) HPM on G0 itself
        for (x,y,O0,cn0,fb0,n0) in QA:
            for D in ext:
                acc["hpm_check"]+=1
                if not is_valid_adjustment_set(D,x,y,set(O0)): acc["HPM_FAIL"]+=1
        for (a0,b0) in NA:
            for (a,b) in ((a0,b0),(b0,a0)):
                H,info=bk_assert(G0,[(a,b)])
                coh = (not info["conflict"]) and (not has_directed_cycle(H)) and pdag_extendable(H)
                pdX0=None
                for (x,y,O0,cn0,fb0,n0) in QA:
                    O1,n1,am1,cn1,fb1=ostar(H,x,y)
                    if coh:
                        if not am1: continue
                        acc["report"]+=1
                        S1=set(O1); moved=(O1!=O0)
                        pdX0=poss_de(G0,{x}); pdX1=poss_de(H,{x})
                        Mfail = not (pdX0<=pdX1)
                        acc["Mfail"]+=Mfail
                        if Mfail:
                            acc["Mfail_moved"]+=moved            # (C)
                        # (B) Lemma U0
                        u0=(parents_of_set(H,cn1)&pdX0)-(cn1|{x})
                        if u0:
                            acc["U0_FAIL"]+=1
                            if len(wit)<10: wit.append(dict(tag="U0",G0=G0.tolist(),H=H.tolist(),x=x,y=y,a=a,b=b,bad=sorted(u0)))
                        acc["U0_nonvac"]+= bool(parents_of_set(H,cn1)&pdX0)
                        nb=0
                        for D in ext:
                            v=is_valid_adjustment_set(D,x,y,S1)
                            acc["checks"]+=1
                            if not v: nb+=1; acc["VIOL"]+=1
                            # (F) back-door form
                            deX=poss_de(D,{x})
                            if S1&deX: acc["deX_hit"]+=1
                        if nb and nb!=len(ext): acc["split"]+=1
                    else:
                        # (D) control arm
                        if not am1: continue
                        acc["ctrl_report"]+=1
                        S1=set(O1); nb=0
                        for D in ext:
                            if not is_valid_adjustment_set(D,x,y,S1): nb+=1
                        acc["ctrl_checks"]+=len(ext)
                        if nb:
                            acc["ctrl_viol"]+=1
                            acc["ctrl_viol_allD" if nb==len(ext) else "ctrl_viol_someD"]+=1
                            if O1!=O0: acc["ctrl_viol_moved"]+=1
                # (E) sharpness: G0 NOT amenable, H amenable & coherent -> report anyway
                if coh:
                    for (x,y,O0,cn0,fb0,n0) in QN:
                        O1,n1,am1,cn1,fb1=ostar(H,x,y)
                        if not am1: continue
                        acc["noamen_report"]+=1
                        S1=set(O1); nb=0
                        for D in ext:
                            if not is_valid_adjustment_set(D,x,y,S1): nb+=1
                        acc["noamen_checks"]+=len(ext)
                        if nb:
                            acc["noamen_viol"]+=1
                            if len(wit)<10 and n0>0: wit.append(dict(tag="NOAMEN",G0=G0.tolist(),H=H.tolist(),x=x,y=y,a=a,b=b,O1=sorted(O1)))
    return dict(acc),wit

def main():
    p=int(sys.argv[1]); nw=int(sys.argv[2])
    Cs=cpdags(p); acc=defaultdict(int); W=[]; t0=time.time()
    with Pool(nw) as pool:
        for a,w in pool.imap_unordered(scan,Cs,chunksize=2):
            for k,v in a.items(): acc[k]+=v
            W.extend(w[:2])
    acc["elapsed"]=round(time.time()-t0,1)
    print("p=%d"%p,json.dumps(dict(sorted(acc.items())),indent=1))
    if W: json.dump(W[:20],open("checkwit_p%d.json"%p,"w"),indent=1)
if __name__=="__main__": main()
