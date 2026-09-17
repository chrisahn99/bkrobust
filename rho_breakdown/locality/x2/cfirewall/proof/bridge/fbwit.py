from common import *
import json,sys
def fmt(G): return " ".join([f"{i}->{j}" for (i,j) in directed_edges(G)]+[f"{i}--{j}" for (i,j) in undirected_edges(G)])
p=5; cp={}
for D in all_dags(p):
    C=dag_to_cpdag(D); cp.setdefault(C.tobytes(),C)
n=0
for kb in cp:
    C=cp[kb]; ref=v_structures(C)
    for G0 in reachable_mpdags(C):
        NA=nonadjacent_pairs(G0)
        if not NA: continue
        ext=consistent_dag_extensions(G0,ref_vstructs=ref)
        Q=[]
        for x in range(p):
            for y in range(p):
                if x!=y:
                    O0,n0,am0=X.ostar_and_paths(G0,x,y)
                    if am0: Q.append((x,y,O0))
        if not Q: continue
        for (a0,b0) in NA:
          for (a,b) in ((a0,b0),(b0,a0)):
            H,info=bk_assert(G0,[(a,b)])
            if info["conflict"] or has_directed_cycle(H) or not pdag_extendable(H): continue
            for (x,y,O0) in Q:
                O1,n1,am1=X.ostar_and_paths(H,x,y)
                if not am1: continue
                fbH=adjust.forb(H,x,y)
                for D in ext:
                    fbD=adjust.forb(D,x,y)
                    if not (fbD<=fbH):
                        print(json.dumps(dict(C=fmt(C),G0=fmt(G0),H=fmt(H),D=fmt(D),stmt=f"{a}->{b}",
                            x=x,y=y,forbD=sorted(map(int,fbD)),forbH=sorted(map(int,fbH)),
                            cnD=sorted(map(int,adjust.causal_nodes(D,x,y))),
                            cnH=sorted(map(int,adjust.causal_nodes(H,x,y))),
                            O1=sorted(map(int,O1)),
                            valid=bool(adjust.is_valid_adjustment_set(D,x,y,set(O1))))))
                        n+=1
                        if n>=3: sys.exit()
