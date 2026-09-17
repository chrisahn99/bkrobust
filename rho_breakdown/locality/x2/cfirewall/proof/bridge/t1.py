"""T1': bk_assert(G0,[(a,b)]) == meek_closure(G0 + a->b) for every non-adjacent pair,
   p=4 census; and conflict never fires."""
from common import *
n=bad=0; conf=0
for p in (3,4):
    cp={}
    for D in all_dags(p):
        C=dag_to_cpdag(D); cp.setdefault(C.tobytes(),C)
    for kb in cp:
        C=cp[kb]
        for G0 in reachable_mpdags(C):
            for (a0,b0) in nonadjacent_pairs(G0):
                for (a,b) in ((a0,b0),(b0,a0)):
                    H,info=bk_assert(G0,[(a,b)])
                    P=G0.copy(); P[a,b]=1; P[b,a]=0
                    Hm=meek_closure(P)
                    n+=1
                    if not np.array_equal(H,Hm): bad+=1
                    conf+=int(info["conflict"])
                    assert not np.array_equal(skeleton(H),skeleton(G0))
print(f"T1': {n} statements, restamp-differs={bad}, conflicts={conf}, skeleton always grew")
