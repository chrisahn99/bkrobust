"""POSITIVE CONTROL: the same detector, with ONE guard removed at a time.
If the detector is capable of seeing O* move, it must fire here.
"""
import sys, json, random
sys.path.insert(0, "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/code")
import numpy as np
import x2lib as X
from graphs import (meek_closure, v_structures, has_directed_cycle,
                    undirected_edges, dag_to_cpdag, random_dag)


def run(mode, p, n, seed):
    rng = np.random.default_rng(seed); prng = random.Random(seed)
    n_set = 0; d = 0; ex = None
    for _ in range(n):
        D = random_dag(p, prng.choice([1.5, 2.0, 2.5, 3.0]), rng)
        C = dag_to_cpdag(D); ref = v_structures(C)
        for (u, v) in undirected_edges(C):
            for (a, b) in ((u, v), (v, u)):
                H = C.copy(); H[b, a] = 0
                if mode != "nomeek":
                    H = meek_closure(H)
                if has_directed_cycle(H):
                    continue
                if mode != "novstruct" and v_structures(H) != ref:
                    continue
                for x in range(p):
                    for y in range(p):
                        if x == y: continue
                        O0, _, a0 = X.ostar_and_paths(C, x, y)
                        O1, _, a1 = X.ostar_and_paths(H, x, y)
                        if not (a0 and a1): continue
                        d += 1
                        if O0 != O1:
                            n_set += 1
                            if ex is None:
                                ex = dict(x=x, y=y, e=[a, b], G0=C.tolist(),
                                          H=H.tolist(), O0=sorted(O0), O1=sorted(O1))
    return dict(mode=mode, p=p, n_dags=n, D_am=d, N_set=n_set, example=ex)


if __name__ == "__main__":
    for mode in ("novstruct", "nomeek"):
        print(json.dumps(run(mode, 6, 400, 7), default=float)[:900])
