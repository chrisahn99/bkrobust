"""Do the three intrinsic coherence checks carry independent information?

Every cfirewall run reported rej_cycle == rej_notext exactly.  Equal counts are
not identity, so this tests the predicates POINTWISE:
   conflict?  cycle?  Dor-Tarsi non-extendable?
over the full p<=5 census (k=1) and a p=6 sample at k=1,2,3.
"""
import sys, importlib.util
from collections import Counter
import numpy as np
BASE = "/Users/josecosta/mugango/output/2026-08-19_x1-x2-locality-spurious"
sys.path.insert(0, BASE + "/x2/code"); sys.path.insert(0, BASE + "/x2/cfirewall/code")
from graphs import dag_to_cpdag, has_directed_cycle, random_dag, v_structures
from run_lemma import all_dags, reachable_mpdags
import cfirewall as CF
x1 = CF.x1_ops

def tally(G0, k, rng, n):
    c = Counter()
    NA = x1.nonadjacent_pairs(G0)
    if len(NA) < k: return c
    seqs = []
    if k == 1:
        seqs = [[(a,b)] for (a0,b0) in NA for (a,b) in ((a0,b0),(b0,a0))]
    else:
        for _ in range(n):
            idx = rng.choice(len(NA), size=k, replace=False)
            seqs.append([(NA[int(t)][0],NA[int(t)][1]) if rng.random()<.5
                         else (NA[int(t)][1],NA[int(t)][0]) for t in idx])
    for s in seqs:
        H, info = x1.bk_assert(G0, s)
        c[(bool(info["conflict"]), bool(has_directed_cycle(H)),
           not x1.pdag_extendable(H))] += 1
    return c

def census(p, k=1):
    c = Counter(); cp = {}
    for D in all_dags(p): C = dag_to_cpdag(D); cp.setdefault(C.tobytes(), C)
    for C in cp.values():
        for G0 in reachable_mpdags(C): c += tally(G0, k, None, 0)
    return c

def sample(p, k, n, seed):
    rng = np.random.default_rng(seed); c = Counter()
    for _ in range(n):
        C = dag_to_cpdag(random_dag(p, float(rng.choice([1.,1.5,2.,2.5,3.,3.5,4.,5.])), rng))
        G0 = CF.random_mpdag(C, v_structures(C), rng)
        c += tally(G0, k, rng, 25)
    return c

def show(lbl, c):
    tot = sum(c.values())
    print(f"{lbl}: n={tot}")
    for kk in sorted(c): print(f"   conflict={kk[0]} cycle={kk[1]} notext={kk[2]} : {c[kk]}")
    dis = sum(v for kk, v in c.items() if kk[1] != kk[2])
    print(f"   cycle XOR notext (disagreements) = {dis}")

if __name__ == "__main__":
    for p in (3,4,5): show(f"census p={p} k=1", census(p))
    for k in (1,2,3): show(f"sample p=6 k={k}", sample(6,k,4000,900+k))
    show("sample p=7 k=4", sample(7,4,2500,905))
