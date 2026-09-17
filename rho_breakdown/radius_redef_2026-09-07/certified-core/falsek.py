"""The realistic setting: the elicited K contains FALSE claims consistent with the CPDAG."""
from __future__ import annotations
import os, sys, json, time, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import (random_dag, true_orientations, canonicalise, Problem,
                  minimal_unsafe, min_hitting_sets, powerset)
from bkrobust.demo.graph import MPDAG, canon
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.evaluate import (random_sem, adjusted_estimand, optimal_adjustment_set_dag,
                                    optimal_adjustment_set_mpdag, is_valid_adjustment_set_mpdag,
                                    is_valid_adjustment_set_dag)

def make(rng, kmin, kmax, p_flip):
    n = int(rng.integers(7, 10)); p = float(rng.uniform(0.22, 0.40))
    dag = random_dag(rng, n, p)
    try: cpdag = dag_to_cpdag(dag)
    except ValueError: return None
    A = true_orientations(dag, cpdag)
    if not A: return None
    size = int(rng.integers(1, len(A)+1))
    sub = [A[i] for i in sorted(rng.choice(len(A), size=size, replace=False))]
    flip = [i for i in range(len(sub)) if rng.random() < p_flip]
    K = [ (sub[i][1], sub[i][0]) if i in flip else sub[i] for i in range(len(sub)) ]
    G0 = apply_orientations(cpdag, K)
    if G0 is None: return None
    nodes = list(dag.nodes)
    pairs = [(a,b) for a in nodes for b in nodes if a != b and b in dag.descendants(a)]
    if not pairs: return None
    x, y = pairs[int(rng.integers(len(pairs)))]
    Z = optimal_adjustment_set_mpdag(G0, x, y)
    if Z is None or not is_valid_adjustment_set_mpdag(G0, x, y, Z): return None
    Kc = canonicalise(cpdag, K, G0)
    if not (kmin <= len(Kc) <= kmax): return None
    return dict(dag=dag, cpdag=cpdag, Kc=Kc, G0=G0, x=x, y=y, Z=Z, n_flip=len(flip))

def run(seed, n_problems, kmin, kmax, p_flip, tag):
    rng = np.random.default_rng(seed); rows=[]; t0=time.time(); tries=0
    while len(rows) < n_problems and tries < 600000 and time.time()-t0 < 1500:
        tries += 1
        P = make(rng, kmin, kmax, p_flip)
        if P is None: continue
        cpdag, Kc, x, y, Z, dag, G0 = P["cpdag"], P["Kc"], P["x"], P["y"], P["Z"], P["dag"], P["G0"]
        m = len(Kc); prob = Problem(cpdag, Kc, x, y, Z)
        mins = minimal_unsafe(prob, m); c, hs = min_hitting_sets(mins, m); core = set(hs[0])
        sem = random_sem(dag, rng)
        th = sem.true_total_effect(x, y)
        pt = adjusted_estimand(sem, x, y, Z)
        cache = {}
        def vals_of(g):
            key = (g.directed_edges, g.undirected_edges)
            if key not in cache:
                cache[key] = [adjusted_estimand(sem, x, y, optimal_adjustment_set_dag(d, x, y))
                              for d in enumerate_dag_extensions(g)]
            return cache[key]
        M=[pt]; W=[pt]; B1=[pt]; F=[]
        for J in powerset(m):
            g = prob.world(J)
            if g is None: continue
            v = vals_of(g)
            W.extend(v)
            if len(J) <= 1: B1.extend(v)
            if (not J) or (set(J) & core): M.extend(v)
            elif J: F.extend(v)
        z_valid_truth = is_valid_adjustment_set_dag(dag, x, y, Z)
        tol = 1e-8 + 1e-6*abs(th)
        def cov(L): return bool(min(L)-tol <= th <= max(L)+tol)
        def wid(L): return max(L)-min(L)
        rows.append(dict(m=m, c=c, nflip=P["n_flip"], zvalid=bool(z_valid_truth),
                         wM=wid(M), wW=wid(W), wB1=wid(B1),
                         covM=cov(M), covW=cov(W), covB1=cov(B1),
                         covPt=bool(abs(pt-th)<=tol),
                         bias=abs(pt-th), scale=abs(th),
                         Fsame=bool(not F or max(abs(v-pt) for v in F) < 1e-9),
                         nM=sum(1 for J in powerset(m) if prob.world(J) is not None and ((not J) or set(J)&core)),
                         nW=sum(1 for J in powerset(m) if prob.world(J) is not None)))
    json.dump(dict(tag=tag,rows=rows), open(f"falsek_{tag}.json","w"))
    R=rows; n=len(R); g=lambda k: np.array([r[k] for r in R])
    print(f"--- falsek {tag}: n={n} |K*|={g('m').mean():.2f} c={g('c').mean():.3f} P(c=0)={np.mean(g('c')==0):.4f} mean flips/problem={g('nflip').mean():.2f}")
    print(f"  P(Z valid in the TRUE dag) = {g('zvalid').mean():.4f}   (so the point estimate is wrong on {1-g('zvalid').mean():.4f})")
    print(f"  structural COVERAGE: point={g('covPt').mean():.4f}  radius-1 hedge={g('covB1').mean():.4f}  certified-core hedge={g('covM').mean():.4f}  blanket/IDA={g('covW').mean():.4f}")
    print(f"  structural WIDTH   : point=0.0000  radius-1={g('wB1').mean():.4f}  certified-core={g('wM').mean():.4f}  blanket/IDA={g('wW').mean():.4f}")
    print(f"  certified worlds give the SAME estimand as the point estimate on {g('Fsame').mean():.4f}")
    print(f"  wM == wW on {np.mean(np.abs(g('wM')-g('wW'))<1e-9):.4f};  support |W|={g('nW').mean():.2f} |M|={g('nM').mean():.2f} ratio={np.mean(g('nW')/g('nM')):.3f}")
    z=g('c')==0
    if z.sum(): print(f"  c=0 stratum (n={z.sum()}): P(Z valid in truth)={g('zvalid')[z].mean():.4f} mean wW there={g('wW')[z].mean():.4f}")

if __name__=="__main__":
    run(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), float(sys.argv[6]), sys.argv[1])
