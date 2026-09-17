"""Certified-core sweep: reversal-only corruption model, canonicalised K."""
from __future__ import annotations
import itertools, json, sys, time
from collections import Counter
import numpy as np

from bkrobust.demo.graph import MPDAG, canon
from bkrobust.demo.meek import apply_orientations, enumerate_dag_extensions
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.demo.evaluate import (
    is_valid_adjustment_set_mpdag, optimal_adjustment_set_mpdag,
)

# ---------------------------------------------------------------- sampling
def random_dag(rng, n, p):
    nodes = [f"n{i}" for i in range(n)]
    order = list(nodes); rng.shuffle(order)
    edges = [(order[i], order[j]) for i in range(n) for j in range(i+1, n) if rng.random() < p]
    return MPDAG(nodes, directed=edges)

def true_orientations(dag, cpdag):
    """The true direction of every undirected edge of the CPDAG."""
    out = []
    for a, b in sorted(cpdag.undirected_edges):
        out.append((a, b) if dag.is_directed_edge(a, b) else (b, a))
    return out

# ---------------------------------------------------------------- machinery
class Problem:
    def __init__(self, cpdag, K, x, y, Z):
        self.cpdag, self.K, self.x, self.y, self.Z = cpdag, list(K), x, y, frozenset(Z)
        self._w = {}; self._v = {}
        self.n_checks = 0

    def world(self, J):
        """MPDAG after reversing the claims indexed by J (or None if inadmissible)."""
        J = frozenset(J)
        if J not in self._w:
            orients = [(h, t) if i in J else (t, h) for i, (t, h) in enumerate(self.K)]
            self._w[J] = apply_orientations(self.cpdag, orients)
        return self._w[J]

    def retract_world(self, J):
        J = frozenset(J)
        orients = [e for i, e in enumerate(self.K) if i not in J]
        return apply_orientations(self.cpdag, orients)

    def safe(self, J):
        J = frozenset(J)
        if J not in self._v:
            g = self.world(J)
            self.n_checks += 1
            self._v[J] = True if g is None else is_valid_adjustment_set_mpdag(g, self.x, self.y, self.Z)
        return self._v[J]

    def safe_retract(self, J):
        g = self.retract_world(J)
        return True if g is None else is_valid_adjustment_set_mpdag(g, self.x, self.y, self.Z)

def minimal_unsafe(prob, m, arm="reversal"):
    """Antichain of minimal unsafe subsets, with downward-closure pruning."""
    mins = []
    for size in range(1, m + 1):
        if size > 1 and not any(len(u) < size for u in mins) and size > m:
            break
        for J in itertools.combinations(range(m), size):
            Js = set(J)
            if any(set(u) <= Js for u in mins):
                continue
            ok = prob.safe(J) if arm == "reversal" else (prob.safe(J) and prob.safe_retract(J))
            if not ok:
                mins.append(tuple(J))
    return mins

def min_hitting_sets(mins, m):
    """All minimum-cardinality hitting sets of the antichain."""
    if not mins:
        return 0, [frozenset()]
    ground = sorted({i for u in mins for i in u})
    for size in range(1, len(ground) + 1):
        hs = [frozenset(S) for S in itertools.combinations(ground, size)
              if all(set(u) & set(S) for u in mins)]
        if hs:
            return size, hs
    return len(ground), [frozenset(ground)]

def canonicalise(cpdag, K, G0):
    """Meek-irredundant generating set of K: drop claims G0 re-derives without them."""
    cur = list(K)
    for k in list(K):
        trial = [e for e in cur if e != k]
        g = apply_orientations(cpdag, trial)
        if g is not None and g.directed_edges == G0.directed_edges:
            cur = trial
    return cur

# ---------------------------------------------------------------- one problem
def make_problem(rng, kmin, kmax):
    n = int(rng.integers(7, 10)); p = float(rng.uniform(0.22, 0.40))
    dag = random_dag(rng, n, p)
    try:
        cpdag = dag_to_cpdag(dag)
    except ValueError:
        return None
    A = true_orientations(dag, cpdag)
    if not A:
        return None
    size = int(rng.integers(1, len(A) + 1))
    K_raw = [A[i] for i in sorted(rng.choice(len(A), size=size, replace=False))]
    G0 = apply_orientations(cpdag, K_raw)
    if G0 is None:
        return None
    nodes = list(dag.nodes)
    pairs = [(a, b) for a in nodes for b in nodes if a != b and b in dag.descendants(a)]
    if not pairs:
        return None
    x, y = pairs[int(rng.integers(len(pairs)))]
    Z = optimal_adjustment_set_mpdag(G0, x, y)
    if Z is None or not is_valid_adjustment_set_mpdag(G0, x, y, Z):
        return None
    Kc = canonicalise(cpdag, K_raw, G0)
    if not (kmin <= len(Kc) <= kmax):
        return None
    return dict(dag=dag, cpdag=cpdag, K_raw=K_raw, Kc=Kc, G0=G0, x=x, y=y, Z=Z)

def analyse(P, do_retraction=True, do_padding=True):
    cpdag, Kc, x, y, Z, G0 = P["cpdag"], P["Kc"], P["x"], P["y"], P["Z"], P["G0"]
    m = len(Kc)
    prob = Problem(cpdag, Kc, x, y, Z)
    mins = minimal_unsafe(prob, m)
    c, hs = min_hitting_sets(mins, m)
    kappa = 1.0 - c / m
    phi1 = sum(1 for i in range(m) if not prob.safe((i,)))
    core = set(hs[0]); free = set(range(m)) - core
    # hedge support
    W = [J for J in powerset(m) if prob.world(J) is not None]
    M = [J for J in W if set(J) & core] + [()]
    M = {frozenset(j) for j in M}
    # theorem check
    thm = None
    if c == 0:
        thm = bool(is_valid_adjustment_set_mpdag(cpdag, x, y, Z))
    # conditional-vs-joint gap: core corruption safe, but adding a free subset breaks it
    gap = False
    for Jc in powerset_of(sorted(core)):
        if not Jc:
            continue
        if not prob.safe(Jc):
            continue
        if prob.world(Jc) is None:
            continue
        for Jf in powerset_of(sorted(free)):
            if not Jf:
                continue
            J = tuple(sorted(set(Jc) | set(Jf)))
            if prob.world(J) is not None and not prob.safe(J):
                gap = True; break
        if gap:
            break
    out = dict(m=m, c=c, kappa=kappa, phi1=phi1, n_mins=len(mins),
               max_min_size=max((len(u) for u in mins), default=0),
               n_checks=prob.n_checks, n_hs=len(hs), unique_core=len(hs) == 1,
               W=len(W), M=len(M), pow_m=2**m, pow_c=2**c, thm=thm, gap=gap,
               c_eq_phi1=(c == phi1))
    if do_retraction:
        mins_r = minimal_unsafe(Problem(cpdag, Kc, x, y, Z), m, arm="both")
        c_r, _ = min_hitting_sets(mins_r, m)
        out["c_retract"] = c_r
        out["retract_changed"] = (c_r != c)
    if do_padding:
        pads = [(t, h) for (t, h) in sorted(G0.directed_edges)
                if canon(t, h) in cpdag.undirected_edges and (t, h) not in Kc]
        out["n_pads"] = len(pads)
        if pads:
            Kp = list(Kc) + pads
            if len(Kp) <= 11:
                gp = apply_orientations(cpdag, Kp)
                same = gp is not None and gp.directed_edges == G0.directed_edges
                pp = Problem(cpdag, Kp, x, y, Z)
                mp = len(Kp)
                minsp = minimal_unsafe(pp, mp)
                cp, _ = min_hitting_sets(minsp, mp)
                out["pad_same_G0"] = same
                out["pad_c"] = cp
                out["pad_kappa"] = 1.0 - cp / mp
                out["pad_m"] = mp
    return out

def powerset(m):
    for size in range(m + 1):
        for J in itertools.combinations(range(m), size):
            yield J

def powerset_of(elems):
    for size in range(len(elems) + 1):
        for J in itertools.combinations(elems, size):
            yield J

# ---------------------------------------------------------------- sweep
def sweep(seed, n_problems, kmin, kmax, tag):
    rng = np.random.default_rng(seed)
    rows = []
    t0 = time.time()
    tries = 0
    while len(rows) < n_problems and tries < 200000 and time.time() - t0 < 900:
        tries += 1
        P = make_problem(rng, kmin, kmax)
        if P is None:
            continue
        rows.append(analyse(P))
    return dict(tag=tag, seed=seed, secs=round(time.time() - t0, 1), tries=tries, rows=rows)

if __name__ == "__main__":
    tag = sys.argv[1]; seed = int(sys.argv[2]); n = int(sys.argv[3])
    kmin, kmax = int(sys.argv[4]), int(sys.argv[5])
    res = sweep(seed, n, kmin, kmax, tag)
    with open(f"sweep_{tag}.json", "w") as f:
        json.dump(res, f)
    print(tag, "n=", len(res["rows"]), "secs=", res["secs"], "tries=", res["tries"])
