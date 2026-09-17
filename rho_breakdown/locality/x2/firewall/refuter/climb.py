"""Hill-climbing adversarial search.

Objective = where the theorem is LOAD-BEARING, not where it is easy:
  crit  : (D,w) with w in O*_G0 \\ O*_H and O*_G0\\{w} INVALID in D
          -> the stamp removed a node that was doing essential blocking work,
             so something else must have replaced it, or the theorem breaks.
  tight : (D,w) with w in O*_H \\ O*_G0 and O*_H\\{w} INVALID in D
          -> the replacement node is essential: the self-repair is exercised.
Any actual violation is reported immediately.
"""
import sys, time, json, random
sys.path.insert(0, "/private/tmp/claude-501/-Users-josecosta-mugango/bceac0a4-93a7-4757-b3b3-a3d4e6c42d55/scratchpad/cfw_ref_A")
from cfwlib import *
import cfwlib as _CF


def evaluate(D0, ufull=14, want_hits=None):
    p = D0.shape[0]
    C = dag_to_cpdag(D0); ref = v_structures(C)
    G0 = C
    U0 = undirected_edges(G0); S = skeleton(G0)
    NA = [(a, b) for a in range(p) for b in range(a + 1, p) if S[a, b] == 0]
    if not NA or len(U0) > ufull or len(U0) == 0:
        return -1.0, 0, None
    full = consistent_dag_extensions(G0, ref_vstructs=ref)
    if not full: return -1.0, 0, None
    crit = tight = 0; nmoved = 0
    for x in range(p):
        for y in range(p):
            if x == y: continue
            O0, np0, amen0 = X.ostar_and_paths(G0, x, y)
            if not amen0: continue
            for (a0, b0) in NA:
                for (a, b) in ((a0, b0), (b0, a0)):
                    H, info = x1_ops.bk_assert(G0, [(a, b)])
                    if info["conflict"] or has_directed_cycle(H) or not x1_ops.pdag_extendable(H):
                        continue
                    O1, np1, amen1 = X.ostar_and_paths(H, x, y)
                    if not amen1 or O1 == O0: continue
                    nmoved += 1
                    for Dc in full:
                        if not valid_in(Dc, x, y, O1):
                            if want_hits is not None:
                                want_hits.append(dict(G0=G0.tolist(), H=H.tolist(),
                                                      x=int(x), y=int(y), edge=[int(a), int(b)],
                                                      O0=sorted(map(int, O0)), O1=sorted(map(int, O1)),
                                                      Dbad=Dc.tolist()))
                            return 1e9, nmoved, "VIOLATION"
                        for w in (set(O0) - set(O1)):
                            if not valid_in(Dc, x, y, frozenset(set(O0) - {w})):
                                crit += 1
                        for w in (set(O1) - set(O0)):
                            if not valid_in(Dc, x, y, frozenset(set(O1) - {w})):
                                tight += 1
    return float(2 * crit + tight), nmoved, None


def perturb(D, rng):
    p = D.shape[0]
    E = D.copy()
    for _ in range(rng.integers(1, 3)):
        i, j = rng.integers(0, p, 2)
        if i == j: continue
        if E[i, j] or E[j, i]:
            E[i, j] = E[j, i] = 0
        else:
            E[i, j] = 1
            A = E.astype(bool); R = A.copy()
            for _ in range(p): R = R | (R @ A)
            if np.any(np.diag(R)): E[i, j] = 0; E[j, i] = 1
            A = E.astype(bool); R = A.copy()
            for _ in range(p): R = R | (R @ A)
            if np.any(np.diag(R)): E[j, i] = 0
    return E


def climb(p, seed, steps, tlimit):
    rng = np.random.default_rng(seed)
    best = random_dag(p, 2.5, rng)
    bs, _, v = evaluate(best)
    hits = []
    t0 = time.time(); nev = 1
    for s in range(steps):
        if time.time() - t0 > tlimit: break
        cand = perturb(best, rng)
        sc, nm, v = evaluate(cand, want_hits=hits); nev += 1
        if v == "VIOLATION":
            return "VIOLATION", cand, hits, nev
        if sc >= bs:
            best, bs = cand, sc
        if s % 400 == 0:
            best = random_dag(p, float(rng.choice([2.0, 2.5, 3.0, 3.5])), rng)
            bs, _, _ = evaluate(best)
        if len(_CF._vcache) > 300000: _CF._vcache.clear()
    return None, best, hits, nev


if __name__ == "__main__":
    p = int(sys.argv[1]); seed = int(sys.argv[2]); steps = int(sys.argv[3]); tl = float(sys.argv[4])
    r, D, hits, nev = climb(p, seed, steps, tl)
    print("p", p, "seed", seed, "result", r, "n_eval", nev)
    if r: json.dump(hits, open(f"climb_hit_{p}_{seed}.json", "w"))
