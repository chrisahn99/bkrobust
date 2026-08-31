"""
X2 core library -- locality of O* under Meek-consistent background knowledge.

Everything graph-theoretic is IMPORTED UNMODIFIED from the e1prime-se modules
(`graphs.py`, `adjust.py`, `scm.py`, `se.py`) copied into this directory.  The
only new code here is

  * `hop_dist_inf`        -- BFS distance with a real `inf` instead of the
                             1<<20 sentinel (PREREG 1.2, gate G7)
  * `ostar_and_paths`     -- O* computed from a SINGLE path enumeration (the
                             library recomputes the same DFS 3-4x); verified
                             exhaustively equal to `optimal_adjustment_set`
                             in test_machinery.py
  * `pcp_capped`          -- the same DFS with the PREREG 3.1 hard cap
  * `classify_trial`      -- the six-event taxonomy of PREREG 1.3 with the M15
                             nopath/unamenable splits
  * `draw_scm`            -- SCM draw replicating run_arm1.analyse_one's RNG
                             stream exactly (gate G3), plus ARM C's query rule
  * `prune_partition`     -- ARM D's B^r_rho (PREREG 2.5)

NOTHING under ~/latent-causal/e1prime-se/ or rho-breakdown-knowledge/ is edited.
"""
from itertools import combinations

import numpy as np

from adjust import (causal_nodes, cov_linear, forb, is_amenable,
                    is_valid_adjustment_set, optimal_adjustment_set,
                    ols_coefficient, parents_of_set, poss_de,
                    possibly_causal_paths, possibly_directed_neighbors,
                    total_effect_linear)
from graphs import (MeekFail, apply_background_knowledge, dag_agrees_with,
                    dag_to_cpdag, directed_edges, is_directed, is_undirected,
                    meek_closure, random_dag, skeleton, undirected_edges,
                    v_structures)
from scm import make_linear_iscm
from se import ols_with_se

INF = float("inf")
SENTINEL = 1 << 20
PATH_CAP = 200_000            # PREREG 3.1, registered
SCM_CAP = 20.0                # seconds, PREREG 3.1
TAU_FLOOR_DRAW = 1e-6         # M16: fixed in writing
TAU_FLOOR_SENS = 1e-3         # M16: reported as a sensitivity
N_REPORT = 20000              # the reported finite n; n = inf is decisive
Z95 = 1.959963984540054
EPS_BETA = 1e-12              # |beta' - beta0| > EPS  <=>  E_est(inf)


class PathCapExceeded(Exception):
    pass


# ------------------------------------------------------------------ distance
def hop_dist_inf(C, srcs):
    """BFS on skeleton(C) from `srcs`; unreachable -> float('inf') (NOT 1<<20)."""
    p = C.shape[0]
    S = skeleton(C)
    dist = np.full(p, np.inf, dtype=float)
    frontier = []
    for v in srcs:
        dist[v] = 0.0
        frontier.append(int(v))
    while frontier:
        nxt = []
        for v in frontier:
            for w in np.flatnonzero(S[v]):
                w = int(w)
                if dist[w] > dist[v] + 1:
                    dist[w] = dist[v] + 1
                    nxt.append(w)
        frontier = nxt
    return dist


def stmt_dist(dist, s):
    """(dmin, dmax) of a statement s=(u,v).  PREREG 1.2."""
    du, dv = float(dist[s[0]]), float(dist[s[1]])
    return (du if du <= dv else dv), (du if du >= dv else dv)


def dbucket(d):
    """Distance stratum label.  `disc` is NEVER merged into any '>= k' bucket (G7)."""
    if not np.isfinite(d):
        return "disc"
    d = int(d)
    return str(d) if d <= 3 else "ge4"


def far_finite(d):
    """PREREG 1.4 as amended by M8: Omega(r>=1) is on FINITE dmin only."""
    return np.isfinite(d) and d >= 1


# ------------------------------------------------------------------ O* in one pass
def pcp_capped(G, x, y, cap=PATH_CAP):
    """`possibly_causal_paths` with a hard cap.  Raises PathCapExceeded."""
    paths = []

    def dfs(path, visited):
        v = path[-1]
        if v == y:
            paths.append(list(path))
            if len(paths) > cap:
                raise PathCapExceeded()
            return
        for w in possibly_directed_neighbors(G, v):
            if w in visited:
                continue
            visited.add(w)
            path.append(w)
            dfs(path, visited)
            path.pop()
            visited.remove(w)

    dfs([x], {x})
    return paths


def ostar_and_paths(G, x, y, cap=PATH_CAP):
    """
    Returns (O*, n_paths, amenable).

    O* is `pa(cn) \\ forb` exactly as adjust.optimal_adjustment_set, but the
    path DFS is run ONCE instead of once per call inside is_amenable /
    causal_nodes / forb.  Equality with the library function is a machinery
    test, not an assumption (test_machinery.T5).

    `n_paths == 0` distinguishes M15's `nopath` from `unamenable`: the library
    collapses both into `amenable = False` (adjust.py:85-86).
    """
    paths = pcp_capped(G, x, y, cap)
    npaths = len(paths)
    if npaths == 0:
        return None, 0, False
    if not all(is_directed(G, pth[0], pth[1]) for pth in paths):
        return None, npaths, False
    cn = set()
    for pth in paths:
        cn.update(pth[1:])
    fb = (poss_de(G, cn) | {x}) if cn else {x}
    O = parents_of_set(G, cn) - fb
    return frozenset(O), npaths, True


# ------------------------------------------------------------------ the six events
EVENT_FIELDS = (
    "meek_fail", "e_ident", "e_ident_nopath", "e_ident_unamen",
    "e_gain", "e_gain_paths", "e_gain_orient",
    "e_set", "e_bias", "e_est_inf", "e_est_n", "cascade_pos",
)


def classify_trial(C, D, Sigma, x, y, Kp, base, cap=PATH_CAP):
    """
    Classify one perturbed knowledge set K' against the base.

    `base` is the dict returned by `make_base`.  Returns a flat dict with the
    six events of PREREG 1.3, the M15 splits, and the bookkeeping needed for
    the four (five, with M8) denominators.
    """
    out = {f: False for f in EVENT_FIELDS}
    out["consistent"] = False
    out["amenable0"] = base["amen0"]
    out["amenable1"] = False
    out["cascade"] = None
    out["est"] = None
    out["ostar"] = None
    try:
        G = apply_background_knowledge(C, Kp)
    except MeekFail:
        out["meek_fail"] = True
        return out
    out["consistent"] = True
    O, npaths, amen1 = ostar_and_paths(G, x, y, cap)
    out["amenable1"] = bool(amen1)
    out["cascade"] = int(len(directed_edges(G)) - base["ndir_C"] - len(Kp))
    out["cascade_pos"] = out["cascade"] > 0
    out["expels_true_dag"] = not dag_agrees_with(D, G)

    a0, a1 = base["amen0"], amen1
    if a0 and not a1:
        out["e_ident"] = True
        out["e_ident_nopath"] = (npaths == 0)
        out["e_ident_unamen"] = (npaths > 0)
    elif (not a0) and a1:
        out["e_gain"] = True
        out["e_gain_paths"] = (base["npaths0"] == 0)      # paths appeared
        out["e_gain_orient"] = (base["npaths0"] > 0)      # paths were oriented
    elif a0 and a1:
        out["ostar"] = O
        if O != base["O0"]:
            out["e_set"] = True
            out["e_bias"] = not bool(is_valid_adjustment_set(D, x, y, O))
        est, _, _ = ols_with_se(Sigma, x, y, O, np.inf)
        out["est"] = float(est)
        d = abs(est - base["est0"])
        out["e_est_inf"] = bool(d > EPS_BETA)
        out["e_est_n"] = bool(d > Z95 * base["se_n"])
    return out


def make_base(C, D, Sigma, x, y, K, cap=PATH_CAP):
    """
    The baseline G0 = M(K), O0 = O*(G0).  K = [] gives ARM A (G0 = C).
    Returns None if M(K) itself is MeekFail (cannot happen for an all-true K --
    machinery test T7a).
    """
    try:
        G0 = apply_background_knowledge(C, K)
    except MeekFail:
        return None
    O0, npaths0, amen0 = ostar_and_paths(G0, x, y, cap)
    base = dict(K=list(K), G0=G0, O0=O0, npaths0=npaths0, amen0=bool(amen0),
                ndir_C=len(directed_edges(C)), ndir_G0=len(directed_edges(G0)),
                est0=None, se_n=np.inf, k_reg=None, se_factor=None)
    if amen0:
        est0, _, sig2 = ols_with_se(Sigma, x, y, O0, np.inf)
        S0 = [x] + sorted(O0)
        Sinv = np.linalg.inv(Sigma[np.ix_(S0, S0)])
        se_factor = float(np.sqrt(max(sig2, 0.0) * Sinv[0, 0]))
        base["est0"] = float(est0)
        base["se_factor"] = se_factor
        base["k_reg"] = len(S0)
        base["se_n"] = se_factor / np.sqrt(N_REPORT - len(S0) - 1)
    return base


# ------------------------------------------------------------------ SCM draw
def draw_scm(seed, p, deg, need_u, query_rule="uniform", tau_floor=TAU_FLOOR_DRAW):
    """
    Replicates run_arm1.analyse_one's RNG stream EXACTLY for
    query_rule='uniform' (gate G3): random_dag -> pairs -> rng.integers ->
    make_linear_iscm -> rng.permutation(len(K_all)).

    query_rule='argmax_delta_C' is ARM C after M16: the candidate set is
    defined on C (`possibly_causal_paths(C,a,b)`), not on the true DAG, and the
    query maximises skeleton distance delta(X,Y).  The overlap with the D-based
    candidate set is returned as a diagnostic.
    """
    rng = np.random.default_rng(seed)
    D = random_dag(p, deg, rng)
    C = dag_to_cpdag(D)
    U = undirected_edges(C)
    if len(U) < need_u:
        return None
    pairs_D = [(a, b) for a in range(p) for b in range(p)
               if a != b and len(possibly_causal_paths(D, a, b)) > 0]
    if not pairs_D:
        return None
    if query_rule == "uniform":
        x, y = pairs_D[int(rng.integers(len(pairs_D)))]
        overlap = None
    elif query_rule == "argmax_delta_C":
        pairs_C = [(a, b) for a in range(p) for b in range(p)
                   if a != b and len(possibly_causal_paths(C, a, b)) > 0]
        if not pairs_C:
            return None
        sD, sC = set(pairs_D), set(pairs_C)
        overlap = dict(n_D=len(sD), n_C=len(sC), n_both=len(sD & sC))
        S = skeleton(C)
        best, bestd = [], -1.0
        for (a, b) in pairs_C:
            d = float(hop_dist_inf(C, [a])[b])
            if not np.isfinite(d):
                continue
            if d > bestd:
                best, bestd = [(a, b)], d
            elif d == bestd:
                best.append((a, b))
        if not best:
            return None
        x, y = best[int(rng.integers(len(best)))]
    else:
        raise ValueError(query_rule)
    scm = make_linear_iscm(D, rng)
    Sigma = cov_linear(scm["A"], scm["omega"])
    tau = total_effect_linear(scm["A"], x, y)
    if abs(tau) < tau_floor:
        return None
    K_all = [(i, j) if D[i, j] == 1 else (j, i) for (i, j) in U]
    order = rng.permutation(len(K_all))
    K_all = [K_all[i] for i in order]
    return dict(D=D, C=C, U=U, x=int(x), y=int(y), Sigma=Sigma, tau=float(tau),
                K_all=K_all, rng=rng, p=int(p), deg=float(deg), seed=int(seed),
                overlap=overlap, n_edges=int(D.sum()))


# ------------------------------------------------------------------ ARM D partition
def prune_partition(K, dmins, members, max_rho, r):
    """
    PREREG 2.5.  `members` is the E1'-shaped list (index <-> flip bijection over
    `for rho: for flip in combinations(range(k), rho)`), `dmins[i]` the dmin of
    statement K[i].

    Returns (I_full, I_prune, n_full, n_prune) where the intervals are
    [min,max] over consistent+amenable members with |flip| <= max_rho, plus the
    base K itself (flip = ()).  r = -1 prunes everything.
    """
    k = len(K)
    near = set(i for i in range(k) if np.isfinite(dmins[i]) and dmins[i] <= r) if r >= 0 else set()
    idx = 0
    vals_full, vals_prune = [], []
    n_full = n_prune = 1                     # K itself
    for rho in range(1, min(max_rho, k) + 1):
        for flip in combinations(range(k), rho):
            m = members[idx]
            idx += 1
            inside = set(flip) <= near
            n_full += 1
            if inside:
                n_prune += 1
            if not (m.get("consistent") and m.get("amenable")):
                continue
            e = m.get("est")
            if e is None:
                continue
            vals_full.append(e)
            if inside:
                vals_prune.append(e)
    return vals_full, vals_prune, n_full, n_prune


def interval(vals, est0):
    v = list(vals) + [est0]
    return (min(v), max(v))


# ------------------------------------------------------------------ statistics
def clopper_pearson(k, n, alpha=0.05):
    """Two-sided Clopper-Pearson interval.  M19."""
    from scipy.stats import beta as _beta
    if n == 0:
        return (0.0, 1.0)
    lo = 0.0 if k == 0 else float(_beta.ppf(alpha / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(_beta.ppf(1 - alpha / 2, k + 1, n - k))
    return (lo, hi)


def rule_of_three(n):
    """One-sided 95% upper bound on a rate given 0 events in n trials."""
    return 3.0 / n if n > 0 else float("inf")
