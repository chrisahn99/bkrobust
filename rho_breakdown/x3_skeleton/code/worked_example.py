"""
A CONCRETE USE OF THE FRAMEWORK, end to end, on one small problem.

The agenda item for the Thursday meeting: "preparer un cas concret d'utilisation du framework
pour mieux comprendre les mecanismes".

What it walks through, in the order an analyst meets them:
    the CPDAG the discovery method returned  ->  what it leaves undetermined
    the expert's knowledge K                 ->  the MPDAG G0 = Meek(Chat, K)
    O*(G0) and the estimate                  ->  what is being claimed
    r_val, with a WITNESS                    ->  what would have to be wrong
    r_eps and the ignorance interval         ->  how much it would move the number

Two radii, never one (framing section 4.2), and the witness matters more than the number
(section 5).  Atomic moves are the framing's own (section 4.1): ORIENT an undirected edge, or
UN-ORIENT a knowledge-oriented one.  A FLIP is a retraction plus an assertion and costs two.

SYNTHETIC.  The graph is drawn from the same generator as every other run in this line.  The
variable names are a reading aid so the witness is contestable in words; no empirical claim
about smoking, statins or anything else is made or implied.

Usage: python worked_example.py [--search] > ../results/worked_example.txt
"""
import sys
from itertools import combinations

import numpy as np

from adjust import (cov_linear, is_valid_adjustment_set,
                    optimal_adjustment_set, ols_coefficient,
                    possibly_causal_paths, total_effect_linear)
from graphs import (MeekFail, apply_background_knowledge,
                    consistent_dag_extensions, dag_to_cpdag, directed_edges,
                    random_dag, topological_order, undirected_edges)
from scm import make_linear_iscm

# Names are attached BY TOPOLOGICAL POSITION in the true DAG, so every edge runs
# from an earlier name to a later one and the skin is at least internally coherent.
# The structure is whatever the generator produced; the names carry no claim.
ORDER = ["Age", "Smoking", "Exercise", "BMI", "BloodPressure", "Statin", "CardiacEvent"]
P = 7
NAMES = list(ORDER)          # rebound per case in build()


# --------------------------------------------------------------------- the poset
def neighbours(Chat, K):
    """One atomic move from Meek(Chat, K), returned as the resulting K'.

    DOWN: orient one still-undirected edge of G0, either way  (assert something new)
    UP:   drop one statement of K                             (retract something)
    Meek closure is applied by apply_background_knowledge; MeekFail is the bottom
    element of section 4.1 and the move is discarded.
    """
    out = []
    try:
        G0 = apply_background_knowledge(Chat, K)
    except MeekFail:
        return out
    for (a, b) in undirected_edges(G0):
        for (u, v) in ((a, b), (b, a)):
            Kp = list(K) + [(u, v)]
            try:
                apply_background_knowledge(Chat, Kp)
            except MeekFail:
                continue
            out.append((tuple(sorted(Kp)), ("orient", u, v)))
    for i in range(len(K)):
        Kp = list(K[:i]) + list(K[i + 1:])
        try:
            apply_background_knowledge(Chat, Kp)
        except MeekFail:
            continue
        out.append((tuple(sorted(Kp)), ("un-orient", K[i][0], K[i][1])))
    return out


def valid_in_mpdag(G, x, y, Z, limit=4000):
    """Z is valid in an MPDAG iff it is valid in EVERY DAG the MPDAG represents."""
    dags = consistent_dag_extensions(G, limit=limit)
    if dags is None:
        return None
    for D in dags:
        if not is_valid_adjustment_set(D, x, y, Z):
            return False
    return True


def grow_ball(Chat, K0, x, y, Z, Sigma, tau0, max_r=6, eps=None):
    """Shell-by-shell BFS.  Returns r_val with its witness, and the ignorance interval
    per shell (which gives r_eps)."""
    seen = {tuple(sorted(K0))}
    frontier = [(tuple(sorted(K0)), [])]
    shells = []
    r_val, witness = None, None
    for r in range(1, max_r + 1):
        nxt, ests = [], []
        for (K, path) in frontier:
            for (Kp, move) in neighbours(Chat, list(K)):
                if Kp in seen:
                    continue
                seen.add(Kp)
                G = apply_background_knowledge(Chat, list(Kp))
                nxt.append((Kp, path + [move]))
                # certificate side: Z is FIXED, validity re-evaluated
                if r_val is None:
                    ok = valid_in_mpdag(G, x, y, Z)
                    if ok is False:
                        r_val, witness = r, path + [move]
                # ignorance side: O* recomputed, estimate moves
                O = optimal_adjustment_set(G, x, y)
                if O is not None:
                    try:
                        ests.append(float(ols_coefficient(Sigma, x, y, O)))
                    except Exception:
                        pass
        shells.append(dict(r=r, n=len(nxt), lo=(min(ests) if ests else None),
                           hi=(max(ests) if ests else None)))
        frontier = nxt
        if not frontier:
            break
    return r_val, witness, shells


# --------------------------------------------------------------------- selection
def build(seed):
    rng = np.random.default_rng(seed)
    D = random_dag(P, 2.0, rng)
    C = dag_to_cpdag(D)
    U = undirected_edges(C)
    if not (3 <= len(U) <= 6):
        return None
    pairs = [(a, b) for a in range(P) for b in range(P)
             if a != b and len(possibly_causal_paths(D, a, b)) > 0]
    if not pairs:
        return None
    x, y = pairs[int(rng.integers(len(pairs)))]
    scm = make_linear_iscm(D, rng)
    Sigma = cov_linear(scm["A"], scm["omega"])
    tau = total_effect_linear(scm["A"], x, y)
    if abs(tau) < 0.2:
        return None
    K_all = [(i, j) if D[i, j] == 1 else (j, i) for (i, j) in U]
    K = [tuple(int(v) for v in e) for e in K_all[:2]]
    topo = topological_order(D)
    names = [None] * P
    for pos, node in enumerate(topo):
        names[node] = ORDER[pos]
    try:
        G0 = apply_background_knowledge(C, K)
    except MeekFail:
        return None
    O0 = optimal_adjustment_set(G0, x, y)
    if O0 is None or not is_valid_adjustment_set(D, x, y, O0):
        return None
    est0 = float(ols_coefficient(Sigma, x, y, O0))
    return dict(seed=seed, D=D, C=C, K=K, G0=G0, x=x, y=y, O0=O0,
                Sigma=Sigma, tau=float(tau), est0=est0, U=U, names=names)


def fmt(v):
    return NAMES[v]


def edges_str(G, kind):
    if kind == "dir":
        return ", ".join(f"{fmt(i)} -> {fmt(j)}" for (i, j) in directed_edges(G))
    return ", ".join(f"{fmt(i)} - {fmt(j)}" for (i, j) in undirected_edges(G))


STAT = __import__("collections").Counter()


def main():
    if "--search" in sys.argv:
        for s in range(6000):
            c = build(s)
            if c is None:
                continue
            rv, w, sh = grow_ball(c["C"], c["K"], c["x"], c["y"], c["O0"],
                                  c["Sigma"], c["est0"], max_r=3)
            nm = c["names"]
            STAT["built"] += 1
            if rv in (2, 3) and w: STAT["rval"] += 1
            if len(c["O0"]) >= 2: STAT["O2"] += 1
            if len(undirected_edges(c["G0"])) >= 1: STAT["undet"] += 1
            if w and any(k == "orient" for (k, _, _) in w): STAT["orient"] += 1
            _w = max(((x["hi"] - x["lo"]) if x["lo"] is not None else 0.0) for x in sh) if sh else 0.0
            if _w > 1e-9: STAT["width"] += 1
            if rv not in (2, 3) or not w:
                continue
            if len(c["O0"]) < 2:
                continue
            wid = [(x["r"], (x["hi"] - x["lo"]) if x["lo"] is not None else 0.0)
                   for x in sh]
            wmax = max((v for _, v in wid), default=0.0)
            if wmax <= 1e-9:
                continue
            STAT['PASS'] += 1
            print(f"seed={s} r_val={rv} |O*|={len(c['O0'])} |U|={len(c['U'])} "
                  f"undet(G0)={len(undirected_edges(c['G0']))} "
                  f"x={nm[c['x']]} y={nm[c['y']]} "
                  f"O*={{{', '.join(nm[v] for v in sorted(c['O0']))}}} "
                  f"widths={[f'{r}:{v:.3f}' for r, v in wid]}")
        print("\nATTRITION:", dict(STAT))
        return

    SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    c = build(SEED)
    assert c is not None, "seed does not build a case"
    global NAMES
    NAMES = c["names"]
    x, y, C, K, G0, O0 = c["x"], c["y"], c["C"], c["K"], c["G0"], c["O0"]

    print("=" * 78)
    print("A CONCRETE USE OF THE FRAMEWORK".center(78))
    print("=" * 78)
    print(f"""
SYNTHETIC.  The graph comes from the same generator as every other run in this
line (random_dag, p=7, expected degree 2, iSCM linear-Gaussian, population Sigma
so Monte-Carlo error is exactly zero).  The variable names are a reading aid, so
that the witness at the end is a sentence a domain expert can argue with.  No
empirical claim about any of these variables is made.  Seed {SEED}.
""")
    print("-" * 78)
    print("STEP 1.  What the discovery method returned")
    print("-" * 78)
    print(f"  CPDAG, oriented   : {edges_str(C, 'dir')}")
    print(f"  CPDAG, UNDETERMINED: {edges_str(C, 'und')}")
    print(f"""
  The undirected edges are what the observational distribution cannot resolve at
  ANY sample size.  They are the whole arena: background knowledge can only ever
  orient these.  Query: effect of {fmt(x)} on {fmt(y)}.""")

    print("\n" + "-" * 78)
    print("STEP 2.  What the expert asserted, and what it bought")
    print("-" * 78)
    for (a, b) in K:
        print(f"  K: {fmt(a)} -> {fmt(b)}")
    print(f"\n  MPDAG G0 = Meek(Chat, K)")
    print(f"    oriented    : {edges_str(G0, 'dir')}")
    print(f"    still undet.: {edges_str(G0, 'und') or '(none)'}")
    n_prop = len(directed_edges(G0)) - len(directed_edges(C)) - len(K)
    print(f"""
  {len(K)} statements were asserted; Meek's rules PROPAGATED {n_prop} further
  orientations for free.  That propagation is why counting wrong assertions is
  the wrong unit: one statement in a dense region can collapse a whole component,
  the same statement elsewhere propagates nothing.""")

    print("\n" + "-" * 78)
    print("STEP 3.  The claim being made")
    print("-" * 78)
    print(f"  O*(G0) = {{{', '.join(fmt(v) for v in sorted(O0))}}}")
    print(f"  estimate, adjusting for O*      : {c['est0']:+.4f}")
    print(f"  true total effect (hidden)      : {c['tau']:+.4f}")
    print(f"""
  The analyst never sees the second line.  Everything from here answers the only
  question available at analysis time: how wrong would K have to be?""")

    rv, wit, shells = grow_ball(C, K, x, y, O0, c["Sigma"], c["est0"], max_r=5)

    print("\n" + "-" * 78)
    print("STEP 4.  r_val, and the witness")
    print("-" * 78)
    print(f"""  The ball is centred on the graph the ANALYST HOLDS, not on the truth --
  that is what makes it computable in deployment.  Z is fixed once at O*(G0);
  validity is then re-evaluated in every graph of every shell, where "valid in an
  MPDAG" means valid in every DAG that MPDAG represents.
""")
    for s in shells:
        print(f"    shell r={s['r']}: {s['n']:4} realisable revisions")
    if rv is None:
        print(f"\n  r_val > {len(shells)}  (no revision within the budget breaks O*(G0))")
    else:
        print(f"\n  r_val = {rv}")
        print("\n  WITNESS -- the specific sequence that breaks it:")
        for i, (kind, a, b) in enumerate(wit, 1):
            print(f"    {i}. {kind:>9}  {fmt(a)} -> {fmt(b)}")
        print(f"""
  In words, and this is the output that matters:

     "The estimate holds unless {' and '.join(f'{fmt(a)} -> {fmt(b)} is ' + ('retracted' if k=='un-orient' else 'added') for (k,a,b) in wit)}."

  A domain expert answers that directly, with no metric literacy at all.  The
  radius says how alarmed to be; the witness says what to go and check.""")
    print(f"""
  A FLIP costs TWO (retract, then assert), so r_val = {rv if rv else '>'+str(len(shells))} reads to a
  practitioner as: {'one flip, or two omissions' if rv==2 else 'see the shell table'}.""")

    print("\n" + "-" * 78)
    print("STEP 5.  r_eps -- what it would do to the NUMBER, against the standard error")
    print("-" * 78)
    print(f"""  Here O* is RE-derived in each member, so the estimate moves.  The spread is
  the ignorance interval: what the data plus this much knowledge-doubt allow.
""")
    print(f"    {'radius':>7}  {'ignorance interval':>28}   width")
    print(f"    {0:>7}  [{c['est0']:+.4f}, {c['est0']:+.4f}]   0.0000")
    for s in shells:
        if s["lo"] is None:
            continue
        print(f"    {s['r']:>7}  [{s['lo']:+.4f}, {s['hi']:+.4f}]   {s['hi']-s['lo']:.4f}")
    first_move = next((s["r"] for s in shells
                       if s["lo"] is not None and (s["hi"] - s["lo"]) > 1e-9), None)
    print(f"""
  Read this against the standard error, which is the comparison a practitioner
  actually makes.  If the interval at r=1 sits inside the confidence interval,
  structural doubt is not the binding constraint.  If it is wider, no amount of
  data helps: SE shrinks like n^-1/2 and this width does not shrink at all.""")
    if first_move is not None and rv is not None and first_move > rv:
        w = [s for s in shells if s["r"] == first_move][0]
        print(f"""
  AND THIS CASE IS THE REASON THE FRAMEWORK REPORTS TWO RADII.

    r_val = {rv}          the adjustment set stops being provably valid here
    r_eps > {first_move - 1}          the estimate does not move at all until r = {first_move}

  Between them the analysis is in the middle row of the framing's own table:
  FRAGILE IDENTIFICATION, STABLE ESTIMATE.  A validity-only analysis raises the
  alarm at r = {rv}; the number it is alarmed about has not moved by so much as
  1e-9.  At r = {first_move} it finally moves, by {w['hi'] - w['lo']:.4f} against an effect of
  {c['est0']:+.4f}, i.e. {abs((w['hi'] - w['lo']) / c['est0']) * 100:.0f}% of the estimate.

  Collapsing the two into one number would have reported this analysis as
  fragile at r = {rv}.  It is not: it is fragile in IDENTIFICATION and stable in
  ESTIMATION, and those are different things to tell a practitioner.""")
    print(f"""
------------------------------------------------------------------------------
A STRUCTURAL REMARK THAT FELL OUT OF BUILDING THIS
------------------------------------------------------------------------------
  Across every case searched, the last move of the minimal breaking path was
  ALWAYS a retraction, never an assertion.  That is not luck, it is forced:

    if G' is a REFINEMENT of G, then [G'] is a SUBSET of [G], so a Z valid in
    every DAG of [G] is valid in every DAG of [G'].

  Orienting one more edge can therefore never invalidate a fixed valid Z.  For
  r_val -- where Z is fixed at O*(G0) -- only the UPWARD moves bind.  The DOWNWARD
  move earns its place in the framework through r_eps, where O* is re-derived and
  the estimate does move.  Worth stating in the paper: it removes a whole
  direction from the r_val search.
""")
    print("=" * 78)


if __name__ == "__main__":
    main()
