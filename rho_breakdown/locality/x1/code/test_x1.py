"""
X1 MACHINERY TESTS — run BEFORE the experiment, and reported with it.

"A run whose machinery is untested is not a result."  Every new operator is
validated against its DEFINITION by brute force on small graphs, never against a
reimplementation of itself.

Usage:  python test_x1.py [seed]
Exit 0 iff every test passes.  Prints a table; T5/T6 also print MEASUREMENTS
(what the resulting object IS) which are part of the deliverable, not just gates.
"""
import sys
from itertools import combinations, product

import numpy as np

from adjust import (PathBlowup, is_valid_adjustment_set, optimal_adjustment_set,
                    possibly_causal_paths)
from graphs import (MeekFail, adjacent, apply_background_knowledge,
                    common_orientation_graph, consistent_dag_extensions,
                    dag_to_cpdag, directed_edges, has_directed_cycle,
                    is_consistent, is_directed, is_undirected, meek_closure,
                    random_dag, skeleton, undirected_edges, v_structures)
from x1_ops import (UNREACH, bk_assert, bk_assert_stampall, bload_pool,
                    draw_sloc, draw_suni, hop_dist_from, label_statement,
                    nonadjacent_pairs, pdag_extendable, vstruct_decomposition)

RES = []


def report(name, ok, detail=""):
    RES.append((name, bool(ok), detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}  {detail}", flush=True)


# ---------------------------------------------------------------------------
def T1_dortarsi_vs_bruteforce(rng, n=4000):
    """Dor-Tarsi `pdag_extendable` against the brute-force definition.

    Definition (graphs.consistent_dag_extensions): orient the undirected edges,
    keep only orientations that are acyclic and introduce no unshielded collider
    beyond those already in G.  Extendable iff at least one survives.
    Random PDAGs, INCLUDING objects that are not MPDAGs of any MEC -- which is
    precisely the population arm S produces.
    """
    bad = 0
    tot = 0
    n_ext = 0
    for _ in range(n):
        p = int(rng.integers(3, 7))
        # a deliberately unstructured PDAG: random skeleton, random orientation
        # of each edge into {->, <-, --}
        G = np.zeros((p, p), dtype=np.int8)
        for a, b in combinations(range(p), 2):
            r = rng.random()
            if r < 0.35:
                continue
            k = rng.integers(0, 3)
            if k == 0:
                G[a, b] = 1
            elif k == 1:
                G[b, a] = 1
            else:
                G[a, b] = G[b, a] = 1
        if len(undirected_edges(G)) > 10:
            continue
        tot += 1
        fast = pdag_extendable(G)
        slow = len(consistent_dag_extensions(G, ref_vstructs=v_structures(G),
                                             limit=1)) > 0
        n_ext += int(slow)
        if fast != slow:
            bad += 1
    report("T1 pdag_extendable == brute-force consistent extension exists",
           bad == 0, f"(disagreements {bad}/{tot}; extendable {n_ext}/{tot})")
    return bad == 0


# ---------------------------------------------------------------------------
def T2_bk_assert_reproduces_meek(rng, n=1500):
    """On knowledge that lies on U(C) with TRUE orientations -- the unperturbed
    anchor -- bk_assert must return exactly apply_background_knowledge's graph.
    If it does not, arm S is not comparable to arm R at rho = 0."""
    bad = 0
    tot = 0
    for _ in range(n):
        p = int(rng.integers(4, 9))
        D = random_dag(p, float(rng.choice([1.5, 2.0, 2.5])), rng)
        C = dag_to_cpdag(D)
        U = undirected_edges(C)
        if len(U) < 3:
            continue
        K = [(i, j) if D[i, j] == 1 else (j, i) for (i, j) in U][:4]
        try:
            G_ref = apply_background_knowledge(C, K)
        except MeekFail:
            continue
        tot += 1
        G_new, info = bk_assert(C, K)
        if not np.array_equal(G_ref, G_new):
            bad += 1
        if info["n_added"] or info["conflict"]:
            bad += 1
    report("T2 bk_assert(C, true K on U(C)) == apply_background_knowledge",
           bad == 0, f"(mismatches {bad}/{tot})")
    return bad == 0


# ---------------------------------------------------------------------------
def T3_meek_consistency_is_trivially_false(rng, n=400):
    """THE DESIGN QUESTION, answered by brute force rather than by argument.

    Claim: for a required edge on a pair non-adjacent in C, "is K' consistent
    with C?" is FALSE for every such statement, and it is false for a
    DEFINITIONAL reason -- no DAG in MEC(C) contains the edge at all.  So the
    consistency check carries no information about this class.

    Verified two ways on the same graphs:
      (i)  is_consistent(C, [(a,b)]) is False   for every non-adjacent (a,b)
      (ii) the brute-force MEC of C (all consistent extensions of C) contains no
           DAG with that edge -- because every member has skeleton(C).
    """
    bad_i = bad_ii = 0
    tot = 0
    n_mec = 0
    for _ in range(n):
        p = int(rng.integers(4, 7))
        D = random_dag(p, float(rng.choice([1.5, 2.0, 2.5])), rng)
        C = dag_to_cpdag(D)
        if len(undirected_edges(C)) > 10:
            continue
        N = nonadjacent_pairs(C)
        if not N:
            continue
        mec = consistent_dag_extensions(C)
        n_mec += len(mec)
        for (a, b) in N:
            for (i, j) in ((a, b), (b, a)):
                tot += 1
                if is_consistent(C, [(i, j)]):
                    bad_i += 1
                if any(H[i, j] == 1 for H in mec):
                    bad_ii += 1
    ok = (bad_i == 0 and bad_ii == 0)
    report("T3 Meek-consistency is trivially FALSE for the spurious class",
           ok, f"(is_consistent True {bad_i}/{tot}; MEC members carrying the "
               f"edge {bad_ii}/{tot}; MEC size total {n_mec})")
    return ok


# ---------------------------------------------------------------------------
def T4_operator_does_what_it_says(rng, n=1200):
    """bk_assert must (i) leave the asserted edge present and oriented as
    asserted, (ii) grow the skeleton by exactly the non-adjacent pairs asserted,
    (iii) never touch any other adjacency."""
    bad = 0
    tot = 0
    for _ in range(n):
        p = int(rng.integers(4, 9))
        D = random_dag(p, float(rng.choice([1.5, 2.0, 2.5])), rng)
        C = dag_to_cpdag(D)
        U = undirected_edges(C)
        N = nonadjacent_pairs(C)
        if len(U) < 4 or len(N) < 1:
            continue
        K = [(i, j) if D[i, j] == 1 else (j, i) for (i, j) in U][:4]
        a, b = N[int(rng.integers(len(N)))]
        if rng.random() < 0.5:
            a, b = b, a
        Kp = [(a, b)] + K[1:]
        G, info = bk_assert(C, Kp)
        tot += 1
        if not is_directed(G, a, b):
            bad += 1
        want = skeleton(C).copy()
        want[a, b] = want[b, a] = 1
        if not np.array_equal(skeleton(G), want):
            bad += 1
        if info["n_added"] != 1:
            bad += 1
    report("T4 bk_assert adds exactly the asserted edge, oriented as asserted",
           bad == 0, f"(violations {bad}/{tot})")
    return bad == 0


# ---------------------------------------------------------------------------
def T5_what_is_the_resulting_object(rng, n=1500):
    """MEASUREMENT, not a gate: WHAT is the object bk_assert returns?

    An MPDAG is a PDAG that is (a) extendable and (b) MAXIMALLY oriented -- equal
    to the common orientation of all its consistent extensions.  bk_assert runs
    Meek closure, so it is Meek-closed; but Meek closure is complete only
    relative to a CPDAG plus background knowledge, and the augmented skeleton has
    its OWN v-structures that C never encoded.  So the output need not be the
    MPDAG of any MEC.

    Measured by brute force over the extensions (|U(G)| <= 12):
      extendable        -- admits any consistent DAG extension at all
      sound             -- every directed edge of G is oriented that way in EVERY
                           extension (Meek never over-orients)
      maximally_oriented-- G equals the common orientation graph of its extensions
    """
    tot = ext = sound = maximal = 0
    for _ in range(n):
        p = int(rng.integers(4, 8))
        D = random_dag(p, float(rng.choice([1.5, 2.0, 2.5])), rng)
        C = dag_to_cpdag(D)
        U = undirected_edges(C)
        N = nonadjacent_pairs(C)
        if len(U) < 4 or len(N) < 1:
            continue
        K = [(i, j) if D[i, j] == 1 else (j, i) for (i, j) in U][:4]
        a, b = N[int(rng.integers(len(N)))]
        if rng.random() < 0.5:
            a, b = b, a
        G, info = bk_assert(C, [(a, b)] + K[1:])
        if len(undirected_edges(G)) > 12:
            continue
        tot += 1
        E = consistent_dag_extensions(G, ref_vstructs=v_structures(G))
        if not E:
            if pdag_extendable(G):
                report("T5 INTERNAL: extendable disagrees with enumeration", False)
                return False
            continue
        ext += 1
        com = common_orientation_graph(E, skeleton(G))
        ok_sound = all(is_directed(com, i, j) for (i, j) in directed_edges(G))
        sound += int(ok_sound)
        maximal += int(np.array_equal(com, G))
    print(f"       T5 MEASUREMENT  n={tot}  extendable={ext}/{tot}={ext/max(tot,1):.4f}  "
          f"sound={sound}/{ext}={sound/max(ext,1):.4f}  "
          f"maximally_oriented={maximal}/{ext}={maximal/max(ext,1):.4f}", flush=True)
    report("T5 the resulting object is SOUND (Meek never over-orients)",
           sound == ext, f"(unsound {ext - sound}/{ext})")
    return sound == ext


# ---------------------------------------------------------------------------
def T6_restamp_and_order(rng, n=1500):
    """G7: the re-stamp must be a no-op (meek_closure only orients UNDIRECTED
    pairs, so it cannot flip a stamped edge).  G8: close-after-each vs
    stamp-all-then-close is a MEASUREMENT, printed with its n."""
    tot = pre = cur = dis = 0
    for _ in range(n):
        p = int(rng.integers(4, 9))
        D = random_dag(p, float(rng.choice([1.5, 2.0, 2.5])), rng)
        C = dag_to_cpdag(D)
        U = undirected_edges(C)
        N = nonadjacent_pairs(C)
        if len(U) < 4 or len(N) < 4:
            continue
        K = [(i, j) if D[i, j] == 1 else (j, i) for (i, j) in U][:4]
        sp = draw_suni(C, 4, rng)
        nrep = int(rng.integers(1, 5))
        flip = set(rng.choice(4, size=nrep, replace=False).tolist())
        Kp = [sp[t] if t in flip else K[t] for t in range(4)]
        tot += 1
        G, info = bk_assert(C, Kp)
        pre += int(info["restamp_changed_prefix"])
        cur += int(info["restamp_changed_current"])
        dis += int(not np.array_equal(G, bk_assert_stampall(C, Kp)))
    print(f"       T6 MEASUREMENT  n={tot}  restamp_changed(prefix)={pre}  "
          f"restamp_changed(current)={cur}  "
          f"close-after-each != stamp-all = {dis}/{tot} = {dis/max(tot,1):.4f}", flush=True)
    report("T6/G7 re-stamp is a no-op under both scopes", pre == 0 and cur == 0,
           f"(prefix {pre}, current {cur}, n={tot})")
    return pre == 0 and cur == 0


# ---------------------------------------------------------------------------
def T7_cycles_and_illposed_are_reachable(rng, n=3000):
    """The operator MUST be able to build objects O* is not defined on -- that is
    b-LOAD's real behaviour and A5's finding.  A hand-built witness plus a
    prevalence count on random draws."""
    # hand-built witness: C = 0 -- 1 -- 2 chain, assert 2 -> 0 after 0 -> 1 -> 2
    D = np.zeros((3, 3), dtype=np.int8)
    D[0, 1] = 1
    D[1, 2] = 1
    C = dag_to_cpdag(D)
    G, info = bk_assert(C, [(0, 1), (1, 2), (2, 0)])
    witness = info["cycle"] and has_directed_cycle(G)
    cyc = nonext = tot = 0
    for _ in range(n):
        p = int(rng.integers(5, 9))
        Dr = random_dag(p, float(rng.choice([1.5, 2.0, 2.5])), rng)
        Cr = dag_to_cpdag(Dr)
        U = undirected_edges(Cr)
        N = nonadjacent_pairs(Cr)
        if len(U) < 4 or len(N) < 4:
            continue
        K = [(i, j) if Dr[i, j] == 1 else (j, i) for (i, j) in U][:4]
        sp = draw_suni(Cr, 4, rng)
        Kp = [sp[0]] + K[1:]
        Gr, inf = bk_assert(Cr, Kp)
        tot += 1
        cyc += int(inf["cycle"])
        nonext += int(not pdag_extendable(Gr))
    print(f"       T7 MEASUREMENT  n={tot}  cycle={cyc}={cyc/max(tot,1):.4f}  "
          f"non-extendable={nonext}={nonext/max(tot,1):.4f}  (rho=1, S-uni)", flush=True)
    report("T7 the operator can build cyclic / non-extendable objects", witness,
           f"(witness cycle={witness}; random rho=1 cycles {cyc}/{tot})")
    return witness


# ---------------------------------------------------------------------------
def T8_draw_laws(rng, n=1200):
    """S-uni draws lie in N(C) and are distinct as unordered pairs (M9.1).
    S-loc's pool is b-LOAD's verbatim and its label partition is exhaustive."""
    bad = 0
    tot = 0
    labels = {}
    for _ in range(n):
        p = int(rng.integers(5, 10))
        D = random_dag(p, float(rng.choice([1.5, 2.0, 2.5])), rng)
        C = dag_to_cpdag(D)
        N = nonadjacent_pairs(C)
        if len(N) < 4:
            continue
        pairs = [(a, b) for a in range(p) for b in range(a + 1, p)
                 if len(undirected_edges(C)) >= 0]
        x, y = 0, 1
        tot += 1
        sp = draw_suni(C, 4, rng)
        S = skeleton(C)
        if any(S[a, b] != 0 for (a, b) in sp):
            bad += 1
        if len({frozenset(e) for e in sp}) != 4:
            bad += 1
        pool = bload_pool(D, x, y)
        if any(D[a, b] == 1 for (a, b) in pool):
            bad += 1
        if any(a not in (x, y) for (a, b) in pool):
            bad += 1
        exp = [(a, b) for a in (x, y) for b in range(p) if b != a and D[a, b] != 1]
        if sorted(pool) != sorted(exp):
            bad += 1
        for (a, b) in pool:
            k, q = label_statement(C, D, a, b, x, y)
            labels[k] = labels.get(k, 0) + 1
    print(f"       T8 MEASUREMENT  S-loc pool label counts {labels}", flush=True)
    report("T8 draw laws match their definitions", bad == 0,
           f"(violations {bad}/{tot})")
    return bad == 0


# ---------------------------------------------------------------------------
def T9_hop(rng, n=600):
    """hop_dist_from reproduces run_arm1.hop_dist exactly (so hop_C is E1'-
    comparable), and the sentinel is UNREACH, not inf and not -1."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "e1run", "/home/costaj/latent-causal/e1prime-se/code/run_arm1.py")
    e1 = importlib.util.module_from_spec(spec)
    sys.modules["e1run"] = e1
    spec.loader.exec_module(e1)
    bad = 0
    tot = 0
    seen_sentinel = False
    for _ in range(n):
        p = int(rng.integers(4, 12))
        D = random_dag(p, float(rng.choice([1.0, 1.5, 2.5])), rng)
        C = dag_to_cpdag(D)
        ref = e1.hop_dist(C, [0, 1], p)
        mine = hop_dist_from(skeleton(C), [0, 1], p)
        tot += 1
        if not np.array_equal(ref, mine):
            bad += 1
        if (mine == UNREACH).any():
            seen_sentinel = True
    report("T9 hop_dist_from == E1' run_arm1.hop_dist, sentinel = 1048576",
           bad == 0 and seen_sentinel and UNREACH == 1048576,
           f"(mismatches {bad}/{tot}; sentinel observed {seen_sentinel})")
    return bad == 0 and seen_sentinel


# ---------------------------------------------------------------------------
def T10_unreachable_is_a_structural_zero(rng, n=2500):
    """AUDIT M11 / gate G12, checked as a THEOREM on small graphs BEFORE it is
    used as a gate -- and the check refutes the mitigation as written.

    M11 registers: "unreachable statements have damage rate 0.000 by
    construction", because adding an edge between two vertices in a component
    containing neither X nor Y cannot change cn(X,Y,.), forb or O*.

    That is true OF THE ADDED EDGE and FALSE of the registered operator, because
    K'(F) REPLACES slot t: the true statement (u_t,v_t) is WITHDRAWN.  So three
    quantities are measured, not one:
        raw    O*(S) != O*(G0)   -- what M11's gate G12 actually tests
        drop   O*(D) != O*(G0)   -- the withdrawal alone, no edge asserted
        edge   O*(S) != O*(D)    -- the added edge alone.  THIS is the theorem.
    The gate is `edge == 0`.  `raw` is reported as a measurement.
    """
    raw = drop = edge = 0
    tot = 0
    for _ in range(n):
        p = int(rng.integers(6, 12))
        D = random_dag(p, float(rng.choice([1.0, 1.5, 2.0])), rng)
        C = dag_to_cpdag(D)
        U = undirected_edges(C)
        if len(U) < 4:
            continue
        pairs = [(a, b) for a in range(p) for b in range(p)
                 if a != b and len(possibly_causal_paths(D, a, b)) > 0]
        if not pairs:
            continue
        x, y = pairs[int(rng.integers(len(pairs)))]
        dist = hop_dist_from(skeleton(C), [x, y], p)
        far = [(a, b) for a in range(p) for b in range(a + 1, p)
               if dist[a] == UNREACH and dist[b] == UNREACH
               and skeleton(C)[a, b] == 0]
        if not far:
            continue
        K = [(i, j) if D[i, j] == 1 else (j, i) for (i, j) in U][:4]
        try:
            G0 = apply_background_knowledge(C, K)
        except MeekFail:
            continue
        O0 = optimal_adjustment_set(G0, x, y)
        if O0 is None:
            continue
        a, b = far[int(rng.integers(len(far)))]
        t = int(rng.integers(4))
        Kp = [(a, b) if s == t else K[s] for s in range(4)]
        Kd = [K[s] for s in range(4) if s != t]
        Gs, _ = bk_assert(C, Kp)
        Gd, _ = bk_assert(C, Kd)
        Os = optimal_adjustment_set(Gs, x, y)
        Od = optimal_adjustment_set(Gd, x, y)
        tot += 1
        raw += int(Os != O0)
        drop += int(Od != O0)
        edge += int(Os != Od)
    print(f"       T10 MEASUREMENT  n={tot}  raw(S vs G0)={raw}={raw/max(tot,1):.4f}  "
          f"drop-alone(D vs G0)={drop}={drop/max(tot,1):.4f}  "
          f"ADDED EDGE alone(S vs D)={edge}={edge/max(tot,1):.4f}", flush=True)
    report("T10/G12 an UNREACHABLE added edge cannot move O* (structural zero)",
           edge == 0, f"(edge-attributable moves {edge}/{tot}; "
                      f"raw replacement moves {raw}/{tot} -- the difference is "
                      f"the WITHDRAWN true statement, not the class)")
    return edge == 0


# ---------------------------------------------------------------------------
def T10b_drop_confound_is_general(rng, n=2000):
    """The withdrawal confound is not special to unreachable pairs: it is a
    property of the registered replacement law and it applies to EVERY arm-S
    rate.  Arm R reverses a slot (all 4 pairs still asserted); arm S replaces one
    (3 remain).  Measured at rho = 1 with a UNIFORM S-uni draw."""
    tot = raw = drop = edge = 0
    for _ in range(n):
        p = int(rng.integers(5, 9))
        D = random_dag(p, float(rng.choice([1.5, 2.0, 2.5])), rng)
        C = dag_to_cpdag(D)
        U = undirected_edges(C)
        N = nonadjacent_pairs(C)
        if len(U) < 4 or len(N) < 4:
            continue
        pairs = [(a, b) for a in range(p) for b in range(p)
                 if a != b and len(possibly_causal_paths(D, a, b)) > 0]
        if not pairs:
            continue
        x, y = pairs[int(rng.integers(len(pairs)))]
        K = [(i, j) if D[i, j] == 1 else (j, i) for (i, j) in U][:4]
        try:
            G0 = apply_background_knowledge(C, K)
        except MeekFail:
            continue
        O0 = optimal_adjustment_set(G0, x, y)
        if O0 is None:
            continue
        sp = draw_suni(C, 4, rng)
        t = int(rng.integers(4))
        Gs, _ = bk_assert(C, [sp[t] if s == t else K[s] for s in range(4)])
        Gd, _ = bk_assert(C, [K[s] for s in range(4) if s != t])
        Os = optimal_adjustment_set(Gs, x, y)
        Od = optimal_adjustment_set(Gd, x, y)
        tot += 1
        raw += int(Os != O0)
        drop += int(Od != O0)
        edge += int(Os != Od)
    print(f"       T10b MEASUREMENT  n={tot} (rho=1, S-uni, ALL hops)  "
          f"raw={raw}={raw/max(tot,1):.4f}  drop-alone={drop}={drop/max(tot,1):.4f}  "
          f"edge-alone={edge}={edge/max(tot,1):.4f}", flush=True)
    report("T10b the withdrawal confound is present at every hop, so ARM D is "
           "required", drop > 0, f"(drop-alone O* moves {drop}/{tot})")
    return drop > 0


# ---------------------------------------------------------------------------
def T11_no_tau_leak():
    """G9: tau must not appear on the instrument path."""
    import x1_ops
    import inspect
    src = inspect.getsource(x1_ops)
    ok = "tau" not in src
    report("T11/G9 no tau on the instrument path (x1_ops)", ok)
    return ok


# ---------------------------------------------------------------------------
def T12_path_cap():
    """M16: the cap fires and is catchable, rather than hanging."""
    import adjust
    old = adjust.PATH_CAP
    adjust.PATH_CAP = 5
    p = 8
    G = np.ones((p, p), dtype=np.int8)
    np.fill_diagonal(G, 0)
    fired = False
    try:
        adjust.possibly_causal_paths(G, 0, 7)
    except PathBlowup:
        fired = True
    adjust.PATH_CAP = old
    report("T12/M16 possibly_causal_paths cap raises PathBlowup", fired)
    return fired


# ---------------------------------------------------------------------------
def T13_wilson_and_bootstrap(rng):
    """The two inferential tools, checked against their definitions."""
    from stats_x1 import wilson, cluster_bootstrap_rr
    # Wilson's DEFINING property: the endpoints are the two roots of
    #     (phat - p)^2 = z^2 p(1-p)/n
    # Checked as a root condition rather than against constants I typed in.
    z = 1.959963984540054
    ok1 = True
    for (k, n_) in [(14, 100), (0, 100), (100, 100), (1, 7), (37, 1234)]:
        lo_, hi_ = wilson(k, n_)
        ph = k / n_
        for e in (lo_, hi_):
            if e in (0.0, 1.0):
                continue
            r = (ph - e) ** 2 - z * z * e * (1 - e) / n_
            if abs(r) > 1e-12:
                ok1 = False
    lo, hi = wilson(14, 100)
    lo0, hi0 = wilson(0, 100)
    ok2 = (lo0 < 1e-12) and (0.0 < hi0 < 0.05) and (wilson(0, 0) == (0.0, 0.0))
    # bootstrap: identical arms must give RR ~ 1 with an interval covering 1
    m = 800
    num_r = rng.binomial(4, 0.15, size=m)
    den_r = np.full(m, 4)
    num_s = rng.binomial(4, 0.15, size=m)
    den_s = np.full(m, 4)
    rr, lo_b, hi_b = cluster_bootstrap_rr(num_s, den_s, num_r, den_r, B=2000,
                                          seed=7)
    ok3 = lo_b < 1.0 < hi_b and 0.7 < rr < 1.4
    report("T13 Wilson CI and paired cluster bootstrap match their definitions",
           ok1 and ok2 and ok3,
           f"(wilson(14,100)=[{lo:.5f},{hi:.5f}]; wilson(0,100)=[{lo0:.5f},{hi0:.5f}]; "
           f"null RR={rr:.3f} [{lo_b:.3f},{hi_b:.3f}])")
    return ok1 and ok2 and ok3


# ---------------------------------------------------------------------------
def T14_vstruct_artefact_is_real(rng, n=2000):
    """A4's skeletal artefact, demonstrated rather than asserted: adding {a,b}
    DELETES every unshielded collider of C whose parent pair is exactly {a,b}."""
    hits = tot = 0
    pure = 0
    for _ in range(n):
        p = int(rng.integers(5, 9))
        D = random_dag(p, float(rng.choice([1.5, 2.0, 2.5])), rng)
        C = dag_to_cpdag(D)
        U = undirected_edges(C)
        N = nonadjacent_pairs(C)
        if len(U) < 4 or len(N) < 1:
            continue
        K = [(i, j) if D[i, j] == 1 else (j, i) for (i, j) in U][:4]
        a, b = N[int(rng.integers(len(N)))]
        G, info = bk_assert(C, [(a, b)] + K[1:])
        tot += 1
        dec = vstruct_decomposition(C, G, [(a, b)])
        if info["vstruct_vs_C"]:
            hits += 1
            if dec["n_gained"] == 0 and dec["n_lost_other"] == 0 \
               and dec["n_lost_shielded_by_added"] > 0:
                pure += 1
    print(f"       T14 MEASUREMENT  vstruct_vs_C fired {hits}/{tot}="
          f"{hits/max(tot,1):.4f}; of those, PURELY the shielding artefact "
          f"{pure}/{hits}={pure/max(hits,1):.4f}", flush=True)
    report("T14 vstruct_vs_C is contaminated by the shielding artefact",
           pure > 0, f"(pure-artefact firings {pure})")
    return pure > 0


def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 20260819
    rng = np.random.default_rng(seed)
    sys.setrecursionlimit(20000)
    fns = [T1_dortarsi_vs_bruteforce, T2_bk_assert_reproduces_meek,
           T3_meek_consistency_is_trivially_false, T4_operator_does_what_it_says,
           T5_what_is_the_resulting_object, T6_restamp_and_order,
           T7_cycles_and_illposed_are_reachable, T8_draw_laws, T9_hop,
           T10_unreachable_is_a_structural_zero, T10b_drop_confound_is_general]
    for f in fns:
        f(rng)
    T11_no_tau_leak()
    T12_path_cap()
    T13_wilson_and_bootstrap(rng)
    T14_vstruct_artefact_is_real(rng)
    npass = sum(1 for _, ok, _ in RES if ok)
    print(f"\n==== {npass}/{len(RES)} machinery tests passed ====")
    sys.exit(0 if npass == len(RES) else 1)


if __name__ == "__main__":
    main()
