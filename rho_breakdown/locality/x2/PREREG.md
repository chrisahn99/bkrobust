# X2 PRE-REGISTRATION — locality: can a statement far from the query move `O*`?

**Written 2026-08-19, BEFORE writing or running any X2 code and before any X2 number exists.**

The claim under test is asserted in `~/rho-breakdown-iclr2027/abstract/round2/C-triage.tex` as
though proved:

> **S5** "The optimal adjustment set is a local functional of the graph, and no misstatement at
> graph-distance at least 1 from the query moved it (0 of 791)."
> **S8** "Locality restricts the enumeration to statements near the query, so cost does not scale
> with |K|."
> **S9** "Knowledge outside that neighbourhood need not be elicited."

**Meek closure is the mechanism that propagates orientations across a graph.** "A statement far
from (X,Y) cannot reach `pa(cn(X,Y)) \ forb(X,Y)`" is exactly the kind of statement closure exists
to violate. R1 (theorist, `REVIEWS-round2.md:22-25`) named the proof obligation and it is unmet.
This run does not discharge it. It decides the empirical question that must be settled first, and
it decides **separately** for S5, S8 and S9, because they have different truth conditions and can
have different fates.

---

## 0. 🔴 DISCLOSURE — numbers already seen before this file was written

E1′'s PREREG discloses a 24-SCM smoke run rather than back-dating around it. Same discipline here.
X2 is **not** being designed blind. Everything below was on the record, or was produced by the
scouting phase earlier today, before this file existed. Registering predictions against numbers
one has already seen is worthless; so each item is marked with what it forecloses.

| # | already seen | source | what it forecloses |
|---|---|---|---|
| P1 | E1′ censoring of ρ\*_any at dist-1: `original` 1.000 (n=93), `licensed` 0.935 (n=229), `large` 0.977 (n=638), `k8` \|K\|=8 0.955 (n=22), **`k8` \|K\|=4 0.778 (n=27)** | `../../2026-08-19_e1prime-se-criterion/results/arm1_analysis_*.json` | the *apparent* contradiction. Strata are **per-SCM** (`HOPMIN = min` over all 4 statements, `analyse_arm1.py:87`); S5 is **per-statement**. Not comparable. §1.6 |
| P2 | **0 of 15,578** single flips at `dmin ≥ 1` moved the estimate or changed `O*`, recomputed per-statement from E1′'s raw `arm1_*.json` (original 0/3,908 · licensed 0/4,009 · large 0/7,661) | scouting phase, betelgeuse, today | **the ρ=1 result on E1′'s own K-sampling is already known.** ARM B below is therefore *replication with a corrected denominator*, not discovery. The unrun parts are ARM A, ARM C, ARM D |
| P3 | non-amenability (`E_ident`) fires at `dmin` **exactly 1** and never at ≥ 2: original 109, licensed 129, large 24 | scouting phase, today | the direction of the S9 verdict. §4.3 registers the **bar**, not the direction |
| P4 | all 15 licensed non-censored dist-≥1 SCMs have `d(ρ) ≡ [0,0,0,0,0]` and `ρ*_se = 5` — every one is `ident`-driven, none `se`-driven | scouting phase, today | the resolution of the contradiction. §1.6 |
| P5 | a ρ=2 counterexample exists: seed 508, p=7, deg=2.0, x=4, y=0, τ=−0.138462, K=[(6,5),(1,3),(4,6),(5,2)], flip statements 0 and 1 (both `dmin`=1), O\*={1,5}, invalid in D, \|bias\|/\|τ\|=0.318 | `../../2026-08-14_latent-causal-iclr2027-rerank/RANKING.md:272-278` | **S5-without-a-ρ-qualifier is already false.** The ρ≥2 leg is confirmatory replication + rate measurement, and §4.2 says so |
| P6 | "0/791" does not reproduce from any archived artefact; on-disk counterparts are 654 (E3 binning) / 794 (unreachable folded in) / 1,002 (both arms); retirement ordered 2026-08-14, twice, as mandatory | `RANKING.md:270-272,:387,:825` | X2 **will not print 791**. Whatever X2 measures replaces it, with its own denominators (§1.4) |
| P7 | 258 of the 794 (32.5 %) were Meek-**rejected** and could not have moved `O*` under any circumstance | scouting phase, recomputed from `e3_rows.json` | the denominator discipline of §1.4. A rejected perturbation is not a trial |

**Consequence for the design.** The pre-registered content of X2 is concentrated in the four things
P1–P7 do *not* answer:

1. **ARM A** — an *exhaustive* single-statement search (every undirected edge, **both** orientations,
   no `|K|=4` sampling), which E1′ never ran and which is the only honest form of "does there exist".
2. **ARM C** — an adversarial ensemble that deliberately maximises the distance from the query to
   the causal-path set, i.e. the regime where locality *should* break if it breaks at all.
3. **ARM D** — the **pruning-equivalence** test, which is what S8 actually needs and which no
   experiment in the campaign has ever run.
4. **`E_gain`** — the possibility that a far statement *creates* identification (§2.4). Never
   measured, and if it fires it refutes S9 with a *positive* result rather than a negative one.

---

## 1. The estimand

### 1.1 Objects and notation

All graphs are `p × p` int8 adjacency matrices on `V = {0,…,p−1}`, the codebase's own convention
(`graphs.py:6-9`):

```
G[i,j]=1, G[j,i]=0   :  i -> j        (directed)
G[i,j]=1, G[j,i]=1   :  i -- j        (undirected)
G[i,j]=0, G[j,i]=0   :  no edge

D            the true DAG;  C = dag_to_cpdag(D)  its CPDAG
U(C)       = { {u,v} : C[u,v]=C[v,u]=1 }                     the undirected edges of C
s = (u,v)    a STATEMENT, read "u -> v is required".
             ADMISSIBLE iff {u,v} in U(C).  (graphs.py:179-180 raises MeekFail otherwise —
             see §6.3: the non-adjacent class is X1's, not X2's.)
K            a set of statements on pairwise-distinct pairs, |K| = k
M(K)       = apply_background_knowledge(C, K)   -- orient each s in K, close under Meek R1-R4
             after EACH, then reject if a directed cycle or a NEW v-structure appears.
             M(K) = FAIL denotes MeekFail.
```

For a query `(X,Y)` with `X ≠ Y`, on an MPDAG `G`:

```
pcp(G)     = all simple proper possibly-causal paths X -> … -> Y             (adjust.py:39)
cn(G)      = union over pi in pcp(G) of pi[1:]                               (adjust.py:60)
forb(G)    = poss_de(G, cn(G)) u {X}    ,  = {X} if cn(G) = {}               (adjust.py:64)
pa(S,G)    = { a : a -> b directed, for some b in S }                        (adjust.py:72)
amen(G)    = pcp(G) is non-empty AND every pi in pcp(G) has X -> pi[1] DIRECTED  (adjust.py:82)
O*(G)      = pa(cn(G), G) \ forb(G)      if amen(G)
           = UNDEFINED                    otherwise                          (adjust.py:90)
beta(O)    = population OLS coefficient of X in  Y ~ X + O                   (se.py:23)
SE_n(O)    = sqrt( sigma2(O) * [Sigma_{S S}^{-1}]_{XX} / (n - |S| - 1) ) ,  S = {X} u O
tau        = ((I - A)^{-1})[X,Y]         the TRUE total effect. ORACLE. Used for evaluation
             only, never as an input to any statistic X2 reports as analyst-facing.
```

### 1.2 Distance — stated once, unambiguously, because the whole dispute is here

The walk is on the **undirected skeleton of `C`**, orientation ignored. This is deliberate and it is
the *conservative* choice: skeleton distance is a lower bound on any orientation-respecting
distance, so `{s : d_skel ≥ 1}` is a **superset** of `{s : d_directed ≥ 1}`, and a locality claim
proved on the larger set is the stronger claim.

```
S          = skeleton(C)
delta(w)   = BFS distance in S from the source set {X, Y}
           = 0 iff w in {X,Y}
           = INFINITY if w is in a different connected component of S
dmin(s)    = min( delta(u), delta(v) )      for s = (u,v)   -- PRIMARY
dmax(s)    = max( delta(u), delta(v) )                      -- SECONDARY, reported
dmin(K')   = min over s in flipped(K,K') of dmin(s)         -- for rho >= 2
```

Three consequences fixed now, in writing:

- **`dmin(s) ≥ 1` ⟺ neither endpoint of `s` is `X` or `Y`** ⟺ "the statement is not incident to the
  query". These are the same set; the scouting phase confirmed both give n=794 on E3's data. So
  hypothesis (b) of the caller's brief — *"`stmt_hopdist` is a min over endpoints, so a statement
  with `d(u)=1, d(v)=0` is being scored as near"* — is **not a bug, it is the definition**, and it
  is the only definition under which S5 is a non-trivial claim.
- The **`dmax` reading is reported and is expected to be trivially false**: nearly every statement
  has an endpoint outside `{X,Y}`, so `dmax ≥ 1` is close to the whole population. Printing
  `N_set(dmax ≥ 1)` closes the ambiguity in the paper's own words rather than leaving a reviewer to
  find it. §4.1 registers this as a *reported* output with no verdict attached.
- 🔴 **`skeleton(C) == skeleton(D)` always** (`dag_to_cpdag` starts from `skeleton(D)` and only ever
  zeroes one direction), and **`skeleton(M(K')) == skeleton(C)` for every consistent `K'`** for the
  same reason. Distances computed once on `C` are therefore valid for every member of the ball.
  This is asserted as **integrity gate G2** rather than assumed.
- 🔴 **The unreachable sentinel.** `hop_dist` (`run_arm1.py:49`) returns `1 << 20 = 1048576`, not
  `inf` and not `−1`. E1′'s own published `large` row "dist ≥ 2 | 1.000 (n=331)" silently contains
  **57 SCMs at distance 1048576**, i.e. 17 % of that denominator is *infinite* distance, not "two or
  more hops". X2 maps the sentinel to `float('inf')` at ingest and gives it its **own reported row,
  labelled `disconnected`, never merged into any `≥ k` bucket.** Gate **G7**.

### 1.3 The event taxonomy — six mutually informative outcomes, never merged

Given a base `K` with `M(K) ≠ FAIL` and a perturbed `K'`, write `G0 = M(K)`, `G' = M(K')`,
`O0 = O*(G0)`, `O' = O*(G')`:

| id | condition | channel | who sees it |
|---|---|---|---|
| `E_meek` | `M(K') = FAIL` | caught free | the analyst (Meek rejects) |
| `E_ident` | consistent, `amen(G0)=True`, `amen(G')=False` | **LOUD** — identification lost | the analyst |
| `E_gain` | consistent, `amen(G0)=False`, `amen(G')=True` | **LOUD** — identification *gained* | the analyst |
| `E_set` | consistent, both amenable, `O' ≠ O0` | **SILENT** — the S5 notion **(a)** | nobody |
| `E_bias` | `E_set` ∧ `O'` not a valid adjustment set for (X,Y) in `D` | **SILENT** — notion **(b)**, oracle-adjudicated | nobody |
| `E_est(n,z)` | `E_set`-eligible ∧ `\|beta(O') − beta(O0)\| > z·SE_n(O0)` | notion **(c)**, analyst-facing | nobody |

`E_est(∞, z) ⟺ beta(O') ≠ beta(O0)` exactly, because `SE_∞ = 0` (`se.py:47`). **All decisive X2
statistics are computed at `n = ∞`.** Rationale: only at `SE = 0` does the analyst-facing notion (c)
coincide with the oracle notion (a)/(b) — and it does coincide, because `changed_still_valid = 0`
on 27,000+ pilot perturbations and on all six E1′ ensemble/arm cells, which makes
`O' ≠ O0 ⟺ O' invalid` empirically. At finite `n` a counterexample can be dismissed as a
thresholding artefact; at `n = ∞` it cannot. **`n = 20,000` is reported alongside, never instead.**

### 1.4 The estimand proper, with its denominators

For a distance radius `r ≥ 0` and ball radius `ρ ≥ 1`:

```
Omega(r,rho) = { (SCM, K') : |flipped(K,K')| = rho , dmin(K') >= r }

N_set (r,rho) = # { w in Omega(r,rho) : E_set  }
N_ident(r,rho)= # { w in Omega(r,rho) : E_ident}
N_gain (r,rho)= # { w in Omega(r,rho) : E_gain }
N_bias (r,rho)= # { w in Omega(r,rho) : E_bias }

D_all (r,rho) = |Omega(r,rho)|                                     every trial attempted
D_con (r,rho) = # { w : M(K') != FAIL }                            survived Meek
D_am  (r,rho) = # { w : M(K') != FAIL and amen(G0) and amen(G') }  could have moved O*

L_set(r,rho)  = N_set(r,rho) / D_am(r,rho)                         THE reported rate
```

🔴 **Denominator discipline, binding.** `D_am` is the only denominator under which `L_set` is a
rate of the thing being claimed: a Meek-rejected perturbation never produces an MPDAG and cannot
move `O*`, so counting it inflates the denominator and makes a zero look stronger than it is. The
retired "0 of 791" did exactly this — 258 of its 794 (32.5 %) were Meek-rejected. **Every X2 table
prints `N`, `D_am`, `D_con` and `D_all` on the same line. No fraction is emitted without all four.**
Gate **G8**.

### 1.5 The locality radius — the number S9 needs and has never been given

```
r*_set    = max { dmin(K') : E_set   fires }        , -1 if the set is empty
r*_ident  = max { dmin(K') : E_ident fires }        , -1 if empty
r*_gain   = max { dmin(K') : E_gain  fires }        , -1 if empty
r*_elicit = max( r*_set , r*_ident , r*_gain )
```

`r*_elicit` is the **elicitation radius**: statements strictly beyond it never changed the analysis
in any way the analyst or the oracle can see. S9 as written asserts `r*_elicit = 0`. Reported per
ρ, per ensemble, and pooled — never pooled *only*.

### 1.6 The apparent contradiction, and how X2 dissolves it in code rather than in argument

E1′'s 0.935 at "dist 1" and S5's "0 of 791" differ on **four** axes simultaneously:

| axis | E1′'s 0.935 | S5's 0/791 | X2 |
|---|---|---|---|
| unit | **per SCM** (`HOPMIN = min` over all 4 statements) | per statement | **per statement**, and the per-SCM roll-up printed beside it |
| ball | pooled ρ = 1…4 | ρ = 1 | **each ρ separately**, never pooled |
| event | `ρ*_any = min(ρ*_se, ρ*_ident)` — silent ∪ **loud** | `E_set` / `E_bias` — silent only | **six events, never merged** (§1.3) |
| yardstick | `z·SE_n` at n = 2×10⁴ | exact, population | **n = ∞ decisive**, n = 2×10⁴ reported |

**Registered as a step, not as a conclusion:** X2 **Step 0** recomputes, from
`~/latent-causal/e1prime-se/results/arm1_*.json` on betelgeuse, the three-way split of every
non-censored `dmin ≥ 1` SCM into `{ident-only, ρ*_se = 1, ρ*_se ≥ 2}`, at `n = ∞`. Attribution of
a `members[]` index to its flip is exact: members are appended in `for rho: for flip in
combinations(range(k), rho)` order and **both** branches append, so index ↔ flip is a bijection
(verified 0 length-mismatches over 7,167 SCMs in the scouting phase). Step 0's output is a table,
not a verdict. If it disagrees with P2/P4 the discrepancy is the finding and X2 halts to report it.

---

## 2. The operators — unambiguous enough that two people implement the same thing

Base machinery is **copied unmodified** from `~/latent-causal/e1prime-se/code/` into
`~/latent-causal/x2-locality/code/` (`graphs.py`, `adjust.py`, `scm.py`, `se.py`). Neither
`e1prime-se/` nor `rho-breakdown-knowledge/` is edited. The perturbation operator is the pilot's
own and is **unchanged**: reverse the orientation of statements already on `U(C)`.

### 2.1 ARM A — SOLO. *Does there exist a single distant statement that moves `O*`?*

**This is the primary existence test and it is the arm no prior experiment ran.**

```
for each SCM (D, C, X, Y, Sigma):
    require amen(C) is not needed here -- both branches are informative (see E_gain)
    G0 := C            (the EMPTY knowledge set; K = {})
    O0 := O*(C)        (may be UNDEFINED)
    for each {u,v} in U(C):
        for s in [ (u,v), (v,u) ]:            # BOTH orientations, exactly one is true in D
            K' := { s }
            classify (E_meek / E_ident / E_gain / E_set / E_bias / E_est) per §1.3
            record dmin(s), dmax(s), delta(u), delta(v),
                   true_in_D := (D[s[0], s[1]] == 1),
                   cascade   := |directed_edges(M(K'))| - |directed_edges(C)| - 1
```

Three things this buys that ARM B cannot:

- **Exhaustive.** E1′ samples 4 of the `|U(C)|` undirected edges at random. Rare far edges are
  systematically under-sampled precisely where the counterexample would live. ARM A enumerates all.
- **Both orientations.** The pilot and E1′ only ever assert the **true** orientation and then flip
  it. So the campaign has never asked what a *correct* distant statement does. For S8 and S9 the
  correct one matters just as much: if a true distant statement moves `O*`, the pruning is invalid
  and the elicitation advice is wrong, regardless of whether anyone lied. **Reported split by
  `true_in_D`.**
- **No confound from a base `K`.** The baseline is `C` itself, so any observed change is attributable
  to `s` alone, with no interaction with three other statements.

### 2.2 ARM B — CONTEXT. *Replication of the E1′/pilot operator, per statement, with honest denominators.*

`K` = 4 statements drawn from `U(C)` oriented to agree with `D` (`run_arm1.py:139`), permuted by the
same RNG stream, `ρ = 1…min(4, |K|)`, `K'` obtained by reversing the ρ selected statements. Identical
to `build_arm` (`run_arm1.py:63-109`) except that X2 additionally records, per member, the flip index
set, `dmin`, `dmax` and the six events of §1.3. **Job seed `20260819`, identical to E1′**, so the
SCMs are the same objects and G3 can compare them exactly.

### 2.3 ARM C — ADVERSARIAL. *The regime where locality should break if it breaks at all.*

Mechanistic reasoning, registered before the run: `O* = pa(cn) \ forb`, and `cn` is the set of
nodes on proper possibly-causal `X → Y` paths. A node at position `i` on a path of length `L` sits
at skeleton distance `≤ min(i, L−i)` from `{X,Y}` — which is **unbounded in `L`**. So `O*` is *not*
confined to a bounded neighbourhood of the query; it is confined to a neighbourhood of the
**causal-path set**. E3 already found exactly this on a per-member feature ("the predictive content
is causal-path locality, not `pa(cn)\forb` membership; the latter is anti-predictive, AUC 0.409").
**If S5's zero is an artefact of short causal paths in small ER graphs, ARM C is where it dies.**

ARM C differs from ARM A only in the SCM draw:

- `p ∈ 12…20`, `deg ∈ {1.5, 2.0, 2.5}` (high degree is *excluded* — `possibly_causal_paths` is an
  uncapped DFS, `adjust.py:39`, and p=25 × deg=6 does not finish; see §3.1).
- the query `(X,Y)` is chosen **not uniformly** but as `argmax` of `delta_skel(X,Y)` over pairs with
  at least one possibly-causal path in `D`, ties broken by the RNG. Registered rationale: S8's
  pruning claim must hold for **every** query, not on average, and a uniform draw over a small ER
  graph puts `X` and `Y` close together and makes `cn` small.
- `need_u ≥ 4`, and a secondary variant `need_u ≥ 8` (dense equivalence classes — the scouting phase
  flags `k8`/`K4` as the weakest cell in the entire corpus, dist-1 censoring 0.778 at n=27, and it
  is the row E1′'s RESULTS.md omitted).

Both ARM A and ARM C run the ARM-A operator. ARM C is a different ensemble, not a different operator.

### 2.4 `E_gain` — registered as a first-class outcome, not a curiosity

If `amen(C) = False` (the CPDAG alone does not identify the effect) and a **single distant**
statement makes `M({s})` amenable, then that statement is *exactly the knowledge that must be
elicited*, and S9 ("knowledge outside that neighbourhood need not be elicited") is false in the
most damaging possible way: not because far knowledge is dangerous, but because far knowledge is
**valuable and the paper tells the practitioner to skip it.**

Mechanism that would produce it: `s` orients a distant edge, R1 fires along a chain, and the
orientation cascades into `X — v1`, making every proper possibly-causal path start `X → v1`.
`cascade` is recorded per member precisely to exhibit the chain length when it happens.

Note `E_gain` is unreachable in ARM B by construction, because `build_arm` returns early with
`mpdag_amenable=False` and no members whenever `O*(M(K))` is `None` (`run_arm1.py:71-72`). **ARM A
is the only place `E_gain` can be observed**, which is one more reason ARM A is the primary arm.

### 2.5 ARM D — PRUNING EQUIVALENCE. *What S8 actually claims, tested directly.*

No prior experiment tests S8. `build_arm` enumerates `combinations(range(k), rho)` over **all** of
`K`, unconditionally; no distance filter exists in either codebase. So "cost does not scale with
|K|" currently describes an algorithm that has never been written. X2 tests its **validity**
(§6.1 states plainly that it does not test its wall-clock).

The formal property S8 requires:

> **(P-r) — `r`-locality closure.** For every `K'` in the ball, `O*(M(K'))` depends on `K'` only
> through `K' ∩ {s : dmin(s) ≤ r}`.

Operationalised on ARM B's already-enumerated ball (so ARM D costs nothing extra — it is a
re-partition, not a new run):

```
near_r(K)      = { s in K : dmin(s) <= r }
B_rho(K)       = { K' : |flipped(K,K')| <= rho }                            , includes K itself
B^r_rho(K)     = { K' in B_rho(K) : K' agrees with K on every s NOT in near_r(K) }

I_full (rho)   = [ min , max ] of beta(O*(M(K'))) over K' in B_rho(K) that are consistent+amenable
I_prune(r,rho) = the same over B^r_rho(K)

Q_r(rho)       = P( I_prune(r,rho) != I_full(rho) )   , endpoints compared to 1e-9 ABSOLUTE
Cost_r(rho)    = |B^r_rho(K)| / |B_rho(K)|            , the advertised saving
```

`r ∈ {−1, 0, 1, 2}`. `r = −1` prunes everything (`near = {}`, so `I_prune = {beta(O0)}`) and is the
**dead-switch placebo**, gate G5. **S8 is bought only if some `r` satisfies `Q_r ≈ 0` AND
`Cost_r ≪ 1` simultaneously.** A radius large enough to be safe but too large to prune anything buys
nothing, and that trade-off — not either number alone — is the reportable object.

⚠️ Declared **definitional, not evidence**: `I_prune(r,ρ) ⊆ I_full(ρ)` always, since `B^r_ρ ⊆ B_ρ`.
So `Q_r` is one-sided by construction: pruning can only ever make the interval too *narrow*,
never too wide. That is the dangerous direction (an under-wide ignorance interval is a false
reassurance), which is why it is the direction measured.

---

## 3. Ensembles, seeds, counts, power

| id | grid | query draw | `need_u` | K | ρ | seed | target SCMs | speaks to |
|---|---|---|---|---|---|---|---|---|
| A-`licensed` | p 5–10, deg {1.5,2,2.5,4,6} | uniform | ≥ 3 | ∅ (solo) | 1 | 20260820 | 6,000 | S5, S9 |
| A-`large` | p 15–25, deg {1.5,2,2.5} | uniform | ≥ 3 | ∅ | 1 | 20260820 | 4,000 | S5, S9 |
| A-`k8` | p 5–10, deg {1.5,2,2.5,4,6} | uniform | ≥ 8 | ∅ | 1 | 20260820 | 2,000 | S5, S9 |
| B-`licensed` | p 5–10, deg {1.5,2,2.5,4,6} | uniform | ≥ 3 | 4 true | 1–4 | **20260819** | 4,000 | S5, S8, S3/S4 |
| B-`large` | p 15–25, deg {1.5,2,2.5} | uniform | ≥ 3 | 4 true | 1–4 | **20260819** | 4,000 | S5, S8 |
| B-`k8` | p 5–10, deg {1.5,2,2.5,4,6} | uniform | ≥ 8 | 4 true | 1–4 | **20260819** | 1,200 | S5, S8 |
| C-`deep` | p 12–20, deg {1.5,2,2.5} | **argmax δ(X,Y)** | ≥ 4 | ∅ | 1 | 20260821 | 4,000 | S5, S8, S9 |
| C-`deep8` | p 12–20, deg {1.5,2,2.5} | **argmax δ(X,Y)** | ≥ 8 | ∅ | 1 | 20260821 | 2,000 | S5, S9 |

**Seeds.** ARM B reuses **20260819** deliberately: identical job stream ⇒ identical `(p, deg, seed)`
triples ⇒ the same SCMs as E1′, so gate G3 is an exact equality test rather than a distributional
one. ARMs A and C use **20260820 / 20260821** so that a zero from ARM A is *not* a zero on the same
draws that already produced P2. Reusing E1′'s seed for the discovery arm would make the independent
search a re-read of the same 15,578 trials.

**Power.** ARM A yields `2·|U(C)|` trials per SCM. Conservative floors: `|U| ≥ 3` gives ≥ 6 trials at
p 5–10 (observed mean ≈ 5 undirected edges ⇒ ≈ 10 trials); `large` has many more. Assume 40 % of
trials sit at `dmin ≥ 1` at p 5–10 and 75 % at p 15–25 (E1′'s own distance histograms; used for
*sizing* only, not registered as a prediction):

```
A-licensed  6,000 x 10 x 0.40  ~=  24,000 far trials
A-large     4,000 x 24 x 0.75  ~= 72,000
A-k8        2,000 x 16 x 0.40  ~= 12,800
C-deep      4,000 x 20 x 0.85  ~= 68,000
C-deep8     2,000 x 20 x 0.85  ~= 34,000
                                 --------
POOLED ARM A + ARM C  ~= 210,000 far single-statement trials
```

**Rule-of-three bound.** If `N_set(r ≥ 1, ρ=1) = 0` over `D_am` far trials, the exact one-sided
95 % upper bound on the per-trial rate is `3 / D_am`. Registered targets:

| pooled `D_am` at `dmin ≥ 1` | licensed bound |
|---|---|
| ≥ 150,000 | **≤ 2.0 × 10⁻⁵** |
| ≥ 50,000 | ≤ 6.0 × 10⁻⁵ |
| < 20,000 | too weak to license S8 — verdict downgraded to SUPPORTED-WEAK, §4.1 |

Per-ensemble bounds are reported separately; the pooled bound is reported **in addition**, never
instead, because pooling across ensembles hides exactly the p-dependence S5's "⚠️ p ≤ 8 only" caveat
is about.

**Budget.** ≤ 20 min wall-clock per ensemble on 12 workers (`~/e1venv/bin/python`, 16 cores).
E1′'s whole pipeline — 13,554 SCMs × ~15 members plus arm2 — ran in under 15 min, and ARM A's
per-trial cost is ~1 Meek closure + 1 `O*` (cheaper than a `build_arm` member, which also does two
matrix inversions).

### 3.1 The path-enumeration guard — fixed now, because it is a selection rule

`possibly_causal_paths` (`adjust.py:39-59`) is an **uncapped DFS over all simple possibly-directed
paths** and is called several times per member, inside `causal_nodes`, `forb`, `is_amenable` and
`optimal_adjustment_set`. Its cost is exponential in `p` and degree; E1′ keeps it finite only by
never pairing `p ≥ 15` with `deg ≥ 4`. ARM C deliberately pushes `p` up *and* selects the query with
the **longest** `δ(X,Y)`, which is the worst case for path count. A guard is therefore mandatory —
and a guard is a **selection rule on the ensemble**, not an implementation detail, so it is
registered here rather than chosen at the keyboard:

```
PATH_CAP = 200_000     paths enumerated per (G, X, Y) call
SCM_CAP  = 20.0        seconds of wall-clock per SCM
```

An SCM that trips either cap is **aborted whole** — never partially recorded, because a partial
record would silently drop its far statements and bias `L_set` toward zero in exactly the regime
the search is aimed at. Aborted SCMs are counted, and per gate **G9** the count is printed per
ensemble; **> 5 % aborts in any ensemble invalidates that ensemble's rule-of-three bound**, because
the aborted SCMs are precisely the large, dense, long-path ones. If `C-deep` aborts above 5 %, the
registered response is to lower `p` to 12–16 and re-run, reporting both, not to raise the cap.

---

## 4. Pre-registered decision rules — numeric, fixed now

### 4.1 S5 — "no misstatement at graph-distance ≥ 1 from the query moved `O*`"

Primary statistic: `N_set(r ≥ 1, ρ = 1)` at `n = ∞`, pooled over ARM A + ARM C, and reported per
ensemble. Existence questions are asymmetric, so the rule is asymmetric.

- **REFUTED** iff `N_set(r ≥ 1, ρ = 1) ≥ 1`. **One counterexample is enough.** It must be
  regenerable from its `(seed, p, deg)` triple and is reported as §4.5 requires.
- **SUPPORTED** iff `N_set = 0` **and** pooled `D_am(r ≥ 1, ρ = 1) ≥ 150,000`. Reported as
  *"0 of D_am, 95 % upper bound 3/D_am"*, with `D_con` and `D_all` on the same line.
- **SUPPORTED-WEAK** iff `N_set = 0` and `20,000 ≤ D_am < 150,000`. The abstract may then say
  *"empirically zero on N trials at p ≤ 25"* and may not say *"local functional"* without the
  theorem.
- **INCONCLUSIVE** iff `D_am < 20,000` (an ensemble failed to produce far statements) — the honest
  report is that the search was too small, not that locality holds.

Reported alongside, with **no verdict attached** because they are known or definitional:

- `N_set(dmax ≥ 1, ρ = 1)` — the other reading of "graph-distance ≥ 1"; expected to be large, and
  printed so the ambiguity is closed in the paper rather than by a reviewer (§1.2).
- `N_set(r ≥ 1, ρ = 1)` split by `true_in_D` — a *correct* far statement that moves `O*` is not a
  "misstatement" and so does not refute S5, **but it does refute S8 and S9** (§4.2, §4.3).
- The E1′ replication (ARM B), which P2 already places at 0 — printed as replication, labelled as
  such, and **not** counted as independent evidence.

🔴 **Regardless of the verdict, S5's wording must change**, and this is registered now so the run
cannot be read as vindicating the sentence as printed: P5 is a ρ=2 counterexample already on the
record. **The strongest sentence any X2 outcome can license is scoped to `ρ = 1`.** Mandated
wording, pre-committed: *"single orientation errors at skeleton graph-distance ≥ 1 from the query
never changed `O*` (0 of ⟨D_am⟩); joint errors did (⟨N_set(r≥1, ρ≥2)⟩ of ⟨D_am⟩)."*

### 4.2 S8 — "locality restricts the enumeration, so cost does not scale with |K|"

Primary statistic: `Q_r(ρ = 4)` on ARM D, per ensemble, `n_SCM ≥ 2,000` per cell.

- **S8 SUPPORTED at radius `r`** iff `Q_r ≤ 0.001` **and** `Cost_r ≤ 0.5` on **every** ensemble.
  (Both halves required: a radius that is safe but prunes nothing is not a saving.)
- **S8 REFUTED at radius `r`** iff `Q_r ≥ 0.01` on any ensemble.
- **INCONCLUSIVE** in between → the abstract must be rewritten as a checked heuristic with the
  measured `Q_r`, not as a property.

**Registered prediction, direction fixed in advance.** `Q_0 ≥ 0.01` on `large` and on `k8`
(**REFUTED at r = 0**), because P5's seed-508 counterexample has both flipped statements at
`dmin = 1`, so an `r = 0` pruner discards exactly the pair that moves the interval. `Q_1` is
genuinely unknown: predicted `Q_1 ≤ 0.001` on `licensed` and `Q_1 > 0.001` on `large`, on the
causal-path-length mechanism of §2.3.

**Named surprises.**
- `Q_0 < 0.001` on every ensemble ⇒ P5 does not replicate under this operator/denominator and the
  ρ=2 counterexample must be re-audited before anything is written.
- `Q_1 ≥ 0.01` on `licensed` ⇒ locality fails at radius 1 in the paper's own primary box, and S8/S9
  are not repairable by re-scoping; they are cut.

### 4.3 S9 — "knowledge outside that neighbourhood need not be elicited"

S9 is **not** S5. A statement that destroys or creates identification has changed the analysis
materially without ever moving `O*`. Primary statistic: `r*_elicit` (§1.5) and the counts
`N_ident(r ≥ 1, ρ = 1)`, `N_gain(r ≥ 1, ρ = 1)`.

- **S9-as-written REFUTED** iff `N_ident(r ≥ 1, ρ=1) + N_gain(r ≥ 1, ρ=1) ≥ 10` pooled, with
  `D_am ≥ 20,000`. The correction is then quantitative: report `r*_elicit` and rewrite S9 as
  *"knowledge more than `r*_elicit` hops from the query never changed the analysis"*.
- **S9 SUPPORTED** iff `r*_elicit = 0`, i.e. all three counts are zero at `dmin ≥ 1`.
- **S9 PARTIAL** iff `r*_elicit = 1` and the counts at `dmin ≥ 2` are all zero — S9 survives only
  with an explicit radius of 1, and the word "neighbourhood" must be given that number in the text.

⚠️ **Disclosure (P3):** `E_ident` is already known to fire at `dmin = 1` on E1′'s data (109/129/24).
So the *direction* of this verdict is largely foreclosed and X2 claims no discovery for it. What is
**not** foreclosed and is registered here: (i) whether it ever fires at `dmin ≥ 2`, (ii) `E_gain`,
which no experiment has ever measured, (iii) the rate with an honest `D_am` denominator.

### 4.4 The `E_gain` prediction — registered, direction fixed

**PREDICTED:** `N_gain(r ≥ 1, ρ = 1) ≥ 1` on `A-large` or on `C-deep`, i.e. at least one distant
single statement makes a non-identified query identified via a Meek cascade.
**SURPRISE (would refute the mechanism):** `N_gain(r ≥ 1, ρ = 1) = 0` pooled over ≥ 100,000 far
trials on ensembles that contain ≥ 1,000 SCMs with `amen(C) = False`. Meek closure would then be
shown not to propagate amenability across distance at all, which is a *stronger* locality result
than anything the paper currently claims and should be written as such.

### 4.5 What a counterexample must contain to count

Not a row in a table. A counterexample is only accepted if it is exhibited as:

1. `(seed, p, deg, ensemble)` and a one-line regeneration command;
2. the CPDAG `C` as an **edge list** (`i -> j` / `i -- j`), plus `D` as an edge list;
3. `(X, Y)`, `delta(·)` for every vertex, the statement `s`, `dmin(s)`, `dmax(s)`;
4. `O*(G0)` and `O*(G')` as explicit vertex sets, and `beta` for each;
5. **the Meek rule chain**: the ordered list of `(rule, a, b)` orientations fired by the closure,
   showing how the assertion reached `pa(cn) \ forb`;
6. whether `O*(G')` is valid in `D`, and `|beta(O') − tau| / |tau|`;
7. **minimality**: the smallest `p` at which the phenomenon reproduces, found by re-running the
   search restricted to `p ≤ p_found − 1` and reporting whether it recurs.

If a counterexample is found, X2 re-runs a targeted minimal search over `p ∈ 5…8` (cheap,
exhaustive over ARM A's operator) before reporting, so the paper can print the smallest instance.

---

## 5. Declared definitional — NOT evidence, written down before the run

1. `I_prune(r,ρ) ⊆ I_full(ρ)` always. `Q_r` is one-sided by construction (§2.5).
2. `E_est(∞, z) ⟺ beta(O') ≠ beta(O0)`, because `SE_∞ = 0` (`se.py:47`). The `n = ∞` column is not
   a finding about sample size.
3. `skeleton(M(K')) = skeleton(C)` for every consistent `K'`; distances are ball-invariant. Used as
   gate G2, never as a result.
4. `dmin(s) ≥ 1 ⟺ s not incident to {X,Y}`, since `delta(w) = 0 ⟺ w ∈ {X,Y}`. These are the same
   set and only one of them will be reported as a discovery: neither.
5. `E_set ⟹ E_est(∞)` is **not** definitional and must not be asserted — `O*` can change while
   `beta` stays numerically equal. Both are computed and any discrepancy is reported.
6. `E_bias ⟹ E_set` by construction. `E_bias/E_set = 1` empirically in every prior run
   (`changed_still_valid = 0`); X2 recomputes it and, if it stays 1, reports that the three
   phrasings *"O\* changed"*, *"the estimate moved"* and *"silent bias"* are **one measurement**,
   not three pieces of evidence.

---

## 6. Integrity gates — all must pass or no verdict is readable

| id | gate | pass condition |
|---|---|---|
| **G1** | Meek idempotence: `meek_closure(meek_closure(G)) == meek_closure(G)` | 100 % of drawn graphs |
| **G2** | skeleton invariance: `skeleton(M(K')) == skeleton(C)` for every consistent `K'` | 100 % |
| **G3** | **replication**: ARM B on `licensed`/`large`/`k8` at seed 20260819 reproduces E1′'s `stmt_hopdist` arrays and `O0` sets exactly on the overlapping SCMs | **0 mismatches** |
| **G4** | 🔴 **positive control / dynamic range**: the identical detector must find `E_set` events at `dmin = 0` | `L_set(0, ρ=1)` ∈ **[0.05, 0.35]** on `B-licensed`. Outside that band the detector is not demonstrated to work and **no zero at `dmin ≥ 1` is readable** |
| **G5** | **dead-switch placebo**: `Q_{−1}(ρ=4)` must equal the independently computed fraction of SCMs whose `I_full(4)` is non-degenerate | exact equality |
| **G6** | **seed-508 anchor**: regenerate P5 from the *pilot's* generator (`rho-breakdown-knowledge/run_linear.py` stream, **not** `run_arm1.py`'s — different RNG order) and confirm `O* = {1,5}`, invalid, `\|bias\|/\|τ\| = 0.318` | reproduces, or the ρ≥2 leg is reported as **unanchored** |
| **G7** | **sentinel**: no trial with `delta = 1<<20` enters any `≥ k` bucket; `disconnected` is its own reported row | 0 leaks, count printed |
| **G8** | **denominator disclosure**: every emitted fraction carries `N`, `D_am`, `D_con`, `D_all` | enforced in the writer, not by review |
| **G9** | **abort accounting**: SCMs killed by the path-enumeration guard (§6.2) are counted and reported per ensemble | count printed; > 5 % of draws in any ensemble invalidates that ensemble's bound |

**G4 is the gate that decides whether X2 is worth reading.** A search that reports zero
counterexamples is worthless unless the same code path is shown to find them where they are known
to be. The `dmin = 0` band [0.05, 0.35] is deliberately wider than E1′'s own G1 band [0.12, 0.18]
because X2's denominator is `D_am` where E1′'s was the consistent-member count; the band is a
functioning-detector check, not a reproduction check.

---

## 7. What this run does NOT cover — in writing, so a gap cannot later read as coverage

1. 🔴 **It is not a proof.** An empirical zero at `D_am = 200,000` is not the theorem R1 asked for
   (*"no distant statement can enter the neighbourhood through Meek closure"* — `REVIEWS-round2.md:22`).
   S8 and S9 stated as properties require that proof. X2 can only ever license *"checked heuristic,
   0 of N at p ≤ 25"* or *"refuted, here is the instance"*.
2. 🔴 **S8 is only half-tested.** X2 tests the **validity** of pruning (`Q_r`), not its **cost**.
   No pruned enumerator is implemented and no wall-clock benchmark is run. "Cost does not scale
   with |K|" therefore remains unbenchmarked even if `Q_r = 0`. `Cost_r` is a *combinatorial*
   ratio `|B^r_ρ| / |B_ρ|`, not a timing.
3. 🔴 **The X1 class is excluded and it can invalidate this entire analysis.** `graphs.py:179-180`
   refuses a statement on a non-adjacent pair, so "required edges asserted between non-adjacent
   variables" (S7) cannot be perturbed here. This is not merely uncovered: **a spurious edge adds
   to the skeleton and therefore changes `delta(·)` itself**, so every distance in X2 is computed
   under the assumption that the skeleton is correct. If X1 finds that class is *protected* rather
   than caught, the locality analysis must be redone against the perturbed skeleton. Registered
   dependency, not a footnote.
4. **The CPDAG is the oracle's** (`C = dag_to_cpdag(D)`). An analyst runs PC/GES and gets `Ĉ ≠ C`.
   Distances, `U(C)` and admissibility are all defined against the true CPDAG. Same limitation as
   E1′ `AUDIT.md` A2; `causal-learn` is absent from every local env.
5. **Scope:** linear-Gaussian iSCM, ER DAGs `p ∈ [5,25]`, expected degree ≤ 6, single-node `X`,
   single-node `Y`, orientation statements only, population `Σ` (`n = ∞` decisive, `n = 2×10⁴`
   reported). Not licensed: nonlinear mechanisms, latent confounding, multi-node interventions,
   non-Gaussian errors, real or semi-synthetic data.
6. **Tiered knowledge is not tested.** Bang & Didelez tiering exists in the pilot
   (`run_linear.py:57-79`) and was dropped by E1′. X8, not X2.
7. **No p-value, ever.** `Q_r` and `L_set` are worst-case-over-a-ball sensitivity quantities in the
   Rosenbaum-Γ tradition, maximised over up to `2^|K| − 1` members. Read as a test this is
   uncorrected multiplicity. It is not a test.
8. **`ρ > 4` is not explored.** The ball is capped at `|K| = 4` (the pilot's `MAX_K`). A pruning
   failure that first appears at ρ = 5 would not be seen.
9. **ARM C's `argmax δ(X,Y)` query rule makes it non-comparable to ARM A/B by design.** Its rates
   may not be pooled with the uniform-draw ensembles and will not be.

---

## 8. Which abstract sentence each output speaks to

| output | sentence | if it fires | if it does not |
|---|---|---|---|
| `N_set(r ≥ 1, ρ=1)` ≥ 1 | **S5** | S5 **cut**. S8 and S9 fall with it. The paper's locality leg becomes a negative result with an exhibited instance | S5 survives **only** rescoped to ρ=1 and to `p ≤ 25`, with `0 of D_am` replacing 791 (P6) |
| `N_set(r ≥ 1, ρ≥2)` ≥ 1 | **S5** | already expected (P5): S5 **must** gain a ρ=1 qualifier regardless of everything else | P5 fails to replicate → re-audit the 08-14 counterexample before writing anything |
| `Q_0`, `Q_1`, `Cost_r` | **S8** | `Q_r ≥ 0.01` ⇒ S8 **cut** or demoted to a checked heuristic with its measured failure rate | `Q_r ≤ 0.001` with `Cost_r ≤ 0.5` ⇒ S8 licensed as *validity only*; the cost half still unbenchmarked (§7.2) |
| `r*_elicit`, `N_ident`, `N_gain` at `r ≥ 1` | **S9** | S9-as-written **cut**; replaced by *"beyond `r*_elicit` hops"* with the number printed | S9 survives with an explicit radius, never the word "neighbourhood" alone |
| `L_set(0, ρ=1)`, `E_meek` rate | **S3, S4** | X2's `D_am`-denominated rates will differ from 67 % / 14.1 %, which are generic-arm, ρ=1, `p ≤ 8` figures. Reported as a **conditioning disclosure**, and every X2 figure prints its conditioning on the same line | — |
| `E_bias / E_set` ratio | **S13** | if < 1, *"O\* changed" / "the estimate moved" / "silent bias"* stop being one measurement and the campaign's triple-counting is corrected | if = 1, X2 states plainly that they are **one** measurement and must not be cited as three |
| distance table with `disconnected` broken out | **S5, S8** | corrects E1′'s published `large` "dist ≥ 2 (n=331)" row, 17 % of which is infinite distance (§1.2) | — |
| ARM A `true_in_D` split | **S8, S9** | a **correct** far statement that moves `O*` refutes the pruning and the elicitation advice without refuting S5 | — |

### 8b. Vault items this run bears on — flagged, not acted on

- 🔴 `wiki/activities/active-claims.md:570` records b-LOAD's graph-distance falsifier as *"currently
  NOT firing"* on the strength of `0/791` — a number retired 2026-08-14 by the campaign's own order
  (P6), on a day when the campaign's own data produced 1/529 at min-distance 1. **The b-LOAD CIKM
  camera-ready is 2026-08-23.** Under Operating Rule #1 this register entry must be corrected before
  any locality sentence is pasted into that camera-ready. **Zé's decision and Zé's wording; no agent
  edits a live claim.**
- `wiki/projects/rho-breakdown-adjustment.md:204-205, :224` still carries `0/791` twice, once with
  the word *"provably"*, which no measurement of any size can support.
- E1′ `AUDIT.md` A10 remains open: there is **no** ρ-breakdown entry in `active-claims.md`, and
  `grep -rl 'relevance:.*rho-breakdown' wiki/literature/` returns zero notes, so `check_claims.py`
  is structurally blind to Taeb–Guo–Henckel 2511.10625, gadjid 2402.08616, Fang 2207.05067 and
  CausalGuard 2605.21928. X2 produces evidence for a claim that is not registered.

---

## 9. Mandatory citations

- **Meek (1995)** and **Perković et al. (UAI'17) Algorithm 1** — the closure and the background-
  knowledge consistency check being tested. The propagation X2 hunts for is *their* mechanism.
- **Henckel, Perković & Maathuis arXiv:1907.02435** — `O*` and its variance optimality. `O*` is
  their estimand; every locality claim here is a claim about *their* functional, not a new one.
- **Guo & Perković** — one adjustment set per MPDAG member is not a contribution.
- **Taeb, Guo & Henckel arXiv:2511.10625** — nearest competitor, structural error against a
  population proxy. 🔴 **not in the corpus**; ingesting it is a prerequisite for writing, not for
  running.
- **Grinsztajn-style scope honesty** — everything is linear-Gaussian iSCM, ER, `p ≤ 25`, and the
  claim is licensed to that box and says so.

---

## 10. Deliverables

```
~/mugango/output/2026-08-19_x1-x2-locality-spurious/x2/
  PREREG.md      <- this file, frozen
  AUDIT.md       <- adversarial read of this design, written BEFORE implementation
  RESULTS.md     <- verdicts on S5, S8, S9 separately; every fraction with four denominators
  code/          <- copied from ~/latent-causal/x2-locality/code/
  results/       <- arm{A,B,C,D}_*.json + analyses
  logs/
```

Work dir on betelgeuse: `~/latent-causal/x2-locality/{code,results,logs}`. Nothing under
`~/latent-causal/e1prime-se/` or `~/latent-causal/rho-breakdown-knowledge/` is modified. No git
commit is made by any agent.
