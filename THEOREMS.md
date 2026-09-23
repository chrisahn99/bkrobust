# THEOREMS — formal statements and status

Session 2 & Session 7, Axis B. Each item states its hypotheses, and is either **proved** (a
proof a referee can check) or **verified** (with the verification scope stated
explicitly). Nothing is labelled proved on the strength of evidence alone.

---

## 0. Setting and notation

Let `Ĉ` be a CPDAG and `[Ĉ]` its Markov equivalence class — a set of DAGs.

For a set `K` of orientations of edges adjacent in `Ĉ`, write `Meek(Ĉ, K)` for
the result of imposing `K` and closing under Meek's rules R1–R4, or FAIL if `K`
is inconsistent with `Ĉ`.

**Meek's Theorem** (Meek 1995, used throughout as an established input). If `K`
is consistent with `Ĉ` then `Meek(Ĉ, K)` is the *maximally oriented* PDAG
representing exactly

```
[Meek(Ĉ,K)] = { D ∈ [Ĉ] : D ⊨ K } .
```

*Maximally oriented* means: for every edge left undirected, both orientations
occur among the represented DAGs. (Verified computationally for every element of
every space examined this session — see §1.)

**The space.**

```
𝔊_Ĉ  =  { Meek(Ĉ, K) : K consistent with Ĉ }
```

— the set of **knowledge states**. For `G ∈ 𝔊_Ĉ` write `K_G = dir(G) \ dir(Ĉ)`
for the orientations `G` adds to the CPDAG, so `[G] = { D ∈ [Ĉ] : D ⊨ K_G }`.

**Order.** `G ⪯ H` iff `[G] ⊆ [H]`. More knowledge means fewer models means
*lower*; `Ĉ` is the maximum; DAGs are minimal. Moving **up** is **retraction**.

**Distance.** `d(·,·)` is hop count on the undirected covering graph, whose edges
are exactly the covering pairs of `⪯`. Covering is computed from model inclusion
by enumerating represented DAGs — never from a conjectured characterisation.

**Radius convention** (unchanged): `r = min{ d(G₀,G) : the property fails at G }`,
so shells `0 … r−1` are certified clean and the practitioner's safe-move count is
`r − 1`.

---

## 1. Task 0 — the space membership and fixpoint definition  *(settled)*

**Statement.** `𝔊_Ĉ` as defined above is *strictly larger* than the set produced
by the inherited `enumerate_space`, which filters candidates through a chordality
test on the **undirected subgraph**.

**Status: settled, with a minimal witness.**

That chordality condition characterises *CPDAGs* — essential graphs are chain
graphs with chordal chain components. It is **not** a property of an MPDAG built
by adding background knowledge: knowledge can place a directed edge *inside* an
otherwise-undirected component, and the chord that would make the component
chordal is then present but **directed**, so the undirected subgraph reads as
chordless.

*Minimal witness* (`results/axisb2/task0_space_membership.json`):

```
Ĉ     : V0-V1 V0-V2 V0-V3 V1-V2 V1-V3        (K₄ minus V2–V3)
K     : { V0→V1 }                             — a single assertion
state : V0->V1 V0-V2 V0-V3 V1-V2 V1-V3
        Meek-closed ✓   5 DAG extensions ✓   maximally oriented ✓
        old is_valid_mpdag = False  →  excluded
```

Every excluded graph found is reachable, Meek-closed, non-empty, and maximally
oriented — so it is a legitimate MPDAG and the *test* was wrong, not the notion.

**Corrected membership predicate** (fixpoint form):

```
G ∈ 𝔊_Ĉ   ⟺   Meek(Ĉ, dir(G) \ dir(Ĉ)) = G
```

Verified equivalent to reachability on **1,588 states** over all CPDAGs on 3 and
4 nodes, with zero failures.

**Scale.** The previous session reported 0.09%, measured per Meek closure. Per
element of the space the loss is larger and **grows with density** (n = 5):

| undirected edges `k` | elements missing | CPDAGs affected |
|---|---|---|
| ≤ 4 | 0.00% | 0 / 12 |
| 5 | 3.85% | **12 / 12** |
| 6 | 2.14% | 9 / 12 |
| 7 | 8.25% | 12 / 12 |
| 8 | 10.19% | 12 / 12 |

The previous session's sweeps ran to `k = 6`, so its numbers were computed on a
space missing elements for most dense CPDAGs.

**Re-derived on the corrected space** (session 2): all 8,782 CPDAGs on 5 nodes,
1,977,820 radius comparisons, **0 counterexamples** to
Conjecture 2. The correction added 720 space elements and
13,860 comparisons involved a state the old space did not contain.

---

## 2. Theorem 1 — failure is upward-closed  *(inherited, proved)*

**Statement.** Fix `Z`, `X`, `Y`. If `Z` fails in `G` and `[G] ⊆ [G']`, then `Z`
fails in `G'`.

**Proof.** Validity is a for-all over represented DAGs, so `Z` failing in `G`
means some `D ∈ [G]` has `Z` invalid. Since `[G] ⊆ [G']`, that same `D` lies in
`[G']` and witnesses failure there. ∎

**Caveat** (convention, not mathematics). This implementation defines validity to
be *false* when `[G]` is empty. An empty-extension graph would then "fail"
vacuously while `[] ⊆ [G']` for every `G'`, breaking the implication. It does not
bite because no element of `𝔊_Ĉ` has an empty extension set — checked, not
assumed.

---

## 3. Theorem 2 — `𝔊_Ĉ` is a join-semilattice, and the join is explicit  *(proved)*

**Statement.** For `G, H ∈ 𝔊_Ĉ`,

```
G ∨ H  =  Meek(Ĉ, K_G ∩ K_H)
```

is the least upper bound of `G` and `H` in `⪯`. Hence `𝔊_Ĉ` is a join-semilattice
with maximum `Ĉ`.

**Proof.** Write `J = Meek(Ĉ, K_G ∩ K_H)`.

*Upper bound.* Let `D ∈ [G]`. Then `D ⊨ K_G`, hence `D ⊨ K_G ∩ K_H`, hence
`D ∈ [J]` by Meek's Theorem. Symmetrically `[H] ⊆ [J]`. So `[G] ∪ [H] ⊆ [J]`.

*Least.* Let `W ∈ 𝔊_Ĉ` with `[G] ∪ [H] ⊆ [W]`. Take any `e = (a→b) ∈ K_W`. Since
`[G] ⊆ [W]`, every `D ∈ [G]` satisfies `e`; that is, the edge `a—b` is oriented
identically across all of `[G]`. Because `G` is maximally oriented (Meek's
Theorem), `G` therefore contains `a→b` as a directed edge. Since `e ∉ dir(Ĉ)` by
definition of `K_W`, we get `e ∈ K_G`, and symmetrically `e ∈ K_H`. Hence
`K_W ⊆ K_G ∩ K_H`, so every `D ⊨ K_G ∩ K_H` satisfies `K_W`, i.e. `[J] ⊆ [W]`. ∎

**Verified independently.** A least upper bound exists, and the construction
above *is* that least upper bound, in **100%** of pairs at both scopes:

| scope | pairs | join exists | construction is the join |
|---|---|---|---|
| all CPDAGs n = 3, 4 (exhaustive) | 40,240 | 100% | 100% |
| n = 5, sampled across k = 1…7 | 220,476 | 100% | 100% |
| **total** | **260,716** | **100%** | **100%** |

(`results/axisb2/lattice/joins_lemmaL_n34.json`, `.../lattice_n5.json`.)

---

## 4. Anti-Exchange Property  *(PROVED)*

**Statement.** Let `cl(S)` denote the logical closure operator `Meek(Ĉ, S)` over
admissible orientation sets. For any closed set `S` (representing an MPDAG `G ∈ 𝔊_Ĉ`)
and distinct orientations `x, y ∉ S`:

```
y ∈ cl(S ∪ {x})   ⟹   x ∉ cl(S ∪ {y}) .
```

**Semantic translation.** Suppose the implication fails, so both `y ∈ cl(S ∪ {x})`
and `x ∈ cl(S ∪ {y})` hold. By the monotonicity and idempotence of the closure
operator `cl`, this mutual deduction requires:

```
cl(S ∪ {x}) = cl(S ∪ {y}) = cl(S ∪ {x, y}) .
```

Semantically, this equality implies that across the entire set of valid DAG
extensions `[G]`, the constraints `x` and `y` are strictly equivalent:

```
∀ D ∈ [G],   D ⊨ x  ⟺  D ⊨ y .
```

Hence, the anti-exchange property is equivalent to the fundamental structural assertion:
**no two distinct uncompelled edge orientations are perfectly correlated across
the DAG extensions of an MPDAG.**

---

### Proof of the Anti-Exchange Property

We partition all candidate pairs into two exhaustive structural cases:

#### Case A: `x` and `y` orient the same edge in opposite directions.
Let `x = (u → v)` and `y = (v → u)` be opposing orientations of the same undirected
edge `u — v` in the MPDAG represented by `S`.

1. Since `S` is closed, `u — v` is undirected in `G`.
2. By Meek's Theorem, `G` is maximally oriented, meaning both orientations `u → v`
   and `v → u` occur among the DAG extensions in `[G]`.
3. Pick any DAG extension `D₁ ∈ [cl(S ∪ {x})]`. By construction, `D₁ ⊨ (u → v)`.
4. In any valid directed acyclic graph, an edge cannot be oriented in both directions
   simultaneously without forming a directed 2-cycle; hence `D₁ ⊭ (v → u)`.
5. Therefore, `y ∉ cl(S ∪ {x})`, and the antecedent `y ∈ cl(S ∪ {x})` is false.
   The implication holds vacuously. ∎

#### Case B: `x` and `y` orient distinct undirected edges.
Let `x = (u₁ → v₁)` and `y = (u₂ → v₂)` target two distinct undirected edges
`u₁ — v₁` and `u₂ — v₂` in `G`.

**Assumption for contradiction:** Suppose `y ∈ cl(S ∪ {x})` and `x ∈ cl(S ∪ {y})`.

1. **Absence of mixed configurations:**
   By the semantic translation above, for every DAG `D ∈ [G]`, `D ⊨ x ⟺ D ⊨ y`.
   Because `G` is maximally oriented and the edges are undirected in `G`, there
   exists at least one DAG `D_A ∈ [G]` satisfying `x`, and at least one DAG
   `D_B ∈ [G]` satisfying the opposing orientation `¬x = (v₁ → u₁)`.
   By the equivalence condition:
   - `D_A ⊨ (x ∧ y)`
   - `D_B ⊨ (¬x ∧ ¬y)`
   
   Crucially, the mixed configurations `(x ∧ ¬y)` and `(¬x ∧ y)` cannot appear in
   any extension in `[G]`:
   ```
   { D ∈ [G] : D ⊨ (¬x ∧ y) }  =  ∅
   { D ∈ [G] : D ⊨ (x ∧ ¬y) }  =  ∅
   ```

2. **Chickering transformational path:**
   By Chickering's Transformational Theorem (Chickering 1995, 2002), any two DAGs
   belonging to the same Markov equivalence class can be connected by a finite
   sequence of covered edge reversals:
   ```
   D_A = D₀  —(r₁)→  D₁  —(r₂)→  D₂  —(r₃)→  …  —(r_m)→  D_m = D_B
   ```
   where:
   - An edge `a → b` is *covered* in `D_k` if `Pa(b) = Pa(a) ∪ {a}`. Reversing a
     covered edge produces a DAG in the identical Markov equivalence class `[Ĉ]`.
   - Both `D_A` and `D_B` are extensions of the MPDAG `G`, so they share all
     compelled directed edges `dir(G)`.
   - Reversing a compelled directed edge creates a non-equivalent v-structure or a
     directed cycle, which leaves the equivalence class `[Ĉ]`. Therefore, the
     sequence operates exclusively on uncompelled edges.
   - Consequently, every intermediate DAG in the sequence preserves all compelled
     edges of `G`, which means:
     ```
     ∀ k ∈ {0, …, m},   D_k ∈ [G] .
     ```

3. **Step analysis of the reversals:**
   In the sequence `D₀, D₁, …, D_m`:
   - Edge `u₁ — v₁` begins as `x` at `D₀` and ends as `¬x` at `D_m`.
   - Edge `u₂ — v₂` begins as `y` at `D₀` and ends as `¬y` at `D_m`.

   Each individual transition `D_{k-1} → D_k` reverses exactly one covered edge.
   Let `i ∈ {1, …, m}` be the step at which edge `u₁ — v₁` is reversed from `x` to `¬x`.
   Let `j ∈ {1, …, m}` be the step at which edge `u₂ — v₂` is reversed from `y` to `¬y`.

   Because `x` and `y` correspond to distinct edges, and each step reverses exactly
   one edge, `i ≠ j`. There are two exhaustive possibilities:
   - **Subcase B.1 (`i < j`):** At step `i`, edge `u₁ — v₁` is flipped to `¬x`.
     Because `i < j`, step `j` has not yet been reached, so edge `u₂ — v₂`
     still maintains its initial orientation `y`. Thus:
     ```
     D_i ⊨ (¬x ∧ y) .
     ```
   - **Subcase B.2 (`j < i`):** At step `j`, edge `u₂ — v₂` is flipped to `¬y`.
     Because `j < i`, step `i` has not yet been reached, so edge `u₁ — v₁`
     still maintains its initial orientation `x`. Thus:
     ```
     D_j ⊨ (x ∧ ¬y) .
     ```

4. **Contradiction:**
   In both subcases, the transformational sequence produces an intermediate DAG
   `D* ∈ {D_i, D_j} ⊂ [G]` that satisfies a mixed configuration: either
   `D* ⊨ (¬x ∧ y)` or `D* ⊨ (x ∧ ¬y)`.
   
   This directly contradicts Step 1, which established that no extension in `[G]`
   can exhibit a mixed configuration.

Therefore, the assumption of mutual implication is false, proving unconditionally:
```
x ∉ cl(S ∪ {y}) .
```
Both Case A and Case B hold, completing the proof of the Anti-Exchange Property. ∎

---

## 4b. Lemma R — every cover adds exactly one orientation  *(PROVED)*

**Statement.** If `X ⋖ Y` in `𝔊_Ĉ` then `|K_X \ K_Y| = 1`.

**Proof.** 
By Edelman & Jamison (1985), a closure operator satisfies the anti-exchange property
if and only if its closed sets form a **convex geometry**.
In any convex geometry, the lattice of closed sets is graded, and every covering
relation `K_X ⋖ K_Y` in the lattice of closed sets differs by exactly one extreme
element of the closure. In `𝔊_Ĉ`, this translates directly to:
```
|K_X \ K_Y| = 1 .
```
Since the Anti-Exchange Property is proved unconditionally (§4), Lemma R is proved. ∎

**Verified directly:** **5,348 covering pairs** (2,658 exhaustive at n ≤ 4,
2,690 sampled at n = 5), every one of them adding exactly one orientation,
**0 violations**. Equivalently, `rank(G) = |dir(G)| − |dir(Ĉ)|` is a valid rank
function on the poset `(𝔊_Ĉ, ⪯)`.

---

## 4c. Property S — upper semimodularity  *(PROVED)*

**Statement.** For all `X, Y, Z ∈ 𝔊_Ĉ` with `X ⋖ Y`:

```
X ∨ Z = Y ∨ Z      or      X ∨ Z ⋖ Y ∨ Z .
```

**Proof.** Work on the orientation side, where an element `W` is identified with
its closed knowledge set `K_W` and, by Theorem 2, `W ∨ Z` is identified with
`K_W ∩ K_Z`. Two facts are used: the intersection of closed sets is
closed (`cl(A ∩ B) ⊆ cl(A) ∩ cl(B) = A ∩ B`), and by **Lemma R** (§4b), `X ⋖ Y`
gives `K_X = K_Y ⊔ {e}` for a single orientation `e`.

*Case 1: `e ∉ K_Z`.* Then `K_X ∩ K_Z = K_Y ∩ K_Z`, so `X ∨ Z = Y ∨ Z`.

*Case 2: `e ∈ K_Z`.* Then `K_X ∩ K_Z = (K_Y ∩ K_Z) ⊔ {e}`. Both are closed sets,
so they are elements of `𝔊_Ĉ`, and their orientation sets differ by the single
element `e`. No closed set lies strictly between two sets differing by one element,
so `X ∨ Z ⋖ Y ∨ Z`. ∎

**Verified:** **1,386,626 triples** (191,484 exhaustive at n = 3, 4; 1,195,142 sampled
at n = 5) with zero gaps of height ≥ 2.

---

## 5. Lemma L — distance non-expansion under join  *(PROVED)*

**Statement.** For all `G₀, G ∈ 𝔊_Ĉ`: `d(G₀, G ∨ G₀) ≤ d(G₀, G)`.

**Proof from Property S.** Let `G₀ = P₀, P₁, …, P_d = G` be a shortest path in the
covering graph, so each consecutive pair `{P_i, P_{i+1}}` is a covering pair in one
direction or the other. Put `Q_i = P_i ∨ G₀`. Then `Q₀ = G₀ ∨ G₀ = G₀` and
`Q_d = G ∨ G₀`.

For each `i`, `{P_i, P_{i+1}}` is a covering pair, so Property S applied with
`Z = G₀` gives `Q_i = Q_{i+1}` or `Q_i ⋖ Q_{i+1}` or `Q_{i+1} ⋖ Q_i`. In every
case `Q_i` and `Q_{i+1}` are equal or adjacent in the covering graph. Deleting
repetitions leaves a walk from `G₀` to `G ∨ G₀` of length at most `d`. Hence
`d(G₀, G ∨ G₀) ≤ d = d(G₀, G)`. ∎

**Verified independently:** **521,432** `(G₀, G)` pairs, **0 violations** (80,480
exhaustive at n = 3, 4; 440,952 at n = 5 across densities).

---

## 6. Conjecture 2 — Monotone Reachability of Nearest Failure  *(PROVED THEOREM)*

**Statement.** The nearest failure to `G₀` is reachable by a purely monotone upward
path in the covering graph of `𝔊_Ĉ`. Consequently, the breakdown radius `r_val`
equals the minimum number of background knowledge retractions from `G₀` required
to induce failure.

**Theorem.**
```
Anti-Exchange (PROVED) ⟹ Lemma R (PROVED) ⟹ Property S (PROVED) ⟹ Lemma L (PROVED) ⟹ Conjecture 2 (PROVED).
```

**Proof.** Let `G` be a nearest failure to `G₀` in the entire knowledge space `𝔊_Ĉ`,
with graph distance `d(G₀, G) = r_val`. Define the join state `J = G ∨ G₀`.

1. By definition of the join, `G ⪯ J` (i.e., `[G] ⊆ [J]`). By **Theorem 1**
   (upward-closure of adjustment set failure), since `Z` fails in `G`, `Z`
   definitively fails in `J`.
2. By **Lemma L** (§5), `d(G₀, J) ≤ d(G₀, G) = r_val`.
3. Because `r_val` is the minimal distance from `G₀` to *any* failing state in
   `𝔊_Ĉ`, and `J` is a failing state, we must have `d(G₀, J) ≥ r_val`.
   Combining the inequalities yields:
   ```
   d(G₀, J) = r_val .
   ```
4. By Property S, `𝔊_Ĉ` is upper semimodular, hence satisfies the Jordan–Dedekind
   chain condition, ensuring that the poset `(𝔊_Ĉ, ⪯)` is strictly graded with
   rank function `ρ(W) = |K_W|`.
5. For comparable states `G₀ ⪯ J`, every path in the covering graph changes rank by
   at most 1 per step, so any path has length at least `ρ(G₀) − ρ(J)`.
   Furthermore, there exists a saturated chain from `G₀` up to `J` of length
   exactly `ρ(G₀) − ρ(J)`.
   Therefore, `d(G₀, J) = ρ(G₀) − ρ(J)` and this shortest path is realised by a
   **monotone upward path** formed entirely of single-edge retractions.

Hence, the optimal failure witness is reachable purely by monotonic retractions,
and the retraction-only radius equals `r_val`. ∎

---

## 7. Operational guarantee for `radius_local_up`

`radius_local_up` performs upward BFS via knowledge retractions.
With the formal closure of the Anti-Exchange Property, Conjecture 2 is now an
unconditional theorem.

| | Status |
|---|---|
| Exactness of `radius_local_up` | **Mathematically guaranteed exact**; BFS upward search cannot overestimate or underestimate `r_val`. |
| Poset structure | Proved to be a **graded convex geometry** and an **upper semimodular join-semilattice**. |
| Search complexity | Retraction-only search provably suffices; no downward exploration is required. |

---

## 8. Summary Table of Formal Statuses

| # | Statement | Status | Notes |
|---|---|---|---|
| T1 | Failure is upward-closed | **PROVED** | Universal validity over extensions |
| T2 | `G ∨ H = Meek(Ĉ, K_G ∩ K_H)`; `𝔊_Ĉ` is a join-semilattice | **PROVED** | Verified on 260,716 pairs |
| — | Space membership = reachability (fixpoint predicate) | **SETTLED** | Verified on 1,588 states |
| AE-A | Anti-Exchange, opposing orientations of same edge | **PROVED** | Follows from maximal orientation |
| AE-B | Anti-Exchange, distinct edges | **PROVED** | §4's argument has a false premise; repaired in §24.2 by counting Chickering's reversals |
| AE | Anti-Exchange Property (full) | **PROVED** | Case A from maximal orientation; Case B from §24.2. Unconditional |
| R | Every cover adds exactly one orientation (`|K_X \ K_Y| = 1`) | **PROVED** | From AE (§24); verified directly on 5,348 covering pairs |
| S | Upper semimodularity (`X ⋖ Y ⟹ X∨Z = Y∨Z` or `X∨Z ⋖ Y∨Z`) | **PROVED** | Via Lemma R; verified on 1.38M triples |
| L | Distance non-expansion under join (`d(G₀, G∨G₀) ≤ d(G₀, G)`) | **PROVED** | Via Property S; verified on 521k pairs |
| C2 | Monotone reachability of failure (Conjecture 2) | **PROVED THEOREM** | From Lemma L, T1, gradedness; 1,977,820 comparisons, 0 counterexamples |
| Fast-O | Order identification via directed edge containment (`Lemma O`) | **PROVED** | Reduces cover check from exponential to set ops |
| GAC | GAC radius bounds back-door radius (`r_val(GAC) ≥ r_val(BD)`) | **PROVED** | Forb set inclusion |
| Opt | Criteria equivalence on optimal adjustment set `O*` | **PROVED** | Verified across 5,304 instances |

---

# Session 3 addendum — correctness of the declarative encodings

Everything below concerns the CP-SAT encodings in `src/bkrobust/sat/`. The
question these sections answer is not "is the radius right" but "does the
encoding define the object we think it defines", which is prior to it. A wrong
encoding does not crash; it returns plausible numbers.

## 9. The closure encoding defines the corrected space

**Statement.** For a CPDAG `Ĉ`, the satisfying assignments of the constraint
system built by `sat.closure.add_meek_closure` are exactly the orientation sets
`{K_G : G ∈ 𝔊_Ĉ}` of the corrected space of §1.

**Method.** Each Meek rule is encoded as a **forbidden firing configuration**
rather than as an implication: for every tuple of vertices matching a rule's
premises, the conjunction of (premises ∧ ¬consequent) is ruled out by one clause.
Added to this are the no-new-v-structure clauses and acyclicity via position
integer variables.

**Why firing configurations and not implications.** This was wrong once. R1 and R2
survive a weaker encoding that omits the undirectedness premise, because the
alternatives to the consequent are either a new v-structure or a cycle, both
independently forbidden. R3 and R4 do not: dropping their undirectedness premises
makes the system strictly too strong and it loses legitimate elements. The first
version of the encoding did exactly that and reproduced only **66 of 133** CPDAG
spaces.

**Verified.** With the premises restored, the solution set equals
`enumerate_space_correct` on **133 CPDAGs / 1,588 elements**, 0 mismatches
(`tests/sat/test_encodings.py::test_closure_encoding_is_the_corrected_space`).

## 10. The failure predicate agrees with the validity oracle

**Statement.** `sat.failure.add_failure` is satisfiable for a witness DAG `D`
extending `G` exactly when `Z` is an invalid adjustment set in `D`.

**Method.** The two back-door failure modes are encoded directly: (A) some
`z ∈ Z` is a descendant of `X` in `D`, and (B) a d-connecting path from `X` to
`Y` in `D_X̄` given `Z`. Path (B) uses an **explicit bounded simple path** with
position variables, not a Bayes-ball fixpoint; the fixpoint formulation is unsafe
in both directions here. Two fixpoints do remain (`rx`, `ancz`), and they are
safe because they range over an acyclic graph.

**The collider clause was wrong once**, and the error made radii **too small** —
the dangerous direction. A non-collider clause requiring `into_l ∨ into_r` let a
*chain* through a member of `Z` be accepted as a collider, so paths that are
blocked were reported open, so failures were reported that were not failures. A
collider needs **both** arrows pointing in, which is two clauses, not one.

**Verified.** Against the reference oracle on **26,304 cases**, 0 disagreements
after the fix.

## 11. E1, E2 and E3 — what each one assumes

Let `r` be the true BFS radius. All three encodings are **upper bounds** on `r`;
none can return a radius that is too small. That one-sidedness is inherited from
§7 and is the dangerous direction for reporting, so it is stated on every result
object via an `assumes` field.

| | searches | assumes | relation |
|---|---|---|---|
| **E1** | retractions of `K_{G₀}` | Conjecture 2 (Now Proved) | `r ≤ r_E1` (Now `r = r_E1`) |
| **E2** | every closed state, distance via the join | Theorem 2 + Lemma R + up-then-down normalisation | `r ≤ r_E2 ≤ r_E1` |
| **E3** | walks of one-orientation covering steps | **nothing** | `r ≤ r_E3` |

**E2's join needs no closure.** Theorem 2 identifies `G₀ ∨ G` with `K_{G₀} ∩ K_G`,
and the intersection of closed sets is closed, so `K_J = K_{G₀} ∩ K_G` outright.
With `K_{G₀}` constant the objective collapses to the symmetric difference
`|K_{G₀} Δ K_G|` — one linear constraint over a single copy of the variables,
rather than the extra copy the plan budgeted. The assumptions are unchanged.

## 12. E3's soundness, and the direction it runs in

**Statement.** If E3 returns a walk of length `k` ending at a failing state, then
`r ≤ k`, unconditionally.

**Proof.** Consecutive states in an E3 walk are closed and differ by exactly one
orientation. Two sets differing by a single element admit nothing strictly
between them, so the pair is a covering pair by the definition of covering — with
no appeal to Lemma R, to anti-exchange, or to gradedness. Hence the walk is a
walk in the covering graph, and BFS distance is at most its length. ∎

**Consequences:**
Now that Anti-Exchange and Lemma R are proven, every cover in `𝔊_Ĉ` is certified
to differ by exactly one orientation. Thus, E3 searches the full covering space
of paths without loss of generality.

**Certificates, not solver claims.** Every E3 walk reported in this project was
replayed through the ordinary graph code by `sat.verify.verify_walk`, which
re-checks, without consulting the encoding, that each state is a knowledge state
of the CPDAG, that each step changes exactly one orientation, and that the final
state genuinely fails.

---

# Session 4 addendum — the order identification, and the MPDAG criterion

## 13. Lemma O — model inclusion is reverse containment of orientations

**Statement.** For elements `G`, `H` of the corrected space `𝔊_Ĉ`:

```
[G] ⊆ [H]      ⟺      dir(H) ⊆ dir(G) ,
```

and the two are strict together.

**Proof.**

*(⇐)* Suppose `dir(H) ⊆ dir(G)`. Every `D ∈ [G]` is a DAG consistent with the
orientations `dir(G)`, hence with the smaller set `dir(H)`, and it shares the
skeleton and v-structures of `Ĉ`. So `D ∈ [H]`, giving `[G] ⊆ [H]`.

*(⇒)* Suppose `[G] ⊆ [H]` and let `a→b ∈ dir(H)`. Every `D ∈ [H]` orients that
edge `a→b`, and every `D ∈ [G]` lies in `[H]`, so every `D ∈ [G]` orients it
`a→b`. Elements of the space are **maximally oriented** — an edge oriented
identically across the whole of `[G]` is directed in `G`, which is the defining
property of §1 and is what `space_fixed.is_maximally_oriented` checks. Hence
`a→b ∈ dir(G)`. ∎

*Strictness.* An element is determined by its extension set and by its
orientation set alike, so equality on one side forces equality on the other, and
the strict versions correspond.

**Verified.** All ordered pairs of elements, over every CPDAG on 3 and 4 nodes
with at least one undirected edge: **80,480 pairs, 0 disagreements**
(`results/axisb4/order_identification.json`).

**Consequence.** `local_up`'s cover generation decided minimality by comparing
*extension sets*, enumerating `[H]` once per candidate. Profiling put **70.4%**
of `local_up`'s total time in that one test. Lemma O makes it a set comparison on
directed edges, with no enumeration at all. `search/exact_fast.py` implements this;
it returns identical cover sets (1,588 compared, identical including order) and
identical radii (1,044 compared, agreeing with both frozen `local_up` and brute-force BFS).

---

# Session 5 addendum — the two adjustment criteria, and where they coincide

## 14. The GAC accepts a superset, so radii can only grow

**Statement.** For any DAG `D` and any `Z`: back-door-valid ⟹ GAC-valid. Hence
for any MPDAG `G` the GAC failure set is a subset of the back-door failure set,
and for every instance

```
r_val(GAC)  ≥  r_val(back-door) .
```

**Proof.** Back-door requires `Z ∩ de(X) = ∅` and that `Z` block every back-door
path. GAC requires `Z ∩ forb(X,Y) = ∅` and that `Z` block every proper
non-causal path. Since `forb(X,Y) ⊆ de(X) ∪ {X,Y}`, the first condition is
weaker; and given `Z ∩ de(X) = ∅`, blocking every back-door path and blocking
every proper non-causal path coincide. So a back-door-valid set is GAC-valid.
Failure is a for-all over extensions in both cases, so the implication lifts to
MPDAGs. Fewer failures means the nearest one is no nearer. ∎

**Verified.** 4,711,024 DAG-level cases over all 29,824 labelled DAGs on 4 and 5
nodes: **0 implication failures**, and all 227,224 disagreements classified as
"`Z` contains a descendant of `X` off the causal route to `Y`" — 0 unexplainable.
The radius-level invariant was checked on 22,230 instances at n = 3, 4 with
**0 violations**, and is asserted per instance throughout the session-5 census.

## 15. The two criteria coincide on the optimal adjustment set

**Statement.** For the Henckel–Perković–Maathuis optimal set
`O = pa(cn(X,Y)) \ (cn(X,Y) ∪ {X})`, back-door validity and GAC validity are the
same predicate, in every graph.

**Proof.** `O` is disjoint from `de(X)`: every element is a parent of a node in
`cn(X,Y)` and is itself excluded from `cn(X,Y) ∪ {X}`, so no element lies on or
below a causal route from `X`. For any `Z` with `Z ∩ de(X) = ∅`, the two
criteria's first conditions are both satisfied, and their second conditions
coincide (§14). Hence they agree on `O`. ∎

**Verified twice, by different routes.** The criterion sweep found 0
disagreements in 8,154 evaluations (3,360 both-valid, 4,794 both-invalid, so not
vacuous). Independently, radii were recomputed by full BFS over the corrected
space under each definition: **5,304 instances, radii identical in 100.00%**,
and of the 2,790 with back-door radius 1, **0** have a strictly larger GAC
radius.

---

# Session 8 addendum — from a validity radius to a bias radius

`r_val` answers *where does the committed adjustment set stop being valid*. It
does not answer *how wrong the answer becomes*, and the paper's own App. C table
shows the two questions come apart: past `r_val` the bias is a ramp, not a step.

This session defines a bias functional on the same space, shows it is exactly
zero on the ball `r_val` certifies and monotone under retraction, and derives a
radius `r_ε` at which it first exceeds a threshold. Full statements and proofs
are in **`docs/R_EPSILON_THEORY.md`**; this section records the statements and
their status in the format of the rest of this file.

## 16. Setting — the bias functional

The analyst holds one observational law with population covariance `Σ` and
reports one number, `θ_Z = β(Z;Σ)`, the coefficient on `X` in the population
regression of `Y` on `{X} ∪ Z`. It does not depend on which DAG is the truth.
Each candidate truth `D ∈ [Ĉ]` assigns the effect its own value
`τ_D = β(pa_D(X);Σ)` — the IDA quantity, `0` by convention when `Y ∈ pa_D(X)`.
Then

```
T(G) = { τ_D : D ∈ [G] } ,      B(G) = max { |θ_Z − τ| : τ ∈ T(G) } .
```

This replaces the incumbent `oracle.bias_stats(...).mean_abs_bias`, which
redraws a SEM per extension and averages. That quantity is **not monotone** under
`⪯` — a mean over a larger set can be smaller — so `{mean > ε}` is not an up-set,
nothing in this file applies to it, and the upward search is not exact for it.

## 17. Theorem A — exact zero bias inside the certified shell  *(PROVED)*

**Statement.** If `d(G₀,G) ≤ r_val − 1` then `T(G) = {θ_Z}` and `B(G) = 0`.

**Proof.** No element inside the shell fails, so `Z` is valid in every `D ∈ [G]`;
a valid adjustment set identifies the total effect in a linear SCM, so
`τ_D = β(Z;Σ) = θ_Z` for every such `D`. ∎

The claim is *exactly zero at the population level*, for every law compatible
with `Ĉ`, not "no significant bias". Inside the shell the only error is sampling
error — which is the operational content: there is nothing to recompute there.

**Converse fails, in the safe direction.** At distance `r_val` some extension has
`Z` invalid, but `β(Z;Σ)` can still coincide with `τ_D` for the particular `Σ`.
So `r_val ≤ min{d : B > 0}`, with equality for generic `Σ`: **`r_val` is a
universal lower bound on every bias radius**, and understates robustness rather
than overstating it.

## 18. Theorem B — worst-case bias is monotone  *(PROVED)*

**Statement.** `G ⪯ H ⟹ T(G) ⊆ T(H) ⟹ B(G) ≤ B(H)`.

**Proof.** `τ_D` depends on `D` and `Σ` only, so `[G] ⊆ [H]` gives an inclusion
of value sets; a max over a subset is no larger. ∎

**Corollary.** `F_ε = {G : B(G) > ε}` is upward-closed, and `F_{ε'} ⊆ F_ε` for
`ε ≤ ε'`.

*Caveat, stated because the naive reading is tempting and false.* This is
monotonicity **along the order**, not along the distance. The bias of an
individual state is *not* monotone in `d(G₀,·)`; what is monotone is the maximum
over the ball. Counterexamples: `results/epsilon/counterexamples/`.

## 19. Theorem U — the retraction theorem was never about validity  *(PROVED)*

**Statement.** For **any** non-empty up-set `Φ ⊆ 𝔊_Ĉ`,
`min{d(G₀,G) : G ∈ Φ}` is attained at a state comparable to `G₀`, and equals the
least number of retractions from `G₀` that lands in `Φ`.

**Proof.** Verbatim §6, with `J = G ∨ G₀`: `J ∈ Φ` because `Φ` is an up-set,
`d(G₀,J) ≤ d(G₀,G)` by **Lemma L**, and minimality forces equality. ∎

Reading §6 back, the only property of *failure* it consumed is Theorem 1 —
upward-closure. So §6 is a corollary of this, and `r_ε` inherits it for free.

**Proposition R (retraction parameterisation).**
`↑G₀ = { Meek(Ĉ, K_{G₀} \ S) : S ⊆ K_{G₀} }`, and `d(G₀,G) = |K_{G₀}| − |K_G|`
for `G ∈ ↑G₀`. (Lemma O for the containment, Lemma R and Property S for the
distance.) Hence shell `d` is generated **directly** from the `d`-subsets, with
no traversal of any shell below it.

## 20. Theorem C — the staircase  *(PROVED)*

With `β↑(d) = max{B(G) : G ∈ ↑G₀, d(G₀,G) ≤ d}`:

1. `β↑` is non-decreasing and `β↑(d) = 0` for `d < r_val`;
2. `β↑(d)` is attained on the **sphere** at `d`, not merely on the ball;
3. `β↑(d) = max{ B(Meek(Ĉ, K_{G₀}\S)) : |S| = d }`;
4. `r_ε = min{d : β↑(d) > ε}` — the first crossings of `β` and `β↑` coincide;
5. `ε ↦ r_ε` is non-decreasing, and `r_val ≤ r_0 ≤ r_ε`.

Consequence: **one traversal answers every `ε`**, including an `ε` chosen
afterwards, and bisection on `d` is valid.

## 21. Theorem D — the practitioner's certificate  *(PROVED)*

**Statement.** If at most `k` of the orientations in `K_{G₀}` are false of the
truth, then the realised identification error obeys `|θ_Z − τ*| ≤ β↑(k)`. Hence

```
k < r_val              ⟹  θ_Z = τ*  exactly
k < r_ε                ⟹  |θ_Z − τ*| ≤ ε
r_{ε_x} ≤ k < r_{ε_y}  ⟹  ε_x  <  β↑(k)  ≤  ε_y
```

— the bias is pinned into a band: no larger than `ε_y`, and demonstrably able to
exceed `ε_x`.

**Proof.** The state `H = Meek(Ĉ, K_{G₀}\S)` over the false orientations `S`
contains the truth, so `τ* ∈ T(H)` and `|θ_Z − τ*| ≤ B(H) ≤ β↑(|S|)`. ∎

*Units.* `k` counts orientations of the **closure** `K_{G₀}`, the same unit as
`r_val`. It is not the number of *asserted* claims: Meek's rules cascade, so a
smaller number of wrong assertions can produce more wrong closure orientations.
(Getting this wrong makes a sound certificate look violated; it did, once, in
this session.) The claim-level analogue is `r_claim`, paper §6.3.

## 21b. Theorem D′ — the band is tight  *(PROVED)*

**Statement.** For every `k` there is a DAG `D ∈ [Ĉ]` that fits the observed law
exactly, contradicts at most `k` of the analyst's orientations, and realises
`|θ_Z − τ_D| = β↑(k)`.

**Proof.** Take `H ∈ ↑G₀` at distance `≤ k` attaining `β↑(k)` and `D ∈ [H]`
attaining `B(H)`. Then `D ∈ [Ĉ]` fits `Σ`; and `K_H ⊆ K_{G₀}` with
`|K_{G₀}| − |K_H| = d(G₀,H) ≤ k`, so `D` contradicts at most `k` of `K_{G₀}`. ∎

So `β↑(k)` is attained by a world the analyst cannot rule out on the evidence
they hold, and the band cannot be tightened without an assumption beyond the
CPDAG, the data and the error budget. Refining the `ε` grid refines the band;
nothing else does.

## 22. Theorem E — the ambiguity set is semi-local  *(PROVED)*

**Statement.** `{pa_D(X) : D ∈ [G]}` equals the family of `pa_G(X) ∪ S` over
`S ⊆ nb_G(X)` for which
`Meek(Ĉ, K_G ∪ {s→X : s∈S} ∪ {X→t : t ∈ nb_G(X)\S})` does not FAIL.

**Proof.** Both directions from Meek's theorem: a non-FAIL closure has a
non-empty extension set realising exactly that parent set, and any extension
witnesses the consistency of its own local orientation. ∎

Cost: `2^{|nb_G(X)|}` polynomial closures, no DAG enumeration — the same move
Lemma O made for cover generation. The enumeration route is retained as a
differential test only.

## 23. Summary of Session 8 statuses

| # | Statement | Status | Rests on |
|---|---|---|---|
| A | `B ≡ 0` on the ball of radius `r_val − 1` | **PROVED** | definition of `r_val` |
| B | `G ⪯ H ⟹ B(G) ≤ B(H)`; `F_ε` upward-closed | **PROVED** | the order alone |
| U | nearest state of **any** up-set is reachable by retraction | **PROVED** | Lemma L |
| R | `↑G₀ = {Meek(Ĉ,K₀\S)}`; `d = ρ(G₀) − ρ(G)` | **PROVED** | Lemma O, R, Property S |
| C | monotone staircase; sphere sufficiency; direct shells; one pass, all `ε` | **PROVED** | A, B, U, R |
| D | certificate, band `(ε_x, ε_y]` | **PROVED** | C, Meek |
| D′ | the band is attained, hence tight | **PROVED** | R; maximality of `B` |
| E | semi-local ambiguity set | **PROVED** | Meek |
| AE-B.1 | Chickering's sequence touches only *differing* edges | **PROVED** | Chickering strong form (counts the reversals) |
| AE-B.2 | the reversal path between two extensions stays inside `[G]` | **PROVED** | AE-B.1; closes the §4 gap, see §24.2 |

**`r_ε` adds no assumption to `r_val`.** Everything above rests on Lemma L,
Lemma R, Property S and Lemma O exactly as the existing radius does, and the
one-sided error direction is unchanged: an upward search can return a radius
that is too large — overstating robustness — never one that is too small.

## 24. §4's Case B — the gap in the supplied argument, and its repair

### 24.1 The gap

§4 upgrades Anti-Exchange Case B from *verified* to *proved* via Chickering's
covered-edge-reversal theorem. Step 2 asserts that the reversal sequence between
`D_A, D_B ∈ [G]` stays inside `[G]` because it "operates exclusively on
uncompelled edges". **That does not follow.** Chickering's theorem gives a path
inside the Markov equivalence class `[Ĉ]`, and the orientations of `K_G` are
*precisely* edges uncompelled in `[Ĉ]` — they are exactly the ones the theorem
is free to reverse. Nothing in that argument keeps an intermediate DAG inside
`[G]`.

The obvious substitute is unavailable too. The classical route is to reorient
freely inside a chain component using its chordality; but an MPDAG carrying
background knowledge can have a **non-chordal** undirected part, which is the
same phenomenon §1 documents. Re-confirmed computationally this session: for
`Ĉ = K₄ − {V2,V3}` and `K = {V0→V1}`, the state

```
V0→V1   V0−V2   V0−V3   V1−V2   V1−V3
```

is Meek-closed (R1–R4 all fail to fire), has 5 DAG extensions, is maximally
oriented, contains the semi-directed cycle `V0→V1−V2−V0`, and its undirected
part is a **chordless 4-cycle** `V0−V2−V1−V3−V0`. So the graph is a legitimate
MPDAG that is not a chain graph, and the chordality argument does not transfer.

### 24.2 The repair — counting the reversals

The conclusion is nonetheless true. What is needed is not the *existence* of a
Chickering path but its **length**, which the theorem also supplies.

**Chickering's Transformational Theorem, strong form** (Chickering 1995, Thm 2;
2002, Thm 4 — used here as an established input, exactly as §0 uses Meek's
theorem). Let `D` and `D'` be Markov equivalent and let `m` be the number of
edges oriented oppositely in the two. Then there is a sequence of **exactly `m`
distinct covered-edge reversals** in `D`, after each of which the graph remains
in the equivalence class of `D`, and after all of which `D = D'`.

**Lemma AE-B.1 (no agreeing edge is ever touched).** In such a sequence the `m`
reversed edges are exactly the `m` edges on which `D` and `D'` differ, each
reversed once.

*Proof.* Let `δ_k` be the number of edges on which `D_k` and `D'` differ, so
`δ_0 = m` and `δ_m = 0`. Each step reverses exactly one edge, so
`|δ_k − δ_{k−1}| = 1`. A walk of `m` unit steps from `m` to `0` must decrease at
every step; hence every step reverses an edge on which the current graph and
`D'` still differ. The reversals are of *distinct* edges and there are `m` of
them, so they are exactly the initially-differing edges, one each — and no edge
on which `D` and `D'` agree is ever reversed. ∎

**Lemma AE-B.2 (the path stays inside `[G]`).** For `D_A, D_B ∈ [G]` there is a
sequence of covered-edge reversals from `D_A` to `D_B` with **every**
intermediate DAG in `[G]`.

*Proof.* Let `A` be the set of orientations on which `D_A` and `D_B` agree.
Both extend `G`, so `dir(G) ⊆ A`. Apply the strong form to `D_A, D_B`: by Lemma
AE-B.1 no edge of `A` is ever reversed, so `dir(G) ⊆ A ⊆ dir(D_k)` for every
`k`; and each `D_k` lies in `[Ĉ]` by the theorem. A DAG of `[Ĉ]` containing
`dir(G)` is by definition an element of `[G]`. ∎

**Corollary (extension connectivity).** For every `G ∈ 𝔊_Ĉ`, the graph on `[G]`
whose edges are covered-edge reversals that remain in `[G]` is connected. ∎

### 24.3 Case B, proved

Assume for contradiction `y ∈ cl(S ∪ {x})` and `x ∈ cl(S ∪ {y})` for distinct
undirected edges `e₁ = u₁−v₁` and `e₂ = u₂−v₂` of `G`. As in §4, monotonicity
and idempotence of `cl` force
`cl(S∪{x}) = cl(S∪{x,y}) = cl(S∪{y})`, hence

```
∀ D ∈ [G] :   D ⊨ x  ⟺  D ⊨ y ,
```

so no extension exhibits a mixed configuration `(x ∧ ¬y)` or `(¬x ∧ y)`.

`G` is maximally oriented and `e₁` is undirected in `G`, so some `D_A ∈ [G]`
satisfies `x` and some `D_B ∈ [G]` satisfies `¬x`; by the equivalence
`D_A ⊨ x ∧ y` and `D_B ⊨ ¬x ∧ ¬y`. So `D_A` and `D_B` differ on **both** `e₁`
and `e₂`.

Apply Lemma AE-B.2. Every intermediate `D_k` lies in `[G]`, and by Lemma AE-B.1
`e₁` is reversed at exactly one step `i` and `e₂` at exactly one step `j`, with
`i ≠ j` because the reversed edges are distinct. If `i < j` then `D_i` has `e₁`
already flipped and `e₂` not yet, so `D_i ⊨ ¬x ∧ y`; if `j < i` then
`D_j ⊨ x ∧ ¬y`. Either way an element of `[G]` exhibits a mixed configuration,
contradicting the display above. ∎

### 24.4 Status after the repair

Case B is **PROVED**, and with it the full Anti-Exchange Property, Lemma R,
Property S, Lemma L, Conjecture 2 and the exactness of `radius_local_up` —
unconditionally, with no verified-but-unproved link left in the chain. The §8
table is updated accordingly. The corrected argument replaces one false premise
("the sequence touches only uncompelled edges") by a true one that the same
cited theorem already supplies ("the sequence has exactly `m` steps, so it
touches only *differing* edges"); the rest of §4's Case B is unchanged.

The same correction removes the one-sidedness caveat that Sessions 3-7 attached
to every reported radius. Results computed before this session were correct;
what changes is that the caveat "exact **iff** Conjecture 2 holds" becomes
"exact", for `r_val` and for `r_ε` alike.

**Independent verification.** The two load-bearing lemmas are computationally
checkable, and were checked **exhaustively on every CPDAG up to five nodes**,
with no per-CPDAG or per-state sampling anywhere in the run
(`results/epsilon/chickering/`, 241 s):

| check | objects checked | violations | scope |
|---|---|---|---|
| L1 — a covered-edge reversal preserves the equivalence class | 50,758 | **0** | exhaustive, n = 3, 4, 5 (all 29,281 DAGs at n = 5) |
| L1 — same, spot check | 7,034 | **0** | n = 6, 4,000 sampled DAGs |
| L2 — a difference-reducing sequence of exactly `m` reversals exists | 206,878 | **0** | exhaustive: every equivalence class and every ordered pair, n = 3, 4, 5 |
| L3 — every intermediate DAG stays inside `[G]` | 2,907,242 | **0** | exhaustive: all 8,782 CPDAGs and all 116,992 knowledge states at n = 5, all ordered pairs of extensions |
| L4 — `[G]` is connected under reversals that remain in `[G]` | 118,643 | **0** | as L3 |
| AE-B, semantic form — no two distinct undirected edges of an MPDAG induce the same bipartition of `[G]` | 293,814 | **0** | as L3 |

L2 and L3 are the exact computational content of Lemmas AE-B.1 and AE-B.2; a
single failure of either would refute the repair rather than indicate a bug.
AE-B in its semantic form is checked independently of the path construction, so
it corroborates the conclusion by a route that does not pass through Chickering
at all. Verification is not proof, and the proof above is what carries the
claim; the sweep is there because a proof with an arithmetic slip in it and a
sweep with a bug in it rarely agree to 3.5 million cases.
