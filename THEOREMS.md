# THEOREMS — formal statements and status

Session 2, Axis B. Each item states its hypotheses, and is either **proved** (a
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

## 1. Task 0 — the space definition (settled this session)

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

## 3. Theorem 2 — `𝔊_Ĉ` is a join-semilattice, and the join is explicit  *(proved, new)*

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

*Remark.* Meets exist for any pair with a common lower bound
(`[G] ∩ [H] = {D : D ⊨ K_G ∪ K_H}`), but `K_G ∪ K_H` may be inconsistent and the
empty model set is not an element, so `𝔊_Ĉ` is a join-semilattice with top rather
than a lattice, unless a formal bottom is adjoined.

---

## 4. Property S — upper semimodularity  *(PROVED from Lemma R)*

**Statement.** For all `X, Y, Z ∈ 𝔊_Ĉ` with `X ⋖ Y`:

```
X ∨ Z = Y ∨ Z      or      X ∨ Z ⋖ Y ∨ Z .
```

**Proof.** Work on the orientation side, where an element `W` is identified with
its closed knowledge set `K_W` and, by Theorem 2, `W ∨ Z` is identified with
`K_W ∩ K_Z`. Two standard facts are used: the intersection of closed sets is
closed (`cl(A∩B) ⊆ cl(A) ∩ cl(B) = A ∩ B` for closed `A`, `B`), and by
**Lemma R** below `X ⋖ Y` gives `K_X = K_Y ⊔ {e}` for a single orientation `e`.

*Case 1: `e ∉ K_Z`.* Then `K_X ∩ K_Z = K_Y ∩ K_Z`, so `X ∨ Z = Y ∨ Z`.

*Case 2: `e ∈ K_Z`.* Then `K_X ∩ K_Z = (K_Y ∩ K_Z) ⊔ {e}`. Both are closed, so
they are elements, and they differ by the single element `e`. No closed set lies
strictly between two sets differing by one element, so `X ∨ Z ⋖ Y ∨ Z`. ∎

**Verified**, and the case split matches the proof exactly — which is itself a
check that the case analysis is exhaustive:

| scope | triples | Case 1 (collapse) | Case 2 (cover) | gaps of height ≥ 2 |
|---|---|---|---|---|
| all CPDAGs n = 3, 4 (exhaustive) | 191,484 | 129,048 | 62,436 | **0** |
| n = 5, sampled to k = 7 | 1,195,142 | — | — | **0** |
| **total** | **1,386,626** | | | **0** |

The Case-2 key step `cl(S ∪ {e}) = S ∪ {e}` was additionally verified directly on
all 62,436 Case-2 instances, with zero violations.

---

## 4b. Lemma R — every cover adds exactly one orientation  *(reduced to anti-exchange)*

**Statement.** If `X ⋖ Y` in `𝔊_Ĉ` then `|K_X \ K_Y| = 1`.

**Status: follows from Anti-Exchange (§4c) by a classical theorem.**

By Edelman–Jamison, a closure operator satisfies the anti-exchange property iff
its closed sets form a **convex geometry**, and in a convex geometry every cover
in the lattice of closed sets adds exactly one element. So Anti-Exchange ⟹
Lemma R.

**Verified directly** as well: **5,348 covering pairs** (2,658 exhaustive at
n ≤ 4, 2,690 sampled at n = 5), every one of them adding exactly one orientation,
**0 violations**. Equivalently `rank(G) = |dir(G)| − |dir(Ĉ)|` is a rank function,
so the poset is graded — which independently confirms §4.1.

---

## 4c. Anti-Exchange — the single remaining open property

**Statement.** For closed `S` and distinct orientations `x, y ∉ S`:

```
y ∈ cl(S ∪ {x})   ⟹   x ∉ cl(S ∪ {y}) .
```

**Semantic restatement.** Suppose both held. Then
`cl(S∪{x}) = cl(S∪{y}) = cl(S∪{x,y})`, so over the DAGs represented by `S` the
constraints `x` and `y` are *equivalent*. Anti-exchange therefore says exactly:

> **no two distinct orientations are perfectly correlated across the represented
> DAGs.**

**Case A (`x` and `y` orient the same edge): PROVED.** If `x, y ∉ S` and `S` is
closed, that edge is undirected in the graph of `S`, so by maximal orientation
both of its orientations occur among the represented DAGs. Pick `D₁ ⊨ x`. Since
`x` and `y` are opposite orientations of one edge, `D₁ ⊭ y`, so `y ∉ cl(S∪{x})`
and the hypothesis fails vacuously. ∎

**Case B (`x` and `y` orient different edges): OPEN.** Verified with **0
violations** across **5,254** applicable triples (2,178 exhaustive at n ≤ 4,
3,076 sampled at n = 5). A proof would need to rule out two distinct undirected
edges being perfectly correlated over an equivalence class with knowledge
imposed. The natural route is the chordal structure of chain components, which
permits reorienting within a component; the obstruction is that under background
knowledge a component need **not** be chordal — precisely the phenomenon Task 0
uncovered — so the classical argument does not transfer unmodified.

---

## 5. Lemma L  *(proved from Property S)*

**Statement.** For all `G₀, G ∈ 𝔊_Ĉ`: `d(G₀, G ∨ G₀) ≤ d(G₀, G)`.

**Proof from Property S.** Let `G₀ = P₀, P₁, …, P_d = G` be a shortest path in the
covering graph, so each consecutive pair is a covering pair in one direction or
the other. Put `Q_i = P_i ∨ G₀`. Then `Q₀ = G₀ ∨ G₀ = G₀` and `Q_d = G ∨ G₀`.

For each `i`, `{P_i, P_{i+1}}` is a covering pair, so Property S applied with
`Z = G₀` gives `Q_i = Q_{i+1}` or `Q_i ⋖ Q_{i+1}` or `Q_{i+1} ⋖ Q_i`. In every
case `Q_i` and `Q_{i+1}` are equal or adjacent in the covering graph. Deleting
repetitions leaves a walk from `G₀` to `G ∨ G₀` of length at most `d`. Hence
`d(G₀, G ∨ G₀) ≤ d = d(G₀, G)`. ∎

**Verified independently**, without assuming Property S: **521,432** `(G₀, G)`
pairs, **0 violations** — 80,480 over all CPDAGs on 3 and 4 nodes (exhaustive)
and 440,952 at n = 5 across densities.

Slack `d(G₀,G) − d(G₀, G∨G₀)` at n = 5 (440,952 pairs):

| slack | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|---|
| count | 39,637 | 94,498 | 108,655 | 92,864 | 61,076 | 31,048 | 11,058 | 2,116 |

Only **9.0%** of pairs are tight at zero slack (10.5% at n ≤ 4), so the lemma
holds with real margin rather than balancing on the boundary — the more
reassuring of the two possible states, and the distinction the brief asked to be
made explicit. The margin does not shrink with graph size over the range tested.

---

## 6. Conjecture 2 — reduction  *(proved conditional on Property S)*

**Statement.** The nearest failure is reachable from `G₀` by a monotone upward
path, so `r_val` equals the minimum number of retractions from `G₀` that induces
failure.

**Theorem.** Anti-Exchange ⟹ Lemma R ⟹ Property S ⟹ Lemma L ⟹ Conjecture 2.

**Proof.** Let `G` be a nearest failure, `d(G₀,G) = r`. Put `J = G ∨ G₀`.

1. `G ⪯ J`, so by **Theorem 1** `Z` fails at `J`.
2. By **Lemma L** (which follows from Property S), `d(G₀, J) ≤ r`.
3. Property S implies upper semimodularity, hence the Jordan–Dedekind chain
   condition, hence `𝔊_Ĉ` is graded. For comparable `G₀ ⪯ J` every path in the
   covering graph changes rank by `±1` at each step, so has length at least
   `ρ(J) − ρ(G₀)`; and a saturated chain of exactly that length exists. Therefore
   `d(G₀, J) = ρ(J) − ρ(G₀)` and it is realised by a **monotone upward** path.
4. So retraction-only search reaches the failing `J` in `d(G₀,J) ≤ r` steps, and
   it cannot do better than `r` because `r` is the minimum over *all* elements.

Hence the retraction-only radius equals `r`. ∎

**Consequently the session's central open question is no longer Conjecture 2 but
Anti-Exchange Case B** — a single local property of Meek closure, with no
reference to failure, adjustment sets, treatments or outcomes, and with a clean
semantic reading: *no two distinct orientations are perfectly correlated across
the represented DAGs*. Case A is proved. Reducing a conjecture about causal
robustness to that is the main theoretical result of this session.

---

## 7. What this means for `radius_local_up`

`radius_local_up` performs **upward** BFS. Its exactness is therefore *exactly*
Conjecture 2, hence — by §6 — exactly Property S.

| | status |
|---|---|
| `radius_local_up` returns the BFS radius | **Exact if Anti-Exchange Case B holds**; otherwise it can only ever return radii that are **too large**, never too small |
| Property S | **proved**; rests on Anti-Exchange Case B (§4c), which is verified not proved |
| Direction of any error | one-sided: an over-estimate of the radius, i.e. a claim of *more* robustness than is warranted |

The one-sidedness matters for how the method should be reported: a failure of
Property S would make the method **optimistic about robustness**, which is the
dangerous direction, so the conditional must be stated wherever the method's
exactness is claimed.

---

## 8. Summary table

| # | Statement | Status |
|---|---|---|
| T1 | Failure is upward-closed | **Proved** (one line) |
| T2 | `G ∨ H = Meek(Ĉ, K_G ∩ K_H)`; `𝔊_Ĉ` is a join-semilattice | **Proved** (this session), verified on 260,716 pairs |
| — | Space membership = reachability (fixpoint predicate) | **Settled**, verified on 1,588 states |
| — | The chain-space semimodularity refutation | **Refuted** — it conflates an assertion with a covering step |
| — | `𝔊_Ĉ` is graded | Verified on 203 spaces, 0 non-graded; also implied by S |
| S | Upper semimodularity | **Proved from Lemma R**; verified on 1,386,626 triples |
| R | Every cover adds exactly one orientation | **Proved from Anti-Exchange** (Edelman–Jamison); verified on 5,348 covers |
| AE-A | Anti-exchange, same edge | **Proved** |
| AE-B | Anti-exchange, distinct edges | **OPEN**, 0 violations in 5,254 triples |
| L | `d(G₀, G∨G₀) ≤ d(G₀,G)` | **Proved from S**; independently verified on 521,432 pairs |
| C2 | Retraction-optimal witnesses | **Proved from Anti-Exchange**; reduces to AE-B |
