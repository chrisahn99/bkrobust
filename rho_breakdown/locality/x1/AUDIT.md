# X1 DESIGN AUDIT — written against `PREREG.md`, before implementation

Adversarial read of a design I did not write. E1 was refuted 3/3 on *interpretation* while every
number reproduced; the cheapest way to repeat that is to skip this file. X1 is at higher risk than
E1′ was, because X1's two **load-bearing** predictions (P2, P3) are comparisons between an operator
that raises and an operator that is *defined* never to raise, and that asymmetry is worth more than
the effect being hunted.

**Everything below is arithmetic on numbers already on the record** (`arm1_*.json` from E1′,
`RANKING.md:763`, the pilot's `run_linear.py`) plus one **structural** measurement of the ensembles
that computes no X1 outcome: `audit/audit_eligibility.py`, run today on `original` (1 326 SCMs) and
`large` (1 200 SCMs). It records `|N(C)|`, hop compositions, draw-pool compositions and two purely
graph-theoretic bounds. It does **not** implement `bk_assert`, and it measures no catch rate and no
damage rate. P1–P5 remain unspoiled.

---

## A1. 🔴 THE FREE BAR. P2 reads SUPPORTED on the primary ensemble even if the spurious class is *identically* damaging — and even if it is 15 % **less** damaging

This is the E1′-`crit='ident'` failure repeated one level up, and it is computable exactly.

`R_dmg^R` is defined **unconditionally** (§4.2, deliberately). Decompose it:

```
R_dmg^R  =  P(consistent | R)  ×  P(silent | consistent, R)
R_dmg^S  =  P(consistent | S2) × P(silent | consistent, S2)  =  1.000 × P(silent | S2)
```

because S2's `consistent ≡ True` **by construction** (§5.1). The first factor is not a property of
the perturbation class; it is a property of which operator was told to raise. Measured on E1′,
`original`, ρ = 1, `arm1_original.json`, n = 10 548 members over 2 637 SCMs:

| | value | 95 % Wilson |
|---|---|---|
| `P(consistent \| R)` | 6 776/10 548 = **0.6424** | — |
| `R_dmg^R` (unconditional, ρ=1) | 1 007/10 548 = **0.0955** | [0.0900, 0.1012] |
| `P(silent \| consistent, R)` | 1 007/6 776 = **0.1486** | [0.1403, 0.1573] |

Now set the spurious class's **conditional** damage exactly equal to the reversal class's — i.e.
assume X1's claim is *false*, that the two classes are one phenomenon. Then
`R_dmg^S = 0.1486`, and

```
R_dmg^S − R_dmg^R = 0.1486 − 0.0955 = 0.0531   ≥ 0.03    ✅ bar cleared
Wilson at n = 5 304:  R [0.0879, 0.1037]   S [0.1393, 0.1584]   → DISJOINT   ✅
```

**P2 prints SUPPORTED — "the class damages more, S3/S4 do not transfer" — with zero contribution
from the perturbation class.** Solving for the break-even:

| ensemble | `R_dmg^R` | `P(silent\|cons,R)` | SUPPORTED iff `R_dmg^S ≥` | as a ratio to R's conditional rate | REFUTED iff `≤` | ratio |
|---|---|---|---|---|---|---|
| **original** (PRIMARY) | 0.0955 | 0.1486 | 0.1255 | **0.845** | 0.0802 | 0.540 |
| licensed | 0.1085 | 0.1605 | 0.1386 | **0.864** | 0.0903 | 0.563 |
| large | 0.0213 | 0.0257 | 0.0513 | **1.996** | 0.0138 | 0.537 |

Read the `original` row: the spurious class must be **more than 15.5 % relatively *less* damaging**
than the reversal class before P2 stops saying SUPPORTED, and **46 % less damaging** before it says
REFUTED. The INDISTINGUISHABLE window — the one the design says *"is exactly A26's `revive_if`
condition"* and *revives* the b-LOAD transfer — sits at conditional ratios **(0.540, 0.845)**.

> 🔴 **A26's `revive_if` was written to fire when "the two noise classes are one phenomenon". As
> operationalised in §4.2 it fires only when the spurious class is 15–46 % *less* damaging than
> reversals. If the two classes are literally identical, X1 prints the opposite of the truth.**

**Mitigation, binding (M1).**
1. The **primary** contrast for P2 is the **conditional** rate `P(silent | operator did not raise)`
   as a **risk ratio** `RR = P(silent|S) / P(silent|R)`, with a paired cluster-bootstrap CI over
   SCMs (10 000 resamples of SCMs, not statements). Registered bars: SUPPORTED iff `RR ≥ 1.25` and
   `CI_lo > 1`; REFUTED iff `RR ≤ 0.80` and `CI_hi < 1`; INDISTINGUISHABLE otherwise.
2. The unconditional difference stays, **demoted to secondary**, and is printed on one line with the
   two pass rates and this sentence: *"of the difference D, the share attributable to arm S2 never
   raising is `(1 − P(consistent|R)) × P(silent|cons,S2)`; the residual is the class effect."*
   Print that decomposition as two numbers, always.
3. Print `0.845 / 0.540` (the `original` break-even ratios) **in RESULTS**, next to the verdict, as
   the free-bar disclosure. E1′ printed `bar_was_free`; X1 must too.
4. A26's `revive_if` verdict is reported **only** off the conditional RR. Reporting it off the
   unconditional rule would register a false retraction into the negative-results registry.

## A2. 🔴 P3 is free in the same way, and its baseline is not the object it is compared to

`754/8085 = 9.33 % [8.71 %, 9.98 %]` is, per `RANKING.md:763`, one cell of a four-way decomposition:
**47.5 % Meek-FAIL / 23.6 % O\* untouched / 19.5 % loud non-amenable / 9.3 % silent**, over 8 085
**perturbations**. So its denominator is **47.5 % members that raised**. Arm S2 never raises. The
matched level is therefore

```
0.0933 / (1 − 0.475) = 0.1776
```

**A never-raising operator with damage identical to the reversal class prints 0.1776, which is
already outside and above the `revive_if` window [0.087, 0.100].** The registered prediction
(*"does not overlap, lies above, `R^{S-loc} ∈ [0.20, 0.60]`"*) is thus delivered by the operator's
definition. Overlap requires `R^{S-loc} ≤ 0.100`, i.e. a conditional ratio of **0.563** — the same
"must be 44 % milder" structure as A1.

Three further mismatches make "evaluated verbatim" not evaluable as written:

| | the `754/8085` baseline | X1 §4.3 as written |
|---|---|---|
| arms | pilot's **generic + tiered** | one arm |
| ρ | **{1,2,3}** | {1,2,3,4} |
| \|K\| | **{3,4}** (`generic_knowledge` takes `min(\|U\|,4)`, needs ≥3) | 4 exactly (`K_all[:4]`) |
| unit | **member** | §4 says *"all rates below are per statement"* |
| SCMs | `default_rng(20260717)`, `run_linear.py` | `default_rng(20260819)`, `run_arm1.py` — **disjoint draw** |

§8.4 explicitly declares tiered knowledge out of scope while the baseline contains it.

**Mitigation, binding (M2).**
1. Recompute the reversal baseline **in-run from arm R on X1's own SCMs**, per member, pooled
   ρ ≤ 4, `original`. That is the only like-for-like comparator, it costs nothing (arm R is already
   running on the same SCMs), and it removes all five mismatches at once. E1′'s value for it is
   already known: **3 408/39 555 = 0.0862 [0.0834, 0.0890]** unconditional, **0.1861** conditional.
2. The literal `754/8085` figure may be quoted **once**, labelled *"legacy, pilot grid, generic+tiered,
   ρ≤3, |K|∈{3,4}, member unit, disjoint SCM draw"*, and **never** as the thing X1's CI is tested
   against. `RANKING.md:768` names exactly this seam.
3. Report P3 per **member** (the baseline's unit) *and* per statement, both labelled; do not let §4's
   "per statement" default silently redefine the comparison.
4. State in RESULTS that the `revive_if` as written cannot fire under a never-raising operator
   unless the class is ≥ 44 % milder, and propose the amended condition to Zé rather than executing
   a retraction on the old one.

## A3. 🔴 One absolute bar (0.03) across three ensembles whose baselines span 0.021 → 0.109

On `large` the same 0.03 bar demands `R_dmg^S ≥ 0.0513`, a **2.0×** conditional effect; on `original`
it demands **0.845×**. A single registered threshold that is free on one grid and near-unreachable on
another is not a threshold. (E1′'s G1 gate already **failed** on `large` — 0.026 against [0.12, 0.18]
— for the same reason: `large`'s reversal damage rate is a different order of magnitude.)

**Mitigation, binding (M3).** Register the risk-ratio bars of M1 on **all three** grids; keep the
absolute-difference bar only on `original`, with its per-grid free level printed next to it.

## A4. 🔴 "Meek-consistent" is not merely *broken* for arm S — it is **undefined**, so the S3 comparison line is not like-for-like

`apply_background_knowledge` (`graphs.py:179-199`) implements Perković et al. Algorithm 1: its
predicate is *"∃ a DAG in **MEC(C)** satisfying K"*. Every member of MEC(C) has `skeleton(C)`. A
spurious statement asserts an edge **absent from that skeleton**, so **no** DAG in MEC(C) satisfies
K′ — for arm S the answer to "is K′ consistent with C?" is *no*, always, trivially, for every member.
S1's pass rate is therefore **a different predicate**, and §7's S3 row — *"ARM R pass rate vs ARM S1
pass rate vs ARM S2 pass rate ≡ 1, printed on one line"* — puts three incommensurable numbers on one
line. That is the exact failure mode `RANKING.md:768` warns about.

Worse, the S1 predicate as specified is **contaminated by a purely skeletal artefact**. `v_structures`
(`graphs.py:67-76`) requires the two parents to be **non-adjacent**. Adding the edge `{a,b}` therefore
*deletes* from `v_structures(G)` every unshielded collider of `C` whose parent pair is exactly
`{a,b}` — a change with no incoherence in it whatsoever. Measured over the whole `N(C)` pool:

| ensemble | P(a uniform `N(C)` pair is the parent pair of an existing unshielded collider of C) |
|---|---|
| original | **1 750/15 992 = 0.1094** |
| large | **13 400/215 006 = 0.0623** |

So `flags.vstruct` fires for a definitional reason on ≥ 10.9 % of `original` S-uni draws at ρ = 1
before any real incoherence is considered — a floor inside `C_S1`, which is **P1's only empirical
content**. And a second implementer could just as defensibly compare against `v_structures(G0)`, or
against the v-structures of the CPDAG of the *augmented* skeleton, and get three different numbers.

**Mitigation, binding (M4).**
1. **The primary S1 coherence predicate is intrinsic and skeleton-agnostic: does `G` admit a
   consistent DAG extension?** Implement **Dor–Tarsi (1992)** `pdag_extendable(G)` (≈ 25 lines:
   repeatedly find a vertex `x` with no outgoing directed edge whose undirected neighbours are
   adjacent to all its other neighbours; orient into `x`; delete `x`; fail if no such vertex exists).
   `consistent_S1 := pdag_extendable(G) and not has_directed_cycle(G)`.
2. **Gate it** (new **G11**): on every SCM with `|undirected_edges(G)| ≤ 12`, check
   `pdag_extendable(G) == (len(consistent_dag_extensions(G, ref_vstructs=v_structures(G), limit=1)) > 0)`.
   Pass = 100 %. `consistent_dag_extensions` already exists at `graphs.py:211` and `limit=1` makes it
   cheap; the `≤ 12` cap keeps `2^|U|` bounded (|U| mean 4.73 on `original`, 6.11 on `large`).
3. Keep the design's `v_structures(G) != v_structures(C)` test but rename it `vstruct_vs_C` and
   report it as a **diagnostic**, decomposed into four labelled counts, each with its n:
   `conflict` · `cycle` · `collider_of_C_shielded_by_the_new_edge` (the artefact) ·
   `new_unshielded_collider`. `C_S1` is never quoted without that decomposition.
4. §7's S3 row must print R's pass rate and S1's pass rate on **separate lines** with their
   predicates spelled out, and must state that S3's "two thirds" is a statement about
   consistency-with-C that has no arm-S analogue.

## A5. 🔴 `O*` is evaluated on objects for which it is not defined, and the design never says how often

`optimal_adjustment_set` (`adjust.py:90-97`) is Henckel–Perković–Maathuis `O*`, defined for **amenable
MPDAGs**. Under S2 nothing raises, so directed cycles and non-extendable PDAGs enter the estimand
pipeline; `possibly_causal_paths`, `poss_de` and `forb` all terminate on such objects and return
*numbers*, not errors. A purely structural lower bound on the cyclic share — ordered non-adjacent
pairs `(a,b)` for which `C` **already** contains a directed path `b ⇝ a` using directed edges only:

| ensemble | lower bound on P(cycle) per asserted statement |
|---|---|
| original | 375/31 984 = **0.0117** |
| large | 10 037/430 012 = **0.0233** |

This is a *floor*: it ignores the orientations Meek closure adds after `K` is applied, and it ignores
the ρ−1 other added edges at ρ ≥ 2. The true S2 cyclic share at ρ = 4 will be several times larger.

**Mitigation, binding (M5).** Every S2 member carries `mpdag_valid ∈ {True, False}` from the Dor–Tarsi
oracle of M4 and `has_cycle`. Every rate in §4 is reported **twice**, once on `mpdag_valid=True` and
once on the full set, with both n's. Any RESULTS sentence containing the words *optimal adjustment
set*, *variance-optimal*, or *identified* is restricted to the valid stratum. The invalid stratum is
labelled **`ill_posed`** and is a *finding about b-LOAD's semantics* (b-LOAD really does build such
objects) — not a number about `O*`.

## A6. 🔴 The `P-null` dead switch is not a null: it deletes a true statement, so G5 will fail for a non-bug reason

§1.3 defines `P-null`'s replacement as *"a pair already **directed** in `C` with the orientation it
already has"*, and §1.3's `K'(F)` **replaces** slot `t`. So the true statement `(u_t, v_t)` is
**dropped** from `K′`. `G0 = apply_background_knowledge(C, K)` used all four. Dropping one leaves
`{u_t, v_t}` undirected unless Meek recovers it from the others, so `G ≠ G0` routinely and `O*` moves.
**G5's registered pass condition (`G = G0` in 100 %, `ostar_changed = 0`) will fail, and a failing
dead switch reads as "the harness is broken".** The design summary's own words — *"redundant
assertion of an already-directed edge"* — describe a different, correct construction than §1.3's table.

**Mitigation, binding (M6).** Two placebo variants, both run:
- **P-null-identity**: `r_t := (u_t, v_t)` — the same true statement, passed through the full
  `bk_assert` code path. Registered: `G = G0` in **100 %**, `ostar_changed = 0`, `silent = 0`.
- **P-null-append**: `K′ := K ++ [(a,b)]` for `(a,b)` already directed in `C` with its own
  orientation, `|K′| = 5`. Registered: `G = G0` in **100 %**.

A placebo that removes information is a *treatment*. Fix §1.3's table before implementation.

## A7. 🔴 G3 will fail on S-loc, and the reason is a composition fact the design already half-knows

§4.1/G3 register *skeleton-detectability ≡ **1.000**, definitional, any value < 1 is a bug*. That is
true of **S-uni** and **false of S-loc**. b-LOAD's pool `{(a,b) : a ∈ {X,Y}, (a,b) ∉ E(D)}` admits
`(a,b)` whenever `b → a ∈ D` — §0 says so — and such a pair is **adjacent** in `C`, hence not in
`N(C)`, hence not skeleton-detectable. Measured over the pool:

| ensemble | pool size (mean) | `pure_spurious` | `reversal_of_true_edge` | `reorient_cpdag_edge` |
|---|---|---|---|---|
| original | 8.43 | **0.7313** | **0.2687** | 0.0000 |
| large | 35.81 | **0.9207** | **0.0793** | 0.0000 |

So S-loc's skeleton-detectability is ≈ **0.731** on `original`, not 1.000, and G3 as written fires a
false alarm. (The zero in the last column is worth keeping: the pool partitions cleanly into
pure-spurious and reversal-of-a-true-edge, so post-hoc stratification is exact.)

**Mitigation, binding (M7).** G3 is scoped to **S-uni only**. For S-loc, skeleton-detectability is a
**measured** quantity whose expected value is the `pure_spurious` share; register it as a *prediction*
(0.73 ± / 0.92 ±) rather than a gate. Likewise S1's `conflict` branch fires on the reversal stratum
by construction, so S-loc's `C_S1` is a mixture rate and must never be quoted un-stratified.

## A8. 🔴 P3's headline is about reversing the query edge on half the SCMs

`S-loc` anchors every statement at `X` or `Y`, so `(X,Y)` and `(Y,X)` are in the pool whenever they
are not true edges of `D`. Measured:

| ensemble | E[share of one S-loc draw that is the query pair] | P(≥ 1 of the 4 slots is a query-pair statement) |
|---|---|---|
| **original** | **0.1626** | **0.4855** |
| large | 0.0393 | 0.1469 |

Asserting `Y → X` when `D` has `X → Y` is not a spurious edge at all; it is the maximally damaging
single statement available, and it lands in ~49 % of `original` SCMs. P3's predicted band
`[0.20, 0.60]` is substantially a measurement of that. By contrast **S-uni is clean**: only 2.40 % of
draws are the query pair (7.52 % of SCMs get one), because `X,Y` are non-adjacent in only 22.85 % of
`original` SCMs and `|N(C)|` averages 12.06.

**Mitigation, binding (M8).** S-loc's damage rate is reported as a **3-row table** —
`pure_spurious, non-query` · `pure_spurious, query-pair` · `reversal_of_true_edge` — each with n and
Wilson CI, plus the pooled row. **The `pure_spurious, non-query` row is P3's primary.** The pooled
S-loc number may not appear in any sentence without the composition table on the same page.

## A9. Conditioning that masquerades as an effect: the arms do not see the same SCMs, and do not see the same distances

Two distinct eligibility/composition asymmetries, both measured.

**(a) Eligible-SCM set.** §1.3 says S-uni draws four ordered pairs *"uniformly without replacement"*
from `{(a,b) : {a,b} ∈ N(C)}`. Read over **unordered** pairs the arm needs `|N(C)| ≥ 4`; read over
**ordered** pairs it needs `|N(C)| ≥ 2`. On `original`:

```
|N(C)| ≤ 1 :   8/1326 = 0.0060      |N(C)| ≤ 3 : 109/1326 = 0.0822
```

The two readings disagree on **7.6 % of the primary ensemble** (101/1 326) — and the SCMs they
disagree about are the **densest** ones, exactly where `A(C)` is large and the causal neighbourhood is
richest. Arm R has no such exclusion. (`large` is immune: `min |N(C)| = 77`.) The ordered reading has
a second consequence: it permits drawing **both** `(a,b)` and `(b,a)`, in which case the second
assertion overwrites the first, `|K′|` carries only three distinct edges, and the result depends on
`K′` order — which also makes G8 (close-after-each vs stamp-all) order-dependent for a spurious reason.

**(b) Hop composition.** The two arms draw from different pools, and the pools sit at different
distances from the query:

| ensemble | pool | hop 0 | hop 1 | hop 2 | hop ≥3 | unreachable |
|---|---|---|---|---|---|---|
| original | `N(C)` (S-uni) | 0.4920 | 0.3977 | 0.0595 | 0.0080 | **0.0428** |
| original | `K ⊂ U(C)` (R) | 0.6346 | 0.2666 | 0.0402 | 0.0068 | **0.0518** |
| large | `N(C)` | 0.1820 | 0.2939 | 0.2360 | 0.2740 | **0.1141** |
| large | `K ⊂ U(C)` | 0.2108 | 0.2206 | 0.2127 | 0.1870 | **0.1690** |

On `original` arm R's statements are **27 % more likely to be adjacent to the query** than S-uni's.
Since E1′ measured essentially all movement at hop 0, this biases the P2 contrast **against**
SUPPORTED, in the opposite direction to A1's bias — and the design states neither. Two unstated,
opposite-signed confounds do not cancel; they make the number uninterpretable.

**Mitigation, binding (M9).**
1. Fix the draw law explicitly: **without replacement over unordered pairs of `N(C)`, then orient each
   drawn pair by an independent fair coin from the same spawned stream.** This removes the
   both-directions pathology and makes `|N(C)| ≥ 4` the single eligibility rule.
2. **Every arm runs on the same SCM set.** An SCM ineligible for any arm is dropped from **all**
   arms; the dropped count and its `(p, deg, |N(C)|, |O0|)` profile are printed. Otherwise G1's
   "100 %" passes vacuously on an intersection while the §4 rates use different denominators.
3. P2 is reported **three ways**: raw, restricted to hop 0, and **direct-standardised** onto arm R's
   hop distribution (weights = R's hop shares, printed). The standardised contrast is the one that
   answers "is the class worse", the raw one answers "is the b-LOAD regime worse".

## A10. The run will not deliver the registered n, and §2.4's premise is wrong by a factor of three

§2.4 states *"E1′ retained 2 637 / 2 087 / 2 443 analysed SCMs from **4 000 draws**"*. The logs say
otherwise: `[arm1/original] kept 3372 SCMs from <= 12000 draws`, `[arm1/licensed] kept 4000 from
<= 12000`, `[arm1/large] kept 3000 from <= 9000`. `n_scm` is a target on **kept** records (which
include `mpdag_amenable=False`), `OVERSAMPLE` multiplies it into the job count, and `original` was
**oversample-bound**, not `n_scm`-bound.

I ran the exact filter chain at X1's registered setting (`n_scm = 2000`, `OVERSAMPLE = 3`, ⇒ 6 000
jobs) and got **1 326 analysed SCMs**, not 2 000 — yield 0.2210, matching E1′'s 2 637/12 000 = 0.2198.

```
original  n_scm=2000, OVERSAMPLE=3  ->  1326 analysed  ->  5 304 statements/arm at rho=1
licensed  n_scm=2000  ->  ~1044 analysed (kept x 2087/4000)  ->  ~4 176 statements/arm
large     n_scm=1200  ->  1200 analysed within 3 600 jobs    ->   4 800 statements/arm
```

The MDD arithmetic survives (0.0189 at n = 5 304), but the **primary ensemble runs at 50 % of E1′'s
power on the very comparison E1′ was powered for**, and P4's denominator (A12) is bought with no
margin.

**Mitigation, binding (M10).** State every target in **analysed** SCMs. Set `OVERSAMPLE = 5` for
`original`, `6` for `licensed`, `4` for `large`; assert the achieved analysed count at the top of
RESULTS; abort the verdict if `original` analysed < 1 800.

## A11. The sentinel is a stratum, not a curiosity — and it is a structural zero

`hop = 1 048 576` means **both endpoints are unreachable from {X,Y} in `skeleton(C)`**. It is
**4.28 %** of the `N(C)` pool on `original` and **11.41 %** on `large`. §1.6 and G10 already forbid
pooling it — good — but the design does not say what it *is*: adding an edge between two vertices in
a component that does not contain `X` or `Y` **cannot** change `cn(X,Y,·)`, `forb`, or `O*`. These
statements are **guaranteed harmless by construction**, and they dilute every S-uni denominator.

**Mitigation, binding (M11).** Add to §5 (declared definitional): *"unreachable statements have
damage rate 0.000 by construction"*, gate it (**G12**: unreachable-stratum `ostar_changed` = 0, exact),
and report every §4 rate on the **reachable** denominator as primary, with the diluted rate secondary.
Note that the dilution differs by ensemble (4.3 % vs 11.4 %), so cross-ensemble comparisons of
undiluted rates are not comparable.

## A12. P4's rule is a near-certain SUPPORTED, and its statistic answers the wrong question

**Feasibility first** — the denominators exist, but only just:

```
original: 4 x 1326 x 0.4652 (hop>=1, unreachable removed) = 2 467   >= 2 000  ✅ (no margin)
large   : 4 x 1200 x 0.7039                               = 3 379   ✅
```

**But REFUTED is `0 of ≥ 2 000`**, a strict-zero bar. E1′ already found 6.5 % of distance-1
statements uncensored on `licensed`; at n = 2 467 even a rate of 0.0005 yields ~1 event. SUPPORTED
(`≥ 5` events) corresponds to a rate of **0.002** — a count threshold that carries no information
about magnitude, and one that any nonzero mechanism clears. P4 will say SUPPORTED, and the sentence
it licenses ("locality is not a property of knowledge") will rest on 5 events out of 2 467.

Separately: `hop` is computed on `C` **before** the edge is added (§1.6 — correct, it is the
analyst-facing quantity), but the *mechanism* by which a distant assertion damages is that it
**pulls itself into the neighbourhood**. A `hop_C`-only table cannot distinguish "distant knowledge
matters" from "asserted knowledge is never distant once asserted", and S9 turns on exactly that.

**Mitigation, binding (M12).**
1. Record **both** `hop_C` (primary, analyst-facing) and `hop_G` (post-assertion, post-closure), and
   report P4 as a `hop_C × hop_G` cross-tab.
2. Replace the count rule with a **rate + Wilson CI + a materiality ratio**: SUPPORTED iff the
   hop ≥ 1 damage rate is ≥ 0.10 × the hop-0 rate with `CI_lo > 0`; the `≥ 5 events` count survives
   only as an **existence** check, labelled as such.
3. Print the hop ≥ 1 denominator *before* the verdict; if it is < 2 000 the verdict is not readable.

## A13. Ambiguities that would let two implementers produce different numbers

Each of these must be resolved **in the PREREG**, not in the code:

| # | ambiguity | binding resolution |
|---|---|---|
| 1 | "ordered pair without replacement" (A9a) | M9.1 — unordered draw + independent coin |
| 2 | P3's unit: statement or member (A2) | member primary, statement secondary, both labelled |
| 3 | `flags.vstruct` reference graph: `C`, `G0`, or `CPDAG(skeleton(G))` (A4) | M4 — intrinsic extendability is primary; `vs C` is a labelled diagnostic |
| 4 | re-stamp scope: current statement only, or all of `K′` (b-LOAD's `g[mask] = bk[mask]` re-stamps **all** of `K` at every update) | re-stamp the **full `K′` prefix applied so far**; the current-only variant is the G7 comparison |
| 5 | outcome (c) `silent` drops E1′'s `ostar_changed` clause | **checked: empirically null.** On all three E1′ ensembles the two definitions give identical counts (Δ = 0; `changed_still_valid = 0`), because `O0` is always valid. Keep E1′'s definition **verbatim** for comparability and state the equivalence. |
| 6 | `expels_true_dag = not dag_agrees_with(D, G)` (`run_arm1.py:97`) | `dag_agrees_with` first tests `skeleton(D) == skeleton(G)` (`graphs.py:258`), which is **False for every S member with an added edge**. So this field is **≡ True by construction** for arm S. Add it to §5 (declared definitional) or drop it; do not report it. |
| 7 | `cascade` (`run_arm1.py:98`) counts new directed edges minus `|K′|`, mixing *added* edges with *reoriented* ones | redefine for arm S as `len(directed_edges(G)) − n_cpdag_dir − n_added − n_reoriented`, or drop |
| 8 | `n_scm` means kept-records, not analysed SCMs (A10) | M10 |

## A14. Checked and cleared — reported so it is not rediscovered late

- **Intra-SCM clustering is negligible on arm R.** The 4 statements of a ρ = 1 ball share one SCM, so
  Wilson's iid assumption is suspect. Measured design effect (variance of the per-SCM silent count
  against binomial): **0.993** (original), **1.083** (licensed), **1.092** (large) — ICC ≤ 0.031.
  Naive Wilson half-widths are understated by at most **×1.045**. *Binding (M13): re-measure the same
  quantity on arms S1/S2 and print it; if any arm's design effect exceeds 1.25, switch every CI in
  §4 to the SCM-level cluster bootstrap.* The paired structure (same SCMs both arms) still means
  independent-sample CI disjointness is the wrong instrument for the **difference** — M1's paired
  bootstrap fixes that regardless.
- **Oracle leakage: no leak found on the instrument path.** `hop_dist` runs on `C` (`run_arm1.py:46`),
  `N(C)` = `N(D)` because `skeleton(C) ≡ skeleton(D)`, `tau` enters only outcomes (c)/(f), and
  `is_valid_adjustment_set(D, …)` is evaluation-only. **Two labelling caveats, binding (M14):**
  (i) S-loc's pool is defined by `(a,b) ∉ E(D)` — it is an **oracle-guaranteed-false** sampler, which
  is b-LOAD's own design but must be said in §8.8; (ii) the `pure_spurious` / `reversal_of_true_edge`
  labels are **oracle labels**, so a stratified S-loc rate describes a partition the analyst cannot
  perform — it may not be presented as a triage rule.

## A15. The PRIMARY ensemble is the one where the estimand is most degenerate

`original` is primary because 67 % / 14.1 % were measured there (§2.2) — a sound reason. But measured
on the same 1 326 SCMs: **`|O0| = 0` in 65.84 %**, `|O0| ≤ 1` in 87.18 %, mean 0.508. On `large`,
`|O0| = 0` in 23.42 %, mean 1.670. So on the primary grid, two thirds of the queries have an **empty**
optimal adjustment set and `est0` is the unadjusted regression coefficient; "adjustment-set fragility"
is being measured at the boundary of the object. E1′'s G1 gate *failing* on `large` (0.026) is the
same fact seen from the other side.

**Mitigation, binding (M15).** Keep `original` primary for the S3/S4 comparison, and promote `large`
to **co-primary for the mechanism claims (P4, P5)** — it is the grid where `O*` is non-trivial and
where the spurious class's `p²`-vs-`p·deg` prevalence argument (§2.2) actually bites. Print
`P(|O0| = 0)` next to every damage rate; a rate conditioned on `|O0| = 0` and a rate conditioned on
`|O0| ≥ 1` are different quantities.

## A16. Budget: four balls, not three, plus a new oracle

§2.5 says *"three balls per SCM instead of E1′'s one, so ≈ 3× E1′'s per-SCM cost"*. §2.1 lists
**five**: `R`, `S-uni`, `S-loc`, `P-null`, `P-frozen`. `P-frozen` is free; the other four are not — so
**≈ 4×**, plus M4's Dor–Tarsi call per member, plus the fact that adding up to 4 edges to a `large`
graph inflates `possibly_causal_paths` (`adjust.py:39`, an uncapped simple-path DFS) super-linearly.
E1′'s entire arm1+arm2 sweep was ≈ 6 min at 12 workers, so the headroom is real but not what §2.5
claims.

**Mitigation, binding (M16).** (i) `possibly_causal_paths` gets a hard cap of 200 000 enumerated
paths per call; a member that trips it is labelled `path_blowup` and reported with its n, never
silently dropped. (ii) Run the 48-SCM smoke on `large` **first**, not on `original`. (iii) `deg ∈
{4,6}` stays excluded from `large` (§2.5) — keep that line.

## A17. What must be reported even if it embarrasses the design

- The A1 break-even ratios (**0.845 / 0.540** on `original`) printed beside P2's verdict, and the
  decomposition of the unconditional difference into "operator never raises" and "class effect".
- The A2 statement that A26's `revive_if`, as written, cannot fire unless the class is ≥ 44 % milder.
- The achieved **analysed** SCM count and the hop ≥ 1 denominator, **before** any verdict (A10, A12).
- `C_S1` decomposed four ways (M4.3), never as a single number.
- The `mpdag_valid = False` share at every ρ (A5), and every §4 rate reported on both strata.
- The S-loc composition table (A8) on the same page as any S-loc number.
- Gates G1–G12 — including the two that this audit predicts will fail as currently written (**G3 on
  S-loc**, **G5 as specified**) — pass or fail, before any verdict.
- §5.3's S2-vs-S7 testability contradiction. **This audit agrees it is real and agrees it is the most
  valuable thing X1 has to say**, and adds: it is *not* X1's to repair, and RESULTS must not present a
  wording proposal as a measurement.

## A18. 🔴 Carried forward, unclosed, and now dated

§10's two vault items stand and this audit re-verified neither is fixed: `active-claims.md` holds no
ρ-breakdown entry, `grep -rl 'relevance:.*rho-breakdown' wiki/literature/` returns **0**, and
`active-claims.md:570` still records the b-LOAD graph-distance falsifier as *"currently NOT firing"*
on the strength of a number **retired on 2026-08-14**. The b-LOAD CIKM camera-ready is **2026-08-23 —
four days out**. X1 supplies a measurement that bears on that register line; it does not edit it, and
it must not be written into a manuscript before Zé registers the claim in his own words
(CLAUDE.md Operating Rule #1). **Proposed, not committed by an agent.**

---

### Verdict

The design is unusually careful — the definitional declarations in §5, the sentinel discipline in
§1.6, the "NOT RUN" list in §8 and the honest naming of the S2/S7 contradiction are all better than
the run they guard. But **the two predictions the design itself labels *load-bearing* (P2, P3) are
cleared before the experiment starts**, by an operator-definition asymmetry worth 0.053 against a
0.03 bar; two gates (G3 on S-loc, G5 as written) will fail for non-bug reasons; and the S-arm
consistency predicate that S3's comparison line depends on is undefined rather than merely different.
None of this is fatal to the *run* — the same five balls, with M1–M16 applied, answer the question
the abstract actually needs answered. **Implement M1–M16 first; do not implement §4.2 and §4.3 as
currently registered.**

### Audit artefacts

`audit/audit_eligibility.py` · `audit/summ.py` · `audit/audit_elig_original.json` (1 326 SCMs) ·
`audit/audit_elig_large.json` (1 200 SCMs). Run on betelgeuse, `~/e1venv/bin/python`, 12 workers,
`~/latent-causal/x1-spurious/audit/`. Structural only: no `bk_assert`, no catch rate, no damage rate.
