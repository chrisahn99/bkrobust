# VERDICT — X1 (spurious required edges) and X2 (locality)

**Written 2026-08-19 by the synthesiser.** Inputs: `x1/RESULTS.md`, `x2/RESULTS.md`, their
pre-registrations and audits, three independent audits of each, and my own recomputation from raw
data on betelgeuse and from the shipped JSONs.

Everything marked ✅ below I recomputed myself. Everything marked ⚠️ I could not verify and is
flagged as such. No number appears without its denominator.

---

## 1. Headline

**Both experiments ran, both are arithmetically sound, and both are mis-titled: X1 built the
spurious-edge operator the abstract had been promising with no code behind it and found that a
single spurious required edge one hop from the query damages at essentially the adjacent rate
(133 of 3,220 vs 167 of 4,123, `original`) while a single reversal one hop out damages nothing
(0 of 2,560), and X2 found that the pre-registered falsifier it was built to fire returned zero —
0 counterexamples in 12,799 far single-misstatement trials on the archive and 0 in a complete
enumeration of every CPDAG, knowledge state and query at p ≤ 5 (0 of 136,200).**

**Put together, the two runs agree on a single sharp fact that neither headline states: the
elicitation radius for the optimal adjustment set is 1, not 0 and not unbounded — beyond one hop,
O\* did not move in either experiment under either operator (X1: 0 of 931 at hop ≥ 3, pooled;
X2: 0 of 6,726 amenable perturbations whose closest misstatement is at distance ≥ 2) — while
*identification status* changes at distance 2 (360 of 3,600 exhaustively at p = 5), which is the
half of the claim the abstract actually needs and the half both reports blur.**

---

## 2. The fate of each abstract sentence

Abstract under test: `~/rho-breakdown-iclr2027/abstract/round2/C-triage.tex`. **Not edited.**
Replacement text proposed only.

### S3 — WEAKENED

> **Current:** "Two thirds of false knowledge sets pass the consistency check"

| | |
|---|---|
| **Deciding number** | ✅ 0.6441 (k = 5,153, n = 8,000 ρ=1 reversals, X1 `original`, fresh SCMs). Across E1′'s six ensemble/arm cells the pass rate spans **0.5128 to 0.8278** (n_SCM 127–2,637). The pilot's own grid gives 0.6424, not 0.67. |
| **Why not STANDS** | It is a **reversal-only** number, printed in the same paragraph as S7, which promises two error classes. It is also ensemble-dependent by a factor the sentence hides. |

**Replacement:** *"Two thirds of false knowledge sets built from reversed orientations pass the
consistency check (0.64, n = 8,000 single reversals), and the rate ranges from 0.51 to 0.83 across
graph densities."*

---

### S4 — WEAKENED

> **Current:** "and $14.1\%$ of the consistent errors go on to bias the estimated effect"

| | |
|---|---|
| **Deciding number** | ✅ X1 measures 0.1451 [0.1354, 0.1554] (k = 688, n = 4,741 consistent reachable ρ=1 reversals, `original`) — reproducing the pilot's 14.1% (194/1,379). But across E1′'s six cells the same quantity is **0.0257 to 0.2492**, a 9.7× spread (n_SCM 127–2,637). For the *spurious* class it is **0.0398** [0.0357, 0.0444] (k = 309, n = 7,762). |
| **Why not STANDS** | Same reversal-only scope; a 9.7× ensemble spread reported as a point estimate with no CI; and **it is not a property of "knowledge errors" but of reversals specifically** — see the mechanism finding in §5. |

**Replacement:** *"Of the consistent reversals, 14.5% [13.5, 15.5] bias the estimated effect
(n = 4,741); the rate ranges from 2.6% to 24.9% across ensembles. A spurious required edge biases
it in 4.0% of cases (n = 7,762)."*

---

### S5 — WEAKENED (first clause DEAD, second clause survives with a new denominator)

> **Current:** "The optimal adjustment set is a local functional of the graph, and no misstatement
> at graph-distance at least $1$ from the query moved it ($0$ of $791$)."

| | |
|---|---|
| **Clause (a) "local functional"** | **DEAD.** ✅ X1: a spurious required edge at hop ≥ 1 moves O\* and invalidates it in **142 of 3,639** (`original`, ρ=1). ✅ X2: **70 of 10,615** ρ≥2 knowledge sets with *every* misstatement at distance ≥ 1 change O\*, all 70 invalid in the true DAG. |
| **Clause (b) the zero** | **SURVIVES, per statement, for reversals at ρ=1, with a much larger denominator.** ✅ 0 of 12,799 (E1′ archive, far ρ=1, D_am); ✅ 0 of 25,257 (X2 ARM E, reversal sub-operator only — see §4 bug #3); ✅ **0 of 136,200 in a complete enumeration at p ≤ 5** (8,967 CPDAGs, 118,593 MPDAGs), against 184,332 of 307,584 near. |
| **The number 791** | **Has no provenance and must be deleted.** ✅ I verified on betelgeuse: the pilot codebase contains **no graph-distance, hop, or BFS function of any kind**; `791` appears **0 times** in the pilot's `summary.json`, and every match in the raw JSONs is a fragment of a float (`0.2913107791813075`). The only `hop_dist` in the project is `e1prime-se/code/run_arm1.py:46`, written **today**. The number was not mis-stated — the instrument that would have produced it did not exist in the cited artefact. |

**Replacement** (this is also the wording X2's own pre-registration pre-committed at
`x2/PREREG.md` §4.1, before any number existed):

> *"Single reversed orientations at skeleton distance at least 1 from the query never changed the
> optimal adjustment set (0 of 25,257 misstatement trials at p ≤ 25, and 0 in a complete
> enumeration of every knowledge state at p ≤ 5). Joint reversals did: 70 of 10,615 knowledge sets
> whose every misstatement lies at distance at least 1 changed it, each of them with a misstatement
> exactly one hop out. A single spurious required edge one hop out changed it in 142 of 3,639."*

---

### S7 — STANDS, but only as of today, and it forces S3/S4 to be re-scoped

> **Current:** "Statement errors cover reversed orientations and required edges asserted between
> non-adjacent variables."

| | |
|---|---|
| **Deciding number** | Before X1 this sentence had **zero code** behind it in either codebase: `graphs.py:192` raises `MeekFail("i-j not an edge of the CPDAG")` on any non-adjacent assertion, so the class was definitionally impossible. X1 implemented a genuine new operator (`bk_assert`: add the edge to the skeleton, then close) and ran it on 5,200 SCMs. |
| **Caveat that must ship with it** | The two classes are **not interchangeable**, so S3 and S4 may not be printed as if they covered both. ✅ At ρ=1 the spurious class biases at 0.0398 vs the reversal class's 0.1451; ✅ by ρ=4 the two are indistinguishable on the primary grid (RR = 1.119, S 287/1,999 vs R 72/561). |

**Addition required:** *"The two classes behave differently: at a single error the spurious class
biases the estimate roughly a third as often, and by four simultaneous errors the two are
indistinguishable."*

---

### S8 — DEAD

> **Current:** "Locality restricts the enumeration to statements near the query, so cost does not
> scale with $|K|$."

| | |
|---|---|
| **Deciding number** | ✅ Pruning to radius 0 changes the reported ignorance interval in **Q₀ = 0.1244** (n_SCM = 2,637, `original`) and **0.1452** (n_SCM = 2,087, `licensed`) — both with Clopper–Pearson lower bounds above the registered 0.01 REFUTED bar, and both meeting the registered n_SCM ≥ 2,000. At radius 1 the pruner is safe (Q₁ = 0.0011 / 0.0014) but **prunes almost nothing**: mean Cost₁ = 0.827 / 0.808, i.e. it discards under a fifth of the ball. |
| **Second, independent reason** | **No pruned enumerator exists.** Neither codebase has a distance filter; both enumerate `combinations(range(|K|), rho)` over all of K. "Cost does not scale with |K|" describes an algorithm that has never been written or benchmarked. X2 measured pruning *validity* and never measured its *cost* in wall-clock. |

**Replacement:** *"Errors far from the query are rarer but not absent, so restricting the
enumeration is a heuristic rather than a property: pruning to the query's immediate neighbours
changes the reported interval in 12–15% of systems, and pruning to one hop is safe but removes
under a fifth of the ball."*

---

### S9 — DEAD

> **Current:** "Knowledge outside that neighbourhood need not be elicited."

| | |
|---|---|
| **Deciding number** | ✅ Identification status changes at distance **2**, not 1: X2's ARM E records 7 identification losses from single misstatements at dmin = 2 (licensed 3 of 6,485, large 3 of 17,010, k8 1 of 306), and the **exhaustive** p = 5 enumeration records **360 of 3,600** at dmin = 2. ✅ X1 adds that a single *spurious* edge one hop out biases silently in 142 of 3,639. |
| **The other direction** | X2 also measures 29,364 `E_gain` events — a far statement turning a non-identified query into an identified one. ⚠️ These fire only when the base was already non-amenable, i.e. outside the regime where a ρ-ball is being enumerated at all, so they support S9's refutation only weakly. I would not print them. |

**Replacement:** *"Knowledge on the query's immediate neighbours must be elicited: statements one
hop out bias the estimate, and statements two hops out can change whether the effect is identified
at all."*

---

### Summary table

| sentence | verdict | deciding number (with n) |
|---|---|---|
| **S3** | WEAKENED | 0.6441 (n = 8,000); range 0.51–0.83 over 6 cells |
| **S4** | WEAKENED | 0.1451 (k=688, n=4,741); range 0.026–0.249 over 6 cells |
| **S5** | clause (a) **DEAD**, clause (b) WEAKENED, "791" **DELETE** | 142/3,639 and 70/10,615 kill (a); 0/136,200 exhaustive rescues (b); 791 has zero provenance |
| **S7** | **STANDS** (newly earned by X1) | operator built and run on 5,200 SCMs |
| **S8** | **DEAD** | Q₀ = 0.1244 (n=2,637) / 0.1452 (n=2,087); Cost₁ = 0.83 |
| **S9** | **DEAD** | 360/3,600 identification changes at dmin = 2, exhaustive at p = 5 |

---

## 3. The contradiction ledger, resolved

**This is the paragraph that matters most, so it is stated flatly: E1′'s 0.935 and "0 of 791" are
not in contradiction, they never were, and the appearance of contradiction came from comparing
three different things at once.**

E1′ reported that for the `licensed` ensemble, SCMs whose statements all sat at hop-distance ≥ 1
from the query had their breakdown radius censored 93.5% of the time (n = 229 SCMs) — implying
6.5% were *not* censored, i.e. something at distance ≥ 1 moved something. "0 of 791" asserts
nothing at distance ≥ 1 ever moved O\*. The two differ on **three independent axes**:

1. **Unit.** E1′'s strata are **per SCM**: `analyse_arm1.py:87` takes `HOPMIN = min(stmt_hopdist)`
   over *all four* statements, so an SCM enters stratum 1 only if its entire knowledge set is off
   the query. "0 of 791" counts **individual statements**. The n = 229 is a count of graphs, not
   of misstatements.
2. **Ball.** E1′'s censoring statistic is a minimum over the **whole ball ρ = 1…4**, including
   joint flips. "0 of 791" is **ρ = 1** only.
3. **Event.** E1′'s event is non-censoring of `ρ*_any = min(ρ*_se, ρ*_ident)` — the estimate moving
   by more than z·SE at n = 20,000 **or** a consistent member losing amenability. "0 of 791"'s
   event is **silent bias** adjudicated against the true DAG: O\* changed *and* the new set is
   invalid.

X2 decomposed all 38 non-censored HOPMIN ≥ 1 SCMs at n = ∞, and ✅ I reproduced the decomposition
myself from the shipped archive:

| | n |
|---|---|
| SCMs with HOPMIN ≥ 1, all six cells | **1,375** |
| of which not censored | **38** |
| — attributable to `ρ*_ident` alone (**loud** identification failure) | **37** |
| — a genuine ρ ≥ 2 O\*-change | **1** |
| — `ρ*_se = 1`, i.e. a single far statement moving the estimate | **0** |

**Verdict: consistent.** The 6.5% is almost entirely `ρ*_any` folding in *loud* identification
failure — the case where the analyst is told the effect is no longer identified and cannot report
a number at all. That is not silent bias and it is not what S5 denies. The one residual case is a
ρ ≥ 2 joint error, which S5 as written does not quantify over.

**Two consequences that must not be lost.** First, the resolution is a *definitional* one, so it
vindicates S5's zero and simultaneously proves the campaign has been quoting `ρ*_any` in places
where it means `ρ*_se` — those are different instruments and only one of them is silent. Second,
the resolution rescues S5's *number* while X1 and X2 kill the *sentence around it*: nothing at
distance ≥ 1 moves the estimate under a single reversal, and that is precisely the operator the
abstract now promises is only half the story.

---

## 4. Bugs and false claims, ranked

### 🔴 Fatal / must be fixed before anything is published

**1. "0 of 791" has no provenance whatsoever, and it is still in the live abstract.** ✅ Verified:
the pilot codebase (`~/latent-causal/rho-breakdown-knowledge/code/`) contains no BFS, no
`hop_dist`, no graph-distance function; `791` occurs **0 times** in `results/summary.json`; every
occurrence in the raw JSONs is a substring of a float. The number was ordered retired on
2026-08-14 (`output/2026-08-14_latent-causal-iclr2027-rerank/RANKING.md:270-272, :387, :825`) and
`round2/C-triage.tex:13` still prints it. **This is the third audit to say so.**

**2. `wiki/activities/active-claims.md` records a live falsifier as not firing, on the strength of
that number, four days before the b-LOAD CIKM camera-ready (23/08).** ✅ Verified verbatim in the
`b-load-local-knowledge-loop` entry: *"Our own pilot measures 0/791 (rho-breakdown, 2026-07-17), so
this falsifier is currently NOT firing."* It fires for the spurious class (X1: 142 of 3,639 at
hop ≥ 1) and for joint reversals (X2: 70 of 10,615). **The entry cannot be repaired by swapping the
denominator** — the truth value flips with the distance convention and with ρ. It needs new wording,
and per CLAUDE.md Operating Rule #1 that is Zé's to write, not an agent's. **No agent edited it.**

**3. X2's headline "The falsifier fired: S5 as printed is REFUTED" inverts its own
pre-registration.** ✅ Verified at `x2/PREREG.md` §4.1: *"**REFUTED** iff `N_set(r ≥ 1, ρ = 1) ≥ 1`.
One counterexample is enough."* Measured: **N_set = 0** in all six archive cells (D_am = 12,799),
in all five ARM E ensembles, and in the exhaustive p ≤ 5 enumeration. The ρ≥2 refutation that the
headline reports was **pre-disclosed in the same PREREG as already known** (disclosure P5, citing
the 2026-08-14 seed-508 counterexample): *"S5-without-a-ρ-qualifier is already false. The ρ≥2 leg is
confirmatory replication + rate measurement."* `RESULTS.md` never mentions P5, seed 508, or the
pre-disclosure. **The pre-registered falsifier returned zero and the report announces that it
fired.**

### 🟠 Major

**4. X2's S5 denominator pools a structurally inert operator.** ✅ Verified: D_s5 = 54,079 =
**28,822 `add` + 25,257 `rep`**. The `add` sub-operator never moves O\* from an amenable base — 0
events at dmin = 0 across all five ensembles (D_am 2,858 / 1,920 / 328 / 101 / 36) and **0 of
498,792 exhaustively at p ≤ 5, at every distance including 0** — while `rep` fires at
3,355/7,305 near. Including trials that cannot produce the event is the exact error class that
retired the 791. Corrected denominator: **25,257**. Under X2's own M7 thresholds no ensemble then
reaches SUPPORTED (`large` falls 23,056 → 9,930), so the honest grade is SUPPORTED-WEAK.

**5. X1's headline "Locality is a property of the reversal operator, not of knowledge" is
contradicted by X1's own hop table.** ✅ Recomputed from `x1_analysis_*.json`:

| ensemble | S-uni hop 0 | hop ≥ 1 | hop ≥ 2 | hop ≥ 3 | R hop ≥ 1 |
|---|---|---|---|---|---|
| `original` | 167/4,123 | 142/3,639 | **9/419** | **0/44** | 0/2,560 |
| `licensed` | 210/3,663 | 214/4,102 | **2/594** | **0/104** | 0/3,540 |
| `large` | 15/874 | 17/3,403 | **2/1,906** | **0/783** | 0/2,977 |

Pooled: hop ≥ 2 = **13 of 2,919** (0.0045), hop ≥ 3 = **0 of 931**. The materiality ratio of 0.963
that carries the headline is hop **exactly 1**; by hop 2 it is 0.098, below X1's own registered
0.10 bar. Locality does not die — **its radius goes from 0 to 1.**

**6. X1's stated mechanism for the 3.6× is wrong, and the correct one is a better result.**
✅ Recomputed at ρ=1, reachable:

| | O\* moved | P(invalid \| moved) |
|---|---|---|
| `original` R | 688/4,741 = 0.1451 | **1.000** |
| `original` S-uni | 1,140/7,762 = 0.1469 | **0.271** |
| `licensed` R / S | 0.1622 / 0.1700 | **1.000** / 0.321 |
| `large` R / S | 0.0343 / **0.0963** (2.8× more) | **1.000** / 0.078 |

The spurious class disturbs O\* **as often or more often**; the entire ratio is the conditional
factor. 1.012 × 0.271 = 0.274, reproducing the headline exactly. **Every consistent reversal that
moves O\* moves it out of validity — exactly 1.000 on all three grids.** That is the citable
finding and it appears in neither `RESULTS.md` nor the analysis JSON.

**7. X1's "3.6× less damaging" is ρ=1-only and reverses by ρ=4.** ✅ I ran the sweep myself on the
raw data:

| ρ | `original` RR | `licensed` | `large` |
|---|---|---|---|
| 1 | 0.274 | 0.337 | 0.218 |
| 2 | 0.423 | 0.474 | 0.360 |
| 3 | 0.681 | 0.644 | 0.552 |
| **4** | **1.119** (287/1,999 vs 72/561) | 0.805 | 0.793 |

On the primary grid the point estimate is **above 1** at ρ=4 — the registered INDISTINGUISHABLE
branch, which `x1/AUDIT.md` A1 states verbatim "is exactly A26's revive_if condition". `RESULTS.md`
contains no ρ sweep of P2 and concludes A26's revive_if does not fire.

**8. X2's "the radius is sharp — every event involves a misstatement at distance exactly 1" is
contradicted by X2's own primary arm and its own exhaustive lemma.** ✅ ARM E records identification
losses at dmin = 2 (licensed 3/6,485, large 3/17,010, k8 1/306), and `lemma_p5.json` records
**N_ident_rev = 360 of D_am_rev = 3,600 at dmin = 2** — exhaustively, with no sampling error.
X2's own §6 reports r\*_ident = 2, contradicting its own §0 and §4.4. **The radius is 1 for
O\*-movement and ≥ 2 for identification status**; those are two claims and the report merges them.

**9. X1's C_S1 band was applied to a different predicate than the one it was registered for.**
`PREREG.md:288` registers C_S1 ∈ [0.40, 0.95] for "cycle ∨ new-v-structure"; AUDIT M4 replaced the
predicate with intrinsic Dor–Tarsi extendability without re-registering a band; `RESULTS.md:141`
then applies the old band to the new object and declares the "< 0.40 surprise" branch fired.
⚠️ I did not recompute the alternative predicate, but the addendum JSON reports 0.5441 for it —
inside the old band, which would have read CONFIRMED. **X1's surviving claim that the spurious class
is "less catchable" is therefore not established** and should be dropped.

### 🟡 Minor but real

**10. X1's `RESULTS.md` §9.1 and §9.3 rest on a conflation.** ✅ Verified: ARM D (withdrawal)
moves O\* **to a different set in 0 of 7,528 / 7,577 / 3,989** members; the 0.2537 the report
quotes as "withdrawal moves O\* " is the rate at which it makes O\* **undefined** (non-amenable),
counted as a move because `None != frozenset(...)`. The report's claim that an `ostar_changed`
contrast "would not be clean" is therefore wrong, and the clean channel is the stronger result.

**11. E1′'s published locality table selectively omits its weakest cell, and its "13–27×
separation" is contradicted by the table above it.** `k8/K4` (dist-1 censoring 0.7778, n = 27) is
in the JSON and not in `RESULTS.md`; the four printed rows give separations of 26.7×, 13.4×, 59.7×
and **3.0×** (`large`). Two of four fall outside the stated range.

**12. E1′ merges a disconnected-component stratum into a "dist ≥ 2" headline.** 57 of the 331 SCMs
behind `large`'s "1.000 at dist ≥ 2" have all statements in a component not connected to {X,Y}
(`hop_dist`'s sentinel `1<<20 = 1048576`, `run_arm1.py:49`, never converted to infinity). That is
17% of the denominator at infinite distance, not two hops.

**13. `analyse_arm1.py:213` silently drops any stratum with n < 10 and prints "—", which reads as
n = 0.** 6 `licensed`, 5 `original` and 2 `large` SCMs at distance ≥ 2 vanish from the published
table — the most informative cells for the locality law.

**14. X1's G4 "dead switch" gate cannot fail:** `run_x1.py:285` writes `max_abs_dev=0.0` and
`ostar_changed=0` as literal constants, and `analyse_x1.py:122` reads those literals to decide
PASS. Reported as "PASS, n = 30,000", which verified nothing.

**15. X1's `P_consistent_R` divides a restricted numerator by an unrestricted denominator**
(`analyse_x1.py:268`): 4,741/8,000 = 0.5926 where the conditional-on-reachable value is
4,741/7,528 = 0.6298. It propagates into the free/class decomposition printed in §5.

**16. X2's registered gate G6 (the seed-508 cross-codebase anchor for the ρ≥2 leg) was never run
and its omission is not disclosed.** ✅ Verified: `gates` keys in `x2_analysis.json` are
`[machinery_tests, G3_*, G10_*, G4a_*, G4b_*, G7, G9_*, G11]`. `PREREG.md:539` says the ρ≥2 leg is
"reported as **unanchored**" if it does not reproduce; §10's deviations list omits it.

**17. X2's `lemma_p6.json` cannot be regenerated by any committed command.** `run_all3.sh:13`
passes no MPDAG cap; the shipped file records `mpdag_cap: 400` and `n_cpdags_skipped: 37`. The
skipped CPDAGs are the densest, so the p = 6 row is biased **against** finding a counterexample.

**18. X2's pooled 0.66% [0.51%, 0.83%] is composition-driven.** ✅ Per cell the ρ≥2 all-far rate
runs from **1/5,677** (`large`, 0.018%) to **31/1,215** (`k8/K8`, 2.55%) — a 145× spread, with
`large` supplying 53% of the pooled denominator and 1.4% of the events. The pooled CI estimates a
mixing weight, not a population quantity. Report the count and the existence; do not paste the CI.

---

## 5. Where the auditors disagreed with the implementers — adjudicated

I do not average these. Each is decided on the evidence, and where I could recompute, I did.

### X1

| dispute | my adjudication |
|---|---|
| **A3 (FLAWED) vs implementer on "locality is a property of the reversal operator"** | **A3 is right, decisively.** I recomputed the hop table: hop ≥ 2 is 13/2,919 pooled and hop ≥ 3 is 0/931. The counterexamples live at hop exactly 1. The correct statement is that the assertion operator extends the radius from 0 to 1, not that locality dies. |
| **A3 vs implementer on the ρ-sweep** | **A3 is right; I reproduced it independently.** RR is monotone in ρ and reaches 1.119 on the primary grid at ρ=4. The "3.6×" is a ρ=1 statement and must be labelled one. |
| **A3 vs implementer on the mechanism** | **A3 is right, and this is the most valuable single finding in the exercise.** P(invalid \| moved) = **1.000** for reversals on all three ensembles and 0.271/0.321/0.078 for assertions. The spurious class is not inert; its disturbances usually stay valid. |
| **A1 vs implementer on the T10b "withdrawal confound"** | **A1 is right.** ARM D moves O\* to a different set 0 times; the 0.2537 is loss of amenability scored as a move by `None != frozenset`. §9.1 must be retracted, and the machinery-test count is 14/15 once the outcomes are separated. |
| **A2 vs implementer on the C_S1 band** | **A2 is right.** A band registered for one predicate was applied to a different one and the verdict flips under the substitution. "Less catchable" is not established and should be dropped. |
| **A2 vs implementer on the mpdag_valid split (M5)** | **A2 is right that M5 was applied to P2 only** and that `RESULTS.md`'s blanket "all 16 mitigations were implemented" is inaccurate. ⚠️ I did not recompute A2's well-posed rates (0.0200 vs 0.0390), so I flag the magnitudes as unverified — but the direction is corroborated by my hop table, where the far damage collapses beyond hop 1 regardless of stratification. |
| **A2 vs implementer on the moved P2 bar** | **The implementer is fine here.** The AUDIT changed the rule before implementation and before any outcome existed, and A2 verified the PREREG rule as literally written gives the same REFUTED verdict. Disclose both; no claim is at risk. |
| **A1's "four undisclosed scripts on betelgeuse"** | ⚠️ **Unresolved and I flag it.** `refute_x1.py`, `verify_x1.py`, `mech_x1.py`, `harm_x1.py` exist in `x1/audit/refuter/` — they are the *auditor's* files, not withheld implementer output. A1 appears to have mistaken the auditor's working directory for the implementer's. I could not establish that the implementer withheld anything. |

**X1 overall: the measurements are sound and I reproduced every one I checked. The interpretation
is wrong in three specific, fixable places (radius, ρ-scope, mechanism). A3's FLAWED is the right
verdict on the write-up; A1's and A2's SOUND_WITH_CAVEATS is the right verdict on the code.**

### X2

| dispute | my adjudication |
|---|---|
| **A2 (FLAWED) vs implementer on "the falsifier fired"** | **A2 is right, and this is the most serious finding against either report.** ✅ I read the PREREG myself: REFUTED iff `N_set(r ≥ 1, ρ = 1) ≥ 1`; measured 0 everywhere; and the ρ≥2 leg was pre-disclosed as already known in the same document. A report may not announce as a discovery the thing its own pre-registration listed as prior knowledge. |
| **A1 and A2 vs implementer on the `add` operator in D_s5** | **Both are right and they agree.** ✅ Verified: 28,822 of 54,079 trials come from an operator with 0 events in 498,792 exhaustive trials at every distance. The corrected denominator is 25,257 and no ensemble reaches SUPPORTED. |
| **A3 vs implementer on "the radius is sharp = 1"** | **A3 is right.** ✅ ARM E has 7 identification losses at dmin = 2 and the exhaustive p=5 lemma has 360 of 3,600. X2's own §6 says r\*_ident = 2. The report contradicts itself. |
| **A3 vs A2 on whether S5 is refuted at all** | **A3 is closer.** Read per statement — which is how the sentence is written ("no misstatement… moved it") — S5's zero is *confirmed and strengthened by two to three orders of magnitude*. What dies is the clause before it, "the optimal adjustment set is a local functional of the graph", plus the set-level reading at ρ≥2. This distinction matters operationally: the b-LOAD register entry needs **narrowing**, not reversal. |
| **A1 vs implementer on the 6,726 provenance** | **A1 is right that no committed script produces it, and the implementer is right about the number.** ✅ I reproduced it exactly on the third attempt: it is the **ρ≥1**, finite-distance, closest-misstatement ≥ 2 stratum (D_all 10,029 / D_con 6,726 / D_am 6,726, N_set = 0, N_ident = 0). My first two readings (ρ≥2 finite = 2,689; ρ≥2 with sentinel = 5,448) both gave the same zero but different denominators. **A table whose denominator takes three attempts to reproduce is a table that needs its script committed.** |
| **A1's "the five ensembles share seeds"** | ⚠️ **Plausible and unverified by me.** `run_scan.py` uses one JOB_SEED per group; A1 reports 9,153 of k8's 12,000 SCMs are shared with `licensed`. I did not check. If true, "five ensembles" is closer to three independent cells and every pooled figure double-counts. **Flagged as a live question.** |
| **A3 vs implementer on the E_gain events** | **A3 is right.** All 29,364 fire only from a non-amenable base, where no ρ-ball is being enumerated. They should not carry S9's refutation; the 7 in-scope events at dmin = 2 and the exhaustive 360/3,600 should. |
| **A1 vs A3 on the cluster bound** | **A3 is right in principle.** ✅ The 6,726 trials come from **2,746** (ensemble, arm, seed) clusters, so the graph-level rule-of-three bound is ≈ 3/2,746 = 1.1 × 10⁻³, not the trial-level 4.5 × 10⁻⁴. The reported bound is ~2.5× optimistic. |

**X2 overall: A2's FLAWED is correct on the framing and A1's/A3's SOUND_WITH_CAVEATS is correct on
the arithmetic — every number I recomputed reproduced exactly, including the 70, the 0/136,200 and
the full contradiction decomposition. The defect is that the report's title asserts the opposite of
its own pre-registered decision rule.**

### The one place X1 and X2 genuinely disagree with each other

They do not disagree on any measurement. They disagree by **omission**, and it is worth stating
plainly: X1 concludes "locality is a property of the reversal operator, not of knowledge" while X2
concludes "the radius is sharp at 1". ✅ Both are wrong in the same direction and the reconciliation
is arithmetic: X1's own hop ≥ 3 cell is 0 of 931 and X2's own dmin ≥ 2 cell is 0 of 6,726. **The
two runs jointly establish a radius-1 result that neither claims.**

---

## 6. The highest-value next experiment

**Run X0 — the two-line lemma behind S2 — and fold the radius-1 result into it as a proof
obligation. Do not run X3.**

The reasoning, against `~/rho-breakdown-iclr2027/abstract/EXPERIMENTS.md`:

- **X1 and X2 are done and they converged on a theorem-shaped object.** Both runs independently
  produce a clean zero beyond one hop (X1: 0/931 at hop ≥ 3; X2: 0/6,726 at dmin ≥ 2; and the
  exhaustive p ≤ 5 enumeration closes dmin = 1 completely at 0/136,200). ✅ X2's own mechanism
  numbers point at why: O\* is a functional of `cn(X,Y)` and its parents, and in an amenable MPDAG
  nearly every edge touching `cn` is already compelled. **That is a lemma, not a measurement**, and
  it is the exact proof obligation reviewer R1 named ("no distant statement can enter the
  neighbourhood through Meek closure"). It is now the cheapest remaining win because the empirical
  boundary conditions are pinned: it is false at radius 0, true at radius ≥ 2 for O\*-movement, and
  false at radius ≥ 2 for identification. A proof that has its counterexamples already enumerated
  is a much easier proof.

- **It beats X3 (enlarge the ball) decisively.** X3 defends S10 (90.3%), which no reviewer attacked
  and which is a ⚠️ not a 🔴 in the ranking. It requires a new generator and a full rerun ("medium-
  high" cost), and ✅ E1′ already shows the number moves to 88.4% at ρ ≤ 4 — a 1.9-point change that
  buys nothing. Meanwhile S8 and S9 are **dead** and S5's first clause is **dead**; spending the
  next cycle defending an unattacked sentence while three attacked ones are unrepaired is the CIKM
  tutorial failure mode written down in CLAUDE.md.

- **It beats re-running X1 or X2 larger.** Neither is denominator-limited on its core finding. X2's
  zero is now exhaustive at p ≤ 5; X1's radius result is 0/931 with the collapse visible at hop 2.
  More SCMs would tighten intervals nobody disputes and would not touch the three framing defects,
  which are rewrite work, not compute work.

- **The concrete deliverable**, in order: (i) state and prove the radius-1 lemma for the reversal
  operator, using the `cn`-locality mechanism and the enumerated hop-1 counterexamples as the
  sharpness witness; (ii) state the *separate* radius-2 result for identification status, whose
  witness is the exhaustive 360/3,600 at p = 5; (iii) implement the pruned enumerator that S8
  describes and benchmark it, since S8 currently describes an algorithm that has never been written.
  ⚠️ If (i) resists proof, the honest fallback is already measured and the abstract can carry the
  heuristic with its Q₀ = 0.12–0.15 attached.

**The one thing that must happen before any of it:** the b-LOAD `active-claims.md` entry and the
"0 of 791" line, both flagged in §4. The camera-ready is 2026-08-23.

---

## 7. Limitations of this whole exercise

Written so that no gap here can later be read as coverage.

1. **Nothing here is a proof.** X2's zero is exhaustive only at p ≤ 5. Above that everything is
   sampling, and at ρ ≥ 2 the property is *false*, so the theorem R1 asked for is not merely
   unproved — it is refuted as stated and must be restated with a radius before it can be proved.

2. **S7 is now implemented but not validated against b-LOAD.** X1 is not a b-LOAD replication: no
   `mb_by_mb`, no CI tests, no finite samples, no R. X1's own gate G8 measures that the two closure
   conventions (close-after-each vs stamp-all-then-close) disagree on ~28% of members, so
   "b-LOAD's semantics" is itself convention-dependent and every X1 number is on one convention.

3. **X1 never ran the discovery arm.** All X1 graphs use the oracle CPDAG `C = dag_to_cpdag(D)`.
   An analyst running PC or GES would *see* that the asserted pair is absent from the estimated
   skeleton — that is the whole testability asymmetry between the two error classes, and the honest
   answer on it is "not tested".

4. **The spurious/reversal mixture is not covered.** X1 measures two pure strata; b-LOAD's own
   sampler produces a mixture. The pure strata bracket the mixture only if the error types do not
   interact inside Meek closure, which is untested. ✅ X1's own pooled S-loc row (0.1010) inverts its
   pure-stratum verdict (0.0310), so the interaction is not obviously negligible.

5. **"Graph-distance of a statement" remains undefined for an edge, and the answer flips with the
   convention.** ✅ On the same archive, N_set(dmin ≥ 1) = 0 while N_set(dmax ≥ 1) = every event.
   Every sentence in §2 uses the **minimum-endpoint** convention. The defensible phrasing drops
   "graph-distance" and says "not incident to X or Y".

6. **X2's distance measurements assume the skeleton is correct**, and X1's operator changes the
   skeleton. The two experiments therefore do not compose: X1's hop distances are computed on a
   graph whose skeleton already contains the spurious edge. A joint treatment does not exist.

7. **Unverified items, listed so they are not mistaken for checked ones.** ⚠️ A1's seed-sharing
   claim across X2 ensembles; ⚠️ A2's mpdag_valid-stratified X1 rates; ⚠️ A1's design-effect
   recomputations; ⚠️ X2's `deep`/`deep8` ARM C ensembles, which I did not touch; ⚠️ the
   reproducibility of X2's scan, which the implementer discloses is non-deterministic
   (`imap_unordered` with an early break, drift of 3–6 trials in 18,000).

8. **Scope, unchanged from the pilot.** Linear-Gaussian iSCM, Erdős–Rényi DAGs, single-node X and
   Y, population Σ, orientation statements only, expected degree ≤ 6, no latent confounding, no
   real data, no semi-synthetic arm. The abstract's `\todo{Semi-synthetic arm.}` is still open and
   neither X1 nor X2 touched it.

9. **Two experiments in this vault are named "E1" and dated the same day.**
   `output/2026-08-19_e1-uninformative-variates/` is an unrelated EEG experiment. Anything citing
   "E1" must disambiguate.

10. **I did not edit the abstract, the wiki, the claims register, or any prior experiment's
    directory, and I made no commit.** All rewrites in §2 are proposals.

---

## Provenance of the numbers in this document

| verified how | items |
|---|---|
| ✅ recomputed by me from raw data on betelgeuse | X1 ρ-sweep (all 3 ensembles × 4 ρ); X1 mechanism decomposition (moved / P(invalid\|moved)); ARM D channel separation; X2 sharpness stratum (3 conventions); pilot distance-function absence; `791` absence |
| ✅ recomputed by me from shipped JSONs | X1 hop table (3 ensembles × 2 arms); X2's 70 counterexamples and their ρ split; X2 ρ=1 far zero; D_s5 add/rep split; `add` inertness at dmin=0; exhaustive lemma arithmetic; contradiction decomposition (1,375 / 38 / 37 / 1 / 0); E1′ per-ensemble S3/S4 spread |
| ✅ read verbatim | `C-triage.tex`; `x2/PREREG.md` §4.1 and disclosure P5; `active-claims.md` falsifier #3 |
| ⚠️ taken from an audit, not independently checked | items listed in Limitation 7 |
