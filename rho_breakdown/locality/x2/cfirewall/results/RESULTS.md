# C-FIREWALL — exhaustive decision at p ≤ 5, samples at p = 6, 7, 8

**Verdict: NOT REFUTED. 0 violations in 70,375,325 trials / 103,105,398 (trial × D-in-class)
validity checks. The census at p ≤ 5 is complete, so no counterexample exists at p ≤ 5.**

Date 2026-08-23. Code `../code/`, raw counts `SUMMARY.json` and the per-run JSONs beside it.

---

## 1. What was decided, and in which form

The theorem was tested in its **STRONG** form — the one the brief flagged as possibly not
covered by `move2.py`, which draws a single D:

> for **every** D ∈ [G₀], `O*(X,Y,H)` is a valid adjustment set for (X,Y) in D.

Note the restatement this licenses: *"valid in every D of the class"* is precisely *"valid
relative to the MPDAG G₀"*. So the strong C-FIREWALL says the graph **H** — whose skeleton
**no** member of [G₀] possesses — yields a set that is adjustment-valid for the **whole**
equivalence class of the truth.

Sweep, per the brief, with nothing pruned:

```
every labelled DAG on p nodes -> CPDAG (deduped)
  -> every reachable MPDAG G0 = M(C,K)     (run_lemma.reachable_mpdags: exactly {M(C,K)})
    -> every ordered query (X,Y) with G0 amenable
      -> every NON-ADJACENT pair {a,b} of skeleton(G0), BOTH orientations
        -> bk_assert -> (conflict? cycle? Dor-Tarsi?) -> if H amenable:
             O*(X,Y,H) tested against EVERY D in [G0]
```

`[G₀] = consistent_dag_extensions(G₀, v_structures(C))`, verified by brute force (test T2) to
equal `{D ∈ MEC(C) : D agrees with every directed edge of G₀}`.

## 2. Counts

| arm | p | mode | CPDAGs / draws | MPDAGs | trials | aborts | reports | moved | (trial×D) checks | **violations** |
|---|---|---|---|---|---|---|---|---|---|---|
| k=1 | 3 | **CENSUS** | 11 | 50 | 78 | 6 | 72 | 18 | 78 | **0** |
| k=1 | 4 | **CENSUS** | 185 | 1,601 | 10,996 | 720 | 10,276 | 2,232 | 13,536 | **0** |
| k=1 | 5 | **CENSUS** | 8,782 | 116,992 | 2,161,055 | 107,340 | 2,053,715 | 402,840 | 3,267,660 | **0** |
| k=1 | 6 | sample | 400,000 | 400,000 | 22,851,273 | 313,900 | 22,537,373 | 4,400,761 | 31,515,427 | **0** |
| k=1 | 7 | sample | 120,000 | 119,808 | 16,435,189 | 183,098 | 16,252,091 | 2,870,115 | 24,827,139 | **0** |
| k=1 | 8 | sample | 40,000 | 39,983 | 10,839,188 | 103,391 | 10,735,797 | 1,774,275 | 17,684,624 | **0** |
| k=2 | 4 | **CENSUS** | 185 | 1,601 | 25,032 | 2,040 | 22,992 | 8,324 | 29,088 | **0** |
| k=3 | 4 | **CENSUS** | 185 | 1,601 | 39,468 | 2,520 | 36,948 | 16,940 | 44,460 | **0** |
| k=2 | 5 | **CENSUS** | 8,782 | 116,992 | 10,091,230 | 737,215 | 9,354,015 | 3,194,385 | 14,364,835 | **0** |
| k=2 | 6 | sample | 39,659 | 39,659 | 3,666,867 | 66,599 | 3,600,268 | 1,197,631 | 5,165,549 | **0** |
| k=3 | 6 | sample | 39,684 | 39,684 | 2,466,240 | 53,082 | 2,413,158 | 1,078,765 | 3,446,856 | **0** |
| k=4 | 7 | sample | 24,939 | 24,939 | 1,788,709 | 37,207 | 1,751,502 | 876,268 | 2,746,146 | **0** |
| | | | | **TOTAL** | **70,375,325** | | **51,589,324** (k=1) | | **103,105,398** | **0** |

- **p values COMPLETED AS A FULL CENSUS: p = 3, 4, 5** (k = 1; and k = 2 at p = 4, 5; k = 3 at p = 4).
  At p ≤ 5 the enumeration is exhaustive over CPDAGs, over MPDAGs, over queries, over
  non-adjacent pairs, over orientations, **and over [G₀]**. A counterexample at p ≤ 5 does not exist.
- **p = 6, 7, 8 are SAMPLES**, labelled as such. Every sampled G₀ is then swept exhaustively
  (all queries × all non-adjacent pairs × both orientations × all D in [G₀]).
- The predicate is **not vacuous**: 9,450,241 of 51,589,324 single-statement reports (18.3%)
  genuinely **moved** `O*`, and 708,455 trials aborted visibly.
- Rule-of-three 95% upper bound on the violation rate given 0/70,375,325: **< 4.3 × 10⁻⁸**.

### p = 6 as a census: measured NOT tractable

A uniform probe of 40 p=6 CPDAGs cost **5.40 s/CPDAG** ⇒ ≈ **160 core-hours** for the
1,067,825 CPDAGs on 6 nodes. The cost is concentrated in dense chain components and is largely
**wasted**: the K₆ CPDAG alone reaches **130,023** MPDAGs (198 s) and contributes **zero**
trials, because a complete skeleton has no non-adjacent pair. Hence the sample.

## 3. The stronger "all D" version — it holds, and it coincides with the one-D version

This is the point the brief singled out. Two measurements:

1. **The strengthening is genuinely exercised.** 14,483,847 of 51,589,324 single-statement
   reports (28.1%) came from an MPDAG with **|[G₀]| ≥ 2** (max |[G₀]| observed = 156);
   2,954,238 of those also moved `O*`. This is not a class-size-1 census in disguise.
2. **Validity never splits across the class.** `violating_keys_someD = 0` trivially (there are
   no violations at all) — but the same instrumentation on the **control arm**, where 3,005,488
   reported sets *are* invalid, records `ctrl_violating_someD = 0` as well. **Every** invalid
   set was invalid in **every** D of its class simultaneously.

So: the stronger version **holds**, and there is **no** case where it holds for a sampled D and
fails for another. `move2.py`'s single-D draw was not understating anything — but that was an
open possibility until now, and it is now measured rather than assumed.

## 4. The detector has power (control arm, in-census)

Same predicate, same [G₀], on graphs where a coherence check **failed** and the edge was
stamped anyway:

| | trials | reports | moved | invalid | invalid \| moved |
|---|---|---|---|---|---|
| control (p = 3…8) | 13,847,787 | 13,791,476 | 4,739,226 | **3,005,488 (21.8%)** | **63.4%** |

100% of control violations were *moved* reports. A 0/70M treatment result is therefore not the
artefact of a checker that never returns False.

## 5. 🔴 The offered PROOF HANDLE is FALSE — both halves. A replacement is given.

The synthesiser proposed: *"adding a→b can only GROW cn and forb. Show that every parent the
addition contributes to cn_H is already in forb_H."* Measured on the **p ≤ 5 census**
(`structure_p5.json`, n = 2,053,715 coherent amenable trials):

| claim | verdict | counter-count |
|---|---|---|
| cn_H ⊇ cn_G₀ | **FALSE** | 45,720 strict **subsets**, 6,360 **incomparable** |
| forb_H ⊇ forb_G₀ | **FALSE** | 160 strict **subsets** |
| pa_H(cn_H) ⊇ pa_G₀(cn_G₀) | **FALSE** | 4,360 subsets, 2,160 incomparable |
| new parents ⊆ forb_H ("total_loss = 0") | **FALSE** | of 513,560 new parents, **282,480 (55.0%) escape forb_H** and land in `O*_H` |

Mechanism for the shrink: Meek closure inside H can orient `v_{i+1} → v_i` on an edge that was
*undirected* in G₀, **destroying** a possibly-causal path. cn is not monotone in the edge set.
The same pattern reproduces at p = 6 (`structure_p6_sample.json`, n = 6,741,208).

**The surviving invariant — offer this instead (call it F1):**

> **F1.  O*(X,Y,H) ∩ forb(X,Y,G₀) = ∅.**

**0 violations** over the entire p ≤ 5 census (2,064,063 trials) and the p = 6 sample
(6,741,208 trials). F1 is exactly what discharges **clause (1)** of the generalized adjustment
criterion for every D ∈ [G₀] at once, because `forb_D(X,Y) ⊆ forb_{G₀}(X,Y)` for D ∈ [G₀]
(a causal path in D is possibly-causal in G₀, and `de_D ⊆ poss_de_{G₀}`). What a proof still
owes is **clause (2)**: that `O*_H` blocks every proper non-causal path X…Y in D.

F1 is not a trivial re-reading of a containment: `O*_H ⊊ O*_{G₀}` in 123,000 p=5 trials and the
two sets are *incomparable* in 24,480, so no sub/superset argument reaches it.

## 6. 🔴 The three coherence checks collapse to ONE — the hypothesis is over-stated

Tested **pointwise** (`checks_agree.py`, `checks_agree.txt`) over the whole p ≤ 5 census plus
p=6 samples at k = 1, 2, 3 and p=7 at k = 4 — **1,066,670 statement evaluations**:

- `conflict` **never fires** (0 of 1,066,670). A required edge on a non-adjacent pair has no
  existing orientation to contradict; with k ≥ 2 the pairs are distinct, so no self-conflict either.
- `has_directed_cycle(H)` ⟺ `not pdag_extendable(H)`, with **zero disagreements**.

So on this operator the theorem's hypothesis *"passes all three intrinsic coherence checks"* is
equivalent to *"H has no directed cycle"*. The theorem is therefore **true under a weaker
hypothesis than it states**, and the Dor–Tarsi test (x1's M4 amendment) buys nothing here. State
the hypothesis as acyclicity, or the reviewer will ask why three tests are named for one job.

## 7. The abort is always a REFUSAL, never a "zero effect" claim

`abort_nopath = 0` in **all 70,375,325 trials**. Every visible abort is the amenability failure
(`abort_unamen`), never a collapse of the possibly-causal path set. This matters practically:
a `npaths = 0` graph would make the tool *report an effect of zero* — an unsound **claim** —
rather than refuse. That branch is never entered.

## 8. Bonus: the firewall survives k ≥ 2 spurious edges — which does NOT follow by induction

C-FIREWALL does not self-induct: after one assertion H has a skeleton no DAG of the class has,
so H is not an MPDAG of any CPDAG containing D and the theorem cannot be re-applied. Measured
independently anyway (`kfirewall.py`), on ordered sequences of k distinct non-adjacent pairs
with the coherence checks applied to the final graph, as a tool would:
**18,077,546 trials, 6,372,313 of them moved, 0 violations**, including full censuses at
(p=4, k=2), (p=4, k=3), (p=5, k=2).

## 9. Violations

**None.** No witness to report at any p, any k, any arm.

## 10. Machinery tests (`test_cfirewall.py`, all PASS at p = 3 and p = 4)

- **T1** `bk_assert(G₀, a→b)` on a non-adjacent pair == add-the-edge-then-Meek-close (the
  prefix re-stamp is a no-op for a single statement), and skeleton(H) ⊋ skeleton(G₀) always —
  so skeleton(H) ≠ skeleton(D) and [H] ∩ [G₀] = ∅, the gap the theorem must cross.
- **T2** `[G₀]` == brute-forced `{D ∈ MEC(C) : D agrees with G₀}`.
- **T3** `[G₀]` non-empty, and G₀ == the common orientation graph of `[G₀]` (so
  `reachable_mpdags` really enumerates MPDAGs, not arbitrary PDAGs).
- **T4** `ostar_and_paths` == `adjust.optimal_adjustment_set` on every reachable MPDAG.
- **T5** baseline: `O*(X,Y,G₀)` is valid in every D ∈ [G₀] (HPM class invariance) — the checker
  returns True where it must.
- **T6** the checker returns False where it must: 6,640 / 11,892 in the p=4 incoherent control.

## 11. Files

```
code/cfirewall.py       census + sample driver (single statement)
code/kfirewall.py       the k >= 2 extension
code/structure.py       the proof-handle diagnostics (cn / forb / pa / O* / F1)
code/checks_agree.py    pointwise agreement of the three coherence checks
code/test_cfirewall.py  T1-T6 machinery tests
code/summarise.py       rolls results/*.json into results/SUMMARY.json
results/cfirewall_p{3,4,5}.json                     full census, k=1
results/cfirewall_p{6,7,8}_sample.json              samples, k=1
results/kfirewall_p4_k2.json  p4_k3  p5_k2  p6_k2  p6_k3  p7_k4
results/structure_p{3,4,5}.json  structure_p6_sample.json
results/checks_agree.txt        results/SUMMARY.json
```

## 12. What would still refute it

The census closes p ≤ 5 permanently. Any counterexample lives at p ≥ 6, and the p=6/7/8 samples
are uniform over (density, CPDAG, MPDAG) rather than adversarial. The unswept regions are:
(a) p = 6 exhaustively (≈160 core-hours, buys the next closed size), and (b) MPDAGs with a chain
component larger than the |U| ≤ 14–16 cap used to keep `[G₀]` enumerable — 209 of 560,000 draws
at k=1 were skipped for that reason and are recorded as `draw_skipped_bigclass`.
