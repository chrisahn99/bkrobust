# X1 RESULTS — the spurious required-edge arm

**Ran 2026-08-19 on betelgeuse**, 12 workers, `~/e1venv/bin/python`.
Machinery tests **15/15**. Integrity gates **12/12**. Total compute: **~3 min** for the three
sweeps (26.9 s + 75.6 s + 64.7 s), plus gates and analysis.

Pre-registration: `PREREG.md`. Binding design audit: `AUDIT.md` (verdict FLAWED, 16 mitigations).
Every mitigation marked *fatal* or *major* was implemented. Four objections are recorded in §9.

---

## 0. Headline — the falsifier fired on the load-bearing prediction

> **X1's registered claim (PREREG §3): the spurious class is "less catchable and at least as
> damaging" than the reversal class. The first half holds. The second half is REFUTED, on all
> three ensembles, by a margin far outside the audit's free-bar break-even.**

Primary contrast (AUDIT M1: conditional risk ratio, paired SCM-level cluster bootstrap,
10 000 resamples, ρ = 1, reachable statements, ensemble `original`, n = 2 000 analysed SCMs):

```
RR = P(silent bias | S-uni, b-LOAD semantics) / P(silent bias | consistent, R)
   = 0.0398 / 0.1451  =  0.274   [0.242, 0.308]        ->  REFUTED
```

Registered bar: SUPPORTED iff `RR ≥ 1.25` and `CI_lo > 1`; **REFUTED iff `RR ≤ 0.80` and
`CI_hi < 1`**. The audit's free-bar disclosure (A1/M1.3) says the `original` grid only reaches
REFUTED at a conditional ratio of **0.540**. Measured: **0.274**. The refutation is not free.

**And the result that matters more than the headline** — the reversal operator *cannot* change the
skeleton, the assertion operator can, and locality dies with the operator:

| ρ = 1, `original`, hop-distance ≥ 1 from the query, reachable | silent bias | rate |
|---|---|---|
| **ARM R** (reversal — the campaign's only operator to date) | **0 of 1 574** | 0.0000 [0.0000, 0.0024] |
| **ARM S-uni** (spurious required edge, b-LOAD semantics) | **142 of 3 639** | **0.0390 [0.0332, 0.0458]** |

The hop-0 rate for S-uni is 0.0405 [0.0349, 0.0470]. **Materiality ratio 0.963** — a spurious
required edge one or more hops from the query is essentially as damaging as one adjacent to it.
X1-P4 registered SUPPORTED at ≥ 0.10; measured 0.963 with 142 events on a denominator of 3 639.

> **Locality, as the campaign has measured it, is a property of the reversal operator, not a
> property of knowledge.**

---

## 1. What was run, before any verdict (AUDIT A17, M10, M12.3)

| ensemble | analysed SCMs | statements/arm at ρ = 1 | hop ≥ 1 reachable (S-uni) | unreachable | `P(|O0| = 0)` | mean `|O0|` |
|---|---|---|---|---|---|---|
| **`original`** (PRIMARY) | **2 000** ✅ (bar: ≥ 1 800) | 8 000 | **3 639** ✅ (bar: ≥ 2 000) | 238 | 0.6635 [0.6425, 0.6839] | 0.512 |
| `licensed` | 2 000 | 8 000 | 4 102 | 235 | 0.4755 [0.4537, 0.4974] | 1.115 |
| **`large`** (CO-PRIMARY for P4/P5, M15) | 1 200 | 4 800 | 3 403 | 523 | 0.2342 [0.2111, 0.2590] | 1.670 |

`OVERSAMPLE` raised to 6 / 8 / 4 to hit the M10 targets stated in **analysed** SCMs; all three were
hit exactly. Dropped SCMs and their profiles are in `results/x1_analysis_*.json → header.drops`
(`original`: 213 ineligible at mean `|N(C)| = 2.72`, 620 non-amenable, 7 218 pre-filtered).

**Smoke-run disclosure (PREREG §2.5, M16 ii).** A 48-SCM smoke ran on `large` **first**, then on
`original` and `licensed`. It emitted only the analysed-SCM count and wall time; **no outcome was
computed from it** and the three files were deleted before the real sweep.

**Balls per SCM:** `R`, `S-uni`, `S-loc`, **`D` (drop control — see §2)**, `P-null-identity`,
`P-null-append`, `P-frozen`. `S1` is a derived label on the S-arm graphs, not a second run.

---

## 2. 🔴 A confound the audit did not catch, found by the machinery tests before any outcome existed

The registered perturbation law is `K'(F) = (r_t if t ∈ F else (u_t, v_t))` — it **replaces** slot
`t`. Arm R reverses a slot, so all four original pairs are still asserted. Arm S **replaces** one,
so only three are. **Every raw S-vs-R contrast therefore carries a "withdrawn true statement"
effect that has nothing to do with the spurious class.** AUDIT A6 found the special case (P-null)
and prescribed M6; it did not generalise it to the arms.

Found in `test_x1.py::T10`/`T10b` before any damage rate was computed, and answered by adding
**ARM D**: slot `t` removed, nothing asserted in its place. Same SCMs, same slots, same outcomes.

Measured (ρ = 1, `original`, n = 8 000 members):

| | O\* changed | non-amenable (loud) | **silent bias** |
|---|---|---|---|
| **ARM D (withdrawal alone)** | **0 of 8 000** | 0.2388 | **0 of 8 000** |

So the withdrawal channel is **loud or benign, never silently biasing**: dropping a true statement
either leaves `O*` bit-identical or destroys identifiability, and the analyst sees the second.
Consequences, both reported because they cut opposite ways:

- **The P2 contrast on `silent` is clean.** All 309 S-uni silent events at ρ = 1 are attributable to
  the added edge; arm D contributes zero. The headline is not confounded.
- **A P2 contrast on `ostar_changed` would not be clean**, and the campaign's `ostar_changed` and
  `silent` are *not* interchangeable across operators the way they were within arm R.
- `O*(S) ≠ O*(D)` (in set **or** in identifiability) at matched (ρ, flip): **0.2209 [0.2119, 0.2301]**
  (`original`, ρ = 1, n = 8 000). That is the edge-attributable figure; the raw-vs-`G0` figure mixes
  in the withdrawal.

---

## 3. Integrity gates — 12/12, before any verdict (AUDIT A17)

| gate | result | n |
|---|---|---|
| **G1** anchor shared across arms | **STRUCTURAL** — one anchor computed per SCM, passed to every ball. Declared, not claimed as evidence. | — |
| **G2** ARM R reproduces E1′ member-for-member | **PASS**, 0 mismatches on `consistent`/`amenable`/`ostar_changed`/`ostar_valid`, `est` to 1e-12 | **2 000** shared SCMs / 30 000 members (`original`); **1 879** / 28 185 (`licensed`) |
| **G3** S-uni S2 catch rate = 0.000; skeleton-detectability = 1.000 | **PASS** (both definitional) | 30 000 / 8 000 |
| **G3′** S-loc skeleton-detectability — **MEASURED, not a gate** (M7) | 0.7305 [0.7207, 0.7401] | 8 000 |
| **G4** P-frozen dead switch | **PASS** — max abs deviation 0.0, `ostar_changed` 0, `silent` 0 | 30 000 |
| **G5a** P-null-identity (M6) | **PASS** — `G = G0` in 8 000/8 000, `ostar_changed` 0, `silent` 0 | 8 000 |
| **G5b** P-null-append, `|K'| = 5` (M6) | **PASS** — `G = G0` in 2 000/2 000, `ostar_changed` 0 | 2 000 |
| **G6** `d(ρ)` monotone in ρ | **PASS** — 0 violations | 8 000 SCM×arm |
| **G7** re-stamp is a no-op, both scopes | **PASS** — 0 divergences | 60 000 |
| **G8** close-after-each vs stamp-all-then-close | **MEASUREMENT, no bar** — see §6 | 60 000 |
| **G9** no `tau` on the instrument path | **PASS** (source check) | — |
| **G10** sentinel discipline | **PASS** — `unreachable` is a labelled row everywhere, never pooled, never averaged | — |
| **G11** Dor–Tarsi vs brute-force extension oracle, in-run | **PASS** — 0 disagreements | 11 010 (`original`) + 4 767 (`large`) members |
| **G12** unreachable added edge cannot move `O*` | **PASS** — 0 edge-attributable moves | 318 / 351 / 993 |

Machinery tests (`logs/test_x1.log`), **15/15**, all validating the new operator against its
*definition* by brute force, never against a reimplementation: Dor–Tarsi vs enumeration
(0/4 000 disagreements); `bk_assert` on true `K` reproduces `apply_background_knowledge`
(0/744); Meek-consistency trivially false for the spurious class (0/4 616 under **both**
`is_consistent` and MEC-membership brute force); operator adds exactly the asserted edge (0/334);
soundness of the resulting object (0 unsound / 329); re-stamp no-op; cycles reachable; draw laws;
hop identity with E1′; the two structural-zero theorems; `tau` leak; the path cap; Wilson's root
condition; the v-structure artefact.

---

## 4. X1-P1 — catchability. The registered band was missed, downward

Registered prediction: `C_S1` ∈ **[0.40, 0.95]**. Registered surprise branch: *"< 0.40 — even a
coherence-checking tool is nearly blind to the class ⇒ S7 is a live hazard, not a coverage note."*

| ρ = 1, `original` | rate | 95 % Wilson | k / n |
|---|---|---|---|
| ARM R Meek-inconsistency (E1′ reference 0.3576, gate ±0.02) | **0.3559** | [0.3455, 0.3664] | 2 847 / 8 000 |
| ARM S2 (b-LOAD) catch rate | **0.0000** | — | **DEFINITIONAL** (PREREG §5.1) |
| skeleton-detectability of an S-uni statement | **1.0000** | — | **DEFINITIONAL** (PREREG §5.2) |
| **`C_S1` — intrinsic Dor–Tarsi predicate (M4)** | **0.2464** | **[0.2371, 0.2559]** | 1 971 / 8 000 |

**`C_S1` = 0.246 is below the registered band.** The "careful tool" branch fires: an operator that
runs a genuine coherence check still accepts **75.4 %** of single spurious required edges.

`C_S1` decomposed four ways, never quoted as one number (M4.3), ρ = 1, `original`, n = 8 000:

| branch | rate | 95 % Wilson | k |
|---|---|---|---|
| conflict (asserted edge clashes with a Meek-derived orientation) | 0.1650 | [0.1570, 0.1733] | 1 320 |
| directed cycle | 0.0814 | [0.0756, 0.0876] | 651 |
| non-extendable but acyclic | **0.0000** | [0.0000, 0.0005] | **0** |

**Non-extendability is entirely cycles** — 0 of 8 000 acyclic-but-non-extendable objects on
`original`, and the same zero on `licensed` and `large`. Dor–Tarsi buys nothing beyond a cycle
check on this population; that is a measurement, not an assumption, and it is why the intrinsic
predicate is so much weaker than the design expected.

**The `vstruct_vs_C` diagnostic, renamed and decomposed as M4.3 requires** (`original`, ρ = 1):
fires 0.4617 [0.4508, 0.4727]; of the whole population, **0.0836 [0.0778, 0.0899] is *purely* the
skeletal artefact** (an unshielded collider of `C` deleted because the new edge shields it — no
incoherence in it at all), and 0.3781 involves a genuinely new unshielded collider. AUDIT A4
predicted a ≥ 0.1094 artefact floor over the whole `N(C)` pool; measured 0.0836 over the drawn
statements. **For the record only**, the PREREG's original contaminated predicate
(`conflict ∨ cycle ∨ vstruct_vs_C`) gives 0.5441 [0.5332, 0.5550] — more than double the intrinsic
one, which is exactly the substitution M4 was written to prevent.

**`mpdag_valid` (M5).** The b-LOAD operator builds objects `O*` is not defined on at a rate that
grows fast with ρ:

| ρ | 1 | 2 | 3 | 4 |
|---|---|---|---|---|
| `original` S-uni | 0.8751 | 0.7127 | 0.5312 | **0.3640** |
| `large` S-uni | 0.9531 | 0.8897 | 0.8106 | 0.7192 |

At ρ = 4 on the primary grid, **64 % of what b-LOAD's semantics builds is not a PDAG of anything**.
Every rate in §5 is also reported restricted to `mpdag_valid = True`; the `ill_posed` stratum is a
finding about b-LOAD's semantics, not a number about `O*`. `path_blowup`: **0 / 30 000** everywhere
(the M16 cap never fired).

### What the resulting object *is* — the design question, answered by brute force

Measured in `test_x1.py::T5` over 396 small graphs (|U(G)| ≤ 12), by enumerating all consistent
extensions:

```
extendable            329/396 = 0.8308
sound                 329/329 = 1.0000   (Meek never over-orients)
maximally oriented    329/329 = 1.0000   (G == common orientation of its extensions)
```

**When `bk_assert`'s output is extendable at all, it is a genuine MPDAG — of a *different* MEC than
`C`'s.** When it is not (≈ 17 % at ρ = 1 on small graphs, 12.5 % on the run's `original`), it is a
PDAG of nothing and `O*` is ill-posed on it. So the design brief's questions (a) and (b) have
crisp answers: (a) the object is never an MPDAG of MEC(`C`), and is an MPDAG of another MEC exactly
when it is extendable; (b) "Meek-consistent" is **trivially false** for the class — verified two
ways on 4 616 statements: `is_consistent(C, [(a,b)])` is False for every non-adjacent pair, and
**no** member of the brute-forced MEC(`C`) carries the edge, because every member has `skeleton(C)`.
A consistency check that answers "no" for definitional reasons carries no information about this
class, which is why arm S1's predicate had to be made intrinsic.

---

## 5. X1-P2 — the headline. REFUTED on all three grids

ρ = 1, reachable denominator primary (M11), paired SCM-level cluster bootstrap, B = 10 000.

| ensemble | `P(silent \| S-uni)` | `P(silent \| consistent, R)` | **RR** | verdict |
|---|---|---|---|---|
| **`original`** | 0.0398 [0.0357, 0.0444] k=309 n=7 762 | 0.1451 [0.1354, 0.1554] k=688 n=4 741 | **0.274 [0.242, 0.308]** | **REFUTED** |
| `licensed` | 0.0546 [0.0498, 0.0599] k=424 n=7 765 | 0.1622 [0.1522, 0.1727] k=807 n=4 976 | **0.337 [0.299, 0.378]** | **REFUTED** |
| `large` | 0.0075 [0.0053, 0.0105] k=32 n=4 277 | 0.0343 [0.0286, 0.0411] k=111 n=3 236 | **0.218 [0.145, 0.304]** | **REFUTED** |

M3 satisfied: the risk-ratio bar is registered and applied on all three grids, so the three
verdicts are the same test.

**Free-bar disclosure (M1.3), printed beside the verdict as required.** On `original` the PREREG's
unconditional 0.03 bar would print SUPPORTED for an *identically* damaging class; it stops saying
SUPPORTED only below a conditional ratio of **0.845** and says REFUTED only below **0.540**.
Measured ratio **0.274**. The demoted absolute rule agrees: diff = −0.0516 [−0.0586, −0.0446],
intervals disjoint → REFUTED. **Decomposition of the unconditional difference (M1.2):**
free component (arm S2 never raising) = **+0.0162**; residual class effect = **−0.0678**.

**Robustness (M9.3, M5, M11), `original`:**

| view | S-uni rate | RR | verdict |
|---|---|---|---|
| reachable (primary) | 0.0398 | 0.274 [0.242, 0.308] | REFUTED |
| diluted (incl. unreachable) | 0.0386 | 0.289 [0.255, 0.325] | REFUTED |
| hop 0 only | 0.0405 | 0.186 [0.157, 0.217] | REFUTED |
| **direct-standardised onto arm R's hop distribution** | **0.0363** | — | weights printed in `results/x1_analysis_original.json` (hop 0 = 0.6146, hop 1 = 0.2503, hop 2 = 0.0462, hop 3 = 0.0083, hop 4 = 0.0006, unreachable = 0.0800) |
| `mpdag_valid = True` only | 0.0194 | 0.133 [0.111, 0.157] | REFUTED |

Every view points the same way, and the two opposite-signed confounds A9 named (the free bar, and
S-uni's flatter hop distribution) do not rescue SUPPORTED in either direction.

**The like-for-like operator comparison, added because M1's "not-raised" is not symmetric** (§9.2).
Conditioning *both* operators on their own coherence check passing:

```
P(silent | consistent_S1, S-uni) = 0.0110 [0.0086, 0.0140]  k=64  n=5 815
P(silent | consistent,    R    ) = 0.1451 [0.1354, 0.1554]  k=688 n=4 741
RR = 0.076 [0.058, 0.095]          licensed 0.062 [0.047, 0.079]   large 0.031 [0.007, 0.066]
```

A tool that runs a genuine coherence check and refuses when it fails is left with **1.1 %** silent
bias on the spurious class against **14.5 %** on reversals.

**Design effects (M13).** `R` 1.004 · `S-uni` 1.176 · `S-loc` **1.643** · (`licensed` S-uni 1.291,
S-loc 2.104; `large` S-loc 1.759). **Three exceed 1.25, so M13 is binding**: every headline rate is
reported with an SCM-level cluster-bootstrap CI beside its Wilson interval in
`results/x1_addendum.json`. The cluster interval is wider on S-loc (pooled ρ ≤ 4: Wilson
[0.0976, 0.1045] vs cluster **[0.0922, 0.1105]**) and changes no verdict. All primary contrasts
already used the paired cluster bootstrap.

---

## 6. X1-P3 — A26's `revive_if`, and why it inverts

**Composition first (M8) — the pooled S-loc number may not appear without this table.**
`original`, statements, n = 8 000:

| stratum | share | 95 % Wilson |
|---|---|---|
| `pure_spurious` | 0.7305 | [0.7207, 0.7401] |
| `reversal_of_true_edge` | **0.2695** | [0.2599, 0.2793] |
| query pair `{X,Y}` | 0.1495 | [0.1419, 0.1575] |
| SCMs with ≥ 1 query-pair statement | **0.5575** | [0.5356, 0.5791] |

AUDIT A8 measured 0.4855 for the last row from the pool; drawing 4 without replacement gives 0.5575.
**b-LOAD's own sampler is 27 % reversals, and puts the query pair itself into more than half the
`original` SCMs.** These are ORACLE labels (M14 ii): the analyst cannot perform this partition, so
no stratified rate here is a triage rule.

**S-loc damage, three rows plus the pooled row** (member unit — the baseline's unit — ρ ≤ 4,
`original`):

| stratum | rate | 95 % Wilson | k / n |
|---|---|---|---|
| **`pure_spurious`, non-query — P3's PRIMARY** | **0.0310** | [0.0282, 0.0340] · cluster **[0.0252, 0.0374]** | 428 / 13 805 |
| `pure_spurious`, query pair | 0.0538 | [0.0370, 0.0777] | 26 / 483 |
| `reversal_of_true_edge` | **0.2295** | [0.2150, 0.2446] | 712 / 3 103 |
| pooled (all of b-LOAD's sampler) | 0.1010 | [0.0976, 0.1045] · cluster [0.0922, 0.1105] | 3 030 / 30 000 |

**The in-run reversal baseline (M2.1) — the only like-for-like comparator**, arm R on X1's own
SCMs, per member, ρ ≤ 4, same grid, same `|K|`, same draw: unconditional **0.0784
[0.0754, 0.0815]** (2 352 / 30 000), conditional **0.1708 [0.1646, 0.1772]** (2 352 / 13 768).

The legacy figure is quoted **once**, as M2.2 requires and never as the comparator:
*754/8 085 = 0.0933 [0.0871, 0.0998] — pilot grid, generic+tiered arms, ρ ≤ 3, |K| ∈ {3,4}, member
unit, disjoint SCM draw (`default_rng(20260717)`).*

**Verdict on A26's `revive_if`, evaluated as literally as the data permit:**

- On the **primary** stratum, the CI **[0.0252, 0.0374] does not overlap [0.087, 0.100]** and lies
  **below** it. The PREREG predicted **above**, in [0.20, 0.60]. This is the registered
  **SURPRISE-LOW** branch (< 0.15). **A26's `revive_if` does not fire.**
- On the **pooled** row it *would* fire (cluster CI [0.0922, 0.1105] overlaps [0.087, 0.100]) —
  but that number is 27 % reversals-of-true-edges carrying 7.4× the damage of the pure class.
  **Reviving the b-LOAD transfer off the pooled row would revive it on the strength of the very
  class it was meant to be distinguished from.**
- **As AUDIT M2.4 predicted, the rule is not evaluable as written.** Under a never-raising operator
  the matched level for identical damage is 0.0933/(1 − 0.475) = **0.1776**, already above the
  window, so the literal rule can only fire when the class is ≥ 44 % *milder*. **No retraction is
  executed here.** The amended condition proposed for Zé: state the `revive_if` as a **conditional
  risk ratio against an in-run reversal baseline on the same SCMs**, with the S-loc composition
  stratified — not as an overlap with a fixed absolute window from a disjoint draw.

**X1-P5 (severity, A26's second condition).** Cliff's δ on `|est − τ|/|τ|` conditional on silent,
S-uni vs R, ρ = 1, `original`: **δ = −0.0788 [−0.1556, −0.0003]**, n_S = 309, n_R = 688, medians
0.5981 vs 0.6158. Registered rule: dominance declared only at δ ≥ 0.15 with CI excluding 0.
→ **NOT SEPARATED.** (`licensed` δ = −0.0538 [−0.1210, 0.0147]; `large` δ = −0.1633
[−0.3880, 0.0738].) The second half of the `revive_if` is met: conditional damage magnitude does
not stochastically dominate. **The class is rarer, not milder when it bites.** No p-value on ρ\*
appears anywhere in this run.

---

## 7. X1-P4 — locality under a skeleton-changing perturbation. SUPPORTED

Statement of the registered rule (M12.2): **SUPPORTED iff the hop ≥ 1 damage rate ≥ 0.10 × the
hop-0 rate with `CI_lo > 0`**; the "≥ 5 events" count survives only as an existence check.
Denominator printed before the verdict (M12.3).

| ρ = 1, S-uni | hop 0 | hop ≥ 1 reachable | unreachable | materiality | verdict |
|---|---|---|---|---|---|
| **`original`** | 0.0405 [0.0349, 0.0470] k=167 n=4 123 | **0.0390 [0.0332, 0.0458] k=142 n=3 639** | 0.0000 [0.0000, 0.0159] n=238 | **0.963** | **SUPPORTED** |
| `licensed` | 0.0573 [0.0503, 0.0653] k=210 n=3 663 | **0.0522 [0.0458, 0.0594] k=214 n=4 102** | 0.0000 n=235 | **0.910** | **SUPPORTED** |
| **`large`** (co-primary) | 0.0172 [0.0104, 0.0281] k=15 n=874 | **0.0050 [0.0031, 0.0080] k=17 n=3 403** | 0.0000 n=523 | **0.291** | **SUPPORTED** |

Full strata, every one printed with its `n`, none suppressed, the sentinel on its own row
(`original`, S-uni, ρ = 1):

```
hop 0            0.0405 [0.0349,0.0470]  k=167  n=4123
hop 1            0.0413 [0.0350,0.0487]  k=133  n=3220
hop 2            0.0240 [0.0127,0.0450]  k=  9  n= 375
hop 3            0.0000 [0.0000,0.0897]  k=  0  n=  39
hop 4            0.0000 [0.0000,0.4899]  k=  0  n=   4
hop 5            0.0000 [0.0000,0.7935]  k=  0  n=   1
unreachable      0.0000 [0.0000,0.0159]  k=  0  n= 238
```

and the same strata for **ARM R**, ρ = 1, conditional on consistency, on the *same* SCMs:

```
hop 0            0.2172 [0.2032,0.2319]  k=688  n=3167
hop 1            0.0000 [0.0000,0.0030]  k=  0  n=1290
hop 2            0.0000 [0.0000,0.0159]  k=  0  n= 238
hop 3            0.0000 [0.0000,0.0820]  k=  0  n=  43
hop 4            0.0000 [0.0000,0.5615]  k=  0  n=   3
unreachable      0.0000 [0.0000,0.0092]  k=  0  n= 412
```

Arm R's hop ≥ 1 zero replicates on `licensed` (0 of 2 297 reachable) and `large` (0 of 2 406).
**These are X1's own fresh, code-generated denominators for the reversal class at ρ = 1. X1 prints
no `0 of 791` and no derivative of it; X2 owns S5's wording, the ρ ≥ 2 counterexample and that
retirement.**

**Edge-attributable, drop-controlled** (`O*(S) ≠ O*(D)` in set or identifiability at matched
(ρ, flip)), hop ≥ 1: `original` **0.1657 [0.1540, 0.1781]** k=603 n=3 639; `licensed` 0.1643;
`large` 0.0644. Arm D contributes **zero** silent events anywhere, so all 142/214/17 silent events
at hop ≥ 1 are attributable to the added edge, not to the withdrawal.

**`hop_C × hop_G` (M12.1) is diagonal, and that is a proof, not a finding.** For a single asserted
pair `(a,b)`, adding the edge gives `dist_G[a] ≥ min(dist_C[a], dist_C[b]+1)` and symmetrically, so
`min(dist_G[a], dist_G[b]) = min(dist_C[a], dist_C[b]) = hop_C` exactly. Every off-diagonal cell of
the cross-tab is empty in all three ensembles. **The mechanism S9 turns on is therefore *not* "the
assertion pulls itself into the neighbourhood"** — its own hop is invariant. The cross-tab was
implemented as M12 registered and it carries no information at ρ = 1; §9.4 records the objection.

Two candidate mechanisms remain and **X1 does not discriminate between them** (NOT RUN): (i) the
added edge lies on a *new* proper possibly-causal path from `X` to `Y`, changing `cn` and `forb`
directly; (ii) Meek closure propagates the new orientation into the query neighbourhood — which is
precisely the proof obligation reviewer R1 named for S8 (*"no distant statement can enter the
neighbourhood through Meek closure"*). Either way, **a locality-based pruning heuristic is unsafe
for this class**: 142 counterexamples on `original` at ρ = 1 alone.

---

## 8. X1-P6 — the movement channel, three columns, never merged

Censoring at n = 20 000, z = 1.960, `original`, n = 2 000 SCMs. Reported as **three separate
columns**; `ρ*_any = min(ρ*_se, ρ*_ident)` is printed but never quoted alone (PREREG §11.3).

| arm | `ρ*_se` censored | `ρ*_ident` censored | `ρ*_any` censored |
|---|---|---|---|
| R | 0.5520 [0.5301, 0.5737] | 0.2270 [0.2092, 0.2459] | 0.0795 [0.0684, 0.0922] |
| S-uni | 0.6855 [0.6648, 0.7055] | 0.1900 [0.1734, 0.2078] | 0.1380 [0.1236, 0.1538] |
| S-loc | 0.7380 [0.7183, 0.7568] | 0.0910 [0.0792, 0.1044] | 0.0595 [0.0500, 0.0707] |
| **D (drop)** | **1.0000 [0.9981, 1.0000]** | 0.0785 [0.0675, 0.0911] | 0.0785 [0.0675, 0.0911] |

Arm D's `ρ*_se` censoring is **exactly 1.000**, and the reason is an identity worth writing down:
**all valid adjustment sets give the same population estimate**, so `est0 = τ` whenever `O0` is
valid — verified 0 violations / 5 200 SCMs, max |diff| 2.6e-14, and `O0` was valid in every single
SCM of all three ensembles. Hence `|est_m − est0| > 0 ⟺ O_m invalid ⟺ silent`. This is the same
fact E1′ recorded as `changed_still_valid = 0`, seen from the other side, and it means the movement
channel and the silent-bias channel are **the same event at n = ∞** in this population.

**Two other measurements that belong to the operator, reported because they are surprising:**

- **G8, closure order.** `close-after-each` and `stamp-all-then-close` produce a **different graph**
  in **0.2818 [0.2782, 0.2854]** of members on `original` (n = 60 000), 0.2487 on `licensed`,
  0.0954 on `large`. b-LOAD's `initialize_background_knowledge` builds the whole matrix first;
  Perković Alg. 1 closes after each. **These are not the same algorithm on this class**, and every
  number in this report is on the pilot-matched close-after-each convention.
- **G7, re-stamp.** 0 divergences in 60 000 members under both scopes (full-prefix, per M13 #4, and
  current-only). `meek_closure` only orients undirected pairs, so it cannot flip a stamped edge;
  b-LOAD's `g[mask] = bk[mask]` is a genuine no-op *here* — but only because nothing else writes to
  those cells, which is not true inside `mb_by_mb`.

---

## 9. Objections recorded, mitigations implemented anyway

All 16 mitigations were implemented. Four are recorded as disputed, per the brief.

1. **AUDIT A6/M6 caught the special case and missed the general one.** The registered replacement
   law confounds every arm-S rate with a withdrawn true statement. M6 fixed it for `P-null` only.
   Fixed here by adding **arm D**, which is not in the PREREG and not in the audit. Measured
   consequence: the confound is **null on `silent`** and **total on `ostar_changed`**, so the
   headline survives and any future contrast on outcome (b) does not.
2. **M1's "conditional" RR is not symmetric.** For arm S2, "operator did not raise" ≡ all members,
   so M1's conditional rate *is* the unconditional rate, and the RR compares everything b-LOAD
   builds against only what a Meek-checking tool accepts. Implemented as registered (RR = 0.274) and
   supplemented with the genuinely symmetric S1-vs-R contrast (RR = 0.076). Both REFUTE; the
   conclusion does not depend on the choice.
3. **M11's justification generalises only to the edge.** "Unreachable statements have damage rate
   0.000 by construction" is true of the *added edge* — 0 of 335 in machinery test T10, 0 of 1 662
   across the three ensembles in the run (318 + 351 + 993) — and is **not** a construction argument
   for the *replacement*. G12 is reported in the corrected edge-attributable form as well as the
   registered raw form; both pass, but only the first is entailed. T10b shows the withdrawal moves
   `O*` at 0.2539 over all hops, so the registered form passes here by a property of the data
   (in the unreachable stratum the withdrawal happens to change neither `O*` nor amenability),
   not by the argument M11 gives for it.
4. **M12.1's `hop_C × hop_G` cross-tab is degenerate at ρ = 1** by the two-line identity in §7.
   Implemented and reported; it carries no information, and the question it was meant to answer
   ("does the assertion pull itself into the neighbourhood?") is answered **no, provably**, which
   makes the mechanism question harder rather than easier.

---

## 10. Which abstract sentence this speaks to

| sentence | X1's finding |
|---|---|
| **S7** *"errors cover … required edges asserted between non-adjacent variables"* | **Now has evidence, and it is a two-strata table, not a blend.** The class is real, is never caught by b-LOAD's operator (definitional), is caught by an intrinsically coherent tool only 24.6 % of the time, and is **~3.6× less likely to bias the estimate** than a reversal. The sentence can stand only alongside both rows. |
| **S3** *"two thirds pass the consistency check"* | **A reversal-only figure.** ARM R pass rate `original` ρ=1 = 0.6441 [0.6336, 0.6545]. Arm S2 pass rate ≡ 1.000 **by definition**; arm S1's 0.7536 is a **different predicate** (intrinsic extendability), because "consistent with C" is trivially false for this class. The three must be printed on **separate lines** with their predicates spelled out; they are not one row. |
| **S4** *"14.1 % of the consistent errors go on to bias the estimate"* | Reproduced for reversals on X1's own SCMs: **0.1451 [0.1354, 0.1554]** (`original`, ρ=1, reachable, n=4 741). **It does not transfer to the spurious class: 0.0398 [0.0357, 0.0444].** S4 must be re-scoped to reversals and a second row added. |
| **S5 / S8 / S9** | **X1 supplies the spurious-class half only.** At ρ = 1 on `original`, arm R damages 0 of 1 574 statements at hop ≥ 1 while arm S-uni damages 142 of 3 639 at essentially the hop-0 rate (materiality 0.963). **Locality is a property of the reversal operator.** A pruning heuristic (S8) is unsafe for the class S7 names, and "knowledge outside the neighbourhood need not be elicited" (S9) is false for it. X2 owns S5's wording and the reversal class. |
| **S2 vs S7** | PREREG §5.3's analytic point stands and X1 **reports** it rather than repairing it: within an MEC all members entail the same CIs, so orientations are untestable — but **adjacency is not in that equivalence class**, so a spurious required edge *is* testable with power → 1. As both are currently written, S2 and S7 are inconsistent. **This is a wording decision for Zé and Chris, not a measurement, and it must not be presented as one.** |

---

## 11. Limitations — in writing, so a gap cannot later read as coverage

Everything in PREREG §8 stands unchanged, plus what the run added:

1. **An estimated CPDAG is NOT RUN.** `C = dag_to_cpdag(D)` is the oracle's. An analyst running
   PC/GES would *see* that the asserted pair is absent from `Ĉ`'s skeleton — the whole point of the
   testability asymmetry. X1 measures damage **conditional on the tool not performing that check**.
   For the spurious class specifically this is the first question a reviewer asks, and the honest
   answer is "not tested", not "robust".
2. **This is not a b-LOAD replication.** It emulates the *semantics* of
   `initialize_background_knowledge` + `update_graph_mpdag` on the pilot's population harness. No
   `mb_by_mb`, no CI tests, no finite samples, no R. And G8 shows the two closure conventions
   disagree on 28 % of members, so "b-LOAD's semantics" is itself convention-dependent.
3. **The two damage mechanisms at hop ≥ 1 are not discriminated** (§7). New-causal-path vs
   Meek-propagation is the discrimination S8's proof obligation needs, and X1 does not supply it.
4. **`|K| = 4` only; the `k8` ensemble is NOT RUN** (X3's axis). **Mixed knowledge sets are NOT
   covered**: X1 measures the two pure strata and they bracket b-LOAD's mixture only if the two
   error types do not interact inside Meek closure, which is untested.
5. **Forbidden-edge knowledge, tiered knowledge, nonlinear mechanisms, latent confounding,
   non-Gaussian errors, multi-node X or Y, real data: none of it.** Linear-Gaussian iSCM, ER DAGs,
   single X, single Y, population Σ.
6. **`P(|O0| = 0) = 0.66` on the primary grid** (A15). Two thirds of `original`'s queries have an
   empty optimal adjustment set, so "adjustment-set fragility" is measured at the boundary of the
   object; `large` (`P = 0.23`) is co-primary for P4/P5 for that reason and agrees on every verdict.
7. **`P-null-identity` was run at ρ = 1 only** (4 independent passes through the full `bk_assert`
   path per SCM, 8 000 total), not over the full ball, because all 15 members are identical by
   construction. Stated rather than silently trimmed.
8. **Eligibility drops the densest small SCMs.** M9.2's single rule (`|N(C)| ≥ 4` and
   `|b-LOAD pool| ≥ 4`, applied to **all** arms) removed 213 of `original`'s SCMs at mean
   `|N(C)| = 2.72`, mean `p = 5.04`. Arm R had no such exclusion in E1′. `large` is immune.
9. **S-loc's sampler is oracle-guaranteed-false** (M14 i) — its pool is defined by `(a,b) ∉ E(D)`.
   That is b-LOAD's own design, and it is stated rather than assumed away.
10. **Carried forward, unclosed, and now dated (AUDIT A18).** `wiki/activities/active-claims.md`
    still holds **no** ρ-breakdown entry and `grep -rl 'relevance:.*rho-breakdown' wiki/literature/`
    returns **0**. `active-claims.md:570` records the b-LOAD graph-distance falsifier as *"currently
    NOT firing"*. §7 measures the spurious-class half of exactly that falsifier and it **fires** for
    this class. The b-LOAD CIKM camera-ready is **2026-08-23**. X1 supplies the measurement; it does
    not edit the register, and none of this may enter a manuscript before Zé registers the claim in
    his own words (CLAUDE.md Operating Rule #1). **Proposed, not committed by an agent.**

---

## 12. Artefacts

| what | where |
|---|---|
| code (12 modules) | `code/` here, and `~/latent-causal/x1-spurious/code/` on betelgeuse |
| machinery tests (15/15) | `code/test_x1.py`, output `logs/test_x1.log` |
| analysis JSON, 3 ensembles | `results/x1_analysis_{original,licensed,large}.json` |
| gates G2 / G11 | `results/x1_gates_g2_g11.json` |
| M13 cluster CIs, drop-arm decomposition, `est0 = τ` check | `results/x1_addendum.json` |
| human-readable summaries | `logs/summary_P1.txt`, `logs/summary_P2_P6.txt`, `logs/summary_addendum.txt`, `logs/summary_extra.txt` |
| run log | `logs/run_x1.log` |
| **raw per-SCM JSON (108 MB — deliberately NOT copied into the vault)** | `~/latent-causal/x1-spurious/results/x1_{original,licensed,large}.json` on betelgeuse |

Reproduce: `python test_x1.py` → `python run_x1.py <n> <ensemble> ../results/x1_<ensemble>.json 12`
→ `python analyse_x1.py <ensemble>` → `python gate_g2_g11.py` → `python addendum_x1.py`, all with
cwd = `code/` and `~/e1venv/bin/python`.
