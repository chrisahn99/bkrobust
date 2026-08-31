# E1′ RESULTS — the `se` criterion delivers the analyst-facing property, and **refutes** the premise that it fixes degeneracy

**Run 2026-08-19 on betelgeuse** (16 cores, `~/e1venv` Python 3.12.5, numpy 2.4.4, scipy 1.18.0,
networkx 3.6.1). Wall clock ≈ 6 min, 12 workers, no GPU, no network.
Pre-registration: [`PREREG.md`](PREREG.md) + Addendum 1, both written before the run.
Design audit: [`AUDIT.md`](AUDIT.md).

**13,554 SCMs across four ensembles · 2,500 SCMs × 20 datasets × 4 sample sizes in ARM 2 ·
0 regeneration-gate failures.**

---

## 0. Verdict in one paragraph

The `se` criterion **does** what `RANKING.md:839` said it would: it removes `tau` from the
instrument, so ρ\* becomes computable by an analyst who does not have the answer key. That leg of
the paper is now demonstrated rather than asserted. **It does not fix the degeneracy the reviewer
objected to, and the reason is that degeneracy was never about the yardstick.** On the licensed
ensemble the censored mass moves from **0.554 (n=200) to 0.515 (n=∞)** — 3.9 points across five
orders of magnitude of sample size, against a hard floor — and **94.6 % of that censored mass is an
exact zero**: no combination of up to four expert reversals moves the estimate by any amount at all.
The composite `ρ*_any = min(ρ*_se, ρ*_ident)` cuts censoring to **0.166**, but converts
top-degeneracy into bottom-degeneracy (69 % of the mass at ρ\* = 1) and so fails the same
pre-registered non-degeneracy profile. **What the run does produce is a single law that explains E1,
E3 and the large-graph inversion at once, and it is an elicitation-design law, not a bound.**

---

## 1. Integrity gates

| gate | licensed | original | large | k8/K4 | k8/K6 | k8/K8 |
|---|---|---|---|---|---|---|
| **G1** silent bias @ρ=1 (pilot denominator) | **0.161** ✅ | **0.149** ✅ | 0.026 ❌ | **0.173** ✅ | 0.223 ❌ | 0.249 ❌ |
| Meek-inconsistent @ρ=1 (pilot 0.330) | 0.324 | 0.358 | 0.172 | 0.215 | 0.363 | 0.487 |
| **G3** ρ\*_se monotone in n | ✅ 0 viol. | ✅ 0 | ✅ 0 | ✅ 0 | ✅ 0 | ✅ 0 |
| **G4** d(ρ) monotone in ρ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ |
| **G5** placebo (O\* frozen at O₀) censored | **1.000** ✅ | 1.000 ✅ | 1.000 ✅ | 1.000 ✅ | 1.000 ✅ | 1.000 ✅ |
| **G6** floor at SE = 0 equals frac(d ≡ 0) | ✅ exact | ✅ | ✅ | ✅ | ✅ | ✅ |
| `changed_still_valid` | **0** | 0 | 0 | 0 | 0 | 0 |

**G5 is the gate that matters.** The frozen-O₀ ball returns censoring **exactly 1.000 at every n and
every z ∈ {1, 1.96, 2.576}** — 39 cells per ensemble. The instrument is not a dead switch and its
dynamic range is demonstrated *inside this run*, not borrowed from E1's cross-criterion spread.

### 🔴 Two gate corrections, both disclosed rather than quietly fixed

1. **G1's denominator was wrong in my first cut.** I divided by *all* ρ=1 members; E1's published
   0.141 uses **consistent members only** (`RESULTS.md` §6: *"on the pilot's own denominator"*). The
   first reading was low by exactly the Meek-inconsistency rate (0.1085 = 0.161 × (1 − 0.324)).
   Corrected in `analyse_arm1.py`, both denominators now reported, and the run re-analysed. No
   simulation was repeated and no verdict depended on it.
2. **G1 fails on `large`, `k8/K6`, `k8/K8` — and the failures are the finding, not a bug.** The
   [0.12, 0.18] band was calibrated on the pilot's |K| = 4, p ≤ 10 measurement, and the gate passes
   on **every** ensemble in that box (0.161 / 0.149 / 0.173). It fails where silent bias *should*
   move: down on `large` (0.026 — statements are far from the query, §4) and up with |K| (0.223 →
   0.249 — more statements, more chances one lands near the query). The gate is scope-limited; it is
   not evidence of an unfaithful re-implementation.
3. **G2's constant was mis-specified in `PREREG.md`.** I wrote *"within ±0.03 of 0.336"*, which is
   E1's **`original`**-grid number, while the `k8` ensemble here uses the **licensed** grid, whose
   E1 number is **0.271**. Measured: `rel50` @ |K|=8 = **0.284**. Against the correct comparator
   |0.284 − 0.271| = 0.013 → **PASS**; against the constant as literally written → FAIL. Both stated.

---

## 2. ARM 1 — the sweep. The n-axis is nearly inert.

**`licensed`, |K| ≤ 4, unconditioned, full ball ρ ≤ 4, n = 2,087 amenable SCMs, z = 1.96.**

| n | ρ\*_se = 1 | 2 | 3 | 4 | **CENS** | ρ\*_any = 1 | 2 | 3 | 4 | **CENS** |
|---|---|---|---|---|---|---|---|---|---|---|
| 50 | 0.265 | 0.102 | 0.017 | 0.001 | **0.614** | 0.622 | 0.139 | 0.024 | 0.001 | **0.215** |
| 200 | 0.321 | 0.108 | 0.017 | 0.001 | **0.554** | 0.665 | 0.130 | 0.021 | 0.001 | **0.183** |
| 2 000 | 0.347 | 0.110 | 0.016 | 0.001 | 0.526 | 0.685 | 0.125 | 0.021 | 0.001 | 0.170 |
| 20 000 | 0.352 | 0.112 | 0.016 | 0.001 | **0.519** | 0.690 | 0.123 | 0.020 | 0.001 | **0.167** |
| 10⁶ | 0.354 | 0.114 | 0.017 | 0.001 | 0.515 | 0.691 | 0.123 | 0.020 | 0.001 | 0.166 |
| ∞ | 0.354 | 0.113 | 0.017 | 0.001 | **0.515** | 0.691 | 0.123 | 0.020 | 0.001 | **0.166** |

**Five orders of magnitude of sample size buy 9.9 points of censoring (ρ\*_se) and 4.9 points
(ρ\*_any), against floors of 0.515 and 0.166.** Per `AUDIT.md` A1's binding mitigation:
`censored(200) = 0.183` for the primary statistic, i.e. **≤ 0.50 already at the bottom of the grid —
the primary bar was cleared without the n-axis doing any work.** That is stated, not smoothed.

### Pre-registered decision

| statistic | n50 | nND | verdict |
|---|---|---|---|
| **ρ\*_any** (Addendum 1 PRIMARY) | 50 (the grid floor) | **does not exist** | **REFUTED** |
| ρ\*_se (secondary, isolating) | does not exist | 100 | **REFUTED** |
| ρ\*_ident alone (free) | — | — | censored 0.406 |

Both fail, for **opposite** reasons, and that is the result:

- **ρ\*_se is censored-degenerate** — 51.5 % of problems never break.
- **ρ\*_any is one-degenerate** — 69.0 % of problems break at ρ\* = 1, so `max_bin = 0.690 > 0.60`.

> **The reviewer's objection survives both repairs.** `sign` at |K| ≤ 4 says *"safe"* for 88.5 % of
> problems; `any` says *"one error is enough"* for 69.0 %. Two-valued either way.
> ⚠️ **Interpretive note, flagged as such:** these two failures are not equally bad. `CENSORED` means
> *the diagnostic declines to speak*; `ρ\* = 1` means *the diagnostic says something actionable*.
> The pre-registered profile rule does not distinguish them, and it is honoured as written.

### The censored bin was never one thing (Addendum 1's mandated output)

`licensed`, ρ\*_se, n = ∞, 1,074 censored SCMs:

| what actually happened | share of censored mass |
|---|---|
| all perturbations Meek-inconsistent (**caught free**) | **0.000** |
| all perturbations destroy amenability (**loud**) | 0.054 |
| consistent, amenable, and the estimate **never moves** | **0.946** |
| — *of which `d ≡ 0` exactly* | **0.946** |

**94.6 % of the censored mass is an exact zero, not an unresolved value.** This is the real reply to
*"your diagnostic returns 'no breakdown' nine problems in ten"*: for two-thirds of problems the
ignorance interval is **a point**, and "no breakdown" is the correct answer rather than a failure to
resolve. The paper should make that argument instead of chasing a smaller censoring number.

### The continuous statistic — bimodal, and it is what should be printed

`R(ρ=1, n) = d(1) / (z·SE_n)`, licensed:

| n | frac exactly 0 | frac > 1 | mean | q90 |
|---|---|---|---|---|
| 200 | **0.646** | 0.321 | 1.42 | 5.08 |
| 2 000 | 0.646 | 0.347 | 4.54 | 16.18 |
| 20 000 | **0.646** | 0.352 | 14.35 | **51.18** |

A point mass at exactly zero (64.6 %, **n-invariant**) plus a heavy right tail. `R` has no bins and
no ceiling, and it shows the bimodality every binned version hides. **Report `R`; report ρ\* as its
thresholding.**

---

## 3. 🟢 THE RESULT — one locality law explains E1, E3 and the large-graph inversion

Censoring of ρ\*_any at n = 20 000, stratified by the **minimum hop-distance from any asserted
statement to {X, Y}** on the CPDAG skeleton:

| ensemble | dist 0 | dist 1 | dist ≥ 2 |
|---|---|---|---|
| `original` (p 5–8) | **0.037** (n=2 539) | **1.000** (n=93) | — |
| `licensed` (p 5–10) | **0.070** (n=1 852) | **0.935** (n=229) | — |
| `k8` \|K\|=8 | **0.016** (n=563) | 0.955 (n=22) | — |
| `large` (p 15–25) | **0.330** (n=1 472) | 0.977 (n=638) | **1.000** (n=331) |

**A 13–27× separation, invariant to n to four decimals** (0.9345 at n = 2 000, 20 000 and ∞ alike),
and **exactly 1.000** for every one of the 331 `large` SCMs whose statements sit two or more hops
from the query.

> **The breakdown radius is informative if and only if the elicited knowledge touches the query
> neighbourhood.** ρ\* is a functional of `O*`, `O*` is local, so a statement that does not reach the
> query cannot move it — and no sample size repairs that, because zero times anything is zero.

**This unifies three separate findings the campaign carried as unrelated:**

1. **E1's large-graph inversion is not a new phenomenon.** `large` fails because only **60.3 %** of
   its SCMs have a statement at distance 0, against **88.7 %** (licensed) and **96.3 %** (original).
   Condition on distance 0 and `large` behaves like the others (0.330, still worse — larger `O*`,
   more room to be inert — but not a different regime). E1 named O\*-locality as the cause; this
   measures it as the *governing variable*, on the licensed box as well.
2. **E3's conclusion is promoted.** E3 found the predictive content of the O\*-relevance metric is
   *causal-path locality*, not `pa(cn) \ forb` membership (which is anti-predictive, AUC 0.409).
   That was a statement about a per-member feature. Here the same variable governs the **instrument
   itself**.
3. **E2 already pointed the same way.** Its surviving finding was an **elicitation-design** claim
   (*"it does not matter that you elicited a tiering; it matters that the expert's mistake is a
   tier-level mistake"*). Two of four experiments now converge: **the paper's thesis is about how to
   elicit, not about how to bound.**

### The registered `large` prediction (`PREREG.md` §4)

Predicted: `large` censoring at n = 2×10⁴ strictly inside (0.30, 0.85).

| statistic | measured | against the prediction |
|---|---|---|
| ρ\*_se | **0.9034** | ❌ **the "n buys nothing" surprise fired**, exactly as pre-named |
| ρ\*_any | **0.5903** | ✅ inside |

**The registered falsifier of `RANKING.md:346`, evaluated literally** (censoring > 0.755 at every
n ∈ {200 … 2×10⁴, ∞}): **FIRES on `large` for ρ\*_se**, does not fire on `licensed`, `original` or
`k8`. The licensed-scope boundary therefore hardens again: **licensed to p ≤ 10 because REFUTED at
p ≥ 15, with the mechanism now measured rather than named.**

### Honest-risk comparators (`PREREG.md` §3c + Addendum 1) — both legs survive

| comparator | Spearman | bar | outcome |
|---|---|---|---|
| ρ\*_se vs the free t-statistic \|est₀\|/SE (n=2 000) | **+0.199** | > 0.90 ⇒ one sentence | ✅ **not a relabelled t** |
| ρ\*_any vs the free ρ\*_ident (n=2 000) | +0.571 | > 0.95 ⇒ `se` half is nothing | ✅ `se` strictly lowers ρ\* on **30.0 %** of problems |

---

## 4. ARM 2 — the analyst simulation. What a wrong expert actually costs.

2,500 SCMs × 20 datasets × n ∈ {200, 1 000, 5 000, 20 000} × r ∈ {0,1,2}. The analyst is handed only
`(C, K_asserted, Σ̂, n)`. `tau` adjudicates coverage and enters nowhere else.
**Regeneration gate: 0 failures / 2,500** (E3's G1–G5 pattern: every SCM rebuilt from `(seed,p,deg)`
and `|τ − τ_stored| < 1e-9`).

**G7 PASSES exactly:** at r = 0 (expert correct), n = 20 000, naive coverage = **0.9509
[0.9481, 0.9533]** against nominal 0.95. The finite-sample machinery is calibrated.

**Capability check PASSES:** a wrong expert is expensive.

| n = 20 000 | r = 0 | r = 1 | r = 2 |
|---|---|---|---|
| naive CI coverage | **0.951** | **0.720** [0.713, 0.727] | **0.568** [0.559, 0.577] |

> **One reversed background-knowledge statement costs 23 points of coverage; two cost 38.**
> Meek catches **23.6 %** of single errors and **38.6 %** of double errors *for free* (raises
> `MeekFail`); those are excluded from the table above, so 0.720 is the cost **among errors Meek
> does not catch**.

### The coverage–width frontier, r = 1, n = 20 000

| interval | coverage | width ratio median / mean / q90 | ball inert |
|---|---|---|---|
| `naive` (ρ=0) | 0.720 | 1.00 / 1.00 / 1.00 | — |
| `robust(1)` | **0.884** | **1.00** / 6.84 / 23.3 | 69.3 % |
| `robust(2)` | 0.925 | 1.00 / 11.62 / 34.4 | 52.7 % |
| `robust(4)` = full ball | 0.931 | 1.00 / 12.98 / 36.2 | 50.1 % |

**Pre-registered rule → INCONCLUSIVE.** Coverage 0.8835 misses the ≥ 0.90 bar; median width ratio
1.00 clears the < 3× bar. Neither REFUTED threshold (coverage < 0.80, ratio > 10×) is hit *on the
median*.

🔴 **Disclosure, and it matters: I pre-registered the MEDIAN width ratio, and the median flatters the
method.** The ball is inert for 69 % of problems, so the median is 1.00 by arithmetic. On the
**mean** the ratio is **6.84×** at n = 20 000, which would have triggered the > 3× failure. The
pre-registered statistic is honoured; the mean is reported beside it and the paper must quote both.

### 🟢 The structural result of ARM 2 — ignorance does not shrink with data

At ρ = 1, r = 1, the mean width ratio against sample size:

| n | 200 | 1 000 | 5 000 | 20 000 |
|---|---|---|---|---|
| mean width ratio | 1.61 | 2.32 | 3.93 | **6.84** |
| q90 | 3.36 | 6.12 | 12.26 | **23.30** |

`SE_n ∝ n^{−1/2}` while `d(ρ)` is **n-free**. So **background-knowledge ignorance is a fixed-width
component that sampling error eventually stops hiding.** By n = 20 000 the honest interval is ~7×
the confidence interval on average and 23× at the 90th percentile. This is the Γ↔ρ sensitivity
statement in its natural units and it is the one number in this run that is neither definitional nor
already on the record.

### The non-definitional cell (`AUDIT.md` A6): the ball misses the truth

r = 2, ρ = 1 — the analyst budgets for one error and the expert made two:

| n | 200 | 1 000 | 5 000 | 20 000 |
|---|---|---|---|---|
| coverage | 0.809 | 0.785 | 0.772 | **0.772** |

Degrades gracefully (naive would be 0.568) and **stops improving past n = 5 000** — the same
ignorance-floor signature.

⚠️ **Even the full ball never restores nominal coverage**: 0.931 at r=1, 0.912 at r=2. The residue is
the ~7 % of cells where the truth-restoring ball member is itself non-amenable and therefore carries
no estimate to widen with. That is a real limitation of the construction, not noise.

---

## 5. What this changes for the ICLR 2027 paper

1. 🔴 **`RANKING.md:839` prices `crit='se'` as *"the single highest-value remaining experiment"*
   because it converts an oracle statistic into an analyst-facing one. Half of that is confirmed and
   half is refuted.** It *is* analyst-facing now — demonstrated, not asserted, and the `tau`-free
   path is enforced in code. It is **not** the fix to degeneracy, and nothing indexed by `n` can be:
   64.6 % of licensed problems have `d(1) = 0` **exactly**, at every sample size.
2. 🟢 **The paper's §3 should lead with locality, not with a radius.** One measured law
   (0.037–0.070 vs 0.935–1.000 censoring by hop-distance) explains E1's inversion, absorbs E3, and
   agrees with E2's surviving finding. It is a prescription an analyst can act on — *elicit near the
   query* — and it is the only thing here that is neither definitional nor already published.
3. 🟢 **ARM 2 supplies the utility number the abstract needs**: one wrong statement costs 23 points
   of coverage; the ρ=1 robust interval buys 16 of them back; the price is a mean 6.8× interval at
   n = 20 000 and **nothing at all** for the 69 % of problems whose ball is inert.
4. ⚠️ **Never write "ρ\* is a well-behaved 4-valued diagnostic".** Both repairs fail the
   pre-registered profile, in opposite directions. Write the **continuous** `R(ρ,n)` with its point
   mass at zero, and quote the width frontier with the **mean** beside the median.
5. ⚠️ **`sign` stays dropped** (E1's adversarial verdict), and `rel50` is now visibly redundant:
   `ident` alone (0.406 censored) beats `rel50` (0.677) for free on the licensed box.
6. 🔴 **Unchanged and still blocking (`AUDIT.md` A10):** there is **no ρ-breakdown entry** in
   `wiki/activities/active-claims.md`, and `grep -rl 'relevance:.*rho-breakdown' wiki/literature/`
   returns zero notes, so `check_claims.py` is blind to Taeb–Guo–Henckel 2511.10625, gadjid
   2402.08616, Fang 2207.05067 and CausalGuard 2605.21928. **Three claim ids need Zé's own wording
   before any of the above is written into a manuscript.**
7. ⚠️ **NOT RUN and not claimed** (`AUDIT.md` A2): the CPDAG is the oracle's throughout. ρ\*_se is
   *"computable from (C, K, Σ̂, n)"* — never *"robust under a discovered CPDAG"*. `causal-learn` is
   absent from every local env; ≈ 1 day plus an install.

---

## 6. Reproduction

```bash
ssh betelgeuse && cd ~/latent-causal/e1prime-se/code
PY=$HOME/e1venv/bin/python NW=12 ./run_all.sh      # ~6 min, 16 cores, no GPU, no network
```

Seeds fixed (`20260819`); ordered `imap` so the kept set is the first *n* passing seeds.
`graphs.py`, `adjust.py`, `scm.py` are the pilot's, imported **unmodified**; `se.py`,
`run_arm1.py`, `run_arm2.py`, `analyse_arm*.py` are new. Nothing under
`output/2026-08-14_latent-causal-iclr2027-rerank/` was modified.

Kept ensembles: `licensed` 4 000 SCMs (2 087 amenable) · `original` 3 372 (2 637) ·
`large` 3 000 (2 443) · `k8` 1 200 (127/329/585 at |K| = 4/6/8) · ARM 2 2 500 × 20 × 4 × 3.
