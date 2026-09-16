# E1′ PRE-REGISTRATION — the `se` criterion: ρ* stops being an oracle statistic

**Written 2026-08-19, BEFORE any `se`-criterion number was computed anywhere.**
Baselines quoted below are already on the record (`../2026-08-14_latent-causal-iclr2027-rerank/RESULTS.md`,
`results/ksweep_*.json`, `results/summary.json`). Nothing under `crit='se'` has ever been run:
`grep -rn "crit.*se" ../2026-08-14_latent-causal-iclr2027-rerank/code/` returns no such branch, and
`RANKING.md:839` lists it verbatim under **"Experiments not run"**.

---

## 0. What E1 left broken, and which of those this run fixes

E1 (14/08) answered *"is ρ\* degenerate?"* by raising |K| from 4 to 8 and switching the criterion
from `sign` to `rel50`. It returned SUPPORTED on 3 of 4 ensembles, and 3 of 3 adversarial verifiers
refuted the **interpretation**. Six defects are on the record:

| # | defect | fixed here? |
|---|---|---|
| D1 | **ρ\* is an ORACLE statistic** — both criteria compare the perturbed estimate to `tau`, the true total effect (`analyse.py:57,59`). No analyst has `tau`. | ✅ **this is the whole run** |
| D2 | the pass bar (censored ≤ 0.60) was **already cleared for free** at \|K\|≤4 by `crit='ident'` (0.306) | ✅ bar is now a **location**, §3 |
| D3 | `sign`→`rel50` is a **nesting** relation (sign flip ⇒ \|est−τ\| ≥ \|τ\| > 0.5\|τ\|), so the headline movement is partly definitional | ✅ declared definitional ex ante, §5 |
| D4 | \|K\|=8 is definable on **1.5 %** of the licensed ensemble; ~half the gain was **graph conditioning**, not \|K\| | ✅ primary arm is **\|K\| ≤ 4, unconditioned** |
| D5 | top-censoring is **budget-imposed** (`MAX_RHO = 3`), not measured | ✅ **full ball**, ρ = 1…\|K\| |
| D6 | the instrument **inverts** at p ∈ [15,25] (censored 0.877), cause = O\*-locality | ⚠️ measured again under `se`; **prediction registered**, §4 |

---

## 1. The estimand

For a query (X,Y), a CPDAG `C`, an asserted background knowledge set `K` (|K| ≤ 4), and a sample
size `n`:

```
G0   = MeekClosure(C, K)                      O0 = O*(X,Y,G0)         est0 = beta_X( Y ~ X + O0 )
SE_n = sqrt( sigma0^2 * [Sigma_{S0 S0}^{-1}]_{XX} / n ),  S0 = {X} u O0,
       sigma0^2 = Sigma_YY - Sigma_{Y S0} Sigma_{S0 S0}^{-1} Sigma_{S0 Y}
d(rho) = max over Meek-consistent, amenable K' at Hamming distance <= rho of | est(K') - est0 |

rho*_se(n, z) = min { rho : d(rho) > z * SE_n }        z = 1.96 primary
              = |K| + 1  (CENSORED) if no rho <= |K| overturns
```

**Every input is in the analyst's hands.** `tau` appears nowhere. This is the single sentence that
separates ρ\* from gadjid / SID / the Taeb–Guo–Henckel radius, all of which are computed around a
**known** G₀ (`RANKING.md:325-331`) — and it is the sentence the pilot asserts and never demonstrates.

**Secondary, continuous, never censored:**  `R(rho, n) = d(rho) / (z * SE_n)` — the
**ignorance-to-sampling ratio**. Reported as a distribution. A four-valued censored statistic is what
the reviewer objected to; `R` has no bins and no ceiling.

⚠️ **Semantics, to be printed in the paper and not softened:** ρ\*_se measures **ignorance, not
error**. It answers *"how few expert mistakes would move my estimate by more than my own standard
error"*, not *"how wrong am I"*. `changed_still_valid = 0` across 27,000+ pilot perturbations means
the two coincide empirically in this ensemble — **empirically, not by theorem**.

---

## 2. Arms

### ARM 1 — population Σ, exact, zero Monte Carlo. *Does `se` un-degenerate the instrument?*

Pilot machinery imported unmodified (`graphs.py`, `adjust.py`, `scm.py`). `MAX_K = 4` — **the
pilot's own cap, deliberately kept**: the claim under test is that the `se` criterion removes the
degeneracy *without* the |K|=8 ensemble that does not exist. Full ball, ρ = 1…|K|.

| ensemble | grid | filter |
|---|---|---|
| `licensed` | p ∈ 5..10, deg ∈ {1.5,2,2.5,4,6} | **none** (`n_und ≥ 3`, the pilot's own) — **PRIMARY** |
| `original` | p ∈ 5..8, deg ∈ {1.5,2,2.5} | none |
| `large` | p ∈ 15..25, deg ∈ {1.5,2,2.5} | none |
| `k8` | licensed grid | `n_und ≥ 8`, nested K4 ⊂ K6 ⊂ K8 — **only to evaluate the registered falsifier as literally written** |

n-grid: **exactly the registered falsifier's** — {200, 500, 1e3, 2e3, 5e3, 1e4, 2e4, ∞} — extended
downward/upward {50, 100, 5e4, 1e5, 1e6} for diagnosis only. z ∈ {1.0, 1.96, 2.576}, 1.96 primary.
Legacy criteria `sign` / `rel50` / `ident` recomputed on the same objects for continuity with E1.

### ARM 2 — finite-sample analyst simulation. *Is the number usable by someone without the answer key?*

The pilot's deepest gap (`CHRIS-MEETING.md:97`): **K is true by construction**, so the
analyst-facing quantity is never exercised. Here the expert is genuinely wrong.

For each SCM: draw `n` rows from the SCM (`scm.sample_linear`), form `Σ̂`, and hand the analyst
**only** `(C, K_asserted, Σ̂, n)`, where `K_asserted` is `K_true` with `r ∈ {0,1,2}` statements
reversed (Meek-consistent ones only; inconsistency is recorded as *caught free*, not dropped).
The analyst reports:

| interval | construction |
|---|---|
| `naive` | `est0 ± z·SE_n` — trusts K completely |
| `robust(ρ)` | `[ min_m (est_m − z·SE_m) , max_m (est_m + z·SE_m) ]` over consistent+amenable members at distance ≤ ρ |

Outcome: **does the interval cover `tau`** (adjudicated with the oracle, as *evaluation*, never as
*input*), and at what **width**. n ∈ {200, 1000, 5000, 20000}; 20 datasets per (SCM, n).

---

## 3. Decision rules — fixed now

### 3a. The registered falsifier, evaluated as literally written

`RANKING.md:346` registers: *"ρ\*(n) rebuilt at |K|=8 against a matched |K|=4 control remains
top-censored for > 75.5 % of SCMs at every n ∈ {200, 500, 1e3, 2e3, 5e3, 1e4, 2e4, ∞}."*
Evaluated verbatim on the `k8` ensemble. **FIRES** iff censoring > 0.755 at every n.

### 3b. The primary rule — a LOCATION, not an attainment

D2's lesson: a bar that can be cleared without the intervention is not a bar. Censoring under `se`
is **monotone non-increasing in n by construction** (§5), so "it goes down" is worthless and
"it reaches 0.60 somewhere" is guaranteed. What is **not** guaranteed is *where*. On the
**`licensed` ensemble, |K| ≤ 4, unconditioned**, define

```
n50  = min { n in grid : censored(n) <= 0.50 }
nND  = min { n in grid : max-bin(n) <= 0.60  AND  >= 3 of the 5 bins hold >= 0.05 }
```

- **SUPPORTED** iff `n50 <= 20000` **and** `nND <= 20000`.
- **REFUTED** iff `n50 > 1e6` **or** `nND` does not exist for any n ≤ 1e6.
- **INCONCLUSIVE** otherwise.

Grounding for 20,000: the top of the registered falsifier's own grid, and an `n` an observational
analyst plausibly has. An instrument that only becomes informative at n = 10⁶ is not an instrument —
which is exactly the objection E1 was answering, displaced onto the new axis.

### 3c. Honest-risk comparator, named in advance (the E3 lesson, applied ex ante)

E3's pre-registered comparator was **mis-operationalised** and the contribution collapsed from a
section to a sentence. So, named now: the free alternative to ρ\*_se is the analyst's own
**t-statistic** `t = |est0| / SE_n`, available at zero cost.

> If `|Spearman( rho*_se(n=2000) , t )| > 0.90`, ρ\*_se is a relabelled t-statistic and the
> contribution is **one sentence, not a leg of the paper**. Reported whatever it comes out at.

### 3d. ARM 2 — the vacuity falsifier

A wide enough interval always covers. So:

- **Capability check (must pass or ARM 2 is unreadable):** at `r = 1`, `naive` coverage of `tau`
  must fall **below 0.90** at n = 20000. If a wrong expert costs nothing, there is nothing to measure.
- **SUPPORTED** iff at `r = 1, ρ = 1`: coverage ≥ 0.90 **and** median width ratio vs `naive` **< 3×**.
- **REFUTED** iff coverage < 0.80 **or** median width ratio > 10×.
- Reported regardless: the full coverage–width frontier over ρ ∈ {0,1,2,…,|K|} × r ∈ {0,1,2}.

---

## 4. Registered prediction on the large-graph inversion (D6)

E1 killed the instrument at p ∈ [15,25] (censored 0.877) and named the cause: **O\*-locality** — a
randomly drawn statement is rarely near (X,Y), so O\* does not move and `d(ρ) = 0` exactly.

**Registered prediction, direction fixed in advance:** the inversion **weakens but does not vanish**
under `se`. Mechanism: exactly-zero deviations stay zero at every n (locality is exact at ρ=1 —
0/654 silent at graph-distance ≥ 1), so they can never exceed `z·SE_n`; but members whose O\* moves
*slightly* now count once `n` is large. Quantitatively:

- **PREDICTED:** `large` censoring at n = 2e4 lands strictly between 0.30 and 0.85, i.e. strictly
  better than E1's 0.877 and strictly worse than `licensed` at the same n.
- **SURPRISE (would refute the mechanism):** `large` censoring at n = 2e4 ≤ 0.30 — locality would
  then not be protecting anything and E1's stated cause is wrong.
- **SURPRISE (the other way):** ≥ 0.85 — `n` buys nothing and the licensed scope boundary hardens
  further.

---

## 5. Declared definitional — NOT evidence, written down before the run

D3's lesson. These are true by construction and will not be quoted as findings:

1. **`rho*_se(n)` is monotone non-increasing in `n`.** `SE_n ∝ n^{-1/2}`, `d(ρ)` is n-free.
   Used as an **integrity gate** (§6 G3), never as a result.
2. **`d(ρ)` is monotone non-decreasing in ρ.** Larger ball, max over a superset.
3. **`sign` ⇒ `rel50`.** A sign flip implies |est−τ| ≥ |τ| > 0.5|τ|, so ρ\*_sign ≥ ρ\*_rel50 pointwise.
4. **ARM 2 at ρ ≥ r is near-guaranteed to cover**: the ball around `K_asserted` then contains
   `K_true`, whose O\* is valid. The content of ARM 2 is therefore the **width cost** and the
   behaviour at **ρ < r**, not the coverage at ρ ≥ r.
5. **`ident` is not nested with either.** Reported separately.

---

## 6. Integrity gates — all must pass or no verdict is readable

| id | gate | pass condition |
|---|---|---|
| G1 | pilot faithfulness: silent-bias rate at ρ=1, |K|≤4, pilot grid | within [0.12, 0.18] (pilot 0.141, E1 re-impl 0.153) |
| G2 | legacy criteria reproduce E1 on the `k8` ensemble | `rel50` K8 censoring within ±0.03 of 0.336 |
| G3 | **monotonicity in n**: `rho*_se(n1) >= rho*_se(n2)` for every SCM, every n1<n2 | **100 %** |
| G4 | **monotonicity in ρ**: `d(ρ)` non-decreasing | **100 %** |
| G5 | **placebo P1 (dead-switch, one direction)**: recompute the ball with `O` frozen at `O0` | censoring = **1.000** at every n |
| G6 | **placebo P2 (dead-switch, other direction)**: `SE_n` forced to 0 (n → ∞) | censoring = fraction of SCMs with `d(|K|) = 0` exactly |
| G7 | ARM 2 sanity: at `r = 0`, `n → large`, `naive` coverage → nominal 0.95 ± 0.03 | pass |

**C13 compliance.** The instrument's dynamic range is demonstrated *within this run*, not borrowed:
G5 forces 100 % censoring and G6 forces the floor, on the identical code path. A null anywhere
between them is a property of the estimand.

---

## 7. Mandatory citations — without these, three over-claims become four

- **Taeb, Guo & Henckel arXiv:2511.10625 Eq. (8)** already bounds structural error against a
  population proxy (`d_L(Ĝ_n,G) ≤ d_L(Ĝ_n,G_∞) + d_L(G_∞,G)`, *variability* + *approximation error*)
  — the **same form**, in graph-distance units. Without it, *"nobody compares structural error to
  sampling error"* is the third avoidable over-claim of this project (`CHRIS-MEETING.md:169-172`).
  🔴 **This paper is not in the corpus** — `check_claims.py` is blind to the nearest competitor.
- **Henckel, Perković & Maathuis arXiv:1907.02435** — O\* and its variance-optimality. The variance
  formula this run uses is *their* estimand; `se` inherits their optimality, it does not add to it.
- **Guo & Perković Corollary 8** — one adjustment set per MPDAG member is **not** a contribution.
- **Grinsztajn-style scope honesty:** everything here is linear-Gaussian iSCM, single X, single Y,
  ER graphs, `p ≤ 25`. The claim is licensed to that box and says so.

## 8. Scope statement, fixed now

**Licensed:** linear-Gaussian iSCM, ER DAGs `p ∈ [5,25]`, expected degree ≤ 6, single-node X and Y,
|K| ≤ 4 true-or-perturbed orientation statements on undirected CPDAG edges, homoskedastic OLS SE,
**CPDAG given** (not discovered). **Not licensed and not claimed:** nonlinear mechanisms, latent
confounding, multi-node interventions, a CPDAG estimated from the same data (`arm C`, NOT RUN — see
`AUDIT.md` A2), non-Gaussian errors, real data.

---

# ADDENDUM 1 — 2026-08-19, written after a 24-SCM smoke run, BEFORE the real ensembles

**Disclosure first.** This amendment was made after looking at a **24-SCM smoke run**
(`results/_smoke.json`, discarded) whose purpose was the integrity gates. All seven gates passed
(G1 silent-bias 0.1354 vs pilot 0.141; G3/G4 monotonicity 0 violations; G5 placebo 1.000 censored at
every n and every z; G6 floor exact; `changed_still_valid = 0`). **One number was also seen**:
censoring under `se` at `n = ∞` was 0.5417, i.e. a hard floor no sample size can beat. That number
is what exposed the flaw below. It is **not** a result and the smoke file is deleted; the amendment
is disclosed rather than back-dated.

## The flaw: the censored bin conflates three different events

`rho*_se` counts only **Meek-consistent AND amenable** members, because only those carry an
estimate to compare. So `CENSORED` silently merges:

| what actually happened | is it "no breakdown"? |
|---|---|
| **(a) every perturbation is Meek-inconsistent** — the expert's error is *caught for free* | ✅ genuinely safe |
| **(b) every perturbation destroys amenability** — identification *fails, loudly* | ❌ **this IS a breakdown**, merely a visible one |
| **(c) consistent, amenable, and `O*` never moves far enough** | ✅ genuinely robust |

Case (b) sits in the "instrument says nothing happened" bin while being the loudest possible
failure. E1's `sign` and `rel50` have the **same defect** — which is precisely why `ident` scores
0.306 against `sign`'s 0.903 on identical data. The 60 pp "dynamic range" E1 cites as C13 compliance
is, read correctly, **evidence that a third of the censored mass is misfiled**.

## The fix — a composite that is still analyst-facing

`ident` needs no `tau`: amenability is read off the MPDAG. So the union is computable by the analyst:

```
rho*_any(n) = min( rho*_se(n) , rho*_ident )
```

*"How few expert mistakes before the analysis stops being the analysis I think it is"* — either the
effect stops being identified (loud) or the estimate moves past my own standard error (silent).

## Amended decision rules

- **PRIMARY statistic → `rho*_any`.** The rules of §3b (`n50 ≤ 20000` and `nND ≤ 20000`) now attach
  to `rho*_any` on the `licensed`, |K| ≤ 4, unconditioned ensemble.
- **`rho*_se` alone is reported as the SECONDARY, isolating statistic**, with its own `n50`/`nND`
  under the identical rule. Both verdicts are printed. If they disagree, **both are printed and the
  disagreement is the finding** — no post-hoc choice of the flattering one.
- **`rho*_ident` is reported alone**, as it is free and already on the record (0.306).
- **New mandatory output — the censored-bin decomposition**, per ensemble and per n:
  the censored mass split into (a) all-caught-free, (b) all-loud, (c) genuinely-robust.
  This is the direct answer to *"your diagnostic returns 'no breakdown' nine times in ten"*:
  the correct reply is not a smaller number, it is **that bin was never one thing**.
- **New honest-risk comparator, added to §3c.** `rho*_any` must not be a relabelling of
  `rho*_ident`, which is free: if `Spearman(rho*_any, rho*_ident) > 0.95` at n = 2000, the `se`
  half contributes nothing and the leg is one sentence. Reported at whatever value it takes.

## What does NOT change

Gates G1–G7, the placebos, the `large`-ensemble prediction of §4, the definitional declarations of
§5, the mandatory citations of §7, the scope of §8, and the literal evaluation of the registered
falsifier in §3a — all stand exactly as written before the smoke run.
