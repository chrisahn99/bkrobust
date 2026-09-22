# The breakdown radius ranks survival within a knowledge state, and barely across networks

*Session 9 report. `[RE-11]`: survival and the paired cross-arm design on the
committed real-network corpus. Branch `experiments/todo_1`. Every number traces
to a committed file under `results/axis_robustness_real/` and is re-asserted by
`src/bkrobust/analysis/session9_verify.py`.*

**Pre-registration:** `results/axis_robustness_real/PREREGISTRATION.md` — §§1–10
written before any survival curve, AUC, τ or cross-arm bin existed, with
directional predictions for seven hypotheses. **Appendices A–J record every
amendment, falsification, defect and self-correction in order**, including three
of my own errors.

**Gate:** `benchmarks.measure.fast_gate` exclusively, on every row.
**Assumption on every radius:** *Conjecture 2 (hence Anti-Exchange Case B,
verified not proved)*; the error is one-sided, so radii can only be too **large**
(`THEOREMS.md` §4c, §6, §8).

---

## 1. The headline, in the order a reviewer should read it

**Anchor table: [`table_tau_real.md`](table_tau_real.md).**

| | result |
|---|---|
| **Pre-registered, primary** | τ_b(`r_hop`, survival) is **positive in 9 of 9** strata. Its interval excludes zero in **3 of 9** — against **8 of 9** synthetically. |
| **After a stronger leave-one-network-out check** | **Exactly 1 of 9** strata (`tiered, n_tiers=4`) keeps an interval excluding zero *and* survives dropping any single network (Appendix H). |
| **Post-hoc, and the one that means something** | **Within a fixed network and knowledge state**, the radius ranks survival with an interval excluding zero in **6 of the 9** strata, at mean τ **+0.30 to +0.84** (Appendix J). |
| **Against the baselines** | `r_hop` outranks `\|K\|` in only **5 of 9**. `SHD(G₀, truth)` beats it in **5 of the 8** strata where SHD is defined. No dominance claim is made or supported. |

**The pre-registered result is materially weaker than the synthetic one, and that
is the honest answer to "does this transfer to practice".** The post-hoc
decomposition explains why, and is the finding this campaign actually contributes:
**the pooled τ mixes a within-network comparison, which is what the radius is for,
with a between-network one in which the radius is confounded with `|K|`,
component structure and a measurement artifact.**

Two further results are sharper than the τ table and cost the paper nothing.

**The paper's Result D does not transfer, and the cause is the knowledge model,
not the graphs.** The single-claim contradiction rate on this corpus is **0.021**
against a synthetic 0.627. Holding the same 25 networks fixed and changing only
how the analyst's knowledge is written down moves it to **0.432** — a factor of
twenty (§5).

**`paths` is a clean, exhaustively enumerated demonstration of the two-radii
over-promise.** One statement, fourteen committed orientations, `SHD = 0`, radii
spanning 1–14, survival **exactly 0.000** on all twenty pairs (§6).

---

## 2. What would have falsified this, and what did

Seven directional predictions were recorded before any curve existed. **Four are
falsified**, one exactly right.

| | prediction | outcome |
|---|---|---|
| **P1** | τ>0 in 6–9/9; CI excludes 0 in **4–7/9** | 9/9 positive ✔; **3/9** ✘ |
| **P2** | `r_hop` beats `\|K\|` in **6–9/9** | **5/9** ✘ |
| **P3** | `\|K\|` CI excludes 0 in **≥3/9** | **2/9** ✘ |
| **P4** | SHD undefined in **exactly 1** stratum | **exactly 1** ✔ |
| **P5** | contradiction **0.65–0.90** at depth 1 | **0.021** ✘✘ |
| **P6** | cross-arm ordering replicates | direction ✔, band far narrower |
| **P7** | undefined-separation: higher AUC, **lower** r=1 share | AUC ≈ equal; r=1 share **0.82 vs 0.43** ✘ |

**Falsification triggers.** T1 (≥5 of 9 positive) passed. T3 passed cleanly — see
§7. T4, T6 and T7 passed. T5 passed on its letter only — see §8. No trigger was
weakened after a result was seen.

---

## 3. Design

### 3.1 The frame

The **primary frame** is the committed real-network corpus: every row of
`results/axisa3/instances.jsonl` with `admissible == true` — **831 rows, 25
networks, 543 distinct `(network, X, Y)` pairs**, at coverage 1.0 / 0.5 / 0.25
(543 / 182 / 106). It was built, frozen and hashed **before any survival curve
existed** (`results/axis_robustness_real/frame.jsonl`, sha256
`73881798319206be…`; source `instances.jsonl` sha256 `fc10ad9707fe1229…`), and
every later stage re-checks both hashes before writing a row.

**Falsification trigger T7 passed**: all **831 of 831** radii recompute to exactly
the committed values, with the same dispatch leg, the same commitment size and
the same separation status. The orchestrator spot-checked six random rows against
an independent recomputation rather than taking the builder's word.

Every knowledge-dependent condition is a **status column on that frame, never a
filter**. Nothing is dropped; `optimal_set_undefined`, `n_k_zero`,
`o_g0_extensions_intractable`, `k_contradictory` and `censored_wall_cap` are all
statuses and all counted.

**The two arms have structurally different denominators, and this is stated
rather than smoothed over.** The flip arm is a function of coverage, so its
denominator is the 831 frame rows. The tiered arm generates its own knowledge
from the network's temporal order and has no coverage parameter at all, so its
denominator is the 543 distinct pairs.

### 3.2 The unit of analysis is the network

The intraclass correlation of the `r = 1` indicator is 0.40 over the 25 networks
— a design effect of about 12.7 against 831 rows. Accordingly **every interval in
this report is a cluster bootstrap over networks**, 10,000 resamples, seed 0,
each resample drawing 25 networks with replacement and taking all rows of each.
**No per-row interval is printed anywhere.** Every aggregate is given both
pair-weighted and network-weighted, and the per-network table is printed beside
it.

Because within-network radius variation on this corpus is thin — at coverage 1.0,
ten of the twenty-five networks have a single distinct radius across all their
rows — every headline τ also carries a **leave-one-network-out** sensitivity, and
any cell whose verdict does not survive dropping one network is marked
`⚠ single-network leverage` and is not quoted alone.

### 3.3 The SHD baseline, and the base wrongness it needs

`select_knowledge` draws the analyst's claims from the true orientations, so at
coverage 1.0 `G₀` **is** the ground-truth DAG and `SHD(G₀, truth) ≡ 0`; at partial
coverage it equals exactly the number of undirected edges `G₀` still carries — it
measures *incompleteness*, never *wrongness*. Without intervention the baseline
would be rigged by construction.

The synthetic design's **base wrongness** is therefore reproduced:
`b ∈ {0.00, 0.10, 0.25}` applied by `synth.knowledge.flip` before the sweep, with
**three analyst replicates** at each non-zero level (one at `b = 0`, where the
corruption is the identity). Seeds derive by SHA-256 from
`(network, coverage, "base", b, replicate)`.

**Appendix A of the pre-registration records a defect in that mechanism on this
corpus**, found before any curve existed: because `|K|` is 1 or 2 on most of the
50 `(network, coverage)` cells, `round(b · |K|) = 0` and the rate-based
corruption reverses **nothing at all** — in 38 of 50 cells at `b = 0.10` and 24 of
50 at `b = 0.25`. A fourth, **absolute** level (`bw_abs = 1`, exactly one claim
reversed, expressed through the same `flip` call at rate `1/|K|`) was added as a
**supplementary** stratum so the SHD baseline is non-degenerate on every network.
The nine pre-registered strata are unchanged. The inert cells were run in full
rather than skipped, and every row carries `n_claims_actually_wrong` and
`bw_is_inert`.

### 3.4 The two arms

**Flip** — axis is the **claim depth**. `K_b` is Meek-closed to `G₀`,
`Z* = optimal_adjustment_set_mpdag(G₀, X, Y)` is read **once** and held fixed, and
`flip` reverses exactly `d` claims. The realised reversal count is **asserted**
equal to `d`, never assumed. The depth grid is exactly the set the pre-registered
`AUC_frac` endpoint averages over.

**Tiered** — axis is the **corruption rate**, from the start, never the `d_claims`
column. That column counts assertions added or reversed but **not removed**, and
relocating a node into a shared tier silently deletes a cross-tier assertion, so
its zero bin is contaminated (synthetic `PREREGISTRATION.md` Appendix E). It is
still recorded, under the name `d_claims_defective_do_not_bin`, purely so the
defect can be re-measured; no analysis bins on it.

**The two arms are never merged row to row.**

### 3.5 Endpoint, draws and censoring

`Z*` is read once at depth 0 and held fixed; survival is decided by the
**polynomial GAC predicate** `is_gac_valid_mpdag`, and
`all_valid_adjustment_sets_mpdag` is never on a decision path anywhere in this
campaign.

**N = 1000 draws per grid point, everywhere**, with no cell at lower precision.
Both the raw endpoint and the conservative `_usable` variant (grid points with
`n_eval ≥ 30`) are computed, both are reported, and the conservative one is
quoted in summary claims.

`UNREACHED` is a status, never a radius — there are none in this data. A
wall-capped unit carries `wall_until_timeout_s` and **no row this campaign writes
has a key named `seconds`**, which is checked mechanically by the writer itself.
A contradictory corrupted `K` is `corrupted_k_contradictory` and is not a
survival-0 outcome; `S` is computed over non-contradictory draws and the
contradiction rate is its own series.


---

## 4. Why the pooled figure is weak, and what is underneath it

Three facts about this corpus, all measured, account for the gap between the
synthetic result and this one. **None of them is a property of the radius.**

### 4.1 `|K|` and `k_g0` are network-level constants

`select_knowledge` returns one claim set per `(network, coverage)`, and `flip`
preserves list length, so no base wrongness can make `|K|` vary within a network.
The `|K|` and `k_g0` baselines can therefore rank *networks* but not *queries* on
this corpus — a weaker thing than they could do synthetically. A win over them
here is worth correspondingly less, and the table says so.

### 4.2 The depth grid pins the endpoint at zero on small-`|K|` networks

The flip arm's depth grid is **fractional**. At `|K| = 1` the only corruption it
can draw is "reverse everything", so those networks' survival endpoint is
**exactly 0 by construction**, whatever their radius. Five of the 24 networks at
coverage 1.0 are pinned — `Acid_1996`, `Didelez_2010`, `hailfinder`, `mediator`,
`paths` — and **all five have `|K| ≤ 3`**.

Removing them raises τ_b(`r_hop`) in the three coverage-1.0 flip strata to
**+0.286**, **+0.376** and **+0.583**, two of three intervals excluding zero.
This is a **design–corpus interaction**, not a property of the radius, and it is
reported as such (Appendix F).

At the extreme, `pathfinder` (`|K| = 79`) has a grid whose *smallest* point is
`d = 8`; eight simultaneous reversals on its 85-vertex chain component always
contradict, so it yields **no survival measurement at all** at any depth — a
measurement limit recorded as a status (Appendix I).

### 4.3 The pooled τ is the wrong statistic for the stated unit of analysis

This is my design error, not a discovery. §1.5 fixed the network as the unit of
analysis and bootstrapped over networks, but computed τ pooled over rows from
different networks. Decomposed (Appendix J):

| stratum | pooled τ | **within-network** mean τ [95% CI] | nets pos. | between-network τ |
|---|---|---|---|---|
| flip cov=1.0 bw=0.00 | +0.064 | **+0.296 [+0.016, +0.561]** | 7 / 11 | −0.253 |
| flip cov=1.0 bw=0.10 | +0.126 | **+0.303 [+0.185, +0.412]** | 14 / 15 | −0.144 |
| flip cov=1.0 bw=0.25 | +0.292 | **+0.403 [+0.295, +0.511]** | 17 / 18 | −0.071 |
| flip cov=0.5 bw=0.00 | +0.191 | +0.303 [−0.361, +0.865] | 4 / 5 | −0.036 |
| tiered n_tiers=2 | +0.580 | **+0.836 [+0.615, +1.000]** | 3 / 3 | +0.526 |
| tiered n_tiers=3 | +0.378 | **+0.654 [+0.375, +0.863]** | 7 / 8 | +0.324 |
| tiered n_tiers=4 | +0.483 | **+0.631 [+0.534, +0.717]** | 9 / 9 | +0.452 |

> **Within a fixed network and knowledge state — the comparison a practitioner
> actually faces — the radius ranks survival positively and consistently. Pooled
> across networks it is confounded, and the pooled figure is a mixture.**

This also disposes of an apparent reversal. The supplementary coverage-0.25
strata show a *significant negative* pooled τ (−0.527, −0.513, −0.521, intervals
excluding zero). Decomposed that is **−0.775 between networks and +0.158 within
them**, produced by three `|K| = 1` networks whose endpoint is pinned at 0 and
which happen to carry radius 2. **It is not evidence that the radius ranks
survival backwards** and is not reported as such.

---

## 5. Errors on real structure are not self-revealing — and that is about the knowledge, not the graphs

The paper's Result D says most orientation errors announce themselves: Meek
closure fails about six times in ten. **On this corpus it essentially never
does.**

Measured **exhaustively** — every claim reversed in turn, no sampling, so these
are counts and not estimates:

| | single reversals | contradictory | rate | cells at exactly 0 |
|---|---|---|---|---|
| all 50 `(network, coverage)` cells | 298 | **5** | **0.017** | 46 of 50 |
| coverage 1.0, minimal generator | 238 | **5** | **0.021** | 21 of 25 |
| coverage 1.0, **all undirected edges asserted** | 410 | **177** | **0.432** | 5 of 25 |

The third row is the synthetic arm's knowledge model, on the **same graphs**.
Holding structure fixed and changing only how the knowledge is written down moves
the rate by a factor of twenty. The internal control is clean: on the five
networks where the minimal generator *is* the full assertion set, both models give
**exactly 0.000**, as they must.

> **Result D is a property of how the analyst's knowledge is phrased, not of the
> graphs.** The paper must say so.

This was pre-registered as a diagnostic **with its directional prediction, before
it ran** (Appendix D.3 predicted "≥ 0.30 for `K_full` against ≈ 0.01 for `K_min`";
the outcome is 0.432 against 0.021). It raises `[RE-1]` from a tidying-up item to
a load-bearing one: whether elicited knowledge is Meek-closed decides which of two
numbers differing twentyfold the paper may print.

The rate does climb steeply with depth — 0.027, 0.173, 0.431, 0.565, 0.662, 0.800
at `d = 1…6` at coverage 1.0 — so the arm is not degenerate; it starts from a
floor near zero.

---

## 6. The spotlight: `paths`

Twenty admissible pairs, flip arm, truthful knowledge. Checked **before** it was
written up, per the session-8 lesson.

| quantity | value |
|---|---|
| stated knowledge `\|K\|` | **1** |
| commitment size `k_g0` | **14** (leverage 14.0) |
| `SHD(G₀, truth)` | **0** — `G₀` *is* the ground-truth DAG |
| `r_hop` across the twenty pairs | **1 … 14** |
| survival AUC | **0.000 on every one** |
| evaluable draws per row | 1000, with **zero** contradictory |
| distinct corrupted states | **1** — the corruption space at `\|K\|=1, d=1` has one element |

That last row matters and is stated rather than buried: this is an **exhaustive
enumeration**, not a sample. The result is a deterministic fact.

The analyst made one statement. Meek closure turned it into fourteen committed
orientations. SHD says they are perfectly correct. The hop radius certifies up to
thirteen clean shells. **Reversing the single statement they actually made
destroys the adjustment set every time, at any radius.**

The radius is not wrong: the nearest failing knowledge state really is fourteen
hops away, because reversing a claim is not a one-hop move. What over-promises is
the **practitioner reading** — the certificate denominated in stated claims,
exactly the defect `[RE-16]` names and `docs/PAPER_NARRATIVE.md` §4 warns about.

**A coda worth the paper's attention.** The same network behaves in *opposite*
directions in the two arms: survival **0.000** at radii up to 14 in the flip arm,
**1.000** at high radii in the tiered arm. Not a contradiction — the two arms hand
`paths` completely different knowledge (`|K| = 1` versus a generated tiering).
That is this campaign's central finding restated on a single network.

---

## 7. Precision: N = 1000 is exhaustive for most of this corpus

**Every cell was run at N = 1000 draws per grid point**, with none at lower
precision. But a draw count is not an amount of evidence, and on this corpus the
two differ sharply.

The closure cache records, for free, how many **distinct** corrupted claim sets
each shard actually saw. Against the combinatorial size of the flip arm's
corruption space, `Σ_d C(|K|, d)` over the swept depths:

> **409 of 450 scored flip shards are *exhaustive*** — they enumerated the entire
> corruption space, so their survival numbers carry **no Monte-Carlo error at
> all**. Only 41 shards, on seven networks, are genuinely sampling.

`insurance` at `|K| = 3` sees all 7 possible corruptions; `barley` at `|K| = 7`
sees all 127. Every shard with `|K| ≤ 9` is exhaustive and every sampled shard has
`|K| ≥ 9`. Only `arth150`, `diabetes`, `ecoli70`, `hepar2`, `magic-irri`,
`pathfinder` and `win95pts` sample, the thinnest — `pathfinder`, at `|K| = 79` —
covering 1.1e-19 of its space.

**Both endpoints agree, which discharges the concern behind `[RE-12]` on this
corpus.** Across 45 (stratum, predictor) cells the raw and `n_eval ≥ 30`
endpoints disagree on a verdict **zero times**, with deltas at or below 0.02 in
the flip arm and a maximum of 0.108 in the thinnest tiered stratum. At N = 1000
the filter is **inert** — exactly the diagnostic the synthetic Appendix J
identified and could not run. Falsification trigger **T3 passes cleanly**.

`[RE-12]` itself is **not** discharged: it asks for a re-run of the *synthetic*
nine strata, which this is not.

---

## 8. The paired cross-arm: it works, and it barely reaches

The design is session 8's (synthetic Appendix I), ported intact: every instance
built **once** the tiered way, corrupted **both** ways, scored on the shared
intensity axis `|dir(G₀) Δ dir(G)| / |dir(G₀)|`, paired within instance.

**The pairing itself is sound.** All **501** instances common to both arms have
**identical** `instance_sha256` over `(K_ref, dir(G₀), Z*)` — zero mismatches — so
every pair really is one instance corrupted two ways. The two arms were sampled
in separate shards and the hash is what guarantees it, checked rather than
assumed.

**But it barely reaches.** Only **three** intensity bins carry ≥ 15 matched
instances in both arms, and they are **0.00, 0.05 and 0.10** — against the
synthetic band of 0.15–0.65. The reason is structural: real `|dir(G₀)|` is large,
so the normalised symmetric difference is small for any realistic corruption and
both arms pile up near zero. Those three bins rest on **8, 8 and 3** networks.

Falsification trigger **T5 passes on its letter** (≥ 3 bins) and the result is
reported with that qualification rather than as a clean port.

In the three bins that carry data, the direction of session 8's finding
replicates: tiered contradiction exceeds flip in all three (+0.233, +0.248,
+0.141), and the silent-failure rate is lower for tiered in the two
best-supported bins (−0.062, −0.102). **But under the network-clustered bootstrap
this design mandates, only one of those differences has an interval excluding
zero.** The per-instance sign tests are far more confident than the clustered
intervals — p-values to 1e-26 against intervals spanning zero — and the clustered
interval is the one that respects the stated unit of analysis.

> **The cross-arm ordering is directionally consistent with the synthetic result
> and is not established at the network level on real structure.**

---

## 9. The matched-population coverage contrast ([RE-6])

Conditioned on the **105 `(network, X, Y)` triples admissible at every coverage
level**, on 8 networks, so the composition effect that makes the unmatched
543 → 182 → 106 comparison meaningless is removed. Non-nesting of the coverage
sweep ([RE-4]) remains and is reported, not fixed.

On the same triples, raising coverage 0.25 → 1.0:

| | 0.25 | 0.5 | 1.0 |
|---|---|---|---|
| median survival AUC | 0.347 | 0.428 | **0.569** |
| median `\|K\|` | 2 | 4 | 7 |
| median `k_g0` | 4 | 6 | 10 |
| median leverage | 2.0 | 1.5 | 1.29 |
| median `r_hop` | 2 | 2 | 2 |

Survival rises on **8 of 8 networks** (network-level sign test p = 0.0078, and
underpowered at 8 clusters). The radius is **unchanged on 73 of the 105 triples**
and rises by exactly one, 2 → 3, on the other 32; its median rises on **1 of 8**
networks.

> **More knowledge coverage buys measurably more average-case survival. The radius
> moves in the right direction where it moves at all, but is far less responsive
> to coverage than survival is.**

Note this runs **opposite** to the synthetic association between `|K|` and
survival, which was negative. The two are not comparable designs — synthetic
varied `|K|` across different instances, this varies coverage within the same
instance — and both are stated rather than reconciled.

---

## 10. The undefined-separation stratum (P7)

368 of the 831 frame rows have no member of `Z*` in `X`'s component, so separation
is undefined. Reported on its own, never as a number.

Pooled over the six primary flip strata:

| | n | networks | median AUC (pair) | median AUC (network) | share at `r_hop = 1` |
|---|---|---|---|---|---|
| separation measured | 2599 | 20 | 0.405 | 0.408 | **0.430** |
| separation undefined | 2149 | 14 | 0.428 | 0.411 | **0.820** |

**P7 is half falsified.** The predicted higher AUC holds only pair-weighted and
vanishes network-weighted (0.411 against 0.408); in the tiered arm it *reverses*
(0.777 against 0.889 at `n_tiers=4`). The predicted **lower** share at `r_hop = 1`
is wrong in the opposite direction, and decisively: **0.820 against 0.430**.
Instances whose adjustment set lies outside the treatment's component are far
*more* often at radius 1, not less.
## 11. What is not claimed

- **Nothing is pooled with the synthetic corpus.** Not in a table, not in a
  figure, not in a sentence. Where a synthetic figure is named for contrast, the
  two are separately computed aggregates on different corpora with different
  knowledge models, and the comparison is qualitative.
- **No claim of perfect ranking**, and no claim of dominance. `r_hop` does not
  strictly dominate the baselines on this corpus and this report does not say it
  does.
- **The leverage analysis of §6 is post-hoc.** It was opened after seeing the
  primary result, it is labelled as such, and the nine pre-registered strata are
  not restated in its light.
- **Oracle CPDAG throughout.** Every CPDAG here is computed from the true DAG.
  Real discovery on finite data returns a sparser skeleton, and missing weak
  edges is exactly what breaks back-door blocking. This is the corpus's largest
  scope limit and nothing here touches it.
- **One corpus.** Everything traces to a single pgmpy 1.0.0 sdist.
- **No comparison to the truth on the computation path.** The endpoint is
  structural survival of the committed set across knowledge states, not whether
  `Z*` is correct in the true DAG. The true DAG enters only to define the
  corpus's knowledge and to score the `shd_truth` baseline.
- **The cross-arm comparison covers only the intensity range where both arms
  carry data**, which on real structure is far narrower than it was
  synthetically. Outside it no claim is made.
- **`[RE-12]` is not discharged.** That item asks for a re-run of the *synthetic*
  nine strata at N = 1000; this campaign is on real structure. What it does
  supply is independent evidence for Appendix J's diagnostic — see §5.


---

## 12. The bugs, reported at the same weight as the results

Four defects and one operator error. **Every one of them returned plausible
numbers rather than crashing** — which is this project's standing lesson and the
reason the cross-checks exist. Three are mine. All are in dated appendices, in
order, rather than quietly fixed.

**1. `g0_undirected_edges` was never a measurement.** `benchmarks/measure.py::evaluate`
assigns it only on the `o_g0_extensions_intractable` rejection path, so **all 831
committed rows carry the dataclass default `0`**. Recomputation gives non-zero on
**288** of them. This is `[RE-5b]`, now confirmed with a number. *Caught by*
recomputing the column while building the frame instead of copying it. The
committed file is **not** modified; the frame carries the recomputed column under
its own name.

**2. A measurement limit was reported as a structural property.** My cross-arm
driver stamped `dir_g0_empty` — a fact about a graph that exists and has no
directed edge — when the real reason was `o_g0_extensions_intractable`, a limit on
what we could afford to evaluate. No measurement was affected; only the label.
*Caught by* reading the `pathfinder` markers by hand after those shards finished
in under a second against a budget of hours. Six shards re-run (Appendix B.1).

**3. The contradiction table counted the same draw once per pair.** The corrupted
claim set is drawn once per `(shard, grid point)` and scored against every pair of
a network, so pooling over cells inflated `n_draws` from a true **22,000 to
445,000** and moved the rate from **0.0270 to 0.0117**. **The defect was in my
specification, not the implementation** — the brief said "pooled over units",
which is right about the arithmetic and wrong about the unit. No τ, endpoint or
survival number is affected. *Caught by* an independent recomputation that had
deduplicated by draw-set from the start and therefore disagreed; each number was
individually plausible (Appendix G).

**4. My own leverage explanation, withdrawn.** I explained the weak flip-arm
result by leverage `k_g0/|K|`. Controlling for `|K|` the pattern **reverses**, the
high-leverage bands rest on 3–4 networks, and τ_b(`|K|`, leverage) = −0.395 across
networks, so the two were entangled from the start. The explanation that survives
is the small-`|K|` pinning of §4.2. *Caught by* running the control (Appendix F).
It was wrong in the direction that made the story tidier, which is when a claim
needs the most scepticism.

**5. My pre-registered leave-one-network-out flag is too weak.** It was defined on
point estimates, so it cannot see a stratum whose interval excludes zero with all
networks and includes zero without one while every point estimate stays positive.
All three tiered strata pass it; two fail the bootstrap version. Reported both
ways, the pre-registered column unchanged (Appendix H).

**Incident.** An orchestrator error ran cross-arm shards in a serial loop
alongside the eight-worker pool, so one shard could have been written twice.
Caught within four minutes; **all 64 shards touched in that window were deleted
and requeued without anyone judging whether a file was salvageable**, because the
incrementality contract removes that judgement. Cost: four minutes of one
network's compute (Appendix B.2).

**Not a bug, but recorded:** the per-shard wall cap that §5.4 pre-registered was
implemented an hour into the run. The eight shards already in flight ran uncapped
and the four `pathfinder` ones exceeded the 5-hour figure; they ran to completion
rather than being censored, which yields *more* data than the cap would have, and
their realised elapsed times are in their completion markers (Appendix C).

---

## 13. Figures

All in `figures/`, PDF and PNG, regenerable by
`PYTHONPATH=src /usr/bin/python3 -m bkrobust.analysis.session9_figures`. Each
carries a one-line "take from this" inside the figure, and every figure showing a
radius carries the Conjecture 2 note.

| figure | take from this |
|---|---|
| `s9_f1_tau_forest` | Three predictors' rank correlation per primary stratum. Open markers are leave-one-network-out unstable; predictor-constant cells are **text, not a point at zero**. |
| `s9_f2_survival_by_radius_bucket` | Survival against the grid point by radius bucket, **flip and tiered in separate panels that never share an axis** — the two arms are never merged. |
| `s9_f3_per_network` | Networks differ from each other more than the corpus differs from any average. `pathfinder` is excluded and the exclusion is counted in the footnote, because every one of its units is unresolved. |
| `s9_f4_leverage` | Per-network leverage `k_g0/\|K\|` against median endpoint, log-x, with the `leverage = 1` line marked "K is Meek-closed". `paths` and `hailfinder` highlighted. |
| `s9_f5_knowledge_model` | The same networks under two knowledge models, paired. The five networks where the minimal generator *is* the full set sit on top of each other at zero — the internal control. |
| `s9_f6_crossarm` | Where each corruption process actually puts its mass, and the paired difference per intensity bin with its **bootstrap CI over networks**. Bins with fewer than 15 matched instances are marked as not carrying the comparison. |

**Read `s9_f4` with §4.2 in hand.** It invites a conflation the report does not
make: a median AUC of exactly 0 tracks **small `|K|`**, not leverage — all five
pinned networks have `|K| ≤ 3` while their leverage spans 1.0 to 17.0, and `child`
at leverage 6.0 is not pinned. `results/axis_robustness_real/SMALL_K_VS_LEVERAGE.txt`
keeps the two apart.

---

## 14. What the paper can now say

See [`docs/PAPER_NOTE_RE11.md`](docs/PAPER_NOTE_RE11.md) for the sentences and the
TODO this discharges.

## 15. Reproducing this

```bash
cd bkrobust
export PYTHONPATH=src

# 1. the frozen frame (and falsification trigger T7)
/usr/bin/python3 -m bkrobust.robustness.real_frame

# 2. the sweep -- 725 shards, resumable, safe to re-run
/usr/bin/python3 -m bkrobust.robustness.run_real_survival pool --workers 8 --n-draws 1000
/usr/bin/python3 -m bkrobust.robustness.run_real_survival status   # must read 725

# 3. analysis, the anchor table, and the supplementary measurements
/usr/bin/python3 -m bkrobust.robustness.real_analyse --n-boot 10000
/usr/bin/python3 -m bkrobust.robustness.build_table_tau_real
/usr/bin/python3 -m bkrobust.robustness.real_knowledge_model
/usr/bin/python3 -m bkrobust.robustness.real_matched_coverage
/usr/bin/python3 -m bkrobust.robustness.real_spotlight

# 4. figures
/usr/bin/python3 -m bkrobust.analysis.session9_figures

# 5. independent verification of every number in this report
/usr/bin/python3 -m bkrobust.analysis.session9_verify
```

The sweep is **idempotent**: a shard with a completion marker is skipped, and a
shard without one is redone from scratch. Re-running `pool` after it exits is
required once, and harmless thereafter.

## 16. Environment and seeds

`/usr/bin/python3`, **Python 3.9.6**, macOS-26.6.2-arm64, numpy 2.0.2,
scipy 1.13.1, pandas 2.3.3, matplotlib 3.9.4, ortools 9.15.6755 — the same
interpreter and library set recorded in `results/axisa3/manifest.json`, which
produced the corpus this campaign runs on. `docs/REMAINING_EXPERIMENTS.md`
recommends standing up Python 3.12; that was **declined deliberately**, because
moving off the interpreter that produced the committed corpus would trade
reproducibility for nothing. 269 of 269 of the tests that matter (`tests/core`,
`tests/search`, `tests/criterion`, `tests/hybrid`, `tests/demo`) pass on it.

`num_workers: 1`, `random_seed: 0`, CP-SAT single-threaded, **no global RNG
anywhere**. Every draw goes through an explicit `numpy.random.Generator` seeded
by a SHA-256 derivation from `(instance_id, process, grid_point, rep)` — never
Python's salted built-in `hash()`. Bit-identical output across `PYTHONHASHSEED`
0 and 12345 is **checked, not assumed**, on every shard kind, excluding the three
wall-clock timing fields; the `*.cells.jsonl` files carry no timing at all and
are byte-identical.

**Assumption carried by every radius in this report:** *Conjecture 2 (hence
Anti-Exchange Case B, verified not proved)*. All upward searches are exact iff it
holds, and the error is **one-sided** — radii can only be too **large**, never
too small (`THEOREMS.md` §4c, §6, §8).

**Gate:** `benchmarks.measure.fast_gate` exclusively, stamped on every row.
`synth.runner.gate` is never imported; it calls the exponential
`all_valid_adjustment_sets_mpdag` and is not equivalent (16 disagreements in
46,800 cases, all one direction).
