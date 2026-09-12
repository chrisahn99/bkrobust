# Paper narrative — the v0 brief

**Audience:** a Claude session (or a person) sitting down to draft the ICLR 2027
paper with full access to this repository.

**What this document is for.** The repository contains eight sessions of work.
Most of it must not go in the paper. The failure mode this document exists to
prevent is an exhaustive recital of everything measured, with no argument. This
file fixes **the argument**, names **the results that serve it**, and names
**the results that must be set aside** even though they were expensive.

Read this before writing prose. Then read
[`docs/REMAINING_EXPERIMENTS.md`](REMAINING_EXPERIMENTS.md) for what is not yet
measured — several sections below are marked `TODO [RE-n]` and must be drafted
as placeholders, not invented.

**Two hard rules.**
1. **Never invent a number.** Every quantity in the paper traces to a file under
   `results/`. The provenance map is §10. If a number is not there, it is a
   `TODO`, not a guess.
2. **Read [`src/bkrobust/core/conventions.py`](../src/bkrobust/core/conventions.py)
   before quoting any radius.** `r` is the distance to the *nearest* failure;
   the safe-moves number is `r − 1`. Anyone quoting safe moves must subtract one
   and say so.

---

## 1. The thesis, in one sentence

> When an analyst orients a CPDAG's undirected edges with background knowledge,
> the adjustment set they read off the result inherits a quantifiable amount of
> structural risk — and that risk can be priced, in the analyst's own units,
> without ever consulting the truth.

Everything in the paper is either the construction of that price, evidence that
it is meaningful, or an honest account of what it does not cover.

---

## 2. Scope — write this early and defend it

This is the single most important framing decision in the paper, and the one most
likely to draw fire if left implicit. State it in the introduction and restate it
in the experimental section.

**We are in one specific, common, under-examined setting:**

- A CPDAG `Ĉ` has already been obtained (by discovery, or assumed).
- An analyst supplies **background knowledge** `K` — orientations of edges the
  data left undirected.
- Meek closure gives the analyst's graph `G₀ = Meek(Ĉ, K)`.
- From `G₀` they read an adjustment set `Z` and estimate a total effect.

**Background knowledge is used *post-CPDAG*, to orient undirected edges.** We are
not studying knowledge injected into discovery, not studying constraint-based
discovery itself, not studying latent confounding or selection. Say this
plainly; it is a feature, not a limitation, because it isolates the one place
where unverifiable human input enters a pipeline that is otherwise data-driven.

**We do not compare to the ground truth.** This must be unmistakable. The
quantity we compute is a function of `(Ĉ, K, X, Y, Σ̂, n)` alone. The true DAG
is used in exactly two places, both after the fact:
1. to draw samples in the simulation panels, and
2. to *audit*, once every reported quantity is frozen and hashed, whether the
   certificate held.

It is **never on the computation path**. Enforce this mechanically: carry a hash
of the CPDAG source and of the generating DAG on every row, and open the audit
only after the deployment rows are frozen (`[RE-10]`).

**What we diagnose.** Not whether the analyst is right. Not whether the estimate
is close to the truth. We diagnose **whether the conclusion survives stress on
the knowledge that produced it.**

### The definition of robustness, stated for our context

Write this as a numbered definition, early:

> **Robustness.** An adjustment set `Z`, read once from `G₀` and then held fixed,
> is *robust to degree r* for the query `(X, Y)` if every knowledge state within
> distance `r − 1` of `G₀` still admits `Z` as a valid adjustment set.

Three things this definition deliberately is **not**:
- It is not accuracy. `Z` may be robust and the analyst still wrong about the world.
- It is not efficiency. We show (§5, Result C) that efficiency is *invariant*
  under truthful knowledge, so robustness is the only axis knowledge moves.
- It is not a property of the graph. It is a property of **the conclusion**, under
  perturbation of its input.

---

## 3. The introduction — specification and draft

The introduction must do four things, in this order. A draft follows; treat it as
a starting point with the right shape and the right commitments, not as final prose.

**(i)** Establish the setting: BK applied post-CPDAG to orient undirected edges.
**(ii)** Establish the problem: consistency is checkable, truth is not.
**(iii)** Establish what we offer: a *diagnostic*, not an optimisation target.
**(iv)** Establish what we do not claim: no comparison to ground truth.

### The sentence that fixes the register

The framing the whole paper hangs on, to be stated explicitly near the end of the
introduction:

> Our radius is not an optimisation target. It is a **diagnostic pricing
> mechanism**: it tells a practitioner exactly how much structural risk they take
> on by using the adjustment set they have chosen.

Do not let the paper drift into "maximise the radius." An analyst who maximises
the radius by asserting more is, by Result C, buying nothing but risk. The radius
is a *price tag*, read before purchase.

### Draft introduction

> Causal discovery returns an equivalence class, not a graph. A CPDAG leaves some
> edges undirected, and in applied work the gap is closed by an expert: the
> analyst asserts orientations for edges the data cannot orient, Meek's rules
> propagate those assertions, and an adjustment set is read off the resulting
> MPDAG.
>
> That input is checked for **consistency** and never for **truth**. Meek's
> algorithm either terminates with a valid MPDAG or reports failure; it has no
> way to ask whether the orientation asserted is the one that holds in the world,
> and it cannot, since answering that would require the graph discovery was run
> to find. Knowledge that is *consistent but false* therefore passes every check
> in current practice. The three-node chain is the whole problem in miniature:
> `A − B − C` is equally the CPDAG of `A → B → C` and of `C → B → A`, so an
> expert asserting the reverse of the truth is undetectably, completely wrong.
>
> We do not propose to detect such errors — in this setting that is impossible.
> We propose to **price** them. Fix the analyst's graph `G₀` and the adjustment
> set `Z` they read from it. Perturb the knowledge, one elementary revision at a
> time, in every direction, and record the distance at which `Z` first stops
> being a valid adjustment set for the query. That distance is a **breakdown
> radius**, and it is computed from the analyst's own inputs — the estimated
> CPDAG, their stated knowledge, the query, and the sample — with no appeal to
> the true graph.
>
> This makes the quantity a diagnostic rather than an objective. It does not tell
> the analyst they are right, and it offers nothing to maximise: we show that
> under truthful knowledge, statistical efficiency is invariant to what is
> asserted, so every marginal claim beyond the minimum needed for
> identifiability buys precision of exactly zero while measurably reducing
> robustness. What the radius supplies is the price: *how much structural risk am
> I taking on to use this adjustment set?*
>
> We report the price in two units, because they differ and the difference
> matters in practice. The **geometric radius** is a distance in the space of
> knowledge states — the natural object, on which the quantity is a genuine
> metric [TODO: cite Guo et al., `[RE-18]`]. The **claim radius** is denominated
> in the sentences the analyst actually uttered. One asserted sentence can
> commit an analyst to many orientations through Meek propagation, so the two
> readings can differ by more than an order of magnitude; we characterise when
> they coincide and report both throughout.
>
> Our claims are about the behaviour of adjustment sets under stressed
> background knowledge, not about recovering the truth. The generating graph is
> used to draw samples and, once every reported quantity is fixed, to check
> whether the certificate held. It never enters the computation.

**Notes for whoever drafts this.**
- The last paragraph is not optional politeness. It is the reviewer defence
  (§8) and must survive to the final version.
- Do not oversell the two-radii point in the abstract as a "correction." It is a
  *measurement result*: the two units differ, here is when and by how much.

---

## 4. The two radii — how to present them

This is the paper's conceptual centrepiece and must be presented as a **finding**,
not as an erratum.

| | `r_hop` (geometric) | `r_claim` (practitioner) |
|---|---|---|
| One step is | one covering step in MPDAG space | one revision of one asserted sentence |
| Defined on | `G₀` alone | `(K, G₀)` |
| Mathematical status | a **metric** — symmetry and triangle inequality verified over all triples (`demo/space.py`, `check_metric_axioms`) | not a metric; depends on how `K` was phrased |
| Answers | how much structure hangs on this knowledge state | how many of my statements can be wrong |
| Size on the committed corpus | 1 to 14 | 1 on all 540 rows (see caveat) |

**How to frame it.** The geometric radius is the right object *mathematically*:
the space of knowledge states is a lattice under model inclusion, every vertex is
a reachable knowledge state (verified equivalent to `{Meek(Ĉ,K)}` on 1,588 states,
commit `748f2e1`), and BFS distance on its covering graph is a genuine metric.
This is the lineage to anchor in prior work on distances between causal graphs
and MPDAGs [TODO: Guo et al., `[RE-18]`].

The claim radius is the right object *operationally*, because an analyst's error
budget is denominated in what they said, not in what they were thereby committed
to. The exchange rate between the two is `k_g0 / |K|` — the **leverage** of a
claim, how many data-compelled orientations hang on it.

**The honest statement, which must appear:** the two coincide exactly when the
analyst's knowledge is already Meek-closed (`|K| == k_g0`). Verified with no
exceptions on the committed corpus: all 114 closed-knowledge rows have
`r_hop = 1`, and all 163 rows where the geometric reading over-states the claim
budget have `|K| < k_g0`. Report `|K| == k_g0` as a per-row status column.

**Do not** present the committed `axisa3` numbers as evidence that `r_claim` is
always 1. That corpus uses a **greedy-minimal** generator set by construction
(`demo/example.py`, `knowledge_to_recover`), which maximises the gap, and its
admission gate selects for `r_claim = 1` (§7, "set aside"). Elicited knowledge
may behave entirely differently — that is `[RE-1]`, and it is the first thing to
run.

---

## 5. The results to foreground

Five. In this order. Everything else is appendix or cut.

### Result A — the instrument, and that it is well-founded
The breakdown radius, the space it lives on, and the fact that it is a metric.
Include the lattice structure (join-semilattice, Theorem 2, proved), and state
the one open property honestly: **Anti-Exchange Case B is verified, not proved**,
and every retraction-only radius is exact only if Conjecture 2 holds, with a
**one-sided** error — radii can only be too large. `THEOREMS.md` §4c, §6, §8.

*Evidence:* `THEOREMS.md`; `demo/space.py`; 1.985M exhaustive radius comparisons
with no Conjecture-2 counterexample (commit `4cf4020`).

### Result B — the worst case predicts the average case
The radius is a **worst-case** distance to the nearest failure. Average-case
survival under randomly corrupted knowledge is a different object. That the first
ranks the second is the finding, and it is not a tautology.

*Evidence:* `table_tau_comparisons.md` (the **session-8 anchor table** — use this,
not the earlier session-7 table). τ_b(`r_val`, survival AUC) positive in 9/9
strata; CI excludes zero in 8/9.

**Report this honestly.** `r_val` is the *sole* predictor whose CI excludes zero
in only one of nine strata; elsewhere it shares significance with the SHD and
`|K|` baselines or is beaten by them. Two cells carry an `⚠ unstable` marker at
N = 200 and **must not be quoted alone** (`[RE-12]`). One stratum
(`flip, cov=0.5, bw=0.25`) is a confirmed null at N = 1000. Say all of this. A
reviewer who finds it themselves will discard the paper; a reviewer who is told
will read it as calibration.

### Result C — knowledge buys identifiability, then only risk
Under truthful knowledge, **statistical efficiency is an invariant of `(Ĉ, X, Y)`** —
not a function of what the analyst asserts. Every identifying truthful proposal
yields the same `Z`, hence the same asymptotic variance. Falsified twice under
pre-registration (0 of 24 strata, then 0 of 32 varied in utility, while `r_val`
varied in 18 and 23 respectively), and Phase 2 was abandoned per the
pre-registered trigger rather than redefined a third time.

So there is **no precision/robustness Pareto frontier**. The operative trade-off
is **identifiability against robustness, and it is a step, not a frontier**:
assert enough to identify (26.9% of proposals could justify no valid set at all),
and every further claim is pure radius cost at zero precision gain.

**This is the result that makes the diagnostic framing necessary**, and it should
be cited in the introduction as the reason the radius is not an optimisation
target. Converges with the independent finding that more claims measurably
*reduce* survival (`|K|` negatively associated, τ ≈ −0.38 to −0.42 in every
tiered stratum).

*Evidence:* `report_fragility_and_pareto.md` §3; `results/axis_robustness/pareto*`.
*Not claimed:* that `r_val` falls monotonically with proposal size beyond the
minimum identifying set — `[RE-13]`.

### Result D — most orientation errors announce themselves; the dangerous ones do not
A two-stage filter. Reversing a single truthful claim renders the knowledge set
**inconsistent** with the CPDAG about 6 times in 10 — Meek closure fails and the
analyst finds out. Among corrupted sets that stay consistent, median survival is
0.952: most survivors are harmless.

The dangerous case is **consistent *and* invalidating** — and that is precisely
the case the radius is defined over, since its space contains only consistent
MPDAGs.

This also dissolves an apparent tension worth pre-empting: `r_val = 1` is common
while average survival at depth 1 is high. Both are true; one is a worst case
over a shell, the other an average over it.

*Evidence:* `report_fragility_and_pareto.md` §2.7. Note the corrected table
(2026-09-10) — use the flip-arm figures, not the retracted pooled ones.

### Result E — correlated errors are worse than uniform ones, at matched intensity
Session 7 claimed this, then **retracted** it (different units, different
populations). Session 8 rebuilt it properly: every instance constructed once,
corrupted both ways, scored on a shared intensity axis, paired within instance.

Result: tiered (correlated) corruption degrades validity faster than uniform in
the intensity range where both arms have real data (≈0.15–0.55), significant in
6 of 7 bins.

**Tell the retraction story.** It is short, it is to the project's credit, and it
inoculates against the obvious methodological objection.

*Evidence:* `results/axis_robustness/XARM_RUN_NOTES.md`, `xarm_paired_test.csv`.

---

## 6. Suggested structure

| § | Content | Status |
|---|---|---|
| 1 | Introduction — §3 above | draftable now |
| 2 | Setting and related work — post-CPDAG BK; adjustment criteria; graph metrics (Guo et al.) | draftable; `[RE-18]` for the cite |
| 3 | The breakdown radius — space, metric, the two units | draftable now (Result A, §4) |
| 4 | Theory — lattice, upward-closure, what is proved vs verified | draftable now |
| 5 | Computation — the up-set reduction, the SAT ladder, one-sided error | draftable now |
| 6 | Evaluation protocol — frozen frame, two panels, audit-after-freeze | **TODO `[RE-7]`, `[RE-15]`** |
| 7 | Results: the instrument ranks survival | draftable now (Result B) |
| 8 | Results: knowledge buys identifiability, then risk | draftable now (Result C) |
| 9 | Results: error visibility and correlation structure | draftable now (Results D, E) |
| 10 | Results: the deployment arm | **TODO `[RE-8]`, `[RE-9]`, `[RE-10]`** |
| 11 | Limitations | draftable now — §8 |
| 12 | Conclusion | last |

**Figure and table budget for nine pages.** Everything else to the appendix.
- **Figure 1.** One real instance. The first-failing retraction named in domain
  terms, and the hop count that missed it. This figure carries the two-radii
  point better than any table.
- **Table 1.** The applicability funnel and the ladder.
- **Table 2.** The certificate ledger — deployment sources above the line,
  controls below. **TODO `[RE-8]`.**
- **Table 3.** The decision table (τ anchor table, condensed).
- **Table 4.** Worst-case radius against average-case survival, same instances.
  **TODO `[RE-11]`** for the real-network version.

---

## 7. Set these aside

Expensive, real, and **not part of this paper's argument**. Cutting them is the
main thing that turns a recital into a paper.

- **The free rule / separation heuristic.** A session-local calibration exercise
  (predict the separation where measured, 1 where not; 636/831 exact, 16
  over-certifications). It compares the instrument to a cheap proxy. It does not
  advance the thesis and invites a baseline argument the paper does not need.
  **Omit entirely** — do not even put it in the appendix unless a reviewer asks.
- **The `r_val = min(s, |K_{G₀}|)` law.** Session 5 established it on a designed
  family (1,572/1,572); session 6 found it is an approximation off that family
  that errs low, with 16 counterexamples clustered on six treatments. Interesting
  internally; a distraction in a nine-page paper. One sentence in the appendix at
  most.
- **The saturation-at-1 story.** The arc "saturation was an artefact of
  Erdős–Rényi, real graphs escape it" is session history, and it is partly
  reversed by the two-radii result. Do not narrate the arc; state the current
  position.
- **Per-session narration.** No "session 5 found…", no "the brief asked for…".
  The pre-registration record is a strength if a reviewer asks for it, not a
  structure for the paper.
- **The stage-0 gate tautology, as a headline.** It belongs in the *methods*
  justification for why the deployment arm exists (see §9), not as a result.
- **Anything marked retracted.** Session 7's cross-arm detectability claim
  (0.354 vs 0.653) is superseded by session 8. The pooled contradiction-rate table
  is corrected. Use the current numbers only.

---

## 8. Reviewer defences to write *into* the paper

Do not wait to be asked.

**"You never validate against the truth."**
Correct, and deliberate — see §2. The audit exists and is opened after freezing.
Report commission and omission separately, and print the rule-of-three floor
beside the dangerous-error rate.

**"`r_val = 1` most of the time — the metric is degenerate."**
Two answers. (1) Worst case ≠ average case; the radius still ranks survival
(Result B). (2) Part of the flatness is a **selection effect of the measurement
gate**, not a property of the metric, and the fix is `[RE-3]`. Be candid about
both.

**"Your certificate over-promises."**
Only under the claim-unit reading, only when `|K| < k_g0`, and we report both
units with the coincidence condition stated. This is Result §4, presented up
front. A reviewer cannot land a blow the paper has already thrown.

**"Oracle CPDAG, so it says nothing about practice."**
Two panels, labelled. **Only the estimated-CPDAG panel carries a deployment
claim.** Do not blur them. `[RE-7]`.

**"An LLM as a knowledge source is a gimmick / it is just reciting the benchmark."**
Answered by design, not by argument: scrambled variable names make recitation
show up as invariance; two model families and two presentation orders separate
instrument variance from source variance; a matched random control distinguishes
a wide radius earned by correct knowledge from one earned by volume and
placement. `[RE-8]`.

**"Synthetic only."**
True of the survival phase today. Say so, and cite `[RE-11]` as the completion.
Real networks currently appear only in the aggressive-pair census (28 networks,
396 pairs).

**Standing scope limits, to state once and plainly:** oracle conditional
independence throughout, causal sufficiency assumed, no finite-sample discovery
outside the estimated-CPDAG panel, one benchmark corpus (`[RE-Tier 4]`).

---

## 9. Where the stage-0 audit belongs

It is **methods justification**, not a result. Two sentences, in §6 or §10:

> Measured on the committed oracle corpus, the admission gate's own criterion
> coincides with the single-claim-retraction test, so the claim radius is 1 on
> all 540 rows by construction. A claim radius above 1 can only be measured on
> knowledge that did not come from the truth, which is why the elicited arm
> carries the deployment claims and the oracle rows serve as controls.

That is all. It explains the design; it is not a finding about the world.

---

## 10. Provenance map

Every number in the paper must resolve here.

| Claim | File |
|---|---|
| Radius convention (`r` vs `r − 1`) | `src/bkrobust/core/conventions.py` |
| Metric axioms, space, covering relation | `src/bkrobust/demo/space.py` |
| Theory status, proved vs verified | `THEOREMS.md` |
| Real-network corpus, 831 admissible rows, funnel | `results/axisa3/instances.jsonl` |
| Corpus manifest, checksums, limit semantics | `results/axisa3/manifest.json` |
| Claim radius, `φ₁`, nesting, amenability | `results/stage0/claim_radius_*.{csv,json}` |
| τ anchor table (**use this one**) | `table_tau_comparisons.md`; `results/axis_robustness/table_tau_source.csv` |
| Survival curves, AUC endpoints | `results/axis_robustness/survival_*.csv` |
| Pareto falsification, both rounds | `results/axis_robustness/pareto*` |
| Matched cross-arm | `results/axis_robustness/xarm_*`, `XARM_RUN_NOTES.md` |
| Null-cell N = 1000 re-run | `results/axis_robustness/null_*` |
| Pre-registration and every retraction, in order | `results/axis_robustness/PREREGISTRATION.md` |
| Evaluation design | `docs/EVALUATION_PLAN.md` |
| What is still unrun | `docs/REMAINING_EXPERIMENTS.md` |

**Session reports** (`report_*.md`, `REPORT_*.pdf` at the root) are the long-form
record. They are sources, not drafts: their framing is session-local and several
of their headlines are superseded. Mine them for numbers, not for narrative.

---

## 11. Conventions to hold throughout

- **Naming.** `r_hop` for the geometric radius, `r_claim` for the practitioner
  radius. The legacy name `r_val` appears throughout the code and older reports
  and means `r_hop`; do not use `r_val` in the paper. `r_opt` and `r_eps` keep
  their meanings (first failure of *optimality*, first exceedance of a bias
  tolerance ε).
- **"Two radii" is overloaded — never use the phrase unqualified.** There are two
  independent pairs. One is indexed by *which property breaks*
  (`r_val` / `r_opt` / `r_eps`); the other by *what one move costs*
  (`r_hop` / `r_claim`). The README and the older reports use the phrase for the
  first; this document uses it for the second. In the paper, say
  "validity and optimality radii" or "geometric and claim units" explicitly.
- **`k_g0` is not "the knowledge size."** It is the number of orientations in
  `G₀` — the closure. Older reports call it knowledge size; that is the
  conflation the two-radii result is about. Call it the **commitment size**, and
  call `|K|` the **stated size**.
- **Never aggregate `UNREACHED`.** It means no element of the space fails. It is
  not a large radius, must never be averaged or plotted on a numeric axis.
- **Censoring is a status, never a drop.** Wall-capped instances are reported as
  censored.
- **Denominators.** Every knowledge-dependent condition is a status column on the
  frozen frame, never a filter, so every knowledge source reports the same
  denominator.
- **Networks, not pairs, are the unit of analysis.** The intraclass correlation of
  the `r = 1` indicator is 0.40 over 25 networks (design effect ≈ 12.7, effective
  n ≈ 65 against 831 rows). Bootstrap over networks; print no per-row interval.
- **Assumption carried on every radius:** Conjecture 2, hence Anti-Exchange
  Case B — verified, not proved. One-sided error. State it once prominently, then
  cite the theorem file.

---

## 12. First actions for the next session

1. Read `core/conventions.py`, `THEOREMS.md` §8 (summary table), and
   `demo/space.py`'s module docstring. Twenty minutes, and it prevents the two
   most likely errors.
2. Run the cheap proxy for `[RE-1]` — simulate total-order responses on the
   `axisa3` chain components and compute `|K|` vs `k_g0`. It costs little and it
   decides how §4 is written.
3. Draft §§1–5 and 7–9. They are fully supported by committed data today.
4. Leave §6 and §10 as structured `TODO` stubs naming the `[RE-n]` items they
   wait on. **Do not write results that do not exist yet**, not even as
   placeholders with plausible numbers.
