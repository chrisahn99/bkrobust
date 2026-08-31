# Breakdown radius for background knowledge: what we have measured

**To:** Seong Woo (Chris) Ahn
**From:** José Lucas De Melo Costa
**Date:** 2026-08-31
**On:** *A Breakdown-Radius Theory for Background Knowledge in Causal Effect Estimation* (your rewrite)
**With:** the code, in this directory, and the raw run registers it was computed from

---

## 0. Where we agree, which is most of it

Your rewrite settles four things that were open since the 14 August meeting, and it settles them the
right way.

The object is now a breakdown radius with `r_opt` and `r_val` kept apart, and section 5 states that
the radius is not a probability. That closes the "bound" line from the meeting notes. A continuity
bound was never available to us: conditional on `O*` moving, our damage distribution is heavy-tailed
(median 0.59, up to 128× on scale-free graphs), so any sentence of the form "small knowledge error
implies small bias" would have been falsified by our own data.

Taeb, Guo and Henckel is now the methodological pivot, and section 10 is right that the version has
to be pinned. It is **v3, dated 2026-06-30**. Section 6 demoting metricity to a corollary is also
right, and so is the decision in section 10 not to lead with their open problem.

You also bring a literature we did not have. Masten and Poirier on breakdown frontiers, and Horowitz
and Manski, are the correct econometric anchor for what we are doing, and neither appeared anywhere
in our corpus before your document.

The rest of this report is the three places where our measurements change your plan, plus two work
packages you have scheduled that are already done.

---

## 1. The unit changes every number we hold, by a factor of two

Section 4.1 prices a flip at two: a retraction followed by an assertion. Our pilot's perturbation
operator **is** a flip. It reverses the orientation of an edge that exists in the CPDAG skeleton
(`pilot/code/run_linear.py:132`), and the budget is global at `MAX_K = 4`
(`pilot/code/run_linear.py:34`).

In your units the pilot has been running at radius 8, not radius 4.

Nothing should cross between the two documents until this is written down, because it cuts both
ways. It says our feasibility evidence is stronger than section 8.2 assumes, since we reach radius 8
on `p <= 10` in about twenty minutes of CPU. It also says every censoring figure below sits on a
different abscissa from the one your figures will use.

| quantity | our unit | your unit |
|---|---|---|
| one expert error | 1 (a flip) | 2 (retract, then assert) |
| pilot budget `MAX_K = 4` | 4 | 8 |
| `rho* = 1` (breaks at the first error) | 1 | 2 |

---

## 2. The locality law is the result, not a pruning heuristic

Section 8.2 lists "prune by relevance" among the computational mitigations: orientations that cannot
affect a relevant path near `X`, `Y` and the covariates can be quotiented out. That quotient is not a
speed-up. It is the finding, and it is the one asset we have that is measured, sample-size invariant
and directly actionable for a practitioner.

Censoring of `rho*_any` at `n = 20 000`, stratified by the minimum hop distance from an elicited
statement to `{X, Y}` on the CPDAG skeleton (`locality/x1/`, anchor defined at
`locality/x1/code/x1_ops.py:271`, minimum taken over both endpoints):

| ensemble | dist 0 | dist 1 | dist >= 2 |
|---|---|---|---|
| `original` (p 5-8) | **0.037** (n = 2 539) | **1.000** (n = 93) | |
| `licensed` (p 5-10) | **0.070** (n = 1 852) | **0.935** (n = 229) | |
| `k8`, \|K\| = 8 | **0.016** (n = 563) | 0.955 (n = 22) | |
| `large` (p 15-25) | **0.330** (n = 1 472) | 0.977 (n = 638) | **1.000** (n = 331) |

The separation is 13 to 27 fold. It is invariant in `n` to four decimals (0.9345 at n = 2 000,
20 000 and infinity alike), and it is exactly 1.000 for all 331 `large` SCMs whose statements sit two
or more hops from the query.

The reading is mechanical. `rho*` is a functional of `O*`, `O*` is local, so a statement that does
not reach the query cannot move it, and no sample size repairs that because zero times anything is
zero.

This also dissolves an anomaly we had both been treating as a separate regime. The `large` ensemble
behaves differently because only **60.3 %** of its SCMs have a statement at distance 0, against
88.7 % for `licensed` and 96.3 % for `original`. Condition on distance 0 and it stops being a
different regime.

**Suggestion.** This is what section 3 of the paper should carry, and it gives the abstract a second
currency that costs a reader nothing to understand: elicited knowledge far from the query does not
move the estimate, so stop paying to elicit it.

---

## 3. The principal risk in section 10 is real at both ends, and the switching variable is missing

Section 10 fears that `r_val(O) = 1` almost always, and makes the choice among three papers depend on
WP0. Both extremes are true at once in our runs, and they are separated by hop distance rather than
by budget.

The radius saturates at the **top** for the population. Censoring of `rho*_any` moves from **0.554**
at n = 200 to **0.515** at n = infinity, that is 3.9 points across five orders of magnitude of sample
size, against a hard floor. The cause is that **64.6 %** of licensed problems have `d(1) = 0`
exactly, at every sample size (`e1prime/RESULTS.md`, section on the width ratio).

The radius saturates at the **bottom** as soon as we repair the other end. The composite
`min(rho*_se, rho*_ident)` cuts censoring to **0.166** and puts **69.0 %** of the mass at
`rho* = 1`, which is literally your principal risk.

So WP0 does not have one answer. It has one answer per stratum, and the stratifying variable is not
in the document. Your branch "the O-set always breaks at 1" is true near the query; your branch "the
O-set survives" is true away from it; an unstratified distribution of `r_val` averages two regimes
and describes neither.

One consequence for section 10's venue paragraph. The empirical ICLR hook you propose,
*knowledge-informed methods break at radius 1-2 on standard benchmarks*, is contradicted by our data
as written, because more than half of problems never break inside the budget. The conditional version
survives and is stronger, because it tells the practitioner what to do.

Your own assessment of 22 August arrived at the same place from the other side: if influence is
local, one false statement touches `O(1)` of `O(d^2)` queries, so the mean damage tends to zero as
the graph grows. That is an argument for worst-case or distributional aggregation and against
averages, and I think you are right.

---

## 4. Section 2 puts b-LOAD outside the framework, and the exclusion is not clean

Section 2 argues that adjacency knowledge is vacuous or refuted, so the perturbation space is forced
to be orientations. The argument is good. Its consequence is uncomfortable.

b-LOAD's noise operator asserts a required edge on pairs that are mostly **non-adjacent**
(`src/background_knowledge.py:166-169` in this repository). Under section 2 those assertions are
refuted, so they fall outside `G_Chat` and outside the theory. We measured that disjointness on
14 August; section 2 does not repair it, it formalises it. That matters because the empirical demand
that motivates this whole line sits on that side: all four CIKM reviewers and the metareview raised
the perfect-knowledge assumption, and R4 asked specifically for a robustness analysis on Intervention
Distance in the main body.

The exclusion is also not clean. On the `licensed` ensemble, b-LOAD's pool is only **69.4 %**
`pure_spurious` (5 555 of 8 000). **30.6 %** is `reversal_of_true_edge`, which is inside the space,
and **12.9 %** is the query pair itself. The split moves with the ensemble: 73.1 / 27.0 on
`original`, 91.4 / 8.6 on `large`, so the number to quote is the one for the graphs you benchmark on.
"Disjoint by construction" holds for the operators and fails for the sampler. One caveat that our own
audit attaches to this partition and that I am passing on rather than dropping: it uses **oracle**
labels, so no rate stratified by it can be turned into a triage rule for an analyst.

There is a good way out, and it is yours. Section 2 offers a free consistency check: any adjacency
assertion conflicting with the skeleton can be flagged before any radius is computed. Applied to
b-LOAD that check catches **all** of its spurious noise, which is a publishable result and a direct
answer to R4. It has to be written honestly, though, because it reverses the story: b-LOAD's
robustness under that noise stops being a surprise and becomes a consistency test that was never
run. Note that the current code cannot perform it, since the assertion is written into the
initialisation matrix and then protected by `update_graph_mpdag`.

Whichever of the two papers absorbs that sentence should decide it before either is submitted.

---

## 5. Amenability needs to be a fourth exit beside bottom

`cl(G, a->b)` returns bottom for a directed cycle or a forbidden structure. It says nothing about the
case where `O(G)` is undefined because the MPDAG is not amenable, and in an amenable MPDAG the edges
into `cn(X, Y)` are already forced.

This is not a corner case in our runs. On 2 058 single-error perturbations the pilot classifies
**22.9 %** of Meek-consistent outcomes as loud non-amenable failures (15.3 % of all single errors),
which is a different event from silent bias because the analyst sees it. Non-amenability is also a
filter on which problems enter the analysis at all, and we have never written it down as a scope
condition. We should (`pilot/RESULTS.md`, the outcome table under NUMBER 1).

Pricing a loud failure as a breakdown inflates fragility, which is the objection section 4.2 already
makes against collapsing `r_opt` and `r_val` into one number.

---

## 6. Two work packages are already done

**WP0(d), the decisive statistic.** Section 3 says the argument against naive assertion counting
turns on the fraction of single-assertion flips that are logically inconsistent with the skeleton and
therefore generate no valid MPDAG, and that a two-paragraph invariance argument loses to a one-line
alternative unless that fraction is substantial. It is substantial, and we have it three
times from two independent runs: **33.0 %** of all single-orientation errors are Meek-inconsistent in
the pilot, **32.4 %** on the G1 denominator in E1-prime, and **23.6 %** of single errors are caught
free by Meek in ARM 2 (38.6 % for double errors). Your section 3 wins its argument with a number that
already exists, and the pilot and E1-prime agree to within half a point on it.

**The section 4.3 calibration experiment.** Corrupt the knowledge, compute `r` around the corrupted
`G_0`, and check whether `d(G_0, G_true) < r` predicted validity. That is ARM 2 of E1-prime, run on
19 August over 2 500 SCMs times 20 datasets. Coverage at r = 0 is **0.9509** against 0.95 nominal,
with zero failures of the regeneration gate.

It also produced the number section 5 asks for, the one that is directly comparable to a standard
error (n = 20 000):

| | r = 0 | r = 1 | r = 2 |
|---|---|---|---|
| naive CI coverage | **0.951** | **0.720** [0.713, 0.727] | **0.568** [0.559, 0.577] |

**One reversed background-knowledge statement costs 23 points of coverage, and two cost 38**, among
the errors that Meek does not catch for free.

And the structural result, which is the best answer to "collect more data": `SE_n` shrinks like
`n^{-1/2}` while `d(rho)` is free of `n`, so the ignorance does not shrink with data. The mean width
ratio moves from **1.61** at n = 200 to **6.84** at n = 20 000.

---

## 7. Two numbers of ours that you should stop using

I would rather you heard these from me.

`0/791` is dead. Our own audit on 19 August marked it as having no provenance, and it is the figure
your 22 August assessment quoted back at us. The locality claim is carried by the censoring table in
section 2 above, which is a frequency, not a zero-existence result.

The **median** width ratio of 1.00 flatters the method. The ball is inert for 69 % of problems, so
the median is 1.00 by arithmetic. I pre-registered the median; the mean is 6.84×, which would have
failed the pre-registered 3× bar. Both belong in the paper.

---

## 8. What I need from you, and what I owe you

From you, in order of value before the 18 September abstract:

1. **Read Taeb, Guo and Henckel v3 in full and tell me what section 6 actually does.** Our corpus
   note is explicit that everything it says about section 6, Lemma S6 and the open covering relation
   is second-hand from your framing. If their adjustment-set robustness application already defines a
   tolerance radius centred on the analyst's graph, the novelty surface shrinks a lot, and I would
   rather know now than in a review.
2. **The covering enumeration on four and five vertices** (your WP2). It is also the independent
   correctness check on a preprint that no referee has read, which section 10 correctly treats as a
   small credibility gain rather than a liability.
3. **Guo and Perković 2010.08611 in full.** The framing rests on "they never perturb K" and nobody on
   either side has actually read it.

From me: this directory. It is the code you asked for on 21 August, and I am sorry it took ten days.

---

## What is in this directory

| path | what it is |
|---|---|
| `pilot/` | the 17 July run: 600 primary plus 1 650 sweep plus 60 nonlinear SCMs, iSCM parameterisation frozen to closed-form linear SEMs, Meek closure brute-force validated on 7 688 knowledge sets |
| `e1prime/` | the 19 August run: the `se` criterion, the censoring floor, and ARM 2 (calibration and coverage). `PREREG.md` and `AUDIT.md` were written before the first number |
| `locality/x1/` | the hop-distance stratification, the operator arms, and the pool composition |
| `locality/x2/` | the exhaustive enumeration: events as functions of `(G_0, x, y, s)`, no SCM and no truth, complete at `p <= 5` |
| `NUMBERS.md` | every number in this report, with the file it comes from |

Start with `NUMBERS.md` if you want to check a figure, and with `e1prime/PREREG.md` if you want to
see what was declared before the run.
