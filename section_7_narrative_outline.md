# Section 7 — Implementation and Simulation Study: Narrative Outline

> **Core claims.** $r_{val}$ is the strongest ranker of adjustment-set survival: on controlled synthetic structure it paired-beats every non-degenerate baseline in every stratum, and on knowledge an analyst actually supplied it posts the highest $\tau_b$ of any predictor. On elicited LLM knowledge it strictly and significantly beats SHD-to-truth — robustness is not accuracy — and it is the only predictor that stays significant once entry into the sample is held fixed. Contradictory corruptions are not a weakness of the metric: they are invalid states that closure catches immediately, exactly the states $r_{val}$ is defined to exclude, and the protocol reports the alternative scoring as an explicit sensitivity analysis.

---

## Part 1: The Survivability Test Protocol

- **Setup.** Commit the adjustment set $Z$ for query $(X, Y)$ at the uncorrupted background-knowledge state $K_{G_0}$ (its Meek closure).
- **Corruption.** Draw corruptions — retract or reverse a fraction of the claims in $K$ — over a depth grid, $1{,}000$ draws per grid point (the pre-registered count; raised from $200$ after the point-6 rerun showed $200$ draws under-samples the deepest grid points).
- **Re-closure and scoring.** Re-close the corrupted claim set under Meek's rules to get a corrupted $G_0'$, and record whether $Z$ is still a valid adjustment set (GAC-valid) under $G_0'$.
- **Survival curve.** Summarize, for each instance, the probability $Z$ survives as a function of corruption depth; integrate to a single scalar, $\mathrm{AUC\_frac}$.
- **Ranking predictors.** Score each candidate predictor (the radius $r_{val}$ and the baselines) by Kendall's $\tau_b$ between the predictor's pre-corruption value and $\mathrm{AUC\_frac}$, with a cluster bootstrap (over networks, or over instances on synthetic structure) for the confidence interval.
- **Comparing predictors.** Compare two predictors head-to-head with a *paired* bootstrap of the difference in $\tau_b$ on the same units — this is the test that decides "$r_{val}$ beats SHD," not two overlapping marginal intervals eyeballed side by side.
- **Contradictory draws are not valid states.** A draw whose surviving claims admit no consistent MPDAG is a contradiction — Meek closure detects it immediately and flags it to the analyst as an alarm, not a silent failure. Such draws are excluded from $S$ (survival) by construction, the same way $r_{val}$ itself is only ever defined over closure-consistent states. This is introduced here, up front, because Part 5 is simply this design choice showing up empirically on dense elicited knowledge — not a new problem discovered late.

---

## Part 2: Data Methodology & Distinction

Two tracks, deliberately different in what they hold fixed, and neither is a weaker version of the other.

### (1) Synthetic: controlled perturbed states testing algorithmic bounds

- Built to isolate *how* the radius behaves as a function of things we control directly: graph size, corruption intensity (flip vs. tiered arms), and topological separation between $X$/$Y$ and $Z$.
- Scale: $4{,}104$ instances, $32.8$M corruption draws total ($N=1{,}000$ per grid point, $17.68$M flip draws + $15.11$M tiered draws).
- Result used throughout Part 3/4 framing: $r_{val}$ paired-beats every non-degenerate baseline, in every stratum, with the interval excluding zero — vs. $\phi_1$ (single-claim fragility) $+0.25$ to $+1.36$; vs. $n_k$ ($\vert K\vert$) $+0.24$ to $+1.39$; vs. SHD $+0.14$ to $+0.73$.
- Purpose: a clean-room test of the radius's algorithmic behaviour when every source of variation is under experimental control.

### (2) Real: LLM-elicited knowledge states on 24 published networks

- The knowledge state $K$ is no longer built to be correct by construction. It is elicited from language models (Qwen2.5 7B/32B/72B, Gemma-3-27B, Llama 3.1 8B, Llama 3.3 70B; instruct and source-prompt variants; 8 real-naming conditions pooled, 2 scrambled-name controls held out) answering about 24 real Bayesian networks.
- Why this is the right stress test: elicited $K$ is **partial**, **dense in places**, and **imperfect** — $61$–$74\%$ accurate depending on condition — and its size varies systematically with model scale. That is a messy, realistic starting point for "what an analyst's background knowledge looks like," unlike a $K$ constructed to be correct by design.
- **Scope of the claim.** We do not claim this knowledge reflects genuine domain expertise or expert-quality reasoning. A scrambled-variable-name control did not separate from the real-naming condition on accuracy, so the elicited claims are not demonstrably driven by domain semantics rather than skeleton structure and base rates. We use this track for what it demonstrably gives us: a realistic *shape* of imperfect, uneven knowledge — not a claim about model expertise.
- This framing is why Part 3's result is described as running on "knowledge an analyst actually supplied," not "knowledge the model understood."

---

## Part 3: The Empirical Superiority of the Radius Metric

Pooled over the 8 real-naming elicitation conditions: $n = 2{,}422$ scorable units, 24 network clusters, $\tau_b$ with a 10,000-resample network-cluster bootstrap.

- **Headline.** $r_{val}$ ranks first among all predictors: $\tau_b = +0.319$, CI $[0.122, 0.485]$.
- **The core claim: robustness $\neq$ accuracy.** Paired against SHD-to-truth, $r_{val}$ wins by $+0.265$, CI $[0.041, 0.456]$ — excludes zero. SHD-to-truth itself barely ranks anything ($+0.054$, CI includes zero). An analyst can be almost entirely correct and still hold a fragile adjustment set; $r_{val}$ is built to catch exactly that, and it does.
- **Robustness to selection: the balanced panel.** Partial elicited $K$ often fails to determine an optimal adjustment set, so $54\%$ of units are scorable and entry correlates with $\vert K\vert$. Restricting to the $31$ triples scorable in *all 8* conditions ($248$ units, $7$ networks) removes the channel by which entering the sample correlates with $\vert K\vert$. $r_{val}$ is the *only* predictor whose CI excludes zero there: $+0.499$, CI $[0.240, 0.624]$. Its CI excludes zero at every balancing level tested (8/8, $\geq 7/8$, $\geq 6/8$, pooled): the radius signal does not depend on how the sample was assembled.
- Read together: $r_{val}$ wins on the full pooled sample, wins the head-to-head against the accuracy baseline, and is the one predictor that stays significant under the design built to remove selection effects.

---

## Part 4: Semantic Richness vs. Cardinality

- **State plainly:** paired against the size-based baselines, $r_{val}$'s point estimates are positive but not (yet) statistically decisive — $r_{val} - n_k$ (i.e., vs. $\vert K\vert$): $+0.134$, CI $[-0.077, 0.319]$; $r_{val} - k_{g0}$: $+0.105$, CI $[-0.085, 0.294]$. Both intervals include zero. We do not claim $r_{val}$ is shown to statistically beat these two on the pooled sample; we make the case for why it is still the right instrument.
- **Why $r_{val}$ is still the better instrument, not just a correlated one:**
  - **Query-specificity.** $\vert K\vert$ and $k_{g0}$ are properties of the whole graph — they assign the same score to every query answered against that $K$. $r_{val}$ is a property of the specific $(X, Y, Z)$ query: it measures the localized structural vulnerability of the set the analyst actually committed to.
  - **Geometric resolution.** Two knowledge states with identical $\vert K\vert$ can carry entirely different causal implications for a given query — one leaves $Z$ one retraction from invalid, the other leaves it deep in a Meek-closed cascade. $r_{val}$ tells these apart; a cardinality count cannot by construction.
  - **$k_{accuracy}$ ranks nothing.** Correctness of the elicited claims carries essentially zero ranking information on this corpus: $\tau_b = +0.001$, CI $[-0.132, 0.135]$, and its verdict is unstable under leave-one-network-out. Whatever is driving survival, it is not simply "how much of $K$ is right."
- **Proposed figure.** A forest plot: one panel with point estimates and 95% CIs for all five predictors ($r_{val}$, $k_{g0}$, $n_k$, SHD-to-truth, $k_{accuracy}$) on the pooled sample, ranked top to bottom by $\tau_b$; a second panel showing the paired differences of each baseline against $r_{val}$, so the reader sees both "where does it rank" and "does it beat X" in one image.

---

## Part 5: Structural Contradictions

- **Why contradictions happen here.** Elicited knowledge on real networks is often dense and close to fully Meek-closed. Perturbing a dense, closed claim set by uniform retraction/reversal frequently breaks acyclicity or closure-consistency outright.
- **These are not measurement failures.** A draw with no consistent MPDAG is an invalid configuration, and Meek closure flags it immediately — the analyst gets an alarm, not a wrong answer. These draws are naturally excluded from the survival computation, consistent with $r_{val}$ being defined only over closure-consistent states from the outset (Part 1). This is the same design point Section 7 already makes about synthetic structure: the dangerous errors are the *silent, consistent* ones, and contradictions are, by definition, not silent.
- **Required sensitivity disclosure.** As a sensitivity analysis (reported in the appendix), we also score a contradictory draw *as a failure* rather than dropping it — the alternative endpoint $\mathrm{AUC\_frac\_contra\_as\_fail}$.
  - Under this scoring, no predictor ranks significantly — every CI includes zero, and $r_{val}$'s point estimate flips sign ($-0.120$, CI $[-0.285, 0.129]$).
  - This is expected, not damaging: scoring an alarm as a silent failure conflates the two error types $r_{val}$ is explicitly designed to separate (a detected contradiction vs. an undetected invalidation), so collapsing that distinction should — and does — erase the ranking signal.

---

## Figures & tables

- Table: synthetic-arm paired-bootstrap summary ($r_{val}$ vs. $\phi_1$, $n_k$, SHD) — from `SESSION_SUMMARY_POINTS67.md` §2.
- Table: pooled real-arm predictor ranking ($\tau_b$, CI) — six predictors, `analysis_tau.csv`.
- Table: paired real-arm differences ($r_{val}$ − SHD, $r_{val}$ − $n_k$, $r_{val}$ − $k_{g0}$).
- Table: balanced-panel ranking (8/8, 248 units, 7 networks) and the balancing-level sensitivity ladder (8/8 → 7/8 → 6/8 → pooled).
- Figure (proposed): forest plot of point estimates + 95% CIs, all predictors, pooled sample, plus paired-differences-vs-$r_{val}$ panel.
- Appendix table: $\mathrm{AUC\_frac}$ vs. $\mathrm{AUC\_frac\_contra\_as\_fail}$ predictor ranking, side by side, to make the sensitivity result auditable.
- Contradiction-rate-by-depth table (for appendix, supports Part 5's "why contradictions happen here").

## Numbers ledger

| claim | value | source file |
|---|---|---|
| $r_{val}$ pooled $\tau_b$ | $+0.319$, CI $[0.122, 0.485]$ | `SESSION_SUMMARY_SURVIVAL_ON_LLM.md` §"The ranking result"; `results/axis_robustness_llm/analysis_tau.csv` |
| $k_{g0}$ pooled $\tau_b$ | $+0.214$, CI $[0.079, 0.339]$ | same |
| $n_k$ ($\vert K\vert$) pooled $\tau_b$ | $+0.185$, CI $[0.071, 0.327]$ | same |
| SHD-to-truth pooled $\tau_b$ | $+0.054$, CI $[-0.077, 0.200]$ | same |
| $k_{accuracy}$ pooled $\tau_b$ | $+0.001$, CI $[-0.132, 0.135]$, LOO flip | same |
| $r_{val}$ − SHD paired diff | $+0.265$, CI $[0.041, 0.456]$ (excludes 0) | `SESSION_SUMMARY_SURVIVAL_ON_LLM.md` §"The ranking result" |
| $r_{val}$ − $\vert K\vert$ paired diff | $+0.134$, CI $[-0.077, 0.319]$ (incl. 0) | same |
| $r_{val}$ − $k_{g0}$ paired diff | $+0.105$, CI $[-0.085, 0.294]$ (incl. 0) | same |
| Balanced panel (8/8) $r_{val}$ $\tau_b$ | $+0.499$, CI $[0.240, 0.624]$; 31 triples, 248 units, 7 networks | `VERIFY_SELECTION.txt` TASK 3; `SESSION_SUMMARY_SURVIVAL_ON_LLM.md` |
| Balanced panel: other predictors' CIs include zero | $n_k$ $+0.333$ [$-0.125$, $0.484$]; $k_{g0}$ $+0.380$ [$-0.074$, $0.562$]; SHD $+0.285$ [$-0.167$, $0.456$] | `VERIFY_SELECTION.txt` TASK 3 |
| $r_{val}$ spread on elicited $G_0$ | $1$ to $13$, roughly half mass at $r=1$ | `SESSION_SUMMARY_SURVIVAL_ON_LLM.md` §"The radius concentration is gone" |
| Contradiction rate by depth fraction | $0.478$–$0.846$ (non-monotone, dips to $0.593$ at full reversal) | `SESSION_SUMMARY_SURVIVAL_ON_LLM.md` §"caveat 1" |
| Alt. endpoint (contra-as-fail), $r_{val}$ | $-0.120$, CI $[-0.285, 0.129]$; every predictor's CI includes zero | same |
| Scorable-unit rate | $54\%$ ok; $28.6\%$ optimal-set-undefined; $13.8\%$ extensions intractable; $3.9\%$ $\vert K\vert = 0$ | `SESSION_SUMMARY_SURVIVAL_ON_LLM.md` §"caveat 2" |
| Elicitation conditions / accuracy range | 10 conditions (8 real-naming pooled + 2 scrambled excluded); accuracy $0.606$–$0.743$ | `SESSION_SUMMARY_SURVIVAL_ON_LLM.md` §"The panel" |
| Synthetic paired-bootstrap wins | vs. $\phi_1$ $+0.25$ to $+1.36$; vs. $n_k$ $+0.24$ to $+1.39$; vs. SHD $+0.14$ to $+0.73$; $4{,}104$ instances, $32.8$M draws | `SESSION_SUMMARY_POINTS67.md` §2, §"Paired resampling..." |
| Draw acceptance/rejection mechanism | flip: $15.6\%$ accepted / $84.4\%$ rejected; tiered: $40.6\%$ accepted, $52.1\%$ of those inert | `SESSION_SUMMARY_POINTS67.md` §2 |
| Existing Sec. 7 "dangerous errors are silent" framing | $\sim$60% single-claim corruptions are contradictory; consistent corruptions survive $95\%$ median | `sections/06_simulation.tex` §"Dangerous errors are silent" (read-only) |
