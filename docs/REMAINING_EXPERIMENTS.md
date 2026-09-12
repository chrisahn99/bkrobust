# Remaining experiments

Written 2026-09-12, immediately after `experiments/survival-and-pareto` and
`iclr2027-eval-plan` were merged into `main`. Nothing in this file has been run.

This is the ledger of what must still be measured before the paper's claims are
complete. It is ordered by **what blocks what**, not by effort. Each item states
the question, the artifact it produces, what it would change if it came out the
other way, and where the code already lives.

Companion document: [`docs/PAPER_NARRATIVE.md`](PAPER_NARRATIVE.md) — what the
paper argues and which of these results it needs. Items below are cross-referenced
from there as `[RE-n]`.

Two notes that apply to everything here:

- **Environment.** This machine runs Python 3.9.6. The scaffold stubs need 3.11+
  (`TypeAlias` / `StrEnum`), so 23 tests fail and 3 files fail to collect. The SAT
  legs need `ortools`. Stand up Python 3.12 with `ortools` before running anything
  below; do not interpret the current test failures as regressions.
- **Every committed radius assumes Conjecture 2** (hence Anti-Exchange Case B,
  verified not proved). The error is one-sided: radii can only be too **large**.
  See `THEOREMS.md` §4c and §6. `[RE-9]` is the only item that attacks this.

---

## Tier 0 — blocking. These decide what the paper can claim.

### [RE-1] Is elicited background knowledge Meek-closed?

**Question.** For knowledge as a real analyst actually supplies it, how large is
`|K|` (sentences asserted) against `k_g0` (orientations committed to)?

**Why it blocks everything.** The certificate sentence — *"`r_val − 1` revisions
of your stated knowledge are safe"* — is exactly sound when `K` is Meek-closed
(`|K| == k_g0`) and unsound otherwise. On the committed `axisa3` corpus it is
sound on 114 of 540 rows and violated on 163; the corpus is the worst case by
construction, because `knowledge_to_recover` is greedy-minimal and deliberately
picks the highest-leverage generator at each step (`demo/example.py`).

A real analyst is not a minimal generator set. The questionnaire in the
evaluation plan asks for **a causal order over each chain component**, which
yields every orientation in that component rather than a generator for it. If
that produces `|K| ≈ k_g0`, the defect is confined to the synthetic oracle rows
and `r_hop` keeps its practitioner reading in the deployment arm.

**Artifact.** A two-column table, `|K|` against `k_g0`, per elicited instance,
plus the fraction with `|K| == k_g0`.

**Decides.** Whether the paper leads with one radius plus a caveat, or two radii
as co-equal quantities. Run this *before* rewriting any framing.

**Depends on.** `[RE-5]` (elicitation must exist first). A cheap proxy that can
run today: simulate a total-order response on each chain component of the
`axisa3` networks and compute `|K|` vs `k_g0` without any model in the loop.
**Do the proxy first.**

---

### [RE-2] `r_claim` at partial coverage

**Question.** What is the claim-unit radius at coverage 0.5 and 0.25?

**Status.** `experiments/stage0_claim_radius.py` ran at coverage 1.0 only — 540
rows. There are no `r_claim` numbers anywhere at partial coverage.

**Why it matters.** At coverage 1.0 the claim-radius measurement is degenerate:
`fast_gate`'s admission test *is* the single-claim-retraction test, so 540/540
is the entry requirement read back out, not a measurement (see `[RE-3]`). Partial
coverage is the first setting in the existing corpus where `r_claim` can vary.

**Artifact.** Extend `stage0_claim_radius.py` to loop coverage; same CSV schema.

**Watch out.** The coverage sweep is not nested (`[RE-4]`). Fix that first or the
comparison is confounded.

---

### [RE-3] Record the atomic-perturbation clause as a status, not a filter

**Question.** What do the 1,659 instances rejected as
`no_atomic_perturbation_changes_validity` actually look like?

**Why.** `fast_gate` (`benchmarks/measure.py:62`) admits a pair only if some
single retraction breaks validity. That is precisely the condition
`r_claim == 1`. The 1,659 rows it rejects are therefore *exactly* the
`r_claim ≥ 2` population — the only rows that could exhibit variation in the
quantity the paper wants to report — and they are discarded before measurement.

**This is the cheapest high-value item in the file.** It needs no elicitation, no
new data, and no model: change one clause from a filter to a status column and
re-run the sweep over data already on disk.

**Artifact.** A re-run `results/axisa3/instances.jsonl` in which
`atomic_perturbation_changes_validity` is a boolean column and admissibility no
longer depends on it. Expected yield: 831 → ~2,490 admissible rows, with a
genuine `r_claim ≥ 2` stratum.

**Keep the old file.** Commit the re-run alongside, do not overwrite; prior
published numbers must stay reproducible from the code that produced them (the
precedent is commit `748f2e1`, which left `enumerate_space` unmodified).

---

## Tier 1 — pipeline corrections. Small, each needs a differential test.

These are the five patches listed in `docs/EVALUATION_PLAN.md`. Each should be
validated by re-running against the committed 831 rows and diffing.

### [RE-4] Nested coverage sweep
`select_knowledge` (`benchmarks/measure.py:103`) uses a stride
`⌊i·|K|/keep⌋`, which selects different index sets at different coverages.
`K(0.25) ⊄ K(0.5)` on `diabetes`, `ecoli70`, `magic-irri`, `munin1`. The sweep is
therefore not a monotone ablation, and one pair (`munin1`,
`R_MED_RD_EW → R_MED_CV_EW`) is admissible at 1.0 and 0.25 but not 0.5.

**Fix.** Prefixes of one seeded permutation. **Caution:** `experiments/stage0_claim_radius.py`
re-implements this function locally rather than importing it. Import it, or the
nesting result becomes unreproducible the moment this lands.

### [RE-5a] Three-valued verdict for the optimal set
Distinguish *valid*, *invalid*, and *undefined* rather than collapsing the last
two.

### [RE-5b] `g0_undirected_edges` set in every branch
Currently unset on some rejection paths, so the intractability limit cannot be
separated from structure post hoc.

### [RE-5c] Dispatch leg on every row
Record whether `local_up_fast` or `e1_ladder` produced each radius. Currently
inferable from `method` but not guaranteed present.

### [RE-5d] Split `not_amenable` from `o_g0_not_identified`
All 518 `o_g0_not_identified` rejections were verified non-amenable at `G0`
(stage 0, check 3). They are a structural exclusion, not a definitional artifact,
and should be labelled as such.

### [RE-6] Gate the corpus on matched populations
Per `NEXT.md` item 5: admissible instances collapse 543 → 106 as coverage falls
and survivors have larger separation, so marginal rates across coverage are a
composition effect. Condition on a fixed instance set — the same `(network, X, Y)`
at every coverage — for any coverage comparison.

---

## Tier 2 — the deployment arm. This is the paper's main new experiment.

### [RE-7] The frozen frame
One candidate set per network, **frozen before any knowledge is collected**.
`X` in or adjacent to a chain component of size ≥ 2; `Y` a possible descendant.
Every knowledge-dependent condition is a **status column on that frozen list,
never a filter**, so every knowledge source reports the same denominator. Carry a
hash of the CPDAG source and of the generating DAG on every row.

Two panels, labelled as such:
- **Oracle CPDAG** on the 25 yielding networks.
- **Estimated CPDAG** (PC and GES, three sample sizes, ten seeds) on the four
  Gaussian networks. *This is the only panel that carries a deployment claim.*

### [RE-8] Elicitation
Questionnaire over each chain component containing a candidate treatment, asking
for **a causal order over the component**, not pairwise edge directions (pairwise
questions produce cycles at scale).

Sources above the line (deployment):
- A language model at temperature 0, prompt hashed and cached, in **two model
  families** and **two presentation orders**, so instrument variance and source
  variance are separate columns.
- The same with **scrambled variable names**, so recitation shows up as
  invariance.
- A **random control** matched to the model's answers on size, granularity, and
  distance from the treatment.
- The **empty set** as a status line.

Below the line (controls): the truthful recovering set, its corrupted versions,
and the committed `axisa3` rows.

**Run `[RE-1]` on the output the moment it exists.**

### [RE-9] Radius sweep on elicited knowledge
Per instance, report all four:
- `r_claim` — exhaustive retraction subsets to depth 3, with a `>3` sentinel.
- `r_hop` — the existing search, dispatch leg recorded.
- `φ₁` — the fraction of single retractions that break validity. This is the one
  quantity that varies where `r_claim` is flat (0.034 to 1.0, median 0.30 on the
  committed rows).
- **Identifiability radius** — retractions until the effect stops being
  identifiable by adjustment. Monotone, and defined on rows the current pipeline
  rejects outright.

Polynomial GAC gate only. Wall cap per instance; censoring recorded as a status,
never silently dropped.

### [RE-10] The audit
**Opened only after every deployment row is frozen and hashed.** Commission and
omission counted separately. The consistent-but-false fraction is the reporting
object. Four buckets (dangerous / held / safe / slack) with the rule-of-three
floor printed beside the dangerous rate. Test clustering of wrong assertions in
graph space and in prompt position. Report invariance under scrambled names.

---

## Tier 3 — completing the robustness axis

### [RE-11] Survival and cross-arm on the 831 real rows
Phase 1 of the robustness axis is **synthetic structure only**; real networks
appear only in the P7′ census. Running survival and the matched cross-arm design
on the real corpus is what makes the survival claim a claim about practice.
Code: `src/bkrobust/robustness/run_survival.py`, `run_crossarm.py`.

### [RE-12] N = 1000 re-run of the nine strata
Only the null cell (`flip, coverage=0.5, base_wrongness=0.25`) has been re-run at
N = 1000. Two strata carry the `⚠ unstable` marker in
`table_tau_comparisons.md` and must not be quoted alone until the full re-run
exists. The `n_eval ≥ 30` filter is not a conservative control — it conditions on
a quantity correlated with `r_val` and flips the verdict in opposite directions in
those two cells.

### [RE-13] The step hypothesis
Phase 2 established that statistical efficiency is *invariant* under truthful
knowledge, so the trade-off is identifiability against robustness, **a step, not a
frontier**. What is *not* established: does `r_val` fall monotonically with
proposal size beyond the minimum identifying set? Explicitly not claimed in the
session 7 report. It is the natural completion and it is cheap.

### [RE-14] Validate Conjecture 2 on real structure with E3
Every committed radius comes from `local_up_fast` or `e1_ladder` — both
**retraction-only**, both exact only if Conjecture 2 holds, both one-sided (too
large, never too small). **E2 and E3 have never been run on the real corpus.** E3
permits down-moves and assumes nothing; `r_E3 < r_E1` would refute Conjecture 2
with an explicit witness walk. Given that `r_hop` reaches 14 on real structure,
a spot-check on the largest radii is worth having before print.

---

## Tier 4 — deferred past the deadline

Recorded so they are not silently forgotten, and so reviewers' likely asks have a
stated answer.

- Knowledge harvested from source publications, with a verbatim quotation
  required for every assertion. **Note:** `paper/sections/evaluation-protocol.tex`
  currently promises this arm. Either remove it from that paragraph or move it
  above the line — it cannot stay as written. See `[RE-15]`.
- A human analyst round.
- The estimated-CPDAG panel on the **discrete** networks.
- A second, independent corpus. Everything currently traces to one pgmpy sdist.
- Sortability probes.
- The M-bias expansion.
- The RoCA port.
- `r_ε` on real structures (`NEXT.md` item 3).
- Structural characterisation of the law's exceptions (`NEXT.md` item 6).

---

## Documentation corrections that are not experiments

### [RE-15] Reconcile `paper/sections/evaluation-protocol.tex` with the plan
Two conflicts, both live:
1. It states that a CPDAG is estimated *for each network*; the design restricts
   estimated CPDAGs to the four Gaussian networks. Its own header `TODO` flags this.
2. It lists knowledge read from source publications as one of four sources; that
   arm is deferred past the deadline.

### [RE-16] Correct the practitioner sentence
[`hybrid.py:225`](../src/bkrobust/hybrid.py) `describe()` renders a radius as
*"survives any `{radius-1}` of the analyst's orientation claims being wrong."*
The search counts orientations in `G₀`, not asserted claims. On `paths` this
prints "13 claims" for an analyst who made one. Restate in commitment units.
Same conflation in the reports, which call `k_g0` "the knowledge size."

### [RE-17] Correct the `pathfinder` skip rationale
`stage0_claim_radius.py` says `pathfinder` yields "three censored rows." The
three `censored: true` records are `_meta` **run-level** markers, one per
coverage, each recording a 1,200 s wall timeout. Its three admissible rows are
`"ok"` and uncensored, and they do appear in the ladder. The real caveat — the
sweep was wall-capped, so those rows are a truncated sample — is recorded
nowhere.

### [RE-18] Confirm the Guo et al. citation
`paper/references.bib` has no Guo entry. The narrative anchors the choice of a
graph metric in prior work on distances between causal graphs and MPDAGs; the
exact reference must be supplied and verified by a human before it goes in. Do
not invent it.

---

## Housekeeping

- `results/axis_robustness/` carries ~115 MB of gzipped sample CSVs
  (`survival_samples.csv.gz` 48 MB, `xarm_samples.csv.gz` 49 MB,
  `null_samples.csv.gz` 14 MB). They are whitelisted deliverables and were merged
  as-is. If the repository must be shipped for anonymous review, decide whether
  these travel with it or move to an archive with checksums in the manifest.
- `run_analyse.py` and `table_tau_comparisons.md` sit at the repository root
  deliberately; `results/axis_robustness/PREREGISTRATION.md` explains why, and
  `build_tau_comparisons.py` writes to the root path. Do not tidy them without
  updating both.

---

## Round-1 review responses (branch `experiments/review-round-1`)

Four items below were run in response to an external review of the draft. Each has a
guard that reproduces published numbers before reporting new ones; where a guard would
fail, the script aborts rather than report.

- **[RE-12] closed for the second unstable cell.** `flip, coverage=1.0,
  base_wrongness=0.00` re-run at N = 1000 (`experiments/review1_unstable_stratum_n1000.py`,
  `results/axis_robustness/review1_unstable_cell/`). The published restricted-endpoint
  value of **-0.091** was a low-draw-count filter artefact: at N = 1000 the same endpoint
  gives **+0.537** [+0.484, +0.585]. The stratum is resolved positive; the reproduction
  check returns the published N = 200 value exactly. The remaining marked cell
  (coverage 0.5, base wrongness 0.25) stays a confirmed null.
- **Structural baselines added** (`experiments/review1_structural_baselines.py`,
  `review1_anchor_with_separation.py`). `analyse.PREDICTORS` and
  `build_tau_comparisons.PREDICTORS` never contained `separation`, so the published
  comparison never tested the radius against a truth-free structural statistic. On the
  anchor endpoints the radius beats separation in all six flip strata and **loses in all
  three tiered strata**; the two are equal on 875 of 1,978 instances and on all 240 of
  the `cov 1.0, bw 0.00` stratum. This is now reported in the paper.
- **The efficiency claim re-measured under the standard definition**
  (`src/bkrobust/demo/optimal_mpdag.py`, `experiments/review1_efficiency_definition.py`).
  `optimal_adjustment_set_mpdag` returns a set only when all DAG extensions agree, which
  makes the Phase-2 invariance result near-definitional. Implementing
  `O(x,y,G) = pa(cn) \ forb` (agreeing with the DAG-level optimal set on 3,428 cases,
  0 mismatches) and re-running the same bases and proposals: variance still does not move
  among proposals whose optimal set is GAC-valid (0 of 192 strata), and the two
  definitions agree on the set in 100.00% of proposals where both apply. The apparent
  variation in 24 of 192 strata comes entirely from 16 proposals where the graphical
  formula returns a GAC-invalid set.
- **Covering-graph connectivity** (`results/axisb2/review1_connectivity.json`). "Hop
  distance is a metric" is trivial unless the covering graph is connected. Exhaustive
  over all CPDAGs on 3 and 4 nodes: 196 spaces, 1,651 elements, 0 disconnected, 0
  undefined pairwise distances.

Still open and unchanged: [RE-1], [RE-3], [RE-7] to [RE-11], [RE-13], [RE-14].
