# Session summary — reviewer points 6 and 7

Branch `revision/experiments_points67`. One night session, 2026-09-22.
Orchestrated by an Opus 5 session with four Sonnet workers on disjoint files.

This is the authors' briefing: what was run, how long it took, what it found,
and — importantly — which parts of the reviewer's request could not be answered
as asked and why.

---

## 0. Three premises in the brief were already out of date

Checked against the repository before any compute was allocated:

1. **The 831-instance survival ranking was already done.** `experiments/todo_1`
   completed it at N = 1000 (725/725 shards). It was merged here and
   `make re11-verify` passes **26 of 26** checks. Re-running it would have spent
   roughly 80 CPU-hours reproducing hash-verified results. What was genuinely
   missing on that corpus was the *new baselines*, and those were run.
2. **The "33 cases of `r_val = 0`" are not in the 831-row corpus.** That corpus
   has zero rows with `radius == 0` and zero with an empty `Z`; its minimum
   radius is 1. The 33 live in `results/final_table/` and were *already*
   classified `degenerate` and already excluded from every median.
3. **The matched control already existed as elicited data.** `knowledge.json`
   carries ten conditions including `D_SCRAMBLED`; no API calls were needed.

---

## 1. What was run, and for how long

| pipeline | budget | wall time |
|---|---|---|
| Synthetic survival rerun, 9 shards (`results/axis_robustness_p6/`) | N = 1000 draws, 40 instances/cell, full 72 + 36 cell grid | **4.62 CPU-h** (≤5 concurrent) |
| Shard merge | 4,104 instances, 29,872 curves, 32.8M draws | 7.4 min |
| Paired resampling, real corpus (`results/axis_robustness_real_p6/`) | 10,000 resamples, clustered by network | ~25 min |
| Paired resampling, synthetic (`results/axis_robustness_p6_paired/`) | 10,000 resamples, clustered by instance, 54 comparisons | 2.8 min |
| Zero-radius audit (`results/zero_radius/`) | all 160 final-table rows recomputed | 7.4 min |
| Matched control, semantic arms | `D_SCRAMBLED`, `D_SCRAMBLED_72B` | 4,973 s + 4,956 s |
| Matched control, quantity arm | 20 seeds/network, 3 shards, 6,320 draw attempts | 5,103 + 19,183 + 7,231 s = **8.76 h** |

Nothing was mocked. Every reported number derives from committed result files.
Raw per-draw samples (1.4 GB) stay untracked and regenerate from the seeded
driver; curves and per-instance rows are committed.

**Reproducibility.** Because `derived_seed` is indexed by `rep`, the N = 1000
run's first 200 draws are the committed N = 200 run's draws *exactly*.
Rebuilding the old curves from them reproduced `survival_curves.csv` on **400 of
400** instance/depth cells checked, with `n_eval`, `n_contradictory` and `S` all
identical. The two budgets are therefore paired, and the instability the
reviewer saw is a precision effect rather than an artefact of reseeding.

---

## 2. Point 6 — the headline is the denominators

Every draw now logs attempted / rejected / accepted / actual-closure-change
separately. Over the full 32.8M-draw sweep:

| | flip | tiered |
|---|---|---|
| draws | 17,680,000 | 15,114,000 |
| accepted by closure | **15.6%** | 40.6% |
| rejected by closure | **84.4%** | 59.4% |
| accepted but **inert** (no closure orientation changed) | 0.0% | **52.1%** |

**This explains the reviewer's instability mechanically.** Flip acceptance falls
monotonically with component size — 31.8% (size 6) → 17.7% → 12.6% → **8.3%**
(size 12). So a nominal 200-draw grid point at size 12 rested on roughly **17
evaluable draws**. The committed analysis sensed this and guarded with
`n_eval ≥ 30`, but that guard conditions on `n_eval`, which is itself correlated
with `r_val` (τ from +0.318 to −0.438 across strata) — so the repair introduced
a selection effect that moved τ by up to 0.6 and flipped its sign in one
stratum. At N = 1000 the raw counts are five times higher, the guard excludes
far less, and τ reverts to its unfiltered value. **N = 1000 does not merely add
precision; it defuses the selection effect the N = 200 guard created.**

In the tiered arm, the `corruption_rate` grid's low end is largely a no-op:
100% inert at rate 0.00 and **84% inert at 0.05**, because `round(rate × n_k)`
rounds to zero claims. Survival is therefore pinned at 1.000 at those grid
points for every instance regardless of its radius, which compresses
between-instance variation in `AUC_rate`. The pre-registration's trigger T4
asserts exactly this at rate 0 and passes; nobody extended the check to 0.05.

### The new baselines: one of four is real

| requested baseline | verdict |
|---|---|
| single-claim deletion fragility (`phi_1`) | **genuinely new and informative** — 19 distinct values at coverage 1.0 |
| treatment-to-adjustment `separation` | usable but **undefined on 368/831** real rows; and on synthetic it *equals* `r_val` exactly wherever the generator realises its target separation |
| number of closure orientations (`k_g0`) | **a relabelling** — definitionally `shd_cpdag`, and equal to `n_k` wherever the claims are Meek-closed |
| claim-level radius (`r_claim`) | **uncomputable on this corpus** — constant at 1 on all 831 rows |

`r_claim ≡ 1` is not an empirical finding but a consequence of the gate that
built the corpus: `fast_gate` rejects with `no_atomic_perturbation_changes_validity`,
i.e. it admits a pair only when a single retraction breaks validity, which *is*
the condition `r_claim == 1`. Making it informative requires rebuilding the
corpus with that clause demoted from a filter to a status column — already
item **[RE-3]** in `docs/REMAINING_EXPERIMENTS.md`, expected yield 831 → ~2,490
rows. That is a separate experiment, not a bolt-on.

Degeneracies are **measured per stratum**, not asserted: `k_g0 == n_k` holds at
coverage 1.0 but not at 0.5, where the comparison is real.

### Paired resampling replaces eyeballing two intervals

On **synthetic** structure `r_val` paired-beats every non-degenerate baseline in
every stratum with the interval excluding zero (vs `phi_1` +0.25 to +1.36, vs
`n_k` +0.24 to +1.39, vs `shd_truth` +0.14 to +0.73).

On the **real** corpus the picture is weaker and must be reported as such:

- vs `phi_1`: radius wins at all three coverages, +0.298 / +0.392 / +0.491, all
  intervals excluding zero and all stable under leave-one-network-out.
- vs `n_k`: radius **loses**, robustly, at coverage 0.25 across all three
  base-wrongness replicates (Δ ≈ −1.3, interval entirely negative, LOO-stable).
- vs `k_g0`, `separation`: no reliable difference anywhere.
- **33 of 90 paired cells do not survive leave-one-network-out**, consistent
  with the corpus's known thin-strata leverage.

---

## 3. Point 7 — the LLM attribution does not survive its control

All arms share the questionnaire, the query seed and one interpreter
(`.venv/bin/python` 3.14.7).

| arm | informative | share at `r_val == 1` |
|---|---|---|
| `D_LLM` (treatment) | 53 | **86.8%** |
| `D_SCRAMBLED` (content control) | 42 | 76.2% |
| `D_SCRAMBLED_72B` | 31 | 83.9% |
| `D_RANDOM_MATCHED` (quantity control) | 700 | **74.9%** |

Paired over matched queries, clustered by network, 10,000 resamples:

| comparison | n | paired Δ | 95% CI |
|---|---|---|---|
| `D_LLM` − `D_RANDOM_MATCHED` | 53 | **+0.038** | **[−0.079, +0.134]** |
| `D_LLM` − `D_SCRAMBLED` | 29 | +0.103 | [0.000, +0.333] |

**Random orientations at the same \|K\| concentrate at radius 1 essentially as
hard as elicited LLM knowledge.** The paired difference contains zero. On this
evidence the concentration is a property of the real CPDAGs and the size of the
asserted knowledge set, not of what the model knows, and the paper's causal
attribution should be withdrawn or restated as an observation about the corpus.

**A caveat that cuts both ways.** The scrambled arm is not measurably less
accurate than the LLM arm: `compelled_control` accuracy 77.4% vs 67.3%, paired
per-network difference **+0.092, CI [−0.058, +0.259]**, containing zero. So a
null *semantic* contrast cannot be read as "knowledge does not matter" — the
manipulation may not have removed the knowledge. The quantity control does not
have this problem, since a uniform random draw carries no knowledge by
construction, which is why it carries the conclusion.

### Zero-radius: all 33 are genuine precondition violations

Recomputing every non-informative row of the final table:

- `VALID_EMPTY` = 0, `NO_SET_EXISTS` = 27, `INVALID_RETURNED` = **33**.
- The 33 split into **25** with no possibly-causal path from X to Y at all and
  **8** where `G₀` is not amenable. `optimal_adjustment_set_mpdag` returns `{}`
  without validating its own answer, and the pipeline records that as
  `r_val = 0`.
- They belong with the unidentified rows: **60 unidentified, not 27.** They were
  never "radius zero".

The distinction matters because **27 of the 66 informative rows also carry an
empty adjustment set** and a perfectly real radius. Verified directly:
`child CO2→DuctFlow` and `asia asia→tub` have `Z = {}` and are GAC-**valid** at
`G₀`, while `child Age→BirthAsphyxia` has `Z = {}` and is **invalid**. Valid and
invalid are indistinguishable from `|Z|` alone and separate only under the GAC
check — which is precisely why the three-way split was worth making.

---

## 4. Defects found and fixed this session

1. `cmd_compare` built the paired scrambled statistic with `NaN` for any
   instance informative but not at radius 1, then included it — every such
   instance poisoned the mean and the statistic returned `NaN`. Now `0.0`.
2. `control_accuracy` / `per_network_accuracy_diffs` assumed every network has a
   `compelled_control` block; networks that never reached a compelled question
   carry `None` and raised `TypeError`. Now skipped and counted, never imputed.
3. `paired_resample` parsed right-censored `r_claim` (written as an empty
   string) as a float and crashed. Censored rows are dropped and reported as a
   denominator, never imputed.
4. The shard merge wrote a 496 MB `survival_samples.csv.gz` outside the ignore
   rule; excluded before it reached git.
5. Two `munin4` rows disagree between `final_table` (`ok`, `r_val = 0`) and
   `final_table_eps_gt1` (`timeout`) — a provenance inconsistency between two
   committed tables, recorded but not resolved.

An upstream defect is **not** fixed here and is left for a scoped session:
`optimal_adjustment_set_mpdag` (`src/bkrobust/demo/evaluate.py`) should validate
its return value rather than emitting an invalid `{}`. That file produced
committed results and was out of this session's ownership.

---

## 5. What an author should take to the rebuttal

- **Point 6 is answerable and the answer is strong**, but the strongest material
  is not "we ran more draws" — it is the denominators. 84.4% of flip corruptions
  are rejected, effective sample size collapses to 8.3% at size 12, and 52.1% of
  accepted tiered corruptions change nothing. That is a concrete mechanism for
  the instability, and the N = 1000 rerun demonstrably removes it.
- **Do not claim four new baselines.** One is new and informative, one is a
  relabelling, one is undefined on 44% of rows, and one is constant by
  construction. Saying so plainly, with the gate argument attached, is stronger
  than presenting four comparators of which three are artefacts.
- **Point 7's control returns a negative**, and it should be reported at full
  weight. The radius-1 concentration is not shown to come from LLM knowledge.
- **The real-corpus ranking claim is weaker than the synthetic one**: the radius
  loses to `|K|` at coverage 0.25, and a third of paired cells do not survive
  dropping one network.

Every radius above assumes Conjecture 2 (hence Anti-Exchange Case B, verified
not proved); the error is one-sided, so radii can only be too large.
