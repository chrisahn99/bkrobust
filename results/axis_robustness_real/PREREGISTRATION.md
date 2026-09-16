# Pre-registration — [RE-11]: survival and the paired cross-arm on real graph structure

**Branch:** `experiments/todo_1`
**Written:** 2026-09-16, after a read-only orientation pass over the repository and a
cost probe, and **before any survival curve, AUC, rank correlation or cross-arm bin
existed**. No sweep code had been written when this file was committed.
**Corpus:** the committed real-network corpus, `results/axisa3/instances.jsonl`.
**Interpreter:** `/usr/bin/python3`, Python 3.9.6 — the interpreter recorded in
`results/axisa3/manifest.json`, which produced that corpus.

This file follows the project convention of prediction-before-looking. It exists so
that the τ values in `table_tau_real.md` cannot be read as post-hoc. **Nothing below
may be edited after the first result file lands in `results/axis_robustness_real/`;
corrections go in a dated appendix at the foot of this file**, in order, including
corrections that are embarrassing.

Companion pre-registration for the synthetic experiment this one ports:
`results/axis_robustness/PREREGISTRATION.md`, including its Appendices A–J, which are
binding on this design wherever they record a defect that could recur here.

---

## 0. The one-paragraph statement

An analyst holds a real network's oracle CPDAG `Ĉ`, asserts background knowledge `K`,
Meek-closes to `G₀`, and reads an adjustment set `Z* = optimal_adjustment_set_mpdag(G₀, X, Y)`
**once, at depth 0, and then holds it fixed**. We corrupt `K` and ask whether the set
they already committed to is still valid, decided by the **polynomial GAC predicate**
`is_gac_valid_mpdag(G, X, Y, Z*)` — no enumeration of `[G]` anywhere on the path. We ask
whether the **worst-case** breakdown radius `r_hop` ranks the resulting **average-case**
survival, and whether the baselines `|K|` and `SHD(G₀, truth)` do so. The synthetic
answer is `table_tau_comparisons.md`; this is its real-structure counterpart. **The two
corpora are never pooled and no real row is ever compared with a synthetic row.**

---

## 1. The frame

### 1.1 Definition

The **primary frame** is every row of `results/axisa3/instances.jsonl` with
`admissible == true`: **831 rows, 25 networks, 543 distinct `(network, X, Y)` pairs**,
spanning coverage ∈ {1.0, 0.5, 0.25} (543 / 182 / 106). Using this committed corpus,
rather than re-selecting pairs, is what makes the result speak to §6.2 of the paper.

The frame is built, frozen and **hashed before any survival curve exists**. It is
written to `results/axis_robustness_real/frame.jsonl`, its SHA-256 recorded in
`results/axis_robustness_real/frame_hash.json` together with the SHA-256 of the source
`instances.jsonl`, and every later stage re-checks both hashes before writing a row.

### 1.2 Columns carried on every frame row

`network`, `x`, `y`, `coverage`, `largest_component_size`, `component_size`,
`n_k` (**stated size** `|K|`), `k_g0` (**commitment size**), `separation` and
`separation_status`, `dispatch_leg` (`local_up_fast` or `e1_ladder`), `gate`
(`fast_gate`, always), `radius_committed` and `r_status_committed`,
`g0_undirected_edges_recomputed`, `z_star`, `z_size`, `assumes`, and the SHA-256 of the
source row.

Naming follows `docs/PAPER_NARRATIVE.md` §11: `k_g0` is the **commitment size**, not
"the knowledge size"; `|K|` is the **stated size**. The legacy column name `r_val` is
retained in machine-readable files because the existing analysis code reads it, and is
rendered `r_hop` in every table and figure.

### 1.3 Status, never filter

Every knowledge-dependent condition is a **status column on the frame, never a filter**,
so every arm reports the same denominator. In particular `optimal_set_undefined`,
`o_g0_extensions_intractable`, `k_contradictory`, `censored_wall_cap` and
`unresolved_all_contradictory` are statuses. Nothing is dropped.

**The two arms have structurally different denominators and this is stated, not
smoothed over.** The flip arm is a function of `coverage`, so its denominator is the
**831 frame rows**. The tiered arm generates its own knowledge from the network's
temporal order and has no coverage parameter at all, so its denominator is the **543
distinct `(network, X, Y)` pairs** the frame covers. Both are printed beside every
aggregate.

### 1.4 A defect in the committed corpus, recorded here because the frame must work around it

`benchmarks/measure.py::evaluate` assigns `g0_undirected_edges` **only** on the
`o_g0_extensions_intractable` rejection path. On the admissible path it is never
assigned, so all 831 admissible rows carry the dataclass default `0`, which is not a
measurement. Recomputation gives non-zero values at partial coverage (e.g. `insurance`
at coverage 0.5 leaves 6 undirected edges in `G₀`). This is [RE-5b]. The frame carries
the **recomputed** column under a distinct name, `g0_undirected_edges_recomputed`; the
committed file is not modified, and no previously published number changes, since no
published number used that column.

### 1.5 Unit of analysis

**The network, not the pair.** The intraclass correlation of the `r = 1` indicator is
0.40 over the 25 networks (design effect ≈ 12.7, effective n ≈ 65 against 831 rows).
Accordingly:

- every bootstrap is a **cluster bootstrap over networks** (§6.2);
- **no per-row interval is printed anywhere**;
- every aggregate is given **both pair-weighted and network-weighted**;
- the **per-network table is printed beside every aggregate**, because networks in this
  corpus differ from each other more than the corpus differs from Erdős–Rényi
  (`report_real_graphs.md` §5.2).

---

## 2. The SHD baseline, and why it needs base wrongness

`select_knowledge(dag, cpdag, coverage)` draws the analyst's claims from
`knowledge_to_recover(dag, cpdag)`, which reads the **true** orientation off the
ground-truth DAG. Measured in the orientation pass, on this corpus:

- at coverage **1.0**, `G₀` **is** the ground-truth DAG, so `SHD(G₀, truth) ≡ 0` on
  every admissible row;
- at coverage 0.5 and 0.25, `SHD(G₀, truth)` equals **exactly** the number of undirected
  edges `G₀` still carries — it measures *incompleteness*, never *wrongness*.

So without intervention the SHD baseline is rigged by construction and a reviewer would
say so. The synthetic design solved this with a **base wrongness** `b ∈ {0.00, 0.10,
0.25}` applied by `synth.knowledge.flip` before the sweep
(`results/axis_robustness/PREREGISTRATION.md` §5.1), and the real corpus already carries
the same machinery in `results/axisa3/wrong_knowledge.jsonl` (three draws per rate).
**That is the mechanism used here**, and it is recorded on every row as
`base_wrongness` and `analyst_replicate`.

**Verification gate, binding:** before a single τ is computed, the analysis asserts that
`shd_truth` actually **varies within each stratum** where it is claimed to be defined,
and reports `undefined (predictor constant)` — never `τ = 0` — where it does not.

### 2.1 Analyst replicates, and why they exist

`select_knowledge` returns a **network-level** claim set: on this corpus `K` is a
property of `(network, coverage)`, not of the query. Since `flip` preserves list length,
**`|K|` cannot vary within a `(network, coverage)` cell under any base wrongness.**
That is a fact about the corpus and is reported as one, not engineered away.

`SHD(G₀, truth)` would be equally constant if a single corrupted `K` were shared by
every pair of a network. To give it within-network variation without pretending that
each query has its own analyst, the design draws **`R` independent base-wrongness
corruptions per `(network, coverage, base_wrongness)`** — "three analysts, same coverage,
same error rate" — each shared across that network's pairs:

- `R = 1` for `b = 0.00` (the corruption is the identity; a second draw would be a
  duplicate, not a replicate);
- `R = 3` for `b ∈ {0.10, 0.25}`, mirroring `wrong_knowledge.jsonl`'s three draws.

Seeds are derived by SHA-256 from `(network, coverage, "base", b, replicate)` — never
Python's salted `hash()`.

**Consequence, stated in advance:** in the `b = 0.00` strata `shd_truth` is constant
within each network by construction, and at coverage 1.0 it is constant at 0 across the
whole corpus. That stratum will therefore report `undefined (predictor constant)` for
`shd_truth`, exactly as the synthetic `flip cov=1.0 bw=0.00` stratum does. This is
predicted, not discovered.

---

## 3. Strata

Nine, mirroring the nine pre-registered synthetic strata so the two tables can be read
side by side — **and never pooled with them**.

| # | arm | stratum | synthetic counterpart |
|---|---|---|---|
| 1 | flip | coverage = 1.0, base_wrongness = 0.00 | flip cov=1.0 bw=0.00 |
| 2 | flip | coverage = 1.0, base_wrongness = 0.10 | flip cov=1.0 bw=0.10 |
| 3 | flip | coverage = 1.0, base_wrongness = 0.25 | flip cov=1.0 bw=0.25 |
| 4 | flip | coverage = 0.5, base_wrongness = 0.00 | flip cov=0.5 bw=0.00 |
| 5 | flip | coverage = 0.5, base_wrongness = 0.10 | flip cov=0.5 bw=0.10 |
| 6 | flip | coverage = 0.5, base_wrongness = 0.25 | flip cov=0.5 bw=0.25 |
| 7 | tiered | n_tiers = 2 | tiered n_tiers=2 |
| 8 | tiered | n_tiers = 3 | tiered n_tiers=3 |
| 9 | tiered | n_tiers = 4 | tiered n_tiers=4 |

**What real structure cannot support, named rather than substituted:**

- The synthetic strata are also crossed with *component size* (6, 8, 10, 12) and
  *separation* (min, median, max), which are **design parameters of a generator**. Real
  structure does not let those be set; they are **observed** and are carried as columns
  and reported in the per-network table, not used as strata.
- **Coverage 0.25 has no synthetic counterpart** and is therefore **not** one of the
  nine. Its 106 rows are run and reported as a **supplementary stratum**, clearly
  labelled.
- **No marginal rate is ever compared across coverage levels.** Admissible instances
  collapse 543 → 182 → 106 as coverage falls and the survivors have larger separation,
  so any such comparison is a composition effect ([RE-4], [RE-6],
  `report_real_graphs.md` §5.3). Any coverage contrast is computed **only** on the
  **105 `(network, X, Y)` triples admissible at all three coverages** (8 networks), is
  labelled as the matched-population contrast, and is supplementary.

---

## 4. The two arms, and the axis each lives on

**The two arms are never merged row to row.** Their intensity parameters count different
things, and `results/axis_robustness/PREREGISTRATION.md` Appendix G records what happens
when they are pooled anyway.

### 4.1 Flip arm — axis is the claim depth

- `K = select_knowledge(dag, cpdag, coverage)`;
  `K_b = flip(K, rng(network, coverage, b, replicate), b)`;
  `G₀ = apply_orientations(cpdag, K_b)`; `Z* = optimal_adjustment_set_mpdag(G₀, X, Y)`,
  read **once** and held fixed.
- Corruption: `flip(K_b, rng(instance, "flip", d, rep), d / |K_b|)`, which reverses
  exactly `d` claims. The realised count is **asserted equal to `d`**, as the synthetic
  driver does; a mismatch is a crash, not a silent row.
- **Depth grid.** `D = sorted(set(clamp(round(f · |K|), 1, |K|)) for f in 0.1 … 1.0)`.
  This is exactly the set of depths the pre-registered `AUC_frac` endpoint averages
  over, so nothing is measured that the endpoint discards and nothing the endpoint needs
  is missing. For `|K| ≤ 10` it is the full dense sweep `1 … |K|`; for `|K| > 10` it is
  the 10 fractional grid points. The deviation from the synthetic sweep (which was dense
  because its `|K|` was small) is recorded per row as `depth_grid = "frac10"`.
- **Endpoint.** `AUC_frac` — `survival.auc_frac` **imported, not reimplemented**, so the
  nearest-defined-depth and tie-toward-smaller-depth rules cannot drift.

### 4.2 Tiered arm — axis is the corruption rate, never the depth column

- `K_ref = tiered(dag, cpdag, rng(network, n_tiers, "ref"), n_tiers, 0.0)`;
  `G₀ = apply_orientations(cpdag, K_ref)`; `Z*` read off it and held fixed. This is a
  **separate instance population** from the flip arm.
- Corruption: `tiered(dag, cpdag, rng(instance, "tiered", rate, rep), n_tiers, rate)`
  over the 11-point grid `0.00, 0.05, …, 0.50` — the same grid as the synthetic arm.
- **The axis is `corruption_rate`, fixed now, from the start.** The tiered `d_claims`
  column counts assertions added or reversed but **not removed**, and relocating a node
  into a shared tier silently deletes a cross-tier assertion, so the depth column's zero
  bin is contaminated (`results/axis_robustness/PREREGISTRATION.md` Appendix E). The
  column is still recorded, marked `d_claims_defective_do_not_bin`, purely so the defect
  can be re-measured on real structure; **no analysis bins on it**.
- **Endpoint.** `AUC_rate` — the mean of `S` over the 11 rates.
- **Binding sanity check (falsification trigger T4).** `S = 1.000` **exactly**, for
  **every instance individually**, at `corruption_rate = 0`. Checked before any τ is
  computed. A single violation halts the arm and is reported.

### 4.3 What is recorded on every sample

`instance_id`, `arm`, `gate`, the stratum keys, the grid point, `rep`, `seed`, `status`
(`ok` / `corrupted_k_contradictory` / `censored_wall_cap`), `survived`, and
`symdiff_proxy_not_distance` — labelled a proxy everywhere it appears, never a distance
(`THEOREMS.md` §10–§11).

---

## 5. Endpoints, draw count and censoring

### 5.1 Both endpoints, always

For each curve, both the **raw** endpoint (`AUC_frac` / `AUC_rate`) and the
**conservative `_usable`** variant restricted to grid points with `n_eval ≥ 30` are
computed and **both are reported**. Following the brief, the `_usable` figure is the one
quoted in any summary claim.

**This is done with Appendix J in view, not in ignorance of it.** Appendix J withdrew
the claim that `_usable` is a *conservative* control: it conditions on `n_eval`, which
is correlated with `r_hop`, and at N = 200 it flipped the verdict in two synthetic cells
in opposite directions. Appendix J also established the clean diagnostic: **at N = 1000
the two endpoints converged to 0.0002 on the same 220 instances**, so divergence is a
symptom of too few draws rather than evidence that either is safer. Therefore:

- the **gap** `|AUC_raw − AUC_usable|` and the two τ verdicts are reported side by side
  for every stratum, as a diagnostic;
- **falsification trigger T3**: if the two endpoints disagree in sign or in
  CI-excludes-zero verdict in more than 2 of the 9 strata *at N = 1000*, those strata are
  reported **unresolved**, not resolved by choosing an endpoint.

### 5.2 Draw count

**N = 1000 per grid point, everywhere, from the start.** [RE-12] established that N = 200
is not enough for the `n_eval ≥ 30` filter to be inert, and real structure has thinner
strata than synthetic, so it is more exposed, not less. The cost probe puts the campaign
at ≈ 80 CPU-hours on a 10-core machine, which is affordable, so no cell is run at lower
precision. `n_draws` is recorded **per row** regardless, and any cell that ends up below
1000 for any reason is named explicitly in the report.

`n_eval` is recorded **per grid point**. The **contradiction rate is recorded per grid
point** — the fraction of corruptions Meek closure rejects outright. On synthetic
structure that was 0.627 at depth 1; whether real structure is more or less
self-revealing is a finding either way.

### 5.3 Censoring and sentinels

- `UNREACHED = −1` is a **status**, never a radius. It is never averaged, never plotted
  on a numeric axis, never fed to a correlation, and instances carrying it are **counted
  out loud**.
- A wall-capped unit of work is status `censored_wall_cap` and carries
  `wall_until_timeout_s`. **It never carries a key that a real measurement uses** — in
  particular never `seconds`, which is the exact shape of the defect that cost session 4
  a day.
- A corrupted `K` admitting no consistent MPDAG is `corrupted_k_contradictory`. It is
  **not** a survival-0 outcome; `S` is computed over non-contradictory samples, and the
  contradiction rate is its own reported series. The pre-specified sensitivity
  `S_contra_as_fail` (contradictions scored as failures) is computed everywhere. If the
  two disagree materially that is reported prominently rather than resolved by choosing.
- An undefined separation is `None` with a status string, never `0` and never `−1`.

### 5.4 Wall caps

Per-shard wall cap **5 hours**; per-radius-computation cap **300 s** (as the committed
corpus used). Exceeding either produces a `censored_wall_cap` status with
`wall_until_timeout_s` and never a silent drop. `pathfinder` is the network expected to
approach the shard cap: its radius computes in 0.16 s, but a single corrupted Meek
closure on its 85-vertex component costs ≈ 1 s.

---

## 6. Predictors, statistics and the analysis plan

### 6.1 Predictors under test

| predictor | definition | unit | known structure on this corpus |
|---|---|---|---|
| `r_hop` (`r_val`) | `hybrid.breakdown_radius(...).radius` at the instance's own `G₀` | hops | varies within network on 15 of 25 networks at coverage 1.0 |
| `shd_truth` | `SHD(G₀, true DAG)` | edges | constant within `(network, coverage, b, replicate)`; ≡ 0 when `b = 0` and coverage = 1.0 |
| `n_k` | `|K|`, the **stated** size | claims | **constant within `(network, coverage)` by construction** |
| `k_g0` | the **commitment** size | orientations | constant within `(network, coverage)`; carried as a secondary baseline |

### 6.2 Statistics

- **Kendall τ-b** (`scipy.stats.kendalltau(variant="b")`). Ties are material — `r_hop` is
  small-integer valued — so never τ-a.
- **Cluster bootstrap over networks.** 10,000 resamples, seed 0, explicit
  `numpy.random.Generator`. Each resample draws **25 networks with replacement** and
  takes **all** rows of each drawn network; τ is recomputed on the resampled rows.
  Percentile interval at 2.5 / 97.5. Resamples that leave a predictor constant yield NaN
  and are **counted and reported**, never silently dropped.
- **Stratified, never pooled.** No pooled τ is a headline.
- **Leave-one-network-out sensitivity on every headline τ**, because within-network
  radius variation is thin: 10 of 25 networks are constant in `r_hop` at coverage 1.0,
  while `paths` alone spans `r = 1…14` in 20 rows and therefore carries large leverage.
  **Any stratum whose verdict flips when a single network is removed is marked
  `⚠ single-network leverage` in the anchor table and is not quoted alone.**
- **Effect sizes by radius bucket** (`r = 1, 2, 3, 4, 5+`), median and IQR of the
  endpoint, with per-bucket instance counts printed. `UNREACHED` is its own bucket and
  is never placed on the numeric axis.
- **Decomposition identity check**, on every curve row:
  `S_contra_as_fail(d) = (1 − contradiction_rate(d)) · S(d)`, to 1e-12
  (falsification trigger T6).
- **The undefined-separation stratum** (368 of 831 rows) is analysed on its own and
  reported as its own row, never pooled into a separation number.
- Any **discordance spotlight** must have its per-grid-point `n_eval` printed beside it
  before it is written up. Session 8 withdrew a spotlight that rested on depths with
  `n_eval` of 35, 10, 1 and 1; a spotlight here requires **every** contributing grid
  point to have `n_eval ≥ 30`.

### 6.3 Gate discipline

`benchmarks.measure.fast_gate` **exclusively**, imported and never reimplemented, with
the literal string `fast_gate` on every row. `synth.runner.gate` is **not used**: it
calls the exponential `all_valid_adjustment_sets_mpdag`, and it is not equivalent —
16 disagreements in 46,800 cases, all one direction, biasing **against** the hypothesis
under test. A number without its gate label is not quotable.

### 6.4 The assumption string

Every radius row carries `HybridResult.assumes` verbatim: *"Conjecture 2 (hence
Anti-Exchange Case B, verified not proved)"*. All upward searches are exact iff
Conjecture 2 holds and the error is **one-sided** — radii can only be too **large**,
never too small (`THEOREMS.md` §4c, §6, §8). The string travels with every number
reported, including into `table_tau_real.md`.

`all_valid_adjustment_sets_mpdag` is never on a decision path.

---

## 7. The paired cross-arm design

Session 7 claimed correlated errors are harder to detect than uniform ones and
**retracted it** (Appendix E.3): the two arms relocate nodes versus reverse claims, so an
unpaired comparison at nominal intensity is meaningless. Session 8 rebuilt it properly
(Appendix I) and **that** is the design ported here. The retracted version is not
reconstructed under any circumstances.

- Every instance is constructed **once**, the tiered way:
  `K_ref = tiered(dag, cpdag, rng, n_tiers, 0.0)`, `G₀`, `Z*` read off it.
- That same `(Ĉ, K_ref, G₀, Z*)` is corrupted **both ways** — `tiered(…, rate)` and
  `flip(K_ref, …, rate)`.
- Both are scored on the **shared intensity axis**
  `intensity = |dir(G₀) Δ dir(G)| / |dir(G₀)|`, in bins of width 0.05. The brief's
  "fraction of the skeleton altered" is identically zero, since corruption moves
  orientations and never adjacency.
- Comparison is **paired within instance**, per bin, requiring **≥ 30 non-contradictory
  samples in both arms** for that instance and bin.
- Reported per bin: matched instance count, mean and median paired difference, sign
  test, Wilcoxon signed-rank, and a **paired bootstrap over networks**.
- Reported **only over the intensity range where both arms actually carry data**. Flip
  has a floor near 0.10 because reversing one claim changes at least two orientations;
  tiered concentrates below it. The achieved-intensity histogram of each arm is printed.
- Also reported, as the endpoint that answers the practitioner's question (Appendix I.3):
  the **silent-failure rate** `(1 − contradiction_rate) · (1 − S)`.

**For engineering reasons the two arms are sampled in separate shards.** The instance is
still constructed once in the sense that matters: both shards recompute it from the same
deterministic derivation, and **each shard writes the SHA-256 of
`(K_ref, dir(G₀), Z*)`; the pairing step refuses to pair any instance whose two hashes
differ.** That is a stronger guarantee than shared process memory, and it is checked,
not assumed.

**Falsification trigger T5.** If fewer than 3 intensity bins carry ≥ 15 matched
instances in both arms, the pairing is declared **not workable on real structure** and
reported as such, with the achieved-intensity histograms as evidence. The unpaired
comparison is **not** substituted.

---

## 8. Predictions — the part that must not move

Recorded before any survival curve, AUC, τ, contradiction rate or cross-arm bin existed
on real structure. Written by the orchestrator, not by a subagent.

- **P1 — `r_hop` ranks survival on real structure.** τ_b(`r_hop`, `AUC_*_usable`) > 0 in
  the **majority** of the nine strata. Point predictions: **positive in 6–9 of 9**, with
  point estimates **+0.15 to +0.55 in the flip strata** and **+0.30 to +0.70 in the
  tiered strata** — that is, *weaker than synthetic* (which gave +0.63 to +0.79). Because
  the bootstrap is now clustered over 25 networks rather than over ~240 instances, the
  intervals will be **materially wider**, and I predict the CI **excludes zero in only
  4–7 of 9 strata**, against 8 of 9 synthetic. **I expect the real-structure result to
  be weaker than the synthetic one, and I am saying so before looking.**
- **P2 — `r_hop` outranks `|K|`.** τ_b(`r_hop`) > τ_b(`n_k`) in the majority of strata
  where both are defined: predicted **6–9 of 9**. Caveat recorded in advance: `n_k` is a
  network-level constant on this corpus, so it can rank networks but not queries, and a
  win over it is worth less here than it was synthetically. The report must say this.
- **P3 — `|K|` is not inert.** H4's "inert baseline" framing was wrong synthetically;
  `n_k` was *negatively* associated (τ ≈ −0.38 to −0.42 in every tiered stratum). I
  predict the same sign on real structure in the tiered strata, and an **unstable sign**
  in the flip strata. Predicted: `n_k`'s CI excludes zero in **at least 3 of 9**.
- **P4 — `SHD` is undefined by construction in exactly one stratum.** `shd_truth` is a
  constant predictor in **flip coverage = 1.0, base_wrongness = 0.00** and nowhere else
  among the nine. Elsewhere I predict it is defined but has **no consistent sign** across
  strata, as synthetically.
- **P5 — real structure is MORE self-revealing than synthetic.** The flip-arm
  contradiction rate at depth 1 was **0.627** synthetically. Real CPDAGs are far more
  compelled (median undirected fraction 10.7%), so a reversed claim more often collides
  with a compelled edge. I predict **0.65–0.90 at depth 1 on real structure**, and
  directionally **higher than synthetic**. If it comes out *lower*, that is the more
  interesting result and it is reported as the headline of that section.
- **P6 — the cross-arm ordering survives the port.** At matched intensity, tiered
  contradiction rate exceeds flip in most populated bins, and the **silent-failure rate
  is lower for tiered than flip** — i.e. uniform, independent errors remain the more
  insidious. I predict the same direction as Appendix I, over a **narrower** usable
  intensity band (≈0.10–0.50 rather than ≈0.15–0.65), because real `|dir(G₀)|` is much
  larger, which compresses the normalised intensity.
- **P7 — the undefined-separation stratum is more robust.** The 368 rows with no member
  of `Z*` in `X`'s component will show **higher** median survival AUC than the 463 rows
  with a measured separation, and a **lower** share at `r_hop = 1`. Held weakly.

### 8.1 Falsification triggers

Restated compactly; these are the conditions under which the campaign **stops and
reports a null** rather than redefining the analysis.

| id | trigger | consequence |
|---|---|---|
| **T1** | τ_b(`r_hop`, `AUC_*_usable`) not positive in ≥ 5 of the 9 strata | the real-structure ranking claim is **not established**; the report leads with the null |
| **T2** | τ_b(`r_hop`) ≤ τ_b(`n_k`) in ≥ 5 of the strata where both are defined | report that the radius does not outrank `|K|` on real structure |
| **T3** | raw and `_usable` endpoints disagree in sign or verdict in > 2 of 9 strata **at N = 1000** | those strata are reported **unresolved**; no endpoint is chosen to break the tie |
| **T4** | tiered `S ≠ 1.000` exactly at rate 0 for any instance | tiered arm halts; defect reported before any τ |
| **T5** | < 3 intensity bins with ≥ 15 matched instances in both arms | pairing declared **not workable**; unpaired comparison **not** substituted |
| **T6** | decomposition identity fails beyond 1e-12 on any row | sweep halts |
| **T7** | base-wrongness-zero radii do not reproduce `results/axisa3/instances.jsonl` exactly | the frame is wrong; campaign halts |

No trigger may be weakened after a result is seen. If one fires, the appendix records it
and the report leads with it.

---

## 9. Determinism, provenance and incrementality

- `num_workers: 1`, `random_seed: 0`, CP-SAT single-threaded, **no global RNG anywhere**.
  Every draw goes through an explicit `numpy.random.Generator` seeded by a SHA-256
  derivation from `(instance_id, process, grid_point, rep)` — never Python's salted
  built-in `hash()`.
- **Bit-identical output across `PYTHONHASHSEED` 0 and 12345 is checked, not assumed**,
  for every stage, and the check is committed.
- Every results subtree carries `manifest.json` with git SHA, dirty flag, environment,
  solver parameters, `RADIUS_CONVENTION` and a determinism block.
- **Append-only JSONL, one row per unit of work**, keyed by a deterministic id derived
  from `(network, X, Y, coverage, arm, intensity, seed)`. **No row is ever rewritten in
  place** — session 8 lost an entire survival sweep to a worker that relaunched and
  overwrote its own outputs.
- **The unit of work is one `(instance, arm, grid point)` cell, and that is the row
  written.** Each row carries `n_draws`, `n_eval`, `n_contradictory`, `n_survived`, `S`,
  `S_contra_as_fail`, and — for the cross-arm — the median measured intensity of the
  cell's non-contradictory draws and its bin. **Individual Monte-Carlo draws are not
  persisted.** The synthetic sweep persisted 3.14M sample rows because it later had to
  re-bin the tiered arm off a defective depth column (Appendix E); this design bins on
  `corruption_rate` from the start, so no analysis here needs a raw draw, and every
  quantity the analysis plan of §6 uses — including the decomposition identity, the
  per-grid-point contradiction rate and the `n_eval ≥ 30` filter — is computable from the
  cell row. Persisting ≈ 40M draw rows would cost several gigabytes and buy nothing that
  is used. This is a deliberate deviation from the synthetic sweep's file layout and is
  recorded here rather than discovered later.
- **Work is sharded** at `(arm, network, cell, replicate)` granularity. A shard owns
  exactly one file and writes `_done/<shard_id>.json` — with row count and file
  SHA-256 — **only when it is fully finished**.
- **Partial shard output is cleanly discardable by construction.** A shard with no
  completion marker has its file truncated and is redone from scratch on re-invocation.
  Nothing is resumed mid-file. Re-invoking any stage is therefore safe to repeat and
  skips completed work.
- **Existing modules that produced committed results are not modified.** All new code
  lands in new modules alongside them; the precedent is commit `748f2e1`.

---

## 10. What this design does not claim, fixed in advance

- **Oracle CPDAG throughout.** Every CPDAG here is computed from the true DAG. Real
  discovery on finite data returns a sparser skeleton, and missing weak edges is exactly
  what breaks back-door blocking. This is the corpus's largest scope limit
  (`report_real_graphs.md` §7) and nothing here touches it.
- **One corpus.** Everything traces to a single pgmpy 1.0.0 sdist.
- **No comparison to the truth.** The endpoint is structural survival of the committed
  set across knowledge states, not whether `Z*` is correct in the true DAG.
- **No pooling with synthetic.** Not in a table, not in a figure, not in a sentence.
- **No claim of perfect ranking.** The endpoint is a rank correlation with a bootstrap
  interval. The word "perfectly" appears nowhere in this analysis plan, and "dominates"
  will not be claimed unless it is measured — Appendix J.6 records what happened last
  time that word was used ahead of the data.

---

# Appendix A — 2026-09-16, the rate-based base wrongness is inert on most of this corpus

Written after the frozen frame landed and **before any survival curve, AUC or τ
existed**. It amends §2 by *adding* a level, and changes nothing that was already
fixed. Recorded as a dated appendix rather than an edit to §2, per the rule at the
head of this file.

## A.1 The measurement

`select_knowledge` returns a **greedy-minimal** generator set
(`demo/example.py::knowledge_to_recover`), so `|K|` is small on most of this
corpus: 1 or 2 on 24 of the 50 `(network, coverage)` cells, and 3 or less on 32 of
them. `synth.knowledge.flip` reverses exactly `round(rate · |K|)` claims. Counting
the claims actually reversed, per cell, across all three analyst replicates:

| level | cells where **nothing at all** is reversed |
|---|---|
| base wrongness 0.10 | **38 of 50** |
| base wrongness 0.25 | **24 of 50** |

Worked examples: `|K| = 1` (Acid_1996, Didelez_2010, hailfinder, paths, and six
partial-coverage cells) gives `round(0.25 · 1) = 0`; `|K| = 2` gives
`round(0.5) = 0` under Python's banker's rounding. Only `arth150` (29),
`diabetes` (26), `pathfinder` (79), and the cells with `|K| ≥ 6` are reached by
`b = 0.10` at all.

**Consequence if left alone.** On those cells `G₀` is bit-identical to the
`b = 0.00` graph, so `SHD(G₀, truth)` stays identically zero and the baseline
comparison is rigged by construction on exactly the cells §2 set out to protect.
The `b = 0.10` and `b = 0.25` strata would also contain large numbers of rows that
are exact duplicates of the `b = 0.00` rows.

## A.2 The amendment

**A fourth base-wrongness level is added: `bw_abs = 1`, "exactly one claim
reversed", with 3 analyst replicates.** It is expressed as
`flip(K, rng, 1/|K|)`, so **no new corruption operator enters the design** — it is
the same operator at a rate chosen to reverse one claim on every cell, including
the `|K| = 1` and `|K| = 2` cells the rate-based levels cannot touch.

- It is a **supplementary stratum**, not one of the nine. The nine strata of §3
  are unchanged, so the real and synthetic anchor tables still line up row for row.
- It exists for one stated purpose: to give `SHD(G₀, truth)` a non-degenerate
  comparison on every network, so §2's verification gate can be met corpus-wide.
- **Prediction, recorded now.** On the `bw_abs = 1` stratum `shd_truth` will vary
  across networks and across analyst replicates, and τ_b(`shd_truth`, AUC) will be
  **defined** (non-constant predictor) where the rate-based strata leave it
  undefined. No directional prediction is made for its **sign** — the synthetic
  result gave `shd_truth` no consistent direction across strata, and I have no
  reason to expect one here.

## A.3 The inert cells are run anyway, not skipped

Where `round(b · |K|) = 0` the cell is a bit-identical duplicate of the `b = 0.00`
cell. It would be cheaper to detect that and reuse the `b = 0.00` result. **It is
run in full instead.** The duplicated cells are precisely the small-`|K|` ones,
which are also the cheapest — every expensive network (`pathfinder`, `diabetes`,
`arth150`) has a non-zero realised wrongness at both rates — so the saving would be
negligible, while a reuse-by-join would introduce exactly the kind of silent
encoding bug this project has been caught by before. Every row carries
`n_claims_actually_wrong` and a boolean `bw_is_inert`, and the report states which
strata are near-duplicates of `b = 0.00` and why.

## A.4 Scope of the determinism claim

Bit-identical output across `PYTHONHASHSEED` 0 and 12345 is verified on every shard
kind **excluding the three wall-clock fields** `r_search_seconds`,
`r_total_seconds` and `wall_until_timeout_s`, which are timings and cannot be
reproducible. The `*.cells.jsonl` files carry no timing at all and are verified
**byte-identical**. Verified on five shards spanning all three arms before the full
sweep started; the check is repeated over the whole campaign at the end.

---

# Appendix B — 2026-09-16, two defects during the sweep, and what caught them

Recorded in order, while the sweep was running, rather than tidied up afterwards.

## B.1 A measurement limit was being reported as a structural property

`run_real_survival.run_xarm_shard` computed the shard's status as

```python
n_dir_g0 = len(g0.directed_edges) if g0 is not None else 0
g0_status = g0_reason if n_dir_g0 else "dir_g0_empty"
```

When `build_g0` returned `None` — because `G₀` still carried more than
`MAX_G0_UNDIRECTED_FOR_EXTENSIONS = 12` undirected edges, i.e. the
**`o_g0_extensions_intractable` measurement limit** — `n_dir_g0` fell to 0 and the
shard was stamped `dir_g0_empty` instead. Those are different facts:
`dir_g0_empty` says a graph exists and has no directed edge, so the shared
intensity axis has no denominator; `o_g0_extensions_intractable` says we could not
afford to evaluate the pair. Keeping the two apart is the whole point of §1.3, and
`benchmarks/measure.py` says so in its own module docstring.

**No measurement was wrong** — the shard was correctly blocked from sampling either
way, and no survival number came from a mislabelled row. Only the *label* was
wrong, which is precisely the failure mode this project keeps catching: a wrong
answer that looks entirely plausible in the output file.

**What caught it.** Reading the completion markers of the six `pathfinder`
cross-arm shards by hand because they finished in under a second when the budget
said hours, and noticing that `dir_g0_empty` was an odd thing to say about a
network whose `G₀` has 195 directed edges.

**Fix.** The two conditions are now tested separately, with the comment that says
why. Every shard carrying the wrong label — 6 of them, all blocked, all
zero-cell — had its marker removed and was re-run.

## B.2 Incident: an operator error, and the property that made it free

While re-running those six shards the orchestrator invoked
`run_real_survival list --kind xarm`, which returns **every** pending cross-arm
shard, not the six intended, and piped all of them into a serial loop. That loop
ran alongside the eight-worker pool, so the same shard could in principle have been
written by two processes at once.

It was noticed within four minutes, the stray process was killed, and **every
cross-arm shard touched in that window — 64 of them — had its marker and both its
files deleted and was returned to the pending list.** Nothing was inspected to
decide whether a given file was salvageable, because the incrementality contract of
§9 says a shard is either complete-with-a-marker or discardable, and that removes
the judgement call entirely.

This is the first time the discardable-by-construction property has been exercised
for real on this campaign, and it cost four minutes of one network's compute. It is
recorded because session 8 lost an entire survival sweep to the same class of event
under a design that allowed resumption mid-file.

**Consequence to carry forward.** A pool worker that dequeued one of those 64
shards *before* the markers were deleted will have found it complete, skipped it,
and exited cleanly, so the pool's own queue may not revisit it. **The pool must
therefore be re-invoked after it finishes**, which is safe and idempotent by
construction, and the final shard count must be checked against 725 before any
analysis is quoted.

---

# Appendix C — 2026-09-16, the per-shard wall cap was implemented after the sweep started

§5.4 pre-registered a **5-hour per-shard wall cap** with a `censored_wall_cap`
status. When the sweep was launched, the per-*radius* cap (300 s, matching the
committed corpus) was implemented but the per-*shard* cap was not: it existed in
the pre-registration and not in the code. It was implemented and tested roughly
an hour into the run.

**What this means for the data, stated exactly.**

- The cap is a **safety net on wall clock**, not a scientific parameter, and it
  is not a function of any knowledge-dependent or radius-dependent quantity.
- Grid points are swept in ascending order, so a cap that *fires* truncates the
  **high-intensity end** of a curve. That is a bias and is flagged wherever it
  occurs: a censored grid point is written as a cell with
  `status == "censored_wall_cap"`, `n_draws == 0` and `wall_until_timeout_s`
  set, it is excluded from every endpoint, and any unit containing one carries
  `endpoint_censored = True`.
- The eight shards already in flight when the cap landed ran, and will finish,
  **uncapped**: a running worker holds the code it was launched with, so the cap
  cannot apply to them retroactively. Four of them — the `pathfinder` flip
  shards — were still running when this appendix was written, at about 1 h 50 m
  each and roughly half way through their grids, which projects to between
  3 and 6 hours. **One or more of them may therefore exceed the 5-hour cap that
  §5.4 pre-registered**, and if so it will have run to completion rather than
  being censored. The realised elapsed time of **every** shard is recorded in its
  completion marker; the final report states the maximum observed against the
  18,000 s cap and names any shard that exceeded it.
- Every other shard in the campaign ran under the cap as implemented.
- Uncapped is the scientifically cleaner outcome here, since a fired cap would
  have removed data rather than added it: an over-running shard yields a
  *complete* curve, a censored one yields a truncated curve at exactly the
  high-intensity end. The gap is recorded because the pre-registration claimed
  an enforcement that did not exist for the first hour of the run, and that
  claim has to be made true or withdrawn — here it is narrowed to "enforced for
  every shard launched after the first hour", which is what actually happened.

The cap's behaviour was verified before it was relied on, by running a real
shard with `--wall-cap-s 0.5`: 16 of 22 cells came back `censored_wall_cap` with
`n_draws == 0`, `S == None` and `wall_until_timeout_s` set, the censored grid
points were exactly the tail of the ascending grid, and no censored row carried
a key named `seconds`.

---

# Appendix D — 2026-09-16, P5 is falsified, and a post-hoc diagnostic registered before it is run

## D.1 P5 is falsified, decisively and in the direction the pre-registration named as
## the more interesting one

§8 P5 predicted the flip-arm contradiction rate at depth 1 on real structure at
**0.65–0.90**, "directionally **higher** than synthetic", against the synthetic
figure of 0.627 — on the reasoning that real CPDAGs are far more compelled
(median undirected fraction 10.7%), so a reversed claim should more often collide
with a compelled edge.

**Measured, at truthful base knowledge: 0.0127 pooled, 0.0135 network-weighted.**
Forty-three of the 47 `(network, coverage)` draw-sets show a contradiction rate of
**exactly 0.000**; only `mediator` (0.343), `sachs` (0.137) and `ecoli70` (0.115)
are non-zero at all. The prediction is wrong by roughly a factor of fifty, and it
is wrong in the direction §8 pre-committed to treating as the headline: *"If it
comes out lower, that is the more interesting result and it is reported as the
headline of that section."*

The rate does rise with base wrongness — 0.013 / 0.040 / 0.082 / 0.096 at
`b = 0.00 / 0.10 / 0.25 / abs 1` — which is the **opposite** of the synthetic
direction (0.653 / 0.647 / 0.563, falling). It also rises steeply with depth
(0.013 → 0.163 → 0.406 → 0.544 → 0.690 → 0.800 at `d = 1…6`), so the arm is not
degenerate; it simply starts from a floor near zero.

**Why this matters beyond the prediction.** The paper's Result D says most
orientation errors are self-revealing: Meek closure fails and the analyst finds
out about six times in ten. On this corpus it is essentially never. That does not
weaken the paper's motivation — it strengthens it, since the dangerous case is
precisely *consistent and invalidating* — but Result D's **number** does not
transfer to real structure and the paper must not imply that it does.

## D.2 The confound, named before it is resolved

The comparison real-versus-synthetic at depth 1 confounds two things:

1. **structure** — real CPDAGs versus `component_generator`'s small dense
   components; and
2. **the knowledge model** — the real corpus's `K` comes from
   `demo/example.py::knowledge_to_recover`, a **greedy-minimal generator set**,
   while the synthetic arm's comes from `synth/knowledge.py::draw_k_true`, a
   random subset of the CPDAG's undirected edges. Those are different objects. A
   minimal generator asserts one claim per "independent" decision, so reversing
   one of them need conflict with nothing; a full assertion set states many
   orientations that Meek would otherwise derive from each other, so reversing one
   of them contradicts its neighbours immediately.

A difference in contradiction rate could be caused by either. As it stands the
comparison cannot separate them, and **no causal reading of the real-versus-
synthetic gap may be published until it can**.

## D.3 The diagnostic, and its directional prediction, recorded before running it

**Design.** On the *same* real networks, at coverage 1.0, compare two knowledge
models by **exhaustive** single-claim reversal — every claim reversed in turn, no
sampling:

- `K_min = select_knowledge(dag, cpdag, 1.0)` — the corpus's greedy-minimal
  generator, exactly what the main sweep uses;
- `K_full = draw_k_true(dag, cpdag, rng, 1.0)` — every undirected CPDAG edge
  asserted truthfully, which is the synthetic arm's knowledge model at full
  coverage.

Structure is held fixed; only the knowledge model moves. This is a **post-hoc
diagnostic**, opened in response to a falsified prediction, and it is labelled as
one everywhere it is reported. It is **not** one of the nine strata and it
produces no τ.

**Prediction, recorded now, before the diagnostic is run.** `K_full` will show a
**much higher** single-reversal contradiction rate than `K_min` — I predict
**≥ 0.30 pooled for `K_full` against ≈ 0.01 for `K_min`** — because with every
orientation asserted, reversing one puts it in direct conflict with the claims
Meek would derive from its neighbours.

**What each outcome licenses.**

- If `K_full` is high and `K_min` near zero: the real-versus-synthetic gap is
  driven by the **knowledge model**, not by real structure, and the paper must
  say that Result D's figure is a property of *how the knowledge was written
  down*, not of the graphs. This would be the most useful outcome, and it would
  also bear directly on `[RE-1]`, which asks exactly how real elicited knowledge
  compares with a minimal generator.
- If **both** are near zero: real structure genuinely is not self-revealing, and
  Result D does not transfer for structural reasons.
- If `K_full` is *lower* than `K_min`: the mechanism above is wrong and nothing is
  claimed until it is understood.

---

# Appendix E — 2026-09-16, the Appendix D diagnostic, and the mechanism it exposes

## E.1 The D.3 prediction is confirmed, and by a factor of twenty

Appendix D.3 predicted, before running it, that the all-undirected-edges knowledge
model would show a **much higher** single-reversal contradiction rate than the
greedy-minimal generator — *"≥ 0.30 pooled for `K_full` against ≈ 0.01 for
`K_min`"*. **Exhaustive** measurement on the same 25 networks, every claim
reversed in turn, no sampling:

| knowledge model | single reversals | contradictory | rate | networks at exactly 0 |
|---|---|---|---|---|
| `minimal_generator` (coverage 1.0) | 238 | **5** | **0.021** | 21 of 25 |
| `all_undirected_edges` (coverage 1.0) | 410 | **177** | **0.432** | 5 of 25 |

The prediction is confirmed. The internal control is clean: on the five networks
where the minimal generator *is* the full assertion set (`Acid_1996`,
`Sebastiani_2005`, `asia`, `diabetes`, `hepar2`, all with `|K_min| = |K_full|`),
both models give **exactly 0.000** — as they must, since they are the same set.

**So the real-versus-synthetic contradiction gap is driven by the knowledge
model, not by real structure.** Holding the graphs fixed and changing only how
the knowledge is written down moves the rate from 0.021 to 0.432, a factor of
twenty. Real structure with the synthetic knowledge model (0.432) is much closer
to the synthetic figure (0.627) than to the real corpus's own (0.021).

The paper's Result D must therefore be restated as a property of **how the
analyst's knowledge is phrased**, not of the graphs. This bears directly on
`[RE-1]`, which asks exactly how real elicited knowledge compares with a minimal
generator: the answer here is that the difference is not cosmetic, it is the
whole effect.

## E.2 A post-hoc leverage analysis, labelled post-hoc

**This was not pre-registered.** It was opened after seeing that τ_b(`r_hop`,
AUC) is weak in the flip arm, and it is reported as a post-hoc explanation. The
nine pre-registered strata of §3 remain the primary result and are not restated
in its light.

The mechanism. The flip arm corrupts in **claim** units; `r_hop` is a distance in
**hop** units. `docs/PAPER_NARRATIVE.md` §4 already names the exchange rate — the
**leverage** `k_g0 / |K|`, how many data-compelled orientations hang on one
statement — and states that the two radii coincide exactly when `|K| == k_g0`. On
this corpus the leverage runs from **1.0 to 17.0**:

| network | `\|K\|` | `k_g0` | leverage | `r_hop` range | median AUC |
|---|---|---|---|---|---|
| hailfinder | 1 | 17 | 17.0 | 2 | 0.000 |
| paths | 1 | 14 | 14.0 | 1–14 | 0.000 |
| child | 2 | 12 | 6.0 | 1–2 | 0.358 |
| insurance | 3 | 18 | 6.0 | 1–3 | 0.569 |
| … | | | | | |
| diabetes | 26 | 26 | 1.0 | 1 | 0.445 |
| hepar2 | 9 | 9 | 1.0 | 1 | 0.247 |

`paths` is the clearest case and it is worth stating on its own: **an analyst who
made one statement has a hop radius of up to 14 and a survival AUC of exactly
0.000.** The hop radius certifies thirteen clean shells; reversing the single
statement they actually made destroys the adjustment set every time. That is
precisely the over-promise §4 of the narrative warns about, measured on real
structure against an average-case endpoint.

Splitting the flip strata by leverage, with the same cluster bootstrap over
networks:

| stratum | leverage = 1 (`\|K\| == k_g0`) | leverage ≤ 2 | leverage > 2 | all |
|---|---|---|---|---|
| flip cov=1.0 bw=0.00 | undefined (predictor constant) | +0.069 [−0.116, +0.285] | +0.167 [−0.535, +0.579] | +0.064 [−0.351, +0.379] |
| flip cov=1.0 bw=0.10 | **+0.230 [+0.184, +0.316]** | **+0.239 [+0.084, +0.351]** | +0.146 [−0.590, +0.536] | +0.126 [−0.299, +0.412] |
| flip cov=1.0 bw=0.25 | **+0.309 [+0.205, +0.478]** | **+0.402 [+0.267, +0.515]** | +0.296 [−0.533, +0.694] | +0.292 [−0.240, +0.603] |

**Where the analyst's knowledge is Meek-closed — the exact condition under which
the paper says the certificate is sound — the radius ranks survival with an
interval that excludes zero. Where leverage is large it does not, and the
interval is enormous.** The pooled figure sits between the two because it mixes
them.

This also explains why the **tiered** arm is the one where the radius performs
well (+0.378 to +0.580, all three intervals excluding zero): tiered knowledge
asserts every cross-tier edge, so it is far closer to Meek-closed than a greedy
minimal generator is.

## E.3 What is and is not claimed from E.2

- **Claimed:** on the committed real corpus, τ_b(`r_hop`, survival) is materially
  larger and its interval materially tighter in the sub-population where
  `|K| == k_g0` than where `|K| ≪ k_g0`, in the two flip strata where the
  predictor is not constant.
- **Not claimed:** any pre-registered status for that comparison. It is post-hoc.
- **Not claimed:** that leverage *causes* the difference. The leverage-1 networks
  differ from the others in other ways too, and with 5–9 networks per band the
  comparison is not powered to isolate a cause.
- **Not claimed:** that the primary table should be read through this lens. The
  nine strata are the pre-registered result and they are weaker than the
  synthetic ones; that is the headline, and this appendix is the explanation
  offered for it, not a replacement.

---

# Appendix F — 2026-09-16, CORRECTING Appendix E.2: the leverage explanation does not survive a control

## F.1 What E.2 claimed, and why it is withdrawn

Appendix E.2 offered **leverage** `k_g0 / |K|` as the explanation for the weak
flip-arm result, on the strength of a split showing τ_b(`r_hop`, AUC) larger and
tighter where `|K| == k_g0` than where leverage is large. **That explanation is
withdrawn.** It is confounded with `|K|`, and when `|K|` is controlled the pattern
does not merely weaken — it **reverses**.

Restricting to `|K| ≥ 3` and then splitting by leverage, same cluster bootstrap
over networks:

| stratum | `\|K\| ≥ 3`, leverage ≤ 2 | `\|K\| ≥ 3`, leverage > 2 | `\|K\| ≥ 3`, all |
|---|---|---|---|
| flip cov=1.0 bw=0.00 | +0.058 [−0.147, +0.290] | +0.596 [+0.000, +0.673] | +0.286 [−0.021, +0.513] |
| flip cov=1.0 bw=0.10 | +0.235 [+0.075, +0.353] | +0.582 [+0.241, +0.639] | +0.376 [+0.157, +0.528] |
| flip cov=1.0 bw=0.25 | +0.407 [+0.287, +0.511] | +0.674 [−0.344, +0.877] | +0.583 [+0.327, +0.684] |

High leverage now looks **better**, not worse, in all three strata. The
high-leverage bands rest on **3 to 4 networks** apiece, so their intervals are
close to meaningless either way — which is itself reason enough not to have built
an explanation on them. Across the 24 networks, τ_b(`|K|`, leverage) = **−0.395**:
the two are entangled, because a large chain component with a single minimal
generator produces high leverage and small `|K|` at once.

**This is a correction to the orchestrator's own post-hoc analysis, recorded in
order rather than edited into E.2.** E.2's numbers are correct as computed; its
*interpretation* was wrong, and it was wrong in the direction that made the
campaign's story tidier, which is exactly when a claim needs the most scepticism.

## F.2 The explanation that does survive

What actually distinguishes the two bands in E.2 is not leverage but **`|K|`**,
through a design–corpus interaction that has nothing to do with the radius:

> The flip arm's depth grid is **fractional**. At `|K| = 1` the only corruption
> the arm can draw is "reverse everything". So a network whose analyst asserted
> one claim has its survival endpoint pinned at **exactly 0**, by construction,
> whatever its radius.

Five of the 24 networks at coverage 1.0 are pinned at a median AUC of exactly 0 —
`Acid_1996`, `Didelez_2010`, `hailfinder`, `mediator`, `paths` — and **all five
have `|K| ≤ 3`**, while their leverage spans 1.0 to 17.0. Conversely `child`, at
leverage 6.0, is not pinned. Leverage does not sort them; `|K|` does.

Removing them — a **post-hoc sensitivity, labelled as one** — raises τ_b(`r_hop`)
in all three coverage-1.0 flip strata, to **+0.286**, **+0.376** and **+0.583**,
with two of the three intervals excluding zero.

## F.3 What is claimed after the correction

- **Claimed:** on this corpus the flip arm's fractional depth grid interacts with
  a greedy-minimal generator set to pin the endpoint at zero on the small-`|K|`
  networks, and those networks sit at the bottom of the endpoint ranking
  irrespective of their radius. Removing them substantially strengthens the
  radius's rank correlation. This is a property of the **measurement design
  applied to this corpus**, not of the radius and not of real structure.
- **Claimed:** the `paths` spotlight stands unchanged. It does not depend on the
  leverage explanation; it depends only on `|K| = 1`, `k_g0 = 14`, radii spanning
  1–14 and survival exactly 0, all of which are exhaustively enumerated facts.
- **Withdrawn:** E.2's claim that the radius ranks survival *because* the
  analyst's knowledge is Meek-closed. The data do not support it once `|K|` is
  controlled.
- **Unchanged:** the nine pre-registered strata. None of this restates them, and
  the primary result remains the weaker one.

---

# Appendix G — 2026-09-16, the contradiction table was counting the same draw once per pair

## G.1 The defect

`real_analyse.build_contradiction_table` pooled the contradiction rate over
**cells**: it summed `n_contradictory` and `n_draws` across every
`(instance, grid point)` row in a stratum.

On this corpus that double-counts, heavily. The analyst's knowledge is a property
of the `(network, coverage)` cell, so the sweep draws **one** corrupted claim set
per `(shard, grid point)` and scores it against **every** admissible pair of that
network. Whether a corruption contradicts the CPDAG is a property of that claim
set alone — it has nothing to do with which `(X, Y)` query is being asked — so a
network with 95 admissible pairs contributed the same draw 95 times.

The consequence, at flip coverage 1.0, base wrongness 0.00, depth 1:

| | cells pooled | reported `n_draws` | contradiction rate |
|---|---|---|---|
| as computed (pair-weighted) | 445 | **445,000** | **0.0117** |
| correct (draw-set weighted) | 445 | **22,000** | **0.0270** |

The draw total was inflated by a factor of twenty, and the rate itself moved by
more than a factor of two, because networks were weighted by their pair count.

**No τ, no endpoint and no survival number is affected** — those are per-pair
quantities and are correctly per-pair. Only the contradiction series was wrong.

## G.2 Where the defect came from, and what caught it

**The defect originated in the orchestrator's specification, not in the
implementation.** The written brief for the analysis harness said, verbatim,
*"pooled contradiction rate over units (sum of contradictory draws over sum of
draws — NOT a mean of rates)"*. That instruction is correct about the arithmetic
— a pooled rate, not a mean of rates — and wrong about the unit, because on this
corpus the draw is shared across a network's pairs and "over units" therefore
counts it many times. The harness implemented exactly what it was told. Recorded
this way round because the alternative — letting it read as a subagent's error —
would misplace the lesson: the specification is where this had to be caught, and
it was not.

**What caught it** was the orchestrator's own independent recomputation, which
had deduplicated by `(shard_id, grid_point)` from the start precisely because the
draw-set is shared, and which therefore disagreed with the harness. The
disagreement was noticed only because the two were compared; each number was
individually plausible, and a contradiction rate of 0.0117 would have been
reported without anything looking wrong.

## G.3 The fix

The primary `contradiction_rate` is now computed over **distinct draw-sets**, and
the pair-weighted figure is reported beside it as
`contradiction_rate_pair_weighted` rather than discarded, so the two can never be
silently confused again. `n_draw_sets` is printed so the real amount of evidence
is visible.

**Every contradiction-rate number in the report is the draw-set-weighted one**,
and the P5 falsification of Appendix D is unaffected: it was computed by the
orchestrator's deduplicating script and independently confirmed by an exhaustive
enumeration that involves no sampling at all.
