# Pareto Phase 2 continuation — P6' and P7' (Appendix A.3)

Run: 2026-09-09. Code: `src/bkrobust/robustness/pareto2.py`,
`src/bkrobust/robustness/run_pareto2.py`. Reads (read-only) round 1's
`src/bkrobust/robustness/pareto.py` for base-CPDAG selection, the truthful
proposal enumerator, and SEM handling. Outputs, all `pareto2_`-prefixed under
this directory: `pareto2_proposals.csv` (1431 rows), `pareto2_menus.csv` (243
rows), `pareto2_aggressive_census.csv` (39 rows), `pareto2_manifest.json`.

Command: `PYTHONPATH=src python3 -m bkrobust.robustness.run_pareto2 --mode full`.
Total wall time: **97.6 s**. `git_sha=d3cd71f`, `git_dirty=True` (recorded in
the manifest; other workers' files were dirty at run time, mine are additive).

## Headline result, stated first: P6' is falsified too, and Phase 2 is abandoned

Per the pre-registered trigger in Appendix A.3 ("If achievable utility is also
invariant, the Phase 2 claim is abandoned entirely rather than redefined a
third time"):

- **Task A (synthetic, 24 strata, same grid as round 1): strata with >1
  distinct `achievable_utility_median`: 0/24.**
- **Task B follow-up (8 real-network `(network, x, y)` bases that survived
  `fast_gate`): bases with >1 distinct `achievable_utility_median`: 0/8.**

Not "approximately zero" — exactly zero, in every one of 32 measured strata,
across two structurally different instance families (synthetic components and
real benchmark networks). **P6' is falsified. Phase 2 (the trade-off between
`r_val` and adjustment-set efficiency) is dropped entirely, per the
pre-registered abandonment clause. No further redefinition of "utility" is
attempted.**

### Why this happened (mechanistic, not just "invariant again")

Inspecting `n_valid_menu_members` and `winning_member_index` alongside
`achievable_utility_median` shows something stronger than round 1's finding.
Round 1 found the *optimal* set was proposal-invariant. Here, the *number* of
GAC-valid menu members among "ok"-status proposals is itself exactly constant
within every one of the 24 synthetic strata (e.g. every "ok" proposal in
`c06_s04_cov100_seed566415042` has exactly 8 valid menu members out of the
8-member menu; every "ok" proposal in `c06_s01_cov100_seed152819226` has
exactly 2), even though proposal size and which specific truthful edges are
asserted vary substantially within that stratum. And `winning_member_index`
(the argmin-variance member) is the same single index every time.

The proposed mechanism, offered as a plausible explanation and not proved:
`is_gac_valid_mpdag` on this small, fixed, true-DAG-derived menu appears to be
governed by a small number of "gating" undirected edges near `(x, y)` rather
than by the full assortment of asserted edges — once enough of `Ĉ` is
truthfully oriented to admit *any* menu member, it already suffices to admit
the same full subset. Combined with the Henckel/Perkovic/Maathuis result that
the population-optimal adjustment set is uniformly at least as efficient as
every other valid set (so `O_true`, when GAC-valid on `G0p`, always wins the
arg-min-variance competition against the rest of the menu — confirmed
empirically: `winning_member_source == "O_true"` for 653/884 "ok" synthetic
rows, with the fallback fixed within a stratum whenever `O_true` itself isn't
GAC-valid there), the achievable-efficiency floor turns out to be an
all-or-nothing property of the stratum, not a continuum the analyst's
specific proposal moves along.

## Acceptance criteria evidence

**C1 — smoke, timed, not estimated.** `run_pareto2.py --mode smoke` on one
synthetic stratum: 0.08 s (24-base projection: 2.0 s). Loading and building
CPDAGs for all 39 real networks: 63.8 s in the smoke run (dominated by the
four `munin*` networks and `link`, matching `results/axisa3/descriptive_structure.csv`'s
own `seconds_to_cpdag`). The **actual full run** (24 synthetic strata + 39-network
census + 8-base real-network follow-up) completed in **97.6 s wall time**,
far under the 45-minute budget — no cuts to base CPDAGs, networks, SEM draws,
or menu size were needed or made.

**C2 — menu sanity.** Every menu member found GAC-valid on a proposal's `G0`
is bias-checked against every SEM draw for that stratum (1046 (member × SEM)
groups checked across all 1431 proposal rows that reached a validity
decision). Global max `|bias|` = **1.51e-14**, comfortably under 1e-8, for
both the synthetic and real-network rows. No STOP fired. Reasoning: every
proposal here (Task A and the Task B follow-up alike) is built by round 1's
`sample_proposals` from `truthful_orientation_menu`, so it only ever asserts
a subset of the *true* orientations; the true DAG is therefore always a
member of `G0`'s extension set, and GAC-validity (checked over *every*
extension) on `G0` implies validity against the true DAG in particular, hence
zero asymptotic bias.

**C3 — determinism.** `run_pareto2.py --mode single --base-index 2` under
`PYTHONHASHSEED=0` and `PYTHONHASHSEED=12345`: identical
`sha256=d4b665c190dd269e5f13063b44d9c3de6702b471890fffdb8be8a3fa4a4c7125`
for both.

**C4 — purity.**
```
grep -n "all_valid_adjustment_sets_mpdag\|synth.runner\|np.random.seed\|random.seed\|import random" \
  src/bkrobust/robustness/pareto2.py src/bkrobust/robustness/run_pareto2.py
```
returns nothing (exit code 1, no matches). The module docstring was worded to
avoid these literal strings (paraphrased, per round 1's own note about this
exact false-positive risk).

**C5 — exclusions, counted out loud.**

Task A (synthetic, 1087 proposal rows, identical count to round 1's run since
the base grid and proposal enumerator are reused unmodified):
- `ok`: 884 (81.3%)
- `no_menu_member_valid`: 203 (18.7%) — matches round 1's 203/1087
  `optimal_set_undefined` count almost exactly, as expected since both are
  driven by the same underlying non-identification proposals.
- `avar_undefined`: 0
- `proposal_contradictory`: 0
- `r_val == UNREACHED (-1)` among `ok` rows: 0
- `radius_timed_out`: 0

Task B follow-up (real networks, 344 proposal rows across 8 gated bases):
- `ok`: 162 (47.1%)
- `no_menu_member_valid`: 182 (52.9%)
- `r_val == UNREACHED (-1)` among `ok` rows: 0
- `radius_timed_out`: 0

Task B network loading: 40 files found, **1 excluded** (`M-bias.txt`,
"excluded_admg_not_dag_by_design", per the task brief), **39 loaded**, 0
further parse/cycle/CPDAG failures.

Task B follow-up gating: 6 networks selected (smallest with >0 on-path
pairs: `mediator`, `Didelez_2010`, `Schipf_2010`, `asia`, `sachs`,
`Kampen_2014`), up to 3 on-path pairs each considered (18 candidate pairs),
**10 gate-excluded** (9 `empty_set_trivially_valid`, 1
`no_atomic_perturbation_changes_validity` — `sachs Erk->Akt`), **8
gate-admitted** and measured: `mediator{X,I,Y}` (all 3 pairs),
`sachs Mek->Akt`, `sachs Mek->Erk`, `Kampen_2014 AFF->{ALN,APA,CDR}`.
`Didelez_2010` and `Schipf_2010` had zero pairs survive the gate at all.

**C6 — census honesty.** `pareto2_aggressive_census.csv` reports all 39
networks, including the 6 with `n_undirected_edges = 0` (`cancer`,
`confounding`, `earthquake`, `mildew`, `pigs`, `survey`) at `n_pairs_sampled =
0`/`n_pairs_with_on_path_undirected_edge = 0` and the four `munin*` networks
plus `win95pts` where sampling caught few or no descendant pairs at all
(e.g. `munin`: 1 candidate pair had a causal path, out of 200 drawn from
~1.08M possible ordered pairs on 1041 nodes — small-sample noise, reported as
is, not resampled to manufacture a nonzero row). No row was dropped or
filtered.

## Task B / P7': structurally available, unlike round 1 — but untestable on the efficiency axis

Round 1 found **zero** aggressive proposals on the synthetic generator by
construction. On the real benchmark corpus: **396 on-path pairs** across the
39-network census (sum of `n_pairs_with_on_path_undirected_edge`), present in
23 of 39 networks. Notably `paths` (96/98 sampled pairs on-path), `sachs`
(22/22), `Polzer_2012` (64/77), `Kampen_2014` (33/47), `mediator` (6/6) — real
networks routinely have undirected CPDAG edges sitting on treatment-outcome
causal paths, confirming the round-1 diagnosis that this was a property of
the synthetic generator's component-construction rule, not a general fact
about CPDAGs.

**P7' is therefore structurally available**, and the follow-up stage produced
both aggressive (215 rows) and conservative (129 rows) labeled proposals
across the 8 gated real-network bases. But since `achievable_utility_median`
is invariant within every one of those 8 bases too (see headline result),
**the intended contrast — "aggressive concentrates in the high-utility /
low-`r_val` region" — cannot be tested on the utility axis**: utility never
moves, aggressive or not.

On the surviving axis (`r_val` alone, which does vary: 2 distinct values in
7/8 bases), the aggressive/conservative contrast is **mixed, not
directional**: aggressive proposals had a lower mean `r_val` than
conservative ones in 2 bases (`Kampen_2014 AFF->ALN`: 1.53 vs 1.78;
`AFF->APA`: 1.70 vs 2.0), an equal mean in 2 bases (`Kampen_2014 AFF->CDR`:
1.5 vs 1.5; `mediator I->Y`: 1.0 vs 1.0), a *higher* mean in 1 base (`sachs
Mek->Akt`: 1.70 vs 1.33 — the wrong direction for the P7 prediction), and no
conservative proposals survived gating at all in 2 bases (`mediator X->I`,
`mediator X->Y`) to compare against. This is not reported as support for or
against P7's directional claim — 8 bases with this much internal noise is not
enough signal either way, and the claim was already moot once the utility
axis (the half of the prediction that mattered) was found invariant.

## Design choices worth flagging (not tuning, documented as they were made)

- **Menu-member `O_true ± v` filtering order** (Appendix A.3 step 5/6): the
  spec's "`v` in `sorted(...)[:4]`, `v ∉ O_true ∪ {x,y}`" is read as *filter
  first, then take the first 4 of the filtered list* (not "slice first, then
  filter, possibly yielding fewer than 4"), so the menu actually gets up to 4
  new sets per rule rather than being at the mercy of alphabetic coincidence.
  Implemented in `pareto2.build_menu`.
- **Representative SEM draw.** Menu-validity (`n_valid_menu_members`) is
  structural — independent of any SEM draw, since `is_gac_valid_mpdag` never
  looks at SEM parameters. Only *which* valid member has lowest asymptotic
  variance depends on the draw's random weights. `achievable_utility_median`
  is computed as specified (min-over-valid per draw, median over the 20
  draws), but `winning_member_index` and the committed `Z` fed to
  `breakdown_radius` for `r_val` are taken from the single draw whose utility
  is closest to that reported median (ties broken by smallest `sem_k`) — one
  fixed, auditable, reproducible committed set per proposal, rather than 20
  different `r_val`s.
- **Census pair sampling.** "Up to 200 ordered pairs" is read as: draw up to
  200 *candidate* ordered pairs (all of them, if the node-pair space is
  ≤200), then `n_pairs_sampled` in the CSV is the count of those candidates
  that satisfy the causal-path precondition (`y` a descendant of `x`) — the
  only ones a background-knowledge claim on an edge could ever be
  "aggressive" against. `n_candidate_pairs_drawn` (the pre-filter count) is
  also recorded for transparency.
- **Real-network `BaseCpdag.component_size` / `.separation`.** Set to `-1`
  (not a synthetic-generator concept) rather than repurposing the field
  incorrectly; `origin="real_network"` and `network=<name>` columns
  distinguish these rows unambiguously in `pareto2_proposals.csv` and
  `pareto2_menus.csv`.
- **Follow-up scope.** 6 networks (the maximum of the 3-6 range), 3 pairs per
  network considered, chosen once by ascending `n_nodes` before any gating or
  measurement was run — not re-picked after seeing which pairs survived the
  gate.
