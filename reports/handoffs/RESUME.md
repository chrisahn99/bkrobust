# RESUME — [RE-11]: survival and the paired cross-arm on real graph structure

**Branch:** `experiments/todo_1` · **Campaign owner:** an Opus 5 orchestrator session
**Started:** 2026-09-16 · **Status: COMPLETE — all ten phases run, all deliverables committed**

Read this file first. It is written for someone who has never seen this campaign.
It is updated **before** anything else whenever state changes.

---

## 1. What this campaign is

The ICLR 2027 draft carries a TODO in its results section:

> *The ranking evaluation of §6.1 has not yet been run on these 831 real
> instances; survival and the paired cross-arm design exist only on synthetic
> structure. This is the experiment that turns the ranking claim into a claim
> about practice, and Table 1 must not be read as covering it.*

That is **[RE-11]** in `docs/REMAINING_EXPERIMENTS.md`. This campaign runs it.

**The claim under test.** The breakdown radius `r_hop` (legacy name `r_val`) is a
**worst-case** distance: the number of elementary knowledge revisions from the
analyst's graph `G₀` to the *nearest* knowledge state at which the committed
adjustment set `Z*` stops being valid. Average-case **survival** — the probability
that `Z*` is still valid after a randomly corrupted knowledge state at some
corruption intensity — is a different object. The paper's Result B says the first
**ranks** the second, and that the baselines `|K|` (stated knowledge size) and
`SHD(G₀, truth)` do not do so consistently.

That is currently established on **synthetic structure only** (1,978 generated
instances, `table_tau_comparisons.md`). Real networks enter the paper only through
a structural census and a radius distribution (`report_real_graphs.md`). This
campaign establishes, or fails to establish, the same ranking claim on the
**committed real-network corpus**: the 831 admissible rows on 25 networks in
`results/axisa3/instances.jsonl`.

**What would falsify it.** Stated before any survival curve exists, and binding:

- **T1.** If τ_b(`r_hop`, survival AUC) is not positive in at least 5 of the 9
  pre-registered real strata, the real-structure ranking claim is **not
  established** and the report leads with a null. No re-stratification, no
  endpoint substitution, no dropping of an inconvenient stratum.
- **T2.** If τ_b(`r_hop`) fails to exceed τ_b(`n_k`) in the majority of strata
  where both are defined, the report says the radius does not outrank the `|K|`
  baseline on real structure.
- **T3.** If the raw and the `n_eval ≥ 30` endpoints disagree in sign or in
  CI-excludes-zero verdict in more than 2 of the 9 strata *at N = 1000*, those
  strata are reported **unresolved**, not resolved by choosing an endpoint
  (`PREREGISTRATION.md` Appendix J).
- **T4.** If tiered survival is not exactly `1.000` at corruption rate 0 for
  every instance individually, the tiered arm halts and the defect is reported
  before any τ is computed.
- **T5.** If the cross-arm pairing yields fewer than 3 intensity bins with ≥ 15
  matched instances in *both* arms, the pairing is declared **not workable on
  real structure** and reported as such. The retracted unpaired comparison is
  **not** substituted.
- **T6.** If the decomposition identity
  `S_contra_as_fail(d) = (1 − contradiction_rate(d)) · S(d)` fails beyond 1e-12
  on any row, the sweep halts.
- **T7.** If the base-wrongness-zero radii do not reproduce
  `results/axisa3/instances.jsonl` **exactly** on every row, the frame is wrong
  and the campaign halts.

A null is a real result and is reported at the same weight as a positive.

**The assumption every radius carries:** *Conjecture 2 (hence Anti-Exchange
Case B, verified not proved)*. The error is one-sided — a radius can be too
large, never too small (`THEOREMS.md` §4c, §6, §8).

---

## 2. Current state

**THE CAMPAIGN IS COMPLETE.** 725 of 725 shards, every derived table rebuilt from
them, the report and the anchor table written, and the verifier passing.

| | |
|---|---|
| **Shards** | **725 / 725**, none failed, none censored, no shard over the 18,000 s cap (max 3.01 h, `pathfinder`) |
| **Draws** | N = 1000 per grid point everywhere; 409 of 450 scored flip shards are *exhaustive* |
| **Deliverables** | `table_tau_real.md`, `report_real_survival.md`, `docs/PAPER_NOTE_RE11.md`, six figures in `figures/s9_*`, `src/bkrobust/analysis/session9_verify.py` |
| **Pre-registration** | `results/axis_robustness_real/PREREGISTRATION.md`, §§1–10 plus **Appendices A–J** |
| **Docs** | `docs/REMAINING_EXPERIMENTS.md`: [RE-11] marked done; [RE-12] explicitly **not** discharged; [RE-1] and [RE-16] raised |

### To reproduce or re-check anything

```
make re11-verify     # re-assert every number in the report
make re11-all        # analysis + figures + verify, from the committed shards
make re11-sweep      # resume the sweep itself; idempotent, skips completed shards
```

The sweep is idempotent: a shard with a completion marker is skipped, one without
is redone from scratch. Nothing ever needs cleaning up after an interruption.

### The result, in one line each

- **Pre-registered:** τ_b(`r_hop`, survival) positive in **9 of 9** strata, interval
  excludes zero in **3 of 9** — against 8 of 9 synthetically.
- **After a stronger leave-one-network-out check (Appendix H):** **1 of 9**
  survives dropping any single network.
- **Post-hoc (Appendix J):** **within** a fixed network and knowledge state the
  radius ranks survival with an interval excluding zero in **6 of 9**, mean τ
  **+0.30 to +0.84**. The pooled figure mixes that with a confounded
  between-network comparison.
- **Result D does not transfer**, and the cause is the knowledge model, not the
  graphs: **5 of 298** single reversals contradict under a minimal generator
  against **177 of 410** under a full assertion set, on the same networks.
- **`paths`**: `|K| = 1`, `k_g0 = 14`, `SHD = 0`, radii 1–14, survival **exactly
  0.000** — the two-radii over-promise, by exhaustive enumeration.

---

## 3. Facts established in Phase 0 (all verified against the repo, not assumed)

1. **The interpreter is `/usr/bin/python3`, Python 3.9.6**, with numpy 2.0.2,
   scipy 1.13.1, pandas 2.3.3, matplotlib 3.9.4 and a working `ortools`
   (`from ortools.sat.python import cp_model` succeeds). This is *exactly* the
   environment recorded in `results/axisa3/manifest.json`, which produced the
   committed corpus. `docs/REMAINING_EXPERIMENTS.md` recommends standing up
   Python 3.12; that recommendation is **declined, deliberately**, because moving
   off the interpreter that produced the committed corpus would trade
   reproducibility for nothing — the live tree runs on 3.9 and only the scaffold
   stubs (which this campaign does not touch) need 3.11+. Run everything with
   `PYTHONPATH=src /usr/bin/python3`.
2. **The committed corpus reproduces.** Re-running `benchmarks.measure.evaluate`
   on `asia`'s three admissible rows returns the committed radius, method,
   `k_g0` and separation exactly. This is the basis of falsification trigger T7.
3. **The corpus is 831 admissible rows** = 543 at coverage 1.0 + 182 at 0.5 +
   106 at 0.25, over 25 networks and **543 distinct `(network, X, Y)` pairs**.
   Only **105** of those triples are admissible at all three coverages, on 8
   networks — so any cross-coverage comparison must be conditioned on that
   matched set ([RE-6]) and is supplementary, never a headline.
4. **368 of 831 rows have undefined separation** (`no_z_member_in_component`).
   A status, never a number, and a stratum reported on its own.
5. **`SHD(G₀, truth) ≡ 0` on every admissible row at coverage 1.0**, and at
   partial coverage it equals exactly the number of undirected edges `G₀` still
   carries — i.e. it measures *incompleteness*, never *wrongness*, because
   `select_knowledge` draws the analyst's claims from the true orientations. The
   SHD baseline is therefore degenerate by construction unless a **base
   wrongness** is applied before the sweep, exactly as the synthetic design did
   and as `results/axisa3/wrong_knowledge.jsonl` does. This is design
   requirement §3.2 of the brief and it is honoured.
6. **`|K|` (stated size) and `k_g0` (commitment size) are constants within a
   `(network, coverage)` cell.** `select_knowledge` returns a network-level
   claim set, and `flip` preserves list length, so no base wrongness can make
   `|K|` vary within a network. On this corpus the `|K|` baseline can therefore
   only rank *networks*, not queries. That is a finding about the corpus, not a
   defect to be engineered away, and it is reported as one.
7. **A defect in the committed corpus, found in Phase 0.**
   `benchmarks.measure.evaluate` sets `g0_undirected_edges` only on the
   `o_g0_extensions_intractable` rejection path. On the admissible path it is
   never assigned, so **all 831 admissible rows carry the dataclass default
   `0`** — which is *not* a measurement. Recomputing it directly gives non-zero
   values at coverage 0.5 and 0.25 (e.g. `insurance` at 0.5 leaves 6 undirected
   edges). This is exactly [RE-5b]. The frame recomputes the column; the
   committed file is **not** modified.
8. **Within-network radius variation is thin.** At coverage 1.0, 10 of 25
   networks have a single distinct radius across all their rows (`diabetes`: 89
   rows, all `r = 1`), while `paths` alone spans `r = 1…14` in 20 rows. A
   network-clustered analysis is therefore dominated by between-network
   variation and `paths` carries large leverage. A **leave-one-network-out**
   sensitivity is pre-registered for every headline τ.
9. **Cost is dominated by Meek closure on two networks.** Measured per corrupted
   draw: `pathfinder` (109 nodes, |K| = 79) ≈ 985–1,600 ms; `diabetes` (413
   nodes, |K| = 26) ≈ 42 ms closure + 9.4 ms per row for the GAC check; every
   other network is under 6 ms. Projected total at N = 1000 with closures shared
   across the pairs of a network: ≈ 37 CPU-h flip + ≈ 14 CPU-h tiered + ≈ 28
   CPU-h cross-arm. The machine has 10 cores. **N = 1000 everywhere is
   affordable** and is what will be run.

---

## 4. How the campaign is built to survive being killed

- Every sweep writes **append-only JSONL**, one row per unit of work, under
  `results/axis_robustness_real/`.
- **Work is sharded** at `(arm, network, cell, replicate)` granularity, sized so
  no shard exceeds ≈ 3 h.
- A shard writes `_done/<shard_id>.json` **only when it is fully finished**,
  carrying its row count and the SHA-256 of its JSONL.
- **Partial shard output is cleanly discardable, by construction.** A shard owns
  exactly one file. On re-invocation, a shard with no completion marker has its
  file truncated and is redone from scratch. Nothing is ever resumed mid-file
  and nothing is ever rewritten in place — this is the defect that destroyed
  session 8's first survival sweep.
- Re-invoking any stage skips completed shards by reading the markers, and is
  safe to run repeatedly.
- An interrupt costs at most one shard per running worker.

---

## 5. Where things live

| | |
|---|---|
| Raw and derived results | `results/axis_robustness_real/` |
| Pre-registration | `results/axis_robustness_real/PREREGISTRATION.md` |
| Anchor table (deliverable) | `table_tau_real.md` (repo root) |
| Report (deliverable) | `report_real_survival.md` (repo root) |
| Figures | `figures/s9_*.{pdf,png}` |
| Verifier | `src/bkrobust/analysis/session9_verify.py` |
| New sweep code | `src/bkrobust/robustness/real_*.py` (added alongside; **no existing module that produced a committed result is modified** — precedent `748f2e1`) |

---

## 6. Log

- **2026-09-16** — Phase 0 complete. Repository read; environment settled on
  `/usr/bin/python3` 3.9.6; committed `asia` rows reproduced exactly; corpus
  profiled; cost budget probed; the `g0_undirected_edges` defect found.
- **2026-09-16** — Phase 1. Pre-registration committed with nine strata, seven
  falsification triggers and directional predictions for all seven hypotheses.
- **2026-09-16** — Phase 2. Frame frozen and hashed. **T7 passed**: all 831 radii
  reproduce `results/axisa3/instances.jsonl` exactly. Appendix A added: the
  rate-based base wrongness is inert on most of the corpus, so an absolute
  one-claim level was added as a supplementary stratum.
- **2026-09-16** — Phase 3. Smoke run at N = 1000 on three cheap networks; every
  pre-registered check passed, including the polynomial GAC predicate agreeing
  with the enumeration oracle on 254 of 254 comparisons.
- **2026-09-16** — Phases 4-6. Sweep launched, 8 workers. Appendix B records two
  defects: a measurement limit mislabelled as a structural property, and an
  orchestrator error that ran shards outside the pool (64 discarded and
  requeued). Appendix C records that the per-shard wall cap was implemented an
  hour into the run.
- **2026-09-16** — **P5 falsified.** The flip-arm contradiction rate at depth 1
  on real structure is **0.0127**, not the predicted 0.65-0.90. Confirmed
  exhaustively: **5 of 298** single-claim reversals across the whole corpus are
  rejected by Meek closure, and 46 of 50 cells are exactly zero. The paper's
  Result D number does not transfer to real structure. Appendix D records the
  falsification and registers a post-hoc knowledge-model diagnostic, with its
  directional prediction, before running it.
- **2026-09-16** — **T5 preview, computed independently of the analysis harness**
  (`ORCHESTRATOR_T5_CHECK.txt`). The cross-arm pairing is workable but only just:
  exactly **3** intensity bins carry >= 15 matched instances in both arms, and
  they are 0.00, 0.05 and 0.10 — a far narrower band than the synthetic
  0.15-0.65, because real `|dir(G0)|` is large and compresses the normalised
  intensity. Those bins rest on 8, 8 and 3 networks. The instance-identity hash
  agrees across the two arms on all 501 instances, so the pairing itself is
  sound.
- **2026-09-16** — Scored-unit denominators measured. Flip strata are healthy
  (540/540, 1579/1629, 1386/1629 at coverage 1.0). **Tiered strata are thin**:
  60, 176 and 265 scored units on 11, 18 and 19 networks at `n_tiers` 2, 3, 4,
  the rest lost to `n_k_zero` and `optimal_set_undefined` — the same structural
  property of the tiered generator that the synthetic Appendix I counted aloud.
