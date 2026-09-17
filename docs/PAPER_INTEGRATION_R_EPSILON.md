# Paper integration brief — Session 8 (the ε-bias radius, and an unconditional `r_val`)

**Audience.** A session sitting down to update `ICLR_2027___bkrobust_v0` with
this session's work. You do not need to read any code. This file tells you what
changed, what to claim, where every number lives, and what must not be claimed.

**The two hard rules of `docs/PAPER_NARRATIVE.md` still bind.** Never invent a
number — §4 below maps every claim to the file it comes from. And read
`src/bkrobust/core/conventions.py` before quoting any radius: `r` is the
distance to the *nearest* failure, so the practitioner-facing safe-move count is
`r − 1`.

---

## 1. The two things that change the paper

### 1.1 A structural claim the paper currently hedges is now a theorem

The paper attaches a caveat to every radius it reports: the upward search is
exact *iff* Conjecture 2 holds, which rests on Anti-Exchange Case B, which
§A.3 states is **open** ("verified with no violation on 36,094 applicable
triples … the classical argument does not transfer"). Section 7 and the
reproducibility statement both repeat it.

**That caveat is now removable.** Case B is proved. The argument the session was
handed had a false premise; the corrected proof is in `THEOREMS.md` §24 and is
summarised in §3.1 below. No previously reported number changes — what changes
is that "exact if the conjecture holds" becomes "exact", for `r_val` and for
`r_ε` alike.

This is the single highest-value edit available. It removes a hedge from the
abstract-level claim of Section 5, and it retires the open problem in §A.3.

### 1.2 The paper's "what does failure cost?" gap is now answered

Appendix C's shell table has a bias column that is a *sampled mean* with an
explicit caveat that per-element maxima are "a sampled proxy for a worst case
and not a supremum". This session replaces the proxy with the supremum, proves
it has the properties a certificate needs, and derives `r_ε` from it.

The practitioner statement the paper can now make: **if at most `k` of your
asserted orientations are wrong, your reported effect is off by at most `ε_y`,
and a bias above `ε_x` is already reachable.** The bound is attained, so it
cannot be tightened without a further assumption.

---

## 2. Files to read, in priority order

Only high-level artefacts are listed. Everything under `src/` and
`experiments/` is implementation and can be ignored for writing.

### Tier 1 — read these

| File | What it is |
|---|---|
| **`report_r_epsilon.md`** | The session report. Findings, scope, and an explicit "what is not established" section. Start here. |
| **`docs/R_EPSILON_THEORY.md`** | Every statement and proof, with a per-item status table (§6) saying exactly what each one assumes. The source for all theorem text. |
| **`paper/sections/epsilon-radius.tex`** | **A drafted section, already in the paper's prose register, with no `TODO` numbers left.** Written to drop in as a new section after Section 5. Contains the theorem statements, the algorithm discussion, the assumption paragraph (§1.1 above) and the full evaluation write-up. |
| **`THEOREMS.md` §§16–24** | The formal record in house format. §24 is the Anti-Exchange gap and its repair; §§16–23 are the new results with a status table at §23. |

### Tier 2 — evidence, when you need to cite or check a number

| File | What it holds |
|---|---|
| `results/epsilon/study/SUMMARY.md` | The 900-instance sweep, analyses 1–6, human-readable. |
| `results/epsilon/study/analysis.json` | The same, machine-readable, with full distributions. |
| `results/epsilon/chickering/SUMMARY.md` | The exhaustive n ≤ 5 verification of the Anti-Exchange repair. Scope table. |
| `results/epsilon/counterexamples/SUMMARY.md` | The two negative results, with minimal witnesses and stated search scope. |
| `results/epsilon/worked/SUMMARY.md` | The running example, shell by shell, for all three knowledge states — the direct replacement for Appendix C's table. |
| `results/epsilon/finite_sample/SUMMARY.md` | Plug-in vs conservative `r_ε` under an estimated covariance. |
| `results/epsilon/verification/manifest.json` | The P1–P7 falsification counts (436,798 checks, 0 violations). |

### Tier 3 — figures, ready to use

| Figure | Claim it carries | Suggested caption |
|---|---|---|
| `figures/fige1_staircase.{png,pdf}` | The headline object. Three panels, one per knowledge state of the running example. | "The bias staircase on the running example. Worst-case bias is exactly zero out to `r_val` (shaded), then steps to the ceiling. The dotted rules are the ε grid, so each radius is readable off the plot. Contrast Appendix C, whose gradually-rising column is a sampled mean." |
| `figures/fige2_what_epsilon_buys.{png,pdf}` | What the refinement actually buys. **The `unreached` bar is the story, not the `0`/`1` bars.** | "Distribution of `r_ε − r_val` over 900 instances. Tolerating ε rarely buys an extra certified shell; what it buys is a ceiling — at ε = 100% the threshold is unreachable on 70.1% of instances, a statement `r_val` cannot make." |
| `figures/fige3_cost.{png,pdf}` | Cost of the two exact strategies. | "Median states evaluated by the two exact strategies against `|K_G0|`. Both return the same radius; incremental dominates because `r_ε` is small almost everywhere." |

`replay.py` sits in `results/epsilon/counterexamples/` and re-derives every
headline counterexample number from its stored witness, if a referee asks.

---

## 3. The findings, with their evidence

### 3.1 Anti-Exchange Case B is proved

**The gap.** The supplied proof claimed Chickering's covered-edge-reversal path
between two extensions `D_A, D_B` of an MPDAG `G` stays inside `[G]` because it
"operates exclusively on uncompelled edges". False: the orientations of `K_G`
are *precisely* the edges uncompelled in `[Ĉ]`, so they are exactly what the
theorem is free to reverse. The classical substitute — reorienting freely inside
a chordal chain component — is unavailable, because an MPDAG carrying background
knowledge can have a non-chordal one. That is the same phenomenon §A.1 of the
paper already documents, and the `K₄ − {V2,V3}` witness was re-confirmed.

**The repair.** Use the *strong* form of Chickering's theorem, which counts its
reversals: between Markov-equivalent graphs differing on `m` edges there is a
sequence of **exactly `m`** distinct covered-edge reversals. A walk of `m` unit
steps taking the difference count from `m` to `0` must decrease at every step,
so every reversal flips a *currently differing* edge and never touches one on
which `D_A` and `D_B` agree. Both extend `G`, so they agree on all of `dir(G)`;
every intermediate therefore retains it and lies in `[G]`. Case B follows as
written. → `THEOREMS.md` §24.2–24.3.

**Verified**, exhaustively over every CPDAG on ≤ 5 nodes, no sampling:

| check | objects | violations |
|---|---|---|
| L1 — covered reversal preserves the class | 50,758 | 0 |
| L2 — an exactly-`m` difference-reducing sequence exists | 206,878 | 0 |
| L3 — every intermediate stays in `[G]` | 2,907,242 | 0 |
| L4 — `[G]` connected under in-class reversals | 118,643 | 0 |
| AE-B, semantic form | 293,814 | 0 |

**3,584,369 checks, 0 violations.** → `results/epsilon/chickering/SUMMARY.md`.

### 3.2 The construction

Fix `Σ`. The analyst reports `θ_Z = β(Z;Σ)`; each candidate truth `D` assigns
`τ_D = β(pa_D(X);Σ)`, the IDA quantity. Then `T(G) = {τ_D : D ∈ [G]}` and
`B(G) = max{|θ_Z − τ| : τ ∈ T(G)}`.

- **A** — `d(G₀,G) < r_val ⟹ B(G) = 0`, exactly, at the population level, for
  every law compatible with `Ĉ`. Not "insignificant": zero.
- **B** — `B` is monotone under the order, so `{B > ε}` is an up-set.
- **U** — the paper's Proposition 1 only ever consumed upward-closure of
  failure, so it is a theorem about *arbitrary* up-sets and `r_ε` inherits it.
  **This is the reason `r_ε` costs no more than `r_val`.**
- **C** — the staircase `β↑` is non-decreasing, zero below `r_val`, attained on
  the sphere, and shell `d` is built directly from the `d`-subsets of `K_{G₀}`.
  One traversal answers every ε.
- **D / D′** — the certificate, and its attainment.
- **E** — `T(G)` is semi-local: `2^{deg(X)}` Meek closures, no DAG enumeration.

`r_ε` **adds no assumption to `r_val`**, which after §3.1 means both are
unconditional. → `docs/R_EPSILON_THEORY.md`, table at §6.

**Falsification harness**: 400 instances × 3 parameter draws, tolerance `1e-9`.
P1 zero-shell 1,758; P2 monotonicity 372,099 ordered pairs; P3 staircase 17,160;
P4 three algorithms vs brute-force BFS 5,844; P5 nesting 5,844; P6 certificate
audit 22,680; P7 semi-local vs enumerated 11,413. **436,798 checks, 0
violations.** → `results/epsilon/verification/manifest.json`.

### 3.3 The sweep — 900 instances

→ `results/epsilon/study/SUMMARY.md`

**What ε buys is a ceiling, not a longer runway.** This is the finding that
should shape how the section is written; the obvious framing ("tolerating bias
extends the safe region") is not what the data says.

| ε | `r_ε > r_val` | ε never reached | extra shells when both finite |
|---|---|---|---|
| 1% | 1.0% | 1.4% | 1.00 |
| 5% | 3.9% | 3.8% | 1.00 |
| 10% | 8.9% | 7.4% | 1.05 |
| 25% | 24.5% | 20.2% | 1.09 |
| 50% | 48.5% | 40.1% | 1.14 |
| 100% | 75.6% | **70.1%** | 1.12 |

`r_val` itself: n = 893, mean 1.09, median 1, max 3; unreachable on 0.8%.

Other results from the same sweep:

- **`r_0 = r_val` on all 900 instances, 0 exceptions.** The converse to Theorem A
  fails only non-generically and the non-generic case was never hit.
- **Staircase shape: step 757 (84.1%), ramp 136 (15.1%), flat 7 (0.8%).** The
  gradual rise in Appendix C's column is, in most instances, an artefact of
  averaging over a shell whose contaminated fraction is growing — not of the
  worst case rising.
- **Certificate: 0 violations over the 884 instances where `β↑(k)` was actually
  computed.** 16 further instances had the traversal stop before depth `k`, so
  the recorded value is a lower bound and the comparison is not a test of the
  theorem; recomputing their full staircases, the certificate holds there too.
  Among non-degenerate bounds the conservativeness ratio has **median 1.000 and
  maximum 1.000** — the bound is routinely *attained*. (This is the bias
  analogue of the paper's 0.43–0.47 coverage-conservativeness figure, and it
  points the opposite way: this bound is tight, not slack.)
- **Cost**: both exact strategies agree on all 900 (0 disagreements).
  Incremental median 1 state / mean 2.09 / max 65; bisection median 1 / mean
  5.14 / max 172. By `|K_{G₀}|` the gap opens up: at 9, mean 24.3 vs 129.3. The
  anytime greedy chain, which enumerates no shell, returns the **exact** radius
  on **98.4%** of instances.
- **Against the incumbent** (matched subsample, n = 50): disagrees on ≥ 1
  threshold in **84%** of instances, both directions (23 larger, 14 smaller, 8
  equal, 5 mixed). On the same instances the incumbent's mean was non-monotone
  under the order in **49 of 50**.

### 3.4 Two negative results

→ `results/epsilon/counterexamples/SUMMARY.md`

**C1 — bias is not monotone in distance.** Minimal at four nodes: a state at
distance 1 with `B = 1.3405`, a state at distance 2 with `B = 0.0000`. They are
incomparable in the order, so monotonicity-under-⪯ is untouched; what fails is
the reading that confuses distance with knowledge. **This is why `β↑` maximises
over a shell rather than following a perturbation path**, and it is worth a
sentence in the paper because the naive reading is the tempting one.

**C2 — the mean-based alternative is not merely unjustified, it is wrong.**
Flagship witness: `[G] ⊊ [H]` with mean bias 0.96–1.03 at `G` against 0.43–0.46
at `H`, reproduced 8/8 at a minimum of **11.4 combined standard errors**. Rate
on identical pairs, Monte-Carlo-aware (3 s.e. in each of two repetitions, 300
draws, so a **lower bound**): **34.7%** for the mean, **0.0%** for `B`.
Consequence: the retraction-only search disagrees with brute force on **269 of
644** threshold probes, in the recorded cases by reporting that nothing reaches
the threshold anywhere when something does at distance 2 — an error in the
direction that overstates robustness.

### 3.5 Finite samples

→ `results/epsilon/finite_sample/SUMMARY.md`

`r_val` reads no numbers and has no sampling distribution. `r_ε` reads a
covariance and does. Over 120 instances, n = 100 … 10,000: the **plug-in**
radius is too large — the direction that overstates robustness — in up to
**5.0%** of cases at n = 100, falling to 0–0.8% at n = 10,000. Thresholding a
bootstrap upper bound and reporting the 5% quantile of the radius gives
**0.0% in every cell**, at the cost of understating on 1–17%. The conservative
radius is the one to report.

This connects directly to the evaluation protocol's stated input
`(Ĉ, K, Σ̂, n)` — the `Σ̂` and `n` in that tuple only matter for `r_ε`, and now
there is something to say about them.

---

## 4. Concrete edits to the existing paper

| Where | Edit |
|---|---|
| **Abstract** | The radius currently "prices structural risk". It can now price it *in effect units*: add the band statement. |
| **§5.1, after Prop. 1** | Note that the proposition's proof uses only upward-closure, so it holds for any up-set — this is what licenses reusing the search for `r_ε`. |
| **§5.3 / Table 3** | Unchanged. `r_ε`'s cost discussion is in the new section. |
| **§7 "Assumptions inherited by every number"** | Delete the one-sided-error sentence about `??`. It is now a theorem. |
| **§A.3 "Status of ??"** | Replace wholesale with the proof from `THEOREMS.md` §24.2–24.3 plus the verification table of §3.1 above. The section currently says Case B "is open". |
| **App. C, Table 2** | Keep, but re-caption: the bias column is a sampled mean, and the worst case behaves differently (step, not ramp). Put `results/epsilon/worked/SUMMARY.md`'s table beside it or in place of it. |
| **Reproducibility statement** | "with the one property that is verified rather than proved isolated in Section A.3" — no longer true; update. |
| **New section, after §5** | Drop in `paper/sections/epsilon-radius.tex`. |
| **References** | Chickering (1995/2002) is already cited in spirit; add the explicit entry, since the repair now leans on the *strong* form of the theorem. Maathuis–Kalisch–Bühlmann (IDA) needs adding — `τ_D` is their quantity. |

---

## 5. What must NOT be claimed

- **Verification is not proof.** §3.1's sweep corroborates the repair and does
  not extend it; the proof in `THEOREMS.md` §24 is what carries the claim. The
  sweep is exhaustive only through n = 5 (n = 6 is an L1 spot check only).
- **The 34.7% non-monotonicity rate is a lower bound** by construction. Do not
  round it up or describe it as "about a third of pairs" without the word
  *at least*. The naive rate on the same data was 43.8%; the gap is
  Monte-Carlo noise that the conservative rule deliberately discards.
- **Linear-Gaussian is where the *cost* claims hold, not the theorems.**
  Theorems A–E need only that the estimand is a functional of `P` and `Z` and
  that each DAG identifies a number from `P`; the g-formula under positivity
  satisfies both. Linearity only makes `B` a closed-form evaluation. **Nothing
  beyond linear-Gaussian was measured**, so do not claim empirical support there.
- **Scope of the sweep**: n ≤ 12, `|undirected| ≤ 10`, synthetic generators
  only. **The real-graph corpus was not touched** — `r_ε` on the 85-vertex
  component is untested, and §6.2's scalability claim does not transfer to it.
- **`Z` is the optimal adjustment set throughout.** The robust-selection
  analysis §C runs for `r_val` across all valid sets; nothing equivalent was
  done for `r_ε`.
- **The 84% step / 15% ramp split is a property of this corpus, not a theorem.**
  Nothing forbids a long ramp.
- **`r_ε` is not a replacement for `r_val`** and must not be presented as one.
  `r_val ≤ r_ε` always, `r_val` is parameter-free while `r_ε` reads data, and on
  most instances they coincide. The refinement adds a magnitude, not a bigger
  number.

---

## 6. One open decision, deliberately not taken

`src/bkrobust/hybrid.py` still carries the old caveat on every result object
(`HybridResult.assumes` = "Conjecture 2 (hence Anti-Exchange Case B, verified
not proved)") and in its module docstring. After §3.1 that string is wrong, in
the conservative direction.

Updating it changes **no computed number**, but it changes the caveat attached
to every radius this repository has produced, and §A.3, §7 and the
reproducibility statement have to move with it. It was left unchanged so that
the paper and the code change together rather than drifting. Decide this before
editing §7.
