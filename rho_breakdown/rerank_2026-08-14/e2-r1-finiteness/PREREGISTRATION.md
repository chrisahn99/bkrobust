# E2 — R1-closed finiteness: PRE-REGISTRATION

**Written 2026-08-14, BEFORE any code was run against the tier-perturbation ball.**
Nothing below was chosen after seeing a result. The only numbers quoted here are the
pilot's already-published ones (2026-07-17), used to set the direction of the test.

---

## 0. What is being tested

The ρ-breakdown paper's *stated mathematical core* is an **R1-closed finiteness theorem**:
tiered background knowledge (Bang & Didelez, UAI 2023, arXiv:2306.01638) needs only Meek's
**R1**, and R1-closure is therefore the lever that **controls the Meek cascade**, making a
breakdown radius finite/tractable on the tiered sublattice.

The pilot (`RESULTS.md` §5.2) states this is **UNTESTED**, because its tiered arm perturbed
a tiered `K` by *flipping statements*, which destroys the tiering and frees R2–R4. Its
proper test — **perturb the tier assignment, keeping `K'` a tiering** — has never been run.

This experiment runs the empirical half **before** anyone attempts the proof.

## 1. Design

Ensemble: the pilot's own 600-SCM primary arm, regenerated bitwise from the stored
`(seed, p, deg)` triples (`p ∈ {5..8}`, expected degree ∈ {1.5, 2.0, 2.5}, linear-Gaussian
iSCM, single X, single Y with a directed path, causal sufficiency, faithfulness).
Verification requirement: regenerated `τ` must match `results/linear_raw.json` to < 1e-12
on every SCM, or the run is void.

Three arms, all on the **same** SCMs:

| arm | knowledge `K` | perturbation | keeps `K'` a tiering? |
|---|---|---|---|
| **G-flip** | generic: orientations of undirected edges of `C`, true, capped at \|K_T\| | reverse ρ ∈ {1,2,3} statements | n/a (never was one) |
| **T-flip** | tiered: **full** cross-tier orientation set, no cap | reverse ρ ∈ {1,2,3} statements | **NO** — this is the pilot's arm |
| **T-move** | tiered: **full** cross-tier set, no cap | move ρ ∈ {1,2} nodes to a different tier | **YES** — this is the new experiment |

Two code defects in the pilot are fixed first, and both are load-bearing:

1. `run_linear.py:55` — `K = [(u,v) for (u,v) in K if D[u,v]==1]` silently **drops every
   statement that disagrees with the true DAG**. Under a perturbed tiering those statements
   are exactly the signal. **Removed.**
2. `run_linear.py:56-58` — subsamples tiered `K` to `MAX_K = 4`. Bang & Didelez's R1-only
   property is about the **full** set of cross-tier constraints; a random size-4 subsample
   is not that object. **Cap removed for both tiered arms;** the generic arm is matched to
   `|K_T|` per SCM so the knowledge budget is equal.

New machinery written for this experiment (absent from the pilot):
`meek_closure_rules(G, rules)` with a per-rule firing counter, and
`apply_background_knowledge_rules(C, K, rules)`.

## 2. Metrics

Let `U` = undirected edges of `C`; `G0` = MPDAG under the true `K`; `G'` = MPDAG under `K'`.
"Orientation status" of an edge ∈ {i→j, j→i, undirected}.

- **n_changed(K→K')** = `|K Δ K'|` over *oriented pairs*, counting a reversed edge as **1**,
  an added statement as 1, a removed statement as 1. This is the perturbation size, and it
  is what makes a tier move and a statement flip comparable.
- **Δorient** = # edges of `U` whose orientation status differs between `G0` and `G'`.
  *This is literally the quantity the theorem is about: "edges whose orientation changes
  per unit perturbation."*
- **rate r = Δorient / n_changed** — the primary statistic.
- **cascade** = # edges of `U` oriented in `G'` **beyond** those directly named in `K'`
  (rule-propagated orientations). Secondary/mechanistic.
- **R1-sufficiency** = 1 if `apply_bk(C, K', rules={R1})` equals `apply_bk(C, K', rules={R1..R4})`.

Only **Meek-consistent** members enter the rate comparison (inconsistent ones are "caught
free" and have no MPDAG).

## 3. PRE-REGISTERED DECISION RULE

### GATE 1 — the premise (R1-closure must survive the perturbation)

Let `s` = R1-sufficiency rate over all consistent members of **T-move**.

- `s ≥ 0.99` → premise holds; proceed to Gate 2.
- `s < 0.99` → **REFUTED.** A perturbed tiering is no longer R1-closed, so the "R1-closed
  sublattice" is not closed under the perturbation the theorem would range over, and the
  theorem has no object.

### GATE 2 — containment (the primary falsifier)

Compare `r` between **T-move** (R1-closed) and **G-flip** (generic, R1–R4 free).

- **REFUTED** if `median(r_T-move) > median(r_G-flip)` **OR** `q90(r_T-move) > q90(r_G-flip)`,
  **and** a one-sided Mann–Whitney U test of `H1: r_T-move > r_G-flip` gives `p < 0.05`.
  *(R1-closed knowledge cascades MORE ⇒ "R1 controls the cascade" is false ⇒ the theorem's
  premise is backwards.)*
- **SUPPORTED** if `median(r_T-move) ≤ median(r_G-flip)` **AND** `q90(r_T-move) ≤ q90(r_G-flip)`
  **AND** `max(r_T-move) ≤ max(r_G-flip)` **AND** one-sided MWU of `H1: r_T-move < r_G-flip`
  gives `p < 0.05`.
- **INCONCLUSIVE** otherwise.

### GATE 3 — finiteness flavour (secondary, reported, NOT decisive)

Distribution of `|{distinct O* over the ball}|` per SCM, and the tail of `|bias|/|τ|`
conditional on `O*` changing, per arm. No threshold: the campaign was burned once by a
heavy tail, so the tail is reported and never summarised by a mean alone.

### Overall

**REFUTED** if Gate 1 or Gate 2 returns REFUTED. **SUPPORTED** only if both pass.
Otherwise **INCONCLUSIVE**.

## 4. C13 compliance (vault registry rule)

"Any null in this vault must state how its threshold parameters were chosen and show they
were capable of producing a non-null."

- The thresholds are **not** free constants: Gate 2 is a *matched two-sample comparison
  against an in-house control measured on the same SCMs with the same estimator*. The only
  free constant is the MWU α = 0.05 and the 0.99 in Gate 1.
- Capability check, mandatory, reported in RESULTS.md: the instrument must be shown able to
  return **both** verdicts. Specifically, `r` must have non-degenerate spread in **both**
  arms (≥ 3 distinct values with > 5% mass each), and the **T-flip** arm is included
  precisely as a positive control — it is *known* from the pilot to cascade more than
  generic (0.760 vs 0.278 extra orientations), so if the pipeline cannot reproduce
  `T-flip > G-flip` it is broken and no verdict is readable.

## 5. Prior expectation, recorded so it cannot be retrofitted

The pilot's published tiered-vs-generic cascade (0.760 vs 0.278) and the ground-truth note
that *"the tiered arm measured R1 as the LONG-RANGE propagator (50.1% vs 19.9% of flips
cascading)"* both point the **wrong way for the theorem**. R1 is the rule that chains along
chordless undirected paths of arbitrary length; R2–R4 are local shielded-triple rules. So
the honest prior is that **R1-closure buys tractability of the RULE SET, not containment of
the CASCADE**, and Gate 2 is expected to fire REFUTED.

Recording this in advance is the point: if Gate 2 fires, it is a pre-registered kill, not a
post-hoc story.
