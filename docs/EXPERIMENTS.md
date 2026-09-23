# Experiments

What each experiment establishes, what would falsify it, and what it cannot
show. Ordered by dependency: later experiments assume earlier ones came out.

Every experiment is driven by a Hydra config, takes a `seed`, and writes to
`results/<experiment>/<timestamp>_<git-sha>/`.

---

## Control arms

Three arms run at every radius in every sweep. Omitting any of them makes the
others uninterpretable.

| Arm | Knowledge | Answers |
|---|---|---|
| `no_knowledge` | none | What does the CPDAG alone give? |
| `true_knowledge` | consistent **and true**, size-matched | What does correct knowledge buy? |
| `false_knowledge` | consistent but false | The object of study. |

The size match matters. Knowledge shrinks the equivalence class whether or not
it is correct, and some of the resulting change in `O*` has nothing to do with
falsity. Comparing false knowledge against *no* knowledge conflates the two
effects; comparing it against size-matched true knowledge does not.

---

## E1 — Synthetic sweep (`run_synthetic_sweep.py`)

**Establishes.** The two breakdown radii exist, are computable, and are ordered
`δ_opt ≤ δ_valid` (theorem target T1). The primary structural result.

**Protocol.** Sample a DAG → attach an SCM → sample data → compute the exact
CPDAG from the true DAG → select a target pair → compute both radii exactly →
sweep the radius grid over all three arms → record structural and effect metrics
together.

The CPDAG comes from the *true DAG*, not from running discovery on the sampled
data. That isolates the effect of the knowledge from the effect of discovery
error; a radius measured against an estimated CPDAG confounds the two and is not
comparable. E4 uses the estimated CPDAG, which is the practitioner's situation.

**Falsified by.** An exact-method violation of the ordering. That is a
counterexample to T1 and the single most valuable output of this codebase — the
witness is persisted in full rather than aggregated into a rate.

**Cannot show.** Anything about real data. Synthetic graphs and SCMs are drawn
from families chosen for tractability, and the radii are properties of those
families until E3 and E4 say otherwise.

**Watch for.** High varsortability. If the data leaks its causal order through
marginal variances, discovery recovers the order regardless of the knowledge,
the knowledge never gets to matter, and the radii look large for a reason that
has nothing to do with robustness. Recorded per replicate; see `data/leakage.py`.

---

## E2 — Cascade study (`run_cascade_study.py`)

**Establishes.** The mechanism. `k` imposed constraints force `m ≥ k`
orientations, and the amplification factor `m / k` explains why small knowledge
errors have large structural consequences.

**Protocol.** For each graph family and density, impose `k = 1, 2, ...`
consistent constraints of each kind and record the forced count, the cascade
reach, the per-rule attribution, and whether `O*` moved.

**Report separately.** R4's contribution. R4 only fires once background
knowledge is present, so its share is the part of the cascade that exists
*because* the analyst supplied knowledge rather than because the data determined
something.

**Cannot show.** That amplification causes harm. A cascade that never reaches
the parents of the causal nodes costs nothing measurable, which is why the hit
rate on `O*` is recorded alongside the raw factor. Amplification is necessary
for harm, not sufficient.

---

## E3 — Semi-synthetic (`run_semisynthetic.py`)

**Establishes.** That the radii survive realistic covariate distributions —
i.e. that E1 is not an artefact of synthetic marginals.

**Protocol.** Real covariates, simulated outcomes (DREAM, RealCause). True graph
known, so the radii are still computable; heuristic path where networks are too
large for exact enumeration, with those radii labelled as upper bounds.

**Falsified by.** Radii that differ systematically from E1 at matched graph size
and density. That would mean the synthetic families are not representative and
E1's numbers do not transfer.

**Cannot show.** Robustness on the real process the covariates came from. The
true graph here is the *simulator's*, and it is true of the simulated outcome,
not of whatever generated the covariates.

---

## E4 — Real data (`run_real_data.py`)

**Establishes.** That the certificate is computable and informative on real
datasets, and that real elicited knowledge sits at radii where it matters.

**Protocol.** Load the dataset → read or elicit the analyst's knowledge →
bootstrap the discovery step → check consistency on each replicate → certify on
each → report the spread of radii, verdicts, and reachable conclusions.

The bootstrap is not optional. The CPDAG is an estimate, and a certificate
computed against a single point estimate understates the uncertainty by exactly
the discovery error. A knowledge set that FAILs on some bootstrap replicates and
not others is itself worth reporting, since the analyst would have seen one draw.

**Cannot show.** Bias, or either radius as a measured quantity. There is no
ground-truth DAG. What is measurable: the spread of effect estimates across the
ball of knowledge that all passes the same consistency check — the range of
conclusions an analyst could have reached while doing everything correctly.

---

## E5 — CFM audit (`run_cfm_audit.py`)

**Establishes.** Whether amortized causal models inherit the same breakdown
structure, are more robust, or are less.

**Protocol.** Data fixed; knowledge varied over radius, conditioning mode and
bias scale, across the registered checkpoints, with a classical OLS-on-`O*` arm
on identical data.

**Two gates run first, and abort the sweep rather than warn.**

1. **Zero-scale identity.** `bias_scale = 0.0` must reproduce `mode="none"` bit
   for bit. Otherwise the injection is changing the forward pass through some
   route other than the intended bias, and every downstream number measures that
   bug.
2. **True knowledge helps.** The true-knowledge arm must beat the unconditioned
   arm. Otherwise the model is not reading the conditioning input at all, and a
   flat bias-versus-δ curve looks like robustness while meaning nothing.

Report the gate results whatever they say. A checkpoint that ignores its
conditioning input is a finding about that checkpoint.

**Cannot show.** That attention-bias injection *is* knowledge conditioning.
Imposing knowledge on a CPDAG has a semantics — the forced orientations are
entailed. An attention bias has none; it is a nudge, and the model may follow
it, ignore it, or overreact. The audit measures what these models do with a
conditioning signal, and the paper must say so in those terms.

---

## E6 — Knowledge embedding (`run_knowledge_embedding.py`) — CONDITIONAL

**Establishes.** Whether the constraint poset's learned geometry predicts `δ` or
separates true from consistent-but-false knowledge.

**Protocol.** Build the poset per knowledge set → embed → probe → compare
against the structural baselines and the shuffle control. Splits are **by
graph**: knowledge sets from one graph share its structure, so a random split
leaks and the score is meaningless.

**Two reasons for scepticism, both of which shape the protocol.**

The radius is exactly computable at the graph sizes that matter, so a learned
predictor of it is a claim about speed or about generalisation to larger graphs
— not about access to something otherwise unavailable. Say which is being
claimed.

Separating true from consistent-but-false knowledge would be a strong result,
and consistent-but-false knowledge is *by construction* indistinguishable from
true knowledge given the CPDAG. So a probe that separates them cannot be
detecting truth; it is detecting how the two arms were generated. The shuffle
control (`representations/probes.py`) is reported next to every score, always.

**Cut if.** The probes do not beat constraint count, cascade size, and SHD.

---

## Reporting rules

These apply to every experiment and exist because each names a specific way a
number can be quietly wrong.

1. **Label heuristic radii as upper bounds.** In the `exact` flag *and* in the
   rendered text. The difference between a guarantee and a guess lives there.
2. **`UNREACHED` is not a large radius.** It means no breaking perturbation was
   found within budget. Aggregating it as a number turns a budget limit into a
   robustness claim.
3. **Report signed and absolute bias.** Signed bias averages toward zero under
   symmetric perturbations, which reads as "no bias".
4. **Report bias in standard-error units too.** That is what determines whether
   anyone would notice.
5. **Never count a missing interval as non-covering.** Most CFM checkpoints
   report no interval; counting those as failures fabricates a finding.
6. **Carry leakage diagnostics on every synthetic and semi-synthetic row.**
7. **Record which CPDAG source produced each row** — true-DAG or discovery.
   The two are not comparable.
8. **Persist witnesses for ordering violations individually.** Never aggregate a
   potential counterexample into a rate.
