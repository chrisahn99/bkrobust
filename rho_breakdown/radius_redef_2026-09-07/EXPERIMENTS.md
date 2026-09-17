# Radius redefinition — design, audit, and the experiments actually run

**Design phase:** 67 agents, 9.1M tokens, 3.6 h. 18 candidate definitions → 6, after two lemmas
collapsed most of them into one object. **Execution phase:** an independent harness written here,
verified against `~/bkrobust`, three arms, ~2.6 s.

---

## 1. The artefact, reproduced independently and sharpened

`harness.py`, 400 random problems, 2.9 s:

| `r_val` | share |
|---|---|
| **1** | **47.5 %** |
| 2 | 1.5 % |
| 3 | 0.8 % |
| **UNREACHED** | **50.2 %** |

**97.7 % of the mass sits on {1, ∞}.** The statistic is not merely biased toward 1, it is very
nearly **binary**: it answers *does some single retraction break the adjustment set*, and almost
nothing else. A quantile or a mean of a near-binary variable is still near-binary, which is a
prediction against the whole "less strict radius" family, made before the design phase reported.

---

## 2. Two lemmas the design phase proved and measured

**L1 — retraction equals reversal.** The union of retractions of size ≤ b equals the union of
reversals of size ≤ b, because a retraction is the disjoint union of the reversals it permits.
Measured: 0 differences over 573 problems. Several apparently distinct proposals were the same
object.

**L2 — screens are width-neutral.** For any claim whose retraction leaves the adjustment set
valid, the hull is unchanged. Hence the hull over the *load-bearing* claims equals the hull over
*all* claims. Measured 250/250 and 300/300.

🔴 **Consequence: direction (4), identifying the exact false claims, has no width content at
all.** It cannot enter the width table. It is a **cost** result (the same interval from fewer
enumerations) and a **retrieval** result, and it must be presented as one.

---

## 3. The decisive experiment, run here

`confirm.py`. Random DAGs on 7 nodes, edge probability 0.3, CPDAGs with 2–5 undirected edges,
up to 4 elicited claims, population covariance, effects by the IDA convention over DAG
extensions. 300 problems per arm.

🐛 **A bug found by reading the output.** The first run sampled the expert's claims *truthfully*,
so there were no wrong claims and coverage was 1.000 at every budget. Coverage that is perfect
everywhere is a tell, not a result. Fixed by reversing `n_false` of the claims, keeping only
draws that survive Meek closure, so the population is the errors the consistency check does not
catch for free.

### One elicited claim reversed

| hedge | covers | width / blanket | equals blanket |
|---|---|---|---|
| none (report as it stands) | 0.937 | 0.000 | 0 % |
| **one revision** | **1.000** | **0.710** | 67.2 % |
| two revisions | 1.000 | 0.945 | 91.8 % |
| blanket over the class | 1.000 | 1.000 | 100 % |

### Two elicited claims reversed

| hedge | covers | width / blanket | equals blanket |
|---|---|---|---|
| none | 0.850 | 0.000 | 0 % |
| one revision | 0.927 | 0.450 | 43.1 % |
| **two revisions** | **1.000** | **0.903** | 90.3 % |
| blanket | 1.000 | 1.000 | 100 % |

**The headline.** The first budget reaching 95 % coverage costs **0.710** of the blanket width
with one claim wrong, and **0.903** with two. So selective hedging saves **29 %** of width at one
error and **10 %** at two.

⚖️ **Against the threshold fixed in the plan (ρ★ ≤ 0.85):** it **passes at one error and fails at
two.** The design phase's own pilot reported 0.936 and read it as the width bet dying; this run
says that verdict belongs to the two-error arm, and the one-error arm is a genuine 29 % saving.
Two different generators, same ordering, and the disagreement is about *which arm was quoted*.

⚠️ **And the saving comes from a minority of problems.** At one error the selective hull *is* the
blanket on 67.2 % of them. The mean is carried by the remaining third.

---

## 4. 🔴 The finding that reframes the paper

Measured three times, on two independent generators:

| | design phase | one error (here) | two errors (here) |
|---|---|---|---|
| informative problems where the class disputes whether X causes Y | 73.1 % | **86.9 %** | 84.7 % |
| width restricted to worlds where Y descends from X, as a share of blanket | 0.335 | **0.314** | 0.359 |
| that width is exactly zero | 57.9 % | **60.7 %** | 58.3 % |

**Roughly two thirds of the ignorance this construction hedges is the equivalence class
disagreeing about whether the treatment affects the outcome at all, not about which confounders
to adjust for.** On about 60 % of informative problems, once the worlds that agree on the causal
direction are isolated, the interval collapses to a point.

If this is not printed, a width table over the class is a table about **causal-order disputes
wearing an adjustment-set label**. It is the single most consequential thing either phase found,
and it was found by an audit, not by a proposal.

---

## 5. A scope fact neither phase had

Only **20–24 %** of drawn problems are informative at all: on the rest, the equivalence class
leaves no ignorance about the effect and there is nothing to hedge. Every width number above is
conditional on that stratum, and the paper has to say so.

---

## 6. What this changes

- **Direction (4) leaves the width table.** L2 is decisive. It becomes a cost and retrieval result.
- **Direction (2) and (3), the softer radius, are in trouble** before they start: the underlying
  variable is near-binary, so quantiles and means of it inherit that.
- **Direction (1), the interval on the treatment effect, is the survivor**, and it is the one that
  feeds the width-at-95 %-coverage metric directly.
- **The paper's headline width claim should be stated at one reversed claim**, where it is 0.710,
  with the two-error arm printed beside it at 0.903 rather than omitted.
- 🔴 **The direction-dispute stratum must be a co-primary column**, not an appendix note.

## Files

`harness.py` problem generator and `r_val` reference · `confirm.py` the three-arm experiment ·
`plan.json` the full design-phase plan and audits · ~99 audit scripts the design agents wrote
and ran while reasoning.
