# Prior art on the retraction operator

Working note, 12 September. The plan gates the build on one question, asked
before any elicitation budget is spent: does the object we are about to claim
as new already have a name somewhere else?

**It does.** The general object is a minimal correction set. The causal
instantiation is where whatever is new has to live.

## The object, stated so it can be recognised under another name

Given a set `K` of asserted constraints, here orientations of edges in a
partially directed graph closed under Meek's rules, and a property that
currently holds, here the validity of a chosen adjustment set, we report the
cardinality of the smallest subset of `K` whose removal makes the property
fail. Failure is upward-closed: adding constraints cannot break the property,
only removing them can.

## Verdict by area

Eight areas were searched on 12 September, none of which had ever been queried
on this project.

| area | the object there | ours? |
|---|---|---|
| Stability radius in scheduling and combinatorial optimisation | largest data perturbation preserving optimality | no, a continuous magnitude in a norm |
| Radius of robust feasibility | smallest perturbation losing feasibility | no, continuous |
| Renegar's distance to ill-posedness | norm distance to the nearest ill-posed instance | no, continuous |
| AGM contraction and epistemic entrenchment | which beliefs to give up, ordered by entrenchment | close in spirit, but it prescribes an order rather than reporting a cardinality |
| **Minimal correction sets, MUS and MaxSAT diagnosis** | **smallest subset of constraints whose removal restores satisfiability** | **yes, this is the object** |
| Inconsistency measures for knowledge bases | degree of inconsistency, often via minimal inconsistent subsets | related structure, different question |
| Specification curve and multiverse analysis | the distribution of results across defensible specifications | no, a distribution rather than a minimum |
| Breakdown frontiers in econometrics | the boundary in a continuous assumption space | no, continuous, and it is the analogy we already cite |

A ninth query asked whether any of the eight has been applied to causal graphs,
background knowledge, Meek closure or adjustment-set validity. Nothing was
found. The nearest results are the knowledge-informed local discovery line, the
tiered-knowledge equivalence work, and the imperfect-constraints discovery paper
of Wang et al., and none of them computes a minimal removal cardinality or
reports one as a certificate.

## What this changes

**The operator is not ours to name.** "The smallest subset of asserted
constraints whose retraction breaks a property" is a minimal correction set,
and the paper should say so in the sentence that introduces it rather than let
a reviewer say it first. The narrow thing that survives is the instantiation:
the constraint set is a knowledge set closed under Meek's rules, the property
is adjustment-set validity in a maximal PDAG, and the certificate is reported
against an analyst's own claims.

**It buys a literature rather than costing one.** MCS enumeration has thirty
years of solver work behind it, including the hitting-set duality with minimal
unsatisfiable subsets and block-based enumeration. Our depth-3 exhaustive
enumeration is the naive algorithm for it. Where the enumeration is the
bottleneck, that literature is the fix, and the repository already has a CP-SAT
encoding to hang it on.

**It sharpens what has to be measured.** If the framework is standard, the
contribution is the measurement: what the certificate is worth against the free
statistics, and whether it holds when the truth is revealed. That is the
argument the evaluation plan is built around, and this verdict strengthens it
by removing a claim that would not have survived review.

## References, and their status

Returned by a literature sweep on 12 September. **They have not been checked
against the papers themselves**, and one is already known to be wrong: the
GRASP paper is 1999, not 2005, and it is a satisfiability solver rather than the
source of the correction-set definition. Verify every line below before any of
it enters a bibliography.

- Liffiton and Sakallah, "Algorithms for computing minimal unsatisfiable subsets
  of constraints", Journal of Automated Reasoning 40(1), 2008.
  doi:10.1007/s10817-007-9084-z
- Bacchus and Katsirelos, "Using minimal correction sets to more efficiently
  compute minimal unsatisfiable sets", SAT 2015.
  doi:10.1007/978-3-319-21668-3_5
- Hunter and Konieczny, "Measuring inconsistency through minimal unsatisfiable
  subsets", KR 2008.
- Gardenfors, *Knowledge in Flux*, MIT Press 1988, for contraction and
  entrenchment.
- Sotskov, Tanaev and Werner, "Stability radius of an optimal schedule", 1996.
  doi:10.1007/3-540-61310-2_26
- Goberna, Jeyakumar and Li, "Radius of robust feasibility for semi-infinite
  convex programs", Optimization 60(10-11), 2011.
  doi:10.1080/02331934.2010.546581
- Renegar, "Incorporating condition measures into the complexity theory of
  linear programming", Journal of Complexity 11(1), 1995.
  doi:10.1006/jcom.1995.1005
- Simonsohn, Simmons and Nelson, "Specification curve analysis", Nature Human
  Behaviour 4(11), 2020. doi:10.1038/s41562-020-0912-z
- Masten and Poirier, "Inference on breakdown frontiers", Quantitative Economics
  11(1), 2020. doi:10.3982/QE1288

## Still open

Whether the fixed-skeleton sub-poset walk and the correction-set enumeration are
the same search seen twice. They are defined on different objects, retraction
subsets against covering steps, and the measurements now on disk show they give
different numbers on the same row, so they are at least not interchangeable.
Whether one reduces to the other is a paragraph somebody has to write.

Whether an MCS solver beats the depth-3 enumeration at the sizes we care about.
On the committed corpus the median row is settled by closing a single subset, so
the answer today is no, and the question returns only if elicited knowledge
produces larger asserted sets.
