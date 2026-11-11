# E21 preregistration: compute the walkthrough toy exactly before drawing it

Dated 2026-11-06. Not to be edited after the first run. Code: submitted branch, commit 5421cabc.

**Question.** Does the planned four-node walkthrough toy (CPDAG U-X, X-Y, X-W, W-Y, all undirected; K = {U->X, X->Y, X->W}; query X->Y) have the optimal set, up-set, radius, witnesses and claim radius we intend to draw? What is the radius of the Figure 2(a) toy (CPDAG A-W, W-X, W-Y, X-Y; K = {A->W, X->Y}; Z = {W})?

**Registered statistic.** Per toy, the breakdown radius r_val of Z = O(G0) under the generalized adjustment criterion (GAC), returned by `bkrobust.hybrid.breakdown_radius`, accepted only if five independent computations agree: (1) the hybrid with the GAC criterion oracle, (2) the hybrid with the enumeration oracle, (3) the hybrid forced onto the SAT ladder (search budget 0), (4) BFS distance from G0 to the nearest failing state on the covering graph of the fully enumerated space (`space_fixed.build_space_correct`), (5) the closed form min over D in [C] with Z invalid of |K_G0 \ dir(D)| (Cor. distance (c)). The values are exact, so the "interval" is the degenerate [r, r]; there is no sampling. Alongside: |space|, the up-set of G0 with each state's validity, the nearest failing states (witnesses), the DAGs covered by G0, and the claim radius r_claim (all subsets of K retracted, cross-checked against its closed form, rem:rclaim-form).

**Prediction (hand computation, model-only).**
- Walkthrough: Z = {} (empty). Up-set of G0 has 5 states: G0 = {U->X, X->Y, X->W} valid, {X->Y, X->W} valid, {X->Y} invalid, {X->W} invalid, the CPDAG invalid. Two DAGs directly below G0, both valid. r_val = 2, witnesses {X->Y} and {X->W}. r_claim = 2. The subset {U->X, X->W} is not a state (closure re-adds X->Y). |space| not predicted.
- Figure 2(a) toy: Z = {W}, r_val = 1 (retracting X->Y leaves X-Y undirected). Panel (b) of Figure 2 draws r_val = 3, which is schematic.
- Running example, scenario B (added, see below): Z = {Age, Smoke}, r_val = 3, 48 states.

**Outcome against the authors.** Any disagreement between the code and the hand values: the figure and the response are drawn from what the code prints, and the hand values are withdrawn. If the Figure 2(a) toy has r_val = 1, the response says so and Figure 2(b) is relabelled or redrawn with the computed value. If the running example does not give 48 states and r_val = 3, Figure 1's caption is corrected.

**Populations, seeds, exclusions.** Three fixed toys, no randomness, no seeds (the SAT ladder runs single-threaded with seed 0). Nothing is excluded. The smoke run computes the two four-node toys; the full run adds the running example.

**Consistency failures (exit non-zero, row still written).** The five radii disagree; the hybrid's witness is not among the nearest failing states; the criterion oracle and the enumeration oracle disagree on any state; Lemma O (inclusion of models = reverse inclusion of directed edges) fails on any pair; Lemma 2 (join = Meek(C, K_G intersect K_H)) fails on any pair; the CPDAG is not the CPDAG of its own extensions; r_claim from subsets differs from its closed form, or r_claim > r_val; Z is not valid at G0.

**Stopping rule.** One pass over the toys. No rerun after seeing values, except to fix a crash, recorded in the run log.

**Deviations from the plan.** (1) The running example (scenario B of `demo/scenario.py`) is added as a third unit because the drafted reply to C2Bd-W1/Q2 promises a figure on it; the plan named only the two four-node toys. (2) The claim radius is computed over all subsets of K, not only single retractions as in `experiments/stage0_claim_radius.py`, since the toy's predicted r_claim is 2. (3) The back-door radius is recorded beside the GAC one, as information only.
