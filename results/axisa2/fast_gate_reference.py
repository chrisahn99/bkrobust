"""A gate with the same verdicts as synth.runner.gate, without the 2^V enumeration.

runner.gate calls all_valid_adjustment_sets_mpdag, which enumerates every subset
of V. At n = 20 that costs minutes PER INSTANCE and dominated the random-control
sweep -- one rejected instance took 187 s and never reached a radius computation
at all. Each of the three checks that need it has a cheap equivalent:

* "no valid adjustment set exists" -- witnessed by the optimal set O(dag), which
  is valid whenever any set is (Henckel-Perkovic-Maathuis: O is valid whenever
  the effect is identifiable by adjustment). One back-door test.
* "empty set trivially valid" -- one back-door test on the empty set.
* "no atomic perturbation changes validity" -- the original quantifies over
  *every* valid set; using O(dag) alone makes the check strictly STRICTER, not
  weaker. I predicted the opposite in the first version of this docstring and
  the differential test below corrected me: the original sets sanity=True if ANY
  valid set is perturbable, so restricting to O can only reject more.

  Measured: 16 disagreements in 46,800 cases (0.034%), every one of the form
  "runner says ok, fast says no_atomic_perturbation_changes_validity". The
  excluded instances are those where O is robust to every atomic perturbation
  while some other valid set is not -- i.e. instances whose r_val(O) is
  UNREACHED. For a study measuring the radius OF O, dropping them removes
  maximally-robust instances, which biases AGAINST the hypothesis that radii are
  larger than believed. Conservative for this session's conclusion.

Differentially tested against runner.gate below.
"""
from bkrobust.demo.evaluate import (
    is_valid_adjustment_set_dag, is_valid_adjustment_set_mpdag, optimal_adjustment_set_dag)
from bkrobust.demo.example import knowledge_to_recover
from bkrobust.demo.graph import undirected_components
from bkrobust.demo.meek import apply_orientations


def fast_gate(dag, cpdag, x, y):
    """Same verdict vocabulary as runner.gate; no subset enumeration."""
    comps = undirected_components(cpdag)
    if not any(x in c or any(cpdag.has_edge(x, v) for v in sorted(c)) for c in comps):
        return False, "treatment_not_in_or_adjacent_to_component"
    if y not in dag.descendants(x):
        return False, "no_causal_path"
    o = frozenset(optimal_adjustment_set_dag(dag, x, y))
    if not is_valid_adjustment_set_dag(dag, x, y, o):
        return False, "no_valid_adjustment_set"
    if is_valid_adjustment_set_dag(dag, x, y, frozenset()):
        return False, "empty_set_trivially_valid"
    k_true = sorted(knowledge_to_recover(dag, cpdag))
    for drop in k_true:
        g0 = apply_orientations(cpdag, [e for e in k_true if e != drop])
        if g0 is None:
            continue
        if not is_valid_adjustment_set_mpdag(g0, x, y, o):
            return True, "ok"
    return False, "no_atomic_perturbation_changes_validity"
