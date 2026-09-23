# Results

One directory per sweep. Each carries a `manifest.json` (or per-shard
manifests) recording the resolved configuration, git SHA, interpreter and
package versions, solver seed and a determinism block, and most carry notes
(`*NOTES.md`, `FINDINGS.md`, `PREREGISTRATION.md`, `VERIFY*.txt`) written when
the sweep ran. Paths are relative to the repository root, which is where every
driver must be run from.

A radius of `-1` (`UNREACHED`) in any file is a status, not a number: no
reachable knowledge state invalidates the adjustment set. It is never averaged.

## Behind the paper

| Directory | Content | Written by |
|---|---|---|
| `axis_robustness_p6/` | Synthetic survival under flip and tiered corruption: 4,104 instances, N=1000 draws per grid point | `bkrobust.robustness.run_survival_p6` |
| `axis_robustness_p6_paired/` | Paired bootstrap differences in τ_b on that sweep | `experiments/paired_synth_p6.py` |
| `elicit/` | Orientations elicited from six LLMs under ten conditions (`knowledge.json`) | `experiments/questionnaire_build.py`, `elicit_run.py`, `elicit_import.py` |
| `axis_robustness_real/` | Real-network query frame (`frame.jsonl`, 831 rows) and survival with generated knowledge | `bkrobust.robustness.real_frame`, `run_real_survival` |
| `axis_robustness_llm/` | Real-network survival with LLM-elicited knowledge, all conditions | `scripts/run_llm_survival_panel.sh` |
| `axis_robustness_real_p6/` | Extra baselines and paired differences on real structure | `bkrobust.robustness.real_baselines_p6` |
| `axisa3/` | Real-network corpus (`networks/`, checksummed) and its structural measurements | `bkrobust.benchmarks` |
| `e2e_speedup_gated/` | End-to-end timing of space, search and hybrid on 256 instances | `experiments/e2e_speedup.py --gate` |
| `axisb3/`, `axisb4/`, `search/` | CP-SAT encodings, oracles, differential validation, search speed-ups | `bkrobust.sat`, `bkrobust.search`, `bkrobust.hybrid` |
| `epsilon/` | Bias radius: exhaustive anti-exchange check (`chickering/`), verification, counterexample search, synthetic study, finite-sample check, worked example | `bkrobust.epsilon.*`, `experiments/run_r_epsilon*.py` |
| `final_table/` | Real-network radius panel, `r_val` and `r_eps` per query | `experiments/final_table.py` |
| `final_table_eps_gt1/` | The same panel at bias tolerances above 1 | `experiments/final_table.py` |
| `mean_table_fullspace/` | Mean bias radius `r_mu` per query, on the full knowledge-state space | `experiments/mean_table_fullspace.py` |
| `breakdown_radius_demo/` | The running example, three knowledge scenarios | `python -m bkrobust.demo.run_all` |

## Supporting and earlier sweeps

| Directory | Content | Written by |
|---|---|---|
| `axis_robustness/` | First survival sweep (N=200), claims-vs-hops units, null cell, cross-arm, Pareto frontier | `bkrobust.robustness.run_*`, `experiments/run_analyse.py` |
| `axisa2/` | Synthetic census, separation frontier, gate differential | `bkrobust.synth`, `bkrobust.benchmarks` |
| `axisb2/` | Lattice structure, Conjecture 2 sweep, bounds, scaling | `bkrobust.search`, `bkrobust.core` |
| `synth/` | Synthetic ensembles: pre-registration, pilot, calibration, census, main run | `bkrobust.synth.runner` |
| `frame/` | Real-network candidate frame for elicitation | `experiments/frame_build.py` |
| `stage0/` | Claim-radius and ladder audits on the real corpus | `experiments/stage0_*.py` |
| `zero_radius/` | Audit of degenerate (`r_val = 0`) queries in `final_table/` | `experiments/zero_radius_audit.py` |
| `matched_control/` | Random knowledge matched to the LLM arm, as a control | `experiments/matched_control.py` |
| `mean_vs_max/`, `mean_fullspace/`, `mean_table/`, `mean_epsilon/`, `mean_bounds/`, `mean_monotonicity/`, `mean_adaptive/`, `mean_sampling/` | Design studies for the mean bias profile (see `docs/`) | `experiments/mean_*.py` |
| `e2e_speedup/` | Timing sweep before the hybrid dispatch fix; superseded by `e2e_speedup_gated/` | `experiments/e2e_speedup.py` |
