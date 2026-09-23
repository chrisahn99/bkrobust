# Reports

Development record of the project: one report per working session, the
hand-off notes written between sessions, and the planning documents they
were written against. They are the sources for the numbers in the paper, not
drafts of it, and several of their early headlines were later superseded;
the later session reports say which.

Nothing in this folder is read by the code as input. The files live here only
so the repository root stays readable. Three kinds of code do refer to them by
path, and those paths were updated when the files moved:

- the verifiers `bkrobust.analysis.session{3,4,5,6}_verify`, which re-assert
  every number quoted in a report against the committed files under `results/`;
- the PDF builders `bkrobust.analysis.session*_pdf`, which write the PDF
  renderings in `sessions/`;
- the anchor-table builders `bkrobust.robustness.build_tau_comparisons` and
  `bkrobust.robustness.build_table_tau_real`, which write `tables/`.

Run all of them from the repository root, e.g.

```bash
PYTHONPATH=src python -m bkrobust.analysis.session6_verify
```

## `sessions/`: session reports, in order

| # | Date | Report | PDF | Verifier | Topic |
|---|---|---|---|---|---|
| 0 | 2026-09-03 | [`report.md`](sessions/report.md) | | | Fully worked example of the breakdown radius on a small clinical graph (`bkrobust.demo`) |
| 1 | 2026-09-03 | [`report_synth_and_search.md`](sessions/report_synth_and_search.md) | [`SESSION_REPORT.pdf`](sessions/SESSION_REPORT.pdf) | | Synthetic-graph generalisation and search heuristics |
| 2 | 2026-09-05 | [`report_axisb_deep.md`](sessions/report_axisb_deep.md) | [`REPORT_AXISB_DEEP.pdf`](sessions/REPORT_AXISB_DEEP.pdf) | | Axis B in depth: proof, scale and bounds |
| 3 | 2026-09-05 | [`report_axisb_sat.md`](sessions/report_axisb_sat.md) | [`REPORT_AXISB_SAT.pdf`](sessions/REPORT_AXISB_SAT.pdf) | `session3_verify` | Declarative (CP-SAT) encodings of the radius |
| 4 | 2026-09-05 | [`report_axisb_oracle.md`](sessions/report_axisb_oracle.md) | [`REPORT_AXISB_ORACLE.pdf`](sessions/REPORT_AXISB_ORACLE.pdf) | `session4_verify` | The validity oracle, enumeration overhead, reach of exact radii |
| 5 | 2026-09-06 | [`report_saturation_gac.md`](sessions/report_saturation_gac.md) | [`REPORT_SATURATION_GAC.pdf`](sessions/REPORT_SATURATION_GAC.pdf) | `session5_verify` | Does the radius saturate at 1? It equals a structural separation |
| 6 | 2026-09-06 | [`report_real_graphs.md`](sessions/report_real_graphs.md) | | `session6_verify` | Separation on real benchmark networks vs. Erdős–Rényi |
| 7 | 2026-09-10 | [`report_fragility_and_pareto.md`](sessions/report_fragility_and_pareto.md) | [`REPORT_FRAGILITY_AND_PARETO.pdf`](sessions/REPORT_FRAGILITY_AND_PARETO.pdf) | | Survival under corrupted knowledge; Pareto frontier of assumptions |
| 8 | 2026-09-17 | [`report_r_epsilon.md`](sessions/report_r_epsilon.md) | | | The ε-bias radius and an unconditional `r_val` |
| 9 | 2026-09-16 | [`report_real_survival.md`](sessions/report_real_survival.md) | | `session9_verify` | [RE-11] Survival and the paired cross-arm on real graph structure |

Session numbers are the ones the reports use for themselves and that the
verifier and builder modules are named after.

## `handoffs/`: summaries written for the next session

| Date | File | Content |
|---|---|---|
| 2026-09-16 | [`RESUME.md`](handoffs/RESUME.md) | State of the [RE-11] campaign (session 9), phase by phase |
| 2026-09-22 | [`SESSION_SUMMARY_POINTS67.md`](handoffs/SESSION_SUMMARY_POINTS67.md) | Reviewer points 6 and 7: N=1000 synthetic survival rerun and new baselines |
| 2026-09-22 | [`SESSION_SUMMARY_SURVIVAL_ON_LLM.md`](handoffs/SESSION_SUMMARY_SURVIVAL_ON_LLM.md) | Survival and ranking with LLM-elicited knowledge |
| 2026-09-23 | [`SESSION_SUMMARY_E2E_SPEEDUP.md`](handoffs/SESSION_SUMMARY_E2E_SPEEDUP.md) | End-to-end speedup sweep and the hybrid dispatch fix |

## `planning/`: plans and logs

| File | Content |
|---|---|
| [`PLAN.md`](planning/PLAN.md) | The living plan from session 1 on, with its recorded revisions |
| [`NEXT.md`](planning/NEXT.md) | Status after session 2: what was safe to claim, what came next |
| [`LOG.md`](planning/LOG.md) | Audit trail of sessions 1–6 |
| [`section_7_narrative_outline.md`](planning/section_7_narrative_outline.md) | Outline of the simulation-study section of the paper |

## `tables/`: publication anchor tables

| File | Built by | Content |
|---|---|---|
| [`table_tau_comparisons.md`](tables/table_tau_comparisons.md) | `bkrobust.robustness.build_tau_comparisons` | Kendall τ_b(predictor, survival AUC) by stratum, synthetic structure (N=200 run) |
| [`table_tau_real.md`](tables/table_tau_real.md) | `bkrobust.robustness.build_table_tau_real` | The same table on real graph structure |

The numbers in the submitted paper come from the later runs in
`results/axis_robustness_p6*/` and `results/axis_robustness_llm/`, not from
these two tables.

Other working notes (experiment briefs, paper-integration briefs, result
notes) are in [`../docs/`](../docs/README.md).
