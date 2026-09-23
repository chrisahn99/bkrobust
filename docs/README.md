# Docs

Files stay at these paths because docstrings in `src/` and `experiments/`
cite them by name. Proved statements for the validity radius are in
[`../THEOREMS.md`](../THEOREMS.md); session reports are in
[`../reports/`](../reports/README.md).

## Reference

| File | Content |
|---|---|
| [`R_EPSILON_THEORY.md`](R_EPSILON_THEORY.md) | The ε-bias radius: statements, proofs, and what each one assumes |
| [`THEORY.md`](THEORY.md) | Original definitions and theorem *targets*, written before any result; see `THEOREMS.md` for what was proved |
| [`EXPERIMENTS.md`](EXPERIMENTS.md) | Planned protocols for the Hydra entry points in `experiments/run_*.py` (scaffold tree) |
| [`DATA.md`](DATA.md) | Dataset provenance, licences and manual download instructions |
| [`CHECKPOINTS.md`](CHECKPOINTS.md) | Causal-foundation-model checkpoint sources and licences |

## Result notes (bias radius and its mean companion)

| File | Content | Results |
|---|---|---|
| [`EXPERIMENT_MEAN_VS_MAX.md`](EXPERIMENT_MEAN_VS_MAX.md) | Exploratory: mean vs. max bias per shell | `results/mean_vs_max/` |
| [`RESULT_MEAN_SOUNDNESS_GATE.md`](RESULT_MEAN_SOUNDNESS_GATE.md) | The mean over the retraction up-set is unsound | `results/mean_fullspace/` |
| [`MEAN_REPORTING.md`](MEAN_REPORTING.md) | Mean bias profile on the up-set (superseded for reporting) | `results/mean_*` |
| [`MEAN_TABLE_FULLSPACE.md`](MEAN_TABLE_FULLSPACE.md) | `r_mu` on the full space, the table used in the paper | `results/mean_table_fullspace/` |
| [`EXPERIMENT_EPSILON_ABOVE_ONE.md`](EXPERIMENT_EPSILON_ABOVE_ONE.md) | Brief: the ε-radius at tolerances above 1 | |
| [`RESULT_EPSILON_ABOVE_ONE.md`](RESULT_EPSILON_ABOVE_ONE.md) | Result for that brief | `results/final_table_eps_gt1/` |

## Paper-writing briefs and plans

| File | Content |
|---|---|
| [`PAPER_NARRATIVE.md`](PAPER_NARRATIVE.md) | The argument, the results that serve it, the results to set aside |
| [`PAPER_INTEGRATION_FINAL_TABLE.md`](PAPER_INTEGRATION_FINAL_TABLE.md) | Real-network radius table at LLM-elicited knowledge |
| [`PAPER_INTEGRATION_R_EPSILON.md`](PAPER_INTEGRATION_R_EPSILON.md) | Session 8: the ε-bias radius |
| [`PAPER_INTEGRATION_MEAN_RADIUS.md`](PAPER_INTEGRATION_MEAN_RADIUS.md) | The mean bias radius as a companion to `r_val` |
| [`PAPER_NOTE_RE11.md`](PAPER_NOTE_RE11.md) | What [RE-11] changes in the paper |
| [`EVALUATION_PLAN.md`](EVALUATION_PLAN.md) | Evaluation plan, 11 September |
| [`REMAINING_EXPERIMENTS.md`](REMAINING_EXPERIMENTS.md) | Ledger of what still had to be measured, 12 September |
| [`ANONYMITY_CHECKLIST.md`](ANONYMITY_CHECKLIST.md) | Checklist for the double-blind release |
