# Docs

Working notes behind the paper. Proofs as submitted are in the paper's
appendix; [`../THEOREMS.md`](../THEOREMS.md) and
[`R_EPSILON_THEORY.md`](R_EPSILON_THEORY.md) are the working versions, with the
status of each statement and the computational checks behind it.

## Reference

| File | Content |
|---|---|
| [`R_EPSILON_THEORY.md`](R_EPSILON_THEORY.md) | The ε-bias radius: statements, proofs, and what each one assumes |
| [`THEORY.md`](THEORY.md) | Original definitions and theorem *targets*, written before any result; see `THEOREMS.md` for what was proved |
| [`EXPERIMENTS.md`](EXPERIMENTS.md) | Planned protocols for the Hydra entry points of the scaffold tree |
| [`DATA.md`](DATA.md) | Dataset provenance and licences |
| [`CHECKPOINTS.md`](CHECKPOINTS.md) | Causal-foundation-model checkpoint sources (scaffold tree) |

## Result notes: the bias radius and its mean companion

| File | Content | Results |
|---|---|---|
| [`EXPERIMENT_MEAN_VS_MAX.md`](EXPERIMENT_MEAN_VS_MAX.md) | Mean vs. max bias per shell | `results/mean_vs_max/` |
| [`RESULT_MEAN_SOUNDNESS_GATE.md`](RESULT_MEAN_SOUNDNESS_GATE.md) | The mean over the retraction up-set is unsound | `results/mean_fullspace/` |
| [`MEAN_REPORTING.md`](MEAN_REPORTING.md) | Mean bias profile on the up-set (superseded for reporting) | `results/mean_*` |
| [`MEAN_TABLE_FULLSPACE.md`](MEAN_TABLE_FULLSPACE.md) | `r_mu` on the full space, as reported in the paper | `results/mean_table_fullspace/` |
| [`EXPERIMENT_EPSILON_ABOVE_ONE.md`](EXPERIMENT_EPSILON_ABOVE_ONE.md) | Protocol: the ε-radius at tolerances above 1 | |
| [`RESULT_EPSILON_ABOVE_ONE.md`](RESULT_EPSILON_ABOVE_ONE.md) | Result of that protocol | `results/final_table_eps_gt1/` |
