# A Diagnostic Radius for Knowledge-Induced Causal Effect Estimation

Code and committed results for the ICLR 2027 submission of the same title.

Observational causal discovery returns an equivalence class (a CPDAG). An
analyst orients the remaining edges with background knowledge, Meek's rules
propagate it, and a covariate adjustment set `Z` is read off the resulting
MPDAG `G0`. The **breakdown radius** `r_val` is the fewest orientation claims
of `G0` that, if wrong, can make `Z` an invalid adjustment set. It is computed
from the analyst's own inputs only, never from the true graph.

```python
from bkrobust.demo.graph import MPDAG
from bkrobust.hybrid import breakdown_radius, describe

# Estimated CPDAG: the triangle C - X - Y - C is fully undirected.
cpdag = MPDAG(nodes="CXY", undirected=[("C", "X"), ("C", "Y"), ("X", "Y")])
# The analyst asserts C -> X and X -> Y; Meek closure does the rest.
knowledge = [("C", "X"), ("X", "Y")]

result = breakdown_radius(cpdag, knowledge, x="X", y="Y", z=frozenset({"C"}))
print(result.radius)       # 1
print(describe(result))
```

`result.radius` is `UNREACHED = -1` when no reachable knowledge state
invalidates `Z`. That value is a status, not a number: never average it.

---

## Install

```bash
python -m venv .venv
.venv/bin/pip install -e ".[dev]"
```

Every command below runs from the repository root with `PYTHONPATH=src`.
Nothing downloads at run time: the benchmark networks, the elicited knowledge
and every sweep output are committed under `results/`.

The committed results were produced by two interpreters, and each sweep's
`manifest.json` records which one, with package versions:

- Python 3.9.6 (numpy 2.0, networkx 3.2, pandas 2.3) for the survival sweeps
  (`results/axis_robustness*`) and the real-graph corpus (`results/axisa3`);
- Python 3.14 (numpy 2.5, networkx 3.6, OR-Tools 9.15) for the end-to-end
  timing sweep, the bias-radius work and the real-network tables.

Sweeps are bit-identical across `PYTHONHASHSEED` values, and the CP-SAT solver
runs single-threaded with a fixed seed. Manifests record the git SHA of the
development history, which is not part of this anonymised snapshot.

## Tests

```bash
PYTHONPATH=src .venv/bin/python -m pytest
```

Six test files (`tests/test_{adjustment,consistency,distances,meek,perturb,radius}.py`)
target the unimplemented scaffold packages listed under *Repository map* and
fail with `NotImplementedError` by design; every other test passes.

---

## Reproducing the paper

Each row names the driver that wrote the committed files. The drivers' module
docstrings give the full set of flags. Tables and figures are rebuilt from the
committed files without re-running any sweep.

### Figures

```bash
PYTHONHASHSEED=0 PYTHONPATH=src .venv/bin/python scripts/make_paper_figures.py
```

writes `figures/paper/fig_running_example.pdf` (Figure 1),
`fig_survival_ranking.pdf` (Figure 2, with its table of τ_b values as
`fig_survival_ranking_tau.csv`) and `fig_e2e_speedup.pdf` (Figure 4), in
about ten seconds. `PYTHONHASHSEED=0` only fixes the node layout of Figure 1.

The remaining figures are in `figures/`:

| Paper | File | Rebuilt by |
|---|---|---|
| Fig. 5 (bias staircase) | `fige1_staircase.pdf` | `python -m bkrobust.epsilon.figures` |
| Fig. 6 (what ε buys) | `fige2_what_epsilon_buys.pdf` | `python -m bkrobust.epsilon.figures` |
| Fig. 7 (separation on real networks) | `s6_f3_separation.pdf` | `python -m bkrobust.analysis.session6_figures` |

Figures 5 and 6 are rebuilt identically under Python 3.14; Figure 7
reproduces byte for byte under the Python 3.9 / Matplotlib 3.9 stack that
drew it, and differs only in rendering under newer Matplotlib. The three
`mean_*.png` files illustrate notes in `docs/`.

### Section 5 and Appendix D: computation

| Claim | Results | Driver |
|---|---|---|
| End-to-end timing on 256 instances (Fig. 4, Table 5) | `results/e2e_speedup_gated/` | `experiments/e2e_speedup.py --gate --out results/e2e_speedup_gated` (serial, ~80 min); plots: `experiments/plot_e2e.py` |
| Encodings E1/E2/E3 and their differential validation | `results/axisb3/`, `results/axisb4/`, `results/search/` | `bkrobust.sat`, `bkrobust.search.differential` |
| Anti-exchange property, exhaustive check (Table 4) | `results/epsilon/chickering/` | `python -m bkrobust.epsilon.chickering` |

### Section 6 and Appendix G: does the radius track survival?

| Claim | Results | Driver |
|---|---|---|
| Synthetic survival, 4,104 instances, 32.8M draws | `results/axis_robustness_p6/` | `python -m bkrobust.robustness.run_survival_p6` (nine shards, then `--merge`; see `RUN_NOTES.md` there) |
| Paired τ_b differences (Table 10) | `results/axis_robustness_p6_paired/` | `experiments/paired_synth_p6.py` |
| Knowledge elicited from six LLMs, ten conditions | `results/elicit/knowledge.json` | `experiments/questionnaire_build.py`, `experiments/elicit_run.py`, `experiments/remote/` (vLLM on a GPU node) |
| Real-network query frame (831 rows, 25 networks) | `results/axis_robustness_real/frame.jsonl`, `results/frame/` | `python -m bkrobust.robustness.real_frame`, `experiments/frame_build.py` |
| Real-structure survival with LLM knowledge (Tables 1, 11–14) | `results/axis_robustness_llm/` | `scripts/run_llm_survival_panel.sh` (resumable, skips finished shards) |
| Real-network corpus | `results/axisa3/` | `bkrobust.benchmarks`; parsed from one checksummed archive |
| Claims-vs-hops units of the corruption axis | `results/axis_robustness/hops_*` | `python -m bkrobust.robustness.run_hops` |

The elicitation step needs model weights and a GPU. Everything downstream of
it reads the committed `results/elicit/knowledge.json`, so no sweep needs
network access or a model API.

### Section 7 and Appendix E: the radius on the effect scale

| Claim | Results | Driver |
|---|---|---|
| Real-network radius panel, `r_val` and `r_eps` (Table 6) | `results/final_table/` | `experiments/final_table.py --queries-per-network 5 --time-limit 120 --out results/final_table` (~75 min); `--report-only` re-renders the table |
| Tolerances above 1 | `results/final_table_eps_gt1/` | `experiments/final_table.py` (see `docs/RESULT_EPSILON_ABOVE_ONE.md`) |
| Mean bias radius `r_mu` (Tables 2 and 7) | `results/mean_table_fullspace/` | `experiments/mean_table_fullspace.py --budget 900 --max-undirected 10 --out results/mean_table_fullspace` (~35 min) |
| Theorem verification | `results/epsilon/verification/` | `python -m bkrobust.epsilon.verify --target 400 --sem-draws 3` |
| Counterexample search | `results/epsilon/counterexamples/` | `python -m bkrobust.epsilon.counterexamples` |
| Synthetic study, 900 instances | `results/epsilon/study/` | `experiments/run_r_epsilon.py` |
| Finite-sample check | `results/epsilon/finite_sample/` | `experiments/run_r_epsilon_finite_sample.py` |
| Worked example (Tables 8 and 9) | `results/breakdown_radius_demo/`, `results/epsilon/worked/` | `python -m bkrobust.demo.run_all`, `experiments/run_r_epsilon_worked.py` |

Proofs of the order-theoretic results are in the paper's appendix;
[`THEOREMS.md`](THEOREMS.md) and [`docs/R_EPSILON_THEORY.md`](docs/R_EPSILON_THEORY.md)
are the working versions, with the status of every statement and the
computational checks behind them.

`results/` holds a few further sweeps (`axisa2`, `axisb2`, `synth`, `stage0`,
`mean_*`, `matched_control`, `zero_radius`, …) that support statements in the
appendix or were used to design the experiments above. See
[`results/README.md`](results/README.md).

---

## Repository map

```
src/bkrobust/
  hybrid.py          THE ENTRY POINT: breakdown_radius(cpdag, K, x, y, z)
  core/              the radius convention (conventions.py), space construction,
                     BFS distance, enumeration oracle, results IO
  demo/              MPDAG type, Meek rules R1-R4, DAG extensions, adjustment
                     sets and linear-SEM bias; the running example
  search/            exact upward (retraction) search, its accelerated variant,
                     and the differential test harness
  sat/               CP-SAT encodings E1/E2/E3; E1 is on the default path
  gac/, mpdag_criterion/
                     generalized adjustment criterion, DAG and MPDAG level
  synth/             instance generators and knowledge corruptions
  benchmarks/        benchmark-network acquisition, parsers, measurement driver
  robustness/        survival under corrupted knowledge, baselines, bootstraps
  epsilon/           bias functional B, eps-bias radius, mean profile, checks
  analysis/          figure scripts and a verifier for the real-structure run
experiments/         sweep drivers for the results above
scripts/             make_paper_figures.py, run_llm_survival_panel.sh
results/             one directory per sweep, each with a manifest
figures/             paper figures
tests/               tests for the implemented tree
docs/                theory notes, data provenance, result notes
THEOREMS.md          formal statements, proofs and their status
```

**Scaffold, not implemented.** `graphs/`, `knowledge/`, `theory/`,
`estimation/`, `cfm/`, `representations/`, `data/`, `metrics/` and `utils/`
under `src/bkrobust/`, together with `configs/`, `notebooks/` and the Hydra
entry points `experiments/run_{synthetic_sweep,cascade_study,semisynthetic,real_data,cfm_audit,knowledge_embedding}.py`
and `experiments/make_figures.py`, are an interface that was planned and not
built. Every function there raises `NotImplementedError`, and none of it is on
the path of any result. `graphs/mpdag.py` and `demo/graph.py` both define an
`MPDAG`; the implemented code uses `demo/` throughout.

## License

MIT. See [LICENSE](LICENSE).
