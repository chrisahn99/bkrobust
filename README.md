# bkrobust

**Robustness of causal effect estimation to consistent-but-false background knowledge.**

Background knowledge injected into causal discovery — required and forbidden
edges, tiers, ancestral constraints — is checked only for *consistency*: Meek's
Algorithm 1 either terminates with an MPDAG or reports FAIL. It is never checked
for *truth*, and it cannot be, since checking truth would require the graph that
discovery was run to find. Knowledge that is **consistent but false** therefore
passes every check in current practice. The three-node chain is the whole
problem in miniature: `A - B - C` is equally the CPDAG of `A → B → C` and of
`C → B → A`, so an expert asserting the reverse of the truth is undetectably,
completely wrong. This repository builds the quantitative theory of what such
knowledge costs — two breakdown radii for the optimal adjustment set, a bias
bound past the second, an efficiency gap between them, and an audit of whether
amortized causal foundation models inherit the same failure.

> **Status: scaffold.** Every function raises `NotImplementedError`. Structure,
> interfaces, and documented intent are in place; no algorithm is implemented.

---

## The two radii

Fix a CPDAG `C`, a true DAG `G` in its equivalence class, a target pair
`(X, Y)`, and a distance on knowledge sets. Write `O*` for the
Henckel–Perković–Maathuis optimal adjustment set.

| Radius | Definition | What breaks | Detectable? |
|---|---|---|---|
| `δ_opt` | smallest perturbation at which `O*` stops being **optimal** | efficiency — unbiased, larger variance | **No** |
| `δ_valid` | smallest perturbation at which `O*` stops being **valid** | identification — biased at any `n` | **No** |

The conjecture is `δ_opt ≤ δ_valid`. The band between them is the interesting
regime: the analyst is wrong, is paying for it in effective sample size, and has
no diagnostic that would tell them so.

Full definitions and the three theorem targets: [docs/THEORY.md](docs/THEORY.md).

---

## Install

Python 3.11+. [`uv`](https://github.com/astral-sh/uv) recommended.

```bash
uv venv && source .venv/bin/activate
uv pip install -e ".[dev]"
```

Or with pip:

```bash
python -m pip install -e ".[dev]"
```

The heavy stack — torch, transformers, the foundation-model checkpoints — is a
separate extra, so the core theory and simulation code installs without a GPU
environment:

```bash
python -m pip install -e ".[cfm]"
```

Development hooks:

```bash
pre-commit install
```

One dependency is deliberately unresolved: the causal-discovery backend
(`causal-learn` vs `gcastle`). See the commented decision block in
`pyproject.toml`.

---

## Repository map

```
configs/          Hydra configs -- every experiment is config-driven
  experiment/     one per experiment; graph/ perturbation/ estimator/ model/
src/bkrobust/
  graphs/         MPDAG, Meek rules, consistency, adjustment, optimality, distances
  knowledge/      BackgroundKnowledge, misspecification taxonomy, perturbation
                  samplers, Meek cascade measurement, elicitation adapters
  theory/         the two radii, bias bound, efficiency gap, certificate
  estimation/     OLS adjustment, AIPW/DML, conformal wrapper
  cfm/            foundation-model registry, knowledge conditioning, adapters, audit
  representations/ constraint-poset embeddings and probes (conditional component)
  data/           synthetic and semi-synthetic generators, real loaders, leakage
  metrics/        effect metrics and structural metrics
  utils/          seeding, results IO, logging, shared plotting style
experiments/      thin Hydra entrypoints, one per experiment
tests/            canonical graphs and the contracts everything depends on
docs/             theory, experiment protocols, data and checkpoint provenance
paper/            references and generated figures
```

Two files carry more weight than the rest and are worth reading first:
`knowledge/perturb.py`, which documents the sampling contract every experiment
rests on, and `theory/certificate.py`, which is the practitioner-facing
deliverable.

---

## Running an experiment

Every entrypoint is a Hydra application. Nothing is hardcoded; if an experiment
needs a number, it comes from a config.

```bash
python experiments/run_synthetic_sweep.py
```

Override config groups or individual values:

```bash
python experiments/run_synthetic_sweep.py graph=scale_free perturbation=tier seed=7
```

Sweep:

```bash
python experiments/run_synthetic_sweep.py -m seed=1,2,3 perturbation=orientation_flip,ancestral
```

| Entrypoint | What it establishes |
|---|---|
| `run_synthetic_sweep.py` | both radii and their ordering (T1) |
| `run_cascade_study.py` | Meek amplification — the mechanism |
| `run_semisynthetic.py` | the radii survive realistic marginals |
| `run_real_data.py` | certificates on real datasets |
| `run_cfm_audit.py` | do foundation models inherit the failure? |
| `run_knowledge_embedding.py` | conditional component — may be cut |
| `make_figures.py` | regenerate every paper figure |

Results are written to `results/<experiment>/<timestamp>_<git-sha>/` with a
manifest recording the resolved config, git SHA, dirty-tree flag, and package
versions. Protocols, control arms, and reporting rules:
[docs/EXPERIMENTS.md](docs/EXPERIMENTS.md).

---

## Development

```bash
make test        # pytest
make lint        # ruff check + format --check
make typecheck   # mypy
make format      # ruff format + fix
make figures     # rebuild paper figures from the latest results
```

---

## Data and checkpoints

Nothing downloads automatically. Several datasets carry redistribution terms and
the model checkpoints have varying licences, so `scripts/download_data.sh` and
`scripts/download_checkpoints.sh` are deliberate stubs. Provenance, licences and
manual instructions: [docs/DATA.md](docs/DATA.md) and
[docs/CHECKPOINTS.md](docs/CHECKPOINTS.md).

---

## Anonymity

This repository is linked from a double-blind submission. Before making it
public, work through [docs/ANONYMITY_CHECKLIST.md](docs/ANONYMITY_CHECKLIST.md)
— including the item most easily missed, which is citing the authors' own prior
work in the third person like any other prior work.

```bash
make anonymity-check
```

---

## Citation

```bibtex
@inproceedings{anonymous2027bkrobust,
  title     = {TODO},
  author    = {Anonymous},
  booktitle = {Submitted to the International Conference on Learning Representations},
  year      = {2027},
  note      = {Under review}
}
```

## License

MIT. See [LICENSE](LICENSE).
