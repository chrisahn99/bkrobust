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

> **Status: two trees, and the difference matters.** The package contains an
> *implemented* tree that produced every committed result, and an *unimplemented*
> scaffold tree in which every function still raises `NotImplementedError`. The
> repository map below marks which is which. Read
> [`src/bkrobust/core/conventions.py`](src/bkrobust/core/conventions.py) before
> quoting any radius.

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

The `δ` names above are the original framing. The **normative** definition that
all implemented code and every committed number obeys lives in
[`src/bkrobust/core/conventions.py`](src/bkrobust/core/conventions.py), which
writes them `r_val` and `r_opt`:

    r_P = min { d(G0, G) : G in the space, P fails at G }

with `d` the BFS hop count on the covering-relation neighbour graph. Shells
`0..r-1` are certified clean, so the practitioner-facing "how many claims may I
get wrong" number is **`r - 1`**, not `r`. `UNREACHED = -1` is a status, not a
radius: never average it or plot it on a numeric axis.

Full definitions and the three theorem targets: [docs/THEORY.md](docs/THEORY.md);
proved and verified statements, with their assumption chains, in
[THEOREMS.md](THEOREMS.md).

### Two senses of "two radii" -- do not conflate them

The pair above (`r_val` / `r_opt`) is indexed by **which property breaks**:
validity or optimality. A second, independent pair is indexed by **what a step
costs**:

| Pair | Axis | Members |
|---|---|---|
| property | what breaks | `r_val`, `r_opt`, `r_eps` |
| unit | what one move is | `r_hop` (a covering step in MPDAG space), `r_claim` (one revision of one asserted sentence) |

`r_val` in the code and in every committed result is a **`r_hop`-unit** radius:
its move set is the orientations of `G0`, i.e. the Meek closure of `K`, not the
sentences the analyst uttered. The two units coincide exactly when the analyst's
knowledge is already Meek-closed (`|K| == k_g0`). They can differ by more than an
order of magnitude when it is not. See
[docs/PAPER_NARRATIVE.md](docs/PAPER_NARRATIVE.md) section 4.

### Writing the paper

Start with [docs/PAPER_NARRATIVE.md](docs/PAPER_NARRATIVE.md) -- the argument,
the results that serve it, and the results to set aside. What is still unmeasured
is in [docs/REMAINING_EXPERIMENTS.md](docs/REMAINING_EXPERIMENTS.md). The
session reports at the root are sources for numbers, not drafts: several of their
headlines are superseded.

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

**Implemented.** Everything below runs, is tested, and produced the committed
results.

```
src/bkrobust/
  core/           THE RADIUS CONVENTION (conventions.py), space construction and
                  BFS distance (spacelib.py), enumeration oracle, results IO
  demo/           the live MPDAG type (graph.py), Meek rules R1-R4 and DAG
                  extensions (meek.py), brute-force space (space.py), adjustment
                  sets and linear-SEM bias (evaluate.py), CPDAG and knowledge
                  construction (example.py)
  hybrid.py       THE ENTRY POINT -- breakdown_radius(cpdag, K, x, y, z)
  search/         exact upward BFS (exact.py, frozen reference), the Lemma O
                  acceleration (exact_fast.py), differential test harness
  sat/            CP-SAT encodings E1/E2/E3; e1.py is on the default path
  gac/            generalized adjustment criterion, DAG and MPDAG level
  mpdag_criterion/ polynomial MPDAG-level validity via edge-state reachability
  synth/          instance generators, knowledge corruptions (omit/flip/
                  compound/tiered), sweep runners
  benchmarks/     benchmark-network acquisition, parsers, and the real-network
                  measurement driver (measure.py)
  analysis/       per-session figures, PDF builders, and the verifiers that
                  re-assert every reported number against committed files
results/          axisa*/ and axisb*/, each with a manifest recording git SHA,
                  environment, solver seed, and a determinism block
tests/            canonical graphs and the contracts everything depends on
docs/             theory, experiment protocols, data and checkpoint provenance
paper/            references and generated figures
```

**Scaffold — not implemented.** Every function in these packages raises
`NotImplementedError`. Their docstrings are the clearest prose statements of the
definitions in the project and are worth reading; the code is not callable, and
none of it is on the path of any committed result.

```
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
configs/          Hydra configs for the scaffold entrypoints
experiments/      thin Hydra entrypoints -- these drive the scaffold tree
```

Note the collisions: `graphs/mpdag.py` and `demo/graph.py` both define an
`MPDAG`, and `graphs/meek.py` and `demo/meek.py` both define Meek's rules. The
implemented code uses the `demo/` versions throughout.

Worth reading first, in order: `core/conventions.py` for what a radius means,
`hybrid.py` for how one is computed, `gac/mpdag_level.py` for what "fails"
means, and `search/differential.py` for the harness that every accelerated
method has to pass. In the scaffold tree, `theory/radius.py` and
`graphs/distances.py` carry the definitional docstrings.

---

## Running an experiment

> **These entrypoints drive the scaffold tree and therefore do not run yet.**
> They are the planned interface, kept because the protocol they encode is still
> the intended one. The results committed under `results/` were produced by the
> implemented tree instead -- see "Reproducing the committed results" below.

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

## Worked example: the breakdown radius

A fully enumerated, exactly computed demonstration of the breakdown radius on a
small clinical example lives in [`report.md`](report.md). It walks through what
one unit of perturbation is, builds the whole perturbation space layer by layer,
and reports where an adjustment set first fails and what that means for a
practitioner.

```bash
PYTHONPATH=src python3 -m bkrobust.demo.run_all
```

Writes `results/breakdown_radius_demo/` and 20 figures to `figures/`; runtime
about 7 seconds. Code in `src/bkrobust/demo/`, tests in `tests/demo/`.

Note that the demonstration runs on Python 3.9 as well as the 3.11+ this package
targets, because it was written against the interpreter available at the time.

---

## Reproducing the committed results

The implemented tree is imported as a library rather than driven by Hydra. The
single documented entry point for the metric is:

```python
from bkrobust.hybrid import breakdown_radius, describe

result = breakdown_radius(cpdag, knowledge, x, y, z)   # never sees the true DAG
print(result.radius, result.method, result.assumes)
```

`result.assumes` records the assumption chain the number inherits: every upward
search is exact **iff Conjecture 2 holds**, which rests on Anti-Exchange Case B
(verified, not proved -- [THEOREMS.md](THEOREMS.md) sections 4c and 6). The error
is one-sided, so a radius can be too large but never too small. Do not drop that
string when a number is copied into a table.

Two sweep drivers sit above it: `benchmarks/measure.py` for the benchmark
networks and `synth/runner.py` for generated instances. **They apply different
admissibility gates, and the gates are not equivalent** --
`benchmarks.measure.fast_gate` is strictly stricter (16 disagreements in 46,800
cases, all one direction). Always record which one produced a number.

Each `results/` subtree carries a `manifest.json` with the git SHA, environment,
solver parameters and a determinism block, and the session reports are checked by
`analysis/session*_verify.py`, which re-asserts every quoted number against the
committed files.

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
