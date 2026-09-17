# Environment, 12 September

The machine every Stage 0 number was produced on, and the two things that had
to be fixed before anything would run.

## Interpreter and packages

| | |
|---|---|
| python | 3.12.7 (`uv venv --python 3.12 .venv`) |
| ortools | 9.15.6755 |
| numpy | 2.5.3 |
| networkx | 3.6.1 |
| platform | macOS 26.5.2, arm64 |

Installed with `uv pip install -e ".[dev]"` followed by `uv pip install "ortools>=9.15"`.

## Two things that do not install on their own

**`ortools` is not declared anywhere.** It is imported at module level by
`hybrid.py`, so `from bkrobust.hybrid import breakdown_radius` fails on a clean
install, and every radius in the project goes through that import. It is
absent from `dependencies` in `pyproject.toml`, from the optional extras, and
from `requirements.txt`. Pinning it resolves to the same 9.15.6755 the axisa3
manifest records, so nothing about the committed results is in question; the
install instructions are simply incomplete. It belongs in `dependencies`.

**Python 3.9 cannot collect the suite.** `pyproject.toml` asks for 3.11 or
later and the system interpreter here is 3.9.6, which is what produced the
"23 failures and 3 files that will not collect" recorded in `NEXT.md`. On 3.12
the three files collect, and the failure count rises to 50 for that reason
rather than a new defect.

## Test suite

`python -m pytest -q`: **50 failed, 443 passed** in 95 s.

Forty-nine of the fifty are `NotImplementedError` raised by the scaffold stubs
at the top level of `tests/`, in six files that test an interface nobody has
built: `test_perturb` 10, `test_meek` 9, `test_adjustment` 9, `test_radius` 8,
`test_consistency` 7, `test_distances` 6. The modules the project actually uses
live under `tests/demo`, `tests/gac`, `tests/search`, `tests/core`,
`tests/criterion`, `tests/hybrid`, `tests/sat` and `tests/benchmarks`, and they
pass.

The fiftieth is real but is not about this code:
`tests/demo/test_evaluate.py::test_dseparation_matches_networkx_on_random_dags`
calls `networkx.algorithms.d_separated`, removed in networkx 3.5 and replaced
by `is_d_separator`. The test's cross-check against networkx cannot run on a
current networkx; the d-separation implementation under test is untouched.

`NEXT.md` says a suite in this state will eventually be mistaken for a clean
signal. Until the stubs are filled or deselected, the signal to read is the
subdirectories.

## Reproduction gate

Nothing downstream was run until the committed radii came back from the
committed inputs. `experiments/stage0_reproduce.py` rebuilds the inputs the way
`benchmarks.measure.evaluate` does, calls `hybrid.breakdown_radius`, and
compares the radius, the dispatch leg, the exactness flag and the adjustment
set size against the committed row.

**30 rows, stratified across both dispatch legs, 3.7 s, 0 disagreements.**

`pathfinder` is excluded from every Stage 0 run: an 85-node chain component,
three rows, all three censored in the committed sweep.

## Working tree

Twenty figure files and `results/axisa3/networks/acquisition_manifest.json`
carry uncommitted modifications that predate this work; the manifest one is a
local absolute path replacing the path on the machine that built it. They are
left alone. The 8 MB `pgmpy-1.0.0.tar.gz` under `results/axisa3/networks/` is
untracked and stays untracked.
