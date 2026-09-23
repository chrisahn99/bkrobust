# bkrobust -- development tasks.
# Requires `uv` (preferred) or a Python 3.11+ venv with pip.

PYTHON ?= python3
UV     := $(shell command -v uv 2>/dev/null)

.PHONY: help install test lint format typecheck clean figures hooks \
        re11-sweep re11-analyse re11-figures re11-verify re11-all

help:  ## Show this help.
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) \
		| awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-18s\033[0m %s\n", $$1, $$2}'

install:  ## Install the package with the dev extra (editable).
ifneq ($(UV),)
	uv pip install -e ".[dev]"
else
	$(PYTHON) -m pip install -e ".[dev]"
endif

hooks:  ## Install pre-commit git hooks.
	pre-commit install

test:  ## Run the test suite.
	pytest

lint:  ## Lint without modifying files.
	ruff check .
	ruff format --check .

format:  ## Auto-format and auto-fix.
	ruff format .
	ruff check --fix .

typecheck:  ## Static type check.
	mypy

figures:  ## Regenerate the paper's figures 1, 2 and 4 into figures/paper/.
	PYTHONHASHSEED=0 PYTHONPATH=src $(PYTHON) scripts/make_paper_figures.py

# --- [RE-11]: survival and the paired cross-arm on the 831 real rows ---------
# Everything below runs on /usr/bin/python3 (3.9.6), which is the interpreter
# that produced the committed corpus. Each target is idempotent.

RE11_PY ?= PYTHONPATH=src /usr/bin/python3

re11-sweep:  ## [RE-11] Run or resume the 725-shard sweep. Safe to re-run; skips completed shards.
	$(RE11_PY) -m bkrobust.robustness.real_frame
	$(RE11_PY) -m bkrobust.robustness.run_real_survival pool --workers 8 --n-draws 1000
	$(RE11_PY) -m bkrobust.robustness.run_real_survival status

re11-analyse:  ## [RE-11] Rebuild every derived table and the anchor table from the committed shards.
	$(RE11_PY) -m bkrobust.robustness.real_analyse --n-boot 10000
	$(RE11_PY) -m bkrobust.robustness.real_within_between
	$(RE11_PY) -m bkrobust.robustness.real_loo_bootstrap
	$(RE11_PY) -m bkrobust.robustness.real_knowledge_model
	$(RE11_PY) -m bkrobust.robustness.real_matched_coverage
	$(RE11_PY) -m bkrobust.robustness.real_spotlight
	$(RE11_PY) -m bkrobust.robustness.real_precision
	$(RE11_PY) -m bkrobust.robustness.build_table_tau_real

re11-figures:  ## [RE-11] Regenerate every session-9 figure (PDF + PNG) from committed results.
	$(RE11_PY) -m bkrobust.analysis.session9_figures

re11-verify:  ## [RE-11] Re-assert every number in report_real_survival.md against the committed files.
	$(RE11_PY) -m bkrobust.analysis.session9_verify

re11-all: re11-analyse re11-figures re11-verify  ## [RE-11] Everything downstream of the sweep.

clean:  ## Remove caches and build artefacts (leaves results/ and data/ alone).
	rm -rf build dist *.egg-info .pytest_cache .mypy_cache .ruff_cache .coverage htmlcov
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
	find . -type d -name .ipynb_checkpoints -prune -exec rm -rf {} +
