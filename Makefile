# bkrobust -- development tasks.
# Requires `uv` (preferred) or a Python 3.11+ venv with pip.

PYTHON ?= python3
UV     := $(shell command -v uv 2>/dev/null)

.PHONY: help install test lint format typecheck clean figures hooks anonymity-check

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

figures:  ## Regenerate every paper figure from the latest results.
	$(PYTHON) experiments/make_figures.py

anonymity-check:  ## Grep tracked files for de-anonymising strings.
	@echo "== git author strings =="
	@git log --format='%an <%ae>' | sort -u
	@echo "== absolute home/cluster paths in tracked files =="
	@git grep -nE '/(home|Users|scratch|nfs|mnt)/[A-Za-z0-9._-]+' -- . || echo "  none"
	@echo "== notebooks with output cells =="
	@git grep -lE '"(outputs|execution_count)": *\[?[^]n]' -- '*.ipynb' || echo "  none"
	@echo
	@echo "Work through docs/ANONYMITY_CHECKLIST.md -- these greps are a"
	@echo "reminder, not a substitute for it."

clean:  ## Remove caches and build artefacts (leaves results/ and data/ alone).
	rm -rf build dist *.egg-info .pytest_cache .mypy_cache .ruff_cache .coverage htmlcov
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
	find . -type d -name .ipynb_checkpoints -prune -exec rm -rf {} +
