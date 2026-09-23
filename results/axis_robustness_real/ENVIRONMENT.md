# Environment Verification Report

Repository: `<repo>`
Branch at time of verification: `experiments/todo_1`
Date: 2026-09-16

## 1. Interpreter and package versions

`/usr/bin/python3 -VV`:
```
Python 3.9.6 (default, Apr 17 2026, 18:15:52)
[Clang 21.0.0 (clang-2100.1.1.101)]
```

Platform string (`platform.platform()`):
```
macOS-26.6.2-arm64-arm-64bit
```

Package versions:
```
numpy       2.0.2
scipy       1.13.1
pandas      2.3.3
matplotlib  3.9.4
ortools     9.15.6755
```
`from ortools.sat.python import cp_model` succeeded (no exception).

## 2. "Tests that matter"

Command:
```
PYTHONPATH=src /usr/bin/python3 -m pytest tests/core tests/search tests/criterion tests/hybrid tests/demo -q
```

Result (verbatim tail):
```
........................................................................ [ 26%]
........................................................................ [ 53%]
........................................................................ [ 80%]
.....................................................                    [100%]
=============================== warnings summary ===============================
(14 matplotlib PyparsingDeprecationWarning warnings, unrelated to repo code)
-- Docs: https://docs.pytest.org/en/stable/how-to/capture-warnings.html
269 passed, 14 warnings in 46.01s
```

**Exact counts: 269 passed, 0 failed, 0 errors.** No failures to quote — all tests in `tests/core tests/search tests/criterion tests/hybrid tests/demo` passed.

## 3. Full suite baseline

### 3a. Collect-only

Command:
```
PYTHONPATH=src /usr/bin/python3 -m pytest tests -q --collect-only 2>&1 | tail -30
```

Result (verbatim tail):
```
=========================== short test summary info ============================
ERROR tests/test_consistency.py
ERROR tests/test_meek.py
ERROR tests/test_perturb.py
!!!!!!!!!!!!!!!!!!! Interrupted: 3 errors during collection !!!!!!!!!!!!!!!!!!!!
467 tests collected, 3 errors in 0.50s
```

**Exact counts: 467 tests collected, 3 collection errors.**

### 3b. Full run

Command:
```
PYTHONPATH=src /usr/bin/python3 -m pytest tests -q 2>&1 | tail -40
```

Result (verbatim tail):
```
=========================== short test summary info ============================
ERROR tests/test_consistency.py
ERROR tests/test_meek.py
ERROR tests/test_perturb.py
!!!!!!!!!!!!!!!!!!! Interrupted: 3 errors during collection !!!!!!!!!!!!!!!!!!!!
14 warnings, 3 errors in 0.51s
```

**Important factual note:** pytest's default behavior is to interrupt the entire session when there are collection errors (it does not fall back to `--continue-on-collection-errors`). This means `pytest tests -q` reports **0 passed, 0 failed, 3 errors** — none of the 464 collectible-and-runnable tests were actually executed in this invocation, because pytest aborted before running anything. This is a literal, faithful report of the exact command specified; the 269-passed count above (from the "tests that matter" subset, which does not touch the 3 broken files) is the only "ran and passed" evidence collected in this task. No `--continue-on-collection-errors` or other flag was added, per instructions to run the exact command given.

## 4. Confirmation that pre-existing failures are unrelated to the working tree

```
$ git status --porcelain
(no output — clean tree)
$ git stash list
(no output — empty stash)
```

Since `git status --porcelain` produced no output, the working tree is clean and already identical to HEAD. **No stashing is needed or was performed** — the 3 collection errors are therefore attributable to the committed state of the repository on this Python 3.9.6 interpreter, not to any uncommitted change. No `git stash` command was run, per instructions.

## 5. Test files that fail to collect, with exact ImportError

### `tests/test_consistency.py`
```
ImportError while importing test module '<repo>/tests/test_consistency.py'.
Hint: make sure your test modules/packages have valid Python names.
Traceback:
.../importlib/__init__.py:127: in import_module
    return _bootstrap._gcd_import(name[level:], package, level)
tests/test_consistency.py:13: in <module>
    from bkrobust.knowledge.base import BackgroundKnowledge
src/bkrobust/knowledge/base.py:29: in <module>
    from typing import TYPE_CHECKING, TypeAlias
E   ImportError: cannot import name 'TypeAlias' from 'typing' (/Library/Developer/CommandLineTools/Library/Frameworks/Python3.framework/Versions/3.9/lib/python3.9/typing.py)
```

### `tests/test_meek.py`
```
ImportError while importing test module '<repo>/tests/test_meek.py'.
Hint: make sure your test modules/packages have valid Python names.
Traceback:
.../importlib/__init__.py:127: in import_module
    return _bootstrap._gcd_import(name[level:], package, level)
tests/test_meek.py:12: in <module>
    from bkrobust.graphs import meek
src/bkrobust/graphs/meek.py:33: in <module>
    from bkrobust.graphs.mpdag import MPDAG
src/bkrobust/graphs/mpdag.py:24: in <module>
    from typing import TYPE_CHECKING, TypeAlias
E   ImportError: cannot import name 'TypeAlias' from 'typing' (/Library/Developer/CommandLineTools/Library/Frameworks/Python3.framework/Versions/3.9/lib/python3.9/typing.py)
```

### `tests/test_perturb.py`
```
ImportError while importing test module '<repo>/tests/test_perturb.py'.
Hint: make sure your test modules/packages have valid Python names.
Traceback:
.../importlib/__init__.py:127: in import_module
    return _bootstrap._gcd_import(name[level:], package, level)
tests/test_perturb.py:16: in <module>
    from bkrobust.knowledge.taxonomy import MisspecificationKind
src/bkrobust/knowledge/taxonomy.py:46: in <module>
    from enum import StrEnum
E   ImportError: cannot import name 'StrEnum' from 'enum' (/Library/Developer/CommandLineTools/Library/Frameworks/Python3.framework/Versions/3.9/lib/python3.9/enum.py)
```

All three ImportErrors match `docs/REMAINING_EXPERIMENTS.md`'s stated cause: `TypeAlias` (typing, needs 3.10+) and `StrEnum` (enum, needs 3.11+) are not available in Python 3.9.6. This matches the documented count of 3 files failing to collect. Note: the doc also states "23 tests fail"; that count was not separately reproduced here because the full run interrupts entirely on collection errors before executing any test (see section 3b) — with `PYTHONPATH=src /usr/bin/python3 -m pytest tests -q` as literally specified, the observed result is 0 ran / 3 collection errors, not 23 failed. This discrepancy is reported as a fact, not resolved.

## Raw log files saved alongside this report
- `full_suite_raw.txt` — full untruncated output of `pytest tests -q`
- `collect_only_raw.txt` — full untruncated output of `pytest tests -q --collect-only`
- `det_seed0.txt`, `det_seed12345.txt` — determinism check outputs
