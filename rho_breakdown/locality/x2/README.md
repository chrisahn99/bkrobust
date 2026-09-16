# X2 — locality: can a knowledge statement far from the query move `O*`?

Ran 2026-08-19 on **betelgeuse** (16 cores, 12 workers, `~/e1venv/bin/python`, no GPU).
Work dir: `~/latent-causal/x2-locality/{code,results,logs}`.

**Read in this order:** `PREREG.md` (frozen before implementation) → `AUDIT.md` (adversarial read of
that design, verdict FLAWED, 19 binding mitigations) → `RESULTS.md`.

## The one-line answer

The falsifier fired. **S5 as printed is refuted** (70 knowledge sets whose *every* misstatement lies
at distance ≥ 1 change `O*`, all 70 silently biasing, minimal `p = 5`, exhibited with the Meek
chain). **S8 is not supported at any radius.** **S9 is refuted**, `r*_elicit ≥ 4`, by far knowledge
*creating* identification rather than destroying it. Only **S5 restricted to `ρ = 1`** survives — and
there it is far stronger than the campaign ever showed, being *exhaustively* true over every CPDAG,
every knowledge state and every query at `p ≤ 5`.

The apparent E1′ contradiction (0.935 censoring at distance 1) is not a refutation: 37 of the 38
non-censored SCMs are loud identification failure folded in by `ρ*_any`, 1 is a genuine ρ ≥ 2
`O*`-change.

## Artefacts

| file | what it is |
|---|---|
| `RESULTS.md` | verdicts on S5 / S8 / S9 separately, every fraction with five denominators |
| `code/test_machinery.py` | **17 brute-force machinery tests, run before any experiment; 17/17 pass** |
| `code/x2lib.py` | the only new machinery: `inf`-safe BFS distance, one-pass `O*`, the six-event taxonomy, ARM D's prune partition |
| `code/step0_archive.py` | reproduces E1′'s locality table from raw and decomposes the contradiction |
| `code/verify_armB.py` | gate G3: re-derives every archived member independently (**0 mismatches / 176,592**) |
| `code/run_scan.py` | ARM A (mechanism / `E_gain`) and ARM E (the primary S5 arm after M3) |
| `code/run_lemma.py` | the exhaustive small-`p` existence search — the proof-substitute |
| `code/exhibit.py` | the counterexample dossiers (PREREG 4.5) with traced Meek rule chains |
| `code/analyse_x2.py` | all tables, gates and pre-registered verdicts |
| `results/` | analyses and aggregates. The raw per-SCM scan JSONs (~100 MB) and `results/pass1/` stay on betelgeuse at `~/latent-causal/x2-locality/results/`, as with the E1′ mirror |
| `logs/` | run logs |

`graphs.py`, `adjust.py`, `scm.py`, `se.py` are **copied unmodified** from
`~/latent-causal/e1prime-se/code/`. Nothing under `e1prime-se/` or `rho-breakdown-knowledge/` was
touched. No git commit was made by any agent.
