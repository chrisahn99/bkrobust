# E1 — Rebuild the ρ* instrument: rel50 @ |K|=8 vs sign @ |K|≤4

Run 2026-08-14, local Mac (8 cores, no GPU), `/Users/josecosta/miniforge3/bin/python3`
(3.12.7, numpy 2.3.3, networkx 3.5, scipy 1.16.3). Total compute ≈ 55 s.
Pre-registration: `PREREG.md`, written before any |K|=8 simulation was run.

**Verdict: SUPPORTED on the pre-registered rule — but the fix is not the re-analysis
the vault believes it is, and it inverts outside a narrow graph regime.**

---

## 0. Answer to task question 1: can rel50 @ |K|=8 be recomputed from data on disk?

**No.** `MAX_K = 4` (`code/run_linear.py:34`) is global — `sweep_density.py` imports
`analyse_one` from `run_linear`, so it inherits the cap. Measured on disk:

| file | records | |K| observed |
|---|---|---|
| `results/linear_raw.json` | 600 SCMs | generic 258×|K|=3, 321×|K|=4 |
| `results/density_sweep.json` | 1,650 SCMs | 3 or 4 only, all 11 cells |

`rel50` **is** free at |K|≤4 — `breakdown_radius(..., crit='rel50')` already exists at
`code/analyse.py:42,61`. Re-running it reproduces the published table verbatim
(§1). |K|=8 required **re-simulation**, reported as such throughout.

---

## 1. Baseline, from data on disk, zero simulation (`code/disk_reanalysis.py`)

Generic arm, n = 579 SCMs, |K| ∈ {3,4}, ball radius ρ ≤ 3. Bins are ρ* = 1 / 2 / 3 / censored.

| crit | 1 | 2 | 3 | censored | censored 95% CI | max bin | bins ≥5% |
|---|---|---|---|---|---|---|---|
| sign  | 0.057 | 0.028 | 0.012 | **0.903** | [0.877, 0.926] | 0.903 | 2 |
| rel50 | 0.168 | 0.066 | 0.012 | **0.755** | [0.718, 0.788] | 0.755 | 3 |
| ident | 0.510 | 0.162 | 0.023 | 0.306 | [0.268, 0.344] | 0.510 | 3 |

Counts: sign 33/16/7/523, rel50 97/38/7/437, ident 295/94/13/177. Reproduces
`results/summary.json` exactly. Tiered arm (n=180): sign 0.900, rel50 0.733, ident 0.239.

Ignorance-interval width / |τ|, generic arm: ρ=1 mean 0.351, median 0.000, **72.5% exactly
zero-width**; ρ=2 mean 0.480, 61.1% zero; ρ=3 mean 0.514, 59.6% zero.

**The 0.755 in row 2 is the pre-registered REFUTED threshold.** It is also the vault's
own registered falsifier threshold (75.5%) — this run identifies where that number
came from: it is the rel50 censoring rate at |K|≤4.

---

## 2. 🔴 The blocker the campaign never checked: |K|=8 is undefinable on the licensed ensemble

|K|=8 requires a CPDAG with ≥ 8 undirected edges. Measured on the pilot's **own kept
ensemble** (`linear_raw.json`, n=600), `n_undirected` distribution is
{3:258, 4:183, 5:91, 6:30, 7:29, 8:5, 9:3, 10:1}:

| |K| target | definable for | share |
|---|---|---|---|
| 4 | 342/600 | 57.0% |
| 5 | 159/600 | 26.5% |
| 6 | 68/600 | 11.3% |
| 7 | 38/600 | 6.3% |
| **8** | **9/600** | **1.5%** |

Density sweep (1,650 SCMs, p ∈ {6,8,10} × deg ∈ {1.5,2.5,4.0,6.0}): 127/1650 = **7.7%**.
Densest cell p8_deg6.0 tops out at 26%. A fresh p × deg scan (200 graphs/cell) shows
P(n_undirected ≥ 8) never exceeds 0.07 anywhere in p ≤ 10, and reaches 0.43 only at
p = 25, deg = 1.5 — because denser graphs manufacture v-structures, which *remove*
undirected edges.

**Consequence.** "Rebuild §3 on rel50-ρ* at |K|=8" is **not** "a re-plot of data already
in hand" (`rho-breakdown-adjustment.md:69-75`). It requires a new ensemble, because on
98.5% of the ensemble the paper is licensed over, |K|=8 does not exist.

---

## 3. The |K| sweep (`code/run_ksweep.py` → `code/analyse_ksweep.py`, `code/paired_test.py`)

**Design.** One uniform random ordering of the orientable statements per SCM, then
K4 = order[:4] ⊂ K6 = order[:6] ⊂ K8 = order[:8]. Identical graph, (x,y), Σ, τ across
arms. K4 remains a uniform random 4-subset, i.e. distributionally identical to the
original protocol, but now **paired** — this sidesteps the RNG-stream shift that makes a
naive `MAX_K` change produce a different ensemble (`run_linear.py:65-67`). All SCMs
rejection-sampled to `n_undirected ≥ 8`. ρ ≤ 3 throughout (same radius, bigger ball).

Three ensembles, 600 SCMs each:

| id | grid | note |
|---|---|---|
| `original` | p ∈ 5..8, deg ∈ {1.5,2,2.5} | **the pilot's exact grid** — scope-faithful |
| `licensed` | p ∈ 5..10, deg ∈ {1.5,2,2.5,4,6} | claim's full licensed box |
| `large` | p ∈ 15..25, deg ∈ {1.5,2,2.5} | outside scope; where |K|=8 is common |

### 3a. Head-to-head, rel50, scope-faithful `original` ensemble

PAIRED subset (all three arms amenable under true K), n = 137:

| arm | 1 | 2 | 3 | censored | max bin | bins ≥5% |
|---|---|---|---|---|---|---|
| K4 | 0.292 | 0.110 | 0.000 | 0.599 | 0.599 | 3 |
| K6 | 0.299 | 0.161 | 0.036 | 0.504 | 0.504 | 3 |
| **K8** | **0.292** | **0.168** | **0.080** | **0.460** | **0.460** | **4** |

Paired difference K8 − K4 = **−0.1387 [−0.2044, −0.0730]**; McNemar b=21, c=2, **p = 7e-05**.

UNPAIRED (every SCM where that arm is amenable under true K — the deployment view):

| arm | n | 1 | 2 | 3 | censored | censored 95% CI |
|---|---|---|---|---|---|---|
| K4 | 137 | 0.292 | 0.110 | 0.000 | 0.599 | [0.518, 0.679] |
| K6 | 329 | 0.313 | 0.195 | 0.052 | 0.441 | [0.386, 0.495] |
| **K8** | **524** | **0.424** | **0.170** | **0.071** | **0.336** | **[0.294, 0.374]** |

Widths / |τ|, K8 unpaired (n=524): ρ=1 mean **1.448 [0.981, 2.059]**, median 0.389,
q90 2.257, zero-width 33.4%; ρ=2 mean 2.114, zero 13.5%; ρ=3 mean 2.372, zero 10.7%.
(Compare |K|≤4 on disk: ρ=1 mean 0.351, zero-width 72.5%.)

`sign` at |K|=8 on the same ensemble is still bad: censored 0.744 [0.708, 0.781].
**The fix requires rel50 AND |K|=8 jointly; neither alone suffices.**

### 3b. `licensed` ensemble (p ≤ 10, deg ≤ 6)

rel50 UNPAIRED: K4 0.574 (n=61) → K6 0.341 (n=173) → **K8 0.271 [0.221, 0.320] (n=303)**,
bins 0.449/0.188/0.092/0.271, 4 bins ≥5%. PAIRED (n=61): 0.574 → 0.459, diff
−0.1148 [−0.2131, −0.0164], McNemar p = 0.065. Widths ρ=1 mean 2.814, zero 29.7%.

### 3c. 🔴 `large` ensemble (p 15–25) — the fix **inverts**

| crit | K4 (n=384) | K6 (n=474) | K8 (n=560) |
|---|---|---|---|
| sign censored | 0.990 | 0.975 | **0.970** [0.955, 0.982] |
| rel50 censored | 0.935 | 0.907 | **0.877** [0.848, 0.904] |
| ident censored | 0.542 | 0.477 | 0.443 |

rel50 @ K8 bins: 0.084/0.039/0.000/0.877 — **only 2 bins ≥5%, ρ*=3 is empty**. Widths
ρ=1 mean 0.186, **86.6% zero-width**. Paired diff K8−K4 = −0.0312 [−0.0495, −0.0156]:
statistically real, practically nil.

**0.877 ≥ 0.755 ⇒ REFUTED on this ensemble by the pre-registered rule.** Mechanism: O*
is a local functional, so on large graphs a randomly-drawn statement is almost never
adjacent to the (X,Y) neighbourhood and the perturbation does nothing. This is the *same*
locality result the paper wants to sell, acting as the instrument's own destroyer.

---

## 4. Decomposition — how much of the "fix" is actually |K|?

| step | rel50 censored | Δ |
|---|---|---|
| disk, unconditioned, |K| ≤ 4 (n=579) | 0.755 | — |
| conditioned on n_undirected ≥ 8, |K| = 4 (n=137) | 0.599 | **−15.6 pp (graph conditioning)** |
| conditioned, |K| = 8, paired (n=137) | 0.460 | **−13.9 pp (|K|, paired, p=7e-05)** |
| conditioned, |K| = 8, all K8-amenable (n=524) | 0.336 | −12.4 pp (amenability selection) |

**Roughly half the headline improvement is not |K| at all — it is restricting to graphs
with ≥ 8 undirected edges**, a 1.5% slice of the pilot's own ensemble. This must be
reported; presenting 0.336 against 0.755 as a |K| effect would be a confounded comparison.

## 5. Two side results worth carrying

- **`ident` gets WORSE with more knowledge, paired**: +0.124 [+0.073, +0.183] (original,
  McNemar p=2e-05), +0.312 [+0.197, +0.426] (licensed, p≈0). More true knowledge → more
  MPDAGs amenable → fewer identification failures reachable inside the ball. The
  unpaired trend runs the other way (0.029→0.111 original), so this is a pairing/selection
  effect and must not be quoted un-stratified.
- **Utility of knowledge, cleanly measured**: fraction of problems amenable under TRUE K
  rises 0.228 → 0.548 → 0.873 (original, |K|=4/6/8) and 0.102 → 0.288 → 0.505 (licensed).
  Independent of the instrument question and a genuine positive for §4.
- **`changed_still_valid = 0` replicates**: 0 across all three ensembles, all three |K|
  arms, 27,000+ perturbations. Still empirical, still not a theorem.

## 6. Faithfulness of the re-implementation (`code/validate_runner.py`)

Same runner, pilot grid, |K|=4, no n_undirected≥8 filter, 600 SCMs / 2,400 ρ=1 members.
On the pilot's own denominator (consistent members):

| quantity | this run | pilot published |
|---|---|---|
| Meek-inconsistent (of all) | 0.350 | 0.330 |
| loud non-amenable | 0.262 | 0.23 |
| O* untouched | 0.584 | 0.63 |
| **silent bias** | **0.153** | **0.141** |
| changed_still_valid | **0** | **0** |

Residual gap is the |K|=3 exclusion (pilot admits |K|=3; this arm does not) plus a
different seed stream. The machinery is the pilot's, imported unmodified.

---

## 7. Verdict against the pre-registration

| ensemble | rel50 @ |K|=8 censored | max bin | bins ≥5% | verdict |
|---|---|---|---|---|
| `original` unpaired (n=524) | **0.336** | 0.424 | 4 | **SUPPORTED** |
| `original` paired (n=137) | 0.460 | 0.460 | 4 | **SUPPORTED** |
| `licensed` unpaired (n=303) | 0.271 | 0.449 | 4 | **SUPPORTED** |
| `large` unpaired (n=560) | 0.877 | 0.877 | 2 | **REFUTED** (≥ 0.755) |

Rule (pre-registered): SUPPORTED iff censored ≤ 0.60 **and** max bin ≤ 0.60 **and** ≥3 of
4 bins ≥ 0.05; REFUTED iff censored ≥ 0.755.

**On the pilot's own graph regime the declared fix works and is better than the campaign's
own un-evidenced prose figure** (claimed 26/14/4/56 with mean width 1.0|τ|; measured
paired 29/17/8/46, unpaired 42/17/7/34, mean width 1.448|τ| [0.981, 2.059]). The
reviewer objection "your diagnostic returns no-breakdown nine problems in ten" is
answerable: 0.903 → 0.336, with all four bins populated.

**Three qualifications the campaign record does not contain**, each of which a reviewer
will find:
1. |K|=8 is definable for only **1.5%** of the licensed ensemble. §3 cannot be re-plotted;
   it must be re-simulated on a different graph distribution, and the claim's scope
   sentence changes with it.
2. Roughly **half** the improvement is the graph conditioning, not |K|.
3. On p ≥ 15 the instrument **inverts** (0.877 censored, 86.6% zero-width intervals),
   and the cause is O*-locality — the paper's own headline mechanism.

## 8. C13 compliance

The registry rule (a null must show its threshold was capable of producing a non-null) is
satisfied in both directions here: the same `breakdown_radius` returns 0.306 under
`ident` and 0.903 under `sign` on identical data (60 pp of measured dynamic range), and
the |K|=8 arm produces 0.336 on one ensemble and 0.877 on another. The instrument is not
a dead switch, and the `large`-ensemble REFUTED is a property of the estimand.

## 9. Reproduction

```
cd /Users/josecosta/mugango/output/2026-08-14_latent-causal-iclr2027-rerank/code
PY=/Users/josecosta/miniforge3/bin/python3
$PY disk_reanalysis.py                                            # ~2 s, no simulation
$PY run_ksweep.py 600 original ../results/ksweep_original.json    # ~9 s
$PY run_ksweep.py 600 licensed ../results/ksweep_licensed.json    # ~7 s
$PY run_ksweep.py 600 large    ../results/ksweep_large.json       # ~33 s
$PY analyse_ksweep.py ../results/ksweep_original.json original    # (+ licensed, large)
$PY paired_test.py                                                # ~3 s
$PY validate_runner.py                                            # ~3 s
```

Seeds fixed (`20260814`); `Pool(8)` with **ordered** `imap` so the kept set is the first
n passing seeds and is reproducible. Nothing under
`/Users/josecosta/research-pilots/latent-causal-rho-breakdown-knowledge` was modified;
`code/*.py` there was copied, and `graphs.py` / `adjust.py` / `scm.py` / `analyse.py` are
imported unmodified.
