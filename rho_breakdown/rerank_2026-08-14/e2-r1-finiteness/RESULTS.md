# E2 — R1-closed finiteness: the empirical half

**Run 2026-08-14, local Mac, 8 cores, ~6 s of compute. Verdict: SUPPORTED** on the
pre-registered rule in [`PREREGISTRATION.md`](PREREGISTRATION.md) — **against the prior
recorded in that document**, which expected REFUTED.

> This experiment does **not** prove the R1-closed finiteness theorem. It runs the
> empirical half the pilot never ran, and its result is that the theorem's *premise*
> survives the perturbation it must range over, and that the cascade *is* contained —
> by a factor of ~2 in the tail, not by an order of magnitude. The theorem is worth
> attempting. It is not free.

---

## 0. Commands

```bash
cd output/2026-08-14_latent-causal-iclr2027-rerank/e2-r1-finiteness/code
python test_meek_rules.py                        # machinery validation, ~2 s
python run_e2.py ../results/e2_raw.json          # 600 SCMs, 3 arms, ~6 s on 8 cores
python analyse_e2.py  ../results/e2_raw.json     # PRE-REGISTERED gates
python analyse_e2b.py ../results/e2_raw.json     # post-hoc attacks A-D
python analyse_e2c.py ../results/e2_raw.json     # post-hoc attack E (decisive)
python make_fig.py
```
Interpreter `/Users/josecosta/miniforge3/bin/python3` (3.12.7, numpy 2.3.3, networkx 3.5,
scipy 1.16.3). No GPU, no betelgeuse, no network.

## 1. What was built (it did not exist)

| piece | file | why |
|---|---|---|
| rule-selectable Meek closure + per-rule counter | `code/meek_rules.py` | `graphs.py:111` hard-codes R1→R2→R3→R4 behind one boolean; it can size a cascade but never attribute it |
| `apply_background_knowledge_rules(C,K,rules)` | `code/meek_rules.py` | Perković UAI'17 Alg. 1 with a rule filter |
| tier-perturbation operator + tiering ball | `code/run_e2.py:96` `move_ball` | absent from the pilot entirely |
| tier-space distance | `n_changed = \|K Δ K'\|` | absent; this was the load-bearing design decision |
| non-negative cascade metric | `cascade_size` | the pilot's `\|dir(G)\|−\|dir(C)\|−\|K\|` goes negative when a statement of a *full* tiered `K` names an edge an earlier closure already oriented |

**Two pilot defects fixed, both load-bearing.**
`run_linear.py:55` (`if D[u,v]==1`) silently drops every statement disagreeing with the
true DAG — a no-op today, but under a perturbed tiering it deletes exactly the signal.
**Removed.** `MAX_K = 4` subsampling of tiered `K` — Bang & Didelez's property is about the
**full** cross-tier set. **Removed**; the generic arm is matched to `|K_T|` per SCM instead.

**Machinery validation** (`results/test_machinery_out.txt`), against the pilot's own
already-validated code, not a re-implementation:

```
T-A closure equality   : 1166 graphs        OK      (rules=ALL == graphs.meek_closure)
T-B alg-1 equality     : 3498 knowledge sets OK     (680 matched MeekFails)
T-C R1 subset-of-full  : 2818 MPDAGs         OK     (R1-only never over-orients)
T-D rule firings       : {'R1': 1669, 'R2': 713, 'R4': 169}
```
⚠️ R3 is credited 0 times because the counter is **first-rule-wins** (it mirrors the
pilot's evaluation order), so it *under*-attributes R3/R4. The counter is descriptive only.
**Gate 1 does not use it** — R1-sufficiency is measured by graph equality, which is
attribution-free.

## 2. Design

Three arms, same 600 pilot SCMs, regenerated **bitwise** from the stored `(seed,p,deg)`
(`worst |τ − τ_pilot| = 3.331e-16`, threshold 1e-12, else the run is void).

| arm | `K` | perturbation | `K'` still a tiering? |
|---|---|---|---|
| **G-flip** | generic, true, `\|K_G\| = \|K_T\|` | reverse ρ ∈ {1,2,3} statements | n/a |
| **T-flip** | tiered, true, **full** cross-tier set | reverse ρ ∈ {1,2,3} statements | **NO** ← the pilot's arm |
| **T-move** | tiered, true, **full** cross-tier set | move ρ ∈ {1,2} nodes to another tier | **YES** ← the new experiment |

523/600 SCMs retained, **355 with all three arms**. `|K|` = {2:159, 3:124, 4:47, 5:20, 6:4, 7:1}.
Balls enumerated exhaustively (no subsampling occurred at the caps used).
Members: G-flip 2,730 · T-flip 2,730 · T-move 48,911 (of which 9,863 = 20.2% are **no-ops**,
`n_changed = 0`, excluded from every rate below).

## 3. GATE 1 — the premise. **PASS**, and it is exact.

Pre-registered: `s ≥ 0.99` → pass; `s < 0.99` → REFUTED.

```
UNPERTURBED tiered K,  R1-only closure == R1-R4 closure : 355/355       = 1.0000
UNPERTURBED generic K, same check                       : 319/355       = 0.8986
PERTURBED  G-flip                                       : 1342/1528     = 0.8783  [0.8609,0.8937]
PERTURBED  T-flip                                       : 1296/1540     = 0.8416  [0.8225,0.8589]
PERTURBED  T-move                                       : 39267/39267   = 1.0000  [0.9999,1.0000]
>>> GATE 1 = PASS (s = 1.0000)
```
Holds at **1.0000 in every `|K|` stratum** (≥2, ≥3, ≥4) and **1.0000 on the subset that
asserts a false statement** (13,584/13,584).

🔴 **Honesty flag, and it matters — this is limitation 5.6 repeating.** `s = 1.0000` for
T-move is **entailed**, not discovered: a perturbed tiering *is* a tiering, so if Bang &
Didelez (arXiv:2306.01638) hold, R1-sufficiency follows by their lemma. Do **not** report
1.0000 as a finding. What is *not* entailed, and is the reportable number, is **0.8416** —
the rate at which the pilot's own flip-perturbation **broke** the premise. That quantifies
limitation 5.2: the pilot's tiered arm violated the lemma's hypothesis on ~16% of members,
which is why it saw R1 behaving as a long-range propagator.

**Consequence for the paper: the ball must be over TIERINGS, not over flip-sets.** That is
a design fact the theorem needs, and it was not established before today.

## 4. GATE 2 — containment. **SUPPORTED.**

`r = Δorient / n_changed`, consistent members, `n_changed > 0`.

```
arm            n     mean   median    q75    q90    q99    max
G-flip      1528    1.380    1.000  2.000  2.000  3.000  5.000
T-flip      1540    1.407    1.000  2.000  2.000  3.610  6.000
T-move     29404    0.899    1.000  1.000  1.000  2.000  4.000

Mann-Whitney U  H1: r(T-move) > r(G-flip)   p = 1.000e+00
Mann-Whitney U  H1: r(T-move) < r(G-flip)   p = 0.000e+00
median  T-move 1.000 vs G-flip 1.000  -> not greater
q90     T-move 1.000 vs G-flip 2.000  -> not greater
max     T-move 4.000 vs G-flip 5.000
>>> GATE 2 = SUPPORTED
```

Raw `Δorient` at **matched** `n_changed` (the rate normalisation is a design choice; this
is the same comparison without it):

| `n_changed` | G-flip n / med / q90 / max | T-move n / med / q90 / max |
|---|---|---|
| 1 | 750 / 1.00 / 2.00 / 5 | 15268 / 1.00 / **1.00** / 4 |
| 2 | 575 / 2.00 / 3.00 / 6 | 10183 / 2.00 / **2.00** / 5 |
| 3 | 203 / 3.00 / 4.00 / 8 | 3321 / **2.00** / **3.00** / 5 |

**The effect is entirely in the tail.** The median is 1.0 in all three arms; what R1-closure
buys is q90 **1.0 vs 2.0** and a shorter maximum. Containment is real and it is a factor of
~2, not an order of magnitude.

Cascade (rule-propagated orientations beyond those named in `K'`):
`G-flip 0.524 mean / 0.420 frac>0` · `T-flip 0.577 / 0.428` · `T-move 0.303 / 0.260`.

## 5. Post-hoc attacks — all five survived

Everything in this section was written **after** Gate 2 returned SUPPORTED, in order to
attack it, and is labelled post-hoc.

**A. No-op contamination (`n_changed = 0` removed).** Decomposition over real perturbations:

```
arm            N   caught     loud    inert   silent chg_valid
G-flip      2730    0.440    0.192    0.285    0.082     0.000
T-flip      2730    0.436    0.187    0.286    0.091     0.000
T-move     39048    0.247    0.198    0.532    0.023     0.000
```
`changed_still_valid = 0` reproduces in all three arms (still an **empirical** zero — do not
print as a theorem).

**B. The better-controlled contrast — T-move vs T-flip, *same* `K`, *same* SCM, only the
error model differs.** This isolates R1-closure from elicitation form and should have been
the primary comparison:
`MWU T-move < T-flip p = 0.000e+00`. And the unplanned **null control behaves**:
`MWU T-flip vs G-flip p(less) = 0.491, p(greater) = 0.509` — i.e. the *form* of `K`
(tiered vs generic) does **nothing** under the same error model. **The entire effect is the
error model, not the elicitation.** That is the sharpest sentence this experiment produces.

**C. `|K|` robustness.** `|K|≥2 / ≥3 / ≥4`: T-move mean `0.899 / 0.878 / 0.867`,
G-flip `1.380 / 1.281 / 1.243`; MWU T-move<G-flip `p = 0 / 1.1e-261 / 1.3e-124`.
Gate 1 = 1.0000 in every stratum.

**D. Capability (C13).** The Gate-1 statistic returns 1.0000 for T-move and 0.8416 for
T-flip, so a sub-0.99 value was demonstrably reachable. The Gate-2 statistic has range
(44% of G-flip members exceed the median; max 5–6 vs median 1) and correctly returns a
**null** between the two flip arms — the pipeline does not manufacture effects.

**E. The decisive one — is a tier move a real error, or just *less* knowledge?**
A tier move can reverse, remove *or* add statements. If it mostly removed true ones, the
arm would measure nothing.

```
n_rev   mean 0.506  frac>0 0.395
n_rem   mean 0.668  frac>0 0.553
n_add   mean 0.576  frac>0 0.510
n_false mean 0.795  frac>0 0.595     <- asserts >=1 statement CONTRADICTING D*

subset                        n   caught   silent  d_orient med   q90   R1-suff
asserts >=1 FALSE stmt    23228    0.415    0.039          2.00  3.00    1.0000
only true add/remove      15820    0.000    0.000          1.00  2.00    1.0000
```
**59.5% of tier moves assert something false**, so the arm measures real error. Restricted
to that subset:

```
arm                               n     mean   median      q90      max
T-move (FALSE-asserting)      13584    0.978    1.000    1.000    4.000
T-flip (all, same K)           1540    1.407    1.000    2.000    6.000
G-flip (all)                   1528    1.380    1.000    2.000    5.000
  MWU T-move(false) < T-flip : p = 6.637e-287
  MWU T-move(false) < G-flip : p = 2.427e-296
```
This also **retires the apparent cost** noted in attack A: caught-free is 0.415 on the
false-asserting subset vs 0.440 for G-flip — comparable. The low 0.247 was the 40.5% of
pure-omission moves, which are correctly harmless (caught 0.000, silent 0.000).

## 6. What did NOT improve — the constraint on the theorem's *form*

Damage `|bias|/|τ|` conditional on `O*` changing, `n_changed > 0`:

```
arm            n     mean   median      q90      q99      max
G-flip       224    1.123    0.575    2.762    6.210   15.334
T-flip       249    1.057    0.509    2.816    5.977   15.334
T-move       897    1.049    0.463    2.762    5.492    6.425
```
**Heavy-tailed in every arm, including the R1-closed one.** q90 ≈ 2.76 |τ| and q99 ≈ 5.5 |τ|
in T-move. The pilot's §3 conclusion is unchanged and now holds *inside* the tiered
sublattice: **the theorem must be a breakdown-radius / probability statement and can never
be a continuity or Lipschitz bound.** `"small knowledge error ⇒ small bias"` remains false
and remains falsified by our own data.

Utility ⊥ safety reproduces exactly, in all three arms: **0 damaging perturbations of
10,896 where the CPDAG was already amenable** (T-move), vs 897/38,015 where `K` is
load-bearing.

## 7. Limitations

1. **Gate 1's 1.0000 is entailed by Bang & Didelez, not discovered** (§3). The non-entailed
   number is T-flip's 0.8416.
2. **Gate 3 (finiteness flavour) is uninformative as measured.** Distinct `O*` per ball:
   G-flip 1.245, T-flip 1.310, T-move 1.315 — but the T-move ball averages ~138 members vs
   ~7.7 for the flip balls, an 18× difference. Ball sizes must be matched before any
   "reachable-`O*` is bounded" statement. **Not done here.**
3. **No enumerator for the MPDAG sublattice reachable under tier-consistent `K`.** This
   experiment perturbs a *given* tiering; it does not enumerate the sublattice, which is
   what a *finiteness theorem* is ultimately about. That remains the largest missing piece
   and is likely a proof, not a script.
4. **Scope, unchanged and not widened:** ER graphs, `p ∈ {5..8}`, expected degree
   ∈ {1.5,2.0,2.5}, linear-Gaussian iSCM, single X, single Y, causal sufficiency,
   faithfulness, population (no sampling error). Scale-free/hub graphs, `p > 10`,
   degree > 6 and every real dataset remain untested.
5. **77/600 SCMs dropped** (no tiering with ≥2 cross-tier orientable edges in 20 draws);
   another 168 lacked all three arms. Selection is on the *elicitation draw*, never on an
   outcome, and both arms of a retained SCM share `|K|`.
6. **The tiering is derived from a topological order of `D*`**, so the base `K` is true by
   construction. An expert whose *whole* mental model is wrong (correlated error — pilot
   limitation 5.3) is still not modelled. A tier move is a structured, correlated error,
   which is closer than an independent flip, but it is not the adversarial case.
7. `changed_still_valid = 0` is an **empirical** zero in all arms, not a theorem.

## 8. What this hands the paper

- The R1-closed finiteness theorem's **premise survives the right perturbation**, and the
  right perturbation is **over tierings**. The pilot's flip-ball was the wrong ball, and
  0.8416 measures how wrong.
- **Containment is real but modest**: tail halved (q90 1.0 vs 2.0), median unchanged.
  A theorem should be stated over the tail/probability, consistent with §6.
- The **cleanest sentence** available: *the form of the elicited knowledge buys nothing
  (T-flip vs G-flip, p = 0.49); what buys containment is that a tier-level error stays
  inside the R1-closed sublattice* (T-move vs T-flip, p ≈ 1e-287, same `K`).
- Per the ICLR-anatomy constraints, **do not lead with the theorem** (adding a theorem to a
  paper is a measured null, OR 0.96 [0.80,1.15], bought at reviewer confidence −0.22).
  This result supports keeping it as a §4 lemma, which is where the refinement phase
  already put it.

## 9. Files

```
PREREGISTRATION.md            thresholds, fixed before the run
code/meek_rules.py            rule-selectable Meek closure (new)
code/test_meek_rules.py       validation against the pilot's machinery
code/run_e2.py                three arms, 600 SCMs
code/analyse_e2.py            PRE-REGISTERED gates
code/analyse_e2b.py           post-hoc attacks A-D
code/analyse_e2c.py           post-hoc attack E
code/make_fig.py              fig_e2.png/pdf
results/e2_raw.json.gz        48,911 + 2,730 + 2,730 member records
results/e2_summary.json       machine-readable verdict
results/*_out.txt             verbatim console output of every command above
```
