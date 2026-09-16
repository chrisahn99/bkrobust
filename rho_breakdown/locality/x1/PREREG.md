# X1 PRE-REGISTRATION — the spurious required-edge arm: the perturbation class the abstract promises and no run has ever performed

**Written 2026-08-19, BEFORE any line of X1 code and BEFORE any X1 number exists anywhere.**

Everything quoted below as a baseline is already on the record and was produced by earlier runs:
`../../2026-08-19_e1prime-se-criterion/` (E1′), `../../2026-08-14_latent-causal-iclr2027-rerank/RANKING.md`
(the 08-14 rerank), and the pilot at `~/research-pilots/latent-causal-rho-breakdown-knowledge/`.
Nothing under an *assertion* operator has ever been run: both codebases contain exactly one
perturbation operator, `Kp = [(v,u) if idx in flip else (u,v) …]`
(`run_linear.py:132`, `run_arm1.py:84`), and `graphs.py:191-192` **raises** on any statement whose
pair is not an edge of the CPDAG. `RANKING.md:729-742` records the same three lines and calls the
two classes disjoint.

---

## 0. What X1 is for, and the six defects it inherits

The ICLR abstract (`~/rho-breakdown-iclr2027/abstract/round2/C-triage.tex:19-21`) says:

> *Statement errors cover reversed orientations **and required edges asserted between non-adjacent
> variables.***

`EXPERIMENTS.md` scores the evidence behind that sentence as **none**. It is the only sentence in
the abstract whose evidence column is empty, and round-2 review credit was earned partly on it
(`REVIEWS-round2.md:41-42`, R4: *"All three now name spurious required edges between non-adjacent
variables, which was my central complaint"*; C's mean rose 6 → 7.25). X1 either buys that sentence
or removes it.

| # | defect on the record | addressed here? |
|---|---|---|
| D1 | **S7 has zero code behind it.** Neither codebase ever asserts a non-adjacent required edge. | ✅ this is the whole run |
| D2 | **The two taxonomies are inverted.** `graphs.py:191` files a b-LOAD-style false edge in the *most reassuring* bucket ("caught free by Meek-inconsistency"); b-LOAD's `update_graph_mpdag` (`mb_by_mb.py:495`, `g[mask] = background_knowledge[mask]`) re-stamps it after every local update, so it is never caught at all. | ✅ both semantics run as named arms, §2 |
| D3 | **S3 (67 %) and S4 (14.1 %) are scoped to reversals only** and the abstract prints them in the same paragraph that promises both classes (`EXPERIMENTS.md:13-14`). | ✅ matched-SCM two-strata table, §4 |
| D4 | **A26's `revive_if` is unevaluated.** `RANKING.md:783-787` registers a runnable condition and it has never been run. | ✅ evaluated verbatim, §4.3 |
| D5 | **`0 of 791` was ordered retired on 2026-08-14 and is still in the live abstract.** | ⛔ **not X1's to fix** — X2 owns S5. X1 prints no 0/791 and no derivative of it. |
| D6 | The locality property (S5/S8/S9) has only ever been tested against a **skeleton-preserving** operator. | ✅ X1 supplies the skeleton-*changing* stress test, §4.4 |

**A correction to the vault's own wording, found by reading b-LOAD's source today and registered
here before the run.** `RANKING.md` and the proposed registry entry A26 say the two classes are
*disjoint*. That is true of **our** operator and **not** true of b-LOAD's sampler.
`sample_local_background_knowledge_noised` (`~/LOAD_in_MPDAG/src/background_knowledge.py:166-170`)
draws false required edges from

```
false_req_candidates = [(u,v) for u in target_set for v in nodes
                        if u != v and (u,v) not in true_edges_set]
```

which **admits `(u,v)` whenever `v -> u` is a true edge** — i.e. b-LOAD's false-required class is a
**mixture** of reversals and spurious edges, not a pure second class. The reversal share of the
candidate pool for a target `u` is exactly

```
|pa_D(u)| / ( (p-1) - |ch_D(u)| )        (a definition, not a measurement)
```

so it is O(deg / p) and small but nonzero. **X1 measures the two pure strata and declares the
mixture NOT COVERED (§8.6).** The registry entry A26 should be amended to "disjoint *as our
operator implements them*"; that is a proposal for Zé, not an agent edit.

---

## 1. The estimand

### 1.1 Objects, all defined

| symbol | definition |
|---|---|
| `V = {0,…,p-1}` | vertices |
| `D` | true DAG, `p x p` int8 amat; `D[i,j]=1` ⇔ `i -> j` |
| `C = CPDAG(D)` | `dag_to_cpdag(D)`; `G[i,j]=G[j,i]=1` ⇔ `i -- j` |
| `S = skeleton(C)` | `((C+C.T)>0)`. **`S = skeleton(D)` identically** (`dag_to_cpdag` starts from `skeleton(D)` and only ever zeroes one direction) |
| `A(C) = {{a,b} : S[a,b]=1}` | the **adjacent** pairs |
| `N(C) = {{a,b} : a≠b, S[a,b]=0}` | the **non-adjacent** pairs — **the class S7 names, and the class no run has ever touched** |
| `U(C)` | `undirected_edges(C)`, the pairs `a--b` |
| `(X,Y)` | query, single-node, drawn among pairs with a possibly-causal path in `D` |
| `K = ((u_1,v_1),…,(u_k,v_k))` | asserted knowledge, `k = |K| = 4`, `{u_t,v_t} ∈ U(C)`, oriented to agree with `D`. **Every statement in `K` is TRUE.** |
| `Meek(·)` | `meek_closure`, R1–R4 to fixpoint |
| `O*(x,y,G)` | `parents_of_set(G, cn(x,y,G)) \ forb(x,y,G)`, or `None` if `G` is not amenable |

### 1.2 The unperturbed anchor — identical in every arm, by construction

```
G0    = apply_background_knowledge(C, K)          # K is true and lies on U(C): no arm differs here
O0    = O*(X,Y,G0)                                # arm is skipped if O0 is None
est0  = beta_X( Y ~ X + O0 )   under the population Sigma
S0    = {X} u O0 ,  k_reg = |S0|
sigma0^2 = Sigma_YY - Sigma_{Y S0} Sigma_{S0 S0}^{-1} Sigma_{S0 Y}
SE_n  = sqrt( sigma0^2 * [Sigma_{S0 S0}^{-1}]_{XX} / (n - k_reg - 1) ) ,   SE_inf = 0
tau   = ((I-A)^{-1})[X,Y]                          # ORACLE. Evaluation only. Never an input.
```

`O0`, `est0`, `se_factor`, `k_reg`, `sigma0^2` are **shared across all arms of the same SCM**. This
is the like-for-like anchor and it is enforced as gate **G1**.

### 1.3 The perturbed knowledge sets

Fix `rho ∈ {1,…,k}` and a slot set `F ⊆ {1,…,k}`, `|F| = rho` (`itertools.combinations` order,
`run_arm1.py:83`). Each arm defines a replacement `r_t` for slot `t`, **pre-drawn once per SCM and
held fixed across every `rho`** (so the balls nest and `d(rho)` is monotone — gate G6):

```
K'(F) = ( r_t if t in F else (u_t, v_t) )_{t=1..k}
```

| arm | replacement `r_t` | draw law |
|---|---|---|
| **R** *(reversal — the control)* | `(v_t, u_t)` | deterministic |
| **S-uni** *(spurious, uniform — PRIMARY)* | `(a_t, b_t)` | ordered pair drawn uniformly without replacement from `{(a,b) : {a,b} ∈ N(C)}` |
| **S-loc** *(spurious, b-LOAD-matched — SECONDARY)* | `(a_t, b_t)` | drawn uniformly without replacement from `{(a,b) : a ∈ {X,Y}, b ∈ V\{a}, (a,b) ∉ E(D)}` — **b-LOAD's own candidate pool verbatim**, mixture included, with each drawn statement labelled `spurious` / `reversal-of-true-edge` / `reorientation-of-directed-CPDAG-edge` so the strata can be separated post hoc |
| **P-null** *(placebo, gate G5)* | a pair already **directed** in `C` with the orientation it already has | deterministic |

### 1.4 The three operators — unambiguous enough that two people implement the same thing

**Operator R — `bk_reverse` — the status quo, imported unmodified.** `apply_background_knowledge`
(`graphs.py:179-199`), verbatim. Raises `MeekFail` on: conflict with an existing orientation,
**non-adjacent pair**, directed cycle, new v-structure vs `C`.

**Operator S — `bk_assert(C, K', mode)`.** The only new code. Pseudocode, normative:

```
G = C.copy()
for (i, j) in K':                                  # in K' order; closure after EACH statement,
    if G[i,j] == 0 and G[j,i] == 0:                #   matching apply_background_knowledge
        G[i,j] = 1 ; G[j,i] = 0                    # ADD the edge. skeleton(G) grows. flag: added
    elif is_directed(G, j, i):
        if mode == 'strict':  raise MeekFail('conflict')
        else:                 G[i,j] = 1 ; G[j,i] = 0      # b-LOAD overwrites; flag: overwrote
    else:
        G[j,i] = 0                                  # ordinary orientation of an undirected edge
    G = meek_closure(G)
    G[i,j] = 1 ; G[j,i] = 0                         # RE-STAMP (b-LOAD's g[mask] = bk[mask])
flags.cycle   = has_directed_cycle(G)
flags.vstruct = (v_structures(G) != v_structures(C))
flags.skel    = (skeleton(G) != skeleton(C))
if mode == 'strict' and (flags.cycle or flags.vstruct):  raise MeekFail(...)
return G, flags
```

- **ARM S1 = `mode='strict'`** — operator (i) of the design brief: *add the edge to the skeleton,
  re-run Meek closure, and still apply the pilot's own coherence checks*. The one branch that is
  **deliberately disabled** is `"not an edge of the CPDAG"`; everything else the pilot checks is
  kept. This is "a careful tool that does not audit its knowledge against the estimated skeleton".
- **ARM S2 = `mode='protect'`** — operator (ii): **b-LOAD's semantics**. Nothing ever raises.
  Cycles and new v-structures are *recorded*, not rejected, because `update_graph_mpdag` records
  nothing and rejects nothing.

**S1 is a derived view of S2, not a second run.** The two constructions coincide whenever S1 does
not raise (declared definitional, §5.5), so the implementation performs **one** ball per arm and
`strict` is reconstructed at analysis time from `flags`. This is a cost decision and a correctness
one: it makes S1 ⊆ S2 exact rather than approximately so.

**Two operator facts that must be gated, not assumed** (§6, G7/G8):
- `meek_closure` only ever orients *undirected* pairs (`graphs.py:118-160`), so it cannot flip a
  stamped directed edge and the re-stamp should be a **no-op**. Gate: re-stamp changes `G` in
  **0.000** of members. If it ever fires, that is a real Meek/b-LOAD divergence and is reported.
- `apply_background_knowledge` closes after **each** statement; b-LOAD stamps **all** then closes.
  Gate G8 reports the rate at which `close-after-each` and `stamp-all-then-close` disagree on `G`.
  Primary is close-after-each (pilot-matched). A nonzero disagreement rate is a **finding**, printed
  with its `n`, not a bug to be quietly resolved.

### 1.5 The five outcomes — reported separately, never merged

For each member `m = K'(F)` of each arm:

```
(a) caught      ∈ {none, nonadjacent, conflict, cycle, vstruct}   # WHICH check fired, not just that one did
(a') skel_detect = 1 iff any asserted pair of m lies in N(C)      # free, data-side detectability
(b) ostar_changed = ( O*(X,Y,G_m) != O0 )                          # frozensets, set()-normalised out of JSON
(c) silent       = consistent AND amenable AND NOT is_valid_adjustment_set(D, X, Y, O_m)
(d) loud         = consistent AND (O*(X,Y,G_m) is None)
(e) moved(n,z)   = |est_m - est0| > z * SE_n ,  est_m = beta_X( Y ~ X + O_m ) under Sigma
(f) relbias      = |est_m - tau| / |tau|            # reported ONLY conditional on (c)
```

with `consistent` arm-specific: **R** = `apply_background_knowledge` did not raise; **S1** = `bk_assert(strict)`
did not raise; **S2** ≡ `True` **by construction** (declared definitional, §5.1).

Derived, as in E1′: `d(rho) = max_{m : |F| ≤ rho, consistent, amenable} |est_m − est0|`;
`rho*_se(n,z) = min{rho : d(rho) > z·SE_n}`, `= k+1` if CENSORED; `rho*_ident`;
`rho*_any = min(rho*_se, rho*_ident)`. **Reported per channel, never only as `rho*_any`** — the
E1′ post-mortem showed that folding identification failure into the movement statistic is what
manufactured a false contradiction with S5.

### 1.6 The hop-distance convention, fixed now

`dist = hop_dist(C, [X,Y], p)` — BFS on `skeleton(C)`, **computed before any edge is added**. A
statement `(a,b)` carries `hop = min(dist[a], dist[b])`.

🔴 **`hop_dist` returns the sentinel `1 << 20 = 1048576` for unreachable nodes** (`run_arm1.py:49`).
X1 converts it to the string label `"unreachable"` at record time. It is **never** pooled into
`>= 2`, **never** averaged, and always printed as its own row with its own `n`. (E1′'s
`large` "dist ≥ 2 | 1.000 (n=331)" cell silently contained 57 unreachable SCMs; X1 does not repeat
that.) **Every stratum is printed with its `n`, however small.** No `< 10` suppression, no em-dash
for an empty cell — `analyse_arm1.py:213-214` suppresses and that is how the most informative
strata vanished.

---

## 2. Arms, ensembles, seeds, power

### 2.1 Arms per SCM (all inside one worker call, on the same `C, D, Sigma, X, Y, K`)

`R`, `S-uni`, `S-loc`, plus two placebo balls `P-null` and `P-frozen`. Five balls of ≤ 15 members
each; `S1` is derived from `S-uni`/`S-loc` flags at analysis time.

### 2.2 Ensembles

Grids imported **unchanged** from `run_arm1.py:35-42` so that ARM R reproduces E1′ exactly (gate G2).

| ensemble | grid | target `n_scm` | why |
|---|---|---|---|
| `original` | p ∈ 5..8, deg ∈ {1.5,2,2.5}, `need_u ≥ 3` | **2000** | **PRIMARY.** This is the pilot's own grid — the *only* grid on which 67 % and 14.1 % were measured, so it is the only grid on which S3/S4 can be compared like for like. |
| `licensed` | p ∈ 5..10, deg ∈ {1.5,2,2.5,4,6} | 2000 | E1′'s primary; keeps the two runs commensurable |
| `large` | p ∈ 15..25, deg ∈ {1.5,2,2.5} | 1200 | `|N(C)|` grows like `p^2` while `|U(C)|` grows like `p·deg`, so the spurious class is *the* class whose prevalence and locality are p-dependent. Also the ensemble where E1′'s locality law was weakest. |
| `k8` | — | **NOT RUN** | 127 analysed SCMs at K4 ⇒ Wilson half-width ≈ 0.06 at p ≈ 0.14, four times too wide to decide §4.2. `|K|` is X3's axis, not X1's. Stated so the gap cannot read as coverage. |

`OVERSAMPLE = 3`, `MAX_K = 4`, `max_rho = 4` (**full ball**, `rho = 1..|K|`), `TAU_FLOOR = 1e-3`
at analysis with a 1e-6 sensitivity reported (AUDIT A8 of E1′).

### 2.3 Seeds — and the reason for a spawned stream

The `(p, deg)` job stream is `np.random.default_rng(20260819)` and the per-SCM stream is
`np.random.default_rng(seed)`, **both byte-identical to `run_arm1.py:171-174`**, and the draw order
inside `analyse_one` up to and including `K_all` is untouched. Spurious pairs are drawn from a
**separate spawned stream** `np.random.default_rng([seed, 0xA55E27])`. Without this, adding one
`rng` call re-rolls every SCM and gate G2 becomes untestable — the exact failure mode
`run_arm2.regenerate` was built to catch.

### 2.4 Power, stated as a number before the run

At `rho = 1` there are `4 × n_analysed` statements per arm. E1′ retained 2 637 / 2 087 / 2 443
analysed SCMs from 4 000 draws on `original` / `licensed` / `large`, so X1 at `n_scm = 2000`
expects **≈ 4 000–8 000 statements per arm per ensemble**.

Two-proportion, `alpha = 0.05`, `n = 4 000` per arm, baseline `p ≈ 0.14`:

```
se_diff = sqrt( 2 * 0.14 * 0.86 / 4000 ) = 0.00776
MDD at 80 % power = (1.96 + 0.84) * 0.00776 = 0.0217
```

**The registered bar of 0.03 in §4.2 is therefore attainable at this `n`, with margin.** Wilson
95 % half-width at `p = 0.14, n = 4 000` is `0.0107`. Every rate in the RESULTS carries `n` and a
Wilson interval; **no fraction is printed without its denominator.**

### 2.5 Budget

12 workers on betelgeuse, `~/e1venv/bin/python`. Three balls per SCM instead of E1′'s one, so
≈ 3× E1′'s per-SCM cost; E1′'s four `arm1` sweeps plus `arm2` finished in under 15 min.
**Hard cap: 20 min per ensemble**, `large` capped by its grid (`deg ≤ 2.5` at `p ≤ 25`, never
`deg ∈ {4,6}` — `possibly_causal_paths` is an uncapped DFS and that pairing does not finish).
A **48-SCM smoke run** runs the gates first; if it produces a number, that number is disclosed in an
addendum and the smoke file is deleted, exactly as E1′ ADDENDUM 1 did.

---

## 3. What is being claimed

> **X1's claim.** The spurious required-edge class — required edges asserted between variables
> non-adjacent in the CPDAG — has a **materially different risk profile** from the reversed-orientation
> class the campaign has measured, in the direction of **less catchable and at least as damaging**.
> Consequently the abstract's 67 % / 14.1 % figures do **not** transfer to it, and S7 cannot stand on
> the reversal arm's evidence.

Two things about this claim, fixed now:

- It is a claim about **two measured strata**, not about a single blended number. If X1 fires, C
  survives by printing a two-row table, which is a better paper than the one that quietly let one
  class stand for both.
- It is **partly definitional and X1 says so up front** (§5): under b-LOAD semantics the catch rate
  is zero by construction, and every spurious statement is skeleton-detectable by construction. The
  *empirical* content is entirely in the **damage** rates and their **locality**.

---

## 4. Pre-registered predictions, with numeric bars, before the answer is known

All rates below are **per statement**, at `rho = 1`, on `original`, unless stated. Every one is
reported with `n` and a 95 % Wilson interval whichever way it comes out.

### 4.1 X1-P1 — catchability

| quantity | prediction | what would surprise |
|---|---|---|
| ARM S2 Meek-catch rate | **0.000 exactly** — definitional, gate G3 | any nonzero value is a bug |
| skeleton-detectability of a spurious statement | **1.000 exactly** — definitional, gate G3 | any value < 1 is a bug |
| ARM S1 catch rate `C_S1` (cycle ∨ new-v-structure only) | **PREDICTED `C_S1` ∈ [0.40, 0.95]** | **< 0.40**: even a coherence-checking tool is nearly blind to the class ⇒ S7 is a live hazard, not a coverage note. **> 0.95**: the class is caught essentially always by v-structure bookkeeping ⇒ S7 is safe but *empty*, and the honest edit is to keep the sentence only alongside the two-strata table. |
| ARM R Meek-inconsistency rate | reproduce E1′ to ±0.02 (`original` 0.3576) | gate G2 |

Mechanism behind the band: adding `a -> b` creates a new unshielded collider `(a,b,c)` for every
`c -> b` with `c` non-adjacent to `a`, and destroys any existing v-structure `(a,m,b)` — so the
`v_structures(G) != v_structures(C)` test has real power, but only when `b` already has parents.

### 4.2 X1-P2 — the headline: does the class damage more? *(load-bearing)*

`R_dmg^R` = P(silent bias | ARM R), `R_dmg^S` = P(silent bias | ARM S2, S-uni), same SCMs, same
slots, same `rho = 1`, **unconditional** (not conditioned on passing the check — the two arms'
pass rates differ *by construction*, so conditioning on passing would compare different populations).

- **SUPPORTED** — the class is worse, S3/S4 do not transfer:
  `R_dmg^S − R_dmg^R ≥ 0.03` **and** the two 95 % Wilson intervals are disjoint.
- **REFUTED** — the class is *less* damaging: `R_dmg^S − R_dmg^R ≤ 0` with disjoint intervals.
  Then S7 is a coverage note, S3/S4 hold a fortiori, and the sentence still needs X1's table printed
  because "we checked and it is milder" is a result and "we never checked" is not.
- **INDISTINGUISHABLE** — the intervals overlap. **This is exactly A26's `revive_if` condition** and
  it *revives* the b-LOAD transfer that 08-14 retracted. It must be written that way, not smoothed
  into "no difference found".

The same triple is reported on `licensed` and `large` as secondaries, with the caveat that S3/S4
were never measured on those grids.

### 4.3 X1-P3 — A26's `revive_if`, evaluated verbatim *(load-bearing)*

`RANKING.md:783-787` registers:

> *"a matched-budget run of an assertion operator … returns a silent-bias rate whose 95 % Wilson CI
> **overlaps** the reversal arm's in-house baseline of 754/8085 = 9.3 % [8.7 %, 10.0 %], and whose
> conditional |bias|/|τ| does not stochastically dominate it."*

Evaluated on **ARM S2 / S-loc**, the b-LOAD-matched sampler, pooled over the cell that matches
`754/8085` as literally as the new data permit (all `rho ∈ {1,…,4}`, `original` grid), and reported
with the exact cell definition printed on the same line as the number — `RANKING.md:768` already
names *"two pooled figures for the same object"* as the seam a reviewer finds.

- **PREDICTED**: the CI does **not** overlap `[0.087, 0.100]` and lies **above** it,
  `R^{S-loc} ∈ [0.20, 0.60]`. Mechanism: b-LOAD anchors every false edge at a target node, so every
  S-loc statement is at hop-distance **0 by construction**, and hop-0 is where every O\*-change in
  E1′ lives; and `X -> b` manufactures a possibly-causal path while `Y -> b` puts `b` into
  `forb`, so both draw directions move `O*`.
- **SURPRISE-LOW** (`< 0.15`): query-anchored false edges are nearly harmless. That would
  *strengthen* S9 sharply and would mean b-LOAD's observed robustness has a different explanation
  than locality — a result worth its own paragraph.
- **SURPRISE-HIGH** (`> 0.60`): the class is catastrophic and the abstract's "narrow damage"
  framing is a property of the reversal operator alone.

### 4.4 X1-P4 — locality under a skeleton-changing perturbation *(feeds S5/S8/S9; coordinates with X2)*

The reversal operator **cannot** change `skeleton(C)`, so it cannot create a possibly-causal path
that did not already exist. The assertion operator can. This is the first stress test locality has
ever had against that mechanism.

- **Statistic**: among ARM S2 / S-uni statements at `hop >= 1` (unreachable reported separately),
  at `rho = 1`, the count and rate of `silent` and of `ostar_changed`.
- **SUPPORTED (locality is *not* a property of knowledge, only of reversals)**: `>= 5` damaging
  statements at `hop >= 1`, denominator printed.
- **REFUTED (locality survives a skeleton-changing perturbation)**: `0` damaging statements out of
  `>= 2 000` at `hop >= 1`. That is a **stronger** result than the abstract currently claims and
  S8/S9 gain — though S8 still needs the proof R1 named (*"no distant statement enters the
  neighbourhood through Meek closure"*), which X1 does not supply.
- Anything between: reported as a rate with `n` and Wilson CI, no verdict word.

**Division of labour, fixed now:** X1 measures locality for the **spurious** class; X2 owns the
**reversal** class, the per-statement/per-SCM unit mismatch, and the retirement of `0/791`. Neither
prints the other's denominator.

### 4.5 X1-P5 — severity, the second half of A26's `revive_if`

Conditional on `silent`, compare `relbias = |est − tau|/|tau|` between ARM R and ARM S2 on
`original`, `rho = 1`.

- **Stochastic dominance declared only if** Cliff's δ ≥ 0.15 with a 95 % bootstrap CI excluding 0.
  Otherwise: **"not separated"**, which is the `revive_if`'s second condition met.
- Medians, IQRs and the full ECDFs printed either way.
- ⚠️ A two-sample test on a *bias magnitude* distribution is legitimate. A p-value on `rho*` is
  **not** and never appears (E1′ AUDIT A5: `rho*` is a worst-case-over-the-ball sensitivity radius
  in the Rosenbaum-Γ tradition, not a test).

### 4.6 X1-P6 — the movement channel, on E1′'s grid

Outcome (e) reported over `n ∈ {200, 500, 1e3, 2e3, 5e3, 1e4, 2e4, ∞}`, `z ∈ {1.0, 1.96, 2.576}`,
`z = 1.96` and `n = 20 000` headline, `n = ∞` printed alongside. `rho*_se`, `rho*_ident` and
`rho*_any` reported **as three separate columns**. No prediction registered: this channel exists so
that X1's arms are commensurable with E1′, not to decide anything.

---

## 5. Declared definitional — NOT evidence, written down before the run

1. **ARM S2's Meek-catch rate is 0.000.** The operator is *defined* not to raise. Quoting it as
   "spurious edges are never caught" would be quoting the operator's definition back as a finding.
   The reportable statement is the **S1** rate (§4.1) and the **skeleton** check (item 2).
2. **Every spurious statement is skeleton-detectable, rate 1.000.** `{a,b} ∈ N(C)` is the definition
   of the class and `skeleton(C) = skeleton(D)`, so the assertion contradicts an adjacency the data
   determine. Again: definitional.
3. 🔴 **The testability asymmetry is analytic, not measured, and it is the most important thing X1
   has to say about the abstract.** Within an MEC, all members entail the same conditional
   independences (Verma & Pearl 1990; Chickering 2002), so the *orientation* of an edge is not a
   testable hypothesis — this is the abstract's S2 and it is correct. **Adjacency is not in that
   equivalence class.** `a` and `b` non-adjacent in `D` ⇒ they are d-separated by some `Z ⊆ V\{a,b}`
   (e.g. `pa_D(a)` when `a` is not a descendant of `b`), so a spurious required edge **is** a
   testable hypothesis, with power → 1. Therefore:
   > The abstract asserts *"no test of it has power above its size at any sample size"* (S2) and then
   > asserts that statement errors include a class for which a test has power → 1 (S7). **As written,
   > S2 and S7 are inconsistent with each other.** That is an internal contradiction a theorist
   > reviewer finds in one pass, and it exists in the current draft independently of anything X1
   > measures.
   X1 does not resolve it by measurement. It reports it, and the repair is a wording decision for Zé
   and Chris: either S2 is narrowed to orientations, or S7 is dropped, or S7 is kept **with** the
   sentence that this class is testable and the paper measures what happens when the tool does not
   test it.
4. `d(rho)` is monotone non-decreasing in `rho` (fixed per-slot replacements, max over a superset)
   and `rho*_se(n)` is monotone non-increasing in `n` (`SE_n ∝ n^{-1/2}`, `d` is n-free).
5. **S1 ⊆ S2**: whenever `bk_assert(strict)` does not raise, its graph is identical to
   `bk_assert(protect)`'s. S1 is therefore a *view* of S2 and the two are not independent evidence.
6. `O0`, `est0`, `SE_n` are identical across arms (§1.2). Any cross-arm difference is attributable
   to the perturbation and nothing else — that is the point of the design, not a result.
7. **ARM S-loc statements are at hop-distance 0 by construction** (`a ∈ {X,Y}`). Its damage rate must
   never be compared to a distance-mixed rate without saying so on the same line.

---

## 6. Integrity gates — all must pass or no verdict is readable

| id | gate | pass condition |
|---|---|---|
| **G1** | **like-for-like anchor**: `O0`, `est0`, `se_factor`, `k_reg`, `sigma2` identical across R / S-uni / S-loc / P-null for every SCM | exact equality, **100 %** |
| **G2** | **reversal-arm reproduction**: on SCM seeds shared with `~/latent-causal/e1prime-se/results/arm1_licensed.json` and `arm1_original.json`, ARM R's members match E1′ member-for-member on `consistent`, `amenable`, `ostar_changed`, `ostar_valid`, and `est` to 1e-12 | **100 %**, ≥ 500 shared SCMs |
| **G3** | **definitional gates fire exactly** (§5.1, §5.2): S2 catch rate `= 0.000`; skeleton-detectability `= 1.000` | exact |
| **G4** | **placebo P-frozen** (dead switch, one direction): recompute every ball with `O` frozen at `O0` | `ostar_changed = 0`, `silent = 0`, censoring `= 1.000` at every `n`, every arm |
| **G5** | **placebo P-null** (dead switch, other direction): assert edges already directed in `C` with the orientation they already have, through the **full new code path** | `G = G0` in 100 %; `ostar_changed = 0`; `silent = 0` |
| **G6** | monotonicity of `d(rho)` in `rho`, per SCM per arm | **100 %**, 0 violations |
| **G7** | **re-stamp is a no-op**: `G` after re-stamp equals `G` after closure | rate reported; expected **0.000** divergence. A nonzero rate is printed as a finding |
| **G8** | **closure-order diagnostic**: `close-after-each` vs `stamp-all-then-close` produce the same `G` | rate reported with `n`; no bar, this is a measurement |
| **G9** | **no `tau` leak**: `tau` appears in no function on the instrument path; only in evaluation of (c) and (f) | source check, as E1′ `se.assert_no_tau_leak` |
| **G10** | **sentinel discipline**: no record carries `1048576` into a numeric stratum; `unreachable` is a labelled row | grep of the results JSON returns 0 hits |

**C13 compliance, demonstrated inside this run and not borrowed.** G4 forces every damage rate to
its floor and G5 forces the operator itself to the null on the identical code path. A rate anywhere
between them is a property of the perturbation class, not of the harness.

---

## 7. Which abstract sentence each output speaks to

| sentence | X1 output | if the output goes one way | if it goes the other |
|---|---|---|---|
| **S2** *"no test has power above its size at any sample size"* | §5.3, analytic | — | **Already inconsistent with S7 as both are written.** X1 reports it; the repair is a wording decision, not a measurement. |
| **S3** *"two thirds of false knowledge sets pass the consistency check"* | ARM R pass rate vs ARM S1 pass rate vs ARM S2 pass rate `≡ 1`, `original` grid, `rho = 1`, all three printed on one line with `n` | if S1's pass rate is high, "two thirds" is a **reversal-only** figure and the abstract must say `rho = 1`, generic arm, reversals | if S1's pass rate is low, the strict operator saves the sentence, but only for tools that check — and b-LOAD does not |
| **S4** *"14.1 % of the consistent errors go on to bias the estimate"* | §4.2 `R_dmg^S` vs `R_dmg^R`, unconditional and conditional-on-pass both printed | SUPPORTED ⇒ S4 is re-scoped to reversals and a second row is added | INDISTINGUISHABLE ⇒ A26 revives and the numbers transfer, which must be stated as a *revival*, not as reassurance |
| **S5** *"no misstatement at graph-distance ≥ 1 moved it"* | §4.4, spurious class only | ≥ 5 damaging at hop ≥ 1 ⇒ S5 is false for the class S7 names, and cannot be stated over "statement errors" as a whole | 0 of ≥ 2 000 ⇒ locality survives a skeleton-changing perturbation, a stronger claim than the abstract makes. **X2 owns S5's wording, the ρ ≥ 2 counterexample, and the retirement of `0/791`; X1 prints none of them.** |
| **S7** *"errors cover … required edges asserted between non-adjacent variables"* | the whole run | **buys the sentence** if the two-strata table is printed | if X1 is not run, the sentence has no evidence and R4's round-2 credit was earned on a promise (`REVIEWS-round2.md:41-42`) |
| **S8** *"locality restricts the enumeration, cost does not scale with |K|"* | §4.4 provides the *empirical* half for the spurious class | a hop ≥ 1 counterexample **kills** the pruning heuristic outright | zero counterexamples still leave S8 **unproved and unimplemented** — no pruning step exists in either codebase (`run_arm1.py:83` enumerates all of `K`), so S8 describes an algorithm nobody has written or benchmarked. X1 does not change that. |
| **S9** *"knowledge outside that neighbourhood need not be elicited"* | §4.3 vs §4.4: S-loc (hop 0 by construction) against S-uni stratified by hop | if S-uni at hop ≥ 1 is harmless while S-loc is harmful, S9 **survives and gains a mechanism** | if S-uni damages at hop ≥ 1, S9 falls with S5 |
| **S10–S14** | untouched | — | — |

---

## 8. What this run does NOT cover — in writing, so a gap cannot later read as coverage

1. **An estimated CPDAG.** `C = dag_to_cpdag(D)` is the oracle's. The analyst who runs PC/GES on
   their own data gets `Ĉ`, and — the point of §5.3 — would **see** that the asserted pair is not in
   `Ĉ`'s skeleton. X1 measures damage **conditional on the tool not performing that check**, which
   is b-LOAD's regime and not every regime. **ARM C (discovery on the same rows) is NOT RUN.** This
   is the first thing a reviewer asks about the spurious class specifically, and the honest answer
   is "not tested", not "robust".
2. **b-LOAD itself.** X1 emulates the *semantics* of `initialize_background_knowledge` +
   `update_graph_mpdag` on the pilot's population harness. It is **not** a b-LOAD replication: no
   local Markov-blanket-by-Markov-blanket discovery, no CI tests, no finite samples, no `mb_by_mb`.
   R and Rscript are absent locally and `rpy2` imports in none of the local envs; the b-LOAD-native
   run is Chris's machine.
3. **Forbidden-edge knowledge.** Not perturbed, not asserted. `--frac_forb` is a measured dead
   parameter in b-LOAD (`initialize_background_knowledge` starts from an all-zero matrix), so
   sweeping it returns a guaranteed null for an implementation reason and C13 forbids reporting that.
4. **Tiered knowledge** (Bang & Didelez). The pilot's `tiered_knowledge` is not used here.
5. **`|K| > 4` and the `k8` ensemble** (§2.2). X3's axis.
6. **Mixed knowledge sets.** b-LOAD's sampler produces a mixture of reversals and spurious edges
   (§0). X1 measures the two **pure** strata. Whether the mixture's damage rate is the convex
   combination of the two depends on the two error types not interacting inside Meek closure, and
   **that non-interaction is untested here.** The pure strata bracket the mixture only under that
   assumption.
7. **Nonlinear mechanisms, latent confounding, non-Gaussian errors, multi-node X or Y, real data.**
   Everything is linear-Gaussian iSCM, ER DAGs, single X, single Y, population `Sigma`.
8. **Whether an analyst would ever state such a knowledge set.** X1 measures consequences, not
   plausibility. The uniform draw is a stress test; the b-LOAD draw is one tool's actual sampler;
   neither is elicited from a domain expert.

---

## 9. Scope statement, fixed now

**Licensed:** linear-Gaussian iSCM, ER DAGs `p ∈ [5,25]`, expected degree ≤ 6 (≤ 2.5 for `p ≥ 15`),
single-node `X` and `Y`, `|K| = 4` statements of which `rho ≤ 4` are replaced by either a reversal or
a required edge on a non-adjacent pair, **CPDAG given**, homoskedastic OLS SE, population `Sigma`.
**Not licensed and not claimed:** everything in §8.

---

## 10. Mandatory citations — without these the run produces an over-claim rather than a result

- **Meek 1995** and **Perković, Kalisch & Maathuis, UAI 2017, Algorithm 1** — the background-knowledge
  contract that operator S deliberately breaks. The paper must say *which* clause it drops
  (`"not an edge of the CPDAG"`) and why that is the practitioner-relevant case.
- **Verma & Pearl 1990**; **Chickering 2002** — all MEC members entail identical conditional
  independences. This is what makes orientations untestable **and adjacencies testable**, i.e. §5.3.
- **Henckel, Perković & Maathuis, arXiv:1907.02435** — `O*` and its variance optimality. The
  estimand is theirs.
- **Chen et al., arXiv:2306.07032 (IEEE TPAMI 2025)** — robustness of discovery to *incorrect*
  background knowledge. Nearest prior art for X1 specifically; `RANKING.md` records it as a
  scoop-in-shape for the b-LOAD follow-up. X1's differentiator must be stated against it explicitly,
  not against a generic literature.
- **Constantinou et al., arXiv:2102.00473 (KAIS)** — impact of incorrect knowledge on structure
  learning. Same obligation.
- **gadjid, arXiv:2402.08616 (UAI 2024)** — owns the outcome variable for "how wrong is this graph".
- **Fang et al., arXiv:2207.05067 (JMLR 26)** — identifiability in MPDAGs, polynomial-time.
- **Taeb, Guo & Henckel, arXiv:2511.10625** — the nearest competitor for the whole ρ-breakdown frame.

🔴 **Vault-level gap this run does not close, carried forward verbatim from E1′ AUDIT A10.**
`wiki/activities/active-claims.md` still holds **no** entry for ρ-breakdown, and
`grep -rl 'relevance:.*rho-breakdown' wiki/literature/` returns **0** notes — re-verified today.
`check_claims.py` is therefore structurally blind to every paper listed above. Under CLAUDE.md
Operating Rule #1 the register entry must exist before any of this is written into a manuscript.
Registering it needs Zé's own claim wording: **proposed, not committed by an agent.**

🔴 **A second live vault item X1's result bears on directly.** `active-claims.md:570`, claim
`b-load-local-knowledge-loop`, records its graph-distance falsifier as *"currently NOT firing"* on
the strength of `0/791`. §4.4 measures the spurious-class half of that falsifier. Whichever way it
comes out, the register line is anchored on a retired number and needs Zé's decision **before the
b-LOAD CIKM camera-ready on 2026-08-23**. X1 supplies the measurement; it does not edit the register.

---

## 11. Standing prohibitions for the RESULTS write-up

Carried from the 08-14 rerank and the E1′ post-mortem, so they cannot be rediscovered late:

1. **Never print `0 of 791`**, or any number derived from it. Retired 2026-08-14, twice, as mandatory.
2. **Every fraction carries its denominator, its `n`, and its conditioning on the same line** —
   arm, `rho`, ensemble, hop stratum. `RANKING.md:768`: two pooled figures for the same object is
   how a reviewer finds the seam.
3. **Never merge the five outcomes.** In particular `rho*_any = min(rho*_se, rho*_ident)` folds
   silent movement into loud identification failure; report the two channels separately or the
   locality answer will be wrong in the same way E1′'s was.
4. **Never suppress a stratum** for small `n`, and never print an em-dash where a suppressed cell was.
5. **No p-value on `rho*`, ever.** §4.5's Cliff's δ is on bias magnitudes and is the only inferential
   statistic in the run.
6. **A declared placebo is never quoted as a finding** (G4, G5). E1′'s RESULTS quoted its own G6
   dead-switch as the headline reply to a reviewer; X1 does not repeat that.
7. **Naming hazard:** `output/2026-08-19_e1-uninformative-variates/` is an unrelated EEG experiment
   also called "E1". X1's outputs are labelled `x1_*` throughout.
