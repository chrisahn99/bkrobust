# Data

Provenance, licences, and download instructions for every dataset.

**Nothing in this repository downloads data automatically.** Loaders read from
`data/` (gitignored) and raise `FileNotFoundError` pointing here when a file is
absent. Several of these datasets have redistribution terms that a convenience
auto-download would quietly violate, and a registration wall is not an obstacle
to route around.

`scripts/download_data.sh` is a stub. Fill it in per dataset as each is
actually used, and record the retrieval date in the table below.

---

## Status

| Dataset | Tier | Ground truth | Licence | Obtained | Notes |
|---|---|---|---|---|---|
| Sachs | real | consensus network | TODO | ☐ | 11 nodes; observational + interventional |
| DREAM | semi-synthetic | gold-standard network | TODO | ☐ | sizes 10/50/100 |
| RealCause | semi-synthetic | known simulated effect | TODO | ☐ | fitted generators |
| IHDP | real covariates | simulated outcome | TODO | ☐ | 1000 replications |
| ACIC | real covariates | simulated outcome | TODO | ☐ | 2016/2017/2018 settings |
| Twins | real | both potential outcomes | TODO | ☐ | true ITE observed |
| Jobs | real | randomised arm | TODO | ☐ | LaLonde-derived |
| LaLonde | real | experimental benchmark | TODO | ☐ | PSID/CPS controls |
| Perturb-seq | real | partial (interventions) | TODO | ☐ | very large; subset required |

Fill in each licence **before** any result from that dataset appears in the
paper. "TODO" in a licence column at submission time is a defect.

---

## Layout

```
data/
├── sachs/
├── dream/
├── realcause/
├── ihdp/
├── acic/
├── twins/
├── jobs/
├── lalonde/
└── perturb_seq/
```

`data/` is gitignored in full. No dataset is ever committed, whatever its size
or licence.

---

## Per-dataset notes

### Sachs

Protein-signalling network, 11 phosphoproteins, flow cytometry. The standard
real-data structure benchmark.

**The consensus network is not ground truth.** It is expert belief that has held
up under scrutiny, which is a different epistemic object, and this project of all
projects should not blur that line — the consensus network *is itself* background
knowledge, elicited from experts and never verified against a known truth.
Results should say "consensus network", not "true DAG".

Preprocessing: log-transform by default (raw abundances are heavily skewed).
Standardisation is off by default — it destroys the variance information that
varsortability diagnostics read, hiding a leak rather than removing it.

### DREAM

Gene-regulatory network challenges with published gold standards and simulators.
DREAM4 size-10 networks are small enough for exact radius computation; size-50
and above need the heuristic path.

### RealCause

Generative models fitted to real causal benchmarks, sampled to give realistic
marginals with a known effect. The caveat that must ride along: the "true"
mechanism is the fitted generator's, not the real process's.

### IHDP / ACIC

Real covariates, simulated outcomes. Effect known by construction. IHDP is
1000 replications of the same covariates with different outcome draws — treat
replications as one dataset with 1000 seeds, not as 1000 datasets, or the
effective sample size of any aggregate is overstated.

### Twins

Both potential outcomes observed (one twin each), so the individual effect is
known without simulation. Rare, and the best PEHE benchmark available here.

### Jobs / LaLonde

Randomised and observational arms, so a benchmark effect exists. The
observational control choice (PSID vs CPS) changes the difficulty substantially
and must be reported.

### Perturb-seq

Single-cell expression with genetic perturbations. The interventional arm gives
partial structural ground truth — knocking out a gene reveals its downstream
targets — which is the closest thing here to a verified graph. Far too large for
exact radius computation; a gene subset is mandatory and the selection criterion
must be reported, since selecting genes by a criterion related to the graph
would bias the result.

---

## Leakage

Every synthetic and semi-synthetic result carries its varsortability and
R²-sortability. See `bkrobust/data/leakage.py` and Reisach et al. on why
standard synthetic benchmarks make a trivial variance-sorting baseline
competitive with published methods.

The concern specific to this project: if data is highly varsortable, discovery
recovers the causal order from marginal variances regardless of what the
background knowledge says. The knowledge never gets to matter, the radii look
large, and the result is a benchmark artefact reported as robustness.

Responses, in decreasing order of preference:

1. Sample noise scales to break the pattern — this is why `noise_scale` is a
   range and not a constant.
2. Report results split by varsortability.
3. Standardise. Worst option: it hides the leak from the diagnostic without
   removing the information, and changes the estimand's scale.

---

## Variable glossaries

`data/loaders.py::variable_glossary` returns human-readable variable
descriptions. These are not documentation — they are an experimental condition.
An LLM asked to orient `praf` against `pmek` performs differently depending on
what it was told those are, so the glossary used goes into the run manifest
alongside the prompt template.
