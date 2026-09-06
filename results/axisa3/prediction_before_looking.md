# Directional predictions, recorded before seeing any descriptive result

**Timestamped by its commit, which precedes the commit of
`descriptive_structure.csv`.** The acquisition subagent is running as this is
written; I have seen the *file inventory* of the pgmpy sdist (24 `.bif.gz`
benchmark networks and 9 applied-paper `.txt` DAGs) and nothing else. No node
count, no edge count, no CPDAG, no component size, no separation.

The full pre-registration (`preregistration.md`) is written after Phase 2, since
it must fix the network list and that list depends on what parsed. This file
exists so that the *directional* predictions cannot be tuned to the data, which
is the part a post-Phase-2 pre-registration could not protect.

## What I predict, and why

**On undirectedness (H13).** Benchmark Bayesian networks are largely
expert-elicited and I expect them to be collider-rich, which compels edges. So I
predict the **undirected fraction is low to moderate** — median across networks
somewhere in **15–45%** — rather than the near-total undirectedness a sparse
random DAG gives. I expect wide variance across networks and I expect the large
sparse ones (MUNIN, LINK, PIGS, DIABETES) to differ from the small dense ones.

**On component size.** I predict real networks have **larger maximum chain
components than the synthetic families**, where session 4's largest at n = 24 was
6 vertices and session 5's random ensembles had 92.8% of instances at sizes 2–6.
Concretely: **more than half the networks will have a maximum component larger
than 6.** The reason is scale — these networks have 8 to 1,041 nodes, and a
component only needs a v-structure-free neighbourhood to grow.

**On separation (H11) — the session's actual question.** I predict separation
≥ 2 is **more common in real networks than in Erdős–Rényi, but still a
minority**. Erdős–Rényi gave 85.8% at s = 1 with a maximum of 4. I predict
**between 15% and 40% of measurable pairs at s ≥ 2**, and a maximum separation
above 4 in at least one network.

**On the law (H12).** `r_val = min(s, |K_{G₀}|)` was a law of a designed family,
close to true by construction. I predict it **does not hold universally** on real
structure — I expect agreement **above 70% but below 100%**, with exceptions
arising where the adjustment set has more than one member in the component, or
where a second route to failure exists that the spine construction excluded.
I expect the exceptions to be the interesting part.

## What would embarrass each prediction

- Undirected fraction near 0% everywhere: the framework would apply narrowly and
  the project's premise would be in trouble.
- Undirected fraction near 100% everywhere: my collider-rich reasoning is wrong.
- Separation ≥ 2 at Erdős–Rényi rates (~14%): the factual question resolves
  negative, and the honest claim changes again.
- The law holding at exactly 100%: I would suspect the measurement of inheriting
  the construction rather than testing it, and would look for why.
