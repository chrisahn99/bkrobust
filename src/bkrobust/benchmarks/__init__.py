"""Real-world benchmark Bayesian networks: acquisition, parsing, description.

Every experiment so far has run on *synthetic* structure -- Erdos-Renyi
ensembles and hand-designed spine families -- where the chain-component
structure is a consequence of the generator's parameters. The question this
package exists to answer is whether the structural facts those experiments
established survive contact with structure nobody designed for us: the
expert-elicited Bayesian networks (ALARM, MUNIN, HEPAR2, ...) that the causal
discovery literature has used as ground truth for thirty years.

The package is deliberately split by concern:

* :mod:`.acquire` -- get the files, and pin exactly which bytes we used.
* :mod:`.bif`, :mod:`.dagitty`, :mod:`.bnjson` -- parse the three on-disk
  formats into a graph structure. Structure only; no CPDs.
* :mod:`.describe` -- the descriptive sweep over the parsed networks.

The governing rule, stated here because it is the whole point of
:mod:`.acquire` recording hashes: **no network is ever reconstructed from
memory**. Every node and every edge comes from a fetched file whose sha256 is
recorded next to the numbers derived from it. A hand-typed adjacency list that
resembles ALARM would be indistinguishable, in a results table, from the real
thing -- and so would poison every downstream claim.

Runs on Python 3.9 as well as the repository's 3.11 target.
"""

from __future__ import annotations

__all__ = ["acquire", "bif", "bnjson", "dagitty", "describe"]
