"""Graph representations and operations on DAGs, CPDAGs and MPDAGs.

An MPDAG (maximally oriented PDAG) is what you get by taking the CPDAG output
by a discovery algorithm, imposing background knowledge, and closing under
Meek's rules. Every structural claim in this project is a statement about how
the MPDAG -- and the optimal adjustment set read off it -- moves when the
imposed knowledge is wrong.

Module map:

``mpdag``
    The :class:`~bkrobust.graphs.mpdag.MPDAG` container and the ``Node`` /
    ``Edge`` type aliases used package-wide.
``meek``
    Meek rules R1-R4, the closure, and Meek's Algorithm 1 with FAIL detection.
``consistency``
    Whether a body of background knowledge is consistent with a CPDAG.
``adjustment``
    Valid adjustment sets, b-adjustment, and the construction of ``O*``.
``optimality``
    Asymptotic variance of an adjustment set and the optimality test.
``distances``
    Metrics on graphs and on knowledge: SHD, orientation-flip count, and the
    model-oriented poset distance.
"""

from __future__ import annotations

__all__: list[str] = []
