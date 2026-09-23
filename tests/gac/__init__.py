"""Tests for :mod:`bkrobust.gac`, the generalised adjustment criterion.

The package makes two claims that a unit test cannot check by inspection, so
both are checked differentially instead:

* at the DAG level, the single-d-separation decision equals the
  path-enumerating reference (which uses none of the reduction the fast one is
  derived from), and it accepts a strict superset of the back-door criterion;
* at the MPDAG level, the polynomial criterion equals its own semantics --
  "GAC-valid in every DAG extension" -- exhaustively over the same graph scope
  session 4's back-door criterion was verified on.

There is also a performance test, because the criterion is on a hot path and
correctness alone was what let session 4's first criterion ship exponential.
"""
