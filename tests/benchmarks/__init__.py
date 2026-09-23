"""Tests for the benchmark-network acquisition and descriptive-structure module.

These tests guard the one property that cannot be recovered after the fact: that
the networks in the results table are the networks in the files. They therefore
cross-check the BIF parser against an independent counting path, assert every
parsed benchmark really is a DAG, and assert that parsing is deterministic.

They run offline against the cache written by
:func:`bkrobust.benchmarks.acquire.acquire`, and skip cleanly when that cache
has not been produced.
"""
