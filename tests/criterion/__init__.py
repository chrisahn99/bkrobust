"""Tests for :mod:`bkrobust.mpdag_criterion`.

The package under test decides adjustment-set validity on an MPDAG graphically,
so its only real specification is the enumerating oracle it replaces. These
tests are therefore mostly differential: exhaustive at n=3, deterministically
sampled at n=4, plus hand-worked unit cases for each component (possibly causal
paths, definite status, possible descendants, amenability, blocking) so a
failure points at a component rather than only at the end-to-end predicate.
"""
