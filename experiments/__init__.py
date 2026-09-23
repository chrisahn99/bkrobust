"""Hydra entrypoints, one per experiment.

Each script does four things and no more: parse its config, set up logging and
seeding, resolve the results directory, and call into :mod:`bkrobust`. No
experiment logic lives here -- if a script grows a loop that computes something,
that computation belongs in the package where it can be tested.
"""

from __future__ import annotations
