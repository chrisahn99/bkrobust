"""THE NORMATIVE RADIUS CONVENTION. Read this before quoting any number.

The report and the project's framing document disagreed by one. This module
settles it; everything downstream obeys the definition here, and
``tests/core/test_conventions.py`` asserts it.

Normative definition
--------------------

For an adjustment set ``Z`` fixed once from the analyst's graph ``G0``, and a
property ``P`` (validity, optimality, or bias-below-epsilon):

    r_P  =  min { d(G0, G) : G in the space, P fails at G }

where ``d`` is BFS hop count on the covering-relation neighbour graph. That is:

    **the radius is the DISTANCE TO THE NEAREST FAILURE.**

Consequences, stated explicitly because the alternative reading is natural:

* Shells ``0, 1, ..., r-1`` are certified clean. Shell ``r`` contains at least
  one failure.
* The practitioner's phrasing -- "how many atomic moves may I make and still be
  safe?" -- is therefore ``r - 1``, **not** ``r``. Anyone quoting the
  safe-moves number must subtract one and say so.
* ``r >= 1`` whenever ``P`` holds at ``G0``. Since ``Z`` is read off ``G0``,
  validity holds there by construction, so ``r_val >= 1`` on any non-degenerate
  instance.
* ``r = 0`` means ``P`` already fails at ``G0`` itself. For validity this is a
  degenerate instance and is gated out at generation time, not reported as a
  radius of zero.
* ``r = UNREACHED`` means no element of the fully enumerated space fails. This
  is NOT a large radius and must never be aggregated as a number -- see
  :data:`UNREACHED`.

Why this convention and not the other
-------------------------------------

1. It is a minimum over a set, which is a clean mathematical object; the
   safe-moves reading is that object minus one, and is easy to derive from it.
2. Axis B's admissible lower bounds are naturally of the form ``r >= L``,
   certifying shells ``0 .. L-1`` clean without visiting them. That composes
   directly with this convention and awkwardly with the other.
3. It is what the existing implementation and ``report.md`` already do, so no
   previously published number needs restating. Had the other convention been
   chosen, every number in the prior report would have shifted by one.

Prior numbers therefore carry over unchanged: ``r_val = 3, 3, 2`` for the three
scenarios of the worked example means shells 0-2 clean with the first failure at
shell 3 (scenarios A and B), and shells 0-1 clean with the first failure at
shell 2 (scenario C).
"""

from __future__ import annotations

#: No element of the enumerated space fails the property.
#:
#: Not a radius. Never average it, never plot it on a numeric axis, never
#: compare it with ``<`` against a real radius. Downstream code must branch on
#: it explicitly.
UNREACHED: int = -1

#: The property already fails at ``G0``. Degenerate; gated out at generation.
DEGENERATE: int = 0

#: Human-readable statement of the convention, embedded in every manifest so a
#: results file is self-describing.
RADIUS_CONVENTION: str = (
    "r = min distance from G0 to a graph where the property fails; "
    "shells 0..r-1 are certified clean; safe-moves = r - 1"
)


def safe_moves(radius: int) -> int:
    """Convert a radius to the practitioner's "how many moves am I safe for".

    Args:
        radius: A radius under the normative convention.

    Returns:
        ``radius - 1``.

    Raises:
        ValueError: If ``radius`` is :data:`UNREACHED`, where the question has
            no finite answer, or negative for any other reason.
    """
    if radius == UNREACHED:
        raise ValueError("safe_moves is undefined for UNREACHED")
    if radius < 0:
        raise ValueError(f"not a radius: {radius}")
    return radius - 1


def is_certified_clean(radius: int, shell: int) -> bool:
    """Whether the convention certifies ``shell`` to contain no failure.

    Args:
        radius: A radius under the normative convention.
        shell: A shell index.

    Returns:
        ``True`` if every element at ``shell`` is known clean. With
        :data:`UNREACHED` every shell is clean, since nothing in the space
        fails.
    """
    if radius == UNREACHED:
        return True
    return shell < radius
