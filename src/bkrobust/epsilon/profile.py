"""The bias profile ``beta_up(d)``, the radius ``r_eps``, and the bounds on it.

Everything here is a consequence of two facts proved in
``docs/R_EPSILON_THEORY.md`` and of nothing else:

* ``B`` is monotone under model inclusion (Theorem B), so ``{B > eps}`` is an
  up-set and the universal retraction theorem (Theorem U) applies to it exactly
  as Proposition 1 of the paper applies to the validity failure set. The search
  is therefore confined to retractions of the analyst's own commitments.
* The up-set of ``G0`` is parameterised by subsets of ``K_{G0}`` (Proposition R),
  and shell ``d`` is generated **directly** from the ``d``-subsets, with no
  traversal of shells below it (Theorem C.4).

Four consequences shape the code:

1. **Shells below ``r_val`` are never evaluated.** Theorem A puts ``B`` at
   exactly zero there. The expensive part of an ``r_eps`` computation is the
   ``B`` evaluation, and it is skipped on the whole certified ball -- which is
   the operational content of "no need to compute estimates inside this shell".
2. **Only the frontier shell is evaluated.** ``beta_up(d)`` is the maximum over
   the sphere at ``d``, not over the ball (Theorem C.3).
3. **One traversal answers every ``eps``.** The staircase is computed once and
   every radius is a lookup on it -- including for an ``eps`` chosen afterwards.
   There is no per-``eps`` search.
4. **Monotonicity licenses bisection.** ``beta_up`` is non-decreasing, so the
   first crossing can be bracketed by binary search over ``d``, evaluating
   ``O(log |K_{G0}|)`` shells instead of all of them.

The two strategies are exact and must agree; they differ only in which shells
they touch, and :func:`r_epsilon` dispatches between them the way
:mod:`bkrobust.hybrid` dispatches between the search and the ladder -- on a
discriminator that is free at runtime.

**Every radius here inherits the assumption chain of ``r_val`` and adds nothing
to it.** The retraction reduction rests on Lemma L of ``THEOREMS.md``, whose
Anti-Exchange Case B premise is verified rather than proved; the error is
one-sided in the same direction (a radius can come back too large, i.e.
overstating robustness, never too small). That caveat travels on
:attr:`EpsilonResult.assumes`.
"""

from __future__ import annotations

import itertools
import time
from collections.abc import Sequence
from dataclasses import dataclass, field

from bkrobust.core.conventions import UNREACHED
from bkrobust.demo.graph import MPDAG
from bkrobust.demo.meek import apply_orientations
from bkrobust.epsilon.bias import BiasAt, BiasContext, bias_at, knowledge_of

Edge = tuple[str, str]

#: Refuse to build a shell with more subsets than this. ``C(|K_{G0}|, d)`` peaks
#: at ``d = |K_{G0}| / 2``; at ``|K_{G0}| = 24`` -- the largest the paper's
#: envelope reports (Table 3) -- that peak is 2.7M closures. The cap turns a
#: hang into a reported status, which is this project's house rule for anything
#: that could otherwise be silently truncated.
DEFAULT_SHELL_CAP: int = 200_000

#: Shells the incremental strategy explores before handing over to bisection.
#: The paper measures ``r_val`` at 1 for ~80% of synthetic instances and <= 3
#: almost always, and ``r_eps >= r_val``, so a small budget catches the regime
#: the incremental pass is good at while bounding the work it can be forced into.
DEFAULT_SHELL_BUDGET: int = 4


# --- shells ---------------------------------------------------------------


def retraction_shell(
    cpdag: MPDAG,
    g0: MPDAG,
    d: int,
    *,
    cap: int = DEFAULT_SHELL_CAP,
) -> list[MPDAG]:
    """The states reachable by retracting exactly ``d`` of ``G0``'s orientations.

    By Corollary R1 this family sits between the sphere at ``d`` and the ball at
    ``d``, which is exactly what makes its maximum of ``B`` equal ``beta_up(d)``
    (Theorem C.4). It is built straight from the ``d``-subsets of ``K_{G0}``:
    **no shell below ``d`` is visited.**

    Distinct subsets frequently share a closure -- Meek's rules re-derive some of
    what was withdrawn -- so the result is deduplicated, and is typically far
    smaller than ``C(|K_{G0}|, d)``.

    Args:
        cpdag: The CPDAG.
        g0: The analyst's state.
        d: How many orientations to retract.
        cap: Refuse to enumerate more than this many subsets.

    Returns:
        The distinct closures, ordered by ``(number of directed edges, edge
        string)`` for determinism.

    Raises:
        ValueError: If ``d`` is negative or exceeds ``|K_{G0}|``.
        MemoryError: If the subset count exceeds ``cap``. Raised rather than
            truncated: a silently partial shell would produce a ``beta_up``
            that is too small and hence a radius that is too large, in the
            direction that overstates robustness.
    """
    k0 = knowledge_of(cpdag, g0)
    if d < 0 or d > len(k0):
        raise ValueError(f"d must be in 0..{len(k0)}, got {d}")

    n_subsets = _binomial(len(k0), d)
    if n_subsets > cap:
        raise MemoryError(
            f"shell {d} of a {len(k0)}-orientation knowledge set needs {n_subsets} "
            f"subsets, above the cap of {cap}"
        )

    seen: dict[str, MPDAG] = {}
    for drop in itertools.combinations(range(len(k0)), d):
        keep = [k0[i] for i in range(len(k0)) if i not in set(drop)]
        h = apply_orientations(cpdag, keep)
        if h is None:  # unreachable: a subset of a consistent set is consistent
            continue
        seen.setdefault(h.edge_string(), h)
    return sorted(seen.values(), key=lambda g: (len(g.directed_edges), g.edge_string()))


def _binomial(n: int, k: int) -> int:
    """``C(n, k)`` without importing math.comb's overflow behaviour into a guard."""
    if k < 0 or k > n:
        return 0
    result = 1
    for i in range(k):
        result = result * (n - i) // (i + 1)
    return result


# --- the staircase --------------------------------------------------------


@dataclass(frozen=True)
class ShellBias:
    """One step of the bias staircase.

    Attributes:
        d: Retraction depth.
        n_subsets: ``C(|K_{G0}|, d)``, the subsets enumerated.
        n_states: Distinct closures among them -- the real width of the shell.
        shell_max: ``max B`` over this shell alone, **before** the running
            maximum is taken. Kept separately so that Theorem C.3 (the shell
            maximum already is the ball maximum, hence is itself non-decreasing)
            is falsifiable from the recorded output rather than enforced by it.
        beta_up: ``beta_up(d)``, the running maximum through depth ``d``.
        witness: Edge string of a state attaining ``shell_max``.
        detail: The full :class:`~bkrobust.epsilon.bias.BiasAt` at the witness.
        seconds: Wall time for this shell.
    """

    d: int
    n_subsets: int
    n_states: int
    shell_max: float
    beta_up: float
    witness: str | None
    detail: BiasAt | None
    seconds: float


def shell_bias(
    ctx: BiasContext,
    g0: MPDAG,
    d: int,
    *,
    running_max: float = 0.0,
    method: str = "semilocal",
    cap: int = DEFAULT_SHELL_CAP,
) -> ShellBias:
    """Evaluate ``B`` across one shell and return the step of the staircase.

    Args:
        ctx: The per-instance bias context.
        g0: The analyst's state.
        d: Retraction depth.
        running_max: ``beta_up(d - 1)``, folded in so the returned ``beta_up``
            is the ball maximum even if a shell were to come back lower (which
            Theorem C.3 says cannot happen, and which ``shell_max`` records
            independently so that the prediction stays checkable).
        method: Passed through to :func:`~bkrobust.epsilon.bias.bias_at`.
        cap: Passed to :func:`retraction_shell`.

    Returns:
        A :class:`ShellBias`.
    """
    t0 = time.perf_counter()
    states = retraction_shell(ctx.cpdag, g0, d, cap=cap)

    best = -1.0
    best_state: MPDAG | None = None
    best_detail: BiasAt | None = None
    for h in states:
        value = bias_at(ctx, h, method=method)
        if value.worst > best:
            best, best_state, best_detail = value.worst, h, value

    shell_max = best if best_state is not None else 0.0
    return ShellBias(
        d=d,
        n_subsets=_binomial(len(knowledge_of(ctx.cpdag, g0)), d),
        n_states=len(states),
        shell_max=float(shell_max),
        beta_up=float(max(running_max, shell_max)),
        witness=best_state.edge_string() if best_state is not None else None,
        detail=best_detail,
        seconds=time.perf_counter() - t0,
    )


def bias_profile(
    ctx: BiasContext,
    g0: MPDAG,
    *,
    start: int = 0,
    max_d: int | None = None,
    stop_above: float | None = None,
    method: str = "semilocal",
    cap: int = DEFAULT_SHELL_CAP,
) -> list[ShellBias]:
    """The staircase ``beta_up(start), ..., beta_up(max_d)``.

    Shells strictly below ``start`` are reported as exactly zero without any
    ``B`` evaluation. Callers pass ``start = r_val``, which Theorem A licenses:
    the whole certified ball is algebraically clean, so evaluating it would be
    computing a number already known to be zero.

    Args:
        ctx: The per-instance bias context.
        g0: The analyst's state.
        start: First depth actually evaluated. Depths below it are recorded with
            ``beta_up = 0.0``.
        max_d: Last depth to evaluate; defaults to ``|K_{G0}|``, the top of the
            space (``Meek(Chat, {}) = Chat``).
        stop_above: Stop once ``beta_up`` exceeds this. Passing the largest
            ``eps`` of interest makes the traversal stop at ``r_eps`` for that
            ``eps``, and hence at or before ``r_eps`` for every smaller one.
        method: Passed through.
        cap: Passed through.

    Returns:
        One :class:`ShellBias` per depth from 0 up to where the traversal
        stopped, in increasing depth order.

    Raises:
        ValueError: If ``start`` is negative.
    """
    if start < 0:
        raise ValueError(f"start must be non-negative, got {start}")
    k0 = knowledge_of(ctx.cpdag, g0)
    top = len(k0) if max_d is None else min(max_d, len(k0))

    out: list[ShellBias] = []
    for d in range(0, min(start, top + 1)):
        out.append(
            ShellBias(
                d=d,
                n_subsets=_binomial(len(k0), d),
                n_states=0,
                shell_max=0.0,
                beta_up=0.0,
                witness=None,
                detail=None,
                seconds=0.0,
            )
        )

    running = 0.0
    for d in range(start, top + 1):
        step = shell_bias(ctx, g0, d, running_max=running, method=method, cap=cap)
        out.append(step)
        running = step.beta_up
        if stop_above is not None and running > stop_above:
            break
    return out


def r_epsilon_from_profile(profile: Sequence[ShellBias], eps: float) -> int:
    """Read ``r_eps`` off a computed staircase.

    The reason the staircase is the right object: this is a lookup, so a whole
    grid of ``eps`` values -- including ones chosen after the traversal -- costs
    nothing beyond the single traversal.

    Args:
        profile: Output of :func:`bias_profile`, in increasing depth order.
        eps: The bias threshold, in whatever units the profile was built in.

    Returns:
        The smallest depth whose ``beta_up`` exceeds ``eps``, or
        :data:`~bkrobust.core.conventions.UNREACHED` if none does. ``UNREACHED``
        is a sentinel, not a number: it means "not found within the traversal",
        and is only a proof of ``r_eps = infinity`` when the traversal reached
        the top of the space.
    """
    for step in profile:
        if step.beta_up > eps:
            return step.d
    return UNREACHED


# --- the radius, and the bounds on it -------------------------------------


@dataclass
class EpsilonResult:
    """A radius ``r_eps``, how it was obtained, and what it assumes.

    Attributes:
        radius: The depth at which the worst-case bias first exceeds ``eps``, or
            :data:`~bkrobust.core.conventions.UNREACHED`.
        eps: The threshold this radius is for.
        units: ``"absolute"``, ``"relative"`` or ``"standardised"``.
        exact: False when the traversal stopped on a budget rather than on an
            answer, in which case ``lower_bound``/``upper_bound`` are what is
            known.
        lower_bound: A certified lower bound on ``r_eps`` (no state at a smaller
            depth exceeds ``eps``).
        upper_bound: A certified upper bound, or ``UNREACHED`` if none was
            established.
        strategy: ``"incremental"``, ``"bisection"``, ``"greedy_chain"`` or
            ``"top_only"``.
        beta_at_radius: ``beta_up`` at ``radius``, i.e. how far past ``eps`` the
            first crossing actually goes -- the staircase is a step function and
            the overshoot is worth reporting next to the radius.
        witness: Edge string of a state attaining the crossing.
        shells_evaluated: How many shells had ``B`` computed on them.
        states_evaluated: How many ``B`` evaluations were performed.
        seconds: Wall time.
        assumes: The inherited assumption chain.
    """

    radius: int
    eps: float
    units: str = "absolute"
    exact: bool = True
    lower_bound: int = 0
    upper_bound: int = UNREACHED
    strategy: str = "incremental"
    beta_at_radius: float | None = None
    witness: str | None = None
    shells_evaluated: int = 0
    states_evaluated: int = 0
    seconds: float = 0.0
    assumes: str = (
        "Lemma L of THEOREMS.md (Anti-Exchange Case B, proved in THEOREMS.md "
        "section 4), exactly as r_val does"
    )
    profile: list[ShellBias] = field(default_factory=list)


def top_bias(ctx: BiasContext, *, method: str = "semilocal") -> BiasAt:
    """``B(Chat)``: the largest bias any perturbation whatsoever can produce.

    ``Chat`` is the maximum of the order, so by Theorem B this dominates ``B``
    everywhere. One evaluation therefore settles ``r_eps = infinity`` for every
    ``eps >= B(Chat)`` -- the bias analogue of the paper's "no failure anywhere"
    case, and like it the case an exhaustive search handles worst.

    Args:
        ctx: The per-instance bias context.
        method: Passed through.

    Returns:
        The :class:`~bkrobust.epsilon.bias.BiasAt` at the CPDAG.
    """
    return bias_at(ctx, ctx.cpdag, method=method)


def greedy_chain_bound(
    ctx: BiasContext,
    g0: MPDAG,
    eps: float,
    *,
    method: str = "semilocal",
) -> EpsilonResult:
    """An anytime **upper** bound on ``r_eps`` from a single monotone chain.

    Walks one retraction chain, at each step withdrawing whichever single
    orientation maximises ``B``, and stops at the first state exceeding ``eps``.
    Such a state is a witness, so its depth bounds ``r_eps`` from above by
    Definition 5 -- with **no shell enumerated at all**, at a cost of
    ``O(|K_{G0}|^2)`` closures rather than ``2 ** |K_{G0}|``.

    Depth is measured as the rank drop ``|K_{G0}| - |K_G|``, not as the number of
    greedy steps: one retraction can cascade through Meek's rules and withdraw
    several orientations at once, and the distance is the former.

    Args:
        ctx: The per-instance bias context.
        g0: The analyst's state.
        eps: The threshold.
        method: Passed through.

    Returns:
        An :class:`EpsilonResult` with ``strategy="greedy_chain"``,
        ``exact=False``, and ``upper_bound`` set. ``radius`` carries the same
        value as ``upper_bound``: it is a bound, not a computed radius, and the
        ``exact`` flag is what separates them.
    """
    t0 = time.perf_counter()
    k0 = knowledge_of(ctx.cpdag, g0)
    cur, cur_k = g0, list(k0)
    evaluated = 0

    while cur_k:
        best_value, best_state, best_k = -1.0, None, None
        for e in cur_k:
            keep = [o for o in cur_k if o != e]
            h = apply_orientations(ctx.cpdag, keep)
            if h is None:
                continue
            evaluated += 1
            value = bias_at(ctx, h, method=method).worst
            if value > best_value:
                best_value, best_state, best_k = value, h, knowledge_of(ctx.cpdag, h)
        if best_state is None:
            break
        cur, cur_k = best_state, best_k or []
        depth = len(k0) - len(cur_k)
        if best_value > eps:
            return EpsilonResult(
                radius=depth,
                eps=eps,
                exact=False,
                lower_bound=0,
                upper_bound=depth,
                strategy="greedy_chain",
                beta_at_radius=best_value,
                witness=cur.edge_string(),
                shells_evaluated=0,
                states_evaluated=evaluated,
                seconds=time.perf_counter() - t0,
            )

    return EpsilonResult(
        radius=UNREACHED,
        eps=eps,
        exact=False,
        lower_bound=0,
        upper_bound=UNREACHED,
        strategy="greedy_chain",
        states_evaluated=evaluated,
        seconds=time.perf_counter() - t0,
    )


def r_epsilon(
    ctx: BiasContext,
    g0: MPDAG,
    eps: float,
    *,
    r_val: int = 0,
    strategy: str = "hybrid",
    shell_budget: int = DEFAULT_SHELL_BUDGET,
    method: str = "semilocal",
    cap: int = DEFAULT_SHELL_CAP,
) -> EpsilonResult:
    """The exact ``r_eps``, by incremental traversal, bisection, or a dispatch.

    Both strategies are exact and must return the same radius; they differ only
    in which shells they touch.

    * ``"incremental"`` evaluates shells ``r_val, r_val + 1, ...`` and stops at
      the first crossing. Cost is the sum of the shells up to ``r_eps``: cheap
      when the crossing is near, which the paper's measurements say is usual.
    * ``"bisection"`` binary-searches ``d`` in ``[r_val, |K_{G0}|]``, building
      each probed shell directly by Theorem C.4. ``O(log |K_{G0}|)`` shells, but
      a probe near the middle is the widest shell there is.
    * ``"hybrid"`` runs the incremental pass under ``shell_budget`` and switches
      to bisection if it has not crossed -- the same discriminator
      :mod:`bkrobust.hybrid` uses, and free at runtime.

    Args:
        ctx: The per-instance bias context.
        g0: The analyst's state.
        eps: The threshold, in ``ctx``'s effect units.
        r_val: The validity radius. Shells below it are skipped outright
            (Theorem A) and it is reported as the certified ``lower_bound``.
            Passing 0 is always safe and merely forgoes the saving.
        strategy: As above.
        shell_budget: Shells the incremental leg explores before handing over.
        method: Passed through.
        cap: Passed to :func:`retraction_shell`.

    Returns:
        An :class:`EpsilonResult`.

    Raises:
        ValueError: On an unknown strategy.
    """
    if strategy not in {"incremental", "bisection", "hybrid"}:
        raise ValueError(f"unknown strategy {strategy!r}")

    t0 = time.perf_counter()
    k0 = knowledge_of(ctx.cpdag, g0)
    start = max(0, r_val)

    # The infinity test (Bound U2): Chat dominates B everywhere, so one
    # evaluation at the top settles the whole question for large eps.
    top = top_bias(ctx, method=method)
    if top.worst <= eps:
        return EpsilonResult(
            radius=UNREACHED,
            eps=eps,
            exact=True,
            lower_bound=len(k0) + 1,
            upper_bound=UNREACHED,
            strategy="top_only",
            beta_at_radius=top.worst,
            shells_evaluated=1,
            states_evaluated=1,
            seconds=time.perf_counter() - t0,
        )

    if strategy in {"incremental", "hybrid"}:
        limit = len(k0) if strategy == "incremental" else min(len(k0), start + shell_budget)
        profile = bias_profile(
            ctx, g0, start=start, max_d=limit, stop_above=eps, method=method, cap=cap
        )
        hit = r_epsilon_from_profile(profile, eps)
        evaluated_shells = sum(1 for s in profile if s.n_states > 0)
        evaluated_states = sum(s.n_states for s in profile)
        if hit != UNREACHED:
            step = next(s for s in profile if s.d == hit)
            return EpsilonResult(
                radius=hit,
                eps=eps,
                exact=True,
                lower_bound=hit,
                upper_bound=hit,
                strategy="incremental",
                beta_at_radius=step.beta_up,
                witness=step.witness,
                shells_evaluated=evaluated_shells,
                states_evaluated=evaluated_states,
                seconds=time.perf_counter() - t0,
                profile=list(profile),
            )
        if strategy == "incremental":
            # The traversal reached the top without crossing, yet top_bias said
            # it would; that is a contradiction and must not be reported as a
            # radius.
            raise AssertionError(
                "the incremental traversal exhausted the up-set without crossing eps, "
                "but B(Chat) > eps -- the profile and the top evaluation disagree"
            )
        start = limit + 1  # hybrid falls through with a tightened lower bound

    # Bisection. Invariant: beta_up(lo - 1) <= eps < beta_up(hi).
    lo, hi = max(start, 0), len(k0)
    shells, states = 1, 1  # the top evaluation already counted
    best_step: ShellBias | None = None
    probes: list[ShellBias] = []
    while lo < hi:
        mid = (lo + hi) // 2
        step = shell_bias(ctx, g0, mid, method=method, cap=cap)
        probes.append(step)
        shells += 1
        states += step.n_states
        if step.shell_max > eps:
            hi, best_step = mid, step
        else:
            lo = mid + 1

    if best_step is None or best_step.d != lo:
        best_step = shell_bias(ctx, g0, lo, method=method, cap=cap)
        probes.append(best_step)
        shells += 1
        states += best_step.n_states

    return EpsilonResult(
        radius=lo,
        eps=eps,
        exact=True,
        lower_bound=lo,
        upper_bound=lo,
        strategy="bisection" if strategy != "hybrid" else "hybrid_bisection",
        beta_at_radius=best_step.shell_max,
        witness=best_step.witness,
        shells_evaluated=shells,
        states_evaluated=states,
        seconds=time.perf_counter() - t0,
        profile=sorted(probes, key=lambda s: s.d),
    )


def epsilon_grid(
    ctx: BiasContext,
    g0: MPDAG,
    epsilons: Sequence[float],
    *,
    r_val: int = 0,
    method: str = "semilocal",
    cap: int = DEFAULT_SHELL_CAP,
) -> tuple[dict[float, int], list[ShellBias]]:
    """Every ``r_eps`` on a grid, from **one** traversal.

    This is the operational form of Theorem C.4: the staircase is built once, up
    to the largest ``eps`` on the grid, and each radius is then a lookup. The
    marginal cost of an extra ``eps`` is zero, which is what makes the interval
    certificate of Theorem D -- which needs at least two radii -- free.

    Args:
        ctx: The per-instance bias context.
        g0: The analyst's state.
        epsilons: The thresholds, in ``ctx``'s effect units. Order is irrelevant.
        r_val: The validity radius; shells below it are skipped.
        method: Passed through.
        cap: Passed through.

    Returns:
        ``(radii, profile)``. ``radii`` maps each ``eps`` to its radius or to
        :data:`~bkrobust.core.conventions.UNREACHED`; ``profile`` is the
        staircase the lookups were made on, so a reader can check them.
    """
    if not epsilons:
        return {}, []
    largest = max(epsilons)
    profile = bias_profile(ctx, g0, start=max(0, r_val), stop_above=largest, method=method, cap=cap)
    return {e: r_epsilon_from_profile(profile, e) for e in epsilons}, profile


def certified_band(radii: dict[float, int], k: int) -> tuple[float | None, float | None]:
    """The band ``(eps_x, eps_y]`` that Theorem D pins the bias into at budget ``k``.

    Args:
        radii: Output of :func:`epsilon_grid`.
        k: The analyst's error budget -- how many asserted orientations they are
            willing to concede might be false.

    Returns:
        ``(lower, upper)``. ``lower`` is the largest ``eps`` whose radius is at
        most ``k`` (so a bias above it is demonstrably reachable), or ``None``
        if no threshold is reached. ``upper`` is the smallest ``eps`` whose
        radius exceeds ``k`` (so the bias is certified at or below it), or
        ``None`` if every threshold on the grid is already reachable and the
        grid therefore does not bound the bias from above.
    """
    reached = [e for e, r in sorted(radii.items()) if r != UNREACHED and r <= k]
    unreached = [e for e, r in sorted(radii.items()) if r == UNREACHED or r > k]
    return (max(reached) if reached else None, min(unreached) if unreached else None)
