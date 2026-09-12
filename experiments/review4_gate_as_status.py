"""[RE-3]: measure the population the admission gate throws away.

`fast_gate` admits a pair only if retracting one claim from the recovering set
breaks validity. At coverage 1.0 that test *is* the definition of `r_claim == 1`,
so the 1,659 pairs it rejects as `no_atomic_perturbation_changes_validity` are
exactly the `r_claim >= 2` population -- the only rows on which the claim radius
could vary -- and every published claim-radius number is a constant by
construction because of it.

Those rows already passed the four structural conditions before being rejected, so
completing them needs no re-screening of the 115,974-pair corpus: we read them back
out of the committed `instances.jsonl` and run the rest of `measure_pair` on each.

Writes `results/axisa3/instances_gate_as_status.jsonl` **alongside** the committed
file, never over it, so previously published numbers stay reproducible from the code
that produced them.

Guard: for a sample of admitted pairs, `gate_without_atomic_clause` must return the
same verdict as `fast_gate` and must set the status flag True, or the script aborts.
"""
from __future__ import annotations

import hashlib
import json
import pathlib
import signal
import sys
import time
from collections import Counter
from contextlib import contextmanager

#: Hard per-instance timeout. The network wall budget is checked between instances,
#: which cannot stop a single one: `optimal_adjustment_set_mpdag` enumerates DAG
#: extensions with no time limit of its own, and on a 200-node network one instance
#: can run unbounded. These rows are the expensive case by construction, so the guard
#: has to be able to interrupt.
PER_INSTANCE_S = 15


class _Timeout(Exception):
    pass


@contextmanager
def time_budget(seconds: int):
    def _fire(signum, frame):
        raise _Timeout()

    old = signal.signal(signal.SIGALRM, _fire)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from bkrobust.benchmarks.describe import parse_file  # noqa: E402
from bkrobust.benchmarks.measure import (  # noqa: E402
    MAX_G0_UNDIRECTED_FOR_EXTENSIONS,
    O_INTRACTABLE,
    component_of,
    fast_gate,
    gate_without_atomic_clause,
    select_knowledge,
    separation,
)
from bkrobust.demo.example import dag_to_cpdag  # noqa: E402
from bkrobust.demo.evaluate import optimal_adjustment_set_mpdag  # noqa: E402
from bkrobust.demo.graph import MPDAG  # noqa: E402
from bkrobust.demo.meek import apply_orientations  # noqa: E402
from bkrobust.gac import is_gac_valid_mpdag  # noqa: E402
from bkrobust.hybrid import breakdown_radius  # noqa: E402

A3 = ROOT / "results" / "axisa3"
NETDIR = A3 / "networks" / "example_models"
OUT = A3 / "instances_gate_as_status.jsonl"
REASON = "no_atomic_perturbation_changes_validity"
TIME_LIMIT = 20.0  # per instance; exhaustion is recorded as censored, never dropped
#: Per-network cap, by deterministic stride over the sorted rows.
#:
#: These are the rows the gate rejected because no single retraction breaks
#: validity, so they are exactly the instances on which the upward search must
#: exhaust the whole up-set -- the expensive case by construction. Measuring all
#: 1,659 would finish only on the small networks, biasing the result toward them.
#: A stride cap covers all 27 networks instead, and the sampling is reported.
PER_NETWORK_CAP = 15
#: Wall budget per network, mirroring the original sweep's 1,200 s cap. A network that
#: exhausts it contributes the rows it finished and is recorded as censored, never
#: dropped and never silently truncated -- the project's convention for a measurement
#: limit as opposed to a structural rejection.
PER_NETWORK_WALL_S = 60.0
CLAIM_DEPTH = 3  # exhaustive retraction subsets to this depth, then a sentinel
UNREACHED_SENTINEL = -1  # matches the project convention: a status, never a number


def load_network(name: str) -> tuple[MPDAG, MPDAG] | None:
    for path in sorted(NETDIR.iterdir()):
        stem = path.name.split(".")[0]
        if stem == name:
            pn = parse_file(path, hashlib.sha256(path.read_bytes()).hexdigest())
            dag = MPDAG(nodes=set(pn.nodes), directed=set(pn.edges), undirected=set())
            return dag, dag_to_cpdag(dag)
    return None


def claim_radius(cpdag, k, x, y, z, max_depth=CLAIM_DEPTH, budget=20_000):
    """Smallest number of asserted sentences whose retraction breaks validity.

    Returns a negative sentinel when no subset up to ``max_depth`` breaks validity
    (``-max_depth - 1``) or when the subset budget is spent (``-1``, i.e. the same
    UNREACHED status the hop radius uses, since neither is a number).

    The caller skips this entirely when the hop radius is UNREACHED: no element of
    the space fails there, and every claim retraction lands in the space, so the
    claim radius is UNREACHED too and no search can say otherwise.
    """
    import itertools

    spent = 0
    for d in range(1, min(max_depth, len(k)) + 1):
        for drop in itertools.combinations(range(len(k)), d):
            spent += 1
            if spent > budget:
                return UNREACHED_SENTINEL
            kept = [e for i, e in enumerate(k) if i not in drop]
            g = apply_orientations(cpdag, kept)
            if g is None:
                continue
            if not is_gac_valid_mpdag(g, x, y, z):
                return d
    return -max_depth - 1  # sentinel: > max_depth


def main() -> None:
    rows = [json.loads(l) for l in (A3 / "instances.jsonl").open() if '"_meta"' not in l]
    targets = [r for r in rows if r.get("reject_reason") == REASON]
    print(f"{len(targets)} rows rejected as {REASON}")

    by_net: dict[str, list[dict]] = {}
    for r in targets:
        by_net.setdefault(r["network"], []).append(r)

    out, stats, guard_checked = [], Counter(), 0
    censored_networks: dict[str, int] = {}
    t0 = time.perf_counter()
    sink = OUT.open("w", encoding="utf-8")  # incremental: a long run stays inspectable
    for ni, (net, group) in enumerate(sorted(by_net.items()), 1):
        loaded = load_network(net)
        if loaded is None:
            stats["network_not_found"] += len(group)
            print(f"  [{ni}/{len(by_net)}] {net}: NETWORK NOT FOUND, {len(group)} rows skipped")
            continue
        dag, cpdag = loaded

        # --- guard: on pairs fast_gate admits, the statused gate must agree -----
        admitted = ([r for r in rows if r["network"] == net and r.get("admissible")][:2]
                    if guard_checked < 20 else [])
        for r in admitted:
            ok_f, reason_f = fast_gate(dag, cpdag, r["x"], r["y"])
            ok_s, reason_s, flag = gate_without_atomic_clause(dag, cpdag, r["x"], r["y"])
            if not (ok_f and ok_s and flag and reason_f == reason_s):
                raise SystemExit(
                    f"ABORT: gate disagreement on {net} {r['x']}->{r['y']}: "
                    f"fast_gate=({ok_f},{reason_f}) statused=({ok_s},{reason_s},{flag})"
                )
            guard_checked += 1

        group = sorted(group, key=lambda d: (d["coverage"], d["x"], d["y"]))
        if len(group) > PER_NETWORK_CAP:
            step = len(group) / PER_NETWORK_CAP
            group = [group[min(len(group) - 1, int(i * step))] for i in range(PER_NETWORK_CAP)]
            stats["sampled_networks"] += 1
        t_net = time.perf_counter()
        for r in group:
            if time.perf_counter() - t_net > PER_NETWORK_WALL_S:
                stats["censored_network_wall"] += 1
                censored_networks.setdefault(net, 0)
                censored_networks[net] += 1
                continue
            x, y, cov = r["x"], r["y"], r["coverage"]
            try:
                with time_budget(PER_INSTANCE_S):
                    row = measure_one(dag, cpdag, net, x, y, cov, stats)
            except _Timeout:
                stats["censored_instance_wall"] += 1
                continue
            if row is None:
                continue
            stats["ok"] += 1
            out.append(row)
            sink.write(json.dumps(row) + "\n")
            sink.flush()
        print(f"  [{ni}/{len(by_net)}] {net}: {len(group)} rows -> {stats['ok']} total ok "
              f"({time.perf_counter() - t0:.0f}s)", flush=True)

    sink.close()
    (OUT.with_suffix(".meta.json")).write_text(json.dumps({
        "source": "instances.jsonl rows rejected as " + REASON,
        "guard_pairs_checked": guard_checked, "claim_radius_depth_cap": CLAIM_DEPTH,
        "time_limit_s": TIME_LIMIT, "per_network_cap": PER_NETWORK_CAP,
        "per_network_wall_s": PER_NETWORK_WALL_S, "per_instance_s": PER_INSTANCE_S,
        "censored_networks": censored_networks, "outcomes": dict(stats),
        "n_rows": len(out)}, indent=2))
    print(f"\nguard: {guard_checked} admitted pairs agreed between the two gates")
    print(f"outcomes: {dict(stats)}")
    print(f"wrote {OUT} ({len(out)} newly measured rows)")


def measure_one(dag, cpdag, net, x, y, cov, stats):
    """Complete one rejected row: gate status, G0, optimal set, radius, claim radius."""
    ok, reason, flag = gate_without_atomic_clause(dag, cpdag, x, y)
    if not ok or flag:
        stats["unexpected_gate_verdict"] += 1
        return None
    k = select_knowledge(dag, cpdag, cov)
    g0 = apply_orientations(cpdag, k)
    if g0 is None:
        stats["knowledge_inconsistent"] += 1
        return None
    if len(g0.undirected_edges) > MAX_G0_UNDIRECTED_FOR_EXTENSIONS:
        stats[O_INTRACTABLE] += 1
        return None
    o = optimal_adjustment_set_mpdag(g0, x, y)
    if not o:
        stats["o_g0_not_identified"] += 1
        return None
    z = frozenset(o)
    if not is_gac_valid_mpdag(g0, x, y, z):
        stats["z_invalid_at_g0"] += 1
        return None
    comp = component_of(cpdag, x)
    sep, sep_status = separation(cpdag, x, z)
    t = time.perf_counter()
    # These rows were rejected precisely because shell 1 holds no failure, so the
    # bounded upward search has nothing to find and would exhaust the whole up-set --
    # which it does without consulting the time limit, since only the ladder is
    # time-bounded. Budget the search at 1 and let the CP-SAT ladder answer.
    res = breakdown_radius(cpdag, None, x, y, z, g0=g0, search_budget=1,
                           time_limit_s=TIME_LIMIT)
    if not res.exact:
        stats["censored"] += 1
    secs = round(time.perf_counter() - t, 5)
    # If nothing in the space fails, no claim retraction can fail either -- every
    # retraction lands in the space. Skipping the search here is not an approximation,
    # and it is what makes this sweep affordable.
    rc = UNREACHED_SENTINEL if res.radius == UNREACHED_SENTINEL else claim_radius(
        cpdag, k, x, y, z)
    return {
        "network": net, "x": x, "y": y, "coverage": cov, "admissible": True,
        "reject_reason": "ok", "atomic_perturbation_changes_validity": False,
        "component_size": len(comp) if comp else 0,
        "separation": sep, "separation_status": sep_status,
        "n_asserted": len(k),
        "k_g0": sum(1 for (a, b) in g0.directed_edges
                    if tuple(sorted((a, b))) in cpdag.undirected_edges),
        "z_size": len(z), "radius": res.radius, "method": res.method,
        "oracle": res.oracle, "exact": res.exact, "seconds": secs,
        "claim_radius": rc, "claim_radius_depth_cap": CLAIM_DEPTH,
        "g0_undirected_edges": len(g0.undirected_edges),
    }


if __name__ == "__main__":
    main()
