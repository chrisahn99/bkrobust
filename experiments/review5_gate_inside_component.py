"""[RE-3], properly: the discarded rows on which the question can be asked at all.

Round 4 sampled the rejected population by a stride over all 1,659 rows and reported
that every measured row was UNREACHED. That conclusion was an artefact of the sample.
`fast_gate`'s atomic clause is near-collinear with a trivial structural predicate on
this corpus: 1,602 of the 1,659 rejected rows have the treatment OUTSIDE every
undirected chain component, while all 831 admitted rows have it inside one. A stride
sample therefore lands almost surely in the degenerate slice, where no adjustment set
can break because there is nothing to perturb near X, and where the separation baseline
is undefined too.

The rows that carry the question are the 57 with X inside a component. This measures
every one of them -- no sampling, no stride -- and reports the claim-radius distribution
the round-4 write-up wrongly declared absent.
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

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from bkrobust.benchmarks.describe import parse_file  # noqa: E402
from bkrobust.benchmarks.measure import component_of  # noqa: E402
from bkrobust.demo.example import dag_to_cpdag  # noqa: E402
from bkrobust.demo.graph import MPDAG  # noqa: E402

import experiments.review4_gate_as_status as g4  # noqa: E402

A3 = ROOT / "results" / "axisa3"
ND = A3 / "networks" / "example_models"
OUT = A3 / "instances_gate_inside_component.jsonl"
REASON = "no_atomic_perturbation_changes_validity"
PER_INSTANCE_S = 120


class _Timeout(Exception):
    pass


@contextmanager
def budget(seconds):
    def fire(s, f):
        raise _Timeout()
    old = signal.signal(signal.SIGALRM, fire)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)


_cache: dict[str, tuple] = {}


def load(name):
    if name not in _cache:
        _cache[name] = None
        for p in sorted(ND.iterdir()):
            if p.name.split(".")[0] == name:
                pn = parse_file(p, hashlib.sha256(p.read_bytes()).hexdigest())
                dag = MPDAG(nodes=set(pn.nodes), directed=set(pn.edges), undirected=set())
                _cache[name] = (dag, dag_to_cpdag(dag))
                break
    return _cache[name]


def main() -> None:
    rows = [json.loads(l) for l in (A3 / "instances.jsonl").open() if '"_meta"' not in l]
    rej = [r for r in rows if r.get("reject_reason") == REASON]
    adm = [r for r in rows if r.get("admissible")]

    inside, outside = [], 0
    for r in rej:
        g = load(r["network"])
        if g is None:
            continue
        (inside.append(r) if component_of(g[1], r["x"]) is not None else None)
        outside += component_of(g[1], r["x"]) is None
    adm_inside = sum(1 for r in adm if load(r["network"])
                     and component_of(load(r["network"])[1], r["x"]) is not None)
    print(f"rejected: {len(rej)}; X inside a component: {len(inside)}; outside: {outside}")
    print(f"admitted: {len(adm)}; X inside a component: {adm_inside}")

    out, stats = [], Counter()
    t0 = time.perf_counter()
    for i, r in enumerate(sorted(inside, key=lambda d: (d["network"], d["coverage"],
                                                        d["x"], d["y"])), 1):
        dag, cpdag = load(r["network"])
        try:
            with budget(PER_INSTANCE_S):
                row = g4.measure_one(dag, cpdag, r["network"], r["x"], r["y"],
                                     r["coverage"], stats)
        except _Timeout:
            stats["censored_instance_wall"] += 1
            print(f"  [{i}/{len(inside)}] {r['network']} {r['x']}->{r['y']} "
                  f"cov{r['coverage']}: CENSORED", flush=True)
            continue
        if row is None:
            continue
        stats["ok"] += 1
        out.append(row)
        print(f"  [{i}/{len(inside)}] {r['network']} {r['x']}->{r['y']} cov{r['coverage']}: "
              f"comp={row['component_size']} r_hop={row['radius']} "
              f"r_claim={row['claim_radius']} ({time.perf_counter() - t0:.0f}s)", flush=True)

    with OUT.open("w", encoding="utf-8") as f:
        for r in out:
            f.write(json.dumps(r) + "\n")
    meta = {
        "question": "[RE-3] on the rejected rows where the question is well posed",
        "n_rejected": len(rej), "n_inside_component": len(inside),
        "n_outside_component": outside, "n_admitted": len(adm),
        "n_admitted_inside_component": adm_inside, "n_measured": len(out),
        "per_instance_s": PER_INSTANCE_S, "outcomes": dict(stats),
        "radius": dict(sorted(Counter(r["radius"] for r in out).items())),
        "claim_radius": dict(sorted(Counter(r["claim_radius"] for r in out).items())),
        "by_network": dict(Counter(r["network"] for r in out)),
    }
    (OUT.with_suffix(".meta.json")).write_text(json.dumps(meta, indent=2))
    print("\n" + json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
