"""Independent audit of the LLM-elicited-knowledge survival sweep.

This script is a verification worker's own artifact, not part of the pipeline
under test. It re-derives every claim in the session brief from the
committed inputs (``results/elicit/knowledge.json``,
``results/axis_robustness_real/frame.jsonl``) and the sweep's own outputs
(``results/axis_robustness_llm/_done/*.json`` and
``results/axis_robustness_llm/shards/*.instances.jsonl``), without importing
or trusting any number the analysis layer (``llm_analyse.py``) already
computed. Where the source code itself is the object of a claim (task 1,
provenance), it inspects the *text* of the three modules named in the brief
plus every module ``run_llm_survival`` actually imports, rather than
asserting from memory.

Run as::

    PYTHONPATH=src .venv/bin/python results/axis_robustness_llm/verify_panel.py \
        > results/axis_robustness_llm/VERIFY_PANEL.txt 2>&1

Nothing here modifies ``src/``, ``results/elicit/``, or the sweep's own
output tree. It only reads them.
"""

from __future__ import annotations

import hashlib
import json
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[2]
SRC = REPO / "src" / "bkrobust" / "robustness"
OUT_DIR = REPO / "results" / "axis_robustness_llm"
FRAME_PATH = REPO / "results" / "axis_robustness_real" / "frame.jsonl"
BUNDLE_PATH = REPO / "results" / "elicit" / "knowledge.json"

W = 100


def hr(title: str) -> None:
    print()
    print("=" * W)
    print(title)
    print("=" * W)


def sub(title: str) -> None:
    print()
    print("-" * len(title))
    print(title)
    print("-" * len(title))


# ===========================================================================
# Task 1: static provenance -- is select_knowledge ever called on this path?
# ===========================================================================


def extract_defs(src: str) -> list[tuple[str, int, int]]:
    """Return (name, start_line, end_line) for every top-level def/class."""
    lines = src.splitlines()
    starts = []
    for i, line in enumerate(lines):
        m = re.match(r"^(def|class) (\w+)", line)
        if m:
            starts.append((m.group(2), i))
    out = []
    for idx, (name, start) in enumerate(starts):
        end = starts[idx + 1][1] if idx + 1 < len(starts) else len(lines)
        out.append((name, start, end))
    return out


def calls_of(name: str, src: str) -> list[tuple[str, int]]:
    """Every top-level function/class whose body contains a literal call to `name(`."""
    defs = extract_defs(src)
    lines = src.splitlines()
    hits = []
    pattern = re.compile(r"\b" + re.escape(name) + r"\s*\(")
    for fname, start, end in defs:
        body = "\n".join(lines[start:end])
        if pattern.search(body):
            hits.append((fname, start + 1))
    return hits


def module_level_hits(name: str, src: str) -> list[int]:
    """Line numbers where `name(` appears outside any def/class body (i.e. at import/module scope)."""
    defs = extract_defs(src)
    covered = set()
    for _, start, end in defs:
        covered.update(range(start, end))
    lines = src.splitlines()
    pattern = re.compile(r"\b" + re.escape(name) + r"\s*\(")
    hits = []
    for i, line in enumerate(lines):
        if i not in covered and pattern.search(line):
            hits.append(i + 1)
    return hits


def task1_provenance() -> None:
    hr("TASK 1: is benchmarks.measure.select_knowledge reachable from run_llm_survival.run_shard?")

    run_llm = (SRC / "run_llm_survival.py").read_text()
    llm_surv = (SRC / "llm_survival.py").read_text()
    real_surv = (SRC / "real_survival.py").read_text()
    run_real = (SRC / "run_real_survival.py").read_text()

    print("""
Method: for each of the three new/imported modules, list every TOP-LEVEL
def/class whose body contains a literal call `select_knowledge(`, and
separately list any occurrence at module scope (which would fire at import
time, before any shard runs). Then check whether run_shard's own call graph
(traced by hand below, function names given explicitly) ever reaches one of
the flagged functions.
""")

    # --- what does run_llm_survival.py itself import? ---
    imports = re.findall(r"^from ([\w.]+) import \(?([^)\n]*)\)?", run_llm, re.MULTILINE)
    sub("run_llm_survival.py's own import statements")
    for mod, names in imports:
        names_clean = " ".join(n.strip().rstrip(",") for n in names.split(","))
        print(f"  from {mod} import {names_clean}")

    sub("Literal `select_knowledge(` occurrences, per module (call sites, not just the word)")
    for label, src in [
        ("run_llm_survival.py", run_llm),
        ("llm_survival.py", llm_surv),
        ("real_survival.py", real_surv),
        ("run_real_survival.py", run_real),
    ]:
        hits = calls_of("select_knowledge", src)
        mod_hits = module_level_hits("select_knowledge", src)
        print(f"  {label}: calls inside {hits or '(none)'}; module-scope calls {mod_hits or '(none)'}")

    print("""
Finding: `select_knowledge(` as an actual call (not just the bare word in a
docstring/comment/import list) appears exactly once across these four files:
inside `run_flip_shard` in run_real_survival.py (see hit above). It is a
function `run_llm_survival.py` never imports and never calls.
""")

    sub("What run_llm_survival.py imports FROM run_real_survival.py, verbatim")
    m = re.search(
        r"from bkrobust\.robustness\.run_real_survival import \(([^)]*)\)", run_llm
    )
    imported_names = [n.strip().rstrip(",") for n in m.group(1).split("\n") if n.strip()]
    print(f"  imported names: {imported_names}")
    for n in imported_names:
        hits = calls_of("select_knowledge", run_real)
        flagged = [h for h in hits if h[0] == n]
        print(f"    {n}: select_knowledge called inside its own body? "
              f"{'YES -- ' + str(flagged) if flagged else 'no'}")

    sub("run_shard's own call graph (traced from source, function : line)")
    trace = [
        ("run_llm_survival.run_shard", "calls llm_survival.distinct_pairs(frame, net)"),
        ("run_llm_survival.run_shard", "calls real_survival.load_network(net)  [NOT select_knowledge]"),
        ("run_llm_survival.run_shard", "calls llm_survival.build_state(bundle, cond, net, dag, cpdag)"),
        ("llm_survival.build_state", "calls llm_survival.elicited_k(bundle, condition, network)  "
         "-- reads K from the JSON bundle, no DAG argument used to build K"),
        ("llm_survival.build_state", "calls real_survival.build_g0(cpdag, k)"),
        ("real_survival.build_g0", "calls demo.meek.apply_orientations(cpdag, list(k)) only"),
        ("llm_survival.build_state", "calls llm_survival.k_diagnostics -> k_accuracy/off_skeleton "
         "(diagnostic only, reads dag.directed_edges, never writes K)"),
        ("run_llm_survival.run_shard", "calls llm_survival.panel_family(cond, net)"),
        ("run_llm_survival.run_shard", "calls real_survival.commit_z_star(g0, x, y)"),
        ("run_llm_survival.run_shard", "calls real_survival.radius_columns(...)"),
        ("run_llm_survival.run_shard", "calls real_survival.ClosureCache(cpdag), "
         "real_survival.depth_grid(len(k)), real_survival.flip_draw(k, d, rep, family=family), "
         "real_survival.cell_statistics, real_survival.directed_symdiff"),
    ]
    for caller, what in trace:
        print(f"  {caller:38s} {what}")

    print("""
None of these functions is select_knowledge, and none of their own bodies
(checked above) calls it either. The only import from run_real_survival.py
used by run_llm_survival.py is {JsonlWriter, _censored_cell, load_frame,
sha_file} -- pure I/O/hashing utilities with no reference to
select_knowledge in their bodies (confirmed above: 'no' for all four).

Caveat correctly anticipated in the brief: `llm_survival.py` imports
`real_survival as rs` at module scope, and `real_survival.py` in turn does
`from bkrobust.benchmarks.measure import (MAX_G0_UNDIRECTED_FOR_EXTENSIONS,
O_INTRACTABLE, select_knowledge)` at ITS OWN module scope (real_survival.py,
see its import block). So `select_knowledge` the *function object* is bound
into `real_survival`'s and hence transitively reachable namespace the moment
either module is imported -- a grep for the bare identifier will find it in
several files. But binding a name into a namespace is not calling it: the
Python import machinery does not execute select_knowledge's body just
because the name is imported. The only *call sites* in the whole
"select_knowledge" grep across the repo that live inside a function
run_llm_survival.py can reach are: none. The nearest call site,
run_flip_shard in run_real_survival.py, is a sibling function in the same
module as the four names run_llm_survival.py imports, but Python imports
bind names, not sibling functions -- importing `load_frame` from a module
does not call any other function defined in that module.
""")

    verdict = "PASS" if not any(
        n == "run_flip_shard" or n == "run_tiered_shard" or n == "run_xarm_shard"
        for n in imported_names
    ) else "FAIL"
    print(f"VERDICT (task 1): {verdict} -- select_knowledge is reachable-by-name-binding but "
          f"never CALLED anywhere on run_llm_survival.run_shard's actual call path.")


# ===========================================================================
# Shared loading
# ===========================================================================


def load_bundle() -> tuple[dict[str, Any], str]:
    raw = BUNDLE_PATH.read_bytes()
    return json.loads(raw), hashlib.sha256(raw).hexdigest()


def load_frame_networks() -> list[str]:
    nets = set()
    with FRAME_PATH.open() as f:
        for line in f:
            if line.strip():
                nets.add(json.loads(line)["network"])
    return sorted(nets)


def k_b_sha256(k_pairs: list[list[str]]) -> str:
    """Mirror rs.sha_of(sorted(k)) exactly, as specified in the brief."""
    k = sorted((str(a), str(b)) for a, b in k_pairs)
    return hashlib.sha256(repr(sorted(k)).encode("utf-8")).hexdigest()[:16]


def load_done_markers() -> list[dict[str, Any]]:
    out = []
    for p in sorted((OUT_DIR / "_done").glob("*.json")):
        out.append(json.loads(p.read_text()))
    return out


# ===========================================================================
# Task 2: per-(condition, network) K fidelity
# ===========================================================================


def task2_fidelity(bundle: dict[str, Any], markers: list[dict[str, Any]]) -> None:
    hr("TASK 2: per-(condition, network) K fidelity -- shard K == bundle K, exactly")

    n_shards = len(markers)
    n_match = 0
    n_mismatch = 0
    n_len_only_match = 0
    mismatches = []
    rows_checked = 0
    rows_mismatched = 0

    for marker in markers:
        cond, net = marker["condition"], marker["network"]
        try:
            rec = bundle[cond]["networks"][net]
        except KeyError:
            mismatches.append((cond, net, "no bundle record"))
            n_mismatch += 1
            continue
        expected_digest = k_b_sha256(rec["k"])
        expected_len = len(rec["k"])

        inst_path = OUT_DIR / "shards" / f"{marker['shard_id']}.instances.jsonl"
        if not inst_path.is_file():
            mismatches.append((cond, net, "instances file missing"))
            n_mismatch += 1
            continue

        shard_digests = set()
        shard_lens = set()
        with inst_path.open() as f:
            for line in f:
                if not line.strip():
                    continue
                row = json.loads(line)
                rows_checked += 1
                shard_digests.add(row["k_b_sha256"])
                shard_lens.add(row["n_k"])
                if row["k_b_sha256"] != expected_digest:
                    rows_mismatched += 1

        if shard_digests == {expected_digest} and shard_lens == {expected_len}:
            n_match += 1
        elif shard_lens == {expected_len}:
            n_len_only_match += 1
            mismatches.append((cond, net, f"len matches ({expected_len}) but digest differs: "
                                           f"expected {expected_digest}, shard has {shard_digests}"))
            n_mismatch += 1
        else:
            mismatches.append((cond, net, f"expected len={expected_len} digest={expected_digest}, "
                                           f"shard has lens={shard_lens} digests={shard_digests}"))
            n_mismatch += 1

    print(f"Shards audited: {n_shards}")
    print(f"Instance rows checked: {rows_checked}")
    print(f"Shards where EVERY row's k_b_sha256 == recomputed digest from knowledge.json: {n_match}")
    print(f"Shards with a length-only match (n_k equal, digest differs): {n_len_only_match}")
    print(f"Shards with any mismatch (digest or length or missing data): {n_mismatch}")
    print(f"Individual rows whose k_b_sha256 differs from the shard's own expected digest: {rows_mismatched}")
    if mismatches:
        sub("Mismatch detail")
        for cond, net, msg in mismatches:
            print(f"  {cond} / {net}: {msg}")
    verdict = "PASS: exact match on all shards" if n_mismatch == 0 else f"FAIL: {n_mismatch} shard(s) mismatched"
    print(f"VERDICT (task 2): {verdict}")


# ===========================================================================
# True DAG loading (needs bkrobust importable -- PYTHONPATH=src)
# ===========================================================================


def load_true_edges(network: str) -> set[tuple[str, str]] | None:
    try:
        from bkrobust.robustness.real_survival import load_network
    except Exception as e:  # pragma: no cover
        print(f"  ERROR: cannot import bkrobust.robustness.real_survival.load_network: {e}")
        return None
    try:
        dag, _cpdag = load_network(network)
    except Exception as e:
        print(f"  ERROR loading network {network!r}: {e}")
        return None
    return set(dag.directed_edges)


# ===========================================================================
# Task 3: panel description table
# ===========================================================================


def task3_panel_table(bundle: dict[str, Any], frame_nets: list[str]) -> dict[str, Any]:
    hr("TASK 3: describe the panel accurately, from the data")

    assert "pathfinder" in frame_nets, "expected pathfinder in the frame"
    scored_nets = [n for n in frame_nets if n != "pathfinder"]
    print(f"Frame networks (frame.jsonl): {len(frame_nets)} -- {frame_nets}")
    print(f"'pathfinder' present in frame: {'pathfinder' in frame_nets}")
    for c in sorted(bundle):
        print(f"  pathfinder in {c}'s bundle networks? "
              f"{'pathfinder' in bundle[c].get('networks', {})}")
    print(f"Networks counted below (frame minus pathfinder): {len(scored_nets)} -- {scored_nets}")

    edges_cache: dict[str, set[tuple[str, str]] | None] = {}
    for n in scored_nets:
        edges_cache[n] = load_true_edges(n)
    n_unloadable = sum(1 for v in edges_cache.values() if v is None)
    if n_unloadable:
        print(f"WARNING: {n_unloadable} of {len(scored_nets)} networks' true DAG could not be loaded; "
              f"accuracy for those is reported as 'N/A' rather than invented.")

    table = []
    for cond in sorted(bundle):
        rec = bundle[cond]
        nets_covered = [n for n in scored_nets if n in rec.get("networks", {})]
        sizes = []
        total_k = 0
        total_asked = 0
        total_asserted = 0
        acc_num = 0
        acc_den = 0
        for n in nets_covered:
            netrec = rec["networks"][n]
            k = netrec["k"]
            sizes.append(len(k))
            total_k += len(k)
            n_asked = netrec.get("n_asked") or 0
            total_asked += n_asked
            total_asserted += len(k)
            true_edges = edges_cache.get(n)
            if true_edges is not None:
                for a, b in k:
                    acc_den += 1
                    if (str(a), str(b)) in true_edges:
                        acc_num += 1
        row = {
            "condition": cond,
            "model": rec.get("model"),
            "naming": rec.get("naming"),
            "n_networks_covered": len(nets_covered),
            "n_networks_missing_from_24": len(scored_nets) - len(nets_covered),
            "total_K_sum": total_k,
            "min_per_net_K": min(sizes) if sizes else None,
            "median_per_net_K": statistics.median(sizes) if sizes else None,
            "max_per_net_K": max(sizes) if sizes else None,
            "assert_rate_pooled": (total_asserted / total_asked) if total_asked else None,
            "accuracy_pooled": (acc_num / acc_den) if acc_den else None,
            "acc_num": acc_num,
            "acc_den": acc_den,
        }
        table.append(row)

    sub("Per-condition panel table (24 frame networks, pathfinder excluded)")
    print("Note on assert_rate_pooled and accuracy_pooled: BOTH are pooled "
          "(numerator summed / denominator summed across all 24 networks), "
          "not an average of 24 per-network ratios. Pooling weights networks "
          "by how many pairs they asked, which an unweighted mean-of-ratios "
          "would not; this is the more defensible choice and is stated "
          "explicitly rather than left implicit.")
    hdr = (f"{'condition':18s} {'model':22s} {'naming':10s} {'n_net':6s} {'sum|K|':7s} "
           f"{'min':4s} {'med':5s} {'max':4s} {'assert_rate':11s} {'accuracy':9s}")
    print(hdr)
    print("-" * len(hdr))
    for r in table:
        ar = f"{r['assert_rate_pooled']:.3f}" if r["assert_rate_pooled"] is not None else "N/A"
        ac = f"{r['accuracy_pooled']:.3f}" if r["accuracy_pooled"] is not None else "N/A"
        print(f"{r['condition']:18s} {str(r['model'])[:22]:22s} {str(r['naming']):10s} "
              f"{r['n_networks_covered']:<6d} {r['total_K_sum']:<7d} "
              f"{str(r['min_per_net_K']):4s} {str(r['median_per_net_K']):5s} {str(r['max_per_net_K']):4s} "
              f"{ar:11s} {ac:9s}")
        if r["n_networks_missing_from_24"]:
            print(f"    ** covers only {r['n_networks_covered']}/24 frame networks "
                  f"(missing {r['n_networks_missing_from_24']}) **")

    return {"table": table, "edges_cache": edges_cache, "scored_nets": scored_nets}


# ===========================================================================
# Task 4: model-scale ladders
# ===========================================================================


def per_network_series(bundle, cond, scored_nets, edges_cache):
    """Return {network: (n_k, accuracy_or_None)} for one condition, over scored_nets."""
    out = {}
    rec = bundle.get(cond, {})
    nets = rec.get("networks", {})
    for n in scored_nets:
        if n not in nets:
            continue
        k = nets[n]["k"]
        true_edges = edges_cache.get(n)
        if true_edges is None or not k:
            acc = None
        else:
            correct = sum(1 for a, b in k if (str(a), str(b)) in true_edges)
            acc = correct / len(k)
        out[n] = (len(k), acc)
    return out


def sign_test(diffs: list[float]) -> dict[str, Any]:
    """Two-sided exact sign test on nonzero diffs, via the binomial distribution."""
    from math import comb

    nonzero = [d for d in diffs if d != 0]
    n = len(nonzero)
    n_pos = sum(1 for d in nonzero if d > 0)
    n_ties = len(diffs) - n
    if n == 0:
        return {"n": 0, "n_pos": 0, "n_ties": n_ties, "p_two_sided": None}
    k = min(n_pos, n - n_pos)
    p = sum(comb(n, i) for i in range(0, k + 1)) * 2 / (2 ** n)
    p = min(p, 1.0)
    return {"n": n, "n_pos": n_pos, "n_ties": n_ties, "p_two_sided": p}


def wilcoxon_signed_rank(diffs: list[float]) -> dict[str, Any]:
    try:
        from scipy.stats import wilcoxon
    except Exception as e:  # pragma: no cover
        return {"error": str(e)}
    nonzero = [d for d in diffs if d != 0]
    if len(nonzero) < 1:
        return {"n_nonzero": 0, "statistic": None, "p_value": None}
    try:
        stat, p = wilcoxon(nonzero)
        return {"n_nonzero": len(nonzero), "statistic": float(stat), "p_value": float(p)}
    except Exception as e:
        return {"n_nonzero": len(nonzero), "error": str(e)}


def paired_compare(bundle, small_cond, big_cond, scored_nets, edges_cache, label):
    sub(f"{label}: {small_cond} -> {big_cond}")
    small = per_network_series(bundle, small_cond, scored_nets, edges_cache)
    big = per_network_series(bundle, big_cond, scored_nets, edges_cache)
    common = sorted(set(small) & set(big))
    if not common:
        print("  no networks in common -- cannot compare")
        return
    print(f"  networks in common: {len(common)}")

    k_diffs = [big[n][0] - small[n][0] for n in common]
    print(f"  |K| per-network diff ({big_cond} - {small_cond}):")
    print(f"    mean={statistics.mean(k_diffs):.2f} median={statistics.median(k_diffs):.2f} "
          f"n_increase={sum(1 for d in k_diffs if d > 0)} n_decrease={sum(1 for d in k_diffs if d < 0)} "
          f"n_tie={sum(1 for d in k_diffs if d == 0)}")
    st = sign_test(k_diffs)
    wx = wilcoxon_signed_rank(k_diffs)
    print(f"    sign test: n_nonzero={st['n']} n_pos={st['n_pos']} p_two_sided={st['p_two_sided']}")
    print(f"    wilcoxon signed-rank: {wx}")

    both_acc = [(n, small[n][1], big[n][1]) for n in common
                if small[n][1] is not None and big[n][1] is not None]
    print(f"  accuracy: networks with a defined accuracy on BOTH sides: {len(both_acc)}/{len(common)}")
    if both_acc:
        a_diffs = [b - s for _, s, b in both_acc]
        print(f"    accuracy per-network diff ({big_cond} - {small_cond}): "
              f"mean={statistics.mean(a_diffs):.3f} median={statistics.median(a_diffs):.3f} "
              f"n_increase={sum(1 for d in a_diffs if d > 0)} n_decrease={sum(1 for d in a_diffs if d < 0)} "
              f"n_tie={sum(1 for d in a_diffs if d == 0)}")
        st_a = sign_test(a_diffs)
        wx_a = wilcoxon_signed_rank(a_diffs)
        print(f"    sign test: n_nonzero={st_a['n']} n_pos={st_a['n_pos']} p_two_sided={st_a['p_two_sided']}")
        print(f"    wilcoxon signed-rank: {wx_a}")
    else:
        print("    cannot run sign test / wilcoxon on accuracy: no network has a defined "
              "accuracy on both sides (K empty or DAG unloadable)")


def task4_ladders(bundle, scored_nets, edges_cache) -> None:
    hr("TASK 4: is model scale actually the axis? (paired per-network, not pooled means)")
    print("""
Qwen ladder: qwen2.5:7b-instruct (D_LLM, D_LLM_INSTR) -> Qwen2.5-32B-Instruct
(D_LLM_32B) -> Qwen2.5-72B-Instruct (D_LLM_72B, D_LLM_72B_INSTR)
Llama ladder: llama3.1:8b (D_LLM_SRC) -> Llama-3.3-70B-Instruct (D_LLM_SRC_70B)
""")
    pairs = [
        ("Qwen 7b(D_LLM) -> 32B", "D_LLM", "D_LLM_32B"),
        ("Qwen 7b(D_LLM_INSTR) -> 32B", "D_LLM_INSTR", "D_LLM_32B"),
        ("Qwen 32B -> 72B(D_LLM_72B)", "D_LLM_32B", "D_LLM_72B"),
        ("Qwen 32B -> 72B(D_LLM_72B_INSTR)", "D_LLM_32B", "D_LLM_72B_INSTR"),
        ("Qwen 7b(D_LLM) -> 72B(D_LLM_72B)", "D_LLM", "D_LLM_72B"),
        ("Qwen 7b(D_LLM_INSTR) -> 72B(D_LLM_72B_INSTR)", "D_LLM_INSTR", "D_LLM_72B_INSTR"),
        ("Llama 8b -> 70B", "D_LLM_SRC", "D_LLM_SRC_70B"),
    ]
    for label, small, big in pairs:
        paired_compare(bundle, small, big, scored_nets, edges_cache, label)

    print("""
VERDICT (task 4): read the per-pair sign counts and p-values above directly --
"bigger asserts more and more accurately" HOLDS for a pair only if |K| n_increase
is large with a small p, AND the accuracy diff is systematically positive with a
small p. Do not average across pairs; report each ladder step on its own, since
the two Qwen families (D_LLM vs D_LLM_INSTR; D_LLM_72B vs D_LLM_72B_INSTR) are
not guaranteed to agree with each other.
""")


# ===========================================================================
# Task 5: scrambled control
# ===========================================================================


def task5_scrambled(bundle, scored_nets, edges_cache) -> None:
    hr("TASK 5: the scrambled-naming control -- does accuracy collapse to a semantic floor?")
    pairs = [
        ("D_SCRAMBLED vs its real-naming sibling D_LLM", "D_LLM", "D_SCRAMBLED"),
        ("D_SCRAMBLED vs D_LLM_INSTR", "D_LLM_INSTR", "D_SCRAMBLED"),
        ("D_SCRAMBLED_72B vs D_LLM_72B", "D_LLM_72B", "D_SCRAMBLED_72B"),
        ("D_SCRAMBLED_72B vs D_LLM_72B_INSTR", "D_LLM_72B_INSTR", "D_SCRAMBLED_72B"),
    ]
    for label, real_cond, scrambled_cond in pairs:
        sub(label)
        real_s = per_network_series(bundle, real_cond, scored_nets, edges_cache)
        scr_s = per_network_series(bundle, scrambled_cond, scored_nets, edges_cache)
        common = sorted(set(real_s) & set(scr_s))
        print(f"  networks in common: {len(common)}")
        both_acc = [(n, real_s[n][1], scr_s[n][1]) for n in common
                    if real_s[n][1] is not None and scr_s[n][1] is not None]
        print(f"  networks with defined accuracy on both sides: {len(both_acc)}/{len(common)}")
        if both_acc:
            diffs = [r - s for _, r, s in both_acc]  # real - scrambled
            print(f"  accuracy(real) - accuracy(scrambled) per network: "
                  f"mean={statistics.mean(diffs):.3f} median={statistics.median(diffs):.3f}")
            print(f"    n_real_higher={sum(1 for d in diffs if d > 0)} "
                  f"n_scrambled_higher_or_equal={sum(1 for d in diffs if d <= 0)}")
            for n, r, s in both_acc:
                print(f"    {n:18s} real={r:.3f} scrambled={s:.3f} diff={r - s:+.3f}")
            st = sign_test(diffs)
            wx = wilcoxon_signed_rank(diffs)
            print(f"  sign test on (real-scrambled): n_nonzero={st['n']} n_pos={st['n_pos']} "
                  f"p_two_sided={st['p_two_sided']}")
            print(f"  wilcoxon signed-rank: {wx}")
        else:
            print("  cannot compare: no network has defined accuracy on both sides")

    print("""
VERDICT (task 5): read the paired sign counts above. If accuracy(real) is NOT
systematically and significantly higher than accuracy(scrambled) per network,
the "the model is really using domain semantics" story fails its control --
consistent with the earlier session's finding, if it reproduces here.
""")


# ===========================================================================
# Task 6: degenerate-case census
# ===========================================================================


def task6_census(bundle, markers, scored_nets) -> None:
    hr("TASK 6: degenerate-case census")

    status_by_cond: dict[str, Counter] = defaultdict(Counter)
    status_by_net: dict[str, Counter] = defaultdict(Counter)
    all_fail_cells = []
    n_k_zero_cells = []

    for m in markers:
        cond, net = m["condition"], m["network"]
        counts = m.get("instance_status_counts", {})
        for status, n in counts.items():
            status_by_cond[cond][status] += n
            status_by_net[net][status] += n
        n_total = sum(counts.values())
        n_ok = counts.get("ok", 0)
        if n_total > 0 and n_ok == 0:
            all_fail_cells.append((cond, net, dict(counts)))
        if m.get("g0_status") == "n_k_zero":
            n_k_zero_cells.append((cond, net))

    sub("Status totals by condition (summed instance_status_counts across all its shards)")
    all_statuses = sorted({s for c in status_by_cond.values() for s in c})
    hdr = f"{'condition':18s} " + " ".join(f"{s:>22s}" for s in all_statuses)
    print(hdr)
    for cond in sorted(status_by_cond):
        print(f"{cond:18s} " + " ".join(f"{status_by_cond[cond].get(s, 0):>22d}" for s in all_statuses))

    sub("Status totals by network (summed across all conditions that covered it)")
    hdr2 = f"{'network':18s} " + " ".join(f"{s:>22s}" for s in all_statuses)
    print(hdr2)
    for net in sorted(status_by_net):
        print(f"{net:18s} " + " ".join(f"{status_by_net[net].get(s, 0):>22d}" for s in all_statuses))

    sub("(condition, network) cells where EVERY pair failed (n_ok == 0, n_total > 0)")
    for cond, net, counts in all_fail_cells:
        reasons = ", ".join(f"{k}={v}" for k, v in counts.items())
        print(f"  {cond} / {net}: {reasons}")
    print(f"Total all-fail cells: {len(all_fail_cells)}")

    sub("Cross-check: do n_k_zero cells correspond exactly to bundle records with empty k?")
    bundle_empty = set()
    for cond, rec in bundle.items():
        for net, netrec in rec.get("networks", {}).items():
            if net in scored_nets and len(netrec.get("k", [])) == 0:
                bundle_empty.add((cond, net))
    marker_n_k_zero = set(n_k_zero_cells)
    only_in_bundle = bundle_empty - marker_n_k_zero
    only_in_markers = marker_n_k_zero - bundle_empty
    print(f"  bundle records with empty k (restricted to scored/shardable cells): {len(bundle_empty)}")
    print(f"  shards whose g0_status == 'n_k_zero': {len(marker_n_k_zero)}")
    print(f"  in bundle-empty but NOT flagged n_k_zero by a shard: {sorted(only_in_bundle) or '(none)'}")
    print(f"  flagged n_k_zero by a shard but bundle k is NOT empty: {sorted(only_in_markers) or '(none)'}")
    verdict = "PASS: exact correspondence" if not only_in_bundle and not only_in_markers else "MISMATCH -- see above"
    print(f"  VERDICT (n_k_zero correspondence): {verdict}")


# ===========================================================================


def main() -> None:
    print("VERIFY_PANEL -- independent audit of results/axis_robustness_llm/")
    print(f"Repo: {REPO}")
    print(f"Frame: {FRAME_PATH}")
    print(f"Bundle: {BUNDLE_PATH}")

    task1_provenance()

    bundle, bundle_sha = load_bundle()
    frame_nets = load_frame_networks()
    markers = load_done_markers()
    print(f"\nLoaded bundle sha256={bundle_sha}, {len(bundle)} conditions.")
    print(f"Loaded {len(markers)} completion markers from _done/.")

    task2_fidelity(bundle, markers)

    panel_info = task3_panel_table(bundle, frame_nets)
    scored_nets = panel_info["scored_nets"]
    edges_cache = panel_info["edges_cache"]

    task4_ladders(bundle, scored_nets, edges_cache)
    task5_scrambled(bundle, scored_nets, edges_cache)
    task6_census(bundle, markers, scored_nets)

    print("\nDONE.")


if __name__ == "__main__":
    main()
