"""State builders for the E11 frames and the synthetic survival instances, vendored for E25.

Sources, copied on 2026-11-07 from this repository's working tree (untracked rebuttal folders,
HEAD 5421cabc23b28af8f7c6e4124a698d66a76e0a08):

* experiments/rebuttal/E11_survival_on_cpdag_only_queries/run.py: FRAMES, frame_units,
  load_network, get_bundle, and the state lines of score_item (ls.build_state, rs.commit_z_star).
  Code unchanged except that the survival draws are dropped (E25 needs the state and Z* only).
* experiments/rebuttal/E15_aid_as_survival_baseline/run.py: UNTARGETED_CSV, TARGETED_CSV, rebuild.
  Unchanged.

Every bkrobust name imported here exists in the submitted tree.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from bkrobust.benchmarks import describe as _describe
from bkrobust.benchmarks.acquire import load_cached
from bkrobust.demo.example import dag_to_cpdag
from bkrobust.robustness import llm_survival as ls
from bkrobust.robustness import real_survival as rs
from bkrobust.robustness import survival as sv
from bkrobust.robustness import unconfounded_generator as ug

# ---------------------------------------------------------------- E11 (frames and states)

FRAMES = {
    "cpdag_only": ("results/frame/frame.jsonl",
                   "6f05f044b13117d731aa345d4d549eebfa0586d7da11ceff953502003c7843fc"),
    "table1": ("results/axis_robustness_real/frame.jsonl",
               "73881798319206be0a657c849431ad35c4cbddc5c4c3eecc4790d18c84981f38"),
}
BUNDLE_SHA = "45119393225fee57ed1dd4c002130f805d6a655c2ee7648006470918c2dfc7bb"


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def frame_units(name: str, rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """One record per distinct (network, x, y), as E11's frame_units (columns trimmed)."""
    out: list[dict[str, Any]] = []
    if name == "table1":
        for net in sorted({r["network"] for r in rows}):
            for r in ls.distinct_pairs(rows, net):
                out.append({"network": net, "x": r["x"], "y": r["y"], "frame_row_id": r["row_id"]})
        return out
    seen: set[tuple[str, str, str]] = set()
    for r in rows:
        key = (r["network"], r["x"], r["y"])
        if key in seen:
            raise SystemExit(f"cpdag_only frame repeats {key}")
        seen.add(key)
        out.append({"network": r["network"], "x": r["x"], "y": r["y"],
                    "frame_row_id": hashlib.sha256(f"cpdag_only|{r['network']}|{r['x']}|{r['y']}"
                                                   .encode()).hexdigest()[:16]})
    return out


_NET: dict[str, Any] = {}
_BUNDLE: dict[str, Any] = {}
_STATE: dict[tuple[str, str], Any] = {}


def load_network(name: str, networks_dir: str):
    if name not in _NET:
        acq = load_cached(Path(networks_dir))
        recs = {r.name: r for r in acq.files}
        fname = rs.network_filename(name, recs)
        parsed = _describe.parse_file(recs[fname].path, recs[fname].sha256)
        dag = _describe.to_mpdag(parsed)
        if len(_NET) > 4:
            _NET.clear()
        _NET[name] = (dag, dag_to_cpdag(dag))
    return _NET[name]


def get_bundle(path: str) -> dict[str, Any]:
    if path not in _BUNDLE:
        _BUNDLE[path] = ls.load_bundle(Path(path))[0]
    return _BUNDLE[path]


def e11_state(elicit: str, networks_dir: str, cond: str, net: str):
    """(dag, cpdag, g0 or None, status) for one (condition, network), as E11's score_item."""
    key = (cond, net)
    if key not in _STATE:
        bundle = get_bundle(elicit)
        if not ls.has_network(bundle, cond, net):
            _STATE[key] = (None, None, None, "network_absent_from_bundle")
        else:
            dag, cpdag = load_network(net, networks_dir)
            _k, g0, reason, _shared = ls.build_state(bundle, cond, net, dag, cpdag)
            if len(_STATE) > 8:
                _STATE.clear()
            _STATE[key] = (dag, cpdag, g0, "ok" if g0 is not None else reason)
    return _STATE[key]


def e11_commit(g0, x: str, y: str):
    return rs.commit_z_star(g0, x, y)


def e11_conditions(elicit: str) -> list[str]:
    return ls.panel_conditions(get_bundle(elicit), naming=ls.REAL_NAMING)


# ---------------------------------------------------------------- E15 (synthetic instances)

UNTARGETED_CSV = "results/axis_robustness_unconfounded_v2/survival_instances.csv"
TARGETED_CSV = "results/axis_robustness_p6/survival_instances.csv"


def rebuild(row: dict):
    arm = row["arm"]
    c, seed = int(row["component_size"]), int(row["seed"])
    if arm == "tiered_unconfounded":
        return ug.build_unconfounded_tiered_instance(c, seed, int(row["n_tiers"]))
    if arm == "flip":
        return sv.build_flip_instance(c, int(row["separation"]), float(row["coverage"]),
                                      float(row["base_wrongness"]), seed)
    if arm == "tiered":
        return sv.build_tiered_instance(c, int(row["separation"]), seed, int(row["n_tiers"]))
    return None, f"unknown_arm:{arm}"
