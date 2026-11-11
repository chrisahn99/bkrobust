"""E11 driver: the submitted survival shard runner on the CPDAG-only frame.

``python -m bkrobust.robustness.run_llm_survival pool`` launches each shard as
``run --shard`` in a subprocess without forwarding ``--frame-dir``, so every
subprocess silently reads the default frame (the paper's Table 1 frame)
instead of the one given to ``pool``; shards of networks absent from that
frame fail as ``unknown shard``. This driver calls the same
``run_llm_survival.run_one`` in-process, with the CPDAG-only frame, one shard
per pool task. Nothing in src/ is changed. 1,000 draws per grid point, as in
the paper.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing as mp
from pathlib import Path

from common import ROOT

from bkrobust.robustness import llm_survival as ls
from bkrobust.robustness import real_survival as rs
from bkrobust.robustness import run_llm_survival as rl

FRAME_DIR = ROOT / "results/rebuttal/e11_frame"
OUT_DIR = ROOT / "results/rebuttal/e11"


def _one(spec):
    frame = rl.load_frame(FRAME_DIR)
    frame_sha = json.loads((FRAME_DIR / "frame_hash.json").read_text())["frame_sha256"]
    bundle, bundle_sha = ls.load_bundle()
    sid = rl.shard_id(spec)
    if rl.is_complete(OUT_DIR, sid):
        return sid, "skip"
    try:
        p = rl.run_one(spec, frame, bundle, bundle_sha, frame_sha, OUT_DIR, rs.N_DRAWS_DEFAULT)
        return sid, f"ok scored={p['n_scored']} {p['elapsed_s']}s"
    except Exception as e:  # noqa: BLE001
        return sid, f"FAIL {type(e).__name__}: {e}"


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    frame = rl.load_frame(FRAME_DIR)
    bundle, _ = ls.load_bundle()
    specs = rl.enumerate_shards(frame, bundle, naming="real")
    specs.sort(key=lambda s: -rl.estimated_cost_s(s, frame, bundle, rs.N_DRAWS_DEFAULT))
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(len(specs), "shards", flush=True)
    with mp.Pool(a.workers, maxtasksperchild=1) as P:
        for i, (sid, msg) in enumerate(P.imap_unordered(_one, specs)):
            print(i, sid, msg, flush=True)
