#!/usr/bin/env python3
"""Fidelity gate: did a run implement the experiment it is named after? Decided apart from the score.

    python experiments/rebuttal/fidelity_gate.py experiments/rebuttal/E1_cpdag_check_on_lower_bound_rows/FIDELITY.yaml
    cd experiments/rebuttal/E4_class_enumeration_baseline && python ../fidelity_gate.py FIDELITY.yaml

Each experiment folder holds a FIDELITY.yaml rubric, written before the run. Its keys:

    arm:   short name of the experiment
    says:  what the experiment is, in one paragraph
    repo:  base directory for `check` commands and `code` paths, relative to the rubric
    code:  the files that implement the experiment
    runs:  glob (relative to the rubric) of the run configs the run wrote
    rubric: list of items, each with an `id`, a `says` line and one or more of
      config: {key: <name>, equals | in | absent | present | absent_or: <value>}
              checked against every run config matched by `runs`
      check:  a read-only shell command run in `repo`; exit 0 = pass
      judge:  true marks a line that needs a person to read the code; this script reports it
              as "to read" and does not decide it

A bare key is searched anywhere in the nested config and must have a single value; a dotted key
(`a.b.0`) is a path from the root.

Writes FIDELITY.verdict.json beside the rubric (sha256 of every file read, one verdict per item)
and prints a table. Exit 0 = every decided item passes, 1 = at least one FAIL, 2 = no FAIL but
something unclear, 3 = rubric unreadable. Needs PyYAML.
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

_MISSING = object()


def expand(p, base: Path) -> Path:
    p = Path(os.path.expanduser(str(p)))
    return p if p.is_absolute() else base / p


def load_any(path: Path):
    text = path.read_text()
    if path.suffix.lower() == ".json":
        return json.loads(text)
    try:
        import yaml
    except ImportError:
        sys.exit("PyYAML is needed to read FIDELITY.yaml: pip install pyyaml")
    return yaml.safe_load(text)


def sha256_path(path: Path) -> str:
    h = hashlib.sha256()
    if path.is_dir():
        for f in sorted(q for q in path.rglob("*") if q.is_file()):
            h.update(str(f.relative_to(path)).encode() + b"\0" + sha256_path(f).encode() + b"\n")
        return "dir:" + h.hexdigest()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def find_key(obj, key: str):
    if "." in key:
        cur = obj
        for part in key.split("."):
            if isinstance(cur, dict) and part in cur:
                cur = cur[part]
            elif isinstance(cur, list) and part.lstrip("-").isdigit() and -len(cur) <= int(part) < len(cur):
                cur = cur[int(part)]
            else:
                return _MISSING
        return cur
    hits = []

    def walk(o):
        if isinstance(o, dict):
            for k, v in o.items():
                if k == key:
                    hits.append(v)
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)

    walk(obj)
    if not hits:
        return _MISSING
    first = json.dumps(hits[0], sort_keys=True, default=str)
    if any(json.dumps(h, sort_keys=True, default=str) != first for h in hits[1:]):
        raise ValueError(f"key {key!r} occurs with different values; use a dotted path")
    return hits[0]


def _eq(a, b) -> bool:
    if isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool):
        return abs(float(a) - float(b)) <= 1e-9 * max(1.0, abs(float(b)))
    return a == b


def eval_config(spec: dict, cfg) -> tuple[str, str]:
    key = spec.get("key")
    if not key:
        return "unclear", "config item has no key"
    try:
        v = find_key(cfg, key)
    except ValueError as e:
        return "unclear", str(e)
    missing = v is _MISSING
    if spec.get("absent"):
        return ("pass", f"{key} absent") if missing else ("fail", f"{key} = {v!r}, must be absent")
    if spec.get("present"):
        return ("fail", f"{key} missing") if missing else ("pass", f"{key} = {v!r}")
    if missing:
        if "absent_or" in spec:
            v = spec["absent_or"]
        else:
            return "fail", f"{key} missing from the run config"
    if "equals" in spec:
        return ("pass" if _eq(v, spec["equals"]) else "fail"), f"{key} = {v!r} (want {spec['equals']!r})"
    if "absent_or" in spec:
        return ("pass" if _eq(v, spec["absent_or"]) else "fail"), f"{key} = {v!r} (want {spec['absent_or']!r} or absent)"
    if "in" in spec:
        return ("pass" if any(_eq(v, w) for w in spec["in"]) else "fail"), f"{key} = {v!r} (want one of {spec['in']!r})"
    return "unclear", "config item has no operator (equals / in / absent / present / absent_or)"


def eval_check(cmd: str, repo: Path) -> tuple[str, str]:
    try:
        p = subprocess.run(cmd, shell=True, cwd=str(repo), capture_output=True, text=True, timeout=120)
    except subprocess.TimeoutExpired:
        return "unclear", "check timed out (120 s)"
    tail = (p.stdout + p.stderr).strip().splitlines()[-1:] or [""]
    return ("pass" if p.returncode == 0 else "fail"), f"exit {p.returncode} {tail[0][:160]}"


def run_rubric(path: Path, write: bool = True) -> dict:
    rub = load_any(path)
    if not isinstance(rub, dict) or not isinstance(rub.get("rubric"), list):
        raise ValueError("rubric file needs a top-level `rubric:` list")
    base = path.parent
    repo = expand(rub.get("repo", "."), base)
    pats = [rub["runs"]] if isinstance(rub.get("runs"), str) else (rub.get("runs") or [])
    run_files = [Path(p) for pat in pats for p in sorted(glob.glob(str(expand(pat, base))))]
    hashes = {str(f): sha256_path(f) for f in run_files}
    for rel in rub.get("code") or []:
        p = expand(rel, repo)
        if p.is_file():
            hashes[str(p)] = sha256_path(p)
    results = []
    for item in rub["rubric"]:
        rid = item.get("id", "?")
        if "config" in item:
            if not run_files:
                results.append({"id": rid, "kind": "config", "verdict": "unclear", "detail": "no run config matched `runs`"})
            for f in run_files:
                try:
                    v, d = eval_config(item["config"], load_any(f))
                except Exception as e:  # noqa: BLE001
                    v, d = "unclear", f"cannot read: {e}"
                results.append({"id": rid, "kind": "config", "run": str(f), "verdict": v, "detail": d})
        if "check" in item:
            v, d = eval_check(item["check"], repo)
            results.append({"id": rid, "kind": "check", "verdict": v, "detail": d})
        if item.get("judge"):
            results.append({"id": rid, "kind": "read", "verdict": "to read",
                            "detail": "decided by reading the code listed in `code`"})
        if not any(k in item for k in ("config", "check")) and not item.get("judge"):
            results.append({"id": rid, "kind": "none", "verdict": "unclear", "detail": "item has no config, check or judge"})
    verdicts = [r["verdict"] for r in results]
    overall = "FAIL" if "fail" in verdicts else ("UNCLEAR" if "unclear" in verdicts else "PASS")
    out = {"rubric": str(path), "arm": rub.get("arm"), "overall": overall,
           "checked_at": time.strftime("%Y-%m-%dT%H:%M:%S"), "files_sha256": hashes, "results": results}
    if write:
        path.with_name(path.stem + ".verdict.json").write_text(json.dumps(out, indent=2, ensure_ascii=False))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("rubric", type=Path)
    ap.add_argument("--no-judge", action="store_true", help="accepted for compatibility; code-reading items are never decided here")
    ap.add_argument("--no-write", action="store_true", help="do not write FIDELITY.verdict.json")
    a = ap.parse_args()
    try:
        out = run_rubric(a.rubric, write=not a.no_write)
    except Exception as e:  # noqa: BLE001
        print(f"rubric unreadable: {e}", file=sys.stderr)
        return 3
    print(f"{out['overall']:7s} {out['arm']}  ({out['rubric']})")
    for r in out["results"]:
        run = f" [{Path(r['run']).name}]" if r.get("run") else ""
        print(f"  {r['verdict']:8s} {r['id']:<18s} {r['kind']:<6s}{run} {r['detail']}")
    return {"PASS": 0, "FAIL": 1, "UNCLEAR": 2}[out["overall"]]


if __name__ == "__main__":
    sys.exit(main())
