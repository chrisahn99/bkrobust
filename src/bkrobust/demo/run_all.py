"""Entrypoint: run every scenario, write results, build figures, print a summary.

    PYTHONPATH=src python3 -m bkrobust.demo.run_all

Seeded throughout; see ``report.md`` for reproduction details and runtime.
"""

from __future__ import annotations

import time


def main() -> None:
    """Run the whole demonstration end to end."""
    from bkrobust.demo.export import write_all

    start = time.time()
    results = write_all()
    for label, r in results.items():
        print(f"[{label}] |space|={len(r['space'])}  radii={r['radii']}")
    print(f"results written to results/breakdown_radius_demo/ in {time.time() - start:.1f}s")

    try:
        from bkrobust.demo.figures import build_all_figures
    except ImportError:
        print("figures module not available; skipping figures")
        return
    paths = build_all_figures(results)
    print(f"{len(paths)} figures written to figures/")


if __name__ == "__main__":
    main()
