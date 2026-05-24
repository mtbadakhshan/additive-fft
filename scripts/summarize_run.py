#!/usr/bin/env python3
"""Top-level index for a combined measurement run (bench + perf subdirs)."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from summary_common import emit_header, emit_run_context, parse_run_meta


def main() -> int:
    if len(sys.argv) != 2:
        print(f"usage: {sys.argv[0]} <measure_results_dir>", file=sys.stderr)
        return 2
    root = Path(sys.argv[1])
    if not root.is_dir():
        print(f"error: not a directory: {root}", file=sys.stderr)
        return 1

    meta = parse_run_meta(root / "RUN_META.txt")
    emit_header(
        "measurement run",
        root,
        "Index for wall-clock benchmarks (`bench/`) and hardware counters (`perf/`).",
    )
    emit_run_context(
        meta,
        prefer_keys=[
            "date_utc",
            "repo",
            "build",
            "quick",
            "MIN_RANGE",
            "MAX_RANGE",
            "STEP",
            "THREAD_LIST",
            "bench_status",
            "perf_status",
            "config_file",
        ],
        only_preferred=True,
    )

    print("## Reports\n")
    print("| Tool | Status | Report |")
    print("|------|--------|--------|")
    for sub, label in (("bench", "Google Benchmark"), ("perf", "perf stat")):
        status = meta.get(f"{sub}_status", meta.get(f"{sub}_dir", "not run"))
        if (root / sub / "SUMMARY.md").is_file():
            link = f"[{sub}/SUMMARY.md]({sub}/SUMMARY.md)"
        else:
            link = "—"
        print(f"| {label} | {status} | {link} |")
    print()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
