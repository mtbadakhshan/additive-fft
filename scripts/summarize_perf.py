#!/usr/bin/env python3
"""
Build markdown tables from perf_stat.sh artifacts (*.perf.txt).
"""

from __future__ import annotations

import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from summary_common import (
    NON_MEASURED,
    emit_header,
    emit_run_context,
    expected_thread_counts,
    parse_run_meta,
    perf_expected_m_values,
    perf_expected_variants,
)


def parse_perf_txt(path: Path) -> dict:
    text = path.read_text()
    out: dict[str, float | str] = {"file": path.name}

    m = re.search(
        r"\[perf_driver\] variant=(\S+) m=(\d+) iters=(\d+) total=([\d.]+)s per_call=([\d.]+)ms",
        text,
    )
    if m:
        out["variant"] = m.group(1)
        out["m"] = int(m.group(2))
        out["iters"] = int(m.group(3))
        out["total_s"] = float(m.group(4))
        out["per_call_ms"] = float(m.group(5))

    pref = "cpu_core/"
    for line in text.splitlines():
        line = line.strip()
        if "insn per cycle" in line and pref in line:
            ipc_m = re.search(r"#\s*([\d.]+)\s+insn per cycle", line)
            if ipc_m:
                out["ipc"] = float(ipc_m.group(1))
        for key, ev in (
            ("cycles", "cycles"),
            ("instructions", "instructions"),
            ("cache_references", "cache-references"),
            ("cache_misses", "cache-misses"),
        ):
            mm = re.match(rf"^\s*([\d,]+)\s+{re.escape(pref)}{ev}/", line)
            if mm:
                out[key] = int(mm.group(1).replace(",", ""))

    if "cycles" not in out:
        for line in text.splitlines():
            line = line.strip()
            if "insn per cycle" in line and "ipc" not in out:
                ipc_m = re.search(r"#\s*([\d.]+)\s+insn per cycle", line)
                if ipc_m:
                    out["ipc"] = float(ipc_m.group(1))
            for key, pat in (
                ("cycles", r"^\s*([\d,]+)\s+cycles\b"),
                ("instructions", r"^\s*([\d,]+)\s+instructions\b"),
                ("cache_references", r"^\s*([\d,]+)\s+cache-references\b"),
                ("cache_misses", r"^\s*([\d,]+)\s+cache-misses\b"),
            ):
                if key in out:
                    continue
                mm = re.match(pat, line)
                if mm:
                    out[key] = int(mm.group(1).replace(",", ""))

    refs = out.get("cache_references")
    misses = out.get("cache_misses")
    if isinstance(refs, int) and refs > 0 and isinstance(misses, int):
        out["miss_frac_pct"] = 100.0 * misses / refs

    return out


def fmt_opt(x, fmt: str) -> str:
    if x is None:
        return NON_MEASURED
    if isinstance(x, float):
        return format(x, fmt)
    return str(x)


def main() -> int:
    if len(sys.argv) != 2:
        print(f"usage: {sys.argv[0]} <perf_results_dir>", file=sys.stderr)
        return 2
    root = Path(sys.argv[1])
    if not root.is_dir():
        print(f"error: not a directory: {root}", file=sys.stderr)
        return 1

    meta = parse_run_meta(root / "RUN_META.txt")
    rows: list[dict] = []
    for path in sorted(root.glob("*.perf.txt")):
        stem = path.name[: -len(".perf.txt")]
        tm = re.match(r"^(.+)_m(\d+)_t(\d+)$", stem)
        row = parse_perf_txt(path)
        if tm:
            row.setdefault("variant", tm.group(1))
            row.setdefault("m", int(tm.group(2)))
            row["threads"] = int(tm.group(3))
        rows.append(row)

    emit_header(
        "perf stat",
        root,
        "Per-call wall time from `perf_driver`; hardware counters from `perf stat`. "
        "Catalog grid uses the suite from `PERF_SUITE`; unrun cells are **non-measured**.",
    )
    emit_run_context(
        meta,
        prefer_keys=[
            "config_file",
            "date_utc",
            "quick",
            "PERF_MIN_RANGE",
            "PERF_MAX_RANGE",
            "PERF_STEP",
            "PERF_ITERS",
            "PERF_REPETITIONS",
            "THREAD_LIST",
            "PERF_SUITE",
        ],
        only_preferred=True,
    )

    merged: dict[tuple[str, int], dict[int, dict]] = defaultdict(dict)
    for r in rows:
        if "variant" not in r or "m" not in r:
            continue
        variant = str(r["variant"])
        m = int(r["m"])
        p = int(r["threads"])
        merged[(variant, m)][p] = r

    measured_threads = sorted({p for d in merged.values() for p in d})
    measured_variants = sorted({k[0] for k in merged})
    measured_ms = sorted({k[1] for k in merged})

    exp_ms = perf_expected_m_values(meta, measured_ms)
    exp_variants = perf_expected_variants(meta, measured_variants)
    cfg_threads = expected_thread_counts(meta, measured_threads)

    print("## Overview\n")
    print(f"- **Measured variants:** {len(measured_variants)}")
    print(f"- **Catalog variants ({meta.get('PERF_SUITE', 'cantor')}):** {len(exp_variants)}")
    print(f"- **Measured m:** {', '.join(str(x) for x in measured_ms) if measured_ms else '—'}")
    print(f"- **Configured m sweep:** {', '.join(str(x) for x in exp_ms)}")
    print(f"- **Measured thread counts:** {', '.join(str(p) for p in measured_threads) if measured_threads else '—'}")
    print()

    if not rows:
        print(f"## Results\n\n_{NON_MEASURED} — no `*.perf.txt` files found._\n")
        return 0

    print("## Per-call wall time (ms) — measured\n")
    hdr = "| Variant | m | " + " | ".join(f"P={p}" for p in measured_threads) + " |"
    sep = "|---------|---:|" + "|".join("------:" for _ in measured_threads) + "|"
    print(hdr)
    print(sep)
    for variant in measured_variants:
        for m in measured_ms:
            cells = []
            for p in measured_threads:
                hit = merged.get((variant, m), {}).get(p)
                cells.append(fmt_opt(hit.get("per_call_ms") if hit else None, ".3f"))
            print(f"| `{variant}` | {m} | " + " | ".join(cells) + " |")
    print()

    print("## Per-call wall time (ms) — catalog grid\n")
    print("Full suite × configured m × THREAD_LIST; missing cells are **non-measured**.\n")
    hdr2 = "| Variant | m | " + " | ".join(f"P={p}" for p in cfg_threads) + " |"
    sep2 = "|---------|---:|" + "|".join("------:" for _ in cfg_threads) + "|"
    print(hdr2)
    print(sep2)
    for variant in exp_variants:
        for m in exp_ms:
            cells = []
            for p in cfg_threads:
                if p not in measured_threads:
                    cells.append(NON_MEASURED)
                    continue
                hit = merged.get((variant, m), {}).get(p)
                cells.append(fmt_opt(hit.get("per_call_ms") if hit else None, ".3f"))
            print(f"| `{variant}` | {m} | " + " | ".join(cells) + " |")
    print()

    print("## Hardware counters — measured\n")
    print("| Variant | m | P | IPC | miss % | cycles | instructions |")
    print("|--------|---:|---:|---:|---:|---:|---:|")
    for variant in measured_variants:
        for m in measured_ms:
            for p in measured_threads:
                hit = merged.get((variant, m), {}).get(p)
                if not hit:
                    continue
                print(
                    "| `{v}` | {m} | {p} | {ipc} | {miss} | {cyc} | {insn} |".format(
                        v=variant,
                        m=m,
                        p=p,
                        ipc=fmt_opt(hit.get("ipc"), ".2f"),
                        miss=fmt_opt(hit.get("miss_frac_pct"), ".1f"),
                        cyc=fmt_opt(hit.get("cycles"), ",.0f"),
                        insn=fmt_opt(hit.get("instructions"), ",.0f"),
                    )
                )
    print()

    print("## Hardware counters — catalog grid\n")
    print("| Variant | m | P | IPC | miss % |")
    print("|--------|---:|---:|---:|---:|")
    for variant in exp_variants:
        for m in exp_ms:
            for p in cfg_threads:
                if p not in measured_threads:
                    ipc = miss = NON_MEASURED
                else:
                    hit = merged.get((variant, m), {}).get(p)
                    ipc = fmt_opt(hit.get("ipc") if hit else None, ".2f")
                    miss = fmt_opt(hit.get("miss_frac_pct") if hit else None, ".1f")
                print(f"| `{variant}` | {m} | {p} | {ipc} | {miss} |")
    print()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
