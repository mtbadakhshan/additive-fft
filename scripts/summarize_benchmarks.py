#!/usr/bin/env python3
"""
Merge Google Benchmark JSON files produced by bench.sh into markdown tables.
"""

from __future__ import annotations

import json
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from summary_common import (
    NON_MEASURED,
    emit_header,
    emit_run_context,
    expected_m_values,
    expected_thread_counts,
    family_label,
    parse_run_meta,
    radix_label,
)

SERIAL_PARALLEL_CATALOG: list[tuple[str, str, str]] = [
    ("BM_cantor_additive_fft", "BM_cantor_additive_fft_parallel", family_label("Cantor", "affine", 1)),
    (
        "BM_cantor_additive_fft_radix2k<2>",
        "BM_cantor_additive_fft_radix2k_parallel<2>",
        family_label("Cantor", "affine", 2),
    ),
    (
        "BM_cantor_additive_fft_radix2k<3>",
        "BM_cantor_additive_fft_radix2k_parallel<3>",
        family_label("Cantor", "affine", 3),
    ),
    (
        "BM_cantor_additive_fft_radix2k<4>",
        "BM_cantor_additive_fft_radix2k_parallel<4>",
        family_label("Cantor", "affine", 4),
    ),
    (
        "BM_cantor_additive_fft_precmp_basis",
        "BM_cantor_additive_fft_precmp_basis_parallel",
        family_label("Cantor", "table", 1),
    ),
    (
        "BM_cantor_additive_fft_precmp_basis_radix2k<2>",
        "BM_cantor_additive_fft_precmp_basis_radix2k_parallel<2>",
        family_label("Cantor", "table", 2),
    ),
    (
        "BM_cantor_additive_fft_precmp_basis_radix2k<3>",
        "BM_cantor_additive_fft_precmp_basis_radix2k_parallel<3>",
        family_label("Cantor", "table", 3),
    ),
    (
        "BM_cantor_additive_fft_precmp_basis_radix2k<4>",
        "BM_cantor_additive_fft_precmp_basis_radix2k_parallel<4>",
        family_label("Cantor", "table", 4),
    ),
    (
        "BM_lch_additive_fft_precmp_basis",
        "BM_lch_additive_fft_parallel_precmp_basis",
        family_label("LCH", "table", 1),
    ),
    (
        "BM_lch_additive_fft_radix2k_precmp_basis<2>",
        "BM_lch_additive_fft_radix2k_parallel_precmp_basis<2>",
        family_label("LCH", "table", 2),
    ),
    (
        "BM_lch_additive_fft_radix2k_precmp_basis<3>",
        "BM_lch_additive_fft_radix2k_parallel_precmp_basis<3>",
        family_label("LCH", "table", 3),
    ),
    (
        "BM_lch_additive_fft_radix2k_precmp_basis<4>",
        "BM_lch_additive_fft_radix2k_parallel_precmp_basis<4>",
        family_label("LCH", "table", 4),
    ),
    (
        "BM_lch_additive_fft_radix2k_precmp_basis<5>",
        "BM_lch_additive_fft_radix2k_parallel_precmp_basis<5>",
        family_label("LCH", "table", 5),
    ),
]

FFT_IFFT_CATALOG: list[tuple[str, str, str]] = [
    ("BM_cantor_additive_fft", "BM_cantor_additive_ifft", family_label("Cantor", "affine", 1)),
    (
        "BM_cantor_additive_fft_radix2k<2>",
        "BM_cantor_additive_ifft_radix2k<2>",
        family_label("Cantor", "affine", 2),
    ),
    (
        "BM_cantor_additive_fft_radix2k<3>",
        "BM_cantor_additive_ifft_radix2k<3>",
        family_label("Cantor", "affine", 3),
    ),
    (
        "BM_cantor_additive_fft_radix2k<4>",
        "BM_cantor_additive_ifft_radix2k<4>",
        family_label("Cantor", "affine", 4),
    ),
    (
        "BM_cantor_additive_fft_precmp_basis",
        "BM_cantor_additive_ifft_precmp_basis",
        family_label("Cantor", "table", 1),
    ),
    (
        "BM_cantor_additive_fft_precmp_basis_radix2k<2>",
        "BM_cantor_additive_ifft_precmp_basis_radix2k<2>",
        family_label("Cantor", "table", 2),
    ),
    (
        "BM_cantor_additive_fft_precmp_basis_radix2k<3>",
        "BM_cantor_additive_ifft_precmp_basis_radix2k<3>",
        family_label("Cantor", "table", 3),
    ),
    (
        "BM_cantor_additive_fft_precmp_basis_radix2k<4>",
        "BM_cantor_additive_ifft_precmp_basis_radix2k<4>",
        family_label("Cantor", "table", 4),
    ),
    ("BM_lch_additive_fft_precmp_basis", "BM_lch_additive_ifft_precmp_basis", family_label("LCH", "table", 1)),
    (
        "BM_lch_additive_fft_radix2k_precmp_basis<2>",
        "BM_lch_additive_ifft_radix2k_precmp_basis<2>",
        family_label("LCH", "table", 2),
    ),
    (
        "BM_lch_additive_fft_radix2k_precmp_basis<3>",
        "BM_lch_additive_ifft_radix2k_precmp_basis<3>",
        family_label("LCH", "table", 3),
    ),
    (
        "BM_lch_additive_fft_radix2k_precmp_basis<4>",
        "BM_lch_additive_ifft_radix2k_precmp_basis<4>",
        family_label("LCH", "table", 4),
    ),
    (
        "BM_lch_additive_fft_radix2k_precmp_basis<5>",
        "BM_lch_additive_ifft_radix2k_precmp_basis<5>",
        family_label("LCH", "table", 5),
    ),
]

COMPACT_BLOCKS: list[tuple[str, list[tuple[str, str, str]]]] = [
    (
        "Cantor — affine chart",
        [
            (radix_label(1), "BM_cantor_additive_fft", "BM_cantor_additive_fft_parallel"),
            (radix_label(2), "BM_cantor_additive_fft_radix2k<2>", "BM_cantor_additive_fft_radix2k_parallel<2>"),
        ],
    ),
    (
        "Cantor — table path",
        [
            (radix_label(1), "BM_cantor_additive_fft_precmp_basis", "BM_cantor_additive_fft_precmp_basis_parallel"),
            (
                radix_label(2),
                "BM_cantor_additive_fft_precmp_basis_radix2k<2>",
                "BM_cantor_additive_fft_precmp_basis_radix2k_parallel<2>",
            ),
        ],
    ),
    (
        "LCH — table path",
        [
            (radix_label(1), "BM_lch_additive_fft_precmp_basis", "BM_lch_additive_fft_parallel_precmp_basis"),
            (
                radix_label(2),
                "BM_lch_additive_fft_radix2k_precmp_basis<2>",
                "BM_lch_additive_fft_radix2k_parallel_precmp_basis<2>",
            ),
        ],
    ),
]


def load_medians(path: Path) -> tuple[dict, dict]:
    data = json.loads(path.read_text())
    ctx = data.get("context", {})
    medians: dict[str, dict[str, float]] = {}
    for b in data.get("benchmarks", []):
        if b.get("run_type") != "aggregate":
            continue
        if b.get("aggregate_name") != "median":
            continue
        if b.get("aggregate_unit") != "time":
            continue
        name = b.get("run_name")
        if not name:
            continue
        medians[name] = {
            "real_us": float(b["real_time"]),
            "cpu_us": float(b["cpu_time"]),
        }
    if not medians:
        for b in data.get("benchmarks", []):
            if b.get("run_type") != "iteration":
                continue
            name = b.get("run_name")
            if not name:
                continue
            medians[name] = {
                "real_us": float(b["real_time"]),
                "cpu_us": float(b["cpu_time"]),
            }
    return ctx, medians


def omp_threads_from_context(ctx: dict) -> int:
    v = ctx.get("OMP_NUM_THREADS_env", "(unset)")
    if v == "(unset)":
        return 1
    try:
        return int(v)
    except ValueError:
        return -1


def m_values(merged: dict) -> list[int]:
    out: set[int] = set()
    for nm in merged:
        m = re.search(r"/(\d+)$", nm)
        if m:
            out.add(int(m.group(1)))
    return sorted(out)


def wall_us(merged: dict, run_base: str, m_val: int, p_val: int) -> str:
    key = f"{run_base}/{m_val}"
    t = merged.get(key, {}).get(p_val)
    if not t or t["real_us"] <= 0:
        return NON_MEASURED
    return f"{t['real_us']:.0f}"


def emit_table(headers: list[str], rows: list[list[str]]) -> None:
    if not rows:
        return
    align = [":---" if i == 0 else "---:" for i in range(len(headers))]
    print("| " + " | ".join(headers) + " |")
    print("| " + " | ".join(align) + " |")
    for row in rows:
        print("| " + " | ".join(row) + " |")
    print()


def emit_time_tables(merged: dict, names: list[str], thread_counts: list[int]) -> None:
    hdr = ["Benchmark"] + [f"P={p}" for p in thread_counts]
    wall_rows = []
    cpu_rows = []
    for name in names:
        wcells = []
        ccells = []
        for p in thread_counts:
            t = merged[name].get(p)
            wcells.append(f"{t['real_us']:.0f}" if t else NON_MEASURED)
            ccells.append(f"{t['cpu_us']:.0f}" if t else NON_MEASURED)
        wall_rows.append([f"`{name}`"] + wcells)
        cpu_rows.append([f"`{name}`"] + ccells)

    print("## Median wall time (µs) — measured\n")
    emit_table(hdr, wall_rows)

    print("## Median CPU time (µs) — measured\n")
    emit_table(hdr, cpu_rows)


def emit_speedup_table(
    merged: dict,
    names: list[str],
    measured_threads: list[int],
    configured_threads: list[int],
) -> None:
    others = [p for p in configured_threads if p != 1]
    print("## Speedup vs P=1 (wall time)\n")
    if not others:
        print(f"_{NON_MEASURED} — `THREAD_LIST` has only P=1._\n")
        return
    print("`speedup = real_us(P=1) / real_us(P)`; missing runs marked **non-measured**.\n")
    rows = []
    for name in names:
        base = merged[name].get(1)
        cells = []
        for p in others:
            if p not in measured_threads:
                cells.append(NON_MEASURED)
                continue
            t = merged[name].get(p)
            if not base or base["real_us"] <= 0 or not t or t["real_us"] <= 0:
                cells.append(NON_MEASURED)
            else:
                cells.append(f"{base['real_us'] / t['real_us']:.2f}×")
        rows.append([f"`{name}`"] + cells)
    emit_table(["Benchmark"] + [f"P={p}" for p in others], rows)


def emit_fft_ifft_catalog(
    merged: dict,
    exp_ms: list[int],
    thread_counts: list[int],
    measured_threads: list[int],
) -> None:
    print("## FFT vs IFFT (catalog)\n")
    print("Pairs from the default suite; missing runs marked **non-measured**.\n")
    rows = []
    for fft_base, ifft_base, label in FFT_IFFT_CATALOG:
        for m in exp_ms:
            for p in thread_counts:
                if p not in measured_threads:
                    fft_v = ifft_v = ratio = NON_MEASURED
                else:
                    tf = merged.get(f"{fft_base}/{m}", {}).get(p)
                    ti = merged.get(f"{ifft_base}/{m}", {}).get(p)
                    if tf and ti and tf["real_us"] > 0:
                        fft_v = f"{tf['real_us']:.0f}"
                        ifft_v = f"{ti['real_us']:.0f}"
                        ratio = f"{ti['real_us'] / tf['real_us']:.2f}×"
                    else:
                        fft_v = f"{tf['real_us']:.0f}" if tf else NON_MEASURED
                        ifft_v = f"{ti['real_us']:.0f}" if ti else NON_MEASURED
                        ratio = NON_MEASURED
                rows.append([label, str(m), str(p), fft_v, ifft_v, ratio])
    emit_table(["Family", "m", "P", "FFT µs", "IFFT µs", "IFFT/FFT"], rows)


def emit_serial_parallel_catalog(
    merged: dict,
    exp_ms: list[int],
    max_p: int,
    measured_threads: list[int],
) -> None:
    print(f"## Serial vs parallel (catalog, parallel P={max_p})\n")
    print("Default suite pairs; missing runs marked **non-measured**.\n")
    rows = []
    for ser_base, par_base, label in SERIAL_PARALLEL_CATALOG:
        for m in exp_ms:
            ts = merged.get(f"{ser_base}/{m}", {}).get(1)
            if max_p not in measured_threads:
                tp = None
            else:
                tp = merged.get(f"{par_base}/{m}", {}).get(max_p)
            ser_v = f"{ts['real_us']:.0f}" if ts else NON_MEASURED
            par_v = f"{tp['real_us']:.0f}" if tp else NON_MEASURED
            if ts and tp and tp["real_us"] > 0:
                spd = f"{ts['real_us'] / tp['real_us']:.2f}×"
            else:
                spd = NON_MEASURED
            rows.append([label, str(m), ser_v, par_v, spd])
    emit_table(["Family", "m", "serial µs (P=1)", f"parallel µs (P={max_p})", "speedup"], rows)


def emit_compact_catalog(
    merged: dict,
    exp_ms: list[int],
    p_lo: int,
    p_hi: int,
    measured_threads: list[int],
) -> None:
    print("## Compact comparison (catalog)\n")
    print("Default suite families; missing runs marked **non-measured**.\n")
    for title, pairs in COMPACT_BLOCKS:
        print(f"### {title}\n")
        block_rows = []
        for label, left_bm, right_bm in pairs:
            for m in exp_ms:
                a = wall_us(merged, left_bm, m, p_lo) if p_lo in measured_threads else NON_MEASURED
                b = wall_us(merged, right_bm, m, p_hi) if p_hi in measured_threads else NON_MEASURED
                block_rows.append([label, str(m), a, b])
        emit_table([f"variant (P={p_lo})", "m", "serial µs", f"parallel µs (P={p_hi})"], block_rows)


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: summarize_benchmarks.py <results_dir>", file=sys.stderr)
        return 2
    out_dir = Path(sys.argv[1])
    files = sorted(out_dir.glob("threads_*.json"))
    if not files:
        print(f"no threads_*.json under {out_dir}", file=sys.stderr)
        return 1

    merged: dict[str, dict[int, dict[str, float]]] = defaultdict(dict)
    contexts: dict[int, dict] = {}

    for fp in files:
        ctx, med = load_medians(fp)
        p = omp_threads_from_context(ctx)
        if p < 0:
            m = re.search(r"threads_(\d+)\.json$", fp.name)
            p = int(m.group(1)) if m else -1
        contexts[p] = ctx
        for run_name, times in med.items():
            merged[run_name][p] = times

    measured_threads = sorted({p for d in merged.values() for p in d})
    if not measured_threads:
        print(f"no benchmark rows found under {out_dir}", file=sys.stderr)
        return 1

    names = sorted(merged.keys())
    measured_ms = m_values(merged)
    meta = parse_run_meta(out_dir / "RUN_META.txt")
    exp_ms = expected_m_values(meta, measured_ms)
    cfg_threads = expected_thread_counts(meta, measured_threads)
    max_p = max(cfg_threads) if cfg_threads else max(measured_threads)
    p_lo = 1 if 1 in cfg_threads else min(cfg_threads)

    emit_header(
        "google-benchmark",
        out_dir,
        "Median wall/CPU time per benchmark iteration (microseconds). "
        "Catalog sections list the default suite; unrun cells are **non-measured**.",
    )

    emit_run_context(
        meta,
        prefer_keys=[
            "config_file",
            "date_utc",
            "quick",
            "MIN_RANGE",
            "MAX_RANGE",
            "STEP",
            "BM_MIN_RANGE",
            "BM_MAX_RANGE",
            "BM_STEP",
            "THREAD_LIST",
            "BENCH_REPETITIONS",
            "BENCH_FILTER",
        ],
        only_preferred=True,
    )

    if thread_counts := measured_threads:
        ctx = contexts.get(thread_counts[0], {})
        host_bits = []
        for k in ("host_name", "num_cpus", "mhz_per_cpu", "OMP_NUM_THREADS_env"):
            if k in ctx:
                host_bits.append(f"{k}={ctx[k]}")
        if host_bits:
            print("## Machine (from benchmark JSON)\n")
            print("- " + "\n- ".join(host_bits))
            print()

    print("## Overview\n")
    print(f"- **Measured benchmarks:** {len(names)}")
    print(f"- **Measured m:** {', '.join(str(x) for x in measured_ms) if measured_ms else '—'}")
    print(f"- **Configured m sweep:** {', '.join(str(x) for x in exp_ms)}")
    print(f"- **Measured thread counts:** {', '.join(str(p) for p in measured_threads)}")
    print(f"- **Configured THREAD_LIST:** {meta.get('THREAD_LIST', '—')}")
    print()

    emit_time_tables(merged, names, measured_threads)
    emit_fft_ifft_catalog(merged, exp_ms, cfg_threads, measured_threads)
    emit_speedup_table(merged, names, measured_threads, cfg_threads)
    emit_serial_parallel_catalog(merged, exp_ms, max_p, measured_threads)
    emit_compact_catalog(merged, exp_ms, p_lo, max_p, measured_threads)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
