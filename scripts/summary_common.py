#!/usr/bin/env python3
"""Shared markdown helpers for benchmark and perf summaries."""

from __future__ import annotations

from pathlib import Path

NON_MEASURED = "non-measured"


def radix_label(k: int) -> str:
    """radix2k<K> uses butterfly radix 2^K (K=2 → radix-4, K=3 → radix-8, …)."""
    return f"radix-{2**k}"


def family_label(scope: str, path: str, k: int) -> str:
    """e.g. family_label('Cantor', 'affine', 2) → 'Cantor affine radix-4'."""
    return f"{scope} {path} {radix_label(k)}"


def parse_run_meta(path: Path) -> dict[str, str]:
    meta: dict[str, str] = {}
    if not path.is_file():
        return meta
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            continue
        key, _, val = line.partition("=")
        meta[key.strip()] = val.strip()
    return meta


def emit_header(tool: str, results_dir: Path, blurb: str) -> None:
    print(f"# Summary — {tool}\n")
    print(f"Directory: `{results_dir}`\n")
    print(blurb.rstrip())
    print("\n")


def emit_run_context(
    meta: dict[str, str],
    *,
    prefer_keys: list[str] | None = None,
    only_preferred: bool = False,
) -> None:
    print("## Run context\n")
    print("| Key | Value |")
    print("|-----|-------|")
    keys = prefer_keys or []
    seen: set[str] = set()
    for k in keys:
        if k in meta:
            print(f"| `{k}` | {meta[k]} |")
            seen.add(k)
    if not only_preferred:
        for k in sorted(meta):
            if k in seen:
                continue
            print(f"| `{k}` | {meta[k]} |")
    print()


def expected_m_values(meta: dict[str, str], measured: list[int]) -> list[int]:
    try:
        lo = int(meta.get("BM_MIN_RANGE", meta.get("MIN_RANGE", measured[0] if measured else 14)))
        hi = int(meta.get("BM_MAX_RANGE", meta.get("MAX_RANGE", lo)))
        step = int(meta.get("BM_STEP", meta.get("STEP", 2)))
        if step <= 0:
            step = 2
        return list(range(lo, hi + 1, step))
    except (ValueError, TypeError):
        return measured


def expected_thread_counts(meta: dict[str, str], measured: list[int]) -> list[int]:
    tl = meta.get("THREAD_LIST", "").strip()
    if not tl:
        return measured
    try:
        return sorted(int(x) for x in tl.split())
    except ValueError:
        return measured


def perf_expected_m_values(meta: dict[str, str], measured: list[int]) -> list[int]:
    try:
        lo = int(meta.get("PERF_MIN_RANGE", meta.get("MIN_RANGE", measured[0] if measured else 14)))
        hi = int(meta.get("PERF_MAX_RANGE", meta.get("MAX_RANGE", lo)))
        step = int(meta.get("PERF_STEP", meta.get("STEP", 2)))
        if step <= 0:
            step = 2
        return list(range(lo, hi + 1, step))
    except (ValueError, TypeError):
        return measured


PERF_VARIANTS_CANTOR = [
    "cantor_r2",
    "cantor_r2_par",
    "cantor_r2k2",
    "cantor_r2k2_par",
    "cantor_r2k3",
    "cantor_r2k3_par",
    "cantor_r2k4",
    "cantor_r2k4_par",
]

PERF_VARIANTS_LCH = [
    "r2",
    "r2_par",
    "r2k2",
    "r2k2_par",
    "r2k3",
    "r2k3_par",
    "r2k4",
    "r2k4_par",
    "r2k5",
    "r2k5_par",
]


def perf_expected_variants(meta: dict[str, str], measured: list[str]) -> list[str]:
    suite = meta.get("PERF_SUITE", "cantor")
    if suite == "lch":
        catalog = PERF_VARIANTS_LCH
    elif suite == "all":
        catalog = PERF_VARIANTS_CANTOR + PERF_VARIANTS_LCH
    else:
        catalog = PERF_VARIANTS_CANTOR
    # Preserve catalog order; append any measured variants not in catalog
    seen = set(catalog)
    extra = [v for v in measured if v not in seen]
    return catalog + extra
