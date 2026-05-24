#!/usr/bin/env bash
# Run bench.sh and perf_stat.sh in one unattended sweep.
#
# Usage: ./scripts/run_all.sh [--config PATH]
# Guide: scripts/README.md

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
# shellcheck source=load_config.sh
source "${SCRIPT_DIR}/load_config.sh"

run_all_usage() {
  cat <<EOF
usage: run_all.sh [--config PATH]

Run bench/ and perf/ into measure_results/run_<timestamp>/.
Config: scripts/measure.conf (optional). See scripts/README.md.
EOF
}

__measure_rest=()
if ! measure_parse_config_args "$@"; then
  run_all_usage
  exit 0
fi
[[ ${#__measure_rest[@]} -eq 0 ]] || {
  echo "error: unknown argument(s): ${__measure_rest[*]}" >&2
  run_all_usage >&2
  exit 2
}

measure_load_config
measure_apply_defaults
measure_sync_tool_ranges
export MEASURE_CONFIG_LOADED=1

OUT_ROOT="${OUT:-${REPO_ROOT}/measure_results/run_$(date +%Y%m%d_%H%M%S)}"
BENCH_OUT="${OUT_ROOT}/bench"
PERF_OUT="${OUT_ROOT}/perf"
mkdir -p "${OUT_ROOT}"

bench_status="pending"
perf_status="pending"

{
  measure_write_run_meta_header "measurement-run"
  echo "SKIP_BENCH=${SKIP_BENCH:-0}"
  echo "SKIP_PERF=${SKIP_PERF:-0}"
  echo "PERF_SUITE=${PERF_SUITE:-cantor}"
  echo "TASKSET_CPUS=${TASKSET_CPUS:-}"
  echo "bench_dir=bench"
  echo "perf_dir=perf"
  echo "bench_status=${bench_status}"
  echo "perf_status=${perf_status}"
  measure_write_cpu_info
} | tee "${OUT_ROOT}/RUN_META.txt"

log() { echo "$*" | tee -a "${OUT_ROOT}/console.log"; }
log "=== run_all: ${OUT_ROOT} ==="

if [[ "${SKIP_BENCH:-0}" != "1" ]]; then
  log "--- bench ---"
  OUT="${BENCH_OUT}" MEASURE_CONFIG_LOADED=1 "${SCRIPT_DIR}/bench.sh" 2>&1 | tee -a "${OUT_ROOT}/console.log"
  bench_status="ok"
else
  bench_status="skipped (SKIP_BENCH=1)"
  log "bench skipped"
fi

if [[ "${SKIP_PERF:-0}" != "1" ]]; then
  if command -v perf >/dev/null 2>&1; then
    log "--- perf ---"
    OUT="${PERF_OUT}" MEASURE_CONFIG_LOADED=1 "${SCRIPT_DIR}/perf_stat.sh" 2>&1 | tee -a "${OUT_ROOT}/console.log"
    perf_status="ok"
  else
    perf_status="skipped (perf not in PATH)"
    log "perf skipped: perf not found"
  fi
else
  perf_status="skipped (SKIP_PERF=1)"
  log "perf skipped"
fi

sed -i '/^bench_status=/d;/^perf_status=/d' "${OUT_ROOT}/RUN_META.txt" 2>/dev/null || true
{
  echo "bench_status=${bench_status}"
  echo "perf_status=${perf_status}"
} >>"${OUT_ROOT}/RUN_META.txt"

command -v python3 >/dev/null 2>&1 || {
  echo "error: python3 required for summarize_run.py" >&2
  exit 1
}
python3 "${SCRIPT_DIR}/summarize_run.py" "${OUT_ROOT}" >"${OUT_ROOT}/SUMMARY.md"
log "Wrote ${OUT_ROOT}/SUMMARY.md"
log "Done. Artifacts: ${OUT_ROOT}"
