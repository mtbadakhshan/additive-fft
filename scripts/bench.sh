#!/usr/bin/env bash
# Run additive-FFT Google Benchmarks (Cantor + LCH, serial and OpenMP parallel).
#
# Usage: ./scripts/bench.sh [--config PATH]
# Guide: scripts/README.md

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
# shellcheck source=load_config.sh
source "${SCRIPT_DIR}/load_config.sh"

bench_usage() {
  cat <<EOF
usage: bench.sh [--config PATH]

Wall-clock timings via Google Benchmark. Config: scripts/measure.conf (optional).
See scripts/README.md for defaults and measure.conf.example.
EOF
}

__measure_rest=()
if ! measure_parse_config_args "$@"; then
  bench_usage
  exit 0
fi
[[ ${#__measure_rest[@]} -eq 0 ]] || {
  echo "error: unknown argument(s): ${__measure_rest[*]}" >&2
  bench_usage >&2
  exit 2
}

measure_load_config
measure_apply_defaults
measure_sync_tool_ranges
measure_bench_timing

BUILD="${BUILD}"
RUN_BENCH="${BUILD}/run_benchmark"
[[ -x "${RUN_BENCH}" ]] || {
  echo "error: ${RUN_BENCH} not found (cmake --build ${BUILD} --target run_benchmark)" >&2
  exit 1
}

export BENCH_FILTER="${BENCH_FILTER:-$(measure_bench_filter_default)}"
OUT="${OUT:-${REPO_ROOT}/benchmark_results/run_$(date +%Y%m%d_%H%M%S)}"
mkdir -p "${OUT}"

{
  measure_write_run_meta_header "google-benchmark"
  echo "BM_MIN_RANGE=${BM_MIN_RANGE}"
  echo "BM_MAX_RANGE=${BM_MAX_RANGE}"
  echo "BM_STEP=${BM_STEP}"
  echo "BENCH_REPETITIONS=${BENCH_REP}"
  echo "BENCH_MIN_WARMUP_TIME=${BENCH_WARM}"
  echo "BENCH_WARMUP_ROUNDS=${BENCH_EXTRA_WARM}"
  printf 'BENCH_FILTER=%q\n' "${BENCH_FILTER}"
  measure_write_cpu_info
} | tee "${OUT}/RUN_META.txt"

for P in ${THREAD_LIST}; do
  export OMP_NUM_THREADS="${P}"
  echo "=== OMP_NUM_THREADS=${P} ===" | tee -a "${OUT}/console.log"

  if [[ "${BENCH_EXTRA_WARM}" != "0" ]]; then
    for _w in $(seq 1 "${BENCH_EXTRA_WARM}"); do
      echo "--- warmup ${_w}/${BENCH_EXTRA_WARM} ---" | tee -a "${OUT}/console.log"
      "${RUN_BENCH}" \
        --benchmark_filter="${BENCH_FILTER}" \
        --benchmark_min_time=0.1s \
        --benchmark_repetitions=1 \
        --benchmark_report_aggregates_only=true \
        >>"${OUT}/console.log" 2>&1
    done
  fi

  if ! "${RUN_BENCH}" \
    --benchmark_filter="${BENCH_FILTER}" \
    --benchmark_repetitions="${BENCH_REP}" \
    --benchmark_min_warmup_time="${BENCH_WARM}" \
    --benchmark_report_aggregates_only=true \
    --benchmark_format=json \
    --benchmark_out="${OUT}/threads_${P}.json" \
    >/dev/null 2>>"${OUT}/console.log"; then
    echo "error: run_benchmark failed (OMP_NUM_THREADS=${P})" >&2
    exit 1
  fi
  echo "Wrote ${OUT}/threads_${P}.json" | tee -a "${OUT}/console.log"
done

command -v python3 >/dev/null 2>&1 || {
  echo "error: python3 required for summarize_benchmarks.py" >&2
  exit 1
}
python3 "${SCRIPT_DIR}/summarize_benchmarks.py" "${OUT}" >"${OUT}/SUMMARY.md"
echo "Wrote ${OUT}/SUMMARY.md" | tee -a "${OUT}/console.log"
echo "Done. Artifacts: ${OUT}"
