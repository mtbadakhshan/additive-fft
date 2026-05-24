#!/usr/bin/env bash
# Shared config and defaults for bench.sh, perf_stat.sh, and run_all.sh.
#
# Precedence (lowest → highest):
#   1. Script defaults (measure_apply_defaults)
#   2. Config file (scripts/measure.conf or MEASURE_CONFIG)
#   3. Environment variables already set before the script runs
#
# --config PATH only selects which file to load; it does not override env vars.

measure_config_path() {
  echo "${MEASURE_CONFIG:-${SCRIPT_DIR}/measure.conf}"
}

measure_load_config() {
  [[ "${MEASURE_CONFIG_LOADED:-0}" == "1" ]] && return 0
  local config
  config="$(measure_config_path)"
  [[ -f "${config}" ]] || return 0

  local tmp line key
  tmp="$(mktemp)"
  while IFS= read -r line || [[ -n "${line}" ]]; do
    line="${line%%#*}"
    line="${line#"${line%%[![:space:]]*}"}"
    line="${line%"${line##*[![:space:]]}"}"
    [[ -z "${line}" ]] && continue
    key="${line%%=*}"
    key="${key#"${key%%[![:space:]]*}"}"
    key="${key%"${key##*[![:space:]]}"}"
    [[ "${key}" =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]] || continue
    [[ -v "${key}" ]] && continue
    printf '%s\n' "${line}" >>"${tmp}"
  done <"${config}"

  if [[ -s "${tmp}" ]]; then
    # shellcheck disable=SC1090
    set -a
    source "${tmp}"
    set +a
  fi
  rm -f "${tmp}"
}

measure_parse_config_args() {
  __measure_rest=()
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --config|-c)
        [[ $# -ge 2 ]] || {
          echo "error: $1 requires a path" >&2
          return 2
        }
        export MEASURE_CONFIG="$2"
        shift 2
        ;;
      -h|--help)
        return 1
        ;;
      --)
        shift
        __measure_rest+=("$@")
        return 0
        ;;
      *)
        __measure_rest+=("$1")
        shift
        ;;
    esac
  done
  return 0
}

measure_print_config_source() {
  local config
  config="$(measure_config_path)"
  if [[ -f "${config}" ]]; then
    echo "config_file=${config}"
  else
    echo "config_file=(none)"
  fi
}

# Shared defaults: problem size, OpenMP, build path.
# Optional overrides: BM_MIN_RANGE / PERF_MIN_RANGE (else copy from MIN_RANGE).
measure_apply_defaults() {
  export BUILD="${BUILD:-${REPO_ROOT}/C++/build}"
  export QUICK="${QUICK:-0}"
  export OMP_PLACES="${OMP_PLACES:-cores}"
  export OMP_PROC_BIND="${OMP_PROC_BIND:-close}"

  export MIN_RANGE="${MIN_RANGE:-14}"
  export MAX_RANGE="${MAX_RANGE:-22}"
  export STEP="${STEP:-2}"

  export BM_MIN_RANGE="${BM_MIN_RANGE:-${MIN_RANGE}}"
  export BM_MAX_RANGE="${BM_MAX_RANGE:-${MAX_RANGE}}"
  export BM_STEP="${BM_STEP:-${STEP}}"

  export PERF_MIN_RANGE="${PERF_MIN_RANGE:-${MIN_RANGE}}"
  export PERF_MAX_RANGE="${PERF_MAX_RANGE:-${MAX_RANGE}}"
  export PERF_STEP="${PERF_STEP:-${STEP}}"

  if [[ "${QUICK}" == "1" ]]; then
    export THREAD_LIST="${THREAD_LIST:-1 8}"
  else
    export THREAD_LIST="${THREAD_LIST:-1 2 4 8}"
  fi
}

measure_bench_filter_default() {
  echo '(BM_cantor_additive_fft|BM_cantor_additive_ifft)(/|_parallel|_precmp_basis|_radix2k)|BM_lch_additive_(fft|ifft).*precmp_basis'
}

measure_bench_timing() {
  if [[ "${QUICK:-0}" == "1" ]]; then
    BENCH_REP="${BENCH_REPETITIONS:-7}"
    BENCH_WARM="${BENCH_MIN_WARMUP_TIME:-0}"
    BENCH_EXTRA_WARM="${BENCH_WARMUP_ROUNDS:-0}"
  else
    BENCH_REP="${BENCH_REPETITIONS:-15}"
    BENCH_WARM="${BENCH_MIN_WARMUP_TIME:-1.0}"
    BENCH_EXTRA_WARM="${BENCH_WARMUP_ROUNDS:-1}"
  fi
}

measure_perf_timing() {
  export PERF_EVENTS="${PERF_EVENTS:-cycles,instructions,cache-references,cache-misses}"
  export PERF_SUITE="${PERF_SUITE:-cantor}"
  if [[ "${QUICK:-0}" == "1" ]]; then
    PERF_REP="${PERF_REPETITIONS:-1}"
    PERF_ITERS_DEFAULT=10
    PERF_WARMUP_ROUNDS="${PERF_WARMUP_ROUNDS:-0}"
    PERF_WARMUP_ITERS="${PERF_WARMUP_ITERS:-2}"
  else
    PERF_REP="${PERF_REPETITIONS:-3}"
    PERF_ITERS_DEFAULT=25
    PERF_WARMUP_ROUNDS="${PERF_WARMUP_ROUNDS:-1}"
    PERF_WARMUP_ITERS="${PERF_WARMUP_ITERS:-5}"
  fi
}

measure_write_run_meta_header() {
  local tool="$1"
  echo "tool=${tool}"
  measure_print_config_source
  echo "date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "repo=${REPO_ROOT}"
  echo "build=${BUILD}"
  echo "quick=${QUICK:-0}"
  echo "MIN_RANGE=${MIN_RANGE}"
  echo "MAX_RANGE=${MAX_RANGE}"
  echo "STEP=${STEP}"
  echo "THREAD_LIST=${THREAD_LIST}"
  echo "OMP_PLACES=${OMP_PLACES}"
  echo "OMP_PROC_BIND=${OMP_PROC_BIND}"
}

measure_write_cpu_info() {
  {
    echo ""
    uname -a || true
    grep -m1 'model name' /proc/cpuinfo 2>/dev/null || true
  }
}

# run_all.sh: keep BM_* and PERF_* aligned with MIN_RANGE (ignore stale shell exports).
measure_sync_tool_ranges() {
  export BM_MIN_RANGE="${MIN_RANGE}"
  export BM_MAX_RANGE="${MAX_RANGE}"
  export BM_STEP="${STEP}"
  export PERF_MIN_RANGE="${MIN_RANGE}"
  export PERF_MAX_RANGE="${MAX_RANGE}"
  export PERF_STEP="${STEP}"
}
