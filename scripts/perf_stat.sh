#!/usr/bin/env bash
# Run C++/perf_driver under Linux `perf stat` (hardware counters + wall time).
#
# Full guide (Linux packages, build, cache counters): scripts/README.md
#
# Usage (from repo root):
#   ./scripts/perf_stat.sh                    # full sweep (variants × m × threads)
#   QUICK=1 ./scripts/perf_stat.sh            # faster full sweep
#   ./scripts/perf_stat.sh <variant> <m> <iters>
#   ./scripts/perf_stat.sh --compare-cantor [m] [iters]
#   THREAD_LIST="1 8" ./scripts/perf_stat.sh cantor_r2k4_par 22 25
#
# Environment (aligned with bench.sh where sensible):
#   BUILD, OUT, QUICK, THREAD_LIST, OMP_PLACES, OMP_PROC_BIND, TASKSET_CPUS
#   PERF_MIN_RANGE / PERF_MAX_RANGE / PERF_STEP — m sweep (default 14 22 2)
#   PERF_ITERS          — FFT calls per perf measurement (default 25; 10 if QUICK)
#   PERF_REPETITIONS    — perf stat -r (default 3; 1 if QUICK)
#   PERF_WARMUP_ROUNDS  — dry perf_driver runs before each measurement (default 1; 0 if QUICK)
#   PERF_WARMUP_ITERS   — iters per warmup dry run (default 5; 2 if QUICK)
#   PERF_EVENTS         — perf stat -e list
#   PERF_SUITE          — cantor (default) | lch | all
#
# Outputs:
#   RUN_META.txt, console.log, <variant>_m<m>_t<threads>.perf.txt, SUMMARY.md

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
# shellcheck source=load_config.sh
source "${SCRIPT_DIR}/load_config.sh"
CPP_ROOT="${REPO_ROOT}/C++"
BUILD="${BUILD:-${CPP_ROOT}/build}"
PERF_DRIVER="${BUILD}/perf_driver"

usage() {
  cat <<'EOF'
usage:
  perf_stat.sh [--config PATH]              # full sweep (default Cantor radix family)
  perf_stat.sh [--config PATH] <variant> <m> <iters>
  perf_stat.sh --variants "<v1> <v2> ...>" [<m> <iters>]
  perf_stat.sh --compare-cantor [<m> <iters>]
  perf_stat.sh --compare-lch [<m> <iters>]

Config: scripts/measure.conf (copy from measure.conf.example) if present.
Precedence: script defaults < config file < environment variables.

Full sweep (no args): all default Cantor variants, m in [PERF_MIN_RANGE..PERF_MAX_RANGE]
step PERF_STEP, THREAD_LIST, with warmup + perf repetitions. See scripts/README.md.

Variants (perf_driver):
  LCH:    r2, r2k1..r2k5, r2_par, r2k1_par..r2k5_par
  Cantor: cantor_r2, cantor_r2_par, cantor_r2k2|3|4, cantor_r2k2_par|3_par|4_par
EOF
}

die() {
  echo "error: $*" >&2
  exit 1
}

have_perf() {
  command -v perf >/dev/null 2>&1
}

variants_cantor() {
  CANTOR_VARIANTS=(
    cantor_r2 cantor_r2_par
    cantor_r2k2 cantor_r2k2_par
    cantor_r2k3 cantor_r2k3_par
    cantor_r2k4 cantor_r2k4_par
  )
}

variants_lch() {
  LCH_VARIANTS=(
    r2 r2_par
    r2k2 r2k2_par r2k3 r2k3_par r2k4 r2k4_par r2k5 r2k5_par
  )
}

variants_all() {
  variants_cantor
  variants_lch
  VARIANTS=("${CANTOR_VARIANTS[@]}" "${LCH_VARIANTS[@]}")
}

select_suite_variants() {
  case "${PERF_SUITE:-cantor}" in
    cantor) variants_cantor; VARIANTS=("${CANTOR_VARIANTS[@]}") ;;
    lch) variants_lch; VARIANTS=("${LCH_VARIANTS[@]}") ;;
    all) variants_all ;;
    *) die "unknown PERF_SUITE=${PERF_SUITE} (use cantor, lch, or all)" ;;
  esac
}

gen_m_list() {
  M_LIST=()
  local m
  for ((m = PERF_MIN_RANGE; m <= PERF_MAX_RANGE; m += PERF_STEP)); do
    M_LIST+=("${m}")
  done
}

run_warmup() {
  local variant="$1" m="$2" p="$3"
  local round
  [[ "${PERF_WARMUP_ROUNDS}" -gt 0 ]] || return 0
  for ((round = 1; round <= PERF_WARMUP_ROUNDS; round++)); do
    echo "--- warmup ${round}/${PERF_WARMUP_ROUNDS} ${variant} m=${m} t=${p} ---" \
      >>"${OUT}/console.log"
    "${PERF_DRIVER}" "${variant}" "${m}" "${PERF_WARMUP_ITERS}" \
      >>"${OUT}/console.log" 2>&1 || true
  done
}

run_measurement() {
  local variant="$1" m="$2" p="$3" iters="$4"
  local tag="${variant}_m${m}_t${p}"
  local out_file="${OUT}/${tag}.perf.txt"
  echo "=== ${tag} (OMP_NUM_THREADS=${p}) ===" | tee -a "${OUT}/console.log"

  run_warmup "${variant}" "${m}" "${p}"

  {
    echo "# variant=${variant} m=${m} iters=${iters} OMP_NUM_THREADS=${p}"
    echo "# perf stat -r ${PERF_REP} -e ${PERF_EVENTS}"
    echo ""
  } >"${out_file}"

  if [[ ${#TASKSET_PREFIX[@]} -gt 0 ]]; then
    "${TASKSET_PREFIX[@]}" perf stat -r "${PERF_REP}" -e "${PERF_EVENTS}" -- \
      "${PERF_DRIVER}" "${variant}" "${m}" "${iters}" >>"${out_file}" 2>&1
  else
    perf stat -r "${PERF_REP}" -e "${PERF_EVENTS}" -- \
      "${PERF_DRIVER}" "${variant}" "${m}" "${iters}" >>"${out_file}" 2>&1
  fi
  tee -a "${OUT}/console.log" <"${out_file}" >/dev/null
}

VARIANTS=()
COMPARE_MODE=""
FULL_SWEEP=0
POSITIONAL=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      usage
      exit 0
      ;;
    --config|-c)
      [[ $# -ge 2 ]] || die "--config requires a path"
      export MEASURE_CONFIG="$2"
      shift 2
      ;;
    --full)
      FULL_SWEEP=1
      shift
      ;;
    --variants)
      shift
      [[ $# -ge 1 ]] || die "--variants requires a quoted list"
      read -r -a VARIANTS <<<"$1"
      shift
      ;;
    --compare-cantor)
      COMPARE_MODE="cantor"
      shift
      ;;
    --compare-lch)
      COMPARE_MODE="lch"
      shift
      ;;
    --)
      shift
      POSITIONAL+=("$@")
      break
      ;;
    -*)
      die "unknown option: $1"
      ;;
    *)
      POSITIONAL+=("$1")
      shift
      ;;
  esac
done

measure_load_config
measure_apply_defaults
measure_sync_tool_ranges
measure_perf_timing

PERF_EVENTS="${PERF_EVENTS}"
OUT="${OUT:-${REPO_ROOT}/perf_results/run_$(date +%Y%m%d_%H%M%S)}"

if [[ ${#VARIANTS[@]} -eq 0 ]]; then
  case "${COMPARE_MODE}" in
    cantor) variants_cantor; VARIANTS=("${CANTOR_VARIANTS[@]}") ;;
    lch) variants_lch; VARIANTS=("${LCH_VARIANTS[@]}") ;;
    "")
      if [[ ${#POSITIONAL[@]} -ge 1 ]]; then
        VARIANTS=("${POSITIONAL[0]}")
        POSITIONAL=("${POSITIONAL[@]:1}")
      else
        FULL_SWEEP=1
        select_suite_variants
      fi
      ;;
  esac
fi

# Single-m vs multi-m mode
M_LIST=()
ITERS="${PERF_ITERS:-${PERF_ITERS_DEFAULT}}"

if [[ ${FULL_SWEEP} -eq 1 ]]; then
  gen_m_list
elif [[ ${#POSITIONAL[@]} -ge 2 ]]; then
  M_LIST=("${POSITIONAL[0]}")
  ITERS="${POSITIONAL[1]}"
elif [[ ${#POSITIONAL[@]} -eq 1 ]]; then
  M_LIST=("${POSITIONAL[0]}")
  # iters: default
elif [[ -n "${COMPARE_MODE}" ]]; then
  gen_m_list
else
  usage >&2
  exit 1
fi

if [[ ${#VARIANTS[@]} -eq 0 ]]; then
  usage >&2
  exit 1
fi

if ! have_perf; then
  die "perf not found in PATH (install linux-perf)"
fi

if [[ ! -x "${PERF_DRIVER}" ]]; then
  die "${PERF_DRIVER} not found — build with: cmake --build ${BUILD} --target perf_driver"
fi

if [[ -r /proc/sys/kernel/perf_event_paranoid ]]; then
  paranoia="$(cat /proc/sys/kernel/perf_event_paranoid)"
  if [[ "${paranoia}" -gt 2 ]]; then
    echo "warning: perf_event_paranoid=${paranoia} may block user perf; try sudo or lower the sysctl" >&2
  fi
fi

mkdir -p "${OUT}"

{
  measure_write_run_meta_header "perf-stat"
  echo "perf_driver=${PERF_DRIVER}"
  echo "full_sweep=${FULL_SWEEP}"
  echo "PERF_MIN_RANGE=${PERF_MIN_RANGE}"
  echo "PERF_MAX_RANGE=${PERF_MAX_RANGE}"
  echo "PERF_STEP=${PERF_STEP}"
  echo "m_list=${M_LIST[*]}"
  echo "PERF_ITERS=${ITERS}"
  echo "PERF_EVENTS=${PERF_EVENTS}"
  echo "PERF_REPETITIONS=${PERF_REP}"
  echo "PERF_WARMUP_ROUNDS=${PERF_WARMUP_ROUNDS}"
  echo "PERF_WARMUP_ITERS=${PERF_WARMUP_ITERS}"
  echo "TASKSET_CPUS=${TASKSET_CPUS:-}"
  echo "PERF_SUITE=${PERF_SUITE}"
  echo "variants=${VARIANTS[*]}"
  perf --version 2>/dev/null || true
  measure_write_cpu_info
} | tee "${OUT}/RUN_META.txt"

TASKSET_PREFIX=()
if [[ -n "${TASKSET_CPUS:-}" ]]; then
  TASKSET_PREFIX=(taskset -c "${TASKSET_CPUS}")
fi

n_m=${#M_LIST[@]}
n_var=${#VARIANTS[@]}
n_runs=$((n_m * n_var))
echo "Planning ${n_runs} variant×m configurations × thread counts in THREAD_LIST=${THREAD_LIST}" \
  | tee -a "${OUT}/console.log"

for P in ${THREAD_LIST}; do
  export OMP_NUM_THREADS="${P}"
  for m in "${M_LIST[@]}"; do
    for variant in "${VARIANTS[@]}"; do
      run_measurement "${variant}" "${m}" "${P}" "${ITERS}"
    done
  done
done

if ! command -v python3 >/dev/null 2>&1; then
  die "python3 required for summarize_perf.py"
fi
python3 "${SCRIPT_DIR}/summarize_perf.py" "${OUT}" >"${OUT}/SUMMARY.md"
echo "Wrote ${OUT}/SUMMARY.md" | tee -a "${OUT}/console.log"
echo "Done. Artifacts: ${OUT}"
