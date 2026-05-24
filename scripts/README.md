# Benchmarking and performance measurement

This directory contains scripts to **compare FFT implementations by wall-clock time**
([`bench.sh`](bench.sh)) and, on Linux, to inspect **CPU hardware counters**
([`perf_stat.sh`](perf_stat.sh) + [`C++/perf_driver.cpp`](../C++/perf_driver.cpp)).

| Tool | Answers |
|------|---------|
| **`bench.sh`** | Which variant is fastest? How does it scale with `m` and thread count? |
| **`perf_stat.sh`** | Why is it faster/slower? Did cache behaviour improve (miss rate, IPC)? |

Use both: Benchmark for speed decisions; `perf` when you care about memory behaviour
(e.g. higher-radix Cantor FFT vs radix-2).

---

## Linux prerequisites

Install **build dependencies** once. Package names differ by distribution; adjust as needed.

### Debian / Ubuntu

```bash
sudo apt update
sudo apt install -y \
  build-essential cmake pkg-config git python3 \
  libgmp-dev libboost-all-dev libssl-dev libsodium-dev \
  libomp-dev
# Hardware counters (optional, for perf_stat.sh):
sudo apt install -y linux-perf   # or: linux-tools-generic / linux-tools-$(uname -r)
```

### Fedora / RHEL / Rocky

```bash
sudo dnf install -y \
  gcc-c++ cmake pkgconfig git python3 \
  gmp-devel boost-devel openssl-devel libsodium-devel \
  libomp-devel
sudo dnf install -y perf
```

### Arch Linux

```bash
sudo pacman -S --needed \
  base-devel cmake pkgconf git python \
  gmp boost openssl libsodium openmp
sudo pacman -S --needed perf
```

### What each dependency is for

| Dependency | Used by |
|------------|---------|
| **g++**, **cmake**, **make** | Building `C++/` (libff, Google Benchmark, executables) |
| **pkg-config**, **libsodium** | Root `C++/CMakeLists.txt` (`pkg_check_modules(LIBSODIUM)`) |
| **gmp**, **boost**, **openssl** | Vendored **libff** (`C++/depends/libff`) |
| **OpenMP** (`libomp`) | Parallel Cantor/LCH FFT (`*_parallel` benchmarks); optional but recommended |
| **python3** | `summarize_benchmarks.py` (markdown tables after `bench.sh`) |
| **perf** | `perf_stat.sh` only (Linux PMU counters) |
| **git submodules** | `libff`, `benchmark` under `C++/depends/` |

**Not required on the host:** separate installs of Google Benchmark or libff — they are built from `C++/depends/` via CMake.

### Clone with submodules

```bash
git clone --recurse-submodules <repo-url>
cd additive-fft
# If you already cloned without submodules:
git submodule update --init --recursive
```

### `perf` permissions (for `perf_stat.sh` only)

If `perf stat` fails with permission errors, check:

```bash
cat /proc/sys/kernel/perf_event_paranoid
```

Values **≤ 2** usually allow normal users to read basic counters. To allow temporarily
(root required):

```bash
sudo sysctl -w kernel.perf_event_paranoid=2
```

Alternatively run the script under a user in group `perf` (distro-dependent), or use
`sudo` only for the `perf stat` line (not ideal for daily use).

Works on **Intel and AMD** Linux. Ignore Intel hybrid-only event names (`cpu_core/…`,
`cpu_atom/…`) unless you are on a P/E-core CPU and know you need them.

---

## One-time build

From the **repository root**:

```bash
cmake -S C++ -B C++/build -DCMAKE_BUILD_TYPE=Release
cmake --build C++/build --target run_benchmark perf_driver -j"$(nproc)"
```

Binaries:

| Target | Path | Role |
|--------|------|------|
| `run_benchmark` | `C++/build/run_benchmark` | Google Benchmark suite |
| `perf_driver` | `C++/build/perf_driver` | Fixed-iteration harness for `perf stat` |
| `main` | `C++/build/main` | Small demo / tests (optional) |

Release builds use `-O3 -march=native`. Rebuild on a new machine so native tuning matches
that CPU.

`run_benchmark` **requires** these environment variables (set automatically by `bench.sh`):

- `BM_MIN_RANGE`, `BM_MAX_RANGE`, `BM_STEP` — problem sizes: `m` runs from min to max
  in steps of `STEP` (each benchmark uses `2^m` field elements).

---

## Configuration — `scripts/measure.conf`

Optional. Copy and edit:

```bash
cp scripts/measure.conf.example scripts/measure.conf
./scripts/run_all.sh
```

**Precedence:** script defaults → `measure.conf` → environment variables.

Use **`MIN_RANGE` / `MAX_RANGE` / `STEP`** for both bench and perf. Only set `BM_*` or `PERF_*` if you need them to differ.

| Variable | Default | `QUICK=1` |
|----------|---------|-----------|
| `MIN_RANGE` / `MAX_RANGE` / `STEP` | 14 / 22 / 2 | same |
| `THREAD_LIST` | `1 2 4 8` | `1 8` |
| `QUICK` | `0` | — |
| `BUILD` | `C++/build` | same |
| `BENCH_REPETITIONS` | 15 | 7 |
| `BENCH_MIN_WARMUP_TIME` | 1.0 s | 0 |
| `BENCH_WARMUP_ROUNDS` | 1 | 0 |
| `PERF_ITERS` | 25 | 10 |
| `PERF_REPETITIONS` | 3 | 1 |
| `PERF_WARMUP_ROUNDS` | 1 | 0 |
| `PERF_SUITE` | `cantor` | same |
| `SKIP_PERF` / `SKIP_BENCH` | 0 | same |

Example smoke config:

```bash
QUICK=1
MIN_RANGE=14
MAX_RANGE=18
STEP=2
SKIP_PERF=1
```

---

## Full measurement run — `run_all.sh`

Run **wall-clock benchmarks and perf counters** in one command and leave:

```bash
./scripts/run_all.sh
```

Faster smoke run — either edit `scripts/measure.conf` (see above) or use env vars:

```bash
QUICK=1 MIN_RANGE=14 MAX_RANGE=18 STEP=2 THREAD_LIST="1 $(nproc)" ./scripts/run_all.sh
```

Bench only (no `perf` required):

```bash
SKIP_PERF=1 ./scripts/run_all.sh
```

Output layout (`measure_results/run_<timestamp>/`):

| Path | Contents |
|------|----------|
| `RUN_META.txt` | Shared env, CPU, bench/perf status |
| `SUMMARY.md` | Index with links to `bench/SUMMARY.md` and `perf/SUMMARY.md` |
| `console.log` | Full log from both tools |
| `bench/` | Same as `benchmark_results/run_*` (`threads_*.json`, `SUMMARY.md`, …) |
| `perf/` | Same as `perf_results/run_*` (`*.perf.txt`, `SUMMARY.md`, …) |

Shared env vars: `MIN_RANGE`, `MAX_RANGE`, `STEP` (propagate to both `BM_*` and `PERF_*`),
plus all variables documented for `bench.sh` and `perf_stat.sh`.

Both sub-tools write **`SUMMARY.md`** in the same format: `# Summary — <tool>`,
`## Run context`, then metric tables with **P=1, P=2, …** columns.

---

## Remote server (copy-paste)

All paths are **relative to the repo root**; no hardcoded host paths in the scripts.

```bash
# 1) Clone (with submodules) on the remote machine
git clone --recurse-submodules <repo-url>
cd additive-fft

# 2) One-time build (re-run after git pull or when moving to a different CPU)
cmake -S C++ -B C++/build -DCMAKE_BUILD_TYPE=Release
cmake --build C++/build --target run_benchmark perf_driver -j"$(nproc)"

# 3) Optional: copy and edit config, then run everything
cp scripts/measure.conf.example scripts/measure.conf
# set QUICK=1, MIN_RANGE, MAX_RANGE, THREAD_LIST, SKIP_PERF, etc.
./scripts/run_all.sh

# 4) Fetch results (from your laptop)
scp -r user@server:additive-fft/measure_results/run_* ./measure_results/
```

**Thread count:** set `THREAD_LIST` to match the machine, e.g. `THREAD_LIST="1 4 8 16"`.
Default is `1 2 4 8`, which is fine on most servers.

**Problem sizes:** narrow while testing, e.g. `BM_MIN_RANGE=14 BM_MAX_RANGE=18 BM_STEP=2`.
Full default (`14..22`) with four thread counts can take **hours**.

**IFFT included:** default `BENCH_FILTER` runs Cantor **affine-chart and precmp-basis**
FFT/IFFT (serial + parallel + radix-`2^K`), and LCH precmp-basis FFT/IFFT.
Override with `BENCH_FILTER='...'` if you want a subset.

**Hardware counters (optional):** only on Linux with `perf` installed and sufficient
permissions — see [`perf_stat.sh`](#hardware-counters--perf_statsh).

**Artifacts per run** (`measure_results/run_<timestamp>/`, or standalone `benchmark_results/` / `perf_results/`):

| File | Role |
|------|------|
| `RUN_META.txt` | Env, CPU model, filter — archive this with results |
| `SUMMARY.md` | Human-readable tables (combined at top level when using `run_all.sh`) |
| `console.log` | Text log (errors + progress) |
| `bench/threads_<P>.json` | Raw Google Benchmark output |
| `perf/*.perf.txt` | Raw `perf stat` output per variant |

---

## Google Benchmark — `bench.sh`

### Quick start

```bash
./scripts/bench.sh
```

Outputs land in `benchmark_results/run_<timestamp>/`:

| File | Contents |
|------|----------|
| `RUN_META.txt` | Git repo path, env, CPU model |
| `threads_<P>.json` | Google Benchmark JSON (one file per `OMP_NUM_THREADS`) |
| `console.log` | Full text log |
| `SUMMARY.md` | Markdown tables (`summarize_benchmarks.py`) |

### Faster smoke run

```bash
QUICK=1 ./scripts/bench.sh
```

Uses fewer repetitions and skips extra warmup (see `bench.sh` header).

### What is measured

- **Wall time** and **CPU time** per benchmark (median over repetitions).
- Default **filter** covers Cantor and LCH **FFT and IFFT** on the precmp-basis /
  radix-`2^K` paths (serial + OpenMP parallel). See `DEFAULT_FILTER` in `bench.sh`.
- Benchmarks are registered in [`C++/benchmark.cpp`](../C++/benchmark.cpp).
  Parameter `m` in the name is the log₂ domain size (`2^m` coefficients).

### Useful environment variables

| Variable | Default | Meaning |
|----------|---------|---------|
| `BUILD` | `C++/build` | CMake build directory |
| `OUT` | `benchmark_results/run_<ts>` | Output directory |
| `THREAD_LIST` | `1 2 4 8` | Values for `OMP_NUM_THREADS` |
| `BM_MIN_RANGE` | `14` | Minimum `m` |
| `BM_MAX_RANGE` | `22` | Maximum `m` |
| `BM_STEP` | `2` | Step for `m` |
| `BENCH_FILTER` | (see script) | Regex passed to `--benchmark_filter` |
| `QUICK` | `0` | `1` = fewer reps, no extra warmup, faster sweep |
| `BENCH_REPETITIONS` | `15` (`7` if `QUICK=1`) | Statistical repetitions |
| `OMP_PLACES` | `cores` | OpenMP affinity |
| `OMP_PROC_BIND` | `close` | OpenMP thread binding |

### Run a subset manually

```bash
export BM_MIN_RANGE=18 BM_MAX_RANGE=22 BM_STEP=2
export OMP_NUM_THREADS=4 OMP_PLACES=cores OMP_PROC_BIND=close
C++/build/run_benchmark \
  --benchmark_filter='BM_cantor_additive_fft_radix2k.*/4' \
  --benchmark_repetitions=10
```

Other families (libiop, Gao, legacy Cantor `precmp`) exist in `benchmark.cpp` but are
commented out in `REGISTER_BENCH`; uncomment to enable.

### Interpreting results

- Compare **median real time** at the same `m` and `OMP_NUM_THREADS`.
- **Parallel** variants need `OMP_NUM_THREADS>1`; serial baselines are usually run at `1`.
- Large `m` (e.g. 22) can take minutes per configuration; use `QUICK=1` or narrow
  `BM_MAX_RANGE` while developing.

---

## Hardware counters — `perf_stat.sh`

Use when Benchmark shows a time difference but you need **evidence about caches or IPC**
(e.g. “did radix-4 reduce cache misses even though wall time mixed?”).

### Full sweep (like `bench.sh`)

Run one command and leave; you get all default **Cantor radix variants** (r2, r2k2–r2k4,
serial + parallel), all **`m`** in `[PERF_MIN_RANGE .. PERF_MAX_RANGE]`, and all
**`THREAD_LIST`** values, with warmup and aggregated `SUMMARY.md` tables.

```bash
cmake --build C++/build --target perf_driver
./scripts/perf_stat.sh
```

Faster smoke sweep (fewer `m`, threads `1 8`, no warmup, 10 FFT iters, 1 perf repeat):

```bash
QUICK=1 ./scripts/perf_stat.sh
```

| | **`bench.sh`** | **`perf_stat.sh`** |
|---|----------------|---------------------|
| One command | `./scripts/bench.sh` | `./scripts/perf_stat.sh` |
| `m` sweep | `BM_MIN_RANGE` … `BM_MAX_RANGE` step `BM_STEP` | `PERF_MIN_RANGE` … `PERF_MAX_RANGE` step `PERF_STEP` |
| Radix / variants | Google Benchmark filter (FFT/IFFT names) | Default **Cantor** suite: `cantor_r2`, `r2k2`–`r2k4`, `_par` |
| Threads | `THREAD_LIST` default `1 2 4 8` | same |
| Warmup | `BENCH_WARMUP_ROUNDS` + Benchmark min warmup | `PERF_WARMUP_ROUNDS` dry `perf_driver` runs |
| Repetitions | `BENCH_REPETITIONS` (timing stats) | `PERF_REPETITIONS` (`perf stat -r`) |
| Work per measurement | Benchmark adapts iterations | `PERF_ITERS` FFT calls (default **25**) |
| Report | `SUMMARY.md` via `summarize_benchmarks.py` | `SUMMARY.md` via `summarize_perf.py` |

**Note:** A full perf sweep is **much slower** than Benchmark (PMU + fixed iters × every
variant × `m` × thread count). Use `QUICK=1` or narrow ranges while iterating.

LCH instead of Cantor: `PERF_SUITE=lch ./scripts/perf_stat.sh`  
Both: `PERF_SUITE=all ./scripts/perf_stat.sh`

### Single configuration

```bash
./scripts/perf_stat.sh cantor_r2k4_par 22 25
```

Arguments: `<variant> <m> <iters>` — one variant, one `m`, fixed FFT iteration count.

### More examples

```bash
# All Cantor radix variants, default m sweep (14..22 step 2)
./scripts/perf_stat.sh --compare-cantor

# Same but only m=22, 25 FFT iters per perf run
./scripts/perf_stat.sh --compare-cantor 22 25

# Custom variant list + m sweep
./scripts/perf_stat.sh --variants "cantor_r2 cantor_r2k4_par"

# Pin CPUs (recommended for fair comparison)
TASKSET_CPUS=0-7 THREAD_LIST="8" ./scripts/perf_stat.sh
```

### `perf_driver` variants

**LCH:** `r2`, `r2k1`…`r2k5`, `r2_par`, `r2k1_par`…`r2k5_par`  
**Cantor:** `cantor_r2`, `cantor_r2_par`, `cantor_r2k2|3|4`, `cantor_r2k2_par|3_par|4_par`

Field is fixed to **`gf256`**.

### Output layout

Results under `perf_results/run_<timestamp>/`:

| File | Contents |
|------|----------|
| `RUN_META.txt` | Events, variants, CPU, perf version |
| `<variant>_m<m>_t<threads>.perf.txt` | Full `perf stat` output + `[perf_driver]` timing |
| `console.log` | All runs concatenated |
| `SUMMARY.md` | Extracted counter lines per run |

### Environment variables

| Variable | Default | Meaning |
|----------|---------|---------|
| `PERF_MIN_RANGE` | `14` (`20` if `QUICK=1`) | Minimum `m` in full sweep |
| `PERF_MAX_RANGE` | `22` | Maximum `m` |
| `PERF_STEP` | `2` | Step for `m` |
| `PERF_ITERS` | `25` (`10` if `QUICK=1`) | FFT calls per `perf stat` run |
| `PERF_REPETITIONS` | `3` (`1` if `QUICK=1`) | Outer `perf stat -r` repeats |
| `PERF_WARMUP_ROUNDS` | `1` (`0` if `QUICK=1`) | Dry runs before each measurement |
| `PERF_WARMUP_ITERS` | `5` (`2` if `QUICK=1`) | FFT iters per warmup dry run |
| `PERF_SUITE` | `cantor` | `cantor`, `lch`, or `all` (no-arg full sweep) |
| `PERF_EVENTS` | `cycles,instructions,cache-references,cache-misses` | `perf stat -e` list |
| `THREAD_LIST` | `1 2 4 8` (`1 8` if `QUICK=1`) | `OMP_NUM_THREADS` sweep |
| `TASKSET_CPUS` | (unset) | If set, `taskset -c …` |
| `OUT`, `BUILD` | (see script) | Same idea as `bench.sh` |

**Cache miss fraction** (informal):  
`cache-misses / cache-references` from the same run — compare variants at fixed `m`,
not across machines.

For last-level cache detail (when supported):

```bash
PERF_EVENTS="cycles,instructions,cache-references,cache-misses,LLC-loads,LLC-load-misses" \
  ./scripts/perf_stat.sh cantor_r2k4_par 22 25
```

Check supported events: `perf list`.

---

## Troubleshooting

| Symptom | Likely fix |
|---------|------------|
| `run_benchmark` not found | `cmake --build C++/build --target run_benchmark` |
| `libsodium not found` | Install `libsodium-dev` / `libsodium` |
| `gmp.h` missing | Install `libgmp-dev` / `gmp-devel` |
| OpenMP parallel same as serial | Install `libomp-dev`; rebuild |
| `BM_MIN_RANGE` / getenv crash | Always set via `bench.sh` or export the three `BM_*` vars |
| `perf: command not found` | Install `linux-perf` / `perf` package |
| Zero or empty perf counters | Lower `perf_event_paranoid` or fix permissions |
| Submodule / missing `depends/` | `git submodule update --init --recursive` |
| Very long run | `QUICK=1`, smaller `BM_MAX_RANGE`, or narrower `BENCH_FILTER` |

---

## File index

| File | Purpose |
|------|---------|
| [`load_config.sh`](load_config.sh) | Shared config loader (`measure.conf`) |
| [`measure.conf.example`](measure.conf.example) | Template config (copy to `measure.conf`) |
| [`run_all.sh`](run_all.sh) | Run `bench.sh` + `perf_stat.sh` into one `measure_results/run_*` tree |
| [`bench.sh`](bench.sh) | Orchestrate Google Benchmark + JSON + summary |
| [`perf_stat.sh`](perf_stat.sh) | Wrap `perf_driver` with `perf stat` |
| [`summarize_benchmarks.py`](summarize_benchmarks.py) | Build `SUMMARY.md` from `threads_*.json` |
| [`summarize_perf.py`](summarize_perf.py) | Build `SUMMARY.md` tables from `*.perf.txt` |
| [`summarize_run.py`](summarize_run.py) | Combined top-level `SUMMARY.md` for `run_all.sh` |

**Note:** `perf_driver` and `perf_stat.sh` measure **FFT only** (no IFFT variants yet).
Use `bench.sh` for FFT vs IFFT wall-clock comparison.

| [`../C++/benchmark.cpp`](../C++/benchmark.cpp) | Benchmark registrations |
| [`../C++/perf_driver.cpp`](../C++/perf_driver.cpp) | Minimal FFT loop for PMU measurement |
