# Parallelizing the Radix-$2^K$ Cantor AFFT — Design Notes

This document is a design write-up *before* implementation. It captures where
parallelism lives in the algorithm, why earlier attempts at parallelizing the
radix-2 version were disappointing, what changes once we move to radix-$2^K$
(`additive_FFT_radix2k<K>` in `C++/Cantor/fft.tcc`), and the staged plan for
adding a parallel variant.

It is the parallelization companion to
[`RADIX4_CANTOR_FFT.md`](RADIX4_CANTOR_FFT.md) and
[`RADIX2_CANTOR_FFT_PASSES.md`](RADIX2_CANTOR_FFT_PASSES.md).

---

## 1. TL;DR

* **Why the original radix-2 parallelization underwhelmed**: that kernel does
  $\sim 1$ field op per byte streamed, so a single core already saturates DRAM
  bandwidth. Adding threads multiplied contention without adding compute
  headroom.
* **What changes with radix-$2^K$**: arithmetic intensity rises by a factor
  $\frac{2K}{K+1}$ ($1.33\times$ for $K{=}2$, $1.6\times$ for $K{=}4$), there
  are $K\times$ fewer barriers, and per-task work is $K\times$ larger. All three
  reasons the radix-2 version scaled poorly are improved.
* **Recommended plan**: build the parallel variant in two commits.
  - **Commit 1 — Strategy A (inter-module).** `#pragma omp parallel for` over
    the `module` loop, threshold-gated. Simple, ~50 lines of change. Expected
    speedup on 8 P-cores at $m{=}22$, $K{=}4$: **~3–4×** (Amdahl-bounded by
    super-round 0).
  - **Commit 2 — Strategy B (hybrid intra/inter-module).** Add intra-module
    parallelism for the rounds where `n_modules < num_threads` (super-rounds
    0/1). Lifts the Amdahl ceiling to **~6×** on 8 cores.
* **Tunables**: a small `parallel_config` struct (`num_threads`,
  `parallel_threshold_log2_elems`, `schedule`), defaulting to
  `parallel_config::from_env()` which reads `CANTOR_FFT_*` environment
  variables, with hard-coded fallbacks.

---

## 2. Where parallelism lives in `additive_FFT_radix2k<K>`

Reproducing the per-super-round shape for reference (see
`C++/Cantor/fft.tcc:421–536`):

```cpp
for (size_t round_j = 0; round_j < num_rounds; ++round_j) {
    const size_t L     = 1ull << (m - K * round_j);   // module size
    const size_t chunk = L >> K;
    const size_t half  = L >> 1;

    size_t offset = 0;
    for (size_t module = 0; module < n_modules; ++module) {
        // 1. Build base_T[1..K] from module index           (O(K * log n_modules))

        // 2. Sweep 1: level-1 division
        //    descending k in [half, L); writes to g[offset + k - nz]
        //    HAS read-after-write recurrence — not splittable per k.

        // 3. Sweeps 2..K: level-ell divisions
        //    for ell = 2..K: 2^{ell-1} INDEPENDENT bands, each does a
        //    descending k-loop inside its band.

        // 4. Sweep K+1: in-register K-stage butterfly closer
        //    chunk = L/2^K INDEPENDENT tuples (each an array of Q values).

        offset += L;
    }

    n_modules <<= K;
}
```

There are three structural axes of parallelism, with very different properties:

| Axis | Granularity | Independent? | Where it's wide |
|------|-------------|--------------|-----------------|
| **(A) outer modules within a super-round** | $n_\text{mod}$ pieces of size $L$ | yes (offsets disjoint) | wide for late rounds, narrow for early rounds |
| **(B) bands inside sweeps 2..K** | $2^{\ell-1}$ pieces in sweep $\ell$ | yes | wide *within* a single module — even for super-round 0 |
| **(C) tuples inside the closer (sweep K+1)** | $\text{chunk} = L/Q$ tuples | yes | wide for early rounds (huge $L$), narrow for late rounds |

Sweep 1 has a chain dependency on $g[k - \text{nz}]$ and is **not** parallelizable
per-$k$ without changing the algorithm. That's our irreducible serial fragment.

Importantly: *axes (A) and (C) are complementary*. (A) is wide exactly when (C)
is narrow, and vice versa. (B) is moderate everywhere. A good design exploits
whichever is wide *at this round*.

---

## 3. Round indexing, in this codebase

A reminder, because the convention matters:

| `round_j` | `n_modules` | $L$ (module size) |
|---:|---:|---:|
| 0           | $1$           | $2^m$              |
| 1           | $Q = 2^K$     | $2^{m-K}$          |
| 2           | $Q^2$         | $2^{m-2K}$         |
| $\dots$     | $\dots$       | $\dots$            |
| $m/K - 1$   | $2^{m-K}$     | $2^K$              |

Higher round ⇒ **smaller** modules, **more** of them. The total work per
super-round is approximately constant ($\sim 2^m$ field ops), so no single
super-round dominates the runtime — the relevant differences across rounds are
**granularity** and **memory locality**, not work.

For $m=22$, $K=4$:

| `round_j` | `n_modules` | $L$ | parallel posture |
|---:|---:|---:|---|
| 0  | $1$              | $2^{22}=4{\rm M}$  | only axis (B/C); axis (A) is dead |
| 1  | $16$             | $2^{18}$            | axis (A) at 16-way (fine for $\le 16$ threads) |
| 2  | $256$            | $2^{14}$            | axis (A) saturates any thread count |
| 3  | $2^{12}$         | $2^{10}$            | axis (A) saturates; modules near L1-resident |
| 4  | $2^{16}$         | $64$ elems          | axis (A) saturates but modules are tiny — keep serial |
| residual radix-2 | $2^{20}$+ | $\le 4$ | tiny; serial |

So the parallelization story per round splits into three regimes:

* **Wide rounds** (rounds 2..3): axis (A), `parallel for` over modules, static
  schedule. This is the bulk of the work.
* **Narrow rounds** (rounds 0..1): need axes (B) and (C) to fill threads.
* **Tiny rounds** (round 4 and the radix-2 residual): serial. Threading
  overhead dominates the work.

---

## 4. Why earlier radix-2 parallelization disappointed

For completeness — your original observation:

A radix-2 round does **one byte-level XOR/multiply per element**, then writes
back. Arithmetic intensity is roughly **1 op per byte**, well below the DRAM
ridge point of any modern CPU (typically 5–10 FLOP/byte). One core already
saturates the memory subsystem, so adding threads only multiplies DRAM
contention. On top of that:

* Round 0 has 1 module → axis (A) gives no parallelism.
* The last round has $2^{m-1}$ pairs → modules are below the cache-line size
  → false sharing.
* Radix-2 has $2m$ barriers ($\sim 44$ for $m{=}22$) → on small rounds the
  barrier cost dominates the work.

Diagnosis "memory-bound" was right at the high-$m$ end and "granularity-bound"
at the low-$m$ end. The radix-$2^K$ rewrite materially improves all three of
these.

---

## 5. The three viable strategies

### Strategy A — pure inter-module (`omp for` over modules)

```cpp
#pragma omp parallel
{
    for (size_t round_j = 0; round_j < num_rounds; ++round_j) {
        // ... per-round setup ...

        if (n_modules >= 2 && L >= L_threshold) {
            #pragma omp for schedule(static)
            for (size_t module = 0; module < n_modules; ++module)
                process_module(module);
        } else {
            #pragma omp single
            for (size_t module = 0; module < n_modules; ++module)
                process_module(module);
        }
        // implicit barrier
    }
}
```

| ✓ | ✗ |
|---|---|
| Smallest diff vs current code (~50 lines) | Super-round 0 is serial → Amdahl floor at $K/m$ |
| Each thread keeps the whole per-module body in i-cache | Super-round 1 has only $Q$ modules — under-fills high thread counts |
| Already cache-friendly because per-module work is fused | |

For $K{=}4, m{=}22$, Amdahl ceiling on $P$ cores:
$\text{speedup} \le 1 / (K/m + (1 - K/m)/P) = 1 / (0.18 + 0.82/P)$.
At $P{=}8$ this is **~3.6×**. At $P{=}16$ it caps at **~4.7×**.

### Strategy B — hybrid (inter-module when wide, intra-module when narrow)

Per super-round decision:

| condition | parallelize over |
|---|---|
| $n_\text{mod} \ge P$ and $L \ge L_\text{threshold}$ | modules — axis (A) |
| $n_\text{mod} < P$ and $L \ge L_\text{threshold}$ | inside the module: axis (B) for sweeps 2..K, axis (C) for sweep K+1; sweep 1 stays serial within the module |
| $L < L_\text{threshold}$ | serial |

The intra-module path looks roughly:

```cpp
for (size_t module = 0; module < n_modules; ++module) {
    // sweep 1 (recurrence) — serial, executed by one thread per module
    #pragma omp single
    sweep1(module);

    // sweeps 2..K — bands are independent
    for (size_t ell = 2; ell <= K; ++ell) {
        #pragma omp for schedule(static)
        for (size_t b = 0; b < (1ull << (ell - 1)); ++b)
            do_band(module, ell, b);
    }

    // sweep K+1 — tuples are independent
    #pragma omp for schedule(static)
    for (size_t jj = 0; jj < chunk; ++jj)
        closer_tuple(module, jj);
}
```

For super-round 0, sweep 1 is the only serial fragment per module. Its work is
$\approx \frac{1}{K+1}$ of the round's work, so the serial fraction of the
*whole FFT* drops to roughly:

$$\sigma \;=\; \frac{K/m}{K+1} \;\approx\; 0.036 \quad \text{for } K=4, m=22.$$

Amdahl ceiling on 8 cores becomes $\sim 6.3\times$; on 16 cores $\sim 9\times$.

| ✓ | ✗ |
|---|---|
| Recovers most of super-round 0/1 | Two parallel code paths — more complex |
| Approaches linear scaling for $m \gg K$ | Need careful collapse-loop semantics; risk of false sharing if axis (B) bands are < cache line — check at threshold-decision time |

### Strategy C — task graph

Decompose the whole FFT into independent typed tasks (band, tuple, full module)
and dispatch via `#pragma omp taskloop` or a custom work-stealing pool.
Theoretically cleanest, but for our workload the gains over Strategy B are
marginal and the implementation cost is high. **Not recommended for this
codebase.**

---

## 6. Recommendation: Strategy A first, Strategy B second

The phased plan:

### Commit 1 — Strategy A

1. **Build system**

   * In `C++/CMakeLists.txt`: `find_package(OpenMP REQUIRED)`.
   * Link `OpenMP::OpenMP_CXX` to the targets that consume `cantor::fft`
     (the test binary and the benchmark binary).
   * Pass `-fopenmp` consistently (Clang/GCC).

2. **Lift per-module body into a helper** in `fft.tcc`:

   ```cpp
   template<size_t K, typename FieldT>
   static inline void radix2k_process_module(
       std::vector<FieldT>& g,
       size_t offset,
       size_t L,
       size_t chunk,
       size_t half,
       const std::array<std::vector<size_t>, K + 1>& nz_S,
       const std::array<FieldT, K + 1>& base_T,
       const std::array<FieldT, /*Q_half*/>& block_offsets);
   ```

   The `module` loop body in the existing `additive_FFT_radix2k<K>` is an
   almost-direct copy of this; we just have to pre-compute `base_T` per-module
   outside the helper.

3. **Add `additive_FFT_radix2k_parallel<K>`** in `fft.hpp` / `fft.tcc`:

   ```cpp
   template<size_t K, typename FieldT>
   std::vector<FieldT> additive_FFT_radix2k_parallel(
       const std::vector<FieldT>& poly_coeffs,
       const libiop::affine_subspace<FieldT>& domain,
       parallel_config cfg = parallel_config::from_env());
   ```

   Body: identical to the serial version, except the main super-round loop and
   the residual radix-2 loop are wrapped in **one** `#pragma omp parallel`
   region with `omp for schedule(static)` on the `module` dimension, gated by
   `if (n_modules >= 2 && L >= cfg.parallel_threshold_elems)`.

4. **`parallel_config`** in a new small header `C++/Cantor/parallel_config.hpp`:

   ```cpp
   namespace cantor {
   struct parallel_config {
       size_t num_threads = 0;                 // 0 = auto (omp_get_max_threads)
       size_t parallel_threshold_log2_elems = 0; // 0 = auto (default 10 → 1024)
       enum class Schedule { Static, Dynamic, Guided } schedule = Schedule::Static;

       static parallel_config detect();        // pure auto (compiled-in defaults)
       static parallel_config from_env();      // detect() + CANTOR_FFT_* env vars
   };
   }
   ```

   Env vars:
   * `CANTOR_FFT_THREADS` — overrides `num_threads`.
   * `CANTOR_FFT_THRESHOLD_LOG2` — overrides
     `parallel_threshold_log2_elems`.
   * `CANTOR_FFT_SCHEDULE` — `"static"|"dynamic"|"guided"`.

   No L2 detection in commit 1. A hard-coded threshold of $2^{10}$ elements is
   above the false-sharing and overhead floors for `gf256`.

5. **Correctness test** in `C++/main.cpp`: extend `test_radix2k_correctness`
   to compare `additive_FFT_radix2k_parallel<K>` against `additive_FFT` for
   $K\in\{2,3,4\}$, $m\in\{4,\dots,16\}$, with several random thread counts.

6. **Benchmarks** in `C++/benchmark.cpp`:

   ```cpp
   template<size_t K> static void BM_cantor_additive_fft_radix2k_parallel(...);
   BENCHMARK_TEMPLATE(BM_cantor_additive_fft_radix2k_parallel, 2)
       ->DenseRange(MIN_RANGE, MAX_RANGE, STEP)->Unit(microsecond);
   // ... 3, 4 ...
   ```

   And a runbook for scaling sweeps:

   ```bash
   for P in 1 2 4 8 16; do
     OMP_NUM_THREADS=$P OMP_PLACES=cores OMP_PROC_BIND=close \
     taskset -c 0-$((P-1)) ./run_benchmark \
       --benchmark_filter='radix2k_parallel<4>/22' \
       --benchmark_repetitions=5
   done
   ```

### Commit 2 — Strategy B (only if Commit 1 hits an Amdahl wall)

Add the intra-module branch sketched in §5. The key implementation choices:

* Keep the **same** `#pragma omp parallel` region wrapping the whole FFT;
  switch internally between the two `omp for` shapes per super-round.
* For axis (B) parallelism, pre-compute per-band `(band_start, band_end)`
  pairs into a small `std::array` so the parallel loop body is a flat
  `for (b ...)` with no inner control flow.
* For axis (C), the closer's per-tuple work is already a tight inner loop —
  parallelize the outer `jj` loop.
* For super-round 0, sweep 1 stays serial. We can issue it from the master
  thread inside `#pragma omp single nowait` and let the others race ahead to
  the band-parallel sweeps 2..K of *the same module* (no, they can't — sweep
  2 reads from sweep 1's outputs). So just `#pragma omp single` with the
  default barrier.

No new public API: just the existing
`additive_FFT_radix2k_parallel<K>(poly, domain, cfg)` gets smarter about
switching between the two paths.

### What we defer

* **L2 / cache-size auto-detection.** Linux sysfs is portable enough but
  adds plumbing; `hwloc` is the right answer if we ever do auto-tune. Not
  needed for first cut.
* **Hybrid CPU (P/E core) handling.** Document `OMP_PLACES=cores
  OMP_PROC_BIND=close taskset -c <P-cores>` as the runbook; don't bake it in.
* **NUMA-aware first-touch allocation.** Only relevant on multi-socket boxes.
* **FFTW-style autotuning / wisdom file.** Only if heuristics misfire on a
  real machine.

---

## 7. Expected speedups (rough model)

Using the simple Amdahl model with serial fraction $\sigma$:

| version | $\sigma$ ($K{=}4, m{=}22$) | $P{=}4$ | $P{=}8$ | $P{=}16$ |
|---|---:|---:|---:|---:|
| Strategy A | $K/m \approx 0.18$       | $2.7\times$ | $3.6\times$ | $4.4\times$ |
| Strategy B | $K/(m(K{+}1)) \approx 0.036$ | $3.6\times$ | $6.3\times$ | $9.4\times$ |

These are upper bounds — they ignore DRAM bandwidth contention, barrier cost,
and per-module setup. Realistic speedups will be ~70–80% of the Amdahl
ceiling on a typical 8-core desktop, i.e. ~$3\times$ for Strategy A and
~$5\times$ for Strategy B. The Strategy A → Strategy B delta of $\sim 2\times$
is the headline reason to do Commit 2.

For larger $m$ (e.g., $m=24$, $K=4$), $\sigma$ shrinks, both ceilings rise,
and the parallel version gets relatively better.

---

## 8. Risks / open questions

1. **DRAM bandwidth saturation at large $m$.** Even with arithmetic intensity
   $1.6\times$ the radix-2 baseline, the closer sweep at super-round 0
   ($L=2^m$) is still streaming the full array. On a single-channel desktop
   8 threads will share that channel. Mitigation: don't expect linear scaling
   at $m=24$+; benchmark to confirm we still get a useful win.

2. **`#pragma omp for` overhead on small sweeps.** Sweeps 2..K of super-round
   0 dispatch only a few independent bands ($2^1, 2^2, \dots, 2^{K-1}$). For
   $K{=}4$ that's 2, 4, 8 bands — sometimes fewer than threads. Use
   `schedule(static, 1)` with `nowait` where safe to minimize idle barriers.

3. **`#pragma omp single` granularity for sweep 1 in narrow rounds.** Whichever
   thread executes sweep 1 holds back the rest. We accept this; it's the
   irreducible serial fragment.

4. **False sharing in axis (C).** The closer's per-tuple writes are at strided
   positions `i * chunk + jj`. If two threads land on adjacent `jj` values and
   `chunk * sizeof(FieldT) <` cache-line size, they ping-pong. Threshold the
   axis (C) parallel path on `chunk * sizeof(FieldT) >= 2 * cache_line`.

5. **Determinism.** Field addition is associative and commutative, so the
   parallel result equals the serial result bitwise; no determinism caveat
   unless we later introduce `omp reduction` (which we don't here).

---

## 9. API and migration

After Commit 1 the public surface looks like:

```cpp
// existing
namespace cantor {
template<typename FieldT>
std::vector<FieldT> additive_FFT(const std::vector<FieldT>&,
                                 const libiop::affine_subspace<FieldT>&);

template<typename FieldT>
std::vector<FieldT> additive_FFT_radix4(const std::vector<FieldT>&,
                                         const libiop::affine_subspace<FieldT>&);

template<size_t K, typename FieldT>
std::vector<FieldT> additive_FFT_radix2k(const std::vector<FieldT>&,
                                          const libiop::affine_subspace<FieldT>&);

// new
struct parallel_config { /* see §6 */ };

template<size_t K, typename FieldT>
std::vector<FieldT> additive_FFT_radix2k_parallel(
    const std::vector<FieldT>&,
    const libiop::affine_subspace<FieldT>&,
    parallel_config = parallel_config::from_env());
}
```

No existing call site changes. Callers that want the parallel version opt
in by switching the function name. Callers that want full control pass an
explicit `parallel_config`.

---

## 10. Validation plan

* **Correctness**: `test_radix2k_correctness` extended to cover the parallel
  variant for $K\in\{2,3,4\}$, $m\in[4,16]$, thread counts in $\{1,2,4,8\}$.
  All outputs must equal the serial radix-2 baseline element-by-element.
* **Single-thread regression**: `OMP_NUM_THREADS=1` performance of
  `additive_FFT_radix2k_parallel<K>` must match `additive_FFT_radix2k<K>`
  within 5% — confirms the OpenMP plumbing has zero overhead in serial mode.
* **Scaling sweep**: $P\in\{1,2,4,8\}$ at $m\in\{18,20,22\}$, $K\in\{2,4\}$.
  Plot speedup curves; expect $\ge 2.5\times$ on 8 cores at $m=22$ for
  Strategy A.
* **`perf stat` cross-check**: confirm LLC-load-misses do **not** explode with
  threads (sign of bandwidth contention) and `cycles ÷ instructions` stays
  roughly flat.

If the scaling plot for Strategy A shows a clear plateau at $\sim 3.5\times$
on $P{=}8$, that's the signal to do Commit 2. If it scales linearly to 8
threads, we ship Strategy A and move on.
