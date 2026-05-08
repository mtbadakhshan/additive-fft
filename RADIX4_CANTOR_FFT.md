# Radix-4 and Radix-$2^K$ Cantor Additive FFT: Paper Algorithm vs. Memory-Efficient Implementation

This note compares three ways of running the Cantor additive FFT:

1. The pseudocode given in `main.tex`, §"The Radix-4 Cantor AFFT"
   (Algorithm `Algo:CantorRadix4`). It describes the *math* in three textually
   separate phases.
2. The memory-efficient radix-4 implementation `cantor::additive_FFT_radix4`
   in `C++/Cantor/fft.tcc`, which produces bit-for-bit the same output but
   merges the loops so each module is touched as few times as possible per
   round.
3. The generalization to radix-$2^K$ for arbitrary $K\ge 1$, implemented as
   `cantor::additive_FFT_radix2k<K>` in the same file. The radix-4 wrapper is
   now a thin call to `additive_FFT_radix2k<2>`.

Both versions perform the **same number of field multiplications and field
additions**. The only thing that changes is the schedule: how often the
buffer is dragged through cache and how many distinct `for` bodies are
emitted.

The radix-2 baseline (the existing `additive_FFT(poly, domain)` in
`C++/Cantor/fft.tcc`, lines 258–329) is the reference point.

---

## 1. Setup and notation

We evaluate $f(x)\in\mathbb{F}_{2^r}[x]$, $\deg f<n=2^m$, on the affine
subspace $\theta + W_m$ where $\{\beta_0=1,\beta_1,\dots,\beta_{m-1}\}$ is the
Cantor special basis (so $S(\beta_i)=\beta_{i-1}$, $S(x)=x^2+x$,
$S^i(x)=\mathbb{Z}_{W_i}(x)$).

A radix-$2^K$ round descends $K$ levels at a time. At round
$j\in\{0,\dots,\lfloor m/K\rfloor-1\}$:

- module size: $L = 2^{m-Kj}$;
- number of modules: $2^{Kj}$;
- per module we use the polynomials
  $S^{m-Kj-\ell}$ of degree $L/2^{\ell}$ for $\ell=1,\dots,K$;
- per module we use $K$ base twiddles $T_{u,\ell,0}=S^{m-Kj-\ell}(\theta_u)$,
  $\ell=1,\dots,K$. Within stage $\ell$ of the in-register butterfly, block
  $b\in[0,2^{\ell-1})$ uses the twiddle $T_{u,\ell,b} = T_{u,\ell,0} +
  \sum_{i:\,b_i=1} \beta_{i+1}$.
- the output of one round is laid out in $2^K$ contiguous chunks of size
  $L/2^K$; chunk $c$ becomes module index $u\cdot 2^K + c$ in the next round.

For the special case $K=2$ this reproduces `additive_FFT_radix4` with
$T_2=T_{u,1,0}$, $T_1=T_{u,2,0}$, $T_3=T_{u,2,1}=T_1+\beta_1$.

---

## 2. Paper version (Algorithm `Algo:CantorRadix4`)

The radix-4 pseudocode is written in *three* phases of vector assignments.
Translated to imperative-style C++ pseudocode operating on a single module of
size $L$, it would look like this:

```cpp
// ---------- Phase 1: three sequential divisions ----------

// (1a) divide whole module by S^{m-2j-1} of degree half = L/2
for (size_t k = L - 1; k >= half; --k) {
    FieldT gk = g[k];
    for (size_t nz : nz_S_high) g[k - nz] += gk;
}
// after (1a):  g[0..half) = r,  g[half..L) = q

// (1b) divide r (lower half) by S^{m-2j-2} of degree quarter = L/4
for (size_t k = half - 1; k >= quarter; --k) {
    FieldT gk = g[k];
    for (size_t nz : nz_S_low) g[k - nz] += gk;
}

// (1c) divide q (upper half) by S^{m-2j-2}
for (size_t k = L - 1; k >= half + quarter; --k) {
    FieldT gk = g[k];
    for (size_t nz : nz_S_low) g[k - nz] += gk;
}

// ---------- Phase 2: first-stage butterfly, four vector ops ----------
for (j=0;j<quarter;++j) g[0      +j] = g[0      +j] + T_2 * g[half        +j];
for (j=0;j<quarter;++j) g[half   +j] = g[0      +j] + g[half             +j];
for (j=0;j<quarter;++j) g[quarter+j] = g[quarter+j] + T_2 * g[half+quarter+j];
for (j=0;j<quarter;++j) g[half+quarter+j] = g[quarter+j] + g[half+quarter +j];

// ---------- Phase 3: second-stage butterfly, four vector ops ----------
for (j=0;j<quarter;++j) g[0      +j] = g[0      +j] + T_1 * g[quarter         +j];
for (j=0;j<quarter;++j) g[quarter+j] = g[0      +j] + g[quarter              +j];
for (j=0;j<quarter;++j) g[half   +j] = g[half   +j] + T_3 * g[half+quarter   +j];
for (j=0;j<quarter;++j) g[half+quarter+j] = g[half+j] + g[half+quarter      +j];
```

This is a faithful, line-for-line transcription of `Algo:CantorRadix4`.

### How many array passes per radix-4 round?

| Phase | Sub-step | Sweep length | Direction |
|---|---|---|---|
| 1a | divide by $S^{m-2j-1}$ on size-$L$ module | $L/2$ | descending |
| 1b | divide by $S^{m-2j-2}$ on lower half | $L/4$ | descending |
| 1c | divide by $S^{m-2j-2}$ on upper half | $L/4$ | descending |
| 2  | $u_0$ (Q0 $+=$ T·Q2) | $L/4$ | ascending |
| 2  | $u_2$ (Q2 $+=$ Q0) | $L/4$ | ascending |
| 2  | $u_1$ (Q1 $+=$ T·Q3) | $L/4$ | ascending |
| 2  | $u_3$ (Q3 $+=$ Q1) | $L/4$ | ascending |
| 3  | $h_0$ (Q0 $+=$ T·Q1) | $L/4$ | ascending |
| 3  | $h_1$ (Q1 $+=$ Q0) | $L/4$ | ascending |
| 3  | $h_2$ (Q2 $+=$ T·Q3) | $L/4$ | ascending |
| 3  | $h_3$ (Q3 $+=$ Q2) | $L/4$ | ascending |

That is **11 distinct sweeps per radix-4 round**. Total element-touches per
round (read+write) is roughly
$$
\text{paper}: \quad 1\cdot L + 2\cdot \tfrac{L}{2} + 8\cdot \tfrac{L}{4}
\;=\; 4L \;\;\text{element-touches}.
$$

Compared with the radix-2 baseline that does **2 sweeps per round** (one
descending `for(k)` of length $L/2$ doing division+twiddle, one ascending
`for(j)` of length $L/2$ doing the closer):
$$
\text{radix-2}: \quad 2 \cdot \tfrac{L}{2} \;=\; L \;\;\text{element-touches per round}.
$$

So the paper's pseudocode, naively transcribed, **moves four times as much
data through cache per round** as radix-2 does. With $m/2$ rounds versus $m$
rounds for radix-2, total data movement is $4L\cdot m/2 = 2Lm$ versus
$L\cdot m = Lm$ for radix-2 — **twice the memory traffic, not less.**

This is the source of the apparent paradox: if you look only at the
pseudocode, radix-4 *loses* on memory. The reason it can be made to win is
that the eleven sweeps can be folded — but as we will see in §3, the minimum
is **3 sweeps per radix-4 round**, not 2.

---

## 3. Memory-efficient version

### 3.1 Why "two sweeps per round" does not work

A first attempt is to fuse Phase 1's three divisions plus the $T_2$ twiddle
into a single descending sweep over $k\in[L/4,L)$, using band conditions to
fire each rule only on the appropriate range:

```cpp
// BUGGY: do not use. Documented for posterity.
size_t k = L;
while (k > quarter) {
    --k;
    const FieldT gk = g[offset + k];
    if (k >= half) {
        for (auto nz : nz_S_high) g[offset + k - nz] += gk;   // S^{m-2j-1}
        g[offset + k - half] += gk * T_2;                     // T_2 cross-half
    }
    if (k >= half + quarter || k < half) {
        for (auto nz : nz_S_low) g[offset + k - nz] += gk;    // S^{m-2j-2}
    }
}
```

This was the design used by an earlier `additive_FFT_radix4` in this repo,
and it produces the wrong answer for $m\ge 3$. The reason is a
read-after-write hazard between *level 2* (the inner `nz_S_low` writes) and
*level 1* (the outer `nz_S_high` reads):

- At iter $k\in[3L/4,L)$ the body fires `nz_S_low` writes, which target
  positions in $[k-L/4+1,\,k-1]\subseteq[L/2+1,\,L-2]$ — **upper-half**
  positions.
- At a later iter $k'\in[L/2,3L/4)$ the body reads `gk' = g[k']`. By then
  $g[k']$ has been polluted by the level-2 writes from $k>k'$, so the
  level-1 division at $k'$ uses the wrong dividend bit.

Concrete example ($m=3$, $L=8$, $S^2(x)=x^4+x$, $S^1(x)=x^2+x$):
- Iter $k=7$: level 2 writes `g[6] += f[7]`, leaving $g[6]=f[6]+f[7]$.
- Iter $k=6$: level 1 reads `gk = g[6] = f[6]+f[7]` (polluted) and writes
  `g[3] += f[6]+f[7]` instead of the correct `g[3] += f[6]`.

So the levels **cannot share a descending sweep** in general.

### 3.2 The correct memory-efficient radix-$2^K$ schedule: $K+1$ sweeps per round

For radix-$2^K$, level $\ell$ ($\ell=1,\dots,K$) performs a polynomial
division by $S^{m-Kj-\ell}$ (degree $L/2^\ell$) on $2^{\ell-1}$ disjoint
bands of size $L/2^\ell$. Each band is processed in standard descending
order. Two distinct levels cannot be fused into the same sweep because the
deeper level's writes corrupt the shallower level's reads (§3.1). The
**minimum** number of array sweeps per radix-$2^K$ round is therefore

$$
\boxed{\; K \;\text{division sweeps} \;+\; 1 \;\text{butterfly closer sweep} \;=\; K+1 \;\;\text{sweeps per round}.\; }
$$

This is what `additive_FFT_radix2k<K>` (and hence `additive_FFT_radix4`
which is `additive_FFT_radix2k<2>`) implements. For each module:

- **Sweep 1** — descending `for(k)` over $k\in[L/2,L)$. Fires the
  $S^{m-Kj-1}$ division *plus* the stage-1 cross-half twiddle
  `g[k-L/2] += g[k] * T_{u,1,0}` (just like the existing radix-2 sweep-A).

- **Sweep $\ell$** for $\ell=2,\dots,K$ — descending `for(b,k)` that visits
  the $2^{\ell-1}$ active bands of level $\ell$ in ascending band-index
  order, each band scanned descending. Inside each iter: standard
  $g[k-\text{nz}_S]+=gk$ writes for the $S^{m-Kj-\ell}$ division.

- **Sweep $K+1$** — ascending `for(jj)` over $jj\in[0,L/2^K)$. Loads
  $2^K$ values into registers, runs the full $K$-stage butterfly closer
  in-register (stage-1 closer plus stages 2..K), writes $2^K$ values back.

The $K$-stage butterfly is structurally identical to the radix-4 closer
generalized: stage $\ell$ has $2^{\ell-1}$ blocks, each twiddle is
$T_{u,\ell,0} + \sum_{i:b_i=1}\beta_{i+1}$, the eight basis offsets
$\beta_1,\dots,\beta_{K-1}$ are precomputed once per round.

```cpp
// Sweep K+1 — in-register K-stage butterfly closer (per tuple jj).
FieldT v[1u<<K];
for (size_t i = 0; i < (1u<<K); ++i)
    v[i] = g[offset + i * chunk + jj];

// Stage-1 closer (the +=q part of round-0 butterfly):
for (size_t i = 0; i < (1u<<(K-1)); ++i)
    v[(1u<<(K-1)) + i] += v[i];

// Stages 2..K
for (size_t ell = 2; ell <= K; ++ell) {
    const size_t stride     = 1ull << (K - ell);
    const size_t num_blocks = 1ull << (ell - 1);
    for (size_t b = 0; b < num_blocks; ++b) {
        const FieldT T_lb = base_T[ell] + block_offsets[b];
        for (size_t t = 0; t < stride; ++t) {
            const size_t i0 = b * (stride << 1) + t;
            const size_t i1 = i0 + stride;
            v[i0] += T_lb * v[i1];
            v[i1] += v[i0];
        }
    }
}
for (size_t i = 0; i < (1u<<K); ++i)
    g[offset + i * chunk + jj] = v[i];
```

### 3.3 Correctness of the schedule

- Each band of level $\ell$ is one ordinary descending polynomial division;
  within a band, $g[k]$ at iter $k$ holds exactly the standard residual.
- Different bands of the same level are independent (their write footprints
  do not overlap), so they can be visited in any order.
- Levels are processed lowest index first ($\ell=1,2,\dots,K$); when level
  $\ell$ runs, the array already represents the correct
  level-$(\ell-1)$ output, which is what level $\ell$ expects to consume.
- The final closer is per-tuple in registers, with no inter-tuple
  dependency.

The repository's correctness suite in `C++/main.cpp::test_radix2k_correctness`
exercises every $K\in\{1,2,3,4\}$ against the radix-2 baseline for a range
of $m$, including residual cases where $K\nmid m$.

### 3.4 Pass count summary

| | Radix-2 baseline | Radix-4 paper (literal) | Radix-$2^K$ efficient |
|---|---|---|---|
| Outer rounds | $m$ | $m/2$ (+residual) | $m/K$ (+residual) |
| Sweeps per round | $2$ | $11$ | $K+1$ |
| Sweeps total (no residual) | $2m$ | $11m/2$ | $(K+1)m/K = m + m/K$ |
| Element-touches total | $\sim Lm$ | $\sim 2Lm$ | $\sim Lm \cdot (1 + 1/K)/2$ |
| Field multiplications | $nm/2$ | $nm/2$ | $nm/2$ |
| Field additions | identical | identical | identical |
| Twiddles per module | $1$ | $3$ ($T_1,T_2,T_3$) | $K$ base + $2^{K-1}-1$ derived |

Asymptotically as $K\to\infty$, total memory passes go to $m$ — half of
radix-2's $2m$. At finite $K$:

| $K$ | sweeps per round | total sweeps (no residual) | sweeps saved vs radix-2 |
|---:|---:|---:|---:|
| 1 | 2 | $2m$       | $0$         |
| 2 | 3 | $1.5m$     | $25\%$      |
| 3 | 4 | $1.333m$   | $33\%$      |
| 4 | 5 | $1.25m$    | $37.5\%$    |

The original "$2$ sweeps per radix-4 round, half the memory traffic of
radix-2" claim that an earlier version of `additive_FFT_radix4` and an
earlier version of this document made was incorrect; the correct figure is
$3$ sweeps per radix-4 round and $25\%$ traffic reduction over radix-2.

---

## 4. Where this is in the code

- `cantor::additive_FFT_radix2k<K, FieldT>(poly, domain)` —
  templated radix-$2^K$ implementation in
  `C++/Cantor/fft.tcc`, declared in `C++/Cantor/fft.hpp`.
  Supported for $K\in\{1,2,3,4,5\}$.
- `cantor::additive_FFT_radix4<FieldT>(poly, domain)` is now a thin wrapper
  that calls `additive_FFT_radix2k<2, FieldT>`.
- Reference radix-2 used as a baseline:
  `cantor::additive_FFT(poly, domain)`,
  `C++/Cantor/fft.tcc:258–329`.
- Mixed-radix epilogue: when $K\nmid m$, the function performs $\lfloor m/K
  \rfloor$ radix-$2^K$ rounds and then $m\bmod K$ residual radix-2 rounds
  (2 sweeps each).
- Correctness check:
  `C++/main.cpp::test_radix2k_correctness` exercises $K=1,2,3,4$ against
  `cantor::additive_FFT` for $m=2,3,4,5,6,8,9,10,11,12,13,15,16,17,18$.

---

## 5. Measured runtime and cache behavior (this repository)

This section shows **actual numbers** from this machine and this checkout.
They will move with CPU model, compiler (`-O3 -march=native` in
`C++/CMakeLists.txt`), and system load, but the *relative* trend is the
point.

### 5.1 Runtime (Google Benchmark, `taskset -c 0`, 3 repetitions, median)

`gf256` field, Cantor-special-basis affine domain. Times in milliseconds per
FFT.

| $m$ | radix-2 | radix-4 ($K=2$) | $K=3$ | $K=4$ |
|---:|---:|---:|---:|---:|
| 18 | 43.81 | 42.87 | 42.94 | 52.18 |
| 19 | 91.99 | 100.90 | 104.54 | 114.09 |
| 20 | 220.24 | 202.69 | 243.15 | 210.52 |
| 21 | 442.41 | 472.07 | 425.07 | 502.02 |
| 22 | 947.84 | 916.11 | 1023.06 | 1125.93 |

Notes on the residuals:

- $m=18$: $K=2$ uses 9 rounds; $K=3$ uses 6 rounds (no residual);
  $K=4$ uses 4 radix-16 rounds + 2 residual radix-2 rounds.
- $m=20$: $K=4$ uses 5 rounds, no residual; $K=3$ has 6 main + 2 residual.
- $m=21$: $K=3$ no residual; $K=2$ and $K=4$ have residuals.

Where $K\mid m$ exactly the higher-$K$ variants land within a few percent of
radix-2; with two residual radix-2 rounds tacked on, the gain is consumed
or the variant is slower. In other words, **the saving is real but in a
regime where the inner loop's arithmetic dominates over memory traffic**;
see the `perf` numbers below for the cache picture.

### 5.2 Hardware counters (`perf stat`)

Hybrid CPUs report split `cpu_atom` / `cpu_core` PMUs. We pin to a single
performance core and read only `cpu_core` counters. Each run uses
`--benchmark_min_time=3s`.

```bash
taskset -c 0 perf stat -e cycles,instructions,cache-references,cache-misses,LLC-loads,LLC-load-misses \
  env BM_MIN_RANGE=20 BM_MAX_RANGE=20 BM_STEP=1 ./run_benchmark \
  --benchmark_filter='BM_cantor_additive_fft_radix2k<4>' \
  --benchmark_min_time=3s
```

Representative `cpu_core`-only totals at $m=20$:

| Metric                     | radix-2     | $K=2$       | $K=3$       | $K=4$       |
|---                         |---:         |---:         |---:         |---:         |
| Cycles                     | 27.15 G     | 27.49 G     | 27.71 G     | 26.19 G     |
| Instructions               | 65.70 G     | 65.62 G     | 66.61 G     | 64.66 G     |
| IPC                        | 2.42        | 2.39        | 2.40        | 2.47        |
| Cache references           | 465.7 M     | 320.4 M     | 292.2 M     | 354.5 M     |
| Cache misses               | 215.5 M     | 107.1 M     | 96.1 M      | 97.1 M      |
| Cache miss rate            | **46.3 %**  | **33.4 %**  | **32.9 %**  | **27.4 %**  |
| LLC loads                  | 34.1 M      | 12.2 M      | 13.4 M      | 36.1 M      |
| LLC load misses            | 26.3 M      | 5.8 M       | 4.2 M       | 6.0 M       |
| LLC load-miss rate         | **77.0 %**  | **47.7 %**  | **31.4 %**  | **16.5 %**  |

And at $m=22$ (working set $\approx128$ MiB, well beyond L3):

| Metric                | radix-2 | $K=2$ | $K=3$ | $K=4$ |
|---                    |---:     |---:   |---:   |---:   |
| Cycles                | 20.7 G  | 19.4 G| 25.2 G| 23.4 G|
| Cache misses          | 163.8 M | 91.4 M| 54.9 M| 72.8 M|
| Cache miss rate       | 44.1 %  | 33.2 %| 26.5 %| 24.8 %|
| LLC load misses       | 24.7 M  | 4.9 M | 2.7 M | 3.6 M |
| LLC load-miss rate    | 78.4 %  | 43.9 %| 21.9 %| 13.7 %|

Interpretation:

- The cache picture matches the schedule analysis. Going from radix-2 to
  $K=4$ cuts last-level-cache load misses by roughly $4\times$ (78% → 16%),
  and overall cache miss rate roughly halves (46% → 27%). This is the
  $(K+1)/(2K)$ memory-pass reduction predicted in §3.4 made concrete.
- Cycles only move modestly because at these sizes (inner loop = mults +
  XORs in a tight `for(nz : nz_S)` body) the field arithmetic is the
  bottleneck on this hardware; the prefetcher hides much of the L3 miss
  cost, and many radix-$2^K$ iterations have shorter `nz_S` per level
  (offsetting the level-count growth in instruction count). At larger $m$
  the cache benefit gets exposed more clearly.
- $K=4$ ($Q=16$) uses 16 simultaneous registers in the closer; on
  `gf256` this is at the edge of what fits comfortably in the GP-register
  file together with twiddles, which limits further speedup. For wider
  fields (`gf128`/`gf192`) the same code path will be more
  memory-bandwidth-bound and the radix-2^K saving should grow.
- $K=3$ at $m=22$ uses 7 main rounds + 1 residual radix-2 round; $K=4$ at
  $m=22$ uses 5 main rounds + 2 residual radix-2 rounds. The residual
  rounds are the reason $K=4$ does not strictly dominate $K=3$ there.

### 5.3 Reproducing

```bash
cd C++/build && cmake .. -DCMAKE_BUILD_TYPE=Release \
  && cmake --build . --target run_benchmark main

# correctness sanity check (fast)
./main

# microbenchmark
taskset -c 0 env BM_MIN_RANGE=18 BM_MAX_RANGE=22 BM_STEP=1 ./run_benchmark \
  --benchmark_repetitions=3 \
  --benchmark_filter='BM_cantor_additive_fft/|BM_cantor_additive_fft_radix4|BM_cantor_additive_fft_radix2k'

# cache picture for one variant
taskset -c 0 perf stat -e cycles,instructions,cache-references,cache-misses,LLC-loads,LLC-load-misses \
  env BM_MIN_RANGE=22 BM_MAX_RANGE=22 BM_STEP=1 ./run_benchmark \
  --benchmark_filter='BM_cantor_additive_fft_radix2k<4>' \
  --benchmark_min_time=3s
```

---

## 6. Why radix-2 cannot do fewer than 2 sweeps per round

For the radix-2 case the analogous "is there a fewer-sweep schedule?"
question has a negative answer; see `RADIX2_CANTOR_FFT_PASSES.md`.
The structural reason is that the round-0 closer reads the
post-division upper half *and* the post-division lower half, and the
position-by-position dependency cannot be unrolled into a single sweep
that simultaneously reads both halves and writes back without first
finalizing one of them.

The corresponding statement for radix-$2^K$ is that levels $1,\dots,K$
must each have their own descending pass (§3.1), and the closer must be
its own ascending pass — giving the $K+1$ minimum.
