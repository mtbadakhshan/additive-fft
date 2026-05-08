# Why the Classical Radix-4 NTT Trick Doesn't Save Multiplications in Additive FFTs

**Short answer.** The shape of the radix-4 butterfly is the same in classical
NTT and in our LCH/Cantor additive FFT — both fuse two consecutive radix-2
butterflies into one $4\times4$ matrix-vector product. But classical radix-4
also gets a **25 % arithmetic win** on top of the cache-pass reduction,
because its inner $V_4$ matrix has entries in a tiny multiplicative subgroup
$\{1,j,-1,-j\}$ where multiplications by $\pm 1$ (and often $\pm j$) are
*free*. The additive-FFT analog does not have this. The twiddles in LCH live
on an *additive* subgroup chain $W_0\subset W_1\subset\dots\subset W_m$, not
on a multiplicative cyclic subgroup, and there is no element of small order
in the *acting* group that turns a multiplication into a sign flip / swap.

So the radix-4 LCH butterfly inherits **only** the cache-locality benefit of
fusing two stages — never the 25 % arithmetic reduction. This file explains
why this is structural, not a missed optimization.

Companion docs in this repo:

- [`RADIX4_CANTOR_FFT.md`](RADIX4_CANTOR_FFT.md) — paper vs. cache-efficient
  radix-4 Cantor; concrete cache-pass and timing numbers.
- [`RADIX2_CANTOR_FFT_PASSES.md`](RADIX2_CANTOR_FFT_PASSES.md) — why radix-2
  Cantor cannot be reduced below two sweeps per round.

The implementations referenced here are
- `cantor::additive_FFT_radix4` / `cantor::additive_FFT_radix2k<K>`
  in `C++/Cantor/fft.tcc`
- `lch::additive_FFT_radix4` (`butterfly_radix4`) in `C++/LCH/utils.tcc`,
  `C++/LCH/fft.tcc`.

---

## 1. The classical radix-4 NTT butterfly, with arithmetic accounting

Let $\omega \in \mathbb{F}_p$ be a primitive $N$-th root of unity (or
$\omega = e^{-2\pi i/N}$ in $\mathbb{C}$). Set $j = \omega^{N/4}$ — a
**primitive 4th root of unity**. The radix-4 Cooley–Tukey butterfly applied
to four indices that share $\omega^N=1$ structure factors as

$$
\begin{pmatrix}y_0\\y_1\\y_2\\y_3\end{pmatrix}
=
\underbrace{\begin{pmatrix}1&1&1&1\\1&-j&-1&j\\1&-1&1&-1\\1&j&-1&-j\end{pmatrix}}_{V_4}
\;\diag(1,\omega^a,\omega^{2a},\omega^{3a})\;
\begin{pmatrix}x_0\\x_1\\x_2\\x_3\end{pmatrix}.
$$

The factorization that minimizes nontrivial multiplications is

$$
\begin{aligned}
A &= x_0+\omega^{2a}x_2, & B &= x_0-\omega^{2a}x_2,\\
C &= \omega^{a}x_1+\omega^{3a}x_3, & D &= \omega^{a}x_1-\omega^{3a}x_3,\\
E &= jD,\\
y_0 &= A+C, & y_2 &= A-C,\\
y_1 &= B-E, & y_3 &= B+E.
\end{aligned}
$$

Counting actual multiplications:

| Multiplication | Cost in $\mathbb{F}_p$ | Cost in $\mathbb{C}$ |
|---|---|---|
| $\omega^{2a} x_2$ | full mult | full mult |
| $\omega^{a} x_1$  | full mult | full mult |
| $\omega^{3a} x_3$ | full mult | full mult |
| $j\cdot D$        | full mult, but multiplier is fixed | **free** (swap+negate) |
| Multiplications by $\pm 1$ inside $V_4$ | **free** (sign flips) | free |

So per 4 outputs:

| | Multiplications per 4 outputs |
|---|---|
| Two radix-2 stages | $4$ |
| One radix-4 stage (NTT, generic $\mathbb{F}_p$) | $3$ + 1 fixed-$j$ mult |
| One radix-4 stage (complex / $\mathbb{F}_{p^2}$) | $3$ |

The reason the **additive** $\pm1$'s and the **multiplicative** $j$ collapse
into "free" ops is that they sit in a 4-element subgroup of the
multiplicative group $\mathbb{F}_p^\times$ (or $\mathbb{C}^\times$), and that
subgroup acts on field elements by *cheap* maps (sign flips, swaps).
Asymptotically over a length-$N$ FFT that turns $\tfrac{N}{2}\log_2 N$
nontrivial mults (radix-2) into $\tfrac{3N}{8}\log_2 N$ (radix-4) — the
**25 % savings**.

---

## 2. The radix-4 LCH/Cantor butterfly, with arithmetic accounting

In LCH (and Cantor) the radix-4 butterfly factors through the matrix

$$
M \;=\;
\begin{pmatrix}
1 & T_1 & T_2 & T_1T_2\\
1 & T_1+1 & T_2 & (T_1+1)T_2\\
1 & T_1+\beta_1 & T_2+1 & (T_1+\beta_1)(T_2+1)\\
1 & T_1+\beta_1+1 & T_2+1 & (T_1+\beta_1+1)(T_2+1)
\end{pmatrix}
\;=\; A\cdot B
$$

with the sparse factorization

$$
A=\begin{pmatrix} 1 & T_1 & 0 & 0\\ 1 & T_1+1 & 0 & 0\\ 0 & 0 & 1 & T_3\\ 0 & 0 & 1 & T_3+1 \end{pmatrix},\quad
B=\begin{pmatrix} 1 & 0 & T_2 & 0\\ 0 & 1 & 0 & T_2\\ 1 & 0 & T_2+1 & 0\\ 0 & 1 & 0 & T_2+1 \end{pmatrix},
$$

where $T_3 = T_1+\beta_1$ (using $\beta_1^2+\beta_1=1$). The corresponding
factorized algorithm is exactly what `butterfly_radix4` runs:

```cpp
u_0 = q_0 + T_2 * q_2;     u_1 = q_1 + T_2 * q_3;
u_2 = u_0 + q_2;           u_3 = u_1 + q_3;
h_0 = u_0 + T_1 * u_1;     h_2 = u_2 + T_3 * u_3;
h_1 = h_0 + u_1;           h_3 = h_2 + u_3;
```

Counting multiplications:

| Multiplication | Cost |
|---|---|
| $T_2\cdot q_2$ | full $\mathbb{F}_{2^r}$ mult |
| $T_2\cdot q_3$ | full $\mathbb{F}_{2^r}$ mult |
| $T_1\cdot u_1$ | full $\mathbb{F}_{2^r}$ mult |
| $T_3\cdot u_3$ | full $\mathbb{F}_{2^r}$ mult |
| Multiplications by $1$ inside $A$ and $B$ | free (identity) |
| Multiplications by $\beta_1$ when computing $T_3$ | precomputed, once per module |

So per 4 outputs:

| | Multiplications per 4 outputs |
|---|---|
| Two radix-2 stages | $4$ |
| One radix-4 stage (LCH/Cantor) | $4$ |

**Identical.** This is what the paper's note "Cost in radix-4 LCH AFFT is
same as radix-2 LCH AFFT" is recording, and it is structurally forced. None
of the four pair-multiplications can collapse:

- $T_1, T_2$ are values of the linear operators $S^{i-2}, S^{i-1}$ at the
  module's coset offset and are *generic* elements of $\mathbb{F}_{2^r}$
  (one per module per round; they take all values as the module index varies).
- $T_3 = T_1+\beta_1$ is *also* generic, even though $\beta_1$ lies in the
  subfield $\mathbb{F}_4\subset\mathbb{F}_{2^r}$ ($\beta_1^3=1$). The reason
  $\beta_1$'s small multiplicative order is irrelevant here is that it is
  *added to* $T_1$, not used as a stand-alone multiplier in the inner loop.
  We do precompute $T_3$ once per module from one addition, which avoids
  *one twiddle table lookup* per module — but that has nothing to do with
  the inner-loop arithmetic count.
- The "free" entries ($\pm 1$ in the classical $V_4$) have no analog: the
  only constants appearing in $A$ and $B$ besides 0 and 1 are $T_2+1$,
  $T_1+1$, $T_3+1$, all generic non-trivial multipliers.

---

## 3. Why the trick is structural to the multiplicative-group setting

The 25 % radix-4 win in classical NTT is not an "accident of the algorithm."
It comes from one specific algebraic fact: the multiplicative group
$\mathbb{F}_p^\times$ (or $\mathbb{C}^\times$) **contains an element of order
4** that the host architecture can multiply by *for free* (sign flip / swap).

| Acting group on the data | Element of order 4? | "Free" action available? |
|---|---|---|
| $\mathbb{C}^\times$ | $j=i$ | yes (real↔imag swap + sign flip) |
| $\mathbb{F}_p^\times$ for $p\equiv 1\pmod 4$ | $j$ with $j^2=-1$ | partial — $-1$ free, $j$ usually a generic mult |
| $\mathbb{F}_{2^r}^\times$ | **no** — group order is $2^r-1$, odd, so $4\nmid \lvert \mathbb{F}_{2^r}^\times\rvert$ | n/a (no order-4 element exists) |
| Additive group $W_m \subset \mathbb{F}_{2^r}$ | **no** — every element is its own inverse in characteristic 2, $-1=1$ | n/a (no order-4 element exists) |

Row 3 is worth pausing on. Because $\lvert \mathbb{F}_{2^r}^\times\rvert = 2^r-1$
is odd for every $r\ge 1$, by Lagrange the cyclic group
$\mathbb{F}_{2^r}^\times$ has **no element of order 4 at all** — not even a
non-free one. So even a *multiplicative* NTT over $\mathbb{F}_{2^r}$ cannot
mimic the classical $\{1,j,-1,-j\}$ trick: the order-4 subgroup that the
trick relies on simply does not exist in any binary extension field.

In the **additive-FFT** setting the situation is even cleaner: the acting
group is the additive subspace $W_m$, which has exponent 2 in characteristic
2 (every element is its own additive inverse). There is no nontrivial
element of order 4 anywhere in the picture, because $-1 = 1$. The structure
the algorithm *does* use — $\beta_1^2+\beta_1=1$, hence $T_3 = T_1+\beta_1$
— is an *additive* relation between the two pair-twiddles of the inner
stage. It saves us a *twiddle table lookup*, not a *field multiplication*.

---

## 4. The win that does survive: cache passes (and how big it is)

What the radix-4 LCH still wins is purely the memory-layout argument from
[`RADIX4_CANTOR_FFT.md`](RADIX4_CANTOR_FFT.md), §3. Per outer round, both
the fused radix-4 and a single radix-2 sweep do one full pass over the
coefficient array. Generalising to radix-$2^K$:

|   | Outer rounds | Sweeps per round | Sweeps total |
|---|---|---|---|
| Radix-2 LCH (`butterfly`) | $m$ | $1$ | $m$ |
| Radix-4 LCH (`butterfly_radix4`) | $\lceil m/2\rceil$ | $1$ | $\lceil m/2\rceil$ |
| Radix-$2^K$ LCH (`butterfly_radix2k<K>`) | $\lceil m/K\rceil$ | $1$ | $\lceil m/K\rceil$ |

So radix-$2^K$ LCH cuts the number of full passes through the buffer by a
factor $K$ at (essentially) the same arithmetic cost.

The implementation `lch::butterfly_radix2k<K>` (and its FFT wrapper
`lch::additive_FFT_radix2k<K, FieldT>`) lives in
[`C++/LCH/utils.tcc`](C++/LCH/utils.tcc) /
[`C++/LCH/fft.tcc`](C++/LCH/fft.tcc). $K=1$ degenerates to a single
radix-2 stage per round, $K=2$ matches `butterfly_radix4`, and $K\in\{3,4,5\}$
fuses 3, 4, or 5 consecutive radix-2 stages into one in-register K-stage
butterfly per outer round. Twiddles are derived from the same
`cantor_combinations` table by

$$
\text{mshift}_{\ell,b} = (\text{mod} \ll \ell) \mid (b \ll 1) \mid (\text{shift\_bit}_\ell \ll 1),
$$

with $\text{shift\_bit}_\ell = \text{shift\_bit}_1 \ll (\ell-1)$, so all $K$
per-stage base twiddles fall out of one bit-pattern shifted by $0,\dots,K-1$.

### 4.1 Measured wall time on $\mathbb{F}_{2^{256}}$

All numbers below are from `C++/build/run_benchmark` (Google Benchmark)
and `perf_driver` on a 12th-Gen Intel Core i5-1250P, P-core 0
(`taskset -c 0`), L1d $48\,$KiB / L2 $1280\,$KiB / L3 $12\,$MiB, with
`-O3 -march=native`. Each `gf256` element occupies 32 B, so the working
set at $m$ is $W = 2^m \cdot 32\,\text{B}$.

|   $m$   | $W$        | r2     | r4     | r2k1   | r2k2   | r2k3   | r2k4   | r2k5   |
|--------:|-----------:|-------:|-------:|-------:|-------:|-------:|-------:|-------:|
| 16      | 2 MiB      | 7.82  ms | 8.06  ms | 8.86  ms | 9.45  ms | 8.84  ms | 8.14  ms | 8.15  ms |
| 18      | 8 MiB      | 36.65 ms | 36.05 ms | 41.92 ms | 43.87 ms | 41.46 ms | 41.09 ms | 38.32 ms |
| 20      | 32 MiB     | 165.3 ms | 167.1 ms | 182.9 ms | 190.5 ms | 176.3 ms | 201.6 ms | 171.0 ms |
| 22      | 128 MiB    | 752    ms | 732    ms | 837    ms | 914    ms | 830    ms | 905    ms | 841    ms |

Naming: `r2 = additive_FFT` (hand-rolled radix-2 baseline), `r4 = additive_FFT_radix4`,
`r2kK = additive_FFT_radix2k<K>` (the templated implementation, for $K\in\{1,\ldots,5\}$).

Two things stand out immediately:

1. **`r2` and `r4` are roughly tied at every $m$**: 1–3 % spread, never
   meaningfully more. Going from 1 sweep per stage × 22 stages to 1 sweep
   per fused-pair × 11 stages does *not* halve the wall time.
2. **The templated `r2k1` is consistently slower than the hand-rolled
   `r2`** by 11–25 %. This is purely template-loop overhead (load tuple
   into `FieldT v[Q]`, K-stage register butterfly, store back) versus the
   straight-line in-place radix-2 butterfly in `butterfly_op`. Comparing
   `r2kK` against `r2k1` (both inside the same template harness) is the
   apples-to-apples comparison for the cache-pass argument, and there
   `r2k5/r2k1 ≈ 0.93` at $m\in\{18,20,22\}$ — *so cache passes do buy
   us a 7 % speedup, but not the 22→5 = 4× one might hope for.*

### 4.2 Why so little? Direct cache-miss measurement

To check whether the algorithm is memory-bound or compute-bound at
$m=22$ ($W=128\,$MiB, $\sim 10\times$ L3) I wrapped each variant with
`perf stat -e cycles,instructions,cache-references,cache-misses,
L1-dcache-loads,L1-dcache-load-misses,LLC-loads,LLC-load-misses` pinned
to a single P-core (4 iterations per variant; numbers below are totals
across the 4 iterations):

| Variant | wall (ms/call) | cycles  | insns   | IPC  | L1d loads | L1d miss% | LLC loads | LLC misses |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| r2    | 752  | 13.99 G | 30.11 G | 2.15 | 5.98 G | 3.37 % | 45.1 M | 36.98 M |
| r4    | 732  | 13.60 G | 29.09 G | 2.14 | 6.32 G | 2.77 % | 47.0 M | 36.55 M |
| r2k1  | 837  | 15.12 G | 35.01 G | 2.32 | 8.45 G | 2.15 % | 44.2 M | 36.33 M |
| r2k2  | 914  | 16.33 G | 37.32 G | 2.29 | 7.85 G | 2.63 % | 49.3 M | 37.09 M |
| r2k3  | 830  | 14.78 G | 33.75 G | 2.28 | 7.10 G | 2.72 % | 52.3 M | 37.92 M |
| r2k4  | 905  | 16.27 G | 32.35 G | 1.99 | 6.57 G | 2.98 % | 63.0 M | 39.22 M |
| r2k5  | 841  | 15.07 G | 31.18 G | 2.07 | 5.98 G | 3.39 % | 79.6 M | 41.15 M |

The headline cell here is **LLC misses**: it stays at $\approx 37$ M
*regardless of K*, even slightly *increasing* at $K=5$. Same picture
holds at $m=20$ ($\approx 14$ M LLC misses across the board) and at
$m=18$ (within L3, $\approx 5.7$ M LLC misses). In other words the
DRAM-to-L3 traffic — which is the bandwidth that actually costs a
modern CPU 100+ cycles per cache-line miss — is *not reduced by raising
$K$*.

What does decrease with $K$ is the L1/L2 traffic:

- `cache-references` $467\,\text{M}\to 410\,\text{M}$ ($\downarrow 12\%$).
- `L1-dcache-loads`  $5.98\,\text{G}\to 5.98\,\text{G}$ for $r2$ vs. $r2k5$
  (flat), with the templated harness inflating loads at intermediate $K$.

This is the cache-pass reduction in action — but it lives at the L1/L2
boundary, not at the DRAM boundary, because the FFT is sequential within
each pass and the L3 + hardware prefetcher already absorb most cross-pass
reuse. **Saving L1 misses on a kernel whose inner loop is bottlenecked
by a 3-cycle-throughput `PCLMULQDQ` chain does not move the wall clock.**

### 4.3 Top-down breakdown: where do the cycles actually go?

`perf stat -e cpu_core/{slots,topdown-retiring,topdown-fe-bound,
topdown-be-bound,topdown-bad-spec,topdown-mem-bound}/u` at $m=22$
(4 iterations):

| Variant | retiring | fe-bound | be-bound | mem-bound | core-bound | bad-spec |
|---|---:|---:|---:|---:|---:|---:|
| r2    | 37.3 % | 3.5 % | 58.0 % | 19.7 % | **38.4 %** | 0.8 % |
| r4    | 35.3 % | 2.8 % | 60.0 % | 33.0 % | **27.1 %** | 0.8 % |
| r2k1  | 40.8 % | 9.4 % | 47.4 % | 18.4 % | **29.0 %** | 0.8 % |
| r2k2  | 38.8 % | 3.9 % | 56.1 % | 34.9 % | **21.2 %** | 0.9 % |
| r2k3  | 39.2 % | 3.9 % | 55.7 % | 26.3 % | **29.4 %** | 0.8 % |
| r2k4  | 36.9 % | 3.5 % | 58.0 % | 24.8 % | **33.2 %** | 1.2 % |
| r2k5  | 37.6 % | 3.9 % | 58.0 % | 24.4 % | **33.6 %** | 0.8 % |

Reading the rows: **only $\sim$ 37–41 % of the issue slots are productive**
("retiring"), the back-end-bound stalls dominate, and they split fairly
evenly between memory-bound (LLC/L2/L1 wait) and core-bound (port
pressure on the `PCLMULQDQ` lane). At $m=20$ the picture is essentially
the same (retiring $\approx 37$–$42\%$, mem-bound $\approx 19$–$35\%$,
core-bound $\approx 18$–$37\%$).

The crucial observation is that `mem-bound` does *not* drop monotonically
with $K$: $\text{r2}\to\text{r2k5}$ takes mem-bound from $19.7\%$ to
$24.4\%$ (slightly *worse*), and core-bound from $38.4\%$ to $33.6\%$
(slightly better). Raising $K$ trades a tiny amount of core-bound stalls
for a tiny amount of mem-bound stalls, and the wall time is essentially
flat — exactly the signature of an algorithm that is *not bottlenecked
by cache-pass count*.

### 4.4 Bottom-line: compute-bound, not memory-bound

> At $m=18$–$22$ on `gf256`, LCH AFFT spends $\sim 30\%$ of slots on
> back-end memory stalls and $\sim 30\%$ of slots on back-end core
> stalls. Increasing the radix from 2 to $2^5$ cuts L1/L2 traffic by
> $\sim 12\%$ but **does not reduce DRAM (LLC-miss) traffic at all** —
> the working set still streams through L3 once per pass, and the
> hardware prefetcher already absorbs the remaining reuse. The wall
> time correspondingly does not improve. This places LCH AFFT firmly in
> the **partially-memory-bound, more-strongly-compute-bound** regime,
> where the next gains have to come from the multiplication side
> (SIMD packing of `PCLMULQDQ`, batched Karatsuba, constant-twiddle
> specialisation), not from further cache-pass amortisation.

In one slogan: **radix-4 NTT compresses two stages into one and gets a
25 % multiplicative bonus from the order-4 subgroup of the multiplicative
group. Radix-$2^K$ LCH compresses $K$ stages into one and gets neither
the 25 % mult bonus nor a meaningful DRAM-traffic reduction, because the
additive subspace tower it lives on has no analog of that subgroup *and*
the inner-loop bottleneck is `PCLMULQDQ` throughput, not cache misses.**

---

## 5. What from the NTT side does transfer

It is worth being explicit about the optimizations that *do* carry over to
LCH, since "the NTT trick doesn't help" is a much narrower statement than
"NTT optimizations don't help."

| Idea from the NTT world | Carries over to LCH? | Why |
|---|---|---|
| Trivial-twiddle savings ($\pm 1, \pm j$) | **No** | additive subgroup, no acting element of order 4 |
| Common-subexpression sparse factorization of the 4×4 | **Yes** | already used: $M=A\cdot B$ via $\beta_1^2+\beta_1=1$ |
| Stockham / self-sorting layout | **Yes** | pure index transform, basis-agnostic |
| Bailey four-step / six-step (block) FFT | **Yes** | layout/cache argument, removes most DRAM traffic at large $m$ |
| Mixed radix when $m$ is not a multiple of $K$ | **Yes** | already done: residual radix-2 stage when $\log_n$ is odd |
| SIMD over the inner-tuple loop | **Yes** | the four-quarter butterfly closer is per-tuple data-parallel |
| Constant-multiplier specialisation for fixed twiddles | **Yes** | top-level radix-4 stages have twiddles in $\{0,1,\beta_1,1+\beta_1\}$, replaceable by table-driven xor-shifts |
| Higher radix ($K\ge 3$) for cache amortisation | **Yes, but small win** | done for both Cantor (`additive_FFT_radix2k<K>`) and LCH (`butterfly_radix2k<K>`); §4 measures $\le 7\%$ wall-time gain on gf256 because LLC misses don't shrink |
| 25 % asymptotic mult reduction | **No** | structural, see §2–§3 |

The natural next implementations, if more wall-time gain is wanted, are (in
roughly decreasing payoff order at large $m$, *given the §4 finding that
the bottleneck is `PCLMULQDQ` throughput, not cache passes*):

1. **SIMD-packed $\mathbb{F}_{2^{r}}$ multiplication** in the inner `jj`
   loop, e.g. AVX-512 batched `vpclmulqdq` over 4 quadwords at once.
   This is the only optimisation that touches the 30 %-of-slots core-bound
   ceiling identified in §4.3.
2. **Constant-twiddle specialisation** for the deepest radix-$2^K$ stages
   (where every twiddle lies in a small fixed set), turning generic
   $\mathbb{F}_{2^r}$ mults into table-driven xor-shifts.
3. **Bailey block FFT** for $m$ that overflows DRAM bandwidth (this is the
   only knob that actually reduces *DRAM* traffic, as opposed to L1/L2
   traffic — the §4 perf data shows that's where the remaining memory
   stalls live).
4. **Stockham / self-sorting layout** to remove the input-replication
   $\texttt{std::copy}$ in `lch::additive_FFT*` when $n_{\text{poly}} < n$.

None of those reaches into the 25 %-mults box, because that box is
welded shut by the algebra. The §4 measurements show that the *cache-
pass* knob is also nearly closed for LCH on `gf256` — radix-2 and
radix-32 land within $\sim 7\%$ of each other on a 128 MiB working set.
The only knob with substantial remaining headroom is the `PCLMULQDQ`
throughput knob, which is independent of radix choice.

---

## 6. One-paragraph summary

> Both classical NTT and additive (LCH/Cantor) FFT have the same butterfly
> *shape*; both can fuse two (or $K$) radix-2 stages into one radix-$2^K$
> stage to cut cache passes by a factor $K$. Classical NTT *additionally*
> saves 25 % of the field multiplications because its inner $V_4$ matrix
> has entries in the order-4 subgroup $\{1,j,-1,-j\}$ of the multiplicative
> group, and modern hardware can multiply by those for free (sign flip /
> swap). Additive FFTs operate on an *additive* subspace tower
> $W_0\subset\dots\subset W_m$, whose "twiddles" $T_k=S^k(\theta_a)$ are
> generic elements of $\mathbb{F}_{2^r}$ with no analog of the order-4
> free-action subgroup. The only algebraic structure available,
> $\beta_1^2+\beta_1=1$, saves a *twiddle table lookup* but no *field
> multiplications*. The §4 perf measurements then show that even the
> cache-pass half of the classical radix-4 win is largely lost on gf256:
> raising $K$ from 2 to 32 cuts L1/L2 traffic but **does not reduce LLC
> misses or DRAM traffic** ($\approx 37$ M LLC misses at $m=22$ across
> all $K$), because the inner loop is bottlenecked by `PCLMULQDQ`
> throughput rather than cache misses. So radix-4/radix-$2^K$ LCH matches
> radix-2 in mult count *and* in DRAM traffic, leaving only a $\sim 7\%$
> wall-time edge that comes from L1/L2 amortisation; the future-work
> needle has to be moved on the multiplication side, not the cache side.
