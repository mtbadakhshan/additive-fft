# Can the Radix-2 Cantor AFFT Be Reduced Below Two Sweeps Per Round?

**Short answer: no.** The existing radix-2 implementation in
`C++/Cantor/fft.tcc` (`additive_FFT(poly, domain)`, lines 258–329) already
performs the minimum number of array sweeps possible for a pure radix-2
decomposition. The only way to reduce sweeps further is to merge multiple
decomposition levels per round, i.e. go to radix-4 (already implemented in
`additive_FFT_radix4`) or higher.

This file is the cache-pass companion to
[`RADIX4_CANTOR_FFT.md`](RADIX4_CANTOR_FFT.md). No new function was added,
because no useful one exists.

---

## 1. What the current radix-2 round looks like

For each round $r\in\{0,\dots,m-1\}$, every module of size $L=2^{m-r}$ is
processed in **two sweeps**:

```cpp
// Sweep A: descending k. Division by S^{m-r-1} + cross-half twiddle by T.
for (size_t k = offset2 + half - 1; k >= offset2; --k) {
    FieldT gk = g[k];
    for (const auto& nz : nz_S)
        g[k - nz] += gk;          // division writes (lower or upper half)
    g[k - half] += gk * T;        // twiddle write (lower half only)
}

// Sweep B: ascending j. Closer h_1[j] = u_0[j] + f_1[j].
for (size_t j = 0; j < half; ++j)
    g[offset2 + j] += g[offset + j];
```

After Sweep A:
- lower half holds $u_0(x) = f_0(x) + T \cdot f_1(x)$;
- upper half holds the quotient $f_1(x)$.

After Sweep B:
- lower half = $h_0 = u_0$;
- upper half = $h_1 = u_0 + f_1$.

This is two array passes per module. With $m$ rounds, total is $2m$ passes.

---

## 2. Why fusion into one sweep is impossible

Sweep B's iteration $j$ reads `g[offset + j]` and expects it to hold the
**finalised** $u_0[j]$ — i.e., after every Sweep A write to position $j$ has
landed. The question is therefore: when, during Sweep A's descending $k$ loop,
does each lower-half index $j$ become final?

### Which Sweep A iterations write to lower-half index $j$?

Iteration $k$ of Sweep A ($k\in[L/2,\,L)$) writes to the following positions:

- `g[k - nz_i]` for each $nz_i \in$ `nz_S`. The set `nz_S` always satisfies
  $0 < nz_i < L/2$ (every entry is strictly below the leading degree of
  $S^{m-r-1}$, see `s_i[]` and `n_terms[]` in `Cantor/cantor_basis.hpp`).
- `g[k - L/2]` (the twiddle write).

Translating: position $j$ in the lower half ($0\le j < L/2$) is written by
iteration $k$ iff
$$
k = j + L/2  \qquad\text{(twiddle)} \qquad\text{or}\qquad
k = j + nz_i \;\text{ for some } nz_i\in\texttt{nz\_S}.
$$

### When does position $j$ become final?

Sweep A processes $k$ in **descending** order. So position $j$ becomes final at
the *smallest* such $k$. Two regimes:

1. **Position $j=0$:**
   - $k = 0 + L/2 = L/2$ from the twiddle.
   - $k = 0 + nz_i$: requires $nz_i \ge L/2$ (so that $k$ is in the loop range
     $[L/2, L)$). But every entry of `nz_S` is strictly less than $L/2$.
     **Hence no `nz` write ever targets $j=0$.**

   The only write to position $0$ is the twiddle at $k=L/2$, the **very last
   iteration of Sweep A.**

2. **General $j$:** position $j$ first becomes final at
   $$k_\text{final}(j) = \min\Bigl(j + L/2,\; \min\{j + nz_i : nz_i\ge L/2 - j\}\Bigr).$$

   For all $j$ with $j < L/2 - \max(\texttt{nz\_S})$ this set is just
   $\{j+L/2\}$, so $k_\text{final}(j) = j+L/2$ — the *latest possible*
   iteration that could touch $j$.

### Why this kills fusion

To fold Sweep B into Sweep A, we would need Sweep B's iteration $j$ to be
schedulable *before* Sweep A finishes. But:

- For $j=0$, finalisation happens at $k=L/2$ — the **final** iteration of
  Sweep A. There is no Sweep A iteration left after that to do the closer in.
- For most $j$ (specifically $j < L/2 - \max \texttt{nz\_S}$), finalisation also
  happens at $k=j+L/2$, which is the latest iteration that could touch them.

So the closer for a typical index $j$ has to wait until effectively the end of
Sweep A. There is no way to interleave Sweep B's reads with Sweep A's writes
that produces fewer than two sweeps over the module.

This contrasts sharply with the radix-4 case, where the equivalent fusion *is*
sound because each radix-4 round contains *two* divisions whose write
patterns can be merged into a single descending sweep, and the pair-of-stages
butterfly closer is performed in one ascending in-register sweep over four
quarters. The win is *not* a reduction below two sweeps per round; it is
performing two decomposition levels per (still two-sweep) round, halving
$m$.

### A purely partial fusion is possible but useless

For the small subset of indices $j \in [L/2 - \min(\texttt{nz\_S}),\, L/2)$,
we can do the closer for index $j$ inside Sweep A immediately after iteration
$k = j + \min(\texttt{nz\_S})$ — once that iteration writes the last `nz`
contribution to position $j$, $u_0[j]$ is final. But this covers at most
$\min(\texttt{nz\_S}) \le L/2 - 1$ of the $L/2$ indices — typically a tiny
fraction (e.g. for `s_i[20] = {1048575, 1048560, 983040}` and $L/2=2^{19}$,
$\min = 65536$, so we'd cover $65536/524288 \approx 12\%$ of indices and
still need a second sweep for the remaining 88%).

It does not turn two sweeps into one; it just shaves a few extra iterations
off the second sweep at the cost of extra branching, and is therefore not
worth implementing.

---

## 3. What the radix-2 round actually optimises

The current implementation already *does* fuse the two operations that *can*
be safely fused inside one sweep:

- The **division by $S^{m-r-1}$** and the **cross-half twiddle** are merged
  into Sweep A's single descending body. They share the loaded $g_k$, share
  the iteration counter $k$, and respect the `writes-go-strictly-below-k`
  invariant. This fusion is identical in spirit to the merge inside the
  radix-4 Sweep 1.
- Sweep B is then the unavoidable closer.

So at the level of "what can be merged in radix-2," the existing code is
already optimal. A naïve transcription of the pseudocode would have *three*
sweeps per round (division, then twiddle, then closer); the current code is
already at the minimum **two** for radix-2 by fusing the first two.

---

## 4. The only way to reduce sweeps further: higher radix

The reason `additive_FFT_radix4` runs fewer sweeps in total is **not** that it
breaks the two-sweeps-per-round floor. Each radix-4 round is also two sweeps.
The win is that one radix-4 round does the work of two radix-2 rounds, so the
total round count drops from $m$ to $\lceil m/2\rceil$, and the total sweep
count from $2m$ to $m + \mathbb{1}[m\text{ odd}]$.

| | Radix-2 (current) | Radix-4 (`additive_FFT_radix4`) |
|---|---|---|
| Decomposition levels per round | 1 | 2 |
| Sweeps per round | 2 (already minimum) | 2 (minimum for radix-4 too) |
| Rounds total | $m$ | $\lceil m/2 \rceil$ |
| Sweeps total | $2m$ | $m + [m\text{ odd}]$ |
| Lower bound from data dependency | $2m$ ✗ | $m + [m\text{ odd}]$ ✗ |

The same argument as in §2 applies recursively: any radix-$2^k$ butterfly will
need exactly two sweeps per round (one descending merged-divisions sweep, one
ascending fused-butterfly sweep), so the path to fewer sweeps is to make $k$
larger, not to break the two-sweeps-per-round barrier.

---

## 5. Conclusion

Summary:

1. The existing `cantor::additive_FFT(poly, domain)` already runs at **two
   sweeps per round per module**, which is the **lower bound** for any
   radix-2 Cantor butterfly, dictated by a hard data dependency between the
   division/twiddle phase and the closer.
2. A new function reducing this to one sweep cannot be written, because for
   the typical index $j$ the closer's input $u_0[j]$ does not become final
   until the *last* iteration of the descending division sweep.
3. The route to fewer cache passes is therefore higher radix, which is
   already provided by `cantor::additive_FFT_radix4(poly, domain)` —
   approximately halving the total number of sweeps without changing the
   field-arithmetic operation count.

So no new radix-2 function was added.
