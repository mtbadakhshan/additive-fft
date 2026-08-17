# Benchmarking

The benchmark driver is implemented in `main.c`, which can be modified to select the FFT algorithms and parameter sizes to benchmark.

The **Dyadic AFFT** includes an optional tuning step that determines suitable implementation parameters for the target machine. The following instructions benchmark the Dyadic AFFT against the LCH AFFT implementation provided by [`bitpolymul`](https://github.com/fast-crypto-lab/bitpolymul).

### 1. Initialize dependencies

From the repository root, initialize `bitpolymul` and the other Git submodules:

```bash
git submodule update --init --recursive
```

### 2. Enter the C implementation directory

```bash
cd C
```

### 3. Tune the Dyadic AFFT (optional)

The tuning step can be run once for the target machine. To reduce scheduling variability, the process can be pinned to a single CPU core:

```bash
taskset -c 0 make tune
```

The resulting tuning parameters are then used by the Dyadic AFFT implementation during benchmarking.

### 4. Build and run the benchmark

Compile the benchmark:

```bash
make
```

Then run it on the same CPU core:

```bash
taskset -c 0 ./main.out
```

Using `taskset` is optional, but pinning both tuning and benchmarking to the same CPU core helps reduce run-to-run variability and makes the measurements more reproducible.

