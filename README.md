# ksieve

A specialized integer factoring program for semiprimes `N = u·v` where both prime factors are small multiples of 6 ± 1. Three complementary sieves:

- **K-sieve** and **S-sieve** — factor N ≡ 1 (mod 6), where both factors satisfy p ≡ −1 (mod 6)
- **M-sieve** — factors N ≡ 5 (mod 6), where one factor is ≡ −1 (mod 6) and the other is ≡ +1 (mod 6)

Together these cover every prime-pair semiprime with factors greater than 3. On a 16-thread build, 100-bit semiprimes factor in roughly 0.24 s (S-sieve) to 0.61 s (M-sieve) mean wall-time; 120-bit in roughly 54 s mean (S-sieve).

The sieves share the core mathematical structure of Hittmeir's hyperbolic sieve for Fermat factorization ([arXiv:2205.10074](https://arxiv.org/abs/2205.10074)) — filter candidates by testing whether a specific quadratic polynomial is a residue mod a composite modulus assembled from small primes via CRT — but apply it to the `6n±1` parameterization and add Hensel lifting at p²/p³/p⁴, a tiled 64-bit-word bitset exploiting 2-adic structure, AVX-512 batch-16 acceleration, and multi-threaded prefix work-stealing. See the whitepaper for the full mathematical treatment and benchmarks.

---

## Quickstart

```bash
# Install dependencies (Debian/Ubuntu)
sudo apt install build-essential libgmp-dev

# Build a 16-thread release binary
gcc -O3 -march=native -DKSIEVE_N_THREADS=16 -o ksieve main.c -lgmp -lm -lpthread

# Factor a specific N (auto-detects old/new form from N mod 6)
KSIEVE_AUTO=1 ./ksieve 941635629865253

# Benchmark 10 random 100-bit new-form semiprimes
KSIEVE_MSIEVE=1 ./ksieve --bench 100 10

# See the built-in performance profile
./ksieve --info
```

---

## Requirements

- GCC or Clang with C99 support and `__int128` (any GCC ≥ 4.6 or Clang ≥ 3.2)
- GMP development headers (`libgmp-dev` on Debian/Ubuntu, `gmp-devel` on Fedora, `gmp` on macOS Homebrew)
- POSIX threads
- x86-64. AVX-512 is strongly recommended but optional — the build falls back to a scalar path when AVX-512 isn't available

No other dependencies. The build is a single `gcc` invocation on one `.c` file.

---

## Building

The number of threads is a **compile-time** setting. Pick the value that matches your hardware:

```bash
# Single-threaded (useful for benchmarking per-core performance)
gcc -O3 -march=native -DKSIEVE_N_THREADS=1 -o ksieve_t1 main.c -lgmp -lm -lpthread

# 4 threads (typical laptop)
gcc -O3 -march=native -DKSIEVE_N_THREADS=4 -o ksieve_t4 main.c -lgmp -lm -lpthread

# 16 threads (desktop / workstation)
gcc -O3 -march=native -DKSIEVE_N_THREADS=16 -o ksieve_t16 main.c -lgmp -lm -lpthread
```

Default if `-DKSIEVE_N_THREADS=...` is omitted: **2**. The hard cap is 16.

### Optional: instrumented build for profiling

Add `-DINSTRUMENT` to enable cycle-accurate per-filter rejection counters:

```bash
gcc -O3 -march=native -DKSIEVE_N_THREADS=16 -DINSTRUMENT \
    -o ksieve_instr main.c -lgmp -lm -lpthread
```

The instrumented binary prints a detailed breakdown after `--bench-instr` runs (see below).

### Build sanity check

```bash
./ksieve --demo
```

This runs correctness tests from 40-bit through 100-bit on known inputs. All lines should end in `OK`.

---

## Runtime options

### CLI subcommands

| Command | Purpose |
|---|---|
| `./ksieve <N>` | Factor a specific N (decimal). |
| `./ksieve --demo` | Run built-in correctness demos across 40–100 bit. |
| `./ksieve --bench <bits> <trials>` | Generate `<trials>` random `<bits>`-bit semiprimes and benchmark. |
| `./ksieve --bench-instr <bits> <trials>` | Same as `--bench` but prints instrumented counters (requires `-DINSTRUMENT` at build time). |
| `./ksieve --info` | Print the built-in performance-profile table. |
| `./ksieve --help` | Print usage. |

### Mode-selection environment variables

These select which sieve to run. At most one should be set; `KSIEVE_AUTO` takes precedence only if none of the explicit flags are set.

| Variable | Values | Effect |
|---|---|---|
| `KSIEVE_MSIEVE=1` | 0 or 1 | Use M-sieve (for N ≡ 5 mod 6, new form with mixed-class factors). |
| `KSIEVE_SSIEVE=1` | 0 or 1 | Use S-sieve (for N ≡ 1 mod 6, old form, ∼4.2× smaller scan range than K-sieve). |
| `KSIEVE_AUTO=1` | 0 or 1 | If neither of the above is set: auto-select based on N mod 6 (M-sieve for N ≡ 5, S-sieve for N ≡ 1). |
| _(none)_ |  | Default: K-sieve. |

Mode examples:

```bash
# New-form semiprime (e.g. u ≡ 5, v ≡ 1 mod 6)
KSIEVE_MSIEVE=1 ./ksieve --bench 100 10

# Old-form semiprime, wider-range K-sieve (legacy, rarely optimal)
./ksieve --bench 80 10

# Old-form semiprime, faster S-sieve
KSIEVE_SSIEVE=1 ./ksieve --bench 80 10

# Factor an unknown N, let the program pick
KSIEVE_AUTO=1 ./ksieve 941635629865253
```

### Tuning environment variables

Most users won't need these. They're provided for benchmarking, debugging, and reproducing paper results.

| Variable | Values | Default | Purpose |
|---|---|---|---|
| `KSIEVE_NP` | integer 4–63 | auto | Override the number of CRT primes. The default is chosen by an empirically validated heuristic based on bit size; override only for sweeps. |
| `KSIEVE_BS_NP` | integer 3–33 | mode-aware | Override the number of bitset primes. Defaults are `crt_np + 2` (K-sieve), `crt_np + 3` (S-sieve), `crt_np + 5` (M-sieve). |
| `KSIEVE_CRT_PRIMES` | comma-separated primes | auto | Explicit list of CRT primes (experimental). Must have exactly `np` entries and all must be odd and > 3. Example: `KSIEVE_CRT_PRIMES=5,7,11,13,17,19,23,29`. |
| `KSIEVE_ADAPTIVE_PRIMES` | 0 or 1 | 0 | Enable adaptive last-prime selection at sieve build time. Mixed empirical results (sometimes 60% faster, sometimes 43% slower); off by default. |
| `KSIEVE_W_SORT` | 0 or 1 | 0 | Enable W-sort DFS branch ordering (Jacobi=+1 residues before Jacobi=0). ~1.18× geometric-mean speedup at 120-bit. |

### Example: reproducing the Rev. 5 bitset-tuning sweep

```bash
# Build once
gcc -O3 -march=native -DKSIEVE_N_THREADS=16 -o ksieve main.c -lgmp -lm -lpthread

# Sweep BS_NP values at 90-bit M-sieve
for np in 9 10 11 12 13 14; do
    mean=$(KSIEVE_MSIEVE=1 KSIEVE_BS_NP=$np ./ksieve --bench 90 20 \
           | awk '/Mean total/ {print $4}')
    echo "BS_NP=$np  mean_total=$mean s"
done
```

Expected output: `BS_NP=13` should be fastest (≈0.11 s).

### Example: instrumented run

```bash
gcc -O3 -march=native -DKSIEVE_N_THREADS=16 -DINSTRUMENT \
    -o ksieve_instr main.c -lgmp -lm -lpthread

KSIEVE_MSIEVE=1 ./ksieve_instr --bench-instr 100 5
```

Prints per-filter rejection counters (`qr31=…, qr37=…, cp=…`) after the timing summary — useful for understanding where work is being done and for verify-chain reordering.

---

## Output format

A typical `--bench` invocation produces:

```
#     N-bit   np    crt_ops    k0s        isqrt    p1(s)    total(s)  OK?
----  ------  ----  ---------  ---------  -------  -------  --------  ---
1     100     8     3124153    2180576    1        0.4392   0.4400    OK
      N=27949867069715778095124598391  u=136707860879549  v=204449597045059
...

--- Summary (10 OK / 0 FAIL / 0 TIMEOUT / 10 total) ---
  Mean p1    : 2.3179 s
  Mean total : 2.3185 s
  Min / Max  : 0.3003 / 5.4795 s
  Mean CRT   : 3057373
  Mean K0s   : 2271752
  Mean isqrt : 1
```

Column glossary:
- **crt_ops** — total CRT residue-combination operations performed
- **k0s** — total K₀ values produced by the CRT DFS (each feeds the bitset scan)
- **isqrt** — integer square-root calls; should be 1 for a successful sieve run (meaning the true factor was found on the first isqrt)
- **p1** — Phase-1 wall time (the sieve itself)
- **total** — Phase-1 + Phase-2 (factor recovery and verification)

`isqrt = 1` across all trials is the correctness signal: the filtering pipeline found the answer without producing false positives.

---

## Bit-size expectations

## Bit-size expectations

Mean wall-time on a 16-thread build across 30 trials per cell (Rev. 6 benchmarks):

| N bits | K-sieve | S-sieve | M-sieve |
|---|---|---|---|
| 80  | 0.012 s  | 0.002 s  | 0.002 s  |
| 90  | 0.248 s  | 0.012 s  | 0.021 s  |
| 100 | 4.781 s  | 0.237 s  | 0.455 s  |
| 110 | —        | 2.724 s  | 9.803 s  |
| 120 | —        | 54.0 s   | —        |

Variance across random semiprimes at fixed bit size is high (10–60× min-to-max spread at 100-bit+, driven by factor-balance and DFS ordering effects). Means stabilize with 20–30 trials. At 100-bit, the S-sieve's single-trial range was 0.016 s to 0.794 s; the M-sieve's was 0.022 s to 1.168 s.

Single-thread performance at 100-bit: M-sieve mean 4.83 s (7.5× slower than 16-thread). Multi-threaded speedup is sub-linear due to memory-bandwidth saturation in the bitset scan and early-termination tail effects; see the whitepaper §5.5.

---

## What it cannot do

- Factor N with a factor ≤ 3 (trivial; use GCD).
- Factor N with more than 2 prime factors.
- Factor N where both factors are within a known small difference (use Fermat directly; this sieve is slower for that case).
- Beat the quadratic sieve or number field sieve on arbitrary numbers. This is a special-purpose algorithm for N ≡ ±1 (mod 6) semiprimes.

At 130+ bits, general-purpose methods (msieve, YAFU, CADO-NFS) start outperforming this implementation.

---

## Repository layout

```
main.c                      Single-file implementation (~3500 lines)
gather_benchmarks.sh        Overnight benchmark harness (outputs CSV + markdown tables)
ksieve_whitepaper_v5.docx   Full mathematical and engineering writeup
README.md                   This file
```

---

## License

MIT License. See `LICENSE`.

---

## Citation

If you use this work, please cite:

> Embry, A. *K-Sieve / S-Sieve / M-Sieve: A Specialized Factoring Algorithm for 6n±1 Semiprimes*, Revision 5, April 2026.

And the foundational paper on which the core sieve structure is based:

> Hittmeir, M. *Integer factorization as subset-sum problem*. Journal of Number Theory **249**, 93–118 (2023). [arXiv:2205.10074](https://arxiv.org/abs/2205.10074).
