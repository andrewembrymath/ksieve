# ksieve

A specialized integer factoring program for semiprimes `N = u·v` where both prime factors are small multiples of 6 ± 1. Three complementary sieves:

- **K-sieve** and **S-sieve** — factor N ≡ 1 (mod 6), where both factors satisfy p ≡ −1 (mod 6)
- **M-sieve** — factors N ≡ 5 (mod 6), where one factor is ≡ −1 (mod 6) and the other is ≡ +1 (mod 6)

Together these cover every prime-pair semiprime with factors greater than 3 **except** the case where both factors are ≡ 1 (mod 6). That remaining case (~25% of random RSA-style semiprimes) is not handled by this algorithm family.

On a 16-thread build, 100-bit semiprimes factor in roughly 0.24 s (S-sieve) to 0.61 s (M-sieve) mean wall-time; 120-bit in roughly 12.7 s mean (S-sieve, with multi-pass capped-range scanning); 130-bit in roughly 205 s mean.

The sieves share the core mathematical structure of Hittmeir's hyperbolic sieve for Fermat factorization ([arXiv:2205.10074](https://arxiv.org/abs/2205.10074)) — filter candidates by testing whether a specific quadratic polynomial is a residue mod a composite modulus assembled from small primes via CRT — but apply it to the `6n±1` parameterization and add Hensel lifting at p²/p³/p⁴, a tiled 64-bit-word bitset exploiting 2-adic structure, AVX-512 batch-16 acceleration, multi-threaded prefix work-stealing, and **capped-range multi-pass scanning** that exploits the skewed distribution of `K_true/K_max` to terminate early on the common case. See the whitepaper for the full mathematical treatment and benchmarks.

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

# Benchmark 10 random 100-bit RSA-style semiprimes (auto-selects sieve per N)
KSIEVE_RSA=1 ./ksieve --bench 100 10

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

These select which sieve to run (when factoring a specific N) or which semiprime form to generate (when benchmarking). At most one should be set; precedence order is **RSA > MSIEVE > SSIEVE > AUTO > default**.

| Variable | Effect |
|---|---|
| `KSIEVE_RSA=1` | **Benchmark only.** Generate RSA-style semiprimes (balanced primes with no mod-6 constraint). Automatically selects S-sieve or M-sieve per trial based on `N mod 6`. Rejects the ~25% of draws where both primes are ≡ 1 (mod 6), since no sieve in this family handles that case. |
| `KSIEVE_MSIEVE=1` | Use M-sieve (for N ≡ 5 mod 6, new form with mixed-class factors). |
| `KSIEVE_SSIEVE=1` | Use S-sieve (for N ≡ 1 mod 6, old form, ~4.2× smaller scan range than K-sieve). |
| `KSIEVE_AUTO=1` | If no explicit mode flag is set: auto-select based on `N mod 6` (M-sieve for N ≡ 5, S-sieve for N ≡ 1). |
| _(none)_ | Default: K-sieve. |

Mode examples:

```bash
# New-form semiprime (u ≡ 5, v ≡ 1 mod 6)
KSIEVE_MSIEVE=1 ./ksieve --bench 100 10

# Old-form semiprime, faster S-sieve
KSIEVE_SSIEVE=1 ./ksieve --bench 80 10

# RSA-style: no mod-6 constraint on factors, mixed sieve usage across trials
KSIEVE_RSA=1 ./ksieve --bench 100 10

# Factor an unknown N, let the program pick
KSIEVE_AUTO=1 ./ksieve 941635629865253
```

#### Note on `KSIEVE_RSA` benchmarks

RSA-style generation draws balanced random primes with no `p mod 6` constraint. Of the four possible pairings of factor-residues mod 6:

| (u mod 6, v mod 6) | Resulting N mod 6 | Sieve used |
|---|---|---|
| (5, 5) | N ≡ 1 | S-sieve |
| (1, 5) or (5, 1) | N ≡ 5 | M-sieve |
| (1, 1) | N ≡ 1 | **rejected — not solvable by any sieve here** |

With uniformly random primes the (1, 1) case arises roughly 25% of the time. The benchmark prints the actual rejection count before the trial table — for example:

```
[RSA mode] generated 10 usable semiprimes; rejected 3 draws where both primes ≡ 1 mod 6 (23.1% of 13 total draws)
```

Factor difference is constrained to `|u - v| > 2^(bits/4)`, matching standard RSA keygen guidance and keeping trivial Fermat factorization from solving the instance.

### Capped-range multi-pass scanning

The algorithm exploits the empirical observation that `K_true / K_max` (the factor-imbalance ratio) is heavily concentrated near zero for random semiprimes — at 130-bit, roughly 73% of trials have `K_true/K_max < 5%` and ~100% are below 15%. Rather than scanning the full `[0, K_max]` range every time, the algorithm runs successive passes with progressively larger word-index caps, terminating as soon as any pass finds the factor.

At 130-bit S-sieve, this typically gives a **4–5× speedup** over the uncapped single-pass scan with no loss of correctness (every word is eventually scanned across the pass sequence).

| Variable | Default | Effect |
|---|---|---|
| `KSIEVE_PASS_FRACS` | unset | Comma-separated list of fraction cutoffs (of `n_words`). Each fraction becomes the upper bound for one pass; a final pass to 100% is always appended. Example: `KSIEVE_PASS_FRACS="0.03,0.08,0.15,0.35"` gives 5 passes at (0–3%, 3–8%, 8–15%, 15–35%, 35–100%). |
| `KSIEVE_PASS0_FRAC` | 0.10 | First-pass cutoff (legacy 3-pass knob, used if `KSIEVE_PASS_FRACS` is unset). |
| `KSIEVE_PASS1_FRAC` | 0.32 | Second-pass cutoff (legacy 3-pass knob, used if `KSIEVE_PASS_FRACS` is unset). |
| `KSIEVE_SINGLEPASS=1` | unset | Disable multi-pass — force a single full-range scan. Useful for regression testing against older builds. |
| `KSIEVE_PASS_TRACE=1` | unset | Print the computed pass bounds to stderr before running. |

Input to `KSIEVE_PASS_FRACS` is robust: values are auto-sorted, deduplicated, and invalid entries (≤0, ≥1, garbage) are silently ignored. Up to 10 passes are supported. Examples:

```bash
# Legacy 3-pass (default)
./ksieve --bench 130 10

# Single aggressive cap: 2 passes at (0–15%, 15–100%).
# Reasonable if your observed K_true/K_max distribution is tight.
KSIEVE_PASS_FRACS="0.15" KSIEVE_SSIEVE=1 ./ksieve --bench 130 10

# Fine-grained 5-pass: better for bit ranges where the distribution
# has long tails beyond the first cutoff.
KSIEVE_PASS_FRACS="0.03,0.08,0.15,0.35" KSIEVE_SSIEVE=1 ./ksieve --bench 130 10

# Regression: single-pass (reproduces pre-multi-pass timings)
KSIEVE_SINGLEPASS=1 KSIEVE_SSIEVE=1 ./ksieve --bench 130 10

# See what the algorithm actually does with your fractions
KSIEVE_PASS_TRACE=1 KSIEVE_PASS_FRACS="0.05,0.15,0.35" ./ksieve --bench 130 3
# Prints to stderr, e.g.:
#   [ksieve] n_words=1024, 4 pass(es): 0 51 153 358 1024
```

#### Tuning caveats

Multi-pass is a **variance-reducing optimization**, not a worst-case improvement. Trials with `K_true/K_max` just above a pass cutoff pay the full cost of the preceding (empty) pass plus the work to find `K_true` in the next pass — occasionally worse than single-pass for that specific instance. Mean and median times improve substantially; pathological cases can regress. Choose pass fractions that don't cluster near the likely `K_true/K_max` values for your bit range.

### Other tuning environment variables

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

### Example: tuning pass fractions

```bash
# Baseline: single-pass (no capped-range optimization)
KSIEVE_SINGLEPASS=1 KSIEVE_SSIEVE=1 ./ksieve --bench 130 30

# Default 3-pass
KSIEVE_SSIEVE=1 ./ksieve --bench 130 30

# Sweep different pass configurations
for fracs in "0.15" "0.10,0.32" "0.05,0.15,0.35" "0.03,0.08,0.15,0.35"; do
    mean=$(KSIEVE_PASS_FRACS="$fracs" KSIEVE_SSIEVE=1 \
           ./ksieve --bench 130 10 | awk '/Mean total/ {print $4}')
    echo "PASS_FRACS='$fracs'  mean_total=$mean s"
done
```

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
- **crt_ops** — total CRT residue-combination operations performed (summed across all passes when multi-pass is active)
- **k0s** — total K₀ values produced by the CRT DFS (also summed across passes)
- **isqrt** — integer square-root calls; should be 1 for a successful sieve run (meaning the true factor was found on the first isqrt)
- **p1** — Phase-1 wall time (the sieve itself)
- **total** — Phase-1 + Phase-2 (factor recovery and verification)

`isqrt = 1` across all trials is the correctness signal: the filtering pipeline found the answer without producing false positives.

---

## Bit-size expectations

Mean wall-time on a 16-thread build across 30 trials per cell (with capped-range multi-pass at default fractions):

| N bits | K-sieve | S-sieve | M-sieve |
|---|---|---|---|
| 80  | 0.012 s  | 0.002 s  | 0.002 s  |
| 90  | 0.248 s  | 0.012 s  | 0.021 s  |
| 100 | 4.781 s  | 0.237 s  | 0.455 s  |
| 110 | —        | 2.724 s  | 9.803 s  |
| 120 | —        | 12.7 s   | —        |
| 130 | —        | 205 s    | —        |

The 120-bit and 130-bit S-sieve numbers reflect the capped-range multi-pass scanning optimization (~4.3× speedup over single-pass at those sizes). Use `KSIEVE_SINGLEPASS=1` to reproduce pre-multi-pass timings (~54 s at 120-bit, ~895 s at 130-bit).

Variance across random semiprimes at fixed bit size is high (10–60× min-to-max spread at 100-bit+, driven by factor-balance and DFS ordering effects — multi-pass reduces this but does not eliminate it). Means stabilize with 20–30 trials.

Single-thread performance at 100-bit: M-sieve mean 4.83 s (7.5× slower than 16-thread). Multi-threaded speedup is sub-linear due to memory-bandwidth saturation in the bitset scan and early-termination tail effects; see the whitepaper §5.5.

---

## What it cannot do

- Factor N with a factor ≤ 3 (trivial; use GCD).
- Factor N with more than 2 prime factors.
- Factor N where **both** factors are ≡ 1 (mod 6). About 25% of arbitrary RSA-style semiprimes fall in this category. `KSIEVE_RSA=1` benchmark mode detects and rejects these at generation time.
- Factor N where both factors are within a known small difference (use Fermat directly; this sieve is slower for that case).
- Beat the quadratic sieve or number field sieve on arbitrary numbers. This is a special-purpose algorithm for N ≡ ±1 (mod 6) semiprimes where factor-class mod 6 is compatible with one of the three sieves.

At 130+ bits, general-purpose methods (msieve, YAFU, CADO-NFS) start outperforming this implementation.

---

## Repository layout

```
main.c                      Single-file implementation (~3700 lines)
gather_benchmarks.sh        Overnight benchmark harness (outputs CSV + markdown tables)
ksieve_whitepaper_v6.docx   Full mathematical and engineering writeup
README.md                   This file
```

---

## License

MIT License. See `LICENSE`.

---

## Citation

If you use this work, please cite:

> Embry, A. *K-Sieve / S-Sieve / M-Sieve: A Specialized Factoring Algorithm for 6n±1 Semiprimes*, Revision 6, April 2026.

And the foundational paper on which the core sieve structure is based:

> Hittmeir, M. *Integer factorization as subset-sum problem*. Journal of Number Theory **249**, 93–118 (2023). [arXiv:2205.10074](https://arxiv.org/abs/2205.10074).
