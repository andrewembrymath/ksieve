# ksieve

A specialized integer factoring program for semiprimes `N = u·v` where both prime factors are small multiples of 6 ± 1. Four complementary sieves:

- **K-sieve** — general-purpose, factors any N ≡ 1 (mod 6). Class-agnostic polynomial but widest scan range; slowest of the four.
- **S-sieve** — factors N ≡ 1 (mod 6) where both factors are ≡ −1 (mod 6). ~4.24× smaller scan range than K-sieve plus mod-3 root fold.
- **T-sieve** — factors N ≡ 1 (mod 6) where both factors are ≡ +1 (mod 6). Mirror of the S-sieve; identical performance profile on its input class. *Added in Rev. 8.*
- **M-sieve** — factors N ≡ 5 (mod 6), where one factor is ≡ −1 (mod 6) and the other is ≡ +1 (mod 6).

Together these cover **every** prime-pair semiprime with factors greater than 3. Rev. 8's T-sieve closed the coverage gap — previous revisions did not have a fast sieve for the (+1, +1) factor-pairing class (though the K-sieve always handled it, slowly).

On a 16-thread build, 100-bit semiprimes factor in roughly 0.044–0.059 s mean wall-time (S- and T-sieves) to 0.10 s (M-sieve) to 0.76 s (K-sieve); 120-bit S- and T-sieves each factor in roughly 13.6 s mean (with capped-range multi-pass scanning).

The sieves share the core mathematical structure of Hittmeir's hyperbolic sieve for Fermat factorization ([arXiv:2205.10074](https://arxiv.org/abs/2205.10074)) — filter candidates by testing whether a specific quadratic polynomial is a residue mod a composite modulus assembled from small primes via CRT — but apply it to the `6n±1` parameterization and add Hensel lifting at p²/p³/p⁴, a tiled 64-bit-word bitset exploiting 2-adic structure, AVX-512 batch-16 acceleration, multi-threaded prefix work-stealing, and **capped-range multi-pass scanning** that exploits the skewed distribution of `K_true/K_max` to terminate early on the common case. See the whitepaper for the full mathematical treatment and benchmarks.

---

## Quickstart

```bash
# Install dependencies (Debian/Ubuntu)
sudo apt install build-essential libgmp-dev

# Build a 16-thread release binary
gcc -O3 -march=native -DKSIEVE_N_THREADS=16 -o ksieve main.c -lgmp -lm -lpthread

# Factor a specific N (auto-detects sieve based on N mod 6; tries S- then T-sieve for N ≡ 1)
# This N has factors (5 mod 6) × (1 mod 6) → N ≡ 5 mod 6 → routed to M-sieve directly.
KSIEVE_AUTO=1 ./ksieve 941635629865253

# This N has factors (1 mod 6) × (1 mod 6) → N ≡ 1 mod 6 → S-sieve first, then T-sieve on exhaustion.
KSIEVE_AUTO=1 ./ksieve 665040100429

# Benchmark 10 random 100-bit new-form (N ≡ 5 mod 6) semiprimes
KSIEVE_MSIEVE=1 ./ksieve --bench 100 10

# Benchmark 10 random 100-bit (+1,+1)-form (N ≡ 1 mod 6, factors ≡ +1 mod 6) semiprimes
KSIEVE_TSIEVE=1 ./ksieve --bench 100 10

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

These select which sieve to run (when factoring a specific N) or which semiprime form to generate (when benchmarking). At most one should be set; precedence order is **RSA > MSIEVE > SSIEVE > TSIEVE > AUTO > default**.

| Variable | Effect |
|---|---|
| `KSIEVE_RSA=1` | **Benchmark only.** Generate RSA-style semiprimes (balanced primes with no mod-6 constraint). Automatically selects S-, T-, or M-sieve per trial based on `N mod 6`; for N ≡ 1 (mod 6) it tries S-sieve first and falls back to T-sieve if S-sieve's scan exhausts (since the factor class can't be distinguished from N alone). **As of Rev. 8 all four factor-class pairings are handled** — no draws are rejected. |
| `KSIEVE_MSIEVE=1` | Use M-sieve (for N ≡ 5 mod 6, new form with mixed-class factors). |
| `KSIEVE_SSIEVE=1` | Use S-sieve (for N ≡ 1 mod 6 with both factors ≡ −1 mod 6, "old form", ~4.24× smaller scan range than K-sieve). |
| `KSIEVE_TSIEVE=1` | Use T-sieve (for N ≡ 1 mod 6 with both factors ≡ +1 mod 6, "(+1,+1) form"). Mirror of S-sieve with matched performance. *Rev. 8+.* |
| `KSIEVE_AUTO=1` | If no explicit mode flag is set: auto-select based on `N mod 6`. For N ≡ 5 → M-sieve. For N ≡ 1 → try S-sieve first; if its scan exhausts (i.e. the factors are actually in the (+1,+1) class), retry with T-sieve. Worst-case cost ≈2× S-sieve time; expected cost ≈1.5× on uniform random input. |
| _(none)_ | Default: K-sieve (class-agnostic, handles both (−1,−1) and (+1,+1) old-form inputs transparently at the cost of its wider scan range). |

Mode examples:

```bash
# New-form semiprime (u ≡ 5, v ≡ 1 mod 6)
KSIEVE_MSIEVE=1 ./ksieve --bench 100 10

# Old-form semiprime, factors ≡ −1 mod 6, faster S-sieve
KSIEVE_SSIEVE=1 ./ksieve --bench 80 10

# (+1,+1)-form semiprime, factors ≡ +1 mod 6, faster T-sieve
KSIEVE_TSIEVE=1 ./ksieve --bench 80 10

# RSA-style: no mod-6 constraint on factors, mixed sieve usage across trials
KSIEVE_RSA=1 ./ksieve --bench 100 10

# Factor an unknown N, let the program pick
KSIEVE_AUTO=1 ./ksieve 941635629865253
```

#### Note on `KSIEVE_RSA` benchmarks

RSA-style generation draws balanced random primes with no `p mod 6` constraint. All four possible pairings of factor-residues mod 6 are solvable as of Rev. 8:

| (u mod 6, v mod 6) | Resulting N mod 6 | Sieve used |
|---|---|---|
| (5, 5) | N ≡ 1 | S-sieve |
| (1, 1) | N ≡ 1 | T-sieve (via S-first-T-fallback) |
| (1, 5) or (5, 1) | N ≡ 5 | M-sieve |

For N ≡ 1 (mod 6), the factor class cannot be determined from N alone, so auto-dispatch tries S-sieve first and retries with T-sieve on exhaustion. This typically adds ≈0.5× S-sieve time on average (expected cost ≈1.5× when the factor distribution is uniform over the two sub-classes).

A typical RSA-mode output preamble now looks like:

```
[RSA mode] all 10 draws usable (Rev. 8: T-sieve handles (+1,+1) class; the AUTO/RSA dispatcher retries S→T on N ≡ 1 mod 6)
```

Pre-Rev. 8 users: earlier builds rejected the (1,1) case at generation time and printed a rejection-count line. That behavior is gone — no draws are rejected anymore. If an older codebase somehow exposes the legacy rejection counter (it stays 0 in Rev. 8+), the line prints as `[RSA mode] rejected 0 / N draws` and can be ignored.

**Priority-override behavior.** `KSIEVE_RSA=1` takes precedence over `KSIEVE_SSIEVE=1`, `KSIEVE_TSIEVE=1`, and `KSIEVE_MSIEVE=1`. This is deliberate: RSA generation produces a mix of N ≡ 1 and N ≡ 5 (mod 6) values with varying factor classes, and forcing a single sieve across all of them causes silent correctness failures on the mismatched subsets. If you set both `KSIEVE_RSA=1` and an explicit sieve flag, a one-line warning prints to stderr and the explicit flag is ignored.

Factor difference is constrained to `|u - v| > 2^(bits/4)`, matching standard RSA keygen guidance and keeping trivial Fermat factorization from solving the instance.

### Capped-range multi-pass scanning

The algorithm exploits the empirical observation that `K_true / K_max` (the factor-imbalance ratio) is heavily concentrated near zero for random semiprimes — at 130-bit, roughly 73% of trials have `K_true/K_max < 5%` and ~100% are below 15%. Rather than scanning the full `[0, K_max]` range every time, the algorithm runs successive passes with progressively larger word-index caps, terminating as soon as any pass finds the factor.

At 130-bit S-sieve, this typically gives a **4–5× speedup** over the uncapped single-pass scan with no loss of correctness (every word is eventually scanned across the pass sequence).

| Variable | Default | Effect |
|---|---|---|
| `KSIEVE_PASS_FRACS` | `"0.03,0.08,0.15,0.35"` | Comma-separated list of fraction cutoffs (of `n_words`). Each fraction becomes the upper bound for one pass; a final pass to 100% is always appended. The Rev. 8 default `"0.03,0.08,0.15,0.35"` gives 5 passes at (0–3%, 3–8%, 8–15%, 15–35%, 35–100%). Set explicitly to override. |
| `KSIEVE_PASS0_FRAC` | 0.10 | Legacy 2-knob 3-pass knob. If either `KSIEVE_PASS0_FRAC` or `KSIEVE_PASS1_FRAC` is set, the legacy 3-pass behaviour is used (the unset one defaults to its historical value) and the Rev. 8 5-pass default is overridden. Kept for benchmark reproduction against pre-Rev. 8 builds. |
| `KSIEVE_PASS1_FRAC` | 0.32 | Legacy 2-knob 3-pass knob. See `KSIEVE_PASS0_FRAC`. |
| `KSIEVE_SINGLEPASS=1` | unset | Disable multi-pass — force a single full-range scan. Useful for regression testing against older builds. |
| `KSIEVE_PASS_TRACE=1` | unset | Print the computed pass bounds to stderr before running. Also printed under verbose mode (single-factor invocations via `./ksieve <N>`). |

Input to `KSIEVE_PASS_FRACS` is robust: values are auto-sorted, deduplicated, and invalid entries (≤0, ≥1, garbage) are silently ignored. Up to 10 passes are supported. Examples:

```bash
# Default: 5-pass fine-grained split (Rev. 8). No env var needed.
./ksieve --bench 130 10

# Single aggressive cap: 2 passes at (0–15%, 15–100%).
# Reasonable if your observed K_true/K_max distribution is tight.
KSIEVE_PASS_FRACS="0.15" KSIEVE_SSIEVE=1 ./ksieve --bench 130 10

# Reproduce the Rev. 6/7 legacy 3-pass default (for regression comparison)
KSIEVE_PASS0_FRAC=0.10 KSIEVE_SSIEVE=1 ./ksieve --bench 130 10

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
| `KSIEVE_BS_NP` | integer 3–33 | mode-aware | Override the number of bitset primes. Defaults are `crt_np + 2` (K-sieve), `crt_np + 3` (S- and T-sieves), `crt_np + 5` (M-sieve). |
| `KSIEVE_CRT_PRIMES` | comma-separated primes | auto | Explicit list of CRT primes (experimental). Must have exactly `np` entries and all must be odd and > 3. Example: `KSIEVE_CRT_PRIMES=5,7,11,13,17,19,23,29`. |
| `KSIEVE_ADAPTIVE_PRIMES` | 0 or 1 | 0 | Enable adaptive last-prime selection at sieve build time. Mixed empirical results (sometimes 60% faster, sometimes 43% slower); off by default. |
| `KSIEVE_W_SORT` | 0 or 1 | 0 | Enable W-sort DFS branch ordering (Jacobi=+1 residues before Jacobi=0). ~1.18× geometric-mean speedup at 120-bit. |

#### Internal flags (not for direct use)

`KSIEVE_RETRY_AS_T` is set internally by the `KSIEVE_AUTO=1` / `KSIEVE_RSA=1` dispatch wrapper when the initial S-sieve pass on an N ≡ 1 (mod 6) input exhausts without finding a factor, signalling the fallback to T-sieve. Users should not set this manually; doing so on a specific-N invocation would force T-sieve mode and skip the S-first pass (equivalent to `KSIEVE_TSIEVE=1` with slightly different priority semantics). If you see `KSIEVE_RETRY_AS_T=1` in process environment during an auto-dispatch retry — for instance via `strace` or `ltrace` — that's the wrapper doing its job.

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

Expected output: `BS_NP ∈ {11, 12, 13}` cluster within trial-to-trial variance. Rev. 5/6 measured BS_NP=13 fastest at ≈0.11 s; Rev. 8 on different hardware measured BS_NP=12 fastest at ≈0.012 s. The `crt_np + 5` default (BS_NP=14 at 90-bit) is within ~30% of optimal across hardware configurations. For production benchmarking, sweep on your target hardware.

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

# Default 5-pass (Rev. 8 default: 0.03,0.08,0.15,0.35)
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

Mean wall-time on a 16-thread build across 30 trials per cell (20 trials at 120-bit), Rev. 8 benchmark data with the default 5-pass multi-pass configuration:

| N bits | K-sieve | S-sieve | T-sieve | M-sieve |
|---|---|---|---|---|
| 60  | —        | —        | —        | 0.001 s  |
| 70  | —        | —        | —        | 0.001 s  |
| 80  | 0.003 s  | 0.001 s  | 0.001 s  | 0.001 s  |
| 90  | 0.046 s  | 0.003 s  | 0.003 s  | 0.010 s  |
| 100 | 0.863 s  | 0.059 s  | 0.044 s  | 0.104 s  |
| 110 | —        | 0.623 s  | 0.640 s  | 2.911 s  |
| 120 | —        | 13.586 s | 13.806 s | —        |

Dashes denote bit-sizes not measured in the Rev. 8 sweep — typically because the sieve runs well beyond practical wall-time at that size (K-sieve above 100-bit, M-sieve above 110-bit). S- and T-sieve have indistinguishable performance at matched bit size (within trial-to-trial variance at ≤3%), as predicted algebraically: the two sieves have identical search-range width, root-modulus density, discriminant, and QR144 tables; only the polynomial's linear-term sign differs.

Variance across random semiprimes at fixed bit size is high (10–90× min-to-max spread at 110–120 bit, driven by factor-balance and DFS ordering effects — multi-pass reduces this but does not eliminate it). Means stabilize with 20–30 trials.

Single-thread performance at 100-bit M-sieve: 1.157 s mean (5.3× slower than 16-thread, Rev. 8 10-trial run; Rev. 6 measured 7.5× on different hardware). Multi-threaded speedup is sub-linear due to memory-bandwidth saturation in the bitset scan and early-termination tail effects; see the whitepaper §6.6.

Scaling exponents from linear regression on log₂(mean total time) vs bit size over 80–120 bit: S-sieve 2^(0.355·bits), T-sieve 2^(0.356·bits) — indistinguishable — M-sieve 2^(0.367·bits), K-sieve 2^(0.402·bits) (over 80–100 bit).

---

## What it cannot do

- Factor N with a factor ≤ 3 (trivial; use GCD).
- Factor N with more than 2 prime factors.
- Factor N where both factors are within a known small difference (use Fermat directly; this sieve is slower for that case).
- Factor N whose factors fall outside the 6n±1 family (e.g. N with a factor equal to 2 or 3; handled by trial division before invoking ksieve).
- Beat the quadratic sieve or number field sieve on arbitrary numbers. This is a special-purpose algorithm for N whose factors are in the 6n±1 family and fall into one of the four residue-class pairings handled by the K/S/T/M-sieves.

As of Rev. 8, all four prime-pair residue-class pairings mod 6 are solvable: (−1,−1) by S-sieve, (+1,+1) by T-sieve, and the two mixed pairings by M-sieve. Previous revisions could not factor the (+1,+1) case quickly (only the K-sieve handled it, slowly); Rev. 8 closes that gap.

At 130+ bits, general-purpose methods (msieve, YAFU, CADO-NFS) start outperforming this implementation.

---

## Repository layout

```
main.c                      Single-file implementation (~4240 lines)
gather_benchmarks.sh        Overnight benchmark harness (outputs CSV + markdown tables)
ksieve_whitepaper_v8.docx   Full mathematical and engineering writeup (Rev. 8)
README.md                   This file
```

---

## License

MIT License. See `LICENSE`.

---

## Citation

If you use this work, please cite:

> Embry, A. *K-Sieve / S-Sieve / T-Sieve / M-Sieve: A Specialized Factoring Algorithm for 6n±1 Semiprimes, with Multi-Layer Quadratic-Residue Filtering, Sum-Variable Range Reduction, Complete (6x±1)-Class Coverage, and AVX-512 Acceleration*, Revision 8, April 2026.

And the foundational paper on which the core sieve structure is based:

> Hittmeir, M. *Integer factorization as subset-sum problem*. Journal of Number Theory **249**, 93–118 (2023). [arXiv:2205.10074](https://arxiv.org/abs/2205.10074).

The 6n±1 parameterization underlying all four sieves is due to:

> Creft, R. *Hexile Sieve Analysis of Prime and Composite Integers*, [arXiv:1202.5948](https://arxiv.org/abs/1202.5948) [math.NT], February 2012.
