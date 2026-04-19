/*
 * ksieve.c  —  K-Sieve integer factorization
 *
 * Factors N = u·v where u ≡ v ≡ 5 mod 6  (primes of form 6k−1).
 *
 * IDENTITY:  K = (v−u)/6,  T = (u+v)/2,  9K² + N = T²
 * GOAL:      Find K ∈ [1, K_max = floor(sqrt(N)/3)] such that 9K²+N is a perfect square.
 *
 * ── ALGORITHM ─────────────────────────────────────────────────────────────
 *
 * For each sieve prime p, the valid residue set is:
 *   valid_p = {k mod p : Jacobi(9k²+N, p) >= 0}   (|valid_p| ~= p/2)
 *
 * PHASE 1 — CRT DFS enumeration:
 *   Iterative depth-first search through the CRT residue tree.
 *   Primes sorted by density (nvp/p) ascending: best filter first.
 *   Leaves give K0 values satisfying all Jacobi filters simultaneously.
 *   Uses __int128 arithmetic so M = prod(p) can exceed 2^64 (np >= 15 safe).
 *   Prune: if K_new > K_max, skip entire subtree.
 *
 * PHASE 2 — Stride scan + isqrt test:
 *   For each K0: scan K = K0, K0+M, K0+2M, ... <= K_max.
 *   Hard filters (parity, mod-8) skip ~81% before isqrt.
 *   Native u128 path for N<=64-bit (~8 ns/test); GMP for N>64-bit (~300 ns/test).
 *
 * Hard filters (0% false-negative rate):
 *   Parity: (K & 1) == k_parity  where k_parity = (N%4==1 ? 0 : 1)
 *   Mod-8:  (n^2 - K^2) == 0 mod 8  where n = (N-1)/6
 *
 * ── RANGE ─────────────────────────────────────────────────────────────────
 *   Supported: N up to ~100 bits, both factors must be == 5 mod 6.
 *   For N > 64-bit phase 2 uses GMP isqrt (~300 ns/test).
 *
 * ── BUILD ─────────────────────────────────────────────────────────────────
 *   gcc -O3 -march=native -o ksieve ksieve.c -lgmp -lm
 *
 * ── USAGE ─────────────────────────────────────────────────────────────────
 *   ./ksieve <N>                    Factor decimal N
 *   ./ksieve --demo                 Run built-in examples (40-100 bit)
 *   ./ksieve --bench <bits> <n>     Benchmark n random b-bit semiprimes
 *   ./ksieve --info                 Show performance profile table
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <gmp.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <signal.h>
#include <unistd.h>
#include <pthread.h>
#include <stdatomic.h>
#ifdef __AVX2__
#include <immintrin.h>
#endif

/* Threading config. N_THREADS=1 gives single-threaded (no threading overhead). */
#ifndef KSIEVE_N_THREADS
#define KSIEVE_N_THREADS 2
#endif
#define SPLIT_DEPTH_MAX 4    /* max DFS depth at which we split into work items */

/* ── Integer types ────────────────────────────────────────────────────── */
typedef unsigned __int128   u128;
typedef unsigned long long  u64;
typedef __int128            i128;
typedef long long           i64;

/* ── Sieve mode ─────────────────────────────────────────────────────────
 * MODE_K_SIEVE: old form (6x-1)(6y-1)=6n+1, scan K=y-x.
 *   g_K(K) = 9K² + N = T² where T=3(x+y)-1.
 *
 * MODE_S_SIEVE: old form (6x-1)(6y-1)=6n+1, scan S=x+y.
 *   g_S(S) = 9S² - 6S - 6n = (3K)². Range 4.24× smaller than K-sieve.
 *   Root mod folds mod-3 (S ≡ -n mod 3 from 6xy=n+S), giving +1.585 bits.
 *
 * MODE_M_SIEVE: NEW form (6x-1)(6y+1)=6n-1, scan M=x+y.
 *   g_M(M) = 9M² - N = L² where L = 3K_Δ+1, K_Δ = y-x.
 *   Range same as S-sieve. NO mod-3 fold available (6xy = n+K_Δ folds to K_Δ, not M).
 *   Target: N ≡ 5 (mod 6) semiprimes (which K/S-sieve cannot factor).
 *
 * The enum replaces the earlier boolean `s_sieve` field. Old code "if (s_sieve)"
 * becomes either "if (mode != MODE_K_SIEVE)" (for sum-var behavior shared by S and M)
 * or "if (mode == MODE_M_SIEVE)" (for the new-form-specific paths). */
enum sieve_mode {
    MODE_K_SIEVE = 0,
    MODE_S_SIEVE = 1,
    MODE_M_SIEVE = 2,
};

/* ── Sieve prime table (starting at 5; 2 and 3 excluded by 6k-1 structure) */
static const int SP[] = {
     5,  7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
    67, 71, 73, 79, 83, 89, 97,101,103,107,109,113,127,131,137,139,
   149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,
   233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317
};
#define NSP      (int)(sizeof(SP)/sizeof(SP[0]))
#define MAX_NP   64
#define MAX_NVP  320   /* nvp <= ceil(p/2)+1; p <= 317 */

/* ── Extra primes for post-DFS stride pre-filtering ─────────────────────
 * These primes are NOT in the CRT modulus M. After Phase 1 produces K0
 * residues, for each K0 we build a bitset over the stride indices j in
 * [0, K_max/M) marking which j values are consistent with all extra primes.
 * We then test only those j values with isqrt, skipping the rest.
 * This has 0% FNR: j_true is always included in the bitset.
 * For 80-bit N with np=8 CRT primes (M=2^30, ~256 strides):
 *   8 extra primes reduce survivors from 256 to ~1, saving ~80% of isqrt.   */
#define EXTRA_NP      8     /* max extra filter primes (prod fits in long)    */
#define EXTRA_MAXQ    80    /* largest extra prime supported (< 80)           */
#define EXTRA_MAXNVP  42    /* max valid K residues per prime (ceil(79/2)+1)  */

/* Extra prime data precomputed per N */
typedef struct {
    int     np;                            /* number of extra primes used      */
    int     q   [EXTRA_NP];               /* the primes                       */
    long    Minv[EXTRA_NP];               /* modinv(M, q[i])                  */
    long    inv_cache[EXTRA_NP];          /* modinv(prod(q[0..d-1]), q[d])    */
    int     nvp [EXTRA_NP];               /* count of valid K residues        */
    uint8_t valid_v[EXTRA_NP][EXTRA_MAXNVP]; /* packed list of valid residues */
} ExtraFilter;

/* ── Timing ───────────────────────────────────────────────────────────── */
static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* ── Cycle-accurate instrumentation ───────────────────────────────────────
 * Controlled by INSTRUMENT macro. When INSTRUMENT is defined at compile time
 * the hot paths emit rdtsc counters into global accumulators. When not
 * defined, all measurement macros are no-ops and the compiler deletes them.
 * We use rdtsc (not rdtscp) since we want unordered-timestamp reads; the
 * ordering is established by volatile memory barriers on the accumulators. */
#ifdef INSTRUMENT
static inline u64 rdtsc_now(void) {
    unsigned int lo, hi;
    __asm__ __volatile__("rdtsc" : "=a"(lo), "=d"(hi));
    return ((u64)hi << 32) | lo;
}

/* Accumulators for each measured region. Counts are in cycles.
 * _n counts number of times the region fired; _c accumulates cycles. */
typedef struct {
    u64 setup_c,       setup_n;      /* batch16: pre-loop setup */
    u64 word_fast_c,   word_fast_n;  /* batch16: per-word on fast path (nz==0) */
    u64 word_slow_c,   word_slow_n;  /* batch16: per-word with nz!=0 */
    u64 testK_c,       testK_n;      /* fused_test_K total (incl. isqrt) */
    u64 testK_fail_c,  testK_fail_n; /* fused_test_K rejections (no isqrt) */
    u64 crt_node_c,    crt_node_n;   /* per DFS node in fused loop */
    u64 total_c,       total_n;      /* ksieve_factor wall-cycles */
    u64 dead_lanes,    live_lanes;   /* batch16 lane stats (qr64 mask == 0 vs not) */
    u64 skipped_batches;              /* batches where all 16 lanes were dead */
    u64 survivors_found;              /* total K candidates yielded by bitset to testK */
    u64 np_bitset_sum;                /* sum of ctx->n_primes across batch16 calls */
    u64 reject_qr256, reject_qr65536, reject_qr144, reject_mod8;
    u64 reject_qr31, reject_qr37, reject_qr41, reject_qr49;
    u64 reject_qr27, reject_qr25, reject_qr121, reject_qr169;
    u64 reject_cp, reject_jac, reject_isqrt;
} InstrCounters;

static InstrCounters g_instr;

static void instr_reset(void) { memset(&g_instr, 0, sizeof(g_instr)); }

static void instr_print(FILE *f) {
    #define _PR(name, label) do { \
        double avg = g_instr.name##_n ? (double)g_instr.name##_c / g_instr.name##_n : 0; \
        fprintf(f, "  %-18s: n=%-12llu  cycles=%-15llu  avg=%.1f\n", \
            label, (unsigned long long)g_instr.name##_n, \
            (unsigned long long)g_instr.name##_c, avg); \
    } while (0)
    fprintf(f, "--- INSTRUMENT counters ---\n");
    _PR(setup,      "batch16_setup");
    _PR(word_fast,  "word_fast_path");
    _PR(word_slow,  "word_slow_path");
    _PR(testK_fail, "testK (reject)");
    _PR(testK,      "testK (all)");
    _PR(crt_node,   "crt_dfs_node");
    _PR(total,      "ksieve_total");
    #undef _PR
    fprintf(f, "  dead_lanes=%llu  live_lanes=%llu  dead_frac=%.3f  skipped_batches=%llu\n",
            (unsigned long long)g_instr.dead_lanes,
            (unsigned long long)g_instr.live_lanes,
            (double)g_instr.dead_lanes / (double)(g_instr.dead_lanes + g_instr.live_lanes + 1),
            (unsigned long long)g_instr.skipped_batches);
    if (g_instr.setup_n) {
        fprintf(f, "  survivors_total=%llu  avg_np_bitset=%.2f  survivors_per_slow_word=%.3f\n",
                (unsigned long long)g_instr.survivors_found,
                (double)g_instr.np_bitset_sum / g_instr.setup_n,
                (double)g_instr.survivors_found / (double)(g_instr.word_slow_n + 1));
    }
    fprintf(f, "  reject distribution (how many testK calls stopped at each filter):\n");
    fprintf(f, "    qr256=%llu qr65536=%llu qr144=%llu mod8=%llu\n",
        (unsigned long long)g_instr.reject_qr256, (unsigned long long)g_instr.reject_qr65536,
        (unsigned long long)g_instr.reject_qr144, (unsigned long long)g_instr.reject_mod8);
    fprintf(f, "    qr31=%llu qr37=%llu qr41=%llu qr49=%llu\n",
        (unsigned long long)g_instr.reject_qr31, (unsigned long long)g_instr.reject_qr37,
        (unsigned long long)g_instr.reject_qr41, (unsigned long long)g_instr.reject_qr49);
    fprintf(f, "    qr27=%llu qr25=%llu qr121=%llu qr169=%llu\n",
        (unsigned long long)g_instr.reject_qr27, (unsigned long long)g_instr.reject_qr25,
        (unsigned long long)g_instr.reject_qr121, (unsigned long long)g_instr.reject_qr169);
    fprintf(f, "    cp=%llu jac=%llu isqrt=%llu\n",
        (unsigned long long)g_instr.reject_cp, (unsigned long long)g_instr.reject_jac,
        (unsigned long long)g_instr.reject_isqrt);
    /* Derived percentages */
    if (g_instr.total_c) {
        double tot = (double)g_instr.total_c;
        double fast_pct = 100.0 * g_instr.word_fast_c / tot;
        double slow_pct = 100.0 * g_instr.word_slow_c / tot;
        double setup_pct = 100.0 * g_instr.setup_c / tot;
        double crt_pct = 100.0 * g_instr.crt_node_c / tot;
        double testK_pct = 100.0 * g_instr.testK_c / tot;
        fprintf(f, "  %% of total: fast=%.1f  slow=%.1f  setup=%.1f  crt=%.1f  testK=%.1f\n",
                fast_pct, slow_pct, setup_pct, crt_pct, testK_pct);
    }
}

#define INSTR_START(var) u64 _t_##var = rdtsc_now()
#define INSTR_END(name)  do { u64 _e = rdtsc_now(); g_instr.name##_c += _e - _t_##name; g_instr.name##_n++; } while (0)
/* Use INSTR_END_VAR when the accumulator name differs from the START var name.
 * e.g., INSTR_START(word_fast) then on slow path: INSTR_END_VAR(word_slow, word_fast) */
#define INSTR_END_VAR(name, var)  do { u64 _e = rdtsc_now(); g_instr.name##_c += _e - _t_##var; g_instr.name##_n++; } while (0)
#else
#define INSTR_START(var) ((void)0)
#define INSTR_END(name)  ((void)0)
#define INSTR_END_VAR(name, var) ((void)0)
static inline void instr_reset(void) {}
static inline void instr_print(FILE *f) { (void)f; }
#endif

/* ── Jacobi symbol (odd n > 0) ────────────────────────────────────────── */
static int jacobi(long a, long n) {
    int r = 1;
    a = ((a % n) + n) % n;
    while (a) {
        while (!(a & 1)) {
            a >>= 1;
            if ((n & 7) == 3 || (n & 7) == 5) r = -r;
        }
        { long t = n; n = a; a = t; }
        if ((a & 3) == 3 && (n & 3) == 3) r = -r;
        a %= n;
    }
    return (n == 1) ? r : 0;
}

/* ── Modular inverse (gcd(a,m) must be 1) ────────────────────────────── */
static long modinv(long a, long m) {
    long g = m, x = 0, y = 1, ta = a;
    while (ta) {
        long q = g / ta;
        g -= q * ta; { long t = g; g = ta; ta = t; }
        x -= q * y;  { long t = x; x = y;  y = t; }
    }
    return (x % m + m) % m;
}

/* ── Per-prime data ───────────────────────────────────────────────────── */
typedef struct {
    int    p;
    int    nvp;
    int    vp[MAX_NVP];
    double density;
} PrimeSieve;

typedef struct {
    int              np;
    PrimeSieve       ps[MAX_NP];
    int              k_parity;    /* (N%4==1) ? 0 : 1              */
    int              nval_mod8;   /* ((N-1)/6) % 8                 */
    int              N_mod8;      /* N % 8                         */
    int              T12_mask;    /* bitmask of valid T%12 values  */
    const uint8_t   *QR144;       /* points to QR144_A or QR144_B  */
    long             N_mod_p2[MAX_NP];  /* N mod p² for each CRT prime  */
    enum sieve_mode  mode;        /* MODE_K_SIEVE / MODE_S_SIEVE / MODE_M_SIEVE */
    long             n_mod_p[MAX_NP+1]; /* n=(N-1)/6 mod p, for S/M-sieve poly */
} Sieve;

static int cmp_density(const void *a, const void *b) {
    double da = ((const PrimeSieve *)a)->density;
    double db = ((const PrimeSieve *)b)->density;
    return (da < db) ? -1 : (da > db) ? 1 : 0;
}

/* Forward declarations for QR144 tables (defined later with other QR tables) */
static const uint8_t QR144_N1[144];
static const uint8_t QR144_N3[144];
static const uint8_t QR144_N5[144];
static const uint8_t QR144_N7[144];
/* QR144 tables for S-sieve (see definitions for derivation). */
static const uint8_t QR144S_N1[144];
static const uint8_t QR144S_N3[144];
static const uint8_t QR144S_N5[144];
static const uint8_t QR144S_N7[144];
/* QR144 tables for M-sieve (new form, g_M = 9M²-N).
 * Partitioned by N mod 8: determines L parity and L² mod 16.
 * g mod 9 ∈ {1,4,7} always (from L ≡ 1 mod 3 → L² mod 9 ∈ {1,4,7}).
 * Combined with N mod 8 → 3 or 6 valid residues out of 144.
 * Verified zero false negatives on 10000 new-form semiprimes. */
static const uint8_t QR144M_N1[144];
static const uint8_t QR144M_N3[144];
static const uint8_t QR144M_N5[144];
static const uint8_t QR144M_N7[144];

static void sieve_build(Sieve *sv, const mpz_t N, int np, enum sieve_mode mode, u64 S_min_l) {
    sv->np        = np;
    sv->mode      = mode;
    sv->k_parity  = (mpz_fdiv_ui(N, 4) == 1) ? 0 : 1;

    mpz_t n_val; mpz_init(n_val);
    /* For M-sieve: n = (N+1)/6 since N = 6n-1.
     * For K/S-sieve: n = (N-1)/6 since N = 6n+1.
     * Keep the same variable name n_val for parallel code structure. */
    if (mode == MODE_M_SIEVE) {
        mpz_add_ui(n_val, N, 1);
        mpz_fdiv_q_ui(n_val, n_val, 6);
    } else {
        mpz_sub_ui(n_val, N, 1);
        mpz_fdiv_q_ui(n_val, n_val, 6);
    }
    sv->nval_mod8 = (int)mpz_fdiv_ui(n_val, 8);

    sv->N_mod8 = (int)mpz_fdiv_ui(N, 8);

    /* T mod 12 valid bitmask from N mod 8. For K-sieve this was derived for 6k-1 pairs.
     * For M-sieve the T constraints are different (T = 3M not 3S-1), so we set a
     * permissive mask and rely on other filters. T12_mask is primarily used in the
     * K-sieve verify chain; the M-sieve uses its own L ≡ 1 mod 3 check. */
    switch (sv->N_mod8) {
        case 3:  sv->T12_mask = (1<<2);          break;
        case 7:  sv->T12_mask = (1<<8);          break;
        default: sv->T12_mask = (1<<5)|(1<<11);  break;
    }
    switch (mode) {
        case MODE_K_SIEVE:
            switch (sv->N_mod8) {
                case 1:  sv->QR144 = QR144_N1; break;
                case 3:  sv->QR144 = QR144_N3; break;
                case 5:  sv->QR144 = QR144_N5; break;
                default: sv->QR144 = QR144_N7; break;
            }
            break;
        case MODE_S_SIEVE:
            switch (sv->N_mod8) {
                case 1:  sv->QR144 = QR144S_N1; break;
                case 3:  sv->QR144 = QR144S_N3; break;
                case 5:  sv->QR144 = QR144S_N5; break;
                default: sv->QR144 = QR144S_N7; break;
            }
            break;
        case MODE_M_SIEVE:
            switch (sv->N_mod8) {
                case 1:  sv->QR144 = QR144M_N1; break;
                case 3:  sv->QR144 = QR144M_N3; break;
                case 5:  sv->QR144 = QR144M_N5; break;
                default: sv->QR144 = QR144M_N7; break;
            }
            break;
    }

    /* ── Root modulus ──────────────────────────────────────────────────────
     * K-sieve: mod 4 (N%8 ∈ {1,5}) or mod 2 (N%8 ∈ {3,7}).
     * S-sieve: root modulus is MULTIPLIED BY 3 to enforce S ≡ -n (mod 3),
     *   which follows from 6xy = n + S (so mod 3: 0 ≡ n + S).
     *   Gives p=12 (N%8 ∈ {3,7}) or p=6 (N%8 ∈ {1,5}) with nvp=1.
     *   This triples M (3× fewer K0s, 3× smaller n_words) for free.
     * M-sieve: mod 2 or mod 4 only. NO mod-3 fold — 6xy = n + K_Δ folds to K_Δ
     *   (which M-sieve does not scan), and M is free mod 3 (verified empirically).
     *   Root ps: N≡1,5 mod 8 → mod 2 (M odd); N≡3 → mod 4, target 2; N≡7 → mod 4, target 0. */
    if (mode == MODE_K_SIEVE) {
        /* K-sieve: original mod-4/mod-2 logic */
        if (sv->N_mod8 == 1 || sv->N_mod8 == 5) {
            sv->ps[0].p      = 4;
            sv->ps[0].nvp    = 1;
            sv->ps[0].vp[0]  = (sv->N_mod8 == 1) ? 0 : 2;
            sv->ps[0].density = 0.25;
        } else {
            sv->ps[0].p      = 2;
            sv->ps[0].nvp    = 1;
            sv->ps[0].vp[0]  = sv->k_parity;
            sv->ps[0].density = 0.5;
        }
    } else if (mode == MODE_S_SIEVE) {
        /* S-sieve with mod-3 folded: target for (S - S_min) mod (root_p).
         * Compute S_target mod root_p from:
         *   S ≡ S_mod_base (mod base)  where base ∈ {2, 4}
         *   S ≡ -n        (mod 3)
         * then ps[0].vp[0] = (S_target - S_min) mod root_p. */
        long n3 = (long)mpz_fdiv_ui(n_val, 3);         /* n mod 3 */
        int s3 = (int)((3 - n3) % 3);                  /* S ≡ -n mod 3 */
        int base, s_base_target;
        if (sv->N_mod8 == 3) { base = 4; s_base_target = 1; }
        else if (sv->N_mod8 == 7) { base = 4; s_base_target = 3; }
        else                 { base = 2; s_base_target = 0; }  /* N%8 ∈ {1,5}: S even */

        int root_p = base * 3;                          /* 12 or 6 */

        /* CRT: find S_target mod root_p such that
         *   S_target ≡ s_base_target (mod base)
         *   S_target ≡ s3           (mod 3)         */
        int s_target = -1;
        for (int t = 0; t < root_p; t++) {
            if ((t % base) == s_base_target && (t % 3) == s3) { s_target = t; break; }
        }
        /* s_target always exists since gcd(base, 3) = 1. */

        int s_min_mod = (int)(S_min_l % (u64)root_p);
        int vp0 = ((s_target - s_min_mod) % root_p + root_p) % root_p;

        sv->ps[0].p       = root_p;
        sv->ps[0].nvp     = 1;
        sv->ps[0].vp[0]   = vp0;
        sv->ps[0].density = 1.0 / (double)root_p;
    } else {
        /* M-sieve root modulus: mod 2 or mod 4, no mod-3 fold.
         * Compute target for (M - M_min) mod root_p. Here S_min_l = M_min. */
        int root_p, m_target;
        if (sv->N_mod8 == 3)      { root_p = 4; m_target = 2; }
        else if (sv->N_mod8 == 7) { root_p = 4; m_target = 0; }
        else                      { root_p = 2; m_target = 1; }  /* N%8 ∈ {1,5}: M odd */

        int m_min_mod = (int)(S_min_l % (u64)root_p);
        int vp0 = ((m_target - m_min_mod) % root_p + root_p) % root_p;

        sv->ps[0].p       = root_p;
        sv->ps[0].nvp     = 1;
        sv->ps[0].vp[0]   = vp0;
        sv->ps[0].density = 1.0 / (double)root_p;
    }

    /* ── Per-prime Jacobi + Hensel enumeration ─────────────────────────── */
    /* Prime selection for CRT — three modes, in priority order:
     *   (1) KSIEVE_CRT_PRIMES env: explicit prime list (experimental override).
     *   (2) Default: adaptive selection — pick np primes with smallest nvp/p
     *       from candidates < 41 (avoid BS_POOL overlap). Reduces K0-count
     *       variance across N's at same bit size.
     *   (3) KSIEVE_ADAPTIVE_PRIMES=0: legacy fixed SP[0..np-1]. */
    int local_primes[MAX_NP];
    int use_custom_primes = 0;
    {
        const char *env = getenv("KSIEVE_CRT_PRIMES");
        if (env && *env) {
            int cnt = 0;
            const char *s = env;
            while (*s && cnt < np) {
                while (*s == ',' || *s == ' ') s++;
                if (!*s) break;
                int v = 0;
                while (*s >= '0' && *s <= '9') { v = v*10 + (*s - '0'); s++; }
                if (v > 3) local_primes[cnt++] = v;
                while (*s && *s != ',' && *s != ' ') s++;
            }
            if (cnt == np) use_custom_primes = 1;
        }
    }

    if (!use_custom_primes) {
        /* Adaptive selection is OFF by default. Set KSIEVE_ADAPTIVE_PRIMES=1 to opt in.
         * Experimental results were mixed: sometimes 60% faster, sometimes 43% slower.
         * The DFS leaf count grows with prime size, so swapping small primes hurts. */
        int adaptive_enabled = 0;
        const char *env2 = getenv("KSIEVE_ADAPTIVE_PRIMES");
        if (env2 && *env2 == '1') adaptive_enabled = 1;

        if (adaptive_enabled) {
            /* Adaptive selection: fix the (np-1) smallest CRT primes (SP[0..np-2]) and
             * only choose the LAST one from a small pool of candidates at the tail.
             * Why: replacing small primes explodes DFS leaf count (nvp grows with p,
             * so product(nvp[i]) explodes if we swap small for medium primes).
             * Swapping only the LAST prime keeps DFS cost stable while giving
             * 10-15% filter improvement on lucky N's — the regime verified empirically. */
            const int tail_slack = 4;  /* candidates for the last prime: SP[np-1..np-1+slack-1] */
            int pool[16];
            int pool_size = 0;
            for (int i = np - 1; i < NSP && i < np - 1 + tail_slack; i++) {
                if (SP[i] >= 41) break;
                pool[pool_size++] = SP[i];
            }

            if (pool_size <= 1) {
                /* No choice available — just use the default. */
                adaptive_enabled = 0;
            } else {
                /* Fix the first np-1 primes as SP[0..np-2]. */
                for (int i = 0; i < np - 1; i++) local_primes[i] = SP[i];

                /* Pick the last prime: smallest nvp/p among candidates. */
                int best_p = pool[0];
                double best_density = 2.0;
                for (int i = 0; i < pool_size; i++) {
                    int p = pool[i];
                    long Np = (long)mpz_fdiv_ui(N, (unsigned long)p);
                    long np_val_p = (long)mpz_fdiv_ui(n_val, (unsigned long)p);
                    int cnt = 0;
                    for (int k = 0; k < p; k++) {
                        long val;
                        if (mode == MODE_K_SIEVE) {
                            val = (9LL * k * k + Np) % p;
                        } else if (mode == MODE_S_SIEVE) {
                            long Sm = S_min_l % p;
                            long sk = (Sm + k) % p;
                            val = ((9LL * sk % p * sk % p + p - 6LL * sk % p + p
                                    - 6LL * np_val_p % p + p) % p + p) % p;
                        } else {
                            /* M-sieve: g_M(M) = 9M² - N where M = M_min + k */
                            long Mm = S_min_l % p;
                            long mk = (Mm + k) % p;
                            val = ((9LL * mk % p * mk % p + p - Np % p) % p + p) % p;
                            (void)np_val_p;
                        }
                        int j = jacobi(val, p);
                        if (j >= 0) cnt++;
                    }
                    double d = (double)cnt / (double)p;
                    if (d < best_density) {
                        best_density = d;
                        best_p = p;
                    }
                }
                local_primes[np - 1] = best_p;

                /* Re-sort ascending so prime order matches existing code assumptions. */
                for (int i = 0; i < np - 1; i++) {
                    for (int j = i + 1; j < np; j++) {
                        if (local_primes[j] < local_primes[i]) {
                            int t = local_primes[i]; local_primes[i] = local_primes[j]; local_primes[j] = t;
                        }
                    }
                }
                use_custom_primes = 1;
            }
        }
    }

    for (int i = 0; i < np; i++) {
        int  p   = use_custom_primes ? local_primes[i] : SP[i];
        long Np  = (long)mpz_fdiv_ui(N, (unsigned long)p);
        long np_val = (long)mpz_fdiv_ui(n_val, (unsigned long)p);
        long Np2 = (long)mpz_fdiv_ui(N, (unsigned long)(p * p));
        long Np4 = (long)mpz_fdiv_ui(N, (unsigned long)(p * p) * (unsigned long)(p * p));
        sv->n_mod_p[i + 1] = np_val;
        sv->ps[i + 1].p   = p;
        sv->ps[i + 1].nvp = 0;

        /* KSIEVE_W_SORT env flag: when set, puts Jac=+1 residues before Jac=0 residues
         * in vp[]. This implements "W-ranking" — the CRT DFS visits high-W subtrees
         * (where most primes contribute +1/p) before low-W subtrees (where some primes
         * have p|T, contributing 0). Theoretical expected speedup ~1.35× from simulation. */
        int w_sort_enabled = 0;
        {
            const char *env = getenv("KSIEVE_W_SORT");
            if (env && *env == '1') w_sort_enabled = 1;
        }

        /* Temporary buffers for two-pass collection when w_sort_enabled */
        int vp_jac_plus[MAX_NVP];   int n_plus = 0;
        int vp_jac_zero[MAX_NVP];   int n_zero = 0;

        for (int k = 0; k < p; k++) {
            long val;
            if (mode == MODE_K_SIEVE) {
                val = (9LL * k * k + Np) % p;               /* g_K = 9K² + N */
            } else if (mode == MODE_S_SIEVE) {
                /* g'(k) = g_S(S_min + k) = 9(S_min+k)² - 6(S_min+k) - 6n
                 * Precompute S_min mod p once per prime. */
                long Sm = S_min_l % p;
                long sk = (Sm + k) % p;
                val = ((9LL * sk % p * sk % p + p - 6LL * sk % p + p
                        - 6LL * np_val % p + p) % p + p) % p;
            } else {
                /* M-sieve: g_M(M_min + k) = 9(M_min+k)² - N */
                long Mm = S_min_l % p;
                long mk = (Mm + k) % p;
                val = ((9LL * mk % p * mk % p + p - Np % p) % p + p) % p;
            }
            int  j   = jacobi(val, p);
            if (j > 0) {
                if (w_sort_enabled) vp_jac_plus[n_plus++] = k;
                else sv->ps[i + 1].vp[sv->ps[i + 1].nvp++] = k;
            } else if (j == 0) {
                /* Hensel lift — same structure, different polynomial evaluation */
                long p2 = (long)p * p;
                long p3 = p2 * (long)p;
                long p4 = p2 * p2;
                int keep = 0;
                long Np3 = (long)mpz_fdiv_ui(N, (unsigned long)p3);
                long np3 = (long)mpz_fdiv_ui(n_val, (unsigned long)p3);
                long Np4v = (long)mpz_fdiv_ui(N, (unsigned long)p4);
                long np4 = (long)mpz_fdiv_ui(n_val, (unsigned long)p4);
                (void)np3; (void)np4;  /* unused in M-sieve branch */
                for (int jj = 0; jj < p && !keep; jj++) {
                    long kp = k + (long)jj * p;
                    long v2;
                    if (mode == MODE_K_SIEVE) {
                        v2 = (9LL * kp * kp + Np2) % p2;
                    } else if (mode == MODE_S_SIEVE) {
                        long np2 = (long)mpz_fdiv_ui(n_val, (unsigned long)p2);
                        long sk2 = (S_min_l % p2 + kp) % p2;
                        v2 = ((9LL * (sk2 * sk2 % p2) % p2 + p2 - 6LL * sk2 % p2 + p2
                               - 6LL * np2 % p2 + p2) % p2 + p2) % p2;
                    } else {
                        /* M-sieve mod p² */
                        long Np2m = Np2;
                        long mk2 = (S_min_l % p2 + kp) % p2;
                        v2 = ((9LL * (mk2 * mk2 % p2) % p2 + p2 - Np2m % p2) % p2 + p2) % p2;
                    }
                    if (v2 != 0) continue;

                    for (int cc = 0; cc < p && !keep; cc++) {
                        long kp3 = kp + (long)cc * p2;
                        long v3;
                        if (mode == MODE_K_SIEVE) {
                            v3 = (9LL * kp3 * kp3 + Np3) % p3;
                        } else if (mode == MODE_S_SIEVE) {
                            long sk3 = (S_min_l % p3 + kp3) % p3;
                            v3 = ((9LL * (sk3 * sk3 % p3) % p3 + p3 - 6LL * sk3 % p3 + p3
                                   - 6LL * np3 % p3 + p3) % p3 + p3) % p3;
                        } else {
                            long mk3 = (S_min_l % p3 + kp3) % p3;
                            v3 = ((9LL * (mk3 * mk3 % p3) % p3 + p3 - Np3 % p3) % p3 + p3) % p3;
                        }
                        if (v3 != 0) { keep = 1; break; }

                        for (int dd = 0; dd < p && !keep; dd++) {
                            long kp4 = kp3 + (long)dd * p3;
                            long v4;
                            if (mode == MODE_K_SIEVE) {
                                v4 = (9LL * kp4 * kp4 + Np4v) % p4;
                            } else if (mode == MODE_S_SIEVE) {
                                long sk4 = (S_min_l % p4 + kp4) % p4;
                                v4 = ((9LL * (sk4 * sk4 % p4) % p4 + p4 - 6LL * sk4 % p4 + p4
                                       - 6LL * np4 % p4 + p4) % p4 + p4) % p4;
                            } else {
                                long mk4 = (S_min_l % p4 + kp4) % p4;
                                v4 = ((9LL * (mk4 * mk4 % p4) % p4 + p4 - Np4v % p4) % p4 + p4) % p4;
                            }
                            if (v4 != 0) continue;
                            keep = 1;
                        }
                    }
                }
                if (keep) {
                    if (w_sort_enabled) vp_jac_zero[n_zero++] = k;
                    else sv->ps[i + 1].vp[sv->ps[i + 1].nvp++] = k;
                }
            }
        }

        /* If W-sort was requested, commit the two-pass ordering: Jac=+1 first, Jac=0 last. */
        if (w_sort_enabled) {
            for (int a = 0; a < n_plus; a++)
                sv->ps[i + 1].vp[sv->ps[i + 1].nvp++] = vp_jac_plus[a];
            for (int a = 0; a < n_zero; a++)
                sv->ps[i + 1].vp[sv->ps[i + 1].nvp++] = vp_jac_zero[a];
        }

        sv->ps[i + 1].density = (double)sv->ps[i + 1].nvp / p;
    }

    qsort(sv->ps + 1, (size_t)np, sizeof(PrimeSieve), cmp_density);
    sv->np = np + 1;

    sv->N_mod_p2[0] = 0;
    for (int i = 1; i < sv->np; i++) {
        long p = sv->ps[i].p;
        sv->N_mod_p2[i] = (long)mpz_fdiv_ui(N, (unsigned long)(p * p));
        sv->n_mod_p[i] = (long)mpz_fdiv_ui(n_val, (unsigned long)p);
    }
    mpz_clear(n_val);
}

/* ── Result ───────────────────────────────────────────────────────────── */
typedef struct {
    int    success;
    int    np, N_bits, K_bits;
    long   n_crt_ops, n_k0s, n_isqrt;
    double t_total, t_precomp, t_p1, t_p2;
} Res;

/* ─────────────────────────────────────────────────────────────────────────
 * np SELECTION
 * Minimise estimated wall-clock = phase1_cost + phase2_cost.
 * Works for any bit size; M stored as double for estimation.
 * ───────────────────────────────────────────────────────────────────────── */
static int choose_np(int N_bits, double K_max_d, enum sieve_mode mode) {
    /* Simple and empirically validated: add CRT primes until n_words
     * (= J_max / 64) falls into the sweet spot [~100, TARGET_NW_MAX].
     * This balances bitset scan cost (proportional to n_words) against
     * CRT tree cost (grows with np). Validated against empirical np sweeps
     * at 80-130 bit on 16-core AVX-512 hardware.
     *
     * Root modulus:
     *   K-sieve: 4 (mod-4) for half of Ns, 2 (mod-2) for the other half.
     *   S-sieve: 12 (N%8 ∈ {3,7}) or 6 (N%8 ∈ {1,5})  — mod-3 folded into root.
     *   M-sieve: 4 (mod-4) for half of Ns, 2 (mod-2) for the other half.
     *            No mod-3 fold (M is uniform mod 3, verified empirically).
     * Use the conservative (larger) estimate so we don't over-add CRT primes. */
    const int TARGET_NW_MAX = 600;
    double M;
    switch (mode) {
        case MODE_K_SIEVE: M = 4.0;  break;  /* mod 4 conservative */
        case MODE_S_SIEVE: M = 12.0; break;  /* mod 12 with mod-3 fold */
        case MODE_M_SIEVE: M = 4.0;  break;  /* mod 4 conservative, no mod-3 */
        default:           M = 4.0;  break;
    }
    int best_np = 4;

    for (int i = 0; i < NSP && i < MAX_NP; i++) {
        M *= SP[i];
        double J_max_f = (M >= K_max_d) ? 1.0 : K_max_d / M;
        long nw = (long)((J_max_f + 63) / 64);
        if (nw < 1) nw = 1;
        best_np = i + 1;
        if (nw <= TARGET_NW_MAX) break;
    }

    if (best_np < 4) best_np = 4;
    if (best_np > NSP) best_np = NSP;

    /* K-sieve has 8.5× larger range than S-sieve. The same np gives
     * 8.5× more n_words. Subtract 1 to compensate (empirically validated
     * at 90-100 bit K-sieve).
     * M-sieve has the same tight range as S-sieve, but 3× less root fold,
     * so it sees 3× more n_words than S-sieve at matched np. Don't subtract. */
    if (mode == MODE_K_SIEVE && N_bits >= 85 && best_np > 4) best_np--;

    return best_np;
}

/* ─────────────────────────────────────────────────────────────────────────
 * PHASE 1: CRT DFS
 * K_acc and mod_acc use __int128 so M can exceed 2^64 (needed for np >= 15).
 * Collected K0 values always fit in long: K0 <= K_max <= sqrt(N)/3 <= 2^49.
 * K0s array grows dynamically via realloc.
 * ───────────────────────────────────────────────────────────────────────── */
static long *phase1_dfs(const Sieve *sv, long K_max_l,
                        long *nK0_out, long *ncrt_out) {
    int np = sv->np;

    /* inv_cache[d] = modinv(prod(p[0..d-1]), p[d]) */
    long inv_cache[MAX_NP];
    {
        i128 running = 1;
        for (int d = 0; d < np; d++) {
            int p = sv->ps[d].p;
            inv_cache[d] = (d == 0) ? 1 : modinv((long)(running % (i128)p), p);
            running *= (i128)p;
        }
    }

    long  cap = 1L << 20;
    long *K0s = (long *)malloc((size_t)cap * sizeof(long));
    if (!K0s) { perror("malloc K0s"); exit(1); }
    long nK0 = 0, ncrt = 0;

    int  vi     [MAX_NP];
    i128 K_acc  [MAX_NP + 1];
    i128 mod_acc[MAX_NP + 1];

    K_acc[0]   = 0;
    mod_acc[0] = 1;
    for (int i = 0; i < np; i++) vi[i] = 0;

    int depth = 0;
    while (depth >= 0) {
        if (vi[depth] >= sv->ps[depth].nvp) {
            vi[depth] = 0;
            depth--;
            if (depth >= 0) vi[depth]++;
            continue;
        }

        int  p   = sv->ps[depth].p;
        int  r   = sv->ps[depth].vp[vi[depth]];
        i128 Ka  = K_acc[depth];
        i128 ma  = mod_acc[depth];

        long Ka_p  = (long)(Ka % (i128)p);
        long diff  = ((long)(r - Ka_p) % p + p) % p;
        long t     = diff * inv_cache[depth] % p;
        i128 K_new = Ka + ma * (i128)t;
        ncrt++;

        if (K_new > (i128)K_max_l) {
            vi[depth]++;
            continue;
        }

        K_acc[depth + 1]   = K_new;
        mod_acc[depth + 1] = ma * (i128)p;

        if (depth + 1 == np) {
            if (K_new > 0) {
                if (nK0 >= cap) {
                    cap *= 2;
                    K0s  = (long *)realloc(K0s, (size_t)cap * sizeof(long));
                    if (!K0s) { perror("realloc K0s"); exit(1); }
                }
                K0s[nK0++] = (long)K_new;
            }
            vi[depth]++;
            if (vi[depth] >= sv->ps[depth].nvp) {
                vi[depth] = 0;
                depth--;
                if (depth >= 0) vi[depth]++;
            }
        } else {
            depth++;
        }
    }

    *nK0_out  = nK0;
    *ncrt_out = ncrt;
    return K0s;
}

/* ─────────────────────────────────────────────────────────────────────────
 * QR TABLES for mod pre-screens
 *
 * g = 9K²+N must be a perfect square. These tables pre-reject candidates
 * where g is provably a non-square mod m — machine-word only, no GMP.
 * Combined pass rate after all screens: ~0.05% reach GMP isqrt.
 * ───────────────────────────────────────────────────────────────────────── */

/* QR mod 256: 44/256 pass (~17.2%).  Subsumes mod 64 (strictly stronger). */
static const uint8_t QR256[256] = {
    1,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0, 1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,
    0,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0, 0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,
    1,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0, 0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,
    0,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0, 0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,
    0,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0, 1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,
    0,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0, 0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,
    0,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0, 0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,
    0,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0, 0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,
};

/* QR mod 2^16: 10924/65536 pass (~16.7%).  Generated at startup. */
static uint8_t QR65536[65536];

/* QR mod 45 = 9×5: 12/45 pass (~26.7%) */
/* (QR45 removed: mod-5 part redundant since 5|M; mod-9 part dominated by QR27) */

/* QR mod 49 = 7²: 22/49 pass (~44.9%) — p²-valuation of 7, varies since M%49≠0 */
static const uint8_t QR49[49] = {
    1,1,1,0,1,0,0,0,1,1, 0,1,0,0,0,1,1,0,1,0,
    0,0,1,1,0,1,0,0,0,1, 1,0,1,0,0,0,1,1,0,1,
    0,0,0,1,1,0,1,0,0
};

/* QR mod 31: 16/31 pass (51.6%) — extra prime, M%31≠0 so varies with stride */
static const uint8_t QR31[31] = {
    1,1,1,0,1,1,0,1,1,1,
    1,0,0,0,1,0,1,0,1,1,
    1,0,0,0,0,1,0,0,1,0,
    0
};

/* QR mod 37: 19/37 pass (51.4%) */
static const uint8_t QR37[37] = {
    1,1,0,1,1,0,0,1,0,1,
    1,1,1,0,0,0,1,0,0,0,
    0,1,0,0,0,1,1,1,1,0,
    1,0,0,1,1,0,1
};

/* QR mod 41: 21/41 pass (51.2%) */
static const uint8_t QR41[41] = {
    1,1,1,0,1,1,0,0,1,1,
    1,0,0,0,0,0,1,0,1,0,
    1,1,0,1,0,1,0,0,0,0,
    0,1,1,1,0,0,1,1,0,1,
    1
};

/* QR mod 27 = 3³: 11/27 pass (~40.7%) — best 3-adic screen */
static const uint8_t QR27[27] = {
    1,1,0,0,1,0,0,1,0,1,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0
};

/* QR mod 25 = 5²: 11/25 pass (~44%) */
static const uint8_t QR25[25] = {
    1,1,0,0,1,0,1,0,0,1,0,1,0,0,1,0,1,0,0,1,0,1,0,0,1
};

/* QR mod 121 = 11^2: 56/121 pass (~46.3%) — p² filter, independent of QR64 */
static const uint8_t QR121[121] = {
    1,1,0,1,1,1,0,0,0,1,0,0,1,0,1,1,1,0,0,0,1,0,0,1,0,1,1,1,0,0,0,1,
    0,0,1,0,1,1,1,0,0,0,1,0,0,1,0,1,1,1,0,0,0,1,0,0,1,0,1,1,1,0,0,0,
    1,0,0,1,0,1,1,1,0,0,0,1,0,0,1,0,1,1,1,0,0,0,1,0,0,1,0,1,1,1,0,0,
    0,1,0,0,1,0,1,1,1,0,0,0,1,0,0,1,0,1,1,1,0,0,0,1,0
};

/* QR mod 169 = 13^2: 79/169 pass (~46.7%) — p² filter */
static const uint8_t QR169[169] = {
    1,1,0,1,1,0,0,0,0,1,1,0,1,0,1,0,1,1,0,0,0,0,1,1,0,1,0,1,0,1,1,0,
    0,0,0,1,1,0,1,0,1,0,1,1,0,0,0,0,1,1,0,1,0,1,0,1,1,0,0,0,0,1,1,0,
    1,0,1,0,1,1,0,0,0,0,1,1,0,1,0,1,0,1,1,0,0,0,0,1,1,0,1,0,1,0,1,1,
    0,0,0,0,1,1,0,1,0,1,0,1,1,0,0,0,0,1,1,0,1,0,1,0,1,1,0,0,0,0,1,1,
    0,1,0,1,0,1,1,0,0,0,0,1,1,0,1,0,1,0,1,1,0,0,0,0,1,1,0,1,0,1,0,1,
    1,0,0,0,0,1,1,0,1
};

/* QR144: encodes T%12 constraint as a pre-isqrt screen on g mod 144.
 * Since T²=g, g mod 144 determines which T mod 12 values are possible.
 * Four tables — one per N%8 class — selected at sieve_build time.
 * This moves the T%12 check BEFORE isqrt, saving the isqrt cost entirely
 * for candidates that would have failed the T%12 post-check.
 * Pass rates: 2.1% for N%8∈{3,7}, 4.2% for N%8∈{1,5}.                 */

/* N%8=1: T%12 in {5,11}, 6/144 pass (4.2%) */
static const uint8_t QR144_N1[144] = {
    0,1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,
    0,1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,
    0,1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,
    0,1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,
    0,1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,
    0,1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,
};

/* N%8=3: T%12 in {2}, 3/144 pass (2.1%) */
static const uint8_t QR144_N3[144] = {
    0,0,0,0,1,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,1,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,1,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,
};

/* N%8=5: T%12 in {5,11}, 6/144 pass (4.2%) */
static const uint8_t QR144_N5[144] = {
    0,1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,
    0,1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,
    0,1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,
    0,1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,
    0,1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,
    0,1,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,
};

/* N%8=7: T%12 in {8}, 3/144 pass (2.1%) */
static const uint8_t QR144_N7[144] = {
    0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,1,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,1,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,1,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,
};

/* ─── S-sieve QR144 tables ───────────────────────────────────────────────
 * At the true S: g_S = 9S²-6S-6n = (3K)². So valid g_S values mod 144 are
 * restricted to squares (3K mod 144)² whose sqrt 3K mod 12 satisfies the
 * T-parity constraint derived from N mod 8 (via T = 3S-1):
 *   N%8 = 1  → 3K ≡ 0  mod 12  →  good g_S mod 144 = {0}
 *   N%8 = 3  → 3K ∈ {3,9} mod 12 → good g_S mod 144 = {9, 81}
 *   N%8 = 5  → 3K ≡ 6  mod 12  →  good g_S mod 144 = {36}
 *   N%8 = 7  → 3K ∈ {3,9} mod 12 → good g_S mod 144 = {9, 81}
 * Each S-sieve table is indexed by g_S mod 144.
 *
 * Empirical filter strength: ~12/144 = 3.58 raw bits (measured over all
 * admissible n residues).  Marginal bits in the verify chain are smaller
 * due to correlation with QR256 / QR27, but non-negligible.
 * True-S pass rate: 100% (verified over 400 random 40-bit semiprimes). */
static const uint8_t QR144S_N1[144] = {
    [0]=1
};
static const uint8_t QR144S_N3[144] = {
    [9]=1, [81]=1
};
static const uint8_t QR144S_N5[144] = {
    [36]=1
};
static const uint8_t QR144S_N7[144] = {
    [9]=1, [81]=1
};

/* ── M-sieve QR144 tables ────────────────────────────────────────────────
 * g_M(M) = 9M² - N at true M equals L² where L = 3K_Δ + 1.
 * Constraint: L ≡ 1 (mod 3), so L² ≡ 1, 4, or 7 (mod 9).
 * Combined with L parity (opposite of T parity, determined by N mod 8):
 *   N≡1 mod 8: L ≡ 0 mod 4,  so L² ≡ 0 mod 16 → 3 values mod 144: {16, 64, 112}
 *   N≡3 mod 8: L odd,         so L² ∈ {1,9} mod 16 → 6 values: {1,25,49,73,97,121}
 *   N≡5 mod 8: L ≡ 2 mod 4,  so L² ≡ 4 mod 16 → 3 values: {4, 52, 100}
 *   N≡7 mod 8: L odd,         so L² ∈ {1,9} mod 16 → 6 values: {1,25,49,73,97,121}
 * Raw filter strength: ~4.6-5.6 bits (3 or 6 of 144 pass).
 * Generated and verified zero false negatives over 10000 random new-form semiprimes. */
static const uint8_t QR144M_N1[144] = {
    [16]=1, [64]=1, [112]=1
};
static const uint8_t QR144M_N3[144] = {
    [1]=1, [25]=1, [49]=1, [73]=1, [97]=1, [121]=1
};
static const uint8_t QR144M_N5[144] = {
    [4]=1, [52]=1, [100]=1
};
static const uint8_t QR144M_N7[144] = {
    [1]=1, [25]=1, [49]=1, [73]=1, [97]=1, [121]=1
};

static void init_qr_tables(void) {
    for (int k = 0; k < 65536; k++)
        QR65536[((unsigned)k * (unsigned)k) & 0xFFFF] = 1;
}

/* ─────────────────────────────────────────────────────────────────────────
 * PHASE 2a: Native u128 path  (N <= 64-bit)
 * ───────────────────────────────────────────────────────────────────────── */
static inline u64 isqrt_u128(u128 g) {
    u64 t = (u64)sqrtl((long double)g);
    if ((u128)(t + 1) * (t + 1) <= g) t++;
    else if ((u128)t * t > g)          t--;
    return t;
}

/* Thread-local scratch mpz_t's for the cross-prime filter's mpz path.
 * Used at N ≥ 127-bit where g = 9K² + N overflows u128. Lazy-init on first use
 * per thread to avoid per-call mpz_inits/clears overhead. */
static __thread mpz_t cp_scratch_g, cp_scratch_T, cp_scratch_K;
static __thread int   cp_scratch_initialized = 0;

static inline void cp_scratch_ensure(void) {
    if (!cp_scratch_initialized) {
        mpz_inits(cp_scratch_g, cp_scratch_T, cp_scratch_K, NULL);
        cp_scratch_initialized = 1;
    }
}

static int phase2_native(const long *K0s, long nK0,
                         u64 N64, u64 M64, u64 K_max,
                         int kpar, int nval8, const uint8_t *QR144_tab,
                         long *n_isqrt_out, mpz_t uo, mpz_t vo) {
    unsigned long N_mod144 = (unsigned long)(N64 % 144);
    unsigned long M_mod144 = (unsigned long)(M64 % 144);
    unsigned long twoA_144 = (18UL * M_mod144 % 144 * M_mod144) % 144;

    for (long ki = 0; ki < nK0; ki++) {
        u64 Kv = (u64)K0s[ki];
        unsigned long K0_144 = (unsigned long)(K0s[ki] % 144);
        unsigned long g_144  = (9UL*(K0_144*K0_144%144)%144 + N_mod144) % 144;
        unsigned long d_144  = (18UL*(K0_144*M_mod144%144)%144 + 9UL*(M_mod144*M_mod144%144)) % 144;

        while (Kv <= K_max) {
            if ((int)(Kv & 1) != kpar) goto next_K_native;
            {
                int K8 = (int)(Kv & 7);
                if (((nval8*nval8 - K8*K8 + 64) & 7) != 0) goto next_K_native;
            }

            /* QR144 pre-screen: encodes T%12 constraint, 2-4% pass */
            if (!QR144_tab[g_144]) goto next_K_native;

            {
                u128 gv = (u128)9 * Kv * Kv + (u128)N64;
                if (!QR256[gv & 255]) goto next_K_native;

                (*n_isqrt_out)++;
                u64 Tv = isqrt_u128(gv);
                if ((u128)Tv * Tv == gv) {
                    mpz_set_ui(uo, (unsigned long)(Tv - 3 * Kv));
                    mpz_set_ui(vo, (unsigned long)(Tv + 3 * Kv));
                    return 1;
                }
            }

        next_K_native:
            g_144 = (g_144 + d_144) % 144;
            d_144 = (d_144 + twoA_144) % 144;
            Kv += M64;
        }
    }
    return 0;
}

/* ─────────────────────────────────────────────────────────────────────────
 * Build the ExtraFilter: precompute valid-K tables for extra primes.
 * Extra primes are chosen from SP[] beyond index np (not already in M).
 * ───────────────────────────────────────────────────────────────────────── */
static void extra_filter_build(ExtraFilter *ef, const mpz_t N,
                                const mpz_t M_mpz, int crt_np, long J_max) {
    memset(ef, 0, sizeof(*ef));
    if (J_max <= 1) return;

    int n = 0;
    long running = 1;
    for (int si = crt_np; si < NSP && n < EXTRA_NP; si++) {
        int q = SP[si];
        if (q >= EXTRA_MAXQ) break;
        /* Stop if adding this prime would overflow long m_acc in DFS */
        if (running > (long)(9e18) / q) break;

        long Nq  = (long)mpz_fdiv_ui(N, (unsigned long)q);
        long Mq  = (long)mpz_fdiv_ui(M_mpz, (unsigned long)q);
        long inv = modinv(Mq, q);

        int nvp = 0;
        for (int k = 0; k < q; k++) {
            long val = (9LL * k * k + Nq) % q;
            if (jacobi(val, q) >= 0)
                ef->valid_v[n][nvp++] = (uint8_t)k;
        }
        ef->q[n]    = q;
        ef->nvp[n]  = nvp;
        ef->Minv[n] = inv;
        running    *= q;
        n++;
    }
    ef->np = n;

    /* Precompute inv_cache[d] = modinv(prod(q[0..d-1]), q[d]) */
    {
        long prod = 1;
        for (int d = 0; d < n; d++) {
            ef->inv_cache[d] = (d == 0) ? 1 : modinv(prod % ef->q[d], ef->q[d]);
            prod *= ef->q[d];
        }
    }
}

/* Mini-DFS on stride index j for one K0.
 * For each extra prime q: valid j values satisfy (K0 + j*M) mod q ∈ valid_K[q],
 * i.e. j ≡ (v - K0) * Minv (mod q) for each valid residue v.
 * This is the same CRT-DFS structure as Phase 1 but over [0, J_max).
 * Collects all valid j into j_out[]. Returns count.
 * With enough primes (prod(q) > J_max): count is 0 or 1 per K0.       */
static int extra_filter_dfs(const ExtraFilter *ef, long K0, long J_max,
                             long *j_out) {
    int np = ef->np;
    if (np == 0) return 0;

    /* Translate valid K residues → valid j residues mod q[d] for this K0.
     * Uses packed valid_v: O(nvp) per prime instead of O(q).              */
    int  nvj[EXTRA_NP];
    int  vj [EXTRA_NP][EXTRA_MAXNVP];

    for (int d = 0; d < np; d++) {
        int  q    = ef->q[d];
        long K0_q = K0 % q;
        long Minv = ef->Minv[d];
        int  nvp  = ef->nvp[d];
        int  n    = 0;
        for (int vi = 0; vi < nvp; vi++) {
            int  v  = ef->valid_v[d][vi];
            long jr = ((long)(v - K0_q) % q + q) % q * Minv % q;
            vj[d][n++] = (int)jr;
        }
        nvj[d] = n;
        if (n == 0) return 0;
    }

    /* Iterative DFS over j in [0, J_max) */
    int  vi   [EXTRA_NP];
    long j_acc[EXTRA_NP + 1];
    long m_acc[EXTRA_NP + 1];
    for (int i = 0; i < np; i++) vi[i] = 0;
    j_acc[0] = 0;
    m_acc[0] = 1;

    int count = 0;
    int depth = 0;
    while (depth >= 0) {
        if (vi[depth] >= nvj[depth]) {
            vi[depth] = 0; depth--;
            if (depth >= 0) vi[depth]++;
            continue;
        }

        int  q     = ef->q[depth];
        int  r     = vj[depth][vi[depth]];
        long ja    = j_acc[depth];
        long ma    = m_acc[depth];
        long diff  = ((long)(r - ja % q) % q + q) % q;
        long t     = diff * ef->inv_cache[depth] % q;
        long j_new = ja + ma * t;

        if (j_new >= J_max) { vi[depth]++; continue; }

        j_acc[depth + 1] = j_new;
        m_acc[depth + 1] = ma * q;

        if (depth + 1 == np) {
            j_out[count++] = j_new;
            vi[depth]++;
            if (vi[depth] >= nvj[depth]) {
                vi[depth] = 0; depth--;
                if (depth >= 0) vi[depth]++;
            }
        } else {
            depth++;
        }
    }
    return count;
}
static int phase2_gmp(const long *K0s, long nK0,
                      const mpz_t N, const mpz_t M_mpz, const mpz_t K_max,
                      int kpar, int nval8, const uint8_t *QR144_tab,
                      const ExtraFilter *ef,
                      long *n_isqrt_out, mpz_t uo, mpz_t vo) {

    /* Only g_65 uses a recurrence. All other screens compute on-demand. */
    unsigned long M_mod65 = (unsigned long)mpz_fdiv_ui(M_mpz, 65536);
    unsigned long twoA_65 = (18UL * M_mod65 % 65536 * M_mod65) % 65536;
    unsigned long N_mod65 = (unsigned long)mpz_fdiv_ui(N, 65536);

    /* N mod small moduli for on-demand g_mod checks */
    unsigned long N_mod144= (unsigned long)mpz_fdiv_ui(N, 144);
    unsigned long N_mod31 = (unsigned long)mpz_fdiv_ui(N, 31);
    unsigned long N_mod37 = (unsigned long)mpz_fdiv_ui(N, 37);
    unsigned long N_mod41 = (unsigned long)mpz_fdiv_ui(N, 41);
    unsigned long N_mod49 = (unsigned long)mpz_fdiv_ui(N, 49);
    unsigned long N_mod27 = (unsigned long)mpz_fdiv_ui(N, 27);
    unsigned long N_mod25 = (unsigned long)mpz_fdiv_ui(N, 25);

    mpz_t K, g, T, T2, threeK;
    mpz_inits(K, g, T, T2, threeK, NULL);

    long M_long  = (long)mpz_get_ui(M_mpz);
    long K_max_l = mpz_get_si(K_max);
    long J_max   = (M_long > 0) ? (K_max_l / M_long + 1) : 1;

    /* Cross-prime: leaf verifier only, no per-step tracking. */
    u64 cp_p = 0, cp_N_mod_p = 0, cp_T_init = 0;
    u128 cp_N_128 = 0;
    {
        mpz_t p_mpz, sqrtN_mpz;
        mpz_inits(p_mpz, sqrtN_mpz, NULL);
        mpz_sqrt(sqrtN_mpz, N);
        cp_T_init = mpz_get_ui(sqrtN_mpz);
        mpz_add_ui(p_mpz, sqrtN_mpz, 1);
        mpz_nextprime(p_mpz, p_mpz);
        if (mpz_fits_ulong_p(p_mpz)) {
            cp_p = mpz_get_ui(p_mpz);
            cp_N_mod_p = (u64)mpz_fdiv_ui(N, (unsigned long)cp_p);
        }
        mpz_clears(p_mpz, sqrtN_mpz, NULL);
        mpz_t nlo, nhi; mpz_inits(nlo, nhi, NULL);
        mpz_tdiv_r_2exp(nlo, N, 64); mpz_tdiv_q_2exp(nhi, N, 64);
        cp_N_128 = ((u128)mpz_get_ui(nhi) << 64) | mpz_get_ui(nlo);
        mpz_clears(nlo, nhi, NULL);
    }

    /* On-demand g mod m check: (9*K² + N) mod m, table lookup. */
    #define G_MOD_CHECK(K_, N_mod_, table_, m_) ({                         \
        unsigned long _k = (unsigned long)(((K_) % (m_) + (m_)) % (m_));   \
        unsigned long _g = (9UL * (_k * _k % (m_)) % (m_) + (N_mod_)) % (m_); \
        (table_)[_g]; })

    /* Bitset path is now handled by the fused CRT+bitset in ksieve_factor.
     * This fallback only handles: DFS pre-filter (J_max > 2000) or lean stride. */
    int use_dfs_bitset = (ef && ef->np > 0 && J_max > 2000);
    int found = 0;

    for (long ki = 0; ki < nK0 && !found; ki++) {
        long K0 = K0s[ki];

        /* Mini-DFS pre-filter path (large J_max only) */
        if (use_dfs_bitset) {
            long j_out[64];
            int  nj = extra_filter_dfs(ef, K0, J_max, j_out);
            if (nj == 0) continue;
            for (int ji = 0; ji < nj && !found; ji++) {
                long K_long = K0 + j_out[ji] * M_long;
                if (K_long > K_max_l) continue;
                if ((K_long & 1) != (unsigned long)kpar) continue;
                { int K8=(int)(K_long&7); if(((nval8*nval8-K8*K8+64)&7)!=0) continue; }
                mpz_set_si(K, K_long);
                if (kpar == 1 && mpz_jacobi(N, K) < 0) continue;
                (*n_isqrt_out)++;
                mpz_mul(g, K, K); mpz_mul_ui(g, g, 9); mpz_add(g, g, N);
                mpz_sqrt(T, g); mpz_mul(T2, T, T);
                if (mpz_cmp(T2, g) != 0) continue;
                mpz_mul_ui(threeK, K, 3);
                mpz_sub(uo, T, threeK); mpz_add(vo, T, threeK);
                found = 1;
            }
            continue;
        }

        /* ── Lean stride loop ─────────────────────────────────────────────
         * Hot path: ONLY g_65 recurrence + QR256/QR65536.
         * Everything else computed on-demand from K_long.                  */
        unsigned long K0_65 = (unsigned long)((K0%65536+65536)%65536);
        unsigned long g_65  = (9UL*(K0_65*K0_65%65536)%65536+N_mod65)%65536;
        unsigned long d_65  = (18UL*(K0_65*M_mod65%65536)%65536+9UL*(M_mod65*M_mod65%65536))%65536;
        long K_long = K0;

#ifdef __AVX2__
        static const int16_t J_arr[16]   = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
        static const int16_t JJ1_arr[16] = {0,0,1,3,6,10,15,21,28,36,45,55,66,78,91,105};
        const __m256i Jv   = _mm256_loadu_si256((const __m256i*)J_arr);
        const __m256i JJ1v = _mm256_loadu_si256((const __m256i*)JJ1_arr);
        const __m128i qr_lo  = _mm_loadu_si128((const __m128i*)(QR256 +   0));
        const __m128i qr_hi  = _mm_loadu_si128((const __m128i*)(QR256 + 128));
        const __m128i mask80 = _mm_set1_epi8((char)0x80);
#endif
        while (K_long <= K_max_l) {
#ifdef __AVX2__
            if (__builtin_expect(K_long + 16 * M_long <= K_max_l, 1)) {
                __m256i g_vec = _mm256_add_epi16(
                    _mm256_set1_epi16((short)g_65),
                    _mm256_add_epi16(
                        _mm256_mullo_epi16(_mm256_set1_epi16((short)d_65), Jv),
                        _mm256_mullo_epi16(_mm256_set1_epi16((short)twoA_65), JJ1v)));
                __m128i lo = _mm256_castsi256_si128(g_vec);
                __m128i hi = _mm256_extracti128_si256(g_vec, 1);
                __m128i bytes = _mm_packus_epi16(
                    _mm_and_si128(lo, _mm_set1_epi16(0xFF)),
                    _mm_and_si128(hi, _mm_set1_epi16(0xFF)));
                __m128i lo_r = _mm_shuffle_epi8(qr_lo, bytes);
                __m128i hi_r = _mm_shuffle_epi8(qr_hi, _mm_xor_si128(bytes, mask80));
                __m128i qr_r = _mm_or_si128(lo_r, hi_r);

                if (_mm_testz_si128(qr_r, qr_r)) {
                    g_65  = (g_65 + 16*d_65 + 120*twoA_65) & 0xFFFF;
                    d_65  = (d_65 + 16*twoA_65) & 0xFFFF;
                    K_long += 16 * M_long;
                    continue;
                }
                int pass_mask = _mm_movemask_epi8(
                    _mm_cmpgt_epi8(qr_r, _mm_setzero_si128()));
                int fj = __builtin_ctz((unsigned)pass_mask);
                if (fj > 0) {
                    unsigned long sf = (unsigned long)(fj - 1);
                    if (sf > 0) {
                        unsigned long _qf = sf*(sf-1)/2;
                        g_65 = (g_65 + sf*d_65 + _qf*twoA_65) & 0xFFFF;
                        d_65 = (d_65 + sf*twoA_65) & 0xFFFF;
                        K_long += (long)sf * M_long;
                    }
                    goto adv_and_check;
                }
            }
#endif
            /* Scalar tier-1 */
            if (!QR256  [g_65 & 0xFF]) goto adv_only;
            if (!QR65536[g_65])        goto adv_only;

            /* On-demand tier-2 (~3% of steps) */
            if (!G_MOD_CHECK(K_long, N_mod144, QR144_tab, 144)) goto adv_only;
            { int K8=(int)(K_long&7); if(((nval8*nval8-K8*K8+64)&7)!=0) goto adv_only; }

            /* On-demand tier-3 (~0.1%) */
            if (!G_MOD_CHECK(K_long, N_mod31, QR31, 31)) goto adv_only;
            if (!G_MOD_CHECK(K_long, N_mod37, QR37, 37)) goto adv_only;
            if (!G_MOD_CHECK(K_long, N_mod41, QR41, 41)) goto adv_only;
            if (!G_MOD_CHECK(K_long, N_mod49, QR49, 49)) goto adv_only;
            if (!G_MOD_CHECK(K_long, N_mod27, QR27, 27)) goto adv_only;
            if (!G_MOD_CHECK(K_long, N_mod25, QR25, 25)) goto adv_only;

            /* On-demand cross-prime (~0.06%) — compute from K_long, no tracking */
            if (cp_p) {
                u64 Ku = (u64)(K_long < 0 ? -K_long : K_long);
                u128 cpg = (u128)9 * Ku * Ku + cp_N_128;
                u64 cpT = cp_T_init + (u64)((u128)9 * Ku * Ku / (2 * (u128)cp_T_init));
                for (int _i = 0; _i < 10; _i++) {
                    u128 sq = (u128)cpT * cpT;
                    u64 nxt;
                    if (sq > cpg) nxt = cpT - (u64)((sq - cpg) / (2*cpT)) - 1;
                    else          nxt = cpT + (u64)((cpg - sq) / (2*cpT));
                    if (nxt == cpT) break;
                    cpT = nxt;
                }
                if ((u128)(cpT+1)*(cpT+1) <= cpg) cpT++;
                if ((u128)cpT*cpT > cpg) cpT--;
                u64 Tm = cpT % cp_p;
                u64 Km = Ku % cp_p;
                u64 gm = (u64)(((u128)9 * Km % cp_p * Km + cp_N_mod_p) % cp_p);
                if ((u64)(((u128)Tm * Tm) % cp_p) != gm) goto adv_only;
            }

            /* Jacobi(N,K) filter */
            mpz_set_si(K, K_long);
            if (kpar == 1 && mpz_jacobi(N, K) < 0) goto adv_only;

            /* GMP isqrt (essentially never reached) */
            (*n_isqrt_out)++;
            mpz_mul(g, K, K); mpz_mul_ui(g, g, 9); mpz_add(g, g, N);
            mpz_sqrt(T, g); mpz_mul(T2, T, T);
            if (mpz_cmp(T2, g) != 0) goto adv_only;

            mpz_mul_ui(threeK, K, 3);
            mpz_sub(uo, T, threeK); mpz_add(vo, T, threeK);
            found = 1; goto done;

        adv_and_check:
        adv_only:
            g_65 = (g_65 + d_65) & 0xFFFF; d_65 = (d_65 + twoA_65) & 0xFFFF;
            K_long += M_long;
        }
    }

done:
    mpz_clears(K, g, T, T2, threeK, NULL);
    return found;
}


/* ─────────────────────────────────────────────────────────────────────────
 * FUSED CRT DFS + BITSET SCAN (no K0 array, u64 arithmetic)
 *
 * Combines phase 1 (CRT DFS) and phase 2 (bitset mini-CRT) into a single
 * pass. At each DFS leaf, immediately runs the bitset check and verification.
 * Eliminates multi-GB K0 array allocation. Uses u64 instead of __int128 for
 * all CRT arithmetic (valid for N ≤ ~126-bit where K_max < 2^62).
 * ───────────────────────────────────────────────────────────────────────── */

/* Adaptive bitset mini-CRT: extra primes and word count chosen at runtime.
 * LAYERED BITSET: each prime stores q entries × 1 u64 (tiled 64-bit pattern).
 * Total table fits L1/L2 regardless of J_max. Per word: advance offset, load, AND.
 *
 * Pool extended beyond 131 to give more filter bits at 130+ bit N (each prime adds
 * ~1 bit of filter). Only the first ctx->n_primes entries are actually used per run,
 * so storing up to 30 primes doesn't cost runtime unless ctx->n_primes is large. */
#define BS_MAX_NP    30                        /* max extra primes           */
#define BS_MAX_Q     199                       /* largest extra prime        */
static const int BS_POOL[] = {
    41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,101,103,
   107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,
   191,193,197,199
};
#define BS_POOL_N (int)(sizeof(BS_POOL)/sizeof(BS_POOL[0]))

typedef struct {
    /* Layered tiled tables: tile[qi][s] = 64-bit mask for starting offset s.
     * Bit j is set iff (s + j) % q is in the valid set for prime qi.
     * Each prime: q entries × 8 bytes. All primes: ~2.4KB total → L1.        */
    u64  tile[BS_MAX_NP][BS_MAX_Q];
    long Minv[BS_MAX_NP];
    int  primes[BS_MAX_NP];
    int  n_primes;
    int  n_words;                              /* ceil(J_max / 64)           */
    int  J_max;
    /* Capped-range multi-pass scan: only process words in [w_start, w_end).
     * When w_end >= n_words (the default), this scans the entire range.
     * Multi-pass scheme sets these to sub-ranges for early-out on the common
     * case where true K is small. Deterministic: every word is eventually
     * scanned across the sequence of passes. */
    int  w_start, w_end;
    u64  last_word_mask;                       /* mask for bits in final word */
    /* Cross-prime.
     *   cp_uses_mpz = 1 for N ≥ 127-bit (g = 9K²+N overflows u128), else 0 (u128 fast path).
     *   cp_T_init: u128 seed for Newton in u128 mode.
     *   N_mpz_ptr: pointer to the full N mpz used in mpz mode (not owned).
     */
    u64  cp_p, cp_N_mod_p;
    u128 cp_T_init;
    u128 cp_N_128;
    int  cp_uses_mpz;
    const __mpz_struct *cp_N_mpz_ptr;
    /* QR64 tiled table: 64 entries × 1 u64 = 512 bytes.
     * tile_qr64[K0 mod 64] = mask where bit j set iff g(K0+j*M) mod 64 in QR64.
     * Same mask for ALL words (period 64 = word size). ~2.3 extra bits.     */
    u64  tile_qr64[64];
    unsigned long M_mod64;
    /* On-demand screen constants */
    unsigned long N_mod65, N_mod144, N_mod256, N_mod31, N_mod37, N_mod41, N_mod49, N_mod27, N_mod25;
    unsigned long N_mod121, N_mod169;
    /* S-sieve / M-sieve support (new-form M-sieve added) */
    enum sieve_mode mode;
    u64 S_min;  /* offset: actual S (or M) = S_min + CRT residue + j*M */
    u128 cp_n_128;    /* n as u128 for S-sieve cross-prime */
    u64  cp_n_mod_p;  /* n mod cp_p (S-sieve) */
    unsigned long n_mod31, n_mod37, n_mod49, n_mod25, n_mod121, n_mod169;
    unsigned long n_mod256, n_mod65, n_mod41, n_mod27;
    unsigned long n_mod144;   /* n = (N-1)/6 mod 144, for S-sieve QR144 screen */
    /* M-sieve QR144: built at runtime in fused_ctx_init.
     * qr144m[m] = 1 iff g_M(m) = (9m² - N) mod 144 is a square mod 144 AND
     *             satisfies the L ≡ 1 (mod 3) constraint at the true solution.
     * Filter density: exactly 36/144 = 25% pass (2 bits of filtering). */
    uint8_t qr144m[144];
} FusedCtx;

/* Choose optimal number of extra primes for a given J_max.
 * Jacobi primes only in bitset. Survivors go to verify chain (~50ns each). */
static int bs_choose_nprimes(int J_max) {
    int n_words = (J_max + 63) / 64;
    if (n_words < 1) n_words = 1;
    int best = 5;
    double best_cost = 1e99;
    for (int np = 3; np <= BS_POOL_N && np <= 12; np++) {
        double density = 1.0;
        for (int i = 0; i < np; i++)
            density *= ((BS_POOL[i]+1)/2.0) / BS_POOL[i];
        /* QR64 provides additional ~0.188 density factor */
        double combined_density = density * 0.188;
        double survivors_per_word = 64.0 * combined_density;
        /* Per-word: np gathers (~2 cycles each) + AND + QR64 AND + survivor extraction.
         * Survivor handling: ~120 cycles per testK call. */
        double per_word = np * 2.0 + 3.0 + survivors_per_word * 120.0;
        /* Fixed per-K0: batch setup, shift computation per prime */
        double fixed = np * 5.0 + 20.0;
        double cost = fixed + (double)n_words * per_word;
        if (cost < best_cost) { best_cost = cost; best = np; }
    }
    return best;
}

static void fused_ctx_init(FusedCtx *ctx, const mpz_t N, const mpz_t M_mpz,
                           long M_long, long K_max_l, int kpar, int nval8,
                           const uint8_t *QR144_tab, int crt_np, enum sieve_mode mode,
                           u64 S_min_l) {
    long J_max = (M_long > 0) ? (K_max_l / M_long + 1) : 1;
    ctx->J_max   = (int)J_max;
    ctx->n_words = (int)((J_max + 63) / 64);
    /* Default: scan entire range. Multi-pass dispatcher may override. */
    ctx->w_start = 0;
    ctx->w_end   = ctx->n_words;
    /* Bitset prime count: mode-aware default.
     *   K-sieve: crt_np + 2 (legacy default).
     *   S-sieve: crt_np + 3 (optimal at 90-bit; within noise of +4 at 100-bit).
     *   M-sieve: crt_np + 5 (optimal across 90-100 bit sweeps).
     * The mod-3 root fold in S-sieve means its BS_POOL tiles already enjoy
     * ~1.5 bits of extra root-modulus filtering for free — M-sieve has no
     * such fold, so it needs 2-3 more bitset primes to reach equivalent
     * survivor-per-word density. Ratios measured via BS_NP sweep at 90 & 100-bit.
     * Net speedup from this tuning: M-sieve 1.58× faster at 90-bit, S-sieve 1.23×. */
    int bs_np_default;
    switch (mode) {
        case MODE_S_SIEVE: bs_np_default = crt_np + 3; break;
        case MODE_M_SIEVE: bs_np_default = crt_np + 5; break;
        default:           bs_np_default = crt_np + 2; break;  /* K-sieve */
    }
    ctx->n_primes = bs_np_default;
    if (ctx->n_primes < 5) ctx->n_primes = 5;
    if (ctx->n_primes > BS_POOL_N) ctx->n_primes = BS_POOL_N;
    /* KSIEVE_BS_NP env override for empirical tuning / regression testing */
    {
        const char *env = getenv("KSIEVE_BS_NP");
        if (env) {
            int override_np = atoi(env);
            if (override_np >= 3 && override_np <= BS_POOL_N) {
                ctx->n_primes = override_np;
            }
        }
    }

    int tail = (int)(J_max % 64);
    ctx->last_word_mask = (tail == 0) ? ~0ULL : ((1ULL << tail) - 1);

    /* Precompute n mod q for S-sieve and M-sieve polynomials.
     * S-sieve: n = (N-1)/6 (since N = 6n+1).
     * M-sieve: n = (N+1)/6 (since N = 6n-1). M-sieve doesn't use n in the
     *   polynomial g_M(M) = 9M²-N, but we still keep n_val computed for the
     *   Jacobi signal / QR144 precompute path if extended later. */
    mpz_t n_val;
    mpz_init(n_val);
    if (mode == MODE_S_SIEVE) {
        mpz_sub_ui(n_val, N, 1);
        mpz_fdiv_q_ui(n_val, n_val, 6);
    } else if (mode == MODE_M_SIEVE) {
        mpz_add_ui(n_val, N, 1);
        mpz_fdiv_q_ui(n_val, n_val, 6);
    }

    for (int qi = 0; qi < ctx->n_primes; qi++) {
        int q = BS_POOL[qi];
        long Mq = (long)mpz_fdiv_ui(M_mpz, (unsigned long)q);

        if (Mq == 0) {
            ctx->primes[qi] = q;
            ctx->Minv[qi] = 0;
            for (int s = 0; s < q; s++)
                ctx->tile[qi][s] = ~0ULL;
            continue;
        }

        ctx->primes[qi] = q;
        long Nq = (long)mpz_fdiv_ui(N, (unsigned long)q);
        long nq = (mode == MODE_S_SIEVE) ? (long)mpz_fdiv_ui(n_val, (unsigned long)q) : 0;
        ctx->Minv[qi] = modinv(Mq, q);
        (void)nq;  /* used only in S-sieve branch below */

        uint8_t valid_j[BS_MAX_Q + 1];
        memset(valid_j, 0, sizeof(valid_j));
        for (int k = 0; k < q; k++) {
            long val;
            if (mode == MODE_K_SIEVE) {
                val = (9LL * k * k + Nq) % q;
            } else if (mode == MODE_S_SIEVE) {
                long Smq = S_min_l % q;
                long sk = (Smq + k) % q;
                val = ((9LL * sk % q * sk % q + q - 6LL * sk % q + q
                        - 6LL * nq % q + q) % q + q) % q;
            } else {
                /* M-sieve: g_M(M_min + k) = 9(M_min+k)² - N (mod q) */
                long Mmq = S_min_l % q;
                long mk = (Mmq + k) % q;
                val = ((9LL * mk % q * mk % q + q - Nq % q) % q + q) % q;
            }
            if (jacobi(val, q) >= 0) {
                int jb = (int)(((long)k * ctx->Minv[qi] % q + q) % q);
                valid_j[jb] = 1;
            }
        }

        for (int s = 0; s < q; s++) {
            u64 mask = 0;
            for (int j = 0; j < 64; j++) {
                if (valid_j[(s + j) % q])
                    mask |= 1ULL << j;
            }
            ctx->tile[qi][s] = mask;
        }
    }

    /* ── QR64 tiled table ────────────────────────────────────────────────
     * For each K0/S0/M0 mod 64: precompute which of 64 stride positions have
     * g mod 64 in QR64. Same mask applies to ALL words (period = word size). */
    {
        static const uint8_t _QR64[64] = {
            1,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,
            0,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0};
        ctx->M_mod64 = (unsigned long)mpz_fdiv_ui(M_mpz, 64);
        unsigned long Nm64 = (unsigned long)mpz_fdiv_ui(N, 64);
        unsigned long nm64 = (mode == MODE_S_SIEVE)
                              ? (unsigned long)mpz_fdiv_ui(n_val, 64) : 0;
        unsigned long Sm64 = (unsigned long)(S_min_l & 63);
        for (int k0m = 0; k0m < 64; k0m++) {
            u64 mask = 0;
            for (int j = 0; j < 64; j++) {
                unsigned long K64 = ((unsigned long)k0m + (unsigned long)j * ctx->M_mod64) & 63;
                unsigned long g64;
                if (mode == MODE_K_SIEVE) {
                    g64 = (9UL * K64 * K64 + Nm64) & 63;
                } else if (mode == MODE_S_SIEVE) {
                    /* actual S mod 64 = (S_min + K64) mod 64 */
                    unsigned long S64 = (Sm64 + K64) & 63;
                    g64 = (9UL * S64 * S64 + 64 - (6UL * S64 & 63) + 64
                           - (6UL * nm64 & 63)) & 63;
                } else {
                    /* M-sieve: g_M(M) = 9M² - N mod 64, with M = M_min + K64 */
                    unsigned long M64 = (Sm64 + K64) & 63;
                    g64 = (9UL * M64 * M64 + 64 - Nm64) & 63;
                }
                if (_QR64[g64])
                    mask |= 1ULL << j;
            }
            ctx->tile_qr64[k0m] = mask;
        }
    }

    mpz_clear(n_val);

    /* Cross-prime setup.
     *   cp_p: for N < 128-bit, use nextprime(√N); for N ≥ 128-bit, cap at nextprime(2^62)
     *         so cp_p fits u64. Loses ~2 bits vs full √N prime but preserves the filter.
     *   cp_uses_mpz: true for N ≥ 127-bit, where g = 9K²+N overflows u128 and the isqrt
     *                must be computed with mpz. False otherwise (u128 fast path).
     *   cp_T_init: u128 seed for Newton in u128 mode. Ignored in mpz mode. */
    ctx->cp_p = 0;
    ctx->cp_T_init = 0;
    ctx->cp_uses_mpz = (mpz_sizeinbase(N, 2) >= 127) ? 1 : 0;
    ctx->cp_N_mpz_ptr = (const __mpz_struct *)N;  /* borrow; N outlives ctx */
    {
        mpz_t p_mpz, sqrtN_mpz;
        mpz_inits(p_mpz, sqrtN_mpz, NULL);
        mpz_sqrt(sqrtN_mpz, N);
        /* Seed cp_T_init as u128 (used only in u128 mode) */
        if (!ctx->cp_uses_mpz && mpz_fits_ulong_p(sqrtN_mpz)) {
            ctx->cp_T_init = (u128)mpz_get_ui(sqrtN_mpz);
        }
        mpz_add_ui(p_mpz, sqrtN_mpz, 1);
        mpz_nextprime(p_mpz, p_mpz);
        if (!mpz_fits_ulong_p(p_mpz)) {
            /* √N exceeds u64 — use nextprime(2^62) for u64-fitting filter */
            mpz_set_ui(p_mpz, 1UL << 62);
            mpz_nextprime(p_mpz, p_mpz);
        }
        if (mpz_fits_ulong_p(p_mpz)) {
            ctx->cp_p = mpz_get_ui(p_mpz);
            ctx->cp_N_mod_p = (u64)mpz_fdiv_ui(N, (unsigned long)ctx->cp_p);
        }
        mpz_clears(p_mpz, sqrtN_mpz, NULL);
        mpz_t nlo, nhi; mpz_inits(nlo, nhi, NULL);
        mpz_tdiv_r_2exp(nlo, N, 64); mpz_tdiv_q_2exp(nhi, N, 64);
        ctx->cp_N_128 = ((u128)mpz_get_ui(nhi) << 64) | mpz_get_ui(nlo);
        mpz_clears(nlo, nhi, NULL);
    }

    /* S-sieve cross-prime: precompute n = (N-1)/6 as u128 and n mod cp_p.
     * M-sieve: n = (N+1)/6, but the M-sieve cross-prime uses g_M = 9M²-N,
     * which doesn't need n. So we don't precompute n for M-sieve here. */
    if (mode == MODE_S_SIEVE) {
        mpz_t nv2, nlo2, nhi2;
        mpz_inits(nv2, nlo2, nhi2, NULL);
        mpz_sub_ui(nv2, N, 1); mpz_fdiv_q_ui(nv2, nv2, 6);
        mpz_tdiv_r_2exp(nlo2, nv2, 64); mpz_tdiv_q_2exp(nhi2, nv2, 64);
        ctx->cp_n_128 = ((u128)mpz_get_ui(nhi2) << 64) | mpz_get_ui(nlo2);
        if (ctx->cp_p)
            ctx->cp_n_mod_p = (u64)mpz_fdiv_ui(nv2, (unsigned long)ctx->cp_p);
        mpz_clears(nv2, nlo2, nhi2, NULL);
    }

    /* On-demand screen constants */
    ctx->N_mod65  = (unsigned long)mpz_fdiv_ui(N, 65536);
    ctx->N_mod256 = (unsigned long)mpz_fdiv_ui(N, 256);
    ctx->N_mod144 = (unsigned long)mpz_fdiv_ui(N, 144);
    ctx->N_mod31  = (unsigned long)mpz_fdiv_ui(N, 31);
    ctx->N_mod37  = (unsigned long)mpz_fdiv_ui(N, 37);
    ctx->N_mod41  = (unsigned long)mpz_fdiv_ui(N, 41);
    ctx->N_mod49  = (unsigned long)mpz_fdiv_ui(N, 49);
    ctx->N_mod27  = (unsigned long)mpz_fdiv_ui(N, 27);
    ctx->N_mod25  = (unsigned long)mpz_fdiv_ui(N, 25);
    ctx->N_mod121 = (unsigned long)mpz_fdiv_ui(N, 121);
    ctx->N_mod169 = (unsigned long)mpz_fdiv_ui(N, 169);

    /* S-sieve: precompute n mod m for verify chain polynomial.
     * M-sieve: doesn't need n_mod_* in the test_K chain (polynomial is 9M²-N, no n). */
    ctx->mode = mode;
    if (mode == MODE_S_SIEVE) {
        mpz_t nv; mpz_init(nv);
        mpz_sub_ui(nv, N, 1); mpz_fdiv_q_ui(nv, nv, 6);
        ctx->n_mod31  = (unsigned long)mpz_fdiv_ui(nv, 31);
        ctx->n_mod37  = (unsigned long)mpz_fdiv_ui(nv, 37);
        ctx->n_mod49  = (unsigned long)mpz_fdiv_ui(nv, 49);
        ctx->n_mod25  = (unsigned long)mpz_fdiv_ui(nv, 25);
        ctx->n_mod121 = (unsigned long)mpz_fdiv_ui(nv, 121);
        ctx->n_mod169 = (unsigned long)mpz_fdiv_ui(nv, 169);
        ctx->n_mod256 = (unsigned long)mpz_fdiv_ui(nv, 256);
        ctx->n_mod65  = (unsigned long)mpz_fdiv_ui(nv, 65536);
        ctx->n_mod41  = (unsigned long)mpz_fdiv_ui(nv, 41);
        ctx->n_mod27  = (unsigned long)mpz_fdiv_ui(nv, 27);
        ctx->n_mod144 = (unsigned long)mpz_fdiv_ui(nv, 144);
        mpz_clear(nv);
    } else {
        /* K-sieve and M-sieve: zero these out to be safe; they're unused in those paths. */
        ctx->n_mod31 = ctx->n_mod37 = ctx->n_mod49 = ctx->n_mod25 = 0;
        ctx->n_mod121 = ctx->n_mod169 = ctx->n_mod256 = ctx->n_mod65 = 0;
        ctx->n_mod41 = ctx->n_mod27 = ctx->n_mod144 = 0;
    }
    /* K-sieve resets S_min to 0 (K scanned from 0). S-sieve and M-sieve use
     * the scan variable offset from S_min_l (or M_min). */
    if (mode == MODE_K_SIEVE) ctx->S_min = 0;
}

/* Test a single K or S candidate: on-demand screens → CP → isqrt.
 * In S-sieve or M-sieve mode, K_long actually holds S or M, and the
 * polynomial is g_S = 9S²-6S-6n (S) or g_M = 9M²-N (M) instead of g_K = 9K²+N. */
static int fused_test_K(u64 K_long, const FusedCtx *ctx, const mpz_t N,
                        int kpar, int nval8, const uint8_t *QR144_tab,
                        long *n_isqrt, mpz_t uo, mpz_t vo) {
    INSTR_START(testK);
    #define _TESTK_REJECT() do { INSTR_END_VAR(testK_fail, testK); return 0; } while (0)
#ifdef INSTRUMENT
    #define _TESTK_REJECT_AT(ctr) do { g_instr.reject_##ctr++; _TESTK_REJECT(); } while (0)
#else
    #define _TESTK_REJECT_AT(ctr) _TESTK_REJECT()
#endif

    /* K-sieve polynomial: g = (9k² + N) mod m */
    #define _G_CHK_K(K_, Nm_, tbl_, m_) ({                                   \
        unsigned long _k = (unsigned long)(((K_) % (m_) + (m_)) % (m_));     \
        unsigned long _g = (9UL * (_k * _k % (m_)) % (m_) + (Nm_)) % (m_); \
        (tbl_)[_g]; })

    /* S-sieve polynomial: g = (9s² - 6s - 6n) mod m
     * Rewritten as (9s² + m - 6s%m + m - 6n%m) % m to stay positive. */
    #define _G_CHK_S(S_, nm_, tbl_, m_) ({                                           \
        unsigned long _s = (unsigned long)(((S_) % (m_) + (m_)) % (m_));             \
        unsigned long _g = (9UL * (_s * _s % (m_)) % (m_)                            \
                            + (m_) - (6UL * _s) % (m_)                               \
                            + (m_) - (6UL * (nm_)) % (m_)) % (m_);                   \
        (tbl_)[_g]; })

    /* M-sieve polynomial: g = (9m² - N) mod `mod_`. Kept positive via +m_ offset. */
    #define _G_CHK_M(M_, Nm_, tbl_, m_) ({                                   \
        unsigned long _m = (unsigned long)(((M_) % (m_) + (m_)) % (m_));     \
        unsigned long _g = (9UL * (_m * _m % (m_)) % (m_)                    \
                            + (m_) - (Nm_) % (m_)) % (m_);                   \
        (tbl_)[_g]; })

    if (ctx->mode == MODE_K_SIEVE) {
        /* --- K-sieve filter chain --- */
        { unsigned long _K256 = (unsigned long)(((K_long) % 256 + 256) % 256);
          unsigned long _g256 = (9UL * (_K256 * _K256 % 256) % 256 + ctx->N_mod256) % 256;
          if (!QR256[_g256]) _TESTK_REJECT_AT(qr256); }
        if (!_G_CHK_K(K_long, ctx->N_mod31, QR31, 31)) _TESTK_REJECT_AT(qr31);
        if (!_G_CHK_K(K_long, ctx->N_mod37, QR37, 37)) _TESTK_REJECT_AT(qr37);
        if (!_G_CHK_K(K_long, ctx->N_mod65, QR65536, 65536)) _TESTK_REJECT_AT(qr65536);
        if (!_G_CHK_K(K_long, ctx->N_mod49, QR49, 49)) _TESTK_REJECT_AT(qr49);
        if (!_G_CHK_K(K_long, ctx->N_mod25, QR25, 25)) _TESTK_REJECT_AT(qr25);
        if (!_G_CHK_K(K_long, ctx->N_mod121, QR121, 121)) _TESTK_REJECT_AT(qr121);
        if (!_G_CHK_K(K_long, ctx->N_mod169, QR169, 169)) _TESTK_REJECT_AT(qr169);
        if (!_G_CHK_K(K_long, ctx->N_mod144, QR144_tab, 144)) _TESTK_REJECT_AT(qr144);
        { int K8=(int)(K_long&7); if(((nval8*nval8-K8*K8+64)&7)!=0) _TESTK_REJECT_AT(mod8); }
        if (!_G_CHK_K(K_long, ctx->N_mod41, QR41, 41)) _TESTK_REJECT_AT(qr41);
        if (!_G_CHK_K(K_long, ctx->N_mod27, QR27, 27)) _TESTK_REJECT_AT(qr27);
    } else if (ctx->mode == MODE_S_SIEVE) {
        /* --- S-sieve filter chain --- */
        { unsigned long _S256 = (unsigned long)(((K_long) % 256 + 256) % 256);
          unsigned long _g256 = (9UL * (_S256 * _S256 % 256) % 256
                                 + 256 - (6UL * _S256) % 256
                                 + 256 - (6UL * ctx->n_mod256) % 256) % 256;
          if (!QR256[_g256]) _TESTK_REJECT_AT(qr256); }
        if (!_G_CHK_S(K_long, ctx->n_mod31, QR31, 31)) _TESTK_REJECT_AT(qr31);
        if (!_G_CHK_S(K_long, ctx->n_mod37, QR37, 37)) _TESTK_REJECT_AT(qr37);
        if (!_G_CHK_S(K_long, ctx->n_mod65, QR65536, 65536)) _TESTK_REJECT_AT(qr65536);
        if (!_G_CHK_S(K_long, ctx->n_mod49, QR49, 49)) _TESTK_REJECT_AT(qr49);
        if (!_G_CHK_S(K_long, ctx->n_mod25, QR25, 25)) _TESTK_REJECT_AT(qr25);
        if (!_G_CHK_S(K_long, ctx->n_mod121, QR121, 121)) _TESTK_REJECT_AT(qr121);
        if (!_G_CHK_S(K_long, ctx->n_mod169, QR169, 169)) _TESTK_REJECT_AT(qr169);
        /* S-sieve QR144: pre-screen g_S mod 144. Uses QR144S_tab selected by
         * N%8 in sieve_build. At true S, g_S = (3K)² with 3K mod 12 matching
         * the N%8 constraint → only {0}, {9,81}, {36}, or {9,81} are valid.
         * Density of passing S mod 144: 12/144 (~3.6 raw bits). */
        if (!_G_CHK_S(K_long, ctx->n_mod144, QR144_tab, 144)) _TESTK_REJECT_AT(qr144);
        if (!_G_CHK_S(K_long, ctx->n_mod41, QR41, 41)) _TESTK_REJECT_AT(qr41);
        if (!_G_CHK_S(K_long, ctx->n_mod27, QR27, 27)) _TESTK_REJECT_AT(qr27);
    } else {
        /* --- M-sieve filter chain ---
         * Polynomial g_M(M) = 9M² - N where K_long holds the actual M value.
         *
         * Order determined by instrumented rejection counts at 100-bit M-sieve:
         *   QR31 (~23M rejections, strongest), QR37 (~11M), QR256 (~8.6M),
         *   QR25 (~2.7M), QR49 (~2.3M), QR169 (~0.7M), QR65536 (~0.7M).
         * Filters with ~0 rejections at this bit size have been dropped from the
         * hot path: QR41, QR27, QR121, QR144. (QR144M is provably correct but
         * fully subsumed by QR9 ∩ QR16 factors of QR49, QR25, QR256.) */
        if (!_G_CHK_M(K_long, ctx->N_mod31, QR31, 31)) _TESTK_REJECT_AT(qr31);
        if (!_G_CHK_M(K_long, ctx->N_mod37, QR37, 37)) _TESTK_REJECT_AT(qr37);
        { unsigned long _M256 = (unsigned long)(((K_long) % 256 + 256) % 256);
          unsigned long _g256 = (9UL * (_M256 * _M256 % 256) % 256
                                 + 256 - ctx->N_mod256) % 256;
          if (!QR256[_g256]) _TESTK_REJECT_AT(qr256); }
        if (!_G_CHK_M(K_long, ctx->N_mod25, QR25, 25)) _TESTK_REJECT_AT(qr25);
        if (!_G_CHK_M(K_long, ctx->N_mod49, QR49, 49)) _TESTK_REJECT_AT(qr49);
        if (!_G_CHK_M(K_long, ctx->N_mod169, QR169, 169)) _TESTK_REJECT_AT(qr169);
        if (!_G_CHK_M(K_long, ctx->N_mod65, QR65536, 65536)) _TESTK_REJECT_AT(qr65536);
        (void)nval8;  /* unused in M-sieve (no mod-8 identity on n²-K²) */
        (void)QR144_tab;  /* M-sieve QR144 tables generated but redundant in hot path */
    }
    #undef _G_CHK_K
    #undef _G_CHK_S
    #undef _G_CHK_M
    #undef _TESTK_REJECT_AT
    #undef _TESTK_REJECT

    /* Cross-prime filter: T² ≡ g (mod cp_p) where T = isqrt(g) in integers.
     * At N < 127-bit, g fits u128 → fast u128 Newton path.
     * At N ≥ 127-bit, g overflows u128 → mpz path (slower per call but preserves the filter).
     * Either way, gives ~log₂(cp_p) bits of filtering (62+ bits). */
    if (ctx->cp_p) {
        if (!ctx->cp_uses_mpz) {
            /* ── u128 fast path (N ≤ 126-bit) ─────────────────────────────── */
            if (ctx->mode == MODE_K_SIEVE) {
                /* K-sieve: g = 9K² + N. Check √g mod p. */
                u64 Ku = K_long;
                u128 cpg = (u128)9 * Ku * Ku + ctx->cp_N_128;
                u64 cpT = (u64)ctx->cp_T_init
                        + (u64)((u128)9 * Ku * Ku / (2 * ctx->cp_T_init));
                for (int i = 0; i < 10; i++) {
                    u128 sq = (u128)cpT * cpT;
                    u64 nxt = (sq > cpg) ? cpT - (u64)((sq-cpg)/(2*cpT)) - 1
                                         : cpT + (u64)((cpg-sq)/(2*cpT));
                    if (nxt == cpT) break;
                    cpT = nxt;
                }
                if ((u128)(cpT+1)*(cpT+1) <= cpg) cpT++;
                if ((u128)cpT*cpT > cpg) cpT--;
                u64 Tm = cpT % ctx->cp_p, Km = Ku % ctx->cp_p;
                u64 gm = (u64)(((u128)9*Km%ctx->cp_p*Km+ctx->cp_N_mod_p)%ctx->cp_p);
                if ((u64)(((u128)Tm*Tm)%ctx->cp_p) != gm) {
#ifdef INSTRUMENT
                    g_instr.reject_cp++;
#endif
                    INSTR_END_VAR(testK_fail, testK); return 0;
                }
            } else if (ctx->mode == MODE_S_SIEVE) {
                /* S-sieve: g_S = 9S² - 6S - 6n = (3K)². Check √(g_S) mod p. */
                u64 Su = (u64)K_long;  /* K_long holds actual S (always positive) */
                u128 cpg = (u128)9 * Su * Su - (u128)6 * Su - (u128)6 * ctx->cp_n_128;
                int bits = 0;
                { u128 tmp = cpg; while (tmp) { bits++; tmp >>= 1; } }
                u64 cpR = (bits <= 1) ? 1 : (u64)1 << ((bits + 1) / 2);
                for (int i = 0; i < 20; i++) {
                    if (cpR == 0) { cpR = 1; break; }
                    u128 sq = (u128)cpR * cpR;
                    u64 nxt = (u64)((sq + cpg) / (2 * (u128)cpR));
                    if (nxt >= cpR) nxt = cpR - 1;
                    if (nxt == cpR || nxt == cpR - 1) break;
                    cpR = nxt;
                }
                while ((u128)(cpR+1)*(cpR+1) <= cpg) cpR++;
                while ((u128)cpR*cpR > cpg) cpR--;
                u64 Rm = cpR % ctx->cp_p;
                u64 Sm = Su % ctx->cp_p;
                u64 nm = ctx->cp_n_mod_p;
                u64 gm = (u64)(((u128)9*Sm%ctx->cp_p*Sm%ctx->cp_p
                                + ctx->cp_p - (u128)6*Sm%ctx->cp_p
                                + ctx->cp_p - (u128)6*nm%ctx->cp_p) % ctx->cp_p);
                if ((u64)(((u128)Rm*Rm)%ctx->cp_p) != gm) {
#ifdef INSTRUMENT
                    g_instr.reject_cp++;
#endif
                    INSTR_END_VAR(testK_fail, testK); return 0;
                }
            } else {
                /* M-sieve: g_M = 9M² - N = L². Check √(g_M) mod p.
                 * K_long holds actual M. Need isqrt in u128 like the S-sieve branch
                 * (g_M is positive unsigned within range). */
                u64 Mu = (u64)K_long;
                u128 ninemsq = (u128)9 * Mu * Mu;
                if (ninemsq < ctx->cp_N_128) {
                    /* g would be negative — M below M_min or floating-point slop. Reject. */
#ifdef INSTRUMENT
                    g_instr.reject_cp++;
#endif
                    INSTR_END_VAR(testK_fail, testK); return 0;
                }
                u128 cpg = ninemsq - ctx->cp_N_128;
                int bits = 0;
                { u128 tmp = cpg; while (tmp) { bits++; tmp >>= 1; } }
                u64 cpR = (bits <= 1) ? 1 : (u64)1 << ((bits + 1) / 2);
                for (int i = 0; i < 20; i++) {
                    if (cpR == 0) { cpR = 1; break; }
                    u128 sq = (u128)cpR * cpR;
                    u64 nxt = (u64)((sq + cpg) / (2 * (u128)cpR));
                    if (nxt >= cpR) nxt = cpR - 1;
                    if (nxt == cpR || nxt == cpR - 1) break;
                    cpR = nxt;
                }
                while ((u128)(cpR+1)*(cpR+1) <= cpg) cpR++;
                while ((u128)cpR*cpR > cpg) cpR--;
                u64 Lm = cpR % ctx->cp_p;
                u64 Mm = Mu % ctx->cp_p;
                /* gm = (9*Mm² - N_mod_p) mod cp_p, kept positive */
                u64 gm = (u64)(((u128)9*Mm%ctx->cp_p*Mm%ctx->cp_p
                                + ctx->cp_p - ctx->cp_N_mod_p % ctx->cp_p) % ctx->cp_p);
                if ((u64)(((u128)Lm*Lm)%ctx->cp_p) != gm) {
#ifdef INSTRUMENT
                    g_instr.reject_cp++;
#endif
                    INSTR_END_VAR(testK_fail, testK); return 0;
                }
            }
        } else {
            /* ── mpz path (N ≥ 127-bit, g overflows u128) ──────────────────
             * Uses thread-local scratch mpz_t's (cp_scratch_*) to avoid heap
             * allocations per call. */
            cp_scratch_ensure();
            mpz_ptr g = cp_scratch_g, T = cp_scratch_T, Kmpz = cp_scratch_K;
            if (ctx->mode == MODE_K_SIEVE) {
                /* K-sieve: g = 9K² + N, T = isqrt(g), check T² ≡ g (mod cp_p). */
                mpz_set_ui(Kmpz, (unsigned long)K_long);
                mpz_mul(g, Kmpz, Kmpz);          /* g = K² */
                mpz_mul_ui(g, g, 9);             /* g = 9K² */
                mpz_add(g, g, ctx->cp_N_mpz_ptr); /* g = 9K² + N */
                mpz_sqrt(T, g);                  /* T = isqrt(g) */
                u64 Tm = (u64)mpz_fdiv_ui(T, (unsigned long)ctx->cp_p);
                u64 Km = (u64)K_long % ctx->cp_p;
                u64 gm = (u64)(((u128)9*Km%ctx->cp_p*Km+ctx->cp_N_mod_p)%ctx->cp_p);
                if ((u64)(((u128)Tm*Tm)%ctx->cp_p) != gm) {
#ifdef INSTRUMENT
                    g_instr.reject_cp++;
#endif
                    INSTR_END_VAR(testK_fail, testK); return 0;
                }
            } else if (ctx->mode == MODE_S_SIEVE) {
                /* S-sieve: g_S = 9S² - 6S - 6n = 9S² - 6S - (N-1).
                 * T = isqrt(g_S), check T² ≡ g_S (mod cp_p). */
                mpz_set_ui(Kmpz, (unsigned long)K_long);  /* Kmpz = S */
                mpz_mul(g, Kmpz, Kmpz);                    /* g = S² */
                mpz_mul_ui(g, g, 9);                       /* g = 9S² */
                mpz_submul_ui(g, Kmpz, 6);                 /* g -= 6S */
                mpz_sub(g, g, ctx->cp_N_mpz_ptr);          /* g -= N */
                mpz_add_ui(g, g, 1);                       /* g += 1, so g -= (N-1) = 6n */
                mpz_sqrt(T, g);
                u64 Tm = (u64)mpz_fdiv_ui(T, (unsigned long)ctx->cp_p);
                u64 Sm = (u64)K_long % ctx->cp_p;
                u64 nm = ctx->cp_n_mod_p;
                u64 gm = (u64)(((u128)9*Sm%ctx->cp_p*Sm%ctx->cp_p
                                + ctx->cp_p - (u128)6*Sm%ctx->cp_p
                                + ctx->cp_p - (u128)6*nm%ctx->cp_p) % ctx->cp_p);
                if ((u64)(((u128)Tm*Tm)%ctx->cp_p) != gm) {
#ifdef INSTRUMENT
                    g_instr.reject_cp++;
#endif
                    INSTR_END_VAR(testK_fail, testK); return 0;
                }
            } else {
                /* M-sieve: g_M = 9M² - N, L = isqrt(g_M), check L² ≡ g_M (mod cp_p). */
                mpz_set_ui(Kmpz, (unsigned long)K_long);  /* Kmpz = M */
                mpz_mul(g, Kmpz, Kmpz);                    /* g = M² */
                mpz_mul_ui(g, g, 9);                       /* g = 9M² */
                mpz_sub(g, g, ctx->cp_N_mpz_ptr);          /* g -= N */
                if (mpz_sgn(g) < 0) {
                    /* M below M_min or bad residue — g would be negative; reject */
#ifdef INSTRUMENT
                    g_instr.reject_cp++;
#endif
                    INSTR_END_VAR(testK_fail, testK); return 0;
                }
                mpz_sqrt(T, g);                             /* T = isqrt(g_M) = L */
                u64 Lm = (u64)mpz_fdiv_ui(T, (unsigned long)ctx->cp_p);
                u64 Mm = (u64)K_long % ctx->cp_p;
                u64 gm = (u64)(((u128)9*Mm%ctx->cp_p*Mm%ctx->cp_p
                                + ctx->cp_p - ctx->cp_N_mod_p % ctx->cp_p) % ctx->cp_p);
                if ((u64)(((u128)Lm*Lm)%ctx->cp_p) != gm) {
#ifdef INSTRUMENT
                    g_instr.reject_cp++;
#endif
                    INSTR_END_VAR(testK_fail, testK); return 0;
                }
            }
        }
    }

    /* Jacobi(N,K) — only K-sieve has this filter (Jac(N, K) ≥ 0 at true K). */
    if (ctx->mode == MODE_K_SIEVE) {
        mpz_t K_tmp; mpz_init(K_tmp);
        mpz_set_ui(K_tmp, K_long);
        if (kpar == 1 && mpz_jacobi(N, K_tmp) < 0) {
            mpz_clear(K_tmp);
#ifdef INSTRUMENT
            g_instr.reject_jac++;
#endif
            INSTR_END_VAR(testK_fail, testK); return 0;
        }
        mpz_clear(K_tmp);
    }

    /* GMP isqrt + factor recovery */
    (*n_isqrt)++;
    mpz_t g, T, T2, threeK, K;
    mpz_inits(g, T, T2, threeK, K, NULL);

    if (ctx->mode == MODE_K_SIEVE) {
        /* K-sieve: g = 9K² + N */
        mpz_set_ui(K, K_long);
        mpz_mul(g, K, K); mpz_mul_ui(g, g, 9); mpz_add(g, g, N);
        mpz_sqrt(T, g); mpz_mul(T2, T, T);
        int ok = (mpz_cmp(T2, g) == 0);
#ifdef INSTRUMENT
        if (!ok) g_instr.reject_isqrt++;
#endif
        if (ok) {
            mpz_mul_ui(threeK, K, 3);
            mpz_sub(uo, T, threeK); mpz_add(vo, T, threeK);
        }
        mpz_clears(g, T, T2, threeK, K, NULL);
        INSTR_END(testK);
        return ok;
    } else if (ctx->mode == MODE_S_SIEVE) {
        /* S-sieve: g_S = 9S² - 6S - 6n.  sqrt(g_S) = 3K. */
        mpz_t S, n_mpz;
        mpz_inits(S, n_mpz, NULL);
        mpz_set_ui(S, K_long);
        mpz_sub_ui(n_mpz, N, 1); mpz_fdiv_q_ui(n_mpz, n_mpz, 6);

        mpz_mul(g, S, S); mpz_mul_ui(g, g, 9);   /* 9S² */
        mpz_submul_ui(g, S, 6);                    /* -6S */
        mpz_submul_ui(g, n_mpz, 6);                /* -6n */

        mpz_sqrt(T, g); mpz_mul(T2, T, T);         /* T = isqrt(g) = 3K candidate */
        int ok = (mpz_cmp(T2, g) == 0);
        if (ok && mpz_fdiv_ui(T, 3) != 0) ok = 0;  /* 3K must be divisible by 3 */
        if (ok) {
            mpz_fdiv_q_ui(K, T, 3);                 /* K = sqrt(g) / 3 */
            mpz_mul_ui(T, S, 3); mpz_sub_ui(T, T, 1); /* T_real = 3S - 1 */
            mpz_mul_ui(threeK, K, 3);
            mpz_sub(uo, T, threeK);
            mpz_add(vo, T, threeK);
            /* Verify u*v == N */
            mpz_mul(T2, uo, vo);
            if (mpz_cmp(T2, N) != 0) ok = 0;
        }
#ifdef INSTRUMENT
        if (!ok) g_instr.reject_isqrt++;
#endif
        mpz_clears(S, n_mpz, g, T, T2, threeK, K, NULL);
        INSTR_END(testK);
        return ok;
    } else {
        /* M-sieve: g_M = 9M² - N.  isqrt(g_M) = L = |(v-u)/2|.
         *   T = 3M = (u+v)/2 (always positive; K_long holds M)
         *   u = T - L, v = T + L  (or swap; we take the ordering where u < v) */
        mpz_t Mmpz;
        mpz_init(Mmpz);
        mpz_set_ui(Mmpz, K_long);

        mpz_mul(g, Mmpz, Mmpz); mpz_mul_ui(g, g, 9);  /* g = 9M² */
        mpz_sub(g, g, N);                              /* g = 9M² - N */

        int ok = 0;
        if (mpz_sgn(g) >= 0) {
            mpz_sqrt(T, g); mpz_mul(T2, T, T);         /* T = isqrt(g) = L candidate */
            ok = (mpz_cmp(T2, g) == 0);
            if (ok) {
                /* T_real = 3M; L_abs = T (the isqrt). Factor recovery: try u=T_real-L. */
                mpz_t T_real;
                mpz_init(T_real);
                mpz_mul_ui(T_real, Mmpz, 3);            /* T_real = 3M */

                mpz_sub(uo, T_real, T);                 /* u = 3M - L */
                mpz_add(vo, T_real, T);                 /* v = 3M + L */

                /* Order u < v. If L < 0 semantically (v < u in class), swap via sign. */
                if (mpz_sgn(uo) > 0 && mpz_sgn(vo) > 0) {
                    mpz_mul(T2, uo, vo);
                    if (mpz_cmp(T2, N) != 0) ok = 0;
                } else {
                    ok = 0;
                }
                mpz_clear(T_real);
            }
        }
#ifdef INSTRUMENT
        if (!ok) g_instr.reject_isqrt++;
#endif
        mpz_clear(Mmpz);
        mpz_clears(g, T, T2, threeK, K, NULL);
        INSTR_END(testK);
        return ok;
    }
}

/* Process one K0 with the bitset mini-CRT */
/* Layered bitset scan: process one 64-bit word at a time.
 * Per word: advance each prime's offset, load tiled mask, AND.
 * All table accesses hit L1 (~2.4KB total for 5 primes + 2KB QR256). */
static int fused_bitset_scan(long K0, long M_long, long K_max_l,
                             const FusedCtx *ctx, const mpz_t N,
                             int kpar, int nval8, const uint8_t *QR144_tab,
                             long *n_isqrt, mpz_t uo, mpz_t vo) {
    int nw = ctx->n_words;
    int np = ctx->n_primes;
    /* Capped-range: only scan words in [w_start, w_end). Defaults to full range. */
    int w_start = ctx->w_start;
    int w_end   = ctx->w_end;
    if (w_end > nw) w_end = nw;
    if (w_start >= w_end) return 0;

    /* Compute initial shift for each Jacobi prime */
    int shift[BS_MAX_NP];
    for (int qi = 0; qi < np; qi++) {
        int q = ctx->primes[qi];
        shift[qi] = (int)((((long long)K0 % q + q) % q * ctx->Minv[qi]) % q);
        if (shift[qi] < 0) shift[qi] += q;
    }

    /* QR256 initial offset: K0 mod 256 */
    int qr_off = (int)(((K0 % 256) + 256) % 256);

    /* Precompute advance-by-64 offset for each prime (branchless) */
    int adv[BS_MAX_NP];
    for (int qi = 0; qi < np; qi++)
        adv[qi] = 64 % ctx->primes[qi];

    /* QR64 mask: same for all words (period 64 = word size) */
    u64 qr64_mask = ctx->tile_qr64[((unsigned long)K0) & 63];

    /* Advance shifts to w_start. For w_start==0 this is a no-op. */
    for (int w = 0; w < w_start; w++) {
        for (int qi = 0; qi < np; qi++) {
            shift[qi] += adv[qi];
            shift[qi] -= ctx->primes[qi] & -(shift[qi] >= ctx->primes[qi]);
        }
    }

    /* Process one word at a time, only within [w_start, w_end) */
    for (int w = w_start; w < w_end; w++) {
        /* Build combined mask from all Jacobi primes + QR64 */
        u64 mask = qr64_mask;
        for (int qi = 0; qi < np; qi++)
            mask &= ctx->tile[qi][shift[qi]];

        /* Mask off invalid bits in last word */
        if (w == nw - 1)
            mask &= ctx->last_word_mask;

        /* Scan survivors */
        while (mask) {
            int bit = __builtin_ctzll(mask);
            mask &= mask - 1;
            u64 K_long = (u64)K0 + (u64)(w * 64 + bit) * (u64)M_long + ctx->S_min;
            if (K_long == 0) continue;

            if (fused_test_K(K_long, ctx, N, kpar, nval8, QR144_tab, n_isqrt, uo, vo))
                return 1;
        }

        /* Branchless advance offsets for next word (+64 stride positions) */
        for (int qi = 0; qi < np; qi++) {
            shift[qi] += adv[qi];
            shift[qi] -= ctx->primes[qi] & -(shift[qi] >= ctx->primes[qi]);
        }
    }
    return 0;
}

#ifdef __AVX512F__
/* ─────────────────────────────────────────────────────────────────────────
 * AVX-512 K0-BATCHED BITSET SCAN
 * Process 8 K0s in parallel at each word index. Uses _mm512_i64gather_epi64
 * to load 8 tile values (one per K0) per prime per word.
 *
 * Returns index (0..count-1) of K0 that found a factor, or -1.
 * ───────────────────────────────────────────────────────────────────────── */
static int fused_bitset_scan_batch8(const long K0s[8], int count, long M_long,
                                     long K_max_l, const FusedCtx *ctx, const mpz_t N,
                                     int kpar, int nval8, const uint8_t *QR144_tab,
                                     long *n_isqrt, mpz_t uo, mpz_t vo) {
    int nw = ctx->n_words;
    int np = ctx->n_primes;

    /* Per-K0 state: shift[BS_MAX_NP], qr64_mask, active flag */
    __attribute__((aligned(64))) long long shifts[BS_MAX_NP][8];
    __attribute__((aligned(64))) long long adv_ll[BS_MAX_NP];
    u64 qr64_masks[8];
    int active_mask = 0;  /* bit i set if K0[i] is alive */

    for (int i = 0; i < count; i++) {
        long K0 = K0s[i];
        qr64_masks[i] = ctx->tile_qr64[((unsigned long)K0) & 63];
        if (qr64_masks[i]) active_mask |= (1 << i);
        for (int qi = 0; qi < np; qi++) {
            int q = ctx->primes[qi];
            long s = ((((long long)K0 % q + q) % q * ctx->Minv[qi]) % q);
            if (s < 0) s += q;
            shifts[qi][i] = s;
        }
    }
    /* Fill unused lanes with 0 shifts to avoid garbage gathers */
    for (int i = count; i < 8; i++) {
        qr64_masks[i] = 0;
        for (int qi = 0; qi < np; qi++) shifts[qi][i] = 0;
    }
    for (int qi = 0; qi < np; qi++)
        adv_ll[qi] = 64 % ctx->primes[qi];

    if (!active_mask) return -1;

    /* Load qr64_masks into __m512i */
    __m512i qr64_vec = _mm512_loadu_si512((const __m512i*)qr64_masks);

    /* Precompute (adv, q) vectors once for the inner loop */
    __m512i adv_vec[BS_MAX_NP], q_vec_arr[BS_MAX_NP];
    for (int qi = 0; qi < np; qi++) {
        adv_vec[qi]   = _mm512_set1_epi64(adv_ll[qi]);
        q_vec_arr[qi] = _mm512_set1_epi64((long long)ctx->primes[qi]);
    }

    for (int w = 0; w < nw; w++) {
        __m512i mask = qr64_vec;

        /* FUSED loop: per prime, load shift -> gather tile -> AND into mask,
         * AND advance shift -> store. Keeping the current shift in register
         * between gather and advance removes a redundant reload and gives
         * the compiler one long dataflow chain per prime to schedule. */
        for (int qi = 0; qi < np; qi++) {
            __m512i sh   = _mm512_load_si512((const __m512i*)shifts[qi]);
            __m512i vals = _mm512_i64gather_epi64(sh, ctx->tile[qi], 8);
            mask         = _mm512_and_si512(mask, vals);

            /* Advance this prime's shift by 64 mod q (branchless wrap). */
            __m512i advanced = _mm512_add_epi64(sh, adv_vec[qi]);
            __mmask8 overflow = _mm512_cmpge_epi64_mask(advanced, q_vec_arr[qi]);
            __m512i advanced_minus_q = _mm512_sub_epi64(advanced, q_vec_arr[qi]);
            __m512i next_sh = _mm512_mask_blend_epi64(overflow, advanced, advanced_minus_q);
            _mm512_store_si512((__m512i*)shifts[qi], next_sh);
        }

        /* Mask off invalid bits in last word (broadcast last_word_mask) */
        if (w == nw - 1) {
            __m512i lwm = _mm512_set1_epi64((long long)ctx->last_word_mask);
            mask = _mm512_and_si512(mask, lwm);
        }

        /* Fast path: skip extraction if all lanes are zero */
        __mmask8 nz = _mm512_test_epi64_mask(mask, mask) & active_mask;
        if (nz) {
            /* Extract 8 u64 masks and scan each non-zero one */
            u64 masks[8] __attribute__((aligned(64)));
            _mm512_store_si512((__m512i*)masks, mask);

            while (nz) {
                int i = __builtin_ctz(nz);
                nz &= nz - 1;
                u64 wmask = masks[i];
                while (wmask) {
                    int bit = __builtin_ctzll(wmask);
                    wmask &= wmask - 1;
                    u64 K_long = (u64)K0s[i] + (u64)(w * 64 + bit) * (u64)M_long + ctx->S_min;
                    if (K_long == 0) continue;

                    if (fused_test_K(K_long, ctx, N, kpar, nval8, QR144_tab,
                                    n_isqrt, uo, vo))
                        return i;
                }
            }
        }
    }
    return -1;
}

/* ─────────────────────────────────────────────────────────────────────────
 * AVX-512 16-WIDE K0-BATCHED SCAN  (opt #10)
 *
 * Doubles batch width to 16 K0s by holding two ZMM per prime (lanes 0..7
 * and 8..15). Two gathers per prime per word. Throughput rationale:
 * _mm512_i64gather_epi64 has latency ~10-15c but throughput ~5c on SPR;
 * issuing two back-to-back gathers on the same tile keeps the gather
 * port better utilised and halves the per-K0 word-scanning overhead.
 *
 * Register budget: per prime we hold (sh_lo, sh_hi, adv, q) = 4 ZMM.
 * For np=5 that's 20 ZMM + mask + qr64_vec + temps = ~26, fits in 32.
 * For np>6 we may spill; batch8 is kept as the fallback path.
 * ───────────────────────────────────────────────────────────────────────── */
static int fused_bitset_scan_batch16(const long K0s[16], int count, long M_long,
                                      long K_max_l, const FusedCtx *ctx, const mpz_t N,
                                      int kpar, int nval8, const uint8_t *QR144_tab,
                                      long *n_isqrt, mpz_t uo, mpz_t vo) {
    int nw = ctx->n_words;
    int np = ctx->n_primes;

    INSTR_START(setup);

    __attribute__((aligned(64))) long long shifts_lo[BS_MAX_NP][8];
    __attribute__((aligned(64))) long long shifts_hi[BS_MAX_NP][8];
    __attribute__((aligned(64))) long long adv_ll[BS_MAX_NP];
    u64 qr64_masks_lo[8], qr64_masks_hi[8];
    int active_lo = 0, active_hi = 0;

    for (int i = 0; i < count; i++) {
        long K0 = K0s[i];
        u64 qm = ctx->tile_qr64[((unsigned long)K0) & 63];
        if (i < 8) {
            qr64_masks_lo[i] = qm;
            if (qm) active_lo |= (1 << i);
        } else {
            qr64_masks_hi[i - 8] = qm;
            if (qm) active_hi |= (1 << (i - 8));
        }
#ifdef INSTRUMENT
        if (qm) g_instr.live_lanes++; else g_instr.dead_lanes++;
#endif
        for (int qi = 0; qi < np; qi++) {
            int q = ctx->primes[qi];
            /* K0 >= 0 in all call sites (CRT DFS only emits positive K_new),
             * so K0 % q >= 0 directly. One expensive mod (K0 range ≤ 2^49),
             * one cheap mod (inputs ≤ q < 131). */
            long k_mod = (long)(K0 % q);
            long s = (k_mod * ctx->Minv[qi]) % q;
            if (i < 8) shifts_lo[qi][i] = s;
            else       shifts_hi[qi][i - 8] = s;
        }
    }
    /* Pad unused lanes */
    for (int i = count; i < 16; i++) {
        if (i < 8) {
            qr64_masks_lo[i] = 0;
            for (int qi = 0; qi < np; qi++) shifts_lo[qi][i] = 0;
        } else {
            qr64_masks_hi[i - 8] = 0;
            for (int qi = 0; qi < np; qi++) shifts_hi[qi][i - 8] = 0;
        }
    }
    for (int qi = 0; qi < np; qi++)
        adv_ll[qi] = 64 % ctx->primes[qi];

    if (!active_lo && !active_hi) {
#ifdef INSTRUMENT
        g_instr.skipped_batches++;
#endif
        return -1;
    }
#ifdef INSTRUMENT
    g_instr.np_bitset_sum += np;
#endif

    __m512i qr64_vec_lo = _mm512_loadu_si512((const __m512i*)qr64_masks_lo);
    __m512i qr64_vec_hi = _mm512_loadu_si512((const __m512i*)qr64_masks_hi);

    __m512i adv_vec[BS_MAX_NP], q_vec_arr[BS_MAX_NP];
    for (int qi = 0; qi < np; qi++) {
        adv_vec[qi]   = _mm512_set1_epi64(adv_ll[qi]);
        q_vec_arr[qi] = _mm512_set1_epi64((long long)ctx->primes[qi]);
    }

    INSTR_END(setup);

    for (int w = 0; w < nw; w++) {
        INSTR_START(word_fast);
        __m512i mask_lo = qr64_vec_lo;
        __m512i mask_hi = qr64_vec_hi;

        /* Fused gather+advance for both halves */
        for (int qi = 0; qi < np; qi++) {
            __m512i sh_lo = _mm512_load_si512((const __m512i*)shifts_lo[qi]);
            __m512i sh_hi = _mm512_load_si512((const __m512i*)shifts_hi[qi]);

            __m512i vals_lo = _mm512_i64gather_epi64(sh_lo, ctx->tile[qi], 8);
            __m512i vals_hi = _mm512_i64gather_epi64(sh_hi, ctx->tile[qi], 8);

            mask_lo = _mm512_and_si512(mask_lo, vals_lo);
            mask_hi = _mm512_and_si512(mask_hi, vals_hi);

            /* Advance both halves */
            __m512i adv_vv = adv_vec[qi];
            __m512i q_vv   = q_vec_arr[qi];

            __m512i a_lo = _mm512_add_epi64(sh_lo, adv_vv);
            __m512i a_hi = _mm512_add_epi64(sh_hi, adv_vv);

            __mmask8 o_lo = _mm512_cmpge_epi64_mask(a_lo, q_vv);
            __mmask8 o_hi = _mm512_cmpge_epi64_mask(a_hi, q_vv);

            __m512i amq_lo = _mm512_sub_epi64(a_lo, q_vv);
            __m512i amq_hi = _mm512_sub_epi64(a_hi, q_vv);

            __m512i n_lo = _mm512_mask_blend_epi64(o_lo, a_lo, amq_lo);
            __m512i n_hi = _mm512_mask_blend_epi64(o_hi, a_hi, amq_hi);

            _mm512_store_si512((__m512i*)shifts_lo[qi], n_lo);
            _mm512_store_si512((__m512i*)shifts_hi[qi], n_hi);
        }

        if (w == nw - 1) {
            __m512i lwm = _mm512_set1_epi64((long long)ctx->last_word_mask);
            mask_lo = _mm512_and_si512(mask_lo, lwm);
            mask_hi = _mm512_and_si512(mask_hi, lwm);
        }

        __mmask8 nz_lo = _mm512_test_epi64_mask(mask_lo, mask_lo) & active_lo;
        __mmask8 nz_hi = _mm512_test_epi64_mask(mask_hi, mask_hi) & active_hi;

        if (nz_lo | nz_hi) {
            u64 masks_lo[8] __attribute__((aligned(64)));
            u64 masks_hi[8] __attribute__((aligned(64)));
            _mm512_store_si512((__m512i*)masks_lo, mask_lo);
            _mm512_store_si512((__m512i*)masks_hi, mask_hi);

            while (nz_lo) {
                int i = __builtin_ctz(nz_lo);
                nz_lo &= nz_lo - 1;
                u64 wmask = masks_lo[i];
                while (wmask) {
                    int bit = __builtin_ctzll(wmask);
                    wmask &= wmask - 1;
                    u64 K_long = (u64)K0s[i] + (u64)(w * 64 + bit) * (u64)M_long + ctx->S_min;
                    if (K_long == 0) continue;
#ifdef INSTRUMENT
                    g_instr.survivors_found++;
#endif
                    if (fused_test_K(K_long, ctx, N, kpar, nval8, QR144_tab,
                                    n_isqrt, uo, vo))
                        { INSTR_END_VAR(word_slow, word_fast); return i; }
                }
            }
            while (nz_hi) {
                int i = __builtin_ctz(nz_hi);
                nz_hi &= nz_hi - 1;
                u64 wmask = masks_hi[i];
                while (wmask) {
                    int bit = __builtin_ctzll(wmask);
                    wmask &= wmask - 1;
                    int lane = i + 8;
                    u64 K_long = (u64)K0s[lane] + (u64)(w * 64 + bit) * (u64)M_long + ctx->S_min;
                    if (K_long == 0) continue;
#ifdef INSTRUMENT
                    g_instr.survivors_found++;
#endif
                    if (fused_test_K(K_long, ctx, N, kpar, nval8, QR144_tab,
                                    n_isqrt, uo, vo))
                        { INSTR_END_VAR(word_slow, word_fast); return lane; }
                }
            }
            INSTR_END_VAR(word_slow, word_fast);
        } else {
            INSTR_END(word_fast);
        }
    }
    return -1;
}
#endif

/* ══════════════════════════════════════════════════════════════════════════
 * THREADING: DFS subtree workers
 *
 * Strategy: run DFS serially to split_depth collecting partial-prefix work
 * items, then distribute prefixes across worker threads. Each worker runs
 * DFS from its prefix down to leaves and batches K0s into batch16 scans.
 * First worker to find the factor sets an atomic flag; others exit on the
 * next batch boundary.
 *
 * Shared read-only state: N, sv (sieve), fctx (FusedCtx), M_long, K_max_l,
 *   inv_cache, QR144_tab, kpar, nval8, split_depth, fnp.
 * Per-worker mutable: K0_batch, n_isqrt counter, u/v result.
 * Cross-thread: atomic `found` flag, one mutex-protected result slot.
 * ══════════════════════════════════════════════════════════════════════════ */
typedef struct {
    /* Prefix state at split_depth */
    u64 K_acc;       /* K ≡ K_acc (mod mod_acc) */
    u64 mod_acc;     /* CRT modulus at split_depth */
} PrefixWork;

typedef struct {
    /* Shared inputs (read-only) */
    mpz_srcptr N;
    const Sieve *sv;
    const FusedCtx *fctx;
    const uint8_t *QR144_tab;
    long M_long, K_max_l;
    int fnp;
    int split_depth;
    int kpar, nval8;
    const long *inv_cache;
    /* Shared mutable */
    const PrefixWork *work;
    _Atomic int *next_work;      /* index of next prefix to grab */
    int n_work;
    _Atomic int *found_flag;
    pthread_mutex_t *result_mtx;
    mpz_ptr result_u;
    mpz_ptr result_v;
    /* Per-worker outputs */
    long n_isqrt;
    long nK0_local;
    long ncrt_local;
} WorkerCtx;

static void *worker_run(void *arg) {
    WorkerCtx *wc = (WorkerCtx *)arg;
    const Sieve *sv   = wc->sv;
    const FusedCtx *fctx = wc->fctx;
    const long *inv_cache = wc->inv_cache;
    int fnp = wc->fnp;
    int split_depth = wc->split_depth;
    long M_long = wc->M_long;
    long K_max_l = wc->K_max_l;

    /* Worker-local result slots. Only copied to shared result on success. */
    mpz_t local_u, local_v;
    mpz_inits(local_u, local_v, NULL);

    /* DFS state - we'll reset for each prefix */
    int  vi      [MAX_NP];
    u64  K_acc   [MAX_NP + 1];
    u64  mod_acc [MAX_NP + 1];
    long K0_batch[16];

    for (;;) {
        /* Early exit if any other worker found the answer */
        if (atomic_load_explicit(wc->found_flag, memory_order_relaxed)) break;

        int idx = atomic_fetch_add_explicit(wc->next_work, 1, memory_order_relaxed);
        if (idx >= wc->n_work) break;

        const PrefixWork *pw = &wc->work[idx];

        /* Set DFS state to the prefix and continue from split_depth */
        K_acc[split_depth]   = pw->K_acc;
        mod_acc[split_depth] = pw->mod_acc;
        for (int i = split_depth; i < fnp; i++) vi[i] = 0;
        int depth = split_depth;
        int batch_count = 0;
        int thread_found = 0;

        while (depth >= split_depth && !thread_found) {
            /* Check global flag periodically (cheap at DFS backtrack) */
            if (vi[depth] >= sv->ps[depth].nvp) {
                vi[depth] = 0;
                if (depth == split_depth) break;  /* done with this prefix */
                depth--;
                vi[depth]++;
                /* Periodic early-exit check */
                if (atomic_load_explicit(wc->found_flag, memory_order_relaxed)) break;
                continue;
            }

            int p  = sv->ps[depth].p;
            int r  = sv->ps[depth].vp[vi[depth]];
            u64 Ka = K_acc[depth];
            u64 ma = mod_acc[depth];

            long Ka_p = (long)(Ka % (u64)p);
            long diff = ((long)(r - Ka_p) % p + p) % p;
            long t    = diff * inv_cache[depth] % p;
            u64 K_new = Ka + ma * (u64)t;
            wc->ncrt_local++;

            if (K_new > (u64)K_max_l) {
                vi[depth]++;
                continue;
            }

            K_acc[depth+1]   = K_new;
            mod_acc[depth+1] = ma * (u64)p;

            if (depth + 1 == fnp) {
                if (K_new > 0) {
                    wc->nK0_local++;
#ifdef __AVX512F__
                    K0_batch[batch_count++] = (long)K_new;
                    if (batch_count == 16) {
                        int hit = fused_bitset_scan_batch16(
                            K0_batch, 16, M_long, K_max_l, fctx, wc->N,
                            wc->kpar, wc->nval8, wc->QR144_tab,
                            &wc->n_isqrt, local_u, local_v);
                        if (hit >= 0) { thread_found = 1; break; }
                        batch_count = 0;
                        /* Check global flag after each batch */
                        if (atomic_load_explicit(wc->found_flag, memory_order_relaxed)) break;
                    }
#else
                    /* No AVX-512: use per-K0 scan. Still thread-safe since
                     * it only reads fctx and writes to worker-local state. */
                    (void)K0_batch; (void)batch_count;
                    if (fused_bitset_scan((long)K_new, M_long, K_max_l, fctx,
                                          wc->N, wc->kpar, wc->nval8,
                                          wc->QR144_tab, &wc->n_isqrt,
                                          local_u, local_v)) {
                        thread_found = 1; break;
                    }
                    if (atomic_load_explicit(wc->found_flag, memory_order_relaxed)) break;
#endif
                }
                vi[depth]++;
                if (vi[depth] >= sv->ps[depth].nvp) {
                    vi[depth] = 0;
                    if (depth == split_depth) break;
                    depth--;
                    vi[depth]++;
                }
            } else {
                depth++;
            }
        }

        /* Flush any leftover K0s for this prefix (only needed for batch16 path) */
#ifdef __AVX512F__
        if (!thread_found && batch_count > 0 &&
            !atomic_load_explicit(wc->found_flag, memory_order_relaxed)) {
            int hit = fused_bitset_scan_batch16(
                K0_batch, batch_count, M_long, K_max_l, fctx, wc->N,
                wc->kpar, wc->nval8, wc->QR144_tab,
                &wc->n_isqrt, local_u, local_v);
            if (hit >= 0) thread_found = 1;
        }
#endif

        if (thread_found) {
            /* Publish result. Use try-lock-first to avoid contention if another
             * thread is already writing. Only first writer wins. */
            int expected = 0;
            if (atomic_compare_exchange_strong_explicit(
                    wc->found_flag, &expected, 1,
                    memory_order_release, memory_order_relaxed)) {
                pthread_mutex_lock(wc->result_mtx);
                mpz_set(wc->result_u, local_u);
                mpz_set(wc->result_v, local_v);
                pthread_mutex_unlock(wc->result_mtx);
            }
            break;
        }
    }

    mpz_clears(local_u, local_v, NULL);
    return NULL;
}

/* Collect all partial-prefix work items at split_depth by running DFS
 * serially from root to split_depth only. Returns count written. */
static int collect_prefixes(const Sieve *sv, const long *inv_cache, int fnp,
                            int split_depth, long K_max_l,
                            PrefixWork *out, int out_cap) {
    int  vi      [MAX_NP];
    u64  K_acc   [MAX_NP + 1];
    u64  mod_acc [MAX_NP + 1];
    K_acc[0] = 0; mod_acc[0] = 1;
    for (int i = 0; i < fnp; i++) vi[i] = 0;

    int depth = 0;
    int count = 0;

    while (depth >= 0) {
        if (vi[depth] >= sv->ps[depth].nvp) {
            vi[depth] = 0;
            depth--;
            if (depth >= 0) vi[depth]++;
            continue;
        }

        int p  = sv->ps[depth].p;
        int r  = sv->ps[depth].vp[vi[depth]];
        u64 Ka = K_acc[depth];
        u64 ma = mod_acc[depth];

        long Ka_p = (long)(Ka % (u64)p);
        long diff = ((long)(r - Ka_p) % p + p) % p;
        long t    = diff * inv_cache[depth] % p;
        u64 K_new = Ka + ma * (u64)t;

        if (K_new > (u64)K_max_l) {
            vi[depth]++;
            continue;
        }

        K_acc[depth+1]   = K_new;
        mod_acc[depth+1] = ma * (u64)p;

        if (depth + 1 == split_depth) {
            /* Record this prefix */
            if (count < out_cap) {
                out[count].K_acc   = K_new;
                out[count].mod_acc = ma * (u64)p;
                count++;
            }
            vi[depth]++;
            if (vi[depth] >= sv->ps[depth].nvp) {
                vi[depth] = 0;
                depth--;
                if (depth >= 0) vi[depth]++;
            }
        } else {
            depth++;
        }
    }

    return count;
}

/* ══════════════════════════════════════════════════════════════════════════
 * TOP-LEVEL FACTORIZATION
 * ══════════════════════════════════════════════════════════════════════════ */
Res ksieve_factor(const mpz_t N, mpz_t uo, mpz_t vo, int verbose) {
    Res res;
    memset(&res, 0, sizeof(res));
    double t0 = now_sec();
    INSTR_START(total);

    res.N_bits = (int)mpz_sizeinbase(N, 2);

    /* ── Mode selection (in priority order) ─────────────────────────────
     *   KSIEVE_RSA=1    → auto-select per N (for RSA-style benchmarks).
     *                      Overrides MSIEVE/SSIEVE if both are set.
     *   KSIEVE_MSIEVE=1 → M-sieve (new form, N ≡ 5 mod 6).
     *   KSIEVE_SSIEVE=1 → S-sieve (old form, N ≡ 1 mod 6, 4.24× range).
     *   KSIEVE_AUTO=1   → auto-select per N (N mod 6 == 5 → M-sieve; else S-sieve).
     *   (none)          → K-sieve (old form, N ≡ 1 mod 6, wide range).
     */
    enum sieve_mode mode = MODE_K_SIEVE;
    u64 S_min_l = 0;
    {
        const char *env_m = getenv("KSIEVE_MSIEVE");
        const char *env_s = getenv("KSIEVE_SSIEVE");
        const char *env_a = getenv("KSIEVE_AUTO");
        const char *env_r = getenv("KSIEVE_RSA");
        int explicit_m = (env_m && atoi(env_m));
        int explicit_s = (env_s && atoi(env_s));
        int explicit_r = (env_r && atoi(env_r));
        /* Priority: RSA > MSIEVE > SSIEVE > AUTO > default (K-sieve).
         *
         * KSIEVE_RSA=1 MUST force per-N sieve selection: RSA-style generators
         * produce a mix of N mod 6 = 1 and N mod 6 = 5 values, and using a
         * fixed sieve for both causes silent factoring failures on half the
         * trials. When KSIEVE_RSA=1 is set alongside MSIEVE or SSIEVE, the
         * RSA per-N behavior takes precedence silently here; run_bench prints
         * a single warning at startup if the combination is detected. */
        if (explicit_r) {
            unsigned long N_mod6 = mpz_fdiv_ui(N, 6);
            mode = (N_mod6 == 5) ? MODE_M_SIEVE : MODE_S_SIEVE;
        }
        else if (explicit_m) mode = MODE_M_SIEVE;
        else if (explicit_s) mode = MODE_S_SIEVE;
        else if (env_a && atoi(env_a)) {
            unsigned long N_mod6 = mpz_fdiv_ui(N, 6);
            mode = (N_mod6 == 5) ? MODE_M_SIEVE : MODE_S_SIEVE;
        }
    }

    /* Sanity-check N mod 6 against the chosen mode. */
    {
        unsigned long N_mod6 = mpz_fdiv_ui(N, 6);
        if (mode == MODE_M_SIEVE && N_mod6 != 5) {
            if (verbose) fprintf(stderr,
                "WARNING: M-sieve requires N ≡ 5 (mod 6); got N mod 6 = %lu. "
                "Factoring will almost certainly fail.\n", N_mod6);
        }
        if (mode != MODE_M_SIEVE && N_mod6 != 1) {
            if (verbose) fprintf(stderr,
                "WARNING: K-sieve / S-sieve requires N ≡ 1 (mod 6); got N mod 6 = %lu. "
                "Factoring will almost certainly fail.\n", N_mod6);
        }
    }

    /* K-sieve: K_max = floor(sqrt(N)/3) + 1 (wide range).
     * S-sieve: [S_min, S_max] where S_min = ceil(sqrt(N)/3), S_max = sqrt(5N/4)/3.
     *   Effective range = S_max - S_min ≈ K_max / 8.5.
     * M-sieve: [M_min, M_max] where M_min = ceil(sqrt(N)/3), M_max = sqrt(5N)/6.
     *   Same range width as S-sieve (verified algebraically).  */
    mpz_t K_max, sqrtN;
    mpz_inits(K_max, sqrtN, NULL);
    mpz_sqrt(sqrtN, N);

    if (mode == MODE_K_SIEVE) {
        mpz_fdiv_q_ui(K_max, sqrtN, 3);
        mpz_add_ui(K_max, K_max, 1);
    } else if (mode == MODE_S_SIEVE) {
        /* S_min = ceil(sqrt(N) / 3) = floor(sqrt(N)/3) + 1 when sqrt(N) not div by 3 */
        mpz_t S_min, S_max, fiveN4;
        mpz_inits(S_min, S_max, fiveN4, NULL);

        mpz_fdiv_q_ui(S_min, sqrtN, 3);
        mpz_add_ui(S_min, S_min, 1);  /* ceil */

        /* S_max = floor(sqrt(5N/4) / 3) */
        mpz_mul_ui(fiveN4, N, 5);
        mpz_fdiv_q_ui(fiveN4, fiveN4, 4);
        mpz_sqrt(S_max, fiveN4);
        mpz_fdiv_q_ui(S_max, S_max, 3);

        S_min_l = (u64)mpz_get_ui(S_min);

        /* Effective search range = S_max - S_min */
        mpz_sub(K_max, S_max, S_min);
        mpz_add_ui(K_max, K_max, 1);

        if (verbose) {
            gmp_printf("S-sieve: S_min=%Zd  S_max=%Zd  range=%Zd\n", S_min, S_max, K_max);
        }

        mpz_clears(S_min, S_max, fiveN4, NULL);
    } else {
        /* M-sieve: M_min = ceil(sqrt(N) / 3), M_max = floor(sqrt(5N) / 6).
         * Range width is (sqrt(5N)/6 - sqrt(N)/3) ≈ 0.039 * sqrt(N), same as S-sieve. */
        mpz_t M_min, M_max, fiveN;
        mpz_inits(M_min, M_max, fiveN, NULL);

        mpz_fdiv_q_ui(M_min, sqrtN, 3);
        mpz_add_ui(M_min, M_min, 1);  /* ceil */

        mpz_mul_ui(fiveN, N, 5);
        mpz_sqrt(M_max, fiveN);
        mpz_fdiv_q_ui(M_max, M_max, 6);

        S_min_l = (u64)mpz_get_ui(M_min);

        mpz_sub(K_max, M_max, M_min);
        mpz_add_ui(K_max, K_max, 1);

        if (verbose) {
            gmp_printf("M-sieve: M_min=%Zd  M_max=%Zd  range=%Zd\n", M_min, M_max, K_max);
        }

        mpz_clears(M_min, M_max, fiveN, NULL);
    }

    res.K_bits = (int)mpz_sizeinbase(K_max, 2);

    long   K_max_l = mpz_get_si(K_max);
    double K_max_d = mpz_get_d(K_max);

    int np = choose_np(res.N_bits, K_max_d, mode);
    /* Allow testing: KSIEVE_NP=NN overrides choose_np decision. */
    {
        const char *env = getenv("KSIEVE_NP");
        if (env) {
            int override_np = atoi(env);
            if (override_np >= 4 && override_np < NSP) {
                np = override_np;
            }
        }
    }
    res.np = np;

    Sieve *sv = (Sieve *)calloc(1, sizeof(Sieve));
    if (!sv) { perror("calloc Sieve"); exit(1); }
    sieve_build(sv, N, np, mode, S_min_l);

    /* Build M from the ACTUAL primes the sieve chose (which may differ from SP[]
     * when adaptive prime selection is enabled). sv->ps[0] is the root; sv->ps[1..np]
     * are the CRT primes after sieve_build sorted them by density. */
    mpz_t M_mpz;
    mpz_init_set_ui(M_mpz, (unsigned long)sv->ps[0].p);
    for (int i = 1; i <= np; i++)
        mpz_mul_ui(M_mpz, M_mpz, (unsigned long)sv->ps[i].p);

    /* Native u128 path is K-sieve-only (polynomial 9K²+N at ≤64-bit N). */
    int use_native = (mode == MODE_K_SIEVE) && (res.N_bits <= 64) && mpz_fits_ulong_p(M_mpz);
    u64  N64    = use_native ? (u64)mpz_get_ui(N)     : 0;
    u64  M64    = use_native ? (u64)mpz_get_ui(M_mpz) : 0;
    long M_long = mpz_fits_ulong_p(M_mpz) ? (long)mpz_get_ui(M_mpz) : 0;

    long J_max = (M_long > 0) ? (K_max_l / M_long + 1) : 1;
    /* S-sieve and M-sieve MUST use fused path — lean stride only supports K-sieve. */
    int  use_fused = !use_native && (J_max > 1) &&
                     (mode != MODE_K_SIEVE || J_max <= (1L << 24));

    if (verbose) {
        int total_np = sv->np;
        const char *range_label =
            (mode == MODE_K_SIEVE) ? "K_max" :
            (mode == MODE_S_SIEVE) ? "S_range" : "M_range";
        const char *mode_label =
            (mode == MODE_K_SIEVE) ? "K-sieve" :
            (mode == MODE_S_SIEVE) ? "S-sieve" : "M-sieve";
        gmp_printf("N   = %d-bit (mode = %s)\n%s = %d-bit  (~2^%.1f)\nnp  = %d (incl. p=2)\nM   = 2^%.1f\nPath: %s\nPrime order (density): ",
               res.N_bits, mode_label, range_label,
               res.K_bits, log2(K_max_d),
               total_np, log2(mpz_get_d(M_mpz)),
               use_native ? "native u128" : (use_fused ? "fused CRT+bitset" : "GMP lean stride"));
        for (int i = 0; i < total_np && i < 12; i++)
            printf("p%d(%.3f)%s", sv->ps[i].p, sv->ps[i].density, i<total_np-1?" ":"");
        if (total_np > 12) printf("...");
        printf("\n");
    }

    res.t_precomp = now_sec() - t0;

    long n_isqrt = 0;

    if (use_fused) {
        /* ── Fused CRT DFS + Bitset scan, threaded ────────────────────── */
        FusedCtx fctx;
        fused_ctx_init(&fctx, N, M_mpz, M_long, K_max_l,
                       sv->k_parity, sv->nval_mod8, sv->QR144, np, mode, S_min_l);
        fctx.S_min = S_min_l;

        double tp = now_sec();

        int fnp = sv->np;
        long inv_cache[MAX_NP];
        {
            u64 running = 1;
            for (int d = 0; d < fnp; d++) {
                int p = sv->ps[d].p;
                inv_cache[d] = (d == 0) ? 1 : modinv((long)(running % (u64)p), p);
                running *= (u64)p;
            }
        }

        /* Choose split_depth: deepest depth such that #prefixes ≥ 2*n_threads
         * for good load balancing. Cap at SPLIT_DEPTH_MAX. At depth d the
         * prefix count is product of nvp[0..d-1] (pruning-dependent, but
         * we use a cheap estimate). */
        int n_threads = KSIEVE_N_THREADS;
        if (n_threads > 16) n_threads = 16;
        if (n_threads < 1)  n_threads = 1;

        int split_depth = 0;
        {
            int estimate = 1;
            int d;
            for (d = 0; d < fnp && d < SPLIT_DEPTH_MAX; d++) {
                estimate *= sv->ps[d].nvp;
                if (estimate >= n_threads * 4) { d++; break; }
            }
            split_depth = d;
            if (split_depth < 1)  split_depth = 1;
            if (split_depth > fnp - 1) split_depth = fnp - 1;
        }

        /* Collect prefixes at split_depth */
        enum { MAX_PREFIXES = 4096 };
        PrefixWork *prefixes = (PrefixWork *)malloc(sizeof(PrefixWork) * MAX_PREFIXES);
        if (!prefixes) { perror("malloc prefixes"); exit(1); }
        int n_prefixes = collect_prefixes(sv, inv_cache, fnp, split_depth,
                                           K_max_l, prefixes, MAX_PREFIXES);

        /* If too few prefixes for threading, collapse to single-threaded */
        if (n_prefixes < 2 || n_threads <= 1) {
            n_threads = 1;
        }

        long total_n_isqrt = 0;
        long total_nK0 = 0;
        long total_ncrt = 0;
        int found_any = 0;

        if (n_prefixes == 0) {
            /* No work after K_max pruning — nothing to do */
            goto no_work;
        }

        _Atomic int found_flag = 0;
        _Atomic int next_work  = 0;
        pthread_mutex_t result_mtx = PTHREAD_MUTEX_INITIALIZER;

        /* Set up worker contexts */
        WorkerCtx wctx[16];
        pthread_t tids[16];
        for (int i = 0; i < n_threads; i++) {
            wctx[i].N          = N;
            wctx[i].sv         = sv;
            wctx[i].fctx       = &fctx;
            wctx[i].QR144_tab  = sv->QR144;
            wctx[i].M_long     = M_long;
            wctx[i].K_max_l    = K_max_l;
            wctx[i].fnp        = fnp;
            wctx[i].split_depth = split_depth;
            wctx[i].kpar       = sv->k_parity;
            wctx[i].nval8      = sv->nval_mod8;
            wctx[i].inv_cache  = inv_cache;
            wctx[i].work       = prefixes;
            wctx[i].next_work  = &next_work;
            wctx[i].n_work     = n_prefixes;
            wctx[i].found_flag = &found_flag;
            wctx[i].result_mtx = &result_mtx;
            wctx[i].result_u   = uo;
            wctx[i].result_v   = vo;
            wctx[i].n_isqrt    = 0;
            wctx[i].nK0_local  = 0;
            wctx[i].ncrt_local = 0;
        }

        /* ── Capped-range multi-pass scan ──────────────────────────────────
         * Empirical observation: at 80-130 bit, K_true (S_true, M_true) lies
         * in the bottom 10% of [0, K_max] with ~50% probability, bottom 32%
         * with ~75%, and bottom 100% always. Scanning the full range every
         * time wastes work on the common case.
         *
         * Strategy: run multiple passes with increasing word-index caps.
         * Each pass scans words [w_{k}, w_{k+1}), and we stop as soon as
         * any worker finds the factor. Every word is eventually scanned
         * across the pass sequence, so correctness is preserved.
         *
         * Pass split points are controlled by fractions f_k of n_words:
         *   KSIEVE_PASS_FRACS="f_0,f_1,...,f_{K-1}"  (K+1 passes, final is 1.0)
         *
         * Default (if not set): 0.10, 0.32 → 3 passes (0, 10%, 32%, 100%).
         * Legacy compatibility: if unset, KSIEVE_PASS0_FRAC and
         * KSIEVE_PASS1_FRAC still work as single-value overrides.
         *
         * Set KSIEVE_SINGLEPASS=1 to force single-pass (regression testing).
         *
         * Example tunings to try:
         *   KSIEVE_PASS_FRACS=0.15                      (2 passes: 0-15%, 15-100%)
         *   KSIEVE_PASS_FRACS=0.05,0.15,0.35            (4 passes)
         *   KSIEVE_PASS_FRACS=0.03,0.08,0.15,0.35       (5 passes)
         */
        const char *env_single = getenv("KSIEVE_SINGLEPASS");
        enum { MAX_PASSES = 10 };
        int n_passes = 0;
        int pass_w_bounds[MAX_PASSES + 1];  /* sentinel at index n_passes */
        if (env_single && atoi(env_single)) {
            /* Regression mode: one pass covering everything. */
            n_passes = 1;
            pass_w_bounds[0] = 0;
            pass_w_bounds[1] = fctx.n_words;
        } else {
            /* Build fraction list from env vars. Priority:
             *   1. KSIEVE_PASS_FRACS (comma-separated list)
             *   2. KSIEVE_PASS0_FRAC + KSIEVE_PASS1_FRAC (legacy 3-pass)
             *   3. Default: 0.10, 0.32
             */
            double fracs[MAX_PASSES];
            int n_fracs = 0;
            const char *env_list = getenv("KSIEVE_PASS_FRACS");
            if (env_list && *env_list) {
                /* Parse comma-separated list. Malformed entries are ignored. */
                const char *p = env_list;
                while (*p && n_fracs < MAX_PASSES) {
                    char *end;
                    double v = strtod(p, &end);
                    if (end == p) break;  /* no number parsed */
                    if (v > 0.0 && v < 1.0) fracs[n_fracs++] = v;
                    if (*end == ',') p = end + 1;
                    else break;
                }
            } else {
                /* Legacy 2-knob defaults */
                double frac0 = 0.10, frac1 = 0.32;
                const char *env0 = getenv("KSIEVE_PASS0_FRAC");
                const char *env1 = getenv("KSIEVE_PASS1_FRAC");
                if (env0) frac0 = atof(env0);
                if (env1) frac1 = atof(env1);
                if (frac0 > 0.0 && frac0 < 1.0) fracs[n_fracs++] = frac0;
                if (frac1 > frac0 && frac1 < 1.0) fracs[n_fracs++] = frac1;
            }
            /* Sort fractions ascending (defensive against misordered input). */
            for (int i = 1; i < n_fracs; i++) {
                double key = fracs[i];
                int j = i - 1;
                while (j >= 0 && fracs[j] > key) { fracs[j+1] = fracs[j]; j--; }
                fracs[j+1] = key;
            }
            /* De-duplicate adjacent entries */
            int n_uniq = 0;
            for (int i = 0; i < n_fracs; i++) {
                if (n_uniq == 0 || fracs[i] > fracs[n_uniq - 1] + 1e-9)
                    fracs[n_uniq++] = fracs[i];
            }
            n_fracs = n_uniq;

            /* Convert fractions to word indices, de-duplicate, and build bounds. */
            int w_cuts[MAX_PASSES];
            int n_cuts = 0;
            for (int i = 0; i < n_fracs; i++) {
                int w = (int)(fctx.n_words * fracs[i]);
                if (w < 1) w = 1;
                if (w >= fctx.n_words) break;  /* subsequent fractions would too */
                if (n_cuts == 0 || w > w_cuts[n_cuts - 1])
                    w_cuts[n_cuts++] = w;
            }
            /* Assemble pass_w_bounds: [0, w_cuts..., n_words] */
            pass_w_bounds[0] = 0;
            for (int i = 0; i < n_cuts; i++) pass_w_bounds[i + 1] = w_cuts[i];
            pass_w_bounds[n_cuts + 1] = fctx.n_words;
            n_passes = n_cuts + 1;

            /* Edge case: too few words or all cuts degenerate → single pass. */
            if (fctx.n_words <= 2 || n_passes < 1) {
                n_passes = 1;
                pass_w_bounds[0] = 0;
                pass_w_bounds[1] = fctx.n_words;
            }
        }

        /* Optional trace of computed pass bounds for debugging */
        if (getenv("KSIEVE_PASS_TRACE")) {
            fprintf(stderr, "[ksieve] n_words=%d, %d pass(es):",
                    fctx.n_words, n_passes);
            for (int i = 0; i <= n_passes; i++)
                fprintf(stderr, " %d", pass_w_bounds[i]);
            fprintf(stderr, "\n");
        }

        for (int pass = 0; pass < n_passes; pass++) {
            /* Set scan bounds on shared fctx. Safe: no workers running yet. */
            fctx.w_start = pass_w_bounds[pass];
            fctx.w_end   = pass_w_bounds[pass + 1];

            /* Reset shared work counter so this pass re-dispatches all prefixes. */
            atomic_store(&next_work, 0);

            if (n_threads == 1) {
                worker_run(&wctx[0]);
            } else {
                for (int i = 1; i < n_threads; i++) {
                    if (pthread_create(&tids[i], NULL, worker_run, &wctx[i]) != 0) {
                        perror("pthread_create"); exit(1);
                    }
                }
                /* Main thread does work too */
                worker_run(&wctx[0]);
                for (int i = 1; i < n_threads; i++)
                    pthread_join(tids[i], NULL);
            }

            /* If any worker found the factor this pass, stop. */
            if (atomic_load(&found_flag)) break;
        }

        found_any = atomic_load(&found_flag);
        for (int i = 0; i < n_threads; i++) {
            total_n_isqrt += wctx[i].n_isqrt;
            total_nK0     += wctx[i].nK0_local;
            total_ncrt    += wctx[i].ncrt_local;
        }

        pthread_mutex_destroy(&result_mtx);

no_work:
        free(prefixes);

        res.success   = found_any;
        res.n_crt_ops = total_ncrt;
        res.n_k0s     = total_nK0;
        res.n_isqrt   = total_n_isqrt;
        res.t_p1      = now_sec() - tp;
        res.t_p2      = 0;
        res.t_total   = now_sec() - t0;

        if (verbose)
            printf("Fused threaded: %d prefixes, %d threads, %ld K0s, %ld isqrt -> %s\n",
                   n_prefixes, n_threads, total_nK0, total_n_isqrt,
                   found_any ? "FOUND" : "not found");

        /* No dynamic table to free — layered tiles are in-struct */

    } else {
        /* ── Legacy path: separate phase 1 + phase 2 ──────────────────── */
        double tp1 = now_sec();
        long nK0 = 0, ncrt = 0;
        long *K0s = phase1_dfs(sv, K_max_l, &nK0, &ncrt);
        res.t_p1      = now_sec() - tp1;
        res.n_crt_ops = ncrt;
        res.n_k0s     = nK0;

        if (verbose)
            printf("Phase 1: %ld CRT ops -> %ld K0 residues\n", ncrt, nK0);

        double tp2 = now_sec();

        ExtraFilter ef;
        memset(&ef, 0, sizeof(ef));
        if (!use_native)
            extra_filter_build(&ef, N, M_mpz, np, J_max);

        if (use_native)
            res.success = phase2_native(K0s, nK0, N64, M64,
                                        (u64)K_max_l, sv->k_parity,
                                        sv->nval_mod8, sv->QR144,
                                        &n_isqrt, uo, vo);
        else
            res.success = phase2_gmp(K0s, nK0, N, M_mpz, K_max,
                                     sv->k_parity, sv->nval_mod8, sv->QR144,
                                     &ef, &n_isqrt, uo, vo);

        res.n_isqrt = n_isqrt;
        res.t_p2    = now_sec() - tp2;
        res.t_total = now_sec() - t0;

        if (verbose)
            printf("Phase 2: %ld isqrt tests -> %s\n",
                   n_isqrt, res.success ? "FOUND" : "not found");

        free(K0s);
    }

    free(sv);
    mpz_clears(K_max, sqrtN, M_mpz, NULL);
    INSTR_END(total);
    return res;
}

/* ══════════════════════════════════════════════════════════════════════════
 * UTILITIES
 * ══════════════════════════════════════════════════════════════════════════ */
static int verify(const mpz_t N, const mpz_t u, const mpz_t v) {
    mpz_t c; mpz_init(c);
    mpz_mul(c, u, v);
    int ok = (mpz_cmp(c, N) == 0);
    mpz_clear(c);
    return ok;
}

static void print_mpz_short(const char *label, const mpz_t x) {
    char *s = mpz_get_str(NULL, 10, x);
    int   n = (int)strlen(s);
    if (n <= 44) printf("%s = %s\n", label, s);
    else         printf("%s = %.*s...%s  (%d digits)\n", label, 14, s, s+n-8, n);
    free(s);
}

static void print_res(const Res *r, int ok) {
    printf("  [%d-bit] np=%d crt=%ld k0=%ld isqrt=%ld"
           " | pre=%.4fs p1=%.4fs p2=%.4fs tot=%.4fs  %s\n",
           r->N_bits, r->np, r->n_crt_ops, r->n_k0s, r->n_isqrt,
           r->t_precomp, r->t_p1, r->t_p2, r->t_total,
           ok ? "OK" : (r->success ? "WRONG!" : "FAIL"));
}

/* ── Semiprime generator ─────────────────────────────────────────────── */
static void gen_semiprime(mpz_t N, mpz_t u, mpz_t v,
                          int bits, gmp_randstate_t rng) {
    int half = bits / 2;
    for (;;) {
        mpz_urandomb(u, rng, half);
        mpz_setbit(u, half - 1);
        mpz_setbit(u, half - 2);
        {
            long adj = ((long)(5 - mpz_fdiv_ui(u, 6)) % 6 + 6) % 6;
            mpz_add_ui(u, u, (unsigned long)adj);
        }
        if (!mpz_probab_prime_p(u, 20)) continue;

        mpz_urandomb(v, rng, half);
        mpz_setbit(v, half - 1);
        {
            long adj = ((long)(5 - mpz_fdiv_ui(v, 6)) % 6 + 6) % 6;
            mpz_add_ui(v, v, (unsigned long)adj);
        }
        if (mpz_cmp(v, u) <= 0) continue;
        if (!mpz_probab_prime_p(v, 20)) continue;

        mpz_mul(N, u, v);
        int nb = (int)mpz_sizeinbase(N, 2);
        if (nb == bits || nb == bits - 1) break;
    }
}

/* Generate a NEW-FORM semiprime: N = u*v with u ≡ 5 (mod 6), v ≡ 1 (mod 6).
 * Produces N ≡ 5 (mod 6), the target for M-sieve. u and v are both large primes
 * of roughly bits/2 bits each. */
static void gen_semiprime_newform(mpz_t N, mpz_t u, mpz_t v,
                                  int bits, gmp_randstate_t rng) {
    int half = bits / 2;
    for (;;) {
        /* u ≡ 5 (mod 6) */
        mpz_urandomb(u, rng, half);
        mpz_setbit(u, half - 1);
        mpz_setbit(u, half - 2);
        {
            long adj = ((long)(5 - mpz_fdiv_ui(u, 6)) % 6 + 6) % 6;
            mpz_add_ui(u, u, (unsigned long)adj);
        }
        if (!mpz_probab_prime_p(u, 20)) continue;

        /* v ≡ 1 (mod 6) */
        mpz_urandomb(v, rng, bits - half);
        mpz_setbit(v, bits - half - 1);
        {
            long adj = ((long)(1 - mpz_fdiv_ui(v, 6)) % 6 + 6) % 6;
            mpz_add_ui(v, v, (unsigned long)adj);
        }
        if (mpz_cmp(v, u) == 0) continue;
        if (!mpz_probab_prime_p(v, 20)) continue;

        mpz_mul(N, u, v);
        int nb = (int)mpz_sizeinbase(N, 2);
        if (nb == bits || nb == bits - 1) break;
    }
}

/* Generate an RSA-STYLE semiprime: N = u*v where u and v are both balanced primes
 * of bits/2 bits each, with NO constraint on their residue mod 6. Standard RSA
 * keygen follows this pattern (modulo additional safe-prime / strong-prime
 * requirements that are orthogonal to this algorithm's behavior).
 *
 * Of the three possible (u mod 6, v mod 6) combinations for primes > 3:
 *   (5, 5)  →  N ≡ 1 (mod 6)  →  S-sieve can factor this
 *   (1, 1)  →  N ≡ 1 (mod 6)  →  NEITHER sieve handles this (both factors ≡ 1)
 *   (1, 5) or (5, 1)  →  N ≡ 5 (mod 6)  →  M-sieve can factor this
 *
 * The (1, 1) case occurs ~25% of the time with truly random primes and is
 * not solvable by any sieve in this family. This function rejects and resamples
 * those cases (reporting counts via *n_rejected_1_1 if non-NULL), so the
 * generator output is always factorable by auto-mode.
 *
 * Also enforces |u - v| > 2^(bits/4) to keep Fermat factorization from
 * trivially solving the instance (matches standard RSA key-gen guidance). */
static void gen_semiprime_rsa(mpz_t N, mpz_t u, mpz_t v,
                              int bits, gmp_randstate_t rng,
                              long *n_rejected_1_1) {
    int half = bits / 2;
    mpz_t diff, min_diff;
    mpz_inits(diff, min_diff, NULL);
    /* Minimum factor difference: 2^(bits/4). Prevents trivial Fermat case. */
    mpz_set_ui(min_diff, 1);
    mpz_mul_2exp(min_diff, min_diff, bits / 4);

    for (;;) {
        /* Random prime u of ~half bits. No mod-6 constraint. */
        mpz_urandomb(u, rng, half);
        mpz_setbit(u, half - 1);        /* full-size */
        mpz_setbit(u, half - 2);        /* high bit cluster, like OpenSSL keygen */
        mpz_setbit(u, 0);               /* odd */
        mpz_nextprime(u, u);
        if ((int)mpz_sizeinbase(u, 2) > half) continue;  /* overflow via nextprime */

        /* Random prime v of ~(bits-half) bits. No mod-6 constraint. */
        mpz_urandomb(v, rng, bits - half);
        mpz_setbit(v, bits - half - 1);
        mpz_setbit(v, bits - half - 2);
        mpz_setbit(v, 0);
        mpz_nextprime(v, v);
        if ((int)mpz_sizeinbase(v, 2) > bits - half) continue;

        if (mpz_cmp(u, v) == 0) continue;
        /* Ensure u < v for consistent output ordering. */
        if (mpz_cmp(u, v) > 0) mpz_swap(u, v);

        /* Reject too-close factors (RSA best practice, and makes the instance
         * non-trivial for our sieve too). */
        mpz_sub(diff, v, u);
        if (mpz_cmp(diff, min_diff) < 0) continue;

        /* Check N mod 6 structure. Both primes ≡ 1 (mod 6) means N ≡ 1 but
         * neither sieve handles it; resample. */
        unsigned long u6 = mpz_fdiv_ui(u, 6);
        unsigned long v6 = mpz_fdiv_ui(v, 6);
        if (u6 == 1 && v6 == 1) {
            if (n_rejected_1_1) (*n_rejected_1_1)++;
            continue;
        }
        /* (u6 == 5 && v6 == 5): N ≡ 1 (mod 6), S-sieve handles it.
         * (u6 == 1 && v6 == 5) or (5 && 1): N ≡ 5 (mod 6), M-sieve handles it. */

        mpz_mul(N, u, v);
        int nb = (int)mpz_sizeinbase(N, 2);
        if (nb == bits || nb == bits - 1) break;
    }
    mpz_clears(diff, min_diff, NULL);
}

/* ══════════════════════════════════════════════════════════════════════════
 * COMMANDS
 * ══════════════════════════════════════════════════════════════════════════ */
static void run_demo(void) {
    printf("=============================================================\n");
    printf("  K-Sieve  --  CRT Enumeration  (up to ~100-bit N)\n");
    printf("=============================================================\n\n");
    printf("  K-Sieve / S-Sieve: factors N = u*v with u == v == 5 (mod 6).\n");
    printf("  M-Sieve (new form): factors N = u*v with u == 5, v == 1 (mod 6).\n");
    printf("  Phase 1: density-sorted CRT DFS enumeration.\n");
    printf("  Phase 2: stride scan + native-u128 or GMP isqrt.\n\n");

    /* ── Old-form (K-sieve / S-sieve) test vectors.
     * u ≡ v ≡ 5 (mod 6), both prime, N ≡ 1 (mod 6). */
    static const struct { const char *lbl, *N, *u, *v; } ex_old[] = {
        {"40-bit", "689323137091",                    "655181",          "1052111"},
        {"48-bit", "210196469228359",                  "11130047",        "18885497"},
        {"56-bit", "60604137223561573",                "235829777",       "256982549"},
        {"64-bit", "10519796552841031447",             "2562877337",      "4104682031"},
        {"72-bit", "3125743374876489738127",           "50035872779",     "62470048013"},
        {"80-bit", "891800461426303612800523",         "857229368789",    "1040328871007"},
        {"90-bit", "816622855523043691765836997",      "26697544623989",  "30587938592273"},
        {"100-bit","771506244974010922570371819361",   "742710680628167", "1038770903794583"},
        {NULL, NULL, NULL, NULL}
    };

    /* ── New-form (M-sieve) test vectors.
     * u ≡ 5 (mod 6), v ≡ 1 (mod 6), both prime, N ≡ 5 (mod 6). */
    static const struct { const char *lbl, *N, *u, *v; } ex_new[] = {
        {"40-bit", "758563019861",                     "764783",          "991867"},
        {"48-bit", "264002805167909",                   "15864473",        "16641133"},
        {"56-bit", "46118103084614651",                 "206651117",       "223168903"},
        {"64-bit", "16096344416567063921",              "3962945243",      "4061712547"},
        {"72-bit", "2474954624063883618389",            "40327988759",     "61370643571"},
        {"80-bit", "1029594954498727924405661",         "941042720573",    "1094100121057"},
        {"90-bit", "945085799437462866985954391",       "32104595210813",  "29437711119907"},
        {"100-bit","896679287854581519369726704303",    "980496649757699", "914515402042597"},
        {NULL, NULL, NULL, NULL}
    };

    mpz_t N, u, v;
    mpz_inits(N, u, v, NULL);

    printf(">>> OLD FORM (K-sieve, default mode; set KSIEVE_SSIEVE=1 for S-sieve)\n\n");
    for (int i = 0; ex_old[i].lbl; i++) {
        printf("--- %s ---\n", ex_old[i].lbl);
        mpz_set_str(N, ex_old[i].N, 10);
        print_mpz_short("  N", N);
        Res r = ksieve_factor(N, u, v, 0);
        int ok = r.success && verify(N, u, v);
        if (ok) { print_mpz_short("  u", u); print_mpz_short("  v", v); }
        print_res(&r, ok);
        putchar('\n');
    }

    /* M-sieve demo: temporarily set KSIEVE_MSIEVE=1 so ksieve_factor dispatches to M-sieve.
     * Save/restore the caller's environment so --demo leaves no side effects. */
    printf(">>> NEW FORM (M-sieve, auto-enabled for this section)\n\n");
    const char *prev_m = getenv("KSIEVE_MSIEVE");
    char prev_m_buf[64] = {0};
    if (prev_m) {
        strncpy(prev_m_buf, prev_m, sizeof(prev_m_buf) - 1);
    }
    setenv("KSIEVE_MSIEVE", "1", 1);

    for (int i = 0; ex_new[i].lbl; i++) {
        printf("--- %s ---\n", ex_new[i].lbl);
        mpz_set_str(N, ex_new[i].N, 10);
        print_mpz_short("  N", N);
        Res r = ksieve_factor(N, u, v, 0);
        int ok = r.success && verify(N, u, v);
        if (ok) { print_mpz_short("  u", u); print_mpz_short("  v", v); }
        print_res(&r, ok);
        putchar('\n');
    }

    /* Restore env */
    if (prev_m) setenv("KSIEVE_MSIEVE", prev_m_buf, 1);
    else        unsetenv("KSIEVE_MSIEVE");

    mpz_clears(N, u, v, NULL);
}

static void run_info(void) {
    printf("=============================================================\n");
    printf("  K-Sieve  --  Performance Profile\n");
    printf("=============================================================\n\n");
    printf("  %-8s  %-4s  %-7s  %-12s  %-12s  %-10s  %s\n",
           "N bits","np","log2(M)","est CRT ops","est isqrt","est time","path");
    printf("  %s\n","-----------------------------------------------------------------------");

    for (int bits = 40; bits <= 100; bits += 4) {
        double K_max_d = ldexp(1.0, bits/2 - 2);
        int    np      = choose_np(bits, K_max_d, 0);

        double M_log2 = 0, M = 1.0;
        for (int j = 0; j < np; j++) { M_log2 += log2(SP[j]); M *= SP[j]; }

        double prev=1.0, total_crt=0.0, nvp_prod=1.0, M_run=1.0;
        for (int j=0; j<np; j++) {
            double nvp = SP[j]/2.0;
            total_crt += prev*nvp; M_run *= SP[j];
            double frac = (M_run>K_max_d) ? K_max_d/M_run : 1.0;
            prev = prev*nvp*frac; nvp_prod *= nvp;
        }
        double ext   = (M>=K_max_d) ? 1.0 : K_max_d/M;
        double isqrt = nvp_prod*ext*0.19;
        const char *path = (bits<=64) ? "native" : "GMP+screens";
        double ins_ns    = (bits<=64) ?   8.0 :   5.0;  /* incremental+screens */
        double crt_ns    = (bits<=64) ?  20.0 :  35.0;
        double est_s     = (total_crt*crt_ns + isqrt*ins_ns)*1e-9;

        printf("  %-8d  %-4d  %-7.1f  %-12.2e  %-12.2e  %-10.3f  %s\n",
               bits, np, M_log2, total_crt, isqrt, est_s, path);
    }
    printf("\n  Note: estimates use avg-case density (nvp ~= p/2).\n");
}

/* Per-trial result communicated from child to parent via a pipe. */
typedef struct {
    Res    r;
    int    ok;
    char   N_str[128];   /* decimal string of N for display */
    char   u_str[64];    /* decimal string of u factor */
    char   v_str[64];    /* decimal string of v factor */
} BenchTrialResult;

/* Instrumented bench: runs in-process (no fork), accumulates INSTRUMENT
 * counters across all trials, prints them at end. Build with -DINSTRUMENT
 * to get actual numbers; otherwise counters are zero. */
static void run_bench_instr(int bits, int trials) {
    printf("=============================================================\n");
    printf("  Instrumented Benchmark: %d-bit N, %d trials (in-process)\n",
           bits, trials);
    printf("=============================================================\n\n");

    gmp_randstate_t rng;
    gmp_randinit_mt(rng);
    gmp_randseed_ui(rng, 0xC0FFEE42UL);

    mpz_t N, u, v, ug, vg;
    mpz_inits(N, u, v, ug, vg, NULL);

    instr_reset();

    double sum_total = 0;
    int n_ok = 0;

    int gen_mode_instr = 0;  /* 0 = old-form, 1 = new-form, 2 = RSA-style */
    {
        const char *env_m = getenv("KSIEVE_MSIEVE");
        const char *env_r = getenv("KSIEVE_RSA");
        if (env_r && atoi(env_r)) gen_mode_instr = 2;
        else if (env_m && atoi(env_m)) gen_mode_instr = 1;
    }
    /* Note: in RSA mode, ksieve_factor auto-selects sieve per N based on the
     * KSIEVE_RSA env var; no additional setup is needed here. */
    long n_rsa_rejected_instr = 0;

    for (int t = 0; t < trials; t++) {
        if (gen_mode_instr == 2)      gen_semiprime_rsa(N, ug, vg, bits, rng, &n_rsa_rejected_instr);
        else if (gen_mode_instr == 1) gen_semiprime_newform(N, ug, vg, bits, rng);
        else                          gen_semiprime(N, ug, vg, bits, rng);
        Res r = ksieve_factor(N, u, v, 0);
        int ok = r.success && verify(N, u, v);
        if (ok) n_ok++;
        sum_total += r.t_total;
        printf("  %3d: total=%7.4fs  crt_ops=%9ld  isqrt=%6ld  %s\n",
               t+1, r.t_total, r.n_crt_ops, r.n_isqrt, ok ? "OK" : "FAIL");
    }

    printf("\n--- Summary ---\n");
    printf("  OK: %d/%d  mean_total=%.4fs\n", n_ok, trials, sum_total/trials);
    printf("\n");
    instr_print(stdout);

    mpz_clears(N, u, v, ug, vg, NULL);
    gmp_randclear(rng);
}

static void run_bench(int bits, int trials) {
    /* Timeout per trial: scale with bits.  A healthy 80-bit run takes ~2s;
     * give 5× headroom so only genuinely hung trials are killed.           */
    int timeout_s = (bits <= 72) ? 10 : (bits <= 80) ? 30 : (bits <= 100) ? 120
                  : (bits <= 110) ? 300 : (bits <= 120) ? 3600 : 20000;

    printf("=============================================================\n");
    printf("  Benchmark: %d-bit N, %d trials  (timeout %ds/trial)\n",
           bits, trials, timeout_s);
    printf("=============================================================\n\n");

    /* Generate all N values upfront so the RNG sequence is deterministic
     * regardless of fork/pipe timing.                                      */
    gmp_randstate_t rng;
    gmp_randinit_mt(rng);
    gmp_randseed_ui(rng, 0xC0FFEE42UL);

    mpz_t *Ns = malloc(trials * sizeof(mpz_t));
    mpz_t ug, vg;
    mpz_inits(ug, vg, NULL);
    /* Mode selection:
     *   KSIEVE_RSA=1    → RSA-style (balanced primes, no mod-6 constraint).
     *                      Auto-select sieve per-N based on N mod 6.
     *   KSIEVE_MSIEVE=1 → new-form generator (N ≡ 5 mod 6).
     *   otherwise       → old-form generator (N ≡ 1 mod 6). */
    int gen_mode = 0;   /* 0 = old-form, 1 = new-form, 2 = RSA-style */
    {
        const char *env_m = getenv("KSIEVE_MSIEVE");
        const char *env_r = getenv("KSIEVE_RSA");
        if (env_r && atoi(env_r)) gen_mode = 2;
        else if (env_m && atoi(env_m)) gen_mode = 1;
    }
    long n_rsa_rejected_1_1 = 0;
    for (int t = 0; t < trials; t++) {
        mpz_init(Ns[t]);
        if (gen_mode == 2)      gen_semiprime_rsa(Ns[t], ug, vg, bits, rng, &n_rsa_rejected_1_1);
        else if (gen_mode == 1) gen_semiprime_newform(Ns[t], ug, vg, bits, rng);
        else                    gen_semiprime(Ns[t], ug, vg, bits, rng);
    }
    mpz_clears(ug, vg, NULL);
    gmp_randclear(rng);

    if (gen_mode == 2) {
        /* RSA mode: sieve mode is auto-selected per-N inside ksieve_factor() based
         * on the KSIEVE_RSA=1 flag. Warn once (before forking child processes)
         * if the user also set an explicit sieve flag — RSA takes precedence. */
        const char *env_m_warn = getenv("KSIEVE_MSIEVE");
        const char *env_s_warn = getenv("KSIEVE_SSIEVE");
        int has_m = (env_m_warn && atoi(env_m_warn));
        int has_s = (env_s_warn && atoi(env_s_warn));
        if (has_m || has_s) {
            fprintf(stderr,
                "WARNING: KSIEVE_RSA=1 is set together with KSIEVE_%sSIEVE=1; "
                "RSA mode overrides the explicit sieve choice and auto-selects "
                "per N based on N mod 6.\n\n",
                has_m ? "M" : "S");
        }
        long total_draws = (long)trials + n_rsa_rejected_1_1;
        printf("[RSA mode] generated %d usable semiprimes; rejected %ld draws "
               "where both primes ≡ 1 mod 6 (%.1f%% of %ld total draws; "
               "theoretical expectation ~25%% for uniform random primes)\n\n",
               trials, n_rsa_rejected_1_1,
               total_draws > 0 ? (100.0 * n_rsa_rejected_1_1 / total_draws) : 0.0,
               total_draws);
    }

    printf("%-4s  %-6s  %-4s  %-10s  %-10s  %-9s  %-10s  %-10s  %s\n",
           "#","N-bit","np","crt_ops","k0s","isqrt","p1(s)","total(s)","OK?");
    printf("%-4s  %-6s  %-4s  %-10s  %-10s  %-9s  %-10s  %-10s  %s\n",
           "----","------","----","----------","----------","---------",
           "----------","----------","---");

    double s_tot=0, s_p1=0;
    long   s_isq=0, s_crt=0, s_k0=0;
    int    n_ok=0, n_fail=0, n_timeout=0;
    double t_min=1e18, t_max=0;

    for (int t = 0; t < trials; t++) {
        /* Get N string for display before forking */
        char N_str[128];
        gmp_snprintf(N_str, sizeof(N_str), "%Zd", Ns[t]);

        /* Pipe: child writes BenchTrialResult, parent reads it */
        int pfd[2];
        if (pipe(pfd) != 0) { perror("pipe"); exit(1); }

        pid_t child = fork();
        if (child < 0) { perror("fork"); exit(1); }

        if (child == 0) {
            /* ── CHILD ── */
            close(pfd[0]);
            alarm(timeout_s);   /* SIGALRM kills child if it hangs */

            mpz_t uf, vf;
            mpz_inits(uf, vf, NULL);

            BenchTrialResult tr;
            tr.r  = ksieve_factor(Ns[t], uf, vf, 0);
            tr.ok = tr.r.success && verify(Ns[t], uf, vf);
            strncpy(tr.N_str, N_str, sizeof(tr.N_str)-1);
            tr.N_str[sizeof(tr.N_str)-1] = '\0';
            if (tr.ok) {
                gmp_snprintf(tr.u_str, sizeof(tr.u_str), "%Zd", uf);
                gmp_snprintf(tr.v_str, sizeof(tr.v_str), "%Zd", vf);
            } else {
                tr.u_str[0] = '\0';
                tr.v_str[0] = '\0';
            }

            /* Write result to pipe (best-effort; if alarm fires we never get here) */
            write(pfd[1], &tr, sizeof(tr));
            close(pfd[1]);
            mpz_clears(uf, vf, NULL);
            _exit(0);
        }

        /* ── PARENT ── */
        close(pfd[1]);

        /* Wait for child with our own wall-clock fallback */
        BenchTrialResult tr;
        memset(&tr, 0, sizeof(tr));
        int timed_out = 0;

        ssize_t got = 0, total = 0;
        char *buf = (char*)&tr;
        /* Blocking read loop — child may write in chunks.
         * Returns 0 (EOF) when child closes pipe (either via _exit or SIGALRM). */
        while (total < (ssize_t)sizeof(tr)) {
            got = read(pfd[0], buf + total, sizeof(tr) - total);
            if (got <= 0) break;
            total += got;
        }
        close(pfd[0]);

        /* Reap child (blocking).  The read above already waited for the child
         * to finish or die, so this returns almost immediately.               */
        int status;
        waitpid(child, &status, 0);

        /* Determine outcome: if we got a short read the child was killed
         * (SIGALRM or other signal) before it could write the result.         */
        if (total < (ssize_t)sizeof(tr) || WIFSIGNALED(status)) {
            timed_out = 1;
        }

        if (timed_out) {
            printf("%-4d  %-6d  %-4s  %-10s  %-10s  %-9s  %-10s  %-10.1f  %s\n",
                   t+1, bits, "?", "?", "?", "?", "?",
                   (double)timeout_s, "TIMEOUT");
            printf("      N=%s\n", N_str);
            fflush(stdout);
            n_timeout++;
        } else {
            int ok = tr.ok;
            printf("%-4d  %-6d  %-4d  %-10ld  %-10ld  %-9ld  %-10.4f  %-10.4f  %s\n",
                   t+1, tr.r.N_bits, tr.r.np, tr.r.n_crt_ops, tr.r.n_k0s, tr.r.n_isqrt,
                   tr.r.t_p1, tr.r.t_total,
                   ok ? "OK" : (tr.r.success ? "WRONG" : "FAIL"));
            printf("      N=%s", tr.N_str);
            if (ok && tr.u_str[0])
                printf("  u=%s  v=%s", tr.u_str, tr.v_str);
            printf("\n");
            fflush(stdout);

            if (ok) {
                s_tot += tr.r.t_total; s_p1 += tr.r.t_p1;
                s_isq += tr.r.n_isqrt; s_crt += tr.r.n_crt_ops; s_k0 += tr.r.n_k0s;
                if (tr.r.t_total < t_min) t_min = tr.r.t_total;
                if (tr.r.t_total > t_max) t_max = tr.r.t_total;
                n_ok++;
            } else { n_fail++; }
        }
    }

    for (int t = 0; t < trials; t++) mpz_clear(Ns[t]);
    free(Ns);

    printf("\n--- Summary (%d OK / %d FAIL / %d TIMEOUT / %d total) ---\n",
           n_ok, n_fail, n_timeout, trials);
    if (n_ok > 0) {
        printf("  Mean p1    : %.4f s\n", s_p1  / n_ok);
        printf("  Mean total : %.4f s\n", s_tot / n_ok);
        printf("  Min / Max  : %.4f / %.4f s\n", t_min, t_max);
        printf("  Mean CRT   : %.0f\n", (double)s_crt / n_ok);
        printf("  Mean K0s   : %.0f\n", (double)s_k0 / n_ok);
        printf("  Mean isqrt : %.0f\n", (double)s_isq  / n_ok);
    }
}

static void run_factor(const char *s) {
    mpz_t N, u, v;
    mpz_inits(N, u, v, NULL);
    if (mpz_set_str(N, s, 10) != 0) {
        fprintf(stderr, "Error: '%s' is not a valid decimal integer\n", s);
        exit(1);
    }
    if (mpz_cmp_ui(N, 1) <= 0 || mpz_probab_prime_p(N, 20)) {
        fprintf(stderr, "Error: N must be composite and > 1\n");
        exit(1);
    }
    {
        /* Warn if N mod 6 doesn't match the selected mode.
         * K/S-sieve require N ≡ 1 (mod 6); M-sieve requires N ≡ 5 (mod 6). */
        unsigned long nm6 = mpz_fdiv_ui(N, 6);
        const char *env_m = getenv("KSIEVE_MSIEVE");
        const char *env_s = getenv("KSIEVE_SSIEVE");
        const char *env_a = getenv("KSIEVE_AUTO");
        int explicit_m = (env_m && atoi(env_m));
        int explicit_s = (env_s && atoi(env_s));
        int auto_mode  = (env_a && atoi(env_a));
        int expect_new = explicit_m || (auto_mode && nm6 == 5);
        int expect_old = explicit_s || (!explicit_m && !auto_mode);
        if (expect_new && nm6 != 5) {
            fprintf(stderr,
                "Warning: N == %lu (mod 6); M-sieve requires N == 5 (mod 6)\n"
                "  (factors must be primes of form (6k-1) and (6k+1))\n", nm6);
        } else if (expect_old && nm6 != 1) {
            fprintf(stderr,
                "Warning: N == %lu (mod 6); K/S-sieve requires N == 1 (mod 6)\n"
                "  (both factors must be primes of the form 6k-1)\n"
                "  Set KSIEVE_MSIEVE=1 (or KSIEVE_AUTO=1) for N == 5 (mod 6).\n", nm6);
        }
    }

    int bits = (int)mpz_sizeinbase(N, 2);
    if (bits > 104)
        fprintf(stderr,
            "Warning: N is %d-bit; optimised for <= 100 bits. May be slow.\n",
            bits);

    printf("Factoring %d-bit N...\n\n", bits);
    Res r = ksieve_factor(N, u, v, 1);
    printf("\n");

    int ok = r.success && verify(N, u, v);
    if (ok) { printf("Factored!\n"); print_mpz_short("  u", u); print_mpz_short("  v", v); }
    else if (r.success) printf("Internal error: result fails verification.\n");
    else printf("Not found. Verify both factors are == 5 (mod 6).\n");
    print_res(&r, ok);

    mpz_clears(N, u, v, NULL);
}

/* ── Main ─────────────────────────────────────────────────────────────── */
int main(int argc, char *argv[]) {
    init_qr_tables();
    if (argc < 2 || !strcmp(argv[1], "--demo")) {
        run_demo();
    } else if (!strcmp(argv[1], "--bench")) {
        if (argc < 4) { fprintf(stderr,"Usage: --bench <bits> <n>\n"); return 1; }
        int bits = atoi(argv[2]), trials = atoi(argv[3]);
        if (bits < 10 || bits > 160) { fprintf(stderr,"bits: 10-160\n"); return 1; }
        if (trials < 1 || trials > 10000) { fprintf(stderr,"trials: 1-10000\n"); return 1; }
        run_bench(bits, trials);
    } else if (!strcmp(argv[1], "--bench-instr")) {
        if (argc < 4) { fprintf(stderr,"Usage: --bench-instr <bits> <n>\n"); return 1; }
        int bits = atoi(argv[2]), trials = atoi(argv[3]);
        if (bits < 10 || bits > 160) { fprintf(stderr,"bits: 10-160\n"); return 1; }
        if (trials < 1 || trials > 10000) { fprintf(stderr,"trials: 1-10000\n"); return 1; }
        run_bench_instr(bits, trials);
    } else if (!strcmp(argv[1], "--info")) {
        run_info();
    } else if (!strcmp(argv[1], "--help") || !strcmp(argv[1], "-h")) {
        printf("Usage:\n");
        printf("  %s <N>               Factor decimal N (u,v must be == 5 mod 6)\n", argv[0]);
        printf("  %s --demo            Examples from 40-100 bit\n", argv[0]);
        printf("  %s --bench <b> <n>   Benchmark n random b-bit semiprimes\n", argv[0]);
        printf("  %s --info            Performance profile table\n", argv[0]);
    } else {
        run_factor(argv[1]);
    }
    return 0;
}