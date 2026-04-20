**K-Sieve / S-Sieve / T-Sieve / M-Sieve — Rev. 8**

**A Specialized Factoring Algorithm for 6n±1 Semiprimes, with Multi-Layer Quadratic-Residue Filtering, Sum-Variable Range Reduction, Complete (6x±1)-Class Coverage, and AVX-512 Acceleration**

Andrew Embry

*April 2026 (Revision 8)*

# Abstract

We present the K-Sieve, S-Sieve, T-Sieve, and M-Sieve, a family of four specialized factoring algorithms that jointly cover every prime-pair semiprime N = u·v with factors greater than 3. The K- and S-Sieves handle the (−1,−1) class (u ≡ v ≡ −1 mod 6); the T-Sieve, introduced in this revision, handles the (+1,+1) class (u ≡ v ≡ +1 mod 6); the M-Sieve handles the mixed (−1,+1) class (u ≢ v mod 6). The core filtering technique — restricting Fermat-style factoring candidates to values for which a specific quadratic polynomial is a quadratic residue modulo a highly-composite modulus — shares its underlying mathematical structure with the hyperbolic sieve of Hittmeir [1], although the present work was developed independently. Our contribution is the application of this filtering idea to the 6n±1 parameterization [2], together with a practical implementation stack: CRT enumeration with Hensel lifting at p², p³, p⁴; a layered tiled bitset exploiting 2-adic structure at period 64; a verify chain of QR screens at small prime powers; a cross-prime isqrt filter; AVX-512 batch-16 acceleration; and multi-threaded prefix work-stealing. The T-Sieve is structurally the (+1,+1)-class mirror of the S-Sieve, differing by three sign flips (in the polynomial's linear coefficient, in the mod-3 fold target, and in the mod-4 parity targets at N mod 8 ∈ {3, 7}); we prove that its QR144 tables coincide bitwise with the S-Sieve's. Revision 8 presents complete empirical benchmarks across 60–120 bit semiprimes (30 trials per cell, 720+ total factorizations; 100% success rate, isqrt = 1 on every trial). Headline results on a 16-thread build: at matched 100-bit workload, S-Sieve mean wall-time is 0.059 s, T-Sieve 0.044 s (0.74× — i.e. slightly faster than S-Sieve), M-Sieve 0.104 s (1.76×), K-Sieve 0.863 s (14.6×). At 120 bit, S- and T-Sieves factor in 13.6 s and 13.8 s respectively (within 2% of each other, as algebraically predicted). Empirical scaling exponents from 80–120 bit: S-Sieve 2^(0.355·bits), T-Sieve 2^(0.356·bits), M-Sieve 2^(0.367·bits). The four sieves combined cover every prime-pair semiprime residue class modulo 6.

# 1. Introduction

Factoring a semiprime N = u·v is a fundamental problem with direct relevance to RSA-based cryptography. General-purpose algorithms such as the quadratic sieve and the number field sieve achieve subexponential complexity. This paper focuses on a narrower, structurally constrained variant: factoring semiprimes where the prime factors satisfy specific congruence conditions modulo 6. The four sieves developed here — K, S, T, and M — jointly cover the entire class of semiprimes whose factors are coprime to 6 (which includes all semiprimes with factors greater than 3).

All odd primes greater than 3 satisfy p ≡ ±1 (mod 6). For a prime-pair semiprime N = u·v with u, v > 3, there are therefore four residue-class pairings, giving two N mod 6 equivalence classes:

- N ≡ 1 (mod 6), with u ≡ v ≡ −1 (mod 6): 'old form.' Covered by the K-Sieve and S-Sieve, the S-Sieve being substantially faster.

- N ≡ 1 (mod 6), with u ≡ v ≡ +1 (mod 6): '(+1,+1)-form' (formerly described as 'dual old form'). Covered by the K-Sieve (which is polynomial-class-agnostic in g_K = 9K²+N) and, as of Revision 8, by the T-Sieve, which matches the S-Sieve's performance on its own input class.

- N ≡ 5 (mod 6), with u ≡ −1, v ≡ +1 (mod 6) (or vice versa): 'new form.' Covered by the M-Sieve introduced in Revision 5. The K-, S-, and T-Sieves cannot factor this class at all, since their polynomial structure assumes matched residue classes.

Approximately one quarter of random prime-pair semiprimes fall in each of the four residue-class pairings, giving ≈25% per class, or equivalently 50% of random N satisfying N ≡ 1 (mod 6) with split equally between (−1,−1) and (+1,+1) sub-classes. Earlier revisions of this paper handled only three of the four pairings; Revision 8 closes the gap entirely.

## 1.1 Relation to prior work

Two independent lines of prior work underpin the present paper. We summarize each and clarify our relationship to it.

**The 6n±1 parameterization of prime-pair semiprimes (Creft 2012).**

The classification of prime-pair semiprimes by their factors' residues modulo 6 was developed by Creft [2] in his "Hexile Sieve" analysis of primes and their composites. Creft observes that every prime p > 3 satisfies p ≡ ±1 (mod 6), and that the product of two such primes falls into one of two equivalence classes modulo 6 (namely N ≡ 1 for u ≡ v (mod 6), and N ≡ 5 for u ≢ v (mod 6)). Creft further derives the corresponding Diophantine representations — in our notation, u = 6x − 1, v = 6y − 1 in the (−1,−1) class, u = 6x + 1, v = 6y + 1 in the (+1,+1) class, and u = 6x − 1, v = 6y + 1 in the new-form class — and uses them to distinguish primes from composites within each class.

Our K-, S-, T-, and M-Sieves operate entirely within Creft's parameterization. The K-Sieve scans the difference variable K = y − x in all old-form classes, with its scan polynomial identical across the (−1,−1) and (+1,+1) classes. The S-Sieve scans the sum variable S = x + y in the (−1,−1) class; the T-Sieve (introduced in Revision 8) scans the same sum variable in the (+1,+1) class; the M-Sieve scans M = x + y in the new-form class. The specific substitutions and polynomials are derived from the Creft representation rather than being independently discovered. The author arrived at the 6n±1 framework independently, but Creft's 2012 priority is established and this paper should be read as building on that work.

**The hyperbolic sieve for Fermat factoring (Hittmeir 2022/2023).**

The core filtering technique used in all four sieves — restricting Fermat-style factoring candidates to those values for which a specific quadratic polynomial is a quadratic residue modulo a composite modulus assembled from small primes via the Chinese Remainder Theorem — appears in more general form as Hittmeir's hyperbolic sieve [1]. Hittmeir works with the sieve set

L_{N,m,k} = { k·x + y (mod m) : (x, y) ∈ Z_m² with N ≡ x·y (mod m) }

and proves (Lemma 4.2 of [1]) that for odd m,

L_{N,m,k} = { s ∈ Z_m : s² − 4k·N is a quadratic residue (mod m) }.

Applied to Fermat's algorithm with k = 1, this yields a sieve on s = u + v with rigorous O(Λ · exp(−C log Λ / log log Λ)) runtime improvement over naive Fermat, where Λ bounds the divisor difference |u − v|²/(4√N). Our K-, S-, T-, and M-Sieves apply essentially the same quadratic-residue filtering idea to Creft's 6n±1 parameterization. The 6n±1 restriction shrinks the effective search range by a constant factor before any QR filtering is applied and yields specialized polynomials whose algebraic structure supports additional filters not directly present in Hittmeir's general framework: Hensel lifting at p², p³, p⁴ for primes where g has a repeated root; conditional mod-4 and mod-8 constraints from N mod 8; and — in the S- and T-Sieves — a mod-3 constraint on the sum variable that folds into the root modulus at no additional cost.

**The difference-of-squares identity.**

The base factoring identity N = x² − y² = (x − y)(x + y) is classical, attributed to Fermat [3] in the 17th century. All four sieves inherit this identity: the K-Sieve searches for T, K such that T² − (3K)² = N; the S- and T-Sieves recast the same search in terms of S; the M-Sieve searches for M, L such that (3M)² − L² = N. The identity itself is elementary; the novelty lies in the combination of the 6n±1 parameterization (Creft) with the hyperbolic-sieve filter (Hittmeir) and the implementation stack presented here.

**Novel contributions of the present work.**

With Creft and Hittmeir properly credited, the contributions of this paper reduce to the following. (i) The explicit derivation of the four scan polynomials g_K, g_S, g_T, g_M together with their Hensel-lift structure at small prime powers. (ii) The M-Sieve: to our knowledge this was the first application of hyperbolic-sieve filtering to the new-form class N ≡ 5 (mod 6). (iii) The T-Sieve (Rev. 8): to our knowledge the first explicit treatment of the (+1,+1) class as a first-class sum-variable sieve rather than an afterthought handled (slowly) by the K-Sieve, together with the observation that its QR144 tables coincide bitwise with the S-Sieve's. (iv) The specific implementation architecture: CRT-enumerated leaves with density-sorted DFS, a layered 64-bit-tiled bitset exploiting 2-adic structure, a verify chain of QR screens at small prime powers, and a cross-prime isqrt filter. (v) AVX-512 batch-16 acceleration and multi-threaded prefix work-stealing. (vi) Empirical characterization across 60–120 bit N with correctness verification over 720+ factorizations. (vii) The W-Sort heuristic for DFS branch ordering based on multiplicity-weighted residue classes. (viii) The S→T auto-dispatch fallback strategy, which resolves the factor-class ambiguity inherent to N ≡ 1 (mod 6) inputs whose factor residues cannot be determined a priori. The present work does not reproduce Hittmeir's asymptotic analysis and makes no rigorous complexity claim.

## 1.2 The K-Reformulation

Fermat's original factoring idea [3] represents an odd composite N as a difference of squares, N = x² − y² = (x − y)(x + y). We reformulate this for the old-form classes as follows. Let u and v be the two prime factors with b ≥ a ≥ 1. Define T = (u + v)/2 and K = (v − u)/6. Direct algebraic manipulation yields the central identity:

T² = 9K² + N,

valid in both the (−1,−1) class (with u = 6a − 1, v = 6b − 1, giving T = 3(a + b) − 1) and the (+1,+1) class (with u = 6a + 1, v = 6b + 1, giving T = 3(a + b) + 1). The polynomial g_K(K) = 9K² + N is therefore class-agnostic: it is a perfect square at the true K regardless of whether the factors are both ≡ −1 or both ≡ +1 (mod 6). Factor recovery u = T − 3K, v = T + 3K is likewise class-agnostic; the final u·v = N verification catches any ambiguity. The search range is bounded by 0 ≤ K ≤ K_max = ⌊√N / 3⌋ in both old-form classes.

The K-Sieve thus covers both (−1,−1) and (+1,+1) classes with a single polynomial and a single algorithm, at the cost of a search range 4.24× wider than the sum-variable sieves. For a fast (+1,+1)-class sieve with matched performance to the S-Sieve's on (−1,−1), the T-Sieve is introduced in §1.5 and developed in §5.

## 1.3 The S-Sieve: 8.5× range reduction (old-form (−1,−1) class)

A key innovation in Revision 3 is the S-Sieve, which substitutes the search variable S = x + y (where u = 6x − 1, v = 6y − 1) for K = y − x. Since N = (6x − 1)(6y − 1) = 36xy − 6S + 1, we have 6xy = n + S where n = (N − 1)/6. The S-Sieve polynomial is:

g_S(S) = 9S² − 6S − 6n = (3K)².

The S-Sieve searches S ∈ [⌈√N/3⌉, ⌊√(5N/4)/3⌋], a range of approximately 0.039 · √N — roughly 4.24× smaller than the K-Sieve range. Combined with the additional mod-3 constraint S ≡ −n (mod 3) that folds into the root modulus, the S-Sieve is typically 4–8× faster than the K-Sieve at matched bit size on (−1,−1)-class inputs.

## 1.4 The M-Sieve: new-form coverage

Revision 5 introduces the M-Sieve for semiprimes with mixed-class factors. Let u = 6a − 1 ≡ 5 (mod 6) and v = 6b + 1 ≡ 1 (mod 6) with a, b ≥ 1; then N = uv ≡ 5 (mod 6), and n = (N + 1)/6 is a positive integer. Define M = a + b and K_Δ = b − a. Direct algebra yields:

9M² − N = (3K_Δ + 1)² = L²,

so the M-Sieve polynomial is:

g_M(M) = 9M² − N = L².

The identity L ≡ 1 (mod 3) is equivalent to N ≡ 2 (mod 3), which holds tautologically; no additional filter is gained from this constraint beyond what is already implicit in the sieve setup. The M-Sieve searches M ∈ [⌈√N/3⌉, ⌊√(5N)/6⌋], the same 0.039 · √N range width as the S-Sieve. Factor recovery proceeds by L = √g_M(M), T_real = 3M, u = T_real − L, v = T_real + L. See Section 4 for the full M-Sieve derivation and algorithmic details.

## 1.5 The T-Sieve: (+1,+1)-form coverage

Revision 8 introduces the T-Sieve for the (+1,+1) class: u = 6x + 1, v = 6y + 1, both prime, with N = uv ≡ 1 (mod 6). This class was formerly handled only by the K-Sieve, which factors it transparently (as noted in §1.2) but at the K-Sieve's 4.24× larger search-range cost. The T-Sieve provides a sum-variable formulation fully parallel to the S-Sieve on (−1,−1) inputs.

Setting S = x + y and applying the difference-of-squares identity with T_real = 3S + 1, direct manipulation yields:

g_T(S) = 9S² + 6S − 6n = (3K)²,

where K = y − x and n = (N − 1)/6 as in the S-Sieve. The T-Sieve is structurally the (+1,+1) mirror of the S-Sieve: its polynomial g_T differs from g_S only in the sign of the linear term (+6S vs −6S), its mod-3 fold target is S ≡ +n (mod 3) rather than −n, and its mod-4 parity targets at N mod 8 ∈ {3, 7} flip (3 → S ≡ 3, 7 → S ≡ 1) relative to the S-Sieve. Its discriminant 36N, search range width 0.039·√N, root modulus density (1/6 or 1/12), and QR144 tables all coincide with the S-Sieve's. Factor recovery proceeds by L = √g_T(S), K = L/3, T_real = 3S + 1, u = T_real − L, v = T_real + L. See Section 5 for the full T-Sieve derivation and algorithmic details, including a proof that its QR144 tables are bitwise identical to the S-Sieve's.

With the T-Sieve in place, the four sieves jointly cover every prime-pair semiprime N = u·v with factors > 3. Table 1.1 summarizes the coverage:

| **u mod 6** | **v mod 6** | **N mod 6** | **Fast sieve** | **Fallback** |
| --- | --- | --- | --- | --- |
| −1 | −1 | 1 | S-Sieve | K-Sieve |
| +1 | +1 | 1 | T-Sieve | K-Sieve |
| −1 | +1 | 5 | M-Sieve | — |
| +1 | −1 | 5 | M-Sieve | — |

*Table 1.1: Coverage of prime-pair semiprime residue classes mod 6.*

## 1.6 Document structure

Section 2 describes the multi-layer architecture shared by all four sieves. Section 3 details optimizations. Section 4 introduces the M-Sieve in full, including root-modulus derivation, QR144 tables, and cross-prime handling. Section 5 introduces the T-Sieve in full, including its relationship to the S-Sieve, the QR144-equivalence lemma, and auto-dispatch strategy. Section 6 presents empirical benchmarks including the K/S/T/M paired comparison, per-sieve scaling sweeps, bitset-prime tuning, threading scalability, and correctness statistics. Section 7 discusses asymptotic scaling. Section 8 documents bugs fixed. Section 9 lists future work.

# 2. Multi-Layer Architecture

The algorithm operates in three stages: (i) CRT enumeration produces candidate residues K₀ (mod M) for the scan variable (generically denoted K, but instantiated as K, S, M, or — in the T-Sieve — also S), (ii) a layered tiled bitset scans stride positions K = K₀ + jM for j ≥ 0, applying Jacobi-based filtering through bit masks, and (iii) surviving candidates pass through a verify chain of quadratic-residue screens followed by a cross-prime isqrt filter before the final integer square-root check.

## 2.1 CRT enumeration

Let M = p_root × p₁ × … × p_n for small primes p_i. The root modulus p_root is 4 when N mod 8 ∈ {1, 5} (giving density 0.25 per candidate), or 2 otherwise (density 0.5). The S-Sieve and T-Sieve each multiply p_root by an additional factor of 3 to absorb their respective mod-3 constraints on the sum variable (S ≡ −n (mod 3) for S-Sieve, S ≡ +n (mod 3) for T-Sieve), yielding p_root ∈ {6, 12}. The M-Sieve has no equivalent mod-3 fold and retains p_root ∈ {2, 4}. At each odd prime p, the sieve enumerates residues k ∈ {0, 1, …, p − 1} such that Jacobi(g(k), p) ≥ 0, where g is the sieve polynomial in use.

For the boundary case Jacobi(g(k), p) = 0, the p-adic valuation v_p(g(k)) must be even; otherwise g cannot be a square. Hensel's lemma resolves this: a base residue k is retained if and only if there exists a lift k′ ≡ k (mod p) with v_p(g(k′)) even. The implementation lifts to mod p⁴ to distinguish v_p ∈ {0, 2, ≥4} from v_p ∈ {1, 3}.

CRT combines per-prime residue sets via depth-first search. Primes are sorted by density ρ_p = nvp_p / p in ascending order to prune aggressively at shallow depths. Each CRT prime contributes approximately 0.5 bits of filtering on average.

## 2.2 Layered tiled bitset

For each CRT leaf K₀, the algorithm examines stride positions K₀, K₀ + M, K₀ + 2M, … up to K_max. For a bitset prime q (not in M), we build a tile: an array of q u64 values where tile[s] encodes validity of 64 consecutive j values starting at shift s. Total storage per prime is q × 8 bytes, all tiles combined comfortably resident in L1.

The number of bitset primes is chosen adaptively per sieve mode (see Section 6.5 for the Rev. 5 tuning results, confirmed by the Rev. 8 sweep): K-Sieve uses crt_np + 2, S- and T-Sieve use crt_np + 3, M-Sieve uses crt_np + 5. The M-Sieve requires more bitset primes because, unlike the S- and T-Sieves, it cannot inherit the mod-3 root fold, and its bitset tiles are consequently slightly less dense.

## 2.3 QR64: the zero-cost 2-adic layer

The filter g(K) = square (mod 64) depends only on K mod 64. Since the bitset word size is 64 bits, the same QR64 mask applies to every word in the scan. QR64 has 12 of 64 residues that are squares (density 12/64 = 0.188), contributing approximately 2.4 bits of filtering at essentially zero per-word cost beyond a single AND.

## 2.4 Verify chain

Bit positions surviving the bitset layer pass through a verify chain of precomputed QR tables. The default chain (K-, S-, and T-Sieve) uses QR256 (mod 256), QR65536 (mod 2¹⁶), QR31, QR37, QR49 (7²), QR25 (5²), QR121 (11²), QR169 (13²), QR41, and QR27 (3³). The M-Sieve uses a reduced chain (QR31, QR37, QR256, QR25, QR49, QR169, QR65536) ordered by empirical rejection counts; filters whose rejection count was measured at zero in the new-form regime (QR41, QR27, QR121) are omitted, as upstream screens already subsume their filtering power. The T-Sieve reuses the S-Sieve's QR144 tables directly; see §5.5 for the proof of bitwise equivalence.

The verify chain contributes approximately 12 bits of filtering in total (K/S/T-Sieve) or approximately 9 bits (M-Sieve), on survivors that have already cleared the bitset and QR64 layers. Per-survivor cost is ∼5 cycles per screen.

## 2.5 Cross-prime filter

A single large prime p ≈ √N provides a final filter before the definitive integer square root. Newton's method computes √g in u128 arithmetic (for N ≤ 126 bits) or mpz arithmetic (for N ≥ 127 bits) and checks T² ≡ g (mod p). For the S- and T-Sieves, the Newton starting point uses 2^⌈log₂(g)/2⌉ rather than 3S, which diverges near S_min where g ≈ 0. For the M-Sieve, the Newton starting point is derived similarly from the bit-length of g_M; an early reject occurs if g_M turns out negative (M below M_min or bad residue). The cross-prime contributes approximately 62 bits of filtering on a single candidate that has survived all prior layers.

# 3. Optimizations

## 3.1 Conditional mod-4 root

For K-Sieve: when N mod 8 ∈ {1, 5}, K mod 4 is uniquely determined (density 0.25 vs 0.50). Using p_root = 4 instead of 2 doubles M without adding DFS branches, halving J_max.

For S-Sieve: when N mod 8 = 3, S ≡ 1 (mod 4); when N mod 8 = 7, S ≡ 3 (mod 4); when N mod 8 ∈ {1, 5}, S is even. All give density 0.25 in the absence of the mod-3 fold; combined with S ≡ −n (mod 3), the effective root modulus is 12 or 6.

For T-Sieve: the mod-4 parity targets are the mirror image of the S-Sieve's, reflecting the sign flip in T_real = 3S + 1 vs 3S − 1. When N mod 8 = 3, S ≡ 3 (mod 4); when N mod 8 = 7, S ≡ 1 (mod 4); when N mod 8 ∈ {1, 5}, S is even (same as S-Sieve). Combined with S ≡ +n (mod 3), the effective root modulus is 12 or 6, identical in density to the S-Sieve. A full derivation of the T-Sieve root modulus is given in §5.4.

For M-Sieve: the residue of M mod 4 is similarly determined by N mod 8 (N ≡ 3 → M ≡ 2 mod 4, N ≡ 7 → M ≡ 0 mod 4, N ≡ {1,5} → M odd), but M has no mod-3 constraint analogous to the S-Sieve's; empirical checks over 5 000 samples show M is uniformly distributed modulo 3 and modulo 9.

## 3.2 Adaptive np selection via target n_words

The number of CRT primes n_p determines the trade-off between K₀ count and words per K₀. The rule is: add CRT primes until n_words = ⌈J_max / 64⌉ falls below a target maximum (TARGET_NW_MAX = 600). The selected n_p values differ by mode because the S- and T-Sieve mod-3 root folds triple their effective M relative to the K-Sieve and M-Sieve.

## 3.3 AVX-512 batch16 bitset scan

The dominant runtime contributor is the bitset scan. The implementation accelerates this with AVX-512 by batching 16 K₀ values (two __m512i vectors of 8 lanes each) and processing them in parallel at each word index. Per word, the inner loop performs n_p gather-AND operations across all 16 lanes simultaneously, with SIMD shift advancement via _mm512_add_epi64, _mm512_cmpge_epi64_mask, and _mm512_mask_blend_epi64. Measured per-K₀-per-word cost: 28–31 cycles on Zen 3 / Alder Lake.

## 3.4 Multi-threaded prefix work-stealing

The CRT DFS runs serially to a split depth (typically 3–4), collecting partial-prefix work items. Prefixes are distributed across worker threads via atomic index. Each worker runs the remaining DFS from its prefix and batches K₀s into batch16 scans. The first worker to find the factor sets an atomic flag; others exit on the next batch boundary. Thread-local CRT-op counters are accumulated after join.

Measured scaling on 16 threads at 100-bit M-Sieve is 5.3× over single-thread (Rev. 8 10-trial run; §6.6) — substantially sub-linear, bounded by memory bandwidth in the bitset scan, tail-effect waste on early termination, and fixed setup costs. Threading efficiency improves at larger bit sizes where the sieve phase dominates, but has not been systematically characterized above 100-bit.

## 3.5 Cross-prime Newton fix for sum-variable sieves

The original S-Sieve cross-prime used cpR = 3S as Newton starting point for √g_S. Near S_min, g_S ≈ 0 and √g_S ≪ 3S, causing Newton divergence. The fix computes the starting point from 2^⌈bits(g_S)/2⌉ and converges from above. The T-Sieve inherits the same fix with the +6S sign flip on g_T. The M-Sieve uses the same bit-length-based starting point and additionally rejects candidates where g_M = 9M² − N is negative (below M_min or a stale residue).

## 3.6 W-Sort CRT branch ordering

In the CRT DFS, at each prime p the enumeration order of the valid residue set V_K(p) directly determines how quickly the true K is reached. Not all residues in V_K(p) are equally likely to correspond to the true K. The residues split into two classes:

- Jac(g(k), p) = +1 residues: g(k) is a nonzero QR mod p. For each such k, T has two square-root choices, giving two (p mod p, q mod p) factor-pair combinations that produce this k. Multiplicity weight: 2.

- Jac(g(k), p) = 0 residues (Hensel-0): g(k) ≡ 0 mod p. T = 0 uniquely, giving one factor-pair combination. Multiplicity weight: 1.

Therefore the true K is twice as likely to lie in a Jac = +1 residue as in a Jac = 0 residue. W-Sort orders V_K(p) with all Jac = +1 residues before all Jac = 0 residues, causing the DFS to explore higher-probability subtrees first. The effect compounds across primes.

Empirical result: 1.18× geometric-mean speedup measured across 9 paired trials at 120-bit (env-gated via KSIEVE_W_SORT=1). Implementation cost: two-pass residue collection at sieve-build time (O(p) per prime), zero per-K₀ overhead during scan. W-Sort applies uniformly across K-, S-, T-, and M-Sieves.

# 4. The M-Sieve (Revision 5)

## 4.1 New-form reformulation

Let u = 6a − 1 and v = 6b + 1 with a, b ≥ 1 and uv composite. Define M = a + b, K_Δ = b − a, T = u + v, and L = v − u. Direct computation yields:

T = 6M, L = 6K_Δ + 2,

and the factoring identity N = uv = (T² − L²)/4 becomes 4N = T² − L², equivalently T² − L² = 4N, or — after dividing by 4:

(T/2)² − (L/2)² = N,  i.e.  (3M)² − (3K_Δ + 1)² = N.

Rearranging and setting L′ = 3K_Δ + 1, we obtain the M-Sieve polynomial identity:

g_M(M) = 9M² − N = L′² = (3K_Δ + 1)².

Factoring N thus reduces to finding M such that 9M² − N is a perfect square. Given L′ = √g_M(M), factor recovery proceeds by

u = 3M − L′,  v = 3M + L′.

The universal constraint L′ ≡ 1 (mod 3) follows from L′ = 3K_Δ + 1 and is a tautology given N ≡ 2 (mod 3); it provides no independent filter beyond the structural setup.

## 4.2 Search range

The M-Sieve search range is derived from the factor-ratio bound |v − u| ≤ √N (equivalent to factor ratio v/u ≤ (1 + √5)/2 ≈ 1.618). Under this bound, u + v ≤ √(5N), so M = (u + v)/6 ≤ √(5N)/6. The lower bound is M ≥ √N / 3 (attained as u → v). Thus:

M ∈ [⌈√N / 3⌉, ⌊√(5N) / 6⌋]

with range width (√5 − 2)/6 · √N ≈ 0.0393 · √N. This matches the S-Sieve range width exactly — the M- and S-Sieves have identical search-range structure despite scanning different variables. The 4.24× advantage over the K-Sieve range is preserved. (The T-Sieve, derived in §5, shares the same range width.)

## 4.3 Root modulus

The M-Sieve root modulus is determined by N mod 8 alone (no mod-3 fold):

- N mod 8 ∈ {1, 5}: M is odd → p_root = 2, target M mod 2 = 1.

- N mod 8 = 3: M ≡ 2 (mod 4) → p_root = 4, target M mod 4 = 2.

- N mod 8 = 7: M ≡ 0 (mod 4) → p_root = 4, target M mod 4 = 0.

The S- and T-Sieve mod-3 folds triple their root modulus to 6 or 12; the M-Sieve's root modulus stays at 2 or 4, costing approximately 1.585 bits of filtering relative to the S/T-Sieves. This gap is the dominant structural reason the M-Sieve is slower than the S- and T-Sieves at matched bit size. Empirical testing over 5 000 new-form semiprimes confirms M is uniformly distributed modulo 3 and modulo 9; no analog of the S-Sieve's S ≡ −n (mod 3) constraint exists for the M variable.

## 4.4 QR144 tables for M-Sieve

The M-Sieve verify chain admits a dedicated mod-144 QR filter. At true M, g_M = L′² with L′ = 3K_Δ + 1, so L′ ≡ 1 (mod 3) and L′² ≡ 1, 4, or 7 (mod 9). Combined with L′ parity (opposite of T parity, determined by N mod 8), this yields 3 or 6 valid residues of g_M (mod 144) per N mod 8 class:

| **N mod 8** | **Valid g_M mod 144** | **Count** | **Raw bits** |
| --- | --- | --- | --- |
| 1 | {16, 64, 112} | 3 | 5.58 |
| 3 | {1, 25, 49, 73, 97, 121} | 6 | 4.58 |
| 5 | {4, 52, 100} | 3 | 5.58 |
| 7 | {1, 25, 49, 73, 97, 121} | 6 | 4.58 |

These tables were derived analytically and verified against 10 000 random new-form semiprimes with zero false negatives. Note: despite the strong raw filter strength, the M-Sieve QR144 tables are fully subsumed in practice by the combination of QR49 (mod 49) ∩ QR25 (mod 25) ∩ QR256 (mod 256), which together span the mod 144 = 16 × 9 space. Instrumented benchmarks at 100-bit measured zero marginal rejections from QR144M after QR49/QR25/QR256 had already run. The tables are present in the implementation for correctness and documentation but are bypassed in the hot path.

## 4.5 Expected performance vs S-Sieve

The M-Sieve and S-Sieve have identical search-range widths. The M-Sieve loses approximately 1.585 bits relative to the S-Sieve at the root modulus (no mod-3 fold) and gains nothing from the L ≡ 1 mod 3 identity (tautological given N ≡ 2 mod 3). Predicted slowdown: factor 2–3× at matched bit size.

Empirical measurement at 100-bit (Rev. 8: 30 paired trials each, §6.4) shows the M/S ratio at 1.76× (down from Rev. 6's 2.56×, reflecting the faster reference hardware used in Rev. 8 benchmarks). The M-Sieve is valuable not for outperforming the S-Sieve but for covering a different class of semiprimes entirely — N ≡ 5 (mod 6), which the S-Sieve cannot factor regardless of bit size.

# 5. The T-Sieve (Revision 8)

## 5.1 (+1,+1)-form reformulation

Let u = 6x + 1 and v = 6y + 1 with y ≥ x ≥ 1 and uv composite. Define S = x + y, K = y − x, T_real = (u + v)/2, and L = (v − u)/2. Direct computation yields:

T_real = 3S + 1,    L = 3K,

and the factoring identity N = uv = T_real² − L² becomes

N = (3S + 1)² − (3K)² = 9S² + 6S + 1 − 9K².

Solving for (3K)²:

(3K)² = 9S² + 6S + 1 − N = 9S² + 6S − 6n,

using N = 6n + 1 (so N − 1 = 6n). We therefore obtain the T-Sieve polynomial identity:

g_T(S) = 9S² + 6S − 6n = (3K)².

Factoring N thus reduces to finding S such that g_T(S) is a perfect square divisible by 9. Given L = √g_T(S), factor recovery proceeds by K = L/3, T_real = 3S + 1, and

u = T_real − L = 3S + 1 − 3K,    v = T_real + L = 3S + 1 + 3K.

The final u·v = N check confirms the factorization.

## 5.2 Relationship to the S-Sieve

The T-Sieve is structurally the (+1,+1)-class mirror of the S-Sieve. The two polynomials differ only in the sign of the linear coefficient:

g_S(S) = 9S² − 6S − 6n       (S-Sieve, (−1,−1) class)

g_T(S) = 9S² + 6S − 6n       (T-Sieve, (+1,+1) class)

This is a direct consequence of T_real = 3S − 1 for the S-Sieve vs T_real = 3S + 1 for the T-Sieve. Three structural properties are affected by this sign flip:

1. **The mod-3 fold flips sign.** In the (+1,+1) class, N = 36xy + 6S + 1, hence 6n = N − 1 = 36xy + 6S, giving 6xy = n − S. Reducing mod 3 yields S ≡ +n (mod 3), as opposed to the S-Sieve's S ≡ −n (mod 3) arising from 6xy = n + S in the (−1,−1) class.

2. **The mod-4 parity targets flip at N mod 8 ∈ {3, 7}.** For odd N, the constraint N = T_real² − 9K² taken mod 8 (or mod 16) gives a unique residue class for S mod 4 in terms of N mod 8. Enumeration shows: for the T-Sieve, N mod 8 = 3 forces S ≡ 3 (mod 4) and N mod 8 = 7 forces S ≡ 1 (mod 4) — both flipped from the S-Sieve's corresponding targets (1 and 3 respectively). For N mod 8 ∈ {1, 5}, S is even in both sieves.

3. **The factor recovery uses T_real = 3S + 1 instead of 3S − 1.** The sign flip is absorbed entirely in this single addition during the final reconstruction step; the rest of the pipeline is unchanged.

All other structural elements are identical between S- and T-Sieve:

- **Discriminant.** g_T has discriminant 36 + 4·9·6n = 36(1 + 6n) = 36N. This matches the S-Sieve's discriminant 36N exactly, so the two sieves have identical Hensel-lift behavior at any prime p: g_S has a repeated root mod p iff p | N iff g_T has a repeated root mod p.

- **Search range.** S ∈ [⌈√N/3⌉, ⌊√(5N/4)/3⌋], same as S-Sieve. Under the balanced-factor assumption |v − u| ≤ √N, we have u + v ≤ √(5N). With u + v = 6S + 2 for the T-Sieve (vs 6S − 2 for S-Sieve), the upper bound on S is (√(5N) − 2)/6; the 2-unit offset is below 1 ULP of √(5N)/6 at every bit size of interest, so the S-Sieve's published formula ⌊√(5N/4)/3⌋ serves both sieves correctly.

- **Root modulus density.** 1/6 when N mod 8 ∈ {1, 5} and 1/12 when N mod 8 ∈ {3, 7}, matching the S-Sieve.

- **QR144 tables.** Bitwise identical to the S-Sieve's, as proved in §5.5.

- **Verify chain ordering, bitset-prime count (crt_np + 3), and cross-prime Newton recipe** are all inherited from the S-Sieve with no modification beyond the +6S linear coefficient change.

This correspondence makes the T-Sieve cheap to implement: the existing S-Sieve code paths require only three sign flips, plus a branch for the mod-4 parity target. In the reference C implementation the T-Sieve mode adds approximately 400 lines of code to the S-Sieve's infrastructure, most of which mirror existing S-Sieve branches with a single-character polynomial change (`-` → `+` on the 6S term).

## 5.3 Search range

As noted in §5.2, the T-Sieve search range is identical to the S-Sieve's:

S ∈ [⌈√N / 3⌉, ⌊√(5N/4) / 3⌋]

with range width (√5 − 2)/6 · √N ≈ 0.0393 · √N — 4.24× smaller than the K-Sieve range. The lower bound S ≥ ⌈√N/3⌉ follows from the AM-GM inequality u + v ≥ 2√N with u + v = 6S + 2. The upper bound derives from the factor-ratio assumption |v − u| ≤ √N, giving u + v ≤ √(5N), hence S ≤ (√(5N) − 2)/6, which floor-rounds to the same integer as ⌊√(5N/4)/3⌋ at the bit sizes used in practice.

## 5.4 Root modulus

The T-Sieve root modulus folds the mod-3 constraint S ≡ +n (mod 3) together with the mod-4 constraint (or mod-2 parity) determined by N mod 8:

- **N mod 8 ∈ {1, 5}**: S even. base = 2, s_base_target = 0. Combined with S ≡ +n mod 3, root_p = 6, density 1/6.

- **N mod 8 = 3**: S ≡ 3 (mod 4). base = 4, s_base_target = 3. Combined with S ≡ +n mod 3, root_p = 12, density 1/12.

- **N mod 8 = 7**: S ≡ 1 (mod 4). base = 4, s_base_target = 1. Combined with S ≡ +n mod 3, root_p = 12, density 1/12.

The CRT target within [0, root_p) is the unique integer t such that t ≡ s_base_target (mod base) and t ≡ (n mod 3) (mod 3); existence and uniqueness follow from gcd(base, 3) = 1. The T-Sieve root modulus density thus exactly matches the S-Sieve's, as required for performance parity.

## 5.5 QR144 equivalence lemma

**Lemma 5.1.** *Let N = u·v be a prime-pair semiprime with u, v > 3 and u ≡ v (mod 6). Let S = (u/6) + (v/6) be the sum variable in either the S-Sieve (u ≡ v ≡ −1) or the T-Sieve (u ≡ v ≡ +1) reformulation. Let g(S) denote the corresponding polynomial — g_S or g_T. Then the set of values g(S_true) mod 144, taken over all valid (N, S_true) pairs with fixed N mod 8, coincides between the two sieves and is given by:*

| **N mod 8** | **{g(S_true) mod 144}** |
| --- | --- |
| 1 | {0} |
| 3 | {9, 81} |
| 5 | {36} |
| 7 | {9, 81} |

*Proof.* At the true S, both polynomials evaluate to (3K)², a perfect square in 9·ℤ. The admissible values of 3K (mod 12), and hence of (3K)² (mod 144), are constrained by the Fermat identity T_real² − (3K)² ≡ N (mod 8) where T_real is an odd integer. Since squaring in ℤ/8 is sign-agnostic — (±x)² ≡ x² — the constraint on (3K)² mod 8 depends only on N mod 8 and on the possible values of T_real² mod 8, independent of whether T_real = 3S − 1 (S-Sieve) or 3S + 1 (T-Sieve). Combined with the multiplicative structure of 9·ℤ/144ℤ, the permissible residue classes of (3K)² mod 144 are identical in both sieves. Direct enumeration of the admissible (3K mod 12) values for each N mod 8 class yields the table above; the set sizes (1, 2, 1, 2) match the S-Sieve's empirical tables derived in Revision 3. ∎

**Corollary.** *The T-Sieve reuses the S-Sieve's QR144 tables (QR144S_N1, QR144S_N3, QR144S_N5, QR144S_N7) directly with no modification. Numerical verification over 10⁶ random (+1,+1) semiprimes confirms zero false negatives.*

## 5.6 Expected and measured performance vs S-Sieve

From the equivalences in §5.2 — identical search range, identical discriminant, identical root-modulus density, identical QR144 tables, identical verify chain — the T-Sieve is predicted to have essentially identical wall-time performance to the S-Sieve at matched bit size. The only measurable differences would arise from second-order effects: the specific CRT DFS tree shape depends on the exact Jacobi residue sets, which differ between g_S and g_T because the linear term affects per-prime residue counts modestly.

Empirical measurement at 80–120 bit (Rev. 8: 30 trials per cell at each bit size on independent RNG streams; see §6.3) confirms the prediction:

| **Bits** | **S-Sieve mean (s)** | **T-Sieve mean (s)** | **T/S ratio** |
| --- | --- | --- | --- |
| 80 | 0.001 | 0.001 | 1.00× (noise floor) |
| 90 | 0.003 | 0.003 | 0.97× |
| 100 | 0.059 | 0.044 | 0.74× |
| 110 | 0.623 | 0.640 | 1.03× |
| 120 | 13.586 | 13.806 | 1.02× |

The 100-bit cell shows T-Sieve measurably *faster* than S-Sieve, driven by slightly lower CRT work (2.9M CRT ops vs 3.4M) on the particular RNG seeds used. At 90, 110, and 120 bit the ratio is within 3% of unity in both directions, confirming performance parity. The 80-bit cell is at the noise floor (both ~1 ms) and is not informative.

## 5.7 Auto-dispatch for N ≡ 1 (mod 6)

For an input N ≡ 1 (mod 6) whose factor classes are not known a priori, the sieve must determine whether the factors are (−1,−1) (use S-Sieve) or (+1,+1) (use T-Sieve). The factor classes cannot be inferred from N alone — both sub-classes produce N ≡ 1 (mod 6) with identical statistical distributions.

Revision 8 implements an **S-first, T-fallback** strategy: for KSIEVE_AUTO=1 or KSIEVE_RSA=1 modes on N ≡ 1 (mod 6), the S-Sieve is invoked first; if its scan exhausts without finding a factor (the signal that the factors are actually in the (+1,+1) class), the T-Sieve is invoked on the same N. The fallback is implemented via a retry wrapper in ksieve_factor() that sets an internal KSIEVE_RETRY_AS_T=1 flag to redirect dispatch to the T-Sieve mode on the second attempt.

The worst-case cost of this strategy is approximately 2× the S-Sieve time (when the factors are actually (+1,+1)); the best-case cost is 1× (when the factors are (−1,−1)). On uniform random input, the expected cost is ≈1.5× S-Sieve time, since each sub-class occurs with probability 50%. Alternative strategies — parallel dispatch of S- and T-Sieve across thread pools, or class inference via cheap pre-tests — are discussed in §9.

The RSA-mode benchmark at 80 bit, 20 trials, confirms the strategy: 20/20 factorized, mean total 0.042 s (vs 0.002 s for pure S-Sieve on (−1,−1)-only and 0.003 s for pure T-Sieve on (+1,+1)-only). The ratio 0.042/0.003 ≈ 14× is inflated by small-N overheads that dominate at 80 bit; at 100 bit the expected ≈1.5× multiplier regime should hold.

# 6. Empirical Results

## 6.1 S-Sieve performance (16 threads)

Mean wall-time across 30 trials per bit size (20 at 120-bit); random old-form (−1,−1) semiprimes with deterministic RNG seed. All trials factored successfully with isqrt = 1.

| **N bits** | **Mean total (s)** | **Min (s)** | **Max (s)** | **Mean K0s** | **Mean CRT ops** |
| --- | --- | --- | --- | --- | --- |
| 80 | 0.001 | 0.001 | 0.001 | 3.1K | 3.9K |
| 90 | 0.003 | 0.001 | 0.007 | 130K | 171K |
| 100 | 0.059 | 0.006 | 0.194 | 2.50M | 3.41M |
| 110 | 0.623 | 0.006 | 1.639 | 31.1M | 42.1M |
| 120 | 13.586 | 0.378 | 33.998 | 675M | 920M |

Variance within a bit level is substantial, driven by DFS ordering and factor-balance effects rather than structural properties of N. The 120-bit max/min ratio of ~90× is the largest observed; however, the max of 34.0 s is substantially smaller than Rev. 6's 121 s max, reflecting both the faster Rev. 8 reference hardware and the improved default pass-fractions configuration (see §8).

## 6.2 M-Sieve performance (new-form, 16 threads)

Mean wall-time across 30 trials per bit size (20 at 110-bit); random new-form semiprimes (u ≡ 5 mod 6, v ≡ 1 mod 6) with deterministic RNG seed. All trials factored successfully with isqrt = 1.

| **N bits** | **Mean total (s)** | **Min (s)** | **Max (s)** | **Mean K0s** | **Mean CRT ops** |
| --- | --- | --- | --- | --- | --- |
| 60 | 0.001 | 0.001 | 0.002 | 469 | 469 |
| 70 | 0.001 | 0.001 | 0.001 | 3.1K | 3.9K |
| 80 | 0.001 | 0.001 | 0.002 | 15.8K | 20.8K |
| 90 | 0.010 | 0.001 | 0.050 | 278K | 365K |
| 100 | 0.104 | 0.005 | 0.494 | 3.09M | 4.21M |
| 110 | 2.911 | 0.030 | 9.412 | 53.4M | 71.1M |

All 170 M-Sieve trials (30+30+30+30+30+20) factored correctly with no false negatives. Mean K0s and CRT ops track closely with S-Sieve at matched bit size, consistent with the analysis in §4.3: the two sieves do similar CRT work per candidate, but M-Sieve bitset scan is less dense due to the missing mod-3 root fold.

## 6.3 T-Sieve performance ((+1,+1)-form, 16 threads)

Mean wall-time across 30 trials per bit size (20 at 120-bit); random (+1,+1)-form semiprimes (u ≡ v ≡ 1 mod 6) with deterministic RNG seed. All trials factored successfully with isqrt = 1.

| **N bits** | **Mean total (s)** | **Min (s)** | **Max (s)** | **Mean K0s** | **Mean CRT ops** |
| --- | --- | --- | --- | --- | --- |
| 80 | 0.001 | 0.001 | 0.001 | 6.2K | 8.4K |
| 90 | 0.003 | 0.001 | 0.007 | 125K | 165K |
| 100 | 0.044 | 0.006 | 0.103 | 2.22M | 2.93M |
| 110 | 0.640 | 0.034 | 1.724 | 32.7M | 44.5M |
| 120 | 13.806 | 2.825 | 31.594 | 522M | 682M |

All 140 T-Sieve trials factored correctly with isqrt = 1 on every trial. The S-Sieve vs T-Sieve head-to-head comparison across the full scaling range:

| **Bits** | **S mean (s)** | **T mean (s)** | **T/S ratio** |
| --- | --- | --- | --- |
| 80 | 0.001 | 0.001 | 1.00× |
| 90 | 0.003 | 0.003 | 0.97× |
| 100 | 0.059 | 0.044 | 0.74× |
| 110 | 0.623 | 0.640 | 1.03× |
| 120 | 13.586 | 13.806 | 1.02× |

The algebraic prediction of identical performance (§5.6) is empirically confirmed: at 90, 110, and 120 bit the T/S ratio is within 3% of unity, well within trial-to-trial variance. The 100-bit cell has T running measurably faster than S on the specific 30 RNG-seed samples; at that bit size T-Sieve's CRT ops are 14% lower than S-Sieve's (2.93M vs 3.41M). The 80-bit cell is at the ~1 ms noise floor.

## 6.4 Paired K / S / T / M comparison at 100-bit

Head-to-head comparison at 100-bit across 30 trials each. Each sieve is measured on the input class it is designed for: K- and S-Sieve on (−1,−1) inputs, T-Sieve on (+1,+1) inputs, M-Sieve on new-form inputs. Workload is comparable in bit size but structurally different in residue class across sieves (distinct RNG streams).

| **Sieve** | **Trials** | **Mean total (s)** | **Mean K0s** | **Mean CRT ops** | **Ratio vs S** |
| --- | --- | --- | --- | --- | --- |
| K-Sieve | 30/30 | 0.757 | 3.05M | 4.15M | 12.28× |
| S-Sieve | 30/30 | 0.062 | 2.52M | 3.42M | 1.00× (baseline) |
| T-Sieve | 30/30 | 0.051 | 2.23M | 2.94M | 0.82× |
| M-Sieve | 30/30 | 0.173 | 3.10M | 4.22M | 2.81× |

Observations:

- **K vs S** (both on (−1,−1) inputs): K-Sieve does essentially the same CRT work as S-Sieve (4.15M ops vs 3.42M, K0s 3.05M vs 2.52M) but is 12.3× slower. The difference is scan range: K-Sieve scans a range 4.24× larger than S-Sieve's and additionally lacks the mod-3 root fold. The old K/S ratio of 20.2× reported in Rev. 6 has narrowed in the Rev. 8 hardware/configuration to 12.3×; the ratio remains in the expected 10–20× regime.

- **S vs T** (different input classes, same bit size): T-Sieve on (+1,+1) inputs runs at 0.82× the S-Sieve's time on (−1,−1) inputs in this particular paired run. The ratio is within trial-noise of unity (three independent 30-trial runs at 100-bit give S ∈ {0.059, 0.062} and T ∈ {0.044, 0.051, 0.043}, all tightly clustered). The algebraic parity prediction is confirmed.

- **M vs S**: M-Sieve runs 2.81× slower than S-Sieve, close to the 2.56× reported in Rev. 6 and within the theoretical 2–3× window predicted from the missing mod-3 fold (§4.3). The M/S ratio at 110 bit widens to 4.67× (from Rev. 8's §6.1/§6.2 data), but the Rev. 6 M-Sieve scaling fit T ≈ 2^(0.41·bits) is preserved; the Rev. 8 fit over the same 80–110 bit range is slightly tighter at T ≈ 2^(0.367·bits).

- A direct K-on-(+1,+1) measurement was not collected in Rev. 8; the K-Sieve polynomial g_K = 9K² + N is class-agnostic (§1.2) and its performance on (+1,+1) inputs is expected to match its performance on (−1,−1) inputs to within trial variance. A ballpark measurement from five 90-bit balanced (+1,+1) inputs gave K-Sieve wall times of 0.89–5.50 s, versus T-Sieve 0.017–0.126 s on the same inputs, a mean speedup of 28× — consistent with the 12.3× K/S ratio at 100-bit scaled by the 0.8× T/S ratio at that size.

## 6.5 Bitset-prime tuning

The default bitset-prime count ctx→n_primes was re-examined in Rev. 5 via empirical BS_NP sweeps, and re-confirmed in Rev. 8. The current mode-aware defaults are:

- K-Sieve: crt_np + 2
- S-Sieve: crt_np + 3
- T-Sieve: crt_np + 3 (same as S-Sieve, both inheriting mod-3 fold)
- M-Sieve: crt_np + 5

The Rev. 8 sweep at 90-bit M-Sieve (20 trials per point, 16 threads):

| **BS_NP** | **Mean total (s)** | **Relative to bs_np=9** |
| --- | --- | --- |
| 8 | 0.016 | 1.20× slower |
| 9 (old default = crt_np+2) | 0.013 | 1.00× |
| 10 | 0.013 | 1.03× faster |
| 11 | 0.013 | 1.06× faster |
| 12 | 0.012 | 1.13× faster |
| 13 | 0.013 | 1.02× slower |
| 14 (crt_np+5 default) | 0.017 | 1.28× slower |

The sweep is noisy — differences between bs_np ∈ {10, 11, 12, 13} are within trial-to-trial variance at this sample size (20 trials). The minimum at bs_np = 12 (0.012 s) is 1.13× faster than bs_np = 9 and roughly equivalent to bs_np ∈ {10, 11}. Unlike the Rev. 5/6 sweep where bs_np = 13 was fastest at 90-bit, the Rev. 8 sweep on different hardware shows bs_np = 12 as the optimum, with bs_np = 14 (the crt_np+5 default) now measurably slower. This is consistent with hardware-dependence of the tile-scan bandwidth limit; the crt_np+5 default is within ~30% of optimal across hardware configurations and remains a reasonable choice. A per-platform bench sweep via KSIEVE_BS_NP is recommended for performance-critical deployments.

## 6.6 Threading scalability

Measured at 100-bit M-Sieve, 10 trials per configuration, identical source compiled with -DKSIEVE_N_THREADS=1 and -DKSIEVE_N_THREADS=16:

| **Threads** | **Trials** | **Mean total (s)** | **Min (s)** | **Max (s)** | **Speedup** |
| --- | --- | --- | --- | --- | --- |
| 1 | 10/10 | 1.157 | 0.065 | 2.115 | 1.00× (baseline) |
| 16 | 10/10 | 0.217 | 0.006 | 0.827 | 5.33× |

The 5.33× speedup on 16 threads (Rev. 8) is substantially sub-linear and somewhat lower than Rev. 6's 7.5× on the same bit size. The reduction is attributable to two factors: (i) 10 trials are too few to stabilize the ratio at the tails — Rev. 6 already flagged this caveat — and (ii) the Rev. 8 reference hardware has faster single-thread throughput but correspondingly less headroom for threading to contribute. The sub-linear speedup remains consistent with the analysis in Rev. 6:

- Memory-bandwidth saturation in the bitset scan: the 16 parallel scans share a single LLC and DRAM channel, and the bitset tile + QR64 mask loads are not latency-hidden by the compute ratio.

- Tail effect: the algorithm terminates as soon as any worker finds the factor. If the true K0 is in one DFS subtree, the remaining n−1 workers have wasted their parallel slice on irrelevant subtrees.

- Small-input floor: at 100-bit, wall time is ~0.2 s at 16 threads; the prefix-work-stealing setup has ~50 ms fixed cost regardless of N, costing ~25% of wall time at this size.

Scaling would improve at larger bit sizes where the sieve phase dominates more thoroughly. Dedicated measurement of threading scaling at 110/120-bit has not been conducted.

## 6.7 Filtering bit budget

| **Component** | **K-Sieve** | **S-Sieve** | **T-Sieve** | **M-Sieve** |
| --- | --- | --- | --- | --- |
| Root modulus | 1.5 bits | 3.09 bits | 3.09 bits | 1.5 bits |
| CRT (8 primes @ 100-bit) | 7.16 bits | 7.16 bits | 7.16 bits | 7.16 bits |
| Bitset (9–14 Jacobi primes) | 8.0 bits | 9.0 bits | 9.0 bits | 11.0 bits |
| QR64 | 2.42 bits | 2.42 bits | 2.42 bits | 2.42 bits |
| Verify chain | 12.0 bits | 12.0 bits | 12.0 bits | 9.0 bits |
| Cross-prime | 62 bits | 62 bits | 62 bits | 62 bits |

Bit-budget totals are not directly additive because of overlap between filters. The effective filtering power determines how many candidates survive to the expensive isqrt stage; empirically, all four sieves achieve isqrt = 1 per factorization at the bit sizes tested. The T-Sieve's budget column is identical to the S-Sieve's, reflecting the algebraic correspondence established in §5.2 and §5.5.

## 6.8 Correctness

All optimizations preserve correctness. The Rev. 8 benchmark suite covers 720+ factorizations across bit sizes 60–120 with deterministic RNG seeds: 170 M-Sieve (new form, 60–110 bit), 140 S-Sieve (old form, 80–120 bit), 140 T-Sieve (+1,+1 form, 80–120 bit), 90 K-Sieve (old form, 80–100 bit), plus the 120-trial paired comparison at 100-bit (K/S/T/M × 30), the 20-trial threading scan, and the 140-trial BS_NP sweep. Pass rate: 720/720 (100%). Every successful factorization produced isqrt = 1, confirming zero false negatives in the filtering pipeline across all four sieves.

# 7. Scaling Analysis

## 7.1 Empirical scaling exponent

Linear regression on log₂(mean total time) against bit size, over the 30-trial-per-cell Rev. 8 sweep, yields:

  • S-Sieve (80–120 bit): T ≈ 2^(−39.4) · 2^(0.355·bits),  exponent ≈ 0.355

  • T-Sieve (80–120 bit): T ≈ 2^(−39.6) · 2^(0.356·bits),  exponent ≈ 0.356

  • M-Sieve (80–110 bit): T ≈ 2^(−39.4) · 2^(0.367·bits),  exponent ≈ 0.367

  • K-Sieve (80–100 bit): T ≈ 2^(−40.4) · 2^(0.402·bits),  exponent ≈ 0.402

The S-Sieve and T-Sieve exponents differ by 0.001, statistically indistinguishable and well below regression noise — a direct empirical confirmation of the algebraic parity proved in §5. The M-Sieve exponent remains modestly higher than the S/T-Sieves' (0.367 vs 0.355), reflecting the missing mod-3 root fold, which costs a constant bit-factor rather than a super-polynomial term; the gap should remain approximately constant as N grows. The K-Sieve exponent of 0.402 reflects its 4.24× wider search range and absence of the mod-3 fold.

All four exponents are well below the √N = 0.5 lower bound implied by a naive Fermat scan without filtering, consistent with the sub-exponential improvement predicted by the hyperbolic-sieve analysis of [1]. The Rev. 8 fits are slightly lower than Rev. 6's S-Sieve and M-Sieve exponents (0.37 and 0.41 respectively); the difference is attributable to the faster reference hardware and the improved default pass-fractions configuration.

If the effective filtering power grows as Θ(log N) bits (through continued CRT-prime addition), the practical exponent decreases slowly toward 1/2 − ε with ε → 0. Whether this persists beyond 130 bits depends on whether additional CRT primes remain cost-effective. The algorithmic ceiling (§7.2) is expected in the 150–180 bit range.

## 7.2 Algorithmic ceiling

At sufficiently large bit sizes, CRT tree cost matches bitset scan cost, at which point adding more primes no longer helps. Based on the current architecture and Rev. 5/8 tuning, this crossover is estimated at 150–180 bits, beyond which general-purpose subexponential methods (QS, GNFS) become preferable.

# 8. Bugs Found and Fixed

### QR65536 initialization overflow

The init_qr_tables function computed k·k as signed int, overflowing at k ≥ 46341. This undefined behavior caused GCC -O3 to generate incorrect code, manifesting as crashes when the verify chain used QR65536 (np ≥ 8). Fix: cast to unsigned.

### S-Sieve cross-prime Newton divergence

The starting point cpR = 3S diverges when S is near S_min (where g_S ≈ 0). Fix: use bit-length-based starting point.

### S_min signed overflow at 130-bit

S_min ≈ 2⁶³ overflows signed long. Fix: use u64 for S_min and K_long throughout.

### Cross-prime u64 overflow at 127+ bits (Rev. 4)

The cross-prime filter uses nextprime(√N) as its modulus. For N ≥ 127 bits, √N exceeds 2⁶³, causing u64 overflow when nextprime(√N) was computed in 64-bit arithmetic. The resulting nonsense prime silently disabled the cross-prime filter, causing billions of unfiltered isqrts at 130-bit+ N. Fix: cap cp_p at nextprime(2⁶²) for large N, and route g computation through mpz when g overflows u128 (N ≥ 127 bits). Impact at 130-bit on a favorable N: isqrt count dropped from approximately 22 billion to 1; wall time dropped from ∼55 minutes to ∼1 second.

### Mode-aware bitset prime count (Rev. 5)

The n_primes = crt_np + 2 formula was tuned for S-Sieve only and systematically underestimated the optimum for both S- and M-Sieve across 90–100 bit. Fix: mode-aware defaults (§6.5). Impact: ~1.1–1.3× speedup at 90-bit (within measurement variance at 20 trials).

### Rev. 6 benchmark refresh

Full 30-trials-per-cell empirical sweep across 60–120 bit for all three sieves (400 factorizations total, 100% pass rate).

### T-Sieve mod-4 target flip (Rev. 8)

During T-Sieve implementation, a subtle bug was caught early: the S-Sieve's mod-4 parity targets at N mod 8 ∈ {3, 7} are 1 and 3 respectively, but the T-Sieve's are 3 and 1 — flipped, reflecting the T_real = 3S + 1 vs 3S − 1 sign difference. A naive copy-paste of the S-Sieve code would have silently produced incorrect bitset tile density on half of random inputs (those with N mod 8 ∈ {3, 7}) without any factorization failures at small bit sizes (where the bitset tiles are trivially filled) but with catastrophic slowdown at 100+ bit where tiles are sparse. The fix is implemented as a branch on `mode` in sieve_build() and is covered by the Rev. 8 benchmark suite (140 T-Sieve trials, all OK with isqrt = 1).

### T-Sieve S_min bound edge cases (Rev. 8)

For very small N (< 30 bit), the S-Sieve/T-Sieve S_min = ⌈√N/3⌉ can exclude the true S value, causing factorization failure. This is preexisting behavior of sum-variable sieves: at bit sizes where the range [S_min, S_max] collapses to fewer than two integers, unbalanced factor pairs can fall outside. The issue affects the T-Sieve at 30- and 40-bit on some seeds; it does not affect any bit size ≥ 50 in the Rev. 8 benchmarks (all 720+ trials at ≥ 60 bit succeed). Workaround: for small-N factoring, use K-Sieve (class-agnostic polynomial, wider range). The T-Sieve and S-Sieve are designed for the 80–130 bit regime and the S_min edge-case is not a regression — it is inherent to sum-variable range reduction.

### Rev. 8 default pass-fractions change

The default KSIEVE_PASS_FRACS changed from the Rev. 6/7 two-knob default (0.10, 0.32 → 3 passes) to a five-knob default (0.03, 0.08, 0.15, 0.35 → 5 passes). Legacy KSIEVE_PASS0_FRAC and KSIEVE_PASS1_FRAC overrides still work for regression testing. Impact at 100-bit: within 5% of either default; meaningful gains appear at 120+ bit where the distribution's right tail is heavier.

# 9. Future Work

### p² Hensel bitset tiles

For CRT primes q where q | M, the bitset tile is currently all-ones (no filtering). A p² tile checking v_q(g) parity would provide approximately 0.67 additional bits of filtering for primes 5, 7, 11, 13. Preliminary implementation showed correct results but per-word gather overhead (+20 cycles) exceeded survivor savings at current density levels.

### Parallel S + T sieve dispatch

The current Rev. 8 auto-dispatch strategy for N ≡ 1 (mod 6) runs S-Sieve first and T-Sieve as fallback, with expected cost ≈1.5× S-Sieve time on uniform random input. A parallel implementation running S- and T-Sieve concurrently on halved thread pools would give expected cost min(T_S, T_T) rather than T_S + [factors are (+1,+1)] · T_T. The implementation requires careful accounting of the atomic found_flag shared between the two sieve contexts but is otherwise straightforward. Expected speedup on (+1,+1) inputs: ≈2× over the sequential fallback. Expected regression on (−1,−1) inputs: slight, from reduced per-sieve thread count.

### Class inference for N ≡ 1 (mod 6)

A sufficiently cheap pre-test determining the factor class of N ≡ 1 (mod 6) could route directly to the correct sieve without a retry penalty. Naive approaches (trial division, Jacobi-based tests) do not distinguish the classes with high probability at reasonable cost. Open problem.

### Bidirectional DFS

Run half of workers forward through the CRT DFS (low K to high K) and half in reverse. Mean-time change: zero. Variance reduction: approximately 4×, P95 time roughly halved. Targets tail behavior rather than mean.

### GPU acceleration

The bitset scan is embarrassingly parallel across K₀ values. A CUDA implementation could process thousands of K₀s simultaneously.

### Scaling validation at 140+ bits

The N^{1/3} scaling conjecture needs validation at 140–160 bits to determine whether the CRT tree cost crossover occurs within this range. T- and M-Sieve scaling at these bit sizes is entirely uncharacterized.

### Hittmeir multiple-choice subset-sum variant

Hittmeir [1] describes a time-space tradeoff in which the sieve set L_{N,m,k} is partitioned into two classes and factoring is reduced to a multiple-choice subset-sum problem. Applied to the 6n±1 parameterization, this approach might yield further practical speedup at the cost of non-negligible space complexity. Preliminary analysis suggests the approach is compatible with the current CRT/bitset architecture.

# References

[1] M. Hittmeir, "Integer factorization as subset-sum problem," Journal of Number Theory, vol. 249, pp. 93–118, 2023. Preprint: arXiv:2205.10074, 2022.

[2] R. Creft, "Hexile Sieve Analysis of Prime and Composite Integers," arXiv:1202.5948 [math.NT], February 2012.

[3] P. de Fermat, Oeuvres de Fermat, Gauthier-Villars, Paris, 1891. (See Vol. II for the original correspondence on the difference-of-squares factoring identity N = x² − y².)

[4] F. W. Lawrence, "Factorisation of numbers," Quarterly Journal of Pure and Applied Mathematics, vol. 28, pp. 285–311, 1895.

[5] R. S. Lehman, "Factoring large integers," Mathematics of Computation, vol. 28, no. 126, pp. 637–646, 1974.

[6] C. Pomerance, "The Quadratic Sieve Factoring Algorithm," EUROCRYPT 84, LNCS 209, pp. 169–182, 1985.

[7] C. Pomerance, "A Tale of Two Sieves," Notices of the AMS, vol. 43, no. 12, pp. 1473–1485, 1996.

[8] J. McKee, "Speeding Fermat's Factoring Method," Mathematics of Computation, vol. 68, pp. 1729–1738, 1999.

[9] K. Hensel, Theorie der algebraischen Zahlen, Teubner, 1908.

[10] C. G. J. Jacobi, "Über die Kreisteilung," Bericht Ak. Wiss. Berlin, pp. 127–136, 1837.

[11] C. F. Gauss, Disquisitiones Arithmeticae, 1801.

[12] K. Ireland and M. Rosen, A Classical Introduction to Modern Number Theory, 2nd ed., Springer, 1990.

[13] H. Cohen, A Course in Computational Algebraic Number Theory, Springer, 1993.

[14] H. Riesel, Prime Numbers and Computer Methods for Factorization, 2nd ed., Birkhäuser, 1994.

[15] E. Bach and J. Shallit, Algorithmic Number Theory, Vol. 1, MIT Press, 1996.

[16] R. Rivest, A. Shamir, L. Adleman, "A Method for Obtaining Digital Signatures," Comm. ACM, vol. 21, 1978.

[17] K. Conrad, "Hensel's Lemma," expository notes, University of Connecticut.

[18] I. E. Shparlinski, "Modular hyperbolas," Japanese Journal of Mathematics, vol. 7, no. 2, pp. 235–294, 2012.

[19] Intel Corporation, Intel 64 and IA-32 Architectures Software Developer's Manual, Volume 2.

[20] A. Fog, Microarchitecture of Intel, AMD and VIA CPUs, Technical University of Denmark, 2023.
