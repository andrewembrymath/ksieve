#!/usr/bin/env bash
# =============================================================================
# ksieve benchmark gathering script
# =============================================================================
# Runs a full publication-quality benchmark suite overnight.
#
# What it does:
#   1. Builds binaries for 1-thread and 16-thread variants (from the SAME source)
#   2. Runs M-sieve scaling sweep (60..110 bit, 30 trials each, 16-thread)
#   3. Runs S-sieve scaling sweep (80..120 bit, 30 trials each, 16-thread)
#   4. Runs T-sieve scaling sweep (80..120 bit, 30 trials each, 16-thread)  [Rev. 8+]
#   5. Runs K-sieve scaling sweep (80..100 bit, 30 trials each, 16-thread)
#   6. Runs paired K vs S vs T vs M comparison at 100-bit (30 trials each)
#   7. Runs single-thread vs 16-thread comparison for M-sieve at 100-bit
#   8. Runs BS_NP tuning sweep at 90-bit (reproduces Rev. 5 result)
#
# Output:
#   - Raw per-trial logs in $OUT/logs/*.log (one file per run)
#   - Parsed CSV in $OUT/results.csv
#   - Markdown tables in $OUT/tables.md (ready to paste into paper)
#   - Progress log in $OUT/progress.log
#
# Resumable: if killed, re-running will skip any run whose log already exists
# with a "--- Summary" line (meaning it completed).
#
# Expected wall time on 16-core 3GHz+ machine: 10-14 hours total (Rev. 8).
#   Most of it is the two 120-bit cells (S-sieve and T-sieve, each ~45 min at
#   ~90s/trial × 30 trials) plus M-sieve 110-bit (6+ hours worst case).
#   Adding T-sieve roughly doubles the S-sieve section's wall time.
# You can kill and resume at any time.
#
# Usage:
#   ./gather_benchmarks.sh [OUTPUT_DIR] [SRC_PATH]
#
#   OUTPUT_DIR  defaults to ./bench_out_<timestamp>
#   SRC_PATH    defaults to ./main.c
#
# Tips:
#   - Run under `nohup ... &` or in a tmux/screen session
#   - `tail -f $OUT/progress.log` to watch progress
#   - To skip a section, create its marker file (e.g. `touch $OUT/.done.scaling-M`)
# =============================================================================

set -u   # but NOT -e — we want to continue past individual trial failures

# ---------- configuration ----------
SRC_PATH="${2:-./main.c}"
OUT_DIR="${1:-./bench_out_$(date +%Y%m%d_%H%M)}"

# Edit these arrays to adjust the sweep.
# Format: "<bits> <trials> <timeout_per_trial_seconds>"
# Timeouts are per-trial within ksieve's own subprocess; the bench wrapper will
# simply record timeouts as non-OK trials. Keep timeouts generous.
M_SIEVE_GRID=(
    "60 30 30"
    "70 30 30"
    "80 30 60"
    "90 30 120"
    "100 30 600"
    "110 20 3600"
)

S_SIEVE_GRID=(
    "80 30 30"
    "90 30 60"
    "100 30 300"
    "110 30 1800"
    "120 20 7200"
)

# T-sieve grid mirrors S-sieve grid. Both have identical search range,
# root modulus density (1/6 or 1/12), discriminant (36N), and filter chain;
# expected performance parity with S-sieve at each bit size. Running the same
# grid gives a direct apples-to-apples comparison for the paper.
T_SIEVE_GRID=(
    "80 30 30"
    "90 30 60"
    "100 30 300"
    "110 30 1800"
    "120 20 7200"
)

K_SIEVE_GRID=(
    "80 30 30"
    "90 30 120"
    "100 30 600"
)

# Bit size for the K-vs-S-vs-M paired comparison (Table 5.3 in paper).
# Single bit for cleanest headline number; paper can pick 100.
PAIRED_BITS=100
PAIRED_TRIALS=30

# Bit size for the single-thread vs 16-thread comparison.
THREADING_BITS=100
THREADING_TRIALS=10

# BS_NP sweep (Rev. 5 reproducibility)
BSNP_BITS=90
BSNP_TRIALS=20
BSNP_MODE="KSIEVE_MSIEVE"   # MSIEVE for new-form; flip to SSIEVE for old-form sweep

# Threading-variant counts to build and compare
THREAD_VARIANTS=(1 16)

# Compile flags
CFLAGS="-O3 -march=native"
LIBS="-lgmp -lm -lpthread"
CC="gcc"

# ---------- setup ----------
mkdir -p "$OUT_DIR/logs" "$OUT_DIR/bin"
PROGRESS="$OUT_DIR/progress.log"
CSV="$OUT_DIR/results.csv"
TABLES="$OUT_DIR/tables.md"

log()  { echo "[$(date +%H:%M:%S)] $*" | tee -a "$PROGRESS"; }
banner() {
    echo ""                                | tee -a "$PROGRESS"
    echo "============================================================" | tee -a "$PROGRESS"
    echo "=== $*"                          | tee -a "$PROGRESS"
    echo "============================================================" | tee -a "$PROGRESS"
}

if [[ ! -f "$SRC_PATH" ]]; then
    log "ERROR: source file not found: $SRC_PATH"
    exit 1
fi

# Initialize CSV header once (safe to repeat — the same header every time)
if [[ ! -f "$CSV" ]]; then
    cat > "$CSV" <<EOF
run_id,sieve,bits,trials_requested,trials_ok,trials_fail,trials_timeout,mean_p1_s,mean_total_s,min_s,max_s,mean_crt_ops,mean_k0s,mean_isqrt,threads,extra_env,timestamp
EOF
fi

# ---------- build variants ----------
build_variant() {
    local threads="$1"
    local out="$OUT_DIR/bin/ksieve_t${threads}"
    if [[ -x "$out" ]]; then
        log "build: t${threads} already exists at $out"
        return 0
    fi
    log "build: t${threads} → $out"
    $CC $CFLAGS -DKSIEVE_N_THREADS="$threads" -o "$out" "$SRC_PATH" $LIBS 2> "$OUT_DIR/bin/ksieve_t${threads}.build.log"
    if [[ ! -x "$out" ]]; then
        log "ERROR: build failed for t${threads}; see $OUT_DIR/bin/ksieve_t${threads}.build.log"
        tail -20 "$OUT_DIR/bin/ksieve_t${threads}.build.log"
        return 1
    fi
    return 0
}

banner "BUILD: compiling threading variants"
for t in "${THREAD_VARIANTS[@]}"; do
    build_variant "$t" || { log "FATAL: cannot proceed without t${t} binary"; exit 2; }
done

# ---------- trial runner ----------
# Args:
#   $1 run_id (used as log filename and CSV key)
#   $2 sieve name (K|S|M) — informational only
#   $3 bits
#   $4 trials
#   $5 threads (must be one of THREAD_VARIANTS)
#   $6 extra env vars to pass (string like "KSIEVE_MSIEVE=1" or "KSIEVE_SSIEVE=1 KSIEVE_BS_NP=11")
run_bench() {
    local run_id="$1" sieve="$2" bits="$3" trials="$4" threads="$5" extra_env="$6"
    local logfile="$OUT_DIR/logs/${run_id}.log"
    local bin="$OUT_DIR/bin/ksieve_t${threads}"

    # Resumability: skip if log already has a valid Summary line
    if [[ -f "$logfile" ]] && grep -q "^--- Summary" "$logfile" 2>/dev/null; then
        log "skip: $run_id (already complete)"
        return 0
    fi

    log "RUN:  $run_id  [sieve=$sieve bits=$bits trials=$trials threads=$threads env='$extra_env']"
    local t0=$(date +%s)

    # shellcheck disable=SC2086   # deliberate word splitting of extra_env
    ( env $extra_env "$bin" --bench "$bits" "$trials" ) > "$logfile" 2>&1
    local rc=$?
    local t1=$(date +%s)
    local dur=$((t1 - t0))

    if ! grep -q "^--- Summary" "$logfile"; then
        log "WARN: $run_id  no Summary line (rc=$rc, $dur s); will not parse"
        return 1
    fi

    # Parse the summary block into CSV
    parse_summary "$run_id" "$sieve" "$bits" "$trials" "$threads" "$extra_env" "$logfile"
    log "DONE: $run_id  in ${dur}s"
    return 0
}

# parse_summary $run_id $sieve $bits $trials_requested $threads $extra_env $logfile
# Extracts the "Summary" block fields and appends one CSV row.
parse_summary() {
    local run_id="$1" sieve="$2" bits="$3" trials_req="$4" threads="$5" extra_env="$6" logfile="$7"

    # Example format we parse:
    #   --- Summary (10 OK / 0 FAIL / 0 TIMEOUT / 10 total) ---
    #     Mean p1    : 2.3179 s
    #     Mean total : 2.3185 s
    #     Min / Max  : 0.3003 / 5.4795 s
    #     Mean CRT   : 3057373
    #     Mean K0s   : 2271752
    #     Mean isqrt : 1
    local line
    line=$(grep "^--- Summary" "$logfile" | head -1)
    local ok fail timeout total
    ok=$(     echo "$line" | sed -n 's/.*(\([0-9]*\) OK.*/\1/p')
    fail=$(   echo "$line" | sed -n 's/.*\/ \([0-9]*\) FAIL.*/\1/p')
    timeout=$(echo "$line" | sed -n 's/.*\/ \([0-9]*\) TIMEOUT.*/\1/p')
    total=$(  echo "$line" | sed -n 's/.*\/ \([0-9]*\) total.*/\1/p')

    local mean_p1    mean_total   min_s      max_s      mean_crt     mean_k0    mean_isqrt
    mean_p1=$(   awk '/Mean p1/    {print $4; exit}' "$logfile")
    mean_total=$(awk '/Mean total/ {print $4; exit}' "$logfile")
    # "Min / Max  : 0.0040 / 0.1808 s"  →  $5 is min, $7 is max
    # (fields are: Min[$1]  /[$2]  Max[$3]  :[$4]  0.0040[$5]  /[$6]  0.1808[$7]  s[$8])
    min_s=$(     awk '/Min \/ Max/ {print $5; exit}' "$logfile")
    max_s=$(     awk '/Min \/ Max/ {print $7; exit}' "$logfile")
    mean_crt=$(  awk '/Mean CRT/   {print $4; exit}' "$logfile")
    mean_k0=$(   awk '/Mean K0s/   {print $4; exit}' "$logfile")
    mean_isqrt=$(awk '/Mean isqrt/ {print $4; exit}' "$logfile")

    local ts; ts=$(date -Iseconds)
    # Escape any commas in extra_env (shouldn't have any, but defensive)
    local env_escaped="${extra_env//,/;}"

    echo "$run_id,$sieve,$bits,$trials_req,${ok:-0},${fail:-0},${timeout:-0},${mean_p1:-NA},${mean_total:-NA},${min_s:-NA},${max_s:-NA},${mean_crt:-NA},${mean_k0:-NA},${mean_isqrt:-NA},$threads,\"$env_escaped\",$ts" >> "$CSV"
}

# ---------- SECTION 1: M-sieve scaling ----------
banner "SECTION 1: M-sieve scaling (new-form, 16-thread)"
if [[ ! -f "$OUT_DIR/.done.scaling-M" ]]; then
    for row in "${M_SIEVE_GRID[@]}"; do
        read -r bits trials _timeout <<< "$row"
        run_bench "scale-M-${bits}b" "M" "$bits" "$trials" 16 "KSIEVE_MSIEVE=1" || true
    done
    touch "$OUT_DIR/.done.scaling-M"
fi

# ---------- SECTION 2: S-sieve scaling ----------
banner "SECTION 2: S-sieve scaling (old-form, 16-thread)"
if [[ ! -f "$OUT_DIR/.done.scaling-S" ]]; then
    for row in "${S_SIEVE_GRID[@]}"; do
        read -r bits trials _timeout <<< "$row"
        run_bench "scale-S-${bits}b" "S" "$bits" "$trials" 16 "KSIEVE_SSIEVE=1" || true
    done
    touch "$OUT_DIR/.done.scaling-S"
fi

# ---------- SECTION 3: T-sieve scaling (Rev. 8+, (+1,+1) form) ----------
banner "SECTION 3: T-sieve scaling ((+1,+1)-form, 16-thread)"
if [[ ! -f "$OUT_DIR/.done.scaling-T" ]]; then
    for row in "${T_SIEVE_GRID[@]}"; do
        read -r bits trials _timeout <<< "$row"
        run_bench "scale-T-${bits}b" "T" "$bits" "$trials" 16 "KSIEVE_TSIEVE=1" || true
    done
    touch "$OUT_DIR/.done.scaling-T"
fi

# ---------- SECTION 4: K-sieve scaling (smaller range, reference only) ----------
banner "SECTION 4: K-sieve scaling (reference, 16-thread)"
if [[ ! -f "$OUT_DIR/.done.scaling-K" ]]; then
    for row in "${K_SIEVE_GRID[@]}"; do
        read -r bits trials _timeout <<< "$row"
        run_bench "scale-K-${bits}b" "K" "$bits" "$trials" 16 "" || true
    done
    touch "$OUT_DIR/.done.scaling-K"
fi

# ---------- SECTION 5: K/S/T/M paired comparison at PAIRED_BITS ----------
banner "SECTION 5: K/S/T/M paired @ ${PAIRED_BITS}-bit (16-thread, $PAIRED_TRIALS trials each)"
if [[ ! -f "$OUT_DIR/.done.paired" ]]; then
    # NOTE: these run on different RNG seed streams because each generator
    # (gen_semiprime, gen_semiprime_tform, gen_semiprime_newform) uses a
    # distinct residue-class filter; 30 trials each stabilizes the mean.
    # K-sieve draws old-form (-1,-1) semiprimes by default, same as S-sieve.
    run_bench "paired-K-${PAIRED_BITS}b" "K" "$PAIRED_BITS" "$PAIRED_TRIALS" 16 ""                  || true
    run_bench "paired-S-${PAIRED_BITS}b" "S" "$PAIRED_BITS" "$PAIRED_TRIALS" 16 "KSIEVE_SSIEVE=1"   || true
    run_bench "paired-T-${PAIRED_BITS}b" "T" "$PAIRED_BITS" "$PAIRED_TRIALS" 16 "KSIEVE_TSIEVE=1"   || true
    run_bench "paired-M-${PAIRED_BITS}b" "M" "$PAIRED_BITS" "$PAIRED_TRIALS" 16 "KSIEVE_MSIEVE=1"   || true
    touch "$OUT_DIR/.done.paired"
fi

# ---------- SECTION 6: single-thread vs 16-thread at THREADING_BITS ----------
banner "SECTION 6: threading scan — M-sieve ${THREADING_BITS}-bit, $THREADING_TRIALS trials"
if [[ ! -f "$OUT_DIR/.done.threads" ]]; then
    for t in "${THREAD_VARIANTS[@]}"; do
        run_bench "thread-M-${THREADING_BITS}b-t${t}" "M" "$THREADING_BITS" "$THREADING_TRIALS" "$t" "KSIEVE_MSIEVE=1" || true
    done
    touch "$OUT_DIR/.done.threads"
fi

# ---------- SECTION 7: BS_NP sweep at BSNP_BITS ----------
banner "SECTION 7: BS_NP sweep — ${BSNP_MODE} ${BSNP_BITS}-bit, $BSNP_TRIALS trials each"
if [[ ! -f "$OUT_DIR/.done.bsnp" ]]; then
    for np in 8 9 10 11 12 13 14; do
        run_bench "bsnp-${BSNP_MODE}-${BSNP_BITS}b-np${np}" "M" "$BSNP_BITS" "$BSNP_TRIALS" 16 "${BSNP_MODE}=1 KSIEVE_BS_NP=$np" || true
    done
    touch "$OUT_DIR/.done.bsnp"
fi

# ---------- build markdown tables from CSV ----------
banner "POST: building markdown tables from CSV"

python3 - "$CSV" > "$TABLES" <<'PYEOF'
import csv, sys, re
from collections import defaultdict

path = sys.argv[1]
rows = []
with open(path, 'r') as f:
    for r in csv.DictReader(f):
        rows.append(r)

def fnum(x, fmt="{:.3f}"):
    try: return fmt.format(float(x))
    except: return "—"

def irow(d):
    """Render one row with normalized columns."""
    try: bits = int(d['bits'])
    except: bits = d['bits']
    return (bits, d)

# Index by run_id prefix for each section.
by_prefix = defaultdict(list)
for r in rows:
    prefix = re.sub(r'-\d+b.*$', '', r['run_id']).rstrip('-')
    by_prefix[prefix].append(r)

out = []

out.append("# ksieve benchmark results\n")
out.append("Generated from `results.csv`. Paste sections into the paper as desired.\n")

# --- Section 1-3: scaling ---
def scaling_section(title, prefix_keys):
    buf = [f"\n## {title}\n"]
    buf.append("| Sieve | Bits | Trials OK/Req | Mean total (s) | Min (s) | Max (s) | Mean K0s | Mean CRT ops |")
    buf.append("|---|---|---|---|---|---|---|---|")
    entries = []
    for p in prefix_keys:
        for r in by_prefix.get(p, []):
            entries.append(r)
    entries.sort(key=lambda r: (r['sieve'], int(r['bits'])))
    for r in entries:
        ok = r['trials_ok']; req = r['trials_requested']
        buf.append(f"| {r['sieve']} | {r['bits']} | {ok}/{req} | {fnum(r['mean_total_s'])} | {fnum(r['min_s'])} | {fnum(r['max_s'])} | {r['mean_k0s']} | {r['mean_crt_ops']} |")
    return "\n".join(buf)

out.append(scaling_section("Scaling — M-sieve (16-thread)", ["scale-M"]))
out.append(scaling_section("Scaling — S-sieve (16-thread)", ["scale-S"]))
out.append(scaling_section("Scaling — T-sieve (16-thread)", ["scale-T"]))
out.append(scaling_section("Scaling — K-sieve (16-thread, reference)", ["scale-K"]))

# --- Combined S vs T scaling comparison (Rev. 8+) ---
# S-sieve and T-sieve should have essentially identical performance at matched
# bit size; this table makes the comparison easy to eyeball.
st_s = {r['bits']: r for r in rows if r['run_id'].startswith('scale-S-')}
st_t = {r['bits']: r for r in rows if r['run_id'].startswith('scale-T-')}
common_bits = sorted(set(st_s.keys()) & set(st_t.keys()), key=int)
if common_bits:
    out.append("\n## S-sieve vs T-sieve head-to-head (matched bit size, 16-thread)\n")
    out.append("| Bits | S mean (s) | T mean (s) | T/S ratio | S OK/Req | T OK/Req |")
    out.append("|---|---|---|---|---|---|")
    for b in common_bits:
        rs, rt = st_s[b], st_t[b]
        ratio = "—"
        try:
            ratio = f"{float(rt['mean_total_s'])/float(rs['mean_total_s']):.2f}×"
        except Exception:
            pass
        out.append(f"| {b} | {fnum(rs['mean_total_s'])} | {fnum(rt['mean_total_s'])} | {ratio} | "
                   f"{rs['trials_ok']}/{rs['trials_requested']} | "
                   f"{rt['trials_ok']}/{rt['trials_requested']} |")

# --- Section 5: paired ---
out.append(f"\n## Paired K / S / T / M comparison (same bit size, 16-thread)\n")
out.append("| Sieve | Bits | Trials | Mean total (s) | Mean K0s | Mean CRT ops | Ratio vs S |")
out.append("|---|---|---|---|---|---|---|")
paired = [r for r in rows if r['run_id'].startswith('paired-')]
# Present in canonical order: K (reference), S (-1,-1 old-form), T (+1,+1 old-form), M (new-form)
order = {'K': 0, 'S': 1, 'T': 2, 'M': 3}
paired.sort(key=lambda r: order.get(r['sieve'], 99))
s_mean = None
for r in paired:
    if r['sieve'] == 'S':
        try: s_mean = float(r['mean_total_s'])
        except: pass
for r in paired:
    ratio = "—"
    if s_mean:
        try: ratio = f"{float(r['mean_total_s'])/s_mean:.2f}×"
        except: pass
    out.append(f"| {r['sieve']} | {r['bits']} | {r['trials_ok']}/{r['trials_requested']} | {fnum(r['mean_total_s'])} | {r['mean_k0s']} | {r['mean_crt_ops']} | {ratio} |")

# --- Section 6: threading ---
out.append(f"\n## Threading scaling — M-sieve\n")
out.append("| Threads | Bits | Trials | Mean total (s) | Min (s) | Max (s) | Speedup vs 1t |")
out.append("|---|---|---|---|---|---|---|")
thread_rows = [r for r in rows if r['run_id'].startswith('thread-')]
thread_rows.sort(key=lambda r: int(r['threads']))
t1_mean = None
for r in thread_rows:
    if r['threads'] == '1':
        try: t1_mean = float(r['mean_total_s'])
        except: pass
for r in thread_rows:
    speedup = "—"
    if t1_mean:
        try: speedup = f"{t1_mean/float(r['mean_total_s']):.2f}×"
        except: pass
    out.append(f"| {r['threads']} | {r['bits']} | {r['trials_ok']}/{r['trials_requested']} | {fnum(r['mean_total_s'])} | {fnum(r['min_s'])} | {fnum(r['max_s'])} | {speedup} |")

# --- Section 7: BS_NP sweep ---
out.append(f"\n## BS_NP sweep (bitset prime count)\n")
out.append("| BS_NP | Bits | Trials | Mean total (s) | vs default |")
out.append("|---|---|---|---|---|")
bsnp_rows = [r for r in rows if r['run_id'].startswith('bsnp-')]
def extract_np(r):
    m = re.search(r'KSIEVE_BS_NP=(\d+)', r['extra_env'])
    return int(m.group(1)) if m else 0
bsnp_rows.sort(key=extract_np)
# Default is listed in code; 9 is historical default (crt_np+2 at 90-bit)
default_np = 9
default_mean = None
for r in bsnp_rows:
    if extract_np(r) == default_np:
        try: default_mean = float(r['mean_total_s'])
        except: pass
for r in bsnp_rows:
    np = extract_np(r)
    mark = ""
    rel = "—"
    if default_mean:
        try:
            ratio = default_mean / float(r['mean_total_s'])
            if ratio > 1.02: rel = f"{ratio:.2f}× faster"
            elif ratio < 0.98: rel = f"{1/ratio:.2f}× slower"
            else: rel = "≈ same"
        except: pass
    if np == default_np: mark = " (old default)"
    out.append(f"| {np}{mark} | {r['bits']} | {r['trials_ok']}/{r['trials_requested']} | {fnum(r['mean_total_s'])} | {rel} |")

print("\n".join(out))
PYEOF

banner "ALL DONE"
log "CSV     : $CSV"
log "Tables  : $TABLES"
log "Progress: $PROGRESS"
log "Per-trial logs in: $OUT_DIR/logs/"
echo
echo "To re-run only a specific section, delete the corresponding marker:"
echo "  rm $OUT_DIR/.done.scaling-M   # M-sieve scaling"
echo "  rm $OUT_DIR/.done.scaling-S   # S-sieve scaling"
echo "  rm $OUT_DIR/.done.scaling-T   # T-sieve scaling"
echo "  rm $OUT_DIR/.done.scaling-K   # K-sieve scaling"
echo "  rm $OUT_DIR/.done.paired      # K/S/T/M paired"
echo "  rm $OUT_DIR/.done.threads     # threading scan"
echo "  rm $OUT_DIR/.done.bsnp        # BS_NP sweep"
echo
echo "To re-run only specific bit sizes, delete the individual log:"
echo "  rm $OUT_DIR/logs/scale-M-100b.log"
echo "  rm $OUT_DIR/logs/scale-T-100b.log"
