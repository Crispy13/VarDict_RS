#!/usr/bin/env bash
set -euo pipefail

SCRIPT_NAME=$(basename "${BASH_SOURCE[0]}")
PROJECT_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$PROJECT_ROOT"

readonly REF_FASTA="testdata/hs37d5.fa"
readonly REF_FASTA_FAI="${REF_FASTA}.fai"
readonly BAM_PATH="testdata/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam"
readonly JAVA_JAR="VarDictJava/build/libs/VarDict-1.8.3.jar"
readonly RUST_BIN="target/debug-release/vardict"
readonly BASE_DIR="tmp/na12878_parity"
readonly DISK_ABORT_GB=5

JAVA_DIR=""
RUST_DIR=""
DIFF_DIR=""

CHR="20"
REQUESTED_CHR_LEN=""
CHR_LEN=""
SHARD_SIZE=1000000
MAX_PARALLEL=10
JAVA_MAX_PARALLEL=10
JAVA_HEAP="2g"
FREQ="0.01"
STOP_ON_FAIL=1
ALL_CHR=0
CLEAN_RUST=0
CLEANUP=1
EXTRA_OPTS=""
OPTS_LABEL="default"
SHARD_TIMEOUT=600
MEM_WARN_GB=4
MEM_ABORT_GB=2
DISK_WARN_GB=10

TOTAL_SHARDS=0
PASS_COUNT=0
FAIL_COUNT=0
EMPTY_COUNT=0
LAST_CHR_ELAPSED="0.0s"
LAST_CHR_ELAPSED_SECONDS="0.0"

GLOBAL_PASS=0
GLOBAL_FAIL=0
GLOBAL_EMPTY=0
GLOBAL_SHARDS=0
GLOBAL_CHR_PASS=0
GLOBAL_CHR_FAIL=0
declare -a GLOBAL_FAIL_DETAILS=()

declare -a SHARD_IDS=()
declare -a SHARD_STARTS=()
declare -a SHARD_ENDS=()
declare -a ACTIVE_PIDS=()
declare -a FAIL_MESSAGES=()
declare -a JSON_FAILURES=()

BATCH_FIRST_FAIL_IDX=""
FIFO_SEM=""
FIFO_OPEN=0
PROGRESS_REPORTER_PID=""
PHASE_PROGRESS_FILE=""
PHASE_ABORT_FLAG=""

usage() {
    cat <<EOF
Usage: ${SCRIPT_NAME} [options]

Options:
  --chr CHR            Chromosome name (default: 20)
  --all-chr            Run all chromosomes from the FAI index sequentially
  --chr-len LEN        Chromosome length; defaults to lookup in ${REF_FASTA_FAI}
  --shard-size N       Override shard size in bases (default: 1000000)
  --parallel N         Override worker count for Java/Rust phases (default: 10)
  --java-parallel N    (Ignored — kept for backwards compatibility)
  --java-heap SIZE     Java heap size per worker (default: 2g)
  --freq F             Override VarDict frequency threshold (default: 0.01)
  --no-stop            Continue past mismatches instead of stopping after the first failed batch
  --clean-rust         Remove only Rust + diff outputs (preserves Java cache)
  --no-cleanup         Keep all output files after comparison (default: cleanup passing outputs)
  --opts "FLAGS"       Extra CLI flags passed to both Java and Rust (default: none)
  --opts-label LBL     Cache subdirectory label for this option set (default: "default")
  --timeout SECONDS    Per-shard timeout for Java and Rust commands (default: 600)
  --mem-warn-gb N      Warn if MemAvailable drops below N GiB (default: 4)
  --mem-abort-gb N     Abort if MemAvailable drops below N GiB (default: 2)
  --disk-warn-gb N     Warn if available disk drops below N GiB (default: 10)
  --help               Show this help

Execution model:
  Phase 1: Generate Java outputs for needed shards only (FIFO worker pool)
  Phase 2: Generate Rust outputs for needed shards only (FIFO worker pool)
  Phase 3: Compare outputs (sequential, lightweight)
  Phase 4: Cleanup passing outputs to save disk (unless --no-cleanup)
  Java and Rust never run simultaneously, preventing memory exhaustion.

Caching:
  Java outputs are cached and reused across runs (deterministic).
  Rust outputs are automatically regenerated when the binary is newer.
  Passing shards receive shard-level .verified markers and are skipped on resume.
  After all shards pass, a chromosome-level .verified marker is created and outputs are deleted.
  Re-runs skip verified configs unless the Rust binary has been rebuilt.
EOF
}

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

require_file() {
    local path=$1

    [[ -f "$path" ]] || fail "Required file not found: $path"
}

format_elapsed() {
    local start_ns=$1
    local end_ns=$2

    awk -v start="$start_ns" -v end="$end_ns" 'BEGIN { printf "%.1fs", (end - start) / 1000000000 }'
}

format_elapsed_seconds() {
    local start_ns=$1
    local end_ns=$2

    awk -v start="$start_ns" -v end="$end_ns" 'BEGIN { printf "%.1f", (end - start) / 1000000000 }'
}

remove_stale_diff() {
    local diff_file=$1

    if [[ -f "$diff_file" ]]; then
        rm -f "$diff_file"
    fi
}

stop_progress_reporter() {
    if [[ -n "$PROGRESS_REPORTER_PID" ]] && kill -0 "$PROGRESS_REPORTER_PID" 2>/dev/null; then
        kill "$PROGRESS_REPORTER_PID" 2>/dev/null || true
        wait "$PROGRESS_REPORTER_PID" 2>/dev/null || true
    fi
    PROGRESS_REPORTER_PID=""
}

close_fifo_semaphore() {
    if (( FIFO_OPEN )); then
        exec 3>&-
        FIFO_OPEN=0
    fi
    if [[ -n "$FIFO_SEM" ]]; then
        rm -f "$FIFO_SEM"
        FIFO_SEM=""
    fi
}

cleanup_temp_files() {
    if [[ -n "$PHASE_PROGRESS_FILE" ]]; then
        rm -f "$PHASE_PROGRESS_FILE"
        PHASE_PROGRESS_FILE=""
    fi
    if [[ -n "$PHASE_ABORT_FLAG" ]]; then
        rm -f "$PHASE_ABORT_FLAG"
        PHASE_ABORT_FLAG=""
    fi
}

kill_active_jobs() {
    local pid
    local child_pid

    for pid in "${ACTIVE_PIDS[@]}"; do
        if [[ -n "$pid" ]] && kill -0 "$pid" 2>/dev/null; then
            for child_pid in $(pgrep -P "$pid" 2>/dev/null); do
                kill "$child_pid" 2>/dev/null || true
            done
            kill "$pid" 2>/dev/null || true
        fi
    done

    for pid in "${ACTIVE_PIDS[@]}"; do
        if [[ -n "$pid" ]]; then
            wait "$pid" 2>/dev/null || true
        fi
    done

    ACTIVE_PIDS=()
}

cleanup_script() {
    stop_progress_reporter
    close_fifo_semaphore
    kill_active_jobs
    cleanup_temp_files
}

trap cleanup_script EXIT

summary_chr_label() {
    if [[ "$CHR" == chr* ]]; then
        printf '%s' "$CHR"
    else
        printf 'chr%s' "$CHR"
    fi
}

display_region() {
    local start=$1
    local end=$2

    printf '%s:%s-%s' "$(summary_chr_label)" "$start" "$end"
}

resolve_chr_length() {
    if [[ -n "$REQUESTED_CHR_LEN" ]]; then
        CHR_LEN=$REQUESTED_CHR_LEN
        return 0
    fi

    require_file "$REF_FASTA_FAI"
    CHR_LEN=$(awk -F '\t' -v chr="$CHR" '$1 == chr { print $2; exit }' "$REF_FASTA_FAI")
    [[ -n "$CHR_LEN" ]] || fail "Chromosome '$CHR' not found in $REF_FASTA_FAI"
}

validate_positive_integer() {
    local label=$1
    local value=$2

    [[ "$value" =~ ^[0-9]+$ ]] || fail "$label must be a positive integer: $value"
    (( value > 0 )) || fail "$label must be greater than zero: $value"
}

generate_shards() {
    local start=1
    local end=0

    SHARD_IDS=()
    SHARD_STARTS=()
    SHARD_ENDS=()
    TOTAL_SHARDS=0

    while (( start <= CHR_LEN )); do
        TOTAL_SHARDS=$((TOTAL_SHARDS + 1))
        end=$((start + SHARD_SIZE - 1))
        if (( end > CHR_LEN )); then
            end=$CHR_LEN
        fi
        SHARD_IDS+=("$(printf '%03d' "$TOTAL_SHARDS")")
        SHARD_STARTS+=("$start")
        SHARD_ENDS+=("$end")
        start=$((end + 1))
    done
}

prepare_dirs() {
    JAVA_DIR="${BASE_DIR}/${OPTS_LABEL}/${CHR}/java"
    RUST_DIR="${BASE_DIR}/${OPTS_LABEL}/${CHR}/rust"
    DIFF_DIR="${BASE_DIR}/${OPTS_LABEL}/${CHR}/diff"

    if (( CLEAN_RUST )); then
        rm -rf "$RUST_DIR" "$DIFF_DIR"
    fi

    mkdir -p "$JAVA_DIR" "$RUST_DIR" "$DIFF_DIR"
}

count_file_lines() {
    local file=$1

    if [[ -f "$file" ]]; then
        wc -l < "$file" | tr -d '[:space:]'
    else
        printf '0'
    fi
}

find_first_difference_line() {
    local left_file=$1
    local right_file=$2

    paste -d$'\x01' "$left_file" "$right_file" | awk -F$'\x01' '{ if ($1 != $2) { print NR; exit } }'
}

write_failure_meta() {
    local meta_file=$1
    local region=$2
    local shard_num=$3
    local java_file=$4
    local rust_file=$5
    local diff_file=$6
    local java_log=$7
    local rust_log=$8
    local java_lines=$9
    local rust_lines=${10}
    local first_diff_line=${11}
    local summary=${12}
    local reason=${13}

    {
        printf 'REGION=%q\n' "$region"
        printf 'SHARD=%q\n' "$shard_num"
        printf 'JAVA_FILE=%q\n' "$java_file"
        printf 'RUST_FILE=%q\n' "$rust_file"
        printf 'DIFF_FILE=%q\n' "$diff_file"
        printf 'JAVA_LOG=%q\n' "$java_log"
        printf 'RUST_LOG=%q\n' "$rust_log"
        printf 'JAVA_LINES=%q\n' "$java_lines"
        printf 'RUST_LINES=%q\n' "$rust_lines"
        printf 'FIRST_DIFF_LINE=%q\n' "$first_diff_line"
        printf 'SUMMARY=%q\n' "$summary"
        printf 'REASON=%q\n' "$reason"
    } >"$meta_file"
}

status_file_for_idx() {
    local idx=$1

    printf '%s/shard_%s.status' "$DIFF_DIR" "${SHARD_IDS[$idx]}"
}

meta_file_for_idx() {
    local idx=$1

    printf '%s/shard_%s.meta' "$DIFF_DIR" "${SHARD_IDS[$idx]}"
}

diff_file_for_idx() {
    local idx=$1

    printf '%s/shard_%s.diff' "$DIFF_DIR" "${SHARD_IDS[$idx]}"
}

failure_summary_for_idx() {
    local idx=$1
    local meta_file
    local SUMMARY=""

    meta_file=$(meta_file_for_idx "$idx")
    if [[ -f "$meta_file" ]]; then
        # shellcheck disable=SC1090
        source "$meta_file"
        if [[ -n "$SUMMARY" ]]; then
            printf '%s' "$SUMMARY"
            return 0
        fi
    fi

    printf 'mismatch'
}

resource_log_file() {
    printf '%s/resource_log.tsv' "$DIFF_DIR"
}

ensure_resource_log() {
    local log_file

    log_file=$(resource_log_file)
    if [[ ! -f "$log_file" ]]; then
        printf 'timestamp\tphase\tstage\tmem_available_kb\tmem_available_gb\tdisk_available_kb\tdisk_available_gb\tstatus\tmessage\n' >"$log_file"
    fi
}

mem_available_kb() {
    awk '/MemAvailable:/ { print $2; exit }' /proc/meminfo
}

disk_available_kb() {
    df --output=avail . | awk 'NR == 2 { print $1 }'
}

kb_to_gb() {
    local kb=$1

    awk -v value="$kb" 'BEGIN { printf "%.1f", value / 1048576 }'
}

log_resource_snapshot() {
    local phase=$1
    local stage=$2
    local status=$3
    local message=$4
    local log_file
    local mem_kb
    local disk_kb
    local safe_message

    log_file=$(resource_log_file)
    mem_kb=$(mem_available_kb)
    disk_kb=$(disk_available_kb)
    safe_message=${message//$'\t'/ }
    safe_message=${safe_message//$'\n'/ }

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$(date -Iseconds)" \
        "$phase" \
        "$stage" \
        "$mem_kb" \
        "$(kb_to_gb "$mem_kb")" \
        "$disk_kb" \
        "$(kb_to_gb "$disk_kb")" \
        "$status" \
        "$safe_message" >>"$log_file"
}

check_preflight_resources() {
    local phase_label=$1
    local mem_kb
    local disk_kb
    local mem_warn_kb=$((MEM_WARN_GB * 1024 * 1024))
    local mem_abort_kb=$((MEM_ABORT_GB * 1024 * 1024))
    local disk_warn_kb=$((DISK_WARN_GB * 1024 * 1024))
    local disk_abort_kb=$((DISK_ABORT_GB * 1024 * 1024))
    local messages=()
    local warn=0

    mem_kb=$(mem_available_kb)
    disk_kb=$(disk_available_kb)

    if (( mem_kb < mem_abort_kb )); then
        log_resource_snapshot "$phase_label" "preflight" "abort" "memory below abort threshold"
        fail "${phase_label}: MemAvailable $(kb_to_gb "$mem_kb")GB is below abort threshold ${MEM_ABORT_GB}GB"
    fi
    if (( disk_kb < disk_abort_kb )); then
        log_resource_snapshot "$phase_label" "preflight" "abort" "disk below abort threshold"
        fail "${phase_label}: available disk $(kb_to_gb "$disk_kb")GB is below abort threshold ${DISK_ABORT_GB}GB"
    fi

    if (( mem_kb < mem_warn_kb )); then
        messages+=("MemAvailable $(kb_to_gb "$mem_kb")GB < ${MEM_WARN_GB}GB")
        warn=1
    fi
    if (( disk_kb < disk_warn_kb )); then
        messages+=("disk $(kb_to_gb "$disk_kb")GB < ${DISK_WARN_GB}GB")
        warn=1
    fi

    if (( warn )); then
        log_resource_snapshot "$phase_label" "preflight" "warn" "${messages[*]}"
        echo "WARNING: ${phase_label}: ${messages[*]}" >&2
    else
        log_resource_snapshot "$phase_label" "preflight" "ok" "resources within thresholds"
    fi
}

monitor_inflight_resources() {
    local phase_label=$1
    local abort_flag=$2
    local mem_kb
    local mem_abort_kb=$((MEM_ABORT_GB * 1024 * 1024))

    mem_kb=$(mem_available_kb)
    if (( mem_kb < mem_abort_kb )); then
        printf 'PHASE=%s\nMEM_AVAILABLE_KB=%s\n' "$phase_label" "$mem_kb" >"$abort_flag"
        log_resource_snapshot "$phase_label" "monitor" "abort" "memory below abort threshold"
        echo "WARNING: ${phase_label}: stopping new shard launches because MemAvailable is $(kb_to_gb "$mem_kb")GB (< ${MEM_ABORT_GB}GB)" >&2
        return 1
    fi

    log_resource_snapshot "$phase_label" "monitor" "ok" "memory within threshold"
    return 0
}

print_phase_progress_line() {
    local phase_label=$1
    local total=$2
    local progress_file=$3
    local start_epoch=$4
    local completed
    local elapsed
    local eta="0.0s"

    completed=$(count_file_lines "$progress_file")
    if (( completed > total )); then
        completed=$total
    fi

    elapsed=$(awk -v start="$start_epoch" -v now="$(date +%s)" 'BEGIN { printf "%.1f", now - start }')
    if (( completed > 0 && completed < total )); then
        eta=$(awk -v elapsed="$elapsed" -v completed="$completed" -v total="$total" 'BEGIN { printf "%.1fs", (elapsed / completed) * (total - completed) }')
    elif (( completed == 0 )); then
        eta="n/a"
    fi

    printf '  %s: [%d/%d shards] %ss elapsed, ~%s ETA\n' "$phase_label" "$completed" "$total" "$elapsed" "$eta" >&2
}

progress_reporter_loop() {
    local phase_label=$1
    local total=$2
    local progress_file=$3
    local start_epoch=$4
    local abort_flag=$5
    local tick=0

    while true; do
        sleep 5
        print_phase_progress_line "$phase_label" "$total" "$progress_file" "$start_epoch"
        tick=$((tick + 1))
        if (( tick % 3 == 0 )); then
            monitor_inflight_resources "$phase_label" "$abort_flag" || true
        fi
    done
}

start_progress_reporter() {
    local phase_label=$1
    local total=$2
    local progress_file=$3
    local abort_flag=$4
    local start_epoch=$5

    if (( total == 0 )); then
        return 0
    fi

    progress_reporter_loop "$phase_label" "$total" "$progress_file" "$start_epoch" "$abort_flag" &
    PROGRESS_REPORTER_PID=$!
}

shard_verified_file_for_idx() {
    local idx=$1

    printf '%s/shard_%s.verified' "$DIFF_DIR" "${SHARD_IDS[$idx]}"
}

is_shard_verified() {
    local idx=$1
    local verified_file
    local stored_mtime
    local current_mtime

    verified_file=$(shard_verified_file_for_idx "$idx")
    if [[ ! -f "$verified_file" ]]; then
        return 1
    fi

    stored_mtime=$(grep '^RUST_BINARY_MTIME=' "$verified_file" | cut -d= -f2)
    current_mtime=$(stat -c %Y "$RUST_BIN" 2>/dev/null) || return 1
    [[ "$stored_mtime" == "$current_mtime" ]]
}

write_shard_verified() {
    local idx=$1
    local verified_file
    local rust_mtime

    verified_file=$(shard_verified_file_for_idx "$idx")
    rust_mtime=$(stat -c %Y "$RUST_BIN" 2>/dev/null || echo unknown)
    {
        printf 'RUST_BINARY_MTIME=%s\n' "$rust_mtime"
        printf 'VERIFIED_AT=%s\n' "$(date -Iseconds)"
    } >"$verified_file"
}

remove_shard_verified() {
    local idx=$1
    local verified_file

    verified_file=$(shard_verified_file_for_idx "$idx")
    rm -f "$verified_file"
}

phase_state_file_for_idx() {
    local mode=$1
    local idx=$2

    printf '%s/shard_%s.%s_phase' "$DIFF_DIR" "${SHARD_IDS[$idx]}" "$mode"
}

record_phase_failure() {
    local mode=$1
    local idx=$2
    local exit_code=$3
    local reason=$4
    local log_file=$5
    local region=$6
    local state_file

    state_file=$(phase_state_file_for_idx "$mode" "$idx")
    {
        printf 'MODE=%q\n' "$mode"
        printf 'EXIT_CODE=%q\n' "$exit_code"
        printf 'REASON=%q\n' "$reason"
        printf 'LOG_FILE=%q\n' "$log_file"
        printf 'REGION=%q\n' "$region"
    } >"$state_file"
}

clear_phase_failure() {
    local mode=$1
    local idx=$2
    local state_file

    state_file=$(phase_state_file_for_idx "$mode" "$idx")
    rm -f "$state_file"
}

read_phase_failure() {
    local mode=$1
    local idx=$2
    local state_file

    state_file=$(phase_state_file_for_idx "$mode" "$idx")
    if [[ ! -f "$state_file" ]]; then
        return 1
    fi

    # shellcheck disable=SC1090
    source "$state_file"
    return 0
}

needs_missing_outputs() {
    local mode=$1
    local -n missing_indices_ref=$2
    local idx
    local out_file

    missing_indices_ref=()

    for ((idx = 0; idx < TOTAL_SHARDS; idx++)); do
        if is_shard_verified "$idx"; then
            continue
        fi

        if [[ "$mode" == "java" ]]; then
            out_file="${JAVA_DIR}/shard_${SHARD_IDS[$idx]}.tsv"
            if [[ ! -f "$out_file" ]]; then
                missing_indices_ref+=("$idx")
            fi
            continue
        fi

        out_file="${RUST_DIR}/shard_${SHARD_IDS[$idx]}.tsv"
        if [[ ! -f "$out_file" || "$RUST_BIN" -nt "$out_file" ]]; then
            missing_indices_ref+=("$idx")
        fi
    done

    (( ${#missing_indices_ref[@]} > 0 ))
}

generate_java_shard() {
    local idx=$1
    local shard_num=${SHARD_IDS[$idx]}
    local start=${SHARD_STARTS[$idx]}
    local end=${SHARD_ENDS[$idx]}
    local region="${CHR}:${start}-${end}"
    local java_file="${JAVA_DIR}/shard_${shard_num}.tsv"
    local java_log="${JAVA_DIR}/shard_${shard_num}.log"
    local java_rc=0
    local freq_flag="-f $FREQ"

    if [[ "$EXTRA_OPTS" =~ (^|[[:space:]])-f[[:space:]] ]]; then
        freq_flag=""
    fi

    if [[ -f "$java_file" ]]; then
        clear_phase_failure "java" "$idx"
        return 0
    fi

    # shellcheck disable=SC2086
    timeout "$SHARD_TIMEOUT" java "-Xmx${JAVA_HEAP}" -classpath "$JAVA_JAR" com.astrazeneca.vardict.Main \
        -G "$REF_FASTA" \
        -b "$BAM_PATH" \
        -N NA12878 $freq_flag -th 1 \
        $EXTRA_OPTS \
        -R "$region" \
        >"$java_file" 2>"$java_log" || java_rc=$?

    if (( java_rc != 0 )); then
        rm -f "$java_file"
        if (( java_rc == 124 )); then
            record_phase_failure "java" "$idx" "$java_rc" "timeout" "$java_log" "$region"
        else
            record_phase_failure "java" "$idx" "$java_rc" "java_command_failed" "$java_log" "$region"
        fi
        echo "JAVA_FAIL shard_${shard_num} (${region}): exit ${java_rc}" >&2
        return 1
    fi

    clear_phase_failure "java" "$idx"
}

generate_rust_shard() {
    local idx=$1
    local shard_num=${SHARD_IDS[$idx]}
    local start=${SHARD_STARTS[$idx]}
    local end=${SHARD_ENDS[$idx]}
    local region="${CHR}:${start}-${end}"
    local rust_file="${RUST_DIR}/shard_${shard_num}.tsv"
    local rust_log="${RUST_DIR}/shard_${shard_num}.log"
    local rust_rc=0
    local freq_flag="-f $FREQ"

    if [[ "$EXTRA_OPTS" =~ (^|[[:space:]])-f[[:space:]] ]]; then
        freq_flag=""
    fi

    if [[ -f "$rust_file" ]] && ! [[ "$RUST_BIN" -nt "$rust_file" ]]; then
        clear_phase_failure "rust" "$idx"
        return 0
    fi

    rm -f "$rust_file"
    # shellcheck disable=SC2086
    timeout "$SHARD_TIMEOUT" "$RUST_BIN" \
        -G "$REF_FASTA" \
        -b "$BAM_PATH" \
        -N NA12878 $freq_flag \
        $EXTRA_OPTS \
        -R "$region" \
        >"$rust_file" 2>"$rust_log" || rust_rc=$?

    if (( rust_rc != 0 )); then
        rm -f "$rust_file"
        if (( rust_rc == 124 )); then
            record_phase_failure "rust" "$idx" "$rust_rc" "timeout" "$rust_log" "$region"
        else
            record_phase_failure "rust" "$idx" "$rust_rc" "rust_command_failed" "$rust_log" "$region"
        fi
        echo "RUST_FAIL shard_${shard_num} (${region}): exit ${rust_rc}" >&2
        return 1
    fi

    clear_phase_failure "rust" "$idx"
}

write_truncated_diff() {
    local java_file=$1
    local rust_file=$2
    local diff_file=$3

    ( diff -u "$java_file" "$rust_file" || true ) | head -n 500 >"$diff_file"
}

compare_shard() {
    local idx=$1
    local shard_num=${SHARD_IDS[$idx]}
    local start=${SHARD_STARTS[$idx]}
    local end=${SHARD_ENDS[$idx]}
    local region="${CHR}:${start}-${end}"
    local java_file="${JAVA_DIR}/shard_${shard_num}.tsv"
    local rust_file="${RUST_DIR}/shard_${shard_num}.tsv"
    local java_log="${JAVA_DIR}/shard_${shard_num}.log"
    local rust_log="${RUST_DIR}/shard_${shard_num}.log"
    local status_file
    local meta_file
    local diff_file
    local java_lines=0
    local rust_lines=0
    local first_diff_line=""
    local summary=""
    local reason=""
    local MODE=""
    local EXIT_CODE=""
    local REASON=""
    local LOG_FILE=""

    status_file=$(status_file_for_idx "$idx")
    meta_file=$(meta_file_for_idx "$idx")
    diff_file=$(diff_file_for_idx "$idx")

    rm -f "$status_file" "$meta_file" "$diff_file"

    if is_shard_verified "$idx"; then
        clear_phase_failure "java" "$idx"
        clear_phase_failure "rust" "$idx"
        printf 'PASS\n' >"$status_file"
        remove_stale_diff "$diff_file"
        return 0
    fi

    if [[ ! -f "$java_file" || ! -f "$rust_file" ]]; then
        summary="missing output"
        reason="missing_output"

        if [[ ! -f "$java_file" ]] && read_phase_failure "java" "$idx"; then
            if [[ "$REASON" == "timeout" ]]; then
                summary="java timeout (${SHARD_TIMEOUT}s)"
                reason="timeout"
            else
                summary="java command failed (exit ${EXIT_CODE})"
                reason="$REASON"
            fi
        fi
        if [[ ! -f "$rust_file" ]] && read_phase_failure "rust" "$idx"; then
            if [[ "$REASON" == "timeout" ]]; then
                summary="rust timeout (${SHARD_TIMEOUT}s)"
                reason="timeout"
            else
                summary="rust command failed (exit ${EXIT_CODE})"
                reason="$REASON"
            fi
        fi

        printf 'FAIL\n' >"$status_file"
        printf '%s\n' "$summary" >"$diff_file"
        remove_shard_verified "$idx"
        write_failure_meta "$meta_file" "$region" "$shard_num" "$java_file" "$rust_file" "$diff_file" "$java_log" "$rust_log" "0" "0" "" "$summary" "$reason"
        return 0
    fi

    if [[ ! -s "$java_file" && ! -s "$rust_file" ]]; then
        clear_phase_failure "java" "$idx"
        clear_phase_failure "rust" "$idx"
        printf 'EMPTY\n' >"$status_file"
        write_shard_verified "$idx"
        return 0
    fi

    if [[ ! -s "$java_file" && -s "$rust_file" ]]; then
        rust_lines=$(count_file_lines "$rust_file")
        printf 'FAIL\n' >"$status_file"
        summary="SUSPECT: java output empty but rust has ${rust_lines} lines (likely transient Java failure — delete java cache and re-run)"
        remove_shard_verified "$idx"
        write_failure_meta "$meta_file" "$region" "$shard_num" "$java_file" "$rust_file" "$diff_file" "$java_log" "$rust_log" "0" "$rust_lines" "" "$summary" "java_empty_suspect"
        return 0
    fi
    if [[ -s "$java_file" && ! -s "$rust_file" ]]; then
        java_lines=$(count_file_lines "$java_file")
        printf 'FAIL\n' >"$status_file"
        summary="SUSPECT: rust output empty but java has ${java_lines} lines (likely transient Rust failure — delete rust cache and re-run)"
        remove_shard_verified "$idx"
        write_failure_meta "$meta_file" "$region" "$shard_num" "$java_file" "$rust_file" "$diff_file" "$java_log" "$rust_log" "$java_lines" "0" "" "$summary" "rust_empty_suspect"
        return 0
    fi

    if cmp -s "$java_file" "$rust_file"; then
        clear_phase_failure "java" "$idx"
        clear_phase_failure "rust" "$idx"
        printf 'PASS\n' >"$status_file"
        remove_stale_diff "$diff_file"
        write_shard_verified "$idx"
        return 0
    fi

    java_lines=$(count_file_lines "$java_file")
    rust_lines=$(count_file_lines "$rust_file")
    first_diff_line=$(find_first_difference_line "$java_file" "$rust_file")
    write_truncated_diff "$java_file" "$rust_file" "$diff_file"

    if [[ "$java_lines" != "$rust_lines" ]]; then
        summary="line_count_mismatch (java=${java_lines}, rust=${rust_lines})"
        reason="line_count_mismatch"
    elif [[ -n "$first_diff_line" ]]; then
        summary="first difference at line ${first_diff_line}"
        reason="content_mismatch"
    else
        summary="byte_mismatch"
        reason="byte_mismatch"
    fi

    printf 'FAIL\n' >"$status_file"
    remove_shard_verified "$idx"
    write_failure_meta "$meta_file" "$region" "$shard_num" "$java_file" "$rust_file" "$diff_file" "$java_log" "$rust_log" "$java_lines" "$rust_lines" "$first_diff_line" "$summary" "$reason"
}

# Legacy process_shard removed in v2 — phases run separately via run_parallel_phase

wait_for_batch_jobs() {
    local pid
    local wait_status=0
    local status=0

    for pid in "${ACTIVE_PIDS[@]}"; do
        wait_status=0
        wait "$pid" || wait_status=$?
        if (( wait_status != 0 )) && (( status == 0 )); then
            status=$wait_status
        fi
    done

    ACTIVE_PIDS=()
    return "$status"
}

append_json_failure() {
    local idx=$1
    local meta_file
    local REGION=""
    local REASON=""
    local FIRST_DIFF_LINE=""

    meta_file=$(meta_file_for_idx "$idx")
    if [[ -f "$meta_file" ]]; then
        # shellcheck disable=SC1090
        source "$meta_file"
    fi

    if [[ -z "$REGION" ]]; then
        REGION=$(display_region "${SHARD_STARTS[$idx]}" "${SHARD_ENDS[$idx]}")
    fi
    if [[ -z "$REASON" ]]; then
        REASON="unknown"
    fi

    JSON_FAILURES+=("${SHARD_IDS[$idx]}|${REGION}|${REASON}|${FIRST_DIFF_LINE}")
}

collect_batch_results() {
    local batch_start=$1
    local batch_end=$2
    local idx
    local shard_num
    local start
    local end
    local status_file
    local status_value=""
    local summary=""

    BATCH_FIRST_FAIL_IDX=""
    JSON_FAILURES=()

    for ((idx = batch_start; idx <= batch_end; idx++)); do
        shard_num=${SHARD_IDS[$idx]}
        start=${SHARD_STARTS[$idx]}
        end=${SHARD_ENDS[$idx]}
        status_file=$(status_file_for_idx "$idx")

        if [[ ! -f "$status_file" ]]; then
            FAIL_COUNT=$((FAIL_COUNT + 1))
            FAIL_MESSAGES+=("  shard_${shard_num} ($(display_region "$start" "$end")): missing status")
            append_json_failure "$idx"
            if [[ -z "$BATCH_FIRST_FAIL_IDX" ]]; then
                BATCH_FIRST_FAIL_IDX=$idx
            fi
            continue
        fi

        status_value=$(<"$status_file")
        case "$status_value" in
            PASS)
                PASS_COUNT=$((PASS_COUNT + 1))
                ;;
            EMPTY)
                PASS_COUNT=$((PASS_COUNT + 1))
                EMPTY_COUNT=$((EMPTY_COUNT + 1))
                ;;
            FAIL)
                FAIL_COUNT=$((FAIL_COUNT + 1))
                summary=$(failure_summary_for_idx "$idx")
                FAIL_MESSAGES+=("  shard_${shard_num} ($(display_region "$start" "$end")): ${summary}")
                append_json_failure "$idx"
                if [[ -z "$BATCH_FIRST_FAIL_IDX" ]]; then
                    BATCH_FIRST_FAIL_IDX=$idx
                fi
                ;;
            *)
                FAIL_COUNT=$((FAIL_COUNT + 1))
                FAIL_MESSAGES+=("  shard_${shard_num} ($(display_region "$start" "$end")): unknown status '${status_value}'")
                append_json_failure "$idx"
                if [[ -z "$BATCH_FIRST_FAIL_IDX" ]]; then
                    BATCH_FIRST_FAIL_IDX=$idx
                fi
                ;;
        esac
    done
}

print_failure_details() {
    local idx=$1
    local shard_num=${SHARD_IDS[$idx]}
    local start=${SHARD_STARTS[$idx]}
    local end=${SHARD_ENDS[$idx]}
    local meta_file
    local REGION=""
    local SHARD=""
    local JAVA_FILE="${JAVA_DIR}/shard_${shard_num}.tsv"
    local RUST_FILE="${RUST_DIR}/shard_${shard_num}.tsv"
    local DIFF_FILE="${DIFF_DIR}/shard_${shard_num}.diff"
    local JAVA_LOG="${JAVA_DIR}/shard_${shard_num}.log"
    local RUST_LOG="${RUST_DIR}/shard_${shard_num}.log"
    local JAVA_LINES=""
    local RUST_LINES=""
    local FIRST_DIFF_LINE=""
    local SUMMARY=""
    local REASON=""
    local java_line_label="missing"
    local rust_line_label="missing"

    meta_file=$(meta_file_for_idx "$idx")
    if [[ -f "$meta_file" ]]; then
        # shellcheck disable=SC1090
        source "$meta_file"
    fi

    if [[ -f "$JAVA_FILE" ]]; then
        if [[ -z "$JAVA_LINES" || "$JAVA_LINES" == "0" ]]; then
            JAVA_LINES=$(count_file_lines "$JAVA_FILE")
        fi
        java_line_label="${JAVA_LINES} lines"
    fi

    if [[ -f "$RUST_FILE" ]]; then
        if [[ -z "$RUST_LINES" || "$RUST_LINES" == "0" ]]; then
            RUST_LINES=$(count_file_lines "$RUST_FILE")
        fi
        rust_line_label="${RUST_LINES} lines"
    fi

    echo
    echo "MISMATCH: $(display_region "$start" "$end") (shard ${shard_num})"
    echo "  Java: ${JAVA_FILE} (${java_line_label})"
    echo "  Rust: ${RUST_FILE} (${rust_line_label})"
    echo "  Diff: ${DIFF_FILE}"
    if [[ -n "$FIRST_DIFF_LINE" ]]; then
        echo "  First difference at line ${FIRST_DIFF_LINE}"
    fi
    if [[ -n "$SUMMARY" && "$SUMMARY" != "first difference at line ${FIRST_DIFF_LINE}" ]]; then
        echo "  Reason: ${SUMMARY}"
    fi
    if [[ "$REASON" == "java_command_failed" || "$REASON" == "timeout" ]]; then
        echo "  Java log: ${JAVA_LOG}"
        echo "  Rust log: ${RUST_LOG}"
    fi
    if [[ "$REASON" == "rust_command_failed" ]]; then
        echo "  Rust log: ${RUST_LOG}"
    fi
}

print_summary() {
    echo
    echo "=== NA12878 Parity Summary ($(summary_chr_label)) ==="
    echo "Shards: ${TOTAL_SHARDS} (pass: ${PASS_COUNT}, fail: ${FAIL_COUNT}, empty: ${EMPTY_COUNT})"
    echo "Time: ${LAST_CHR_ELAPSED}"

    if (( FAIL_COUNT > 0 )); then
        echo
        echo "Failing shards:"
        printf '%s\n' "${FAIL_MESSAGES[@]}"
    fi
}

print_global_summary() {
    local total_chr=$((GLOBAL_CHR_PASS + GLOBAL_CHR_FAIL))

    echo
    echo "==============================="
    echo "=== ALL CHROMOSOMES SUMMARY ==="
    echo "==============================="
    echo "Chromosomes: ${total_chr} (pass: ${GLOBAL_CHR_PASS}, fail: ${GLOBAL_CHR_FAIL})"
    echo "Shards:      ${GLOBAL_SHARDS} (pass: ${GLOBAL_PASS}, fail: ${GLOBAL_FAIL}, empty: ${GLOBAL_EMPTY})"
    if (( GLOBAL_CHR_FAIL > 0 )); then
        echo
        echo "Failing chromosomes:"
        printf '%s\n' "${GLOBAL_FAIL_DETAILS[@]}"
    fi
}

parse_args() {
    while (( $# > 0 )); do
        case "$1" in
            --chr)
                shift
                [[ $# -gt 0 ]] || fail "Missing value for --chr"
                CHR=$1
                ;;
            --all-chr)
                ALL_CHR=1
                ;;
            --chr-len)
                shift
                [[ $# -gt 0 ]] || fail "Missing value for --chr-len"
                REQUESTED_CHR_LEN=$1
                ;;
            --shard-size)
                shift
                [[ $# -gt 0 ]] || fail "Missing value for --shard-size"
                SHARD_SIZE=$1
                ;;
            --parallel)
                shift
                [[ $# -gt 0 ]] || fail "Missing value for --parallel"
                MAX_PARALLEL=$1
                ;;
            --java-parallel)
                shift
                [[ $# -gt 0 ]] || fail "Missing value for --java-parallel"
                JAVA_MAX_PARALLEL=$1
                ;;
            --java-heap)
                shift
                [[ $# -gt 0 ]] || fail "Missing value for --java-heap"
                JAVA_HEAP=$1
                ;;
            --freq)
                shift
                [[ $# -gt 0 ]] || fail "Missing value for --freq"
                FREQ=$1
                ;;
            --no-stop)
                STOP_ON_FAIL=0
                ;;
            --clean-rust)
                CLEAN_RUST=1
                ;;
            --no-cleanup)
                CLEANUP=0
                ;;
            --opts)
                shift
                [[ $# -gt 0 ]] || fail "Missing value for --opts"
                EXTRA_OPTS="$1"
                ;;
            --opts-label)
                shift
                [[ $# -gt 0 ]] || fail "Missing value for --opts-label"
                OPTS_LABEL="$1"
                ;;
            --timeout)
                shift
                [[ $# -gt 0 ]] || fail "Missing value for --timeout"
                SHARD_TIMEOUT=$1
                ;;
            --mem-warn-gb)
                shift
                [[ $# -gt 0 ]] || fail "Missing value for --mem-warn-gb"
                MEM_WARN_GB=$1
                ;;
            --mem-abort-gb)
                shift
                [[ $# -gt 0 ]] || fail "Missing value for --mem-abort-gb"
                MEM_ABORT_GB=$1
                ;;
            --disk-warn-gb)
                shift
                [[ $# -gt 0 ]] || fail "Missing value for --disk-warn-gb"
                DISK_WARN_GB=$1
                ;;
            --help)
                usage
                exit 0
                ;;
            *)
                fail "Unknown option: $1"
                ;;
        esac
        shift
    done
}

is_verified() {
    local verified_file="${BASE_DIR}/${OPTS_LABEL}/${CHR}/.verified"

    if [[ ! -f "$verified_file" ]]; then
        return 1
    fi

    local stored_mtime
    stored_mtime=$(grep '^RUST_BINARY_MTIME=' "$verified_file" | cut -d= -f2)
    local current_mtime
    current_mtime=$(stat -c %Y "$RUST_BIN" 2>/dev/null) || return 1
    [[ "$stored_mtime" == "$current_mtime" ]]
}

write_verified() {
    local verified_file="${BASE_DIR}/${OPTS_LABEL}/${CHR}/.verified"
    local rust_mtime
    rust_mtime=$(stat -c %Y "$RUST_BIN" 2>/dev/null || echo unknown)

    mkdir -p "$(dirname "$verified_file")"
    {
        printf 'RUST_BINARY_MTIME=%s\n' "$rust_mtime"
        printf 'PASS_COUNT=%s\n' "$PASS_COUNT"
        printf 'EMPTY_COUNT=%s\n' "$EMPTY_COUNT"
        printf 'TOTAL_SHARDS=%s\n' "$TOTAL_SHARDS"
        printf 'VERIFIED_AT=%s\n' "$(date -Iseconds)"
    } >"$verified_file"
}

json_escape() {
    local value=$1

    value=${value//\\/\\\\}
    value=${value//\"/\\\"}
    value=${value//$'\n'/\\n}
    value=${value//$'\r'/\\r}
    value=${value//$'\t'/\\t}
    printf '%s' "$value"
}

write_results_json() {
    local elapsed_seconds=$1
    local results_file="${BASE_DIR}/${OPTS_LABEL}/${CHR}/results.json"
    local pass_without_empty=$((PASS_COUNT - EMPTY_COUNT))
    local idx=0
    local entry
    local shard
    local region
    local reason
    local first_diff_line

    mkdir -p "$(dirname "$results_file")"
    {
        echo '{'
        printf '  "chromosome": "%s",\n' "$(json_escape "$CHR")"
        printf '  "opts_label": "%s",\n' "$(json_escape "$OPTS_LABEL")"
        printf '  "total_shards": %d,\n' "$TOTAL_SHARDS"
        printf '  "pass": %d,\n' "$pass_without_empty"
        printf '  "fail": %d,\n' "$FAIL_COUNT"
        printf '  "empty": %d,\n' "$EMPTY_COUNT"
        printf '  "elapsed_seconds": %s,\n' "$elapsed_seconds"
        echo '  "failures": ['
        for entry in "${JSON_FAILURES[@]}"; do
            IFS='|' read -r shard region reason first_diff_line <<<"$entry"
            if (( idx > 0 )); then
                echo ','
            fi
            printf '    {"shard": "%s", "region": "%s", "reason": "%s", "first_diff_line": ' \
                "$(json_escape "$shard")" \
                "$(json_escape "$region")" \
                "$(json_escape "$reason")"
            if [[ -n "$first_diff_line" ]]; then
                printf '%s' "$first_diff_line"
            else
                printf 'null'
            fi
            printf '}'
            idx=$((idx + 1))
        done
        echo
        echo '  ]'
        echo '}'
    } >"$results_file"
}

cleanup_outputs() {
    if (( FAIL_COUNT == 0 )); then
        rm -rf "$JAVA_DIR" "$RUST_DIR" "$DIFF_DIR"
        echo "  Cleaned up all outputs (verified)"
    else
        local idx
        local status_file
        local status_value
        local shard_num
        local cleaned=0

        for ((idx = 0; idx < TOTAL_SHARDS; idx++)); do
            status_file=$(status_file_for_idx "$idx")
            if [[ -f "$status_file" ]]; then
                status_value=$(<"$status_file")
                if [[ "$status_value" == "PASS" || "$status_value" == "EMPTY" ]]; then
                    shard_num=${SHARD_IDS[$idx]}
                    rm -f "${JAVA_DIR}/shard_${shard_num}.tsv" "${JAVA_DIR}/shard_${shard_num}.log"
                    rm -f "${RUST_DIR}/shard_${shard_num}.tsv" "${RUST_DIR}/shard_${shard_num}.log"
                    rm -f "$status_file" "$(meta_file_for_idx "$idx")"
                    rm -f "$(phase_state_file_for_idx java "$idx")" "$(phase_state_file_for_idx rust "$idx")"
                    cleaned=$((cleaned + 1))
                fi
            fi
        done
        echo "  Cleaned ${cleaned} passing shard outputs (kept ${FAIL_COUNT} failing)"
    fi
}

run_parallel_phase() {
    local phase_func=$1
    local phase_label=$2
    local -n phase_indices_ref=$3
    local total_phase_shards=${#phase_indices_ref[@]}
    local idx
    local launch_start_epoch

    if (( total_phase_shards == 0 )); then
        return 0
    fi

    mkdir -p ./tmp
    FIFO_SEM="./tmp/parity_sem_$$"
    PHASE_PROGRESS_FILE="./tmp/parity_progress_$$.log"
    PHASE_ABORT_FLAG="./tmp/parity_abort_$$.flag"
    rm -f "$FIFO_SEM" "$PHASE_PROGRESS_FILE" "$PHASE_ABORT_FLAG"

    mkfifo "$FIFO_SEM"
    exec 3<>"$FIFO_SEM"
    FIFO_OPEN=1
    launch_start_epoch=$(date +%s)

    for ((idx = 0; idx < MAX_PARALLEL; idx++)); do
        printf '\n' >&3
    done

    start_progress_reporter "$phase_label" "$total_phase_shards" "$PHASE_PROGRESS_FILE" "$PHASE_ABORT_FLAG" "$launch_start_epoch"

    ACTIVE_PIDS=()
    for idx in "${phase_indices_ref[@]}"; do
        if [[ -f "$PHASE_ABORT_FLAG" ]]; then
            break
        fi

        read -r -u 3
        (
            local worker_rc=0

            "$phase_func" "$idx" || worker_rc=$?
            printf '%s\n' "${SHARD_IDS[$idx]}" >>"$PHASE_PROGRESS_FILE"
            printf '\n' >&3
            exit "$worker_rc"
        ) &
        ACTIVE_PIDS+=("$!")
    done

    wait_for_batch_jobs || true
    print_phase_progress_line "$phase_label" "$total_phase_shards" "$PHASE_PROGRESS_FILE" "$launch_start_epoch"
    stop_progress_reporter
    close_fifo_semaphore

    if [[ -f "$PHASE_ABORT_FLAG" ]]; then
        cleanup_temp_files
        return 1
    fi

    cleanup_temp_files
    return 0
}

run_chr() {
    local start_ns
    local end_ns
    local idx
    local -a java_needed_indices=()
    local -a rust_needed_indices=()

    if (( ALL_CHR )); then
        CHR_LEN=""
    else
        CHR_LEN="$REQUESTED_CHR_LEN"
    fi
    resolve_chr_length
    generate_shards

    if is_verified; then
        echo "=== $(summary_chr_label): already verified (${TOTAL_SHARDS} shards) ==="
        PASS_COUNT=$TOTAL_SHARDS
        FAIL_COUNT=0
        EMPTY_COUNT=0
        LAST_CHR_ELAPSED="0.0s"
        LAST_CHR_ELAPSED_SECONDS="0.0"
        JSON_FAILURES=()
        write_results_json "$LAST_CHR_ELAPSED_SECONDS"
        GLOBAL_PASS=$((GLOBAL_PASS + TOTAL_SHARDS))
        GLOBAL_SHARDS=$((GLOBAL_SHARDS + TOTAL_SHARDS))
        GLOBAL_CHR_PASS=$((GLOBAL_CHR_PASS + 1))
        return 0
    fi

    prepare_dirs
    ensure_resource_log

    needs_missing_outputs "java" java_needed_indices || true
    needs_missing_outputs "rust" rust_needed_indices || true

    if (( ${#java_needed_indices[@]} > 0 )); then
        require_file "$JAVA_JAR"
    fi
    if (( ${#rust_needed_indices[@]} > 0 )); then
        require_file "$RUST_BIN"
    fi

    PASS_COUNT=0
    FAIL_COUNT=0
    EMPTY_COUNT=0
    FAIL_MESSAGES=()
    JSON_FAILURES=()

    echo "=== Processing $(summary_chr_label) (${TOTAL_SHARDS} shards, ${MAX_PARALLEL} workers) ==="
    start_ns=$(date +%s%N)

    if (( ${#java_needed_indices[@]} > 0 )); then
        check_preflight_resources "Phase 1 (Java)"
        echo "--- Phase 1: Java (${#java_needed_indices[@]} shards, ${MAX_PARALLEL} parallel) ---"
        run_parallel_phase generate_java_shard "Phase 1 (Java)" java_needed_indices || fail "Phase 1 (Java) aborted due to low memory"
    else
        echo "--- Phase 1: Java (all cached) ---"
    fi

    if (( ${#rust_needed_indices[@]} > 0 )); then
        check_preflight_resources "Phase 2 (Rust)"
        echo "--- Phase 2: Rust (${#rust_needed_indices[@]} shards, ${MAX_PARALLEL} parallel) ---"
        run_parallel_phase generate_rust_shard "Phase 2 (Rust)" rust_needed_indices || fail "Phase 2 (Rust) aborted due to low memory"
    else
        echo "--- Phase 2: Rust (all cached) ---"
    fi

    echo "--- Phase 3: Compare ---"
    BATCH_FIRST_FAIL_IDX=""
    for ((idx = 0; idx < TOTAL_SHARDS; idx++)); do
        compare_shard "$idx"
    done
    collect_batch_results 0 $((TOTAL_SHARDS - 1))

    end_ns=$(date +%s%N)
    LAST_CHR_ELAPSED=$(format_elapsed "$start_ns" "$end_ns")
    LAST_CHR_ELAPSED_SECONDS=$(format_elapsed_seconds "$start_ns" "$end_ns")
    write_results_json "$LAST_CHR_ELAPSED_SECONDS"

    if [[ -n "$BATCH_FIRST_FAIL_IDX" ]] && (( STOP_ON_FAIL )); then
        print_failure_details "$BATCH_FIRST_FAIL_IDX"
    fi

    print_summary

    if (( FAIL_COUNT == 0 )); then
        write_verified
    fi

    if (( CLEANUP )); then
        echo "--- Phase 4: Cleanup ---"
        cleanup_outputs
    fi

    GLOBAL_PASS=$((GLOBAL_PASS + PASS_COUNT))
    GLOBAL_FAIL=$((GLOBAL_FAIL + FAIL_COUNT))
    GLOBAL_EMPTY=$((GLOBAL_EMPTY + EMPTY_COUNT))
    GLOBAL_SHARDS=$((GLOBAL_SHARDS + TOTAL_SHARDS))
    if (( FAIL_COUNT > 0 )); then
        GLOBAL_CHR_FAIL=$((GLOBAL_CHR_FAIL + 1))
        GLOBAL_FAIL_DETAILS+=("  $(summary_chr_label): ${FAIL_COUNT}/${TOTAL_SHARDS} failed")
    else
        GLOBAL_CHR_PASS=$((GLOBAL_CHR_PASS + 1))
    fi

    return 0
}

main() {
    parse_args "$@"

    validate_positive_integer "shard size" "$SHARD_SIZE"
    validate_positive_integer "parallel worker count" "$MAX_PARALLEL"
    validate_positive_integer "timeout" "$SHARD_TIMEOUT"
    validate_positive_integer "mem warn threshold" "$MEM_WARN_GB"
    validate_positive_integer "mem abort threshold" "$MEM_ABORT_GB"
    validate_positive_integer "disk warn threshold" "$DISK_WARN_GB"
    if (( MEM_WARN_GB <= MEM_ABORT_GB )); then
        fail "mem warn threshold must be greater than mem abort threshold"
    fi
    if (( DISK_WARN_GB <= DISK_ABORT_GB )); then
        fail "disk warn threshold must be greater than disk abort threshold (${DISK_ABORT_GB}GB)"
    fi
    if [[ -n "$REQUESTED_CHR_LEN" ]]; then
        validate_positive_integer "chromosome length" "$REQUESTED_CHR_LEN"
    fi

    mkdir -p ./tmp
    require_file "$REF_FASTA"
    require_file "$BAM_PATH"

    if (( ALL_CHR )); then
        require_file "$REF_FASTA_FAI"
        while IFS=$'\t' read -r chr_name _; do
            CHR="$chr_name"
            REQUESTED_CHR_LEN=""
            run_chr || exit 1
        done < "$REF_FASTA_FAI"

        print_global_summary
        if (( GLOBAL_FAIL > 0 )); then
            exit 1
        fi
        return 0
    fi

    run_chr
    if (( FAIL_COUNT > 0 )); then
        exit 1
    fi
}

main "$@"