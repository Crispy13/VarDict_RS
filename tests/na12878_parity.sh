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
CLEAN=0
CLEAN_RUST=0

TOTAL_SHARDS=0
PASS_COUNT=0
FAIL_COUNT=0
EMPTY_COUNT=0
LAST_CHR_ELAPSED="0.0s"

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

JAVA_SEM_FD=""
BATCH_FIRST_FAIL_IDX=""

usage() {
    cat <<EOF
Usage: ${SCRIPT_NAME} [options]

Options:
  --chr CHR         Chromosome name (default: 20)
  --all-chr         Run all chromosomes from the FAI index sequentially
  --chr-len LEN     Chromosome length; defaults to lookup in ${REF_FASTA_FAI}
  --shard-size N    Override shard size in bases (default: 1000000)
  --parallel N      Override shard pipeline worker count (default: 10)
  --java-parallel N Override concurrent Java worker count within shard pipelines (default: 10)
  --java-heap SIZE  Java heap size per worker (default: 2g)
  --freq F          Override VarDict frequency threshold (default: 0.01)
  --no-stop         Continue past mismatches instead of stopping after the first failed batch
  --clean           Remove all cached outputs (Java + Rust) before running
  --clean-rust      Remove only Rust + diff outputs (preserves Java cache)
  --help            Show this help

Caching:
  Java outputs are cached and reused across runs (deterministic).
  Rust outputs are automatically regenerated when the binary is newer.
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

remove_stale_diff() {
    local diff_file=$1

    if [[ -f "$diff_file" ]]; then
        rm -f "$diff_file"
    fi
}

kill_active_jobs() {
    local pid

    for pid in "${ACTIVE_PIDS[@]}"; do
        if [[ -n "$pid" ]] && kill -0 "$pid" 2>/dev/null; then
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

cleanup_java_semaphore() {
    if [[ -n "${JAVA_SEM_FD:-}" ]]; then
        exec {JAVA_SEM_FD}>&-
        JAVA_SEM_FD=""
    fi
}

cleanup_script() {
    cleanup_java_semaphore
    kill_active_jobs
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
    JAVA_DIR="${BASE_DIR}/${CHR}/java"
    RUST_DIR="${BASE_DIR}/${CHR}/rust"
    DIFF_DIR="${BASE_DIR}/${CHR}/diff"

    if (( CLEAN )); then
        rm -rf "${BASE_DIR}/${CHR}"
    elif (( CLEAN_RUST )); then
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

    awk '
        NR == FNR { left[FNR] = $0; left_count = FNR; next }
        { right[FNR] = $0; right_count = FNR }
        END {
            max = (left_count > right_count) ? left_count : right_count
            for (i = 1; i <= max; i++) {
                if (left[i] != right[i]) {
                    print i
                    exit
                }
            }
            if (left_count != right_count) {
                print max + 1
            }
        }
    ' "$left_file" "$right_file"
}

init_java_semaphore() {
    local control_dir="${BASE_DIR}/${CHR}/.control"
    local fifo_path="${control_dir}/java_tokens.fifo"
    local idx

    cleanup_java_semaphore
    mkdir -p "$control_dir"
    rm -f "$fifo_path"
    mkfifo "$fifo_path"
    exec {JAVA_SEM_FD}<>"$fifo_path"
    rm -f "$fifo_path"

    for ((idx = 0; idx < JAVA_MAX_PARALLEL; idx++)); do
        printf '.' >&"$JAVA_SEM_FD"
    done
}

acquire_java_slot() {
    local token=''

    read -r -n 1 -u "$JAVA_SEM_FD" token
}

release_java_slot() {
    printf '.' >&"$JAVA_SEM_FD"
}

needs_missing_outputs() {
    local mode=$1
    local idx
    local out_file

    for ((idx = 0; idx < TOTAL_SHARDS; idx++)); do
        if [[ "$mode" == "java" ]]; then
            out_file="${JAVA_DIR}/shard_${SHARD_IDS[$idx]}.tsv"
        else
            out_file="${RUST_DIR}/shard_${SHARD_IDS[$idx]}.tsv"
        fi

        if [[ ! -f "$out_file" ]]; then
            return 0
        fi

        # Rust outputs are stale if the binary is newer
        if [[ "$mode" == "rust" && -f "$out_file" && "$RUST_BIN" -nt "$out_file" ]]; then
            return 0
        fi
    done

    return 1
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

process_shard() {
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
    local java_rc=0
    local rust_rc=0
    local java_lines=0
    local rust_lines=0
    local first_diff_line=""
    local summary=""
    local reason=""

    status_file=$(status_file_for_idx "$idx")
    meta_file=$(meta_file_for_idx "$idx")
    diff_file=$(diff_file_for_idx "$idx")

    rm -f "$status_file" "$meta_file" "$diff_file"

    if [[ ! -f "$java_file" ]]; then
        acquire_java_slot || {
            printf 'FAIL\n' >"$status_file"
            printf 'Unable to acquire Java worker slot\n' >"$diff_file"
            write_failure_meta "$meta_file" "$region" "$shard_num" "$java_file" "$rust_file" "$diff_file" "$java_log" "$rust_log" "0" "0" "" "java worker slot unavailable" "java_slot_unavailable"
            return 0
        }

        java "-Xmx${JAVA_HEAP}" -classpath "$JAVA_JAR" com.astrazeneca.vardict.Main \
            -G "$REF_FASTA" \
            -b "$BAM_PATH" \
            -N NA12878 -f "$FREQ" -th 1 \
            -R "$region" \
            >"$java_file" 2>"$java_log" || java_rc=$?

        release_java_slot

        if (( java_rc != 0 )); then
            rm -f "$java_file"
            printf 'Java command failed with exit code %s\n' "$java_rc" >"$diff_file"
            printf 'FAIL\n' >"$status_file"
            write_failure_meta "$meta_file" "$region" "$shard_num" "$java_file" "$rust_file" "$diff_file" "$java_log" "$rust_log" "0" "0" "" "java command failed (exit ${java_rc})" "java_command_failed"
            return 0
        fi
    fi

    if [[ ! -f "$rust_file" ]] || [[ -f "$rust_file" && "$RUST_BIN" -nt "$rust_file" ]]; then
        rm -f "$rust_file"
        "$RUST_BIN" \
            -G "$REF_FASTA" \
            -b "$BAM_PATH" \
            -N NA12878 -f "$FREQ" \
            -R "$region" \
            >"$rust_file" 2>"$rust_log" || rust_rc=$?

        if (( rust_rc != 0 )); then
            rm -f "$rust_file"
            printf 'Rust command failed with exit code %s\n' "$rust_rc" >"$diff_file"
            printf 'FAIL\n' >"$status_file"
            write_failure_meta "$meta_file" "$region" "$shard_num" "$java_file" "$rust_file" "$diff_file" "$java_log" "$rust_log" "0" "0" "" "rust command failed (exit ${rust_rc})" "rust_command_failed"
            return 0
        fi
    fi

    if [[ ! -f "$java_file" || ! -f "$rust_file" ]]; then
        printf 'FAIL\n' >"$status_file"
        printf 'Missing output for comparison\n' >"$diff_file"
        write_failure_meta "$meta_file" "$region" "$shard_num" "$java_file" "$rust_file" "$diff_file" "$java_log" "$rust_log" "0" "0" "" "missing output" "missing_output"
        return 0
    fi

    if [[ ! -s "$java_file" && ! -s "$rust_file" ]]; then
        printf 'EMPTY\n' >"$status_file"
        return 0
    fi

    # Detect asymmetric empty output — likely a transient failure on one side
    if [[ ! -s "$java_file" && -s "$rust_file" ]]; then
        rust_lines=$(count_file_lines "$rust_file")
        printf 'FAIL\n' >"$status_file"
        summary="SUSPECT: java output empty but rust has ${rust_lines} lines (likely transient Java failure — delete java cache and re-run)"
        write_failure_meta "$meta_file" "$region" "$shard_num" "$java_file" "$rust_file" "$diff_file" "$java_log" "$rust_log" "0" "$rust_lines" "" "$summary" "java_empty_suspect"
        return 0
    fi
    if [[ -s "$java_file" && ! -s "$rust_file" ]]; then
        java_lines=$(count_file_lines "$java_file")
        printf 'FAIL\n' >"$status_file"
        summary="SUSPECT: rust output empty but java has ${java_lines} lines (likely transient Rust failure — delete rust cache and re-run)"
        write_failure_meta "$meta_file" "$region" "$shard_num" "$java_file" "$rust_file" "$diff_file" "$java_log" "$rust_log" "$java_lines" "0" "" "$summary" "rust_empty_suspect"
        return 0
    fi

    if cmp -s "$java_file" "$rust_file"; then
        printf 'PASS\n' >"$status_file"
        remove_stale_diff "$diff_file"
        return 0
    fi

    java_lines=$(count_file_lines "$java_file")
    rust_lines=$(count_file_lines "$rust_file")
    first_diff_line=$(find_first_difference_line "$java_file" "$rust_file")
    diff -u "$java_file" "$rust_file" >"$diff_file" || true

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
    write_failure_meta "$meta_file" "$region" "$shard_num" "$java_file" "$rust_file" "$diff_file" "$java_log" "$rust_log" "$java_lines" "$rust_lines" "$first_diff_line" "$summary" "$reason"
}

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

    for ((idx = batch_start; idx <= batch_end; idx++)); do
        shard_num=${SHARD_IDS[$idx]}
        start=${SHARD_STARTS[$idx]}
        end=${SHARD_ENDS[$idx]}
        status_file=$(status_file_for_idx "$idx")

        if [[ ! -f "$status_file" ]]; then
            FAIL_COUNT=$((FAIL_COUNT + 1))
            FAIL_MESSAGES+=("  shard_${shard_num} ($(display_region "$start" "$end")): missing status")
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
                if [[ -z "$BATCH_FIRST_FAIL_IDX" ]]; then
                    BATCH_FIRST_FAIL_IDX=$idx
                fi
                ;;
            *)
                FAIL_COUNT=$((FAIL_COUNT + 1))
                FAIL_MESSAGES+=("  shard_${shard_num} ($(display_region "$start" "$end")): unknown status '${status_value}'")
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
    if [[ "$REASON" == "java_command_failed" ]]; then
        echo "  Java log: ${JAVA_LOG}"
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
            --clean)
                CLEAN=1
                ;;
            --clean-rust)
                CLEAN_RUST=1
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

run_chr() {
    local start_ns
    local end_ns
    local batch_start=0
    local batch_end=0
    local idx
    if (( ALL_CHR )); then
        CHR_LEN=""
    else
        CHR_LEN="$REQUESTED_CHR_LEN"
    fi
    resolve_chr_length
    generate_shards
    prepare_dirs

    if needs_missing_outputs "java"; then
        require_file "$JAVA_JAR"
    fi
    if needs_missing_outputs "rust"; then
        require_file "$RUST_BIN"
    fi

    PASS_COUNT=0
    FAIL_COUNT=0
    EMPTY_COUNT=0
    FAIL_MESSAGES=()

    init_java_semaphore

    echo "=== Processing $(summary_chr_label) (${TOTAL_SHARDS} shards, ${MAX_PARALLEL} pipelines, ${JAVA_MAX_PARALLEL} Java workers) ==="
    start_ns=$(date +%s%N)

    while (( batch_start < TOTAL_SHARDS )); do
        batch_end=$((batch_start + MAX_PARALLEL - 1))
        if (( batch_end >= TOTAL_SHARDS )); then
            batch_end=$((TOTAL_SHARDS - 1))
        fi

        ACTIVE_PIDS=()
        for ((idx = batch_start; idx <= batch_end; idx++)); do
            process_shard "$idx" &
            ACTIVE_PIDS+=("$!")
        done

        wait_for_batch_jobs || fail "Shard pipeline execution failed unexpectedly"
        collect_batch_results "$batch_start" "$batch_end"

        if [[ -n "$BATCH_FIRST_FAIL_IDX" ]] && (( STOP_ON_FAIL )); then
            end_ns=$(date +%s%N)
            LAST_CHR_ELAPSED=$(format_elapsed "$start_ns" "$end_ns")
            print_failure_details "$BATCH_FIRST_FAIL_IDX"
            cleanup_java_semaphore
            return 1
        fi

        batch_start=$((batch_end + 1))
    done

    end_ns=$(date +%s%N)
    LAST_CHR_ELAPSED=$(format_elapsed "$start_ns" "$end_ns")
    cleanup_java_semaphore
    print_summary

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
    validate_positive_integer "java parallel worker count" "$JAVA_MAX_PARALLEL"
    if [[ -n "$REQUESTED_CHR_LEN" ]]; then
        validate_positive_integer "chromosome length" "$REQUESTED_CHR_LEN"
    fi

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