#!/usr/bin/env bash
set -euo pipefail

SCRIPT_NAME=$(basename "${BASH_SOURCE[0]}")
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR" && git rev-parse --show-toplevel)
cd "$PROJECT_ROOT"

PARITY_SCRIPT="${SCRIPT_DIR}/na12878_parity_v2.sh"
REF_FASTA_FAI="testdata/hs37d5.fa.fai"
RUST_BIN="target/debug-release/vardict"
BASE_DIR="tmp/na12878_parity"
REPORT_FILE="${BASE_DIR}/option_parity_report.json"
CHECKPOINT_FILE="${BASE_DIR}/.sweep_checkpoint"
SHARD_SIZE=1000000

declare -a DEFAULT_CHRS=(20 22 MT)
declare -a TEST_CHRS=("${DEFAULT_CHRS[@]}")
declare -a TIER_FILTERS=()
declare -a CONFIG_FILTERS=()
declare -a BLOCKED_CONFIGS=()
declare -a TEST_MATRIX=(
    # ── Tier 1: Core pipeline (cheap, highest signal) ──────────────────
    "T1-03|nosv|-U|1"
    "T1-05|freq-low|-f 0.001|1"
    "T1-13|fisher|--fisher|1"
    "T1-04|nosv-dedup|-U --deldupvar|1"
    "T1-06|freq-high|-f 0.05|1"
    "T1-07|freq-zero|-f 0.0|1"
    "T1-08|minr-1|-r 1|1"
    "T1-09|minr-5|-r 5|1"
    "T1-11|mapq-10|-Q 10|1"
    "T1-12|mapq-30|-Q 30|1"
    # ── Tier 2: Filter / variant options (mid-cost) ────────────────────
    "T2-09|filter-0x700|-F 0x700|2"
    "T2-16|chimeric|--chimeric|2"
    "T2-01|qual-0|-q 0|2"
    "T2-02|qual-10|-q 10|2"
    "T2-03|qual-30|-q 30|2"
    "T2-04|qual-40|-q 40|2"
    "T2-05|vext-0|-X 0|2"
    "T2-06|vext-5|-X 5|2"
    "T2-07|mismatch-3|-m 3|2"
    "T2-08|mismatch-15|-m 15|2"
    "T2-10|filter-0x100|-F 0x100|2"
    "T2-11|indel-3prime|-3|2"
    "T2-12|readpos-0|-P 0|2"
    "T2-13|readpos-10|-P 10|2"
    "T2-14|qratio-low|-o 0.5|2"
    "T2-15|qratio-high|-o 2.0|2"
    # ── Tier 3: Combo configs ──────────────────────────────────────────
    "T3-01|clinical-wgs|-f 0.001 -Q 10 -F 0x700|3"
    "T4-04|inssize-small|-w 200 -W 50|4"
    "T3-04|fisher-lowfreq|-f 0.001 --fisher|3"
    "T3-05|nosv-tight|-U -f 0.001 -Q 10|3"
    "T3-06|3prime-lowfreq|-3 -f 0.001|3"
    "T3-08|quality-gauntlet|-M 25 -m 5 -Q 10|3"
    "T4-01|extend-150|-x 150|4"
    "T4-02|refext-600|-Y 600|4"
    "T4-03|refext-2000|-Y 2000|4"
    "T4-05|trim-130|-T 130|4"
    # ── Tier 4: Edge-case / expensive (known failure-prone, run last) ──
    "T1-10|no-realign|-k 0|1"
    "T1-14|debug|-D|1"
    "T3-03|no-realign-lowfreq|-k 0 -f 0.001|3"
    # ── Tier 5: Most expensive (pileup — high memory, longest runtime) ─
    "T1-01|pileup|-p|1"
    "T1-02|pileup-max|-p -f 0.0 -r 1|1"
    "T3-02|pileup-strict|-p -q 30|3"
    "T3-07|pileup-relaxed|-p -r 1 -q 10 -Q 0|3"
    "T4-06|include-n-pileup|-K -p -r 1 -f 0.0|4"
)
declare -a RESULTS=()
declare -a COMPLETED_CONFIGS=()
declare -A DONE_MAP=()

USE_ALL_CHR=0
CHR_OVERRIDE=0
NO_STOP=0
NO_CLEANUP=0
FRESH_START=false
RUST_ONLY=0
CLEAN_ALL=0
RELEASE_BUILD=0
RUST_BIN_OVERRIDE=""
NO_BUILD=0
CONFIG_OVERRIDE=0
TIER_USED=0
PRESET=""
REAL_CHRS_ONLY=0
INCLUDE_BLOCKED=0
FULL_GATE=0
DRY_RUN=0
STATUS_ONLY=0
PARALLEL_WORKERS=""
TIMEOUT_SECONDS=""
RETRY_COUNT=0
MEM_WARN_GB=""
MEM_ABORT_GB=""
DISK_WARN_GB=""
_PER_CONFIG_EXTRA_FLAGS=()
TOTAL_CONFIGS=0
PASS_CONFIGS=0
FAIL_CONFIGS=0

usage() {
    cat <<EOF
Usage: ${SCRIPT_NAME} [options]

Options:
  --tier N            Run only tier N (1-4). Can be repeated.
  --config-id ID      Run only a specific config ID. Can be repeated.
    --preset NAME       Apply a standard config+chromosome preset.
  --chr CHR           Override chromosomes (default: 20 22 MT). Can be repeated.
  --all-chr           Run each selected config across all contigs from the FAI index.
  --no-stop           Continue to later configs after a failure.
  --no-cleanup        Keep all output files after comparison (disable disk cleanup).
  --rust-only         Skip Java generation; fail if Java cache is incomplete
    --clean-all         Delete ALL outputs including Java cache on cleanup
    --release           Pass --release through to na12878_parity_v2.sh
    --rust-bin PATH     Pass a specific Rust binary path through to na12878_parity_v2.sh
    --no-build          Pass --no-build through to na12878_parity_v2.sh
    --parallel N        Override parallel worker count (default: 10).
    --dry-run           Print the configs/chromosomes that would run without executing.
    --status            Print the current config x chromosome parity status matrix.
    --fresh             Delete any existing sweep checkpoint before running.
    --timeout SECONDS   Pass through per-shard timeout to na12878_parity_v2.sh.
    --retry N           Pass through transient shard retry count to na12878_parity_v2.sh.
  --mem-warn-gb N     Pass through memory warning threshold.
  --mem-abort-gb N    Pass through memory abort threshold.
  --disk-warn-gb N    Pass through disk warning threshold.
    --include-blocked   Run blocked configs (normally skipped).
  --help              Show this help.

Matrix:
  Tier 1: 14 configs
  Tier 2: 16 configs
  Tier 3: 8 configs
  Tier 4: 6 configs
  Total:  44 configs

Presets:
    smoke:         3 configs x 3 chromosomes  (T1-01, T1-03, T1-13 on 20, 22, MT)
    dev:          10 configs x 10 chromosomes (repr set on 1, 2, 5, 10, 14, 17, 20, 22, X, MT)
    tier1:        14 configs x 3 chromosomes  (all T1-* on 20, 22, MT)
    config-spread: 44 configs x 3 chromosomes (all configs on 20, 22, MT)
    core-wide:    14 configs x 25 chromosomes (all T1-* on all real chromosomes)
    pairwise:     20 configs x 10 chromosomes (PW-000..PW-019 on repr set)
    full-gate:    sequential smoke -> tier1 -> config-spread -> core-wide
    release:      all 44 configs x all real chromosomes from the FAI index
EOF
}

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

validate_positive_integer() {
    local label=$1
    local value=$2

    [[ "$value" =~ ^[0-9]+$ ]] || fail "$label must be a positive integer: $value"
    (( value > 0 )) || fail "$label must be greater than zero: $value"
}

validate_nonnegative_integer() {
    local label=$1
    local value=$2

    [[ "$value" =~ ^[0-9]+$ ]] || fail "$label must be a non-negative integer: $value"
}

array_contains() {
    local needle=$1
    shift
    local item

    for item in "$@"; do
        if [[ "$item" == "$needle" ]]; then
            return 0
        fi
    done

    return 1
}

config_uses_pileup() {
    local cli_flags=$1

    [[ "$cli_flags" =~ (^|[[:space:]])-p($|[[:space:]]) ]]
}

scope_label() {
    if (( USE_ALL_CHR )); then
        printf 'all-chr'
    else
        printf '%s' "${TEST_CHRS[*]}"
    fi
}

render_chr_label() {
    local chr=$1

    if [[ "$chr" == chr* ]]; then
        printf '%s' "$chr"
    else
        printf 'chr%s' "$chr"
    fi
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

display_escape() {
    local value=$1

    value=${value//\\/\\\\}
    value=${value//\"/\\\"}
    printf '%s' "$value"
}

write_checkpoint() {
    local tmp_file="${CHECKPOINT_FILE}.tmp"
    local config_id

    mkdir -p "$BASE_DIR"
    {
        echo '# Sweep checkpoint — auto-generated, do not edit'
        if (( ${#COMPLETED_CONFIGS[@]} == 0 )); then
            echo 'COMPLETED_CONFIGS=()'
        else
            printf 'COMPLETED_CONFIGS=('
            for config_id in "${COMPLETED_CONFIGS[@]}"; do
                printf ' "%s"' "$config_id"
            done
            echo ' )'
        fi
    } >"$tmp_file"
    mv "$tmp_file" "$CHECKPOINT_FILE"
}

load_checkpoint() {
    local config_id
    local loaded_count=0

    if [[ "$FRESH_START" == true && -f "$CHECKPOINT_FILE" ]]; then
        rm -f "$CHECKPOINT_FILE"
    fi

    COMPLETED_CONFIGS=()
    DONE_MAP=()

    if [[ -f "$CHECKPOINT_FILE" ]]; then
        # shellcheck disable=SC1090
        source "$CHECKPOINT_FILE"
        for config_id in "${COMPLETED_CONFIGS[@]}"; do
            DONE_MAP["$config_id"]=1
        done
        loaded_count=${#COMPLETED_CONFIGS[@]}
    fi

    echo "Loaded ${loaded_count} previously completed configs from checkpoint"
}

format_elapsed_seconds() {
    local start_ns=$1
    local end_ns=$2

    awk -v start="$start_ns" -v end="$end_ns" 'BEGIN { printf "%.1f", (end - start) / 1000000000 }'
}

format_duration_human() {
    local raw_seconds=$1
    local total_seconds
    local hours
    local minutes
    local seconds

    total_seconds=$(awk -v value="$raw_seconds" 'BEGIN { if (value < 0) value = 0; printf "%.0f", value }')
    hours=$((total_seconds / 3600))
    minutes=$(((total_seconds % 3600) / 60))
    seconds=$((total_seconds % 60))

    if (( hours > 0 )); then
        printf '%dh%02dm%02ds' "$hours" "$minutes" "$seconds"
    elif (( minutes > 0 )); then
        printf '%dm%02ds' "$minutes" "$seconds"
    else
        printf '%ds' "$seconds"
    fi
}

read_marker_field() {
    local file=$1
    local key=$2

    if [[ ! -f "$file" ]]; then
        return 1
    fi

    awk -F= -v key="$key" '$1 == key { sub(/^[^=]*=/, ""); print; exit }' "$file"
}

current_rust_binary_mtime() {
    if [[ -f "$RUST_BIN" ]]; then
        stat -c %Y "$RUST_BIN" 2>/dev/null || true
    fi
}

configure_rust_bin() {
    if (( RELEASE_BUILD )); then
        RUST_BIN="target/debug-release/vardict"
    fi

    if [[ -n "$RUST_BIN_OVERRIDE" ]]; then
        RUST_BIN="$RUST_BIN_OVERRIDE"
    fi
}

require_fai() {
    [[ -f "$REF_FASTA_FAI" ]] || fail "Reference FAI not found: $REF_FASTA_FAI"
}

resolve_scope_chrs() {
    local -n out_ref=$1

    if (( USE_ALL_CHR )); then
        require_fai
        if (( REAL_CHRS_ONLY )); then
            mapfile -t out_ref < <(awk -F '\t' '$1 ~ /^([1-9]|1[0-9]|2[0-2]|X|Y|MT)$/ { print $1 }' "$REF_FASTA_FAI")
        else
            mapfile -t out_ref < <(awk -F '\t' '{ print $1 }' "$REF_FASTA_FAI")
        fi
    else
        out_ref=("${TEST_CHRS[@]}")
    fi
}

lookup_chr_length() {
    local chr=$1

    require_fai
    awk -F '\t' -v chr="$chr" '$1 == chr { print $2; exit }' "$REF_FASTA_FAI"
}

shard_count_for_chr() {
    local chr=$1
    local chr_len

    chr_len=$(lookup_chr_length "$chr")
    [[ -n "$chr_len" ]] || fail "Chromosome '$chr' not found in $REF_FASTA_FAI"
    printf '%d' "$(((chr_len + SHARD_SIZE - 1) / SHARD_SIZE))"
}

should_run() {
    local config_id=$1
    local tier=$2

    if (( ${#TIER_FILTERS[@]} > 0 )) && ! array_contains "$tier" "${TIER_FILTERS[@]}"; then
        return 1
    fi

    if (( ${#CONFIG_FILTERS[@]} > 0 )) && ! array_contains "$config_id" "${CONFIG_FILTERS[@]}"; then
        return 1
    fi

    if (( INCLUDE_BLOCKED == 0 )) && (( ${#BLOCKED_CONFIGS[@]} > 0 )) && array_contains "$config_id" "${BLOCKED_CONFIGS[@]}"; then
        return 1
    fi

    return 0
}

parse_args() {
    while (( $# > 0 )); do
        case "$1" in
            --tier)
                shift
                [[ $# -gt 0 ]] || fail "Missing value for --tier"
                [[ "$1" =~ ^[1-4]$ ]] || fail "Tier must be 1, 2, 3, or 4: $1"
                TIER_USED=1
                TIER_FILTERS+=("$1")
                ;;
            --config-id)
                shift
                [[ $# -gt 0 ]] || fail "Missing value for --config-id"
                CONFIG_OVERRIDE=1
                CONFIG_FILTERS+=("$1")
                ;;
            --preset)
                shift
                [[ $# -gt 0 ]] || fail "Missing value for --preset"
                case "$1" in
                    smoke|dev|tier1|config-spread|core-wide|pairwise|full-gate|release)
                        PRESET="$1"
                        ;;
                    *)
                        fail "Preset must be one of: smoke, dev, tier1, config-spread, core-wide, pairwise, full-gate, release (got: $1)"
                        ;;
                esac
                ;;
            --chr)
                shift
                [[ $# -gt 0 ]] || fail "Missing value for --chr"
                if (( CHR_OVERRIDE == 0 )); then
                    TEST_CHRS=()
                fi
                CHR_OVERRIDE=1
                USE_ALL_CHR=0
                TEST_CHRS+=("$1")
                ;;
            --all-chr)
                CHR_OVERRIDE=1
                USE_ALL_CHR=1
                ;;
            --parallel)
                shift
                [[ $# -gt 0 ]] || fail "Missing value for --parallel"
                PARALLEL_WORKERS="$1"
                ;;
            --timeout)
                shift
                [[ $# -gt 0 ]] || fail "Missing value for --timeout"
                TIMEOUT_SECONDS="$1"
                ;;
            --retry)
                shift
                [[ $# -gt 0 ]] || fail "Missing value for --retry"
                validate_nonnegative_integer "retry count" "$1"
                RETRY_COUNT="$1"
                ;;
            --mem-warn-gb)
                shift
                [[ $# -gt 0 ]] || fail "Missing value for --mem-warn-gb"
                MEM_WARN_GB="$1"
                ;;
            --mem-abort-gb)
                shift
                [[ $# -gt 0 ]] || fail "Missing value for --mem-abort-gb"
                MEM_ABORT_GB="$1"
                ;;
            --disk-warn-gb)
                shift
                [[ $# -gt 0 ]] || fail "Missing value for --disk-warn-gb"
                DISK_WARN_GB="$1"
                ;;
            --no-stop)
                NO_STOP=1
                ;;
            --no-cleanup)
                NO_CLEANUP=1
                ;;
            --fresh)
                FRESH_START=true
                ;;
            --rust-only)
                RUST_ONLY=1
                ;;
            --clean-all)
                CLEAN_ALL=1
                ;;
            --release)
                RELEASE_BUILD=1
                ;;
            --rust-bin)
                shift
                [[ $# -gt 0 ]] || fail "Missing value for --rust-bin"
                RUST_BIN_OVERRIDE="$1"
                ;;
            --no-build)
                NO_BUILD=1
                ;;
            --dry-run)
                DRY_RUN=1
                ;;
            --status)
                STATUS_ONLY=1
                ;;
            --include-blocked)
                INCLUDE_BLOCKED=1
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

load_pairwise_configs() {
    local tsv_file="${PROJECT_ROOT}/tests/pairwise_configs.tsv"
    local line_num=0
    local index
    local label
    local cli_flags
    local rest

    [[ -f "$tsv_file" ]] || fail "Pairwise config file not found: $tsv_file"

    while IFS=$'\t' read -r index label cli_flags rest; do
        line_num=$((line_num + 1))
        (( line_num == 1 )) && continue
        [[ -z "$index" ]] && continue
        TEST_MATRIX+=("${label}|${label}|${cli_flags}|PW")
    done <"$tsv_file"
}

apply_preset() {
    if [[ -z "$PRESET" ]]; then
        return 0
    fi

    if (( TIER_USED )); then
        fail "--preset cannot be combined with --tier; use one selection mechanism"
    fi

    case "$PRESET" in
        smoke)
            if (( CONFIG_OVERRIDE == 0 )); then
                CONFIG_FILTERS=(T1-01 T1-03 T1-13)
            fi
            if (( CHR_OVERRIDE == 0 )); then
                USE_ALL_CHR=0
                TEST_CHRS=(20 22 MT)
            fi
            ;;
        dev)
            if (( CONFIG_OVERRIDE == 0 )); then
                CONFIG_FILTERS=(T1-01 T1-03 T1-05 T1-10 T1-13 T1-14 T2-09 T2-16 T3-01 T4-04)
            fi
            if (( CHR_OVERRIDE == 0 )); then
                USE_ALL_CHR=0
                TEST_CHRS=(1 2 5 10 14 17 20 22 X MT)
            fi
            ;;
        release)
            if (( CHR_OVERRIDE == 0 )); then
                REAL_CHRS_ONLY=1
                USE_ALL_CHR=1
            fi
            ;;
        tier1)
            if (( CONFIG_OVERRIDE == 0 )); then
                TIER_FILTERS=(1)
            fi
            if (( CHR_OVERRIDE == 0 )); then
                USE_ALL_CHR=0
                TEST_CHRS=(20 22 MT)
            fi
            ;;
        config-spread)
            if (( CHR_OVERRIDE == 0 )); then
                USE_ALL_CHR=0
                TEST_CHRS=(20 22 MT)
            fi
            ;;
        core-wide)
            if (( CONFIG_OVERRIDE == 0 )); then
                TIER_FILTERS=(1)
            fi
            if (( CHR_OVERRIDE == 0 )); then
                REAL_CHRS_ONLY=1
                USE_ALL_CHR=1
            fi
            ;;
        pairwise)
            load_pairwise_configs
            if (( CONFIG_OVERRIDE == 0 )); then
                TIER_FILTERS=(PW)
            fi
            if (( CHR_OVERRIDE == 0 )); then
                USE_ALL_CHR=0
                TEST_CHRS=(1 2 5 10 14 17 20 22 X MT)
            fi
            ;;
        full-gate)
            FULL_GATE=1
            ;;
    esac
}

append_runner_flags() {
    local -n cmd_ref=$1

    if (( NO_STOP )); then
        cmd_ref+=("--no-stop")
    fi
    if (( NO_CLEANUP )); then
        cmd_ref+=("--no-cleanup")
    fi
    if (( RUST_ONLY )); then
        cmd_ref+=("--rust-only")
    fi
    if (( CLEAN_ALL )); then
        cmd_ref+=("--clean-all")
    fi
    if (( RELEASE_BUILD )); then
        cmd_ref+=("--release")
    fi
    if [[ -n "$RUST_BIN_OVERRIDE" ]]; then
        cmd_ref+=("--rust-bin" "$RUST_BIN_OVERRIDE")
    fi
    if (( NO_BUILD )); then
        cmd_ref+=("--no-build")
    fi
    if [[ -n "$PARALLEL_WORKERS" ]]; then
        cmd_ref+=("--parallel" "$PARALLEL_WORKERS")
    fi
    if [[ -n "$TIMEOUT_SECONDS" ]]; then
        cmd_ref+=("--timeout" "$TIMEOUT_SECONDS")
    fi
    cmd_ref+=("--retry" "$RETRY_COUNT")
    if [[ -n "$MEM_WARN_GB" ]]; then
        cmd_ref+=("--mem-warn-gb" "$MEM_WARN_GB")
    fi
    if [[ -n "$MEM_ABORT_GB" ]]; then
        cmd_ref+=("--mem-abort-gb" "$MEM_ABORT_GB")
    fi
    if [[ -n "$DISK_WARN_GB" ]]; then
        cmd_ref+=("--disk-warn-gb" "$DISK_WARN_GB")
    fi
    # Per-config extra flags (e.g., pileup overrides)
    if (( ${#_PER_CONFIG_EXTRA_FLAGS[@]} > 0 )); then
        cmd_ref+=("${_PER_CONFIG_EXTRA_FLAGS[@]}")
    fi
}

run_parity_for_chr() {
    local opts_label=$1
    local cli_flags=$2
    local chr=$3
    local rc=0
    local -a cmd=("bash" "$PARITY_SCRIPT" "--chr" "$chr" "--opts" "$cli_flags" "--opts-label" "$opts_label")

    append_runner_flags cmd
    "${cmd[@]}" || rc=$?
    return "$rc"
}

run_parity_all_chr() {
    local opts_label=$1
    local cli_flags=$2
    local rc=0
    local -a cmd=("bash" "$PARITY_SCRIPT" "--all-chr" "--opts" "$cli_flags" "--opts-label" "$opts_label")

    append_runner_flags cmd
    "${cmd[@]}" || rc=$?
    return "$rc"
}

record_result() {
    local status=$1
    local config_id=$2
    local opts_label=$3
    local tier=$4
    local scope=$5
    local elapsed_seconds=$6
    local cli_flags=$7

    RESULTS+=("${status}|${config_id}|${opts_label}|${tier}|${scope}|${elapsed_seconds}|${cli_flags}")
}

chr_status_info() {
    local opts_label=$1
    local chr=$2
    local chr_dir="${BASE_DIR}/${opts_label}/${chr}"
    local verified_file="${chr_dir}/.verified"
    local diff_dir="${chr_dir}/diff"
    local verified_at=""
    local stored_mtime=""
    local current_mtime=""

    if [[ -f "$verified_file" ]]; then
        verified_at=$(read_marker_field "$verified_file" "VERIFIED_AT" || true)
        stored_mtime=$(read_marker_field "$verified_file" "RUST_BINARY_MTIME" || true)
        current_mtime=$(current_rust_binary_mtime)

        if [[ -n "$current_mtime" && -n "$stored_mtime" && "$current_mtime" != "$stored_mtime" ]]; then
            printf 'STALE|%s\n' "$verified_at"
        else
            printf 'PASS|%s\n' "$verified_at"
        fi
        return 0
    fi

    if compgen -G "${diff_dir}/shard_*.status" >/dev/null; then
        if grep -q '^FAIL$' "${diff_dir}"/shard_*.status 2>/dev/null; then
            printf 'FAIL|\n'
            return 0
        fi
        printf 'PENDING|\n'
        return 0
    fi

    if compgen -G "${diff_dir}/shard_*.diff" >/dev/null || compgen -G "${diff_dir}/shard_*.meta" >/dev/null; then
        printf 'FAIL|\n'
        return 0
    fi

    if [[ -f "${chr_dir}/results.json" || -d "${chr_dir}/java" || -d "${chr_dir}/rust" || -d "$diff_dir" ]]; then
        printf 'PENDING|\n'
        return 0
    fi

    printf 'PENDING|\n'
}

dry_run_status_text() {
    local opts_label=$1
    local chr=$2
    local info
    local status
    local verified_at
    local verified_date

    info=$(chr_status_info "$opts_label" "$chr")
    IFS='|' read -r status verified_at <<<"$info"
    verified_date=${verified_at%%T*}

    case "$status" in
        PASS)
            printf 'verified (%s)' "${verified_date:-unknown}"
            ;;
        STALE)
            printf 'stale (%s)' "${verified_date:-unknown}"
            ;;
        FAIL)
            printf 'failed'
            ;;
        *)
            printf 'unverified'
            ;;
    esac
}

print_dry_run() {
    local entry
    local config_id
    local opts_label
    local cli_flags
    local tier
    local shard_count
    local status_text
    local -a run_chrs=()
    local -a run_configs=()
    local chr
    local cell_count

    resolve_scope_chrs run_chrs

    for entry in "${TEST_MATRIX[@]}"; do
        IFS='|' read -r config_id opts_label cli_flags tier <<<"$entry"
        should_run "$config_id" "$tier" || continue
        run_configs+=("$config_id")
    done

    (( ${#run_configs[@]} > 0 )) || fail "No configs matched the selected filters"
    cell_count=$((${#run_configs[@]} * ${#run_chrs[@]}))

    printf 'DRY-RUN: %d configs x %d chromosomes = %d cells\n' \
        "${#run_configs[@]}" \
        "${#run_chrs[@]}" \
        "$cell_count"

    for entry in "${TEST_MATRIX[@]}"; do
        IFS='|' read -r config_id opts_label cli_flags tier <<<"$entry"
        should_run "$config_id" "$tier" || continue
        printf 'DRY-RUN: %s (%s) tier=%s flags="%s" chrs="%s"\n' \
            "$config_id" "$opts_label" "$tier" "$cli_flags" "${run_chrs[*]}"

        for chr in "${run_chrs[@]}"; do
            shard_count=$(shard_count_for_chr "$chr")
            status_text=$(dry_run_status_text "$opts_label" "$chr")
            printf '  %s: %s shards, status=%s\n' "$(render_chr_label "$chr")" "$shard_count" "$status_text"
        done
    done
}

print_status_matrix() {
    local entry
    local config_id
    local opts_label
    local cli_flags
    local tier
    local info
    local chr_status
    local verified_at
    local -a run_chrs=()
    local chr
    local matched=0
    local chr_width=8
    local label

    resolve_scope_chrs run_chrs
    for chr in "${run_chrs[@]}"; do
        label=$(render_chr_label "$chr")
        if (( ${#label} + 2 > chr_width )); then
            chr_width=$((${#label} + 2))
        fi
    done

    echo "=== PARITY STATUS MATRIX ==="
    printf '%-10s %-16s %-5s' "CONFIG" "LABEL" "TIER"
    for chr in "${run_chrs[@]}"; do
        printf " %-${chr_width}s" "$(render_chr_label "$chr")"
    done
    echo

    for entry in "${TEST_MATRIX[@]}"; do
        IFS='|' read -r config_id opts_label cli_flags tier <<<"$entry"
        should_run "$config_id" "$tier" || continue
        matched=1
        printf '%-10s %-16s %-5s' "$config_id" "$opts_label" "$tier"
        for chr in "${run_chrs[@]}"; do
            info=$(chr_status_info "$opts_label" "$chr")
            IFS='|' read -r chr_status verified_at <<<"$info"
            printf " %-${chr_width}s" "$chr_status"
        done
        echo
    done

    (( matched > 0 )) || fail "No configs matched the selected filters"
}

print_run_banner() {
    local entry
    local config_id
    local opts_label
    local cli_flags
    local tier
    local config_count
    local chr_count
    local cell_count
    local banner_label
    local banner_value
    local -a run_chrs=()
    local -a run_configs=()

    resolve_scope_chrs run_chrs

    for entry in "${TEST_MATRIX[@]}"; do
        IFS='|' read -r config_id opts_label cli_flags tier <<<"$entry"
        should_run "$config_id" "$tier" || continue
        run_configs+=("$config_id")
    done

    config_count=${#run_configs[@]}
    (( config_count > 0 )) || fail "No configs matched the selected filters"

    chr_count=${#run_chrs[@]}
    cell_count=$((config_count * chr_count))

    if [[ -n "$PRESET" ]]; then
        banner_label="Preset"
        banner_value=$PRESET
    else
        banner_label="Selection"
        banner_value="custom"
    fi

    echo
    printf '══════════════════════════════════════════════════════════════\n'
    printf '  %s: %s | %d configs × %d chromosomes = %d cells\n' \
        "$banner_label" \
        "$banner_value" \
        "$config_count" \
        "$chr_count" \
        "$cell_count"
    printf '  Configs: %s\n' "${run_configs[*]}"
    printf '  Chromosomes: %s\n' "${run_chrs[*]}"
    printf '══════════════════════════════════════════════════════════════\n'
}

json_number_field() {
    local file=$1
    local key=$2

    awk -v key="\"${key}\"" '
        index($0, key) {
            sub(/^[^:]*:[[:space:]]*/, "", $0)
            sub(/,[[:space:]]*$/, "", $0)
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", $0)
            print $0
            exit
        }
    ' "$file"
}

write_option_report() {
    local row
    local status
    local config_id
    local opts_label
    local tier
    local scope
    local elapsed_seconds
    local cli_flags
    local idx=0
    local chr_idx=0
    local results_file
    local value_pass
    local value_fail
    local value_empty
    local value_elapsed
    local -a run_chrs=()
    local chr

    mkdir -p "$BASE_DIR"
    resolve_scope_chrs run_chrs

    {
        echo '{'
        printf '  "timestamp": "%s",\n' "$(json_escape "$(date -Iseconds)")"
        printf '  "total_configs": %d,\n' "$TOTAL_CONFIGS"
        printf '  "pass_configs": %d,\n' "$PASS_CONFIGS"
        printf '  "fail_configs": %d,\n' "$FAIL_CONFIGS"
        echo '  "configs": ['

        for row in "${RESULTS[@]}"; do
            IFS='|' read -r status config_id opts_label tier scope elapsed_seconds cli_flags <<<"$row"
            if (( idx > 0 )); then
                echo ','
            fi
            echo '    {'
            printf '      "config_id": "%s",\n' "$(json_escape "$config_id")"
            printf '      "opts_label": "%s",\n' "$(json_escape "$opts_label")"
            if [[ "$tier" =~ ^[0-9]+$ ]]; then
                printf '      "tier": %d,\n' "$tier"
            else
                printf '      "tier": "%s",\n' "$(json_escape "$tier")"
            fi
            printf '      "cli_flags": "%s",\n' "$(json_escape "$cli_flags")"
            printf '      "status": "%s",\n' "$(json_escape "$status")"
            echo '      "chromosomes": {'

            chr_idx=0
            for chr in "${run_chrs[@]}"; do
                results_file="${BASE_DIR}/${opts_label}/${chr}/results.json"
                [[ -f "$results_file" ]] || continue
                value_pass=$(json_number_field "$results_file" "pass")
                value_fail=$(json_number_field "$results_file" "fail")
                value_empty=$(json_number_field "$results_file" "empty")
                value_elapsed=$(json_number_field "$results_file" "elapsed_seconds")
                if (( chr_idx > 0 )); then
                    echo ','
                fi
                printf '        "%s": {"pass": %s, "fail": %s, "empty": %s, "elapsed_seconds": %s}' \
                    "$(json_escape "$chr")" \
                    "${value_pass:-0}" \
                    "${value_fail:-0}" \
                    "${value_empty:-0}" \
                    "${value_elapsed:-0}"
                chr_idx=$((chr_idx + 1))
            done

            echo
            echo '      }'
            echo -n '    }'
            idx=$((idx + 1))
        done

        echo
        echo '  ]'
        echo '}'
    } >"$REPORT_FILE"
}

extract_line_difference() {
    local java_file=$1
    local rust_file=$2
    local target_line=$3

    awk -F '\t' -v target="$target_line" '
        NR == FNR {
            if (FNR == target) {
                jseen = 1
                jn = NF
                for (i = 1; i <= NF; i++) {
                    jf[i] = $i
                }
            }
            next
        }
        FNR == target {
            rseen = 1
            if (!jseen) {
                printf "row-key||%s:%s %s>%s\n", $2, $4, $6, $7
                exit
            }
            jkey = jf[2] ":" jf[4] " " jf[6] ">" jf[7]
            rkey = $2 ":" $4 " " $6 ">" $7
            if (jkey != rkey) {
                printf "row-key|%s|%s\n", jkey, rkey
                exit
            }
            max_fields = (jn > NF ? jn : NF)
            for (i = 1; i <= max_fields; i++) {
                jv = (i in jf ? jf[i] : "")
                rv = (i <= NF ? $i : "")
                if (jv != rv) {
                    printf "col%d|%s|%s\n", i, jv, rv
                    exit
                }
            }
            printf "line|%s|%s\n", jkey, rkey
            exit
        }
        END {
            if (jseen && !rseen) {
                printf "row-key|%s:%s %s>%s|\n", jf[2], jf[4], jf[6], jf[7]
            }
        }
    ' "$java_file" "$rust_file"
}

first_failure_meta_for_chr() {
    local opts_label=$1
    local chr=$2
    local diff_dir="${BASE_DIR}/${opts_label}/${chr}/diff"
    local meta_file
    local status_file

    shopt -s nullglob
    for meta_file in "${diff_dir}"/shard_*.meta; do
        status_file=${meta_file%.meta}.status
        if [[ ! -f "$status_file" || $(<"$status_file") == "FAIL" ]]; then
            printf '%s\n' "$meta_file"
            shopt -u nullglob
            return 0
        fi
    done
    shopt -u nullglob
    return 1
}

print_failure_digest() {
    local row
    local status
    local config_id
    local opts_label
    local tier
    local scope
    local elapsed_seconds
    local cli_flags
    local -a run_chrs=()
    local chr
    local meta_file
    local REGION=""
    local SHARD=""
    local JAVA_FILE=""
    local RUST_FILE=""
    local DIFF_FILE=""
    local JAVA_LOG=""
    local RUST_LOG=""
    local JAVA_LINES=""
    local RUST_LINES=""
    local FIRST_DIFF_LINE=""
    local SUMMARY=""
    local REASON=""
    local diff_info
    local diff_kind
    local java_value
    local rust_value
    local printed=0
    local chr_label

    resolve_scope_chrs run_chrs

    for row in "${RESULTS[@]}"; do
        IFS='|' read -r status config_id opts_label tier scope elapsed_seconds cli_flags <<<"$row"
        [[ "$status" == "FAIL" ]] || continue
        if (( printed == 0 )); then
            echo
            echo 'FAIL DIGEST:'
            printed=1
        fi

        meta_file=""
        for chr in "${run_chrs[@]}"; do
            if meta_file=$(first_failure_meta_for_chr "$opts_label" "$chr"); then
                break
            fi
            meta_file=""
        done

        if [[ -z "$meta_file" ]]; then
            printf '  %s %s: failure artifacts not found\n' "$config_id" "$opts_label"
            continue
        fi

        REGION=""
        SHARD=""
        JAVA_FILE=""
        RUST_FILE=""
        DIFF_FILE=""
        JAVA_LOG=""
        RUST_LOG=""
        JAVA_LINES=""
        RUST_LINES=""
        FIRST_DIFF_LINE=""
        SUMMARY=""
        REASON=""
        # shellcheck disable=SC1090
        source "$meta_file"
        chr=${REGION%%:*}
        chr_label=$(render_chr_label "$chr")

        if [[ -n "$FIRST_DIFF_LINE" && -f "$JAVA_FILE" && -f "$RUST_FILE" ]]; then
            diff_info=$(extract_line_difference "$JAVA_FILE" "$RUST_FILE" "$FIRST_DIFF_LINE")
            IFS='|' read -r diff_kind java_value rust_value <<<"$diff_info"
            printf '  %s %s %s/shard_%s line %s: %s java="%s" rust="%s"\n' \
                "$config_id" \
                "$opts_label" \
                "$chr_label" \
                "$SHARD" \
                "$FIRST_DIFF_LINE" \
                "$diff_kind" \
                "$(display_escape "$java_value")" \
                "$(display_escape "$rust_value")"
        else
            printf '  %s %s %s/shard_%s: %s\n' \
                "$config_id" \
                "$opts_label" \
                "$chr_label" \
                "$SHARD" \
                "$SUMMARY"
        fi
    done
}

print_summary() {
    local row
    local status
    local config_id
    local opts_label
    local tier
    local scope
    local elapsed_seconds
    local cli_flags

    echo
    echo '=== OPTION PARITY SUMMARY ==='
    echo "Configs: ${TOTAL_CONFIGS} total, ${PASS_CONFIGS} pass, ${FAIL_CONFIGS} fail"
    printf '%-6s  %-8s  %-16s  %-4s  %-12s  %-8s  %s\n' "STATUS" "CONFIG" "LABEL" "TIER" "SCOPE" "ELAPSED" "FLAGS"

    for row in "${RESULTS[@]}"; do
        IFS='|' read -r status config_id opts_label tier scope elapsed_seconds cli_flags <<<"$row"
        printf '%-6s  %-8s  %-16s  %-4s  %-12s  %-8s  %s\n' \
            "$status" \
            "$config_id" \
            "$opts_label" \
            "$tier" \
            "$scope" \
            "$(format_duration_human "$elapsed_seconds")" \
            "$cli_flags"
    done

    print_failure_digest
}

main() {
    local entry
    local config_id
    local opts_label
    local cli_flags
    local tier
    local scope
    local config_rc
    local chr
    local start_ns
    local end_ns
    local elapsed_seconds
    local matched_configs=0
    local sweep_complete=true
    local -a run_chrs=()

    parse_args "$@"
    configure_rust_bin
    apply_preset

    if (( DRY_RUN )) && (( STATUS_ONLY )); then
        fail "--dry-run and --status cannot be used together"
    fi

    if [[ -n "$PARALLEL_WORKERS" ]]; then
        validate_positive_integer "parallel worker count" "$PARALLEL_WORKERS"
    fi
    if [[ -n "$TIMEOUT_SECONDS" ]]; then
        validate_positive_integer "timeout" "$TIMEOUT_SECONDS"
    fi
    validate_nonnegative_integer "retry count" "$RETRY_COUNT"
    if [[ -n "$MEM_WARN_GB" ]]; then
        validate_positive_integer "mem warn threshold" "$MEM_WARN_GB"
    fi
    if [[ -n "$MEM_ABORT_GB" ]]; then
        validate_positive_integer "mem abort threshold" "$MEM_ABORT_GB"
    fi
    if [[ -n "$DISK_WARN_GB" ]]; then
        validate_positive_integer "disk warn threshold" "$DISK_WARN_GB"
    fi

    if (( STATUS_ONLY )); then
        print_status_matrix
        return 0
    fi

    if (( FULL_GATE )); then
        local gate_preset
        local gate_rc
        local -a gate_presets=(smoke tier1 config-spread core-wide)
        local -a gate_cmd=()

        for gate_preset in "${gate_presets[@]}"; do
            gate_rc=0

            echo
            echo '═══════════════════════════════════════════════'
            echo "  FULL-GATE: Running ${gate_preset}"
            echo '═══════════════════════════════════════════════'

            gate_cmd=("bash" "${BASH_SOURCE[0]}" "--preset" "$gate_preset" "--fresh")
            (( RUST_ONLY )) && gate_cmd+=("--rust-only")
            (( NO_STOP )) && gate_cmd+=("--no-stop")
            (( NO_BUILD )) && gate_cmd+=("--no-build")
            (( NO_CLEANUP )) && gate_cmd+=("--no-cleanup")
            (( CLEAN_ALL )) && gate_cmd+=("--clean-all")
            (( RELEASE_BUILD )) && gate_cmd+=("--release")
            (( DRY_RUN )) && gate_cmd+=("--dry-run")
            (( STATUS_ONLY )) && gate_cmd+=("--status")
            (( INCLUDE_BLOCKED )) && gate_cmd+=("--include-blocked")
            [[ -n "$RUST_BIN_OVERRIDE" ]] && gate_cmd+=("--rust-bin" "$RUST_BIN_OVERRIDE")
            [[ -n "$PARALLEL_WORKERS" ]] && gate_cmd+=("--parallel" "$PARALLEL_WORKERS")
            [[ -n "$TIMEOUT_SECONDS" ]] && gate_cmd+=("--timeout" "$TIMEOUT_SECONDS")
            gate_cmd+=("--retry" "$RETRY_COUNT")
            [[ -n "$MEM_WARN_GB" ]] && gate_cmd+=("--mem-warn-gb" "$MEM_WARN_GB")
            [[ -n "$MEM_ABORT_GB" ]] && gate_cmd+=("--mem-abort-gb" "$MEM_ABORT_GB")
            [[ -n "$DISK_WARN_GB" ]] && gate_cmd+=("--disk-warn-gb" "$DISK_WARN_GB")

            "${gate_cmd[@]}" || gate_rc=$?

            if (( gate_rc != 0 )); then
                echo "FULL-GATE: ${gate_preset} FAILED (exit ${gate_rc}). Stopping."
                exit "$gate_rc"
            fi

            echo "FULL-GATE: ${gate_preset} PASSED"
        done

        echo
        echo 'FULL-GATE: All tiers passed.'
        exit 0
    fi

    if (( DRY_RUN )); then
        print_dry_run
        return 0
    fi

    load_checkpoint

    [[ -f "$PARITY_SCRIPT" ]] || fail "Parity runner not found: $PARITY_SCRIPT"
    mkdir -p ./tmp
    print_run_banner
    resolve_scope_chrs run_chrs

    for entry in "${TEST_MATRIX[@]}"; do
        IFS='|' read -r config_id opts_label cli_flags tier <<<"$entry"
        should_run "$config_id" "$tier" || continue
        matched_configs=$((matched_configs + 1))

        if [[ -n "${DONE_MAP[$config_id]+isset}" ]]; then
            echo "Skipping $config_id (checkpoint: already completed)"
            continue
        fi

        # Per-config overrides: pileup runs with bounded Java heap to avoid OOM (23G RAM).
        # Java: 5 workers × 4g = 20g max. Actual usage ~13GB for 3 workers.
        # Rust: 10 workers (phases are sequential, so full RAM available when Rust runs).
        # History: 2×8g proven-safe; 3×7g proven-safe; 5×4g+Rust10 OOM on chr5; 5×4g+Rust5 was old.
        # 2026-04-01: Rust back to 10 per user instruction (Java heap freed before Rust phase).
        _saved_parallel="$PARALLEL_WORKERS"
        _PER_CONFIG_EXTRA_FLAGS=()
        if config_uses_pileup "$cli_flags"; then
            _PER_CONFIG_EXTRA_FLAGS=("--java-parallel" "5" "--java-heap" "4g" "--parallel" "10")
        fi

        TOTAL_CONFIGS=$((TOTAL_CONFIGS + 1))
        scope=$(scope_label)

        echo
        echo '================================================================'
        echo "Config: ${config_id} (${opts_label})"
        echo "Tier:   ${tier}"
        echo "Flags:  ${cli_flags}"
        echo "Scope:  ${scope}"
        echo '================================================================'

        start_ns=$(date +%s%N)
        config_rc=0
        if (( USE_ALL_CHR )) && [[ "$PRESET" != "release" ]]; then
            run_parity_all_chr "$opts_label" "$cli_flags" || config_rc=$?
        else
            for chr in "${run_chrs[@]}"; do
                echo
                echo "--- ${config_id}: chromosome ${chr} ---"
                run_parity_for_chr "$opts_label" "$cli_flags" "$chr" || config_rc=$?
                if (( config_rc != 0 )) && (( NO_STOP == 0 )); then
                    break
                fi
            done
        fi
        end_ns=$(date +%s%N)
        elapsed_seconds=$(format_elapsed_seconds "$start_ns" "$end_ns")

        # Restore per-config overrides
        PARALLEL_WORKERS="$_saved_parallel"
        _PER_CONFIG_EXTRA_FLAGS=()

        if (( config_rc == 0 )); then
            PASS_CONFIGS=$((PASS_CONFIGS + 1))
            record_result "PASS" "$config_id" "$opts_label" "$tier" "$scope" "$elapsed_seconds" "$cli_flags"
        else
            FAIL_CONFIGS=$((FAIL_CONFIGS + 1))
            record_result "FAIL" "$config_id" "$opts_label" "$tier" "$scope" "$elapsed_seconds" "$cli_flags"
        fi

        DONE_MAP["$config_id"]=1
        COMPLETED_CONFIGS+=("$config_id")
        write_checkpoint

        if (( config_rc == 0 )); then
            continue
        fi

        if (( NO_STOP == 0 )); then
            echo
            echo "STOPPING: ${config_id} failed. Re-run with --no-stop to continue past failures."
            sweep_complete=false
            break
        fi
    done

    if (( matched_configs == 0 )); then
        fail "No configs matched the selected filters"
    fi

    write_option_report
    print_summary

    if [[ "$sweep_complete" == true ]]; then
        rm -f "$CHECKPOINT_FILE"
        echo 'Sweep complete — checkpoint removed'
    fi

    if (( FAIL_CONFIGS > 0 )); then
        exit 1
    fi
}

main "$@"