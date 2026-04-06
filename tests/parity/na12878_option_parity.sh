#!/usr/bin/env bash
set -euo pipefail

SCRIPT_NAME=$(basename "${BASH_SOURCE[0]}")
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PARITY_SCRIPT="${SCRIPT_DIR}/na12878_parity.sh"

declare -a DEFAULT_CHRS=(20 22 MT)
declare -a TEST_CHRS=("${DEFAULT_CHRS[@]}")
declare -a TIER_FILTERS=()
declare -a CONFIG_FILTERS=()
declare -a TEST_MATRIX=(
    "T1-01|pileup|-p|1"
    "T1-02|pileup-max|-p -f 0.0 -r 1|1"
    "T1-03|nosv|-U|1"
    "T1-04|nosv-dedup|-U --deldupvar|1"
    "T1-05|freq-low|-f 0.001|1"
    "T1-06|freq-high|-f 0.05|1"
    "T1-07|freq-zero|-f 0.0|1"
    "T1-08|minr-1|-r 1|1"
    "T1-09|minr-5|-r 5|1"
    "T1-10|no-realign|-k 0|1"
    "T1-11|mapq-10|-Q 10|1"
    "T1-12|mapq-30|-Q 30|1"
    "T1-13|fisher|--fisher|1"
    "T1-14|debug|-D|1"
    "T2-01|qual-0|-q 0|2"
    "T2-02|qual-10|-q 10|2"
    "T2-03|qual-30|-q 30|2"
    "T2-04|qual-40|-q 40|2"
    "T2-05|vext-0|-X 0|2"
    "T2-06|vext-5|-X 5|2"
    "T2-07|mismatch-3|-m 3|2"
    "T2-08|mismatch-15|-m 15|2"
    "T2-09|filter-0x700|-F 0x700|2"
    "T2-10|filter-0x100|-F 0x100|2"
    "T2-11|indel-3prime|-3|2"
    "T2-12|readpos-0|-P 0|2"
    "T2-13|readpos-10|-P 10|2"
    "T2-14|qratio-low|-o 0.5|2"
    "T2-15|qratio-high|-o 2.0|2"
    "T2-16|chimeric|--chimeric|2"
    "T3-01|clinical-wgs|-f 0.001 -Q 10 -F 0x700|3"
    "T3-02|pileup-strict|-p -q 30|3"
    "T3-03|no-realign-lowfreq|-k 0 -f 0.001|3"
    "T3-04|fisher-lowfreq|-f 0.001 --fisher|3"
    "T3-05|nosv-tight|-U -f 0.001 -Q 10|3"
    "T3-06|3prime-lowfreq|-3 -f 0.001|3"
    "T3-07|pileup-relaxed|-p -r 1 -q 10 -Q 0|3"
    "T3-08|quality-gauntlet|-M 25 -m 5 -Q 10|3"
    "T4-01|extend-150|-x 150|4"
    "T4-02|refext-600|-Y 600|4"
    "T4-03|refext-2000|-Y 2000|4"
    "T4-04|inssize-small|-w 200 -W 50|4"
    "T4-05|trim-130|-T 130|4"
    "T4-06|include-n-pileup|-K -p -r 1 -f 0.0|4"
)
declare -a RESULTS=()

USE_ALL_CHR=0
CUSTOM_CHRS=0
NO_STOP=0
NO_CLEANUP=0
PARALLEL_WORKERS=""
TOTAL_CONFIGS=0
PASS_CONFIGS=0
FAIL_CONFIGS=0

usage() {
    cat <<EOF
Usage: ${SCRIPT_NAME} [options]

Options:
  --tier N          Run only tier N (1-4). Can be repeated.
  --config-id ID    Run only a specific config ID. Can be repeated.
  --chr CHR         Override chromosomes (default: 20 22 MT). Can be repeated.
  --all-chr         Run each selected config across all contigs from the FAI index.
  --no-stop         Continue to later configs after a failure.
  --no-cleanup      Keep all output files after comparison (disable disk cleanup).
  --parallel N      Override parallel worker count (default: 10).
  --help            Show this help.

Matrix:
  Tier 1: 14 configs
  Tier 2: 16 configs
  Tier 3: 8 configs
  Tier 4: 6 configs
  Total:  44 configs
EOF
}

fail() {
    echo "ERROR: $*" >&2
    exit 1
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

parse_args() {
    while (( $# > 0 )); do
        case "$1" in
            --tier)
                shift
                [[ $# -gt 0 ]] || fail "Missing value for --tier"
                [[ "$1" =~ ^[1-4]$ ]] || fail "Tier must be 1, 2, 3, or 4: $1"
                TIER_FILTERS+=("$1")
                ;;
            --config-id)
                shift
                [[ $# -gt 0 ]] || fail "Missing value for --config-id"
                CONFIG_FILTERS+=("$1")
                ;;
            --chr)
                shift
                [[ $# -gt 0 ]] || fail "Missing value for --chr"
                if (( CUSTOM_CHRS == 0 )); then
                    TEST_CHRS=()
                    CUSTOM_CHRS=1
                fi
                USE_ALL_CHR=0
                TEST_CHRS+=("$1")
                ;;
            --all-chr)
                USE_ALL_CHR=1
                ;;
            --parallel)
                shift
                [[ $# -gt 0 ]] || fail "Missing value for --parallel"
                PARALLEL_WORKERS="$1"
                ;;
            --no-stop)
                NO_STOP=1
                ;;
            --no-cleanup)
                NO_CLEANUP=1
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

should_run() {
    local config_id=$1
    local tier=$2

    if (( ${#TIER_FILTERS[@]} > 0 )) && ! array_contains "$tier" "${TIER_FILTERS[@]}"; then
        return 1
    fi

    if (( ${#CONFIG_FILTERS[@]} > 0 )) && ! array_contains "$config_id" "${CONFIG_FILTERS[@]}"; then
        return 1
    fi

    return 0
}

run_parity_for_chr() {
    local opts_label=$1
    local cli_flags=$2
    local chr=$3
    local rc=0
    local -a cmd=("bash" "$PARITY_SCRIPT" "--chr" "$chr" "--opts" "$cli_flags" "--opts-label" "$opts_label")

    if (( NO_STOP )); then
        cmd+=("--no-stop")
    fi
    if (( NO_CLEANUP )); then
        cmd+=("--no-cleanup")
    fi
    if [[ -n "$PARALLEL_WORKERS" ]]; then
        cmd+=("--parallel" "$PARALLEL_WORKERS")
    fi

    "${cmd[@]}" || rc=$?
    return "$rc"
}

run_parity_all_chr() {
    local opts_label=$1
    local cli_flags=$2
    local rc=0
    local -a cmd=("bash" "$PARITY_SCRIPT" "--all-chr" "--opts" "$cli_flags" "--opts-label" "$opts_label")

    if (( NO_STOP )); then
        cmd+=("--no-stop")
    fi
    if (( NO_CLEANUP )); then
        cmd+=("--no-cleanup")
    fi
    if [[ -n "$PARALLEL_WORKERS" ]]; then
        cmd+=("--parallel" "$PARALLEL_WORKERS")
    fi

    "${cmd[@]}" || rc=$?
    return "$rc"
}

record_result() {
    local status=$1
    local config_id=$2
    local opts_label=$3
    local tier=$4
    local scope=$5
    local cli_flags=$6

    RESULTS+=("${status}|${config_id}|${opts_label}|${tier}|${scope}|${cli_flags}")
}

print_summary() {
    local row
    local status
    local config_id
    local opts_label
    local tier
    local scope
    local cli_flags

    echo
    echo "=== OPTION PARITY SUMMARY ==="
    echo "Configs: ${TOTAL_CONFIGS} total, ${PASS_CONFIGS} pass, ${FAIL_CONFIGS} fail"
    printf '%-6s  %-8s  %-16s  %-4s  %-12s  %s\n' "STATUS" "CONFIG" "LABEL" "TIER" "SCOPE" "FLAGS"

    for row in "${RESULTS[@]}"; do
        IFS='|' read -r status config_id opts_label tier scope cli_flags <<<"$row"
        printf '%-6s  %-8s  %-16s  %-4s  %-12s  %s\n' "$status" "$config_id" "$opts_label" "$tier" "$scope" "$cli_flags"
    done
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

    parse_args "$@"

    [[ -f "$PARITY_SCRIPT" ]] || fail "Parity runner not found: $PARITY_SCRIPT"

    for entry in "${TEST_MATRIX[@]}"; do
        IFS='|' read -r config_id opts_label cli_flags tier <<<"$entry"
        should_run "$config_id" "$tier" || continue

        TOTAL_CONFIGS=$((TOTAL_CONFIGS + 1))
        scope="${TEST_CHRS[*]}"
        if (( USE_ALL_CHR )); then
            scope="all-chr"
        fi

        echo
        echo "================================================================"
        echo "Config: ${config_id} (${opts_label})"
        echo "Tier:   ${tier}"
        echo "Flags:  ${cli_flags}"
        echo "Scope:  ${scope}"
        echo "================================================================"

        config_rc=0
        if (( USE_ALL_CHR )); then
            run_parity_all_chr "$opts_label" "$cli_flags" || config_rc=$?
        else
            for chr in "${TEST_CHRS[@]}"; do
                echo
                echo "--- ${config_id}: chromosome ${chr} ---"
                run_parity_for_chr "$opts_label" "$cli_flags" "$chr" || config_rc=$?
                if (( config_rc != 0 )); then
                    break
                fi
            done
        fi

        if (( config_rc == 0 )); then
            PASS_CONFIGS=$((PASS_CONFIGS + 1))
            record_result "PASS" "$config_id" "$opts_label" "$tier" "$scope" "$cli_flags"
            continue
        fi

        FAIL_CONFIGS=$((FAIL_CONFIGS + 1))
        record_result "FAIL" "$config_id" "$opts_label" "$tier" "$scope" "$cli_flags"

        if (( NO_STOP == 0 )); then
            echo
            echo "STOPPING: ${config_id} failed. Re-run with --no-stop to continue past failures."
            break
        fi
    done

    if (( TOTAL_CONFIGS == 0 )); then
        fail "No configs matched the selected filters"
    fi

    print_summary

    if (( FAIL_CONFIGS > 0 )); then
        exit 1
    fi
}

main "$@"
