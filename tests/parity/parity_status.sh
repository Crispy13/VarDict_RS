#!/usr/bin/env bash
set -euo pipefail

SCRIPT_NAME=$(basename "${BASH_SOURCE[0]}")
PROJECT_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && git rev-parse --show-toplevel)
cd "$PROJECT_ROOT"

readonly BASE_DIR="./tmp/na12878_parity"
readonly RUST_BIN="target/debug-release/vardict"
readonly REF_FASTA_FAI="testdata/hs37d5.fa.fai"

declare -a DEFAULT_CHRS=(20 22 MT)
declare -a CONFIG_FILTERS=()
declare -a CHR_FILTERS=()
declare -a TIER_FILTERS=()
declare -a SELECTED_CHRS=()
declare -a SELECTED_CONFIGS=()

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

declare -A STATUS_MAP=()
declare -A DIR_BYTES_MAP=()
declare -A CLEANED_MAP=()
declare -A HAS_DATA_MAP=()
declare -A CONFIG_TOTAL_BYTES_MAP=()
declare -A CONFIG_ALL_CLEANED_MAP=()
declare -A RESULT_EXISTS_MAP=()
declare -A RESULT_TOTAL_MAP=()
declare -A RESULT_PASS_MAP=()
declare -A RESULT_FAIL_MAP=()
declare -A RESULT_EMPTY_MAP=()
declare -A RESULT_ELAPSED_MAP=()

SHOW_DISK=0
SHOW_VERIFIED=0
SHOW_FAILURES=0
SHOW_RESOURCES=0
SHOW_JSON_RESULTS=0
SECTIONS_SPECIFIED=0
OUTPUT_JSON=0

BINARY_EXISTS=0
BINARY_MTIME_EPOCH=""
BINARY_MTIME_TEXT="missing"
BASE_TOTAL_BYTES=0

usage() {
    cat <<EOF
Usage: ${SCRIPT_NAME} [options]

Options:
  --config LABEL    Show only this config label (can repeat)
  --chr CHR         Show only this chromosome (can repeat)
  --tier N          Show only this tier (can repeat)
  --disk            Show only disk usage section
  --verified        Show only verified matrix section
  --failures        Show only failing shards section
  --resources       Show only resource section
    --json            Output all sections as JSON to stdout
  --help            Show help
EOF
}

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

array_contains() {
    local needle=$1
    shift || true
    local item

    for item in "$@"; do
        if [[ "$item" == "$needle" ]]; then
            return 0
        fi
    done

    return 1
}

append_unique() {
    local value=$1
    shift
    local -n array_ref=$1

    if ! array_contains "$value" "${array_ref[@]}"; then
        array_ref+=("$value")
    fi
}

normalize_chr() {
    local chr=$1

    if [[ "$chr" == chr* ]]; then
        chr=${chr#chr}
    fi
    printf '%s' "$chr"
}

chr_label() {
    local chr=$1

    if [[ "$chr" == chr* ]]; then
        printf '%s' "$chr"
    else
        printf 'chr%s' "$chr"
    fi
}

human_bytes() {
    local bytes=${1:-0}

    awk -v bytes="$bytes" 'BEGIN {
        split("B K M G T P", units, " ")
        value = bytes + 0
        unit = 1
        while (value >= 1024 && unit < 6) {
            value /= 1024
            unit++
        }
        if (unit == 1) {
            printf "%d%s", value, units[unit]
        } else if (unit == 2 && value >= 10) {
            printf "%.0f%s", value, units[unit]
        } else {
            printf "%.1f%s", value, units[unit]
        }
    }'
}

percent_text() {
    local numerator=${1:-0}
    local denominator=${2:-0}

    awk -v numerator="$numerator" -v denominator="$denominator" 'BEGIN {
        if (denominator <= 0) {
            printf "0.0%%"
        } else {
            printf "%.1f%%", (numerator / denominator) * 100
        }
    }'
}

json_escape() {
    local value=${1-}

    value=${value//\\/\\\\}
    value=${value//\"/\\\"}
    value=${value//$'\n'/\\n}
    value=${value//$'\r'/\\r}
    value=${value//$'\t'/\\t}
    printf '%s' "$value"
}

marker_field() {
    local file=$1
    local key=$2

    if [[ ! -f "$file" ]]; then
        return 0
    fi

    awk -F= -v key="$key" '$1 == key { print substr($0, index($0, "=") + 1); exit }' "$file"
}

json_scalar() {
    local file=$1
    local key=$2

    if [[ ! -f "$file" ]]; then
        return 0
    fi

    awk -v key="\"${key}\":" '
        $1 == key {
            value = $2
            sub(/,$/, "", value)
            gsub(/^"/, "", value)
            gsub(/"$/, "", value)
            print value
            exit
        }
    ' "$file"
}

dir_has_any_data() {
    local dir=$1
    local first_entry

    if [[ ! -d "$dir" ]]; then
        return 1
    fi

    first_entry=$(find "$dir" -mindepth 1 -print -quit 2>/dev/null || true)
    [[ -n "$first_entry" ]]
}

dir_size_bytes() {
    local dir=$1

    if [[ ! -d "$dir" ]]; then
        printf '0'
        return 0
    fi

    du -sb "$dir" 2>/dev/null | awk 'NR == 1 { print $1 }'
}

is_cleaned_verified_dir() {
    local dir=$1

    [[ -f "$dir/.verified" ]] || return 1
    [[ ! -d "$dir/java" ]] || return 1
    [[ ! -d "$dir/rust" ]] || return 1
    [[ ! -d "$dir/diff" ]] || return 1
    return 0
}

has_fail_status_files() {
    local dir=$1
    local status_file

    shopt -s nullglob
    for status_file in "$dir"/diff/*.status; do
        if [[ -f "$status_file" ]] && grep -qx 'FAIL' "$status_file"; then
            shopt -u nullglob
            return 0
        fi
    done
    shopt -u nullglob

    return 1
}

compute_status() {
    local label=$1
    local chr=$2
    local dir="${BASE_DIR}/${label}/${chr}"
    local verified_file="${dir}/.verified"
    local results_file="${dir}/results.json"
    local stored_mtime
    local result_fail

    if [[ -f "$verified_file" ]]; then
        stored_mtime=$(marker_field "$verified_file" "RUST_BINARY_MTIME")
        if (( BINARY_EXISTS == 1 )) && [[ "$stored_mtime" == "$BINARY_MTIME_EPOCH" ]]; then
            printf 'PASS'
        else
            printf 'STALE'
        fi
        return 0
    fi

    if has_fail_status_files "$dir"; then
        printf 'FAIL'
        return 0
    fi

    if [[ -f "$results_file" ]]; then
        result_fail=$(json_scalar "$results_file" "fail")
        if [[ -n "$result_fail" ]] && (( result_fail > 0 )); then
            printf 'FAIL'
            return 0
        fi
    fi

    if dir_has_any_data "$dir"; then
        printf 'PEND'
    else
        printf 'PEND'
    fi
}

parse_args() {
    while (( $# > 0 )); do
        case "$1" in
            --config)
                shift
                [[ $# -gt 0 ]] || fail "Missing value for --config"
                append_unique "$1" CONFIG_FILTERS
                ;;
            --chr)
                shift
                [[ $# -gt 0 ]] || fail "Missing value for --chr"
                append_unique "$(normalize_chr "$1")" CHR_FILTERS
                ;;
            --tier)
                shift
                [[ $# -gt 0 ]] || fail "Missing value for --tier"
                [[ "$1" =~ ^[1-4]$ ]] || fail "Tier must be 1, 2, 3, or 4: $1"
                append_unique "$1" TIER_FILTERS
                ;;
            --disk)
                SHOW_DISK=1
                SECTIONS_SPECIFIED=1
                ;;
            --verified)
                SHOW_VERIFIED=1
                SECTIONS_SPECIFIED=1
                ;;
            --failures)
                SHOW_FAILURES=1
                SECTIONS_SPECIFIED=1
                ;;
            --resources)
                SHOW_RESOURCES=1
                SECTIONS_SPECIFIED=1
                ;;
            --json)
                OUTPUT_JSON=1
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

    if (( SECTIONS_SPECIFIED == 0 )); then
        SHOW_DISK=1
        SHOW_VERIFIED=1
        SHOW_FAILURES=1
        SHOW_RESOURCES=1
        SHOW_JSON_RESULTS=1
    elif (( OUTPUT_JSON == 0 )); then
        SHOW_JSON_RESULTS=0
    fi

    if (( OUTPUT_JSON == 1 )); then
        SHOW_DISK=1
        SHOW_VERIFIED=1
        SHOW_FAILURES=1
        SHOW_RESOURCES=1
        SHOW_JSON_RESULTS=1
    fi
}

select_chromosomes() {
    local chr

    if (( ${#CHR_FILTERS[@]} > 0 )); then
        SELECTED_CHRS=()
        for chr in "${CHR_FILTERS[@]}"; do
            SELECTED_CHRS+=("$chr")
        done
    else
        SELECTED_CHRS=()
        for chr in "${DEFAULT_CHRS[@]}"; do
            SELECTED_CHRS+=("$chr")
        done
    fi
}

select_configs() {
    local entry
    local config_id
    local label
    local cli_flags
    local tier

    SELECTED_CONFIGS=()
    for entry in "${TEST_MATRIX[@]}"; do
        IFS='|' read -r config_id label cli_flags tier <<<"$entry"
        if (( ${#CONFIG_FILTERS[@]} > 0 )) && ! array_contains "$label" "${CONFIG_FILTERS[@]}"; then
            continue
        fi
        if (( ${#TIER_FILTERS[@]} > 0 )) && ! array_contains "$tier" "${TIER_FILTERS[@]}"; then
            continue
        fi
        SELECTED_CONFIGS+=("$entry")
    done
}

load_binary_state() {
    if [[ -f "$RUST_BIN" ]]; then
        BINARY_EXISTS=1
        BINARY_MTIME_EPOCH=$(stat -c %Y "$RUST_BIN")
        BINARY_MTIME_TEXT=$(date -d "@${BINARY_MTIME_EPOCH}" '+%Y-%m-%d %H:%M:%S')
    else
        BINARY_EXISTS=0
        BINARY_MTIME_EPOCH=""
        BINARY_MTIME_TEXT="missing"
    fi
}

gather_state() {
    local entry
    local config_id
    local label
    local cli_flags
    local tier
    local chr
    local key
    local dir
    local bytes
    local status
    local cleaned
    local all_cleaned_verified
    local results_file
    local result_total
    local result_pass
    local result_fail
    local result_empty
    local result_elapsed

    BASE_TOTAL_BYTES=$(dir_size_bytes "$BASE_DIR")

    for entry in "${SELECTED_CONFIGS[@]}"; do
        IFS='|' read -r config_id label cli_flags tier <<<"$entry"
        CONFIG_TOTAL_BYTES_MAP["$label"]=0
        all_cleaned_verified=1

        for chr in "${SELECTED_CHRS[@]}"; do
            key="${label}|${chr}"
            dir="${BASE_DIR}/${label}/${chr}"
            bytes=$(dir_size_bytes "$dir")
            DIR_BYTES_MAP["$key"]=$bytes
            if dir_has_any_data "$dir"; then
                HAS_DATA_MAP["$key"]=1
            else
                HAS_DATA_MAP["$key"]=0
            fi
            if is_cleaned_verified_dir "$dir"; then
                cleaned=1
            else
                cleaned=0
                all_cleaned_verified=0
            fi
            CLEANED_MAP["$key"]=$cleaned
            status=$(compute_status "$label" "$chr")
            STATUS_MAP["$key"]=$status

            if [[ "$status" != "PASS" || "$cleaned" -ne 1 ]]; then
                all_cleaned_verified=0
            fi

            CONFIG_TOTAL_BYTES_MAP["$label"]=$((CONFIG_TOTAL_BYTES_MAP["$label"] + bytes))

            results_file="${dir}/results.json"
            if [[ -f "$results_file" ]]; then
                RESULT_EXISTS_MAP["$key"]=1
                result_total=$(json_scalar "$results_file" "total_shards")
                result_pass=$(json_scalar "$results_file" "pass")
                result_fail=$(json_scalar "$results_file" "fail")
                result_empty=$(json_scalar "$results_file" "empty")
                result_elapsed=$(json_scalar "$results_file" "elapsed_seconds")
                RESULT_TOTAL_MAP["$key"]=${result_total:-0}
                RESULT_PASS_MAP["$key"]=${result_pass:-0}
                RESULT_FAIL_MAP["$key"]=${result_fail:-0}
                RESULT_EMPTY_MAP["$key"]=${result_empty:-0}
                RESULT_ELAPSED_MAP["$key"]=${result_elapsed:-0}
            else
                RESULT_EXISTS_MAP["$key"]=0
                RESULT_TOTAL_MAP["$key"]=0
                RESULT_PASS_MAP["$key"]=0
                RESULT_FAIL_MAP["$key"]=0
                RESULT_EMPTY_MAP["$key"]=0
                RESULT_ELAPSED_MAP["$key"]=0
            fi
        done

        CONFIG_ALL_CLEANED_MAP["$label"]=$all_cleaned_verified
    done
}

config_label_width() {
    local width=6
    local entry
    local config_id
    local label
    local cli_flags
    local tier

    for entry in "${SELECTED_CONFIGS[@]}"; do
        IFS='|' read -r config_id label cli_flags tier <<<"$entry"
        if (( ${#label} > width )); then
            width=${#label}
        fi
    done
    printf '%s' "$width"
}

verified_label_width() {
    local width=5
    local entry
    local config_id
    local label
    local cli_flags
    local tier

    for entry in "${SELECTED_CONFIGS[@]}"; do
        IFS='|' read -r config_id label cli_flags tier <<<"$entry"
        if (( ${#label} > width )); then
            width=${#label}
        fi
    done
    printf '%s' "$width"
}

print_disk_usage() {
    local width
    local entry
    local config_id
    local label
    local cli_flags
    local tier
    local chr
    local key
    local total_display

    echo "=== DISK USAGE ==="
    echo "Total: $(human_bytes "$BASE_TOTAL_BYTES") (${BASE_DIR}/)"
    echo

    width=$(config_label_width)
    printf "%-${width}s" "Config"
    for chr in "${SELECTED_CHRS[@]}"; do
        printf "  %-10s" "$(chr_label "$chr")"
    done
    printf "  %-20s\n" "Total"

    for entry in "${SELECTED_CONFIGS[@]}"; do
        IFS='|' read -r config_id label cli_flags tier <<<"$entry"
        printf "%-${width}s" "$label"
        for chr in "${SELECTED_CHRS[@]}"; do
            key="${label}|${chr}"
            if [[ ${CLEANED_MAP["$key"]:-0} -eq 1 ]] || [[ ${DIR_BYTES_MAP["$key"]:-0} -eq 0 ]]; then
                printf "  %-10s" "-"
            else
                printf "  %-10s" "$(human_bytes "${DIR_BYTES_MAP["$key"]}")"
            fi
        done

        if [[ ${CONFIG_ALL_CLEANED_MAP["$label"]:-0} -eq 1 ]]; then
            total_display="(verified, cleaned)"
        elif [[ ${CONFIG_TOTAL_BYTES_MAP["$label"]:-0} -eq 0 ]]; then
            total_display="-"
        else
            total_display=$(human_bytes "${CONFIG_TOTAL_BYTES_MAP["$label"]}")
        fi
        printf "  %-20s\n" "$total_display"
    done
}

print_verified_status() {
    local label_width
    local entry
    local config_id
    local label
    local cli_flags
    local tier
    local chr
    local key

    echo "=== VERIFIED STATUS ==="
    echo "Rust binary: ${RUST_BIN} (mtime: ${BINARY_MTIME_TEXT})"
    if [[ -f "$REF_FASTA_FAI" ]]; then
        echo "Reference index: ${REF_FASTA_FAI}"
    fi
    echo

    label_width=$(verified_label_width)
    printf '%-8s  ' "CONFIG"
    printf "%-${label_width}s  " "LABEL"
    printf '%-4s' "TIER"
    for chr in "${SELECTED_CHRS[@]}"; do
        printf '  %-8s' "$(chr_label "$chr")"
    done
    printf '\n'

    for entry in "${SELECTED_CONFIGS[@]}"; do
        IFS='|' read -r config_id label cli_flags tier <<<"$entry"
        printf '%-8s  ' "$config_id"
        printf "%-${label_width}s  " "$label"
        printf '%-4s' "$tier"
        for chr in "${SELECTED_CHRS[@]}"; do
            key="${label}|${chr}"
            printf '  %-8s' "${STATUS_MAP["$key"]}"
        done
        printf '\n'
    done
}

meta_field() {
    local file=$1
    local key=$2

    if [[ ! -f "$file" ]]; then
        return 0
    fi

    awk -F= -v key="$key" '$1 == key { print substr($0, index($0, "=") + 1); exit }' "$file"
}

print_failing_shards() {
    local found=0
    local entry
    local config_id
    local label
    local cli_flags
    local tier
    local chr
    local key
    local status_file
    local shard_name
    local meta_file
    local first_diff_line
    local reason
    local java_lines
    local rust_lines
    local detail

    echo "=== FAILING SHARDS ==="

    for entry in "${SELECTED_CONFIGS[@]}"; do
        IFS='|' read -r config_id label cli_flags tier <<<"$entry"
        for chr in "${SELECTED_CHRS[@]}"; do
            key="${label}|${chr}"
            if [[ "${STATUS_MAP["$key"]}" != "FAIL" ]]; then
                continue
            fi

            found=1
            echo "${label}/$(chr_label "$chr"):"
            shopt -s nullglob
            for status_file in "${BASE_DIR}/${label}/${chr}/diff"/*.status; do
                if ! grep -qx 'FAIL' "$status_file"; then
                    continue
                fi
                shard_name=$(basename "$status_file" .status)
                meta_file="${BASE_DIR}/${label}/${chr}/diff/${shard_name}.meta"
                first_diff_line=$(meta_field "$meta_file" "FIRST_DIFF_LINE")
                reason=$(meta_field "$meta_file" "REASON")
                java_lines=$(meta_field "$meta_file" "JAVA_LINES")
                rust_lines=$(meta_field "$meta_file" "RUST_LINES")
                detail=""
                if [[ -n "$java_lines" || -n "$rust_lines" ]]; then
                    detail=" (java=${java_lines:-?} lines, rust=${rust_lines:-?} lines)"
                fi
                printf '  %s: line %s, %s%s\n' "$shard_name" "${first_diff_line:-?}" "${reason:-unknown}" "$detail"
            done
            shopt -u nullglob
            echo
        done
    done

    if (( found == 0 )); then
        echo "No failing shards found."
        echo
    fi
}

print_resources() {
    local mem_available_kb=0
    local mem_total_kb=0
    local mem_available_bytes=0
    local mem_total_bytes=0
    local disk_available_bytes=0

    if [[ -r /proc/meminfo ]]; then
        mem_available_kb=$(awk '$1 == "MemAvailable:" { print $2; exit }' /proc/meminfo)
        mem_total_kb=$(awk '$1 == "MemTotal:" { print $2; exit }' /proc/meminfo)
    fi
    mem_available_bytes=$((mem_available_kb * 1024))
    mem_total_bytes=$((mem_total_kb * 1024))
    disk_available_bytes=$(df -B1 . | awk 'NR == 2 { print $4 }')

    echo "=== SYSTEM RESOURCES ==="
    echo "Memory:    Available $(human_bytes "$mem_available_bytes") / Total $(human_bytes "$mem_total_bytes") ($(percent_text "$mem_available_bytes" "$mem_total_bytes"))"
    echo "Disk:      Available $(human_bytes "$disk_available_bytes") (${BASE_DIR}/ uses $(human_bytes "$BASE_TOTAL_BYTES"))"
}

print_json_results_summary() {
    local found=0
    local entry
    local config_id
    local label
    local cli_flags
    local tier
    local chr
    local key

    echo "=== JSON RESULTS (from v2 runs) ==="
    for entry in "${SELECTED_CONFIGS[@]}"; do
        IFS='|' read -r config_id label cli_flags tier <<<"$entry"
        for chr in "${SELECTED_CHRS[@]}"; do
            key="${label}|${chr}"
            if [[ ${RESULT_EXISTS_MAP["$key"]:-0} -ne 1 ]]; then
                continue
            fi
            found=1
            printf '%s/%s: %s shards, %s pass, %s fail, %s empty (%ss)\n' \
                "$label" \
                "$(chr_label "$chr")" \
                "${RESULT_TOTAL_MAP["$key"]}" \
                "${RESULT_PASS_MAP["$key"]}" \
                "${RESULT_FAIL_MAP["$key"]}" \
                "${RESULT_EMPTY_MAP["$key"]}" \
                "${RESULT_ELAPSED_MAP["$key"]}"
        done
    done

    if (( found == 0 )); then
        echo "No results.json files found."
    fi
}

print_json_array() {
    local first=1
    local value

    printf '['
    for value in "$@"; do
        if (( first == 0 )); then
            printf ', '
        fi
        printf '"%s"' "$(json_escape "$value")"
        first=0
    done
    printf ']'
}

print_json_disk_section() {
    local first_config=1
    local first_chr
    local entry
    local config_id
    local label
    local cli_flags
    local tier
    local chr
    local key
    local total_display

    printf '  "disk_usage": {\n'
    printf '    "base_dir": "%s",\n' "$(json_escape "$BASE_DIR")"
    printf '    "total_bytes": %s,\n' "$BASE_TOTAL_BYTES"
    printf '    "total_human": "%s",\n' "$(json_escape "$(human_bytes "$BASE_TOTAL_BYTES")")"
    printf '    "configs": [\n'

    for entry in "${SELECTED_CONFIGS[@]}"; do
        IFS='|' read -r config_id label cli_flags tier <<<"$entry"
        if (( first_config == 0 )); then
            printf ',\n'
        fi
        if [[ ${CONFIG_ALL_CLEANED_MAP["$label"]:-0} -eq 1 ]]; then
            total_display="(verified, cleaned)"
        elif [[ ${CONFIG_TOTAL_BYTES_MAP["$label"]:-0} -eq 0 ]]; then
            total_display="-"
        else
            total_display=$(human_bytes "${CONFIG_TOTAL_BYTES_MAP["$label"]}")
        fi
        printf '      {"config_id": "%s", "label": "%s", "tier": %s, "total_bytes": %s, "total_human": "%s", "all_cleaned_verified": %s, "chromosomes": [' \
            "$(json_escape "$config_id")" \
            "$(json_escape "$label")" \
            "$tier" \
            "${CONFIG_TOTAL_BYTES_MAP["$label"]}" \
            "$(json_escape "$total_display")" \
            "$( [[ ${CONFIG_ALL_CLEANED_MAP["$label"]:-0} -eq 1 ]] && printf true || printf false )"
        first_chr=1
        for chr in "${SELECTED_CHRS[@]}"; do
            key="${label}|${chr}"
            if (( first_chr == 0 )); then
                printf ', '
            fi
            printf '{"chromosome": "%s", "display": "%s", "bytes": %s, "human": "%s", "cleaned_verified": %s, "status": "%s"}' \
                "$(json_escape "$chr")" \
                "$(json_escape "$(chr_label "$chr")")" \
                "${DIR_BYTES_MAP["$key"]}" \
                "$(json_escape "$(human_bytes "${DIR_BYTES_MAP["$key"]}")")" \
                "$( [[ ${CLEANED_MAP["$key"]:-0} -eq 1 ]] && printf true || printf false )" \
                "$(json_escape "${STATUS_MAP["$key"]}")"
            first_chr=0
        done
        printf ']}'
        first_config=0
    done
    printf '\n    ]\n'
    printf '  }'
}

print_json_verified_section() {
    local first_config=1
    local first_chr
    local entry
    local config_id
    local label
    local cli_flags
    local tier
    local chr
    local key

    printf '  "verified_status": {\n'
    printf '    "rust_binary": {"path": "%s", "exists": %s, "mtime_epoch": %s, "mtime": "%s"},\n' \
        "$(json_escape "$RUST_BIN")" \
        "$( (( BINARY_EXISTS == 1 )) && printf true || printf false )" \
        "${BINARY_MTIME_EPOCH:-null}" \
        "$(json_escape "$BINARY_MTIME_TEXT")"
    printf '    "reference_fai": "%s",\n' "$(json_escape "$REF_FASTA_FAI")"
    printf '    "rows": [\n'
    for entry in "${SELECTED_CONFIGS[@]}"; do
        IFS='|' read -r config_id label cli_flags tier <<<"$entry"
        if (( first_config == 0 )); then
            printf ',\n'
        fi
        printf '      {"config_id": "%s", "label": "%s", "tier": %s, "chromosomes": [' \
            "$(json_escape "$config_id")" \
            "$(json_escape "$label")" \
            "$tier"
        first_chr=1
        for chr in "${SELECTED_CHRS[@]}"; do
            key="${label}|${chr}"
            if (( first_chr == 0 )); then
                printf ', '
            fi
            printf '{"chromosome": "%s", "display": "%s", "status": "%s"}' \
                "$(json_escape "$chr")" \
                "$(json_escape "$(chr_label "$chr")")" \
                "$(json_escape "${STATUS_MAP["$key"]}")"
            first_chr=0
        done
        printf ']}'
        first_config=0
    done
    printf '\n    ]\n'
    printf '  }'
}

print_json_failures_section() {
    local first_group=1
    local first_shard
    local entry
    local config_id
    local label
    local cli_flags
    local tier
    local chr
    local key
    local status_file
    local shard_name
    local meta_file
    local first_diff_line
    local reason
    local java_lines
    local rust_lines

    printf '  "failing_shards": [\n'
    for entry in "${SELECTED_CONFIGS[@]}"; do
        IFS='|' read -r config_id label cli_flags tier <<<"$entry"
        for chr in "${SELECTED_CHRS[@]}"; do
            key="${label}|${chr}"
            if [[ "${STATUS_MAP["$key"]}" != "FAIL" ]]; then
                continue
            fi
            if (( first_group == 0 )); then
                printf ',\n'
            fi
            printf '    {"config_id": "%s", "label": "%s", "tier": %s, "chromosome": "%s", "display": "%s", "shards": [' \
                "$(json_escape "$config_id")" \
                "$(json_escape "$label")" \
                "$tier" \
                "$(json_escape "$chr")" \
                "$(json_escape "$(chr_label "$chr")")"
            first_shard=1
            shopt -s nullglob
            for status_file in "${BASE_DIR}/${label}/${chr}/diff"/*.status; do
                if ! grep -qx 'FAIL' "$status_file"; then
                    continue
                fi
                shard_name=$(basename "$status_file" .status)
                meta_file="${BASE_DIR}/${label}/${chr}/diff/${shard_name}.meta"
                first_diff_line=$(meta_field "$meta_file" "FIRST_DIFF_LINE")
                reason=$(meta_field "$meta_file" "REASON")
                java_lines=$(meta_field "$meta_file" "JAVA_LINES")
                rust_lines=$(meta_field "$meta_file" "RUST_LINES")
                if (( first_shard == 0 )); then
                    printf ', '
                fi
                printf '{"shard": "%s", "first_diff_line": %s, "reason": "%s", "java_lines": %s, "rust_lines": %s}' \
                    "$(json_escape "$shard_name")" \
                    "${first_diff_line:-null}" \
                    "$(json_escape "${reason:-unknown}")" \
                    "${java_lines:-null}" \
                    "${rust_lines:-null}"
                first_shard=0
            done
            shopt -u nullglob
            printf ']}'
            first_group=0
        done
    done
    printf '\n  ]'
}

print_json_resources_section() {
    local mem_available_kb=0
    local mem_total_kb=0
    local mem_available_bytes=0
    local mem_total_bytes=0
    local disk_available_bytes=0

    if [[ -r /proc/meminfo ]]; then
        mem_available_kb=$(awk '$1 == "MemAvailable:" { print $2; exit }' /proc/meminfo)
        mem_total_kb=$(awk '$1 == "MemTotal:" { print $2; exit }' /proc/meminfo)
    fi
    mem_available_bytes=$((mem_available_kb * 1024))
    mem_total_bytes=$((mem_total_kb * 1024))
    disk_available_bytes=$(df -B1 . | awk 'NR == 2 { print $4 }')

    printf '  "resources": {'
    printf '"memory_available_bytes": %s, ' "$mem_available_bytes"
    printf '"memory_total_bytes": %s, ' "$mem_total_bytes"
    printf '"memory_available_human": "%s", ' "$(json_escape "$(human_bytes "$mem_available_bytes")")"
    printf '"memory_total_human": "%s", ' "$(json_escape "$(human_bytes "$mem_total_bytes")")"
    printf '"memory_available_percent": "%s", ' "$(json_escape "$(percent_text "$mem_available_bytes" "$mem_total_bytes")")"
    printf '"disk_available_bytes": %s, ' "$disk_available_bytes"
    printf '"disk_available_human": "%s", ' "$(json_escape "$(human_bytes "$disk_available_bytes")")"
    printf '"parity_dir_bytes": %s, ' "$BASE_TOTAL_BYTES"
    printf '"parity_dir_human": "%s"' "$(json_escape "$(human_bytes "$BASE_TOTAL_BYTES")")"
    printf '}'
}

print_json_results_section() {
    local first_result=1
    local entry
    local config_id
    local label
    local cli_flags
    local tier
    local chr
    local key

    printf '  "json_results": [\n'
    for entry in "${SELECTED_CONFIGS[@]}"; do
        IFS='|' read -r config_id label cli_flags tier <<<"$entry"
        for chr in "${SELECTED_CHRS[@]}"; do
            key="${label}|${chr}"
            if [[ ${RESULT_EXISTS_MAP["$key"]:-0} -ne 1 ]]; then
                continue
            fi
            if (( first_result == 0 )); then
                printf ',\n'
            fi
            printf '    {"config_id": "%s", "label": "%s", "tier": %s, "chromosome": "%s", "display": "%s", "total_shards": %s, "pass": %s, "fail": %s, "empty": %s, "elapsed_seconds": %s}' \
                "$(json_escape "$config_id")" \
                "$(json_escape "$label")" \
                "$tier" \
                "$(json_escape "$chr")" \
                "$(json_escape "$(chr_label "$chr")")" \
                "${RESULT_TOTAL_MAP["$key"]}" \
                "${RESULT_PASS_MAP["$key"]}" \
                "${RESULT_FAIL_MAP["$key"]}" \
                "${RESULT_EMPTY_MAP["$key"]}" \
                "${RESULT_ELAPSED_MAP["$key"]}"
            first_result=0
        done
    done
    printf '\n  ]'
}

print_json_output() {
    local first_section=1

    printf '{\n'
    printf '  "filters": {"configs": '
    print_json_array "${CONFIG_FILTERS[@]}"
    printf ', "chromosomes": '
    print_json_array "${SELECTED_CHRS[@]}"
    printf ', "tiers": '
    print_json_array "${TIER_FILTERS[@]}"
    printf '},\n'

    if (( SHOW_DISK )); then
        print_json_disk_section
        first_section=0
    fi
    if (( SHOW_VERIFIED )); then
        if (( first_section == 0 )); then
            printf ',\n'
        fi
        print_json_verified_section
        first_section=0
    fi
    if (( SHOW_FAILURES )); then
        if (( first_section == 0 )); then
            printf ',\n'
        fi
        print_json_failures_section
        first_section=0
    fi
    if (( SHOW_RESOURCES )); then
        if (( first_section == 0 )); then
            printf ',\n'
        fi
        print_json_resources_section
        first_section=0
    fi
    if (( SHOW_JSON_RESULTS )); then
        if (( first_section == 0 )); then
            printf ',\n'
        fi
        print_json_results_section
    fi
    printf '\n}\n'
}

main() {
    parse_args "$@"
    select_chromosomes
    select_configs
    load_binary_state
    gather_state

    if (( OUTPUT_JSON == 1 )); then
        print_json_output
        return 0
    fi

    if (( SHOW_DISK )); then
        print_disk_usage
        echo
    fi
    if (( SHOW_VERIFIED )); then
        print_verified_status
        echo
    fi
    if (( SHOW_FAILURES )); then
        print_failing_shards
    fi
    if (( SHOW_RESOURCES )); then
        print_resources
        echo
    fi
    if (( SHOW_JSON_RESULTS )); then
        print_json_results_summary
    fi
}

main "$@"