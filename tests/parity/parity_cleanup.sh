#!/usr/bin/env bash
set -euo pipefail

readonly SCRIPT_NAME=$(basename "${BASH_SOURCE[0]}")
readonly SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
readonly REPO_ROOT=$(cd "$SCRIPT_DIR" && git rev-parse --show-toplevel)

cd "$REPO_ROOT"

readonly BASE_DIR="tmp/na12878_parity"
readonly JAVA_JAR="VarDictJava/build/libs/VarDict-1.8.3.jar"
readonly RUST_BIN="target/debug-release/vardict"

OPERATION=""
CONFIRM=0

declare -a CONFIG_FILTERS=()
declare -a CHR_FILTERS=()
declare -a TARGETS=()
declare -A TARGET_SET=()
declare -A GROUP_BYTES=()
declare -A GROUP_COUNTS=()

BASE_DIR_ABS=""

usage() {
    cat <<EOF
Usage: ${SCRIPT_NAME} <operation> [options]

Operations (exactly one required):
  --stale-java      Remove Java cache older than JAR
  --stale-rust      Remove Rust output older than binary
  --failed-only     Remove only failing shard artifacts
  --verified-only   Remove artifacts for verified configs
  --all             Remove everything under tmp/na12878_parity/

Scope (optional, repeatable):
  --config LABEL    Only affect this config label
  --chr CHR         Only affect this chromosome

Safety:
  --dry-run         List what would be deleted (DEFAULT)
  --confirm         Actually delete files
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
    local -n target_array=$2

    if ! array_contains "$value" "${target_array[@]}"; then
        target_array+=("$value")
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

marker_field() {
    local file=$1
    local key=$2

    [[ -f "$file" ]] || return 1
    awk -F= -v key="$key" '$1 == key { print substr($0, index($0, "=") + 1); exit }' "$file"
}

entry_mtime() {
    local path=$1

    stat -c %Y -- "$path" 2>/dev/null
}

entry_size() {
    local path=$1

    stat -c %s -- "$path" 2>/dev/null || printf '0\n'
}

matches_scope() {
    local label=$1
    local chr=$2

    if (( ${#CONFIG_FILTERS[@]} > 0 )) && ! array_contains "$label" "${CONFIG_FILTERS[@]}"; then
        return 1
    fi

    if (( ${#CHR_FILTERS[@]} > 0 )) && ! array_contains "$chr" "${CHR_FILTERS[@]}"; then
        return 1
    fi

    return 0
}

require_reference_file() {
    local path=$1
    local label=$2

    [[ -f "$path" ]] || fail "$label not found: $path"
}

validate_operation() {
    local count=0

    [[ -n "$OPERATION" ]] && count=1
    (( count == 1 )) || fail "Exactly one operation is required"
}

parse_args() {
    while (( $# > 0 )); do
        case "$1" in
            --stale-java)
                [[ -z "$OPERATION" ]] || fail "Exactly one operation is required"
                OPERATION="stale-java"
                ;;
            --stale-rust)
                [[ -z "$OPERATION" ]] || fail "Exactly one operation is required"
                OPERATION="stale-rust"
                ;;
            --failed-only)
                [[ -z "$OPERATION" ]] || fail "Exactly one operation is required"
                OPERATION="failed-only"
                ;;
            --verified-only)
                [[ -z "$OPERATION" ]] || fail "Exactly one operation is required"
                OPERATION="verified-only"
                ;;
            --all)
                [[ -z "$OPERATION" ]] || fail "Exactly one operation is required"
                OPERATION="all"
                ;;
            --config)
                shift
                [[ $# -gt 0 ]] || fail "Missing value for --config"
                append_unique "$1" CONFIG_FILTERS
                ;;
            --chr)
                shift
                [[ $# -gt 0 ]] || fail "Missing value for --chr"
                append_unique "$1" CHR_FILTERS
                ;;
            --dry-run)
                CONFIRM=0
                ;;
            --confirm)
                CONFIRM=1
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

    validate_operation
}

resolve_base_dir() {
    if [[ ! -d "$BASE_DIR" ]]; then
        echo "No parity data found"
        exit 0
    fi

    BASE_DIR_ABS=$(cd "$BASE_DIR" && pwd -P)
}

assert_safe_target() {
    local path=$1
    local path_abs

    path_abs=$(cd "$(dirname "$path")" && pwd -P)/$(basename "$path")
    case "$path_abs" in
        "$BASE_DIR_ABS"/*) ;;
        *) fail "Refusing to delete outside $BASE_DIR: $path" ;;
    esac
}

add_target() {
    local path=$1

    [[ -e "$path" || -L "$path" ]] || return 0
    assert_safe_target "$path"

    if [[ -z "${TARGET_SET["$path"]+x}" ]]; then
        TARGET_SET["$path"]=1
        TARGETS+=("$path")
    fi
}

add_tree_entries() {
    local dir=$1

    [[ -d "$dir" ]] || return 0
    while IFS= read -r -d '' path; do
        add_target "$path"
    done < <(find -P "$dir" -mindepth 1 ! -type d -print0 2>/dev/null)
}

add_shard_entries() {
    local dir=$1
    local shard_id=$2

    [[ -d "$dir" ]] || return 0
    while IFS= read -r -d '' path; do
        add_target "$path"
    done < <(find -P "$dir" -maxdepth 1 ! -type d -name "shard_${shard_id}.*" -print0 2>/dev/null)
}

iter_chr_dirs() {
    find -P "$BASE_DIR" -mindepth 2 -maxdepth 2 -type d | sort
}

iter_shard_ids() {
    local search_paths=()
    local dir

    for dir in "$@"; do
        if [[ -d "$dir" ]]; then
            search_paths+=("$dir")
        fi
    done

    (( ${#search_paths[@]} > 0 )) || return 0

    find -P "${search_paths[@]}" -maxdepth 1 ! -type d -name 'shard_*' -printf '%f\n' 2>/dev/null \
        | awk '{ sub(/^shard_/, "", $0); sub(/\..*$/, "", $0); print }' \
        | sort -u
}

is_valid_verified_marker() {
    local marker_file=$1
    local reference_file=$2
    local stored_mtime
    local current_mtime

    [[ -f "$marker_file" ]] || return 1
    [[ -f "$reference_file" ]] || return 1

    stored_mtime=$(marker_field "$marker_file" "RUST_BINARY_MTIME" || true)
    current_mtime=$(entry_mtime "$reference_file" || true)

    [[ -n "$stored_mtime" && -n "$current_mtime" && "$stored_mtime" == "$current_mtime" ]]
}

collect_stale_java_targets() {
    local jar_mtime
    local chr_dir
    local label
    local chr
    local shard_id
    local output_file
    local output_mtime

    require_reference_file "$JAVA_JAR" "Java JAR"
    jar_mtime=$(entry_mtime "$JAVA_JAR")

    while IFS= read -r chr_dir; do
        label=$(basename "$(dirname "$chr_dir")")
        chr=$(basename "$chr_dir")
        matches_scope "$label" "$chr" || continue

        while IFS= read -r shard_id; do
            output_file="$chr_dir/java/shard_${shard_id}.tsv"
            if [[ ! -e "$output_file" && ! -L "$output_file" ]]; then
                add_shard_entries "$chr_dir/java" "$shard_id"
                continue
            fi

            output_mtime=$(entry_mtime "$output_file" || true)
            if [[ -z "$output_mtime" || "$jar_mtime" -gt "$output_mtime" ]]; then
                add_shard_entries "$chr_dir/java" "$shard_id"
            fi
        done < <(iter_shard_ids "$chr_dir/java")
    done < <(iter_chr_dirs)
}

collect_stale_rust_targets() {
    local rust_mtime
    local chr_dir
    local label
    local chr
    local shard_id
    local output_file
    local output_mtime

    require_reference_file "$RUST_BIN" "Rust binary"
    rust_mtime=$(entry_mtime "$RUST_BIN")

    while IFS= read -r chr_dir; do
        label=$(basename "$(dirname "$chr_dir")")
        chr=$(basename "$chr_dir")
        matches_scope "$label" "$chr" || continue

        while IFS= read -r shard_id; do
            output_file="$chr_dir/rust/shard_${shard_id}.tsv"
            if [[ ! -e "$output_file" && ! -L "$output_file" ]]; then
                add_shard_entries "$chr_dir/rust" "$shard_id"
                add_shard_entries "$chr_dir/diff" "$shard_id"
                continue
            fi

            output_mtime=$(entry_mtime "$output_file" || true)
            if [[ -z "$output_mtime" || "$rust_mtime" -gt "$output_mtime" ]]; then
                add_shard_entries "$chr_dir/rust" "$shard_id"
                add_shard_entries "$chr_dir/diff" "$shard_id"
            fi
        done < <(iter_shard_ids "$chr_dir/rust" "$chr_dir/diff")
    done < <(iter_chr_dirs)
}

collect_failed_only_targets() {
    local chr_dir
    local label
    local chr
    local shard_id
    local status_file
    local verified_file
    local should_remove

    while IFS= read -r chr_dir; do
        label=$(basename "$(dirname "$chr_dir")")
        chr=$(basename "$chr_dir")
        matches_scope "$label" "$chr" || continue

        while IFS= read -r shard_id; do
            should_remove=0
            status_file="$chr_dir/diff/shard_${shard_id}.status"
            verified_file="$chr_dir/diff/shard_${shard_id}.verified"

            # Only remove shards with explicit FAIL status
            if [[ -f "$status_file" ]] && grep -qx 'FAIL' "$status_file"; then
                should_remove=1
            fi

            if (( should_remove )); then
                add_shard_entries "$chr_dir/java" "$shard_id"
                add_shard_entries "$chr_dir/rust" "$shard_id"
                add_shard_entries "$chr_dir/diff" "$shard_id"
            fi
        done < <(iter_shard_ids "$chr_dir/java" "$chr_dir/rust" "$chr_dir/diff")
    done < <(iter_chr_dirs)
}

collect_verified_only_targets() {
    local chr_dir
    local label
    local chr

    [[ -f "$RUST_BIN" ]] || return 0

    while IFS= read -r chr_dir; do
        label=$(basename "$(dirname "$chr_dir")")
        chr=$(basename "$chr_dir")
        matches_scope "$label" "$chr" || continue

        if is_valid_verified_marker "$chr_dir/.verified" "$RUST_BIN"; then
            add_tree_entries "$chr_dir"
        fi
    done < <(iter_chr_dirs)
}

collect_all_targets() {
    if (( ${#CONFIG_FILTERS[@]} == 0 && ${#CHR_FILTERS[@]} == 0 )); then
        add_tree_entries "$BASE_DIR"
        return 0
    fi

    local chr_dir
    local label
    local chr

    while IFS= read -r chr_dir; do
        label=$(basename "$(dirname "$chr_dir")")
        chr=$(basename "$chr_dir")
        matches_scope "$label" "$chr" || continue
        add_tree_entries "$chr_dir"
    done < <(iter_chr_dirs)
}

collect_targets() {
    case "$OPERATION" in
        stale-java)
            collect_stale_java_targets
            ;;
        stale-rust)
            collect_stale_rust_targets
            ;;
        failed-only)
            collect_failed_only_targets
            ;;
        verified-only)
            collect_verified_only_targets
            ;;
        all)
            collect_all_targets
            ;;
        *)
            fail "Unsupported operation: $OPERATION"
            ;;
    esac
}

build_summary() {
    local path
    local rel_path
    local group
    local size
    local total_bytes=0
    local total_files=0

    GROUP_BYTES=()
    GROUP_COUNTS=()

    for path in "${TARGETS[@]}"; do
        rel_path=${path#${BASE_DIR}/}
        group=$(dirname "$rel_path")
        if [[ "$group" == "." ]]; then
            group="./"
        else
            group="${group}/"
        fi

        size=$(entry_size "$path")
        GROUP_BYTES["$group"]=$(( ${GROUP_BYTES["$group"]:-0} + size ))
        GROUP_COUNTS["$group"]=$(( ${GROUP_COUNTS["$group"]:-0} + 1 ))
        total_bytes=$((total_bytes + size))
        total_files=$((total_files + 1))
    done

    if (( CONFIRM )); then
        printf 'DELETED: %s across %d files\n' "$(human_bytes "$total_bytes")" "$total_files"
    else
        printf 'DRY-RUN: Would remove %s across %d files\n' "$(human_bytes "$total_bytes")" "$total_files"
    fi

    while IFS= read -r group; do
        printf '  %s: %s (%d files)\n' \
            "$group" \
            "$(human_bytes "${GROUP_BYTES["$group"]}")" \
            "${GROUP_COUNTS["$group"]}"
    done < <(printf '%s\n' "${!GROUP_BYTES[@]}" | sort)
}

perform_deletions() {
    local path

    (( CONFIRM )) || return 0

    for path in "${TARGETS[@]}"; do
        rm -f -- "$path"
    done

    find -P "$BASE_DIR" -mindepth 1 -depth -type d -empty -delete 2>/dev/null || true
}

main() {
    parse_args "$@"
    resolve_base_dir
    collect_targets

    if (( ${#TARGETS[@]} == 0 )); then
        echo "Nothing to clean"
        exit 0
    fi

    build_summary
    perform_deletions
}

main "$@"