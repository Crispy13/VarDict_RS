#!/usr/bin/env bash
set -euo pipefail

# Creates a golden Java reference archive from existing parity cache.
# Usage: tests/golden_archive_create.sh [--include-pileup] [--output FILE]

PROJECT_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$PROJECT_ROOT"

CACHE_DIR="tmp/na12878_parity"
INCLUDE_PILEUP=0
OUTPUT=""

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]

Creates a zstd-compressed tarball of Java parity cache for golden reference archiving.

Options:
  --include-pileup  Include pileup configs (pileup, pileup-max, pileup-strict, pileup-relaxed)
                    WARNING: Pileup configs can be 50-80 GB
  --output FILE     Output archive path (default: tmp/golden_java_YYYYMMDD.tar.zst)
  --help            Show this help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --include-pileup)
            INCLUDE_PILEUP=1
            ;;
        --output)
            shift
            [[ $# -gt 0 ]] || { echo "ERROR: Missing value for --output" >&2; exit 1; }
            OUTPUT="$1"
            ;;
        --help)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: Unknown option: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
    shift
done

if [[ -z "$OUTPUT" ]]; then
    OUTPUT="tmp/golden_java_$(date +%Y%m%d).tar.zst"
fi

if [[ ! -d "$CACHE_DIR" ]]; then
    echo "ERROR: Cache directory not found: $CACHE_DIR" >&2
    echo "Run parity tests first to generate Java cache." >&2
    exit 1
fi

# Check for zstd
if ! command -v zstd &>/dev/null; then
    echo "ERROR: zstd not found. Install with: sudo apt install zstd" >&2
    exit 1
fi

# Build exclude patterns for pileup configs
if (( ! INCLUDE_PILEUP )); then
    echo "Excluding pileup configs (use --include-pileup to include)"
fi

# Only archive java/ subdirectories (not rust/, diff/, or markers)
# Find all java/ dirs and create tar from them
echo "Scanning Java cache directories..."
JAVA_DIRS=()
while IFS= read -r -d '' dir; do
    JAVA_DIRS+=("$dir")
done < <(find "$CACHE_DIR" -type d -name "java" -print0 | sort -z)

if (( ${#JAVA_DIRS[@]} == 0 )); then
    echo "ERROR: No Java cache directories found in $CACHE_DIR" >&2
    exit 1
fi

# Filter out pileup if not included
if (( ! INCLUDE_PILEUP )); then
    FILTERED_DIRS=()
    for dir in "${JAVA_DIRS[@]}"; do
        if [[ "$dir" != *"/pileup/"* && "$dir" != *"/pileup-max/"* && "$dir" != *"/pileup-strict/"* && "$dir" != *"/pileup-relaxed/"* ]]; then
            FILTERED_DIRS+=("$dir")
        fi
    done
    JAVA_DIRS=("${FILTERED_DIRS[@]}")
fi

if (( ${#JAVA_DIRS[@]} == 0 )); then
    echo "ERROR: No non-pileup Java cache directories found. Use --include-pileup?" >&2
    exit 1
fi

echo "Found ${#JAVA_DIRS[@]} Java cache directories"

# Count total files
TOTAL_FILES=$(find "${JAVA_DIRS[@]}" -type f | wc -l)
echo "Total files: ${TOTAL_FILES}"

# Calculate uncompressed size
UNCOMPRESSED_SIZE=$(du -shc "${JAVA_DIRS[@]}" | tail -1 | cut -f1)
echo "Uncompressed size: ${UNCOMPRESSED_SIZE}"

# Create archive
mkdir -p "$(dirname "$OUTPUT")"
echo "Creating archive: $OUTPUT"
tar cf - "${JAVA_DIRS[@]}" | zstd -T0 -3 -o "$OUTPUT"

ARCHIVE_SIZE=$(du -h "$OUTPUT" | cut -f1)
echo ""
echo "Archive created successfully:"
echo "  Path: $OUTPUT"
echo "  Size: ${ARCHIVE_SIZE} (uncompressed: ${UNCOMPRESSED_SIZE})"
echo "  Configs: ${#JAVA_DIRS[@]} directories"
echo "  Files: ${TOTAL_FILES}"

# Generate checksum
sha256sum "$OUTPUT" > "${OUTPUT}.sha256"
echo "  Checksum: ${OUTPUT}.sha256"