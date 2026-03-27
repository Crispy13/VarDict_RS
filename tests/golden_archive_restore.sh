#!/usr/bin/env bash
set -euo pipefail

# Restores a golden Java reference archive to the parity cache directory.
# Usage: tests/golden_archive_restore.sh ARCHIVE_FILE

PROJECT_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$PROJECT_ROOT"

CACHE_DIR="tmp/na12878_parity"

usage() {
    cat <<EOF
Usage: $(basename "$0") ARCHIVE_FILE

Restores Java parity cache from a golden reference archive.

Options:
  --help    Show this help

The archive must be a .tar.zst file created by golden_archive_create.sh.
Files are extracted to their original paths under tmp/na12878_parity/.
EOF
}

if [[ $# -lt 1 ]] || [[ "$1" == "--help" ]]; then
    usage
    exit 0
fi

ARCHIVE="$1"

if [[ ! -f "$ARCHIVE" ]]; then
    echo "ERROR: Archive not found: $ARCHIVE" >&2
    exit 1
fi

# Check for zstd
if ! command -v zstd &>/dev/null; then
    echo "ERROR: zstd not found. Install with: sudo apt install zstd" >&2
    exit 1
fi

# Verify checksum if available
CHECKSUM_FILE="${ARCHIVE}.sha256"
if [[ -f "$CHECKSUM_FILE" ]]; then
    echo "Verifying checksum..."
    if sha256sum --check "$CHECKSUM_FILE" --quiet 2>/dev/null; then
        echo "  Checksum OK"
    else
        echo "ERROR: Checksum verification failed!" >&2
        echo "  Expected: $(cat "$CHECKSUM_FILE")" >&2
        echo "  Archive may be corrupted." >&2
        exit 1
    fi
else
    echo "WARNING: No checksum file found (${CHECKSUM_FILE}), skipping verification"
fi

# Extract
echo "Restoring archive: $ARCHIVE"
zstd -d "$ARCHIVE" --stdout | tar xf -

# Count restored directories
RESTORED_DIRS=$(find "$CACHE_DIR" -type d -name "java" | wc -l)
RESTORED_FILES=$(find "$CACHE_DIR" -type d -name "java" -exec find {} -type f \; | wc -l)

echo ""
echo "Archive restored successfully:"
echo "  Directories: ${RESTORED_DIRS}"
echo "  Files: ${RESTORED_FILES}"
echo ""
echo "You can now run parity tests with --rust-only:"
echo "  tests/na12878_option_parity_v2.sh --rust-only"