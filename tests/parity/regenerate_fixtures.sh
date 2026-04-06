#!/usr/bin/env bash
# Regenerate micro smoke test fixtures from VarDictJava.
# Run from project root: bash tests/parity/regenerate_fixtures.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

JAR="$PROJECT_ROOT/VarDictJava/build/libs/VarDict-1.8.3.jar"
REF="$PROJECT_ROOT/testdata/hs37d5.fa"
BAM="$PROJECT_ROOT/testdata/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam"
FIXTURE_DIR="$PROJECT_ROOT/tests/parity/fixtures"
TMP_DIR="$PROJECT_ROOT/tmp"

# Verify prerequisites
for f in "$JAR" "$REF" "$BAM" "${REF}.fai"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: required file not found: $f" >&2
        exit 1
    fi
done

mkdir -p "$TMP_DIR" "$FIXTURE_DIR"

# Test cases: config chrom start end
CASES=(
    "default 20 126269 126333"
    "default 20 168500 168800"
    "default 20 168600 168800"
    "default 20 25456879 25457078"
    "default 20 30000000 30000300"
    "default 22 40000000 40000300"
    "default 22 42522500 42522800"
    "default MT 1 300"
    "default MT 300 600"
    "nosv 20 126269 126333"
    "nosv 20 168500 168800"
    "nosv 20 168600 168800"
    "nosv 20 25456879 25457078"
    "nosv 20 30000000 30000300"
    "nosv 22 40000000 40000300"
    "nosv 22 42522500 42522800"
    "nosv MT 1 300"
    "nosv MT 300 600"
)

TOTAL=${#CASES[@]}
PASS=0
WARN=0

echo "=== Regenerating $TOTAL micro smoke fixtures from VarDictJava ==="
echo ""

for case_str in "${CASES[@]}"; do
    read -r CONFIG CHROM START END <<< "$case_str"
    FIXTURE_NAME="${CONFIG}_chr${CHROM}_${START}_${END}.expected.tsv"
    FIXTURE_PATH="$FIXTURE_DIR/$FIXTURE_NAME"
    REGION="${CHROM}:${START}-${END}"

    # Build Java command
    JAVA_CMD=(java -classpath "$JAR" com.astrazeneca.vardict.Main
        -G "$REF"
        -b "$BAM"
        -N smoke_sample
        -th 1
    )
    if [[ "$CONFIG" == "nosv" ]]; then
        JAVA_CMD+=(-U)
    fi
    JAVA_CMD+=(-R "$REGION")

    # Run and capture
    echo -n "  [$((PASS + WARN + 1))/$TOTAL] $FIXTURE_NAME ... "
    if "${JAVA_CMD[@]}" > "$FIXTURE_PATH" 2>"$TMP_DIR/regen_stderr_${CONFIG}_${CHROM}.txt"; then
        LINES=$(wc -l < "$FIXTURE_PATH")
        if [[ "$LINES" -eq 0 ]]; then
            echo "WARN: empty output (0 lines)"
            WARN=$((WARN + 1))
        else
            echo "OK ($LINES lines)"
            PASS=$((PASS + 1))
        fi
    else
        echo "FAIL: Java exited with non-zero status"
        echo "  stderr: $(cat "$TMP_DIR/regen_stderr_${CONFIG}_${CHROM}.txt")"
        WARN=$((WARN + 1))
    fi
done

echo ""
echo "=== Done: $PASS OK, $WARN warnings out of $TOTAL ==="

# Clean up stderr temp files
rm -f "$TMP_DIR"/regen_stderr_*.txt