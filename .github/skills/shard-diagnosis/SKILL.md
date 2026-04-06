---
name: shard-diagnosis
description: "Diagnose a specific failing shard in the VarDictJava-to-Rust parity workflow. Use when: shard failure, parity mismatch, column diff, output divergence, diagnose failing shard."
argument-hint: "Specify the shard identifier, e.g. 'chr20/shard_042' or '<config>/<chr>/shard_042'."
---

# Shard Diagnosis Workflow

## Purpose

Diagnose why a specific shard produces different output between Java and Rust.

## Scope

This skill is diagnostic only. Do not fix code while using it.

## Operating Rules

- Check for line count differences first, because extra or missing variants can make per-column analysis misleading.
- Always focus on the first differing output line.
- Always report the first divergent column, because that is often the root cause.
- If multiple lines differ, diagnose the first differing line before expanding the scope.

## 6-Step Procedure

### 1. Locate the failing shard

- Cache directory: `tmp/na12878_parity/<opts-label>/<chr>/`
- Check `diff/shard_NNN.diff` for the unified diff.
- Check `diff/shard_NNN.meta` for structured failure metadata.
- If no diff exists, check `diff/shard_NNN.status` for `PASS`, `FAIL`, or `EMPTY`.
- Verify whether Java and Rust TSV files have different line counts before comparing columns.

### 2. Extract divergent rows

Use the shard TSV outputs:

```bash
java_file="tmp/na12878_parity/<label>/<chr>/java/shard_NNN.tsv"
rust_file="tmp/na12878_parity/<label>/<chr>/rust/shard_NNN.tsv"
```

- Take the first differing line from the diff, or identify the first mismatching line by comparing the Java and Rust TSV files directly.
- Grep the position from both Java and Rust TSV files.
- Extract both line versions for the first differing position.

### 3. Column-by-column comparison

Use a tab-split comparison to identify the exact differing columns:

```bash
diff <(echo "$java_line" | tr '\t' '\n' | nl) \
     <(echo "$rust_line" | tr '\t' '\n' | nl)
```

Or use the awk column differ:

```bash
paste <(echo "$java_line") <(echo "$rust_line") | awk -F'\t' '{
    n = NF/2;
    for(i=1; i<=n; i++) {
        if($i != $(i+n)) printf "Col %d: Java=[%s] Rust=[%s]\n", i, $i, $(i+n)
    }
}'
```

- Capture every differing column.
- Mark the first divergent column explicitly.

### 4. Map columns to field names

Use this Simple mode column reference table for 36-column output:

| Col | Field | Col | Field |
|-----|-------|-----|-------|
| 1 | Sample | 19 | MeanQual |
| 2 | Gene | 20 | StdQual |
| 3 | Chr | 21 | MeanMapQual |
| 4 | Start | 22 | MeanMapQualAlt |
| 5 | End | 23 | MapQualMismatch |
| 6 | Ref | 24 | MapQualMismatchRate |
| 7 | Alt | 25 | NM |
| 8 | Depth | 26 | MSI |
| 9 | AltDepth | 27 | MSILen |
| 10 | RefFwdReads | 28 | Shift3 |
| 11 | RefRevReads | 29 | 5pFlankSeq |
| 12 | AltFwdReads | 30 | 3pFlankSeq |
| 13 | AltRevReads | 31 | Segment |
| 14 | Genotype | 32 | VarType |
| 15 | AF | 33 | Duprate |
| 16 | StrandBias | 34 | SplitReads |
| 17 | MeanPosition | 35 | SpanPairs |
| 18 | StdPosition | 36 | Filter |

### 5. Map field to responsible Rust module

Use this field-to-module map to narrow the likely Rust owner:

| Fields (Cols) | Responsible Module | Rust File |
|---------------|-------------------|-----------|
| Sample, Gene, Chr, Segment (1-3, 31) | Region/Config | `conf.rs`, `bin/vardict.rs` |
| Start, End, Ref, Alt (4-7) | CigarParser + Realigner | `cigar_parser.rs`, `variant_realigner.rs` |
| Depth, AltDepth, Reads (8-13) | ToVarsBuilder | `to_vars_builder.rs`, `vardict_pipeline.rs` |
| Genotype (14) | ToVarsBuilder | `to_vars_builder.rs` |
| AF, StrandBias (15-16) | ToVarsBuilder | `to_vars_builder.rs` |
| MeanPosition, StdPosition (17-18) | CigarParser (accumulation) | `cigar_parser.rs` |
| MeanQual, StdQual (19-20) | CigarParser (accumulation) | `cigar_parser.rs` |
| MeanMapQual, MeanMapQualAlt (21-22) | CigarParser | `cigar_parser.rs` |
| MapQualMismatch, Rate (23-24) | CigarParser | `cigar_parser.rs` |
| NM (25) | CigarParser | `cigar_parser.rs` |
| MSI, MSILen (26-27) | ToVarsBuilder | `to_vars_builder.rs` |
| Shift3 (28) | VariantRealigner | `variant_realigner.rs` |
| FlankSeqs (29-30) | VariantRealigner | `variant_realigner.rs` |
| VarType (32) | ToVarsBuilder | `to_vars_builder.rs` |
| Duprate (33) | CigarParser | `cigar_parser.rs` |
| SplitReads, SpanPairs (34-35) | StructuralVariantsProcessor | `structural_variants_processor.rs` |
| Filter (36) | OutputVariant | `output_variant.rs` |

### 6. Produce structured diagnosis report

Output the result in this format:

```text
## Shard Diagnosis: <opts-label>/<chr>/shard_NNN

**Position**: chr:pos
**Config**: <opts>

### Divergent Columns
| Col | Field | Java | Rust | Category |
|-----|-------|------|------|----------|

### Root Cause Module
Module: <module_name>
File: src/mods/<file>.rs

### Recommended Next Step

Follow this sequence — **test-first, Java-fixture-only**:

1. Delegate to java-analyst to trace column {N} logic (if root cause is unclear)
2. **Extract Java fixture**: Save the Java shard output as the expected reference (from `tmp/na12878_parity/<label>/<chr>/java/shard_NNN.tsv.gz` or by running Java). **Never use Rust output as the fixture** — that locks in Rust bugs.
3. **Write a failing test**: Add one new `#[test]` function in `tests/integration_test.rs` named `test_target_bam_{bam_slug}_{chr}_{description}_parity` (or `test_{module}_{description}_parity` for unit-level bugs). The test compares Rust output against the Java fixture. It must **fail before the fix**.
4. **Implement the fix**: Delegate to rust-implementer to fix the specific function.
5. **Verify the test passes**: `cargo test --profile debug-release -- --include-ignored`
6. **Run parity-fix-review**: Load `parity-fix-review` skill and run the review gate before reporting completion.
```

## Diagnostic Checklist

1. Confirm shard status: `PASS`, `FAIL`, or `EMPTY`.
2. Confirm whether the mismatch is a line count issue or a same-line column divergence.
3. Extract the first differing line only.
4. Identify the first divergent column and map it to a field name.
5. Map that field to the likely owning module.
6. Produce the structured report without proposing code changes beyond the next handoff target.