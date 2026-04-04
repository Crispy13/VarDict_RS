# to_vars_builder

**Source**: `src/mods/to_vars_builder.rs`
**LOC**: ~1,414
**Java counterpart**: `ToVarsBuilder.java` → [Java cache](../java/ToVarsBuilder.md)
**Status**: complete
**Last verified**: 2026-04-04

## Overview

The `to_vars_builder` module is the statistics calculation and variant refinement stage. It transforms raw `Variation` objects (from CigarParser and VariationRealigner) into fully-formed `Variant` objects with 20+ statistical properties: allele frequencies, strand bias, mean quality/position/mapping-quality, MSI scores, shift3 values, and flanking sequences. **Important architectural note**: much of the heavy logic from Java's `ToVarsBuilder.process()` — including `collectReferenceVariants()`, MSI/shift3, and CRISPR — lives in `vardict_pipeline.rs::run_to_vars_builder()` in Rust, not in this module.

## Public API

| Function/Method | Purpose |
|----------------|---------|
| `ToVarsBuilder::new()` | Constructor (min_freq=0.02, qual=20, mapq=20) |
| `with_min_frequency()` | Builder: frequency threshold |
| `with_quality_threshold()` | Builder: base quality threshold |
| `with_mapq_threshold()` | Builder: mapping quality threshold |
| `calculate_variant_statistics()` | Core stats engine — accepts reads, outputs populated Variant |
| `build_variants()` | Groups by position+variant, returns position→Vars map |
| `check_strand_bias()` | ≤12 reads: both-strand check; >12: ratio test |
| `determine_genotype()` | "REF/ALT" formatting |
| `validate_ref_allele()` | IUPAC ambiguity → first base |
| `calculate_mean_and_std()` | Statistics helper (mean + stddev) |
| `has_at_least_2_distinct()` | ≥2 unique values predicate |
| `calculate_shift3()` | 3' homopolymer shift calculation |
| `var_type_string()` | SNV/Insertion/Deletion/Complex classifier |
| `create_description_string()` | "A"/"+SEQ"/"-N" formatting |

### Types

| Type | Purpose |
|------|---------|
| `VarType` | enum — SNV, Insertion, Deletion, Complex (Rust-only, no Java counterpart) |
| `StrandBiasValue` | enum — CantAssess(0), HasBias(1), NoBias(2) |
| `StrandBiasFlag` | struct — ref_bias + var_bias; formats "X;Y" |
| `Variant` | struct ~30 fields — fully-formed variant for output |
| `Vars` | variants + reference_variant + sv + sv_flags |

## Java Correspondence

| Rust | Java | Notes |
|------|------|-------|
| `ToVarsBuilder` struct | `ToVarsBuilder` class | Config wrapper; builder pattern vs direct field setting |
| `calculate_variant_statistics()` | `createVariant()` + `createInsertion()` | Partial port; basic stats only |
| `build_variants()` | `process()` top-level loop | In Rust, the main loop lives in `vardict_pipeline.rs::run_to_vars_builder()` |
| `check_strand_bias()` | `VariationUtils.strandBias()` | Direct port |
| `determine_genotype()` | `Variant.genotype()` | Partial port |
| `validate_ref_allele()` | `ToVarsBuilder.validateRefallele()` | Direct port |
| `var_type_string()` | `Variant.varType()` | Direct port |
| `calculate_shift3()` | `findMSI()` (shift3 subset)  | Direct port |
| *(in vardict_pipeline)* | `collectReferenceVariants()` ~340 LOC | Moved to pipeline orchestrator |
| *(in vardict_pipeline)* | `process()` main iteration loop | Moved to pipeline orchestrator |

### Internal Helpers

| Function | Purpose |
|----------|---------|
| `infer_variant_type()` | Parses description string into `VarType` enum |
| `var_key_to_var_type()` | Converts `SimpleVarKey` to `VarType` |
| `Variant::adj_complex()` | Trims common prefix/suffix from ref/var alleles |
| `Variant::is_good_var()` | Simple quality filter (strand bias check) |
| `StrandBiasFlag::is_ref_good_var_biased()` | Checks "2;1" pattern |

## Known Parity Traps

1. **Division-by-Zero / NaN Handling** — Java `mean/count` where count=0 → NaN/Infinity, formatted silently. Variants with varsCount=0 skipped in `createVariant()` but NOT in `createInsertion()`.

2. **High-to-Low Quality Ratio** — Java uses `hicnt / (locnt != 0 ? locnt : 0.5d)` — divisor is 0.5 not 1.0 when locnt=0.

3. **Insertion vs Non-Insertion Different Thresholds** — `extracnt > 0` (strict positive) for non-insertions vs `extracnt != 0` (non-zero) for insertions. Intentional asymmetry.

4. **hicov Local Mutation Across Loop** — `hicov` modified during loop and carries to next iteration. Must not reset by scoping.

5. **Reference Coverage from Next Position** — Insertion with "&" in description: override `totalPosCoverage` to position+1 coverage.

6. **Variant Mutations During Iteration** — `createInsertion()` mutates reference variant at position+1 (decrements strand counts). Required for parity.

7. **LinkedHashMap Insertion Order** — Java uses `LinkedHashMap` for `VariationMap`. Rust must use `IndexMap` or equivalent ordered map.

8. **SV Marker Sentinel** — Key `"SV"` is a sentinel for structural variant metadata. Must handle specially during iteration.

9. **Float Formatting (HALF_EVEN)** — Java `DecimalFormat("0.0000")` with HALF_EVEN rounding. Rust `format!("{:.4}")` may differ; explicit `round_half_even()` needed.

10. **shift3/MSI Floor Interaction** — If `msi <= shift3 / dellen`, then `msi = shift3 / dellen`. Must match exact double division.

11. **chrLengths Side Effect** — Java `getOrElse(chrLengths, chr, 0)` **mutates** the map by inserting 0 if missing.

12. **non_insertion_vars `.remove()` Safety** (from repo memory `non_insertion_vars_current_position_remove_safe_in_tovars_20260403`) — Java does NOT remove entries during `nonInsertionVariants` iteration. Using `.remove(&position)` in Rust is unsafe because position+1 lookups in `collectReferenceVariants()` and insertion coverage adjustment depend on entries remaining alive. If hash ordering places P+1 before P, the P+1 entry is removed before P processes it, causing alt-depth adjustment to be skipped → parity mismatch.

13. **Insertion consume-and-drop** (from repo memory `tovars_insertion_consume_drop_20260324`) — `insertion_vars` CAN be consumed with `.remove(&position)` because the loop has no position+1 lookahead into `insertion_vars`. Safe pattern: `let ins_at_pos = insertion_vars.remove(&position);`.

## Divergences from Java

### Module Split (Major)
Java's `process()` main loop and `collectReferenceVariants()` (~340 LOC) are in `vardict_pipeline.rs::run_to_vars_builder()`, not in this module. This module contains only the stateless helper functions and type definitions. ~60% of Java implementation lives in this module; remaining ~40% in vardict_pipeline.

### Type System
- Rust adds `VarType` enum (not present in Java)
- Rust separates `threshold_frequency` from `frequency` (Java uses a single field)

### Builder Pattern
Rust uses builder pattern (`with_min_frequency()`, etc.) vs Java's direct field assignment.

### Null Safety
`Option<T>` throughout vs Java null checks. Every null-check branch in Java must have an `Option` match in Rust.

### Float Rounding
Explicit `round_half_even()` needed for exact parity — implemented in the pipeline, not in this module.

## Cross-Module Dependencies

### Called By
- `vardict_pipeline.rs::run_to_vars_builder()` — main consumer of all public functions
- Unit tests

### Calls Into
- `simple_variant_caller` — `SimpleVarKey` type
- `scopedata/global_read_only_scope` — configuration access
- `variants/variants` — `StructuralVariantCounts`

### Data Flow
```
RealignedOutput
  → vardict_pipeline::run_to_vars_builder()
    → check_strand_bias(), determine_genotype(), calculate_shift3(), etc.
    → Variant + Vars
  → AlignedVarsData
  → OutputVariant
```

## Completeness Notes

~60% of Java `ToVarsBuilder` logic resides in this Rust module; the remaining ~40% (`collectReferenceVariants`, full MSI computation, CRISPR, main iteration loop) is in `vardict_pipeline.rs`. When diagnosing parity issues related to ToVarsBuilder, both files must be consulted.
