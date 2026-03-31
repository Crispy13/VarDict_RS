---
name: mismatch-triage
description: "Classify parity mismatches by severity and responsible module. Use when: classify mismatches, prioritize parity failures, triage output differences, batch mismatch analysis, severity classification."
argument-hint: "Specify the config label, chromosome, or 'batch' for all current failures."
---

# Mismatch Triage Workflow

## Purpose

Classify parity mismatches by severity and map them to responsible modules for prioritized fixing.

## Scope

This skill is for classification and prioritization only. Do not fix code while using it.

## Operating Rules

- Always process all mismatches in scope, not just the first one.
- Group related mismatches when they likely come from the same root cause.
- Call out cascading mismatches when one upstream bug likely caused multiple downstream diffs.
- Treat line count and column count differences as higher priority than value-level differences.

## Procedure

### Step 1: Gather mismatches

- Input a list of `(line, column, java_value, rust_value)` tuples from the `shard-diagnosis` skill or from diff files.
- Or aggregate mismatches from multiple shard diffs in a single config x chromosome cell.
- Record shard identifier, position, and field name when available so the final report stays actionable.

### Step 2: Apply Decision Tree

For each mismatch, walk this tree in order:

```text
1. Is it a LINE COUNT difference?
   |-- Rust has MORE lines -> CRITICAL: Extra variant (incorrect filter logic)
   |-- Rust has FEWER lines -> CRITICAL: Missing variant (unimplemented branch)
   `-- Same line count -> go to 2

2. Is it a COLUMN COUNT difference?
   |-- Different column counts -> CRITICAL: Missing/extra output field
   `-- Same column count -> go to 3

3. Is the value NUMERIC?
   |-- No -> go to 5 (string comparison)
   `-- Yes -> go to 4

4. NUMERIC comparison:
   |-- Difference < 1e-10 -> LOW: Floating-point noise (likely harmless)
   |-- Values differ only in trailing zeros (0.125 vs 0.1250) -> HIGH: Float formatting
   |-- Difference is exactly in last decimal place -> HIGH: Rounding method
   |-- Values are integers that differ -> go to 4a
   `-- Large numeric difference -> HIGH: Computation error

   4a. INTEGER differences:
   |-- Depth/coverage columns (8-13) -> HIGH: Read counting error
   |-- SV columns (34-35) -> HIGH: SV detection logic
   `-- Other -> HIGH: Accumulation error

5. STRING comparison:
   |-- One is empty, other is not -> CRITICAL: Missing data
   |-- Genotype format differs (e.g., "A/T" vs "T/A") -> HIGH: Allele ordering
   |-- VarType differs -> CRITICAL: Variant classification error
   |-- Variant description differs (+seq, -N, #seq) -> HIGH: CigarParser
   `-- Minor formatting -> MEDIUM: Output formatting
```

### Step 3: Classify Severity

| Severity | Definition | Action |
|----------|------------|--------|
| CRITICAL | Missing/extra variants, wrong variant type, empty fields | Fix immediately - blocks parity |
| HIGH | Wrong numeric values, float formatting, allele ordering | Fix in current stage |
| MEDIUM | Output formatting, whitespace, cosmetic | Fix before milestone close |
| LOW | Floating-point noise below measurable threshold | Document and defer |

### Step 4: Map to Responsible Module

Reference the column-to-module mapping from Step 5 of the `shard-diagnosis` skill instead of duplicating the table here.

Summarize the classified mismatches by responsible module using this format:

```text
## Triage Summary

### By Module
| Module | CRITICAL | HIGH | MEDIUM | LOW | Total |
|--------|----------|------|--------|-----|-------|
| CigarParser | 0 | 2 | 0 | 1 | 3 |
| ToVarsBuilder | 1 | 0 | 0 | 0 | 1 |
| ... | | | | | |

### Fix Priority (ordered)
1. [CRITICAL] ToVarsBuilder - Missing Complex variant at chr3:60830764
2. [HIGH] CigarParser - Read count mismatch at chr1:5934302
...
```

When several mismatches land in the same module and share the same category, collapse them into one root-cause group before ranking fixes.

### Step 5: Produce Triage Report

Output the result in this format:

```text
## Mismatch Triage Report: <scope>

**Input**: {N} mismatches across {M} shards
**Date**: YYYY-MM-DD

### Classification Summary
| Severity | Count | % |
|----------|-------|---|
| CRITICAL | | |
| HIGH | | |
| MEDIUM | | |
| LOW | | |

### Mismatches by Severity

#### CRITICAL
| Shard | Position | Col | Field | Java | Rust | Module |
|-------|----------|-----|-------|------|------|--------|

#### HIGH
...

### Module Impact
| Module | Issues | Recommended Action |
|--------|--------|--------------------|

### Recommended Fix Order
1. {First fix - highest impact}
2. {Second fix}
...
```

## Triage Heuristics

- If one mismatch changes `Ref`, `Alt`, or `VarType`, treat downstream count and formatting mismatches on that row as likely cascading.
- If many rows fail in the same field and module, prioritize the earliest root-cause example and mark the rest as duplicates of the same bug.
- If a batch contains both CRITICAL and HIGH issues in the same module, fix the CRITICAL path first before investigating the value-level mismatches.
- If the same mismatch pattern repeats across shards in one config, report it once as a grouped issue and include the shard count.

## Constraints

- Do not fix code - only classify and prioritize.
- Always process all mismatches, not just the first one.
- Group related mismatches together when they likely share the same root cause.
- Note when multiple mismatches are likely caused by the same bug.