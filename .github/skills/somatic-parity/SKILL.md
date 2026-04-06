---
name: somatic-parity
description: "Scaffold and eventually run somatic-mode parity checks for VarDict-rs against VarDictJava. Use when working with tumor/normal paired mode, somatic fixture planning, or somatic parity smoke tests."
argument-hint: "Describe the somatic parity task, e.g. 'add tumor/normal smoke fixture' or 'diagnose somatic output mismatch'."
deprecated: true
---

# Somatic Parity

Deprecated: This skill has been superseded by `pending-modes`. Use `pending-modes` for new somatic-mode parity work.

Placeholder skill for future somatic-mode parity work.

## Current Status

- Somatic smoke tests exist only as ignored placeholders.
- No committed tumor/normal fixture set is available yet.
- The parity harness is currently centered on Simple mode.

## Intended Scope

Use this skill when work involves:
- tumor/normal paired BAM inputs
- Java vs Rust somatic output comparison
- somatic fixture selection and minimization
- extending parity runners for `-b` and `-b2`

## Required Inputs

- tumor BAM and BAI
- normal BAM and BAI
- reference FASTA and FAI
- target regions or BED
- Java reference output

## Initial Workflow

1. Confirm tumor and normal inputs are aligned to the same reference build.
2. Reproduce Java output with a fixed command line.
3. Run Rust with the exact same flags and regions.
4. Compare output byte-for-byte.
5. Add or refine ignored smoke tests until stable fixture coverage exists.

See `tests/parity/somatic_amplicon_README.md` for the current fixture requirements.