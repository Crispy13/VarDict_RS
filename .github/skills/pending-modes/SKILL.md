---
name: pending-modes
description: "Track placeholder parity workflows for modes without committed fixtures. Use when working with amplicon-mode parity, somatic-mode parity, fixture planning for tumor/normal or targeted inputs, or staging future parity smoke tests."
argument-hint: "Describe the pending mode task, e.g. 'plan amplicon BED fixture' or 'diagnose somatic tumor/normal mismatch'."
---

# Pending Modes Parity

Placeholder skill for parity work that is not fully enabled yet because fixture coverage is incomplete.

## Current Status

- Amplicon smoke tests exist only as ignored placeholders.
- Somatic smoke tests exist only as ignored placeholders.
- No committed amplicon BAM plus BED fixture set is available yet.
- No committed tumor/normal fixture set is available yet.
- The parity harness is currently centered on Simple mode.

## Amplicon Mode Track

Use this track when work involves:
- amplicon BED inputs via `-a`
- targeted-region Java vs Rust parity checks
- amplicon boundary and overlap behavior
- extending parity runners for amplicon-mode workflows

### Required Inputs

- sample BAM and BAI
- reference FASTA and FAI
- amplicon BED
- target regions or fixture manifest
- Java reference output

### Initial Workflow

1. Validate that BAM, BED, and reference naming are consistent.
2. Reproduce Java output with a fixed amplicon command line.
3. Run Rust with the same inputs and flags.
4. Compare output byte-for-byte.
5. Expand from small targeted smoke cases to broader coverage once fixtures exist.

## Somatic Mode Track

Use this track when work involves:
- tumor/normal paired BAM inputs
- Java vs Rust somatic output comparison
- somatic fixture selection and minimization
- extending parity runners for `-b` and `-b2`

### Required Inputs

- tumor BAM and BAI
- normal BAM and BAI
- reference FASTA and FAI
- target regions or BED
- Java reference output

### Initial Workflow

1. Confirm tumor and normal inputs are aligned to the same reference build.
2. Reproduce Java output with a fixed command line.
3. Run Rust with the exact same flags and regions.
4. Compare output byte-for-byte.
5. Add or refine ignored smoke tests until stable fixture coverage exists.

See `tests/parity/somatic_amplicon_README.md` for the current fixture requirements.
