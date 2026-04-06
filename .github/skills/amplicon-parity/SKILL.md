---
name: amplicon-parity
description: "Scaffold and eventually run amplicon-mode parity checks for VarDict-rs against VarDictJava. Use when working with amplicon BED inputs, targeted-mode fixtures, or amplicon parity smoke tests."
argument-hint: "Describe the amplicon parity task, e.g. 'add amplicon BED fixture' or 'investigate amplicon boundary mismatch'."
deprecated: true
---

# Amplicon Parity

Deprecated: This skill has been superseded by `pending-modes`. Use `pending-modes` for new amplicon-mode parity work.

Placeholder skill for future amplicon-mode parity work.

## Current Status

- Amplicon smoke tests exist only as ignored placeholders.
- No committed amplicon BAM plus BED fixture set is available yet.
- The parity harness does not yet have dedicated amplicon-mode coverage.

## Intended Scope

Use this skill when work involves:
- amplicon BED inputs via `-a`
- targeted-region Java vs Rust parity checks
- amplicon boundary and overlap behavior
- extending parity runners for amplicon-mode workflows

## Required Inputs

- sample BAM and BAI
- reference FASTA and FAI
- amplicon BED
- target regions or fixture manifest
- Java reference output

## Initial Workflow

1. Validate that BAM, BED, and reference naming are consistent.
2. Reproduce Java output with a fixed amplicon command line.
3. Run Rust with the same inputs and flags.
4. Compare output byte-for-byte.
5. Expand from small targeted smoke cases to broader coverage once fixtures exist.

See `tests/parity/somatic_amplicon_README.md` for the current fixture requirements.