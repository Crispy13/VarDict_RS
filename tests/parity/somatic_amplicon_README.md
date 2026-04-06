# Somatic and Amplicon Parity Test Data Requirements

This directory documents the placeholder test infrastructure for VarDict somatic and amplicon parity.

Current status:
- Simple-mode parity infrastructure exists and is actively used.
- Somatic-mode parity tests are scaffolded but blocked on tumor/normal fixture data.
- Amplicon-mode parity tests are scaffolded but blocked on amplicon BAM plus BED fixtures.

## Goal

The long-term goal is the same as Simple mode: byte-identical Rust output compared to VarDictJava for the same inputs, options, and reference data.

Mode-specific output expectations:
- Somatic mode: 55 output columns
- Amplicon mode: 38 output columns

## Somatic Mode Requirements

VarDict somatic mode requires paired inputs:
- Tumor BAM and index
- Normal BAM and index
- Shared reference FASTA and FAI
- Region BED or a fixed `-R` region
- Exact Java command line used to generate the expected output

Minimum dataset requirements:
- Tumor and normal BAMs aligned to the same reference build
- Identical chromosome naming between BAMs, BED, and FASTA
- BAMs small enough to run in local development, ideally a few targeted regions first
- At least one region with a clear somatic SNV
- At least one region with an indel
- At least one region that exercises somatic filtering behavior

Preferred additional coverage:
- Low-VAF somatic event
- Region with no somatic call to confirm empty-output parity
- Structural-variant-sensitive region if somatic mode is used with SV-relevant options

Expected command shape:

```bash
java -jar VarDict.jar \
  -G ref.fa \
  -b tumor.bam \
  -b2 normal.bam \
  -N tumor_sample|normal_sample \
  [other flags] \
  regions.bed
```

Rust should be run with the same reference, BAM pair, sample names, and flags. Output comparison must be byte-for-byte.

## Amplicon Mode Requirements

VarDict amplicon mode requires:
- BAM and index
- Reference FASTA and FAI
- Amplicon BED describing target intervals
- Exact Java command line used to generate the expected output

Minimum dataset requirements:
- BED file with valid amplicon intervals for the same reference build as the BAM
- At least one amplicon with a clear SNV
- At least one amplicon with an indel
- At least one boundary case near amplicon edges or overlapping amplicons

Preferred additional coverage:
- Overlapping amplicons
- Off-target or low-coverage amplicon behavior
- Regions that exercise amplicon-specific output columns

Expected command shape:

```bash
java -jar VarDict.jar \
  -G ref.fa \
  -b sample.bam \
  -a amplicons.bed \
  -N sample \
  [other flags] \
  regions.bed
```

Rust should be run with the same BAM, amplicon BED, target regions, and flags. Output comparison must be byte-for-byte.

## Suggested Dataset Strategy

Prefer staged rollout rather than large whole-BAM parity immediately.

Recommended order:
1. One tiny somatic tumor/normal region with a known SNV
2. One tiny somatic region with a known indel
3. One tiny amplicon region with a known SNV
4. One tiny amplicon region with a boundary-sensitive case
5. Broader smoke coverage after the first four cases are stable

If public datasets are used, keep the first fixtures narrow and reproducible. The main requirement is deterministic Java output that can be checked into fixture form or regenerated reliably.

## Harness Extension Notes

The current parity harness is Simple-mode-oriented. Somatic and amplicon coverage will need dedicated runner support.

Planned extensions:
- Somatic runner support for paired `-b` and `-b2` inputs
- Amplicon runner support for `-a amplicon.bed`
- Mode-specific fixture manifests
- Small smoke tests first, then shard or fixture expansion as data becomes available

The placeholder Rust tests in:
- `tests/somatic_parity_smoke.rs`
- `tests/amplicon_parity_smoke.rs`

exist only to reserve names, compile cleanly, and document the intended future coverage.

## Blocking Inputs Checklist

Somatic:
- Tumor BAM
- Tumor BAI
- Normal BAM
- Normal BAI
- Reference FASTA
- Reference FAI
- Region BED or exact region list
- Java reference output

Amplicon:
- Sample BAM
- Sample BAI
- Reference FASTA
- Reference FAI
- Amplicon BED
- Region BED or exact region list
- Java reference output

Until those inputs exist in the repository or an agreed external fixture source, the smoke tests remain intentionally ignored.