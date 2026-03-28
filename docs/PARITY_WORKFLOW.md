# Parity Workflow

This repository uses VarDictJava as the reference implementation and checks Rust output against it byte-for-byte on the same reference genome, BAM, options, and genomic regions.

## Current workflow

1. Run the same configuration against VarDictJava and VarDict-rs.
2. Compare the output shard-by-shard.
3. Investigate the first mismatch and trace it back to the Java and Rust code paths.

## BAM parity status

A BAM is considered complete only when both the default-options baseline and the full option-matrix sweep pass with zero failures. If either is still failing, the BAM is treated as in progress.

| Dataset | Status | Notes |
| --- | --- | --- |
| `testdata/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam` | In progress | Default-options baseline passed. Option-matrix sweep (44 configs × all real contigs) is still in progress. |