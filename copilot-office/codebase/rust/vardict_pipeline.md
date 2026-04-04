# VarDictPipeline

**Source**: `src/mods/vardict_pipeline.rs`
**LOC**: 7,447
**Java counterpart**: `modes/*`, `SAMFileParser`
**Status**: partial

## Overview
Owns the Rust simple/somatic/amplicon orchestration after read collection. It wires `CigarParser`, `VariantRealigner`, `StructuralVariantsProcessor`, and the to-vars / postprocessor stages, then applies Java-compatible output shaping such as genotype selection, reference-call materialization, and simple-mode filtering.

## Public API
| Function/Method | Purpose |
|----------------|---------|
| `process_region_from_bam()` | Main simple-mode region entry point using a shared `BamReader` |
| `process_region_from_bam_streaming()` | Streaming variant emission without collecting full output in memory |
| `process_region_to_aligned_vars_from_bam()` | Produces `AlignedVarsData` before postprocessing |
| `run_to_vars_builder()` | Converts realigned raw variation maps into `Vars` objects per position |
| `run_simple_post_processor()` | Emits final simple-mode tab-delimited rows |

## Java Correspondence
| Rust | Java | Notes |
|------|------|-------|
| `run_to_vars_builder()` | `ToVarsBuilder.process()`-style position materialization | Rust must retain low-frequency simple-mode positions until postprocessing; pruning earlier breaks Java low-AF SNV output |
| `collect_reference_variants()` | `ToVarsBuilder.collectReferenceVariants()` | Reference allele/genotype selection must prefer the reference variant when it exists |
| `run_simple_post_processor()` | `SimplePostProcessModule` output selection | Final simple-mode output filtering is where Java-compatible low-frequency row decisions belong |

## Known Parity Traps
- Do not drop whole single-sample positions in `run_to_vars_builder()` just because `maxfreq <= conf.freq`. Java still emits some low-frequency SNVs from those loci after downstream filtering.
- Simple-mode output cannot reuse the full `is_good_var()` gate unchanged. The output layer must allow Java-retained low-frequency rows that survive the other quality checks.
- Genotype shaping for non-reference variants must keep the reference allele as `genotype1` whenever a reference variant exists, even when the top alternate frequency is below `conf.freq`.

## Divergences from Java
- Rust centralizes simple, somatic, and amplicon orchestration in one large pipeline module instead of Java's separate mode classes.
- Snapshot writers (`VARDICT_*_JSONL`) are Rust-only diagnostics used to localize parity loss between CIGAR parsing, realignment, structural processing, and to-vars materialization.

## Cross-Module Dependencies
- Calls: `cigar_parser.rs`, `variant_realigner.rs`, `structural_variants_processor.rs`, `output_variant.rs`, shared-reference helpers.
- Called by: `parallel_pipeline.rs`, CLI execution paths in `src/bin/vardict.rs`.