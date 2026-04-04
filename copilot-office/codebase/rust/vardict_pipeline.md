# vardict_pipeline

**Source**: `src/mods/vardict_pipeline.rs`
**LOC**: ~7,447
**Java counterpart**: Multiple (Modes, SimpleMode, ToVarsBuilder) → [Java cache](../java/Modes.md)
**Status**: complete
**Last verified**: 2026-04-04

## Overview

The `vardict_pipeline` module is the **main orchestrator and postprocessor** for the VarDict variant-calling pipeline. It wires the four core processing stages (`CigarParser`, `VariantRealigner`, `StructuralVariantsProcessor`, `ToVarsBuilder`) into a coherent workflow, then applies Java-compatible output filtering and formatting. It implements **three distinct postprocessing paths**: simple mode (single-sample), somatic mode (tumor/normal), and amplicon mode (targeted sequencing). The module also owns all Rust-specific diagnostic snapshots (JSONL), pileup handling, and quality-filter logic that determines which variants are emitted to stdout.

## Public API

| Function/Method | Purpose |
|---|---|
| `VarDictPipeline::new(sample_name)` | Constructor; initializes with default settings |
| `with_min_frequency(f64)` | Set minimum variant frequency threshold |
| `with_min_base_quality(f64)` | Set minimum base quality filter |
| `with_min_mapping_quality(u8)` | Set minimum mapping quality filter |
| `with_pileup(bool)` | Enable pileup mode (emits reference call rows) |
| `process_region_from_bam()` | Main entry: BAM → full pipeline → TSV lines |
| `process_region_from_bam_streaming<W: Write>()` | Streaming variant: output directly to writer |
| `process_region_from_cached_records()` | Entry using pre-fetched `Vec<Record>` |
| `process_region_splicing_from_bam()` | Splicing-mode: 4-column intron count table |
| `process_region_to_aligned_vars_from_bam()` | Lower-level: produces `AlignedVarsData` |
| `process_region_to_aligned_vars_from_bam_with_paths()` | With explicit BAM paths for realignment |
| `process_region_to_aligned_vars_from_bam_paths()` | Reads BAM from file paths |
| `process_region<I: Iterator>()` | Generic iterator-based variant |
| `collect_filtered_records()` | Filters BAM records by region, MAPQ, duplicates |
| `run_amplicon_post_processor()` | Amplicon mode postprocessor |
| `run_somatic_post_processor()` | Somatic mode postprocessor |
| `run_somatic_post_processor_with_combine_lookup()` | Somatic with re-entrant combine analysis |

## Java Correspondence

| Rust Function | Java Equivalent | Notes |
|---|---|---|
| `VarDictPipeline` struct | `AbstractMode` + `SimpleMode/SomaticMode/AmpliconMode` | Rust centralizes all modes; Java uses inheritance |
| `process_region_from_bam()` | `SimpleMode.processBamInPipeline()` | Simple mode full pipeline |
| `run_amplicon_post_processor()` | `AmpliconPostProcessModule` | Amplicon output filtering |
| `run_somatic_post_processor()` | `SomaticPostProcessModule` | Tumor/normal comparison |
| `run_to_vars_builder()` | `ToVarsBuilder.process()` | Raw variants → Variant objects |
| `run_simple_post_processor()` | `SimplePostProcessModule` | Quality filtering + output |
| `is_good_var()` | `Variant.isGoodVar()` | Parity filter: freq, minr, meanpos, meanq, etc. |
| `collect_reference_variants()` | `ToVarsBuilder.collectReferenceVariants()` | Reference allele acquisition |
| `create_variant_records()` | `ToVarsBuilder.createVariantRecords()` | Statistics aggregation |
| `sort_variants()` | `Variant` comparators | meanQuality * count desc, then desc asc |
| `convert_desc_string_to_alleles()` | Description string parsing | Maps `-3`, `+AT` to ref/alt sequences |

## Known Parity Traps

### Trap 1: Realigner Must Always Run
CigarParser output is passed through VariantRealigner even for simple mode. Skipping breaks parity.

### Trap 2: Reference Variant Materialization
Reference variants must be created at every position with coverage, even without alternates. Java ToVarsBuilder depends on this.

### Trap 3: Low-Frequency Position Retention
Do NOT drop positions where `frequency < conf.freq`. Java still emits low-frequency SNVs. Frequency check is a postprocessing gate, not early exit.

### Trap 4: SV Marker Position Filtering
Positions outside region bounds are emitted if they carry SV markers. Intentional Java behavior for discordant-pair SVs spanning boundaries.

### Trap 5: Reference Variant Genotyping
When ref variant frequency < conf.freq and alternates exist, genotype1 = reference allele (not top alternate). Simple-mode output cannot reuse the full `is_good_var()` gate unchanged — must allow Java-retained low-frequency rows surviving other quality checks.

### Trap 6: Variant Sort Order
Two-phase sort: (1) meanQuality * variantCount descending, (2) description ascending.

### Trap 7: Large Deletion RightSeq
Uses `historical_del_rightseq_variants` map + `historical_reference_windows` ranges. Empty rightseq correct if variant outside historical windows.

### Trap 8: Insertion-Only Loci
Java does NOT materialize unanchored insertions. Check `is_same_variation_on_ref()` and skip if no non-insertion variant.

### Trap 9: HashMap Iteration Order
Multiple helpers simulate Java 8 HashMap bucket order: `java_hashmap_bucket_index`, `java_hashmap_string_iteration_order`, `java_hashmap_iteration_order_with_capacity`.

## Divergences from Java

### Architecture
1. **Centralized module**: All three modes in one struct vs Java's three mode classes with inheritance.
2. **No global singleton state**: Explicit `SharedReferenceHandle` and `GlobalReadOnlyScope` references instead of Java's `getMode()`.

### Diagnostics
- **JSONL snapshot writers**: Rust-only (env vars: `VARDICT_CIGAR_PARSER_JSONL`, `VARDICT_VARIANT_REALIGNER_JSONL`, `VARDICT_STRUCTURAL_VARIANTS_JSONL`, `VARDICT_TOVARS_JSONL`).

### Reference Handling
- **SharedReferenceHandle**: Arc-wrapped shared reference for cross-thread access vs Java's per-thread `ReferenceResource`.

### Memory
- **`trim_process_allocator()`**: Calls `malloc_trim(0)` on Linux; Java GC handles automatically.
- **Strand bias**: `StrandBiasFlag` struct vs Java's bitwise flags.

## Key Internal Functions

### Record Preprocessing
| Function | Purpose |
|---|---|
| `collect_filtered_records()` | Filter BAM by region, MAPQ, duplicates |
| `passes_preprocess()` | Individual record filter |

### CIGAR Parsing
| Function | Purpose |
|---|---|
| `run_cigar_parser_from_bam()` | BAM → CigarParser |
| `build_cigar_output()` | Package CigarParser state into output struct |

### Realignment & SV
| Function | Purpose |
|---|---|
| `run_variant_realigner_and_sv_processor()` | Chain realigner + SV processor |
| `run_partial_cigar_for_bams()` | Re-entrant partial CIGAR for SV extension |

### To-Vars Building
| Function | Purpose |
|---|---|
| `run_to_vars_builder()` | Main materialization: raw → Variant objects |
| `create_variant_records()` | Aggregate raw counts per description |
| `create_insertion_records()` | Special insertion handling |
| `collect_reference_variants()` | Build reference variant; genotype selection |
| `convert_raw_variant()` | Raw counts → Variant fields |
| `convert_desc_string_to_alleles()` | Short descriptions → ref/alt sequences |
| `detect_microsatellite()` | MSI unit and repeat count detection |

### Postprocessing
| Function | Purpose |
|---|---|
| `run_simple_post_processor()` | Simple mode: filter + emit TSV |
| `run_somatic_post_processor_with_combine_lookup()` | Somatic: tumor/normal comparison |
| `is_good_var_with_type()` | Full Java parity quality gate |
| `determinate_somatic_type()` | Classify: StrongSomatic, LikelySomatic, Germline, etc. |

### Java HashMap Simulation
| Function | Purpose |
|---|---|
| `java_hashmap_bucket_index()` | Simulates Java HashMap bucket for i64 keys |
| `java_string_hash()` | Replicates Java String.hashCode() |
| `java_hashmap_string_iteration_order()` | Reorders entries to match Java iteration |

## Cross-Module Dependencies

### Calls Into
| Module | Purpose |
|--------|---------|
| `cigar_parser` | CigarParser construction, processing, accessor methods |
| `variant_realigner` | VariantRealigner construction, all realignment methods |
| `structural_variants_processor` | SV processor construction, process, accessors |
| `to_vars_builder` | Types: VarType, Variant, Vars, check_strand_bias(), determine_genotype() |
| `output_variant` | SimpleOutputVariant, SomaticOutputVariant, AmpliconOutputVariant |
| `data/bam_reader` | BamReader for BAM I/O |
| `data/reference` | Reference, seed maps |
| `data/shared_reference` | SharedReferenceHandle |
| `scopedata/global_read_only_scope` | Configuration, chromosome info, BAM paths |

### Called By
| Module | Usage |
|--------|-------|
| `parallel_pipeline` | High-level parallel region orchestrator |
| `src/bin/vardict.rs` | CLI entry point |
| Tests | Direct unit test calls |

### Key Data Types
- `CigarParserOutput` — raw variant maps, soft clips, counts
- `RealignedOutput` — realigned variants, historical windows, rightseq
- `AlignedVarsData` — materialized Variant objects per position
- `Vars` — variants at a single position
- `Variant` — fully populated statistics