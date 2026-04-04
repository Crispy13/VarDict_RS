# CLI

**Source**: `src/bin/vardict.rs`
**LOC**: 1,721
**Java counterpart**: `VarDictLauncher.java` + `CmdParser.java` + `Configuration.java` → [Java cache](../java/Configuration.md)
**Status**: complete

## Overview

Single-file Rust CLI entry point replacing VarDictJava's multi-class launcher stack. Delegates CLI parsing to **clap** (typed `Args` struct with ~50 fields), BED region loading/normalization with BAM target awareness, reference genome initialization via `load_shared_reference_chroms()` (loads only needed chromosomes), and pipeline mode dispatch (Simple/Amplicon/Somatic/Splicing). Parallelizes regions using `ParallelPipeline` with ordered output via `OrderedStreamConsumer`.

## Public API

| Function/Method | Purpose |
|----------------|---------|
| `Args` (struct) | ~50 clap-derived fields covering all VarDict CLI params |
| `main()` | Validates inputs, loads reference, determines mode, invokes `run_variant_calling()` |
| `parse_args()` | Legacy flag normalization (`-UN` → long form) then `clap::Parser::parse_from()` |
| `run_variant_calling()` | Init GlobalReadOnlyScope, load SharedReference, select execution mode, dispatch |
| `ExecutionMode` (enum) | `Simple | Amplicon | Somatic | Splicing` |

## Key Internal Functions

| Function | Purpose |
|----------|---------|
| `normalize_legacy_cli_args()` | Rewrite `-UN` and `-fisher` to long-form before clap |
| `resolve_execution_mode()` | Priority: splicing > amplicon > somatic > simple |
| `parse_bam_inputs()` | Split `-b` on `\|` for tumor + normal BAM paths |
| `validate_bam_with_index()` | Check `.bam` and `.bam.bai` exist |
| `infer_sample_name_from_bam()` | Extract sample name via two regex patterns (sorted first, then fallback) |
| `get_regions()` | Load from `-R` string, BED, or default |
| `extend_region()` | Extend by `-x` nucleotides (saturating arithmetic) |
| `parse_region_string()` | Parse `chr:start[-end]`, apply zero-based conversion |
| `parse_bed_file()` | Read BED, detect amplicon (8-col), group regions |
| `parse_amplicon_region_groups()` | Group amplicon regions by chr and insert overlap |
| `normalize_region_chrom()` | Validate chr names against BAM targets |

## Java Correspondence

| Rust | Java | Notes |
|------|------|-------|
| `main()` | `VarDictLauncher.main()` | Inlined into single file |
| `clap` via `Args` derive macro | `CmdParser.parseParams()` + `parseCmd()` | Typed arguments vs. string map |
| `Args` struct `#[arg]` attributes | `CmdParser.buildOptions()` | Clap handles validation |
| `PipelineConfig` (builder) + `GlobalReadOnlyScope.conf` | `Configuration` (POJO) | Split for modularity |
| `INSTANCE.set(scope)` | `GlobalReadOnlyScope.init()` | Thread-safe `OnceLock` init |
| `parse_bed_file()` | `VarDictLauncher.readBedFile()` | Auto-detects amplicon from 8-col BED |
| `infer_sample_name_from_bam()` | `VarDictLauncher.getSampleNames()` | Same regex patterns |

## Known Parity Traps

1. **BED amplicon auto-detection**: 8-col BED auto-sets `zero_based = true`. Must match Java's implicit behavior.
2. **Sample name inference**: Two regex patterns tried in order (sorted pattern first, then fallback). Order matters for parity.
3. **Zero-based coordinate conversion**: `-z 1` flag for single regions. BED 8-col forces zero-based. Conversion: `(start + 1, end)`.
4. **SAM filter hex parsing**: Must handle `0x504`, `0X504`, and decimal formats.
5. **Pileup mode side-effects**: `-p` overrides `freq = -1.0` and `minr = 0` unconditionally (after other flags).
6. **Legacy option normalization**: `-UN` and `-fisher` rewritten to long forms before clap parsing.
7. **Chromosome lengths source**: Rust reads from FASTA `.fai`; Java reads from BAM header. If they disagree, results may differ.
8. **Region extension saturating arithmetic**: `saturating_sub()/saturating_add()` prevents underflow; more conservative than Java's silent underflow.

## Divergences from Java

| Aspect | Java | Rust | Risk |
|--------|------|------|------|
| CLI framework | Apache Commons CLI (string map) | clap (typed struct) | LOW |
| Config storage | Single monolithic `Configuration` | Split: `PipelineConfig` + `GlobalReadOnlyScope.conf` | MED |
| BAM I/O strategy | Upfront full BAM scan | Lazy per-region BAM open | LOW |
| Chromosome lengths | From BAM header | From FASTA `.fai` index | MED |
| Parallelism | Java `ExecutorService` | rayon + crossbeam_channel + `OrderedStreamConsumer` | LOW |
| Allocator | JVM heap | mimalloc (conditional) | LOW |

## Cross-Module Dependencies

**Calls**:
| Module | Used For |
|--------|----------|
| `conf::Configuration` | Built from `Args`, stored in `GlobalReadOnlyScope` |
| `data::bam_reader::BamReader` | BAM validation + target name extraction |
| `data::region::Region` | Parsed regions passed to pipeline |
| `data::shared_reference::load_shared_reference_chroms()` | Reference FASTA loading |
| `mods::parallel_pipeline::ParallelPipeline` | Multi-threaded simple mode |
| `mods::vardict_pipeline::VarDictPipeline` | Amplicon/somatic/splicing modes |
| `mods::output_variant` | Header generation, output formatting |
| `scopedata::global_read_only_scope` | Singleton init |

**Called by**: Top-level entry point — not called by other modules.
