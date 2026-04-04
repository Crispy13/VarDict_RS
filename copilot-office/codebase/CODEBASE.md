# VarDict-rs Codebase Overview

> **This file is a legacy overview.** The structured, per-module codebase cache is now at:
> - **Rust**: [rust/VarDict-rs-CODEBASE.md](rust/VarDict-rs-CODEBASE.md)
> - **Java**: [java/VarDictJava-CODEBASE.md](java/VarDictJava-CODEBASE.md)
>
> Agents should read those indexes instead. This file is kept for backward compatibility.

## Project
Rust port of [VarDictJava](https://github.com/AstraZeneca-NGS/VarDictJava) — a variant discovery tool for next-generation sequencing (NGS) data. Goal: **byte-identical output parity** with the Java implementation.

## Build
- **Toolchain**: Rust 2024 edition, `rust_build_env` conda environment
- **Profiles**: `debug-release` (optimized + debug info) for tests; `release` for production
- **Key deps**: `rust-htslib` (BAM I/O), `indexmap` (insertion-ordered maps), `rayon` (parallelism), `rustc-hash` (fast hashing), `statrs` (statistics), `crackle-kit` (tracing/logging)

```bash
conda activate rust_build_env
export LIBCLANG_PATH=$CONDA_PREFIX/lib
cargo build --profile debug-release
cargo test --profile debug-release
```

## Architecture

### Pipeline Flow
```
CLI (bin/vardict.rs)
  → GlobalReadOnlyScope (scopedata/) — static configuration singleton
  → BamReader (data/bam_reader.rs) — BAM record iteration
  → CigarModifier (mods/cigar_modifier.rs) — CIGAR normalization / realignment
  → CigarParser (mods/cigar_parser.rs) — variant extraction from CIGAR strings
  → VariantRealigner (mods/variant_realigner.rs) — local realignment of variants
  → StructuralVariantsProcessor (mods/structural_variants_processor.rs) — SV detection
  → ToVarsBuilder (mods/to_vars_builder.rs) — variant statistics & filtering
  → OutputVariant (mods/output_variant.rs) — tab-delimited output formatting
  → ParallelPipeline (mods/parallel_pipeline.rs) — rayon-based region parallelism
  → VarDictPipeline (mods/vardict_pipeline.rs) — full pipeline orchestration (7.4K LOC)
```

### Module Size (LOC)
| Module | LOC | Risk | Notes |
|--------|-----|------|-------|
| vardict_pipeline.rs | 7,447 | HIGH | Full pipeline orchestration, all modes |
| variant_realigner.rs | 5,271 | HIGH | Local realignment, position adjustment |
| cigar_parser.rs | 5,004 | HIGH | Core SNP/indel detection from CIGAR |
| structural_variants_processor.rs | 4,702 | HIGH | DEL/DUP/INV from discordant pairs |
| output_variant.rs | 2,677 | MEDIUM | Column formatting, float precision |
| cigar_modifier.rs | 1,984 | MEDIUM | CIGAR normalization edge cases |
| bin/vardict.rs | 1,721 | LOW | CLI arg parsing |
| to_vars_builder.rs | 1,414 | MEDIUM | MSI, genotype, variant filtering |
| parallel_pipeline.rs | 1,057 | LOW | rayon parallelism |
| pipeline.rs | 663 | LOW | Pipeline config/wiring |

### Data Types
- `Configuration` (conf.rs) — all CLI/runtime parameters
- `GlobalReadOnlyScope` (scopedata/) — static singleton for config, chr_lens, BAM paths
- `Variant` (variants/variants.rs) — variant record with all fields
- `VariantMap` = `HashMap<String, Variant, FxBuildHasher>` — position-keyed variant store
- `ModifiedCigar` (data.rs) — normalized CIGAR with adjusted positions
- `Region` (data/region.rs) — genomic interval

### Modes
- **Simple**: Single-sample variant calling (most test coverage)
- **Somatic**: Tumor/Normal paired analysis
- **Amplicon**: Targeted sequencing with amplicon boundaries

### Key Patterns
- `IndexMap` for insertion-ordered maps (Java `LinkedHashMap` parity)
- `FxBuildHasher` (`rustc-hash`) for fast internal hashing
- `SmallVec<[u8; 32]>` for short genomic sequences
- `OnceLock<GlobalReadOnlyScope>` for static configuration (matches Java's `GlobalReadOnlyScope.instance()`)

### Divergences from Java (Optimizations)

#### VecMap (src/utils/vec_map.rs) — INNER VARIANT MAPS
**Java**: `VariationMap<String, Variation>` extends `LinkedHashMap` (inner maps at each genomic position).
**Rust**: `VecMap<VarDesc, RawVariant>` backed by `Vec<(K,V)>` — NOT a faithful port.

This is a **performance/memory optimization**, not Java logic. It assumes inner maps hold 1–4 entries (typical for diploid biology). Rust `hashbrown::HashMap` pre-allocates ~388 bytes per instance vs VecMap's ~24 bytes, saving ~1.9 GB across ~5M maps.

**Risks**: Edge cases (homopolymer runs, >1000x coverage, multi-allelic sites) may produce 10+ entries per map, degrading linear scan performance. **Needs more benchmark tests. Revert to HashMap if production profiling shows regressions.** Type aliases to change back: `RawVarMap` and `CountMap` in `vardict_pipeline.rs`.

## Test Infrastructure

### Parity Tests (tests/integration_test.rs)
- **Manifest-driven**: `tests/parity_case_manifest.csv` lists 99 test cases
- **Java reference data**: `VarDictJava/testdata/integrationtestcases/` contains expected outputs
- **Test groups**: Integration (74 cases), Tier1 (10), SV-Core (38), Somatic (16), Amplicon (3), Unique-Mode (2), Realigner-Complex (4), Splicing (1)
- **Environment vars**: `VARDICT_RUN_NOW_SIMPLE_LIMIT`, `VARDICT_PARITY_VERBOSE`

### Fixture Tests
- `tests/cigar_parser_fixture_test.rs` — JSONL-based CigarParser snapshots
- `tests/variant_realigner_fixture_test.rs` — realigner snapshots
- `tests/structural_variants_fixture_test.rs` — SV processor snapshots
- `tests/tovars_fixture_test.rs` — ToVarsBuilder snapshots

### Large BAM Tests
- NA12878 low-coverage whole-genome (`testdata/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam`)
- hs37d5 reference (`testdata/hs37d5.fa`)
- Test region: chr20 168500-168800 and broader ranges
