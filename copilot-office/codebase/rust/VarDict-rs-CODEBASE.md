# VarDict-rs Codebase Cache

> Progressive knowledge base maintained by agents during Rust implementation work.
> Unlike the Java cache, Rust code evolves — focus on **architecture, Java correspondence, and parity traps** rather than method-level detail.
> Each module has its own file to keep agent context small.

## Build

```bash
conda activate rust_build_env
export LIBCLANG_PATH=$CONDA_PREFIX/lib
cargo build --profile debug-release
cargo test --profile debug-release
```

## Pipeline Flow

```
CLI (bin/vardict.rs)
  → GlobalReadOnlyScope (scopedata/) — static configuration singleton
  → BamReader (data/bam_reader.rs) — BAM record iteration
  → CigarModifier (mods/cigar_modifier.rs) — CIGAR normalization
  → CigarParser (mods/cigar_parser.rs) — variant extraction from CIGAR strings
  → VariantRealigner (mods/variant_realigner.rs) — local realignment
  → StructuralVariantsProcessor (mods/structural_variants_processor.rs) — SV detection
  → ToVarsBuilder (mods/to_vars_builder.rs) — variant statistics & filtering
  → OutputVariant (mods/output_variant.rs) — tab-delimited output formatting
  → ParallelPipeline (mods/parallel_pipeline.rs) — rayon-based region parallelism
  → VarDictPipeline (mods/vardict_pipeline.rs) — full pipeline orchestration
```

## Module Index

| Module | Rust File | LOC | Risk | Java Counterpart | Cache File | Status |
|--------|-----------|-----|------|------------------|------------|--------|
| VarDictPipeline | `src/mods/vardict_pipeline.rs` | 7,447 | HIGH | `modes/*`, `SAMFileParser` | [vardict_pipeline.md](vardict_pipeline.md) | partial |
| VariantRealigner | `src/mods/variant_realigner.rs` | 5,271 | HIGH | `VariationRealigner.java` | [variant_realigner.md](variant_realigner.md) | not started |
| CigarParser | `src/mods/cigar_parser.rs` | 5,004 | HIGH | `CigarParser.java` | [cigar_parser.md](cigar_parser.md) | not started |
| StructuralVariants | `src/mods/structural_variants_processor.rs` | 4,702 | HIGH | `StructuralVariantsProcessor.java` | [structural_variants_processor.md](structural_variants_processor.md) | partial |
| OutputVariant | `src/mods/output_variant.rs` | 2,677 | MEDIUM | `printers/*.java` | [output_variant.md](output_variant.md) | not started |
| CigarModifier | `src/mods/cigar_modifier.rs` | 1,984 | MEDIUM | `CigarModifier.java` | [cigar_modifier.md](cigar_modifier.md) | not started |
| CLI | `src/bin/vardict.rs` | 1,721 | LOW | `CmdParser.java` | [cli.md](cli.md) | not started |
| ToVarsBuilder | `src/mods/to_vars_builder.rs` | 1,414 | MEDIUM | `ToVarsBuilder.java` | [to_vars_builder.md](to_vars_builder.md) | not started |
| ParallelPipeline | `src/mods/parallel_pipeline.rs` | 1,057 | LOW | threading model | [parallel_pipeline.md](parallel_pipeline.md) | not started |
| Configuration | `src/conf.rs` | — | LOW | `Configuration.java` | [configuration.md](configuration.md) | not started |
| Data Types | `src/data.rs`, `src/data/` | — | MEDIUM | `data/*.java` | [data_types.md](data_types.md) | not started |
| ScopeData | `src/scopedata/` | — | LOW | `scopedata/*.java` | [scopedata.md](scopedata.md) | not started |
| Variants | `src/variants/` | — | MEDIUM | `variations/*.java` | [variants.md](variants.md) | not started |
| Utils / VecMap | `src/utils/` | — | MEDIUM | `Utils.java`, `VariationMap.java` | [utils.md](utils.md) | not started |

## Cross-Cutting Concerns

### Key Patterns
- `IndexMap` for insertion-ordered maps (Java `LinkedHashMap` parity)
- `FxBuildHasher` (`rustc-hash`) for fast internal hashing
- `SmallVec<[u8; 32]>` for short genomic sequences
- `OnceLock<GlobalReadOnlyScope>` for static configuration
- `VecMap<K,V>` backed by `Vec<(K,V)>` for small inner variant maps (perf optimization, not Java logic)

### Modes
- **Simple**: Single-sample variant calling (most test coverage)
- **Somatic**: Tumor/Normal paired analysis
- **Amplicon**: Targeted sequencing with amplicon boundaries

### Test Infrastructure
- **Manifest-driven parity tests**: `tests/parity_case_manifest.csv` (99 cases)
- **Fixture tests**: JSONL-based snapshots for CigarParser, VariantRealigner, SV processor, ToVarsBuilder
- **Large BAM tests**: NA12878 low-coverage whole-genome + hs37d5 reference

## Per-Module File Template

When creating a module file for the first time, use this structure:

```markdown
# ModuleName

**Source**: `src/path/to/module.rs`
**LOC**: N
**Java counterpart**: `modules/JavaClass.java` → [Java cache](../java/JavaClass.md)
**Status**: partial | complete

## Overview
One-paragraph summary of the module's role in the pipeline.

## Public API
| Function/Method | Purpose |
|----------------|---------|
| fn_name() | one-liner |

## Java Correspondence
| Rust | Java | Notes |
|------|------|-------|
| rust_fn() | JavaClass.javaMethod() | any divergence |

## Known Parity Traps
- Specific Java behaviors that required non-idiomatic Rust

## Divergences from Java
- Performance optimizations, algorithmic changes, or structural differences
- Why the divergence exists and what to watch for

## Cross-Module Dependencies
- Calls: which modules this depends on
- Called by: which modules invoke this
```
