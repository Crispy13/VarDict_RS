# ScopeData

**Source**: `src/scopedata/` (`global_read_only_scope.rs`, `mod.rs`)
**LOC**: ~150
**Java counterpart**: `data/scopedata/*.java` (7 classes) → [Java cache](../java/ScopeData.md)
**Status**: complete

## Overview

The `scopedata` module holds the **global singleton configuration and sample metadata** shared across all pipeline stages. Unlike Java's 7 scopedata classes (`Scope<T>`, `GlobalReadOnlyScope`, `InitialData`, `VariationData`, `RealignedVariationData`, `AlignedVarsData`, `CombineAnalysisData`), Rust consolidates per-region mutable state into stage-specific structs in `src/mods/` and retains **only `GlobalReadOnlyScope`** here. This is a significant architectural divergence: Rust uses direct data flow with move semantics rather than a generic `Scope<T>` wrapper.

## Public API

| Function/Method | Purpose |
|----------------|---------|
| `GlobalReadOnlyScope` (struct) | Singleton: CLI config, chromosome lengths, sample paths, adapter sequences. Immutable after init. |
| `instance()` | Returns `&'static GlobalReadOnlyScope`. Panics if not initialized. |
| `instance_arc()` | Returns `Arc<GlobalReadOnlyScope>` (created once, memoized). |
| `INSTANCE` (OnceLock) | Underlying once-initialized storage. |

**Key fields**: `chr_lens: HashMap<String, usize>`, `bam_paths: Vec<String>`, `conf: Configuration`, `adaptor_forward/reverse: HashMap<String, usize>`, `amplicon_based_calling: Option<String>`

## Java Correspondence

| Rust | Java | Notes |
|------|------|-------|
| `GlobalReadOnlyScope` | `GlobalReadOnlyScope` | 1:1 port |
| *(not ported)* | `Scope<T>` | Generic wrapper replaced by named stage data structs |
| `CigarParserOutput` (`vardict_pipeline.rs`) | `InitialData` | Merged into pipeline stage output |
| *(implicit in fn params)* | `VariationData` | Fields passed separately, no wrapper |
| `RealignedVariationData` (`structural_variants_processor.rs`) | `RealignedVariationData` | Lives outside scopedata module |
| `AlignedVarsData` (`vardict_pipeline.rs`) | `AlignedVarsData` | Lives outside scopedata module |
| *(not ported)* | `CombineAnalysisData` | Absorbed into pipeline flow |

## Known Parity Traps

1. **No `Scope<T>` wrapper**: Java passes common fields (region, reference, bam path, splice set) through `Scope<T>`. Rust passes these separately. Code relying on `scope.getRegion()` must have an explicit region parameter.

2. **`GlobalReadOnlyScope` test init**: Java's `instance()` never returns null. Rust's panics if not initialized. `OnceLock` doesn't support re-initialization — tests need special handling.

3. **Stage data location divergence**: `RealignedVariationData` lives in `structural_variants_processor.rs`, not `scopedata`. Cross-module imports must reach into the correct module.

4. **`AlignedVarsData` capacity tracking**: Rust includes `aligned_variants_java_capacity: usize` to simulate Java's HashMap iteration order — a Rust-only parity field with no Java counterpart.

5. **SV buckets as separate Vecs**: Rust's `RealignedVariationData` has `svfdel`, `svrdel`, `svfdup`, `svrdup`, etc. as separate fields vs. Java's `SVStructures` maps. Iteration order must be preserved.

6. **SOFTP2SV back-reference removed**: Java's `RealignedVariationData.previousScope` creates a reference cycle. Rust eliminates this — code needing previous-stage data must fetch from external context.

## Divergences from Java

| Aspect | Java | Rust | Impact |
|--------|------|------|--------|
| Generic wrapper | `Scope<T>` carries 7 common fields | Direct stage data structs, params passed separately | Architectural — no output impact |
| Data type location | All in `data/scopedata/` package | Stage structs co-located with processor modules | Import paths differ |
| Module exports | 7 classes | Only `GlobalReadOnlyScope` | Reduced coupling |
| Reference wrapping | Normal class | `Arc<>` wrapped for O(1) clone | Performance optimization |
| Singleton pattern | Java static field | `OnceLock<GlobalReadOnlyScope>` | Thread-safe, one-shot init |

## Cross-Module Dependencies

| Direction | Dependency | Purpose |
|-----------|-----------|---------|
| ← `conf` | `GlobalReadOnlyScope` embeds `Configuration` | CLI config access |
| → All pipeline modules | `instance()` accessor | Global config, chr_lens, adaptors |
| External | `structural_variants_processor.rs` defines `RealignedVariationData` | Stage data lives outside scopedata |
| External | `vardict_pipeline.rs` defines `AlignedVarsData`, `CigarParserOutput` | Stage data lives outside scopedata |
