# Utils / VecMap

**Source**: `src/utils.rs`, `src/utils/aligner.rs`, `src/utils/vec_map.rs`
**Java counterpart**: `Utils.java`, `VariationMap.java` → [Java cache](../java/Utils.md)
**Status**: complete

## Overview

The utils module provides low-level utility functions and data structures used across the entire pipeline. It contains two submodules — **`aligner`** (BAM aligner type selection: BWA vs STAR, controlling which SAM tag is used for edit distance) and **`vec_map`** (a `Vec<(K,V)>`-backed map optimized for small inner variant maps). The module root exports `round_half_even()` for HALF_EVEN banker's rounding (used in 40+ output column computations), `print_exception_and_continue()` for error handling with exception count limits, and extension traits (`SliceExt`, `SliceExt2`, `BytesExt`) for negative-index slicing and UTF-8 conversion. All root functions are faithful Java ports; VecMap is an intentional optimization divergence.

## Public API

| Function/Method | Purpose |
|----------------|---------|
| `round_half_even(pattern, value)` | HALF_EVEN banker's rounding via format→parse round-trip; matches Java's `DecimalFormat` |
| `print_exception_and_continue(...)` | Exception logging with count limits; increments exception counter |
| `SliceExt::get_or_err(index)` | Bounds-checked indexing returning `Result` |
| `SliceExt2::get_with_int(range)` | Signed-integer indexing (Perl-compatible negative index support) |
| `BytesExt::try_as_str()` | UTF-8 validation on byte slices |
| `Aligner::BWA / Star` | Enum selecting NM vs nM SAM tag for edit distance |
| `VecMap<K,V>` | Generic Vec-backed map with Entry API (Occupied/Vacant) |
| `VecMap::entry(key)` | Entry API matching `HashMap::entry()` contract |

## Java Correspondence

| Rust | Java | Notes |
|------|------|-------|
| `round_half_even()` | `Utils.roundHalfEven()` | Exact — format→parse round-trip matches `DecimalFormat("0.0000")` |
| `print_exception_and_continue()` | `Utils.printExceptionAndContinue()` | Exact — logs via crackle_kit tracing; increments counter |
| `SliceExt`, `SliceExt2` | `Utils.substr()`, `charAt()` | Ported — Perl-compatible negative index semantics |
| `VecMap<VarDesc, Variant>` | `VariationMap<String, Variation>` (LinkedHashMap) | **Optimization divergence** — Vec-backed linear scan vs hashing |
| `Aligner` enum | `Utils.getAligner()` | Ported — returns SAM tag name for edit distance |

## Known Parity Traps

1. **Entry Name Collision** (repo: `vecmap_entry_name_collision_and_inner_count_alias_20260327.md`): When switching from HashMap to VecMap, code had to change from `std::collections::hash_map::Entry` to `crate::utils::vec_map::Entry`. Fix: Type alias `InnerMapEntry` in prelude.

2. **Feature Toggle: HashMap Fallback** (repo: `vecmap_feature_toggle_requires_default_ctor_and_cigar_parser_aliases_20260328.md`): VecMap constructors not compatible with `HashMap::new()` when swapping backends. Mitigation: Use `Default::default()` and `InnerMapEntry` trait abstraction.

3. **Size Assumption** (repo: `vecmap_optimization_divergence_not_java_logic_20260327.md`): VecMap assumes 1–4 entries per position (typical diploid). Edge cases with >15 entries would be slower than HashMap. No stress tests exist.

4. **Inner HashMap Replacement Context** (repo: `vecmap_inner_hashmap_replacement_20260327.md`): Memory savings of ~1.9 GB on 5 MB shard. VecMap result: 0.56x Java RSS while maintaining byte-identical output.

## Divergences from Java

| Aspect | Java | Rust | Reason | Impact |
|--------|------|------|--------|--------|
| Inner map backing | `LinkedHashMap<String, Variation>` | `VecMap<VarDesc, Variant>` (Vec-backed) | Memory: 388B HashMap overhead × 5M instances = wasted | Zero output impact; 3.5x memory reduction |
| Lookup algorithm | O(1) hashing | O(n) linear scan (n ∈ [1,4]) | Faster for tiny sets due to cache locality | Faster typical; regresses if >15 entries |
| Entry management | HashMap entry API | VecMap entry API (same contract) | Vec positional lookup | Contract identical; parity-safe |

All root utils functions (`round_half_even`, `print_exception_and_continue`, extension traits, `Aligner`) are Java-faithful with NO divergence.

## Cross-Module Dependencies

**Called by:**
| Caller | Functions/Types Used | Purpose |
|--------|---------------------|---------|
| CigarParser | `VecMap` (via `InnerMap` alias); `round_half_even()` | Per-position variant accumulation |
| VariantRealigner | `VecMap` (via `InnerMap` alias); Entry API | Realignment result storage |
| StructuralVariantsProcessor | `VecMap`; `print_exception_and_continue()` | SV evidence accumulation; error handling |
| ToVarsBuilder | `round_half_even()` (40+ sites) | Output column formatting |
| OutputVariant | `round_half_even()` | TSV column value formatting |
| VarDictPipeline | `round_half_even()`; `InnerMap` alias | Pipeline-level computations |

**Calls:** No significant outbound dependencies (leaf module).
