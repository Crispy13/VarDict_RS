# Variants

**Source**: `src/variants/`, `src/variants/variants.rs`, `src/variants/var_utils.rs`
**LOC**: ~600 (estimated across module files)
**Java counterpart**: `variations/Variation.java`, `variations/Sclip.java`, `variations/Mate.java`, `variations/VariationUtils.java` → [Java cache](../java/Variations.md)
**Status**: complete

## Overview

The `variants` module defines the core data structures flowing through the entire VarDict pipeline: `Variant` (read-level mutation accumulator), `VarDesc` (type-safe variant description keys), `SoftClip` (soft-clipped read data with consensus caching), `Mate` (mate-pair info for SVs), `StructuralVariantCounts` (SV evidence counters), and `InsOrDelLen` (deletion/insertion length discriminator). The sub-module `var_utils` provides type-preserving map accessors and consensus sequence computation. These types are cross-cutting: `Variant` and `VarDesc` flow through CigarParser → VariationRealigner → StructuralVariantsProcessor → ToVarsBuilder, while `SoftClip` carries consensus and SV evidence between stages.

## Public API

| Function/Type | Purpose |
|---------------|---------|
| `Variant` (struct) | Read-level accumulator: depth, strand counts, mean pos/qual/mapq, nm, pstd/qstd flags |
| `VarDesc` (enum: SNV/Ins/Del/Complex/Raw) | Type-safe variant description key; replaces Java's plain String keys |
| `VarDesc::to_key_string()` | Reconstruct Java-compatible key string (byte-for-byte parity required) |
| `VarDesc::cmp_as_key_string()` | Zero-allocation key comparison via KeySegments decomposition (46 test cases) |
| `VarDesc::key_equals(s: &str)` | Zero-alloc string comparison against VarDesc |
| `VarDesc::ref_allele() / alt_allele()` | Output formatting helpers |
| `SoftClip` (struct) | Soft-clip data: contains `Variant`, position maps (`nt`, `seq` as BTreeMap), consensus cache, `used` flag, mate list, SV metadata |
| `Mate` (struct) | Mate-pair info: start/end/len for both read and mate, quality metrics |
| `StructuralVariantCounts` (struct) | SV evidence: splits, pairs, clusters |
| `InsOrDelLen` (enum: None/InsSeq/DelLen) | Explicit insertion/deletion discriminator extracted from Java's implicit string encoding |
| `var_utils::get_variants_from_map()` | Get-or-insert in nested `HashMap<i64, InnerMap<VarDesc, Variant>>`; uses unsafe pointer cast |
| `var_utils::get_variant_from_pos_map()` | Get-or-insert in position-level variant map |
| `var_utils::get_variation_from_seq()` | Get-or-insert in SoftClip seq map (position → base → Variant) |
| `var_utils::find_conseq()` | Compute consensus from soft-clipped reads; caches result in SoftClip |
| `var_utils::HomoPolymerChecker` | State machine for homopolymer run detection in find_conseq() |
| `var_utils::is_has_and_not_equals() / is_has_and_equals()` | Reference boundary check predicates |

## Java Correspondence

| Rust | Java | Notes |
|------|------|-------|
| `Variant` | `Variation.java` | 1:1 field mapping; names normalized to snake_case (varsCount→alt_depth, meanPosition→mean_pos, etc.) |
| `VarDesc` enum | Implicit String keys | **Design improvement**: Java uses plain Strings ("A", "+ACG", "-3#GG^TT&A"); Rust uses discriminated union with SmallVec storage for zero-alloc comparisons |
| `SoftClip` | `Sclip extends Variation` | **Composition over inheritance**: Rust uses `SoftClip { var: Variant, ... }` instead of OOP inheritance |
| `SoftClip.nt/seq` | `TreeMap<Integer, ...>` | BTreeMap preserves sorted iteration order (critical for find_conseq) |
| `SoftClip.soft` | `LinkedHashMap<Integer, Integer>` | `IndexMap<i64, usize>` preserves insertion order |
| `SoftClip.consensus_seq` | `String sequence` (nullable) | `Option<Vec<u8>>`; caches even empty sequences (Some(vec![])) to prevent re-computation |
| `Mate` | `Mate.java` | Direct field mapping; Rust uses full names (mate_start vs ms) |
| `StructuralVariantCounts` | `VariationMap.SV` inner class | Promoted to standalone struct; Java embeds in VariationMap |
| `InsOrDelLen` | Implicit in Java string parsing | Extracted into explicit enum for type safety |
| `var_utils::*` | `VariationUtils.java` static methods | 1:1 mapping of getVariation(), getVariationFromSeq(), findconseq(), helpers |

## Known Parity Traps

1. **Variant refReverseCoverage / refForwardCoverage naming bug**: Java Javadoc says meanings are swapped vs. field names. Rust must follow actual semantics used in ToVarsBuilder.toVars() output, not the misleading Java field names.

2. **SoftClip.nt / seq BTreeMap ordering**: find_conseq() iterates sorted order and exits early at position 3 if evidence < 6. Off-by-one in position comparison causes divergent consensus. Requires extensive fixture testing.

3. **SoftClip.soft insertion order**: IndexMap must match Java LinkedHashMap insertion order. VariationRealigner (line 1231) and StructuralVariantsProcessor iterate entries — order matters for SV cluster detection.

4. **VarDesc::Del parsing separators**: The `#`, `^`, `&` separators in deletion keys determine how CigarParser reconstructs variants. One misplaced separator causes **Silent Key Divergence** — variant maps to wrong VarDesc hash bin. Use to_key_string() extensively in tests.

5. **Variant::pstd/qstd ratchet semantics**: Once set to `true`, these flags **never revert to false**. CigarParser sets on first detection; no later stage may unset them.

6. **find_conseq() poly-A/T export toggle**: The `used` flag is an **export toggle**, not data deletion. The consensus sequence remains stored even when used=true. If Rust clears the sequence instead of just setting the flag, later reads will diverge.

7. **SoftClip.consensus_seq empty caching**: `Some(vec![])` means "already computed, result is empty". find_conseq() returns immediately on `is_some()`. Accidental setting of empty consensus can mask uncomputed state.

8. **get_variants_from_map() unsafe pointer**: Uses `unsafe { &mut *variant_ptr }` to avoid borrow-checker conflict. Safe only because callers don't mutate pos_map while holding the returned reference. All callers are pub(crate).

9. **InsOrDelLen ownership model**: Cloned into VarDesc::Del; mutations happen on VarDesc, not InsOrDelLen. Any code path expecting to mutate InsOrDelLen in-place will hit ownership issues.

10. **Mate field name mapping**: Ensure output columns in OutputVariant/ToVarsBuilder use correct field mappings (mate_start vs Java's mateStart_ms, etc.).

## Divergences from Java

| Area | Java | Rust | Rationale |
|------|------|------|-----------|
| Variant key representation | Plain String HashMap keys | Discriminated `VarDesc` enum with SmallVec | Type safety, zero-alloc comparisons, smaller memory footprint |
| SoftClip class hierarchy | `Sclip extends Variation` | `SoftClip { var: Variant, ... }` composition | Idiomatic Rust; avoids inheritance |
| Base count maps | `TreeMap<Character, Integer>` | `NucBaseMap<usize>` (fixed-size ACGT lookup) | Performance: no tree overhead for 4-element map |
| Map accessor pattern | Java GC handles references | Unsafe pointer cast in get_variants_from_map() | Required to satisfy borrow checker for nested map mutation |
| InsOrDelLen | Implicit in variant string | Explicit enum (None/InsSeq/DelLen) | Type safety; extracted from string parsing logic |

## Cross-Module Dependencies

**Called by (imports types/functions):**
- `mods/cigar_parser.rs` — imports all types; calls get_variants_from_map(), get_variation_from_seq() during CIGAR parsing
- `mods/variant_realigner.rs` — uses SoftClip consensus caching; calls find_conseq() at multiple sites, get_variants_from_map() ~30 times; reads/mutates Variant.pstd, qstd
- `mods/structural_variants_processor.rs` — iterates SoftClip.mates, SoftClip.soft; mutates StructuralVariantCounts
- `mods/to_vars_builder.rs` — reads all Variant accumulator fields; calls VarDesc methods for allele/type extraction
- `scopedata/` — stores Maps of Variants and SoftClips keyed by position and VarDesc

**Depends on:**
- `utils/vec_map.rs` (`VecMap` used as inner map type alias)
- Standard library types, `indexmap`, `smallvec`
