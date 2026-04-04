# variant_realigner

**Source**: `src/mods/variant_realigner.rs`
**LOC**: ~5,271
**Java counterpart**: `VariationRealigner.java` → [Java cache](../java/VariationRealigner.md)
**Status**: complete
**Last verified**: 2026-04-04

## Overview

The `variant_realigner` module is the second stage of the VarDict variant calling pipeline, responsible for **local realignment and refinement** of the raw variant maps produced by CigarParser before structural variant processing. It performs five major tasks: (1) SV filtering and mate clustering (`filter_all_sv_structures`), (2) MNP merging (`adjust_mnp`), (3) short indel realignment by mismatch absorption (`process_deletions`, `process_insertions`), (4) large deletion discovery from unpaired soft clips (`realign_large_deletions`), and (5) large insertion/duplication discovery (`realign_long_insertions_30`, `realign_long_insertions`). All operations mutate shared variant maps in-place—adding counts to target variants, removing consumed variants, and marking soft clips as used.

## Public API

| Function/Method | Purpose |
|----------------|---------|
| `new(reference_seq, reference_seed, ref_start)` | Basic constructor with minimal context |
| `new_with_context(reference_seq, reference_seed, ref_start, chromosome, bam_paths)` | Full constructor with chromosome and BAM paths for optional coverage loading |
| `with_reference_fallback(self, original_reference_seq, original_ref_start)` | Builder: set original reference fallback window for SV realignment |
| `process_deletions(data, position_to_deletions_count)` | Main deletion realignment: absorbs mismatches and soft clips into short deletions |
| `process_insertions(data, position_to_insertion_count)` | Main insertion realignment: absorbs mismatches and soft clips into short insertions |
| `filter_all_sv_structures(data)` | Filters all 8 SV lists + fusions; clusters mates, applies disc/cnt thresholds, sorts SOFTP2SV by varsCount |
| `realign_large_deletions(data)` | Discovers large deletions (>=10bp) from unpaired 5'/3' soft clips via `find_bp` breakpoint search |
| `realign_long_insertions_30(data)` | Pairs 3'/5' soft clips for >=30bp insertions; calls `process_deletions`/`process_insertions` recursively |
| `realign_long_insertions(data)` | Discovers large insertions/DUPs from single soft clips via `find_match` (both 5' and 3' directions) |
| `load_partial_ref_coverage(data, start, end)` | On-demand BAM re-parsing to load additional coverage for a region |
| `adjust_mnp(data, mnp)` | Merges partial MNP sub-variants back into parent MNP; absorbs matching soft clips |
| `is_low_complex_seq(seq) -> bool` | Static: returns true if >75% single base or <3 distinct bases |
| `is_match(seq1, seq2, dir) -> bool` | Static: sequence match check with <=3 mismatches and <15% error rate; direction-aware |
| `is_match_with_threshold(seq1, seq2, dir, mm_threshold) -> bool` | Static: sequence match with configurable mismatch threshold |
| `find_35_match(seq5, seq3) -> Match35` | Static: finds best overlapping match between 5' and 3' sequences |

## Java Correspondence

| Rust Function | Java Method | Notes |
|---|---|---|
| `new`, `new_with_context` | Constructor (L58-L72) | Rust uses `Arc` for shared ownership |
| `with_reference_fallback` | N/A | Rust-specific builder for fallback window |
| `process_deletions` | `realigndel()` (L593-L858) | Two-pass loop structure preserved |
| `process_insertions` | `realignins()` (L864-L1170) | Parallel structure to deletions |
| `filter_all_sv_structures` | `filterAllSVStructures()` (L353-L377) | Filters 8 SV lists + fusions; populates SOFTP2SV |
| `realign_large_deletions` | `realignlgdel()` (L1177-L1594) | Large deletion discovery from soft clips |
| `realign_long_insertions_30` | `realignlgins30()` (L1599-L1799) | Pairs 3'/5' clips for >=30bp |
| `realign_long_insertions` | `realignlgins()` (L1805-L2150) | Large insertion/DUP discovery |
| `load_partial_ref_coverage` | N/A | Rust-specific on-demand BAM re-parsing |
| `adjust_mnp` | `adjustMNP()` (L487-L569) | MNP merging and soft clip absorption |
| `is_low_complex_seq` | `islowcomplexseq()` (L2353-L2385) | Low-complexity filter |
| `is_match` | `ismatch()` (L2313-L2345) | Sequence comparison |
| `find_35_match` | `find35match()` (L2214-L2260) | Overlap matching |

## Known Parity Traps

1. **lgdel stale svcov carry**: Java `realignlgdel` carries `svcov` across loop iterations instead of resetting. Rust must keep `svcov` scoped outside iteration loops to match.

2. **Inversion reference fallback** (row3356): Synthetic inversion alleles from `StructuralVariantsProcessor::find_inv_sub` trigger `process_deletions` on a different reference window. Solved via `with_reference_fallback()` builder providing `original_reference_seq`/`original_ref_start`.

3. **row6306 stale soft-clip consensus**: Fresh NA12878 tail rerun showed content mismatches from stale soft-clip consensus sequences and improper variant state propagation.

4. **Deletion-key micro-optimization regression** (variant_realigner_deletion_lookup_attempts): Fast-path parsing + `get_key_value()` regressed `batch4_t1`. Do not retry without stronger evidence.

5. **Pass 2 loop index** (realigndel): Loops `(1..del_keys.len()).rev()` — element 0 is intentionally skipped. Java does `for (int i = tmp.size()-1; i > 0; i--)`.

6. **Left vs Right MNP merging asymmetry** (adjustMNP): Left check uses `tref.varsCount <= 0` (skip if zero or negative); right check uses `tref.varsCount < 0` (skip only if negative, allow zero). Parity-critical asymmetry.

7. **Reference window fallback**: `get_ref_base()` must fall back to `original_reference_seq`/`original_ref_start` for anchor positions outside the current working window.

8. **Soft-clip consensus caching**: `set_consensus_seq()` should only be called once per soft-clip. Multiple realignment passes must not overwrite existing consensus.

9. **SOFTP2SV sorting stability**: Sort by `varsCount` descending only (no tiebreaker). Non-deterministic order if equal count.

10. **bams shadowing** (realigndel pass 1): `bamsParameter` is shadowed; only null-ness is checked. When null, `noPassingReads()` is skipped.

## Divergences from Java

### Structural
- **Reference window fallback**: Rust adds `original_reference_seq`/`original_ref_start` fields + builder. Java had no fallback.
- **On-demand coverage loading**: Rust `load_partial_ref_coverage()` re-parses BAM for sparse regions. Java parses all data once.
- **Arc-wrapped reference sequences**: Rust uses `Arc<Vec<u8>>` for shared immutable ownership.
- **Public method granularity**: Rust exposes individual realignment methods as public API. Java's `realignIndels()` is a single orchestrator.

### Performance
- **Static utility methods**: `is_match`/`is_match_with_threshold` as free functions, reducing allocation overhead.
- **SmallVecBytes**: Stack-allocated small vector for normalized insertion keys.

### Algorithm
None identified. All algorithms ported faithfully.

## Cross-Module Dependencies

### What This Module Calls
| Module | Purpose |
|--------|---------|
| `data::patterns` | Regex patterns for variant description parsing |
| `data::reference` | Reference sequence and seed map for breakpoint search |
| `data::region` | Region definition for on-demand coverage loading |
| `mods::structural_variants_processor` | `RealignedVariationData` type |
| `mods::vardict_pipeline` | Used in `load_partial_ref_coverage` |
| `prelude` | `InnerMap`, `InnerMapEntry`, `LibDefaultHasher`, `SmallVecBytes` |
| `variants::var_utils` | `find_conseq` — consensus sequence extraction |
| `variants::variants` | `VarDesc`, `Variant`, `SoftClip`, `InsOrDelLen` |

### What Calls This Module
| Module | Usage |
|--------|-------|
| `vardict_pipeline` | Main pipeline orchestrator — calls all public realignment methods in sequence |
| `structural_variants_processor` | Calls `process_deletions()` on synthetic inversion alleles |

### Data Flow
**Input**: `RealignedVariationData` (mutable) — contains non_insertion_variants, insertion_variants, soft_clips_5end/3end, sv_* lists, ref_coverage, sv_counts
**Output**: Modified input (in-place mutations) — variant counts aggregated, consumed variants removed, soft clips marked `used`, SV counts updated

## Key Internal Functions

| Function | Purpose |
|---|---|
| `realign_deletion_mismatches()` | Core deletion realignment: two-pass mismatch absorption and MNP post-processing |
| `realign_with_softclips_5end()` | Absorb 5'-end soft clips matching deletion flanking sequence |
| `realign_with_softclips_3end()` | Absorb 3'-end soft clips matching deletion flanking sequence |
| `find_match()` | Core large insertion/DUP discovery: finds match in reference via seed map or full scan |
| `find_mm5()` | Walks 5' direction, collects mismatches within flanking sequence |
| `find_mm3()` | Walks 3' direction, collects mismatches within flanking sequence |
| `find_bi()` | Finds insertion breakpoint: leftmost match position and adjusted insertion sequence |
| `find_bp()` | Finds breakpoint position: searches for sequence anchor in reference |
| `is_match_ref()` | Checks if sequence matches reference at position with MM tolerance |
| `get_ref_base()` | Gets reference base; checks current window first, then original fallback |
| `merge_duplicate_deletion_keys()` | Merges variant keys representing same deletion with different notation |
| `merge_suffix_shift_insertion_keys()` | Merges insertion keys differing by single-base suffix shift |
| `filter_sv()` | Filters a single SV list: clusters mates, applies disc/cnt filter |
| `check_cluster()` | Clusters mates by proximity; returns dominant cluster with 60% threshold |
| `collect_softp2sv_first_used()` | Populates `softp2sv_first_used` map |
| `ensure_sv_marker()` | Ensures SV placeholder variant exists at position |
| `move_sv_marker()` | Moves SV marker and counts between positions |
| `rm_cnt()` | Subtracts source variant counts from destination (saturating) |
| `extract_inv_flanks()` | Parses inversion flanking sequences from SV description string |
