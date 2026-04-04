# cigar_parser

**Source**: `src/mods/cigar_parser.rs`
**LOC**: ~5,004
**Java counterpart**: `CigarParser.java` → [Java cache](../java/CigarParser.md)
**Status**: complete
**Last verified**: 2026-04-04

## Overview

The `cigar_parser` module is the **core variant detection engine** of VarDict-rs. It processes BAM records to extract all variants (SNVs, MNVs, indels, and complex variants) by interpreting CIGAR strings operation-by-operation. Every variant that VarDict ultimately reports originates from this module. It reads BAM records, applies optional local realignment via `CigarModifier`, then systematically walks each CIGAR element—expanding MNVs on matches, handling indel discovery/combination, processing soft clips for alignment-error evidence, and building structural variant discordant pair clusters. Results accumulate into five core output dictionaries: non-insertion variants, insertion variants, reference coverage, 5'/3' soft-clip observations, and SV structures.

The module is the pipeline bottleneck for variant discovery—every bar of evidence must flow through here. It sits between `CigarModifier` (upstream normalizer) and `VariantRealigner` / `ToVarsBuilder` (downstream consumers).

## Public API

| Function/Method | Purpose |
|----------------|---------|
| `CigarParser::new(region, reference, instance)` | Constructor: allocates buffers, initializes variant maps, seed maps, SV accumulators |
| `process_records(records)` | Entry point: streaming BAM record iteration with error handling |
| `process_record(record)` | Process a single BAM record; dispatches to `parse_cigar()` |
| `get_non_insertion_vars()` / `take_non_insertion_vars()` | Access/extract SNVs, MNVs, deletions |
| `take_non_insertion_vars_insert_index()` | Insertion order metadata (Java LinkedHashMap parity) |
| `get_insertion_vars()` / `take_insertion_vars()` | Access/extract insertion variants |
| `take_mnp()` | Multi-nucleotide polymorphism count map |
| `take_position_to_insertion_count()` | Unique insertion descriptions per position |
| `take_position_to_deletions_count()` | Unique deletion descriptions per position |
| `get_ref_coverage()` / `take_ref_coverage()` | Reference coverage (depth per genomic position) |
| `get_soft_clips_5end()` / `take_soft_clips_5end()` | 5' soft-clip evidence structures |
| `get_soft_clips_3end()` / `take_soft_clips_3end()` | 3' soft-clip evidence structures |
| `get_max_read_len()` | Maximum read length observed (used for MNV/SV window sizing) |
| `get_discordant_count()` | Reads with mate on different chromosome (SV signal) |
| `take_svfdel()` / `take_svrdel()` | Forward/reverse deletion clusters (discordant pairs) |
| `take_svfdup()` / `take_svrdup()` | Forward/reverse duplication clusters |
| `take_svfinv5()` / `take_svrinv5()` / `take_svfinv3()` / `take_svrinv3()` | 5'/3' inversion clusters |
| `take_svffus()` / `take_svrfus()` | Fusion clusters (cross-chromosome, keyed by mate tid) |
| `take_splice_count()` / `take_splice_count_insert_index()` | Splice junction counts and insertion order metadata |

## Java Correspondence

| Rust Function | Java Method (CigarParser.java) | Notes |
|---------------|--------------------------------|-------|
| `process_records()` / `process_record()` | `process(Scope<RecordPreprocessor>)` | Rust is streaming iterator-based vs Java's stateful RecordPreprocessor |
| `parse_cigar()` | `parseCigar(String, SAMRecord)` | Main CIGAR loop: per-record variant extraction |
| `process_soft_clip()` | `processSoftClip(...)` | 5'/3' soft-clip handling with chimeric detection and match-back |
| `process_insertion()` | `processInsertion(...)` | Insertion variant creation with position adjustment and ref count correction |
| `process_deletion()` | `processDeletion(...)` | Deletion variant creation with adjacent indel detection |
| `process_not_matched()` | `processNotMatched(u32)` | N operator (intron/splice junction) handling |
| `add_variation_for_matching_part()` | `addVariationForMatchingPart(...)` | SNV/MNV variant creation from match block |
| `add_variation_for_deletion()` | `addVariationForDeletion(...)` | Deletion variant creation with position/quality metrics |
| `prepare_sv_deletion_structures_for_analysis()` | `prepareSVStructuresForAnalysis(...)` | Discordant pair cluster building for DEL/DUP/INV/FUS |
| `is_trim_at_opt_t_bases()` | `isTrimAtOptTBases(boolean, int)` | Trim read bases at configured distance |
| `is_closer_then_vext_and_good_base()` | `isCloserThenVextAndGoodBase(...)` | Adjacent indel within vext bases (complex variant hint) |
| `find_offset()` | `findOffset(...)` | Scan next match segment for mismatches to append to indel description |
| `skip_overlapping_reads()` | `skipOverlappingReads(...)` | Skip double-counted bases in overlapping paired reads |
| `skip_sites_out_region_of_interest()` | `skipSitesOutRegionOfInterest(CigarOperator[])` | CRISPR mode: skip reads outside cutting site window |
| `clean_up_cigar()` | `cleanupCigar(SAMRecord)` | Convert leading/trailing I→S, remove H |
| `sclip5_high_quality_processing()` | `sclip5HighQualityProcessing(...)` | Populate 5' soft-clip structures |
| `sclip3_high_quality_processing()` | `sclip3HighQualityProcessing(...)` | Populate 3' soft-clip structures |
| `ref_has_and_equals()` / `ref_has_and_not_equals()` | `isHasAndEquals()` / `isHasAndNotEquals()` | Safe reference base comparison (null-guarded) |
| `add_cnt()` / `sub_cnt()` | `addCnt(...)` / `subCnt(...)` | Increment/decrement variant counters |
| `parse_cigar_with_amp_case()` | `parseCigarWithAmpCase(...)` | Amplicon mode: filter overlapping segments |
| `is_followed_by_match_and_indel()` | `isInsertionOrDeletionWithNextMatched()` | Multi-indel pattern recognition |
| `adj_ins_pos()` | `VariationRealigner.adjInsPos()` | Left-align insertion by rotation |

## Known Parity Traps

### Active Traps

1. **Read position parallel counters** — `read_pos_excluding_softclip` vs `read_pos_including_softclip` tracked in parallel. Errors cause false `tp` (distance from read end) and downstream position drift. Java: two parallel counters. Rust: same design. Soft-clip processing resets `start` to original after match-back but creates variants at adjusted positions.

2. **MNV detection while-loop state mutation** — The MNV detection loop in `parse_cigar()` (L2049–2225) mutates `i`, both `read_pos_*` counters, and `start` in-place to absorb consecutive mismatches. If a mismatch is skipped or boundaries miscomputed, all downstream positions diverge. Java: L408–464 CigarParser.java.

3. **Offset statefulness across CIGAR elements** — `offset` persists across the CIGAR loop. Set by deletion/insertion processing, consumed in the next match block's inner loop (`for i = offset...`). After the match loop, if `moffset != 0`, it's applied to reference and read positions. Mishandling resets cause position drift.

4. **Insertion/deletion pattern recognition** — Pattern `D + short_M + I/D` detected by `is_followed_by_match_and_indel()` (L4997+). Requires `ci + 2 < cigar.len()`, next segment is Match with length ≤ `vext`, and segment after is Ins/Del but not the one after that. Return value of `ci` is incremented by 2 to skip consumed CIGAR elements.

5. **Soft-clip match-back position tracking** — 5' soft clip: match-back loop decrements `cigarElementLength` and `start`. After loop, `cigarElementLength` IS restored but `start` remains at first non-matched position. Forgetting to restore `cigarElementLength` shifts all subsequent read position calculations.

6. **Quality accumulation rounding** — Quality accumulated as `f64`, individual scores are `u8` (Phred). Division `q / (qbases + qibases)` must match Java's Phred double accumulation exactly.

7. **Insertion position adjustment (adjInsPos)** — After insertion discovery, if pure ATGC and base before matches reference, insertion left-aligned via `adj_ins_pos()` (L4945+). Rotation logic must match Java's substring rotation exactly.

8. **Reference count subtraction for insertions** — After emitting insertion variant, if base at insertion position matches reference, a non-insertion variant count is SUBTRACTED via `sub_cnt()` to prevent AF > 1.0. The reference variant must already exist or the subtraction is a no-op.

9. **Discordant read direction encoding** — SV structures accumulate discordant clusters by position. New clusters created when distance > `MINSVCDIST * max_read_len`. Direction encoding is inverted: forward = -1, reverse = 1 in Java.

10. **Leading/trailing I→S conversion in `clean_up_cigar()`** — If edge-I conversion is missed, soft clips treated as insertions, causing position arithmetic errors.

11. **Reference window bounds** — `Reference::get(pos)` returns `Option<u8>`. During match-ahead scanning, falling outside the loaded window returns `None`. Java silently returns null; Rust must avoid unwrap.

12. **Splice count key computation** — Keyed by `(start - 1, start + cigar_len - 1)` tuples (intron boundaries). Off-by-one in tuple key computation scatters splice counts incorrectly.

13. **Integer overflow in variant counters** — Variant counters are `u32`; Java `int` silently wraps. Rust uses `saturating_sub()` which may diverge at very high read counts (>2^32).

### Resolved Traps (from Repo Memory)

- **cigar_modifier_edge_order_parity** — Leading/trailing CIGAR element handling (RESOLVED)
- **chr5_softclip_tail_bounds_parity** — Soft-clip quality scanning bounds (RESOLVED)
- **lgdel_stale_svcov_parity_trap** — Large deletion SV coverage tracking (RESOLVED)
- **softp2sv_fusion_used_state_parity** — Soft-clip position in SV tracking (RESOLVED)

## Divergences from Java

1. **Streaming vs. stateful model** — Java: single `process()` receives `Scope<RecordPreprocessor>` with stateful record iterator. Rust: streaming iterator interface (`process_records(I)` where `I: Iterator<Item = &'a mut Record>`). **Parity impact**: NONE.

2. **Owned vs. borrowed data structures** — Java: `LinkedHashMap<Integer, Map<String, Variation>>` with GC cleanup. Rust: `HashMap<i64, InnerMap<VarDesc, Variant>>` with `InnerMap` = `IndexMap`. Variants extracted via `take_*()` methods. **Parity impact**: NONE — ordering preserved via IndexMap.

3. **SmallVec inline storage** — Short sequences (< 32 bytes) stored inline in `SmallVec<[u8; 32]>` to reduce allocations. **Parity impact**: NONE.

4. **Configuration access** — Java: instance methods on `Scope` or static `ScopeData`. Rust: thread-local singleton `GlobalReadOnlyScope::instance()` wrapped in Arc. **Parity impact**: NONE.

5. **Error handling** — Java: exceptions bubble up with try-catch. Rust: `Result<(), Error>` with optional `#[cfg(feature = "catch-panics")]` wrapping. **Parity impact**: NONE.

6. **Reverse complement caching** — Rust uses `crackle_kit::RevComplementor` (threadlocal cache). **Parity impact**: NONE.

## Cross-Module Dependencies

### Upstream (what calls this module)

| Module | Usage |
|--------|-------|
| `vardict_pipeline` | Instantiates `CigarParser::new()`, calls `process_records()`, extracts outputs via `take_*()` |

### Downstream (what this module calls)

| Module | Function/Type | Purpose |
|--------|---------------|---------|
| `CigarModifier` | `modify_cigar()` | Local realignment before CIGAR parsing |
| `Reference` | `get(pos)`, `has_and_equals()`, `has_and_not_equals()` | Reference base access |
| `RefCoverage` | `inc(pos, count)` | Increment coverage |
| `var_utils` | `get_variant_from_pos_map()`, `get_variants_from_map()`, `get_variation_from_seq()` | Variant lookup/creation |
| `GlobalReadOnlyScope` | `conf`, `chr_lens`, `seed` | Configuration access |
| `crackle_kit` | `RevComplementor`, `complement_base()`, `NucBaseMap` | Sequence operations |
| `rust_htslib::bam` | `Record`, `Cigar`, `CigarString` | BAM record parsing |

### Key Data Types

| Type | From Module | Purpose |
|------|------------|---------|
| `Variant` | `variants::variant` | Core variant record |
| `VarDesc` | `variants::variants` | Variant description (SNV/Ins/Del/Raw) |
| `SoftClip` | `variants::variants` | Soft-clip evidence structure |
| `RefCoverage` | `data::data` | Dense coverage array |
| `Reference` | `data::reference` | Reference bases + seed map |
| `Region` | `data::region` | Genomic window definition |
| `InnerMap<K, V>` | `prelude` | Order-preserving HashMap (IndexMap) |

### Data Flow

**Inputs**: Region, Reference, BAM records iterator, GlobalReadOnlyScope

**Outputs** (via `take_*()` methods): non_insertion_vars, insertion_vars, ref_coverage, soft_clips_5end/3end, sv_* clusters, splice_count, metadata (insert indices, max read len, discordant count)

### Key Internal Functions

| Function | Lines | Purpose |
|----------|-------|---------|
| `parse_cigar()` | L1238–2431 | Main per-record CIGAR parsing loop (hot path) |
| `process_soft_clip()` | L2433–2715 | Soft-clip match-back/forward, chimeric detection |
| `process_deletion()` | L2717–3006 | Deletion variant creation, multi-indel pattern detection |
| `process_insertion()` | L2788–3237 | Insertion variant creation, position adjustment |
| `process_not_matched()` | L3509–3531 | Intron (N operator) handling |
| `add_variation_for_matching_part()` | L3077–3366 | SNV/MNV variant creation from match segments |
| `add_variation_for_deletion()` | L2591–2697 | Deletion variant finalization |
| `prepare_sv_deletion_structures_for_analysis()` | L687–1139 | Discordant pair clustering |
| `clean_up_cigar()` | L3337–3384 | Leading I→S, trailing I→S, remove H |
| `sclip5_high_quality_processing()` | L3462–3566 | 5' soft-clip evidence accumulation |
| `sclip3_high_quality_processing()` | L3568–3672 | 3' soft-clip evidence accumulation |
| `add_cnt()` / `sub_cnt()` | free fns | Increment/decrement variant counters |
| `is_followed_by_match_and_indel()` | L4997+ | Multi-indel pattern recognition |
| `adj_ins_pos()` | L4945+ | Left-align insertion by rotation |
