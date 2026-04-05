# structural_variants_processor

**Source**: `src/mods/structural_variants_processor.rs`
**LOC**: ~4,702
**Java counterpart**: `StructuralVariantsProcessor.java` → [Java cache](../java/StructuralVariantsProcessor.md)
**Status**: complete
**Last verified**: 2026-04-05

## Overview

Detects structural variants (deletions, inversions, duplications) from assembled soft-clip consensus sequences and discordant read pair evidence. Positioned after `VariantRealigner` and before `ToVarsBuilder` in the pipeline. Orchestrates six complementary SV finding routines (`find_del`, `find_inv`, `find_svs_del_candidates`/findsv, `find_del_disc`, `find_inv_disc`, `find_dup_disc`) in strict order, each marking processed clusters as used to prevent double-counting. Also provides **always-on** soft-clip rescue logic (`adj_snv`) that adds short soft-clipped reads (<=5 bp) to existing SNV signals, running even when SV detection is disabled. A critical parity challenge is reference window management: Rust maintains a single active contiguous window plus historical snapshots instead of Java's monotonically growing mutable reference map.

## Public API

| Function | Signature | Purpose |
|----------|-----------|---------|
| `new()` | `pub fn new(reference_seq: Vec<u8>, reference_seed: ReferenceSeedMap, ref_start: i64) -> Self` | Create processor with initial region reference |
| `new_with_context()` | `pub fn new_with_context(reference_seq: Arc<Vec<u8>>, reference_seed: Arc<ReferenceSeedMap>, ref_start: i64, chromosome: Option<String>, bam_paths: Vec<String>, shared_reference: Option<SharedReferenceHandle>) -> Self` | Create processor with on-demand reference extension context |
| `process()` | `pub fn process(&mut self, mut data: RealignedVariationData) -> ProcessedVariationData` | Main entry point; runs SV detection + soft-clip rescue |
| `process_with_region()` | `pub fn process_with_region(&mut self, mut data: RealignedVariationData, region: &Region) -> ProcessedVariationData` | Variant of `process()` with region context for optional RSS logging |
| `into_reference()` | `pub fn into_reference(self) -> Reference` | Consume processor; return the original region reference for downstream use |
| `historical_reference_windows()` | `pub fn historical_reference_windows(&self) -> Vec<(i64, i64)>` | Export list of remote windows loaded during SV processing |
| `historical_del_rightseq_variants()` | `pub fn historical_del_rightseq_variants(&self) -> HashMap<i64, HashSet<String>>` | Export DEL variants whose trailing rightseq was loaded from remote windows |

## Java Correspondence

| Rust Function | Java Method | Notes |
|---|---|---|
| `process()` | `process()` | Direct port; orchestrates `find_all_svs()` + `adj_snv()` |
| `find_all_svs()` | `findAllSVs()` | Orchestrates six SV routines in exact order |
| `find_del()` | `findDEL()` | DELs from soft-clip consensus + pair evidence |
| `find_inv()` | `findINV()` | Dispatcher to four `find_inv_sub()` calls |
| `find_inv_sub()` | `findINVsub()` | INVs; triggers secondary realignment via VariantRealigner |
| `find_svs_del_candidates()` | `findsv()` | Hot-path raw soft-clip scanner; **critical findsv-synthetic-probe gate** |
| `find_del_disc()` | `findDELdisc()` | DELs from discordant pairs only |
| `find_inv_disc()` | `findINVdisc()` | INVs from discordant pair combinations |
| `find_dup_disc()` | `findDUPdisc()` | DUPs from discordant pairs |
| `adj_snv()` | `adjSNV()` | Always-on; delegates to `adj_snv_5end()` and `adj_snv_3end()` |
| `find_match()` | `findMatch()` (forward-strand) | Forward seed-based alignment |
| `find_match_rev()` | `findMatchRev()` (reverse-strand) | Reverse-complement seed-based alignment |
| `is_match_ref()` | `ismatchref()` | Validates sequence vs reference with mismatch tolerance |
| `ensure_reference_span()` | `ReferenceResource.getReference()` | **Parity-sensitive**; loads remote reference, replaces active window |
| `join_ref()` | `joinRef()` | Joins reference bases across ranges with fallback |
| `mark_sv()` | `markSV()` | Overlap-based cluster marking |
| `mark_dup_sv()` | `markDUPSV()` | DUP-specific overlap marking |
| `is_overlap()` | `isOverlap()` | Fuzzy interval overlap detection |
| `sort_java_hashmap_order()` | (implicit HashMap iteration) | Sorts positions by Java 8 HashMap bucket order |

## Known Parity Traps

### Trap 1: Forward match in findsv MUST use historical windows (reversed 2026-04-05)
- Java's `findsv()` forward match sees all prior reference loads via its cumulative seed map (Java trap 23)
- Rust forward matches now call `find_match_internal(..., include_historical_windows=true, ...)`
- Previous restriction (`include_historical_windows=false`) caused missed DELs that Java finds
- Spurious INV candidates (e.g. chr13/054) are suppressed by the suspect-window retry mechanism instead
- See also: `findsv_suspect_end3_forward_del_pair_gated_rescue_20260405.md`

### Trap 2: Synthetic pre-INV forward probe gate
- Java accumulates SV buckets into `SOFTP2SV{softp}` map sorted by count, checks `[0].used`
- Rust reconstructs lazily via `is_softp2sv_first_used()`, scanning 8 SV vectors for max-count entry
- Fusion (inter-chromosomal) buckets must be included, or spurious INVs leak through

### Trap 3: 3' inversion historical windows are suspect
- 3' INV clustering is asymmetric to 5'; Rust can create low-support 3' clusters Java doesn't
- `suspect_hist_windows` set tracks these; `extend_seed_positions_from_historical_windows()` skips them

### Trap 4: Reference window restoration must preserve historical coordinates
- Java keeps all prior loads in one growing mutable hash; Rust replaces active window
- Prior windows saved as coordinate snapshots for later shared-reference lookup

### Trap 5: findsv suspect-window retry replaces redundant probe (updated 2026-04-05)
- `should_skip_inv_after_historical_forward_probe()` is removed
- Replaced by `find_match_findsv_allow_suspect_windows()`: a suspect-window retry gated by `check_pairs()`
- When forward match returns no hit, retry allows suspect 3' historical windows; result accepted only if discordant pairs confirm
- This is Rust-only logic not present in Java; validated by 250/250 full parity but should be monitored for novel datasets

### Trap 6: Double-strand-flip INV boundaries
- Left-alignment after INV walks backward checking `ref[softp] == complement(ref[bp])`
- Asymmetry in loop order requires direction/side-dependent control flow

### Associated Repo Memory Files
- `softp2sv_fusion_used_state_parity_20260310.md`
- `sv_reference_window_memory_fix_20260309.md`
- `debug_inv_extra_truncation_20260324.md`
- `row6300_structural_seed_history_fix_20260311.md`
- `row6300_hs37d5_inv_structural_blocker_20260311.md`
- `na12878_18813399_forward_find_match_preempts_inv_20260313.md`
- `na12878_19811894_rightseq_history_window_disproves_cigar_owner_20260313.md`
- `na12878_18807210_del_rightseq_fallback_parity_20260313.md`

## Divergences from Java

### Structural
1. **Reference window lifecycle**: Java = single unbounded growing REF hash; Rust = single active window + historical coordinate snapshots. Saves memory; parity maintained via historical re-fetch from shared reference.
2. **Seed lookups from remote history**: Java has all seeds immediately available; Rust re-loads from shared reference on each seed query via `extend_seed_positions_from_historical_windows()`.
3. **SOFTP2SV map reconstruction**: Java explicitly builds; Rust reconstructs on-demand by scanning 8 SV cluster vectors.
4. **Reference coverage**: Rust uses `RefCoverage` (dense array or sparse map depending on config).
5. **Shared-reference fallback is conditional**: Java never does chromosome-wide linear scan; Rust disables fallback by default (`extend_seed_positions_from_shared_reference()` is intentional no-op).

### Performance
1. Historical window seed scanning guarded by `findsv_active` flag.
2. Shared-reference fallback is no-op to avoid 28x slowdown from chromosome-wide scan.

### Algorithm
None identified. All SV algorithms ported faithfully.

## Cross-Module Dependencies

### Inbound (callers)
| Module | Usage |
|--------|-------|
| `vardict_pipeline` | Main orchestrator; calls `process()` or `process_with_region()` |

### Outbound (callees)
| Module | Purpose |
|--------|---------|
| `VariantRealigner` | `find_inv_sub()` instantiates realigner for secondary INV realignment |
| `Reference`, `ReferenceSeedMap` | Reference sequence and seed storage |
| `SharedReferenceHandle` | On-demand chromosome sequence fetch |
| `var_utils::find_conseq()` | Consensus sequence extraction from soft clips |
| `Configuration` | SEED_1, SEED_2, SVFLANK, MINSVCDIST, DISCPAIRQUAL constants |
| `GlobalReadOnlyScope` | Global config singleton |

### Key Data Types
| Type | Purpose |
|------|---------|
| `RealignedVariationData` | Main input/output structure (defined at top of file) |
| `ProcessedVariationData` | Output wrapper after SV processing |
| `Variant`, `VarDesc`, `SoftClip` | Variant/SV records |
| `StructuralVariantCounts` | Aggregated SV evidence counters |
| `InversionClusterKind`, `InversionSide` | Enums for INV direction/side dispatch |
| `ReferenceProbeState` | Snapshot for reference state save/restore |

### Key Internal Functions
| Function | Purpose |
|---|---|
| `process_internal()` | Core logic: `find_all_svs()` if SV enabled, restore reference, `adj_snv()` always |
| `find_all_svs()` | Orchestrate 6 SV routines in exact order |
| `find_del()` | DELs from soft-clip consensus + pair evidence |
| `find_inv_sub()` | INVs for one direction/side; triggers secondary realignment |
| `find_svs_del_candidates()` | findsv equivalent; raw soft-clip scanner |
| `find_del_disc()`, `find_inv_disc()`, `find_dup_disc()` | Discordant-pair-only SV discovery |
| `adj_snv_5end()`, `adj_snv_3end()` | Short soft-clip SNV rescue |
| `find_match_internal()`, `find_match_rev_internal()` | Core seed-based alignment with configurable window scoping |
| `find_match_findsv_allow_suspect_windows()` | Suspect-window retry for forward match; gated by `check_pairs()` |
| `get_ref_base_from_historical_windows()` | Window-bounds-gated reference base lookup from historical snapshots |
| `is_softp2sv_first_used()` | Reconstruct SOFTP2SV used-state by scanning all SV vectors |
| `check_pairs()`, `peek_pairs()` | Discordant pair overlap scanning |
| `mark_sv()`, `mark_dup_sv()` | Overlap-based cluster marking |
| `ensure_reference_span()` | Load remote reference; replaces active window, saves historical |
| `restore_original_reference_window()` | Restore to pre-SV reference state |
| `seed_positions_with_scope()` | Multi-level seed lookup (current/original/historical/shared) |
| `extend_seed_positions_from_historical_windows()` | Re-load historical windows for seed scanning |
| `get_ref_base()` | Multi-level fallback reference base access |
| `join_ref()` | Join reference bases across window boundaries |
| `sort_java_hashmap_order()` | Replicate Java 8 HashMap iteration order |