# cigar_modifier

**Source**: `src/mods/cigar_modifier.rs`
**LOC**: ~1,984
**Java counterpart**: `CigarModifier.java` → [Java cache](../java/CigarModifier.md)
**Status**: complete
**Last verified**: 2026-04-04

## Overview

The `CigarModifier` module performs iterative CIGAR string normalization between RecordPreprocessor and CigarParser. It transforms raw CIGAR strings through eight phases: strip boundary deletions, convert boundary insertions to soft clips, remove chimeric soft clips, collapse indel+soft-clip edges, and post-loop 5'/3' soft-clip realignment. Unlike Java's regex-based sequential rewrites, Rust uses `VecDeque<Cigar>` structures and pattern-matching helpers. Called only when `performLocalRealignment == true`.

## Public API

| Function/Method | Purpose |
|----------------|---------|
| `new(pos, cigar_str, query_seq, query_qual, ref_data, indel, max_read_len, region, rev_complementor) -> Self` | Constructor: stores per-read alignment data and reference context |
| `modify_cigar(&mut self) -> Result<ModifiedCigar, Error>` | Master transformation: applies all phases, returns adjusted position + CIGAR + trimmed sequences |

Returns `ModifiedCigar`: `align_start_pos`, `cigar` (`VecDeque<Cigar>`), `query_seq`, `query_qual`.

## Java Correspondence

| Rust | Java | Notes |
|------|------|-------|
| `new()` | `CigarModifier()` constructor | Initialize mutable per-read fields |
| `modify_cigar()` | `modifyCigar()` (L57-L248) | Master loop: phases 1-6 |
| `normalize_front_softclip_indel()` | Phase 4a (BEGIN_NUMBER_S_NUMBER_IorD) | Collapse leading S+I/D |
| `normalize_back_indel_softclip()` | Phase 4b (NUMBER_IorD_NUMBER_S_END) | Collapse trailing I/D+S |
| `normalize_front_softclip_match_indel()` | Phase 4c | Collapse leading S+short-M+I/D |
| `normalize_back_indel_match_softclip()` | Phase 4d | Collapse trailing I/D+short-M+S |
| `find_primary_realign_pattern()` | Phase 4g complex patterns | M-D-M-I-M-D, 3-del, 3-indel |
| `combine_to_close_to_correct()` | `combineToCloseToCorrect()` | Merge D-M-D/I (≤15bp) |
| `combine_to_close_to_one()` | `combineToCloseToOne()` | Merge I-M-D/I (≤15bp) |
| `two_dels_ins_to_complex()` | `twoDeletionsInsertionToComplex()` | M-D-M-I-M-D-M → M-D-I-M |
| `three_deletions()` | `threeDeletions()` | M-D-M-D-M-D-M → M-D/I-M |
| `three_indels()` | `threeIndels()` | Most complex: 8 outcome branches |
| `capture_mis_softly3_ms()` | `captureMisSoftlyMS()` | 3' boundary extend/retract |
| `combine_dig_s_dig_m()` | `combineDigSDigM()` | 5' boundary extend/retract |

## Known Parity Traps

### Trap 1: Edge-Normalization Phase Ordering (FIXED)
Java applies trailing I/D-M-S collapse before front-edge short-M-I/D-M. Rust must call back-edge functions before front-edge transforms. Example: `4M59I10M28S` → `101M`. Documented in `cigar_modifier_edge_order_parity_20260310.md`.

### Trap 2: I-M-D/I Pattern Prefix Check
`find_d_i_m_id_i()` returns immediately after first non-D-prefix I-M-[ID] pair — does NOT search later if first is D/H-prefixed. Matches Java's single-match-per-loop behavior.

### Trap 3: Pattern Priority
`find_primary_realign_pattern()` prefers complex > 3-del > 3-indel. Returns first match only.

## Divergences from Java

1. **String vs VecDeque**: Java uses regex patterns + `Pattern.compile()`; Rust uses `VecDeque<Cigar>` with pattern-matching helpers.
2. **Explicit borrowing**: Java mutates `this.position` in-place; Rust passes `&mut u32` through call stack.
3. **Single-pass pattern discovery**: Java evaluates sequential regex matchers; Rust iterates CIGAR vec once with accumulated offsets.
4. **HomoPolymerChecker utility**: Replaces Java's inline heteropolymer detection.
5. **Error propagation**: Java: try/catch wrapper; Rust: `?` operator with `Result`.
6. **Config naming**: `chimeric` (Java) → `chimeric_filter` (Rust).

## Key Internal Functions

| Function | Purpose |
|----------|---------|
| `find_primary_realign_pattern()` | Single-pass scan for complex/3-del/3-indel |
| `find_d_m_di_i()` | Locate D-M-(D\|I) pattern |
| `find_d_i_m_id_i()` | Locate non-D-prefix I-M-(D\|I) |
| `find_d_d()`, `find_i_i()` | Adjacent D-D or I-I for merging |
| `normalize_front_*`, `normalize_back_*` | Collapse leading/trailing indel+softclip |
| `capture_mis_softly3_ms()` | 3' boundary: extend M into S or retract |
| `capture_mis_softly3_mismatches()` | Trailing ≥3 mismatches → soft-clip |
| `combine_dig_s_dig_m()` | 5' boundary: extend M leftward into S |
| `combine_begin_dig_m()` | Leading ≤3 mismatches → soft-clip |
| `two_dels_ins_to_complex()` | M-D-M-I-M-D-M → M-D-I-M |
| `three_deletions()` | M-D-M-D-M-D-M → M-D/I-M |
| `three_indels()` | Most complex: 8 outcome branches |
| `combine_to_close_to_correct/one()` | Merge within 15bp match window |

## Cross-Module Dependencies

### Called By
- `cigar_parser.rs` — instantiated per BAM record when `performLocalRealignment == true`

### Calls Into
- `data/reference.rs` — Reference struct
- `data/region.rs` — Region struct
- `conf.rs` — Configuration (SEED_2, LOW_QUAL, chimeric_filter)
- `scopedata/global_read_only_scope.rs` — global config
- `variants/var_utils.rs` — HomoPolymerChecker, base comparisons
- `rust_htslib::bam` — Cigar types
- `crackle_kit` — RevComplementor for chimeric detection
