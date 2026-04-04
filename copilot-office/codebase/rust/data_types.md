# Data Types

**Source**: `src/data.rs`, `src/data/` (`bam_reader.rs`, `ref_coverage.rs`, `reference.rs`, `shared_reference.rs`, `region.rs`, `patterns.rs`)
**LOC**: ~2,500 (estimated across sub-modules)
**Java counterpart**: `data/Region.java`, `data/Reference.java`, `data/Patterns.java`, `SAMFileParser.java`, `ReferenceResource.java`, InitialData fields
**Status**: complete

## Overview

The data infrastructure module provides five core services plus a regex pattern registry:

1. **BAM I/O** (`bam_reader.rs`) — thread-safe per-region BAM reading with SAM flag filtering via `rust_htslib`
2. **Regional coverage tracking** (`ref_coverage.rs`) — dense array backed by HashMap overflow for positions outside the window
3. **Reference genome + seed maps** (`reference.rs`) — single-chromosome sequence slice with precomputed 17-mer and 12-mer seed indices
4. **Shared multi-chromosome reference** (`shared_reference.rs`) — once-loaded genome shared across threads via `Arc<SharedReference>`
5. **Region specification** (`region.rs`) — genomic region container with display vs. computation coordinates
6. **Compiled regex patterns** (`patterns.rs`) — 26 `LazyLock<Regex>` constants used throughout variant detection

## Public API

| Function/Method | Purpose |
|----------------|---------|
| `BamReader::new()` | Open indexed BAM file |
| `BamReader::fetch_region()` / `BamRecordIter` | Iterator over records in a genomic region |
| `passes_filter()` | SAM flag + MAPQ filtering (subset of Java preprocessRecord) |
| `RefCoverage::new()` | Create coverage tracker for a genomic window |
| `RefCoverage::get()` / `set()` / `inc()` | Position-level coverage access with sentinel-aware zero preservation |
| `RefCoverage::iter_sorted()` | Sorted iteration for post-processing |
| `Reference::from_seq()` / `from_seq_with_start()` | Construct reference from sequence bytes |
| `Reference::build_seed_map()` / `clear_seed_map()` | Build 17/12-mer seed indices; release memory early |
| `Reference::get()` / `has_and_equals()` | Base lookup by 1-based genomic position |
| `FastaReader` / `fetch_seq()` / `get_reference()` | Indexed FASTA I/O and Reference construction |
| `Region::new()` / `new_extended()` / `new_with_insert()` | Construct normal, extended, or amplicon regions |
| `Region::chr()` / `start()` / `end()` / `display_start()` | Coordinate accessors (1-based interface) |
| `SharedReference::load_chromosome()` / `load_all()` | Load genome sequences into shared store |
| `SharedReference::get_base()` / `get_subseq()` | Thread-safe 1-based sequence lookups |
| `get_cap_group!()` macro | Safe regex capture group extraction |

## Java Correspondence

| Rust | Java | Notes |
|------|------|-------|
| `BamReader` + `BamRecordIter` | `SAMFileParser` + htsjdk `SamRecordIterator` | Uses rust_htslib; includes chromosome name normalization |
| `passes_filter()` | `RecordPreprocessor.preprocessRecord()` | Extracted subset — Java version also handles downsampling, dups, mismatch scoring |
| `RefCoverage` | 5 separate `LinkedHashMap<Integer, Integer>` in InitialData | Unified struct with dense array + HashMap overflow |
| `Reference` | `data.Reference` | `Arc<Vec<u8>>` wrapping enables O(1) clone vs Java deep copy |
| `ReferenceSeedMap` | Internal HashMap fields in Reference | `HashMap<[u8; 17], Vec<i64>>` + `HashMap<[u8; 12], Vec<i64>>` using fixed-size byte array keys |
| `Region` | `data.Region` | 1:1 field mapping; Rust uses direct data access, Java uses final fields |
| `SharedReference` | `ReferenceResource` (ThreadLocal) + `GlobalReadOnlyScope` | Explicit `Arc` sharing replaces Java's ThreadLocal singleton pattern |
| `ChromosomeData` | (implicit in ReferenceResource) | Wraps `Arc<Vec<u8>>` + length; enables sequence dedup across regions |
| `Patterns` (26 `LazyLock<Regex>`) | `data.Patterns` (static `Pattern` fields) | 1:1 pattern correspondence |

## Known Parity Traps

1. **RefCoverage Dense Zero Sentinel** — Java `HashMap<i64, usize>` distinguishes `put(pos, 0)` from absent key. Dense Rust array uses `u32::MAX` as sentinel for "never written" to preserve this distinction. Failure would corrupt strand bias and MAF calculations. *(repo: ref_coverage_dense_zero_sentinel_20260403)*

2. **Reference Arc + Allocator Interaction** — Earlier `Arc<Vec<u8>>` deallocation with mimalloc triggered SIGSEGV (C malloc/free mismatch). Resolved after rust-htslib v1.0.0c3 upgrade removed C-malloc workaround. Clone is now O(1) refcount bump. *(repo: reference_arc_internals_mem_opt_20260324, mimalloc_sigsegv_fix_memory_benchmarks_20260324)*

3. **Seed Map Early Release** — Seed maps (~30 MB for 3 MB regions) must be cleared via `clear_seed_map()`. Must check `Arc::strong_count() == 1` before `Arc::make_mut()` to avoid cloning the entire map when shared. If shared, allocates new empty Arc instead. *(repo: reference_seed_map_early_drop_owned_reference_20260324)*

4. **Reference Region Start Is 1-Based** — `Reference.region_start` is a 1-based genomic position. Base at position `p` requires index computation `p - region_start`. Off-by-one errors in coordinate conversion are a frequent parity failure source.

5. **FastaReader 1-Based to 0-Based Conversion** — Input `start`/`end` are 1-based inclusive; `bio::io::fasta::IndexedReader.fetch()` expects 0-based half-open `[start-1, end)`. Misalignment breaks all downstream variant positions.

6. **BamReader 1-Based to 0-Based Conversion** — Same pattern: 1-based inclusive input → 0-based half-open for `rust_htslib::bam::IndexedReader.fetch()`. Uses `start.saturating_sub(1)` for safety.

7. **Chromosome Name Normalization** — BAM files may use "chr1" while FASTA uses "1" or vice versa. Both `BamReader` and `SharedReference` implement independent normalization, storing under both original and normalized names.

## Divergences from Java

| Aspect | Java | Rust | Rationale |
|--------|------|------|-----------|
| Reference cloning | Deep copy of Vec and HashMap | O(1) `Arc` refcount bump | Avoid 30+ MB allocations per thread |
| Coverage storage | 5 separate `LinkedHashMap<i64, usize>` | Single `RefCoverage` with dense array + HashMap overflow | Memory-efficient hybrid; preserves set-zero vs absent semantics via sentinel |
| Chromosome lookup | `SamReader.getFileHeader()` direct name | Bidirectional map (both "chr20" and "20") | CLI may provide unnormalized chromosome names |
| Reference lifetime | `ReferenceResource` ThreadLocal singleton | `Arc<SharedReference>` passed explicitly | No global mutable state; explicit ownership |
| Coordinate system | Region, Reference, BAM all 0-based internal | Region API is 1-based; internal storage varies | Match VarDict CLI conventions; requires careful boundary handling |
| Pattern compilation | Static `Patterns` class initialized once | Per-pattern `LazyLock<Regex>` | Equivalent cost; clearer per-pattern dependency tracking |
| BAM record flow | Direct Record pass-through | Iterator abstraction + `cache_cigar()` per read | Rust ownership model requires explicit lifetime handling |

## Cross-Module Dependencies

**Calls (depended upon by this module)**:
- `rust_htslib` — BAM file I/O in `bam_reader.rs`
- `bio::io::fasta` — Indexed FASTA reading in `reference.rs`
- `regex` — Pattern compilation in `patterns.rs`
- `std::sync::Arc` — Shared ownership throughout

**Called by (modules that depend on this)**:
- `CigarParser` — uses `Reference` (seed map lookups, base access), `RefCoverage` (coverage tracking), `Patterns`
- `VariantRealigner` — uses `Reference` (seed map for realignment)
- `StructuralVariantsProcessor` — uses `RefCoverage`, `Reference`, `Patterns`
- `ToVarsBuilder` — uses `RefCoverage`, `Patterns`
- `VarDictPipeline` — uses `BamReader`, `Region`, `SharedReference`, `Reference`
- `ParallelPipeline` — uses `BamReader`, `SharedReference`, `Region`
- `OutputVariant` — uses `Region` (coordinate display)
- `bin/vardict.rs` — uses `SharedReference` (genome loading), `Region` (BED parsing)
