# Configuration

**Source**: `src/conf.rs`
**Java counterpart**: `Configuration.java`, `CmdParser.java` → [Java cache](../java/Configuration.md)
**Status**: complete

## Overview

The Configuration module is a plain data struct holding all ~35 command-line parameters and algorithmic constants required by the VarDict pipeline. Unlike Java (split across `Configuration.java`, `CmdParser.java`, and `VarDictLauncher.java`), Rust consolidates the struct in `conf.rs` while delegating CLI parsing to `clap` in `src/bin/vardict.rs`. The struct is read-only once constructed and shared across all pipeline modules via `GlobalReadOnlyScope`.

## Public API

| Type | Name | Purpose |
|------|------|---------|
| **Struct** | `Configuration` | All 35+ pub fields: thresholds, filters, SV params, amplicon/CRISPR, MSI, runtime state |
| **Impl** | `Default` | Defaults matching Java `CmdParser` behavior |
| **Constants** | `SEED_1`, `SEED_2`, `ADSEED`, `SVMAXLEN`, `SVFLANK`, `MINSVCDIST`, `MINMAPBASE`, `MINSVPOS`, `DISCPAIRQUAL`, `MAX_EXCEPTION_COUNT`, `LOW_QUAL` | Hardcoded algorithm seeds and thresholds |

### Major Field Categories
- **Calling Thresholds**: `freq` (0.01), `lofreq` (0.05), `goodq` (22.5), `qratio` (1.5), `mapq` (0.0)
- **Depth & Count Limits**: `minr` (2), `min_bias_reads` (2), `min_match` (0), `mismatch` (8)
- **SV Parameters**: `inssize` (300), `insstd` (100), `insstdamt` (4), `sv_min_len` (1000), `disable_sv` (false)
- **Read/Alignment Filters**: `sam_filter` (0x504), `remove_duplicated_reads`, `chimeric_filter`, `perform_local_realignment` (true)
- **Amplicon/CRISPR**: `amplicon_based_calling` (None), `crispr_cutting_site` (0), `crispr_filtering_bp` (0)
- **Output Control**: `debug` (false), `fisher` (false)
- **Runtime**: `exception_counter` (Arc<AtomicUsize>)

## Java Correspondence

All 35 Rust fields have 1:1 Java counterparts with matching defaults.

| Rust | Java | Notes |
|------|------|-------|
| `sam_filter: u32` | `Configuration.samfilter: String` | Parsed at CLI time in Rust; decoded at usage in Java |
| `read_pos_filter: f64` | `Configuration.readPosFilter: int` | Functionally equivalent for comparisons |
| `exception_counter: Arc<AtomicUsize>` | `AtomicInteger` | Both thread-safe counters |
| `downsampling: Option<f64>` | `Double` (boxed) | `None` = Java `null` (not set) |
| `mapping_quality: Option<u8>` | `Integer` (boxed) | `None` = Java `null` (not set) |
| CLI parsing: `clap::Parser` | Apache Commons CLI + `CmdParser` | Structural difference, no output impact |

## Known Parity Traps

1. **`perform_local_realignment` default logic** (HIGH): Java defaults to `true` via `CmdParser` logic (`1 == getIntValue(cmd, "k", 1)`). Rust struct default is `false`. CLI parsing in `src/bin/vardict.rs` must apply the same default-to-1 logic. Realignment ON by default in Java.

2. **`-p` (Pileup) side-effects** (HIGH): Java overwrites `freq = -1.0` and `minr = 0` AFTER parsing other flags when `-p` is present. Rust must replicate this ordering.

3. **`downsampling` tristate** (MEDIUM): `None` vs `Some(0.0)` distinction matters — Java boxed `Double` null = not set.

4. **`mapping_quality` boxed Option** (MEDIUM): Java `Integer` boxed → null if not set. Rust `Option<u8>` must not accidentally use `0` instead of `None`.

5. **`sam_filter` String → u32 parsing** (LOW): Java stores as String, parses with `Integer.decode()` at usage. Rust parses once at CLI time. Must handle `0x`, `0X`, and plain decimal formats.

6. **`read_pos_filter` type mismatch** (LOW): Java `int 5`, Rust `f64 5.0`. Functionally equivalent for comparisons.

## Divergences from Java

| Aspect | Java | Rust | Impact |
|--------|------|------|--------|
| CLI parsing | Apache Commons CLI + CmdParser | `clap::Parser` derive | None (internal) |
| `sam_filter` type | `String` | `u32` | None (bitwise ops identical) |
| `read_pos_filter` type | `int` | `f64` | None (comparison semantics same) |
| Exception counter | `AtomicInteger` | `Arc<AtomicUsize>` | None (both thread-safe) |
| `GlobalReadOnlyScope` layout | Holds sample names, printer type | Different layout | Configuration struct unaffected |

## Cross-Module Dependencies

- **Calls**: None (plain data struct, no method calls)
- **Called by** (consumers of Configuration fields):
  - `cigar_parser` — `goodq`, `vext`, `downsampling`, `include_n_in_total_depth`, `crispr_filtering_bp`
  - `variant_realigner` — `perform_local_realignment`, `move_indels_to_3`, `goodq`, `crispr_cutting_site`
  - `structural_variants_processor` — `inssize`, `insstd`, `insstdamt`, `disable_sv`, `sv_min_len` + SV constants
  - `to_vars_builder` — `freq`, `minr`, `min_bias_reads`, `read_pos_filter`, `qratio`, `mapq`, `mismatch`, MSI frequencies
  - `output_variant` — `debug`, `fisher`, `sv_min_len`
  - `bam_reader` — `mapping_quality`, `remove_duplicated_reads`, `mismatch`, `sam_filter`, `trim_bases_after`
  - `bin/vardict.rs` — All fields (sole writer via clap)
