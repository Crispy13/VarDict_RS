# ParallelPipeline

**Source**: `src/mods/parallel_pipeline.rs`
**LOC**: 1,057
**Java counterpart**: `ExecutorService` + `CompletableFuture` threading model (no single Java class)
**Status**: complete

## Overview

`ParallelPipeline` orchestrates multi-threaded variant calling across sharded genomic regions with **strict input-order preservation** in output. It wraps `rayon::ThreadPool`, replacing Java's `ExecutorService` + `CompletableFuture` pattern. Design: (1) partition regions into chunks, (2) spawn indexed tasks computing `RegionResult` independently, (3) reduce results in manifest order via sorted index reassembly. Four execution paths: batch collect, streaming channels, direct-write to `Write`, and a disabled region-prefetch optimization.

## Public API

| Function/Method | Purpose |
|----------------|---------|
| `ParallelPipeline::new()` / `::with_reference()` / `::with_chromosomes()` | Constructors with reference-loading variants |
| `process_regions_vardict()` | Parallel VarDict pipeline; returns `Vec<RegionResult>` |
| `process_regions_vardict_streaming()` | Parallel with channel sender for real-time consumption |
| `process_regions_vardict_direct_write()` | Single-thread sequential; streams to `Write` sink |
| `process_regions_vardict_to_output()` | Convenience; returns `Vec<String>` |
| `get_header()` | TSV column headers via OutputVariant |
| `RegionResult` (struct) | Output lines + region + optional error from one region |

## Java Correspondence

| Rust | Java | Notes |
|------|------|-------|
| `rayon::ThreadPool` | `ExecutorService` | Auto-scales via work-stealing |
| `(usize, RegionResult)` indexed + sort | `CompletableFuture<RegionResult>` | BED order preserved via enumeration index |
| `crossbeam_channel::bounded()` | `BlockingQueue<Integer>` | Streaming path backpressure |
| `map_init()` per worker scope | Per-thread BAM reuse | Reduces open/close churn |
| `process_regions_vardict_direct_write()` | OutputVariant in reader thread | Matches Java simple-mode I/O |

## Known Parity Traps

- **Output Ordering**: Rayon executes non-deterministically. Must sort by enumeration index before flattening. Forgetting sort = scrambled output. (Repo: `parity_shard_runner_ordered_reduction_20260313.md`)
- **Streaming Channel Ordering**: Consumer must re-sort by index. The streaming path pushes results as they complete; reassembly in BED order is the caller's responsibility. (Repo: `single_thread_streaming_postprocessor_20260324.md`)
- **Thread Determinism**: Regional independence + output sort eliminates reliance on scheduling. Per-region variant calls are deterministic given the same input.
- **BAM Handle per Worker**: `map_init()` ensures one BAM handle per worker lifetime, not per region. No output change but misunderstanding this leads to unnecessary file I/O.
- **Prefetch Path**: Disabled by default (`SIMPLE_MODE_PREFETCH_REGION_GROUPS_ENABLED = false`). Risk of off-by-one in grouping/partitioning if enabled. Defer until parity stable.

## Divergences from Java

| Aspect | Java | Rust | Impact |
|--------|------|------|--------|
| Thread pool | Auto-created ExecutorService | Owned `Arc<ThreadPool>` with explicit `num_threads` | Caller controls parallelism |
| Error handling | CompletableFuture exception chain | Independent per-region; errors in `RegionResult.error` | No cascading failure |
| Memory logging | None | `/proc/self/status` parsing via `VARDICT_RSS_MEMORY_LOG` env var | Rust-only profiling; zero-cost if disabled |
| Region prefetch | No equivalent | Optional group-prefetch + partition (disabled) | Rust-only optimization; not yet proven safe for parity |

## Key Internal Functions

| Function | Purpose |
|----------|---------|
| `process_regions_vardict_one_region_per_task()` | Main worker: enumerate regions, parallel map with BAM via `map_init()`, sort by index |
| `process_regions_vardict_streaming_impl()` | Parallel map sending indexed results to channel |
| `build_prefetch_region_groups()` | Group regions by chr + max span (2.5 Mbp threshold) |
| `collect_prefetched_records()` | BAM fetch for region span, clone records with header |
| `partition_prefetched_records_by_region()` | Distribute records into per-region buckets (multi-overlap) |
| `record_overlaps_region()` | Alignment interval vs. region overlap test |
| `current_process_memory_snapshot()` | Parse `/proc/self/status` for VmRSS/VmSwap (Linux only) |

## Cross-Module Dependencies

**Calls**:
- `data::shared_reference` — `SharedReferenceHandle` (immutable Arc'd reference)
- `data::region` — `Region` struct (genomic locus)
- `data::bam_reader` — `BamReader` (opened per worker thread)
- `mods::vardict_pipeline` — `VarDictPipeline` (full CIGAR→Realign→ToVars pipeline)
- `mods::output_variant` — TSV header generation
- `scopedata::global_read_only_scope` — Global config singleton
- `rayon` — ThreadPool, work-stealing parallelism
- `crossbeam_channel` — Bounded channel for streaming

**Called by**:
- `bin/vardict.rs` — CLI entry point invokes parallel pipeline for all modes
