# Committee Plan: VarDict-rs Performance Optimization

**Date:** 2026-04-02
**Committee:** Opus (Claude Opus 4.6), Sonnet (Claude Sonnet 4.6), Gemini (Gemini 3.1 Pro), GPT (GPT 5.4)
**Rounds:** 1 discussion round (all contested points resolved in round 1)

---

## Plan: Profile-Driven Optimization With Parity Safety

VarDict-rs is already 1.7× faster than Java and uses 0.50–0.61× Java RSS. This plan targets another **15–30% runtime reduction** through allocation elimination in hot loops, data structure improvements, and better parallelism — while preserving byte-identical output parity. The approach is measurement-first: establish profiling baselines before touching any hot paths, then execute high-confidence wins in priority order.

**Steps**

### Phase A: Measurement & Infrastructure (prerequisite for all other phases)

1. **Establish two-track performance scoreboard**
   - Track A (product runtime): wall-clock, instructions, cache-misses, peak RSS at `-th 1/4/8` on 1MB and 5MB regions
   - Track B (developer throughput): `cargo build --profile debug-release --timings`, `cargo test` time, parity cell time
   - Save dated baseline to `./tmp/perf_baseline_YYYYMMDD.md`
   - Commands:
     ```bash
     conda activate rust_build_env && export LIBCLANG_PATH=$CONDA_PREFIX/lib
     cargo build --profile debug-release
     hyperfine --warmup 1 --runs 5 \
       './target/debug-release/vardict -G testdata/hs37d5.fa -b testdata/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam -N NA12878 -f 0.01 --threads 1 -R "1:100000000-101000000" > /dev/null'
     perf stat -e cycles,instructions,cache-references,cache-misses,branches,branch-misses \
       ./target/debug-release/vardict ... -R "1:100000000-101000000" > /dev/null
     ```

2. **Create flamegraph and DHAT profiling scripts** (*parallel with step 1*)
   - `scripts/flamegraph.sh`: `perf record -g --call-graph dwarf` → `inferno-collapse-perf` → `inferno-flamegraph`
   - `scripts/dhat_profile.sh`: `cargo build --features dhat-heap` → run → DHAT viewer
   - `scripts/cachegrind.sh`: `valgrind --tool=cachegrind` for cache-miss hotspot detection
   - Deliverable: top 5 CPU hotspots identified via flamegraph

3. **Expand benchmark coverage** (*parallel with step 1*)
   - Add `benches/region_bench.rs`: 1KB, 50KB, and 100KB regions from NA12878 low-coverage BAM
   - Add structural variants benchmark if SV processor shows up in flamegraph
   - Register in `Cargo.toml` with `harness = false`

4. **Add per-stage timing instrumentation** in `vardict_pipeline.rs`
   - Instrument: `run_cigar_parser_from_bam`, `run_variant_realigner`, `run_structural_variants_processor`, `run_to_vars_builder`
   - Gate behind `Level::INFO` — zero production overhead

### Phase B: P0 Quick Wins — Zero Parity Risk (implement immediately after Phase A baseline)

5. **Remove per-record qname String allocation** in `src/mods/cigar_parser.rs:process_record()`
   - Move `String::from_utf8_lossy(record.qname()).to_string()` into error/trace branches only
   - Impact: eliminates millions of heap allocations per chromosome. Est. 1–3% wall-clock.
   - Risk: Zero — error paths only.

6. **Gate `cigar.to_string()` behind debug check** in `src/mods/cigar_parser.rs`
   - `self.last_modified_cigar = Some(cigar.to_string())` → gate behind `self.needs_current_qname()`
   - Impact: eliminates millions of string allocations in non-debug runs.
   - Risk: Zero — debug-only field.

7. **Remove or gate `catch_unwind` from `process_record()`** in `src/mods/cigar_parser.rs`
   - This wraps the hottest function (4,800 LOC, called millions of times), preventing LLVM cross-boundary inlining and optimization
   - Approach: gate behind `#[cfg(feature = "catch-panics")]` feature flag (default off in release builds), converting panics to `Result::Err` propagation for production
   - Impact: Est. 2–5% CPU improvement.
   - Risk: Low — panics propagate to region-level handler. Add feature flag for optional safety.

8. **Enable incremental compilation for debug-release profile**
   - Add `incremental = true` to `[profile.debug-release]` in `Cargo.toml`
   - Impact: faster dev iteration. Revert if binary bloat >10%.

9. **Add mold linker to `.cargo/config.toml`** (optional, document availability) (*parallel with step 8*)
   ```toml
   [target.x86_64-unknown-linux-gnu]
   linker = "clang"
   rustflags = ["-C", "link-arg=-fuse-ld=mold"]
   ```

### Phase C: P1 Data Structure Improvements — Profile-Validated

10. **Replace `ref_coverage` HashMap with dense Vec** in `src/mods/cigar_parser.rs` (*depends on step 2 flamegraph*)
    - `ref_coverage: HashMap<i64, usize>` → `Vec<usize>` indexed by `(pos - region.start)`
    - Region span is known at construction. Handle out-of-bounds positions (reads extending past region) with bounds checking.
    - Impact: Est. 3–8% CPU. Major cache-line improvement for inner-loop coverage lookups.
    - Risk: Medium — audit all consumers of `ref_coverage` to ensure correct indexing. Requires parity test.

11. **Replace seed map `Vec<u8>` keys with fixed-size representation** in `src/data/reference.rs:build_seed_map()`
    - `HashMap<Vec<u8>, Vec<i64>>` → `HashMap<[u8; 20], Vec<i64>>` or packed `SeedKey { len: u8, bits: u64 }` with 2-bit nucleotide encoding
    - Impact: Est. 2–5% CPU. Eliminates O(ref_len × 2) heap allocations per seed map build.
    - Risk: Low — lookup-only structure, not output-affecting. Must verify hash equality.

12. **CigarParser state pooling across regions** (*depends on step 10*)
    - Add `CigarParser::reset_for_region(&mut self)` that `.clear()`s all maps without deallocating
    - Store parser in thread-local or pass by `&mut` through pipeline
    - Add `#[cfg(test)]` equivalence check: default-constructed vs reset parser produce identical output
    - Impact: Eliminates O(regions × maps) allocations, reuses HashMap capacity.
    - Risk: Low-medium — requires careful audit of all mutable fields. Fixture-backed verification mandatory.

13. **SV seed O(n²) uniqueness → sidecar set** in `src/mods/structural_variants_processor.rs`
    - Replace `positions.contains()` with sidecar `FxHashSet<i64>` for membership + `Vec<i64>` for ordered output in `seed_positions_with_scope`, `seed_positions_with_findsv_scope`, `extend_seed_positions_from_historical_windows`
    - Impact: Est. 3–8% on SV-heavy shards. Classic algorithmic improvement.
    - Risk: Medium — visitation order must be preserved. Parity-test on inversion/deletion shards.

### Phase D: P2 Parallelism & I/O

14. **BAM reader reuse across regions per thread** (*depends on step 2 profiling evidence*)
    - Apply single-threaded reader-reuse pattern from `process_regions_vardict_direct_write` to parallel path via thread-local `BamReader`
    - Impact: Est. 5–10% reduction in region setup overhead for small-region BED files.
    - Risk: Low — `fetch()` between regions is well-understood.

15. **Cost-based region scheduling for skewed BED inputs** (*depends on step 1 scaling data*)
    - Add simple cost heuristic (region span or prefetched record count) to `select_region_batches_for_execution` in `src/bin/vardict.rs`
    - Preserve output order via existing `(usize, RegionResult)` indexing
    - Impact: Est. 5–15% throughput on mixed/skewed inputs. Minor on uniform BED files.
    - Risk: Low-medium — output ordering preserved through existing index-based consumer.

### Phase E: Conditional / Evidence-Gated (only if profiling justifies)

16. **Module splitting** — only if `cargo --timings` shows `cigar_parser.rs` or `vardict_pipeline.rs` dominate rebuild time (>20% of incremental compile). If justified, split `vardict_pipeline.rs` first (cleaner boundaries).

17. **Postprocessor ordering metadata caching** — only if benchmarks show postprocessor is >10% of pipeline time. Cache Java bucket index at variant assembly time.

18. **ToVarsBuilder transient allocation rewrite** — only if DHAT shows `calculate_variant_statistics` is a material allocation hotspot. Must preserve floating-point evaluation order exactly.

**Relevant files**
- `src/mods/cigar_parser.rs` — steps 5, 6, 7, 10, 12 (process_record, parse_cigar, ref_coverage, parser state)
- `src/mods/vardict_pipeline.rs` — step 4, 12 (stage timing, parser threading)
- `src/data/reference.rs` — step 11 (build_seed_map, ReferenceSeedMap type)
- `src/mods/structural_variants_processor.rs` — step 13 (seed_positions_with_scope, ensure_reference_span)
- `src/mods/parallel_pipeline.rs` — step 14 (BamReader init, streaming pipeline)
- `src/bin/vardict.rs` — step 15 (select_region_batches_for_execution)
- `Cargo.toml` — steps 3, 8 (benchmarks, profile config)
- `.cargo/config.toml` — step 9 (mold linker)
- `benches/region_bench.rs` — step 3 (new benchmark)
- `scripts/flamegraph.sh`, `scripts/dhat_profile.sh`, `scripts/cachegrind.sh` — step 2 (new profiling scripts)

**Verification**

For every optimization step:
1. `cargo test` — all 258 tests pass (0 failures)
2. Focused fixture tests for the touched module (cigar_parser, variant_realigner, structural_variants, to_vars_builder)
3. Parity byte-diff on at least 3 chromosomes: `bash copilot-office/m9-option-parity/copilot-desk/scripts/na12878_parity_v2.sh --chr 1 --chr 20 --chr 22 --opts "" --opts-label default --rust-only --no-build --parallel 4 --no-stop`
4. `perf stat -e cycles,instructions,cache-references,cache-misses` before/after with `-th 1`
5. For scheduling changes (step 15): scaling sweep at `-th 1/2/4/8` with `hyperfine`
6. Regression guard: Rust wall-clock must remain ≤0.7× Java on benchmark workloads

**Decisions**
- Seed map optimization uses **local key compaction** (Position B), not global pre-build. Lower risk, preserves algorithm structure. Revisit global pre-build only if local optimization proves insufficient.
- Module splitting is **evidence-gated** (Position D). Run `cargo --timings` before committing to the 4+ hour refactor.
- Thread scheduling is **P2 priority** — after allocation wins but before conditional items. A compromise between GPT's P0 and Opus's P3 positions.
- SmolStr for CountMap is **deferred** — higher-impact zero-dependency wins exist. Revisit if profiling shows CountMap as a hotspot.
- Dedicated BAM I/O decoding thread is **rejected** — insufficient evidence that BAM decompression is a bottleneck vs compute.
- `catch_unwind` removal uses a **feature flag** approach rather than unconditional removal, preserving optional per-record error recovery.

**Further Considerations**
1. **`panic = "abort"` in deploy profile**: Would eliminate all unwind overhead globally without code changes. Tradeoff: process terminates on first panic instead of continuing to next region. For production bioinformatics pipelines this may be acceptable. Recommend investigating after feature-flag approach proves the performance delta.
2. **Pileup mode specifically**: The pileup path (`-p` flag) generates orders of magnitude more output than non-pileup. Phase C/D optimizations may have amplified impact in pileup mode. Consider adding a pileup-specific benchmark in step 3.
3. **Future: SIMD for sequence comparison**: If flamegraphs show a dominant inner comparison loop in realigner or seed matching, consider `memchr` or SIMD-accelerated byte scanning. Only after all safe algorithmic improvements are exhausted.

**Explicit No-Go List**
- Do NOT reintroduce chromosome-wide seed scans in `structural_variants_processor.rs` (caused 28× slowdown)
- Do NOT change `java_hashmap_iteration_order*` or `java_hashmap_bucket_index` functions (parity-critical output ordering)
- Do NOT change float formatting/rounding (`round_half_even`, `java_format_double`)
- Do NOT switch allocators (mimalloc is working well)
- Do NOT add stateful reverse-deletion lookup caches in `variant_realigner.rs` (prior regression)
- Do NOT use `unsafe` code for optimization until all safe alternatives are proven insufficient via profiling

---

## Provenance

How consensus was reached:

| Category | Count |
|----------|-------|
| Full consensus (4/4) | 4 points (profiling first, benchmarks, build time, validation gate) |
| Strong consensus (3/4) | 3 points (qname alloc, seed map keys, no-go list) |
| Resolved via discussion | 2 points (seed map global vs local, module splitting) |
| Resolved by supermajority | 0 points |
| Resolved by Chief decision | 1 point (thread scaling priority → P2 compromise) |
| Escalated to user | 0 points |

Unique contributions adopted: 7 (catch_unwind, ref_coverage Vec, gate cigar.to_string, CigarParser pooling, SV seed sidecar, mold linker, BAM reader reuse)
Unique contributions deferred: 2 (postprocessor cache, SmolStr)
Unique contributions rejected: 1 (dedicated I/O thread)

## Discussion Log

### Round 1
**Contested points resolved:**
- **Seed map strategy**: 3/3 present (Opus, Gemini, GPT) converged on Position B (local key optimization). Sonnet (absent from round) had proposed Position A (global pre-build). Supermajority resolution.
- **Module splitting**: 3/3 present converged on Position D (evidence-first with `cargo --timings`). Gemini added that if evidence warrants it, split both files aggressively.
- **Thread scaling**: 2/3 present (Gemini, GPT) favored Position A (cost-based scheduling as high priority). Opus favored Position B (low priority). Counting absent Sonnet's likely Position B, this was a 2-2 split. **Chief decided**: P2 compromise — implement after allocation wins but before conditional items. GPT's argument about skewed-input throughput was compelling, but Opus's caution about implementation complexity and marginal impact on uniform inputs was valid.

**Unique contributions discussion:**
- Strong agreement (3/3 ADOPT) on: gate `cigar.to_string()`, SV seed sidecar set, mold linker
- Majority agreement (2/3 ADOPT + 1 MODIFY) on: catch_unwind, ref_coverage Vec, CigarParser pooling, BAM reader reuse
- GPT's MODIFY votes added valuable caveats: feature flag for catch_unwind, verify scope for ref_coverage, reset API for CigarParser pooling
- GPT rejected BAM reader reuse (claiming it's already done), but review of the code shows the parallel path does NOT reuse readers (only single-threaded path does). 2/3 majority adopted it.
- GPT rejected SmolStr (type churn concerns) and I/O thread (no evidence). Both decisions accepted.
