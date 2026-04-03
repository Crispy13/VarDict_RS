---
name: mem-optimization
description: "Memory profiling and optimization for VarDict-rs parity. Use when Rust RSS exceeds Java baseline, diagnosing allocation hotspots with DHAT, designing data structure replacements, or measuring before/after memory impact. Covers: baseline measurement, DHAT heap profiling, hotspot analysis, optimization design, benchmark validation, and parity verification."
argument-hint: "Describe the memory issue, e.g. 'Rust uses 7GB vs Java 3GB on 5MB region' or 'profile chr1:100M-105M allocations'"
---

# Memory Profiling & Optimization — VarDict-rs

End-to-end workflow for diagnosing and fixing memory regressions in the Rust VarDict port while preserving byte-identical output parity with Java.

## When to Use

- Rust RSS exceeds Java RSS on a comparable region
- Need to identify allocation hotspots in the pipeline
- Evaluating a data structure replacement (e.g., HashMap → VecMap)
- Verifying a memory optimization didn't break parity

## Prerequisites

- `rust_build_env` conda environment active
- Test BAM + reference FASTA (e.g., `testdata/NA12878.mapped.*.bam`, `testdata/hs37d5.fa`)
- VarDictJava built and runnable for baseline comparison
- DHAT feature available: `Cargo.toml` has `dhat-heap = ["dep:dhat"]`

## Procedure

### Phase 1: Baseline Measurement

Establish the Java and Rust RSS for a specific region before any changes.

1. **Pick a representative region** — use a region that triggers the memory issue.
   Good defaults: `chr1:100000000-101000000` (1MB), `chr1:100000000-105000000` (5MB).

2. **Measure Java RSS**:
   ```bash
   /usr/bin/time -v java -Xmx8g -classpath VarDictJava/build/libs/VarDict-1.8.3.jar \
     com.astrazeneca.vardict.Main \
     -G testdata/hs37d5.fa -b testdata/NA12878.mapped.*.bam \
     -N NA12878 -th 4 -z -c 1 -S 2 -E 3 -g 4 \
     -R "1:100000000-105000000" 2>&1 | grep "Maximum resident"
   ```

3. **Measure Rust RSS** (debug-release build):
   ```bash
   cargo build --profile debug-release
   /usr/bin/time -v ./target/debug-release/vardict \
     -G testdata/hs37d5.fa -b testdata/NA12878.mapped.*.bam \
     -N NA12878 -th 4 -z -c 1 -S 2 -E 3 -g 4 \
     -R "1:100000000-105000000" > /dev/null 2>&1
   ```

4. **Record baseline** in a table:
   ```
   | Region | Java RSS | Rust RSS | Ratio |
   ```

5. **Decision gate**: If Rust RSS ≤ 1.2x Java → STOP, no optimization needed.

### Phase 2: DHAT Heap Profiling

Identify where allocations happen and how much memory they hold at peak.

6. **Build with DHAT** (must disable mimalloc — DHAT needs the system allocator):
   ```bash
   cargo build --profile debug-release --no-default-features --features dhat-heap
   ```

7. **Run and collect DHAT data**:
   ```bash
   ./target/debug-release/vardict \
     -G testdata/hs37d5.fa -b testdata/NA12878.mapped.*.bam \
     -N NA12878 -th 1 -z -c 1 -S 2 -E 3 -g 4 \
     -R "1:100000000-105000000" > /dev/null
   ```
   This produces `dhat-heap.json` in the working directory.
   **IMPORTANT**: Use `-th 1` (single thread) for clean DHAT traces.

8. **Analyze with DHAT viewer**: Open https://nnethercote.github.io/dh_view/dh_view.html
   and load the JSON. Sort by "Total bytes" or "At t-gmax bytes" (bytes at global peak).

9. **Record top-5 allocation sites** at global peak:
   ```
   | Rank | Allocation Site | Bytes at Peak | Block Count | Notes |
   ```

10. **Identify optimization candidates** — look for:
    - Many small allocations (>100k blocks) = per-position overhead → data structure change
    - Large single allocations = bulk buffers → capacity tuning or deferred construction
    - Temporaries that overlap with peak = scheduling opportunity

### Phase 3: Optimization Design

Design the fix. Always prefer the simplest change that preserves Java execution order.

11. **Evaluate options** using this priority:
    1. **Type alias swap** — change inner container type (e.g., HashMap → VecMap). Lowest risk, maximum propagation.
    2. **Capacity tuning** — pre-allocate or shrink. Medium risk.
    3. **Deferred construction** — delay allocation past peak. Changes execution order — **requires user approval**.
    4. **Structural redesign** — new data layout. Highest risk.

12. **Parity constraint check**: Will this change affect output ordering or floating-point paths?
    - If YES → design must preserve Java-equivalent iteration order (use `IndexMap` or ordered alternative).
    - If NO → proceed.

13. **Design the data structure** (if creating a new one):
    - Must implement same API as the type it replaces (get, get_mut, insert, remove, entry, iter, keys, values, len, is_empty, clear, retain)
    - Must implement `Default`, `FromIterator`, `IntoIterator`
    - Entry API must have `Occupied`/`Vacant` with `or_default`, `or_insert`, `or_insert_with`, `and_modify`
    - Write unit tests for each method

14. **Implement via type alias interception**: Minimize call-site changes by changing the type alias, not every usage.
    ```rust
    // In src/prelude.rs — change ONE line to swap backend:
    pub(crate) type InnerMap<K, V> = VecMap<K, V>;       // ← current
    // pub(crate) type InnerMap<K, V> = HashMap<K, V, LibDefaultHasher>; // ← revert
    ```

### Phase 4: Benchmark Validation

Validate the optimization is actually faster/smaller, not just theoretically better.

15. **Create or run microbenchmark** comparing old vs new data structure:
    ```bash
    cargo bench --bench vecmap_vs_hashmap_bench
    ```
    Benchmark must cover: insert, get_hit, get_miss, entry_or_default, iterate, remove, memory size.
    Test at entry counts: 1, 2, 4, 8, 16, 32.

16. **Analyze crossover point**: At what entry count does the new structure lose to the old?
    - If crossover < typical entry count → STOP, optimization is harmful.
    - If crossover > typical entry count → proceed.

17. **Record benchmark results**:
    ```
    | Operation | N=1 | N=2 | N=4 | N=8 | N=16 | N=32 | Crossover |
    ```

### Phase 5: Integration Measurement

Measure actual RSS improvement on the same region as Phase 1.

18. **Rebuild with the optimization** (debug-release, with mimalloc):
    ```bash
    cargo build --profile debug-release
    ```

19. **Measure Rust RSS** with same region and flags as step 3.

20. **Run all cargo tests**:
    ```bash
    cargo test -- --include-ignored
    ```

21. **Run parity check** on the measurement region:
    ```bash
    diff <(./target/debug-release/vardict [flags] -R "1:100000000-105000000") \
         tmp/java_reference_output.tsv
    ```

22. **Record final results**:
    ```
    | Region | Java RSS | Rust Before | Rust After | Ratio | Tests | Parity |
    ```

23. **Decision gate**:
    - Rust After ≤ target AND parity PASS AND tests PASS → **SHIP IT**
    - Rust After > target → loop back to Phase 2 for next hotspot
    - Parity FAIL → revert and diagnose

### Phase 6: Documentation

24. **Update work report** at `copilot-office/missions/mem-opt*/copilot-work-report.md`
25. **Record divergences** — if the optimization deviates from Java logic, document in:
    - Source code doc comment on the new data structure
    - `copilot-office/codebase/CODEBASE.md` under "Divergences from Java"
    - Repo memory: `/memories/repo/<descriptive_name>.md`
26. **Update project plan** status to CLOSED with final measurements.

## Key Constraints

- **Parity is non-negotiable**: Any optimization that changes output is rejected.
- **Follow Java execution order** unless explicitly approved by the user.
- **Use `./tmp`** for all temporary/intermediate files, never `/tmp`.
- **Single-thread for profiling**: Always `-th 1` for DHAT runs to get clean traces.
- **mimalloc vs DHAT**: They are mutually exclusive. DHAT requires `--no-default-features --features dhat-heap`. Production uses `--features mimalloc-global` (default).

## Lessons Learned (from mem-opt2)

| Lesson | Detail |
|--------|--------|
| Inner map overhead dominates | With ~5M positions, even 388 bytes/map overhead = ~1.9 GB waste |
| VecMap wins at 1-4 entries | 1.3-1.7x faster insert, 1.4-1.6x faster get, 2-5x faster iterate, ~2x less memory |
| HashMap wins at 8+ entries for lookups | get_hit crossover at ~6 entries, get_miss at ~3 |
| Type alias swap is safest | Changing 2 lines in `vardict_pipeline.rs` propagated to entire pipeline |
| Deferred construction is risky | May not reduce peak if other allocations fill the gap |
| SmallVec inline storage backfired | For `Variant` (256 bytes value), SmallVec inline made memory worse than HashMap |
| Always benchmark before claiming victory | Theoretical savings don't always materialize (deferred seed_map increased RSS) |
