---
name: perf-optimization
description: "Performance optimization workflow for VarDict-rs. Use when: user requests optimization or runtime is slower than expected, throughput regression, CPU hotspot analysis, cache miss diagnosis, lock contention, I/O bottleneck, flamegraph profiling, benchmarking before/after, or optimizing any Rust code path. Covers: goal setting, baseline measurement, profiling (perf/flamegraph/cachegrind/DHAT), root cause analysis, optimization design, implementation, validation, and regression prevention."
argument-hint: "Describe the performance issue, e.g. 'chr1 5MB region takes 33s vs Java 25s' or 'profile pipeline CPU hotspots'"
---

# Performance Optimization — VarDict-rs

Systematic workflow for diagnosing and resolving performance issues in the Rust VarDict port.
Applies to CPU time, memory, I/O, and concurrency — while preserving byte-identical output parity with Java.

## When to Use

- Rust wall-clock time exceeds Java on comparable input
- Throughput regression after a code change
- Need to identify CPU hotspots, cache misses, or lock contention
- Evaluating an algorithmic or data structure change
- Creating a benchmark suite for a subsystem

## When NOT to Use

- Pure memory (RSS) issues with no runtime impact → use `mem-optimization` skill instead
- Output parity mismatches → use `parity-check` skill instead
- Correctness bugs → standard debugging

---

## Phase 1: Define the Goal

Every optimization starts with a measurable target. Without one, you cannot know when to stop.

### 1.1 — Identify the metric

| Metric | Tool | When to use |
|--------|------|-------------|
| Wall-clock time | `/usr/bin/time -v`, `hyperfine` | End-to-end runtime |
| CPU cycles / instructions | `perf stat` | Micro-optimization |
| Cache miss rate | `perf stat`, `cachegrind` | Data layout changes |
| Branch mispredictions | `perf stat` | Control flow changes |
| Heap allocations | DHAT | Allocation-heavy paths |
| Peak RSS | `/usr/bin/time -v` | Memory-bound workloads |
| Throughput (regions/sec) | Custom harness | Pipeline parallelism |
| Lock contention | `perf lock`, thread traces | Multi-threaded bottlenecks |

### 1.2 — Set the target

State it as: **"Reduce [metric] from [current] to [target] on [workload]."**

Examples:
- "Reduce wall-clock from 33s to ≤25s on chr1:100M-105M with -th 4"
- "Reduce cache miss rate from 12% to <5% in CigarParser hot loop"
- "Achieve ≥0.9x Java throughput on full-chromosome NA12878"

### 1.3 — Identify constraints

- **Parity**: Output must remain byte-identical to Java. Non-negotiable.
- **Correctness**: All `cargo test -- --include-ignored` must pass.
- **Maintenance**: Prefer simple changes. Exotic optimizations need justification.
- **Java execution order**: Preserved unless user explicitly approves deviation.

---

## Phase 2: Baseline Measurement

Establish reproducible measurements before ANY changes.

### 2.1 — Environment preparation

```bash
conda activate rust_build_env
unset CFLAGS CXXFLAGS CPPFLAGS LDFLAGS
export LIBCLANG_PATH=$CONDA_PREFIX/lib

# Build release (with mimalloc, the production allocator)
cargo build --profile debug-release
```

### 2.2 — Measure wall-clock time

Use `hyperfine` for statistical rigor (runs N iterations, reports mean/stddev):

```bash
hyperfine --warmup 1 --runs 5 \
  './target/debug-release/vardict -G testdata/hs37d5.fa \
   -b testdata/NA12878.mapped.*.bam -N NA12878 -th 4 \
   -z -c 1 -S 2 -E 3 -g 4 -R "1:100000000-105000000"' \
  2>&1 | tee ./tmp/baseline_hyperfine.txt
```

If `hyperfine` is unavailable, use `/usr/bin/time -v` with 3 manual runs and take the median.

### 2.3 — Measure Java baseline (same workload)

```bash
hyperfine --warmup 1 --runs 3 \
  'java -Xmx8g -classpath VarDictJava/build/libs/VarDict-1.8.3.jar \
   com.astrazeneca.vardict.Main -G testdata/hs37d5.fa \
   -b testdata/NA12878.mapped.*.bam -N NA12878 -th 4 \
   -z -c 1 -S 2 -E 3 -g 4 -R "1:100000000-105000000"' \
  2>&1 | tee ./tmp/baseline_java_hyperfine.txt
```

### 2.4 — Collect hardware counters

```bash
perf stat -e cycles,instructions,cache-references,cache-misses,branches,branch-misses \
  ./target/debug-release/vardict -G testdata/hs37d5.fa \
  -b testdata/NA12878.mapped.*.bam -N NA12878 -th 1 \
  -z -c 1 -S 2 -E 3 -g 4 -R "1:100000000-105000000" \
  > /dev/null 2>&1
```

**Use `-th 1`** for profiling runs — multi-threading adds noise to per-function attribution.

### 2.5 — Record baseline

```markdown
| Metric | Java | Rust (before) | Target |
|--------|------|---------------|--------|
| Wall-clock (5MB, -th 4) | Xs | Ys | ≤Xs |
| Instructions | — | N | — |
| Cache miss rate | — | X% | — |
| Peak RSS | X MB | Y MB | — |
```

Save to `./tmp/perf_baseline_YYYYMMDD.txt`.

**Decision gate**: If Rust already meets the target → STOP.

---

## Phase 3: Profile

Identify WHERE time is spent before deciding WHAT to optimize.

> **Rule**: Never optimize without a profile. Intuition about hotspots is wrong ~70% of the time.

### 3.1 — Choose the profiling tool

| Symptom | Tool | What it shows |
|---------|------|---------------|
| "It's slow but I don't know where" | `perf record` + flamegraph | CPU time per function, call stacks |
| "Specific function is slow" | `perf annotate`, `cargo-asm` | Instruction-level cost |
| "Too many allocations" | DHAT | Allocation sites, peak heap, block counts |
| "Cache thrashing suspected" | `cachegrind` (valgrind) | L1/LL cache miss per line |
| "Lock contention suspected" | `perf lock record` | Mutex/futex wait times |
| "I/O bound suspected" | `strace -c`, `perf trace` | Syscall counts and latency |
| "Compare two versions" | `criterion` benchmarks | Statistical A/B comparison |

### 3.2 — Flamegraph (recommended first step)

```bash
# Record with DWARF call stacks (better for Rust than frame pointers)
perf record -g --call-graph dwarf -F 997 -- \
  ./target/debug-release/vardict -G testdata/hs37d5.fa \
  -b testdata/NA12878.mapped.*.bam -N NA12878 -th 1 \
  -z -c 1 -S 2 -E 3 -g 4 -R "1:100000000-105000000" \
  > /dev/null

# Generate flamegraph
perf script | inferno-collapse-perf | inferno-flamegraph > ./tmp/flamegraph.svg
```

If `inferno` is not installed:
```bash
cargo install inferno
```

Alternative using `flamegraph` crate:
```bash
cargo install flamegraph
flamegraph -o ./tmp/flamegraph.svg -- \
  ./target/debug-release/vardict [flags] -R "1:100000000-105000000" > /dev/null
```

### 3.3 — Read the flamegraph

1. **Width = time**: Wider frames consume more CPU time.
2. **Top of stack = self time**: Functions at the top are doing actual work (not just calling others).
3. **Look for**:
   - Flat plateaus → tight loops, potential vectorization candidates
   - Deep recursive stacks → potential stack overflow or algorithmic issue
   - Allocator frames (`alloc::`, `malloc`, `mmap`) → allocation pressure
   - Hash computation (`hash`, `fx_hash`) → hashing overhead
   - `memcpy` / `memmove` → unnecessary copies

### 3.4 — Cachegrind (data layout issues)

```bash
valgrind --tool=cachegrind --cachegrind-out-file=./tmp/cachegrind.out \
  ./target/debug-release/vardict -G testdata/hs37d5.fa \
  -b testdata/NA12878.mapped.*.bam -N NA12878 -th 1 \
  -z -c 1 -S 2 -E 3 -g 4 -R "1:100000000-101000000" \
  > /dev/null

cg_annotate ./tmp/cachegrind.out | head -100
```

Use a smaller region (1MB) — cachegrind is 20-50x slower.

### 3.5 — Record hotspot findings

```markdown
| Rank | Function/Site | % CPU (or metric) | Category | Notes |
|------|--------------|-------------------|----------|-------|
| 1 | CigarParser::parse_cigar | 35% | Algorithm | Inner loop |
| 2 | HashMap::insert | 12% | Allocation | Per-position maps |
| 3 | Reference::get_sequence | 8% | I/O | Repeated FASTA fetch |
```

---

## Phase 4: Root Cause Analysis

For each hotspot, determine WHY it's slow before jumping to a fix.

### 4.1 — Classification

| Root Cause | Indicators | Typical Fixes |
|-----------|-----------|---------------|
| **Algorithmic** | O(n²) where O(n) exists; redundant computation | Better algorithm, caching, early exit |
| **Data structure** | Wrong container for access pattern | HashMap↔VecMap, Vec↔BTreeMap, sorted vs unsorted |
| **Allocation pressure** | Many small allocs in DHAT; `malloc` in flamegraph | Pool allocator, arena, pre-allocation, stack allocation |
| **Cache unfriendly** | High cache-miss rate in `perf stat`; scattered access | SoA layout, smaller structs, contiguous storage |
| **Branch misprediction** | High branch-miss rate; conditional-heavy inner loop | Branchless arithmetic, likely/unlikely hints, lookup tables |
| **Unnecessary copies** | `memcpy` in flamegraph; `.clone()` in hot path | Borrow instead of clone, `Cow<>`, move semantics |
| **Lock contention** | `futex` / `pthread_mutex` in flamegraph | Reduce critical section, lock-free structures, per-thread state |
| **I/O bound** | `read`/`write` syscalls dominate; low CPU utilization | Buffered I/O, memory-mapped files, prefetching |
| **String formatting** | `fmt::` in flamegraph on output path | Pre-formatted buffers, `itoa`/`ryu` crates |

### 4.2 — Quantify the opportunity

For each root cause, estimate the theoretical speedup using Amdahl's Law:

$$S = \frac{1}{(1 - p) + \frac{p}{s}}$$

Where $p$ = fraction of time in the hotspot, $s$ = speedup factor of the fix.

**Example**: If a function takes 35% of runtime and you can make it 3x faster:
$S = \frac{1}{0.65 + 0.35/3} = \frac{1}{0.767} = 1.30x$

This tells you the maximum possible improvement — if the effort isn't worth 1.30x, skip it.

### 4.3 — Prioritize

Rank by: (estimated speedup) × (confidence) / (implementation effort).

Only pursue optimizations where the expected gain justifies the risk to parity and maintainability.

---

## Phase 5: Design the Optimization

### 5.1 — Optimization hierarchy (prefer top, avoid bottom)

| Level | Type | Risk | Example |
|-------|------|------|---------|
| 1 | **Algorithm** | Low | Replace O(n²) search with O(1) lookup |
| 2 | **Data structure** | Low-Med | HashMap → VecMap for small maps |
| 3 | **Reduce allocation** | Medium | Pre-allocate, reuse buffers, arena |
| 4 | **Improve data layout** | Medium | SoA, padding removal, cache alignment |
| 5 | **Reduce copies** | Medium | Borrow instead of clone, Cow<> |
| 6 | **Concurrency** | High | Finer parallelism, work stealing |
| 7 | **SIMD / intrinsics** | High | Vectorized sequence comparison |
| 8 | **Unsafe** | Very High | Skip bounds checks, raw pointers — last resort |

### 5.2 — Parity impact assessment

Before implementing, answer:
1. Does this change iteration order of any collection? → May break output ordering
2. Does this change floating-point evaluation order? → May break numeric parity
3. Does this change when/where allocations happen? → Usually safe, but measure
4. Does this remove or reorder side effects? → Must preserve Java behavior

If any answer is "yes", document the risk and get user approval before proceeding.

### 5.3 — Write the benchmark FIRST

Before implementing the optimization, write a targeted benchmark:

```bash
# Create benches/<name>_bench.rs using criterion
```

The benchmark should:
- Test the specific operation being optimized
- Test at realistic input sizes (not just trivial cases)
- Include both the current and proposed implementation
- Report throughput (ops/sec) or latency (ns/op)

This gives you a "before" number and a regression test going forward.

---

## Phase 6: Implement

### 6.1 — Make the smallest possible change

- Prefer type alias swaps over widespread refactoring
- Prefer `#[inline]` hints over manual inlining
- Prefer safe Rust over `unsafe` blocks
- One optimization per commit — never bundle

### 6.2 — Common Rust performance patterns

**Avoid unnecessary allocations in hot loops:**
```rust
// BAD: allocates every iteration
for item in items {
    let key = format!("{}-{}", item.chr, item.pos);
    map.entry(key).or_default();
}

// GOOD: reuse buffer
let mut key_buf = String::with_capacity(64);
for item in items {
    key_buf.clear();
    write!(&mut key_buf, "{}-{}", item.chr, item.pos).unwrap();
    map.entry(key_buf.clone()).or_default(); // clone only when needed
}
```

**Use iterators instead of index loops:**
```rust
// BAD: bounds checks on every access
for i in 0..slice.len() {
    total += slice[i];
}

// GOOD: no bounds checks, auto-vectorizes
let total: i64 = slice.iter().sum();
```

**Pre-size collections:**
```rust
// BAD: repeated reallocations
let mut v = Vec::new();
for item in items { v.push(process(item)); }

// GOOD: single allocation
let mut v = Vec::with_capacity(items.len());
for item in items { v.push(process(item)); }

// BEST: use collect (infers capacity from size_hint)
let v: Vec<_> = items.iter().map(process).collect();
```

**Avoid `.clone()` in hot paths — borrow or move instead.**

**Use `#[inline]` on small functions called across crate boundaries.** (Within a crate, the compiler usually inlines automatically.)

---

## Phase 7: Validate

### 7.1 — Run the benchmark

```bash
cargo bench --bench <name>_bench
```

Compare before/after. Criterion reports:
- `improved` / `regressed` / `no change` with confidence intervals
- Percentage change with statistical significance

If no improvement → revert and re-analyze.

### 7.2 — Measure on real workload

```bash
hyperfine --warmup 1 --runs 5 \
  './target/debug-release/vardict [flags] -R "1:100000000-105000000"' \
  2>&1 | tee ./tmp/after_hyperfine.txt
```

### 7.3 — Verify correctness

```bash
cargo test -- --include-ignored
```

### 7.4 — Verify parity

```bash
diff <(./target/debug-release/vardict [flags] -R "1:100000000-105000000") \
     ./tmp/java_reference_output.tsv
```

### 7.5 — Check for regressions on other workloads

Run at least one different region to ensure the optimization doesn't regress elsewhere:
```bash
/usr/bin/time -v ./target/debug-release/vardict [flags] -R "2:50000000-55000000" > /dev/null
```

### 7.6 — Record results

```markdown
| Metric | Before | After | Change | Target | Status |
|--------|--------|-------|--------|--------|--------|
| Wall-clock (5MB) | Xs | Ys | -Z% | ≤Ws | PASS/FAIL |
| Peak RSS | X MB | Y MB | -Z% | — | — |
| Cargo tests | 258 pass | 258 pass | — | all pass | PASS |
| Parity | ✓ | ✓ | — | identical | PASS |
```

**Decision gate**:
- Target met + parity PASS + tests PASS → proceed to Phase 8
- Target not met → profile again (Phase 3), look for next hotspot
- Parity FAIL or test FAIL → revert immediately, diagnose

---

## Phase 8: Prevent Regression

### 8.1 — Add the benchmark to CI (or at minimum keep it in `benches/`)

Criterion benchmarks in `benches/` serve as regression detectors:
```bash
# Run all benchmarks and save baseline
cargo bench --bench <name>_bench -- --save-baseline before

# After future changes, compare
cargo bench --bench <name>_bench -- --baseline before
```

### 8.2 — Document the optimization

Record in the mission work report:
- What was the bottleneck
- What was changed
- Before/after measurements
- Why this approach was chosen over alternatives
- Known limitations or risks

### 8.3 — Update divergence documentation (if applicable)

If the optimization deviates from Java's algorithm or data structures:
- Add doc comment on the changed code
- Update `copilot-office/codebase/CODEBASE.md` "Divergences from Java" section
- Create repo memory: `/memories/repo/<descriptive_name>.md`

---

## Quick Reference: Profiling Commands

| Task | Command |
|------|---------|
| Wall-clock time (statistical) | `hyperfine --warmup 1 --runs 5 './target/debug-release/vardict ...'` |
| Wall-clock time (single) | `/usr/bin/time -v ./target/debug-release/vardict ...` |
| Hardware counters | `perf stat -e cycles,instructions,cache-references,cache-misses,branches,branch-misses -- ./target/debug-release/vardict ...` |
| CPU flamegraph | `perf record -g --call-graph dwarf -F 997 -- ./target/debug-release/vardict ... && perf script \| inferno-collapse-perf \| inferno-flamegraph > fg.svg` |
| Heap profiling | Build with `--no-default-features --features dhat-heap`, run with `-th 1`, open `dhat-heap.json` in DHAT viewer |
| Cache analysis | `valgrind --tool=cachegrind ./target/debug-release/vardict ...` (use small region) |
| Allocation tracing | `valgrind --tool=massif ./target/debug-release/vardict ...` then `ms_print massif.out.*` |
| Syscall profiling | `strace -c ./target/debug-release/vardict ...` |
| Criterion benchmark | `cargo bench --bench <name>_bench` |
| Criterion A/B compare | `cargo bench -- --save-baseline before` then `cargo bench -- --baseline before` |

## Anti-Patterns

| Don't | Why | Do Instead |
|-------|-----|-----------|
| Optimize without profiling | You'll optimize the wrong thing | Always profile first |
| Chase micro-optimizations before algorithmic fixes | 2% gain vs potential 10x gain | Fix the algorithm first |
| Use `unsafe` to skip bounds checks | Risks UB, breaks safety guarantees | Use iterators (no bounds checks, auto-vectorizes) |
| Add `#[inline(always)]` everywhere | Bloats binary, hurts I-cache | Let compiler decide; use `#[inline]` sparingly |
| Benchmark in debug mode | Debug is 10-50x slower, different hot paths | Always `--profile debug-release` |
| Benchmark with other workloads running | Noisy results, irreproducible | Quiesce the system or use `hyperfine` |
| Bundle multiple optimizations in one change | Can't isolate which helped | One optimization per commit |
| Skip parity check after optimization | Output regression goes undetected | Always `diff` against Java reference |

## Lessons Learned (from this project)

| Lesson | Detail |
|--------|--------|
| Measure, don't guess | Deferred seed_map was theoretically sound but measured WORSE (4,004 MB vs 3,563 MB) |
| Type alias swap is the safest optimization | Changing 2 type aliases propagated VecMap through entire pipeline with minimal risk |
| SmallVec can backfire | Inline storage for large values (256-byte Variant) wastes more than it saves |
| Inner container overhead scales with position count | 5M positions × 388 bytes/HashMap = 1.9 GB of pure overhead |
| Flamegraph before DHAT | CPU profile first — if it's not allocation-bound, DHAT is the wrong tool |
| Single-thread for profiling | Multi-threaded profiles are noisy and hard to attribute |
| `mimalloc` and DHAT are mutually exclusive | Must use `--no-default-features --features dhat-heap` for DHAT |
