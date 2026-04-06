---
name: change-impact-review
description: "Assess performance impact of code changes before approval. Use when: reviewing code changes, pre-merge performance gate, benchmark before/after, regression detection, hot path change, allocator change, collection swap, algorithm change, loop modification. Produces PERF_SAFE / PERF_RISK / PERF_REGRESSION verdict."
argument-hint: "Describe the change, e.g. 'SV processor find_match loop restructured' or 'review PR #42 for perf impact'"
---

# Change Impact Review

Pre-merge performance regression gate for VarDict-rs code changes.
Classifies risk, runs benchmarks when needed, and produces a binding verdict.

> **Origin**: Created after a parity fix in `StructuralVariantsProcessor` caused a ~20x runtime regression. The fix was correct (parity passed) but no performance check existed to catch the slowdown before it was applied.

## When to Use

- Before approving any code change to Rust source files
- After the code reviewer completes correctness and style checks
- When the parity orchestrator requests a performance verdict
- See **Caller Context** below for how the verdict is treated depending on who invokes this skill

## When NOT to Use

- Documentation-only changes (`.md`, comments with no logic change)
- Test-only changes (`#[cfg(test)]` modules, `tests/` directory)
- CI/build configuration changes (`Cargo.toml` dependency bumps without feature changes)
---

## Caller Context

This skill is used by two agents in different modes. The mode affects how the verdict is treated.

### Self-Assessment Mode (rust-implementer)
When the `rust-implementer` loads this skill during Step 5 (Verify Compilation and Performance):
- Classify risk using the decision tree
- Benchmark only for MEDIUM or HIGH risk
- Include the verdict in your implementation report as **advisory**
- The code-reviewer's independent verdict takes precedence

### Independent Gate Mode (code-reviewer)
When the `code-reviewer` loads this skill during Section 3 (Performance Impact):
- Benchmark is **mandatory** for MEDIUM or HIGH risk — do not skip
- Your verdict is **binding** — it determines whether the change is approved
- `PERF_REGRESSION` **blocks approval** and must be escalated via the orchestrator
- If the rust-implementer already included an advisory verdict, review it but produce your own independent classification

## Step 1: Classify Risk

Determine the performance risk level of the change using this decision tree.

```
Is the change in a hot-path module?
├─ YES → Is there an algorithm or data structure change?
│        ├─ YES → HIGH
│        └─ NO  → Is there a new allocation in a loop? New clone()?
│                  ├─ YES → MEDIUM
│                  └─ NO  → LOW
└─ NO  → Is there a new allocation pattern or collection type change?
          ├─ YES → MEDIUM
          └─ NO  → LOW
```

### Hot-Path Modules

These modules execute per-read or per-base — small changes have large impact:

| Module | File | Frequency | Benchmark |
|--------|------|-----------|-----------|
| CigarParser | `src/mods/cigar_parser.rs` | Per-read | `cigar_parser_bench` |
| VariationRealigner | `src/mods/variant_realigner.rs` | Per-region | `pipeline_bench` |
| StructuralVariantsProcessor | `src/mods/structural_variants_processor.rs` | Per-region | `pipeline_bench` |
| ToVarsBuilder | `src/mods/to_vars_builder.rs` | Per-position | `pipeline_bench` |
| Pipeline (record processing) | `src/mods/pipeline.rs` | Per-read | `pipeline_bench` |
| VecMap / InnerMap | `src/data/vecmap.rs` | Per-position | `vecmap_vs_hashmap_bench` |

### Risk Signals

| Signal | Risk | Example |
|--------|------|---------|
| Loop restructure in hot module | HIGH | Changing iteration order in `parse_cigar` |
| New `HashMap`/`IndexMap` in hot loop | HIGH | Adding a lookup table inside per-read processing |
| Algorithm complexity change (O(n)→O(n²)) | HIGH | Nested search replacing direct lookup |
| Collection type swap | MEDIUM | `HashMap` → `BTreeMap`, `Vec` → `LinkedList` |
| New `.clone()` on hot path | MEDIUM | Cloning a `String` or `Vec` per iteration |
| Added `.collect()` mid-pipeline | MEDIUM | Materializing an iterator unnecessarily |
| New branch in cold path | LOW | Additional config check at startup |
| Formatting / output change | LOW | Changing print format in `OutputVariant` |
| Refactor with same logic | LOW | Extracting method, renaming, reordering fields |

---

## Step 2: Act on Risk Level

### LOW Risk

No benchmark required. Proceed with standard code review.

Record in verdict:
```
**Risk**: LOW — {one-line reason}
**Verdict**: PERF_SAFE
```

### MEDIUM Risk

Run the relevant benchmark (from the Hot-Path Modules table above). If no specific benchmark covers the change, use `pipeline_bench` as the integration-level fallback.

```bash
conda activate rust_build_env
export LIBCLANG_PATH=$CONDA_PREFIX/lib

# 1. Benchmark BEFORE the change (stash or checkout base)
git stash  # or: git checkout <base-branch>
cargo bench --profile debug-release --bench <bench_name> -- --save-baseline before

# 2. Benchmark AFTER the change
git stash pop  # or: git checkout <change-branch>
cargo bench --profile debug-release --bench <bench_name> -- --baseline before
```

Evaluate results:
- **≤5% regression**: PERF_SAFE (within noise floor)
- **5–15% regression**: PERF_RISK (document and proceed with acknowledgment)
- **>15% regression**: PERF_REGRESSION (block approval)

### HIGH Risk

Run **both** module benchmark and wall-clock measurement:

```bash
# 1. Module benchmark (same as MEDIUM)
cargo bench --profile debug-release --bench <bench_name> -- --save-baseline before
# ... apply change ...
cargo bench --profile debug-release --bench <bench_name> -- --baseline before

# 2. Wall-clock on representative workload
# BEFORE:
hyperfine --warmup 1 --runs 3 \
  './target/debug-release/vardict -G testdata/hs37d5.fa \
   -b testdata/NA12878.mapped.*.bam -N NA12878 -th 4 \
   -z -c 1 -S 2 -E 3 -g 4 -R "1:100000000-105000000"' \
  2>&1 | tee ./tmp/perf_before.txt

# AFTER:
hyperfine --warmup 1 --runs 3 \
  './target/debug-release/vardict -G testdata/hs37d5.fa \
   -b testdata/NA12878.mapped.*.bam -N NA12878 -th 4 \
   -z -c 1 -S 2 -E 3 -g 4 -R "1:100000000-105000000"' \
  2>&1 | tee ./tmp/perf_after.txt
```

Same thresholds apply (≤5% / 5–15% / >15%).

---

## Step 3: Produce Verdict

Include this block in the code review output:

```
### Performance Verdict

**Risk**: {HIGH | MEDIUM | LOW}
**Modules touched**: {list of hot-path modules, or "none (cold path)"}
**Benchmark**: {bench name and result, or "not required (LOW risk)"}
**Wall-clock**: {before → after, or "not required"}
**Verdict**: {PERF_SAFE | PERF_RISK | PERF_REGRESSION}
**Evidence**: {one-line summary, e.g. "pipeline_bench: +2.3% (within noise)"}
```

### Verdict Definitions

| Verdict | Meaning | Action |
|---------|---------|--------|
| `PERF_SAFE` | No measurable regression (≤5%) or LOW risk change | Approve (performance aspect) |
| `PERF_RISK` | Measurable regression 5–15%, or MEDIUM risk without benchmark | Approve with documented acknowledgment; user should be notified |
| `PERF_REGRESSION` | Regression >15%, or HIGH risk without benchmark | **Block approval**. Escalate via `perf-optimization` skill |

### Escalation

When verdict is `PERF_REGRESSION`:

1. Do NOT approve the change
2. Report to the orchestrator with the benchmark evidence
3. The orchestrator decides next action:
   - **Redesign**: Ask implementer for an alternative approach
   - **Profile**: Invoke `perf-optimization` skill for root cause analysis
   - **User decision**: Present the trade-off (correctness vs performance) to the user

---

## Quick Reference: Benchmark Commands

| Benchmark | Command | Covers |
|-----------|---------|--------|
| CigarParser | `cargo bench --profile debug-release --bench cigar_parser_bench` | Per-read CIGAR parsing |
| Pipeline (integration) | `cargo bench --profile debug-release --bench pipeline_bench` | Full region processing |
| VecMap vs HashMap | `cargo bench --profile debug-release --bench vecmap_vs_hashmap_bench` | Collection operations |
| Wall-clock (5MB) | `hyperfine --warmup 1 --runs 3 './target/debug-release/vardict ... -R "1:100000000-105000000"'` | End-to-end runtime |
| Wall-clock (1MB) | `hyperfine --warmup 1 --runs 5 './target/debug-release/vardict ... -R "1:100000000-101000000"'` | Quick sanity check |

## Anti-Patterns

| Don't | Why | Do Instead |
|-------|-----|------------|
| Skip benchmarks for "trivial" hot-path changes | The SV processor incident was a "trivial" restructure | Always benchmark MEDIUM+ changes |
| Benchmark only in debug mode | 10-50x slower, different hot paths | Always `--profile debug-release` |
| Compare single runs | High variance | Use criterion (`--save-baseline`) or `hyperfine` (≥3 runs) |
| Accept >15% regression "because parity requires it" | There is almost always an alternative implementation | Escalate and redesign |
| Run benchmarks with other heavy processes | Noisy results | Quiesce or use `hyperfine` |
