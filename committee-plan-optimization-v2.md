# Committee Plan: VarDict-rs Data-Driven Optimization v2

**Date:** 2026-04-02
**Committee:** Opus (Claude Opus 4.6), Sonnet (Claude Sonnet 4.6), Gemini (Gemini 3.1 Pro), GPT (GPT 5.4)
**Rounds:** 1 discussion round (all contested points resolved in round 1)
**Basis:** Phase A profiling data — flamegraph, DHAT, perf stat, wall-clock baselines

---

## Plan: Allocation-First Optimization With Parity Safety

VarDict-rs profiling reveals an **allocation-dominated, memory-bound** workload: 2.5 GB total allocations per 1MB region, 21.3% cache miss rate, IPC 1.39. The top optimization lever is the **VarDesc lifecycle** (to_key_string + clone + drop + eq = 14% self-time) and **VecMap allocation patterns** (1.2 GB = 48% of all allocations). Thread scaling failure (6% at 4 threads) is a single-region benchmark artifact — not a code bug.

Target outcome: **20-30% wall-clock reduction** on the 1MB single-threaded benchmark, **~1.5 GB RSS reduction**, zero parity drift.

**Steps**

### Phase B: P0 — Zero-Risk Allocation Elimination

1. **Rewrite `is_same_variation_on_ref` to zero-allocation** (*no dependencies*)
   - Replace `HashSet<String>` + `to_key_string()` per position with direct `VarDesc` pattern matching
   - Check `vars_at_pos.len() == 1 && insertion_vars.is_none()`, then match the single `VarDesc::SNV { ref_base }` against `reference.get(position)`
   - **Impact:** ~1-5% CPU, ~110 MB RSS reduction, eliminates 991K temporary HashSets
   - **Files:** [vardict_pipeline.rs](src/mods/vardict_pipeline.rs) (`is_same_variation_on_ref` ~line 2955)
   - **Risk:** Zero — mathematically equivalent rewrite

2. **Eliminate VecMap clone per position in `run_to_vars_builder`** (*no dependencies*)
   - Replace `non_insertion_vars.get(&position).cloned().unwrap_or_default()` with borrowing: `non_insertion_vars.get(&position).unwrap_or(&EMPTY_RAW_VAR_MAP)` where `EMPTY_RAW_VAR_MAP` is a static default
   - Downstream consumers (`is_same_variation_on_ref`, `calc_hicov`, `create_variant_records`) already take `&RawVarMap`
   - If `create_insertion_records` mutations conflict, use `HashMap::remove_entry` to move out owned data or apply `Cow<RawVarMap>`
   - Treat as a bounded design spike: if borrow checker rejects, abandon and keep clone — bank other wins first
   - **Impact:** ~2-5% CPU, ~244 MB RSS reduction, eliminates 991K VarDesc clones + drops
   - **Files:** [vardict_pipeline.rs](src/mods/vardict_pipeline.rs) (`run_to_vars_builder` ~line 2685)
   - **Risk:** Medium — verify create_insertion_records doesn't mutate current position's entry

3. **Move qname String allocation to error path only** (*no dependencies, parallel with steps 1-2*)
   - In `CigarParser::process_record()`, `String::from_utf8_lossy(record.qname()).to_string()` is called on EVERY record
   - Move into the error/trace branches only — the happy path doesn't need qname
   - **Impact:** ~0.5-1% CPU (millions of eliminated per-record allocations)
   - **Files:** [cigar_parser.rs](src/mods/cigar_parser.rs) (`process_record` ~line 301)
   - **Risk:** Zero — error-path-only change

4. **Pre-size VecMap and outer HashMap** (*no dependencies, parallel with steps 1-2*)
   - Add `VecMap::with_capacity(4)` for inner VecMaps (most positions have 1-4 variants)
   - Add `HashMap::with_capacity(region.len())` for `non_insertion_vars` and `insertion_vars` in CigarParser
   - **Impact:** ~2-4% CPU, reduces 967.8 MB #1 DHAT site by ~60% (fewer reallocations)
   - **Files:** [vec_map.rs](src/utils/vec_map.rs), [cigar_parser.rs](src/mods/cigar_parser.rs) (`get_non_insertion_variant`)
   - **Risk:** Zero — pre-sizing doesn't change behavior

### Phase C: P1 — VarDesc Structural Optimization

5. **Implement `Ord` on VarDesc for zero-allocation sorting** (*depends on step 1 for parity baseline*)
   - Implement `PartialOrd` + `Ord` on `VarDesc` producing the exact same ordering as `sort_by(to_key_string())`
   - Implementation strategy: `key_ordering_tuple()` method returning `(type_ordinal, &payload)` that can be property-tested against `to_key_string()` lexicographic ordering
   - Replace all `sort_by_cached_key(|desc| desc.to_key_string())` with `sort()` — zero allocations
   - **Impact:** ~3-5% CPU (eliminates millions of sort-related String allocations)
   - **Files:** [variants.rs](src/variants/variants.rs) (`VarDesc`), [vardict_pipeline.rs](src/mods/vardict_pipeline.rs) (3 sort sites), [variant_realigner.rs](src/mods/variant_realigner.rs) (6 lookup sites), [structural_variants_processor.rs](src/mods/structural_variants_processor.rs) (3 lookup sites)
   - **Risk:** Medium — must prove sort order matches string sort exactly via property tests + parity

6. **Add non-allocating helper methods on VarDesc** (*parallel with step 5*)
   - `is_deletion_key() -> bool` — replaces `to_key_string().starts_with('-')`
   - `key_equals(s: &str) -> bool` — replaces `to_key_string() == target` comparisons in control-flow paths
   - `is_ref_snv(ref_base: u8) -> bool` — replaces `to_key_string() == ref_base_char` checks
   - Preserve `to_key_string()` for final output (`description_string` in variant records)
   - **Impact:** ~3-5% CPU (eliminates control-flow string allocations across realigner, SV processor, pipeline)
   - **Files:** [variants.rs](src/variants/variants.rs), then update callers in [variant_realigner.rs](src/mods/variant_realigner.rs), [structural_variants_processor.rs](src/mods/structural_variants_processor.rs), [vardict_pipeline.rs](src/mods/vardict_pipeline.rs)
   - **Risk:** Low — behavior-preserving pattern matching

7. **Reduce clone churn in `create_variant_records` / `create_insertion_records`** (*depends on step 5*)
   - The `BTreeMap<String, RawVariant>` aggregation clones RawVariant on first insert
   - Change merge strategy: initialize entry without cloning when possible, accumulate in place
   - Materialize ordered key view once per position, reuse across merge and output preparation
   - **Impact:** ~1-3% CPU, modest additional allocation reduction
   - **Files:** [vardict_pipeline.rs](src/mods/vardict_pipeline.rs) (`create_variant_records`, `create_insertion_records`)
   - **Risk:** Low — internal refactoring of variant assembly

### Phase D: P2 — Seed Map & Locality

8. **Optimize `Reference::build_seed_map` key representation** (*no dependencies*)
   - Replace `HashMap<Vec<u8>, Vec<i64>>` keys with fixed-size `[u8; N]` arrays (SEED_1 and SEED_2 are compile-time constants)
   - Pre-size HashMap from region length estimate
   - Consider `SmallVec<[i64; 4]>` for position values
   - **Impact:** ~2-5% CPU, ~200 MB RSS reduction, improved cache locality
   - **Files:** [reference.rs](src/data/reference.rs) (`build_seed_map`, `ReferenceSeedMap`)
   - **Risk:** Low — lookup-only structure, not output-affecting

9. **REPROFILE GATE** (*depends on steps 1-8*)
   - Re-run flamegraph, perf stat, and DHAT on the same 1MB and 5MB workloads
   - Compare against Phase A baselines
   - Only proceed to Phase E if `ref_coverage`, `parse_cigar`, or `add_variation_for_matching_part` still dominate
   - Gate also determines whether SmallVec backing for VecMap is needed (measure `size_of::<(VarDesc, RawVariant)>()`)

### Phase E: Conditional — Evidence-Gated After Reprofile

10. **Replace `ref_coverage` HashMap with dense Vec** (*first item if cache miss rate still >15%*)
    - `HashMap<i64, usize>` → `Vec<u32>` indexed by `(pos - region.start)`
    - Pad to `region.len() + max_read_len` for reads extending past region
    - **Impact:** ~1-6% CPU depending on post-Phase-D cache pressure
    - **Files:** [cigar_parser.rs](src/mods/cigar_parser.rs), [vardict_pipeline.rs](src/mods/vardict_pipeline.rs)
    - **Risk:** Medium — audit all ref_coverage consumers

11. **SmallVec backing for VecMap** (*only if `size_of::<(VarDesc, RawVariant)>()` ≤ 96 bytes*)
    - Back `VecMap<K,V>` with `SmallVec<[(K,V); 2]>` for inline storage
    - If entry size 96-192 bytes → inline 1; if >192 bytes → skip
    - **Impact:** ~800-950 MB further RSS reduction if applicable
    - **Files:** [vec_map.rs](src/utils/vec_map.rs)
    - **Risk:** Zero (parity-neutral internal structure change)

12. **catch_unwind removal** (*only if parse_cigar still dominates after Phase D*)
    - Gate behind `#[cfg(feature = "catch-panics")]` feature flag (default off)
    - Consider `panic = "abort"` in deployment profile
    - **Impact:** ~2-5% CPU (LLVM optimization unlock)
    - **Files:** [cigar_parser.rs](src/mods/cigar_parser.rs) (`process_record`)

13. **Thread scaling verification with multi-region BED** (*parallel with any phase*)
    - Create 100×10KB BED file from the 1MB region and benchmark 1/2/4/8 threads
    - Verify existing disabled simple-mode prefetch path
    - Only implement cost-based scheduling if scaling proves poor on multi-region workloads
    - **Files:** [parallel_pipeline.rs](src/mods/parallel_pipeline.rs), [vardict.rs](src/bin/vardict.rs)

**Relevant files**
- `src/variants/variants.rs` — VarDesc struct, to_key_string, Ord impl, new helpers
- `src/mods/vardict_pipeline.rs` — run_to_vars_builder, is_same_variation_on_ref, create_variant_records, create_insertion_records, sort sites
- `src/mods/cigar_parser.rs` — process_record (qname, catch_unwind), ref_coverage, VecMap pre-sizing, get_non_insertion_variant
- `src/utils/vec_map.rs` — VecMap with_capacity, potential SmallVec backing
- `src/data/reference.rs` — build_seed_map, seed key representation
- `src/mods/variant_realigner.rs` — find_desc_by_key_string, to_key_string callers
- `src/mods/structural_variants_processor.rs` — get_or_create_variation, to_key_string callers
- `src/mods/parallel_pipeline.rs` — thread scaling, BAM reader reuse

**Verification**

For every optimization step:
1. `cargo test --profile debug-release` — all tests pass (0 failures)
2. Focused fixture tests: `cargo test --profile debug-release cigar_parser_fixture_test variant_realigner_fixture_test structural_variants_fixture_test tovars_fixture_test`
3. Parity byte-diff on chr1, chr20, chr22:
   ```bash
   bash copilot-office/m9-option-parity/copilot-desk/scripts/na12878_parity_v2.sh \
     --chr 1 --chr 20 --chr 22 --opts "" --opts-label default --rust-only --no-build --parallel 4 --no-stop
   ```
4. `perf stat -e cycles,instructions,cache-references,cache-misses,branches,branch-misses` before/after with `-th 1`
5. DHAT before/after for allocation-related changes (steps 1-4, 7-8, 10-11)
6. Flamegraph re-run after completing each Phase (B, C, D)
7. `VarDesc::to_key_string`, `VarDesc::clone`, and `is_same_variation_on_ref` must materially drop in flamegraph before parser-level work (Phase E) is opened
8. Regression guard: Rust wall-clock must remain ≤0.7× Java on benchmark workloads

**Decisions**
- `to_key_string()` is treated as an **output-formatting primitive**, not a control-flow primitive. Goal: stop using it in hot lookups, classifiers, and sorts. Keep it for final `description_string` emission.
- VecMap clone → borrow is treated as a **bounded design spike**: try once, abandon if borrow-checker rejects. Bank lower-risk wins first.
- CigarParser state pooling is **dropped** — pre-sizing outer HashMap captures the same 200 MB benefit with a one-line change.
- SmallVec backing for VecMap is **gated on size_of measurement** — only if inline entries fit within 96 bytes.
- Thread scaling is a **workload-shape issue** — no architecture change until multi-region scaling is verified poor.
- `SharedReference::load_chromosomes` (480 MB) is **cold-start cost** — not addressed unless success metric expands to process startup latency.
- `catch_unwind` removal is **deferred** until post-reprofile evidence confirms it's still material.

**Explicit No-Go List**
- Do NOT change `java_hashmap_iteration_order*` or `java_hashmap_bucket_index` functions
- Do NOT change float formatting/rounding (`round_half_even`, `java_format_double`)
- Do NOT switch allocators (mimalloc is working well)
- Do NOT add stateful reverse-deletion lookup caches in `variant_realigner.rs`
- Do NOT use `unsafe` code until all safe alternatives are proven insufficient
- Do NOT reintroduce chromosome-wide seed scans in `structural_variants_processor.rs`
- Do NOT change VarDesc Hash/Eq semantics (used as map key everywhere)
- Do NOT introduce interior mutability (OnceCell/RefCell) into VarDesc
- Do NOT attempt intra-region parallelism until all allocation wins are shipped and measured

---

## Phase A Profiling Baselines (Reference)

| Metric | Value |
|--------|-------|
| Wall-clock 1MB -th1 | 2.74s median |
| Wall-clock 1MB -th4 | 2.59s median |
| Wall-clock 5MB -th1 | 18.12s median |
| Peak RSS 1MB -th1 | ~1.64 GB |
| Peak RSS 5MB -th1 | ~6.87 GB |
| Total allocations | 2,520 MB in 13M blocks |
| IPC | 1.39 |
| Cache miss rate | 21.3% |
| L1-dcache miss rate | 4.61% |
| Branch miss rate | 0.50% |
| #1 CPU hotspot | VarDesc::to_key_string (9.0% self) |
| #1 DHAT hotspot | VecMap::VacantEntry::insert (967.8 MB) |

---

## Provenance

How consensus was reached:

| Category | Count |
|----------|-------|
| Full consensus (4/4) | 5 points (to_key_string target, is_same_variation rewrite, thread scaling diagnosis, no-go list, verification protocol) |
| Strong consensus (3/4) | 2 points (VecMap clone elimination now, seed map optimization) |
| Resolved via discussion | 5 points (to_key_string approach, VecMap insert opt, CigarParser pooling, ref_coverage priority, catch_unwind/qname) |
| Resolved by supermajority | 3 points (drop pooling, ref_coverage in Phase C, qname now / catch_unwind deferred) |
| Resolved by Chief decision | 0 points |
| Escalated to user | 0 points |

Unique contributions adopted: 6 (VarDesc Ord, static empty RawVarMap, re-profile gate, create_variant_records reduction, verify prefetch path, bounded design spike for VecMap clone)
Unique contributions modified: 1 (SmallVec → gated on size_of)
Unique contributions rejected: 3 (OnceCell caching, pre-compute pairs, region sub-splitting)

## Discussion Log

### Round 1
**Contested point resolution:**
- **CP1 (to_key_string approach):** All 4 members converged on Position A — `Ord` + non-allocating helpers with `key_ordering_tuple()` implementation. Opus conceded Display+scratch is inferior. Gemini dropped pre-compute pairs in favor of Ord. GPT added one-time downstream string materialization.
- **CP2 (VecMap insert):** All 4 converged on `with_capacity(4)` + outer HashMap pre-sizing as immediate action. SmallVec gated on `size_of` measurement — Sonnet proposed conditional capacity thresholds (≤96B→inline 2, 96-192B→inline 1, >192B→skip).
- **CP3 (CigarParser pooling):** 3/4 DROP (Opus, Sonnet conceded, GPT). Gemini dissented (keep at P1). Supermajority resolution: DROP — pre-sizing captures the benefit.
- **CP4 (ref_coverage):** 3/4 agree: Phase C, gated on re-profile (Opus, Gemini, GPT). Sonnet accepted Phase C timing but wanted no gate. Supermajority resolution: Phase C first item but with re-profile gate.
- **CP5 (catch_unwind/qname):** 3/4 agree: qname to error path in Phase B, catch_unwind deferred (Opus, Sonnet, GPT). Gemini would deprioritize both. Supermajority resolution: split — qname now, catch_unwind later.
