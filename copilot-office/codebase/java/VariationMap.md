# VariationMap (collection package)

**Source**: `collection/VariationMap.java`, `collection/ConcurrentHashSet.java`, `collection/DirectThreadExecutor.java`, `collection/Tuple.java`  
**LOC**: 62 (VariationMap) + 88 (ConcurrentHashSet) + 5 (DirectThreadExecutor) + 45 (Tuple) = 200 total  
**Rust counterpart**: `VecMap` (inner map), `StructuralVariantCounts` (SV sidecar), `HashSet` (ConcurrentHashSet), native tuples (Tuple)  
**Status**: complete

## Overview

The `collection` package contains four top-level classes and three inner classes that serve as fundamental building blocks for VarDict's data pipeline:

- **VariationMap<K,V>**: The most parity-critical class. Extends `LinkedHashMap<K,V>` and adds an embedded `SV` (structural variant counters) sidecar field. Used as the **inner map type** for `nonInsertionVariants` and `insertionVariants` at every position. Its insertion-order iteration (from LinkedHashMap) determines the order variants are processed in ToVarsBuilder.
- **VariationMap.SV**: Static inner class storing four structural variant fields (type, pairs, splits, clusters). Attached to a VariationMap at positions where SV evidence exists.
- **ConcurrentHashSet<E>**: A `Set<E>` backed by `ConcurrentHashMap<E, Boolean>`. Used only in SomaticMode and AmpliconMode for the `splice` set.
- **DirectThreadExecutor**: An `Executor` that runs tasks synchronously on the calling thread. Used whenever VarDict needs to run a `CompletableFuture` pipeline stage inline rather than asynchronously.
- **Tuple / Tuple2 / Tuple3**: Generic tuple utility classes with immutable public fields `_1`, `_2`, `_3`. Used for BED file parsing, sample names, amplicon position lists, and variant+description pairing.

## Class Inventory

| Class | File | Lines | Generic Types | Extends/Implements | Parity Risk |
|-------|------|-------|---------------|-------------------|-------------|
| `VariationMap<K,V>` | VariationMap.java:13-63 | 51 | K=String, V=Variation | `LinkedHashMap<K,V>` | **HIGH** |
| `VariationMap.SV` | VariationMap.java:19-24 | 6 | none | none | MEDIUM |
| `ConcurrentHashSet<E>` | ConcurrentHashSet.java:12-99 | 88 | E (typically String) | `Set<E>` | LOW |
| `DirectThreadExecutor` | DirectThreadExecutor.java:5-9 | 5 | none | `Executor` | NONE |
| `Tuple` (outer) | Tuple.java:6-50 | 45 | none (factory) | none | LOW |
| `Tuple.Tuple2<T1,T2>` | Tuple.java:19-32 | 14 | T1, T2 | none | LOW |
| `Tuple.Tuple3<T1,T2,T3>` | Tuple.java:34-48 | 15 | T1, T2, T3 | none | LOW |

## Method Inventory

| Method | Lines | Analyzed? | Summary |
|--------|-------|-----------|---------|
| `VariationMap.getSV()` | L32-L49 | yes | Lazy-init VariationMap + SV sidecar + "SV" sentinel key |
| `VariationMap.removeSV()` | L56-L62 | yes | Null out sv, remove "SV" sentinel key (no null check!) |
| `ConcurrentHashSet.*` | L12-L99 | yes | All methods delegate to ConcurrentHashMap |
| `DirectThreadExecutor.execute()` | L7-L8 | yes | Synchronous Runnable execution on calling thread |
| `Tuple.tuple(f,s)` | L11-L13 | yes | Factory for Tuple2 |
| `Tuple.tuple(f,s,t)` | L15-L17 | yes | Factory for Tuple3 |
| `Tuple2.newTuple(f,s)` | L28-L30 | yes | Redundant factory for Tuple2 |
| `Tuple3.newTuple(f,s,t)` | L43-L45 | yes | Redundant factory for Tuple3 |

## Field-by-Field Analysis

### VariationMap<K,V>

| Field | Type | Default | Writers | Readers |
|-------|------|---------|---------|---------|
| `sv` | `SV` (nullable) | `null` | `getSV()` (lazy-init), `removeSV()` (sets null), StructuralVariantsProcessor, VariationRealigner | ToVarsBuilder (`.sv == null` checks, `.sv.type`), SimplePostProcessModule, SomaticPostProcessModule |

**Inherited from LinkedHashMap**: All standard Map operations. Key behavior: **insertion-order iteration**.

**Critical**: VariationMap has NO additional instance fields beyond `sv`. It is literally LinkedHashMap + one nullable SV pointer. All the map operations (get, put, containsKey, remove, keySet, entrySet) are inherited directly from LinkedHashMap.

### VariationMap.SV

| Field | Type | Default (Java) | Writers | Readers |
|-------|------|-----------------|---------|---------|
| `type` | `String` | `null` | StructuralVariantsProcessor (assigns SV type strings like "DEL", "DUP", "INV", etc.) | ToVarsBuilder line 132 (`varsAtCurPosition.sv.type`), when there's no coverage |
| `pairs` | `int` | `0` | StructuralVariantsProcessor (increment/assign discordant pair counts) | ToVarsBuilder line 443 (serializes to `splits-pairs-clusters`) |
| `splits` | `int` | `0` | StructuralVariantsProcessor (increment/assign split-read counts) | ToVarsBuilder line 443 |
| `clusters` | `int` | `0` | StructuralVariantsProcessor (increment/assign cluster counts) | ToVarsBuilder line 443 |

### ConcurrentHashSet<E>

| Field | Type | Default | Writers/Readers |
|-------|------|---------|-----------------|
| `_map` | `ConcurrentHashMap<E, Boolean>` | new empty map | All Set operations delegate via keySet() |

### DirectThreadExecutor

No fields. Stateless.

### Tuple (outer)

No fields. Private constructor — factory only.

### Tuple.Tuple2<T1,T2>

| Field | Type | Default | Mutability |
|-------|------|---------|------------|
| `_1` | `T1` | constructor arg | `final` (immutable) |
| `_2` | `T2` | constructor arg | `final` (immutable) |

### Tuple.Tuple3<T1,T2,T3>

| Field | Type | Default | Mutability |
|-------|------|---------|------------|
| `_1` | `T1` | constructor arg | `final` (immutable) |
| `_2` | `T2` | constructor arg | `final` (immutable) |
| `_3` | `T3` | constructor arg | `final` (immutable) |

## Method Analyses

### VariationMap.getSV() — PARITY-CRITICAL

**Source**: VariationMap.java:L32-L49  
**Signature**: `public static SV getSV(Map<Integer, VariationMap<String, Variation>> hash, int start)`  
**Called by**: StructuralVariantsProcessor (15 call sites), VariationRealigner (5 call sites)  
**Purpose**: Lazily initialize both the VariationMap entry and its SV sidecar. Returns the SV for mutation.

**Algorithm (Step-by-Step)**:

1. `VariationMap<String, Variation> map = hash.get(start);`
   - Looks up the outer map (typically `nonInsertionVariants`) at the given position.

2. **If map == null** (no entry at this position):
   - Creates `new VariationMap<>()` (empty LinkedHashMap with sv=null)
   - `hash.put(start, map)` — inserts it into the outer map
   - **State change**: new VariationMap entry created at position `start`

3. `SV sv = map.sv;`
   - Reads the SV sidecar from the map.

4. **If sv == null** (no SV data yet at this position):
   - Creates `new VariationMap.SV()` — all fields default (type=null, pairs=0, splits=0, clusters=0)
   - `map.sv = sv;` — attaches SV to the VariationMap
   - `map.put("SV", new Variation());` — **ALSO** puts a sentinel "SV" entry into the map itself
   - **State change**: sv field set, "SV" key added to map with empty Variation

5. **Regardless of whether sv was null**:
   - `if (!map.containsKey("SV"))` → `map.put("SV", new Variation());`
   - This is a **defensive guard**: ensures the "SV" sentinel key exists even if sv was already non-null but someone removed the "SV" key.

6. `return sv;` — returns the (possibly newly created) SV struct for mutation by callers.

**Null/Edge Cases**:
- If `hash` is null → NullPointerException (never happens in practice — hash is always initialized as HashMap)
- If `start` maps to a VariationMap whose sv is already set but "SV" key was removed → step 5 re-adds "SV" key
- The returned SV is mutable; callers increment its fields directly (e.g., `sv.splits++`, `sv.pairs++`)

**Parity-Critical Detail**: The "SV" sentinel Variation is placed into the **LinkedHashMap** at the point of getSV() call. This means iterating the VariationMap's entries will encounter the "SV" key **in the order it was inserted relative to other variants**. In ToVarsBuilder.createVariant(), the "SV" key is explicitly checked: `if (descriptionString.equals("SV")) { ... continue; }` — so it's filtered during variant creation but its **insertion order affects keySet() ordering**.

### VariationMap.removeSV() — PARITY-CRITICAL

**Source**: VariationMap.java:L56-L62  
**Signature**: `public static void removeSV(Map<Integer, VariationMap<String, Variation>> hash, int start)`  
**Called by**: VariationRealigner (3 call sites: lines 1333, 1518, 1524)  
**Purpose**: Remove SV data from a position (used during realignment when SV evidence is invalidated).

**Algorithm (Step-by-Step)**:

1. `hash.get(start).sv = null;`
   - Nulls out the SV sidecar.
   - **WARNING**: No null-check on `hash.get(start)`. If no entry exists at `start`, this throws NullPointerException! In practice, removeSV() is only called at positions where getSV() was previously called, so the entry always exists.

2. `if (hash.get(start).containsKey("SV"))` → `hash.get(start).remove("SV");`
   - Removes the "SV" sentinel key from the LinkedHashMap.
   - **Note**: calls `hash.get(start)` three times — no caching into local. Safe because HashMap.get() is deterministic.

**State Changes**:
- `map.sv` → null
- "SV" key removed from map (if present)
- The empty Variation that was under "SV" is discarded

**Null/Edge Cases**:
- If `hash.get(start)` is null → NullPointerException (caller bug — never happens in practice)
- If "SV" key was already removed → no-op on the containsKey/remove

### ConcurrentHashSet — All Methods

**Source**: ConcurrentHashSet.java:L12-L99

All methods delegate to `_map` (ConcurrentHashMap) or `_map.keySet()`:

| Method | Implementation | Return Value |
|--------|----------------|--------------|
| `size()` | `_map.size()` | int |
| `isEmpty()` | `_map.isEmpty()` | boolean |
| `contains(Object)` | `_map.containsKey(e)` | boolean |
| `iterator()` | `_map.keySet().iterator()` | Iterator<E> |
| `toArray()` | `_map.keySet().toArray()` | Object[] |
| `toArray(T[])` | `_map.keySet().toArray(a)` | T[] |
| `add(E)` | `_map.put(e, Boolean.TRUE) == null` | true if newly added |
| `remove(Object)` | `_map.remove(o) != null` | true if was present |
| `containsAll(Collection)` | `_map.keySet().containsAll(c)` | boolean |
| `addAll(Collection)` | loop + add each | true if any were new |
| `retainAll(Collection)` | `_map.keySet().retainAll(c)` | boolean |
| `removeAll(Collection)` | `_map.keySet().removeAll(c)` | boolean |
| `clear()` | `_map.clear()` | void |
| `toString()` | `_map.keySet().toString()` | String |
| `hashCode()` | `_map.keySet().hashCode()` | int |

**Note**: `equals()` is NOT overridden — relies on Object identity. This is a potential bug but doesn't affect parity since ConcurrentHashSet is only used as a local variable within mode processing.

### DirectThreadExecutor.execute()

**Source**: DirectThreadExecutor.java:L7-L8
```java
public void execute(Runnable command) {
    command.run();
}
```

Runs the Runnable synchronously on the calling thread. No thread creation. This ensures CompletableFuture pipelines execute sequentially when used with this executor.

**Used by**: SimpleMode, SomaticMode, AmpliconMode, SplicingMode, StructuralVariantsProcessor, VariationRealigner, SomaticPostProcessModule — whenever a sub-pipeline is invoked within a pipeline stage ("re-run pipeline on a remote region").

**Rust equivalent**: In Rust, since the pipeline is already synchronous (no CompletableFuture), this pattern is eliminated entirely — just call the function directly.

### Tuple.tuple(f, s)

**Source**: Tuple.java:L11-L13  
Factory method returning `new Tuple2<>(f, s)`.

### Tuple.tuple(f, s, t)

**Source**: Tuple.java:L15-L17  
Factory method returning `new Tuple3<>(f, s, t)`.

### Tuple2.newTuple(f, s)

**Source**: Tuple.java:L28-L30  
Static factory — identical to `new Tuple2<>(f, s)`. Redundant with `Tuple.tuple()`.

### Tuple3.newTuple(f, s, t)

**Source**: Tuple.java:L43-L45  
Static factory — identical to `new Tuple3<>(f, s, t)`. Redundant with `Tuple.tuple()`.

## Data Flow: VariationMap Through the Pipeline

### Creation and Population (CigarParser Stage)

```
InitialData constructor → creates HashMap<Integer, VariationMap<String, Variation>>
    for both nonInsertionVariants and insertionVariants.

CigarParser/RecordPreprocessor → calls VariationUtils.getVariation(hash, position, descString)
    → lazy-creates VariationMap at position if absent
    → lazy-creates Variation at descString if absent
    → returns Variation for mutation (increment varsCount, fwd, rev, quality, etc.)
```

**Key**: `getVariation()` (in VariationUtils.java:L424-L438) is the primary way new VariationMap entries are created. It follows the same lazy-init pattern as `getSV()`.

### SV Tagging (StructuralVariantsProcessor Stage)

```
StructuralVariantsProcessor → calls getSV(nonInsertionVariants, position)
    → lazy-creates VariationMap + SV at position
    → mutates sv.pairs++, sv.splits++, sv.clusters++
    → sets sv.type = "DEL" / "DUP" / "INV" / etc.
```

All 15 getSV() calls in StructuralVariantsProcessor operate on `nonInsertionVariants` only. No SV is ever attached to `insertionVariants`.

### SV Removal (VariationRealigner Stage)

```
VariationRealigner → calls removeSV(nonInsertionVariants, position)
    → nulls out map.sv, removes "SV" key
    → used when realignment invalidates SV evidence
```

### Consumption (ToVarsBuilder Stage)

```
ToVarsBuilder.toVars() iterates nonInsertionVariants.entrySet()
    → for each position, checks varsAtCurPosition.sv == null
    → if sv != null: position may be outside region bounds (SV can extend beyond)
    → if sv != null and no coverage: logs error with sv.type
    → createVariant() iterates varsAtCurPosition.keySet()
        → "SV" key is intercepted, serialized to "splits-pairs-clusters" string
        → other keys (allele descriptions) create Variant objects
```

### Output (PostProcess Stage)

```
SimplePostProcessModule → reads Vars.sv (string), checks isEmpty()
SomaticPostProcessModule → reads v1.sv, v2.sv strings, passes to SomaticOutputVariant
```

The SV data undergoes two type transformations:
1. `VariationMap.SV` (struct with int fields) — **input side** (CigarParser through ToVarsBuilder)
2. `Vars.sv` (String, formatted as "splits-pairs-clusters") — **output side** (ToVarsBuilder through OutputVariant)

## Usage Site Census

### VariationMap — 80+ references across codebase

**Outer map type** (always `Map<Integer, VariationMap<String, Variation>>`):
- `InitialData.nonInsertionVariants` / `.insertionVariants` — **creation point**
- `VariationData.nonInsertionVariants` / `.insertionVariants` — CigarParser → StructuralVariantsProcessor pass-through
- `CigarParser.nonInsertionVariants` / `.insertionVariants` — local references
- `RecordPreprocessor.nonInsertionVariants` / `.insertionVariants` — local references
- `VariationRealigner.nonInsertionVariants` / `.insertionVariants` — local references
- `StructuralVariantsProcessor.nonInsertionVariants` / `.insertionVariants` — local references
- `ToVarsBuilder.insertionVariants` / `.nonInsertionVariants` — final consumer

### ConcurrentHashSet — 5 references

- `SomaticMode` lines 47, 64 — two `splice` sets
- `AmpliconMode` line 87 — one `splice` set

### DirectThreadExecutor — 16 references

Used in all mode classes and post-process modules whenever inline pipeline execution is needed.

### Tuple — 30+ references

- `VarDictLauncher` — BED file parsing (Tuple3), sample names (Tuple2)
- `AmpliconMode` — position→(index, region) lists (Tuple2)
- `AmpliconOutputVariant` — good/bad variant lists (Tuple2<Variant, String>)

## Cross-Module Dependencies

### VariationMap

**Depends on**: `Variation` (value type)

**Depended on by** (direct users):
- `InitialData` — field type
- `VariationData` — field type
- `CigarParser` — field type, iterates
- `RecordPreprocessor` — field type
- `VariationRealigner` — field type, calls getSV/removeSV
- `StructuralVariantsProcessor` — calls getSV, reads sv fields
- `ToVarsBuilder` — iterates, reads sv, serializes to Vars.sv
- `VariationUtils` — getVariation() creates new VariationMap instances

### ConcurrentHashSet

**Depends on**: `ConcurrentHashMap`

**Depended on by**: `SomaticMode`, `AmpliconMode` (for splice sets)

### DirectThreadExecutor

**Depends on**: `java.util.concurrent.Executor`

**Depended on by**: All mode classes, StructuralVariantsProcessor, VariationRealigner, SomaticPostProcessModule

### Tuple

**Depended on by**: `VarDictLauncher`, `AmpliconMode`, `AmpliconOutputVariant`, `AmpliconPostProcessModule`

## Rust Correspondence

| Java Class | Rust Equivalent | Notes |
|-----------|-----------------|-------|
| `VariationMap<String, Variation>` (as inner map) | `InnerMap<VarDesc, Variant>` = `VecMap<VarDesc, Variant>` | VecMap preserves insertion order like LinkedHashMap |
| `VariationMap.SV` | `StructuralVariantCounts` in `src/variants/variants.rs:100` | Fields: splits, pairs, clusters (all usize). No `type` field — type is determined contextually |
| `VariationMap.sv` field (per-position sidecar) | `sv_counts: HashMap<i64, StructuralVariantCounts>` (separate map) | Decoupled from variant map — stored as parallel structure |
| `getSV()` | Custom logic in `structural_variants_processor.rs` | Lazy-creates entry in sv_counts HashMap |
| `removeSV()` | `.sv_counts.remove(&pos)` | Simpler — no sentinel key to manage |
| `"SV" sentinel Variation in map` | **Not present** — sentinel eliminated | Must ensure createVariant() doesn't iterate a phantom "SV" key |
| `ConcurrentHashSet<String>` | `HashSet<String>` | Single-threaded in Rust — no concurrency needed |
| `DirectThreadExecutor` | **Eliminated** | Rust pipeline is synchronous — no executor abstraction |
| `Tuple.Tuple2<T1,T2>` | Rust tuple `(T1, T2)` | Direct mapping |
| `Tuple.Tuple3<T1,T2,T3>` | Rust tuple `(T1, T2, T3)` | Direct mapping |
| `Vars.sv` (String) | `Vars.sv: String` in `to_vars_builder.rs:285` | Same: empty string default, "splits-pairs-clusters" format |

### Key Architectural Divergence

In Java, the SV data is **embedded** in the VariationMap (`.sv` field) — a single object holds both the variant map and SV counters. In Rust, these are **decoupled** into separate data structures:

- Variant data: `HashMap<i64, InnerMap<VarDesc, Variant>>` (by position → by description)
- SV data: `HashMap<i64, StructuralVariantCounts>` (by position → counters)

This divergence eliminates the need for:
1. The "SV" sentinel Variation key in the inner map
2. The `getSV()` double-init pattern
3. The `removeSV()` sentinel cleanup

But it requires careful synchronization: any position that has SV data in Java's VariationMap must have a corresponding entry in Rust's `sv_counts` map.

## Known Parity Traps

### TRAP 1: LinkedHashMap Insertion Order

**VariationMap extends LinkedHashMap**, which preserves **insertion order**. When ToVarsBuilder iterates `varsAtCurPosition.keySet()`, it gets keys in the order they were put(). However, ToVarsBuilder then does:

```java
List<String> keys = new ArrayList<>(varsAtCurPosition.keySet());
Collections.sort(keys);
```

This **sorts the keys lexicographically** before creating variants. So the LinkedHashMap iteration order is actually overridden by the sort at this specific site. BUT insertion order still matters for:
- Whether `varsAtCurPosition.isEmpty()` first encounters the map (no ordering issue)
- The `.sv` sidecar access (not affected by map ordering)
- **Debug output**: `dumpVariationMap()` in CigarParser iterates with `TreeMap` (sorted by position), then inner by natural order of VariationMap (insertion order, not sorted). This affects debug dump output ordering.

**Rust mapping**: The Rust port uses `VecMap` (Vec-based linear scan) as `InnerMap<K,V>`, which naturally preserves insertion order, matching LinkedHashMap. The `sv` sidecar is stored in a **separate** `HashMap<i64, StructuralVariantCounts>` (`sv_counts`).

### TRAP 2: getSV() Creates Both SV AND "SV" Sentinel Variation

When `getSV()` is called, it not only creates the SV sidecar but also inserts a **"SV" key with an empty Variation()** into the LinkedHashMap. This sentinel is visible to all code that iterates the map. In ToVarsBuilder.createVariant(), the "SV" key is explicitly skipped (`if (descriptionString.equals("SV")) { continue; }`). If the Rust port doesn't place an equivalent sentinel, the iteration behavior differs.

**Current Rust approach**: The Rust port stores SV counts in a separate `sv_counts` HashMap, so there is no "SV" sentinel key in the variant map. This is correct as long as the "SV" skip logic in createVariant is also removed (since there's nothing to skip).

### TRAP 3: getSV() Double-Guard on "SV" Key

`getSV()` has a subtle defensive pattern:
```java
if (sv == null) {
    sv = new VariationMap.SV();
    map.sv = sv;
    map.put("SV", new Variation());  // inside sv==null block
}
if (!map.containsKey("SV")) {       // OUTSIDE sv==null block
    map.put("SV", new Variation());
}
```

The outer `containsKey("SV")` check handles the case where sv was already set but someone removed the "SV" key from the map. This is a defensive guard — in normal flow, if sv is non-null, "SV" should already be in the map. But if `removeSV()` was called and then `getSV()` is called again on the same position, `removeSV()` nulls sv AND removes "SV", so the first block handles it. The second guard catches any other edge case.

### TRAP 4: removeSV() Has No Null Check

```java
hash.get(start).sv = null;  // NPE if hash.get(start) is null!
```

`removeSV()` assumes the VariationMap at `start` exists. In Rust, attempting to get from a HashMap and unwrap would panic. The Rust port should either:
- Match the Java behavior (panic on missing key — equivalent to Java NPE)
- Add a defensive check (if the key is always expected to exist from a prior getSV() call)

### TRAP 5: SV Fields Have Two Different Representations

- **Input side**: `VariationMap.SV` — struct with `String type`, `int pairs`, `int splits`, `int clusters`
- **Output side**: `Vars.sv` — String formatted as `"splits-pairs-clusters"` (note: splits FIRST, then pairs, then clusters)

The serialization in ToVarsBuilder line 443:
```java
getOrPutVars(alignedVars, position).sv = sv.splits + "-" + sv.pairs + "-" + sv.clusters;
```

**Parity trap**: The ORDER is `splits-pairs-clusters`, not `pairs-splits-clusters`. Getting the field order wrong changes the output string.

### TRAP 6: Vars.sv Default Is Empty String, Not Null

```java
public String sv = "";  // in Vars.java
```

`Vars.sv` defaults to `""` (empty string), not null. PostProcess modules check `variantsOnPosition.sv.isEmpty()` — this would NPE if sv were null. The empty-string default is important. In Rust, `String::new()` or `String::default()` correctly produces an empty string.

### TRAP 7: ConcurrentHashSet Iteration Order Is Non-Deterministic

ConcurrentHashSet delegates to `ConcurrentHashMap.keySet().iterator()`, which has **no guaranteed order**. The `splice` set is only used for membership testing (`contains()`), never iterated for output, so this is not a parity issue. But if any code path were to iterate it for output purposes, the Rust `HashSet` would also be non-deterministic, matching the Java behavior.

### TRAP 8: getVariation() Lazy-Creates VariationMap Without SV

`VariationUtils.getVariation()` (the primary way variants are added to the map) creates `new VariationMap<>()` with `sv = null`. It does NOT create any SV data. Only `getSV()` creates SV data. This means most VariationMap entries have `sv == null`, and only positions with structural variant evidence have `sv != null`.

### TRAP 9: Tuple Fields Are Immutable (final) But Can Hold Mutable Objects

`Tuple2._1` and `Tuple2._2` are `final`, but if T1/T2 is a mutable type (like `Region`), the referenced object can still be mutated. In Rust, the equivalent would be `(T1, T2)` tuple or a struct with `pub` fields — by default immutable unless wrapped in a mutable reference.

### TRAP 10: VariationMap Is NEVER Used for insertionVariants.sv

Despite `insertionVariants` being typed as `Map<Integer, VariationMap<String, Variation>>`, the `sv` field is **NEVER accessed** on insertion variant maps. All `getSV()` and `removeSV()` calls pass `nonInsertionVariants`. All `.sv` field accesses in ToVarsBuilder operate on `nonInsertionVariants` entries. The insertion map just uses VariationMap as a convenient LinkedHashMap subclass.

### TRAP 11: VariationMap Has No Custom serialVersionUID

VariationMap extends LinkedHashMap which is Serializable. No `serialVersionUID` is declared. This is irrelevant for parity (no serialization in the pipeline) but worth noting for completeness.

### TRAP 12: ConcurrentHashSet Does Not Override equals()

The `equals()` method is inherited from Object (identity equality). This means two ConcurrentHashSets with identical contents are NOT equal. For VarDict's usage this is fine — ConcurrentHashSet is only used as a local accumulator within mode processing, never compared for equality.
