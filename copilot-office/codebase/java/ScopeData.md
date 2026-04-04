# ScopeData

**Source**: `data/scopedata/*.java` + `data/*.java`  
**LOC**: ~950 (7 scopedata files + 13 data files)  
**Rust counterpart**: `src/scopedata/`, `src/data/`  
**Status**: complete  
**Last updated**: 2026-04-04  
**Classes covered**: 20 (7 scopedata + 13 data)

## Overview

The `data/` and `data/scopedata/` packages form the **data layer** of VarDictJava. They contain:

1. **Pipeline state containers** (`Scope<T>`, `GlobalReadOnlyScope`) — carry typed data between pipeline stages
2. **Stage-specific data** (`InitialData`, `VariationData`, `RealignedVariationData`, `AlignedVarsData`, `CombineAnalysisData`) — each holds the outputs of one stage for the next
3. **Reference access** (`Reference`, `ReferenceResource`) — FASTA fetching, sequence caching, seed map construction
4. **Domain value objects** (`Region`, `BaseInsertion`, `Match`, `Match35`, `ModifiedCigar`, `CurrentSegment`, `SVStructures`, `Side`, `SortPositionSclip`)
5. **Compiled regex constants** (`Patterns`) — centralized patterns used by all modules
6. **BAM I/O** (`SamView`) — thread-local BAM reader wrapper

### Pipeline Stage Flow with Scope<T>

```
SAMFileParser.parseSAM()
    → creates Scope<InitialData>
        ↓
CigarParser.cigarParser()
    → consumes Scope<InitialData>, produces Scope<VariationData>
        ↓
VariationRealigner + StructuralVariantsProcessor
    → consumes Scope<VariationData>, produces Scope<RealignedVariationData>
        ↓
ToVarsBuilder.toVars()
    → consumes Scope<RealignedVariationData>, produces Scope<AlignedVarsData>
        ↓
Modes (Simple/Somatic/Amplicon)
    → consumes Scope<AlignedVarsData>
    → (Somatic only) produces Scope<CombineAnalysisData>
        ↓
OutputVariant printers
```

## Class Inventory

| # | Class | Package | LOC | Fields | Constructors | Methods | Status |
|---|-------|---------|-----|--------|--------------|---------|--------|
| 1 | `Scope<T>` | scopedata | 38 | 8 | 2 | 0 | ✅ Analyzed |
| 2 | `GlobalReadOnlyScope` | scopedata | 73 | 8+2 static | 1 | 4 static | ✅ Analyzed |
| 3 | `InitialData` | scopedata | 41 | 5 | 2 | 0 | ✅ Analyzed |
| 4 | `VariationData` | scopedata | 29 | 12 | 0 | 0 | ✅ Analyzed |
| 5 | `RealignedVariationData` | scopedata | 60 | 11 (all final) | 1 | 0 | ✅ Analyzed |
| 6 | `AlignedVarsData` | scopedata | 16 | 2 | 1 | 0 | ✅ Analyzed |
| 7 | `CombineAnalysisData` | scopedata | 11 | 2 | 1 | 0 | ✅ Analyzed |
| 8 | `BaseInsertion` | data | 27 | 3 | 1 | 0 | ✅ Analyzed |
| 9 | `CurrentSegment` | data | 15 | 3 | 1 | 0 | ✅ Analyzed |
| 10 | `Match` | data | 12 | 2 | 1 | 0 | ✅ Analyzed |
| 11 | `Match35` | data | 19 | 3 | 1 | 0 | ✅ Analyzed |
| 12 | `ModifiedCigar` | data | 12 | 4 | 1 | 0 | ✅ Analyzed |
| 13 | `Patterns` | data | 183 | 67 static | 0 | 0 | ✅ Analyzed |
| 14 | `Reference` | data | 81 | 3 | 1 | 0 | ✅ Analyzed |
| 15 | `ReferenceResource` | data | 184 | 1 (ThreadLocal) | 0 | 4 | ✅ Analyzed |
| 16 | `Region` | data | 102 | 6 (all final) | 2 | 4 | ✅ Analyzed |
| 17 | `SVStructures` | data | 56 | 18 | 0 | 0 | ✅ Analyzed |
| 18 | `SamView` | data | 42 | 2+1 static | 1 | 2 | ✅ Analyzed |
| 19 | `Side` | data | 9 | 3 enum | 0 | 1 | ✅ Analyzed |
| 20 | `SortPositionSclip` | data | 19 | 3 | 1 | 0 | ✅ Analyzed |

## Class Analyses

---

### Scope<T>

**Source**: `data/scopedata/Scope.java:L1-L42`  
**Purpose**: Generic pipeline state container carrying common context plus stage-specific typed data.

#### Fields (all `public final`)

| Field | Type | Description |
|-------|------|-------------|
| `bam` | `String` | BAM file path |
| `region` | `Region` | Target genomic region |
| `regionRef` | `Reference` | Reference sequence + seed map for region |
| `referenceResource` | `ReferenceResource` | FASTA accessor for additional reference fetches |
| `maxReadLength` | `int` | Max read length observed during parsing |
| `splice` | `Set<String>` | Splice site identifiers |
| `out` | `VariantPrinter` | Output printer (stdout, file, etc.) |
| `data` | `T` | Stage-specific typed payload |

#### Constructors

1. **Full constructor** (8 params): Sets all fields directly.
2. **Inherit constructor** `Scope(Scope<?> inheritableScope, T data)`: Copies all common fields from a previous scope and replaces only the `data` payload. This is the primary mechanism for pipeline stage transitions.

#### Parity Notes
- `Scope<?>` wildcard usage: Java type erasure means at runtime the generic is erased. Rust must use an enum or separate struct types.
- The inherit constructor preserves the **same** object references for `regionRef`, `referenceResource`, etc. — these are shared between stages, not cloned.

---

### GlobalReadOnlyScope

**Source**: `data/scopedata/GlobalReadOnlyScope.java:L1-L73`  
**Purpose**: Singleton holding global configuration and sample metadata. Accessed everywhere via `instance()`.

#### Static State

| Field | Type | Mutability |
|-------|------|------------|
| `instance` | `volatile GlobalReadOnlyScope` | Write-once via `init()` |
| `mode` | `volatile AbstractMode` | Write-once via `setMode()` |

#### Instance Fields (all `public final`)

| Field | Type | Source |
|-------|------|--------|
| `conf` | `Configuration` | CLI parsed config |
| `chrLengths` | `Map<String, Integer>` | Chromosome name → length from BAM header |
| `sample` | `String` | Sample name (tumor in somatic mode) |
| `samplem` | `String` | Matched sample name (normal in somatic mode) |
| `ampliconBasedCalling` | `String` | Amplicon params string or null |
| `printerTypeOut` | `PrinterType` | Output format (copied from `conf.printerType`) |
| `adaptorForward` | `Map<String, Integer>` | Forward adaptor sequences |
| `adaptorReverse` | `Map<String, Integer>` | Reverse adaptor sequences |

#### Static Methods

| Method | Behavior |
|--------|----------|
| `instance()` | Returns singleton. No null check — callers assume initialized. |
| `init(...)` | Synchronized. Throws `IllegalStateException` if already initialized. Creates instance. |
| `getMode()` | Returns `mode` static field. |
| `setMode(AbstractMode)` | Synchronized. Throws if mode already set. |
| `clear()` | Synchronized. Sets `instance = null; mode = null`. **Test-only.** |

#### Constructor
Takes 7 parameters and assigns all `final` fields. `printerTypeOut` is derived from `conf.printerType`.

#### Parity Notes
- **Singleton → Rust**: The Rust port uses a global static (likely `once_cell::sync::Lazy` or similar). Thread-safety is guaranteed by immutability after init.
- `chrLengths` uses `HashMap<String, Integer>` — ordering does not matter (used only for lookup).
- `adaptorForward` and `adaptorReverse` also `HashMap` — used for lookup-based adaptor trimming.
- `instance()` is called from nearly every module. A `null` return would NPE, but `init()` is always called during startup before any pipeline work.

---

### InitialData

**Source**: `data/scopedata/InitialData.java:L1-L41`  
**Purpose**: Initial empty data containers for variant collection during BAM parsing.

#### Fields (all `public`, mutable)

| Field | Type | Default (no-arg ctor) | Java Map Type |
|-------|------|-----------------------|---------------|
| `nonInsertionVariants` | `Map<Integer, VariationMap<String, Variation>>` | `new HashMap<>()` | **HashMap** outer |
| `insertionVariants` | `Map<Integer, VariationMap<String, Variation>>` | `new HashMap<>()` | **HashMap** outer |
| `refCoverage` | `Map<Integer, Integer>` | `new HashMap<>()` | **HashMap** |
| `softClips5End` | `Map<Integer, Sclip>` | `new HashMap<>()` | **HashMap** |
| `softClips3End` | `Map<Integer, Sclip>` | `new HashMap<>()` | **HashMap** |

Note: `VariationMap<String, Variation>` extends `LinkedHashMap<String, Variation>` — so the **inner** map at each position preserves insertion order, but the **outer** position→map mapping is an unordered HashMap.

#### Constructors

1. **No-arg**: Creates all maps as empty HashMaps.
2. **Full 5-param**: Accepts pre-populated maps (used when CigarParser builds up data, then passes it forward).

#### Parity Notes
- Outer maps: `HashMap<Integer, ...>` → Rust `HashMap<i32, ...>` (order does not affect output since ToVarsBuilder iterates over a sorted key range).
- Inner maps: `VariationMap` = `LinkedHashMap` → Rust `IndexMap` (insertion order preserved). **Critical**: the order of variant description strings at each position affects output row ordering.

---

### VariationData

**Source**: `data/scopedata/VariationData.java:L1-L29`  
**Purpose**: Post-CigarParser data carrying all variation maps, SV structures, and per-shard metadata.

#### Fields (all `public`, mutable, **no constructor**)

| Field | Type | Default | Notes |
|-------|------|---------|-------|
| `nonInsertionVariants` | `Map<Integer, VariationMap<String, Variation>>` | null | Transferred from InitialData |
| `insertionVariants` | `Map<Integer, VariationMap<String, Variation>>` | null | Transferred from InitialData |
| `positionToInsertionCount` | `Map<Integer, Map<String, Integer>>` | null | Position → ins_description → count |
| `positionToDeletionsCount` | `Map<Integer, Map<String, Integer>>` | null | Position → del_description → count |
| `svStructures` | `SVStructures` | null | SV data buckets |
| `refCoverage` | `Map<Integer, Integer>` | null | Position → ref coverage |
| `softClips5End` | `Map<Integer, Sclip>` | null | 5' soft-clip clusters |
| `softClips3End` | `Map<Integer, Sclip>` | null | 3' soft-clip clusters |
| `maxReadLength` | `Integer` | null | Boxed; nullable |
| `splice` | `Set<String>` | null | Splice identifiers |
| `mnp` | `Map<Integer, Map<String, Integer>>` | null | Multi-nucleotide polymorphism counts |
| `spliceCount` | `Map<String, int[]>` | null | Splice → [forward_count, reverse_count] |
| `duprate` | `double` | 0.0 (primitive) | Duplicate rate |

#### Parity Notes
- **No constructor** — all fields are assigned individually by CigarParser after creation via `new VariationData()`. All reference-type fields start as `null`.
- `maxReadLength` is `Integer` (boxed) — can be null. In Rust, `Option<i32>`.
- `duprate` is a primitive `double` — defaults to 0.0, never null. In Rust, `f64`.
- `positionToInsertionCount` and `positionToDeletionsCount`: outer is `HashMap<Integer, ...>`, inner `Map<String, Integer>` — the inner maps are also `HashMap` in Java unless explicitly stated. Check usage: if iteration order matters at consumption, must verify.
- `spliceCount` values are `int[]` of length 2 — `[forward, reverse]`. Rust: `[i32; 2]` or a small struct.

---

### RealignedVariationData

**Source**: `data/scopedata/RealignedVariationData.java:L1-L60`  
**Purpose**: Post-realignment/post-SV data. All fields are `final` — immutable after construction.

#### Fields (all `public final`)

| Field | Type | Notes |
|-------|------|-------|
| `nonInsertionVariants` | `Map<Integer, VariationMap<String, Variation>>` | From VariationData (potentially modified by realigner) |
| `insertionVariants` | `Map<Integer, VariationMap<String, Variation>>` | Same |
| `softClips5End` | `Map<Integer, Sclip>` | May be modified by SV processor |
| `softClips3End` | `Map<Integer, Sclip>` | Same |
| `refCoverage` | `Map<Integer, Integer>` | Position → reference coverage |
| `maxReadLength` | `Integer` | Boxed (nullable) |
| `svStructures` | `SVStructures` | SV analysis results |
| `duprate` | `double` | Duplication rate |
| `CURSEG` | `CurrentSegment` | Current chromosome segment |
| `SOFTP2SV` | `Map<Integer, List<Sclip>>` | Position → list of soft-clips for SV analysis |
| `previousScope` | `Scope<VariationData>` | Reference back to the previous pipeline stage |

#### Constructor
Single constructor taking 11 parameters. Note the parameter ORDER differs from field declaration order (e.g., `softClips3End` before `softClips5End` in params but vice versa in fields).

#### Parity Notes
- `SOFTP2SV` type: `Map<Integer, List<Sclip>>` — a **HashMap** (default). The iteration order of this map may or may not affect output. In StructuralVariantsProcessor the map is iterated to find SV candidates.
- `previousScope` creates a reference cycle: `Scope<RealignedVariationData>` → `RealignedVariationData.previousScope` → `Scope<VariationData>`. In Java this is fine (GC handles it); in Rust may need `Arc` or restructuring.
- `CURSEG` is mutable class `CurrentSegment` stored as final reference — contents can still mutate. Trap for Rust `&` vs `&mut`.

---

### AlignedVarsData

**Source**: `data/scopedata/AlignedVarsData.java:L1-L16`  
**Purpose**: Post-ToVarsBuilder data containing the final variant map.

#### Fields

| Field | Type | Description |
|-------|------|-------------|
| `maxReadLength` | `int` | Primitive (not nullable) |
| `alignedVariants` | `Map<Integer, Vars>` | Position → Vars (ref variant + variant list) |

#### Constructor
Single constructor: `AlignedVarsData(int maxReadLength, Map<Integer, Vars> alignedVariants)`.

#### Parity Notes
- `alignedVariants` map type: the concrete type depends on how ToVarsBuilder creates it. Check ToVarsBuilder — it's typically a `HashMap`, but iteration in Modes may need ordering. This is output-affecting since it determines which positions get printed first.

---

### CombineAnalysisData

**Source**: `data/scopedata/CombineAnalysisData.java:L1-L11`  
**Purpose**: Somatic mode result after combining tumor/normal analysis.

#### Fields

| Field | Type | Description |
|-------|------|-------------|
| `maxReadLength` | `int` | Max read length from both samples |
| `type` | `String` | Variant type string (e.g., "Somatic", "Germline", "LOH", etc.) |

#### Constructor
Single constructor: `CombineAnalysisData(int maxReadLength, String type)`.

---

### BaseInsertion

**Source**: `data/BaseInsertion.java:L1-L27`  
**Purpose**: Holds position and sequence of a base insertion found during CIGAR realignment.

#### Fields

| Field | Type | Java Var | Default | Description |
|-------|------|----------|---------|-------------|
| `baseInsert` | `Integer` (boxed) | `$bi` | ctor param | Starting position of insert |
| `insertionSequence` | `String` | `$ins` | ctor param | Insertion sequence |
| `baseInsert2` | `Integer` (boxed) | `$bi2` | ctor param | Position without extra sequence |

#### Constructor
`BaseInsertion(int baseInsert, String insertionSequence, int baseInsert2)` — note: params are **primitive** `int` but fields are **boxed** `Integer`. Auto-boxing occurs. In practice, these values are never null after construction.

#### Parity Notes
- Fields are `Integer` (boxed), allowing null in theory. In practice always initialized from `int` params. Rust: `i32` should suffice unless the field is later set to null somewhere.

---

### CurrentSegment

**Source**: `data/CurrentSegment.java:L1-L15`  
**Purpose**: Tracks current processing segment (chromosome + start/end).

#### Fields

| Field | Type | Description |
|-------|------|-------------|
| `chr` | `String` | Chromosome name |
| `start` | `int` | Segment start position |
| `end` | `int` | Segment end position |

#### Constructor
`CurrentSegment(String chr, int start, int end)`.

#### Parity Notes
- Fields are mutable (`public`, no `final`). The StructuralVariantsProcessor modifies `CURSEG.start` and `CURSEG.end` in place. Rust needs `&mut` access or interior mutability.

---

### Match

**Source**: `data/Match.java:L1-L12`  
**Purpose**: Result of sequence-to-reference matching.

#### Fields

| Field | Type | Java Var | Description |
|-------|------|----------|-------------|
| `basePosition` | `int` | `$bp` | Matched base position in reference |
| `matchedSequence` | `String` | — | Matched sequence string |

#### Constructor
`Match(int basePosition, String matchedSequence)`.

---

### Match35

**Source**: `data/Match35.java:L1-L19`  
**Purpose**: Match positions at both ends (5' and 3') of a sequence.

#### Fields

| Field | Type | Description |
|-------|------|-------------|
| `matched5end` | `int` | Start matched position of 5' end |
| `matched3End` | `int` | Start matched position of 3' end |
| `maxMatchedLength` | `int` | Maximum matched length |

#### Constructor
`Match35(int matched5end, int matched3End, int maxMatchedLength)`.

---

### ModifiedCigar

**Source**: `data/ModifiedCigar.java:L1-L12`  
**Purpose**: Holds modified CIGAR string and associated query data after CigarModifier processing.

#### Fields

| Field | Type | Description |
|-------|------|-------------|
| `position` | `int` | Adjusted alignment position |
| `cigar` | `String` | Modified CIGAR string |
| `querySequence` | `String` | Possibly trimmed query sequence |
| `queryQuality` | `String` | Possibly trimmed quality string |

#### Constructor
`ModifiedCigar(int position, String cigar, String querySequence, String queryQuality)`.

---

### Patterns (Deep Dive)

**Source**: `data/Patterns.java:L1-L183`  
**Purpose**: Centralized compiled regex pattern constants used across all VarDict modules.

#### Regex Engine Split

VarDictJava uses **two** regex engines:
1. `java.util.regex.Pattern` — standard Java regex (67 patterns)
2. `jregex.Pattern` — third-party library, used for patterns that need group replacement or `globalFind()` (17 patterns)

The `jregex` library provides `jregex.Replacer` and `jregex.Matcher.find()` which differ from `java.util.regex` in some edge cases. In Rust, all patterns will use the `regex` crate.

#### Complete Pattern Inventory

##### SAMRecord Patterns
| # | Name | Engine | Regex | Used In |
|---|------|--------|-------|---------|
| 1 | `MC_Z_NUM_S_ANY_NUM_S` | jregex | `\d+S\S*\d+S` | RecordPreprocessor — MC:Z tag filtering |

##### Variation Patterns
| # | Name | Engine | Regex | Used In |
|---|------|--------|-------|---------|
| 2 | `BEGIN_DIGITS` | java | `^(\d+)` | Variation, CigarParser — extract leading number |
| 3 | `UP_NUMBER_END` | java | `\^(\d+)$` | VariationRealigner — detect `^N` suffix |
| 4 | `BEGIN_MINUS_NUMBER_ANY` | java | `^-\d+(.*)` | Variation — deletion with trailing content |
| 5 | `BEGIN_MINUS_NUMBER_CARET` | java | `^-\d+\^` | VariationRealigner — deletion with caret |
| 6 | `BEGIN_MINUS_NUMBER` | java | `^-(\d+)` | ToVarsBuilder, output — extract deletion length |
| 7 | `MINUS_NUM_NUM` | jregex | `-\d\d` | SomaticPostProcessModule — 2-digit deletion |
| 8 | `HASH_GROUP_CARET_GROUP` | java | `#(.+)\^(.+)` | VariationRealigner — complex variant `#seq^seq` |

##### Sclip Patterns
| # | Name | Engine | Regex | Used In |
|---|------|--------|-------|---------|
| 9 | `B_A7` | java | `^.AAAAAAA` | Sclip filtering — 7-base polyA |
| 10 | `B_T7` | java | `^.TTTTTTT` | Sclip filtering — 7-base polyT |

##### ATGC (Sequence) Patterns
| # | Name | Engine | Regex | Used In |
|---|------|--------|-------|---------|
| 11 | `CARET_ATGNC` | java | `\^([ATGNC]+)` | Multiple — caret-prefixed bases |
| 12 | `CARET_ATGC_END` | java | `\^([ATGC]+)$` | Realigner — terminal caret+bases |
| 13 | `AMP_ATGC` | java | `&([ATGC]+)` | Realigner — ampersand-delimited bases |
| 14 | `BEGIN_PLUS_ATGC` | java | `^\+([ATGC]+)` | ToVarsBuilder — insertion description |
| 15 | `HASH_ATGC` | java | `#([ATGC]+)` | ToVarsBuilder — complex variant hash |
| 16 | `ATGSs_AMP_ATGSs_END` | java | `(\+[ATGC]+)&[ATGC]+$` | Realigner — ins + amp trailing |
| 17 | `MINUS_NUMBER_AMP_ATGCs_END` | java | `(-\d+)&[ATGC]+$` | Realigner — del + amp trailing |
| 18 | `MINUS_NUMBER_ATGNC_SV_ATGNC_END` | java | `^-\d+\^([ATGNC]+)<...\d+>([ATGNC]+)$` | SV — deletion with SV notation |
| 19 | `BEGIN_ATGC_END` | java | `^[ATGC]+$` | Multiple — pure base string check |

##### SV Patterns
| # | Name | Engine | Regex | Used In |
|---|------|--------|-------|---------|
| 20 | `DUP_NUM` | java | `<dup(\d+)` | SV processor — duplication pos |
| 21 | `DUP_NUM_ATGC` | java | `<dup(\d+)>([ATGC]+)$` | SV processor — dup + bases |
| 22 | `INV_NUM` | java | `<inv(\d+)` | SV processor — inversion pos |
| 23 | `SOME_SV_NUMBERS` | java | `<(...)\d+>` | Output — SV type+number |
| 24 | `ANY_SV` | java | `<(...)>` | Output — SV type only |

##### File and Column Patterns
| # | Name | Engine | Regex | Used In |
|---|------|--------|-------|---------|
| 25 | `SAMPLE_PATTERN` | java | `([^\/\._]+).sorted[^\/]*.bam` | Main — extract sample from BAM path |
| 26 | `SAMPLE_PATTERN2` | java | `([^\/]+)[_\.][^\/]*bam` | Main — fallback sample extraction |
| 27 | `INTEGER_ONLY` | java | `^\d+$` | Config — numeric column check |

##### CIGAR Patterns (jregex)
| # | Name | Engine | Regex | Used In |
|---|------|--------|-------|---------|
| 28 | `BEGIN_NUMBER_S_NUMBER_IorD` | jregex | `^(\d+)S(\d+)([ID])` | CigarModifier — leading S+I/D |
| 29 | `NUMBER_IorD_NUMBER_S_END` | jregex | `(\d+)([ID])(\d+)S$` | CigarModifier — trailing I/D+S |
| 30 | `BEGIN_NUMBER_S_NUMBER_M_NUMBER_IorD` | jregex | `^(\d+)S(\d+)M(\d+)([ID])` | CigarModifier — S+M+I/D prefix |
| 31 | `NUMBER_IorD_NUMBER_M_NUMBER_S_END` | jregex | `(\d+)([ID])(\d+)M(\d+)S$` | CigarModifier — I/D+M+S suffix |
| 32 | `BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M` | jregex | `^(\d)M(\d+)([ID])(\d+)M` | CigarModifier — single digit M |
| 33 | `BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M_` | java | `^\dM\d+[ID]\d+M` | CigarModifier — non-capturing variant |
| 34 | `NUMBER_IorD_DIGIT_M_END` | jregex | `(\d+)([ID])(\d)M$` | CigarModifier — trailing single-digit M |
| 35 | `NUMBER_IorD_NUMBER_M_END` | jregex | `(\d+)([ID])(\d+)M$` | CigarModifier — trailing I/D+M |
| 36 | `D_M_D_DD_M_D_I_D_M_D_DD` | jregex | `^(.*?)(\d+)M(\d+)D(\d+)M(\d+)I(\d+)M(\d+)D(\d+)M` | CigarModifier — complex M-D-M-I-M-D-M |
| 37 | `D_M_D_DD_M_D_I_D_M_D_DD_prim` | java | `(\d+)M(\d+)D(\d+)M(\d+)I(\d+)M(\d+)D(\d+)M` | CigarModifier — non-anchored variant |
| 38 | `threeDeletionsPattern` | jregex | `^(.*?)(\d+)M(\d+)D(\d+)M(\d+)D(\d+)M(\d+)D(\d+)M` | CigarModifier — 3 deletions |
| 39 | `threeIndelsPattern` | jregex | `^(.*?)(\d+)M(\d+)([DI])(\d+)M(\d+)([DI])(\d+)M(\d+)([DI])(\d+)M` | CigarModifier — 3 indels |
| 40 | `DIGM_D_DI_DIGM_D_DI_DIGM_DI_DIGM` | java | `\d+M\d+[DI]\d+M\d+[DI]\d+M\d+[DI]\d+M` | CigarModifier — 3-indel detect |
| 41 | `DM_DD_DM_DD_DM_DD_DM` | java | `\d+M\d+D\d+M\d+D\d+M\d+D\d+M` | CigarModifier — 3-del detect |
| 42 | `DIG_D_DIG_M_DIG_DI_DIGI` | java | `(\d+)D(\d+)M(\d+)([DI])(\d+I)?` | CigarModifier — D+M+D/I |
| 43 | `DIG_I_DIG_M_DIG_DI_DIGI` | java | `(\d+)I(\d+)M(\d+)([DI])(\d+I)?` | CigarModifier — I+M+D/I |
| 44 | `NOTDIG_DIG_I_DIG_M_DIG_DI_DIGI` | java | `(\D)(\d+)I(\d+)M(\d+)([DI])(\d+I)?` | CigarModifier — non-digit prefix |
| 45 | `DIG_D_DIG_D` | java | `(\d+)D(\d+)D` | CigarModifier — merge adjacent deletions |
| 46 | `DIG_I_DIG_I` | java | `(\d+)I(\d+)I` | CigarModifier — merge adjacent insertions |
| 47 | `BEGIN_ANY_DIG_M_END` | java | `^(.*?)(\d+)M$` | CigarModifier — trailing M |
| 48 | `DIG_M_END` | java | `\d+M$` | Multiple — has trailing M |
| 49 | `BEGIN_DIG_M` | java | `^(\d+)M` | Multiple — leading M count |
| 50 | `DIG_S_DIG_M` | java | `^(\d+)S(\d+)M` | CigarModifier — leading S+M |
| 51 | `DIG_M_DIG_S_END` | java | `\d+M\d+S$` | CigarModifier — trailing M+S |
| 52 | `ANY_NUMBER_M_NUMBER_S_END` | java | `^(.*?)(\d+)M(\d+)S$` | CigarModifier — capture M+S suffix |
| 53 | `BEGIN_NUMBER_D` | jregex | `^(\d+)D` | CigarParser — leading deletion |
| 54 | `END_NUMBER_D` | jregex | `(\d+)D$` | CigarParser — trailing deletion |
| 55 | `BEGIN_NUMBER_I` | jregex | `^(\d+)I` | CigarParser — leading insertion |
| 56 | `END_NUMBER_I` | jregex | `(\d+)I$` | CigarParser — trailing insertion |
| 57 | `ALIGNED_LENGTH_MND` | jregex | `(\d+)[MND]` | Utils — aligned length (M+N+D) |
| 58 | `ALIGNED_LENGTH_MD` | jregex | `(\d+)[MD=X]` | Utils — aligned length (M+D+=+X) |
| 59 | `SOFT_CLIPPED` | jregex | `(\d+)[MIS]` | RecordPreprocessor — soft-clipped length |
| 60 | `SA_CIGAR_D_S_5clip` | java | `^\d\d+S` | SAMFileParser — SA 5' clip detect (≥2 digits) |
| 61 | `SA_CIGAR_D_S_5clip_GROUP` | java | `^(\d\d+)S` | SAMFileParser — SA 5' clip capture |
| 62 | `SA_CIGAR_D_S_5clip_GROUP_Repl` | jregex | `^\d+S` | SAMFileParser — SA 5' clip replace |
| 63 | `SA_CIGAR_D_S_3clip` | java | `\d\dS$` | SAMFileParser — SA 3' clip detect |
| 64 | `SA_CIGAR_D_S_3clip_GROUP` | java | `(\d\d+)S$` | SAMFileParser — SA 3' clip capture |
| 65 | `SA_CIGAR_D_S_3clip_GROUP_Repl` | jregex | `\d\d+S$` | SAMFileParser — SA 3' clip replace |
| 66 | `BEGIN_dig_dig_S_ANY_dig_dig_S_END` | java | `^\d\dS.*\d\dS$` | SAMFileParser — both clips ≥10bp |
| 67 | `BEGIN_NUM_S_OR_BEGIN_NUM_H` | jregex | `^(\d+)S\|^\d+H` | Utils, CigarParser — leading S or H |
| 68 | `END_NUM_S_OR_NUM_H` | jregex | `(\d+)S$\|H$` | Utils, CigarParser — trailing S or H |

##### Exception Patterns
| # | Name | Engine | Regex | Used In |
|---|------|--------|-------|---------|
| 69 | `UNABLE_FIND_CONTIG` | java | `Unable to find entry for contig` | ReferenceResource — FASTA error |
| 70 | `WRONG_START_OR_END` | java | `Malformed query` | ReferenceResource — FASTA error |

**Total: 70 patterns** (53 java.util.regex, 17 jregex)

#### Parity-Critical Notes for Patterns
1. `jregex` patterns use `jregex.Replacer` for in-place group replacement — Rust `regex` crate `Regex::replace()` semantics differ (especially for group numbering).
2. `jregex.Pattern` alternation `|` may have different precedence than `java.util.regex` for anchored alternation like `^(\d+)S|^\d+H`. In Java `java.util.regex`, this matches either `^(\d+)S` or `^\d+H`. The `jregex` library handles this the same way, but verify in Rust.
3. Pattern `\d\d+S` matches 2-or-more digits followed by S — Rust regex handles this fine but note the minimum quantifier: `\d{2,}S` would be equivalent and clearer.
4. `BEGIN_ANY_DIG_M_END` uses lazy `(.*?)` — Rust regex crate handles lazy quantifiers correctly.
5. All `java.util.regex.Pattern` objects are compiled once (static final) — Rust should use `lazy_static!` or `once_cell::Lazy` with `regex::Regex`.
6. `SOME_SV_NUMBERS` uses `(...)` which means exactly 3 characters — this captures the SV type abbreviation (del, dup, inv, ins).

---

### Reference

**Source**: `data/Reference.java:L1-L81`  
**Purpose**: Container for reference sequence data and seed map for a genomic region.

#### Fields (all `public`)

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `loadedRegions` | `List<LoadedRegion>` | `new ArrayList<>()` | Tracks which regions have been loaded to avoid duplicate fetches |
| `referenceSequences` | `Map<Integer, Character>` | `new HashMap<>()` | Position → reference base. **HashMap** — unordered |
| `seed` | `Map<String, List<Integer>>` | `new HashMap<>()` | K-mer → list of positions. **HashMap** — unordered |

#### Inner Class: LoadedRegion

| Field | Type | Description |
|-------|------|-------------|
| `chr` | `String` | Chromosome |
| `sequenceStart` | `int` | Start position of loaded region |
| `sequenceEnd` | `int` | End position of loaded region |

Implements `equals()` and `hashCode()` (based on all 3 fields).

#### Constructor
No-arg constructor creates empty collections.

#### Seed Map Algorithm (detailed in ReferenceResource)

The seed map indexes the reference by short k-mer subsequences (of lengths `SEED_1=17` and `SEED_2=12`). For each position in the loaded reference, two keys are added to the seed map:
1. The 17-base substring starting at that position
2. The 12-base substring starting at that position

Each key maps to a list of all positions where that k-mer occurs.

#### Parity Notes
- `referenceSequences`: HashMap<Integer, Character>. Key = genomic position (1-based), Value = uppercase base character. In Rust, this is typically a `HashMap<i32, u8>` (since bases are ASCII). **ORDER does NOT matter** — used for point lookups only.
- `seed`: HashMap<String, List<Integer>>. Order of positions within each list IS significant — they're added in left-to-right order and searched sequentially by `findMatch()` in VariationRealigner and StructuralVariantsProcessor. But the map itself (key order) doesn't matter — lookups are by exact k-mer string.
- `loadedRegions` is an `ArrayList` — ordered by insertion order (chronological). Used only for containment checks.

---

### ReferenceResource

**Source**: `data/ReferenceResource.java:L1-L184`  
**Purpose**: FASTA accessor and reference builder. Manages thread-local FASTA file handles, fetches subsequences, and populates `Reference` objects with sequence data and seed maps.

#### Fields

| Field | Type | Description |
|-------|------|-------------|
| `threadLocalFastaFiles` | `ThreadLocal<Map<String, IndexedFastaSequenceFile>>` | Thread-local cache of open FASTA file handles |

#### Methods

##### `fetchFasta(String file)` — L30-L39
- **Synchronized** method
- Lazily opens `IndexedFastaSequenceFile` per FASTA path, caches in ThreadLocal map
- Throws `IllegalArgumentException` if file not found

##### `retrieveSubSeq(String fasta, String chr, int start, int end)` — L48-L63
- Gets FASTA handle via `fetchFasta()`
- Calls `idx.getSubsequenceAt(chr, start, end)` — htsjdk method, 1-based inclusive coordinates
- Returns `String[]` = `[header, bases]`
- **Error handling**:
  - `SAMException` with "Unable to find entry for contig" → `WrongFastaOrBamException`
  - `SAMException` with "Malformed query" → `RegionBoundariesException`
  - Other `SAMException` → rethrown

##### `getReference(Region region)` — L70-L74
- Convenience overload: calls `getReference(region, conf.referenceExtension, new Reference())`

##### `getReference(Region region, int extension, Reference ref)` — L83-L139
**This is the core algorithm. Step by step:**

1. **Calculate sequence boundaries with padding**:
   ```
   sequenceStart = max(1, region.start - conf.numberNucleotideToExtend - extension)
   len = chrLengths[region.chr] or 0
   sequenceEnd = min(len, region.end + conf.numberNucleotideToExtend + extension)
   ```
   Default: `numberNucleotideToExtend` = varies, `referenceExtension` = 1200.

2. **Fetch FASTA subsequence**: `retrieveSubSeq()` → gets uppercase exon string.

3. **Check if already loaded**: `isLoaded(chr, sequenceStart, sequenceEnd, ref)` — if yes, return ref unchanged.

4. **Record loaded region**: Add `LoadedRegion(chr, sequenceStart, sequenceEnd)` to `ref.loadedRegions`.

5. **Determine iteration end**:
   ```
   siteEnd = (len == sequenceEnd) ? exon.length() : exon.length() - SEED_1
   ```
   If at chromosome end, iterate to the very end (don't subtract `SEED_1`). Otherwise stop `SEED_1` (17) bases before the end to ensure seed substrings have enough room.

6. **Populate referenceSequences and seed map** (loop `i = 0..siteEnd`):
   - **Skip if already loaded**: `if ref.referenceSequences.containsKey(i + sequenceStart) continue`
   - **Add base**: `ref.referenceSequences.put(i + sequenceStart, exon.charAt(i))`
   - **Skip seed at chromosome end**: `if len == sequenceEnd && i > exon.length() - SEED_1 continue`
   - **Add SEED_1 key**: `substr(exon, i, SEED_1)` → add `(key, i + sequenceStart)` to seed map
   - **Add SEED_2 key**: `substr(exon, i, SEED_2)` → add `(key, i + sequenceStart)` to seed map

7. **Return** populated ref.

##### `addPositionsToSeedSequence(...)` — L149-L154
- Helper: `ref.seed.getOrDefault(key, new ArrayList<>())` → add position → `ref.seed.put(key, list)`
- Note: creates a new `ArrayList` if key not present. Not using `computeIfAbsent` — slight inefficiency but functionally identical.

##### `isLoaded(chr, start, end, ref)` — L163-L175
- Iterates all `loadedRegions` in ref
- Returns `true` if any loaded region fully contains the requested range
- Linear scan — O(n) where n = number of loaded regions (typically small)

#### Parity Notes
- **Thread-local FASTA handles**: In Rust, this can be per-thread state via `thread_local!` or per-worker state in rayon.
- **1-based inclusive coordinates**: htsjdk `getSubsequenceAt(chr, start, end)` uses 1-based inclusive. Rust's noodles or rust-htslib may differ.
- **Padding calculation**: `conf.numberNucleotideToExtend` and `conf.referenceExtension` are both subtracted/added; verify Rust matches the exact arithmetic.
- **`substr(exon, i, SEED_1)`**: This is `Utils.substr()` which does `string.substring(i, min(i+length, string.length()))` — safe substring that clips to string end. If the substring is shorter than SEED_1/SEED_2, a shorter key is still inserted into the seed map. This can cause false matches in `findMatch()`.
- **`siteEnd` boundary**: When `len == sequenceEnd` (at chromosome end), the loop goes to `exon.length()` but seed generation is skipped for positions where `i > exon.length() - SEED_1`. This means bases near the chromosome end get added to `referenceSequences` but NOT to the seed map. Parity trap: the seed skip condition uses `>` not `>=`.
- **`isLoaded` containment**: Uses `>=` and `<=` for boundaries — fully inclusive check.

---

### Region

**Source**: `data/Region.java:L1-L102`  
**Purpose**: Immutable representation of a BED file region.

#### Fields (all `public final`)

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `chr` | `String` | — | Chromosome name |
| `start` | `int` | — | Start position |
| `end` | `int` | — | End position |
| `gene` | `String` | — | Gene name |
| `insertStart` | `int` | 0 | Amplicon insert start (0 if not amplicon) |
| `insertEnd` | `int` | 0 | Amplicon insert end (0 if not amplicon) |

#### Constructors

1. **4-param** `(chr, start, end, gene)` — delegates to 6-param with `insertStart=0, insertEnd=0`
2. **6-param** `(chr, start, end, gene, insertStart, insertEnd)` — full constructor

#### Static Method
- `newModifiedRegion(Region region, int changedStart, int changedEnd)` — creates a new Region with different start/end but same chr, gene, insertStart, insertEnd. Used by StructuralVariantsProcessor and VariationRealigner for extended regions.

#### Instance Methods
- `equals(Object)` — compares all 6 fields
- `hashCode()` — hash of all 6 fields
- `toString()` — `"Region [chr=X, start=N, end=N, gene=G, istart=N, iend=N]"`
- `printRegion()` — `"chr:start-end"`

#### Coordinate System
- **BED format is 0-based half-open** `[start, end)`, but VarDictJava's `Region.start` and `Region.end` are **1-based inclusive** after BED parsing (the BED parser adds 1 to start). This is a critical parity point.
- `insertStart` and `insertEnd` are also 1-based inclusive when set.
- The value 0 for `insertStart`/`insertEnd` indicates "not set" (non-amplicon mode).

---

### SVStructures

**Source**: `data/SVStructures.java:L1-L56`  
**Purpose**: Container for structural variant candidate collections.

#### Fields

**Deletion SV fields:**
| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `svdelfend` | `int` | 0 | Forward deletion end position |
| `svdelrend` | `int` | 0 | Reverse deletion end position |
| `svfdel` | `List<Sclip>` | `new ArrayList<>()` | Forward deletion soft-clip candidates |
| `svrdel` | `List<Sclip>` | `new ArrayList<>()` | Reverse deletion soft-clip candidates |

**Duplication SV fields:**
| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `svdupfend` | `int` | 0 | Forward duplication end position |
| `svduprend` | `int` | 0 | Reverse duplication end position |
| `svfdup` | `List<Sclip>` | `new ArrayList<>()` | Forward duplication candidates |
| `svrdup` | `List<Sclip>` | `new ArrayList<>()` | Reverse duplication candidates |

**Inversion SV fields:**
| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `svinvfend3` | `int` | 0 | Forward inversion 3' end position |
| `svinvrend3` | `int` | 0 | Reverse inversion 3' end position |
| `svinvfend5` | `int` | 0 | Forward inversion 5' end position |
| `svinvrend5` | `int` | 0 | Reverse inversion 5' end position |
| `svfinv3` | `List<Sclip>` | `new ArrayList<>()` | Forward 3' inversion candidates |
| `svfinv5` | `List<Sclip>` | `new ArrayList<>()` | Forward 5' inversion candidates |
| `svrinv3` | `List<Sclip>` | `new ArrayList<>()` | Reverse 3' inversion candidates |
| `svrinv5` | `List<Sclip>` | `new ArrayList<>()` | Reverse 5' inversion candidates |

**Fusion SV fields:**
| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `svffus` | `Map<String, List<Sclip>>` | `new HashMap<>()` | Forward fusion — keyed by mate chromosome |
| `svrfus` | `Map<String, List<Sclip>>` | `new HashMap<>()` | Reverse fusion |
| `svfusfend` | `Map<String, Integer>` | `new HashMap<>()` | Forward fusion end positions |
| `svfusrend` | `Map<String, Integer>` | `new HashMap<>()` | Reverse fusion end positions |

**Total: 18 fields** (10 ints, 8 lists, 4 maps — all inline-initialized, never null)

#### Parity Notes
- All list fields are initialized to empty `ArrayList` — never null. Rust: `Vec<Sclip>` initialized empty.
- All map fields are initialized to empty `HashMap` — never null. Rust: `HashMap<String, ...>`.
- Integer fields default to 0 (Java primitive default). These are `end` positions — 0 means "not set".
- **No constructor** — Java default no-arg constructor initializes fields with their inline defaults.
- Commented out: `svinsfend`, `svinsrend`, `svfins`, `svrins` — insertion SV was planned but not implemented. Skip in Rust.

---

### SamView

**Source**: `data/SamView.java:L1-L42`  
**Purpose**: Thread-safe BAM reader wrapper with filter support.

#### Fields

| Field | Type | Description |
|-------|------|-------------|
| `threadLocalSAMReaders` | `static ThreadLocal<Map<String, SamReader>>` | Thread-local cache of open BAM readers |
| `iterator` | `SAMRecordIterator` | Current overlapping-interval iterator |
| `filter` | `int` | Hex filter from `-F` option (default `0x504`) |

#### Constructor
`SamView(String file, String samfilter, Region region, ValidationStringency stringency)`:
1. Opens/reuses BAM reader via `fetchReader()`
2. Creates overlapping query iterator for `region.chr:region.start-region.end`
3. Parses `samfilter` string to int via `Integer.decode()`

#### Methods

- `read()`: Returns next SAMRecord that passes filter `(record.getFlags() & filter) == 0`. Returns null at end.
- `close()`: Closes iterator (NOT the reader — reader stays cached).
- `fetchReader(...)`: Synchronized static. ThreadLocal cache of `SamReader` instances keyed by file path.

#### Parity Notes
- ThreadLocal readers may leak file handles — each thread keeps readers open until thread dies. Not a parity issue but affects resource management in Rust.
- `Integer.decode(samfilter)` handles hex (`0x504`) and decimal strings. Rust equivalent: `i32::from_str_radix()` or parse with prefix handling.
- `queryOverlapping` in htsjdk uses 1-based inclusive coordinates for the query interval.

---

### Side

**Source**: `data/Side.java:L1-L9`  
**Purpose**: Enum for SV directionality (3' vs 5' end).

#### Values
- `_3` — 3' end
- `_5` — 5' end
- `UNKNOWN` — fallback

#### Method
`static valueOf(int side)`: returns `_3` if side==3, `_5` if side==5, else `UNKNOWN`.

---

### SortPositionSclip

**Source**: `data/SortPositionSclip.java:L1-L19`  
**Purpose**: Temporary structure for sorting soft-clip candidates by position.

#### Fields

| Field | Type | Description |
|-------|------|-------------|
| `position` | `int` | Genomic position |
| `softClip` | `Sclip` | Associated soft-clip object |
| `count` | `int` | Soft-clip count at this position |

#### Constructor
`SortPositionSclip(int position, Sclip softClip, int count)`.

---

## Pipeline Data Flow

```
                SAMFileParser
                     │
                     ▼
           Scope<InitialData>
           ┌─────────────────────────┐
           │ nonInsertionVariants: HashMap<Integer, VariationMap>  ← outer HashMap, inner LinkedHashMap
           │ insertionVariants:    HashMap<Integer, VariationMap>
           │ refCoverage:          HashMap<Integer, Integer>
           │ softClips5End:        HashMap<Integer, Sclip>
           │ softClips3End:        HashMap<Integer, Sclip>
           └─────────────────────────┘
                     │ CigarParser
                     ▼
           Scope<VariationData>
           ┌─────────────────────────┐
           │ (same 5 maps from InitialData, transferred by reference)
           │ + positionToInsertionCount: HashMap<Integer, HashMap<String, Integer>>
           │ + positionToDeletionsCount: HashMap<Integer, HashMap<String, Integer>>
           │ + svStructures:        SVStructures
           │ + maxReadLength:       Integer (nullable)
           │ + splice:              Set<String>
           │ + mnp:                 HashMap<Integer, HashMap<String, Integer>>
           │ + spliceCount:         HashMap<String, int[2]>
           │ + duprate:             double
           └─────────────────────────┘
                     │ VariationRealigner + StructuralVariantsProcessor
                     ▼
       Scope<RealignedVariationData>
           ┌─────────────────────────┐
           │ (same maps, post-realignment)
           │ + CURSEG:              CurrentSegment (mutable!)
           │ + SOFTP2SV:            HashMap<Integer, List<Sclip>>
           │ + previousScope:       Scope<VariationData> (back-reference)
           └─────────────────────────┘
                     │ ToVarsBuilder
                     ▼
         Scope<AlignedVarsData>
           ┌─────────────────────────┐
           │ maxReadLength:         int (primitive now)
           │ alignedVariants:       Map<Integer, Vars>
           └─────────────────────────┘
                     │ Modes (Somatic only)
                     ▼
       Scope<CombineAnalysisData>
           ┌─────────────────────────┐
           │ maxReadLength:         int
           │ type:                  String
           └─────────────────────────┘
```

### Scope<T> Inheritance Pattern

At each stage transition, the new `Scope` is created via:
```java
new Scope<>(previousScope, newData)
```
This copies `bam`, `region`, `regionRef`, `referenceResource`, `maxReadLength`, `splice`, `out` from the previous scope and replaces only `data`. All common fields are **shared references** (not cloned).

---

## Patterns.java Deep Dive

(See Patterns class analysis above for complete 70-pattern inventory with regex strings, engines, and usage locations.)

### Key groupings:

1. **CIGAR modification patterns** (28-56): Used by CigarModifier to simplify/merge CIGAR operations. Many use jregex for group replacement.
2. **Variant description patterns** (2-8, 11-19): Used by ToVarsBuilder, VariationRealigner, and output printers to parse/construct variant description strings.
3. **SV notation patterns** (20-24): Used by StructuralVariantsProcessor and output formatting.
4. **Soft-clip quality patterns** (9-10): Filter poly-A/T artifacts.
5. **SA tag patterns** (60-68): Process supplementary alignment CIGAR strings.

---

## Reference.java + ReferenceResource.java

### Reference Fetching Algorithm

1. **Padding**: Region is extended by `numberNucleotideToExtend + extension` on both sides. Default extension = 1200bp. Clamped to `[1, chrLength]`.

2. **FASTA fetch**: Uses htsjdk `IndexedFastaSequenceFile.getSubsequenceAt(chr, start, end)` — 1-based inclusive.

3. **Caching**: `isLoaded()` checks if the requested region is fully contained within any previously loaded region. If so, returns immediately without re-fetching.

4. **Seed map construction** (k-mer indexing):
   - For each position `i` in the exon string:
     - Skip if position already in `referenceSequences` (overlap with prior load)
     - Store `referenceSequences[i + sequenceStart] = exon[i]`
     - If NOT at chromosome end special zone: generate two seed keys:
       - SEED_1 (17-mer): `exon[i..i+17]` (clipped if near end)
       - SEED_2 (12-mer): `exon[i..i+12]` (clipped if near end)
     - Each key maps to a list of positions

5. **Boundary behavior**:
   - Loop stops at `exon.length() - SEED_1` (17bp before end) unless at chromosome boundary
   - At chromosome boundary, loop goes to `exon.length()` but skips seed generation for last 17 positions
   - `Utils.substr()` clips to string end — truncated seeds still get indexed

### Seed Map Data Structure
```
seed: HashMap<String, ArrayList<Integer>>
  "ATCGATCGATCGATCGA" → [100234, 100891, 102445]  // SEED_1 = 17
  "ATCGATCGATCG"      → [100234, 100891, 102445, 103001]  // SEED_2 = 12
```

Used by `findMatch()` in VariationRealigner and StructuralVariantsProcessor for rapid breakpoint detection.

---

## Region.java

### Field Mapping

| Java Field | Type | Rust Type | Notes |
|------------|------|-----------|-------|
| `chr` | `String` | `String` | Chromosome name |
| `start` | `int` | `i32` | 1-based inclusive start |
| `end` | `int` | `i32` | 1-based inclusive end |
| `gene` | `String` | `String` | Gene name from BED |
| `insertStart` | `int` | `i32` | 0 = not set |
| `insertEnd` | `int` | `i32` | 0 = not set |

### Coordinate System

- **BED file**: 0-based, half-open `[start, end)`.
- **Region in code**: 1-based, inclusive `[start, end]`. The BED parser (`RegionBuilder`) converts by adding 1 to the start coordinate from BED.
- `insertStart` and `insertEnd`: 1-based inclusive when set (non-zero). Zero means "not applicable" (non-amplicon mode).
- `printRegion()` outputs `chr:start-end` in 1-based inclusive format — this is what appears in log messages.

### `insertionStr` Format
The `insertStart` and `insertEnd` fields are used in amplicon mode to define the "insert" (the region between primers). They come from columns 7-8 in the BED file (after the standard 4 columns + optional columns 5-6).

---

## Cross-Module Dependencies

- **Used by ALL**: `Scope<T>`, `GlobalReadOnlyScope`, `Region`, `Reference`, `Patterns` are used across every pipeline module
- **CigarParser produces**: `InitialData` → `VariationData`
- **VariationRealigner/SV produces**: `RealignedVariationData` (consumes `VariationData`)
- **ToVarsBuilder produces**: `AlignedVarsData` (consumes `RealignedVariationData`)
- **SomaticMode produces**: `CombineAnalysisData` (consumes `AlignedVarsData`)
- **CigarModifier uses**: `ModifiedCigar`
- **VariationRealigner uses**: `BaseInsertion`, `Match`, `Match35`
- **StructuralVariantsProcessor uses**: `SVStructures`, `CurrentSegment`, `Side`, `SortPositionSclip`
- **ReferenceResource**: fetches FASTA data and populates `Reference` objects
- **SamView**: provides BAM iteration to `SAMFileParser`

---

## Known Parity Traps

### Trap 1: GlobalReadOnlyScope Singleton → Rust Static/Thread-Local
Java's `volatile static` singleton is accessed everywhere via `instance()`. In Rust, this maps to either:
- `once_cell::sync::Lazy<GlobalReadOnlyScope>` (if truly immutable after init)
- `Arc<GlobalReadOnlyScope>` passed through function params (more idiomatic)
The `clear()` method for testing requires special handling — Rust's `Lazy` can't be reset.

### Trap 2: Scope<T> Generic Type Erasure → Rust Enum or Separate Types
Java erases generics at runtime. The `Scope<?>` wildcard in the inherit constructor accepts any Scope. In Rust, options:
- Separate structs per stage (most common approach)
- Enum wrapping different data types
- Generic struct `Scope<T>` (works in Rust too, but the inherit constructor needs trait bounds)

### Trap 3: Reference Seed Map Iteration Order
`seed` is `HashMap<String, List<Integer>>` — the map itself is unordered. But within each key, the `List<Integer>` is ordered by insertion (left-to-right position). Consumers like `findMatch()` iterate this list sequentially and take the **first** acceptable match. Rust must preserve list ordering (Vec<i32> — which it does naturally).

### Trap 4: ReferenceResource FASTA Padding and Boundary Handling
- `sequenceStart` calculation: `max(1, region.start - numberNucleotideToExtend - extension)`. Note the `< 1` comparison, not `<= 0`.
- `sequenceEnd` calculation: `min(len, region.end + numberNucleotideToExtend + extension)`. Note `> len` comparison, not `>= len`.
- Chromosome length `len` comes from `chrLengths.get(chr)` which returns `null` (boxed to 0) if chr not found. In Rust, `chrLengths.get(chr).copied().unwrap_or(0)`.
- `siteEnd` when at chromosome end uses `exon.length()` not `exon.length() - 1`.

### Trap 5: Region Coordinate Convention
BED files are 0-based half-open, but Region's `start` is 1-based inclusive (after parsing). This +1 happens in `RegionBuilder`, NOT in Region's constructor. Any code constructing Region directly must use 1-based coordinates. `end` is NOT modified — it's 1-based inclusive (BED's 0-based exclusive end numerically equals 1-based inclusive end).

### Trap 6: Patterns — jregex vs java.util.regex
17 patterns use `jregex.Pattern`, 53 use `java.util.regex.Pattern`. In Rust, all become `regex::Regex`. Key differences:
- `jregex.Replacer` group replacement syntax differs from both Java standard and Rust regex
- `jregex` `Matcher.find()` + group capture works like Java's but `jregex` may handle certain edge cases differently (empty matches, overlapping matches)
- Pattern `BEGIN_NUM_S_OR_BEGIN_NUM_H` uses alternation `^(\d+)S|^\d+H` — both alternatives are anchored at start. Rust handles this correctly.

### Trap 7: SVStructures Default Empty Collections vs Null
All list/map fields in SVStructures are initialized inline (never null). In Rust, these should be initialized as empty `Vec`/`HashMap`. No `Option` wrapping needed. However, `VariationData.svStructures` itself CAN be null — it's assigned after creation. So `Option<SVStructures>` is needed at the VariationData level.

### Trap 8: InitialData — All Maps are HashMap (Outer)
All 5 maps in `InitialData` use `HashMap` for the outer position key. The inner `VariationMap<String, Variation>` is `LinkedHashMap` (insertion-ordered). Rust mapping:
- Outer: `HashMap<i32, VecMap<String, Variation>>` or similar
- Inner: `IndexMap<String, Variation>` or a custom ordered map

### Trap 9: RealignedVariationData — SOFTP2SV Map Type
`SOFTP2SV` is `Map<Integer, List<Sclip>>` — concrete type at construction is whatever the caller creates (typically `HashMap`). Verify in StructuralVariantsProcessor where it's populated. Iteration order may not matter since it's consumed by position-based lookups.

### Trap 10: Field Default Values — Java null vs Rust Option
- `VariationData` has **all** reference fields defaulting to null (no constructor). In Rust, all must be `Option<T>` or the struct must have a builder pattern.
- `RealignedVariationData` has all fields as `final` — set once in constructor, never null. Rust can use plain types (no Option needed).
- `BaseInsertion` fields are boxed `Integer` but always initialized from `int` params — use `i32` in Rust.
- `AlignedVarsData.maxReadLength` is primitive `int` (not nullable) — use `i32` in Rust.
- `VariationData.maxReadLength` is boxed `Integer` (nullable) — use `Option<i32>` in Rust.

### Trap 11: VariationData Has No Constructor
Java's default no-arg constructor initializes all fields to their defaults: `null` for objects, `0.0` for `double`. In Rust, you need a `Default` impl or explicit initialization. Every user of `VariationData` sets fields individually after `new VariationData()`.

### Trap 12: RealignedVariationData.previousScope Back-Reference
Creates `Scope<RealignedVariationData>` → `RealignedVariationData.previousScope: Scope<VariationData>` — a reference from a later stage back to an earlier one. In Rust, this creates ownership complexity. Common solutions:
- `Arc<Scope<VariationData>>` shared ownership
- Extract needed data at construction time instead of holding the reference
- Investigate whether `previousScope` is actually read (it IS — by SomaticMode for combining analyses)

### Trap 13: Utils.substr() Clipping in Seed Map
`Utils.substr(exon, i, SEED_1)` returns a substring of up to SEED_1 characters, but will return fewer if near the end of the string. These truncated seeds are still added to the seed map. This means shorter-than-expected keys exist in the map. The `findMatch()` callers must handle this — they look up full-length seeds, so truncated seeds in the map won't cause false matches (they'll never be looked up). But if any code generates a lookup key from a truncated source, it could match unexpectedly.

### Trap 14: CurrentSegment Mutability Through Final Reference
`RealignedVariationData.CURSEG` is `final` (the reference) but `CurrentSegment` fields (`chr`, `start`, `end`) are mutable. Code in StructuralVariantsProcessor modifies `CURSEG.start` and `CURSEG.end` in place. In Rust, this means the `CurrentSegment` inside `RealignedVariationData` needs interior mutability or must be passed as `&mut`.

---

## Appendix: Configuration Constants Referenced

| Constant | Defined In | Value | Used By |
|----------|-----------|-------|---------|
| `SEED_1` | Configuration | 17 | ReferenceResource — seed map |
| `SEED_2` | Configuration | 12 | ReferenceResource — seed map |
| `ADSEED` | Configuration | 6 | Adaptor detection |
| `EXTENSION` | Configuration | 5000 | SV extended reference |
| `referenceExtension` | Configuration | 1200 (default) | ReferenceResource — padding |
| `numberNucleotideToExtend` | Configuration | CLI param `-x` | ReferenceResource — padding |
| `SVMAXLEN` | Configuration | 150000 | StructuralVariantsProcessor |
| `SVMINLEN` | Configuration | 1000 (default, `-L`) | StructuralVariantsProcessor |
| `SVFLANK` | Configuration | 50 | StructuralVariantsProcessor |
| `LOWQUAL` | Configuration | 10 | Soft-clip quality threshold |
