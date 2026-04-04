# SAMFileParser

**Source**: `modules/SAMFileParser.java` (17 LOC), `modules/RecordPreprocessor.java` (~270 LOC), `data/SamView.java` (55 LOC)  
**LOC**: ~342 total  
**Rust counterpart**: `src/mods/vardict_pipeline.rs` (`passes_preprocess()`, `collect_filtered_records()`, `RecordPreprocessorState`)  
**Risk**: MEDIUM — read filtering affects all downstream variant detection  
**Pipeline Stage**: First module — `AbstractMode → SAMFileParser.process() → CigarParser`  
**Status**: complete

## Overview

SAMFileParser is the entry gate of VarDict's per-region pipeline. It consists of three tightly coupled classes:

1. **SAMFileParser** (17 LOC): A trivial `Module<InitialData, RecordPreprocessor>` implementation that splits the BAM path string by `:`, creates a RecordPreprocessor, and wraps it in a new Scope.

2. **RecordPreprocessor** (~270 LOC): The iterator that CigarParser consumes. It manages a stack of BAM file paths (for multi-BAM tumor/normal support), opens SamView readers one at a time, and applies a cascade of read filters: downsampling → mapping quality → secondary alignment → no-sequence → totalReads increment → duplicate detection.

3. **SamView** (55 LOC): Wraps htsjdk's `SamReader` with thread-local caching, region-based overlapping queries, and SAM flag filtering.

The filter cascade's exact order is parity-critical because `totalReads` is incremented at a specific point in the sequence, and changing the filter order changes the counter value.

## Method Inventory

### SAMFileParser.java

| Method | Lines | Analyzed? | Summary |
|--------|-------|-----------|---------|
| `process()` | L10–L16 | yes | Splits BAM path, creates RecordPreprocessor, wraps in Scope |

### RecordPreprocessor.java

| Method | Lines | Analyzed? | Summary |
|--------|-------|-----------|---------|
| `RecordPreprocessor()` (ctor) | L44–L68 | yes | Init maps, BAM deque (push reverses), JSONL env, calls nextReader() |
| `nextRecord()` | L73–L79 | yes | Public iterator: returns next passing record across BAMs |
| `nextReader()` | L84–L97 | yes | Pops next BAM, opens SamView, resets dup state |
| `nextOnCurrentReader()` | L103–L115 | yes | Reads until a record passes preprocessRecord() |
| `close()` | L117–L125 | yes | Closes reader, writes diagnostic JSONL |
| `preprocessRecord()` | L131–L182 | yes | Core filter cascade (downsampling, mapq, flags, seq, dups) |
| `appendJsonlRecord()` | L184–L213 | yes | Diagnostic JSONL (not parity-relevant) |
| `getChrName()` | L219–L228 | yes | Static: strips "chr" prefix if -C flag set |
| `escapeJson()` | L230–L258 | yes | Static: JSON string escaping (diagnostic only) |

### SamView.java

| Method | Lines | Analyzed? | Summary |
|--------|-------|-----------|---------|
| `SamView()` (ctor) | L20–L24 | yes | Opens queryOverlapping iterator, parses hex filter |
| `read()` | L30–L39 | yes | Returns next record passing SAM flag filter |
| `close()` | L42–L44 | yes | Closes the iterator |
| `fetchReader()` | L46–L51 | yes | Thread-local SamReader cache (computeIfAbsent) |

## Method Analyses

### SAMFileParser.process()

**Source**: `modules/SAMFileParser.java:L10–L16`  
**Purpose**: Pipeline entry point — creates RecordPreprocessor and wraps in Scope.  
**Called By**: `AbstractMode.pipeline()`, `partialPipeline()`, `splicingPipeline()`

#### Parameters
- `scope`: `Scope<InitialData>` — contains `bam` (String), `region` (Region), `data` (InitialData)

#### Algorithm
1. Read `scope.bam` — e.g. `"/path/tumor.bam"` or `"/path/tumor.bam:/path/normal.bam"`
2. Call `scope.bam.split(":")` → `String[]` of BAM file paths
3. Create `new RecordPreprocessor(bams, scope.region, scope.data)`
4. Create `new Scope<>(scope, preprocessor)` — inherits all scope fields (`bam`, `region`, `regionRef`, `referenceResource`, `maxReadLength`, `splice`, `out`), sets `data = preprocessor`
5. Return the new scope

#### Null/Edge Cases
- `scope.bam` null → NullPointerException on `.split(":")`
- `scope.bam` empty `""` → `split(":")` returns `[""]` → htsjdk error opening empty string as BAM path

---

### RecordPreprocessor() constructor

**Source**: `modules/RecordPreprocessor.java:L44–L68`  
**Purpose**: Initialize all fields and open first BAM reader.

#### Parameters
- `bams`: `String[]` — BAM paths from `scope.bam.split(":")`
- `region`: `Region` — genomic region (chr, start, end, gene)
- `data`: `InitialData` — contains 5 mutable variant/coverage maps (shared by reference, NOT copied)

#### Algorithm
1. Assign `this.nonInsertionVariants = data.nonInsertionVariants` — direct reference, not copy
2. Assign `this.insertionVariants = data.insertionVariants`
3. Assign `this.refCoverage = data.refCoverage`
4. Assign `this.softClips5End = data.softClips5End`
5. Assign `this.softClips3End = data.softClips3End`
6. Create `this.bams = new ArrayDeque<>()`
7. For each `bam` in `bams` array:
   - `this.bams.push(bam)` — **push adds to the HEAD (front) of the deque**, reversing array order
   - For array `[A, B]`, after loop: deque is `[B, A]` (B=head, A=tail)
8. Assign `this.region = region`
9. Read environment variable `VARDICT_RECORD_PREPROCESSOR_JSONL`:
   - If non-null and non-blank: `this.jsonlPath = jsonl.trim()`, `this.jsonlEntries = new ArrayList<>()`
   - Else: both set to `null`
10. Call `this.nextReader()` — opens the first BAM (from deque tail = original first element)

#### Mutable State Changes
- `this.bams`: populated with BAM paths in reversed order (push = LIFO)
- `this.currentReader`: set by `nextReader()` call
- `this.duplicates`: set to `new HashSet<>()` by `nextReader()`
- `this.firstMatchingPosition`: set to `-1` by `nextReader()`
- `this.totalReads`: `0` (default int)
- `this.duplicateReads`: `0` (default int)

**Key insight**: `push()` reverses array order into deque, but `nextReader()` uses `pollLast()` (removes from tail), so BAMs are processed in **original array order**.

---

### RecordPreprocessor.nextRecord()

**Source**: `modules/RecordPreprocessor.java:L73–L79`  
**Purpose**: Public iterator API — returns next filtered SAMRecord, spanning across all BAM files.

#### Parameters
- None

#### Algorithm
1. Call `nextOnCurrentReader()` → assign to `record`
2. If `record == null` (current BAM exhausted):
   a. Call `nextReader()` — opens the next BAM file from the deque
   b. Call `nextOnCurrentReader()` again → assign to `record`
3. Return `record` (may be `null` if all BAMs are exhausted)

#### Mutable State Changes
- May change `currentReader`, `duplicates`, `firstMatchingPosition` via `nextReader()`

#### Return Value
- `SAMRecord` that passed all filters, or `null` if all BAMs exhausted

#### Null/Edge Cases
- If `bams` deque is empty when `nextReader()` is called: `nextReader()` returns early, `nextOnCurrentReader()` returns `null` on the already-exhausted reader
- **Limitation**: Only tries ONE additional BAM. If there are 3+ BAMs and the 2nd has zero passing records, the 3rd is never opened. In practice only 1–2 BAMs are used.

---

### RecordPreprocessor.nextReader()

**Source**: `modules/RecordPreprocessor.java:L84–L97`  
**Purpose**: Pop next BAM from deque, open SamView, reset duplicate-detection state.

#### Parameters
- None

#### Algorithm
1. If `bams.isEmpty()`: return immediately (no more BAMs)
2. If `currentReader != null`: `currentReader.close()` — closes htsjdk iterator for previous BAM
3. `bamPath = bams.pollLast()` — removes and returns the **tail** element of the deque
4. Create `currentReader = new SamView(bamPath, instance().conf.samfilter, region, instance().conf.validationStringency)`
   - `conf.samfilter` default: `"0x504"`
   - `conf.validationStringency` default: `LENIENT`
5. `duplicates = new HashSet<>()` — fresh duplicate set for new BAM
6. `firstMatchingPosition = -1` — reset position tracking

#### Mutable State Changes
- `currentReader`: old one closed, replaced with new SamView
- `duplicates`: new empty HashSet
- `firstMatchingPosition`: reset to -1
- `bams`: one element removed from tail

#### Null/Edge Cases
- On first call (from constructor): `currentReader` is `null`, so `close()` is skipped

---

### RecordPreprocessor.nextOnCurrentReader()

**Source**: `modules/RecordPreprocessor.java:L103–L115`  
**Purpose**: Read records from current SamView until one passes the filter cascade, or EOF.

#### Parameters
- None

#### Algorithm
1. Enter infinite loop:
2. Call `currentReader.read()` → assign to `record`
3. If `record == null`: return `null` (BAM exhausted for this region)
4. Call `preprocessRecord(record)` → `passed` (boolean)
5. Call `appendJsonlRecord(record, passed)` — diagnostic logging (if enabled)
6. If `passed == true`: return `record`
7. Else: continue loop (skip record)

#### Mutable State Changes
- `totalReads`, `duplicateReads`: potentially incremented by `preprocessRecord()`
- `duplicates`, `firstMatchingPosition`: potentially modified by `preprocessRecord()`

#### Return Value
- `SAMRecord` that passed all filters, or `null`

---

### RecordPreprocessor.preprocessRecord() ⭐ CRITICAL

**Source**: `modules/RecordPreprocessor.java:L131–L182`  
**Purpose**: Core read filter cascade. Most parity-critical method in this module.

#### Parameters
- `record`: `SAMRecord` — the read to evaluate

#### Algorithm

**Step 1 — Downsampling** (L132–L134):
1. If `instance().conf.isDownsampling()` (i.e., `conf.downsampling != null`):
2. Generate `RND.nextDouble()` — returns value in `[0.0, 1.0)`
3. If `RND.nextDouble() <= instance().conf.downsampling`: return `false` (skip this read)
4. `RND` is `private final static Random RND = new Random(System.currentTimeMillis())` — class-level static, time-seeded, shared across all instances and threads

**Step 2 — Extract read attributes** (L135–L136):
5. `querySequence = record.getReadString()` — read bases as String
6. `mappingQuality = record.getMappingQuality()` — MAPQ integer from SAM

**Step 3 — Mapping quality filter** (L139–L141):
7. If `instance().conf.hasMappingQuality()` (i.e., `conf.mappingQuality != null`):
8. If `mappingQuality < instance().conf.mappingQuality`: return `false`
9. Comparison is **strict less-than** (`<`): reads with MAPQ exactly equal to threshold **pass**
10. `conf.mappingQuality` is boxed `Integer`, checked with `hasMappingQuality()` which returns `mappingQuality != null`

**Step 4 — Secondary alignment filter** (L144–L146):
11. If `record.isSecondaryAlignment()` (SAM flag `0x100`) AND `!instance().conf.samfilter.equals("0")`:
12. Return `false`
13. `isSecondaryAlignment()` = `(flags & 0x100) != 0`
14. The condition `!samfilter.equals("0")` is a **string equality** check against the literal string `"0"` — NOT numeric. `"0x0"`, `"00"`, `"0x000"` would all NOT match.

**Step 5 — No-sequence filter** (L148–L150):
15. If `querySequence.length() == 1 && querySequence.charAt(0) == '*'`: return `false`
16. htsjdk returns `"*"` (length 1) when the SEQ field in SAM is `*`
17. Only checks for single-character `"*"` — a normal read cannot have length 1

**Step 6 — Increment totalReads** (L151):
18. `totalReads++` — counts reads that passed steps 1–5 (including reads that may later be filtered as duplicates)

**Step 7 — Get mate reference name** (L153):
19. `mateReferenceName = getMateReferenceName(record)` — calls `CigarParser.getMateReferenceName()`
20. Returns `"*"` if `record.getMateReferenceName()` is `null`
21. Returns `"="` if `record.getReferenceName().equals(record.getMateReferenceName())` (same chromosome)
22. Otherwise returns the actual mate reference name string

**Step 8 — Duplicate detection** (L156–L181) — only if `instance().conf.removeDuplicatedReads` is `true`:

**8a — Position change check** (L157–L159):
23. If `record.getAlignmentStart() != firstMatchingPosition`:
24. `duplicates.clear()` — resets the dup set when reads move to a new alignment position
25. `getAlignmentStart()` returns 1-based position

**8b — Branch 1: mateAlignmentStart < 10** (L160–L169):
26. Condition: `record.getMateAlignmentStart() < 10`
27. Catches: unmapped mates (start=0), mates near contig start (1–9)
28. Build `dupKey = alignmentStart + "-" + mateReferenceName + "-" + mateAlignmentStart`
29. Example: `"12345-=-0"`, `"12345-*-0"`
30. If `duplicates.contains(dupKey)`: `duplicateReads++`, return `false`
31. Else: `duplicates.add(dupKey)`, `firstMatchingPosition = record.getAlignmentStart()`

**8c — Branch 2: paired + mate unmapped** (L170–L179):
32. Condition: `record.getReadPairedFlag()` (flag `0x1`) AND `record.getMateUnmappedFlag()` (flag `0x8`)
33. Build `dupKey = alignmentStart + "-" + record.getCigarString()`
34. Example: `"12345-100M"`
35. Same contains/add/firstMatchingPosition logic as Branch 1

**8d — Branch 3: neither** (implicit):
36. If `mateAlignmentStart >= 10` AND NOT (paired + mate unmapped): **no duplicate check**
37. Read passes through regardless

**Step 9 — Return true** (L182):
38. Read passed all filters

#### Mutable State Changes
- `totalReads`: incremented at step 6
- `duplicateReads`: potentially incremented at steps 8b or 8c
- `duplicates`: potentially cleared (step 8a), added to (8b/8c)
- `firstMatchingPosition`: potentially updated (8b/8c)

#### Return Value
- `true`: record should be processed by CigarParser
- `false`: record is filtered out

---

### RecordPreprocessor.close()

**Source**: `modules/RecordPreprocessor.java:L117–L125`  
**Purpose**: Close BAM reader and optionally write diagnostic JSONL.

#### Algorithm
1. `currentReader.close()` — closes htsjdk iterator
2. If `jsonlEntries != null && jsonlPath != null`:
   a. Call `RecordPreprocessorJsonlWriter.write(jsonlPath, region, totalReads, duplicateReads, jsonlEntries)`
   b. On exception: print to stderr, do NOT rethrow

#### Parity Impact
None — JSONL is diagnostic only.

---

### RecordPreprocessor.getChrName() (static)

**Source**: `modules/RecordPreprocessor.java:L219–L228`  
**Purpose**: Optionally strip `"chr"` prefix from chromosome name based on `-C` flag.

#### Parameters
- `region`: `Region` — has `.chr` field

#### Algorithm
1. If `instance().conf.chromosomeNameIsNumber` is `true` AND `region.chr.startsWith("chr")`:
   - Return `region.chr.substring(3)` (removes "chr" prefix)
2. Else: return `region.chr` unchanged

#### Edge Cases
- `region.chr == "chr"` (just prefix, no number) → returns `""`
- `region.chr == "chrUn_gl000220"` → returns `"Un_gl000220"`

---

### RecordPreprocessor.appendJsonlRecord() (Diagnostic)

**Source**: `modules/RecordPreprocessor.java:L184–L213`  
**Purpose**: Append diagnostic JSONL entry for each record. Only active when env `VARDICT_RECORD_PREPROCESSOR_JSONL` is set. Not parity-relevant.

---

### RecordPreprocessor.escapeJson() (Diagnostic)

**Source**: `modules/RecordPreprocessor.java:L230–L258`  
**Purpose**: JSON string escaping for diagnostic output. Handles `"`, `\`, `\n`, `\r`, `\t`. Not parity-relevant.

---

### SamView() constructor

**Source**: `data/SamView.java:L20–L24`  
**Purpose**: Open a region-overlapping iterator on BAM file with flag filtering.

#### Parameters
- `file`: `String` — BAM file path
- `samfilter`: `String` — hex string for flag filter (e.g., `"0x504"`)
- `region`: `Region` — genomic region
- `stringency`: `ValidationStringency` — htsjdk validation level

#### Algorithm
1. `reader = fetchReader(file, stringency)` — gets or creates thread-local `SamReader`
2. `iterator = reader.queryOverlapping(region.chr, region.start, region.end)` — 1-based inclusive coordinates
3. `filter = Integer.decode(samfilter)` — handles hex (`0x` prefix), octal (`0` prefix), and decimal
   - `Integer.decode("0x504")` → `1284` decimal
   - Default `0x504` = `0101 0000 0100` binary: bit 2 (unmapped) + bit 8 (secondary) + bit 10 (duplicate)

---

### SamView.read()

**Source**: `data/SamView.java:L30–L39`  
**Purpose**: Return next record passing SAM flag filter.

#### Algorithm
1. While `iterator.hasNext()`:
2. `record = iterator.next()`
3. If `filter != 0` AND `(record.getFlags() & filter) != 0`: continue (skip — any filter bit set in record flags)
4. Return `record`
5. If loop ends (no more records): return `null`

#### Return Value
- `SAMRecord` that has no overlap with filter flags, or `null`

---

### SamView.close()

**Source**: `data/SamView.java:L42–L44`  
**Purpose**: Close the htsjdk iterator. Does NOT close the underlying SamReader (cached per-thread).

#### Algorithm
1. `iterator.close()`

---

### SamView.fetchReader() (static, synchronized)

**Source**: `data/SamView.java:L46–L51`  
**Purpose**: Thread-local `SamReader` cache to avoid reopening BAM files.

#### Algorithm
1. Get thread-local `Map<String, SamReader>` via `threadLocalSAMReaders.get()`
2. `computeIfAbsent(file, f -> SamReaderFactory.makeDefault().validationStringency(stringency).open(SamInputResource.of(f)))`
3. If key `file` exists: return cached reader
4. Else: create new `SamReader`, store in map, return it

#### Parity-Critical Details
- `SamReader` is **never explicitly closed** — persists for thread lifetime
- `synchronized` is redundant given ThreadLocal but prevents races during class initialization
- One reader per file per thread — across multiple regions processed on the same thread, the reader is reused

---

## Cross-Module Dependencies

### Outbound Calls

| Target | Class/File | Method(s) | Purpose |
|--------|-----------|-----------|---------|
| SamView | `data/SamView.java` | constructor, `read()`, `close()` | BAM reading + flag filter |
| CigarParser | `modules/CigarParser.java` | `getMateReferenceName()` (static) | Mate ref normalization: null→"*", same-chr→"=", else name |
| GlobalReadOnlyScope | `data/scopedata/GlobalReadOnlyScope.java` | `instance()` | Access configuration |
| Configuration | `Configuration.java` | `.samfilter`, `.validationStringency`, `.downsampling`, `.mappingQuality`, `.removeDuplicatedReads`, `.chromosomeNameIsNumber` | Config field reads |
| htsjdk SamReaderFactory | external | `makeDefault()`, `.open()` | BAM file opening |
| htsjdk SAMRecord | external | `getReadString()`, `getMappingQuality()`, `isSecondaryAlignment()`, `getAlignmentStart()`, `getMateAlignmentStart()`, `getMateReferenceName()`, `getReadPairedFlag()`, `getMateUnmappedFlag()`, `getCigarString()`, `getFlags()`, `getReadName()`, `getBaseQualityString()`, `getReferenceName()` | Read attribute access |

### Inbound Callers

| Caller | Class/File | Method | How |
|--------|-----------|--------|-----|
| AbstractMode | `modes/AbstractMode.java` | `pipeline()` | `new SAMFileParser().process(scope)` |
| AbstractMode | `modes/AbstractMode.java` | `partialPipeline()` | `new SAMFileParser().process(scope)` |
| AbstractMode | `modes/AbstractMode.java` | `splicingPipeline()` | `new SAMFileParser().process(scope)` |
| CigarParser | `modules/CigarParser.java` | `process()` main loop | Calls `scope.data.nextRecord()` repeatedly until `null` |

## Data Structures Read/Written

| Structure | Java Type | Field Location | R/W | Notes |
|-----------|-----------|----------------|-----|-------|
| `nonInsertionVariants` | `Map<Integer, VariationMap<String, Variation>>` | RecordPreprocessor field | R (pass-through) | Shared ref from InitialData; populated later by CigarParser |
| `insertionVariants` | `Map<Integer, VariationMap<String, Variation>>` | RecordPreprocessor field | R (pass-through) | Same |
| `refCoverage` | `Map<Integer, Integer>` | RecordPreprocessor field | R (pass-through) | Same |
| `softClips5End` | `Map<Integer, Sclip>` | RecordPreprocessor field | R (pass-through) | Same |
| `softClips3End` | `Map<Integer, Sclip>` | RecordPreprocessor field | R (pass-through) | Same |
| `bams` | `Deque<String>` (ArrayDeque) | RecordPreprocessor field | W | push in ctor, pollLast in nextReader |
| `duplicates` | `HashSet<String>` | RecordPreprocessor field | W | Cleared on position change, keys added per read |
| `firstMatchingPosition` | `int` | RecordPreprocessor field | W | Updated when dup key is inserted |
| `totalReads` | `int` | RecordPreprocessor field | W | Counts non-downsampled, mapped, primary, has-seq reads |
| `duplicateReads` | `int` | RecordPreprocessor field | W | Subset of totalReads that are duplicates |
| `currentReader` | `SamView` | RecordPreprocessor field | W | Replaced on each nextReader() call |
| `threadLocalSAMReaders` | `ThreadLocal<Map<String, SamReader>>` | SamView static | W | Cached per-thread SamReader objects |
| `iterator` | `SAMRecordIterator` | SamView field | W | Iterator consumed by read() calls |
| `filter` | `int` | SamView field | R | Parsed from samfilter hex string at construction |

## Configuration Fields Accessed

| Field | Java Type | CLI Flag | Default | Where Used |
|-------|-----------|----------|---------|------------|
| `conf.samfilter` | `String` | `-F` | `"0x504"` | SamView constructor (Integer.decode), preprocessRecord (String.equals) |
| `conf.validationStringency` | `ValidationStringency` | `-VS` | `LENIENT` | SamView constructor |
| `conf.downsampling` | `Double` (boxed) | `-Z` | `null` | preprocessRecord step 1 |
| `conf.mappingQuality` | `Integer` (boxed) | `-Q` | `null` | preprocessRecord step 3 |
| `conf.removeDuplicatedReads` | `boolean` | `-t` | `false` | preprocessRecord step 8 |
| `conf.chromosomeNameIsNumber` | `boolean` | `-C` | `false` | getChrName() |

## Known Parity Traps

1. **BAM path split delimiter**: `scope.bam.split(":")` — BAM paths cannot contain colons. Rust must use the same `:` delimiter.

2. **BAM processing order**: `push()` reverses array order into deque, then `pollLast()` re-reverses = original array order. Rust must process BAMs in the same order since variant maps accumulate across BAMs.

3. **SamView flag filter parsing**: `Integer.decode("0x504")` handles hex/octal/decimal prefixes. Rust must parse the samfilter string with equivalent semantics. Default `0x504` = 1284 = unmapped + secondary + PCR-duplicate.

4. **Flag filter precedes preprocessRecord**: `SamView.read()` applies the `-F` filter BEFORE `preprocessRecord()` sees the record. The Rust equivalent must apply the flag filter at the same pipeline point, not after mapping quality or other filters.

5. **Filter cascade order**: The exact order is: downsampling → mapping quality → secondary alignment → no-sequence → `totalReads++` → duplicate detection. Changing this order changes the `totalReads` counter value, which propagates downstream.

6. **Downsampling non-determinism**: `new Random(System.currentTimeMillis())` — seeded with wall-clock time, non-reproducible. **Downsampling must be disabled** (`-Z` not set) for parity testing.

7. **Secondary alignment samfilter string check**: `!instance().conf.samfilter.equals("0")` is a **string equality** test against the literal `"0"`, NOT a numeric comparison. The strings `"0x0"`, `"00"`, `"000"` would NOT match. Rust must replicate as string comparison.

8. **`getMateAlignmentStart() < 10` heuristic**: The duplicate detection Branch 1 triggers for mate alignment positions 0–9, not a clean unmapped check. Unmapped reads have position 0, but reads mapped to positions 1–9 also trigger this branch. Rust must use the same `< 10` threshold with 1-based coordinates.

9. **Duplicate key format**: The two formats must match exactly:
   - Branch 1: `"{alignmentStart}-{mateRefName}-{mateAlignmentStart}"` (e.g., `"12345-=-0"`)
   - Branch 2: `"{alignmentStart}-{cigarString}"` (e.g., `"12345-100M"`)

10. **Mate reference name normalization**: `CigarParser.getMateReferenceName()` returns `"*"` for null, `"="` for same-chromosome, else the actual reference name. This normalization is used in BOTH RecordPreprocessor (duplicate key) and CigarParser (SV processing). The Rust side must use the same normalization consistently.

11. **Duplicate set clearing on position change**: `duplicates.clear()` when `alignmentStart != firstMatchingPosition`. This assumes reads are sorted by position (standard for BAM). If a BAM is not sorted, duplicate detection breaks silently.

12. **`firstMatchingPosition` only updated in dup branches**: The position is only updated when a duplicate key is actually added (Branch 1 or Branch 2). If a read falls in Branch 3 (neither condition), `firstMatchingPosition` is NOT updated. This means the dup set may persist across reads at multiple positions if intermediate reads don't trigger Branch 1 or 2.

13. **`totalReads` includes then-filtered duplicates**: `totalReads++` happens at step 6, BEFORE duplicate detection at step 8. Reads that are later detected as duplicates are included in `totalReads`. So `totalReads - duplicateReads` = non-duplicate passing reads.

14. **SamReader never explicitly closed**: The `ThreadLocal<Map<String, SamReader>>` cache persists for the thread's lifetime. Only iterators are closed via `SamView.close()`. In Rust with `rust-htslib`, reader lifecycle may differ — ensure BAM index files aren't unnecessarily reopened.

15. **`queryOverlapping` uses 1-based coordinates**: htsjdk's `queryOverlapping(chr, start, end)` uses 1-based inclusive coordinates matching `Region.start` and `Region.end`. In rust-htslib, `fetch()` may use 0-based coordinates. The Rust side must convert appropriately.

16. **`getReadString()` returns `"*"` for no sequence**: htsjdk returns the literal string `"*"` (length 1) when the SEQ field is `*`. The Java check is `length() == 1 && charAt(0) == '*'`. In rust-htslib, `record.seq().len()` returns 0 for reads with no sequence, not 1. The Rust filter must account for this difference.

17. **Static RNG shared across instances**: `RND` is `private final static` — class-level, shared across ALL RecordPreprocessor instances. If multiple regions process concurrently on the same thread, they share the same RNG. The Rust equivalent should use `rand::thread_rng()` or similar.

18. **`getAlignmentStart()` is 1-based**: htsjdk returns 1-based positions. rust-htslib's `record.pos()` returns 0-based. The Rust side must add 1 when computing duplicate keys and when comparing against `firstMatchingPosition`.

19. **Multi-BAM single retry**: `nextRecord()` tries `nextOnCurrentReader()`, and if null calls `nextReader()` then retries ONCE. If there are 3+ BAMs and the second BAM is empty, the third BAM is never opened. In practice VarDict uses at most 2 BAMs (tumor + normal), so this is not a real issue.
