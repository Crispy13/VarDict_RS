# CigarParser

**Source**: `modules/CigarParser.java`
**LOC**: 2,662
**Rust counterpart**: `src/mods/cigar_parser.rs`
**Risk**: HIGH — core variant detection from CIGAR strings
**Pipeline Stage**: After CigarModifier, before VariationRealigner
**Status**: complete

## Overview

CigarParser is the core variant detection engine of VarDictJava. It iterates over BAM records produced by `RecordPreprocessor`, interprets each record's CIGAR string operation-by-operation, and populates variant maps (insertions, non-insertions), reference coverage, soft-clip structures, and structural variant structures. Every SNP, MNV, insertion, deletion, and complex variant that VarDict ultimately calls flows through this module. It implements the `Module<RecordPreprocessor, VariationData>` interface with a single `process()` entry point.

The module processes each CIGAR element in a `switch` on the operator: `M` (match/mismatch) is handled inline in the main loop, while `N` (intron), `S` (soft clip), `H` (hard clip), `I` (insertion), and `D` (deletion) are dispatched to dedicated private methods. During match processing, it detects MNVs by scanning ahead for consecutive mismatches and complex variants by looking for nearby indels within `conf.vext` bases of the end of a matched segment.

## Method Inventory

| Method | Lines | Analyzed? | Summary |
|--------|-------|-----------|---------|
| `process()` | L82–L171 | **yes** | Entry point: loops over records, calls `parseCigar()`, builds `VariationData` output |
| `initFromScope()` | L173–L185 | **yes** | Copies scope fields into instance variables |
| `parseCigar()` | L237–L653 | **yes** | Main per-record CIGAR parsing loop — the core of the module |
| `processNotMatched()` | L1307–L1320 | **yes** | Handles N operator (intron/splice junctions) |
| `processSoftClip()` | L1124–L1304 | **yes** | Handles S operator (5' and 3' soft clips, chimeric detection) |
| `processInsertion()` | L930–L1121 | **yes** | Handles I operator (insertion variants) |
| `processDeletion()` | L672–L883 | **yes** | Handles D operator (deletion variants) |
| `addVariationForMatchingPart()` | L1710–L1776 | **yes** | Creates/updates Variation for M-segment variants (SNPs, MNVs, complex) |
| `addVariationForDeletion()` | L1781–L1835 | **yes** | Creates/updates Variation for deletion variants |
| `prepareSVStructuresForAnalysis()` | L2017–L2310 | **yes** | Fills SV structures (DEL/DUP/INV/FUS) from discordant pairs |
| `addSV()` | L2319–L2377 | **yes** | Adds a single SV observation to an Sclip accumulator |
| `sclip5HighQualityProcessing()` | L1969–L2014 | **yes** | Populates 5' soft-clip Sclip structures |
| `sclip3HighQualityProcessing()` | L1935–L1967 | **yes** | Populates 3' soft-clip Sclip structures |
| `addCnt()` | L1429–L1445 | **yes** | Increments Variation counters |
| `subCnt()` | L1406–L1424 | **yes** | Decrements Variation counters |
| `findOffset()` | L1502–L1547 | **yes** | Scans forward for mismatches to combine with indels |
| `cleanupCigar()` | L1553–L1599 | **yes** | Converts leading/trailing H→remove, I→S in CIGAR |
| `getCigarOperator()` | L1601–L1609 | **yes** | Returns operator, treating edge-I as S |
| `isCloserThenVextAndGoodBase()` | L917–L930 | **yes** | Predicate: near end of M segment + adjacent I/D |
| `isInsertionOrDeletionWithNextMatched()` | L893–L901 | **yes** | Predicate: D/I M D/I pattern for multi-indels |
| `isNextMatched()` | L911–L914 | **yes** | Predicate: next CIGAR element is M |
| `isNextInsertion()` | L907–L910 | **yes** | Predicate: next CIGAR element is I |
| `isTwoInsertionsAhead()` | L903–L905 | **yes** | Predicate: element at ci+1 is I |
| `isNextAfterNumMatched()` | L886–L888 | **yes** | Predicate: element at ci+N is M |
| `parseCigarWithAmpCase()` | L1326–L1400 | **yes** | Amplicon-mode overlap/distance filtering |
| `skipOverlappingReads()` | L1895–L1907 | **yes** | Skips reads to avoid double-counting in overlapping pairs |
| `isReadsOverlap()` | L1914–L1926 | **yes** | Two-case overlap test |
| `isPairedAndSameChromosome()` | L1932–L1934 | **yes** | Paired+same-chr check |
| `isTrimAtOptTBases()` | L655–L663 | **yes** | Trim reads at -T bases |
| `skipSitesOutRegionOfInterest()` | L665–L677 | **yes** | CRISPR mode filtering |
| `skipIndelNextToIntron()` | L1694–L1700 | **yes** | Skip indels adjacent to N (intron) |
| `appendSegments()` | L1878–L1904 | **yes** | Builds complex description string for multi-indel cases |
| `increment()` | L1451–L1460 | **yes** | Increments a count in `Map<Integer, Map<String, Integer>>` |
| `isBEGIN_ATGC_AMP_ATGCs_END()` | L1462–L1477 | **yes** | Checks single-base + '&' + bases pattern (MNV) |
| `isATGC()` | L1479–L1494 | **yes** | Checks if char is A/T/G/C |
| `getAlignedLength()` | L1616–L1623 | **yes** | Sum of M+D lengths |
| `getSoftClippedLength()` | L1630–L1639 | **yes** | Sum of M+I+S lengths |
| `getMatchInsertionLength()` | L1645–L1653 | **yes** | Sum of M+I lengths |
| `getInsertionDeletionLength()` | L1659–L1667 | **yes** | Sum of I+D lengths |
| `getMateReferenceName()` | L1669–L1678 | **yes** | Returns "=", "*", or mate chr name |
| `getBaseQualityString()` | L2383–L2399 | **yes** | Safe quality string extraction with phred capping |
| `isReadChimericWithSA()` | L1841–L1873 | **yes** | Checks SA tag for chimeric reads |
| `adddisccnt()` | L2311–L2313 | **yes** | Increments `.disc` on SV Sclip |
| `getLastSVStructure()` | L2315–L2317 | **yes** | Returns last element of SV list |
| `Offset` (inner class) | L2401–L2413 | **yes** | Data holder for findOffset result |
| `CigarParserJsonlWriter` | L2415–L2662 | no | Debug JSONL writer — not parity-relevant |
| Getters (7 methods) | L187–L235 | no | Trivial accessors |

## Method Analyses

### process()
**Source**: CigarParser.java:L82–L171
**Purpose**: Entry point — iterates all BAM records, dispatches to `parseCigar()`, packages results
**Pipeline Stage**: Called by mode classes (SimpleMode, SomaticMode, AmpliconMode)
**Called By**: Pipeline framework via `Module<RecordPreprocessor, VariationData>` interface

#### Parameters
- `scope`: `Scope<RecordPreprocessor>` — contains BAM path, region, reference, preprocessor with filtered records

#### Algorithm (Step-by-Step)
1. Extract `RecordPreprocessor` from scope, call `initFromScope(scope)` to populate instance fields
2. **Record loop**: While `processor.nextRecord()` returns non-null:
   - Call `parseCigar(chrName, record)` wrapped in try/catch (exceptions print and continue)
3. Close processor
4. **If `conf.outputSplicing`**: Print splice counts to stdout and return early with empty `VariationData`
5. **Compute duprate**:
   - If `svflag`: `duprate = (discordantCount + 1) / (totalReads - duplicateReads + 1) > 0.5 ? 0.0 : 1.0`
   - Else if `removeDuplicatedReads && totalReads != 0`: `String.format("%.3f", (double)duplicateReads / totalReads)` parsed back to double
   - Else: `0.0`
6. Package all maps into `VariationData`, return wrapped in new `Scope`
7. Optional: If `VARDICT_CIGAR_PARSER_JSONL` env var set, write debug JSONL (not parity-relevant)

#### Parity Warnings
- **duprate rounding**: Java formats `%.3f` then parses back to double. Rust must replicate this round-trip to get identical floating-point bits.
- **Integer division in svflag duprate**: `(discordantCount + 1) / (totalReads - duplicateReads + 1)` is integer division, result compared `> 0.5` — always false unless denominator is 1. This means `duprate` is always `1.0` in svflag mode unless `totalReads - duplicateReads == 0`.

---

### parseCigar()
**Source**: CigarParser.java:L237–L653
**Purpose**: Core per-record CIGAR parsing — extracts variants from a single BAM record
**Called By**: `process()`

#### Parameters
- `chrName`: String — chromosome name
- `record`: SAMRecord — the BAM record to parse

#### Algorithm (Step-by-Step)
1. **Set read name** on instance: `this.currentReadName = record.getReadName()`
2. **Extract fields**: `querySequence`, `mappingQuality`, `cigar`, `ref` (from `reference.referenceSequences`)
3. **Compute insertion/deletion length**: `getInsertionDeletionLength(cigar)` sums I+D elements
4. **NM tag mismatch filter** (L251-L266):
   - Try `NM` tag, then `nM` (STAR aligner)
   - If present: `totalNumberOfMismatches = NM - insertionDeletionLength`; if > `conf.mismatch`, **return** (skip record)
   - If absent: if unmapped or `*` cigar, **return**
5. **Get quality string**: `getBaseQualityString(record)` (phred-safe)
6. **Strand direction**: `direction = record.getReadNegativeStrandFlag()` (true = reverse)
7. **Amplicon mode**: If `ampliconBasedCalling != null`, call `parseCigarWithAmpCase()` — if true, **return**
8. **Reset position trackers**: `readPositionIncludingSoftClipped = 0`, `readPositionExcludingSoftClipped = 0`
9. **CigarModifier** (L287-L299): If `performLocalRealignment`, create `CigarModifier` and apply. Otherwise use raw record position/cigar.
10. **Set currentCigarString**, call `cleanupCigar(record)`, set `start = position`, `offset = 0`
11. **Discordant count**: If mate reference != "=", increment `discordantCount`
12. **Double soft-clip filter**: If CIGAR matches `^\d\dS.*\d\dS$` (both ends >= 10bp soft-clipped), **return**
13. **Min match filter**: If `conf.minmatch != 0` and `readLengthIncludeMatchingAndInsertions < conf.minmatch`, **return**
14. **Update maxReadLength**: `totalLengthIncludingSoftClipped = getSoftClippedLength(cigar)`
15. **Supplementary alignment filter**: If not `samfilter "0"` and supplementary flag set, **return**
16. **CRISPR filter**: `skipSitesOutRegionOfInterest()`
17. **Mate direction**: `mateDirection = (flags & 0x20) != 0 ? false : true` — **note inverted logic**: flag bit 0x20 means mate is reverse, so if set, `mateDirection = false` (reverse)
18. **SV preparation**: If paired+mate unmapped -> no-op; else if `mappingQuality > 10 && !disableSV` -> `prepareSVStructuresForAnalysis()`

19. **`processCigar:` labeled loop** over CIGAR elements (L338-L649):
    - Check `skipOverlappingReads()` -> break if true
    - Get `cigarElementLength` and `operator` (via `getCigarOperator` which converts edge-I to S)
    - **Switch on operator**:
      - `N`: `processNotMatched()`, continue
      - `S`: `processSoftClip(...)`, continue
      - `H`: `offset = 0`, continue
      - `I`: `offset = 0`, `ci = processInsertion(...)`, continue
      - `D`: `offset = 0`, `ci = processDeletion(...)`, continue
      - **default (M/=/X)**: fall through to match-processing loop below

20. **Match processing** (the M-segment inner loop, L357-L640):
    - Initialize `nmoff = 0`, `moffset = 0` per CIGAR element
    - **For each base** `i` from `offset` to `cigarElementLength - 1`:
      - **Trim check**: `isTrimAtOptTBases()`
      - Get `ch1 = querySequence.charAt(readPositionIncludingSoftClipped)`, `s = String.valueOf(ch1)`
      - **Skip 'N'**: If `ch1 == 'N'`, optionally inc refCoverage if `includeNInTotalDepth`, advance positions, continue
      - Compute quality `q`, init `qbases = 1`, `qibases = 0`, `ss = new StringBuilder()`

      - **MNV detection while-loop** (L408-L464):
        - While: in region AND `i+1 < cigarElementLength` AND `q >= goodq` AND ref has mismatch at current position AND ref is not 'N':
          - Require next base quality >= `goodq + 5` or break
          - Break if next read base is 'N' or next ref is 'N'
          - If next base also mismatches reference: append to `ss`, advance all positions
          - Else **look-ahead** up to `conf.vext` bases: if another mismatch found within `vext`, grow MNV by including reference-matching bases. If no further mismatch in window, break.

      - If `ss.length() > 0`: `s += "&" + ss` (MNV encoding: first_base & remaining_bases)

      - **Complex variant detection at end of M segment** (L505-L610):
        - **Adjacent deletion case** (`isCloserThenVextAndGoodBase(..., CigarOperator.D)`):
          - Consume remaining M bases, remove `&` from s, prepend `-{delLen}&`
          - Set `startWithDeletion = true`, `ddlen = deletion length`, advance `ci`
          - If two insertions ahead: append `^{insertion_seq}`, adjust positions, advance ci
          - If next-after is M: call `findOffset()` for trailing mismatches
        - **Adjacent insertion case** (`isCloserThenVextAndGoodBase(..., CigarOperator.I)`):
          - Consume remaining M bases, remove `&`, prepend `+`, build insertion+MNV string
          - Advance ci, adjust qibases/qbases

      - **Add variant if not trimmed** (L611-L617): If `pos` (= start - qbases + 1) is in region and `s` has no 'N': `addVariationForMatchingPart(...)`
      - Advance `start`, read positions, check overlap again

    - After inner loop: if `moffset != 0`, apply it to offset/readPos/start
    - Break if `start > region.end`

#### Mutable State Modified
- `this.start` — reference position, advances through each base/operation
- `this.readPositionIncludingSoftClipped` — read index including soft clips
- `this.readPositionExcludingSoftClipped` — read index excluding soft clips
- `this.offset` — carry-over offset between CIGAR elements (for combined indel+mismatch)
- `this.cigarElementLength` — modified during soft-clip match-back (5'/3')
- `this.cigar` — may be replaced by CigarModifier
- `this.currentCigarString` — set to current cigar string
- `this.maxReadLength` — updated if current read is longer
- `this.discordantCount` — incremented for discordant mates
- `nonInsertionVariants`, `insertionVariants`, `refCoverage`, `softClips3End`, `softClips5End` — variant maps populated
- `positionToInsertionCount`, `positionToDeletionCount`, `mnp` — count maps populated
- `splice`, `spliceCount` — splice junction tracking
- `svStructures` — SV accumulators filled

#### Null/Edge Cases
- `numberOfMismatches_NM == null && numberOfMismatches_STAR == null`: Record skipped if unmapped or `*` CIGAR
- `ref.get(start)` may return `null` — all comparison methods (`isHasAndEquals`, `isHasAndNotEquals`) null-guard this
- `querySequence.charAt(idx)` — bounds never explicitly checked; assumes CIGAR arithmetic is correct
- `ci` is mutated inside the loop for D and I operators (they may consume multiple CIGAR elements)
- `offset` carries state between iterations of the main CIGAR loop — after D or I processing, the next M segment starts at `offset` instead of 0

#### Collection Ordering Dependencies
- `nonInsertionVariants` and `insertionVariants` use `VariationMap<String, Variation>` which extends `LinkedHashMap` — insertion order is preserved per position
- Keys are variant description strings (e.g., `"A"`, `"+ATG"`, `"-3"`, `"-3&AT"`, `"A&TG"`)

#### Parity Warnings
1. **`ci` mutation in switch**: `processInsertion()` and `processDeletion()` return a new `ci` value which may skip CIGAR elements. The outer loop uses `ci++` so the return value accounts for this.
2. **`offset` statefulness**: `offset` persists across the CIGAR loop iterations — it's set in D/I processing and consumed in the next M segment's inner loop (`for (int i = offset; ...)`). After the M loop, if `moffset != 0`, it's applied. This is extremely parity-critical.
3. **MNV look-ahead**: The while loop modifies `i`, `readPositionIncludingSoftClipped`, `readPositionExcludingSoftClipped`, and `start` in-place. The variable `nmoff` tracks mismatch count adjustments.
4. **String concatenation order**: `s += "&" + ss` then later `s.replaceFirst("&", "")` — this removes only the FIRST `&`. If `ss` itself contains `&`, subsequent ones survive.
5. **`mateDirection` logic inversion**: `(flags & 0x20) != 0 ? false : true` — SAM flag 0x20 is "mate reverse strand". So `mateDirection = true` means mate is FORWARD. Then `mateDirNum = mateDirection ? 1 : -1`.
6. **Quality as `double`**: `q` accumulates quality as `double` but individual quality scores are `int` (char - 33). Division `q / (qbases + qibases)` produces double — must match Java's double division exactly.
7. **`readPositionExcludingSoftClipped` for tp**: `int tp = readPositionExcludingSoftClipped < rlen1 - readPositionExcludingSoftClipped ? readPositionExcludingSoftClipped + 1 : rlen1 - readPositionExcludingSoftClipped` — min distance from either end.

---

### processDeletion()
**Source**: CigarParser.java:L672–L883
**Purpose**: Creates variant records for deletion CIGAR elements
**Called By**: `parseCigar()` main switch

#### Algorithm (Step-by-Step)
1. **Intron adjacency check**: If adjacent to N operator, skip (adjust `readPositionExcludingSoftClipped`)
2. Init description string: `"-" + cigarElementLength`
3. Get quality of last base before deletion: `queryQuality.charAt(readPositionIncludingSoftClipped - 1)`
4. **Multi-indel within VEXT** (`isInsertionOrDeletionWithNextMatched`):
   - Pattern: `D + short_M + I/D` — combines into complex variant
   - Calls `appendSegments()` to build `#matched_seq^indel_seq` string
   - Adjusts `multoffs`, `multoffp`
   - Can look further with `findOffset()` if another M follows
   - Skips 2 additional CIGAR elements
5. **D + I pattern** (`isNextInsertion`):
   - Appends `^insertion_seq` to deletion description
   - Scans trailing M for mismatches via loop
6. **D + M pattern** (`isNextMatched`):
   - Scans next M segment for nearby mismatches, up to `conf.vext` consecutive matches
   - Sets `offset` if mismatches found, appends those bases
7. **Append trailing mismatches**: If `offset > 0`, append `"&" + sequenceToAppendIfNextSegmentMatched`
8. **Build quality string**: Best of q1 (last before del) and q2 (first after del+offset)
9. **Add variant**: If `start` in region, call `addVariationForDeletion()`
10. **Adjust positions**: `start += cigarElementLength + offset + multoffs`; `readPositionIncludingSoftClipped += offset + multoffp`

#### Parity Warnings
- **`multoffs` vs `multoffp`**: `multoffs` adjusts reference position (includes deletion lengths), `multoffp` adjusts read position (includes insertion lengths). Getting these swapped breaks position tracking.
- **`isInsertionOrDeletionWithNextMatched` accesses `cigar.getCigarElement(ci + 3)`** without bounds check on ci+3 — Java may throw `IndexOutOfBoundsException` (but the caller's conditions prevent it; still, Rust must guard).
- **`qualityOfSegment` char->Phred accumulation**: Quality is built as a `StringBuilder` of raw quality characters. The final quality is computed in `addVariationForDeletion()` by iterating characters and subtracting 33. Rust must replicate this char-based quality encoding.

---

### processInsertion()
**Source**: CigarParser.java:L930–L1121
**Purpose**: Creates variant records for insertion CIGAR elements
**Called By**: `parseCigar()` main switch

#### Algorithm (Step-by-Step)
1. **Intron adjacency check**: Skip if adjacent to N
2. Extract inserted sequence from read: `substr(querySequence, readPositionIncludingSoftClipped, cigarElementLength)`
3. **Multi-indel within VEXT** — same pattern as deletion but mirror: adjusts `multoffs`/`multoffp` inversely
4. **Next matched**: Scan for mismatches in following M segment
5. Append trailing mismatches if `offset > 0`
6. **Insertion position adjustment** (L1023-L1040):
   - If start - 1 is in region and insertion is pure ATGC:
   - Call `VariationRealigner.adjInsPos()` to left-align the insertion
   - Only accept adjustment if `readPositionIncludingSoftClipped - 1 - (start - 1 - adjustedInsertionPosition) > 0`
7. **Add to counts**: `incCnt(positionToInsertionCount, insertionPosition, "+" + descString)`
8. **Create/update Variation** in `insertionVariants` at `insertionPosition`
9. **Reference count adjustment**: If the base at the adjusted insertion position matches the reference base and a non-insertion variant exists there, **subtract** the count from it (via `subCnt`)
10. **Edge insertion adjustment** (L1097-L1114): If insertion is the second CIGAR element and first is S or H, add a reference variant count at `insertionPosition` to avoid AF > 1
11. **Adjust positions**: `readPositionIncludingSoftClipped += cigarElementLength + offset + multoffp`; `start += offset + multoffs`

#### Parity Warnings
- **`adjInsPos` call**: Comes from `VariationRealigner` — cross-module dependency. The returned position and sequence may differ from original. Rust must call the same logic.
- **`subCnt` reference adjustment**: Uses `getVariationMaybe` which may return `null`. Only decrements if non-null. Critical for AF accuracy.
- **Insertion position**: `start - 1` (insertions are placed at the base BEFORE the insertion in VarDict convention)
- **`indexOfInsertionInQuerySequence`**: `readPositionIncludingSoftClipped - 1 - (start - 1 - insertionPosition)` — complex arithmetic prone to off-by-one

---

### processSoftClip()
**Source**: CigarParser.java:L1124–L1304
**Purpose**: Handles soft-clipped segments at 5' and 3' ends — chimeric detection and variant evidence
**Called By**: `parseCigar()` main switch

#### Algorithm (Step-by-Step)
1. **Branch on position**: `ci == 0` -> 5' clip; `ci == cigar.numCigarElements() - 1` -> 3' clip
2. **5' soft clip** (L1129-L1226):
   a. **Chimeric detection**: If not `conf.chimeric`, check SA tag or seed map
      - SA tag: `isReadChimericWithSA()` — if chimeric, advance positions and return
      - Seed: Reverse-complement the clipped sequence, check reference seed map for unique matches within `2 * maxReadLength`
   b. **Match-back loop**: While clipped bases match reference (from end of clip backward), create reference variant + coverage, decrement `cigarElementLength`, decrement `start`
   c. **Remaining clip quality scan**: Loop backward from clip end, count high-quality bases, stop at 'N' or second low-quality base (quality <= 12)
   d. Call `sclip5HighQualityProcessing()` to populate `softClips5End`
   e. **Restore `cigarElementLength`** to original value (it was modified by match-back)
3. **3' soft clip** (L1228-L1296): Mirror of 5' but direction reversed:
   a. Same chimeric detection logic
   b. Match-forward loop: While bases match reference going forward
   c. Forward quality scan
   d. Call `sclip3HighQualityProcessing()` to populate `softClips3End`
4. **Finalize**: `readPositionIncludingSoftClipped += cigarElementLength; offset = 0; start = position` (reset start!)

#### Parity Warnings
- **`start` is reset to `position`** at the end — this is critical. The match-back/forward loop modifies `start`, but it's reset after processing. However, during the match-back, reference variants ARE created at the adjusted positions.
- **`cigarElementLength` mutation**: Both match-back loops decrement `cigarElementLength`. For 5', it's restored from the original CIGAR element; for 3', the remaining value is used for the quality scan, then the original is used for the `readPositionIncludingSoftClipped` advance.
- **5' match-back direction**: Quality indexing uses `queryQuality.charAt(cigarElementLength - 1)` which is the LAST base of the soft clip (closest to alignment). It iterates backward from the alignment edge.
- **3' match-forward direction**: Uses `readPositionIncludingSoftClipped` which is the FIRST base of the soft clip (at the alignment edge). Iterates forward.
- **`getReverseComplementedSequence`**: For chimeric seed check, uses `substr(sequence, -SEED_1, SEED_1)` for 3' — negative index means from end. Rust `substr` equivalent must handle this.
- **`readPositionExcludingSoftClipped` increment**: In the 3' match-forward loop, `readPositionExcludingSoftClipped++` is called but NOT in the 5' match-back loop.

---

### addVariationForMatchingPart()
**Source**: CigarParser.java:L1710–L1776
**Purpose**: Creates or updates Variation records for variants found in M segments
**Called By**: Match processing inner loop in `parseCigar()`

#### Algorithm
1. If `s` starts with `"+"`: insertion — goes to `insertionVariants` + `positionToInsertionCount`
2. Else if `isBEGIN_ATGC_AMP_ATGCs_END(s)`: MNV — goes to `mnp` map AND `nonInsertionVariants`
3. Else: SNP/reference/complex -> `nonInsertionVariants`
4. Increment variant: `varsCount++`, `incDir(dir)`, compute `tp` (min of readPos from either end)
5. Compute average quality: `q = q / (qbases + qibases)`
6. Set pstd/qstd flags, accumulate mean stats
7. **Coverage**: Increase coverage for `qbases - shift` positions backward from `start` (shift = 1 if insertion with `&`)
8. If `startWithDeletion`: add to `positionToDeletionCount`, increase coverage for deletion span

#### Parity Warnings
- **Quality division**: `q / (qbases + qibases)` — integer qibases can be 0, so this is safe, but `qbases` is always >= 1
- **Coverage loop**: `for (int qi = 1; qi <= qbases - shift; qi++) { incCnt(refCoverage, start - qi + 1, 1); }` — covers positions from `start` backward. If insertion with trailing match, `shift = 1` reduces coverage count by 1.

---

### addVariationForDeletion()
**Source**: CigarParser.java:L1781–L1835
**Purpose**: Creates/updates Variation for deletion variants
**Called By**: `processDeletion()`

Similar structure to `addVariationForMatchingPart()` but deletion-specific. Creates or retrieves the Variation from `nonInsertionVariants`, increments counters, computes quality from the raw quality character buffer, and updates coverage for the deletion span.

---

### prepareSVStructuresForAnalysis()
**Source**: CigarParser.java:L2017–L2310
**Purpose**: Populates structural variant accumulators from discordant mate pairs
**Called By**: `parseCigar()` when `mappingQuality > 10 && !disableSV`

#### Algorithm
1. Compute `end` from `start + sum(MND segments)`
2. Detect `soft5` and `soft3` from CIGAR prefix/suffix S elements (quality check: last/first base >= goodq)
3. Compute `readDirNum` (-1 if reverse, 1 if forward) and `mateDirNum` (1 if forward, -1 if reverse)
4. **Intra-chromosomal** (mate ref = "="):
   - Filter: skip MC with both-end softclip, skip MQ < 15
   - **Deletion**: `readDirNum * mateDirNum == -1 && mlen * readDirNum > 0` and mlen exceeds insert size + stddev
     - Forward/reverse branches: possibly create new SV cluster if distance from last > `MINSVCDIST * maxReadLength`
     - Call `addSV()` to accumulate
     - Increment discordant counts on all nearby SV types within MIN_D (75bp)
   - **Duplication**: `readDirNum * mateDirNum == -1 && readDirNum * mlen < 0`
     - Similar clustering and `addSV()` logic
   - **Inversion**: `readDirNum * mateDirNum == 1`
     - Split into 3' vs 5' based on `mlen` sign relative to `3 * maxReadLength`
5. **Inter-chromosomal fusion**: Separate `svffus`/`svrfus` maps by mate chromosome

#### Parity Warnings
- **`mlen` recalculation**: For deletion candidates, `mlen = mateStart > start ? mend - start : end - mateStart` — overwrites the original `record.getInferredInsertSize()` value
- **`MIN_D = 75`**: Hardcoded constant, not configurable
- **SV clustering distance**: Uses `Configuration.MINSVCDIST * maxReadLength` — Rust must use same constants
- **`MINMAPBASE` index**: `queryQuality.charAt(Configuration.MINMAPBASE) - 33` as quality mean — single quality value, not an average

---

### findOffset()
**Source**: CigarParser.java:L1502–L1547
**Purpose**: Scan forward from an indel to find nearby mismatches to combine into a complex variant

#### Algorithm
1. Loop `vi` from 0, while `vsn <= conf.vext` and `vi < cigarLength`:
   - Break on 'N', low quality, null ref
   - If mismatch: `offset = vi + 1`, `tnm++`, `vsn = 0`
   - If match: `vsn++` (consecutive matches counter — breaks when > vext consecutive matches)
2. If `offset > 0`: extract substring, quality, and increment refCoverage for those positions
3. Return `Offset` object

#### Parity Warnings
- The `vsn` reset-on-mismatch logic means it finds the LAST mismatch within a window where no more than `vext` consecutive matches occur. This is NOT simply "mismatches within vext bases".

---

### cleanupCigar()
**Source**: CigarParser.java:L1553–L1599
**Purpose**: Normalize CIGAR: remove H clips, convert leading/trailing I to S

#### Algorithm
1. Copy CIGAR elements to mutable list
2. **Leading elements**: Forward iterate until first consuming-both-bases operator; convert I->S, remove H
3. **Trailing elements**: Reverse iterate same logic
4. Replace CIGAR on record

#### Parity Warnings
- This modifies `record.getCigar()` but `this.cigar` was already set from CigarModifier or raw record. The module uses `this.cigar` for the loop, not the record's cigar. So `cleanupCigar` effectively modifies the record for downstream consumers but does not affect CigarParser's own parsing.

---

### parseCigarWithAmpCase()
**Source**: CigarParser.java:L1326–L1400
**Purpose**: Amplicon-mode overlap/distance filtering

Determines whether a read should be processed based on amplicon boundaries. Returns `true` to skip the record, `false` to continue processing.

---

### isReadChimericWithSA()
**Source**: CigarParser.java:L1841–L1873
**Purpose**: Checks the SA (supplementary alignment) tag for chimeric reads

Parses SA string, extracts alignment position and chromosome, and compares with current read alignment to determine if the read is chimeric.

## Cross-Module Dependencies

### Calls (outbound)

| Module | Method | Purpose |
|--------|--------|--------|
| `RecordPreprocessor` | `nextRecord()`, `getChrName()`, `close()` | Record iteration |
| `CigarModifier` | `modifyCigar()` | CIGAR realignment (if `performLocalRealignment`) |
| `VariationRealigner` | `adjInsPos()` | Left-align insertion positions |
| `VariationUtils` | `getVariation()`, `getVariationMaybe()`, `getVariationFromSeq()`, `isHasAndEquals()`, `isHasAndNotEquals()`, `isNotEquals()`, `isEquals()` | Variation map access and base comparison |
| `Utils` | `substr()`, `getReverseComplementedSequence()`, `globalFind()`, `sum()`, `toInt()`, `incCnt()`, `getOrElse()`, `printExceptionAndContinue()` | String/collection utilities |
| `GlobalReadOnlyScope` | `instance()` | All configuration access (`.conf.*`, `.chrLengths`, `.ampliconBasedCalling`, `.sample`) |
| `Reference` | `.referenceSequences` (Map), `.seed` (seed map) | Reference genome access |
| `Patterns` | Multiple static regex patterns | CIGAR/SA tag matching |

### Called By (inbound)

| Module | Context |
|--------|--------|
| `SimpleMode`, `SomaticMode`, `AmpliconMode` | Via pipeline `Module.process()` interface |

## Data Structures Read/Written

| Structure | Type | Read/Write | Notes |
|-----------|------|------------|-------|
| `nonInsertionVariants` | `Map<Integer, VariationMap<String, Variation>>` | **Write** | SNPs, MNVs, deletions, reference bases |
| `insertionVariants` | `Map<Integer, VariationMap<String, Variation>>` | **Write** | Insertions and complex insertions |
| `refCoverage` | `Map<Integer, Integer>` | **Write** | Per-position coverage counts |
| `softClips5End` | `Map<Integer, Sclip>` | **Write** | 5' soft clip evidence |
| `softClips3End` | `Map<Integer, Sclip>` | **Write** | 3' soft clip evidence |
| `positionToInsertionCount` | `Map<Integer, Map<String, Integer>>` | **Write** | Insertion counts by position+description |
| `positionToDeletionCount` | `Map<Integer, Map<String, Integer>>` | **Write** | Deletion counts by position+description |
| `mnp` | `Map<Integer, Map<String, Integer>>` | **Write** | MNV counts by position+description |
| `spliceCount` | `Map<String, int[]>` | **Write** | Splice junction (intron) counts |
| `splice` | `Set<String>` | **Write** | Splice junction position strings |
| `svStructures` | `SVStructures` | **Write** | 12 SV accumulators (DEL/DUP/INV fwd/rev, FUS) |
| `reference.referenceSequences` | `Map<Integer, Character>` | **Read** | Reference bases |
| `reference.seed` | `Map<String, List<Integer>>` | **Read** | Reference seed map for chimeric detection |

## Known Parity Traps

### 1. Variant Description String Encoding
- `"A"` = SNP (single base)
- `"A&TG"` = MNV (first base + ampersand + remaining mismatched bases)
- `"+ATG"` = insertion
- `"+ATG&CC"` = insertion with trailing mismatches
- `"-3"` = deletion of 3 bases
- `"-3&AT"` = deletion with trailing mismatches
- `"-3^ATG"` = deletion with adjacent insertion
- `"-3#ACG^2"` = complex: deletion + matched + deletion/insertion
- `"A&TG"` detection: `isBEGIN_ATGC_AMP_ATGCs_END()` checks manually — first char ATGC, second char `&`, rest all ATGC

### 2. `offset` Statefulness Across CIGAR Elements
The `offset` variable carries across the CIGAR element loop. After `processDeletion()` or `processInsertion()`, `offset` may be non-zero, causing the next M segment to start at `offset` instead of 0. This is the #1 source of subtle position-tracking bugs.

### 3. `ci` Mutation in Switch Body
`processInsertion()` and `processDeletion()` return modified `ci` values. The for-loop applies `ci++`, so returned values must account for this.

### 4. Mate Direction Logic
`mateDirection = (flags & 0x20) != 0 ? false : true` — flag 0x20 is "mate reverse strand". So `mateDirection = true` means mate is FORWARD. Then `mateDirNum = mateDirection ? 1 : -1`.

### 5. Quality String as char[] Buffer
Quality values are stored as raw quality characters (ASCII offset by 33). Many methods pass `StringBuilder qualityOfSegment` containing these raw chars, then compute averages by iterating and subtracting 33. This pattern must be preserved exactly.

### 6. `replaceFirst("&", "")` in Complex Variant Building
Only the FIRST ampersand is removed. If the string already contains ampersands from MNV detection, the later ones survive. This affects the final description string.

### 7. `LinkedHashMap` via `VariationMap`
All `VariationMap<String, Variation>` instances maintain insertion order. This affects downstream iteration in ToVarsBuilder and output printing. Rust must use `IndexMap` or similar.

### 8. `substr()` with Negative Indices
Java `Utils.substr()` supports negative start indices (counting from string end), used in chimeric seed detection for 3' clips. Rust must replicate this behavior.

### 9. Duprate Double Round-Trip
`String.format("%.3f", (double)duplicateReads / totalReads)` parsed back to `Double.parseDouble()`. This introduces specific rounding that must be matched.

### 10. Integer Division in SV Duprate
`(discordantCount + 1) / (totalReads - duplicateReads + 1)` is integer division (both operands are `int`). Result is compared `> 0.5` which requires the result to be >= 1 (only when denominator = 1). Rust must not accidentally use float division here.

### 11. Soft-Clip Match-Back Modifies `cigarElementLength`
The 5' match-back loop decrements `cigarElementLength` and `start`. At the end, `cigarElementLength` is restored from the original CIGAR element for 5', but for 3', the remaining count is used for the quality scan. The `readPositionIncludingSoftClipped += cigarElementLength` at the end uses the REMAINING count after match-forward for 3'.

### 12. `isInsertionOrDeletionWithNextMatched` Accesses ci+3
The condition checks `cigar.getCigarElement(ci + 3).getOperator()` but only guards `cigar.numCigarElements() > ci + 2`. If the CIGAR has exactly ci+3 elements, accessing ci+3 throws `IndexOutOfBoundsException`. In practice, Java catches exceptions in the record loop and continues. Rust should guard this.

### 13. processNotMatched `splice` vs `spliceCount`
`splice` is a `Set<String>` (key is just added), while `spliceCount` is a `Map<String, int[]>` with count in `[0]`. Both use the same `"start-end"` key format but are separate data structures.

### 14. Configuration Fields Accessed
The following `GlobalReadOnlyScope.instance().conf.*` fields are used:
- `y` (debug), `outputSplicing`, `removeDuplicatedReads`, `mismatch`, `performLocalRealignment`, `chimeric`, `minmatch`, `samfilter`, `crisprCuttingSite`, `crisprFilteringBp`, `uniqueModeAlignmentEnabled`, `uniqueModeSecondInPairEnabled`, `goodq`, `vext`, `trimBasesAfter`, `includeNInTotalDepth`, `disableSV`, `INSSIZE`, `INSSTDAMT`, `INSSTD`
- `instance().ampliconBasedCalling`, `instance().chrLengths`, `instance().sample`
