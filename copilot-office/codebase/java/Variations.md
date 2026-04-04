# Variations

**Source**: `variations/` (7 files: Variation.java, Variant.java, Vars.java, Sclip.java, Mate.java, Cluster.java, VariationUtils.java) + `collection/VariationMap.java`
**LOC**: ~1,125
**Rust counterparts**: `src/variants/variants.rs` (Variant, SoftClip, Mate), `src/mods/to_vars_builder.rs` (output Variant, Vars), `src/variants/var_utils.rs` (utilities)
**Status**: complete

## Overview

The `variations` package contains the core data structures for the entire VarDict pipeline. These are cross-cutting — they flow through ALL pipeline stages from CigarParser through OutputVariant.

- **Variation** — Mutable accumulator for read-level evidence (counts, quality sums, position sums). One Variation exists per unique variant description string per genomic position. Written by CigarParser, read/modified by VariationRealigner, StructuralVariantsProcessor, and ToVarsBuilder.
- **Variant** — Final variant record for output. Built by ToVarsBuilder from Variation data, then formatted by OutputVariant printers. Contains alleles, frequencies, MSI, genotype, and all output-ready fields.
- **Vars** — Per-position container grouping a reference Variant, a list of non-reference Variants (sorted by frequency), and a description-to-Variant lookup map. Built by ToVarsBuilder.toVars().
- **Sclip** — Extends Variation with soft-clip consensus data and structural variant fields. Used by CigarParser (creates/populates), StructuralVariantsProcessor (reads for SV detection), and VariationRealigner (findconseq).
- **Mate** — Stores mate-pair information for SV detection. Created during CigarParser for discordant pairs, consumed by StructuralVariantsProcessor.
- **Cluster** — Extends Mate with a `cnt` field. Produced by checkCluster() during SV processing.
- **VariationUtils** — Static utility methods for getting/creating Variations from maps, adjusting counts (adjCnt), finding soft-clip consensus (findconseq), building reference sequences (joinRef), strand bias calculation, and various helper comparisons.
- **VariationMap** — (in `collection/` package) LinkedHashMap<K,V> with embedded `SV` struct containing SV type/pairs/splits/clusters. Critical for maintaining insertion-order iteration.

## Class Inventory

| # | Class | Extends | Java File | Lines | Role |
|---|-------|---------|-----------|-------|------|
| 1 | `Variation` | Object | Variation.java | 149 | Mutable read-evidence accumulator |
| 2 | `Variant` | Object | Variant.java | ~420 | Final variant record for output |
| 3 | `Vars` | Object | Vars.java | 31 | Per-position variant container |
| 4 | `Sclip` | **Variation** | Sclip.java | 37 | Soft-clip data + SV fields |
| 5 | `Mate` | Object | Mate.java | 44 | Mate-pair info for SV |
| 6 | `Cluster` | **Mate** | Cluster.java | 37 | Mate cluster with count |
| 7 | `VariationUtils` | Object | VariationUtils.java | 522 | Static utility methods |
| 8 | `VariationMap<K,V>` | **LinkedHashMap<K,V>** | collection/VariationMap.java | 62 | Ordered variant map + SV struct |

## Field-by-Field Analysis

### Variation (Variation.java)

The core mutable accumulator. All fields are public with Java default values (0 for int, 0.0 for double, false for boolean).

| # | Field | Type | Default | Perl Var | Writers | Readers | Parity Notes |
|---|-------|------|---------|----------|---------|---------|--------------|
| 1 | `varsCount` | int | 0 | `$cnt` | CigarParser, adjCnt, VariationRealigner | ToVarsBuilder, findconseq, all modules | Core count — Rust must use signed or handle negative from adjCnt |
| 2 | `varsCountOnForward` | int | 0 | `$dirPlus` | CigarParser (incDir), adjCnt (addDir/subDir) | ToVarsBuilder (strandBias), correctCnt | Directional — accessed via incDir/decDir/addDir/subDir/getDir |
| 3 | `varsCountOnReverse` | int | 0 | `$dirMinus` | CigarParser (incDir), adjCnt (addDir/subDir) | ToVarsBuilder (strandBias), correctCnt | Same as above |
| 4 | `meanPosition` | double | 0.0 | `$pmean` | CigarParser += position, adjCnt +=/-= | ToVarsBuilder (divide by count) | Sum, not mean — name is misleading. Divided by varsCount in ToVarsBuilder |
| 5 | `meanQuality` | double | 0.0 | `$qmean` | CigarParser += quality, adjCnt +=/-= | ToVarsBuilder (divide by count) | Sum, not mean — same pattern |
| 6 | `meanMappingQuality` | double | 0.0 | `$Qmean` | CigarParser += mapq, adjCnt +=/-= | ToVarsBuilder (divide by count) | Sum, not mean |
| 7 | `numberOfMismatches` | double | 0.0 | `$nm` | CigarParser += nm, adjCnt +=/-= | ToVarsBuilder (divide by count) | Sum, not mean |
| 8 | `lowQualityReadsCount` | int | 0 | `$locnt` | CigarParser (conditional inc), adjCnt += | ToVarsBuilder (qratio calc) | |
| 9 | `highQualityReadsCount` | int | 0 | `$hicnt` | CigarParser (conditional inc), adjCnt += | ToVarsBuilder (qratio, hicnt, hicov) | |
| 10 | `pstd` | boolean | false | `$pstd` | CigarParser (set true if position differs from pp), adjCnt (forced true) | ToVarsBuilder (output) | Once set true, never goes back to false |
| 11 | `qstd` | boolean | false | `$qstd` | CigarParser (set true if quality differs from pq), adjCnt (forced true) | ToVarsBuilder (output) | Once set true, never goes back to false |
| 12 | `pp` | int | 0 | `$pp` | CigarParser (stores read position) | CigarParser (compares with current for pstd) | Initial 0 means first read at position 0 won't trigger pstd correctly — but this is Java behavior |
| 13 | `pq` | double | 0.0 | `$pq` | CigarParser (stores base quality) | CigarParser (compares with current for qstd) | Same edge case as pp |
| 14 | `extracnt` | int | 0 | `$extracnt` | adjCnt (+= variant.varsCount) | ToVarsBuilder (extraFrequency) | Tracks count added from realignment |

**Methods on Variation:**
- `incDir(boolean dir)` — dir=true increments reverse, dir=false increments forward
- `decDir(boolean dir)` — decrements
- `getDir(boolean dir)` — returns reverse for true, forward for false
- `addDir(boolean dir, int add)` — adds to forward/reverse
- `subDir(boolean dir, int sub)` — subtracts from forward/reverse

**Parity-critical**: The dir boolean convention: `true = reverse, false = forward`. Rust counterpart uses `is_reverse: bool`.

### Variant (Variant.java)

Final variant record. Built by ToVarsBuilder.toVars(). Contains output-ready fields.

| # | Field | Type | Default | Writers | Readers | Parity Notes |
|---|-------|------|---------|---------|---------|--------------|
| 1 | `descriptionString` | String | null | ToVarsBuilder | OutputVariant, isGoodVar, adjComplex | Null until set — Rust must handle Option |
| 2 | `positionCoverage` | int | 0 | ToVarsBuilder | OutputVariant, isNoise | `$cov` |
| 3 | `varsCountOnForward` | int | 0 | ToVarsBuilder | OutputVariant, isNoise, debugVariantsContent | `$fwd` |
| 4 | `varsCountOnReverse` | int | 0 | ToVarsBuilder | OutputVariant, isNoise, debugVariantsContent | `$rev` |
| 5 | `strandBiasFlag` | String | **"0"** | ToVarsBuilder (strandBias result), SV processing | OutputVariant, isGoodVar | Default is "0" not null! Can become "2;1" compound |
| 6 | `frequency` | double | 0.0 | ToVarsBuilder | OutputVariant, isGoodVar, isNoise | `$freq` |
| 7 | `meanPosition` | double | 0.0 | ToVarsBuilder (divided) | OutputVariant, isGoodVar | Now actual mean (divided by count) |
| 8 | `isAtLeastAt2Positions` | boolean | false | ToVarsBuilder (from Variation.pstd) | OutputVariant, isNoise, debugVariantsContent | `$pstd` |
| 9 | `meanQuality` | double | 0.0 | ToVarsBuilder (divided) | OutputVariant, isGoodVar, isNoise | `$qmean` |
| 10 | `hasAtLeast2DiffQualities` | boolean | false | ToVarsBuilder (from Variation.qstd) | OutputVariant, isNoise, debugVariantsContent | `$qstd` |
| 11 | `meanMappingQuality` | double | 0.0 | ToVarsBuilder (divided) | OutputVariant, isGoodVar | `$mapq` |
| 12 | `highQualityToLowQualityRatio` | double | 0.0 | ToVarsBuilder | OutputVariant, isGoodVar, debugVariantsContent | `$qratio` |
| 13 | `highQualityReadsFrequency` | double | 0.0 | ToVarsBuilder | OutputVariant, isNoise, debugVariantsContent | `$hifreq` |
| 14 | `extraFrequency` | double | 0.0 | ToVarsBuilder | OutputVariant | `$extrafreq` |
| 15 | `shift3` | int | 0 | ToVarsBuilder | OutputVariant | |
| 16 | `msi` | double | 0.0 | ToVarsBuilder | OutputVariant, isGoodVar | MSI score |
| 17 | `msint` | int | 0 | ToVarsBuilder | OutputVariant, isGoodVar | MS unit length |
| 18 | `numberOfMismatches` | double | 0.0 | ToVarsBuilder (divided) | OutputVariant | `$nm` |
| 19 | `hicnt` | int | 0 | ToVarsBuilder | OutputVariant, isGoodVar | |
| 20 | `hicov` | int | 0 | ToVarsBuilder | OutputVariant | |
| 21 | `leftseq` | String | **""** | ToVarsBuilder, adjComplex | OutputVariant | Default empty string, NOT null |
| 22 | `rightseq` | String | **""** | ToVarsBuilder, adjComplex | OutputVariant | Default empty string, NOT null |
| 23 | `startPosition` | int | 0 | ToVarsBuilder, adjComplex | OutputVariant, isGoodVar | `$sp` |
| 24 | `endPosition` | int | 0 | ToVarsBuilder, adjComplex | OutputVariant, isGoodVar | `$ep` |
| 25 | `refReverseCoverage` | int | 0 | ToVarsBuilder | OutputVariant | `$rfc` — NOTE: field name says "reverse" but Javadoc says "forward" |
| 26 | `refForwardCoverage` | int | 0 | ToVarsBuilder | OutputVariant | `$rrc` — NOTE: field name says "forward" but Javadoc says "reverse" |
| 27 | `totalPosCoverage` | int | 0 | ToVarsBuilder, isNoise (decrements!) | OutputVariant | `$tcov` — mutated by isNoise() |
| 28 | `duprate` | double | 0.0 | ToVarsBuilder | OutputVariant | |
| 29 | `genotype` | String | null | ToVarsBuilder | OutputVariant | Null by default — must be Option |
| 30 | `varallele` | String | **""** | ToVarsBuilder, adjComplex | OutputVariant, varType, isGoodVar | Default empty string |
| 31 | `refallele` | String | **""** | ToVarsBuilder, adjComplex | OutputVariant, varType, isGoodVar | Default empty string |
| 32 | `vartype` | String | **""** | Post-processing | OutputVariant | Default empty string |
| 33 | `DEBUG` | String | **""** | Debug mode only | OutputVariant | Default empty string |
| 34 | `crispr` | int | 0 | CRISPR mode | OutputVariant | |

**Methods on Variant:**
- `isNoise()` — Checks quality/coverage thresholds; **MUTATES** the variant in-place (zeroes out coverage/counts/freq)
- `adjComplex()` — Trims common prefix/suffix from refallele/varallele, adjusts startPosition/endPosition/leftseq/rightseq
- `varType()` — Derives variant type from refallele/varallele comparison
- `isGoodVar(referenceVar, type, splice)` — Multi-criteria quality filter
- `debugVariantsContent(n)` / `debugVariantsContentSimple` / `debugVariantsContentInsertion` — Debug output formatting

**CRITICAL**: `isNoise()` mutates the variant AND returns a boolean. Side-effect-laden method that must be replicated exactly.

**CRITICAL**: Field names `refReverseCoverage` / `refForwardCoverage` have SWAPPED Javadoc vs names. The Javadoc says `$rfc` = "Reference variant forward strand coverage" but field is called `refReverseCoverage`. This is a naming bug in Java. The OUTPUT column uses `$rfc` / `$rrc` — verify which direction the values actually flow from.

### Vars (Vars.java)

Per-position variant container. Very simple but critical for ordering.

| # | Field | Type | Default | Writers | Readers | Parity Notes |
|---|-------|------|---------|---------|---------|--------------|
| 1 | `referenceVariant` | Variant | **null** | ToVarsBuilder | OutputVariant, isGoodVar, SV processing | Null until set — Rust Option |
| 2 | `variants` | List\<Variant\> | **new ArrayList<>()** | ToVarsBuilder (sorted by freq) | OutputVariant (iterate for output), getVarMaybe(VarsType.var) | Order matters — sorted descending by frequency |
| 3 | `varDescriptionStringToVariants` | Map\<String, Variant\> | **new HashMap<>()** | ToVarsBuilder | getVarMaybe(VarsType.varn), SV processing | HashMap — ordering NOT preserved. Lookup only |
| 4 | `sv` | String | **""** | SV processing | OutputVariant | Default empty string |

### Sclip (Sclip.java) — extends Variation

Inherits ALL 14 fields from Variation plus:

| # | Field | Type | Default | Writers | Readers | Parity Notes |
|---|-------|------|---------|---------|---------|--------------|
| 1-14 | (inherited from Variation) | | | | | See Variation section |
| 15 | `nt` | TreeMap\<Integer, TreeMap\<Character, Integer\>\> | **new TreeMap<>()** | CigarParser | findconseq | Sorted order (TreeMap) is critical for consensus |
| 16 | `seq` | TreeMap\<Integer, Map\<Character, Variation\>\> | **new TreeMap<>()** | CigarParser (via getVariationFromSeq) | findconseq | Outer TreeMap sorted, inner Map is HashMap. The inner Variation tracks per-base quality |
| 17 | `sequence` | String | **null** | findconseq (caches result) | findconseq (short-circuits if not null), SV processing | Null until computed — Rust Option. Once set, even to "", it's never recomputed |
| 18 | `used` | boolean | false | findconseq (if poly-A/T or low complexity), SV processing | SV processing (skips used clips) | Marks clip as already consumed |
| 19 | `start` | int | 0 | SV processing | SV processing | SV region start |
| 20 | `end` | int | 0 | SV processing | SV processing | SV region end |
| 21 | `mstart` | int | 0 | SV processing | SV processing | Mate region start |
| 22 | `mend` | int | 0 | SV processing | SV processing | Mate region end |
| 23 | `mlen` | int | 0 | SV processing | SV processing | Mate length |
| 24 | `disc` | int | 0 | SV processing | SV processing | Discordant count |
| 25 | `softp` | int | 0 | SV processing | SV processing | Soft clip position |
| 26 | `soft` | Map\<Integer, Integer\> | **new LinkedHashMap<>()** | SV processing | SV processing | **LINKED** HashMap — insertion order matters! Rust must use IndexMap |
| 27 | `mates` | List\<Mate\> | **new ArrayList<>()** | SV processing | SV processing (checkCluster) | |

### Mate (Mate.java)

| # | Field | Type | Default | Parity Notes |
|---|-------|------|---------|--------------|
| 1 | `mateStart_ms` | int | 0 | Mate pair mapping start |
| 2 | `mateEnd_me` | int | 0 | Mate pair mapping end |
| 3 | `mateLength_mlen` | int | 0 | Insert size / mate length |
| 4 | `start_s` | int | 0 | Read start |
| 5 | `end_e` | int | 0 | Read end |
| 6 | `pmean_rp` | double | 0.0 | Position in read |
| 7 | `qmean_q` | double | 0.0 | Base quality |
| 8 | `Qmean_Q` | double | 0.0 | Mapping quality — NOTE: capital Q in field name |
| 9 | `nm` | double | 0.0 | Number of mismatches |

Has two constructors: no-arg (all defaults) and all-fields.

### Cluster (Cluster.java) — extends Mate

| # | Field | Type | Default | Parity Notes |
|---|-------|------|---------|--------------|
| 1-9 | (inherited from Mate) | | | See Mate section |
| 10 | `cnt` | int | 0 (set by constructors) | Count of mates in this cluster |

Has two constructors:
1. `Cluster(cnt, ms, me, s, e)` — 5 fields, others default
2. `Cluster(ms, me, cnt, mlen, s, e, rp, q, Q, nm)` — all 10 fields. NOTE: parameter order differs from Mate constructor!

### VariationMap\<K,V\> (collection/VariationMap.java) — extends LinkedHashMap\<K,V\>

| # | Field | Type | Default | Parity Notes |
|---|-------|------|---------|--------------|
| 1 | `sv` | SV | null | Embedded SV data — Rust Option |

**Inner class SV:**

| # | Field | Type | Default |
|---|-------|------|---------|
| 1 | `type` | String | null |
| 2 | `pairs` | int | 0 |
| 3 | `splits` | int | 0 |
| 4 | `clusters` | int | 0 |

**Static methods:**
- `getSV(hash, start)` — Gets or creates SV at position. Side effect: also puts "SV" key with new Variation().
- `removeSV(hash, start)` — Sets sv=null and removes "SV" key from map.

## Method Inventory

| Method | Source | Analyzed? | Summary |
|--------|--------|-----------|---------|
| `VariationUtils.incCnt()` | VariationUtils.java:L32-L41 | yes | Increment counter in map by amount |
| `VariationUtils.strandBias()` | VariationUtils.java:L49-L63 | yes | Calculate strand bias flag (0, 1, or 2) |
| `VariationUtils.findconseq()` | VariationUtils.java:L72-L150 | yes | Find soft-clip consensus sequence |
| `VariationUtils.joinRef()` | VariationUtils.java:L159-L170 | yes | Build reference sequence from position map |
| `VariationUtils.joinRefFor5Lgins()` | VariationUtils.java:L180-L209 | yes | Reference with 5' large insertion consensus |
| `VariationUtils.joinRefFor3Lgins()` | VariationUtils.java:L211-L239 | yes | Reference with 3' large insertion consensus |
| `VariationUtils.adjCnt()` | VariationUtils.java:L253-L296 | yes | Add variant counts to target, subtract from ref |
| `VariationUtils.correctCnt()` | VariationUtils.java:L302-L315 | yes | Clamp negative fields to zero |
| `VariationUtils.getVarMaybe()` | VariationUtils.java:L343-L355 | yes | Get variant from Vars by type selector (has fall-through!) |
| `VariationUtils.getOrPutVars()` | VariationUtils.java:L363-L371 | yes | Get-or-create Vars at position |
| `VariationUtils.getVariationFromSeq()` | VariationUtils.java:L380-L395 | yes | Get-or-create Variation in Sclip.seq |
| `VariationUtils.getVariation()` | VariationUtils.java:L406-L420 | yes | Get-or-create Variation in main hash (mutable aliasing!) |
| `VariationUtils.getVariationMaybe()` | VariationUtils.java:L431-L444 | yes | Nullable lookup by position and ref base |
| `VariationUtils.isHasAndEquals()` | VariationUtils.java:L446-L480 | yes | Null-safe comparison helpers (multiple overloads) |
| `VariationUtils.isHasAndNotEquals()` | VariationUtils.java:L482-L497 | yes | Null-safe not-equal helpers |
| `VariationUtils.isEquals()` | VariationUtils.java:L499-L505 | yes | Null-safe Character equality |
| `VariationUtils.isNotEquals()` | VariationUtils.java:L507-L510 | yes | Null-safe Character inequality |
| `Variant.isNoise()` | Variant.java:L177-L195 | yes | Noise check with state mutation |
| `Variant.adjComplex()` | Variant.java:L201-L241 | yes | Trim shared prefix/suffix from alleles |
| `Variant.varType()` | Variant.java:L249-L274 | yes | Derive variant type from alleles |
| `Variant.isGoodVar()` | Variant.java:L283-L348 | yes | Multi-criteria quality filter |
| `VariationMap.getSV()` | VariationMap.java:L38-L52 | yes | Get/create SV with "SV" sentinel key |
| `VariationMap.removeSV()` | VariationMap.java:L59-L66 | yes | Remove SV data and sentinel key |

## Method Analyses

### VariationUtils.incCnt(Map counts, Object key, int add)
**Source**: VariationUtils.java:L32-L41
**Purpose**: Increment counter in a map by `add` amount.

**Algorithm:**
1. Get current value from map for key → `integer`
2. If `integer == null` → put `(key, add)`
3. Else → put `(key, integer + add)`

**Parity Notes**: Uses raw types (unchecked cast). Key can be Character or Integer. The map is always `Map<K, Integer>` in practice.

### VariationUtils.strandBias(int forwardCount, int reverseCount)
**Source**: VariationUtils.java:L49-L63
**Purpose**: Calculate strand bias flag (0, 1, or 2).

**Algorithm:**
1. If `forwardCount + reverseCount <= 12`:
   - Return `forwardCount * reverseCount > 0 ? 2 : 0`
   - (2 = both strands present, 0 = only one strand or zero)
2. Else:
   - Check if `fwd / (fwd + rev) >= conf.bias` AND `rev / (fwd + rev) >= conf.bias` AND `fwd >= conf.minBiasReads` AND `rev >= conf.minBiasReads`
   - Return 2 if all true, else 1

**Parity Notes**:
- Division is `int / (double)(int + int)` — Java auto-promotes to double
- `conf.bias` default is 0.05, `conf.minBiasReads` default is 2
- Return type is int but stored into Variant.strandBiasFlag which is String — converted by caller

### VariationUtils.findconseq(Sclip softClip, int dir)
**Source**: VariationUtils.java:L72-L150
**Purpose**: Find consensus sequence from soft-clipped reads.

**Algorithm (Step-by-Step):**
1. **Cache check**: If `softClip.sequence != null`, return it immediately (memoized)
2. Initialize: `total = 0, match = 0, seq = StringBuilder(), flag = false`
3. **Iterate** over `softClip.nt` entries (TreeMap — sorted by integer key ascending):
   a. For each position `positionInSclip`, get inner TreeMap `nv` (Character → Integer)
   b. Initialize: `maxCount = 0, maxQuality = 0, chosenBase = null, totalCount = 0`
   c. **Inner loop** over `nv` entries (also sorted — TreeMap):
      - Accumulate `totalCount += currentCount`
      - Update `maxCount`/`chosenBase` if:
        - `currentCount > maxCount` OR
        - `softClip.seq` has a Variation at this position and base with `meanQuality > maxQuality`
      - When updating, read `maxQuality` from `softClip.seq.get(positionInSclip).get(currentBase).meanQuality`
   d. **Early break at position 3**: If `positionInSclip == 3` AND `nt.size() >= 6` AND `totalCount / softClip.varsCount < 0.2` AND `totalCount <= 2` → break
   e. **Low-consensus check**: If `(totalCount - maxCount > 2 || maxCount <= totalCount - maxCount)` AND `maxCount / totalCount < 0.8`:
      - If `flag` is already true → break
      - Set `flag = true` (allow one "bad" position)
   f. Accumulate: `total += totalCount, match += maxCount`
   g. Append `chosenBase` to `seq` if not null

4. **Evaluate consensus quality**:
   - If `total != 0` AND `match / total > 0.9` AND `seq.length() / 1.5 > ntSize - seq.length()` AND (`seq.length() / ntSize > 0.8` OR `ntSize - seq.length() < 10` OR `seq.length() > 25`):
     - `SEQ = seq.toString()`
   - Else: `SEQ = ""`

5. **Poly-A/T and low-complexity check**: If `SEQ` not empty and `length > SEED_2 (12)`:
   - Match against `B_A7` or `B_T7` patterns → mark `softClip.used = true`
   - If `islowcomplexseq(SEQ)` → mark `softClip.used = true`

6. **Adaptor check**: If `SEQ` not empty and `length >= ADSEED (6)`:
   - dir=3: check first ADSEED chars against `adaptorForward` map
   - dir=5: check reversed first ADSEED chars against `adaptorReverse` map
   - If found → `SEQ = ""`

7. **Cache result**: `softClip.sequence = SEQ`
8. Return `SEQ`

**Mutable State Changes:**
- `softClip.sequence` mutated (cache set)
- `softClip.used` may be set to true

**Edge Cases:**
- `softClip.seq` may not have an entry for a position present in `softClip.nt` — the nested `.containsKey()` checks handle this
- `softClip.varsCount` being 0 would cause division by zero at step 3d — but Java integer division doesn't throw, it evaluates to infinity for doubles / NaN
- `nv` (inner TreeMap\<Character, Integer\>) iteration follows character natural order (ASCII)
- **CRITICAL**: `softClip.seq.get(positionInSclip).get(currentBase).meanQuality` — if seq contains the position but NOT the base, this would NPE. But in practice, if `nt` has a base at a position, `seq` should too.

### VariationUtils.joinRef(Map<Integer, Character> baseToPosition, int from, int to)
**Source**: VariationUtils.java:L159-L170
**Purpose**: Build reference sequence string from position-to-base map.

**Algorithm:**
1. For `i = from` to `to` (inclusive):
   - Get `ch = baseToPosition.get(i)`
   - If `ch != null` → append to StringBuilder
2. Return string

**Parity Notes:**
- Missing positions are SKIPPED (not replaced with 'N')
- Inclusive end (`i <= to`)
- Second overload: `joinRef(..., double to)` uses `i < to` (exclusive) — the double is for half-position deletions

### VariationUtils.joinRefFor5Lgins / joinRefFor3Lgins
**Source**: VariationUtils.java:L180-L239
**Purpose**: Build reference sequences incorporating consensus soft-clip sequences for 5'/3' large insertions.

**joinRefFor5Lgins**: For each position from `from` to `to`, if `to - i < seq.length() - EXTRA.length()`, take from seq at index `to - i + EXTRA.length()`; otherwise take from reference map.

**joinRefFor3Lgins**: For each position from `from` to `to`, if `i - from >= shift5` AND `i - from - shift5 < seq.length() - EXTRA.length()`, take from seq at index `i - from - shift5 + EXTRA.length()`; otherwise take from reference map.

**Parity Notes**: These use `charAt()` from Utils which handles negative indexing. Index arithmetic must match exactly.

### VariationUtils.adjCnt(Variation varToAdd, Variation variant, Variation referenceVar)
**Source**: VariationUtils.java:L253-L296
**Purpose**: Add variant's counts to varToAdd; optionally subtract from referenceVar.

**Algorithm:**
1. **Add to varToAdd:**
   - varsCount += variant.varsCount
   - extracnt += variant.varsCount (NOT variant.extracnt!)
   - highQualityReadsCount += variant.highQualityReadsCount
   - lowQualityReadsCount += variant.lowQualityReadsCount
   - meanPosition += variant.meanPosition
   - meanQuality += variant.meanQuality
   - meanMappingQuality += variant.meanMappingQuality
   - numberOfMismatches += variant.numberOfMismatches
   - **pstd = true** (forced!)
   - **qstd = true** (forced!)
   - addDir(true, variant.getDir(true)) — add reverse counts
   - addDir(false, variant.getDir(false)) — add forward counts
2. If `referenceVar == null` → return
3. **Subtract from referenceVar:**
   - varsCount -= variant.varsCount
   - highQualityReadsCount -= variant.highQualityReadsCount
   - lowQualityReadsCount -= variant.lowQualityReadsCount
   - meanPosition -= variant.meanPosition
   - meanQuality -= variant.meanQuality
   - meanMappingQuality -= variant.meanMappingQuality
   - numberOfMismatches -= variant.numberOfMismatches
   - subDir(true, variant.getDir(true))
   - subDir(false, variant.getDir(false))
4. Call `correctCnt(referenceVar)`

**CRITICAL Parity Notes:**
- `extracnt += variant.varsCount` — uses `varsCount`, NOT `extracnt`. This is intentional — extracnt tracks how many counts came from realignment.
- pstd and qstd are ALWAYS set to true after adjCnt, regardless of actual position/quality diversity
- The null check on referenceVar is the 2-arg overload behavior (calls 3-arg with null)

### VariationUtils.correctCnt(Variation varToCorrect)
**Source**: VariationUtils.java:L302-L315
**Purpose**: Clamp negative fields to zero after subtraction.

**Algorithm:**
1. If varsCount < 0 → varsCount = 0
2. If highQualityReadsCount < 0 → 0
3. If lowQualityReadsCount < 0 → 0
4. If meanPosition < 0 → 0
5. If meanQuality < 0 → 0
6. If meanMappingQuality < 0 → 0
7. If getDir(true) < 0 → addDir(true, -getDir(true)) — sets to 0 by adding negation
8. If getDir(false) < 0 → addDir(false, -getDir(false))

**Parity Notes**: The direction clamping uses addDir with negation rather than direct assignment. In Java, `addDir(true, -(-5))` = `addDir(true, 5)`, which adds 5 to the current value. So if `varsCountOnReverse = -5`, it becomes `−5 + 5 = 0`. This is equivalent to setting to 0.

### VariationUtils.getVarMaybe(Vars vars, VarsType type, Object... keys)
**Source**: VariationUtils.java:L343-L355
**Purpose**: Get variant from Vars by type selector.

**Algorithm (SWITCH WITH FALL-THROUGH!):**
```java
switch (type) {
    case var:
        if (vars.variants.size() > (Integer)keys[0]) {
            return vars.variants.get((Integer)keys[0]);
        }
    case varn:     // ← FALLS THROUGH from var!
        return vars.varDescriptionStringToVariants.get(keys[0]);
    case ref:
        return vars.referenceVariant;
}
return null;
```

**CRITICAL BUG/BEHAVIOR**: The `case var:` does NOT have a `break` or `return` after the `if` block fails. If `vars.variants.size() <= keys[0]`, execution falls through to `case varn:` which does `varDescriptionStringToVariants.get(keys[0])`. The key for `varn` is a String, but for `var`, keys[0] is an Integer. So the fallthrough calls `HashMap.get(Integer)` which will always return null (String keys never equal Integer keys).

Net effect: When called with `VarsType.var` and the index is out of bounds, it returns null (after the harmless HashMap miss on Integer key). This is effectively equivalent to having a break, but the fall-through must still be understood.

When called with `VarsType.varn`, keys[0] is a String → normal HashMap lookup.
When called with `VarsType.ref`, returns `referenceVariant` (may be null).

### VariationUtils.getOrPutVars(Map<Integer, Vars> map, int position)
**Source**: VariationUtils.java:L363-L371
**Purpose**: Get-or-create Vars at position.

**Algorithm:**
1. Get Vars from map at position
2. If null → create new Vars(), put in map
3. Return vars

**Rust equivalent**: `map.entry(position).or_insert_with(Vars::default)`

### VariationUtils.getVariationFromSeq(Sclip softClip, int idx, Character ch)
**Source**: VariationUtils.java:L380-L395
**Purpose**: Get-or-create Variation in Sclip's seq double-map.

**Algorithm:**
1. Get outer map `softClip.seq.get(idx)` → if null, create HashMap, put in seq
2. Get variation from inner map at key `ch` → if null, create new Variation(), put in inner map
3. Return variation

**Parity Notes**: Inner map is **HashMap** (not TreeMap/LinkedHashMap), but outer is **TreeMap** (sorted by position).

### VariationUtils.getVariation(hash, start, descriptionString)
**Source**: VariationUtils.java:L406-L420
**Purpose**: Get-or-create Variation in the main variant hash.

**Algorithm:**
1. Get VariationMap from hash at start position → if null, create new VariationMap<>(), put in hash
2. Get Variation from map at descriptionString → if null, create new Variation(), put in map
3. Return variation

**CRITICAL**: Returns a MUTABLE REFERENCE to the Variation object in the map. All subsequent operations modify the same object. This is the primary source of mutable aliasing — the same Variation object is visible through the hash map AND through any local variable holding the returned reference.

### VariationUtils.getVariationMaybe(hash, start, refBase)
**Source**: VariationUtils.java:L431-L444
**Purpose**: Look up Variation by position and ref base (as single-char string).

**Algorithm:**
1. If refBase is null → return null
2. Get map from hash at start → if null → return null
3. Return map.get(refBase.toString())

**Parity Notes**: Converts Character to String for lookup (e.g., 'A' → "A"). This matches the single-letter SNV description string format.

### VariationUtils.isHasAndEquals / isHasAndNotEquals (multiple overloads)
**Source**: VariationUtils.java:L446-L497
**Purpose**: Null-safe comparison helpers for reference base lookups.

All follow pattern:
1. Get value from ref map (may be null)
2. If null → return false
3. Compare with other value

**Overloads:**
- `isHasAndEquals(char, Map<Integer,Character>, int)` — char vs ref[index]
- `isHasAndEquals(int, Map<Integer,Character>, int)` — ref[index1] vs ref[index2]
- `isHasAndEquals(Map<Integer,Character>, int, String, int)` — ref[index1] vs str.charAt(index2)
- `isHasAndNotEquals(Character, Map<Integer,Character>, int)` — char != ref[index]
- `isHasAndNotEquals(Map<Integer,Character>, int, String, int)` — ref[index1] != str.charAt(index2)

### VariationUtils.isEquals / isNotEquals(Character, Character)
**Source**: VariationUtils.java:L499-L510
**Purpose**: Null-safe Character comparison.

**Algorithm**: Both null → true. One null → false. Otherwise `.equals()`.

### Variant.isNoise()
**Source**: Variant.java:L177-L195
**Purpose**: Determine if variant is noise and zero it out if so.

**Algorithm:**
1. `qual = this.meanQuality`
2. If (`qual < 4.5` OR (`qual < 12` AND `!hasAtLeast2DiffQualities`)) AND `positionCoverage <= 3`:
   - OR if `qual < conf.goodq` AND `freq < 2 * conf.lofreq` AND `positionCoverage <= 1`:
   - **MUTATE**: `totalPosCoverage -= positionCoverage; positionCoverage = 0; varsCountOnForward = 0; varsCountOnReverse = 0; frequency = 0; highQualityReadsFrequency = 0`
   - Return **true**
3. Return **false**

**CRITICAL**: The boolean check structure uses short-circuit. The two OR'd conditions are:
```
(((qual < 4.5 || (qual < 12 && !qstd)) && cov <= 3)
 || (qual < goodq && freq < 2*lofreq && cov <= 1))
```

### Variant.adjComplex()
**Source**: Variant.java:L201-L241
**Purpose**: Trim shared prefix/suffix from refallele/varallele.

**Algorithm:**
1. If varAllele starts with '<' → return (SV notation, don't trim)
2. **Prefix trim**: Count matching chars from left while both alleles have length > 1
   - Increment startPosition by n
   - Trim refallele and varallele from left by n
   - Append trimmed chars to leftseq, then trim leftseq from left by n
3. **Suffix trim**: Count matching chars from right while both alleles have length > 0 after removal
   - Decrement endPosition by (n-1)
   - Trim refallele and varallele from right by (n-1)
   - Prepend trimmed suffix to rightseq, trim rightseq from right

**Parity Notes**:
- Uses `substr()` utility (0-based, negative indices count from end)
- Leftseq manipulation: `leftseq += substr(refAllele, 0, n); leftseq = substr(leftseq, n)` — appends prefix then removes same number of chars from the start
- Suffix uses `-n` indexing with `substr(str, -n, 1)` and `substr(str, 0, 1 - n)`

### Variant.varType()
**Source**: Variant.java:L249-L274
**Purpose**: Derive variant type string from alleles.

**Algorithm:**
1. If refallele == varallele AND length == 1 → return ""
2. If both length 1 → "SNV"
3. If varallele matches ANY_SV pattern → return matched group(1)
4. If either is empty → "Complex"
5. If first chars differ → "Complex"
6. If ref length 1 AND var length > 1 AND var starts with ref → "Insertion"
7. If ref length > 1 AND var length 1 AND ref starts with var → "Deletion"
8. Default → "Complex"

### Variant.isGoodVar(referenceVar, type, splice)
**Source**: Variant.java:L283-L348
**Purpose**: Multi-criteria quality filter.

**Algorithm:**
1. If `this == null` OR `refallele` null/empty → false (NOTE: `this == null` can never be true in Java, dead code)
2. If type null/empty → compute via `varType()`
3. If freq < conf.freq OR hicnt < conf.minr OR meanPosition < conf.readPosFilter OR meanQuality < conf.goodq → false
4. If referenceVar != null AND referenceVar.hicnt > conf.minr AND freq < 0.25:
   - Compute `d = meanMappingQuality + refallele.length() + varallele.length()`
   - Compute `f = (1 + d) / (referenceVar.meanMappingQuality + 1)`
   - If `(d - 2 < 5 && referenceVar.meanMappingQuality > 20) || f < 0.25` → false
5. If type == "Deletion" AND splice contains "startPosition-endPosition" → false
6. If highQualityToLowQualityRatio < conf.qratio → false
7. If freq > 0.30 → true (skip remaining checks)
8. If meanMappingQuality < conf.mapq → false
9. MSI filters:
   - If msi >= 15 AND freq <= conf.monomerMsiFrequency AND msint == 1 → false
   - If msi >= 12 AND freq <= conf.nonMonomerMsiFrequency AND msint > 1 → false
10. If strandBiasFlag == "2;1" AND freq < 0.20:
    - If type null/SNV OR both alleles < 3 chars → false
11. Return true

### VariationMap.getSV(hash, start)
**Source**: VariationMap.java:L38-L52
**Purpose**: Get or create SV data at position.

**Algorithm:**
1. Get VariationMap from hash at start → if null, create and put
2. Get sv from map → if null:
   - Create new SV()
   - Set map.sv = sv
   - **Put ("SV", new Variation()) into map** — creates placeholder entry
3. If map doesn't contain key "SV" → put ("SV", new Variation()) — safety net
4. Return sv

**CRITICAL**: The "SV" key in the VariationMap is a REAL entry that participates in iteration. Any code iterating over the VariationMap's entries will see this "SV" → Variation entry. This must be handled in Rust.

### VariationMap.removeSV(hash, start)
**Source**: VariationMap.java:L59-L66
**Purpose**: Remove SV data at position.

**Algorithm:**
1. Set `hash.get(start).sv = null`
2. If contains key "SV" → remove it

**Parity Notes**: No null check on `hash.get(start)` — would NPE if start not in hash. Called only when known to exist.

## Cross-Module Dependencies

**Called by (readers of these data structures):**
- **CigarParser** — Creates Variation, Sclip, Mate objects; populates counts/quality sums; calls getVariation(), getVariationFromSeq(), incCnt()
- **VariationRealigner** — Reads Variation counts; calls adjCnt() to merge variants; calls findconseq() for soft-clip consensus
- **StructuralVariantsProcessor** — Reads Sclip (nt, seq, sequence, SV fields), Mate, Cluster; calls VariationMap.getSV/removeSV; reads VariationMap.sv
- **ToVarsBuilder** — Reads Variation to build Variant; creates Vars; calls strandBias(), isNoise(), isGoodVar(), adjComplex(), varType()
- **OutputVariant** — Reads Variant fields for tab-delimited output; reads Vars.sv
- **Modes (Simple/Somatic/Amplicon)** — Pass data structures between pipeline stages

**Calls (dependencies of this module):**
- **Configuration** — `conf.bias`, `conf.minBiasReads`, `conf.goodq`, `conf.freq`, `conf.lofreq`, `conf.minr`, `conf.readPosFilter`, `conf.mapq`, `conf.qratio`, `conf.monomerMsiFrequency`, `conf.nonMonomerMsiFrequency`
- **Utils** — `substr()`, `charAt()` (with negative index handling)
- **Patterns** — `B_A7`, `B_T7`, `ANY_SV`, `islowcomplexseq()`, adaptor maps

## Known Parity Traps

### Trap 1: Variation fields are SUMS not means
Fields `meanPosition`, `meanQuality`, `meanMappingQuality`, `numberOfMismatches` on Variation are SUMS. They are only divided by count in ToVarsBuilder when creating Variant. Rust field naming must not mislead developers into treating them as averages.

### Trap 2: adjCnt forces pstd/qstd to true
After `adjCnt()`, both `pstd` and `qstd` are unconditionally set to `true` on the target Variation. This means any variant that has been through realignment will always have pstd=true and qstd=true, regardless of actual position/quality diversity.

### Trap 3: extracnt tracks count from adjCnt, not from original extracnt
In `adjCnt()`, `varToAdd.extracnt += variant.varsCount` — it adds the full varsCount, not variant.extracnt. This is intentional: extracnt measures how many reads were added via realignment.

### Trap 4: correctCnt uses addDir to zero out negative direction counts
The `addDir(dir, -getDir(dir))` pattern is equivalent to setting to 0 but goes through the add path. In Rust, this can be simplified to direct assignment, but the arithmetic must produce the same result.

### Trap 5: Sclip.soft uses LinkedHashMap — insertion order matters
`Sclip.soft` is `LinkedHashMap<Integer, Integer>`. Rust must use `IndexMap` to preserve insertion order. Iteration order affects SV detection logic.

### Trap 6: Sclip.nt and Sclip.seq outer maps are TreeMap — sorted by position
Both `nt` and `seq` use `TreeMap<Integer, ...>` at the outer level. This gives ascending integer order during iteration in `findconseq()`. Rust must use `BTreeMap`.

### Trap 7: Sclip.seq inner map is HashMap, not TreeMap
While `nt` uses `TreeMap<Character, Integer>` inner maps, `seq` uses `Map<Character, Variation>` (HashMap in `getVariationFromSeq()`). The inner map ordering doesn't matter for `seq` since it's only used for lookups.

### Trap 8: Sclip.sequence caching — null vs empty string
`Sclip.sequence` starts as null (not computed). After `findconseq()`, it may be set to "" (empty — no consensus) or a non-empty string. The null check `sequence != null` is the cache guard. In Rust: `Option<Vec<u8>>` where `None` = not computed, `Some(vec![])` = computed empty.

### Trap 9: VariationMap is LinkedHashMap — iteration order = insertion order
`VariationMap<String, Variation>` extends `LinkedHashMap`. Iteration over variant description strings follows insertion order. Rust must use `IndexMap`. This affects ToVarsBuilder output ordering.

### Trap 10: "SV" sentinel key in VariationMap
`getSV()` inserts a "SV" key with an empty `Variation()` as a placeholder in the VariationMap. Code iterating over the map must either skip this key or handle it. The associated `SV` struct is a side-car field on the map itself.

### Trap 11: getVarMaybe switch fall-through
The `VarsType.var` case falls through to `VarsType.varn` when the index is out of bounds. The HashMap lookup with an Integer key against String keys returns null, making it effectively a no-op fall-through. But the code structure must be understood so Rust handles the same cases.

### Trap 12: Variant.isNoise() mutates state
`isNoise()` both returns a boolean AND zeros out coverage/count/frequency fields on the Variant. Rust translation must either replicate the mutation or restructure (but output must match).

### Trap 13: Variant field name / Javadoc mismatch (refReverseCoverage / refForwardCoverage)
Javadoc says `refReverseCoverage` is "Reference variant forward strand coverage ($rfc)" — the name and doc are SWAPPED. Verify from ToVarsBuilder which direction each actually represents before translating.

### Trap 14: Vars.varDescriptionStringToVariants is HashMap (unordered)
This is `HashMap<String, Variant>` — NOT LinkedHashMap. It's used only for lookups by description string, never iterated for output. Rust can use regular HashMap.

### Trap 15: Variant default values differ from Variation defaults
- `Variant.strandBiasFlag` defaults to `"0"` (String)
- `Variant.leftseq` defaults to `""` (empty string)
- `Variant.rightseq` defaults to `""` (empty string)
- `Variant.varallele` defaults to `""` (empty string)
- `Variant.refallele` defaults to `""` (empty string)
- `Variant.vartype` defaults to `""` (empty string)
- `Variant.DEBUG` defaults to `""` (empty string)
- `Variant.genotype` defaults to null
- `Variant.descriptionString` defaults to null
These must all be exactly replicated in Rust.

### Trap 16: Variation mutable aliasing through getVariation()
`getVariation()` returns a REFERENCE to the Variation stored in the map. Multiple callers may hold references to the same Variation object, and mutations through any reference affect all. Rust's borrow checker prevents this by default — the Rust port likely uses indices, `RefCell`, or restructuring.

### Trap 17: Mate.Qmean_Q — capital Q field name
The `Qmean_Q` field in Mate uses a capital Q, distinct from `qmean_q`. This matches the Java convention where uppercase Q refers to mapping quality and lowercase q refers to base quality.

### Trap 18: Cluster constructor parameter order differs
`Cluster(ms, me, cnt, mlen, s, e, rp, q, Q, nm)` — note `cnt` is the 3rd parameter but `ms` and `me` come first. This differs from the 5-arg constructor `Cluster(cnt, ms, me, s, e)` where `cnt` is first.

### Trap 19: Variation.pp defaults to 0 — first read at position 0 won't trigger pstd
If the first read's position in the read is 0, and pp's default is 0, then `position != pp` is false, so pstd won't trigger. The first read that sets pstd must have a position different from 0 (or a subsequent read must differ from the first). This is Java's behavior and Rust must match.

### Trap 20: findconseq early break at position 3
The position-3 early break (`positionInSclip == 3 && nt.size() >= 6 && totalCount/varsCount < 0.2 && totalCount <= 2`) — the division is `totalCount/(double)softClip.varsCount`, explicitly cast to double for floating-point division. The `< 0.2` check is against a double.

### Trap 21: findconseq `match / (double)total > 0.9` — double precision
The consensus quality metric uses double division. Edge cases near exactly 0.9 must match Java's double division behavior.

### Trap 22: Variant.adjComplex leftseq manipulation
```java
this.leftseq += substr(refAllele, 0, n);
this.leftseq = substr(this.leftseq, n);
```
This appends the trimmed prefix to leftseq, then removes `n` characters from the start of leftseq. Net effect: leftseq shifts forward by n positions in the reference, gaining the newly exposed reference bases while losing old ones.

### Trap 23: getVariation creates VariationMap (LinkedHashMap) for new positions
When a position doesn't exist in the hash, `getVariation()` creates a `VariationMap<>()` which is a `LinkedHashMap`. First entry inserted becomes first in iteration order. The order in which variants are discovered at a position determines their iteration order forever.

### Trap 24: strandBias returns int but stored as String in Variant
`VariationUtils.strandBias()` returns `int` (0, 1, or 2). But `Variant.strandBiasFlag` is `String` defaulting to `"0"`. The conversion from int to String and the compound values like `"2;1"` are done elsewhere (ToVarsBuilder). The compound `"2;1"` format means the variant-level bias was 2 but the reference-level bias was 1.

### Trap 25: Vars.sv field defaults to empty string ""
Unlike most String fields that default to null, `Vars.sv` defaults to `""`. This affects null vs empty checks downstream.
