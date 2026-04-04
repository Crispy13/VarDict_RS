# ToVarsBuilder

**Source**: `modules/ToVarsBuilder.java`
**LOC**: 1,194
**Rust counterpart**: `src/variants/tovars_builder.rs`
**Status**: complete

## Overview

`ToVarsBuilder` is the pipeline stage that transforms **raw `Variation` maps** (intermediate per-position counters accumulated by CigarParser and refined by VariationRealigner/StructuralVariantsProcessor) into **final `Variant` objects** suitable for output. It implements `Module<RealignedVariationData, AlignedVarsData>`.

**Key responsibilities**:
1. Iterates over all non-insertion variant positions from `nonInsertionVariants`
2. Skips positions outside the region, without coverage, or with only reference matching the ref base
3. Computes per-variant statistics: allele frequency, strand bias, mean quality, mean position, mean mapping quality, hi/lo quality ratio, extra frequency, number of mismatches
4. Computes MSI (Microsatellite Instability) for indels and SNPs
5. Determines ref/var alleles, genotype strings, and left/right flanking sequences
6. Handles CRISPR cut-site adjustment
7. Validates reference alleles for IUPAC ambiguity codes
8. Organizes variants into `Vars` structures keyed by position
9. Assigns reference variant vs. non-reference variants and filters by frequency threshold
10. Produces `AlignedVarsData` containing `Map<Integer, Vars> alignedVariants`

**Input**: `Scope<RealignedVariationData>` containing:
- `nonInsertionVariants`: `Map<Integer, VariationMap<String, Variation>>` — keyed by position
- `insertionVariants`: `Map<Integer, VariationMap<String, Variation>>` — keyed by position
- `refCoverage`: `Map<Integer, Integer>` — total coverage per position
- `duprate`: `double`
- Reference sequence: `Map<Integer, Character>` from `scope.regionRef.referenceSequences`
- `region`: `Region` with `.chr`, `.start`, `.end`

**Output**: `Scope<AlignedVarsData>` containing:
- `maxReadLength`: `int`
- `alignedVariants`: `Map<Integer, Vars>` — position → `Vars` (referenceVariant + variants list + varDescriptionStringToVariants map + sv string)

## Method Inventory

| Method | Lines | Analyzed? | Summary |
|--------|-------|-----------|---------|
| `process()` | L88–L190 | yes | Main entry point — loops positions, calls createVariant/createInsertion, collects and filters |
| `initFromScope()` | L76–L83 | yes | Extracts fields from scope into instance variables |
| `isTheSameVariationOnRef()` | L199–L213 | yes | Skip positions with only reference variant (no pileup/somatic/amplicon) |
| `proceedVrefIsDeletion()` | L221–L247 | yes | Compute MSI/shift3 for deletion variants |
| `proceedVrefIsInsertion()` | L254–L285 | yes | Compute MSI/shift3 for insertion variants |
| `collectVarsAtPosition()` | L293–L314 | yes | Distribute variants into Vars.referenceVariant vs Vars.variants; return maxfreq |
| `sortVariants()` | L320–L331 | yes | Sort variants by quality*coverage descending, then by descriptionString ascending |
| `createInsertion()` | L340–L424 | yes | Convert insertion Variations into Variant objects |
| `createVariant()` | L436–L510 | yes | Convert non-insertion Variations into Variant objects |
| `adjustVariantCounts()` | L517–L540 | yes | Clamp negative fwd/rev/rfc/rrc to zero with stderr warning |
| `calcHicov()` | L548–L564 | yes | Sum highQualityReadsCount across non-insertion variants at position |
| `findMSI()` | L573–L621 | yes | Microsatellite instability detection algorithm |
| `collectReferenceVariants()` | L629–L969 | yes | Core: genotype assignment, ref/var allele construction, MSI, CRISPR, flanking sequences |
| `updateRefVariant()` | L971–L1013 | yes | Fill reference-only variant fields when no non-ref variants exist |
| `constructDebugLines()` | L1018–L1030 | yes | Build DEBUG string from debug lines list |
| `validateRefallele()` | L1039–L1048 | yes | Replace IUPAC ambiguity codes in reference allele |
| Inner class `MSI` | L1053–L1073 | yes | Simple tuple: msi, shift3, msint |
| Inner class `ToVarsJsonlWriter` | L1075–L1194 | yes | JSONL debug output writer (diagnostic, not parity-critical) |
| `getInsertionVariants()` | L66–L68 | yes | Getter for insertionVariants field |
| `getNonInsertionVariants()` | L70–L74 | yes | Getter for nonInsertionVariants field |

## Method Analyses

### `initFromScope(Scope<RealignedVariationData> scope)` — L76–L83

**Purpose**: Extract fields from scope into instance variables for use throughout the class.

**Algorithm**:
1. Set `this.ref` = `scope.regionRef.referenceSequences` (type: `Map<Integer, Character>`)
2. Set `this.region` = `scope.region` (type: `Region`)
3. Set `this.refCoverage` = `scope.data.refCoverage` (type: `Map<Integer, Integer>`)
4. Set `this.insertionVariants` = `scope.data.insertionVariants` (type: `Map<Integer, VariationMap<String, Variation>>`)
5. Set `this.nonInsertionVariants` = `scope.data.nonInsertionVariants` (type: `Map<Integer, VariationMap<String, Variation>>`)
6. Set `this.duprate` = `scope.data.duprate` (type: `double`)

**Null/Edge Cases**: None — all fields come from RealignedVariationData constructor, always populated.

---

### `process(Scope<RealignedVariationData> scope)` — L88–L190

**Purpose**: Main entry point. Iterates over all positions in nonInsertionVariants, builds Variant objects, collects them into alignedVariants map.

**Returns**: `Scope<AlignedVarsData>`

**Algorithm**:
1. Call `initFromScope(scope)`
2. Read `config` = `instance().conf`
3. If `config.y`, print debug segment info to stderr
4. Create `alignedVariants` = new `HashMap<Integer, Vars>()`
5. Set `lastPosition` = 0
6. **Loop** over `getNonInsertionVariants().entrySet()` as `entH`:
   - 6a. Set `position` = `entH.getKey()`
   - 6b. Set `lastPosition` = `position`
   - 6c. Set `varsAtCurPosition` = `entH.getValue()` (type: `VariationMap<String, Variation>`)
   - 6d. **Skip** if `varsAtCurPosition.isEmpty()` AND `!getInsertionVariants().containsKey(position)`
   - 6e. **Skip** if `varsAtCurPosition.sv == null` OR `conf.deleteDuplicateVariants`:
     - If `position < region.start` OR `position > region.end` → `continue`
   - 6f. **Skip** if `varsAtCurPosition.sv == null` AND `!refCoverage.containsKey(position)` → `continue`
   - 6g. If `isTheSameVariationOnRef(position, varsAtCurPosition)` → `continue`
   - 6h. If `!refCoverage.containsKey(position)` OR `refCoverage.get(position) == 0`:
     - Print error to stderr with the SV type
     - `continue`
   - 6i. Set `totalPosCoverage` = `refCoverage.get(position)`
   - 6j. Compute `hicov` = `calcHicov(getInsertionVariants().get(position), varsAtCurPosition)`
   - 6k. Create `var` = new `ArrayList<Variant>()`
   - 6l. Create `keys` = new `ArrayList<>(varsAtCurPosition.keySet())`; `Collections.sort(keys)`
   - 6m. Create `debugLines` = new `ArrayList<String>()`
   - 6n. Call `createVariant(duprate, alignedVariants, position, varsAtCurPosition, totalPosCoverage, var, debugLines, keys, hicov)`
   - 6o. Set `totalPosCoverage` = `createInsertion(duprate, position, totalPosCoverage, var, debugLines, hicov)`
   - 6p. Call `sortVariants(var)`
   - 6q. Set `maxfreq` = `collectVarsAtPosition(alignedVariants, position, var)`
   - 6r. If `!config.doPileup` AND `maxfreq <= config.freq` AND `instance().ampliconBasedCalling == null`:
     - If `!config.bam.hasBam2()` → `alignedVariants.remove(position)`; `continue`
   - 6s. Set `variationsAtPos` = `getOrPutVars(alignedVariants, position)`
   - 6t. Call `collectReferenceVariants(position, totalPosCoverage, variationsAtPos, debugLines)`
7. Catch `Exception` per position: call `printExceptionAndContinue(...)`
8. If `config.y`, print timestamp to stderr
9. Create `alignedData` = new `AlignedVarsData(scope.maxReadLength, alignedVariants)`
10. Check `VARDICT_TO_VARS_JSONL` env var; if set, call `ToVarsJsonlWriter.write(...)`
11. Return new `Scope<>(scope, alignedData)`

**Mutable State Changes**:
- `alignedVariants` map — populated throughout the loop
- `refCoverage` may be mutated indirectly via `getOrElse` in `proceedVrefIsDeletion`/`proceedVrefIsInsertion` (getOrElse puts default value if key missing)
- `nonInsertionVariants` may be mutated in `createInsertion` (adjusts `varsCountOnForward`/`varsCountOnReverse` on next-position reference variation)

**Null/Edge Cases**:
- `getInsertionVariants().get(position)` may return `null` — passed to `calcHicov` which checks for null
- `refCoverage.get(position)` may be null — checked explicitly
- `varsAtCurPosition.sv` may be null — checked before access
- `instance().ampliconBasedCalling` may be null
- `config.bam.hasBam2()` — checks if 2nd BAM configured
- `Exception` caught per position iteration → error printed, loop continues

**Parity-Critical**:
- The iteration order of `getNonInsertionVariants().entrySet()` is `LinkedHashMap` (insertion order) — Rust must use `IndexMap` or BTreeMap
- `keys` sorted with `Collections.sort(keys)` — natural String ordering (lexicographic)
- `totalPosCoverage` can be updated by `createInsertion` and the updated value is used in `collectReferenceVariants`

---

### `isTheSameVariationOnRef(int position, VariationMap<String, Variation> varsAtCurPosition)` — L199–L213

**Purpose**: Determine if position should be skipped because it only has reference variant and no pileup/amplicon/somatic mode.

**Algorithm**:
1. Create `vk` = new `HashSet<>(varsAtCurPosition.keySet())` — copy of keys
2. If `getInsertionVariants().containsKey(position)` → add `"I"` to `vk`
3. If `vk.size() == 1` AND `ref.containsKey(position)` AND `vk.contains(ref.get(position).toString())`:
   - If `!instance().conf.doPileup` AND `!instance().conf.bam.hasBam2()` AND `instance().ampliconBasedCalling == null`:
     - Return `true` (skip this position)
4. Return `false`

**Parity-Critical**:
- The check `ref.get(position).toString()` converts `Character` to `String` — if ref doesn't contain position, this could NPE (but step 3 requires `ref.containsKey(position)`)
- Adding `"I"` as a synthetic key — never matches a single-character reference base

---

### `proceedVrefIsDeletion(int position, int dellen)` — L221–L247

**Purpose**: Compute MSI, shift3, and msint for a deletion variant.

**Algorithm**:
1. Build `leftseq` = `joinRef(ref, max(position - 70, 1), position - 1)` — left 70 bases
2. Get `chr0` = `getOrElse(instance().chrLengths, region.chr, 0)` — chromosome length (defaults to 0 if missing; **and inserts 0 into map**)
3. Build `tseq` = `joinRef(ref, position, min(position + dellen + 70, chr0))` — right dellen+70 bases
4. Call `findMSI(substr(tseq, 0, dellen), substr(tseq, dellen), leftseq)`:
   - First arg: first `dellen` chars of tseq (the deleted bases)
   - Second arg: remaining chars after deletion
   - Third arg: left context
5. Extract `msi`, `shift3`, `msint` from result
6. Call `findMSI(leftseq, substr(tseq, dellen), null)`:
   - Same right context, but with left context as first arg and no third arg
7. Extract `tmsi`, `tmsint`
8. If `msi < tmsi`:
   - Set `msi` = `tmsi`
   - Set `msint` = `tmsint`
   - **Do NOT change shift3**
9. If `msi <= shift3 / (double)dellen`:
   - Set `msi` = `shift3 / (double)dellen`
10. Return new `MSI(msi, shift3, msint)`

**Parity-Critical**:
- `getOrElse` **mutates** `instance().chrLengths` by inserting default 0 if key missing
- `substr(tseq, 0, dellen)` — if `tseq.length() < dellen`, returns what's available
- `substr(tseq, dellen)` — returns empty string if `tseq.length() <= dellen`
- Division `shift3 / (double)dellen` — double division, no integer truncation

---

### `proceedVrefIsInsertion(int position, String vn)` — L254–L285

**Purpose**: Compute MSI, shift3, and msint for an insertion variant.

**Algorithm**:
1. Extract `tseq1` = `vn.substring(1)` — insertion sequence without leading `+`
2. Build `leftseq` = `joinRef(ref, max(position - 50, 1), position)` — left 50 bases (inclusive of position)
3. Get `x` = `getOrElse(instance().chrLengths, region.chr, 0)` — chromosome length
4. Build `tseq2` = `joinRef(ref, position + 1, min(position + 70, x))` — right 70 bases
5. Call `findMSI(tseq1, tseq2, leftseq)` → extract `msi`, `shift3`, `msint`
6. Call `findMSI(leftseq, tseq2, null)` → extract `tmsi`, `tmsint`
7. If `msi < tmsi`:
   - Set `msi` = `tmsi`, `msint` = `tmsint`
   - **Do NOT change shift3**
8. If `msi <= shift3 / (double)tseq1.length()`:
   - Set `msi` = `shift3 / (double)tseq1.length()`
9. Return new `MSI(msi, shift3, msint)`

**Parity-Critical**:
- `leftseq` for insertion includes position itself (`position - 50` to `position`), but for deletion it's `position - 70` to `position - 1` — different ranges
- `vn.substring(1)` vs `descriptionString.substring(1)` — same thing, just different naming

---

### `collectVarsAtPosition(Map<Integer, Vars> alignedVariants, int position, List<Variant> var)` — L293–L314

**Purpose**: Distribute variants into reference vs. non-reference buckets in the Vars structure.

**Algorithm**:
1. Set `maxfreq` = 0
2. **Loop** over `var` as `tvar`:
   - 2a. If `tvar.descriptionString.equals(String.valueOf(ref.get(position)))`:
     - This is the reference variant → set `getOrPutVars(alignedVariants, position).referenceVariant` = `tvar`
   - 2b. Else:
     - Append `tvar` to `getOrPutVars(alignedVariants, position).variants` list
     - Put `tvar` into `getOrPutVars(alignedVariants, position).varDescriptionStringToVariants` map with key `tvar.descriptionString`
     - If `tvar.frequency > maxfreq` → `maxfreq` = `tvar.frequency`
3. Return `maxfreq`

**Parity-Critical**:
- `ref.get(position)` could return `null` if reference doesn't have this position — `String.valueOf(null)` returns `"null"` in Java, which would never match a variant description string. This is a subtle null safety pattern.
- `getOrPutVars` is called multiple times per loop iteration — creates Vars on first call, reuses on subsequent

---

### `sortVariants(List<Variant> var)` — L320–L331

**Purpose**: Sort variant list in-place by quality*coverage descending, with tiebreak by descriptionString ascending.

**Algorithm**:
1. Call `Collections.sort(var, ...)` with comparator:
   - Primary: `Double.compare(o2.meanQuality * o2.positionCoverage, o1.meanQuality * o1.positionCoverage)` — **descending**
   - Secondary (if primary == 0): `o1.descriptionString.compareTo(o2.descriptionString)` — **ascending** (lexicographic)

**Parity-Critical**:
- `Double.compare` handles NaN and -0.0 correctly
- The sort is stable (`Collections.sort` uses TimSort) — equal elements preserve insertion order
- Product `meanQuality * positionCoverage` is computed as `double * int` → double

---

### `createInsertion(double duprate, int position, int totalPosCoverage, List<Variant> var, List<String> debugLines, int hicov)` — L340–L424

**Purpose**: Convert insertion Variations at this position into Variant objects and add to var list. Also update totalPosCoverage.

**Returns**: Updated `totalPosCoverage`

**Algorithm**:
1. Get `insertionVariations` = `getInsertionVariants().get(position)` (may be null)
2. If `insertionVariations != null`:
   - 2a. Create `insertionDescriptionStrings` = new `ArrayList<>(insertionVariations.keySet())`
   - 2b. `Collections.sort(insertionDescriptionStrings)` — lexicographic sort
   - 2c. **Loop** over sorted `insertionDescriptionStrings` as `descriptionString`:
     - 2c.i. If `descriptionString.contains("&")` AND `refCoverage.containsKey(position + 1)`:
       - Set `totalPosCoverage` = `refCoverage.get(position + 1)`
     - 2c.ii. Get `cnt` = `insertionVariations.get(descriptionString)` (type: `Variation`)
     - 2c.iii. Compute `fwd` = `cnt.getDir(false)` (= `cnt.varsCountOnForward`)
     - 2c.iv. Compute `rev` = `cnt.getDir(true)` (= `cnt.varsCountOnReverse`)
     - 2c.v. Compute `bias` = `strandBias(fwd, rev)` → 0, 1, or 2
     - 2c.vi. Compute `vqual` = `roundHalfEven("0.0", cnt.meanQuality / cnt.varsCount)` — base quality
     - 2c.vii. Compute `mq` = `roundHalfEven("0.0", cnt.meanMappingQuality / (double)cnt.varsCount)` — mapping quality
     - 2c.viii. Set `hicnt` = `cnt.highQualityReadsCount`
     - 2c.ix. Set `locnt` = `cnt.lowQualityReadsCount`
     - 2c.x. Set `ttcov` = `totalPosCoverage`
     - 2c.xi. **Adjust ttcov up**: If `cnt.varsCount > totalPosCoverage` AND `cnt.extracnt != 0` AND `cnt.varsCount - totalPosCoverage < cnt.extracnt`:
       - Set `ttcov` = `cnt.varsCount`
     - 2c.xii. If `ttcov < cnt.varsCount`:
       - Set `ttcov` = `cnt.varsCount`
       - If `refCoverage.containsKey(position + 1)` AND `ttcov < refCoverage.get(position + 1) - cnt.varsCount`:
         - Set `ttcov` = `refCoverage.get(position + 1)`
         - Get `variantNextPosition` = `getVariationMaybe(getNonInsertionVariants(), position + 1, ref.get(position + 1))`
         - If `variantNextPosition != null`:
           - `variantNextPosition.varsCountOnForward -= fwd`
           - `variantNextPosition.varsCountOnReverse -= rev`
       - Set `totalPosCoverage` = `ttcov`
     - 2c.xiii. Create new `Variant tvref`
     - 2c.xiv. If `hicov < hicnt` → `hicov` = `hicnt`
     - 2c.xv. Set all Variant fields:
       - `descriptionString` = descriptionString
       - `positionCoverage` = `cnt.varsCount`
       - `varsCountOnForward` = `fwd`
       - `varsCountOnReverse` = `rev`
       - `strandBiasFlag` = `String.valueOf(bias)`
       - `frequency` = `roundHalfEven("0.0000", cnt.varsCount / (double)ttcov)`
       - `meanPosition` = `roundHalfEven("0.0", cnt.meanPosition / (double)cnt.varsCount)`
       - `isAtLeastAt2Positions` = `cnt.pstd`
       - `meanQuality` = `vqual`
       - `hasAtLeast2DiffQualities` = `cnt.qstd`
       - `meanMappingQuality` = `mq`
       - `highQualityToLowQualityRatio` = `hicnt / (locnt != 0 ? locnt : 0.5d)`
       - `highQualityReadsFrequency` = `hicov > 0 ? hicnt / (double)hicov : 0`
       - `extraFrequency` = `cnt.extracnt != 0 ? cnt.extracnt / (double)ttcov : 0`
       - `shift3` = 0
       - `msi` = 0
       - `numberOfMismatches` = `roundHalfEven("0.0", cnt.numberOfMismatches / (double)cnt.varsCount)`
       - `hicnt` = `hicnt`
       - `hicov` = `hicov`
       - `duprate` = `duprate`
     - 2c.xvi. `var.add(tvref)`
     - 2c.xvii. If debug mode, call `tvref.debugVariantsContentInsertion(debugLines, descriptionString)`
3. Return `totalPosCoverage`

**Mutable State Changes**:
- `totalPosCoverage` — may be increased
- `hicov` — may be increased (local copy, but affects subsequent iterations within this loop)
- `variantNextPosition.varsCountOnForward/varsCountOnReverse` — decreased (mutates the Variation in the nonInsertionVariants map at position+1)

**Null/Edge Cases**:
- `insertionVariations` can be null → whole block skipped
- `cnt.varsCount` could be 0 → division by zero in quality calculations produces NaN/Infinity, but note `cnt.varsCount` check is NOT done here (unlike `createVariant`)
- `ref.get(position + 1)` could return null → `getVariationMaybe` handles this
- `hicov` is a local parameter that gets modified — the original `hicov` from the caller is NOT updated (Java passes int by value)

**Parity-Critical**:
- `insertionVariations.keySet()` is from `VariationMap` which extends `LinkedHashMap` — iteration order = insertion order. But then sorted, so insertion order doesn't matter here.
- `highQualityToLowQualityRatio`: when `locnt == 0`, divisor becomes `0.5d`, not `0`. This is a common parity trap.
- `cnt.extracnt != 0` check in line 2c.xi uses `!=` to test for exactly 0 — note in Java `int` comparison
- `hicov` local mutation across loop iterations — the cumulative update within the loop body (`if hicov < hicnt → hicov = hicnt`) carries over to the next iteration. This is parity-critical because Rust could accidentally shadow or reset this variable.

---

### `createVariant(double duprate, Map<Integer, Vars> alignedVars, int position, VariationMap<String, Variation> nonInsertionVariations, int totalPosCoverage, List<Variant> var, List<String> debugLines, List<String> keys, int hicov)` — L436–L510

**Purpose**: Convert non-insertion Variations at this position into Variant objects.

**Algorithm**:
1. **Loop** over `keys` as `descriptionString`:
   - 1a. If `descriptionString.equals("SV")`:
     - Get `sv` = `nonInsertionVariations.sv`
     - Set `getOrPutVars(alignedVars, position).sv` = `sv.splits + "-" + sv.pairs + "-" + sv.clusters`
     - `continue`
   - 1b. Get `cnt` = `nonInsertionVariations.get(descriptionString)` (type: `Variation`)
   - 1c. If `cnt.varsCount == 0` → `continue` (skip zero-count variants)
   - 1d. Compute `fwd` = `cnt.getDir(false)`
   - 1e. Compute `rev` = `cnt.getDir(true)`
   - 1f. Compute `bias` = `strandBias(fwd, rev)`
   - 1g. Compute `baseQuality` = `roundHalfEven("0.0", cnt.meanQuality / cnt.varsCount)`
   - 1h. Compute `mappingQuality` = `roundHalfEven("0.0", cnt.meanMappingQuality / (double)cnt.varsCount)`
   - 1i. Set `hicnt` = `cnt.highQualityReadsCount`
   - 1j. Set `locnt` = `cnt.lowQualityReadsCount`
   - 1k. Set `ttcov` = `totalPosCoverage`
   - 1l. **Adjust ttcov**: If `cnt.varsCount > totalPosCoverage` AND `cnt.extracnt > 0` AND `cnt.varsCount - totalPosCoverage < cnt.extracnt`:
     - Set `ttcov` = `cnt.varsCount`
   - 1m. Create new `Variant tvref` and set fields:
     - `descriptionString` = descriptionString
     - `positionCoverage` = `cnt.varsCount`
     - `varsCountOnForward` = `fwd`
     - `varsCountOnReverse` = `rev`
     - `strandBiasFlag` = `String.valueOf(bias)`
     - `frequency` = `roundHalfEven("0.0000", cnt.varsCount / (double)ttcov)`
     - `meanPosition` = `roundHalfEven("0.0", cnt.meanPosition / (double)cnt.varsCount)`
     - `isAtLeastAt2Positions` = `cnt.pstd`
     - `meanQuality` = `baseQuality`
     - `hasAtLeast2DiffQualities` = `cnt.qstd`
     - `meanMappingQuality` = `mappingQuality`
     - `highQualityToLowQualityRatio` = `hicnt / (locnt != 0 ? locnt : 0.5d)`
     - `highQualityReadsFrequency` = `hicov > 0 ? hicnt / (double)hicov : 0`
     - `extraFrequency` = `cnt.extracnt != 0 ? cnt.extracnt / (double)ttcov : 0`
     - `shift3` = 0
     - `msi` = 0
     - `numberOfMismatches` = `roundHalfEven("0.0", cnt.numberOfMismatches / (double)cnt.varsCount)`
     - `hicnt` = `hicnt`
     - `hicov` = `hicov`
     - `duprate` = `duprate`
   - 1n. `var.add(tvref)`
   - 1o. If debug, call `tvref.debugVariantsContentSimple(debugLines, descriptionString)`

**Key Difference from `createInsertion`**:
- No `totalPosCoverage` update/return — non-insertion doesn't modify it
- `extracnt > 0` check (not `!= 0`) — different from insertion's `!= 0`
- No `hicov` mutation
- `cnt.varsCount == 0` skip — insertions don't have this check
- SV handling: the `"SV"` key gets special treatment — formats as `splits-pairs-clusters`

**Parity-Critical**:
- `cnt.meanQuality / cnt.varsCount` — `meanQuality` is `double`, `varsCount` is `int`, result is double. Same pattern as `meanPosition` and `numberOfMismatches`.
- `cnt.meanMappingQuality / (double)cnt.varsCount` — explicit cast to double for mapping quality
- `extracnt > 0` vs `extracnt != 0` difference between createVariant and createInsertion

---

### `adjustVariantCounts(int p, Variant vref)` — L517–L540

**Purpose**: Clamp negative forward/reverse counts to zero with stderr warning.

**Algorithm**:
1. Build `message` string from position and ref→var alleles
2. If `vref.refForwardCoverage < 0` → set to 0, print warning
3. If `vref.refReverseCoverage < 0` → set to 0, print warning
4. If `vref.varsCountOnForward < 0` → set to 0, print warning
5. If `vref.varsCountOnReverse < 0` → set to 0, print warning

---

### `calcHicov(VariationMap<String, Variation> insertionVariations, VariationMap<String, Variation> nonInsertionVariations)` — L548–L564

**Purpose**: Sum high-quality read counts across all non-insertion variants at a position.

**Algorithm**:
1. Set `hicov` = 0
2. **Loop** over `nonInsertionVariations.entrySet()` as `descVariantEntry`:
   - 2a. If key equals `"SV"` OR key starts with `"+"` → `continue`
   - 2b. `hicov += descVariantEntry.getValue().highQualityReadsCount`
3. If `insertionVariations != null`:
   - **Loop** over `insertionVariations.values()`:
     - **COMMENTED OUT**: `//hicov += variation.highQualityReadsCount;`
4. Return `hicov`

**Parity-Critical**:
- Insertion high-quality counts are **intentionally NOT added** to hicov (commented out in Java)
- Keys starting with `"+"` are skipped — these are insertion description strings that somehow ended up in nonInsertionVariations
- Iteration order of `nonInsertionVariations.entrySet()` is LinkedHashMap order, but since we're summing, order doesn't affect result

---

### `findMSI(String tseq1, String tseq2, String left)` — L573–L621

**Purpose**: Find microsatellite instability by testing repeat patterns of length 1–6 bases.

**Parameters**:
- `tseq1`: variant/deleted sequence (or left context in second call)
- `tseq2`: right reference context
- `left`: left reference context (may be null)

**Returns**: `MSI` (msi count, shift3, msint unit)

**Algorithm**:
1. Set `nmsi` = 1, `shift3` = 0, `maxmsi` = `""`, `msicnt` = 0
2. **While** `nmsi <= tseq1.length()` AND `nmsi <= 6`:
   - 2a. Extract `msint` = `substr(tseq1, -nmsi, nmsi)` — last `nmsi` chars of tseq1
   - 2b. Compile regex `pattern` = `"((" + msint + ")+)$"` — match repeat at end
   - 2c. Match `pattern` against `tseq1`
   - 2d. Set `msimatch` = matched group(1) if found, else `""`
   - 2e. If `left != null && !left.isEmpty()`:
     - Match same pattern against `left + tseq1` (concatenated)
     - If found, update `msimatch` = matched group(1) — **overwrites** previous match
   - 2f. Compute `curmsi` = `msimatch.length() / (double)nmsi`
   - 2g. Match `"^((" + msint + ")+)"` against `tseq2` — repeat at start of right context
   - 2h. If found: `curmsi += mtch.group(1).length() / (double)nmsi`
   - 2i. If `curmsi > msicnt`:
     - `maxmsi` = `msint`, `msicnt` = `curmsi`
   - 2j. `nmsi++`
3. Build `tseq` = `tseq1 + tseq2`
4. Set `shift3` = 0
5. **While** `shift3 < tseq2.length()` AND `tseq.charAt(shift3) == tseq2.charAt(shift3)`:
   - `shift3++`
6. Return new `MSI(msicnt, shift3, maxmsi)`

**Parity-Critical**:
- `substr(tseq1, -nmsi, nmsi)` extracts last `nmsi` chars. If `tseq1.length() < nmsi`, the behavior depends on `substr` clamping.
- The regex is compiled fresh each iteration — no caching
- In step 2e, when `left` is provided, the match against `left + tseq1` **replaces** the match from step 2d entirely (doesn't combine). The longer concatenated string might produce a longer match.
- `shift3` computation at step 4-5: compares `tseq = tseq1 + tseq2` character by character with `tseq2`. Since `tseq[0] = tseq1[0]` and `tseq2[0]` are different strings being compared, this measures how many leading chars of `tseq1+tseq2` match leading chars of `tseq2`, which is effectively how far the variant can be shifted to the right.
- Empty `tseq1` → loop at step 2 doesn't execute → `msicnt` = 0, `shift3` based only on `tseq2` comparison
- Empty `tseq2` → while loop at step 5 doesn't execute → `shift3` = 0

---

### `collectReferenceVariants(int position, int totalPosCoverage, Vars variationsAtPos, List<String> debugLines)` — L629–L969

**Purpose**: The core method that fills in all remaining Variant fields: genotype, refallele, varallele, MSI, CRISPR adjustment, flanking sequences, and reference variant data.

**Algorithm**:
1. Set `referenceForwardCoverage` = 0
2. Set `referenceReverseCoverage` = 0
3. **Determine genotype1**:
   - 3a. If `variationsAtPos.referenceVariant != null` AND `variationsAtPos.referenceVariant.frequency >= instance().conf.freq`:
     - `genotype1` = `variationsAtPos.referenceVariant.descriptionString`
   - 3b. Else if `variationsAtPos.variants.size() > 0`:
     - `genotype1` = `variationsAtPos.variants.get(0).descriptionString` (first/best variant)
   - 3c. Else:
     - `genotype1` = `variationsAtPos.referenceVariant.descriptionString` — **NPE risk if referenceVariant is null AND variants is empty**
4. If `variationsAtPos.referenceVariant != null`:
   - `referenceForwardCoverage` = `variationsAtPos.referenceVariant.varsCountOnForward`
   - `referenceReverseCoverage` = `variationsAtPos.referenceVariant.varsCountOnReverse`
5. **Adjust genotype1 for insertions/duplications**:
   - 5a. If `genotype1.startsWith("+")`:
     - Match `DUP_NUM` (`<dup(\d+)`) against genotype1
     - If match: `genotype1` = `"+" + (SVFLANK + toInt(group(1)))` — adds 50 to dup count
     - Else: `genotype1` = `"+" + (genotype1.length() - 1)` — length of insertion
6. **Determine genotype2** (declared, set per variant below)
7. **Adjust reference coverage from next position**:
   - 7a. If `totalPosCoverage > refCoverage.get(position)` AND `getNonInsertionVariants().containsKey(position + 1)` AND `ref.containsKey(position + 1)` AND `getNonInsertionVariants().get(position + 1).containsKey(ref.get(position + 1).toString())`:
     - Get `tpref` = `getVariationMaybe(getNonInsertionVariants(), position + 1, ref.get(position + 1))`
     - `referenceForwardCoverage` = `tpref.varsCountOnForward`
     - `referenceReverseCoverage` = `tpref.varsCountOnReverse`
8. Create `positionsForChangedRefVariant` = new ArrayList<Integer>()
9. **If non-reference variants exist** (`variationsAtPos.variants.size() > 0`):
   - **Loop** over `variationsAtPos.variants` as `vref`:
     - 9a. Set `genotype1current` = `genotype1` (per-variant copy)
     - 9b. Set `genotype2` = `vref.descriptionString`
     - 9c. If `genotype2.startsWith("+")`:
       - `genotype2` = `"+" + (genotype2.length() - 1)` — length of insertion
     - 9d. Set `descriptionString` = `vref.descriptionString`
     - 9e. Set `deletionLength` = 0
     - 9f. Match `BEGIN_MINUS_NUMBER` (`^-(\d+)`) against descriptionString
     - 9g. If match: `deletionLength` = `toInt(group(1))`
     - 9h. Set `endPosition` = `position`
     - 9i. If `descriptionString.startsWith("-")`:
       - `endPosition` = `position + deletionLength - 1`
     - 9j. Set `refallele` = `""`, declare `varallele`
     - 9k. Set `shift3` = 0, `msi` = 0, `msint` = `""`
     - 9l. Set `startPosition` = `position`
     - **Branch: Insertion** (descriptionString starts with `"+"`):
       - 9m. If `!descriptionString.contains("&")` AND `!descriptionString.contains("#")` AND `!descriptionString.contains("<dup")`:
         - Call `proceedVrefIsInsertion(position, descriptionString)` → get msi, shift3, msint
       - 9n. If `conf.moveIndelsTo3`:
         - `startPosition += shift3`, `endPosition += shift3`
       - 9o. `refallele` = `ref.containsKey(position) ? ref.get(position).toString() : ""`
       - 9p. `varallele` = `refallele + descriptionString.substring(1)` — ref base + insertion sequence
       - 9q. If `varallele.length() > conf.SVMINLEN`:
         - `endPosition += varallele.length()`, `varallele` = `"<DUP>"`
       - 9r. Match `DUP_NUM` against varallele:
         - If match: `dupCount` = toInt(group(1))
         - `endPosition` = `startPosition + (2 * SVFLANK + dupCount) - 1`
         - `genotype2` = `"+" + (2 * SVFLANK + dupCount)`
         - `varallele` = `"<DUP>"`
     - **Branch: Deletion** (descriptionString starts with `"-"`):
       - 9s. Match `INV_NUM` (`<inv(\d+)`) and `BEGIN_MINUS_NUMBER_CARET` (`^-\d+\^`)
       - 9t. If `deletionLength < conf.SVMINLEN`:
         - `varallele` = `descriptionString.replaceFirst("^-\\d+", "")` — strip leading -N
         - Call `proceedVrefIsDeletion(position, deletionLength)` → get msi, shift3, msint
         - If `INV_NUM` matched: `varallele` = `"<INV>"`, `genotype2` = `"<INV" + deletionLength + ">"`
       - 9u. Else if `BEGIN_MINUS_NUMBER_CARET` matched:
         - `varallele` = `"<INV>"`, `genotype2` = `"<INV" + deletionLength + ">"`
       - 9v. Else:
         - `varallele` = `"<DEL>"`
       - 9w. If `!descriptionString.contains("&")` AND `!descriptionString.contains("#")` AND `!descriptionString.contains("^")`:
         - If `conf.moveIndelsTo3`: `startPosition += shift3`
         - If `varallele != "<DEL>"`: `varallele` = `ref.containsKey(position - 1) ? ref.get(position - 1).toString() : ""`
         - `refallele` = `ref.containsKey(position - 1) ? ref.get(position - 1).toString() : ""`
         - `startPosition--`
       - 9x. Match `SOME_SV_NUMBERS` (`<(...)\d+>`) against descriptionString:
         - If match: `refallele` = `ref.containsKey(position) ? ref.get(position).toString() : ""`
         - Else if `deletionLength < conf.SVMINLEN`:
           - `refallele += joinRef(ref, position, position + deletionLength - 1)`
     - **Branch: SNP/MNP** (neither insertion nor deletion):
       - 9y. Build MSI for SNP:
         - `tseq1` = `joinRef(ref, max(position - 30, 1), position + 1)` — 30 left + 2 at position
         - `chr0` = `getOrElse(instance().chrLengths, region.chr, 0)`
         - `tseq2` = `joinRef(ref, position + 2, min(position + 70, chr0))`
         - Call `findMSI(tseq1, tseq2, null)` → msi, shift3, msint
       - 9z. `refallele` = `ref.containsKey(position) ? ref.get(position).toString() : ""`
       - 9aa. `varallele` = `descriptionString`
     - **Handle `&` (followed by matched sequence)**:
       - 9bb. Match `AMP_ATGC` (`&([ATGC]+)`) against descriptionString
       - 9cc. If match:
         - `extra` = group(1) — the matched following sequence
         - Remove `&` from varallele
         - Append `extra.length()` bases from ref to refallele and genotype1current
         - `endPosition += extra.length()`
         - Match `AMP_ATGC` again against varallele for nested `&`:
           - If match: same treatment — remove `&`, extend refallele/genotype1current, extend endPosition
         - If descriptionString starts with `"+"`:
           - `refallele` = `refallele.substring(1)`, `varallele` = `varallele.substring(1)`, `startPosition++`
         - If `varallele.equals("<DEL>")` AND `refallele.length() >= 1`:
           - `refallele` = ref at startPosition
           - If refCoverage has `startPosition - 1`: totalPosCoverage = that coverage
           - If `vref.positionCoverage > totalPosCoverage`: totalPosCoverage = positionCoverage
           - Recalculate `vref.frequency` = `positionCoverage / (double)totalPosCoverage`
     - **Handle `#...^...` (matched sequence + indel tail)**:
       - 9dd. Match `HASH_GROUP_CARET_GROUP` (`#(.+)\^(.+)`)
       - 9ee. If match:
         - `matchedSequence` = group(1), `tail` = group(2)
         - `endPosition += matchedSequence.length()`
         - Append ref bases to refallele for the matched positions
         - Match `BEGIN_DIGITS` against tail:
           - If match: `deletion` = int value → append deletion bases from ref to refallele; `endPosition += deletion`
         - Clean `#` and `^(\d+)?` from varallele
         - Replace `#` → `m` and `^` → `i` in genotype1current and genotype2
     - **Handle `^` (insertion after deletion in novolign)**:
       - 9ff. Match `CARET_ATGNC` (`\^([ATGNC]+)`)
       - 9gg. If match: remove `^` from varallele; replace `^` → `i` in genotypes
     - **CRISPR adjustment** (conf.crisprCuttingSite != 0):
       - 9hh. Complex trimming for 5' (refallele/varallele share prefix → strip prefix)
       - 9ii. Complex shifting towards cut site (shift startPosition/endPosition, rebuild alleles)
       - 9jj. Set `vref.crispr` = n (amount shifted)
     - **Set flanking sequences**:
       - 9kk. `vref.leftseq` = `joinRef(ref, max(startPosition - 20, 1), startPosition - 1)`
       - 9ll. `vref.rightseq` = `joinRef(ref, endPosition + 1, min(endPosition + 20, chr0))`
     - **Build genotype string**:
       - 9mm. `genotype` = `genotype1current + "/" + genotype2`
       - 9nn. Remove `&`, `#` from genotype; replace `^` → `i`
     - **Round and set final fields**:
       - 9oo. `vref.extraFrequency` = `roundHalfEven("0.0000", vref.extraFrequency)`
       - 9pp. `vref.frequency` = `roundHalfEven("0.0000", vref.frequency)`
       - 9qq. `vref.highQualityReadsFrequency` = `roundHalfEven("0.0000", vref.highQualityReadsFrequency)`
       - 9rr. `vref.msi` = `roundHalfEven("0.000", msi)` (3 decimal places)
       - 9ss. `vref.msint` = `msint.length()` — **int**, length of MSI unit string
       - 9tt. `vref.shift3` = `shift3`
       - 9uu. `vref.startPosition` = `startPosition`
       - 9vv. `vref.endPosition` = `endPosition`
       - 9ww. `vref.refallele` = `validateRefallele(refallele)`
       - 9xx. `vref.varallele` = `varallele`
       - 9yy. `vref.genotype` = `genotype`
       - 9zz. `vref.totalPosCoverage` = `totalPosCoverage`
       - 9aaa. `vref.refForwardCoverage` = `referenceForwardCoverage`
       - 9bbb. `vref.refReverseCoverage` = `referenceReverseCoverage`
     - **Set strand bias flag**:
       - 9ccc. If `variationsAtPos.referenceVariant != null`:
         - `vref.strandBiasFlag` = `referenceVariant.strandBiasFlag + ";" + vref.strandBiasFlag`
       - Else:
         - `vref.strandBiasFlag` = `"0;" + vref.strandBiasFlag`
     - 9ddd. Call `adjustVariantCounts(position, vref)`
     - 9eee. If `startPosition != position` AND `conf.doPileup`:
       - Add `position` to `positionsForChangedRefVariant`
     - 9fff. Call `constructDebugLines(debugLines, vref)`
   - **After variant loop**: If `conf.disableSV`:
     - Remove variants where `ANY_SV.matcher(vref.varallele).find()` is true
10. **Else if** `variationsAtPos.referenceVariant != null` (no non-ref variants):
    - Call `updateRefVariant(position, totalPosCoverage, vref, debugLines, referenceForwardCoverage, referenceReverseCoverage)`
11. **Else** (no variants at all):
    - `variationsAtPos.referenceVariant` = new `Variant()`
12. **Pileup ref variant update**: If `referenceVariant != null` AND `conf.doPileup` AND (`positionsForChangedRefVariant.contains(position)` OR `ampliconBasedCalling != null`):
    - Call `updateRefVariant(position, totalPosCoverage, vref, debugLines, referenceForwardCoverage, referenceReverseCoverage)`

**Parity-Critical issues** (see Known Parity Traps for full list):
- `genotype1` determination at step 3 — if referenceVariant is null and variants list is empty, NPE
- The per-variant `genotype1current` copy is important — `genotype1` is reused across variants
- `totalPosCoverage` is a local parameter; it can be mutated in the `&`+`<DEL>` branch (step 9cc)
- The deletion branch's `startPosition--` (step 9w) affects all subsequent calculations
- `SOME_SV_NUMBERS` match at step 9x overrides the refallele built in step 9w
- Multiple regex matches on the same descriptionString — order matters
- `vref.strandBiasFlag` is set in createVariant/createInsertion as just bias value, then **prepended** with reference bias here

---

### `updateRefVariant(int position, int totalPosCoverage, Variant vref, List<String> debugLines, int referenceForwardCoverage, int referenceReverseCoverage)` — L971–L1013

**Purpose**: Fill reference variant fields when only reference reads are observed (or for pileup update).

**Algorithm**:
1. `vref.totalPosCoverage` = totalPosCoverage
2. `vref.positionCoverage` = 0
3. `vref.frequency` = 0
4. `vref.refForwardCoverage` = referenceForwardCoverage
5. `vref.refReverseCoverage` = referenceReverseCoverage
6. `vref.varsCountOnForward` = 0
7. `vref.varsCountOnReverse` = 0
8. `vref.msi` = 0
9. `vref.msint` = 0
10. If `vref.strandBiasFlag.indexOf(';')` == -1:
    - Append `";0"` to strandBiasFlag
11. `vref.shift3` = 0
12. `vref.startPosition` = position
13. `vref.endPosition` = position
14. `vref.highQualityReadsFrequency` = `roundHalfEven("0.0000", vref.highQualityReadsFrequency)`
15. Get `referenceBase` = ref at position or `""` if missing
16. `vref.refallele` = `validateRefallele(referenceBase)`
17. `vref.varallele` = `validateRefallele(referenceBase)`
18. `vref.genotype` = `referenceBase + "/" + referenceBase`
19. `vref.leftseq` = `""`
20. `vref.rightseq` = `""`
21. `vref.duprate` = duprate
22. Call `constructDebugLines(debugLines, vref)`

**Parity-Critical**:
- `strandBiasFlag.indexOf(';')` — checks if semicolon already present. If not (meaning it came from createVariant/createInsertion without being processed by collectReferenceVariants), appends `";0"`
- `highQualityReadsFrequency` is rounded but NOT set to 0 — it retains the value from createVariant/createInsertion, just rounded
- `leftseq` and `rightseq` set to empty strings — NOT computed from reference

---

### `constructDebugLines(List<String> debugLines, Variant vref)` — L1018–L1030

**Purpose**: Build DEBUG field from accumulated debug lines.

**Algorithm**:
1. If `instance().conf.debug`:
   - Create StringBuilder sb
   - For each `str` in debugLines:
     - If `sb.length() > 0`: append `" & "`
     - Append str
   - Set `vref.DEBUG` = `sb.toString()`

---

### `validateRefallele(String refallele)` — L1039–L1048

**Purpose**: Replace IUPAC ambiguity codes with first-alphabetical standard base per VCF 4.3 spec.

**Algorithm**:
1. **Loop** `i` from 0 to `refallele.length() - 1`:
   - 1a. Extract `refBase` = `substr(refallele, i, 1)` — single char at position i
   - 1b. If `IUPAC_AMBIGUITY_CODES.containsKey(refBase)`:
     - `refallele` = `refallele.replaceFirst(refBase, IUPAC_AMBIGUITY_CODES.get(refBase))`
2. Return `refallele`

**IUPAC_AMBIGUITY_CODES map** (field at L54–L65):

| Code | Replacement |
|------|-------------|
| M | A |
| R | A |
| W | A |
| S | C |
| Y | C |
| K | G |
| V | A |
| H | A |
| D | A |
| B | C |

**Parity-Critical**:
- `replaceFirst` replaces **the first occurrence** of the ambiguity code — if the same code appears multiple times, only the first is replaced per iteration. But the loop advances `i`, so subsequent occurrences in later positions would be caught in their own iteration.
- However: `replaceFirst` uses **regex**, and the IUPAC codes (single uppercase letters) are valid regex patterns. No special characters to escape.
- Edge case: if `refallele` is empty, loop doesn't execute.

---

### Inner class `MSI` — L1053–L1073

Simple data holder:
- `double msi` — MSI repeat count
- `int shift3` — bases shiftable to 3'
- `String msint` — MSI unit sequence

Constructor takes all three fields.

---

### Inner class `ToVarsJsonlWriter` — L1075–L1194

**Purpose**: Diagnostic JSONL writer triggered by `VARDICT_TO_VARS_JSONL` env var. NOT parity-critical for output — used only for debugging.

Key methods:
- `write(path, region, data, duprate)` — writes META line, then all positions sorted by TreeMap in ascending order
- `variantJson(v)` — serializes Variant to JSON with specific float formatting
- `fmtDouble(pattern, value)` — uses `DecimalFormat` with `RoundingMode.HALF_EVEN` and `Locale.US`
- `escapeJson(value)` — escapes `"`, `\`, `\n`, `\r`, `\t`

Not analyzed in detail — not parity-critical.

## Cross-Module Dependencies

### Outbound Calls

| Target | Method | Where Used | Purpose |
|--------|--------|------------|---------|
| `VariationUtils` | `strandBias(fwd, rev)` | createVariant, createInsertion | Compute strand bias flag (0/1/2) |
| `VariationUtils` | `getOrPutVars(map, position)` | process, collectVarsAtPosition, createVariant, collectReferenceVariants | Get or create Vars at position |
| `VariationUtils` | `getVariationMaybe(hash, start, refBase)` | createInsertion, collectReferenceVariants | Safe lookup of Variation by position+refBase |
| `VariationUtils` | `joinRef(ref, from, to)` | proceedVrefIsDeletion, proceedVrefIsInsertion, collectReferenceVariants, findMSI | Build reference sequence string |
| `Utils` | `roundHalfEven(pattern, value)` | createVariant, createInsertion, collectReferenceVariants, updateRefVariant | Format double with HALF_EVEN rounding |
| `Utils` | `substr(string, idx)` / `substr(string, idx, len)` | proceedVrefIsDeletion, proceedVrefIsInsertion, findMSI, collectReferenceVariants, validateRefallele | Perl-style substring |
| `Utils` | `getOrElse(map, key, default)` | proceedVrefIsDeletion, proceedVrefIsInsertion, collectReferenceVariants | Map lookup with default (mutates map) |
| `Utils` | `toInt(string)` | collectReferenceVariants | String → int |
| `Utils` | `printExceptionAndContinue(...)` | process | Error handling |
| `GlobalReadOnlyScope` | `instance()` | throughout | Access configuration |
| `Patterns` | Various static Pattern fields | collectReferenceVariants | Regex patterns |
| `Configuration` | `SVFLANK` (=50) | collectReferenceVariants | SV flank size constant |

### Inbound Callers

| Caller | Method | Notes |
|--------|--------|-------|
| `AbstractMode.process()` | `.thenApply(new ToVarsBuilder()::process)` | Pipeline stage in CompletableFuture chain |
| `ToVarsBuilderTest` | Direct calls to `createInsertion`, `createVariant`, `validateRefallele` | Unit tests |

## Data Structures Read/Written

| Structure | Java Type | R/W | Notes |
|-----------|-----------|-----|-------|
| `nonInsertionVariants` | `Map<Integer, VariationMap<String, Variation>>` | R (W on next-pos ref) | Outer map is HashMap; inner is LinkedHashMap (VariationMap). Written: fwd/rev adjusted in createInsertion |
| `insertionVariants` | `Map<Integer, VariationMap<String, Variation>>` | R | Same structure as above |
| `refCoverage` | `Map<Integer, Integer>` | R | Position → total coverage |
| `ref` | `Map<Integer, Character>` | R | Position → reference base |
| `region` | `Region` | R | .chr, .start, .end |
| `duprate` | `double` | R | Duplication rate |
| `alignedVariants` | `Map<Integer, Vars>` | W | Output: position → Vars |
| `Vars.referenceVariant` | `Variant` | W | Reference variant at position |
| `Vars.variants` | `List<Variant>` | W | Non-reference variants at position |
| `Vars.varDescriptionStringToVariants` | `Map<String, Variant>` | W | Key: description string |
| `Vars.sv` | `String` | W | SV info: "splits-pairs-clusters" |
| `instance().chrLengths` | `Map<String, Integer>` | R/W | getOrElse mutates by inserting default |
| `instance().conf` | `Configuration` | R | Various config fields |

## Known Parity Traps

### Float Formatting & Rounding

1. **`roundHalfEven` pattern differences**: `"0.0000"` for frequency/extraFreq/hiFreq, `"0.0"` for quality/position/mismatches, `"0.000"` for MSI. Each pattern has different decimal places. Rust must use the exact same `DecimalFormat` HALF_EVEN rounding as Java. The Java `roundHalfEven` parses the formatted string back to double: `Double.parseDouble(new DecimalFormat(pattern).format(value))`. This double-conversion (double→string→double) introduces specific rounding behavior that must be replicated exactly.

2. **Double division for frequency**: `cnt.varsCount / (double)ttcov` — when `ttcov` is 0, this produces `Infinity` or `NaN`. Java's `DecimalFormat.format(Infinity)` returns `"∞"`, which `Double.parseDouble` would reject. Need to verify this case can't happen (it shouldn't because of coverage checks).

3. **`highQualityToLowQualityRatio` divisor**: `hicnt / (locnt != 0 ? locnt : 0.5d)`. When `locnt == 0`, Java divides by `0.5` (integer hicnt promoted to double). Rust must not accidentally use `0.5f32` or integer division.

### Collection Ordering

4. **`VariationMap` extends `LinkedHashMap`**: Both `nonInsertionVariants` values and `insertionVariants` values use LinkedHashMap iteration order. In `process()`, the keys are explicitly sorted before passing to `createVariant` and `createInsertion`, so **within a single position** the iteration order of keys doesn't matter. However, the **outer loop** over `getNonInsertionVariants().entrySet()` IS in insertion order. Rust must use `IndexMap` for the outer map or process positions in the same order.

5. **`varsAtCurPosition.keySet()` sorted**: Both `createVariant` and `createInsertion` sort keys lexicographically with `Collections.sort`. Rust must sort `String` keys the same way (byte-wise ASCII ordering matches Java's natural String ordering for ASCII data).

6. **`alignedVariants` is `HashMap`**: The output map is a plain HashMap, so its iteration order is non-deterministic. Downstream consumers must not depend on iteration order.

7. **`Vars.variants` list order**: After `sortVariants()`, variants are ordered by `meanQuality * positionCoverage` descending, then `descriptionString` ascending. This order affects `genotype1` selection in `collectReferenceVariants` (index 0 = best variant).

### Null Handling → Option

8. **`ref.get(position)` → `null`**: Multiple places check `ref.containsKey(position)` before access, but some use `ref.get(position)` directly (e.g., `collectVarsAtPosition` does `String.valueOf(ref.get(position))` — Java's `String.valueOf(null)` returns `"null"` string). In Rust, this would be `None` and must be handled differently.

9. **`variationsAtPos.referenceVariant` can be null**: Checked at step 3 of `collectReferenceVariants`. If null and variants list is empty, step 3c accesses `referenceVariant.descriptionString` which would NPE. In practice this shouldn't happen because the variant was already created, but Rust must handle the Option.

10. **`getInsertionVariants().get(position)` returns null**: Checked in `createInsertion` (the whole block is skipped if null). Also checked in `isTheSameVariationOnRef` with `containsKey`.

11. **`refCoverage.get(position)` auto-unboxing**: Java `Map<Integer, Integer>.get()` returns `Integer` (boxed). If null, auto-unboxing to `int` throws NPE. The code checks `containsKey` first in most places but not all (step 7a of collectReferenceVariants does `refCoverage.get(position)` — but `refCoverage` was already verified to contain position at step 6h of process()).

### Integer Arithmetic

12. **`extracnt != 0` vs `extracnt > 0`**: `createInsertion` uses `cnt.extracnt != 0` (catches negative), while `createVariant` uses `cnt.extracnt > 0` (only positive). This is a **real behavioral difference** — if extracnt were negative, createInsertion would trigger the adjustment but createVariant would not.

13. **SV string format**: `sv.splits + "-" + sv.pairs + "-" + sv.clusters` — all ints, concatenated with dashes. Java int-to-String conversion is decimal. Order is splits-pairs-clusters.

14. **`msint.length()`**: The `vref.msint` field is set to `msint.length()` (int), not the string itself. The MSI unit string length is stored, not the sequence.

### String/Regex Operations

15. **`descriptionString.replaceFirst("^-\\d+", "")`**: Uses regex `replaceFirst`. In Rust, this needs equivalent regex behavior. The `^` anchor means it only matches at the start.

16. **`genotype.replace("&", "").replace("#", "").replace("^", "i")`**: `replace` (not `replaceFirst`) — replaces ALL occurrences. After this, the `^` symbols become `i` everywhere.

17. **`varallele.replaceFirst("&", "")`**: Only the FIRST `&` is removed from varallele.

18. **`varallele.replaceFirst("#", "").replaceFirst("\\^(\\d+)?", "")`**: Two sequential replaceFirst calls. The `^` regex is `\\^(\\d+)?` which matches caret optionally followed by digits.

19. **`genotype1current.replaceFirst("#", "m").replaceFirst("\\^", "i")`**: Replace first `#` with `m`, first `^` with `i`.

### MSI/shift3 Computation

20. **Regex compilation in loop**: `findMSI` compiles `Pattern.compile("((" + msint + ")+)$")` fresh each iteration. The `msint` string comes from `substr(tseq1, -nmsi, nmsi)` and contains raw nucleotide characters (ATGC), which are safe in regex. But if the variant contains special regex characters, this could be a problem.

21. **`shift3` not updated by second findMSI call**: In both `proceedVrefIsDeletion` and `proceedVrefIsInsertion`, the second `findMSI` call (with `left=null`) may update `msi` and `msint` but **never** updates `shift3`. The shift3 from the first call is preserved.

22. **`msi <= shift3 / (double)dellen`**: Uses `<=` not `<`. If current msi equals the ratio, it gets replaced.

### Mutation Side Effects

23. **`getOrElse` mutates `chrLengths`**: `Utils.getOrElse(instance().chrLengths, region.chr, 0)` inserts 0 into the shared chrLengths map if the chromosome isn't found. This is a side effect that persists across calls.

24. **`createInsertion` mutates next-position Variation**: In step 2c.xii, if coverage adjustment triggers, `variantNextPosition.varsCountOnForward -= fwd` and `varsCountOnReverse -= rev` are modified on the Variation object in `nonInsertionVariants`. This affects any subsequent processing of position+1.

25. **`totalPosCoverage` can be modified in `collectReferenceVariants`**: In the `&`+`<DEL>` branch (step 9cc), `totalPosCoverage` (a local variable/parameter) is overwritten. This affects `vref.frequency` recalculation.

### Control Flow / Branching

26. **SV position bypass**: At step 6e of `process()`, positions with `varsAtCurPosition.sv != null` AND `!conf.deleteDuplicateVariants` bypass the region boundary check. This means SV positions outside the region can be processed.

27. **`isTheSameVariationOnRef` short-circuit**: Only triggers when exactly one key exists and it matches the reference base, AND no pileup/amplicon/somatic mode. The insertion check adds `"I"` as a synthetic key to prevent skipping positions that have insertions.

28. **Error-and-continue in process loop**: The `try/catch` wraps the entire per-position body. If any exception occurs (e.g., NPE from null ref), the position is skipped and processing continues. Rust doesn't have this implicit error recovery — must be implemented explicitly.

29. **`disableSV` removal**: After the variant loop in `collectReferenceVariants`, variants with SV-type varallele are removed using `removeIf`. This happens after all variant processing, so the variants were fully constructed before removal.

30. **Pileup ref variant double-update**: Step 12 in `collectReferenceVariants` can cause `updateRefVariant` to be called AFTER the variants loop already processed the reference variant. This overwrites fields like `positionCoverage` = 0, `frequency` = 0, etc. This only triggers in pileup mode with either changed start position or amplicon calling.

### CRISPR-Specific

31. **CRISPR adjustment modifies startPosition, endPosition, refallele, varallele**: Only when `conf.crisprCuttingSite != 0`. The adjustment attempts to shift the variant towards the cut site by up to `shift3` positions. The allele sequences are rebuilt from reference after shifting.

32. **CRISPR 5' fix for complex variants**: When both refallele and varallele share a common prefix, the prefix is stripped and startPosition is advanced. This happens before the 3' CRISPR adjustment.

### Output Field Semantics

33. **`strandBiasFlag` format**: For non-reference variants, the final format is `"refBias;varBias"` where each is 0, 1, or 2. If referenceVariant is null, refBias is `"0"`. The varBias was set in createVariant/createInsertion.

34. **`genotype` format**: `"genotype1/genotype2"` where genotype1 is either reference descriptionString, or best variant descriptionString, or dup-adjusted string. genotype2 is the current variant's descriptionString (or adjusted for insertions/dups/INV). Special chars are cleaned: `&` removed, `#` removed, `^` → `i`.

35. **Deletion `startPosition--`**: For simple deletions (no `&`, `#`, `^`), startPosition is decremented to include the preceding base. This means the refallele starts one base before the deletion. The varallele is also set to that preceding base (or `<DEL>` for large deletions).

36. **`updateRefVariant` sets `leftseq`/`rightseq` to empty**: Reference-only positions get empty flanking sequences, NOT the actual reference context. This differs from non-reference variants which get 20bp flanking.
