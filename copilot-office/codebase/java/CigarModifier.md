# CigarModifier

**Source**: `modules/CigarModifier.java` (L1–L787)  
**LOC**: 787  
**Rust counterpart**: `src/mods/cigar_modifier.rs`  
**Status**: complete

## Overview

CigarModifier sits between RecordPreprocessor and CigarParser in the VarDict pipeline. It transforms raw CIGAR strings from BAM records to normalize edge-case alignments before CigarParser processes them for variant calling. Its responsibilities:

1. **Strip leading/trailing deletions** — deletions at read boundaries are artifacts; position is adjusted for leading D  
2. **Convert boundary insertions to soft clips** — insertions at start/end of reads become S operations  
3. **Chimeric soft-clip removal** — large soft clips whose reverse-complement seeds map uniquely and nearby in the reference are removed (chimeric artifact suppression)  
4. **Collapse indel+soft-clip edges** — indels adjacent to soft clips are absorbed into the soft clip  
5. **Collapse short-match+indel edges** — short matched bases (≤9 or ≤10) sandwiched between soft clips and indels get absorbed  
6. **Merge adjacent/nearby indels** — two deletions, two insertions, or D-M-I-M-D complexes within 15bp are merged into a single complex indel  
7. **Soft-clip boundary realignment** — mismatches at read boundary M↔S transitions are reclassified to extend or shrink the soft-clip    
8. **Mismatch-to-softclip conversion** — ≥3 consecutive mismatches at either read end get converted to soft-clips

Called **only** when `performLocalRealignment` is true (the default). Returns a `ModifiedCigar(position, cigarStr, querySequence, queryQuality)`.

## Method Inventory

| Method | Lines | Visibility | Analyzed | Summary |
|--------|-------|------------|----------|---------|
| `CigarModifier()` (constructor) | L39–L51 | package | yes | Initialize all fields from SAM record data |
| `modifyCigar()` | L57–L248 | public | yes | Master method: applies all CIGAR transformations in sequence |
| `combineBeginDigM()` | L255–L278 | private | yes | Convert leading mismatches to soft-clip |
| `combineDigSDigM()` | L284–L358 | private | yes | Extend/shrink 5' soft-clip boundary by matching reference |
| `captureMisSoftly3Mismatches()` | L364–L399 | private | yes | Convert ≥3 trailing mismatches to soft-clip (M-only CIGAR) |
| `captureMisSoftlyMS()` | L405–L486 | private | yes | Extend 3' M into S, or retract M→S boundary for mismatches |
| `combineToCloseToOne()` | L494–L527 | private | yes | Merge I-M-D/I complex (≤15bp match) into D+I |
| `combineToCloseToCorrect()` | L535–L569 | private | yes | Merge D-M-D/I complex (≤15bp match) into D+I |
| `threeIndels()` | L578–L642 | private | yes | Merge 3-indel complex M-D/I-M-D/I-M-D/I-M into simplified form |
| `threeDeletions()` | L651–L701 | private | yes | Merge M-D-M-D-M-D-M into D+optional-I+M |
| `twoDeletionsInsertionToComplex()` | L710–L762 | private | yes | Merge M-D-M-I-M-D-M into D+I+M |
| `beginDigitMNumberIorDNumberM()` | L771–L787 | private | yes | Convert leading short-M+indel into soft-clip |

## Method Analyses

### Constructor: `CigarModifier()`

**Source**: [CigarModifier.java](VarDictJava/src/main/java/com/astrazeneca/vardict/modules/CigarModifier.java#L39-L51)  
**Purpose**: Store all per-read inputs for CIGAR modification.

#### Parameters
- `position`: `int` — 1-based alignment start from `record.getAlignmentStart()`
- `cigarStr`: `String` — CIGAR string from `record.getCigarString()`
- `querySequence`: `String` — read bases
- `queryQuality`: `String` — quality string (ASCII-encoded phred+33)
- `ref`: `Reference` — reference data (contains `referenceSequences` map and `seed` map)
- `indel`: `int` — `insertionDeletionLength` from CigarParser (total indel bases in raw CIGAR)
- `region`: `Region` — BED region for error reporting
- `maxReadLength`: `int` — maximum read length in current batch

#### Algorithm
1. Set `this.position = position`
2. Set `this.cigarStr = cigarStr`
3. Set `this.originalCigar = cigarStr` (saved for error reporting — never mutated)
4. Set `this.querySequence = querySequence`
5. Set `this.queryQuality = queryQuality`
6. Set `this.ref = ref`
7. Set `this.reference = ref.referenceSequences` (shorthand `Map<Integer, Character>`)
8. Set `this.indel = indel`
9. Set `this.region = region`
10. Set `this.maxReadLength = maxReadLength`

---

### `modifyCigar()`

**Source**: [CigarModifier.java](VarDictJava/src/main/java/com/astrazeneca/vardict/modules/CigarModifier.java#L57-L248)  
**Purpose**: Master transformation method. Applies all CIGAR normalizations and returns a `ModifiedCigar`.  
**Return**: `ModifiedCigar(position, cigarStr, querySequence, queryQuality)`

#### Algorithm

1. Set `flag = true` (controls the while-loop for iterative normalization).
2. **Enter try block** — all transformations are wrapped in a try/catch.

**Phase 1: Strip boundary deletions (L63–L72)**

3. Match `cigarStr` against `BEGIN_NUMBER_D` (`^\(\d+\)D`).
4. If match: add the matched number to `position` (skip over the deletion), then remove the matched portion from `cigarStr` using `Replacer("")`.
5. Match `cigarStr` against `END_NUMBER_D` (`\(\d+\)D$`).
6. If match: remove the trailing deletion from `cigarStr` (no position adjustment — trailing D doesn't affect alignment start).

**Phase 2: Convert boundary insertions to soft clips (L73–L82)**

7. Match `cigarStr` against `BEGIN_NUMBER_I` (`^\(\d+\)I`).
8. If match: replace leading `NI` with `NS` where N = same number (insertions at start become soft clips).
9. Match `cigarStr` against `END_NUMBER_I` (`\(\d+\)I$`).
10. If match: replace trailing `NI` with `NS`.

**Phase 3: Chimeric soft-clip detection (L83–L131)**

11. Match `cigarStr` against `SA_CIGAR_D_S_5clip_GROUP` (`^(\d\d+)S`) — 5' soft clip with ≥2 digits (≥10 bases).
12. Match `cigarStr` against `SA_CIGAR_D_S_3clip_GROUP` (`(\d\d+)S$`) — 3' soft clip with ≥2 digits.
13. Get `referenceSeedMap = ref.seed`.
14. **If 5' soft clip matches** (and 3' didn't, due to `else if`):
    - a. Get `cigarElement = toInt(group(1))` — length of the soft clip.
    - b. Check `!instance().conf.chimeric && cigarElement >= Configuration.SEED_2` (SEED_2 = 12).
    - c. Extract `sseq = substr(querySequence, 0, cigarElement)` — soft-clipped bases.
    - d. Compute `sequence = complement(reverse(sseq))` — reverse complement.
    - e. Take `reverseComplementedSeed = sequence.substring(0, Configuration.SEED_2)` — first 12 chars of revcomp.
    - f. Check if `referenceSeedMap.containsKey(reverseComplementedSeed)`.
    - g. If yes, get `positions = referenceSeedMap.get(reverseComplementedSeed)`.
    - h. If `positions.size() == 1` (unique hit) AND `Math.abs(position - positions.get(0)) < 2 * maxReadLength`:
        - Remove leading `\d+S` from `cigarStr` (via `SA_CIGAR_D_S_5clip_GROUP_Repl`).
        - Trim `querySequence = substr(querySequence, cigarElement)` (remove first cigarElement chars).
        - Trim `queryQuality = substr(queryQuality, cigarElement)`.
        - Debug log if `instance().conf.y`.
15. **Else if 3' soft clip matches**:
    - a. Get `cigarElement = toInt(group(1))`.
    - b. Same chimeric check.
    - c. Extract `sseq = substr(querySequence, -cigarElement, cigarElement)` — last cigarElement bases.
    - d. Compute `sequence = complement(reverse(sseq))`.
    - e. Take `reverseComplementedSeed = substr(sequence, -Configuration.SEED_2, Configuration.SEED_2)` — last 12 chars of revcomp.
    - f–h. Same seed lookup logic as 5' case.
        - Remove trailing `\d\d+S$` from `cigarStr`.
        - Trim query strings from the end: `substr(querySequence, 0, querySequence.length() - cigarElement)`.

**Phase 4: Iterative indel normalization loop (L133–L229)**

16. Enter `while (flag && indel > 0)` loop. Set `flag = false` at start of each iteration.

**4a. Leading S+I/D collapse (L136–L145)**
17. Match `BEGIN_NUMBER_S_NUMBER_IorD` (`^(\d+)S(\d+)([ID])`).
18. If match:
    - Compute `tslen = group(1) + (group(3)=="I" ? group(2) : 0)` + "S" — absorb insertion into soft clip, or keep just old soft-clip length for deletion.
    - Adjust `position += (group(3)=="D" ? group(2) : 0)` — skip over deletion.
    - Replace matched portion with `tslen`.
    - Set `flag = true`.

**4b. Trailing I/D+S collapse (L146–L153)**
19. Match `NUMBER_IorD_NUMBER_S_END` (`(\d+)([ID])(\d+)S$`).
20. If match:
    - Compute `tslen = group(3) + (group(2)=="I" ? group(1) : 0)` + "S".
    - Replace matched portion.
    - Set `flag = true`.

**4c. Leading S+short-M+I/D collapse (L154–L165)**
21. Match `BEGIN_NUMBER_S_NUMBER_M_NUMBER_IorD` (`^(\d+)S(\d+)M(\d+)([ID])`).
22. If match and `tmid = group(2) <= 10`:
    - `tslen = group(1) + tmid + (group(4)=="I" ? group(3) : 0)` + "S".
    - `position += tmid + (group(4)=="D" ? group(3) : 0)`.
    - Replace, set `flag = true`.

**4d. Trailing I/D+short-M+S collapse (L166–L176)**
23. Match `NUMBER_IorD_NUMBER_M_NUMBER_S_END` (`(\d+)([ID])(\d+)M(\d+)S$`).
24. If match and `tmid = group(3) <= 10`:
    - `tslen = group(4) + tmid + (group(2)=="I" ? group(1) : 0)` + "S".
    - Replace, set `flag = true`.

**4e. Leading short-M+I/D+M → soft-clip (L180–L184)**
25. Match `BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M` (`^(\d)M(\d+)([ID])(\d+)M`).
26. If match: call `beginDigitMNumberIorDNumberM(mm)`, set `flag` to return value.

**4f. Trailing I/D+short-M → soft-clip (L185–L191)**
27. Match `NUMBER_IorD_DIGIT_M_END` (`(\d+)([ID])(\d)M$`).
28. If match:
    - `tmid = group(3)`, `tslen = tmid + (group(2)=="I" ? group(1) : 0)` + "S".
    - Replace using `NUMBER_IorD_NUMBER_M_END` (note: different pattern from the matcher — uses `\d+` not `\d`).
    - Set `flag = true`.

**4g. Complex indel merging (L194–L210)**
29. Match `D_M_D_DD_M_D_I_D_M_D_DD` (M-D-M-I-M-D pattern).
30. Also match `threeIndelsPattern` (M-D/I-M-D/I-M-D/I-M).
31. Also match `threeDeletionsPattern` (M-D-M-D-M-D-M).
32. If D-M-I-M-D matches: call `twoDeletionsInsertionToComplex(mm, flag)`.
33. Else if threeDeletions matches: call `threeDeletions(matcher, flag)`.
34. Else if threeIndels matches: call `threeIndels(matcher, flag)`.

**4h. D-M-D/I merge (L208–L212)**
35. Match `DIG_D_DIG_M_DIG_DI_DIGI` (`(\d+)D(\d+)M(\d+)([DI])(\d+I)?`).
36. If match: call `combineToCloseToCorrect(cm, flag)`.

**4i. Non-D-prefix I-M-D/I merge (L214–L217)**
37. Match `NOTDIG_DIG_I_DIG_M_DIG_DI_DIGI` (`(\D)(\d+)I(\d+)M(\d+)([DI])(\d+I)?`).
38. If match AND `group(1)` is NOT "D" and NOT "H": call `combineToCloseToOne(cm, flag)`.

**4j. Adjacent D+D or I+I merge (L218–L227)**
39. Match `DIG_D_DIG_D` (`(\d+)D(\d+)D`).
40. If match: sum both numbers, replace with `sum + "D"`, set `flag = true`.
41. Match `DIG_I_DIG_I` (`(\d+)I(\d+)I`).
42. If match: sum both numbers, replace with `sum + "I"`, set `flag = true`.

**End of while loop**.

**Phase 5: Post-loop soft-clip realignment (L231–L247)**

43. Match `ANY_NUMBER_M_NUMBER_S_END` (`^(.*?)(\d+)M(\d+)S$`).
44. If match: call `captureMisSoftlyMS(mtch)`.
45. Else match `BEGIN_ANY_DIG_M_END` (`^(.*?)(\d+)M$`).
46. If match: call `captureMisSoftly3Mismatches(mtch)`.

**Phase 6: 5' soft-clip realignment (L248–L253)**

47. Match `DIG_S_DIG_M` (`^(\d+)S(\d+)M`).
48. If match: call `combineDigSDigM(mtch)`.
49. Else match `BEGIN_DIG_M` (`^(\d+)M`).
50. If match: call `combineBeginDigM(mtch)`.

51. Return `new ModifiedCigar(position, cigarStr, querySequence, queryQuality)`.

**Exception handling (L245–L247)**:
52. On any `Exception`: call `printExceptionAndContinue(exception, "cigar", position + " " + originalCigar, region)`.
53. Still return `new ModifiedCigar(position, cigarStr, querySequence, queryQuality)` (partially-modified state).

#### Mutable State Changes
- `this.position` — adjusted for stripped deletions, absorbed short matches
- `this.cigarStr` — repeatedly rewritten by regex replacements
- `this.querySequence` — trimmed only by chimeric detection (phase 3)
- `this.queryQuality` — trimmed only by chimeric detection (phase 3)

---

### `combineBeginDigM()`

**Source**: [CigarModifier.java](VarDictJava/src/main/java/com/astrazeneca/vardict/modules/CigarModifier.java#L255-L278)  
**Purpose**: If the CIGAR starts with `M` and the first few bases are mismatches, convert them to soft-clip.

#### Parameters
- `matcher`: `Matcher` — matched against `BEGIN_DIG_M` (`^(\d+)M`)

#### Algorithm
1. Extract `mch = toInt(matcher.group(1))` — length of the leading M block.
2. Initialize `rn = 0` (mismatch-tracking counter), `rrn = 0` (position walker), `rmch = 0` (consecutive-match counter).
3. **While** `rrn < mch` AND `rn < mch`:
    - a. If `reference` does not contain key `position + rrn`: **break**.
    - b. If `isHasAndNotEquals(reference, position + rrn, querySequence, rrn)` (mismatch):
        - Set `rn = rrn + 1` (record how many mismatching bases up to and including this one).
        - Reset `rmch = 0`.
    - c. Else if `isHasAndEquals(reference, position + rrn, querySequence, rrn)` (match):
        - Increment `rmch`.
    - d. Increment `rrn`.
    - e. If `rmch >= 3`: **break** (stop scanning once 3 consecutive matches found).
4. If `rn > 0` AND `rn <= 3`:
    - `mch -= rn`.
    - Replace leading M with `rn + "S" + mch + "M"`.
    - `position += rn` (soft-clipped bases no longer align).

#### Edge Cases
- If `mch == 0`: loop doesn't execute, no change.
- If reference is missing at `position`: breaks immediately, no change.
- `rn > 3`: no replacement (only converts ≤3 leading mismatches).

---

### `combineDigSDigM()`

**Source**: [CigarModifier.java](VarDictJava/src/main/java/com/astrazeneca/vardict/modules/CigarModifier.java#L284-L358)  
**Purpose**: Adjust the boundary between a leading soft-clip and matched region by checking whether bases in the soft-clipped region actually match the reference (extend M leftward) or bases in the matched region are mismatches (extend S rightward).

#### Parameters
- `matcher`: `Matcher` — matched against `DIG_S_DIG_M` (`^(\d+)S(\d+)M`)

#### Algorithm

**Part A: Extend M leftward into S (L290–L307)**

1. Extract `mch = toInt(group(2))`, `soft = toInt(group(1))`.
2. Initialize `rn = 0`, `RN = new HashSet<Character>()` (tracks distinct bases for homopolymer detection).
3. **While** `rn < soft` AND `isHasAndEquals(reference, position - rn - 1, querySequence, soft - rn - 1)` AND `queryQuality.charAt(soft - rn - 1) - 33 > Configuration.LOWQUAL` (quality > 10):
    - Increment `rn`.
4. If `rn > 0`:
    - `mch += rn`, `soft -= rn`.
    - If `soft > 0`: rewrite as `soft + "S" + mch + "M"`.
    - Else: rewrite as `mch + "M"`.
    - `position -= rn`.
    - Reset `rn = 0`.

**Part B: Extended leftward extension with gap tolerance (L308–L327)**

5. If `soft > 0` (still have soft-clipped bases):
6. **While** `rn + 1 < soft` AND `isHasAndEquals(reference, position - rn - 2, querySequence, soft - rn - 2)` AND quality check:
    - Increment `rn`.
    - `RN.add(reference.get(position - rn - 2))` — note: uses the **post-increment** `rn` value, so adds the base at `position - rn - 2` where `rn` is already incremented.
7. Compute `rn_nt = RN.size()` — number of distinct reference bases in the extended match region.
8. If `(rn > 4 && rn_nt > 1)` OR `isHasAndEquals(reference, position - 1, querySequence, soft - 1)`:
    - `mch += rn + 1`, `soft -= rn + 1`.
    - Rewrite CIGAR using `DIG_S_DIG_M.matcher(cigarStr).replaceFirst(...)` (re-creates the matcher).
    - `position -= rn + 1`.

**Part C: Extend S rightward for mismatches (L328–L357)**

9. If `rn == 0` (no leftward extension happened in Part B):
10. Initialize `rrn = 0`, `rmch = 0`.
11. **While** `rrn < mch` AND `rn < mch`:
    - a. If reference doesn't contain `position + rrn`: **break**.
    - b. If `isHasAndNotEquals(reference, position + rrn, querySequence, soft + rrn)` (mismatch):
        - `rn = rrn + 1`, `rmch = 0`.
    - c. If `isHasAndEquals(reference, position + rrn, querySequence, soft + rrn)` (match):
        - `rmch++`.
    - d. `rrn++`.
    - e. If `rmch >= 3`: **break**.
12. If `rn > 0` AND `rn < mch`:
    - `soft += rn`, `mch -= rn`.
    - Rewrite using `DIG_S_DIG_M.matcher(cigarStr).replaceFirst(soft + "S" + mch + "M")`.
    - `position += rn`.

#### Mutable State Changes
- `this.cigarStr` — rewritten at most 3 times (Part A, B, C)
- `this.position` — decremented by Part A/B extensions, incremented by Part C

#### Edge Cases
- If `soft == 0` after Part A: Parts B and C are skipped entirely.
- `rn_nt == 1` (homopolymer): Part B only triggers if there's an exact base match at `position - 1`.
- Quality check uses `charAt(index) - 33 > LOWQUAL` where LOWQUAL=10, meaning bases with phred quality ≤10 stop the extension.

---

### `captureMisSoftly3Mismatches()`

**Source**: [CigarModifier.java](VarDictJava/src/main/java/com/astrazeneca/vardict/modules/CigarModifier.java#L364-L399)  
**Purpose**: For CIGARs ending with `M` only (no trailing S), check if the last few bases are mismatches and convert them to soft-clip.

#### Parameters
- `matcher`: `Matcher` — matched against `BEGIN_ANY_DIG_M_END` (`^(.*?)(\d+)M$`)

#### Algorithm
1. Extract `ov5 = group(1)` (prefix before final M), `mch = toInt(group(2))` (matched length).
2. Compute `refoff = position + mch` — reference offset at end of matched block.
3. Compute `rdoff = mch` — read offset at end of matched block.
4. If `ov5` is not empty:
    - `refoff += sum(globalFind(ALIGNED_LENGTH_MND, ov5))` — add all `\d+[MND]` lengths from prefix to get reference position.
    - `rdoff += sum(globalFind(SOFT_CLIPPED, ov5))` — add all `\d+[MIS]` lengths from prefix to get read position.
5. Initialize `rn = 0`, `rrn = 0`, `rmch = 0`.
6. **While** `rrn < mch` AND `rn < mch`:
    - a. If reference doesn't contain `refoff - rrn - 1`: **break**.
    - b. If `rrn < rdoff` AND mismatch at `(refoff - rrn - 1, querySequence[rdoff - rrn - 1])`:
        - `rn = rrn + 1`, `rmch = 0`.
    - c. If `rrn < rdoff` AND match:
        - `rmch++`.
    - d. `rrn++`.
    - e. If `rmch >= 3`: **break**.
7. `mch -= rn`.
8. If `rn > 0` AND `rn <= 3`:
    - Replace `DIG_M_END` with `mch + "M" + rn + "S"`.

#### Edge Cases
- Scans backward from end of M block, so `refoff - rrn - 1` walks leftward.
- Guard `rrn < rdoff` prevents index-out-of-bounds when prefix makes rdoff shorter than mch.
- Only converts ≤3 trailing mismatches (same threshold as `combineBeginDigM`).

---

### `captureMisSoftlyMS()`

**Source**: [CigarModifier.java](VarDictJava/src/main/java/com/astrazeneca/vardict/modules/CigarModifier.java#L405-L486)  
**Purpose**: For CIGARs with trailing `M-S`, extend M rightward into S if soft-clipped bases match reference, or retract M←S boundary if matched bases are actually mismatches. Mirror/complement of `combineDigSDigM` but for the 3' end.

#### Parameters
- `matcher`: `Matcher` — matched against `ANY_NUMBER_M_NUMBER_S_END` (`^(.*?)(\d+)M(\d+)S$`)

#### Algorithm

**Part A: Compute offsets (L406–L416)**

1. Extract `ov5 = group(1)` (prefix), `mch = toInt(group(2))`, `soft = toInt(group(3))`.
2. `refoff = position + mch`, `rdoff = mch`.
3. If `ov5` is not empty:
    - `refoff += sum(globalFind(ALIGNED_LENGTH_MND, ov5))`.
    - `rdoff += sum(globalFind(SOFT_CLIPPED, ov5))`.

**Part B: Extend M rightward (L418–L431)**

4. Initialize `rn = 0`, `RN = new HashSet<Character>()`.
5. **While** `rn < soft` AND `isHasAndEquals(reference, refoff + rn, querySequence, rdoff + rn)` AND quality at `rdoff + rn` > LOWQUAL:
    - Increment `rn`.
6. If `rn > 0`:
    - `mch += rn`, `soft -= rn`.
    - Rewrite using `DIG_M_DIG_S_END` replacer.
    - Reset `rn = 0`.

**Part C: Extended rightward with gap tolerance (L433–L451)**

7. If `soft > 0`:
8. **While** `rn + 1 < soft` AND match at `(refoff + rn + 1, rdoff + rn + 1)` AND quality check:
    - Increment `rn`.
    - `RN.add(reference.get(refoff + rn + 1))` — **post-increment** rn, so this reads `refoff + rn + 1` with already-incremented rn.
9. `rn_nt = RN.size()`.
10. If `rn > 4 && rn_nt > 1` (not homopolymer):
    - `mch += rn + 1`, `soft -= rn + 1`.
    - Rewrite using `DIG_M_DIG_S_END`.

**Part D: Retract M into S for mismatches (L452–L480)**

11. If `rn == 0`:
12. Initialize `rrn = 0`, `rmch = 0`.
13. **While** `rrn < mch` AND `rn < mch`:
    - a. If reference doesn't contain `refoff - rrn - 1`: **break**.
    - b. If `rrn < rdoff` AND mismatch `(refoff - rrn - 1, querySequence[rdoff - rrn - 1])`:
        - `rn = rrn + 1`, `rmch = 0`.
    - c. If `rrn < rdoff` AND match:
        - `rmch++`.
    - d. `rrn++`.
    - e. If `rmch >= 3`: **break**.
14. If `rn > 0` AND `rn < mch`:
    - `soft += rn`, `mch -= rn`.
    - Rewrite using `DIG_M_DIG_S_END`.

#### Edge Cases
- `RN.add(reference.get(refoff + rn + 1))` with post-incremented `rn` is a subtle indexing pattern — the Perl comment says "to my mind in perl value increasing here ? not adding el", confirming this is a ported-from-Perl quirk.
- Part D scans backward from the M-S boundary, not forward from the end of S.

---

### `combineToCloseToOne()`

**Source**: [CigarModifier.java](VarDictJava/src/main/java/com/astrazeneca/vardict/modules/CigarModifier.java#L494-L527)  
**Purpose**: Merge an I-M-D/I complex where the internal M is ≤15bp into a single D+I.

#### Parameters
- `matcher`: `Matcher` — matched against `NOTDIG_DIG_I_DIG_M_DIG_DI_DIGI` (`(\D)(\d+)I(\d+)M(\d+)([DI])(\d+I)?`)
- `flag`: `boolean`

#### Algorithm
1. Extract `op = group(5)`, `g2 = toInt(group(2))`, `g3 = toInt(group(3))`, `g4 = toInt(group(4))`.
2. If `g3 <= 15` (internal match ≤15bp):
    - a. `dlen = g3` (deletion length = matched length).
    - b. `ilen = g2 + g3` (insertion length = first I + matched).
    - c. If `op == "I"`: `ilen += g4`.
    - d. If `op == "D"`: `dlen += g4`.
        - If `groupCount() > 5` and `group(6) != null`: `ilen += toInt(group(6).substring(0, group(6).length() - 1))` (strip trailing "I" char to get the number).
    - e. Re-match using `DIG_I_DIG_M_DIG_DI_DIGI` pattern (without the leading non-digit guard).
    - f. Replace with `dlen + "D" + ilen + "I"`.
    - g. Set `flag = true`.
3. Return `flag`.

#### Edge Cases
- `group(6)` can be null even if `groupCount() > 5` (optional group with `?`).
- The non-digit guard `group(1)` excludes "D" and "H" at the call site (L215), preventing this from matching D-I-M-D and H-I-M-D patterns.
- Replacement uses a **different** pattern (`DIG_I_DIG_M_DIG_DI_DIGI`) than the matcher (`NOTDIG_DIG_I_DIG_M_DIG_DI_DIGI`) — intentionally drops the leading non-digit from the replacement target.

---

### `combineToCloseToCorrect()`

**Source**: [CigarModifier.java](VarDictJava/src/main/java/com/astrazeneca/vardict/modules/CigarModifier.java#L535-L569)  
**Purpose**: Merge a D-M-D/I complex where the internal M is ≤15bp into a single D+I.

#### Parameters
- `matcher`: `Matcher` — matched against `DIG_D_DIG_M_DIG_DI_DIGI` (`(\d+)D(\d+)M(\d+)([DI])(\d+I)?`)
- `flag`: `boolean`

#### Algorithm
1. Extract `g1 = toInt(group(1))`, `g2 = toInt(group(2))`, `g3 = toInt(group(3))`.
2. If `g2 <= 15`:
    - a. Extract `op = group(4)`.
    - b. `dlen = g1 + g2` (both deletions + matched).
    - c. `ilen = g2` (insertion = matched length).
    - d. If `op == "I"`: `ilen += g3`.
    - e. If `op == "D"`: `dlen += g3`.
        - If `groupCount() > 4` and `group(5) != null`: `ilen += toInt(group(5).substring(0, group(5).length() - 1))`.
    - f. Replace first match in `cigarStr` with `dlen + "D" + ilen + "I"`.
    - g. `flag = true`.
3. Return `flag`.

---

### `threeIndels()`

**Source**: [CigarModifier.java](VarDictJava/src/main/java/com/astrazeneca/vardict/modules/CigarModifier.java#L578-L642)  
**Purpose**: Merge a 3-indel complex (M-D/I-M-D/I-M-D/I-M) into a simplified D+I+M form.

#### Parameters
- `matcher`: `jregex.Matcher` — matched against `threeIndelsPattern` (`^(.*?)(\d+)M(\d+)([DI])(\d+)M(\d+)([DI])(\d+)M(\d+)([DI])(\d+)M`)
- `flag`: `boolean`

#### Algorithm
1. Compute `tslen = group(5) + group(8)` (two internal match lengths).
2. If `group(4) == "I"`: `tslen += group(3)`.
3. If `group(7) == "I"`: `tslen += group(6)`.
4. If `group(10) == "I"`: `tslen += group(9)`.
5. Compute `dlen = group(5) + group(8)` (same base as tslen).
6. If `group(4) == "D"`: `dlen += group(3)`.
7. If `group(7) == "D"`: `dlen += group(6)`.
8. If `group(10) == "D"`: `dlen += group(9)`.
9. Compute `mid = group(5) + group(8)` — total internal match length.
10. Extract `ov5 = group(1)` (prefix), `refoff = position + group(2)`, `rdoff = group(2)`, `RDOFF = group(2)`, `rm = group(11)`.
11. If `ov5` not empty:
    - `refoff += sum(globalFind(ALIGNED_LENGTH_MND, ov5))`.
    - `rdoff += sum(globalFind(SOFT_CLIPPED, ov5))`.
12. Initialize `rn = 0`.
13. **While** `rdoff + rn < querySequence.length()` AND `isHasAndEquals(querySequence.charAt(rdoff + rn), reference, refoff + rn)`:
    - Increment `rn`.
14. `RDOFF += rn`, `dlen -= rn`, `tslen -= rn`.
15. Build `newCigarStr = RDOFF + "M"`.
16. **If `tslen <= 0`**:
    - `dlen -= tslen`, `rm += tslen`.
    - If `dlen == 0`: `RDOFF = RDOFF + rm`, `newCigarStr = RDOFF + "M"`.
    - If `dlen < 0`: `tslen = -dlen`, `rm += dlen`.
        - If `rm < 0`: `RDOFF = RDOFF + rm`, `newCigarStr = RDOFF + "M" + tslen + "I"`.
        - Else: `newCigarStr += tslen + "I" + rm + "M"`.
    - Else: `newCigarStr += dlen + "D" + rm + "M"`.
17. **If `tslen > 0`**:
    - If `dlen == 0`: `newCigarStr += tslen + "I" + rm + "M"`.
    - If `dlen < 0`: `rm += dlen`, `newCigarStr += tslen + "I" + rm + "M"`.
    - Else: `newCigarStr += dlen + "D" + tslen + "I" + rm + "M"`.
18. If `mid <= 15`: replace using `DIGM_D_DI_DIGM_D_DI_DIGM_DI_DIGM`, set `flag = true`.
19. Return `flag`.

#### Edge Cases
- `tslen` can go negative after the reference-match extension (step 14). The complex branching in step 16 handles this.
- `dlen` can also go negative, leading to the insertion-only or pure-match forms.
- `rm` can go negative — handled by folding it back into RDOFF.

---

### `threeDeletions()`

**Source**: [CigarModifier.java](VarDictJava/src/main/java/com/astrazeneca/vardict/modules/CigarModifier.java#L651-L701)  
**Purpose**: Merge M-D-M-D-M-D-M into simplified D+optional-I+M form.

#### Parameters
- `matcher`: `jregex.Matcher` — matched against `threeDeletionsPattern` (`^(.*?)(\d+)M(\d+)D(\d+)M(\d+)D(\d+)M(\d+)D(\d+)M`)
- `flag`: `boolean`

#### Algorithm
1. `tslen = group(4) + group(6)` (two internal match lengths = insertion component).
2. `dlen = group(3) + group(4) + group(5) + group(6) + group(7)` (all three deletions + two internal matches).
3. `mid = group(4) + group(6)` (total internal match).
4. Extract `ov5 = group(1)`, `refoff = position + group(2)`, `rdoff = group(2)`, `RDOFF = group(2)`, `rm = group(8)`.
5. If `ov5` not empty: adjust `refoff` and `rdoff` by prefix lengths.
6. `rn = 0`. **While** `rdoff + rn < querySequence.length()` AND match:
    - Increment `rn`.
7. `RDOFF += rn`, `dlen -= rn`, `tslen -= rn`.
8. Build `newCigarStr = RDOFF + "M"`.
9. If `tslen <= 0`:
    - `dlen -= tslen`, `rm += tslen`.
    - `newCigarStr += dlen + "D" + rm + "M"`.
10. Else:
    - `newCigarStr += dlen + "D" + tslen + "I" + rm + "M"`.
11. If `mid <= 15`: replace using `DM_DD_DM_DD_DM_DD_DM`, set `flag = true`.
12. Return `flag`.

#### Edge Cases
- Simpler branching than `threeIndels()` — only two cases (tslen≤0 or >0).
- No `dlen < 0` or `rm < 0` guards (unlike threeIndels), because all three ops are deletions so dlen is always large.

---

### `twoDeletionsInsertionToComplex()`

**Source**: [CigarModifier.java](VarDictJava/src/main/java/com/astrazeneca/vardict/modules/CigarModifier.java#L710-L762)  
**Purpose**: Merge the specific pattern M-D-M-I-M-D-M into a simplified D+I+M form.

#### Parameters
- `matcher`: `jregex.Matcher` — matched against `D_M_D_DD_M_D_I_D_M_D_DD` (`^(.*?)(\d+)M(\d+)D(\d+)M(\d+)I(\d+)M(\d+)D(\d+)M`)
- `flag`: `boolean`

#### Algorithm
1. `tslen = group(4) + group(5) + group(6)` (first internal match + insertion + second internal match).
2. `dlen = group(3) + group(4) + group(6) + group(7)` (both deletions + both internal matches).
3. `mid = group(4) + group(6)` (internal matches only).
4. Extract `ov5 = group(1)`, `refoff = position + group(2)`, `rdoff = group(2)`, `RDOFF = group(2)`, `rm = group(8)`.
5. If `ov5` not empty: adjust offsets.
6. `rn = 0`. Match extension loop.
7. `RDOFF += rn`, `dlen -= rn`, `tslen -= rn`.
8. Build `newCigarStr = RDOFF + "M"`.
9. If `tslen <= 0`:
    - `dlen -= tslen`, `rm += tslen`.
    - `newCigarStr += dlen + "D" + rm + "M"`.
10. Else:
    - `newCigarStr += dlen + "D" + tslen + "I" + rm + "M"`.
11. If `mid <= 15`: replace using `D_M_D_DD_M_D_I_D_M_D_DD_prim`, set `flag = true`.
12. Return `flag`.

---

### `beginDigitMNumberIorDNumberM()`

**Source**: [CigarModifier.java](VarDictJava/src/main/java/com/astrazeneca/vardict/modules/CigarModifier.java#L771-L787)  
**Purpose**: Convert a leading short match followed by an indel into soft-clip. The digit-M (1–9 bases) plus the indel are converted to S, then mismatches at the start of the remaining M are also absorbed.

#### Parameters
- `matcher`: `jregex.Matcher` — matched against `BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M` (`^(\d)M(\d+)([ID])(\d+)M`)

#### Algorithm
1. `tmid = toInt(group(1))` — length of short leading match (single digit = 1–9).
2. `mlen = toInt(group(4))` — length of longer matched block after indel.
3. Initialize `tn = 0`.
4. `tslen = tmid + (group(3) == "I" ? group(2) : 0)` — soft-clip length starts as short match + insertion length (or just short match for deletion).
5. `position += tmid + (group(3) == "D" ? group(2) : 0)` — advance position by short match + deletion length.
6. **While** `tn < mlen` AND `isHasAndNotEquals(querySequence.charAt(tslen + tn), reference, position + tn)` (mismatch):
    - Increment `tn`.
7. `tslen += tn`, `mlen -= tn`, `position += tn`.
8. Replace using `BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M_` (`^\dM\d+[ID]\d+M`) with `tslen + "S" + mlen + "M"`.
9. Return `true` (always sets flag).

#### Edge Cases
- The replacement pattern `BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M_` uses `java.util.regex.Pattern`, while the matcher uses `jregex.Pattern` — different regex engines for match vs. replace.
- The mismatch extension loop (step 6) can consume all of `mlen`, resulting in `0M` in the CIGAR (which is valid but unusual).

---

## Cross-Module Dependencies

### Outbound Calls

| Target | Method/Field | Usage |
|--------|-------------|-------|
| `GlobalReadOnlyScope` | `instance().conf.chimeric` | Controls chimeric soft-clip removal |
| `GlobalReadOnlyScope` | `instance().conf.y` | Debug logging flag |
| `GlobalReadOnlyScope` | `instance().conf.performLocalRealignment` | Checked by caller (CigarParser), not by CigarModifier itself |
| `Configuration` | `SEED_2` (=12) | Minimum soft-clip length for chimeric detection |
| `Configuration` | `LOWQUAL` (=10) | Minimum quality for soft-clip realignment |
| `Reference` | `referenceSequences` (`Map<Integer,Character>`) | Reference base lookups |
| `Reference` | `seed` (`Map<String,List<Integer>>`) | Seed map for chimeric detection |
| `Utils` | `toInt()`, `substr()`, `reverse()`, `complement()` | String/number utilities |
| `Utils` | `sum()`, `globalFind()` | Prefix offset computation |
| `Utils` | `printExceptionAndContinue()` | Error handling |
| `VariationUtils` | `isHasAndEquals(Map,int,String,int)` | Reference-vs-read base comparison |
| `VariationUtils` | `isHasAndEquals(char,Map,int)` | Reference-vs-char comparison |
| `VariationUtils` | `isHasAndNotEquals(Map,int,String,int)` | Reference-vs-read mismatch check |
| `VariationUtils` | `isHasAndNotEquals(char,Map,int)` | Reference-vs-char mismatch check |
| `Patterns` | 30+ regex constants | All CIGAR matching/replacement |

### Inbound Callers

| Caller | Method | Context |
|--------|--------|---------|
| `CigarParser` | `cigarParser()` L286-L300 | Created and called when `performLocalRealignment` is true |

## Data Structures Read/Written

| Structure | Java Type | R/W | Notes |
|-----------|-----------|-----|-------|
| `position` | `int` (field) | R/W | Adjusted by boundary deletion stripping, short-M absorption, soft-clip realignment |
| `cigarStr` | `String` (field) | R/W | Repeatedly rewritten by regex. Final value returned. |
| `querySequence` | `String` (field) | R/W | Trimmed **only** by chimeric detection (phases 3a/3b). Indexed by all soft-clip methods. |
| `queryQuality` | `String` (field) | R/W | Trimmed **only** by chimeric detection. Quality-checked in soft-clip methods. |
| `reference` | `Map<Integer,Character>` (field) | R | Reference base lookups via `isHasAndEquals`/`isHasAndNotEquals` |
| `ref.seed` | `Map<String,List<Integer>>` | R | Chimeric seed lookup |
| `originalCigar` | `String` (field) | R | Only used in error message |
| `indel` | `int` (field) | R | Controls whether the while-loop executes at all |
| `maxReadLength` | `int` (field) | R | Distance threshold for chimeric detection |
| `region` | `Region` (field) | R | Only used in error message |

## Known Parity Traps

1. **Clause ordering in the while-loop is critical.** The matchers are checked sequentially: `BEGIN_NUMBER_S_NUMBER_IorD` before `NUMBER_IorD_NUMBER_S_END` before `BEGIN_NUMBER_S_NUMBER_M_NUMBER_IorD`, etc. A CIGAR can match multiple patterns — the first one wins. Java applies trailing I/D-M-S collapse before front-edge short M-I/D-M (see repo memory `cigar_modifier_edge_order_parity_20260310.md`). Live regression: `4M59I10M28S` must normalize to `101M`.

2. **`jregex` vs `java.util.regex` mixing.** Some patterns use `jregex.Pattern` (for replacement via `Replacer`) and others use `java.util.regex.Pattern` (for `Matcher.replaceFirst`). Critically, `beginDigitMNumberIorDNumberM()` matches with jregex `BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M` but replaces with java.util `BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M_`. Rust must use the same regex for both operations.

3. **`combineToCloseToOne` uses a different replacement pattern than its match pattern.** Matched via `NOTDIG_DIG_I_DIG_M_DIG_DI_DIGI` (includes leading non-digit) but replaced via `DIG_I_DIG_M_DIG_DI_DIGI` (excludes leading non-digit). This means the leading character (e.g., "M", "S") is preserved in the output.

4. **`NUMBER_IorD_DIGIT_M_END` vs `NUMBER_IorD_NUMBER_M_END` mismatch.** In step 4f (L185–L191), the match uses `(\d+)([ID])(\d)M$` (single digit M) but the replacement uses `(\d+)([ID])(\d+)M$` (multi-digit M). The replacement pattern is strictly broader and will always match whatever the match pattern matched. This is correct but must be replicated exactly.

5. **`RN.add(reference.get(pos))` uses post-incremented `rn`.** In `combineDigSDigM` Part B (L314) and `captureMisSoftlyMS` Part C (L440), the `rn` is incremented before the `RN.add()` call. The actual reference position computed uses the new `rn` value. Rust must reproduce this exact indexing.

6. **Homopolymer guard `rn_nt > 1`.** In both `combineDigSDigM` (Part B) and `captureMisSoftlyMS` (Part C), the extension only triggers if `rn > 4 && rn_nt > 1` (more than 4 matches AND more than 1 distinct base). For homopolymer runs of 5+ bases at clip boundaries, the extension is suppressed. In `combineDigSDigM` Part B, there's an additional OR condition `isHasAndEquals(reference, position - 1, querySequence, soft - 1)` that allows extension even for homopolymers if the immediate boundary base matches.

7. **`substr()` with negative indices.** `Utils.substr(string, -n, n)` computes `begin = string.length() + (-n)` and takes `n` chars from there. Used in 3' chimeric seed extraction (L113). Rust equivalent must handle negative-begin semantics identically.

8. **`complement(reverse(sseq))` order.** The 5' chimeric path (L95) takes `substr(querySequence, 0, cigarElement)` → reverse → complement. The 3' path (L113) takes `substr(querySequence, -cigarElement, cigarElement)` → reverse → complement. Then for the 5' seed, `sequence.substring(0, SEED_2)` takes the **first** 12 chars. For 3', `substr(sequence, -SEED_2, SEED_2)` takes the **last** 12 chars. Order matters.

9. **The 5' and 3' chimeric checks are `if/else if`** — only one fires per invocation. If a CIGAR has both ≥10-digit leading and trailing soft clips, only the 5' is checked. (The `else if` at L108.)

10. **`isHasAndEquals()` returns false when key is missing from reference map.** This silently terminates extension loops rather than throwing exceptions. If a Rust implementation uses a different reference lookup that panics on missing keys, parity will break.

11. **Exception swallowing.** On any exception during `modifyCigar()`, the catch block logs the error and returns the **partially-modified** `ModifiedCigar`. This means a crash mid-way through Phase 4 returns a half-transformed CIGAR. Rust must either replicate this behavior or prove each sub-method is infallible.

12. **`flag` initialization and the while-loop.** `flag` starts `true`, but the while condition is `flag && indel > 0`. On first iteration, `flag` is set to `false`. Each sub-check may set it back to `true`. If **no** sub-check matches, the loop exits. If `indel == 0`, the entire Phase 4 is skipped. `indel` is never modified inside the loop — only `flag` controls iteration.

13. **`combineBeginDigM` has range guard `rn <= 3`** while `captureMisSoftly3Mismatches` also uses `rn <= 3`. But `combineDigSDigM` Part C and `captureMisSoftlyMS` Part D use `rn > 0 && rn < mch` (no upper bound). The asymmetry is intentional.

14. **`globalFind(ALIGNED_LENGTH_MND, ov5)` uses `\d+[MND]` pattern** which includes N (skipped reference). `globalFind(SOFT_CLIPPED, ov5)` uses `\d+[MIS]`. These compute reference offset and read offset respectively from the CIGAR prefix. If the Rust implementation uses different patterns for offset computation, all prefix-aware methods (captureMisSoftly3Mismatches, captureMisSoftlyMS, threeIndels, threeDeletions, twoDeletionsInsertionToComplex) will produce wrong offsets.

15. **`DIG_S_DIG_M.matcher(cigarStr).replaceFirst(...)` re-creates the matcher each time** in `combineDigSDigM`. This is not the passed-in `matcher` parameter — it's a fresh match against the current `cigarStr` state, which may have been modified by Part A. Rust must re-apply the regex to the current cigarStr, not cache the original match.

16. **`toInt()` delegates to `Integer.parseInt()`**, which throws `NumberFormatException` on non-numeric strings. All regex groups are expected to be numeric, but if a regex somehow matches incorrectly, the exception is caught by the outer try/catch in `modifyCigar()`. Rust should use infallible parsing or replicate the exception-swallowing.

17. **Quality check boundary**: `queryQuality.charAt(index) - 33 > Configuration.LOWQUAL` means the quality threshold is `> 10` (phred > 10, ASCII > 43). This is `>` not `>=`, so phred 10 (ASCII 43) does NOT pass the quality check.

18. **`referenceSeedMap` is a `HashMap<String, List<Integer>>`** (not LinkedHashMap). The `positions.size() == 1` check requires the seed to be unique. If the seed map has been populated differently in Java vs Rust (e.g., different hashing causing different collision patterns), the `containsKey`/`get` results must still be identical.

19. **The `while` loop in `captureMisSoftlyMS` Part D uses variable `rn` which was last set to 0 by the `if (rn == 0)` guard**, but the loop condition checks `rn < mch` where `rn` is also mutated inside the loop body (set to `rrn + 1`). The dual use of `rn` as both loop guard and mismatch counter is a Perl legacy pattern. `rrn` advances position and `rn` records the last mismatch boundary.

20. **Adjacent D+D and I+I merging (steps 39–42) uses `java.util.regex.Matcher.replaceFirst()`** which only replaces the first occurrence. If the CIGAR contains `5D3D7D`, only the first pair is merged to `8D7D`, requiring another loop iteration for `8D7D` → `15D`. The `flag` mechanism handles this iterative convergence.
