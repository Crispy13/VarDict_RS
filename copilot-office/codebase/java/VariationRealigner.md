# VariationRealigner — Full Module Analysis

**Source**: `VarDictJava/src/main/java/com/astrazeneca/vardict/modules/VariationRealigner.java`
**LOC**: 2,953
**Rust counterpart**: `src/mods/variant_realigner.rs`
**Risk**: HIGH — local realignment, position adjustment, consensus-based variant discovery
**Pipeline Stage**: After CigarParser, before StructuralVariantsProcessor
**Status**: complete

## Overview

VariationRealigner is the second stage of the VarDict variant calling pipeline. It takes the raw variant maps (insertions, non-insertions, soft clips, SV structures) produced by CigarParser and performs local realignment. Its responsibilities are:

1. **SV Filtering** (`filterAllSVStructures`): Clusters SV mates by position proximity, filters false positives by discordant/count ratio, populates `SOFTP2SV` lookup.
2. **MNP Adjustment** (`adjustMNP`): Merges partial MNP sub-variants (left/right fragments) back into their parent MNP, absorbs matching soft clips.
3. **Short Indel Realignment** (`realigndel`, `realignins`): Absorbs nearby mismatches and soft-clipped reads into known short deletions/insertions based on flanking sequence matching.
4. **Large Deletion Discovery** (`realignlgdel`): Discovers large deletions from unpaired 5'/3' soft clips by searching for breakpoints in the reference.
5. **Large Insertion Discovery** (`realignlgins30`, `realignlgins`): Pairs 3'/5' soft clips for ≥30bp insertions, or discovers large insertions/DUPs from single soft clips via `findMatch`.

The module **mutates shared variant maps in place** — adding counts to target variants, removing consumed variants, and marking soft clips as `used`. All methods operate on instance fields set by `initFromScope()`.

## Method Inventory

| Method | Lines | Analyzed? | Summary |
|--------|-------|-----------|---------|
| `process()` | L73–L114 | **yes** | Entry point: orchestrates filter→adjustMNP→realign pipeline |
| `initFromScope()` | L316–L338 | **yes** | Copies all fields from Scope into instance fields |
| `filterAllSVStructures()` | L353–L377 | **yes** | Filters all SV lists + sorts SOFTP2SV by varsCount |
| `filterSV()` | L383–L425 | **yes** | For each SV: clusters mates, applies disc/cnt filter, populates SOFTP2SV |
| `checkCluster()` | L432–L482 | **yes** | Clusters mates by position proximity, returns dominant cluster |
| `adjustMNP()` | L487–L569 | **yes** | Merges partial MNP sub-variants and matching soft clips |
| `realignIndels()` | L571–L586 | **yes** | Dispatcher: calls realigndel→realignins→realignlgdel→realignlgins30→realignlgins |
| `realigndel()` | L593–L858 | **yes** | Absorbs mismatches/soft clips into short deletions (two passes) |
| `realignins()` | L864–L1170 | **yes** | Absorbs mismatches/soft clips into short insertions (two passes) |
| `realignlgdel()` | L1177–L1594 | **yes** | Discovers large deletions from unpaired soft clips (5' pass then 3' pass) |
| `realignlgins30()` | L1599–L1799 | **yes** | Pairs 3'/5' soft clips for ≥30bp insertions |
| `realignlgins()` | L1805–L2150 | **yes** | Discovers large insertions/DUPs from single soft clips (5' then 3') |
| `fillAndSortTmp()` | L2155–L2196 | **yes** | Flattens indel count map into sorted list |
| `SortPositionDescription` (inner class) | L2190–L2207 | **yes** | Data tuple: (position, descriptionString, count) |
| `find35match()` | L2214–L2260 | **yes** | Finds overlapping match between 5'/3' sequences |
| `noPassingReads()` | L2270–L2306 | **yes** | Checks for reads spanning a deletion gap |
| `ismatch()` (2-overload) | L2313–L2345 | **yes** | Compares sequences allowing ≤MM mismatches |
| `islowcomplexseq()` | L2353–L2385 | **yes** | Returns true if >75% single base or <3 distinct bases |
| `count()` | L2392–L2401 | **yes** | Counts char occurrences in a string |
| `adjInsPos()` | L2408–L2423 | **yes** | Left-aligns insertion position |
| `findbi()` | L2431–L2538 | **yes** | Finds insertion breakpoint |
| `findbp()` | L2545–L2600 | **yes** | Finds deletion breakpoint |
| `adjRefCnt()` | L2608–L2647 | **yes** | Adjusts reference counts by position-derived factor |
| `adjRefFactor()` | L2653–L2695 | **yes** | Adjusts reference by multiplicative factor |
| `addVarFactor()` | L2701–L2722 | **yes** | Adds counts by factor |
| `findMM5()` | L2730–L2780 | **yes** | Walks 5' direction for mismatches from variant position |
| `findMM3()` | L2788–L2840 | **yes** | Walks 3' direction for mismatches from variant position |
| `MismatchResult` (inner class) | L2842–L2870 | **yes** | Return type for findMM3/findMM5 |
| `Mismatch` (inner class) | L2872–L2890 | **yes** | Mismatch record: (sequence, position, end) |
| `ismatchref()` (2-overload) | L2890–L2929 | **yes** | Checks sequence matches reference |
| `rmCnt()` | L2936–L2949 | **yes** | Subtracts variant counts |
| `COMP2` (static comparator) | L340–L348 | **yes** | Sort by desc varsCount, asc position |
| `COMP3` (static comparator) | L350–L358 | **yes** | Sort by desc count field, asc position |
| `VariationRealignerJsonlWriter` (inner class) | L123–L310 | **no** | Debug/diagnostic JSONL output — not parity-relevant |

## Method Analyses

---

### `process()`

**Source**: L73–L114
**Purpose**: Entry point implementing `Module<VariationData, RealignedVariationData>`. Orchestrates the entire realignment pipeline.
**Called By**: `SAMFileParser.parseSAM()` via `Module.process()`

**Algorithm (Step-by-Step)**:
1. Call `initFromScope(scope)` to copy all fields.
2. Create `CurrentSegment CURSEG` from `region.chr`, `region.start`, `region.end`.
3. If `!instance().conf.disableSV`: call `filterAllSVStructures()`.
4. Call `adjustMNP()`.
5. If `instance().conf.y`: print timing.
6. If `instance().conf.performLocalRealignment`: call `realignIndels()`.
7. Construct `RealignedVariationData` from all instance fields plus `CURSEG`, `SOFTP2SV`, `scope`.
8. Check env var `VARDICT_VARIANT_REALIGNER_JSONL` — if set, write diagnostic JSONL (non-parity-critical).
9. Return `new Scope<>(scope, realigned)`.

**Parity Warnings**:
- The JSONL writer is diagnostic-only and should NOT be ported for parity. It's gated by an env var.

---

### `initFromScope()`

**Source**: L316–L338
**Purpose**: Copies all fields from the `Scope` into instance variables for convenient access.
**Called By**: `process()`

**Algorithm (Step-by-Step)**:
1. `this.region = scope.region`
2. Copy all variant maps: `nonInsertionVariants`, `insertionVariants`, `positionToInsertionCount`, `positionToDeletionsCount`, `refCoverage`, `softClips5End`, `softClips3End`.
3. `this.reference = scope.regionRef`
4. `this.referenceResource = scope.referenceResource`
5. `this.chr = getChrName(scope.region)` — external call to `RecordPreprocessor.getChrName`
6. `this.maxReadLength = scope.maxReadLength`
7. `this.bams = scope.bam != null ? scope.bam.split(":") : null` — **PARITY CRITICAL**: null check then split on `:`
8. `this.bam = scope.bam`
9. `this.mnp = scope.data.mnp`
10. `this.splice = scope.splice`
11. `this.svStructures = scope.data.svStructures`
12. `this.duprate = scope.data.duprate`
13. `this.variantPrinter = scope.out`

**Parity Warnings**:
- `bams` is set to `null` when `scope.bam` is null. This affects `noPassingReads()` guard checks downstream.
- `bam.split(":")` produces `String[]` — the `:` delimiter is BAM list separator, not path component.

---

### `filterAllSVStructures()`

**Source**: L353–L377
**Purpose**: Filters all SV structure lists to remove false positives, then sorts `SOFTP2SV` entries by descending `varsCount`.
**Called By**: `process()` (when `!disableSV`)

**Algorithm (Step-by-Step)**:
1. Call `filterSV()` on 8 lists: `svfinv3`, `svrinv3`, `svfinv5`, `svrinv5`, `svfdel`, `svrdel`, `svfdup`, `svrdup`.
2. For each entry in `svStructures.svffus` (forward fusions, keyed by chr): call `filterSV()` on value list.
3. For each entry in `svStructures.svrfus` (reverse fusions): call `filterSV()` on value list.
4. For each entry in `SOFTP2SV`:
   - Sort the `List<Sclip>` by **descending** `sclip.varsCount`.
   - Put sorted list back into map.

**Parity Warnings**:
- The SOFTP2SV sort is by `varsCount` only (no tiebreaker). If two Sclips have equal varsCount, the sort is non-deterministic.
- `svffus`/`svrfus` are `Map<String, List<Sclip>>` keyed by chromosome name. Iteration order of the outer map doesn't affect results since `filterSV` processes each list independently.

---

### `filterSV()`

**Source**: L383–L425
**Purpose**: For each SV in a list, clusters its mates, applies disc/cnt filtering, populates `SOFTP2SV`.
**Called By**: `filterAllSVStructures()`

**Algorithm (Step-by-Step)**:
1. For each `sv` in `svList_sva`:
   a. Call `checkCluster(sv.mates, maxReadLength)` → `cluster`.
   b. If `cluster.mateStart_ms != 0`:
      - Copy cluster fields into `sv`: `mstart`, `mend`, `varsCount`, `mlen`, `start`, `end`, `meanPosition`, `meanQuality`, `meanMappingQuality`, `numberOfMismatches`.
   c. Else: `sv.used = true` (mark unusable — zero cluster means no dominant cluster found).
   d. **Disc/cnt filter**: if `sv.disc != 0` AND `sv.varsCount / (double)sv.disc < 0.5`:
      - If NOT (`sv.varsCount / (double)sv.disc >= 0.35 AND sv.varsCount >= 5`): `sv.used = true`.
   e. Sort `sv.soft` entries (a `Map<Integer, Integer>`, soft clip position → count) by descending value.
   f. `sv.softp = soft.size() > 0 ? soft.get(0).getKey() : 0` — take position of highest-count soft clip.
   g. If `sv.softp != 0`: add `sv` to `SOFTP2SV.getOrDefault(sv.softp, new ArrayList<>())`.

**Parity Warnings**:
- The disc/cnt filter has a nested negation: `!(cnt/disc >= 0.35 && cnt >= 5)`. In Rust, ensure the `&&` grouping is preserved.
- `sv.soft` is a `Map<Integer, Integer>` sorted by value desc. If two entries have the same count, the position chosen for `softp` is non-deterministic.

---

### `checkCluster()`

**Source**: L432–L482
**Purpose**: Clusters mates by position proximity. Returns the dominant cluster (≥60% of mates) with averaged statistics.
**Called By**: `filterSV()`

**Algorithm (Step-by-Step)**:
1. Sort `mates` by `mate.mateStart_ms` ascending.
2. Initialize `clusters` list with one cluster from `firstMate`: `(0, ms, me, s, e)`.
3. Set `cur = 0`.
4. For each `mate_m` in `mates`:
   a. Get `currentCluster = clusters.get(cur)`.
   b. If `mate_m.mateStart_ms - currentCluster.mateEnd_me > MINSVCDIST * rlen`:
      - `cur++`; add new cluster at index `cur`.
   c. `currentCluster.cnt++`
   d. Accumulate `mateLength_mlen`, update `mateEnd_me`, `start_s`, `end_e`, `pmean_rp`, `qmean_q`, `Qmean_Q`, `nm`.
5. Sort `clusters` by descending `cnt`.
6. If `firstCluster.cnt / (double)mates.size() >= 0.60`:
   - Return cluster with averaged `mateLength_mlen` (divided by `cnt`), all accumulated stats.
7. Else: return zero cluster `(0,0,0,0,0,0,0,0.0,0,0)`.

**Parity Warnings**:
- The first mate is **both** used to initialize cluster 0 **and** processed in the loop (incrementing cnt). So `cnt` starts at 1 after the first iteration.
- `mateLength_mlen/firstCluster.cnt` — integer division. Rust must use integer division.
- The 0.60 threshold uses `double` division. Rust must use `f64`.
- The zero-cluster return has `mateStart_ms = 0`, which is the sentinel checked by `filterSV`.

---

### `adjustMNP()`

**Source**: L487–L569
**Purpose**: Merges partial MNP sub-variants (left/right fragments) back into their parent MNP, and absorbs matching soft clips.
**Called By**: `process()`

**Algorithm (Step-by-Step)**:
1. `tmp = fillAndSortTmp(mnp)` — flatten MNP count map into sorted list.
2. For each `tpl` in `tmp`:
   a. `position = tpl.position`, `vn = tpl.descriptionString`.
   b. `varsOnPosition = nonInsertionVariants.get(position)` — if null, `continue`.
   c. `vref = varsOnPosition.get(vn)` — if null, `continue` (already consumed).
   d. `mnt = vn.replaceFirst("&", "")` — strip the reference marker.
   e. For `i = 0` to `mnt.length() - 2` (inner loop):
      - **Left fragment**: `left = substr(mnt, 0, i + 1)`. If `left.length() > 1`: insert `&` at index 1.
      - **Right fragment**: `right = substr(mnt, -(mnt.length() - i - 1))`. If `right.length() > 1`: insert `&` at index 1.
      - **Left check**: `tref = varsOnPosition.get(left)`.
        - If `tref != null` AND `tref.varsCount <= 0`: `continue`. (**ASYMMETRY**: uses `<= 0`)
        - If `tref.varsCount < vref.varsCount` AND `tref.meanPosition / tref.varsCount <= i + 1`:
          - `adjCnt(vref, tref)`, `varsOnPosition.remove(left)`.
      - **Right check**: `tref = nonInsertionVariants.get(position + i + 1).get(right)`.
        - If `tref != null` AND `tref.varsCount < 0`: `continue`. (**ASYMMETRY**: uses `< 0`, not `<= 0`)
        - If `tref.varsCount < vref.varsCount`:
          - `adjCnt(vref, tref)`, `incCnt(refCoverage, position, tref.varsCount)`.
          - `nonInsertionVariants.get(position + i + 1).remove(right)`.
   f. **3' soft clip check**: If `softClips3End.containsKey(position)`:
      - Get `sc3v`. If `!sc3v.used`:
        - `seq = findconseq(sc3v, 0)`.
        - If `seq.startsWith(mnt)`:
          - If seq length == mnt length OR `ismatchref(seq.substring(mnt.length()), ref, position + mnt.length(), 1)`:
            - `adjCnt(vref, sc3v)`, `incCnt(refCoverage, position, sc3v.varsCount)`, `sc3v.used = true`.
   g. **5' soft clip check**: If `softClips5End.containsKey(position + mnt.length())`:
      - Get `sc5v`. If `!sc5v.used`:
        - `seq = findconseq(sc5v, 0)`. Reverse it.
        - If `seq.endsWith(mnt)`:
          - If seq length == mnt length OR `ismatchref(seq.substring(0, seq.length() - mnt.length()), ref, position - 1, -1)`:
            - `adjCnt(vref, sc5v)`, `incCnt(refCoverage, position, sc5v.varsCount)`, `sc5v.used = true`.

**Parity Warnings**:
- **CRITICAL ASYMMETRY**: Left uses `tref.varsCount <= 0` (skip if zero or negative). Right uses `tref.varsCount < 0` (skip only if negative, allow zero through). This is a known parity trap (#3).
- Left fragment `meanPosition` check: `tref.meanPosition / tref.varsCount <= i + 1` — integer division. Guards against absorbing left fragments too far from the start.
- Right fragment has NO `meanPosition` check.
- `findconseq()` is from `VariationUtils` — external dependency.
- `ismatchref()` direction: 3' uses `dir=1` (forward), 5' uses `dir=-1` (reverse).

---

### `realignIndels()`

**Source**: L571–L586
**Purpose**: Dispatcher that calls the five realignment sub-procedures in order.
**Called By**: `process()` (when `performLocalRealignment`)

**Algorithm (Step-by-Step)**:
1. `realigndel(bams, positionToDeletionsCount)`
2. `realignins(positionToInsertionCount)`
3. `realignlgdel(svStructures.svfdel, svStructures.svrdel)`
4. `realignlgins30()`
5. `realignlgins(svStructures.svfdup, svStructures.svrdup)`

**Parity Warnings**:
- Order is load-bearing. Each step mutates shared maps that subsequent steps read. Reordering breaks parity.
- `realignlgdel` calls `realigndel` recursively. `realignlgins30` calls `realignins` and `realigndel` recursively.

---

### `realigndel()`

**Source**: L593–L858
**Purpose**: Absorbs nearby mismatches and soft-clipped reads into known short deletions. Two passes: forward (absorb mismatches/clips) and reverse (merge complex deletion variants).
**Called By**: `realignIndels()`, `realignlgdel()` (recursively with singleton map and null bams), `realignlgins30()` (recursively)
**Parameters**:
- `bamsParameter`: `String[]` — BAM file list. **Can be null** when called from `realignlgdel`/`realignlgins30`.
- `positionToDeletionsCount`: `Map<Integer, Map<String, Integer>>`

**Algorithm (Step-by-Step)**:

**Pass 1** (forward, lines L598–L831):
1. `bams = bamsParameter == null ? null : this.bams` — **shadow**: if parameter is null, local `bams` is null. If non-null, uses instance `this.bams`. The parameter value itself is never used (only its null-ness).
2. `tmp = fillAndSortTmp(positionToDeletionsCount)` — sorted by desc count, asc position, desc descriptionString.
3. For each `tpl` in `tmp`:
   a. `p = tpl.position`, `vn = tpl.descriptionString`, `dcnt = tpl.count`.
   b. `vref = getVariation(nonInsertionVariants, p, vn)` — get-or-create.
   c. Parse `dellen` from `vn`: extract digits after leading `-`, add digits before trailing `^N`.
   d. Parse `extra`, `extrains`, `inv5`, `inv3` from `vn` using regex patterns.
   e. Build flanking sequences:
      - `wupseq` = 5' flanking: `joinRef(ref, max(p-200, 1), p-1) + extra`. If `inv3` non-empty: `wupseq = inv3`.
      - `sanpseq` = 3' flanking: `extra + joinRef(ref, p + dellen + extra.length() - extrains.length(), min(p+200, chrLength))`. If `inv5` non-empty: `sanpseq = inv5`.
   f. **Find mismatches**:
      - `r3 = findMM3(ref, p, sanpseq)` — 3' direction.
      - `r5 = findMM5(ref, p + dellen + extra.length() - extrains.length() - 1, wupseq)` — 5' direction.
   g. Combine mismatches: `mmm = mm3 + mm5`.
   h. For each mismatch `(mm, mp, me)` in `mmm`:
      - If `mm.length() > 1`: insert `&` at index 1 (MNP encoding).
      - Look up `tv = nonInsertionVariants[mp][mm]`. Skip if null, zero count, low quality (`meanQuality/varsCount < goodq`), high position (`meanPosition/varsCount > nm+4`), or high count (`>= dcnt + dellen` or ratio `>=8`).
      - If `mp > p && me == 5`: adjust refCoverage by factor, call `adjRefCnt(tv, refVariant, dellen)`.
      - Compute `lref`: reference variant at `p` if applicable.
      - `adjCnt(vref, tv, lref)` — absorb mismatch into deletion.
      - Remove `mm` from `nonInsertionVariants[mp]`.
   i. **Single-mismatch cleanup**: If `misp3 != 0` and only 1 mismatch3 and its count < `dcnt`: remove it.
   j. **5' soft clip absorption**: Match reversed `wupseq` via `ismatch(seq, wupseq, -1)`.
   k. **3' soft clip absorption**: Match `substr(sanpseq, sc3pp - p)` via `ismatch(seq, ..., 1)`.
   l. **noPassingReads check**: If `bams != null && bams.length > 0 && pe - p >= 5 && pe - p < maxReadLength - 10 && ...`:
      - `adjCnt(vref, h, h)` — absorb reference into deletion.

**Pass 2** (reverse, lines L833–L858):
1. Loop `i` from `tmp.size() - 1` down to **`i > 0`** (element 0 is SKIPPED):
   a. Match `MINUS_NUMBER_AMP_ATGCs_END` against `vn` → extract `tn` (base deletion without `&extra`).
   b. Get `tref = nonInsertionVariants[p][tn]`. If `vref.varsCount < tref.varsCount`: `adjCnt(tref, vref)`, remove `vn`.

**Parity Warnings**:
- **bams shadowing** (trap #2): The `bamsParameter` value is never actually read — only its null-ness. When null, `noPassingReads` check is skipped entirely.
- **Pass 2 loop `i > 0`** (trap #12): Index 0 is intentionally skipped. Rust: `(1..tmp.len()).rev()`.
- **`getVariation` vs `getVariationMaybe`**: `getVariation` creates the entry if absent. Critical distinction.
- `fillAndSortTmp` sort order affects which variant gets processed first.
- `adjCnt` with 3 args: the third arg (`lref`) has its counts subtracted to prevent double-counting.

---

### `realignins()`

**Source**: L864–L1170
**Purpose**: Absorbs nearby mismatches and soft clips into known short insertions. Returns `NEWINS` string (modified insertion description).
**Called By**: `realignIndels()`, `realignlgins30()`, `realignlgins()` (recursively with singleton maps)
**Returns**: `String NEWINS` — modified insertion name, or `""` if none.

**Algorithm (Step-by-Step)**:

**Pass 1** (forward, lines L868–L1140):
1. `tmp = fillAndSortTmp(positionToInsertionCount)`.
2. For each `tpl`:
   a. Parse `insert` from `vn` via `BEGIN_PLUS_ATGC`. If no match: `continue`.
   b. Parse `ins3` from `DUP_NUM_ATGC` (dup format). `inslen = insert.length()` (+ dup parts).
   c. Parse `extra`, `compm`, `newins`, `newdel`.
   d. Build flanking sequences:
      - `wupseq` = `joinRef(ref, max(position-150,1), position) + tn`.
      - `sanpseq` depends on `ins3` presence.
   e. Find mismatches via `findMM3`, `findMM5`.
   f. `vref = getVariation(insertionVariants, position, vn)`.
   g. Mismatch absorption loop — same pattern as `realigndel`, targeting `insertionVariants`.
   h. Single-mismatch cleanup.
   i. **5' soft clip absorption**: `ismatch(seq, wupseq, -1)`.
   j. **3' soft clip absorption**: `mseq` differs for `ins3` vs non-`ins3`.
      - **NEWINS mutation**: If `insert.length() + 1 == vn.length()` and `insert.length() > maxReadLength` and `sc3pp >= position + 1 + insert.length()`: walk and correct mismatches, rename variant key, set `NEWINS = tvn`.
   k. **Ref factor adjustment**: If `first3 > first5 + 3` and distance < `maxReadLength * 0.75`:
      - `adjRefFactor(refVariant, (first3 - first5 - 1) / (double)maxReadLength)`.
      - `adjRefFactor(vref, -(first3 - first5 - 1) / (double)maxReadLength)`.

**Pass 2** (reverse, lines L1146–L1170):
1. Same pattern as `realigndel` pass 2 but for insertions.
2. Uses `adjCnt(tref, vref, getVariationMaybe(nonInsertionVariants, p, ref.get(p)))` — note 3rd arg.

**Parity Warnings**:
- **Pass 2 loop `i > 0`** — same as `realigndel`: index 0 is skipped.
- The `NEWINS` mutation logic is extremely rare (requires insertion longer than `maxReadLength`). When it triggers, it renames the insertion key in-place. Rust must handle this HashMap key rename.
- The `mseq` construction for `ins3` vs non-`ins3` cases differs significantly.

---

### `realignlgdel()`

**Source**: L1177–L1594
**Purpose**: Discovers large deletions by finding breakpoints from unpaired 5' and 3' soft clips.
**Called By**: `realignIndels()`

**Algorithm (Step-by-Step)**:

**5' pass** (lines L1183–L1386):
1. Collect `softClips5End` within `region ± EXTENSION` → sort by `COMP2`.
2. For each `(p, sc5v, cnt)`:
   a. Threshold/used checks. `seq = findconseq(sc5v, 5)`. Min 7 chars.
   b. **Direct breakpoint**: `bp = findbp(seq, p - 5, ref, -1, chr)`.
   c. If `bp == 0`: `findMatch(seq, reference, p, -1, SEED_1, 1)`. Require `p - bp > 15 && p - bp < SVMAXLEN`. `bp++`. `markSV`. If `bp < region.start`: extend reference, `partialPipeline`.
   d. `dellen = p - bp`.
   e. Build genotype `gt` with `extra` bases from mismatched boundary.
   f. **Find 3' breakpoint match**: Walk `n` while `ref[bp+n] == ref[bp+dellen+n]` → compute `sc3p`.
   g. Create variant at `bp` with `gt`. `adjCnt(tv, sc5v)`. Mark `sc5v.used = bp != 0`.
   h. Set `refCoverage[bp]` from `refCoverage[p]` if absent.
   i. **Absorb matching 3' clip**: Walk `ip` from `bp+1` to `sc3p`: `rmCnt` mismatched bases.
   j. **Recursive `realigndel`**: `realigndel(bams, singletonMap(bp, singletonMap(gt, tv.varsCount)))`.
   k. SV splits update. `addVarFactor` if `svcov > tv.varsCount`.

**3' pass** (lines L1390–L1594):
Mirror of 5' pass with key differences:
1. `bp = findbp(seq, p + 5, ref, 1, chr)`. If `bp == 0`: `findMatch(...)`.
2. If `bp > region.end`: extend reference, `partialPipeline`.
3. **5' left-alignment**: Walk `bp--` while `ref[bp-1] == ref[bp+dellen-1]`.
4. **`meanPosition` adjustment**: `sc3v.meanPosition += dellen * sc3v.varsCount` before `adjCnt`.
5. `sc3v.used = p + dellen != 0` (different from 5' pass).
6. Recursive `realigndel(bams, dels5)` uses `new HashMap` (not `singletonMap`).

**Parity Warnings**:
- **`sc5v.used = bp != 0`** (trap #8): Boolean from int comparison.
- **`sc3v.used = p + dellen != 0`**: Different formula than 5' pass.
- **`partialPipeline`**: Re-entrant call running entire CigarParser, mutating all shared maps.
- **`singletonMap` vs `new HashMap`**: 5' pass uses immutable `singletonMap`, 3' uses mutable `HashMap`.
- **`joinRef` overloads**: Uses INCLUSIVE endpoints.

---

### `realignlgins30()`

**Source**: L1599–L1799
**Purpose**: Pairs 5'/3' soft clips to discover large insertions (typically ≥30bp) where both clip ends overlap.
**Called By**: `realignIndels()`

**Algorithm (Step-by-Step)**:
1. Collect `softClips5End` → `tmp5`, sort by `COMP3`. Collect `softClips3End` → `tmp3`, sort by `COMP3`.
2. **Nested loop**: For each `(p5, sc5v, cnt5)` in `tmp5`:
   a. For each `(p3, sc3v, cnt3)` in `tmp3`:
      - Distance filters: `p5 - p3 > maxReadLength * 2.5`, `p3 - p5 > maxReadLength - 10`.
      - `seq5 = findconseq(sc5v, 5)`, `seq3 = findconseq(sc3v, 3)`. Both > 10 length.
      - Bias filter: `cnt5 / (double)cnt3 >= 0.08 && <= 12`.
      - **`find35match(seq5, seq3)`** → `(bp5, bp3, score)`. If `score == 0`: `continue`.
      - Build insertion `ins` from seq3/seq5 overlap.
      - If `islowcomplexseq(ins)`: skip.
      - **Branch on `p5 > p3`** (overlap): determine if deletion, insertion, or MNP.
      - **Branch on `p5 <= p3`** (standard): determine if tandem dup or non-tandem.
      - Mark both clips `used`. Set `pstd`, `qstd`.
      - **Post-creation**: `adjCnt` calls, optional `noPassingReads`, recursive `realignins`/`realigndel`.
      - `break` inner loop (one match per 5' clip).

**Parity Warnings**:
- **`COMP3` vs `COMP2`** (trap #11): `COMP3` sorts by `.count` field, not `.softClip.varsCount`.
- **`find35match` early return** (trap #5): First qualifying match, not the best. Loop order matters.
- **Integer division for `smscore`**: `score / 2` is Java int division.
- **Tandem dup `rpt` calculation**: Uses `double` division with implicit `int` cast for `joinRef` endpoint. Truncation must be preserved.
- **`break` on inner loop**: Each 5' clip matches at most one 3' clip.

---

### `realignlgins()`

**Source**: L1805–L2150
**Purpose**: Discovers large insertions/DUPs from single soft clips via `findbi()` or `findMatch()`.
**Called By**: `realignIndels()`

**Algorithm (Step-by-Step)**:

**5' pass** (lines L1811–L1982):
1. Collect `softClips5End` → sort by `COMP2`.
2. For each `(p, sc5v, cnt)`:
   a. `findbi(seq, p, ref, -1, chr)` → `(bi, ins, EXTRA)`.
   b. If `bi == 0`: `findMatch(...)`. Build `<dup>` format or `joinRefFor5Lgins`. `markDUPSV`. Set SV metadata.
   c. Create variant at `bi`: `iref = getVariation(insertionVariants, bi, "+" + ins)`.
   d. **Repeat flag**: Check if ins matches ref at `bi+1`.
   e. **Extend soft clip bases**: From `ii = len + 1` (note: **`len + 1`**) to `sc5v.seq.lastKey()`.
   f. `sc5v.used = bi + len != 0`.
   g. Recursive `realignins(singletonMap(...))`.
   h. **noPassingReads**: If `rpflag && ...`: `adjCnt(kref, mref, mref)`.

**3' pass** (lines L1984–L2150):
Mirror with key differences:
1. `findbi(seq, p, ref, 1, chr)`.
2. If `bi == 0`: `findMatch(...)`. **5' left-shift**: Walk `p--; bi--` while `ref[p-1] == ref[bi-1]`.
3. `lref = getVariationMaybe(...)`. If `p - bi > sc3v.meanPosition / cnt`: `lref = null`.
4. **Extend from `ii = len`** (note: **`len`**, not `len + 1`).
5. `sc3v.used = true` (unconditional, unlike 5' pass).

**Parity Warnings**:
- **Loop offset asymmetry** (trap #15): 5' starts at `ii = len + 1`, 3' starts at `ii = len`.
- **`sc5v.used = bi + len != 0`** vs **`sc3v.used = true`**: Different used-marking logic.
- **`joinRefFor5Lgins` / `joinRefFor3Lgins`**: Different logic between 5' and 3' in `VariationUtils`.
- **`kref` vs `iref`**: After recursive `realignins`, uses `newins` which may differ from original `"+" + ins`.

---

### `fillAndSortTmp()`

**Source**: L2155–L2196
**Purpose**: Flattens a nested `Map<Integer, Map<String, Integer>>` into a sorted `List<SortPositionDescription>`.
**Called By**: `realigndel()`, `realignins()`, `adjustMNP()`

**Algorithm (Step-by-Step)**:
1. Create empty list.
2. For each `(position, v)` in `changes`:
   a. For each `(descriptionString, cnt)` in `v`:
      - Add `SortPositionDescription(position, descriptionString, cnt)`.
3. Sort: **Primary** desc `count`, **Secondary** asc `position`, **Tertiary** desc `descriptionString` (via `s2.compareTo(s1)`).
4. Return sorted list.

**Parity Warnings**:
- **Sort order** (trap #1): Tertiary uses `s2.compareTo(s1)` — reverse lexicographic. Java's `String.compareTo()` matches Rust's `String::cmp()` for ASCII. The reverse is critical.
- The outer map is iterated in arbitrary order (Java HashMap). The sort produces deterministic order.

---

### `find35match()`

**Source**: L2214–L2260
**Purpose**: Finds overlapping match between 5' and 3' sequences by sliding window comparison.
**Called By**: `realignlgins30()`

**Algorithm (Step-by-Step)**:
1. `longMismatch = 2`, `maxMatchedLength = 0`.
2. For `i = 0` to `seq5.length() - 9`:
   a. For `j = 1` to `seq3.length() - 9`:
      - Walk `n` comparing `seq3.charAt(-j - n)` (from end) with `seq5.charAt(i + n)` (from start).
      - If `numberOfMismatch > longMismatch`: break.
      - Qualifying: `n - nm > maxMatchedLength && n - nm > 8 && nm/(double)n < 0.1 && (full length reached)`.
      - If qualifies: **return immediately** (early return, trap #5).
3. Return `Match35(b5, b3, maxMatchedLength)`.

**Parity Warnings**:
- **Early return** (trap #5): Returns FIRST qualifying match. Loop order `(i outer, j inner)` is critical.
- `charAt(seq3, -j - n)` uses negative indexing from end. Rust must implement same semantics.
- Condition `n + j <= seq3.length()` uses `<=` not `<` — intentional boundary.

---

### `noPassingReads()`

**Source**: L2270–L2306
**Purpose**: Checks whether reads span across a deletion gap. Returns `true` if NO spanning reads found.
**Called By**: `realigndel()`, `realignlgins30()`, `realignlgins()`

**Algorithm (Step-by-Step)**:
1. Open `SamView(bam, "0", region, validationStringency)` for each BAM.
2. For each record: skip if CIGAR contains `dlen + "D"`. Count `cnt` (fully spanning) and `midcnt` (ends in middle).
3. Return `cnt <= 0 && midcnt + 1 > 0` (simplifies to `cnt == 0`).

**Parity Warnings**:
- Opens BAM file for each call — expensive. Rust should match this behavior.
- `getAlignedLength` is from `CigarParser` — counts M + D + N bases (not I or S).

---

### `ismatch()` (2 overloads)

**Source**: L2313–L2345
**Purpose**: Tests if two sequences match with ≤MM mismatches and <15% mismatch rate.
**Called By**: `realigndel()`, `realignins()`, `realignlgins30()`

**Algorithm**: Compare `seq1[n]` with `seq2[dir*n - (dir==-1 ? 1 : 0)]`. For `dir=-1`: accesses `seq2[-n-1]` (from end). Return `mm <= MM && mm/(double)seq1.length() < 0.15`.

**Parity Warning**: The 15% threshold uses `seq1.length()` as denominator, not comparison length.

---

### `islowcomplexseq()`

**Source**: L2353–L2385
**Purpose**: Returns true if >75% single base or <3 distinct bases.

**Parity Warning**: `count / (double)len > 0.75` — strict `>` (0.75 exactly is NOT low complexity). `ntcnt < 3` checked after individual checks.

---

### `adjInsPos()`

**Source**: L2408–L2423
**Purpose**: Left-aligns insertion position by rotating the insertion sequence.

**Algorithm**: Walk `bi--` while `ref.get(bi) == ins.charAt(ins.length() - n)` with wraparound. Rotate insertion string.

**Parity Warning**: `ref.get(bi)` may return null below position 1. Java compares `null` with `Character` — Rust should short-circuit.

---

### `findbi()`

**Source**: L2431–L2538
**Purpose**: Finds insertion breakpoint by scanning for indel boundary in soft-clipped sequence.
**Called By**: `realignlgins()`

**Algorithm**: For `n = 6` to `seq.length() - 1`: try aligning `seq[n:]` against reference. Qualifying: ≥3 distinct bases, ≥8 aligned, <15% mismatch. Build insertion from `seq[0:n]` + extra. Direction-specific string construction.

**Parity Warnings**:
- **Direction asymmetry** (trap #14): 5' uses `reverse + append`, 3' uses left-align + rotate.
- `dirExt` is 1 for dir=-1, 0 for dir=1.
- `ept` loop checks TWO consecutive mismatches to continue.

---

### `findbp()`

**Source**: L2545–L2600
**Purpose**: Finds deletion breakpoint by scanning for where soft-clip matches reference at offset.
**Called By**: `realignlgdel()`

**Algorithm**: For `n = 0` to `indelsize - 1`: try aligning `seq` against `ref[startPosition + dir*n + dir*i]`. Qualifying: ≥3 distinct bases, ≥8+n/10 aligned, <12% mismatch.

**Parity Warnings**:
- **`maxmm - n/100`** (trap #9): Integer division. Threshold only changes at n=100, 200, etc.
- Position formula: `direction < 0 ? direction : 0` — adds 1 for dir=-1.

---

### `adjRefCnt()`

**Source**: L2608–L2647
**Purpose**: Adjusts reference variant counts by factor derived from mismatch position relative to deletion.

**Algorithm**: Factor = `(meanPosition/varsCount - len + 1) / (meanPosition/varsCount)`, capped [0, 1]. Subtract `(int)(f * field)` from reference.

**Parity Warning**: `(int)(f * tv.varsCount)` — Java truncates toward zero. Rust `as i32` matches.

---

### `adjRefFactor()`

**Source**: L2653–L2695
**Purpose**: Adjusts reference by multiplicative factor (negative = increase, positive = decrease).

**Algorithm**: Subtract `(int)(factor_f * field)` for counts. Compute `factorCnt = abs(delta) / oldVarsCount`. Negate `factorCnt` when `factor_f < 0`. Apply `factorCnt` to position/quality stats. Apply `factor_f` to nm/direction.

**Parity Warnings**:
- **Sign preservation** (trap #10): `&&` before `||` precedence. `factorCnt` negated when `factor_f < 0`.
- **Inconsistency**: `meanPosition`/`meanQuality`/`meanMappingQuality` use `factorCnt`. `numberOfMismatches`/direction use `factor_f`. Intentional.

---

### `addVarFactor()`

**Source**: L2701–L2722
**Purpose**: Adds counts to a variant by multiplicative factor.

**Parity Warning**: No `correctCnt` call — unlike `adjRefFactor`. Values can go negative.

---

### `findMM5()`

**Source**: L2730–L2780
**Purpose**: Walks 5' direction from variant position, finding mismatches and soft-clip break positions.
**Called By**: `realigndel()`, `realignins()`

**Algorithm**: Walk from end of `wupseq` backward, collecting mismatches. After single mismatch: walk matches, then secondary walk. Mark soft clips `used` as side effect.

**Parity Warnings**:
- **Side effect** (trap #4): Marks soft clips `used` during walk. Load-bearing.
- Uses `mcnt < longmm` (up to 3 iterations). Compare with `findMM3` which uses `mcnt <= longmm`.
- Initial mismatch walk creates cumulative mismatch strings.

---

### `findMM3()`

**Source**: L2788–L2840
**Purpose**: Walks 3' direction from variant position, finding mismatches.
**Called By**: `realigndel()`, `realignins()`

**Algorithm**: Walk matches FIRST, then mismatches (opposite of `findMM5`). Same secondary walk pattern.

**Parity Warnings**:
- **Mismatch loop bound asymmetry** (trap #16): Uses `mcnt <= longmm` (up to 4 iterations). `findMM5` uses `mcnt < longmm` (up to 3).
- Walks matches first, then mismatches. `findMM5` walks mismatches first.

---

### `ismatchref()` (2 overloads)

**Source**: L2890–L2929
**Purpose**: Checks if consensus sequence matches reference at given position/direction.
**Called By**: `adjustMNP()`

**Algorithm**: Compare `seq[dir==1 ? n : -n-1]` against `ref[position + dir*n]`. Return `mm <= MM && mm/(double)len < 0.15`.

**Parity Warning**: When `dir=-1`, the already-reversed sequence is read backward — double-reversal cancels out. `ref.get(...)` returning null → immediate `false`.

---

### `rmCnt()`

**Source**: L2936–L2949
**Purpose**: Subtracts counts of one Variation from another.
**Called By**: `realignlgdel()`

**Parity Warning**: Does NOT subtract `numberOfMismatches` — possible oversight or intentional. `correctCnt` clamps negatives.

---

### COMP2 (Comparator)

**Source**: L340–L348
**Used By**: `realignlgdel()`, `realignlgins()`

Sort by descending `softClip.varsCount`, then ascending `position`.

---

### COMP3 (Comparator)

**Source**: L350–L358
**Used By**: `realignlgins30()`

Sort by descending `.count` field, then ascending `position`.

**Parity Warning** (trap #11): Uses `.count` (explicit field), not `.softClip.varsCount`. Can diverge if `varsCount` is modified after list construction.

## Cross-Module Dependencies

### Outbound Calls

| Target | Method | Called From | Notes |
|--------|--------|-------------|-------|
| `VariationUtils` | `findconseq()` | adjustMNP, realigndel, realignins, realignlgdel, realignlgins30, realignlgins | Consensus sequence from soft clip |
| `VariationUtils` | `adjCnt()` | Many methods | Add counts from one Variation to another |
| `VariationUtils` | `getVariation()` | realigndel, realignins, realignlgdel, realignlgins30, realignlgins | Get-or-create variant at position |
| `VariationUtils` | `getVariationMaybe()` | realigndel, realignins, realignlgdel, realignlgins | Get variant or null |
| `VariationUtils` | `incCnt()` | Many methods | Increment coverage count |
| `VariationUtils` | `correctCnt()` | adjRefCnt, adjRefFactor, rmCnt | Clamp negatives |
| `VariationUtils` | `joinRef()` | realigndel, realignins, realignlgdel, realignlgins30 | Join reference chars into string (INCLUSIVE endpoints) |
| `VariationUtils` | `joinRefFor5Lgins()` | realignlgins (5') | Join ref with seq overlay |
| `VariationUtils` | `joinRefFor3Lgins()` | realignlgins (3') | Join ref with seq overlay |
| `StructuralVariantsProcessor` | `findMatch()` | realignlgdel, realignlgins | Find distant match for soft clip |
| `StructuralVariantsProcessor` | `markSV()` | realignlgdel | Mark DEL SV metadata |
| `StructuralVariantsProcessor` | `markDUPSV()` | realignlgins | Mark DUP SV metadata |
| `CigarParser` | `getAlignedLength()` | noPassingReads | Get aligned length from CIGAR |
| `ReferenceResource` | `getReference()` | realignlgdel, realignlgins | Load reference for extended region |
| `ReferenceResource` | `isLoaded()` | realignlgdel, realignlgins | Check if reference region is loaded |
| `RecordPreprocessor` | `getChrName()` | initFromScope | Get chromosome name |
| `GlobalReadOnlyScope` | `instance().conf.*` | Throughout | Configuration: minr, goodq, y, disableSV, performLocalRealignment, indelsize, chrLengths, SVMINLEN |
| `GlobalReadOnlyScope` | `getMode()` | realignlgdel, realignlgins | Get pipeline mode for partialPipeline |
| `VariationMap` | `getSV()` | realignlgdel, realignlgins | Get/create SV metadata |
| `VariationMap` | `removeSV()` | realignlgdel | Remove SV metadata |
| `Utils` | `substr()`, `charAt()`, `reverse()`, `join()`, `isEquals()`, `isNotEquals()`, `isHasAndEquals()`, `isHasAndNotEquals()` | Throughout | String/char utilities |

### Inbound Callers

| Caller | Calls | Notes |
|--------|-------|-------|
| `SAMFileParser.parseSAM()` | `process()` via Module interface | Main pipeline entry |
| Self (recursive) | `realigndel()` | From `realignlgdel()`, `realignlgins30()` |
| Self (recursive) | `realignins()` | From `realignlgins30()`, `realignlgins()` |

## Data Structures Read/Written

| Data Structure | Type | R/W | Notes |
|----------------|------|-----|-------|
| `nonInsertionVariants` | `Map<Integer, VariationMap<String, Variation>>` | **R/W** | Variants added, removed, modified throughout |
| `insertionVariants` | `Map<Integer, VariationMap<String, Variation>>` | **R/W** | Modified in realignins, realignlgins30, realignlgins |
| `refCoverage` | `Map<Integer, Integer>` | **R/W** | Incremented via incCnt, directly put |
| `softClips5End` | `Map<Integer, Sclip>` | **R/W** | `.used` flag set, counts read via findconseq |
| `softClips3End` | `Map<Integer, Sclip>` | **R/W** | `.used` flag set, counts read via findconseq |
| `svStructures` | `SVStructures` | **R** | SV lists read for filterSV, passed to realignlgdel/realignlgins |
| `SOFTP2SV` | `Map<Integer, List<Sclip>>` | **W** | Populated by filterSV, read for SV lookups |
| `reference.referenceSequences` | `Map<Integer, Character>` | **R** | Reference base lookups throughout |
| `mnp` | `Map<Integer, Map<String, Integer>>` | **R** | Read by adjustMNP |
| `positionToDeletionsCount` | `Map<Integer, Map<String, Integer>>` | **R** | Read by realigndel |
| `positionToInsertionCount` | `Map<Integer, Map<String, Integer>>` | **R** | Read by realignins |
| `Variation` fields | `varsCount, meanPosition, meanQuality, ...` | **R/W** | Mutated via adjCnt, adjRefCnt, adjRefFactor, addVarFactor, rmCnt |
| `Sclip` fields | `used, varsCount, mates, soft, softp, mstart, mend, ...` | **R/W** | Mutated by filterSV, realign methods |

## Known Parity Traps

1. **`fillAndSortTmp` sort order**: Descending count, ascending position, descending descriptionString (Java `String.compareTo` reversed). The tertiary criterion changes which variant is processed first, affecting absorption order. Rust sort must match all three tiebreakers exactly.

2. **`bams` parameter shadowing in `realigndel`**: When `bamsParameter` is null (called from `realignlgdel`/`realignlgins30`), local `bams` is set to null, disabling the `noPassingReads` check. The parameter's actual value is never used — only its null-ness matters.

3. **Asymmetric `varsCount` checks in `adjustMNP`**: Left fragment uses `tref.varsCount <= 0` (skip if zero OR negative). Right fragment uses `tref.varsCount < 0` (skip only if negative, allowing zero through). Off-by-one is intentional.

4. **`findMM3`/`findMM5` mark soft clips `used` as side effect**: During the secondary walk (after finding a single mismatch + matches), these functions set `softClipsXEnd.get(pos).used = true`. Load-bearing — prevents the same soft clip from being consumed again.

5. **`find35match` early return**: Returns the FIRST qualifying match `(b5, b3, score)`, not the best. The nested `(i, j)` loop order means matches near `(i=0, j=1)` are found before `(i=5, j=10)`. Iteration order is parity-critical.

6. **Recursive `realigndel`/`realignins`**: Called with singleton maps and shared mutable variant maps. The outer loop iterates a pre-computed `tmp` list, not the map directly. Rust must ensure the same separation.

7. **`joinRef` overloads**: The module uses `joinRef(ref, from, to)` which is INCLUSIVE on both endpoints. Another overload has exclusive end — using the wrong one silently shifts the window by 1 base.

8. **`sc5v.used = bp != 0`**: Boolean from comparison result. In 5' pass of `realignlgdel`, if `bp` is 0, `used` becomes `false`. Contrast with 3' pass which uses `sc3v.used = p + dellen != 0`.

9. **Integer division thresholds**: `maxmm - n/100` (in `findbp`), `(int)(f * cnt)` (in `adjRefCnt`/`adjRefFactor`), `score/2` (in `realignlgins30`). Java truncation-toward-zero semantics. Rust `as i32` matches for positive values.

10. **`adjRefFactor` sign preservation**: The condition uses Java operator precedence (`&&` before `||`). Since `factorCnt` from `Math.abs()` is always ≥0, the second branch never triggers. Effect: `factorCnt` is negated when `factor_f < 0`.

11. **COMP2 vs COMP3**: COMP2 sorts by `softClip.varsCount`. COMP3 sorts by the `.count` field. These can diverge if `varsCount` is mutated after list construction.

12. **Pass-2 reverse loop `i > 0` (not `>= 0`)**: In both `realigndel` and `realignins`, element 0 is intentionally skipped. Rust must use `(1..tmp.len()).rev()`.

13. **`partialPipeline` re-entrant calls**: `realignlgdel` and `realignlgins` can trigger `partialPipeline`, running the entire CigarParser on an extended region, mutating all shared variant maps. Order-dependent mutation.

14. **`findbi` direction asymmetry**: 5' direction uses `reverse + append`. 3' direction uses left-align + rotate. Different string construction logic.

15. **`realignlgins` loop offset**: 5' pass starts extending soft clip bases at `ii = len + 1`. 3' pass starts at `ii = len`. Off-by-one between directions.

16. **`findMM3` vs `findMM5` mismatch loop bound**: `findMM3` uses `mcnt <= longmm` (iterates up to 4 times). `findMM5` uses `mcnt < longmm` (iterates up to 3 times). Asymmetric by design.

17. **`ismatch` direction formula**: `dir * n - (dir == -1 ? 1 : 0)` — for `dir=-1`, accesses `seq2[-n-1]` (from end). The `substr()` function handles negative indices.

18. **`realignlgdel` 3' `meanPosition` adjustment**: Before calling `adjCnt`, the 3' pass adds `dellen * sc3v.varsCount` to `sc3v.meanPosition`. The 5' pass does NOT do this.
