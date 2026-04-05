# StructuralVariantsProcessor

**Source**: `modules/StructuralVariantsProcessor.java`
**LOC**: 2,307
**Rust counterpart**: `src/mods/structural_variants_processor.rs`
**Risk**: HIGH
**Pipeline Stage**: After VariationRealigner, before ToVarsBuilder
**Status**: complete

## Overview

StructuralVariantsProcessor detects structural variants (DEL, INV, DUP) from two evidence sources: **soft-clipped read consensus sequences** (splits) and **discordant read pairs** (pairs). It sits between VariationRealigner and ToVarsBuilder in the VarDict pipeline.

The module processes SV evidence stored in `SVStructures` (populated by CigarParser during BAM traversal). For each SV cluster, it:
1. Finds the consensus sequence from associated soft-clip positions using `findconseq()`
2. Attempts to align that consensus to the reference using seed-based matching (`findMatch` / `findMatchRev`)
3. Constructs variant description strings (e.g., `-500` for deletions, `-300^<inv250>` for inversions, `+<dup200>` for duplications)
4. Adjusts variant counts, reference coverage, and marks used clusters

After SV processing, `adjSNV()` rescues short soft-clipped reads as SNVs. Debug output of remaining clips is available via `outputClipping()`.

## Method Inventory

| Method | Lines | Analyzed? | Summary |
|--------|-------|-----------|---------|
| `initFromScope()` | L63–L83 | yes | Extract all fields from Scope into instance fields |
| `process()` | L91–L130 | yes | Entry point: findAllSVs → adjSNV → optional debug output |
| `StructuralVariantsJsonlWriter` (inner class) | L131–L319 | yes | Debug JSONL serializer (non-parity-critical) |
| `findAllSVs()` | L324–L349 | yes | Orchestrates all SV finding: DEL → INV → findsv → DELdisc → INVdisc → DUPdisc |
| `findDEL()` | L354–L681 | yes | Find DEL SVs from soft-clip + discordant pair evidence |
| `findINV()` | L686–L690 | yes | Dispatches to findINVsub for all 4 INV lists |
| `findINVsub()` | L700–L892 | yes | Find INV SVs from a single orientation list |
| `findsv()` | L897–L1175 | yes | Find DEL/INV from raw soft-clip scanning (5' and 3') |
| `fillAndSortTmpSV()` | L1177–L1195 | yes | Filter and sort soft-clip entries by count (descending) |
| `findDELdisc()` | L1200–L1430 | yes | Find DEL with discordant pairs only (no soft-clip evidence) |
| `findINVdisc()` | L1435–L1565 | yes | Find INV with discordant pairs only |
| `findDUPdisc()` | L1570–L1825 | yes | Find DUP with discordant pairs only |
| `markSV()` | L1830–L1870 | yes | Mark overlapping SV clusters as used (static) |
| `markDUPSV()` | L1875–L1915 | yes | Mark overlapping DUP SV clusters as used (static) |
| `isOverlap()` | L1920–L1945 | yes | Determine if two SV intervals overlap |
| `checkPairs()` | L1950–L1995 | yes | Check for discordant pair support for a candidate SV |
| `findMatchRev()` (2-arg) | L1997–L2002 | yes | Reverse-strand matching with default MM=3 |
| `findMatchRev()` (6-arg) | L2008–L2075 | yes | Core reverse-strand seed-based alignment |
| `findMatch()` (4-arg) | L2078–L2082 | yes | Forward-strand matching with default MM=3 |
| `findMatch()` (6-arg) | L2088–L2175 | yes | Core forward-strand seed-based alignment (static) |
| `adjSNV()` | L2180–L2240 | yes | Rescue short soft-clipped reads as SNVs |
| `outputClipping()` | L2245–L2285 | yes | Debug output of unused soft-clipped reads |
| `PairsData` (inner class) | L2287–L2310 | yes | Tuple for pair evidence (pairs, pmean, qmean, Qmean, nm) |

## Method Analyses

---

### Analysis: StructuralVariantsProcessor.initFromScope()

**Source**: StructuralVariantsProcessor.java:L63–L83
**Purpose**: Extract all fields from `Scope<RealignedVariationData>` into instance fields for convenience
**Pipeline Stage**: Entry setup
**Called By**: `process()`

#### Parameters
- `scope`: `Scope<RealignedVariationData>` — the pipeline scope containing all data from previous stages

#### Algorithm (Step-by-Step)
1. Copy `scope.regionRef` → `this.reference`
2. Copy `scope.data.CURSEG` → `this.CURSEG`
3. Copy `scope.data.SOFTP2SV` → `this.SOFTP2SV`
4. Copy `scope.referenceResource` → `this.referenceResource`
5. Copy `scope.region` → `this.region`
6. Copy `scope.data.nonInsertionVariants` → `this.nonInsertionVariants`
7. Copy `scope.data.insertionVariants` → `this.insertionVariants`
8. Copy `scope.data.refCoverage` → `this.refCoverage`
9. Copy `scope.data.softClips5End` → `this.softClips5End`
10. Copy `scope.data.softClips3End` → `this.softClips3End`
11. Copy `scope.maxReadLength` → `this.maxReadLength`
12. `this.bams = scope.bam != null ? scope.bam.split(":") : null` — **null-check on bam string**
13. Copy `scope.bam` → `this.bam`
14. Copy `scope.splice` → `this.splice`
15. Copy `scope.data.svStructures` → `this.svStructures`
16. Copy `scope.data.duprate` → `this.duprate`
17. Copy `scope.data.previousScope` → `this.previousScope`
18. Copy `scope.out` → `this.variantPrinter`

#### Parity Warnings
- `bams` can be `null` if `scope.bam` is null; Rust must handle `Option<Vec<String>>`
- All fields are shallow copies (references to same objects) — mutations propagate back to scope data

---

### Analysis: StructuralVariantsProcessor.process()

**Source**: StructuralVariantsProcessor.java:L91–L130
**Purpose**: Main entry point — runs SV detection, SNV rescue, and optional debug output
**Pipeline Stage**: Main orchestrator
**Called By**: Pipeline (Module interface)

#### Parameters
- `scope`: `Scope<RealignedVariationData>` — input scope

#### Algorithm (Step-by-Step)
1. Call `initFromScope(scope)`
2. **If** `!instance().conf.disableSV`:
   - Call `findAllSVs()`
3. Call `adjSNV()` — always runs regardless of disableSV
4. **If** `instance().conf.y` (debug mode):
   - Call `outputClipping()`
   - Print timestamp to stderr
5. Construct new `RealignedVariationData` with all mutated maps
6. **If** env var `VARDICT_STRUCTURAL_VARIANTS_JSONL` is set and non-empty:
   - Write JSONL debug output via `StructuralVariantsJsonlWriter`
7. Return new `Scope<>(scope, realigned)`

#### Parity Warnings
- `adjSNV()` runs even when SV is disabled — must not skip it
- The JSONL writer is debug-only and not part of parity output
- The return creates a new Scope wrapping the same modified data structures

---

### Analysis: StructuralVariantsProcessor.findAllSVs()

**Source**: StructuralVariantsProcessor.java:L324–L349
**Purpose**: Orchestrate all SV finding routines in fixed order
**Pipeline Stage**: SV detection
**Called By**: `process()`

#### Algorithm (Step-by-Step)
1. Call `findDEL()` — deletions from soft-clip + pair evidence
2. Call `findINV()` — inversions from soft-clip + pair evidence
3. Call `findsv()` — general SV from raw soft-clip scanning
4. Call `findDELdisc()` — deletions from discordant pairs only
5. Call `findINVdisc()` — inversions from discordant pairs only
6. Call `findDUPdisc()` — duplications from discordant pairs only

#### Parity Warnings
- **Order matters**: each routine marks clusters as `used`, which prevents later routines from processing them. The order DEL → INV → findsv → DELdisc → INVdisc → DUPdisc must be preserved exactly.

---

### Analysis: StructuralVariantsProcessor.findDEL()

**Source**: StructuralVariantsProcessor.java:L354–L681
**Purpose**: Find DEL SVs using soft-clip consensus + discordant pair evidence
**Pipeline Stage**: SV detection (DEL)
**Called By**: `findAllSVs()`

#### Algorithm (Step-by-Step)

**Part 1: 5' forward deletions (`svfdel` loop, L358–L517)**

For each `del` in `svStructures.svfdel`:
1. Skip if `del.used` or `del.varsCount < instance().conf.minr`
2. Sort `del.soft` entries by value descending → get `softp` (position with most soft-clip support)
3. **If `softp != 0`** (has associated soft-clip position):
   a. Check `softClips3End.containsKey(softp)` → skip if missing
   b. Get `scv = softClips3End.get(softp)` → skip if `scv.used`
   c. `seq = findconseq(scv, 0)` → skip if empty or `< SEED_2` (12)
   d. Ensure reference is loaded for `del.mstart..del.mend`, run partial pipeline if needed
   e. `match = findMatch(seq, reference, softp, 1)` → get `bp` and `extra`
   f. Skip if `bp == 0`
   g. Skip if NOT `(bp - softp > 30 && isOverlap(softp, bp, del.end, del.mstart, maxReadLength))`
   h. `bp--` (convert to 0-based end position)
   i. `dellen = bp - softp + 1`
   j. **Left-align**: while `ref[bp] == ref[softp - 1]`, decrement both `bp` and `softp`
   k. Create variation at `nonInsertionVariants[softp]["-" + dellen]`, set `varsCount = 0`
   l. Create/get SV at position `softp`, set `type = "DEL"`, accumulate `pairs`, `splits`, `clusters`
   m. Update `refCoverage` at `softp` and `bp`
   n. Call `adjCnt(variation, scv, refVariant)` — adds soft-clip counts, subtracts from reference variant
   o. Create temporary `Variation tv` from pair evidence (`del.varsCount`), split counts evenly between fwd/rev
   p. Call `adjCnt(variation, tv)` — adds pair counts
   q. Mark `del.used = true`
   r. Call `markSV(softp, bp, [svrdel], maxReadLength)`

4. **Else** (`softp == 0`, no associated soft-clip):
   a. Ensure reference is loaded
   b. **Iterate all `softClips3End`** entries:
      - Skip if `scv.used`
      - Skip if NOT `(i >= del.end - 3 && i - del.end < 3 * maxReadLength)`
      - Get consensus seq, skip if empty / short
      - `findMatch(seq, reference, softp, 1)` — if bp=0, retry with `SEED_2, 0`
      - Skip if bp=0 or overlap/distance check fails
      - Same DEL construction as step 3(h–r)
      - **Break after first match** — only takes the first soft-clip that works

**Part 2: 3' reverse deletions (`svrdel` loop, L518–L681)**

Mirror logic of Part 1 but reversed:
- Uses `softClips5End` instead of `softClips3End`
- `findMatch` with `dir = -1`
- After match, `bp++; softp--`
- Variation placed at `bp` with `-dellen`
- When no `softp`, iterates `softClips5End` with condition `i <= del.start + 3 && del.start - i < 3 * maxReadLength`
- On 3' no-softp path: calls `incCnt(refCoverage, bp, scv.varsCount)` — **different from 5' no-softp path**

#### Mutable State Modified
- `nonInsertionVariants` — new Variation entries for detected DELs
- `refCoverage` — updated at breakpoint positions
- `svStructures.svfdel[i].used`, `svStructures.svrdel[i].used` — marked as used
- `softClips3End` / `softClips5End` — indirectly through `markSV`
- `reference` — may be extended via `referenceResource.getReference()`

#### Null/Edge Cases
- `soft.isEmpty()` → `softp = 0`, falls to the "no softp" path
- `softClips3End/5End.containsKey(softp)` check before `.get()` — skip if not present
- `ref.containsKey(bp)` and `ref.containsKey(softp - 1)` checked in left-alignment loop
- `nonInsertionVariants.get(softp)` checked with `containsKey` before accessing reference variant in `adjCnt`
- Exception caught per iteration → `printExceptionAndContinue` — continues to next cluster

#### Parity Warnings
- **Integer division for fwd/rev split**: `mcnt / 2` and `mcnt - mcnt / 2` — this is integer division truncating toward zero
- **The 5' happy-reads path breaks after first match** — this is intentional
- **Left-alignment loop** uses `ref.get(bp).equals(ref.get(softp - 1))` with boxed Character comparison (`.equals()`) — `==` would fail for values > 127
- **`dellen` calculated before left-alignment** but the variation string uses the original `dellen` value (not re-calculated after shift)
- **The 3' no-softp path calls `incCnt`** while the 5' path does not — asymmetry that must be preserved
- **`variation.varsCount = 0`** is explicitly set after `getVariation()` creates a default Variation — this resets the count

---

### Analysis: StructuralVariantsProcessor.findINV()

**Source**: StructuralVariantsProcessor.java:L686–L690
**Purpose**: Dispatch INV finding to `findINVsub` for all 4 orientation lists
**Called By**: `findAllSVs()`

#### Algorithm
1. `findINVsub(svStructures.svfinv5, 1, Side._5)`
2. `findINVsub(svStructures.svrinv5, -1, Side._5)`
3. `findINVsub(svStructures.svfinv3, 1, Side._3)`
4. `findINVsub(svStructures.svrinv3, -1, Side._3)`

#### Parity Warnings
- Exactly 4 calls in this exact order. Return value of each is a `Variation` (or null) but is discarded.

---

### Analysis: StructuralVariantsProcessor.findINVsub()

**Source**: StructuralVariantsProcessor.java:L700–L892
**Purpose**: Find INV structural variants from one direction/side combination
**Called By**: `findINV()` (4 calls)

#### Parameters
- `svref`: `Iterable<Sclip>` — one of the 4 INV lists
- `dir`: `int` — 1 for forward (3' soft-clip), -1 for reverse (5' soft-clip)
- `side`: `Side` — `_5` or `_3`

#### Algorithm (Step-by-Step)

For each `inv` in `svref`:
1. Skip if `inv.used` or `inv.varsCount < instance().conf.minr`
2. Sort `inv.soft` by value descending → `softp`
3. Select `sclip = dir == 1 ? softClips3End : softClips5End`
4. Ensure reference loaded for `inv.mstart..inv.mend`, run partial pipeline if needed
5. **If `softp != 0`**:
   a. Skip if `sclip` doesn't contain `softp`, or `scv.used`
   b. `seq = findconseq(scv, 0)` → skip if empty
   c. `findMatchRev(seq, reference, softp, dir)` → if bp=0, retry with `SEED_2, 0`
   d. Skip if bp=0
6. **Else** (no softp — scan within `2 * maxReadLength`):
   a. `sp = dir == 1 ? inv.end : inv.start`
   b. For `i = 1` to `2 * maxReadLength`:
      - `cp = sp + i * dir`
      - Skip if `sclip` doesn't contain `cp`, or `scv.used`
      - `findMatchRev(seq, reference, cp, dir)` — retry with `SEED_2, 0` if bp=0
      - If bp found: `softp = cp`; break if distance check passes (`MINSVCDIST * maxReadLength`)
   c. Skip if bp still 0

7. **Position adjustment** (complex, direction/side dependent):
   - If `side == _5` and `dir == -1`: `bp--`
   - If `side == _3` and `dir == 1`: `bp++; softp--` (if bp != 0)
   - If `side == _3` and `dir == -1`: `softp--`
   - If `side == _3`: swap `bp` and `softp`

8. **Complement-based left-alignment**:
   - If `(dir == -1 && side == _5) || (dir == 1 && side == _3)`:
     While `ref[softp] == complement(ref[bp])`: `softp++; bp--`
   - Then always: While `ref[softp-1] == complement(ref[bp+1])`: `softp--; bp++`

9. **Check viability**: `bp > softp && bp - softp > 150 && (bp - softp) / abs(inv.mlen) < 1.5`
10. Build variant description:
    - `len = bp - softp + 1`
    - `ins5 = reverseComplement(joinRef(ref, bp - SVFLANK + 1, bp))`
    - `ins3 = reverseComplement(joinRef(ref, softp, softp + SVFLANK - 1))`
    - `ins = ins5 + "<inv" + (len - 2*SVFLANK) + ">" + ins3`
    - If `len - 2*SVFLANK <= 0`: `ins = reverseComplement(joinRef(ref, softp, bp))`
    - If `dir == 1` and extra non-empty: `extra = reverseComplement(extra); ins = extra + ins`
    - If `dir == -1` and extra non-empty: `ins = ins + extra`
    - `gt = "-" + len + "^" + ins`

11. Create variation at `nonInsertionVariants[softp][gt]`
12. Set `inv.used = true`, `vref.pstd = true`, `vref.qstd = true`
13. Get/create SV, set `type = "INV"`, accumulate counts
14. Call `adjCnt(vref, scv, vrefSoftp)` — note: `vrefSoftp` is only set when `dir == -1`
15. Build `dels5` map for `realigndel` call
16. Set `refCoverage[softp]` to `refCoverage[softp - 1]` if exists, else `inv.varsCount`
17. Mark `scv.used = true`
18. Call `variationRealigner.realigndel(bams, dels5)` — **triggers a secondary realignment pass**
19. Return `vref`

**After loop**: return `null` if no INV found

#### Mutable State Modified
- `nonInsertionVariants` — new INV entries
- `refCoverage` — updated at softp
- `inv.used`, `scv.used` — marked
- VariationRealigner invoked for secondary realignment of the deletion component

#### Null/Edge Cases
- `scv` initialized as `new Sclip()` before the softp block — if softp=0 path is taken, `scv` remains this empty Sclip (varsCount=0)
- `inv.mlen` used in `abs(inv.mlen)` — could be negative
- `vrefSoftp` is `null` when `dir == 1` — `adjCnt` with null third arg just adds without reference subtraction
- `softp != 0` check in bp adjustment: `if (bp != 0)` → should be `if (softp != 0)` based on Perl original. **This is an idiomatic Java bug that exists in the source** — `bp` was just incremented from 0 potentially.

#### Parity Warnings
- **`complement` returns char (single residue)** — in alignment loops, these are single-char comparisons using `==` (primitive char, safe)
- **`reverseComplement`** is htsjdk's `SequenceUtil.reverseComplement(String)` — complement + reverse
- **The `if (bp != 0)` / `if (softp != 0)` guards** in position adjustment appear to be leftovers from Perl (where 0 is falsy). In Java, `bp` could be 0 after decrement but this path isn't reached because of the earlier `bp == 0` skip. Rust should replicate the exact same conditionals.
- **`realigndel` called here** — this is a cross-module call that can further mutate `nonInsertionVariants` and `refCoverage`
- **Return type**: returns the first successfully created `Variation` or `null` — but callers ignore the return value

---

### Analysis: StructuralVariantsProcessor.findsv()

**Source**: StructuralVariantsProcessor.java:L897–L1175
**Purpose**: Find DEL and INV SVs from raw soft-clip scanning (no pre-associated pair evidence)
**Called By**: `findAllSVs()`

#### Algorithm (Step-by-Step)

**Initialization**:
1. `ref = reference.referenceSequences`
2. `tmp5 = fillAndSortTmpSV(softClips5End.entrySet())` — sorted by count descending
3. `tmp3 = fillAndSortTmpSV(softClips3End.entrySet())` — sorted by count descending

**5' scanning** (L907–L1005):

For each `tuple5` in `tmp5`:
1. Break if `cnt5 < instance().conf.minr` (list is sorted, so all remaining are below threshold)
2. Skip if `sc5v.used`
3. Skip if `SOFTP2SV.containsKey(p5) && SOFTP2SV.get(p5).get(0).used`
4. `seq = findconseq(sc5v, 0)` → skip if empty or `< SEED_2`
5. `match = findMatch(seq, reference, p5, -1)` — **direction -1 for 5' clip**
6. **If `bp != 0`**:
   a. **If `bp < p5`** — candidate deletion:
      - `checkPairs(chr, bp, p5, [svfdel, svrdel], maxReadLength)` → get pairs/stats
      - Skip if `pairs == 0`
      - `p5--; bp++`
      - `dellen = p5 - bp + 1`
      - Create variation `nonInsertionVariants[bp]["-" + dellen]`, set `varsCount = 0`
      - Create SV type "DEL"
      - Update refCoverage, adjCnt with sc5v and pair evidence
   b. **If `bp >= p5`** — candidate duplication: **empty block** (TODO in Java)

7. **Else** (`bp == 0`) — candidate inversion:
   a. `findMatchRev(seq, reference, p5, -1)` → get bp
   b. Skip if `bp == 0` or `abs(bp - p5) <= SVFLANK`
   c. Ensure `bp > p5` (swap if needed, but note: for 5' clip with bp < p5, swap makes bp=p5_orig, p5=bp_orig)
   d. If `bp < p5` before swap → NEVER reached since bp is set to p5 original (no-op swap)
   e. Actually: "bp > p5: bp at 3' side" is the if branch, "bp < p5: bp at 5' side" → swap bp and p5
   f. `bp--`
   g. Left-align: while `ref[bp+1] == complement(ref[p5-1])`: `p5--; bp++`
   h. Build INV variant string with `ins5`, `ins3`, `mid`
   i. Format: `"-" + (bp - p5 + 1) + "^" + ins5 + "<inv" + mid + ">" + ins3 + EXTRA`
   j. If `mid <= 0`: `"-" + (bp - p5 + 1) + "^" + tins + EXTRA` where `tins = reverseComplement(joinRef(ref, p5, bp))`
   k. Create variation and SV, `adjCnt`, `incCnt` refCoverage

**3' scanning** (L1005–L1175):

For each `tuple3` in `tmp3`:
1. Same break/skip logic as 5'
2. `findMatch(seq, reference, p3, 1)` — **direction 1 for 3' clip**
3. **If `bp != 0`**:
   a. **If `bp > p3`** — candidate deletion:
      - `checkPairs`, skip if pairs==0
      - `dellen = bp - p3`
      - `bp--`
      - Left-align: while `ref[bp] == ref[p3 - 1]`: `bp--; p3--`
      - Create DEL variation at `p3`
   b. **If `bp <= p3`** — candidate duplication: **empty block**

4. **Else** — candidate inversion:
   a. `findMatchRev(seq, reference, p3, 1)`
   b. Skip if bp=0 or `abs(bp - p3) <= SVFLANK`
   c. If `bp < p3`: swap p3/bp, then `p3++; bp--`
   d. Left-align with complement
   e. Build variant string: `"-" + (bp - p3 + 1) + "^" + EXTRA + ins5 + "<inv" + mid + ">" + ins3`
   f. **Note**: for 3', EXTRA is prepended to ins5 (before `<inv>`), whereas for 5' it's appended after ins3

#### Parity Warnings
- **Uses `break` on count threshold** — since list is sorted, `cnt < minr` means all remaining are also below
- **SOFTP2SV check**: accesses `get(0)` on the list — assumes non-empty if key exists
- **Empty duplication branches** — these are TODO/not-implemented in Java
- **5' INV EXTRA placement differs from 3'** — `ins5 + "<inv...>" + ins3 + EXTRA` (5') vs `EXTRA + ins5 + "<inv...>" + ins3` (3')
- **3' candidate deletion has left-alignment loop but 5' candidate deletion does not** — asymmetry that must be preserved
- **`int` cast**: in 3' deletion path, `(int)(pairs / 2)` — this is a no-op since `pairs` is already `int`, but the explicit cast is present
- **`mid` computation differs**: 5' mid = `bp - p5 - ins5.length() - ins3.length() + 1`; 3' mid = `bp - p3 - 2 * SVFLANK + 1`

---

### Analysis: StructuralVariantsProcessor.fillAndSortTmpSV()

**Source**: StructuralVariantsProcessor.java:L1177–L1195
**Purpose**: Filter soft-clip entries to current segment and sort by count descending
**Called By**: `findsv()`

#### Parameters
- `entries`: `Set<Map.Entry<Integer, Sclip>>` — soft-clip map entry set

#### Algorithm
1. For each entry:
   - Skip if `sclip.used`
   - Skip if `position < CURSEG.start || position > CURSEG.end`
   - Add `SortPositionSclip(position, sclip, sclip.varsCount)` to list
2. Sort by `count` descending

#### Parity Warnings
- **HashMap iteration order**: `softClips5End`/`softClips3End` are `Map<Integer, Sclip>` — iteration order is non-deterministic for equal-count entries in Java. The sort is stable on `count` only, not on position as tiebreaker. Rust must match this behavior (or ensure identical tie-breaking).
- **CURSEG boundary filtering** — only processes clips within the current segment window

---

### Analysis: StructuralVariantsProcessor.findDELdisc()

**Source**: StructuralVariantsProcessor.java:L1200–L1430
**Purpose**: Find DEL SVs using discordant pairs only (no soft-clip consensus alignment)
**Called By**: `findAllSVs()`

#### Algorithm (Step-by-Step)

**Constants**: `MINDIST = 8 * maxReadLength`

**Forward loop** (`svfdel`, L1207–L1310):

For each `del` in `svStructures.svfdel`:
1. Skip if `del.used`
2. Skip if `!splice.isEmpty() && abs(del.mlen) < 250000` — RNA-Seq stringency
3. Skip if `del.varsCount < instance().conf.minr + 5` — **note: minr + 5**, stricter than soft-clip path
4. Skip if `del.mstart <= del.end + MINDIST`
5. Skip if `del.meanMappingQuality / del.varsCount <= DISCPAIRQUAL` (35)
6. `mlen = del.mstart - del.end - maxReadLength / (del.varsCount + 1)` — **integer division**
7. Skip if NOT `(mlen > 0 && mlen > MINDIST)`
8. `bp = del.end + (maxReadLength / (del.varsCount + 1)) / 2`
9. If `del.softp != 0`: `bp = del.softp`
10. Ensure reference loaded at bp
11. Create variation `nonInsertionVariants[bp]["-" + mlen]`, set varsCount = 0
12. Create SV type "DEL", accumulate splits from soft-clip counts at `del.end + 1` and `del.mstart`
13. Create temporary Variation with `2 * del.varsCount` — **doubled** because discordant pairs represent two reads
14. All stats doubled: `tv.meanQuality = 2 * del.meanQuality`, etc.
15. `adjCnt(vref, tv)`
16. Set `refCoverage[bp]` if not present
17. Mark del.used, call `markSV`

**Reverse loop** (`svrdel`, L1312–L1428):

Similar logic but:
- `mlen = del.start - del.mend - maxReadLength / (del.varsCount + 1)`
- `bp = del.mend + (maxReadLength / (del.varsCount + 1)) / 2`
- If `del.softp != 0 && softClips5End.containsKey(del.softp)`: mark that soft-clip as used
- Reference coverage also checks `refCoverage[del.start]` for max
- Calls `referenceResource.getReference` at `del.mstart - 100..del.mend + 100` before markSV

#### Parity Warnings
- **`minr + 5`** threshold — different from the `minr` threshold used in findDEL
- **Integer division**: `maxReadLength / (del.varsCount + 1)` uses Java integer division (truncate toward zero)
- **Doubled stats**: the `2 *` multiplier is applied to all quality/position sums, not means
- **RNA-Seq check**: `!splice.isEmpty()` gates a minimum abs(mlen) of 250000
- **`DISCPAIRQUAL`** = 35 — this is `meanMappingQuality / varsCount` (sum/count = actual mean), threshold is on mean quality

---

### Analysis: StructuralVariantsProcessor.findINVdisc()

**Source**: StructuralVariantsProcessor.java:L1435–L1565
**Purpose**: Find INV SVs using discordant pairs only (nested loop matching forward+reverse clusters)
**Called By**: `findAllSVs()`

#### Algorithm (Step-by-Step)

**5' INV** (`svfinv5` × `svrinv5`, L1440–L1510):

For each `invf5` in `svStructures.svfinv5`:
1. Skip if used
2. Extract counts and stats: `cnt, me, ms, end, start, nm, pmean, qmean, Qmean`
3. Skip if `Qmean / cnt <= DISCPAIRQUAL`
4. For each `invr5` in `svStructures.svrinv5`:
   a. Skip if used
   b. Skip if `rQmean / rcnt <= DISCPAIRQUAL`
   c. Skip if `cnt + rcnt <= minr + 5`
   d. If `isOverlap(end, me, rstart, rms, maxReadLength)`:
      - `bp = abs((end + rstart) / 2)` — **integer division, then abs**
      - `pe = abs((me + rms) / 2)` — same
      - Ensure reference loaded at pe
      - `len = pe - bp + 1`
      - Build INV string: `ins3 + "<inv" + (len - 2*SVFLANK) + ">" + ins5`
      - **Note**: for 5' disc, `ins5 = revComp(bp..bp+SVFLANK-1)`, `ins3 = revComp(pe-SVFLANK+1..pe)`
      - If `len - 2*SVFLANK <= 0`: full revComp
      - Create variation, mark both used
      - Set `vref.pstd = true`, `vref.qstd = true`
      - Create tmp Variation with combined counts
      - SV splits from `softClips5End[start]` and `softClips5End[ms]`
      - `markSV(bp, pe, [svfinv3, svrinv3], maxReadLength)` — marks _3 clusters

**3' INV** (`svfinv3` × `svrinv3`, L1512–L1563):

Same pattern:
- `pe = abs((end + rstart) / 2)` and `bp = abs((me + rms) / 2)` — **reversed assignment compared to 5'**
- INV string: same `ins3 + "<inv...>" + ins5` pattern
- SV splits from `softClips3End[end + 1]` and `softClips3End[me + 1]`
- `markSV(bp, pe, [svfinv5, svrinv5], maxReadLength)` — marks _5 clusters

#### Parity Warnings
- **`abs()` on integer division result** — `(end + rstart) / 2` is integer division; `abs` redundant if both positive but matches Java exactly
- **5' uses `softClips5End`; 3' uses `softClips3End`** — must not confuse
- **ins5/ins3 naming is confusing**: `ins5` is the revComp of the 5' flank, `ins3` of the 3' flank — their order in the string differs from findINVsub
- **Nested try/catch**: both inner and outer loops have separate exception handlers
- **5' marks 3' clusters; 3' marks 5' clusters** — cross-marking

---

### Analysis: StructuralVariantsProcessor.findDUPdisc()

**Source**: StructuralVariantsProcessor.java:L1570–L1825
**Purpose**: Find DUP SVs using discordant pairs and optional soft-clip refinement
**Called By**: `findAllSVs()`

#### Algorithm (Step-by-Step)

**Forward DUP** (`svfdup`, L1575–L1700):

For each `dup` in `svStructures.svfdup`:
1. Skip if used
2. Skip if `cnt < minr + 5` or `Qmean / cnt <= DISCPAIRQUAL`
3. Calculate: `mlen = end - ms + maxReadLength / cnt`
4. `bp = ms - (maxReadLength / cnt) / 2`
5. `pe = end`
6. Ensure reference loaded
7. Initialize stats: `cntf = cntr = cnt`, etc.
8. **Soft-clip refinement** (if `!dup.soft.isEmpty()`):
   a. Sort `dup.soft` by value descending → `pe = soft[0].getKey()`
   b. Check `softClips3End.containsKey(pe)` and not used
   c. Get consensus, `findMatch(seq, reference, bp, 1)`
   d. If `tbp != 0 && tbp < pe`:
      - Mark sclip3 as used
      - Left-align: while `ref[pe-1] == ref[tbp-1]`: `tbp--; pe--`
      - `mlen = pe - tbp; bp = tbp; pe--; end = pe`
      - If `softClips5End[bp]` exists → extract stats for reverse counts
9. Build DUP string: `ins5 + "<dup" + (mlen - 2*SVFLANK) + ">" + ins3`
10. **DUP goes into `insertionVariants`**: `getVariation(insertionVariants, bp, "+" + ins)` — **different from DEL/INV which use nonInsertionVariants**
11. SV goes into `nonInsertionVariants` (getSV): type "DUP"
12. `tmp.extracnt = tcnt` — **DUP sets `extracnt` explicitly** (unlike DEL/INV)
13. `adjCnt(vref, tmp)`
14. `markDUPSV` — returns `(clusters, pairs)`, clusters added to sv

**Reverse DUP** (`svrdup`, L1703–L1823):

Mirror logic:
- `mlen = me - start + maxReadLength / cnt`
- `bp = start - (maxReadLength / cnt) / 2`
- `pe = mlen + bp - 1`, `tpe = pe`
- Soft-clip refinement uses `softClips5End` for bp position
- `findMatch(seq, reference, pe, -1)` with direction -1
- If match found: `pe = tbp; mlen = pe - bp + 1; tpe = pe + 1`
- Extends `tpe` while `ref[tpe] == ref[bp + (tpe - pe - 1)]`
- If `softClips3End[tpe]` exists → extract stats for forward counts

#### Parity Warnings
- **DUP variant in `insertionVariants`** — DEL and INV go into `nonInsertionVariants`
- **SV struct in `nonInsertionVariants`** — even though the variant is in insertionVariants, the SV metadata is stored in nonInsertionVariants
- **`extracnt` set explicitly** for DUP — this affects downstream frequency calculations in ToVarsBuilder
- **`vref.varsCount = 0`** explicit reset on DUP variants
- **`markDUPSV`** returns `(cnt, pairs)` not `(cov, cnt, pairs)` like `markSV` — different helper
- **Integer division**: `maxReadLength / cnt` — potential divide-by-zero if cnt=0, but guarded by `cnt >= minr + 5`
- **Forward DUP**: `soft` sorted by value then takes `[0]` — **different from findDEL** which streams and sorts tuples. Here it creates an `ArrayList` from `entrySet()` and sorts entries.

---

### Analysis: StructuralVariantsProcessor.markSV()

**Source**: StructuralVariantsProcessor.java:L1830–L1870
**Purpose**: Mark overlapping SV clusters as used (static method)
**Called By**: `findDEL()`, `findINVdisc()`, `findDELdisc()`, `findsv()` (implicitly through flow)

#### Parameters
- `start`, `end`: int — region endpoints
- `structuralVariants_sv`: `List<List<Sclip>>` — lists of SV clusters to check
- `rlen`: int — read length

#### Algorithm
1. For each cluster list in `structuralVariants_sv`:
   a. For each `sv_r`:
      - Determine `start2, end2`:
        If `sv_r.start < sv_r.mstart`: `start2 = sv_r.end; end2 = sv_r.mstart`
        Else: `start2 = sv_r.mend; end2 = sv_r.start`
      - If `isOverlap(start, end, start2, end2, rlen)`:
        Mark `sv_r.used = true`, accumulate cnt, pairs, cov
        `cov += (int)((sv_r.varsCount * rlen) / (sv_r.end - sv_r.start)) + 1`
2. Return `(cov, cnt, pairs)`

#### Null/Edge Cases
- `sv_r.end == sv_r.start` → division by zero in cov calculation — **protected by the if** that checks for overlap, which requires a non-degenerate interval via the inner geometry checks, but **not explicitly guarded**

#### Parity Warnings
- **Integer cast in coverage**: `(int)((sv_r.varsCount * rlen) / (sv_r.end - sv_r.start))` — Java integer division
- **Return value `Tuple3<cov, cnt, pairs>`** but callers mostly ignore it (only findDEL doesn't use the return value at all)

---

### Analysis: StructuralVariantsProcessor.markDUPSV()

**Source**: StructuralVariantsProcessor.java:L1875–L1915
**Purpose**: Mark overlapping DUP SV clusters as used (static method)
**Called By**: `findDUPdisc()`

#### Parameters
Same as `markSV`

#### Algorithm
Same as `markSV` except:
- **start2/end2 assignment differs**: If `sv_r.start < sv_r.mstart`: `start2 = sv_r.start; end2 = sv_r.mend` (uses start/mend, not end/mstart)
  Else: `start2 = sv_r.mstart; end2 = sv_r.end`

#### Parity Warnings
- **Critical difference from `markSV`**: the coordinate extraction logic is different — DUP uses outer boundaries (start, mend) while DEL uses inner boundaries (end, mstart)
- Returns `Tuple2<cnt, pairs>` (no cov return) — though cov is still computed locally

---

### Analysis: StructuralVariantsProcessor.isOverlap()

**Source**: StructuralVariantsProcessor.java:L1920–L1945
**Purpose**: Determine if two SV intervals overlap based on geometric criteria
**Called By**: `findDEL`, `findINVsub`, `findsv`, `findDELdisc`, `findINVdisc`, `markSV`, `markDUPSV`, `checkPairs`

#### Parameters
- `start1, end1`: first interval
- `start2, end2`: second interval
- `rlen`: read length

#### Algorithm
1. If `start1 >= end2 || start2 >= end1`: return `false` (no geometric overlap)
2. Sort all 4 positions: `positions = sort([start1, end1, start2, end2])`
3. `ins = positions[2] - positions[1]` — inner segment overlap
4. **Full overlap check**: if both intervals non-degenerate:
   `ins / (end1 - start1) > 0.75 && ins / (end2 - start2) > 0.75` → return `true`
5. **Proximity check**: `positions[1] - positions[0] + positions[3] - positions[2] < 3 * rlen` → return `true`
6. Return `false`

#### Parity Warnings
- **Double division**: `ins / (double)(end1 - start1)` — explicitly cast to double
- **Degenerate intervals**: `end1 == start1` or `end2 == start2` skips the fraction check (avoids division by zero) and falls through to proximity check
- **Sort is on Integer objects** using `Integer.compareTo` — stable sort but order matches natural integer ordering
- Uses `Arrays.asList()` which creates a fixed-size list that is sortable

---

### Analysis: StructuralVariantsProcessor.checkPairs()

**Source**: StructuralVariantsProcessor.java:L1950–L1995
**Purpose**: Check if discordant pair clusters support a candidate SV
**Called By**: `findsv()` (for candidate deletions)

#### Parameters
- `chr`: String — chromosome name (used only for debug)
- `start, end`: int — candidate SV region
- `sv`: `List<List<Sclip>>` — SV cluster lists to search
- `maxReadLength`: int — read length

#### Algorithm
1. Init `pairs = 0`, stats to 0
2. For each cluster list, for each `svr`:
   a. Skip if `svr.used`
   b. `s = (svr.start + svr.end) / 2; e = (svr.mstart + svr.mend) / 2`
   c. Swap s/e if `s > e`
   d. Skip if `!isOverlap(start, end, s, e, maxReadLength)`
   e. **If `svr.varsCount > pairs`**: update PairsData (take this cluster's stats)
   f. Mark `svr.used = true`
3. Return PairsData with best match

#### Parity Warnings
- **Takes the cluster with highest `varsCount`** as the representative pair data, but **marks ALL overlapping clusters as used** — not just the best one
- Integer division in `(start + end) / 2` — standard truncation

---

### Analysis: StructuralVariantsProcessor.findMatchRev()

**Source**: StructuralVariantsProcessor.java:L1997–L2075 (6-arg version)
**Purpose**: Find matching position on reverse (complemented) strand using seed-based alignment
**Called By**: `findINVsub`, `findsv`

#### Parameters (6-arg version)
- `seq`: String — consensus sequence
- `REF`: Reference — reference data with seed map
- `position`: int — starting position (used only for debug)
- `dir`: int — 1 for 3' clip, -1 for 5' clip
- `SEED`: int — seed length (SEED_1=17 or SEED_2=12)
- `MM`: int — max allowed mismatches (default 3)

#### Algorithm
1. **If `dir == 1`**: `seq = reverse(seq)` (3' soft-clip → reverse first)
2. `seq = complement(seq)` — **complement, NOT reverseComplement** (reverse already done if needed)
3. For `i = seq.length() - SEED` down to `0`:
   a. `seed = substr(seq, i, SEED)`
   b. If `REF.seed.containsKey(seed)` and seed has exactly 1 hit:
      - `firstSeed = seeds[0]`
      - `bp = dir == 1 ? firstSeed + seq.length() - i - 1 : firstSeed - i`
      - If `ismatchref(seq, REF.referenceSequences, bp, -1 * dir, MM)`: return `Match(bp, "")`
      - Else: **complex indel fallback** (up to 15bp trimming):
        For `j = 1` to `15`:
        - `bp -= dir`
        - Trim seq: if `dir == -1`: `sseq = sseq[1:]`; if `dir == 1`: `sseq = sseq[0:-1]`
        - Check 2 consecutive positions match reference
        - Track `eqcnt`; break if `eqcnt >= 3 && eqcnt / j > 0.5`
        - Set `extra` to trimmed portion
        - If `ismatchref(sseq, REF, bp, -1 * dir, 1)`: return `Match(bp, extra)`
4. return `Match(0, "")`

#### Parity Warnings
- **`complement` NOT `reverseComplement`**: The sequence is complemented (A→T, C→G, etc.) but NOT reversed (or reversed separately if dir==1). This is for INV detection where the sequence maps to the complement strand.
- **`-1 * dir`** passed to `ismatchref` — reversed direction for the reference matching
- **Seeds must have exactly 1 hit** (`seeds.size() == 1`) — ambiguous seeds are skipped entirely
- **Extra sequence**: for `dir == -1`, extra is prefix of original seq; for `dir == 1`, extra is suffix
- **`eqcnt` breaking condition**: `eqcnt >= 3 && eqcnt/(double)j > 0.5` — uses double division. If this breaks, no match is returned for this seed.

---

### Analysis: StructuralVariantsProcessor.findMatch()

**Source**: StructuralVariantsProcessor.java:L2088–L2175 (6-arg static version)
**Purpose**: Find matching position on forward strand using seed-based alignment
**Called By**: `findDEL`, `findsv`, `findDUPdisc`

#### Parameters
- `seq`: String — consensus sequence
- `REF`: Reference — reference with seed map
- `position`: int — starting position (debug)
- `dir`: int — 1 for 3' clip, -1 for 5' clip
- `SEED`: int — seed length
- `MM`: int — max mismatches

#### Algorithm
1. **If `dir == -1`**: `seq = reverse(seq)` (5' clip → reverse)
2. For `i = seq.length() - SEED` down to `0`:
   a. `seed = substr(seq, i, SEED)`
   b. If seed found with exactly 1 hit:
      - `bp = dir == 1 ? firstSeed - i : firstSeed + seq.length() - i - 1`
      - **Primary match**: if `ismatchref(seq, REF.referenceSequences, bp, dir, MM)`:
        - **Extra extraction loop**: while `ref[bp] != seq[mm]` (walking backward for mismatches at boundary):
          `mm` starts at `dir == -1 ? -1 : 0`
          `extra += substr(seq, mm, 1); bp += dir; mm += dir`
        - If extra non-empty and `dir == -1`: `extra = reverse(extra)`
        - Return `Match(bp, extra)`
      - **Complex indel fallback** (up to 15bp, same logic as findMatchRev but for forward strand):
        Similar trimming loop with eqcnt breaking condition

#### Key Difference from findMatchRev
- **No complement** — operates on forward strand
- **Extra extraction** happens in successful match path (not just fallback): mismatched bases at the alignment boundary are collected as `extra`
- **bp formula**: `dir==1 ? firstSeed - i` vs findMatchRev's `dir==1 ? firstSeed + seq.length() - i - 1`

#### Parity Warnings
- **`extra` accumulation uses string concatenation**: `extra += substr(seq, mm, 1)` — Java String immutable concat creates new objects. Rust should use `String::push` or `push_str`.
- **`charAt` with negative index**: `charAt(seq, -1)` and `charAt(seq, -2)` — this is a VarDict utility that handles Python-style negative indexing. Must replicate exactly.
- **The `while` loop for extra**: condition is `isHasAndNotEquals(charAt(seq, mm), REF.referenceSequences, bp)` — this loops while there IS a reference base AND it DOESN'T match. Stops on: missing reference base OR matching base.
- **`static` method** — unlike findMatchRev which is an instance method (both call `instance().conf.y` for debug so both need config access)

---

### Analysis: StructuralVariantsProcessor.adjSNV()

**Source**: StructuralVariantsProcessor.java:L2180–L2240
**Purpose**: Rescue short (≤5bp) soft-clipped reads as SNVs if the first base matches an existing variant
**Called By**: `process()` — always runs, even when SV is disabled

#### Algorithm

**5' clips**:
For each `(position, sclip)` in `softClips5End`:
1. Skip if `sclip.used`
2. `seq = findconseq(sclip, 0)` → skip if `seq.length() > 5`
3. `bp = substr(seq, 0, 1)` — first base as string
4. `previousPosition = position - 1`
5. If `nonInsertionVariants[previousPosition]` contains key `bp`:
   a. If `seq.length() > 1`: check `reference[position - 2] != seq.charAt(1)` → skip (second base must match reference)
   b. `adjCnt(nonInsertionVariants[previousPosition][bp], sclip)`
   c. `incCnt(refCoverage, previousPosition, sclip.varsCount)`

**3' clips**:
For each `(position, sclip)` in `softClips3End`:
1. Skip if `sclip.used`
2. Same consensus/length check
3. `bp = substr(seq, 0, 1)`
4. If `nonInsertionVariants[position]` contains key `bp`:
   a. If `seq.length() > 1`: check `reference[position + 1] != seq.charAt(1)` → skip
   b. `adjCnt(nonInsertionVariants[position][bp], sclip)`
   c. `incCnt(refCoverage, position, sclip.varsCount)`

#### Parity Warnings
- **5' uses `previousPosition = position - 1`** and checks `reference[position - 2]`
- **3' uses `position` directly** and checks `reference[position + 1]`
- **`reference.referenceSequences.get(position - 2)`** can return `null` → `isNotEquals(null, seq.charAt(1))` is handled by `isNotEquals` which returns `true` if either is null but they're not equal → would NOT skip, meaning the SNV rescue proceeds
- **adjCnt call without referenceVar** — 2-arg version, no reference adjustment
- Iterates over map entries — **HashMap ordering** in Java is non-deterministic. If multiple soft-clips could match the same variant, all are processed independently.

---

### Analysis: StructuralVariantsProcessor.outputClipping()

**Source**: StructuralVariantsProcessor.java:L2245–L2285
**Purpose**: Debug output — print remaining unused soft-clipped reads to stderr
**Called By**: `process()` (only when `conf.y` is true)

#### Algorithm
1. Print "5' Remaining clipping reads"
2. For each `(position, sclip)` in `softClips5End`:
   - Skip if used or `varsCount < minr`
   - Get consensus, skip if empty or `<= SEED_2`
   - **Reverse the sequence** before printing (5' clips are stored reversed)
   - Print position, count, reversed seq
3. Print "3' Remaining clipping reads"
4. Same for `softClips3End` but **no reversal**

#### Parity Warnings
- Debug output only — does not affect parity of variant output
- The 5' sequence reversal using `StringBuilder.reverse()` is for display only

---

### Analysis: StructuralVariantsProcessor.PairsData (inner class)

**Source**: StructuralVariantsProcessor.java:L2287–L2310
**Purpose**: Simple data holder for pair evidence statistics

#### Fields
- `pairs`: int — number of pairs
- `pmean`: double — mean position
- `qmean`: double — mean base quality
- `Qmean`: double — mean mapping quality
- `nm`: double — number of mismatches

---

## Cross-Module Dependencies

### Outbound Calls

| Target | Method | Purpose |
|--------|--------|---------|
| `VariationUtils` | `findconseq(Sclip, int)` | Find consensus sequence from soft-clips |
| `VariationUtils` | `getVariation(Map, int, String)` | Get/create Variation entry |
| `VariationUtils` | `adjCnt(Variation, Variation, Variation)` | Add variant counts, optionally subtract from reference |
| `VariationUtils` | `incCnt(Map, Object, int)` | Increment coverage counter |
| `VariationUtils` | `isHasAndEquals(...)` | Reference comparison utilities (multiple overloads) |
| `VariationUtils` | `isHasAndNotEquals(...)` | Reference mismatch check |
| `VariationUtils` | `isNotEquals(...)` | Null-safe character inequality |
| `VariationMap` | `getSV(Map, int)` | Get/create SV metadata struct |
| `VariationRealigner` | `ismatchref(String, Map, int, int, int)` | Check if sequence matches reference |
| `VariationRealigner` | `realigndel(String[], Map)` | Secondary realignment of deletion component (from findINVsub) |
| `ReferenceResource` | `getReference(Region, int, Reference)` | Load/extend reference sequence |
| `ReferenceResource` | `isLoaded(String, int, int, Reference)` | Check if reference region is loaded |
| `Utils` | `complement(char)`, `complement(String)` | Nucleotide complement |
| `Utils` | `reverse(String)` | String reversal |
| `Utils` | `substr(String, int, int)` | Perl-style substring |
| `Utils` | `charAt(String, int)` | Perl-style charAt with negative indexing |
| `Utils` | `printExceptionAndContinue(...)` | Error handling |
| `Utils` | `joinRef(Map, int, int)` | Build reference string from position range |
| `SequenceUtil` | `reverseComplement(String)` | htsjdk reverse complement |
| `GlobalReadOnlyScope` | `instance().conf.*` | Configuration access (minr, y, disableSV) |
| `GlobalReadOnlyScope` | `getMode().partialPipeline(...)` | Run partial pipeline for newly loaded regions |

### Inbound Callers

| Source | Context |
|--------|---------|
| Pipeline (Mode) | Calls `process()` as part of SAMFileParser → ... → StructuralVariantsProcessor → ToVarsBuilder chain |

## Data Structures Read/Written

| Structure | Type | Java Type | R/W | Notes |
|-----------|------|-----------|-----|-------|
| `nonInsertionVariants` | Position → (VarDesc → Variation) | `Map<Integer, VariationMap<String, Variation>>` | R/W | DEL/INV variants written here. VariationMap is LinkedHashMap (insertion-ordered). |
| `insertionVariants` | Position → (VarDesc → Variation) | `Map<Integer, VariationMap<String, Variation>>` | R/W | DUP variants written here |
| `refCoverage` | Position → count | `Map<Integer, Integer>` | R/W | Updated at breakpoints |
| `softClips5End` | Position → Sclip | `Map<Integer, Sclip>` | R/W | Read for consensus; `.used` field written |
| `softClips3End` | Position → Sclip | `Map<Integer, Sclip>` | R/W | Read for consensus; `.used` field written |
| `svStructures.svfdel` | List of clusters | `List<Sclip>` | R/W | `.used` field written |
| `svStructures.svrdel` | List of clusters | `List<Sclip>` | R/W | `.used` field written |
| `svStructures.svfinv3/5` | List of clusters | `List<Sclip>` | R/W | `.used` field written |
| `svStructures.svrinv3/5` | List of clusters | `List<Sclip>` | R/W | `.used` field written |
| `svStructures.svfdup` | List of clusters | `List<Sclip>` | R/W | `.used` field written |
| `svStructures.svrdup` | List of clusters | `List<Sclip>` | R/W | `.used` field written |
| `SOFTP2SV` | Position → List\<Sclip\> | `Map<Integer, List<Sclip>>` | R | Checked for `.used` state only |
| `reference` | Reference sequences + seed map | `Reference` | R/W | Extended via `referenceResource.getReference()` |
| `CURSEG` | Current segment bounds | `CurrentSegment` | R | Used to filter clips in `fillAndSortTmpSV` |
| `splice` | Splice junction set | `Set<String>` | R | RNA-Seq check in findDELdisc |

## Known Parity Traps

1. **VariationMap is LinkedHashMap** — `nonInsertionVariants` and `insertionVariants` contain `VariationMap<String, Variation>` which extends `LinkedHashMap`. Rust must use `IndexMap` to preserve insertion ordering, as downstream ToVarsBuilder iterates these maps and order affects output.

2. **DUP variants go into `insertionVariants`, not `nonInsertionVariants`** — All other SV types (DEL, INV) go into `nonInsertionVariants`. DUP variants use key format `"+" + ins` in `insertionVariants`, but SV metadata (via `getSV`) is stored in `nonInsertionVariants` at the same position.

3. **`variation.varsCount = 0` explicit reset** — After `getVariation()` creates a new Variation (with default varsCount=0), several paths explicitly set `varsCount = 0` again. This is redundant for new entries but resets the count if the entry already existed. Rust must preserve this explicit assignment.

4. **Integer division for fwd/rev splitting** — `mcnt / 2` and `mcnt - mcnt / 2` distributes counts between forward and reverse strands. For odd `mcnt`, forward gets the smaller half. This must be exact since it affects reported strand counts.

5. **`extracnt` set explicitly for DUP but not DEL/INV** — In `findDUPdisc`, `tmp.extracnt = tcnt` is set. For DEL/INV, `adjCnt` internally sets `extracnt += varsCount`. The DUP path ends up double-counting extracnt through adjCnt.

6. **`findMatch` is static, `findMatchRev` is not** — Both need `instance().conf.y` access for debug output. This is a Java design inconsistency but doesn't affect parity.

7. **`findMatchRev` complements then matches; `findMatch` operates on forward strand** — The `complement()` call (without reverse) in `findMatchRev` is NOT `reverseComplement`. It only swaps A↔T, C↔G base-by-base. The sequence may be reversed separately if `dir == 1`.

8. **Seeds must have exactly 1 occurrence** — Both `findMatch` and `findMatchRev` skip seeds with `seeds.size() != 1`. This means ambiguous reference regions will fail to find matches even if the correct position is among the hits.

9. **Left-alignment differs between methods** — `findDEL` left-aligns by comparing `ref[bp]` to `ref[softp-1]`; `findINVsub` left-aligns by comparing `ref[softp]` to `complement(ref[bp])`. The `findsv` 3' DEL path left-aligns differently from the 5' path (which doesn't left-align at all for candidate deletions).

10. **`realigndel` called from `findINVsub`** — This cross-module call creates a new `VariationRealigner` and calls `initFromScope(previousScope)` then `realigndel(bams, dels5)`. This can further modify `nonInsertionVariants` and `refCoverage`. The Rust port must ensure this secondary realignment happens with the correct scope data.

11. **`checkPairs` marks ALL overlapping clusters as used but returns only the BEST one's stats** — The function iterates all clusters, marking every overlapping one as used, but only keeps the statistics from the cluster with the highest `varsCount`. This means subsequent calls to `checkPairs` will find fewer available clusters.

12. **Exception handling continues iteration** — Every SV-finding loop wraps each iteration in try/catch and calls `printExceptionAndContinue`. This means any single cluster failure doesn't abort the entire SV search. Rust should handle errors per-iteration similarly (likely with `.ok()` or match on Result).

13. **`SOFTP2SV.get(p).get(0).used`** — Accesses the first element of a list without size check. The code assumes the list is non-empty if the key exists.

14. **`Sclip.soft` is a LinkedHashMap** — The `soft` field maps positions to counts. In `findDEL`, it's streamed/sorted; in `findDUPdisc`, it's dumped to ArrayList and sorted. Both iterate in descending value order. Rust must use `IndexMap` or equivalent for the `soft` field.

15. **Double multiplication before division in stat propagation** — `tv.meanQuality = del.meanQuality * mcnt / del.varsCount` uses Java's left-to-right evaluation: `(del.meanQuality * mcnt) / del.varsCount`. Since `meanQuality` is `double` and `mcnt`/`varsCount` are `int`, this is double×int→double, then double/int→double. No integer truncation risk, but Rust must match the evaluation order.

16. **`markSV` vs `markDUPSV` coordinate logic** — `markSV` extracts inner coordinates (end, mstart) while `markDUPSV` extracts outer coordinates (start, mend). Mixing them up would cause incorrect overlap detection.

17. **`isOverlap` sorts positions as `Integer` objects** — Uses `Arrays.asList()` which creates `List<Integer>` sorted by `Integer.compareTo`. This is natural integer ordering but worth noting for Rust where sorting a `Vec<i32>` is equivalent.

18. **MINDIST in findDELdisc is `8 * maxReadLength`** — This is a local constant, not a Configuration field. Hardcoded multiplier.

19. **`findconseq` with `dir=0`** — Throughout this module, `findconseq` is always called with `dir=0`, which means adaptor checking is never triggered (it only activates for `dir==3` or `dir==5`).

20. **`bp != 0` guards in position adjustment** — In `findINVsub`, conditions like `if (bp != 0) { softp-- }` appear to be vestigial from Perl's truthiness. In Java, `bp` could theoretically be 0 after adjustment, but the earlier `bp == 0` skip should prevent reaching these points with bp=0. Rust should replicate the exact conditional.

21. **`findsv()` pairs==0 triggers `continue`, NOT INV fallthrough** — In both 5' (L943) and 3' (L1055) paths of `findsv()`, when forward match succeeds (`bp != 0`) but `checkPairs()` returns `pairs == 0`: Java executes `continue`, skipping the iteration entirely. No DEL is created, and critically, no INV is attempted. `bp == 0` (forward match failure) is the **ONLY** trigger for the INV path via `findMatchRev()`. Rust must not set `bp = 0` and fall through to the INV path when pairs==0 — that creates rogue INVs not present in Java output.

22. **`findsv()` INV path has NO `checkPairs()` call** — When the forward match fails (`bp == 0`) and `findMatchRev()` succeeds, Java creates the INV unconditionally with `sv.pairs += 0`. There is no discordant pair validation for inversions discovered in `findsv()`. The SV tag is always `"splits-0-0"`. This differs from `findDEL()`/`findINVsub()` which do validate pair evidence.

23. **Seed map is cumulative across `getReference()` calls** — Java's `Reference.seed` HashMap accumulates seeds from every `getReference()` call. Seeds added during `findDEL()` and `findINV()` reference extensions remain visible to `findsv()`'s `findMatch()` calls. This enables cross-region DEL matching where `findsv()` finds a forward match at a distant position seeded by a prior `findDEL()` extension. Rust must ensure seed maps are cumulative (not per-window) to replicate this behavior. The 3-line DEL output for chr2 shard 158 depends on this cumulative seeding.

24. **`findAllSVs()` execution order creates seed dependencies** — Because the seed map is cumulative (trap 23) and `findAllSVs()` runs `findDEL → findINV → findsv → findDELdisc → findINVdisc → findDUPdisc` in fixed order (trap noted in findAllSVs analysis), earlier routines' reference extensions seed the map for later routines. Changing execution order or isolating seed maps per routine would break parity.
