# OutputVariant Printers

**Source**: `printers/*.java` + `postprocessmodules/*.java`  
**Total LOC**: ~904 (printers) + ~600 (post-process modules)  
**Rust counterpart**: `src/mods/output_variant.rs`  
**Status**: complete

## Overview

The OutputVariant printer subsystem is the **final output stage** of the VarDict pipeline. It converts in-memory `Variant` objects into tab-delimited text lines written to stdout. Three layers:

### Layer 1: OutputVariant data classes (data holders)
- **`OutputVariant`** (abstract base) — Common fields: sample, gene, chr, positions, alleles, shift3, msi, msint, leftseq, rightseq, region, vartype, DEBUG. Default delimiter is `\t`.
- **`SimpleOutputVariant extends OutputVariant`** — Single-sample mode. 36 columns (no fisher) or 38 columns (with `--fisher`). Optional crispr column. Optional DEBUG column.
- **`SomaticOutputVariant extends OutputVariant`** — Tumor/normal paired mode. 55 columns (no fisher) or 61 columns (with `--fisher`). Optional DEBUG column.
- **`AmpliconOutputVariant extends OutputVariant`** — Amplicon mode. 38 columns (no fisher) or 40 columns (with `--fisher`). Optional DEBUG column.

### Layer 2: VariantPrinter (output channel)
- **`VariantPrinter`** (abstract) — Holds `PrintStream out`. `print(OutputVariant)` calls `out.println(variant.toString())`.
- **`SystemOutVariantPrinter`** — `out = System.out`.
- **`SystemErrVariantPrinter`** — `out = System.err`.
- **`PrinterType`** — Enum: `OUT`, `ERR`.

### Layer 3: PostProcessModules (callers)
- **`SimplePostProcessModule`** — Iterates `alignedVariants` map, creates `SimpleOutputVariant` per variant.
- **`SomaticPostProcessModule`** — Compares BAM1 vs BAM2 variants, determines labels (StrongSomatic, Germline, etc.), creates `SomaticOutputVariant`.
- **`AmpliconPostProcessModule`** — Aggregates across amplicons, determines good/bad variants + bias flag, creates `AmpliconOutputVariant`.

### Column Count Summary

| Mode | Base Columns | With --fisher | With --crispr | With --debug |
|------|-------------|---------------|---------------|-------------|
| Simple | 36 | 38 (+2: pvalue, oddratio) | +1 (crispr) | +1 (DEBUG) |
| Somatic | 55 | 61 (+6: pvalue1, oddratio1, pvalue2, oddratio2, pvalue, oddratio) | N/A | +1 (DEBUG) |
| Amplicon | 38 | 40 (+2: pvalue, oddratio) | N/A | +1 (DEBUG) |

---

## Method Inventory

### OutputVariant.java
| Method | Lines | Analyzed? | Summary |
|--------|-------|-----------|---------|
| setDelimiter(String) | L31-L33 | yes | Sets delimiter (default `\t`) |

### SimpleOutputVariant.java
| Method | Lines | Analyzed? | Summary |
|--------|-------|-----------|---------|
| constructor(Variant, Region, String, int) | L42-L103 | yes | Populates all fields from Variant, computes Fisher if enabled |
| toString() | L106-L116 | yes | Dispatches to 36 or 38 col formatter, appends crispr/debug |
| create_simple_variant_36columns() | L160-L198 | yes | Non-fisher: zero-check ternary + DecimalFormat |
| create_simple_variant_38columns() | L121-L157 | yes | Fisher: getRoundedValueToPrint + special hifreq/nm |

### SomaticOutputVariant.java
| Method | Lines | Analyzed? | Summary |
|--------|-------|-----------|---------|
| constructor(Variant×4, Region, String×2, String) | L62-L131 | yes | Populates from 4 Variant params, computes Fisher if enabled |
| calculateFisherSomatic(Variant, Variant) | L133-L172 | yes | 3 Fisher tests: tumor, normal, combined somatic p-value |
| toString() | L175-L184 | yes | Dispatches to 55 or 61 col formatter, appends debug |
| create_somatic_variant_55columns() | L222-L280 | yes | Non-fisher format |
| create_somatic_variant_61columns() | L189-L220 | yes | Fisher format with 6 extra columns |

### AmpliconOutputVariant.java
| Method | Lines | Analyzed? | Summary |
|--------|-------|-----------|---------|
| constructor(Variant, Region, List, List, int, int, int, boolean) | L50-L128 | yes | Populates from Variant + amplicon metadata, builds debug string |
| toString() | L131-L141 | yes | Dispatches to 38 or 40 col formatter, appends debug |
| create_amplicon_variant_38columns() | L179-L215 | yes | Non-fisher format |
| create_amplicon_variant_40columns() | L146-L177 | yes | Fisher format |
| debugAmpVariant(String, Variant) | L222-L240 | yes | Debug-mode per-amplicon variant info (space-delimited) |

### VariantPrinter.java
| Method | Lines | Analyzed? | Summary |
|--------|-------|-----------|---------|
| print(OutputVariant) | L17-L19 | yes | `out.println(variant.toString())` |
| print(OutputStream) | L26-L28 | yes | `out.print(outputStream)` |
| setOut(PrintStream) | L34-L36 | yes | Updates output stream |
| getOut() | L38-L40 | yes | Returns output stream |
| createPrinter(PrinterType) | L47-L52 | yes | Factory: creates SystemOut or SystemErr printer |

---

## Method Analyses

### SimpleOutputVariant Constructor

**Source**: `SimpleOutputVariant.java` L42-L103  
**Purpose**: Populate output fields from a single `Variant` object

#### Algorithm (Step-by-Step):

1. Set `sample = instance().sample` (from global config)
2. Set `gene = region.gene`, `chr = region.chr`
3. **If variant == null**:
   - Set `startPosition = position`, `endPosition = position`
   - All other fields remain at defaults: totalCoverage=0, variantCoverage=0, genotype="", bias="0;0", leftSequence="0" (base class), etc.
4. **If variant != null**:
   - Copy all numeric fields directly from Variant
   - `genotype = variant.genotype == null ? "0" : variant.genotype`
   - `pstd = variant.isAtLeastAt2Positions ? 1 : 0` (boolean→int)
   - `qstd = variant.hasAtLeast2DiffQualities ? 1 : 0`
   - `leftSequence = variant.leftseq.isEmpty() ? "0" : variant.leftseq`
   - `rightSequence = variant.rightseq.isEmpty() ? "0" : variant.rightseq`
   - `bias = variant.strandBiasFlag` (NO null check — safe because Variant initializes it to "0")
5. **If conf.fisher**:
   - If variant != null: `FisherExact(refFwd, refRev, varFwd, varRev)`
   - If variant == null: `FisherExact(0, 0, 0, 0)`
   - `pvalue = fisher.getPValue()`, `oddratio = fisher.getOddRatio()`
6. Set `region = region.chr + ":" + region.start + "-" + region.end`
7. Set `sv = sv.equals("") ? "0" : sv`

#### Null/Edge Cases:
- variant == null: skeleton row with position only, all counts zero, genotype=""
- variant.genotype == null → "0"
- variant.leftseq.isEmpty() → "0"
- variant.rightseq.isEmpty() → "0"
- sv == "" → "0"
- variant.strandBiasFlag: NOT null-checked (relies on Variant field initializer "0")

### SimpleOutputVariant.toString()

**Source**: `SimpleOutputVariant.java` L106-L116

1. If `!conf.fisher`: call `create_simple_variant_36columns()`
2. Else: call `create_simple_variant_38columns()`
3. If `conf.crisprCuttingSite != 0`: append `\t` + crispr (int)
4. If `conf.debug`: append `\t` + DEBUG (String)
5. Return final string

### create_simple_variant_36columns()

**Source**: `SimpleOutputVariant.java` L160-L198

**Formatting Rule**: For each double field: `value == 0 ? 0 : new DecimalFormat(pattern).format(value)`. When value is 0, the int literal `0` is appended (produces "0"). When non-zero, the DecimalFormat result is used **with all trailing zeros preserved**.

Special case for nm: `nm > 0 ? new DecimalFormat("0.0").format(nm) : 0` — uses a **positive check** instead of zero check, so negative nm also outputs "0".

### create_simple_variant_38columns()

**Source**: `SimpleOutputVariant.java` L121-L157

**Formatting Rule**: Uses `getRoundedValueToPrint(pattern, value)` for most doubles, EXCEPT:
- **hifreq**: Pre-formatted as `hifreq_f = hifreq == 0 ? "0" : new DecimalFormat("0.0000").format(hifreq)` — NOT using getRoundedValueToPrint!
- **nm**: Pre-clamped to 0 if ≤ 0, then `nm_f = nm == 0 ? "0" : new DecimalFormat("0.0").format(nm)` — NOT using getRoundedValueToPrint!

### SomaticOutputVariant Constructor

**Source**: `SomaticOutputVariant.java` L62-L131

#### Algorithm (Step-by-Step):

1. Set sample, gene, chr from global config / region
2. **From beginVariant** (if != null): startPosition, endPosition, refAllele, varAllele, varType, DEBUG
3. **From endVariant** (if != null): shift3, msi, msint, leftSequence (`empty→"0"`), rightSequence (`empty→"0"`)
4. **From tumorVariant** (if != null): all var1* fields
   - `var1genotype = tumorVariant.genotype == null ? "0" : tumorVariant.genotype`
   - `var1strandBiasFlag = tumorVariant.strandBiasFlag == null ? "0" : tumorVariant.strandBiasFlag`
   - `var1isAtLeastAt2Position = tumorVariant.isAtLeastAt2Positions ? 1 : 0`
   - `var1hasAtLeast2DiffQualities = tumorVariant.hasAtLeast2DiffQualities ? 1 : 0`
5. **From normalVariant** (if != null): all var2* fields (same mapping pattern)
6. Set varLabel, region (`chr:start-end`), var1sv/var2sv (`""→"0"`)
7. If conf.fisher: `calculateFisherSomatic(tumorVariant, normalVariant)`

#### Critical: Constructor Parameter Ordering
```java
SomaticOutputVariant(Variant beginVariant, Variant endVariant, Variant tumorVariant, Variant normalVariant, ...)
```
Callers frequently pass the **SAME Variant object** for multiple parameters:
```java
// StrongSomatic: same variant for begin, end, AND tumor
new SomaticOutputVariant(vref, vref, vref, varForPrint, region, v1.sv, v2.sv, STRONG_SOMATIC)
```

### calculateFisherSomatic(Variant, Variant)

**Source**: `SomaticOutputVariant.java` L133-L172

1. **Fisher test 1 (tumor strand bias)**:
   - If tumorVariant != null: `FisherExact(refFwd, refRev, varFwd, varRev)`
   - Else: `FisherExact(0,0,0,0)`
   - Store `pvalue1`, `oddratio1`
2. **Fisher test 2 (normal strand bias)**:
   - If normalVariant != null: `FisherExact(refFwd, refRev, varFwd, varRev)`
   - Else: `FisherExact(0,0,0,0)`
   - Store `pvalue2`, `oddratio2`
3. **Fisher test 3 (somatic significance)**:
   - `tref = var1totalCoverage - var1variantCoverage` (clamped to 0 if negative)
   - `rref = var2totalCoverage - var2variantCoverage` (clamped to 0 if negative)
   - `FisherExact(var1variantCoverage, tref, var2variantCoverage, rref)`
   - If `pvalue_less < pvalue_greater`: `pvalue = pvalue_less`; else `pvalue = pvalue_greater`
   - Store `oddratio`

### AmpliconOutputVariant Constructor

**Source**: `AmpliconOutputVariant.java` L50-L128

Similar to SimpleOutputVariant but with amplicon-specific fields and different region logic:

1. Same field population from Variant as Simple
2. Additional: `goodVariantsCount`, `totalVariantsCount (= gvscnt + badVariants.size())`, `noCoverage`, `ampliconFlag`
3. **Region handling differs from Simple/Somatic**:
   - variant == null: `region = chr:position-position`
   - variant != null AND goodVariants non-empty: `region = goodVariants.get(0)._2` (first good variant's amplicon region)
   - variant != null AND goodVariants empty: `region = chr:position-position`
4. Debug mode: builds per-amplicon debug output iterating goodVariants/badVariants, each formatted via `debugAmpVariant(" ", variant)`

### debugAmpVariant(String delm, Variant)

**Source**: `AmpliconOutputVariant.java` L222-L240

Formats variant stats with given delimiter (typically space):

| Field | Format | Zero Handling |
|-------|--------|--------------|
| totalPosCoverage | raw int | — |
| positionCoverage | raw int | — |
| refForwardCoverage | raw int | — |
| refReverseCoverage | raw int | — |
| varsCountOnForward | raw int | — |
| varsCountOnReverse | raw int | — |
| genotype | raw | null→"0" |
| frequency | Z "0.0000" | ==0→"0" |
| strandBiasFlag | raw | no null check |
| meanPosition | DecimalFormat("0.0") | **ALWAYS formatted (even 0→"0.0")** |
| pstd | bool→1/0 | — |
| meanQuality | DecimalFormat("0.0") | **ALWAYS formatted** |
| qstd | bool→1/0 | — |
| meanMappingQuality | DecimalFormat("0.0") | **ALWAYS formatted** |
| highQualityToLowQualityRatio | DecimalFormat("0.000") | **ALWAYS formatted** |
| hifreq | Z "0.0000" | ==0→"0" |
| extrafreq | Z "0.0000" | ==0→"0" |

**Parity Trap**: pmean, qual, mapq, qratio are ALWAYS formatted with DecimalFormat in debugAmpVariant (even when 0), producing "0.0"/"0.000". This differs from the main column formatters which output "0" for zero values.

---

## Column Specification

### Key Utility Functions

#### `Utils.getRoundedValueToPrint(String pattern, double value)`
**Source**: `Utils.java` L105-L108
```java
return value == Math.round(value)
        ? new DecimalFormat("0").format(value)
        : new DecimalFormat(pattern).format(value).replaceAll("0+$", "");
```
- Whole number (value == Math.round(value)) → integer string without decimals ("0", "1", "42")
- Non-whole → DecimalFormat with pattern, then **strip trailing zeros** ("0.1234" kept, "0.12340" → "0.1234", "0.10000" → "0.1")

#### `Utils.join(String delim, Object... args)`
**Source**: `Utils.java` L40-L53
```java
for (int i = 0; i < args.length; i++) {
    sb.append(args[i]);
    if (i + 1 != args.length) sb.append(delim);
}
```
Each arg converted via `Object.toString()`. Primitives autoboxed. `int 0` → "0". `String "0"` → "0".

### Format Legend for Column Tables

- **Z** = Zero-check ternary: `value == 0 ? 0 : new DecimalFormat(pattern).format(value)` — preserves trailing zeros when non-zero
- **P** = Positive-check: `nm > 0 ? DecimalFormat(pattern).format(nm) : 0` — negative values also emit "0"
- **R** = `getRoundedValueToPrint(pattern, value)` — strips trailing zeros, whole numbers become integers
- **S** = Special pre-format (documented per field)

---

### Simple Mode — 36 Columns (no --fisher)

| Col# | Name | Java Type | Source Field | Format | Pattern | Zero/Null Output |
|------|------|-----------|-------------|--------|---------|-----------------|
| 1 | sample | String | `instance().sample` | raw | — | — |
| 2 | gene | String | `region.gene` | raw | — | may be null |
| 3 | chr | String | `region.chr` | raw | — | — |
| 4 | startPosition | int | `variant.startPosition` / `position` | raw | — | 0 if variant==null |
| 5 | endPosition | int | `variant.endPosition` / `position` | raw | — | 0 if variant==null |
| 6 | refAllele | String | `variant.refallele` | raw | — | "" if variant==null |
| 7 | varAllele | String | `variant.varallele` | raw | — | "" if variant==null |
| 8 | totalCoverage | int | `variant.totalPosCoverage` | raw | — | 0 |
| 9 | variantCoverage | int | `variant.positionCoverage` | raw | — | 0 |
| 10 | referenceForwardCount | int | `variant.refForwardCoverage` | raw | — | 0 |
| 11 | referenceReverseCount | int | `variant.refReverseCoverage` | raw | — | 0 |
| 12 | variantForwardCount | int | `variant.varsCountOnForward` | raw | — | 0 |
| 13 | variantReverseCount | int | `variant.varsCountOnReverse` | raw | — | 0 |
| 14 | genotype | String | `variant.genotype` | raw | — | null→"0"; variant==null→**""** |
| 15 | frequency | double | `variant.frequency` | **Z** | "0.0000" | ==0→"0" |
| 16 | bias | String | `variant.strandBiasFlag` | raw | — | default "0;0" |
| 17 | pmean | double | `variant.meanPosition` | **Z** | "0.0" | ==0→"0" |
| 18 | pstd | int | `variant.isAtLeastAt2Positions` | bool→1/0 | — | 0 |
| 19 | qual | double | `variant.meanQuality` | **Z** | "0.0" | ==0→"0" |
| 20 | qstd | int | `variant.hasAtLeast2DiffQualities` | bool→1/0 | — | 0 |
| 21 | mapq | double | `variant.meanMappingQuality` | **Z** | "0.0" | ==0→"0" |
| 22 | qratio | double | `variant.highQualityToLowQualityRatio` | **Z** | "0.000" | ==0→"0" |
| 23 | hifreq | double | `variant.highQualityReadsFrequency` | **Z** | "0.0000" | ==0→"0" |
| 24 | extrafreq | double | `variant.extraFrequency` | **Z** | "0.0000" | ==0→"0" |
| 25 | shift3 | int | `variant.shift3` | raw | — | 0 |
| 26 | msi | double | `variant.msi` | **Z** | "0.000" | ==0→"0" |
| 27 | msint | int | `variant.msint` | raw | — | 0 |
| 28 | nm | double | `variant.numberOfMismatches` | **P** | "0.0" | ≤0→"0" |
| 29 | hicnt | int | `variant.hicnt` | raw | — | 0 |
| 30 | hicov | int | `variant.hicov` | raw | — | 0 |
| 31 | leftSequence | String | `variant.leftseq` | raw | — | empty→"0"; variant==null→**""** |
| 32 | rightSequence | String | `variant.rightseq` | raw | — | empty→"0"; variant==null→**""** |
| 33 | region | String | `region.chr:region.start-region.end` | raw | — | — |
| 34 | varType | String | `variant.vartype` | raw | — | "" |
| 35 | duprate | double | `variant.duprate` | **Z** | "0.0" | ==0→"0" |
| 36 | sv | String | sv param | raw | — | ""→"0" |

**Conditional columns** (appended after col 36):

| Condition | Column | Type | Format |
|-----------|--------|------|--------|
| `crisprCuttingSite != 0` | crispr | int | raw |
| `debug` | DEBUG | String | raw |

---

### Simple Mode — 38 Columns (with --fisher)

| Col# | Name | Format | Pattern | Notes |
|------|------|--------|---------|-------|
| 1-7 | (same as 36-col) | | | |
| 8-13 | (same as 36-col) | raw int | | |
| 14 | genotype | raw | | null→"0" |
| 15 | frequency | **R** | "0.0000" | whole→int; else strip trailing 0s |
| 16 | bias | raw | | |
| 17 | pmean | **R** | "0.0" | |
| 18 | pstd | raw int | | |
| 19 | qual | **R** | "0.0" | |
| 20 | qstd | raw int | | |
| **21** | **pvalue** | **R** | "0.00000" | Fisher two-sided p-value |
| **22** | **oddratio** | raw String | | Can be "Inf", "0", or decimal |
| 23 | mapq | **R** | "0.0" | |
| 24 | qratio | **R** | "0.000" | |
| 25 | hifreq | **S** | "0.0000" | `==0 ? "0" : DecimalFormat` (NOT getRoundedValueToPrint!) |
| 26 | extrafreq | **R** | "0.0000" | |
| 27 | shift3 | raw int | | |
| 28 | msi | **R** | "0.000" | |
| 29 | msint | raw int | | |
| 30 | nm | **S** | "0.0" | Clamped to 0 if ≤0; `==0 ? "0" : DecimalFormat` |
| 31 | hicnt | raw int | | |
| 32 | hicov | raw int | | |
| 33 | leftSequence | raw | | |
| 34 | rightSequence | raw | | |
| 35 | region | raw | | |
| 36 | varType | raw | | |
| 37 | duprate | **R** | **"0.00"** | Note: "0.00" not "0.0"! |
| 38 | sv | raw | | |

**Conditional**: crispr appended if `crisprCuttingSite != 0`; DEBUG appended if `debug`.

---

### Somatic Mode — 55 Columns (no --fisher)

| Col# | Name | Source | Format | Pattern | Notes |
|------|------|--------|--------|---------|-------|
| 1 | sample | `instance().sample` | raw | | |
| 2 | gene | `region.gene` | raw | | |
| 3 | chr | `region.chr` | raw | | |
| 4 | startPosition | beginVariant.startPosition | raw | | 0 if null |
| 5 | endPosition | beginVariant.endPosition | raw | | |
| 6 | refAllele | beginVariant.refallele | raw | | "" if null |
| 7 | varAllele | beginVariant.varallele | raw | | "" if null |
| 8 | var1totalCoverage | tumorVariant.totalPosCoverage | raw | | 0 if null |
| 9 | var1variantCoverage | tumorVariant.positionCoverage | raw | | |
| 10 | var1refForwardCoverage | tumorVariant.refForwardCoverage | raw | | |
| 11 | var1refReverseCoverage | tumorVariant.refReverseCoverage | raw | | |
| 12 | var1variantForwardCount | tumorVariant.varsCountOnForward | raw | | |
| 13 | var1variantReverseCount | tumorVariant.varsCountOnReverse | raw | | |
| 14 | var1genotype | tumorVariant.genotype | raw | | null→"0"; tumor==null→"0" |
| 15 | var1frequency | tumorVariant.frequency | **Z** | "0.0000" | ==0→"0" |
| 16 | var1strandBiasFlag | tumorVariant.strandBiasFlag | raw | | null→"0"; tumor==null→"0" |
| 17 | var1meanPosition | tumorVariant.meanPosition | **Z** | "0.0" | |
| 18 | var1isAtLeastAt2Position | bool→1/0 | raw | | |
| 19 | var1meanQuality | tumorVariant.meanQuality | **Z** | "0.0" | |
| 20 | var1hasAtLeast2DiffQualities | bool→1/0 | raw | | |
| 21 | var1meanMappingQuality | tumorVariant.meanMappingQuality | **Z** | "0.0" | |
| 22 | var1highQualityToLowQualityRatio | | **Z** | "0.000" | |
| 23 | var1highQualityReadsFrequency | | **Z** | "0.0000" | |
| 24 | var1extraFrequency | | **Z** | "0.0000" | |
| 25 | var1nm | tumorVariant.numberOfMismatches | **P** | "0.0" | >0→format; else→"0" |
| 26 | var2totalCoverage | normalVariant.totalPosCoverage | raw | | 0 if null |
| 27 | var2variantCoverage | normalVariant.positionCoverage | raw | | |
| 28 | var2refForwardCoverage | normalVariant.refForwardCoverage | raw | | |
| 29 | var2refReverseCoverage | normalVariant.refReverseCoverage | raw | | |
| 30 | var2variantForwardCount | normalVariant.varsCountOnForward | raw | | |
| 31 | var2variantReverseCount | normalVariant.varsCountOnReverse | raw | | |
| 32 | var2genotype | normalVariant.genotype | raw | | null→"0" |
| 33 | var2frequency | normalVariant.frequency | **Z** | "0.0000" | |
| 34 | var2strandBiasFlag | normalVariant.strandBiasFlag | raw | | null→"0" |
| 35 | var2meanPosition | normalVariant.meanPosition | **Z** | "0.0" | |
| 36 | var2isAtLeastAt2Position | bool→1/0 | raw | | |
| 37 | var2meanQuality | normalVariant.meanQuality | **Z** | "0.0" | |
| 38 | var2hasAtLeast2DiffQualities | bool→1/0 | raw | | |
| 39 | var2meanMappingQuality | normalVariant.meanMappingQuality | **Z** | "0.0" | |
| 40 | var2highQualityToLowQualityRatio | | **Z** | "0.000" | |
| 41 | var2highQualityReadsFrequency | | **Z** | "0.0000" | |
| 42 | var2extraFrequency | | **Z** | "0.0000" | |
| 43 | var2nm | normalVariant.numberOfMismatches | **P** | "0.0" | >0→format; else→"0" |
| 44 | shift3 | endVariant.shift3 | raw | | 0 if null |
| 45 | msi | endVariant.msi | **Z** | "0.000" | |
| 46 | msint | endVariant.msint | raw | | |
| 47 | leftSequence | endVariant.leftseq | raw | | empty→"0" |
| 48 | rightSequence | endVariant.rightseq | raw | | empty→"0" |
| 49 | region | `chr:start-end` | raw | | |
| 50 | varLabel | varLabel param | raw | | "StrongSomatic", "Germline", etc. |
| 51 | varType | beginVariant.vartype | raw | | "" if null |
| 52 | var1duprate | tumorVariant.duprate | **Z** | "0.0" | |
| 53 | var1sv | sv1 param | raw | | ""→"0" |
| 54 | var2duprate | normalVariant.duprate | **Z** | "0.0" | |
| 55 | var2sv | sv2 param | raw | | ""→"0" |

---

### Somatic Mode — 61 Columns (with --fisher)

Same base structure but all doubles use **R** (getRoundedValueToPrint) and 6 Fisher columns added:

| Col# | Name | Format | Pattern | Notes |
|------|------|--------|---------|-------|
| 1-7 | (same) | | | |
| 8-14 | tumor stats | (same) | | |
| 15 | var1frequency | **R** | "0.0000" | |
| 16 | var1strandBiasFlag | raw | | |
| 17 | var1meanPosition | **R** | "0.0" | |
| 18 | var1isAtLeastAt2Position | raw | | |
| 19 | var1meanQuality | **R** | "0.0" | |
| 20 | var1hasAtLeast2DiffQualities | raw | | |
| 21 | var1meanMappingQuality | **R** | "0.0" | |
| 22 | var1highQualityToLowQualityRatio | **R** | "0.000" | |
| 23 | var1highQualityReadsFrequency | **R** | "0.0000" | |
| 24 | var1extraFrequency | **R** | "0.0000" | |
| 25 | var1nm | **R** | "0.0" | Clamped to 0 if ≤0, then getRoundedValueToPrint |
| **26** | **pvalue1** | **R** | "0.00000" | **NEW** tumor strand bias Fisher p-value |
| **27** | **oddratio1** | raw String | | **NEW** can be "Inf" |
| 28-38 | normal stats | (same as tumor) | | all **R** format |
| 39 | var2nm | **R** | "0.0" | Clamped to 0 if ≤0 |
| **40** | **pvalue2** | **R** | "0.00000" | **NEW** normal strand bias Fisher p-value |
| **41** | **oddratio2** | raw String | | **NEW** |
| 42 | shift3 | raw | | |
| 43 | msi | **S** | "0.000" | `==0 ? "0" : DecimalFormat("0.000")` — NOT getRoundedValueToPrint! |
| 44 | msint | raw | | |
| 45 | leftSequence | raw | | |
| 46 | rightSequence | raw | | |
| 47 | region | raw | | |
| 48 | varLabel | raw | | |
| 49 | varType | raw | | |
| 50 | var1duprate | **R** | "0.0" | (not "0.00" like Simple fisher) |
| 51 | var1sv | raw | | |
| 52 | var2duprate | **R** | "0.0" | |
| 53 | var2sv | raw | | |
| **54** | **pvalue** | **R** | "0.00000" | **NEW** somatic combined Fisher p-value |
| **55** | **oddratio** | raw String | | **NEW** |

---

### Amplicon Mode — 38 Columns (no --fisher)

| Col# | Name | Format | Pattern | Notes |
|------|------|--------|---------|-------|
| 1-7 | (same as Simple) | | | |
| 8-13 | coverage stats | raw int | | |
| 14 | genotype | raw | | null→"0" |
| 15 | frequency | **Z** | "0.0000" | |
| 16 | bias | raw | | default "0;0" |
| 17 | pmean | **Z** | "0.0" | |
| 18 | pstd | raw int | | bool→1/0 |
| 19 | qual | **Z** | "0.0" | |
| 20 | qstd | raw int | | bool→1/0 |
| 21 | mapq | **Z** | "0.0" | |
| 22 | qratio | **Z** | "0.000" | |
| 23 | hifreq | **Z** | "0.0000" | |
| 24 | extrafreq | **Z** | "0.0000" | |
| 25 | shift3 | raw int | | |
| 26 | msi | **Z** | "0.000" | |
| 27 | msint | raw int | | |
| 28 | nm | **P** | "0.0" | >0→format; else→"0" |
| 29 | hicnt | raw int | | |
| 30 | hicov | raw int | | |
| 31 | leftSequence | raw | | empty→"0" |
| 32 | rightSequence | raw | | empty→"0" |
| 33 | region | raw | | **DIFFERENT**: uses goodVariants region or chr:pos-pos |
| 34 | varType | raw | | |
| 35 | **goodVariantsCount** | raw int | | **AMPLICON-SPECIFIC** |
| 36 | **totalVariantsCount** | raw int | | gvscnt + badVariants.size() |
| 37 | **noCoverage** | raw int | | **AMPLICON-SPECIFIC** |
| 38 | **ampliconFlag** | raw int | | flag ? 1 : 0 |

**No duprate, no sv, no crispr in Amplicon mode.**

---

### Amplicon Mode — 40 Columns (with --fisher)

| Col# | Name | Format | Pattern | Notes |
|------|------|--------|---------|-------|
| 1-14 | (same) | | | |
| 15 | frequency | **R** | "0.0000" | |
| 16 | bias | raw | | |
| 17 | pmean | **R** | "0.0" | |
| 18 | pstd | raw | | |
| 19 | qual | **R** | "0.0" | |
| 20 | qstd | raw | | |
| **21** | **pvalue** | **R** | "0.00000" | **NEW** |
| **22** | **oddratio** | raw String | | **NEW** |
| 23 | mapq | **R** | "0.0" | |
| 24 | qratio | **R** | "0.000" | |
| 25 | hifreq | **S** | "0.0000" | `==0 ? "0" : DecimalFormat` (NOT getRoundedValueToPrint) |
| 26 | extrafreq | **R** | "0.0000" | |
| 27 | shift3 | raw | | |
| 28 | msi | **R** | "0.000" | |
| 29 | msint | raw | | |
| 30 | nm | **S** | "0.0" | Clamped, special pre-format |
| 31-32 | hicnt, hicov | raw | | |
| 33-34 | leftseq, rightseq | raw | | |
| 35 | region | raw | | |
| 36 | varType | raw | | |
| 37 | goodVariantsCount | raw | | |
| 38 | totalVariantsCount | raw | | |
| 39 | noCoverage | raw | | |
| 40 | ampliconFlag | raw | | |

---

## Cross-Module Dependencies

### Outbound Calls (from printers)

| Caller | Called Class | Method | Purpose |
|--------|------------|--------|---------|
| All OutputVariant ctors | `GlobalReadOnlyScope` | `instance().sample`, `instance().conf` | Config access |
| SimpleOutputVariant ctor | `FisherExact` | constructor, `getPValue()`, `getOddRatio()` | Strand bias Fisher test |
| SomaticOutputVariant | `FisherExact` | constructor, `getPValue()`, `getPValueGreater()`, `getPValueLess()`, `getOddRatio()` | 3 Fisher tests |
| AmpliconOutputVariant | `FisherExact` | constructor, `getPValue()`, `getOddRatio()` | Strand bias Fisher test |
| All toString() | `Utils` | `join()`, `getRoundedValueToPrint()` | Formatting |
| VariantPrinter | `PrinterType` | enum values | Factory dispatch |

### Inbound Callers (into printers)

| Caller Class | Method | Creates | When |
|-------------|--------|---------|------|
| `SimplePostProcessModule.accept()` | Per variant per position | `SimpleOutputVariant` | Simple mode |
| `SomaticPostProcessModule.callingForOneSample()` | Per variant | `SomaticOutputVariant` | BAM1 or BAM2 only |
| `SomaticPostProcessModule.printVariationsFromFirstSample()` | Per variant | `SomaticOutputVariant` | Both BAMs, variants from BAM1 |
| `SomaticPostProcessModule.printVariationsFromSecondSample()` | Per variant | `SomaticOutputVariant` | Both BAMs, BAM1 has only ref |
| `AmpliconPostProcessModule.process()` | Per variant per position | `AmpliconOutputVariant` | Amplicon mode |

---

## Data Structures Read/Written

### Variant Fields Accessed

| Variant Field | Java Type | Default | Accessed By | Non-fisher Format | Fisher Format |
|--------------|-----------|---------|-------------|-------------------|--------------|
| `startPosition` | int | 0 | all | raw | raw |
| `endPosition` | int | 0 | all | raw | raw |
| `refallele` | String | "" | all | raw | raw |
| `varallele` | String | "" | all | raw | raw |
| `totalPosCoverage` | int | 0 | all | raw | raw |
| `positionCoverage` | int | 0 | all | raw | raw |
| `refForwardCoverage` | int | 0 | all | raw | raw |
| `refReverseCoverage` | int | 0 | all | raw | raw |
| `varsCountOnForward` | int | 0 | all | raw | raw |
| `varsCountOnReverse` | int | 0 | all | raw | raw |
| `genotype` | String | null | all | null→"0" | null→"0" |
| `frequency` | double | 0.0 | all | Z "0.0000" | R "0.0000" |
| `strandBiasFlag` | String | "0" | all | raw (no null check in Simple!) | raw |
| `meanPosition` | double | 0.0 | all | Z "0.0" | R "0.0" |
| `isAtLeastAt2Positions` | boolean | false | all | bool→1/0 | bool→1/0 |
| `meanQuality` | double | 0.0 | all | Z "0.0" | R "0.0" |
| `hasAtLeast2DiffQualities` | boolean | false | all | bool→1/0 | bool→1/0 |
| `meanMappingQuality` | double | 0.0 | all | Z "0.0" | R "0.0" |
| `highQualityToLowQualityRatio` | double | 0.0 | all | Z "0.000" | R "0.000" |
| `highQualityReadsFrequency` | double | 0.0 | all | Z "0.0000" | S "0.0000" |
| `extraFrequency` | double | 0.0 | all | Z "0.0000" | R "0.0000" |
| `shift3` | int | 0 | all | raw | raw |
| `msi` | double | 0.0 | all | Z "0.000" | R or S (varies) |
| `msint` | int | 0 | all | raw | raw |
| `numberOfMismatches` | double | 0.0 | all | P "0.0" | S "0.0" |
| `hicnt` | int | 0 | Simple, Amplicon | raw | raw |
| `hicov` | int | 0 | Simple, Amplicon | raw | raw |
| `leftseq` | String | "" | all | empty→"0" | empty→"0" |
| `rightseq` | String | "" | all | empty→"0" | empty→"0" |
| `vartype` | String | "" | all | raw | raw |
| `duprate` | double | 0.0 | Simple, Somatic | Z "0.0" | R "0.00"(Simple) / R "0.0"(Somatic) |
| `crispr` | int | 0 | Simple only | raw | raw |
| `DEBUG` | String | "" | all (debug mode) | raw | raw |

### Config Fields Accessed

| Config Field | Where Used | Effect |
|-------------|-----------|--------|
| `instance().sample` | All constructors | Column 1 |
| `conf.fisher` | All `toString()` | 36/38, 55/61, or 38/40 column dispatch |
| `conf.crisprCuttingSite` | `SimpleOutputVariant.toString()` | Appends crispr column if ≠ 0 |
| `conf.debug` | All `toString()` + Amplicon ctor | Appends DEBUG column; builds debug amp string |

---

## Known Parity Traps

### Trap 1: Dual Formatting Scheme — Zero-Check vs getRoundedValueToPrint

The **same field** is formatted with **two different algorithms** depending on `--fisher`:

- **Non-fisher (Z)**: `value == 0 ? 0 : new DecimalFormat(pattern).format(value)`
  - Zero → "0"
  - Non-zero → DecimalFormat with **trailing zeros PRESERVED** (e.g., "12.3400" stays as "12.3400")

- **Fisher (R)**: `getRoundedValueToPrint(pattern, value)`
  - Whole number → integer string without decimals
  - Non-whole → DecimalFormat then **trailing zeros STRIPPED** (e.g., "12.3400" → "12.34")

**Impact**: A value of exactly 1.0 outputs "1.0" in non-fisher (with "0.0" pattern) but "1" in fisher. A value of 0.5000 outputs "0.5000" in non-fisher (with "0.0000") but "0.5" in fisher.

### Trap 2: hifreq Uses Special Formatting in Fisher Mode

In fisher-mode (38-col Simple, 40-col Amplicon), hifreq is NOT formatted with `getRoundedValueToPrint`:
```java
String hifreq_f = hifreq == 0 ? "0" : new DecimalFormat("0.0000").format(hifreq);
```
This **preserves trailing zeros** (unlike getRoundedValueToPrint used for neighboring fields).

### Trap 3: nm Uses Special Formatting in Fisher Mode

Fisher mode:
```java
nm = nm > 0 ? nm : 0;  // clamp first
String nm_f = nm == 0 ? "0" : new DecimalFormat("0.0").format(nm);
```
Non-fisher mode:
```java
nm > 0 ? new DecimalFormat("0.0").format(nm) : 0
```
Both functionally equivalent for positive nm, but the pre-clamp step in fisher must be replicated.

### Trap 4: Somatic 61-col msi Uses Special Formatting

In 61-column somatic (fisher), msi is pre-formatted:
```java
String msi_f = msi == 0 ? "0" : new DecimalFormat("0.000").format(msi);
```
NOT `getRoundedValueToPrint("0.000", msi)`. A whole-number non-zero msi (e.g., 2.0) → "2.000" (preserved) instead of "2" (stripped).

### Trap 5: duprate Pattern Differs Between Modes/Fisher

| Mode | Non-fisher | Fisher |
|------|-----------|--------|
| Simple | Z "0.0" | R **"0.00"** |
| Somatic | Z "0.0" | R "0.0" |
| Amplicon | N/A | N/A |

Simple fisher uses "0.00" (2 decimals) while Somatic uses "0.0" (1 decimal).

### Trap 6: Null Variant Produces Different Defaults Across Modes

| Field | Simple (variant==null) | Somatic (tumorVariant==null) |
|-------|----------------------|---------------------------|
| genotype | **""** (empty) | "0" |
| bias/strandBiasFlag | "0;0" | "0" |
| leftSequence | **""** (base class default!) | "0" (if endVariant sets it) |
| rightSequence | **""** (base class default!) | "0" (if endVariant sets it) |

**CRITICAL**: When variant==null in Simple, `leftSequence` and `rightSequence` stay at base class default `""` because the `empty→"0"` mapping only fires in the `variant != null` branch.

### Trap 7: SomaticOutputVariant Constructor Parameter Reuse

The 4 Variant params (`beginVariant`, `endVariant`, `tumorVariant`, `normalVariant`) serve different roles. Callers pass the SAME object for multiple slots:
```java
// StrongSomatic: vref used 3 times
new SomaticOutputVariant(vref, vref, vref, varForPrint, ...)
// Deletion (isFirstCover=true): variant used for begin, end, and normal (not tumor!)  
new SomaticOutputVariant(variant, variant, null, variant, ...)
```

### Trap 8: Amplicon region String Differs from Simple/Somatic

- **Simple/Somatic**: `region.chr + ":" + region.start + "-" + region.end` (BED region)
- **Amplicon** (variant != null, goodVariants non-empty): `goodVariants.get(0)._2` (first good variant's amplicon region)
- **Amplicon** (variant == null or goodVariants empty): `chr + ":" + position + "-" + position`

### Trap 9: Amplicon totalVariantsCount Calculation

`totalVariantsCount = gvscnt + badVariants.size()` where `gvscnt` is the **parameter** (count of amplicons with this variant), NOT `goodVariants.size()`.

### Trap 10: strandBiasFlag Null Check Inconsistency

- **Simple**: `bias = variant.strandBiasFlag` — NO null check
- **Somatic**: `var1strandBiasFlag = tumorVariant.strandBiasFlag == null ? "0" : ...` — null-checked
- In practice safe (Variant initializes it to "0"), but Rust must match the behavior

### Trap 11: FisherExact.getOddRatio() Returns Special Strings

Returns:
- `"Inf"` — when odds ratio is infinite
- `"0"` — when exactly 0 (formatted via `DecimalFormat("0")`)
- Decimal via `round_as_r()` + `String.valueOf()` — Java's `Double.toString()` formatting

### Trap 12: Somatic Fisher Combined P-value Direction

```java
if (pvalue_less < pvalue_greater) pvalue = pvalue_less;
else pvalue = pvalue_greater;
```
When equal, `pvalue_greater` is chosen (else branch). Must match exactly.

### Trap 13: Amplicon Has No duprate, sv, or crispr Columns

Amplicon replaces Simple's last 2 columns (duprate, sv) with 4 amplicon-specific columns. No crispr column regardless of config.

### Trap 14: println() Adds Newline

`VariantPrinter.print()` uses `out.println()` which adds `\n` on Linux. The `toString()` does NOT include trailing newline.

### Trap 15: DecimalFormat Uses HALF_EVEN Rounding

Java's `DecimalFormat` defaults to `RoundingMode.HALF_EVEN` (banker's rounding). Examples:
- `"0.0".format(0.25)` → "0.2" (rounds to even)
- `"0.0".format(0.35)` → "0.4" (rounds to even)
- `"0.0".format(0.45)` → "0.4" (rounds to even)

### Trap 16: Debug Mode Amplicon Output

Debug appends per-amplicon entries separated by `\t`:
```
{variant.DEBUG}\tGood0 {space-delimited-stats} {region}\tGood1 ...
```
Stats within each entry are SPACE-separated (debugAmpVariant uses `" "`). The pmean/qual/mapq/qratio fields are ALWAYS DecimalFormat-formatted (even when 0 → "0.0"), unlike main columns.

### Trap 17: Somatic varLabel Values

Column 50/48 takes specific strings: `"StrongSomatic"`, `"LikelySomatic"`, `"Germline"`, `"StrongLOH"`, `"LikelyLOH"`, `"AFDiff"`, `"SampleSpecific"`, `"Deletion"`, `""` (empty).

### Trap 18: Null-Variant Row Printing in SimplePostProcessModule

Two paths print null-variant rows:
1. `variants.isEmpty() && doPileup && vref==null` → `SimpleOutputVariant(null, region, sv, position)` → prints skeleton row, continues
2. `vref.startPosition != position && doPileup && vvar.size()==1 && refVar==null` → prints skeleton row, THEN creates `new Variant()` for subsequent use
