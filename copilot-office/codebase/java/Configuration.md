# Configuration

**Source**: `Configuration.java`, `CmdParser.java`, `VarDictLauncher.java`, `data/scopedata/GlobalReadOnlyScope.java`
**LOC**: Configuration.java 374 + CmdParser.java ~640 + VarDictLauncher.java ~290 + GlobalReadOnlyScope.java 76
**Rust counterpart**: `src/conf.rs` (Configuration), CLI handled via `clap` in `src/bin/`
**Status**: complete

## Overview

Configuration is a plain mutable data class (all public fields, no encapsulation) holding every command-line parameter and its default value. It has no behavior beyond 4 boolean query methods and an inner class `BamNames`. `CmdParser` is the sole writer — it creates and populates the object using Apache Commons CLI. `VarDictLauncher` takes the populated `Configuration`, initializes resources (BAM header, sample names, BED regions, adaptor seeds), wraps everything in the `GlobalReadOnlyScope` singleton, then starts the pipeline. Every pipeline module reads configuration via `GlobalReadOnlyScope.instance().conf.*`.

## Method Inventory

| Method | Class | Lines | Analyzed? | Summary |
|--------|-------|-------|-----------|---------|
| `parseParams()` | CmdParser | L25-L54 | yes | Top-level CLI entry: build options, parse, handle missing required |
| `parseCmd()` | CmdParser | L62-L182 | yes | Maps CommandLine values to Configuration fields — sole writer |
| `buildOptions()` | CmdParser | L218-L495 | yes | Registers all CLI options with types and descriptions |
| `setFastaFile()` | CmdParser | L195-L216 | yes | Resolve `-G` value or alias (hg19/hg38/mm10) |
| `getIntValue()` | CmdParser | L600-L603 | yes | Parse int option with default |
| `getColumnValue()` | CmdParser | L605-L608 | yes | Parse column option: subtracts 1 (1-based CLI → 0-based internal) |
| `getDoubleValue()` | CmdParser | L610-L613 | yes | Parse double option with default |
| `readThreadsCount()` | CmdParser | L620-L631 | yes | Thread count: auto-detect CPU count if `-th` without value |
| `start()` | VarDictLauncher | L44-L65 | yes | Init resources, pick mode, run parallel/serial |
| `initResources()` | VarDictLauncher | L73-L133 | yes | Read BAM header, sample names, BED, adaptors → GlobalReadOnlyScope.init() |
| `readBedFile()` | VarDictLauncher | L141-L186 | yes | Read BED, detect amplicon mode from 8-col BED structure |
| `readChr()` | VarDictLauncher | L194-L206 | yes | Read chromosome lengths from BAM header |
| `getSampleNames()` | VarDictLauncher | L213-L245 | yes | Sample name for simple/amplicon mode: -N or regex match |
| `getSampleNamesSomatic()` | VarDictLauncher | L254-L281 | yes | Sample names for somatic mode: per-BAM regex or `-N` split by `\|` |
| `isColumnForChromosomeSet()` | Configuration | L327 | yes | `columnForChromosome >= 0` |
| `isDownsampling()` | Configuration | L331 | yes | `downsampling != null` |
| `hasMappingQuality()` | Configuration | L335 | yes | `mappingQuality != null` |
| `isZeroBasedDefined()` | Configuration | L339 | yes | `zeroBased != null` |
| `init()` | GlobalReadOnlyScope | L23-L30 | yes | Synchronized singleton init, throws if already initialized |
| `setMode()` | GlobalReadOnlyScope | L38-L42 | yes | Sets run mode once |
| `clear()` | GlobalReadOnlyScope | L48-L51 | yes | TEST ONLY — resets singleton |

## Instance Fields (59 Fields)

### CLI Parameter Fields (43 fields)

| # | Java Field | Type | Default (Java) | Default (CmdParser) | CLI Flag | Parsing Method | Modules Reading It |
|---|-----------|------|-----------------|---------------------|----------|----------------|-------------------|
| 1 | `printHeader` | `boolean` | `false` | `false` | `-h` | `cmd.hasOption('h')` | VarDictLauncher (header output) |
| 2 | `delimiter` | `String` | `null` | `"\t"` | `-d` | `cmd.getOptionValue("d", "\t")` | RegionBuilder (BED line splitting) |
| 3 | `bed` | `String` | `null` | `null` | positional arg | `args[0]` from `cmd.getArgs()` | VarDictLauncher (readBedFile) |
| 4 | `numberNucleotideToExtend` | `int` | `0` | `0` | `-x` | `getIntValue(cmd, "x", 0)` | RegionBuilder (region extension) |
| 5 | `zeroBased` | `Boolean` (boxed) | `null` | tristate | `-z` | `1 == getIntValue(cmd, "z", 1)` only if `-z` present | RegionBuilder, VarDictLauncher; `isZeroBasedDefined()` |
| 6 | `ampliconBasedCalling` | `String` | `null` | `null` | `-a` | `cmd.getOptionValue("a")` | VarDictLauncher (amplicon detect), GlobalReadOnlyScope |
| 7 | `columnForChromosome` | `int` | `-1` | from `-c` or `-1` | `-c` | `getColumnValue(cmd, "c", -1)` | `isColumnForChromosomeSet()`, RegionBuilder |
| 8 | `bedRowFormat` | `BedRowFormat` | `null` | constructed | `-c -S -E -s -e -g` | Composite from `getColumnValue()` calls | RegionBuilder |
| 9 | `sampleNameRegexp` | `String` | `null` | `null` | `-n` | Strips leading/trailing `/` | VarDictLauncher (sample name extraction) |
| 10 | `sampleName` | `String` | `null` | `null` | `-N` | `cmd.getOptionValue("N")` | VarDictLauncher (direct sample name) |
| 11 | `fasta` | `String` | `null` | HG19 if unset | `-G` | `setFastaFile(cmd)` — resolves aliases | ReferenceResource |
| 12 | `bam` | `BamNames` | `null` | required | `-b` (required) | `new BamNames(cmd.getOptionValue("b"))` | VarDictLauncher (BAM file access) |
| 13 | `downsampling` | `Double` (boxed) | `null` | `null` | `-Z` | Only set if `-Z` present; `getDoubleValue(cmd, "Z", 0)` | SAMFileParser, CigarParser; `isDownsampling()` |
| 14 | `chromosomeNameIsNumber` | `boolean` | `false` | `false` | `-C` | `cmd.hasOption("C")` | RegionBuilder (deprecated) |
| 15 | `mappingQuality` | `Integer` (boxed) | `null` | `null` | `-Q` | Only if present; `cmd.getParsedOptionValue("Q")` | RecordPreprocessor, CigarParser; `hasMappingQuality()` |
| 16 | `removeDuplicatedReads` | `boolean` | `false` | `false` | `-t` | `cmd.hasOption("t")` | RecordPreprocessor (flag-based filtering) |
| 17 | `mismatch` | `int` | `0` | `8` | `-m` | `getIntValue(cmd, "m", 8)` | RecordPreprocessor (NM-based filtering) |
| 18 | `y` | `boolean` | `false` | `false` | `-y` | `cmd.hasOption("y")` | Various (verbose logging) |
| 19 | `goodq` | `double` | `0.0` | `22.5` | `-q` | `getDoubleValue(cmd, "q", 22.5)` | CigarParser, ToVarsBuilder, VariationRealigner |
| 20 | `vext` | `int` | `2` | `2` | `-X` | `getIntValue(cmd, "X", 2)` | CigarParser (extension around indels) |
| 21 | `trimBasesAfter` | `int` | `0` | `0` | `-T` | `getIntValue(cmd, "T", 0)` | RecordPreprocessor (read trimming) |
| 22 | `performLocalRealignment` | `boolean` | `false` | `true` (!) | `-k` | `1 == getIntValue(cmd, "k", 1)` — default `1` → `true` | VariationRealigner |
| 23 | `indelsize` | `int` | `50` | `50` | `-I` | `getIntValue(cmd, "I", 50)` | ToVarsBuilder (indel size threshold) |
| 24 | `bias` | `double` | `0.05` | N/A | (none) | Hardcoded, no CLI flag | ToVarsBuilder (strand bias cutoff) |
| 25 | `minBiasReads` | `int` | `2` | `2` | `-B` | `getIntValue(cmd, "B", 2)` | ToVarsBuilder (min reads for bias test) |
| 26 | `minr` | `int` | `2` | `2` (→ `0` if `-p`) | `-r` | `getIntValue(cmd, "r", 2)` | ToVarsBuilder (min variant reads) |
| 27 | `debug` | `boolean` | `false` | `false` | `-D` | `cmd.hasOption("D")` | OutputVariant (extra debug columns) |
| 28 | `freq` | `double` | `0.01` | `0.01` (→ `-1` if `-p`) | `-f` | `getDoubleValue(cmd, "f", 0.01d)` | ToVarsBuilder (allele frequency filter) |
| 29 | `moveIndelsTo3` | `boolean` | `false` | `false` | `-3` | `cmd.hasOption("3")` | VariationRealigner |
| 30 | `samfilter` | `String` | `"0x504"` | `"0x504"` | `-F` | `cmd.getOptionValue("F", "0x504")` | RecordPreprocessor (SAM flag filtering) |
| 31 | `regionOfInterest` | `String` | `null` | `null` | `-R` | `cmd.getOptionValue("R")` | VarDictLauncher (region vs BED selection) |
| 32 | `readPosFilter` | `int` | `5` | `5` | `-P` | `getIntValue(cmd, "P", 5)` | ToVarsBuilder (read position filter) |
| 33 | `qratio` | `double` | `1.5` | `1.5` | `-o` | `getDoubleValue(cmd, "o", 1.5d)` | ToVarsBuilder |
| 34 | `mapq` | `double` | `0` | `0` | `-O` | `getDoubleValue(cmd, "O", 0)` | ToVarsBuilder (mean mapping quality filter) |
| 35 | `doPileup` | `boolean` | `false` | `false` | `-p` | Side-effects: sets `freq = -1`, `minr = 0` | ToVarsBuilder, OutputVariant |
| 36 | `lofreq` | `double` | `0.05` | `0.05` | `-V` | `getDoubleValue(cmd, "V", 0.05d)` | Somatic mode (normal freq threshold) |
| 37 | `minmatch` | `int` | `0` | `0` | `-M` | `cmd.getParsedOptionValue("M")` only if present | RecordPreprocessor (min aligned bases) |
| 38 | `outputSplicing` | `boolean` | `false` | `false` | `-i` | `cmd.hasOption('i')` | VarDictLauncher (mode selection), SplicingMode |
| 39 | `validationStringency` | `ValidationStringency` | `LENIENT` | `LENIENT` | `-VS` | `ValidationStringency.valueOf(...)` | BAM reader factory |
| 40 | `includeNInTotalDepth` | `boolean` | `false` | `false` | `-K` | `cmd.hasOption("K")` | CigarParser (N base depth counting) |
| 41 | `uniqueModeAlignmentEnabled` | `boolean` | `false` | `false` | `-u` | `cmd.hasOption("u")` | RecordPreprocessor |
| 42 | `uniqueModeSecondInPairEnabled` | `boolean` | `false` | `false` | `-UN` | `cmd.hasOption("UN")` | RecordPreprocessor |
| 43 | `threads` | `int` | `0` | `max(readThreadsCount, 1)` | `-th` | Auto-detects CPU count if no value | VarDictLauncher (parallel vs serial) |

### SV/Insert-Size Fields (4 CLI fields)

| # | Java Field | Type | Default | CLI Flag | Parsing Method | Modules Reading It |
|---|-----------|------|---------|----------|----------------|-------------------|
| 44 | `INSSIZE` | `int` | `300` | `-w` | `getIntValue(cmd, "w", 300)` | StructuralVariantsProcessor (mean insert size) |
| 45 | `INSSTD` | `int` | `100` | `-W` | `getIntValue(cmd, "W", 100)` | StructuralVariantsProcessor (insert size std dev) |
| 46 | `INSSTDAMT` | `int` | `4` | `-A` | `getIntValue(cmd, "A", 4)` | StructuralVariantsProcessor (std dev multiplier) |
| 47 | `SVMINLEN` | `int` | `1000` | `-L` | `getIntValue(cmd, "L", 1000)` | StructuralVariantsProcessor, ToVarsBuilder |

### Control/Filter Flags (4 CLI fields)

| # | Java Field | Type | Default | CLI Flag | Parsing Method | Modules Reading It |
|---|-----------|------|---------|----------|----------------|-------------------|
| 48 | `chimeric` | `boolean` | `false` | `--chimeric` | `cmd.hasOption("chimeric")` | RecordPreprocessor (chimeric read filter) |
| 49 | `disableSV` | `boolean` | `false` | `-U` / `--nosv` | `cmd.hasOption("U")` | StructuralVariantsProcessor |
| 50 | `deleteDuplicateVariants` | `boolean` | `false` | `--deldupvar` | `cmd.hasOption("deldupvar")` | Post-processing |
| 51 | `fisher` | `boolean` | `false` | `--fisher` | `cmd.hasOption("fisher")` | OutputVariant (Fisher exact test) |

### Extension/Printer Fields (2 CLI fields)

| # | Java Field | Type | Default | CLI Flag | Parsing Method | Modules Reading It |
|---|-----------|------|---------|----------|----------------|-------------------|
| 52 | `referenceExtension` | `int` | `1200` | `-Y` | `getIntValue(cmd, "Y", 1200)` | ReferenceResource (FASTA fetch window) |
| 53 | `printerType` | `PrinterType` | `PrinterType.OUT` | `--DP` | Switch: `"OUT"` or `"ERR"` | GlobalReadOnlyScope → all printers |

### Runtime/Exception Fields (1 non-CLI)

| # | Java Field | Type | Default | CLI Flag | Modules Reading It |
|---|-----------|------|---------|----------|-------------------|
| 54 | `exceptionCounter` | `AtomicInteger` | `0` | (none) | Pipeline (exception counting) |

### Adaptor/CRISPR Fields (3 CLI fields)

| # | Java Field | Type | Default | CLI Flag | Parsing Method | Modules Reading It |
|---|-----------|------|---------|----------|----------------|-------------------|
| 55 | `adaptor` | `List<String>` | `new ArrayList<>()` | `--adaptor` | Split by `,` | VarDictLauncher (adaptor seed build) |
| 56 | `crisprFilteringBp` | `int` | `0` | `-j` | `getIntValue(cmd, "j", 0)` | CigarParser (CRISPR overlap filter) |
| 57 | `crisprCuttingSite` | `int` | `0` | `-J` / `--crispr` | `getIntValue(cmd, "J", 0)` | VariationRealigner (CRISPR position adjust) |

### MSI Frequency Fields (2 CLI fields)

| # | Java Field | Type | Default | CLI Flag | Parsing Method | Modules Reading It |
|---|-----------|------|---------|----------|----------------|-------------------|
| 58 | `monomerMsiFrequency` | `double` | `0.25` | `-mfreq` | `getDoubleValue(cmd, "mfreq", 0.25d)` | ToVarsBuilder (monomer MSI frequency) |
| 59 | `nonMonomerMsiFrequency` | `double` | `0.1` | `-nmfreq` | `getDoubleValue(cmd, "nmfreq", 0.1d)` | ToVarsBuilder (non-monomer MSI frequency) |

## Static Constants (16 Constants)

| # | Name | Type | Value | Usage |
|---|------|------|-------|-------|
| 1 | `HG19` | `String` | `/ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa` | CmdParser default reference |
| 2 | `HG38` | `String` | `/ngs/reference_data/genomes/Hsapiens/hg38/seq/hg38.fa` | CmdParser alias |
| 3 | `MM10` | `String` | `/ngs/reference_data/genomes/Mmusculus/mm10/seq/mm10.fa` | CmdParser alias |
| 4 | `LOWQUAL` | `int` | `10` | CigarParser (soft-clip quality threshold) |
| 5 | `SEED_1` | `int` | `17` | VariationRealigner (large seed) |
| 6 | `SEED_2` | `int` | `12` | VariationRealigner (small seed) |
| 7 | `ADSEED` | `int` | `6` | VarDictLauncher (adaptor seed length) |
| 8 | `MINSVCDIST` | `double` | `1.5` | StructuralVariantsProcessor (min SV cluster distance) |
| 9 | `MINMAPBASE` | `int` | `15` | StructuralVariantsProcessor (min mapping pos for SV) |
| 10 | `MINSVPOS` | `int` | `25` | StructuralVariantsProcessor (inter-chr SV min distance) |
| 11 | `SVMAXLEN` | `int` | `150000` | VariationRealigner (max SV in realignment step) |
| 12 | `SVFLANK` | `int` | `50` | StructuralVariantsProcessor (SV flanking seq length) |
| 13 | `DISCPAIRQUAL` | `int` | `35` | StructuralVariantsProcessor (discordant pair mapq) |
| 14 | `EXTENSION` | `int` | `5000` | ReferenceResource (downstream extension) |
| 15 | `DEFAULT_AMPLICON_PARAMETERS` | `String` | `"10:0.95"` | VarDictLauncher (amplicon default) |
| 16 | `MAX_EXCEPTION_COUNT` | `int` | `10` | Pipeline (exception threshold) — not `static final`, just `static` |

## Inner Class: BamNames

```java
public static class BamNames {
    private final String[] bamNames;  // split by "|" (tumor|normal)
    private final String[] bams;      // split bamNames[0] by ":"
    private final String bamRaw;      // original value

    public String getBam1()   → bamNames[0]          // tumor BAM (or only BAM)
    public String getBam2()   → bamNames[1] or null   // normal BAM (somatic mode)
    public String getBamX()   → bams[0]               // first BAM in colon-separated list
    public boolean hasBam2()  → bamNames.length > 1
    public String getBamRaw() → bamRaw                 // original string
}
```

Parsing: `-b "tumor.bam|normal.bam"` → `bamNames = ["tumor.bam", "normal.bam"]`, `bams = ["tumor.bam"]`.

## Query Methods

| Method | Logic | Used By |
|--------|-------|---------|
| `isColumnForChromosomeSet()` | `columnForChromosome >= 0` | RegionBuilder (custom BED format detection) |
| `isDownsampling()` | `downsampling != null` | Not widely used |
| `hasMappingQuality()` | `mappingQuality != null` | RecordPreprocessor |
| `isZeroBasedDefined()` | `zeroBased != null` | RegionBuilder, VarDictLauncher |

## Method Analyses

### CmdParser.parseParams()
**Source**: CmdParser.java:L25-L54
**Purpose**: Top-level entry point for CLI parsing.

**Algorithm**:
1. Build `Options` via `buildOptions()` — registers all CLI options.
2. Parse with `BasicParser.parse(options, args)`.
3. If no options or `-H` → print help and `System.exit(0)`.
4. Call `new CmdParser().parseCmd(cmd)` (note: creates a NEW CmdParser instance).
5. On `MissingOptionException` → print missing options to stderr, show help, exit.
6. Return `Configuration`.

**Parity Note**: Only `-b` is marked `isRequired(true)`. All other options are optional.

### CmdParser.parseCmd()
**Source**: CmdParser.java:L62-L182
**Purpose**: Maps CommandLine values to Configuration fields — the sole writer.

**Algorithm** (step-by-step):
1. Create `new Configuration()`.
2. Positional args: `cmd.getArgs()` → `config.bed = args[0]` if present.
3. Boolean flags: `printHeader`, `chromosomeNameIsNumber`, `debug`, `removeDuplicatedReads`, `moveIndelsTo3`.
4. `samfilter`: defaults to `"0x504"` — **stored as String in Java**.
5. `zeroBased`: only set if `-z` flag present; `1 == getIntValue(cmd, "z", 1)` → defaults to `true` when flag given without value.
6. `performLocalRealignment`: `1 == getIntValue(cmd, "k", 1)` → **defaults to `true`** when `-k` not given (default value is 1).
7. `fasta`: resolved through `setFastaFile()` which maps `"hg19"` → absolute path. Falls back to HG19 if null.
8. BED row format columns: `-c`, `-S`, `-E`, `-s`, `-e`, `-g` → all parsed via `getColumnValue()` which **subtracts 1** (1-based CLI → 0-based internal).
   - Special logic: if `-S` given but not `-s`, then `s_col = S_col`. Same for `-E`/`-e`.
9. Numeric fields: each parsed with explicit defaults.
10. **`-p` (pileup) side-effects**: `doPileup = true`, `freq = -1`, `minr = 0`. This OVERRIDES prior `-f` and `-r` values.
11. `minmatch`: only set if `-M` present.
12. SV parameters: `INSSIZE`, `INSSTD`, `INSSTDAMT`, `SVMINLEN`.
13. `threads`: `Math.max(readThreadsCount(cmd), 1)`.
14. `adaptor`: split by `,`.
15. `printerType`: `"OUT"` or `"ERR"`.
16. CRISPR and MSI fields last.

### CmdParser.buildOptions()
**Source**: CmdParser.java:L218-L495
**Purpose**: Define all CLI options with their types, descriptions, and constraints.

Key patterns:
- `-z` and `-k` use `hasOptionalArgs(1)` — they can appear without a value.
- `-Q`, `-M` use `.getParsedOptionValue()` which returns `Number`, not raw string.
- `-th` uses `hasOptionalArg()` — if no value, system detects CPU count.
- Type annotations are `Number.class` for numerics, `String.class` for strings.

### CmdParser.setFastaFile()
**Source**: CmdParser.java:L195-L216
**Purpose**: Resolve `-G` value — short aliases `hg19`/`hg38`/`mm10` → absolute paths.

If `-G` is null → print warning to stderr and default to HG19 path.

### CmdParser.getColumnValue()
**Source**: CmdParser.java:L605-L608
**Critical**: Subtracts 1 from CLI value. `((Number)value).intValue() - 1`. This makes CLI 1-based columns into 0-based internal indices.

### CmdParser.readThreadsCount()
**Source**: CmdParser.java:L620-L631
**Purpose**: Thread count from `-th`. If value given → use it. If `-th` without value → `Runtime.getRuntime().availableProcessors()`. If `-th` not present → returns 0 (which `parseCmd` clamps to 1 via `Math.max`).

### VarDictLauncher.start()
**Source**: VarDictLauncher.java:L44-L65
**Purpose**: Init resources, pick mode, run pipeline.

**Mode selection**:
- `-i` → `SplicingMode`
- `-R` or no amplicon → `SomaticMode` (if 2 BAMs) or `SimpleMode` (if 1 BAM)
- Amplicon detected → `AmpliconMode`
- `threads == 1` → `mode.notParallel()`, else `mode.parallel()`

### VarDictLauncher.initResources()
**Source**: VarDictLauncher.java:L73-L133
**Purpose**: Initialize all resources before pipeline execution.

**Algorithm**:
1. Validate region source (`-R` or BED) exists → `RegionMissedSourceException` if not.
2. `readChr(conf.bam.getBamX())` → `Map<String, Integer>` from BAM header.
3. Get sample names: somatic vs simple logic differs.
4. Create `RegionBuilder(chrLengths, conf)`.
5. If `-R` given: `buildRegionFromConfiguration()`.
6. Else: `readBedFile(conf)` → detect amplicon, set `zeroBased`, read lines.
7. Build adaptor seed maps from `conf.adaptor`:
   - For each adaptor, extract seeds of length `ADSEED=6` at offsets 0..6.
   - Forward seed: `substr(sequence, i, ADSEED)`.
   - Reverse seed: `complement(reverse(forwardSeed))`.
   - Store `(seed, i+1)` in maps.
8. `GlobalReadOnlyScope.init(conf, chrLengths, sample, samplem, ampliconBasedCalling, adaptorForward, adaptorReverse)`.

### VarDictLauncher.readBedFile()
**Source**: VarDictLauncher.java:L141-L186
**Purpose**: Read BED file, auto-detect amplicon mode.

**Algorithm**:
1. Open BED file with BufferedReader.
2. Skip lines starting with `#`, `browser`, `track`.
3. For each non-header line:
   - If amplicon not yet detected:
     - Split by `conf.delimiter`.
     - If exactly 8 columns AND columns 7,8 match `INTEGER_ONLY` regex:
       - Parse: `startRegion = col[1]`, `endRegion = col[2]`, `startAmplicon = col[6]`, `endAmplicon = col[7]`.
       - If `startAmplicon >= startRegion && endAmplicon <= endRegion`: amplicon mode ON.
       - Set `ampliconParameters = DEFAULT_AMPLICON_PARAMETERS` ("10:0.95").
       - Set `zeroBased = true` UNLESS user explicitly defined it.
4. Add line to segraw list.
5. Return `(ampliconParameters, zeroBased, segraw)`.

### VarDictLauncher.readChr()
**Source**: VarDictLauncher.java:L194-L206
Opens BAM with `SamReaderFactory.makeDefault()`, reads `SAMSequenceDictionary`, returns `Map<String, Integer>` of chromosome name → length.

### VarDictLauncher.getSampleNames()
**Source**: VarDictLauncher.java:L213-L245
**Purpose**: Extract sample name for simple/amplicon mode.

1. If `-N` set → use directly.
2. Else: compile regex from `-n` (or default `SAMPLE_PATTERN`).
3. Match against `conf.bam.getBamRaw()`. Use `matcher.group(1)`.
4. Fallback: try `SAMPLE_PATTERN2`.
5. Return `(sample, "")`.

### VarDictLauncher.getSampleNamesSomatic()
**Source**: VarDictLauncher.java:L254-L281
**Purpose**: Extract sample names for somatic mode (two BAMs).

1. Start with `getSampleNames()` result.
2. If `-n` regex set: apply to each BAM separately.
3. Else if `-N` set: split by `|`. First → sample, second → samplem. If only one part → samplem = `"${sample}_match"`.

## GlobalReadOnlyScope — Detailed

**Source**: `data/scopedata/GlobalReadOnlyScope.java` (76 LOC)

### Singleton Pattern
```java
private volatile static GlobalReadOnlyScope instance;
private volatile static AbstractMode mode;

public static synchronized void init(...) {
    if (instance != null) throw IllegalStateException("already initialized")
    instance = new GlobalReadOnlyScope(...)
}

public static synchronized void setMode(AbstractMode runMode) {
    if (mode != null) throw IllegalStateException("already initialized")
    mode = runMode;
}

public static synchronized void clear() { // TEST ONLY
    instance = null;  mode = null;
}
```

### Fields (all `public final`)

| Field | Type | Source |
|-------|------|--------|
| `conf` | `Configuration` | Passed from init |
| `chrLengths` | `Map<String, Integer>` | From BAM header |
| `sample` | `String` | Primary sample name |
| `samplem` | `String` | Secondary sample name (somatic) |
| `ampliconBasedCalling` | `String` | Amplicon parameters or null |
| `printerTypeOut` | `PrinterType` | Copied from `conf.printerType` |
| `adaptorForward` | `Map<String, Integer>` | Forward adaptor seed map |
| `adaptorReverse` | `Map<String, Integer>` | Reverse adaptor seed map |

### Rust Equivalent
`GlobalReadOnlyScope` in Rust uses `OnceLock<GlobalReadOnlyScope>` (never re-init). Differences:
- No `sample`/`samplem` — sample names handled differently
- No `printerTypeOut` — always stdout
- `chrLengths` → `chr_lens: HashMap<String, usize>` (usize not i32)
- Added `bam_paths: Vec<String>` — not in Java

## Cross-Module Dependencies

### Configuration is Read by (Instance Fields)
- **RecordPreprocessor**: `mappingQuality`, `removeDuplicatedReads`, `mismatch`, `samfilter`, `trimBasesAfter`, `minmatch`, `uniqueMode*`, `chimeric`
- **CigarParser**: `goodq`, `vext`, `downsampling`, `includeNInTotalDepth`, `crisprFilteringBp`
- **VariationRealigner**: `performLocalRealignment`, `moveIndelsTo3`, `goodq`, `crisprCuttingSite`
- **StructuralVariantsProcessor**: `INSSIZE`, `INSSTD`, `INSSTDAMT`, `disableSV`, `SVMINLEN`
- **ToVarsBuilder**: `freq`, `minr`, `minBiasReads`, `readPosFilter`, `qratio`, `mapq`, `bias`, `indelsize`, `goodq`, `doPileup`, `lofreq`, `monomerMsiFrequency`, `nonMonomerMsiFrequency`
- **OutputVariant**: `debug`, `fisher`, `SVMINLEN`
- **RegionBuilder**: `bedRowFormat`, `delimiter`, `numberNucleotideToExtend`, `zeroBased`, `columnForChromosome`
- **VarDictLauncher**: `regionOfInterest`, `bed`, `bam`, `fasta`, `adaptor`, `threads`, `outputSplicing`, `sampleName`, `sampleNameRegexp`

### Configuration is Written by
- **CmdParser.parseCmd()** — sole writer, creates and populates the object.

### GlobalReadOnlyScope is Read by
- **Every pipeline module** via `instance().conf.*`
- **StructuralVariantsProcessor** via `instance().chrLengths`
- **RecordPreprocessor** via `instance().adaptorForward/Reverse`
- **OutputVariant** via `instance().sample`, `instance().samplem`
- **Mode classes** via `getMode()`

## Known Parity Traps

### Trap 1: `performLocalRealignment` Default is TRUE
Java: `-k` defaults to `1` → `performLocalRealignment = true`. The field initializer is `false`, but `CmdParser` always calls `1 == getIntValue(cmd, "k", 1)` which is TRUE unless user passes `-k 0`.
**Rust risk**: If Rust defaults to `false` and doesn't mirror this `-k 1` default logic, realignment will be silently disabled.

### Trap 2: `samfilter` is a String in Java, u32 in Rust
Java stores `"0x504"` as String. The actual bitwise filtering happens in RecordPreprocessor where it's parsed at usage. Rust parses it to `u32` at CLI time. This is fine for parity but the hex parsing must handle both `0x504` and `0X504` prefixes.

### Trap 3: `zeroBased` is a Boxed Boolean (tristate)
Java `Boolean` type → `null`, `true`, or `false`. The `null` state means "not user-defined" — VarDictLauncher auto-detects from BED format. `isZeroBasedDefined()` checks `zeroBased != null`.
**Rust**: Must use `Option<bool>` and preserve the tristate semantics. Simply defaulting to `false` would break amplicon detection.

### Trap 4: `-p` (pileup) Side-Effects
When `-p` is set, CmdParser ALSO sets `freq = -1` and `minr = 0`. These overwrite values from `-f` and `-r` flags regardless of parse order (because `-p` logic runs after those fields are set).
**Rust**: Must apply the same side-effects after all args are parsed.

### Trap 5: `readPosFilter` is `int` in Java, `f64` in Rust
Java: `int readPosFilter = 5` / `-P` parsed as `int`. Rust: `pub read_pos_filter: f64` with default `5.0`. This shouldn't cause output differences since it's compared with `>` in ToVarsBuilder, but downstream formatting may differ.

### Trap 6: BED Column Indexing — CLI is 1-based, Internal is 0-based
`getColumnValue()` subtracts 1: `((Number)value).intValue() - 1`. The default `BedRowFormat` is `(2, 6, 7, 9, 10, 12)` — already 0-based for refGene format. If `-S` is given as `1` on CLI, internal becomes `0`. Rust CLI uses 1-based defaults in `Args` and must apply the same subtraction.

### Trap 7: `mappingQuality` is Boxed Integer (nullable)
Java: `Integer mappingQuality = null`. Only set if `-Q` present. `hasMappingQuality()` checks `!= null`. "Not set" means NO filtering. If Rust defaults to `0`, that would incorrectly filter reads with mapq 0.
**Rust**: Maps to `Option<u8>` — correct. Java type is `Integer` (32-bit), Rust uses `u8` (0-255). Mapq >255 would overflow but SAM mapq is 0-255 so this is safe.

### Trap 8: `downsampling` is Boxed Double (nullable)
Java: `Double downsampling = null`. Only set if `-Z` present. `isDownsampling()` checks `!= null`.
**Rust**: `Option<f64>` — correct.

### Trap 9: `bias` Has No CLI Flag
`double bias = 0.05d` — hardcoded, never changeable from command line. But ToVarsBuilder reads it.
**Rust**: Ensure the value `0.05` is used identically.

### Trap 10: `indelsize` Not in Rust Configuration
Java: `int indelsize = 50` stored in Configuration, used by ToVarsBuilder. Rust `conf.rs` has no `indelsize` field. May be hardcoded elsewhere — needs verification.

### Trap 11: `doPileup` Not in Rust Configuration
Java: `boolean doPileup = false` — controls whether variant frequency threshold is ignored. Rust has it in `PipelineConfig.pileup` instead. Need to verify parity of all `doPileup` usage paths.

### Trap 12: `outputSplicing` / `y` (verbose) / `printHeader` Not in Rust Configuration
Several Java Configuration fields are handled at the CLI level in Rust and not persisted in Configuration. Architecturally different but not a parity issue as long as logic reaches the same code paths.

### Trap 13: `BedRowFormat` Column Defaults Differ
Java DEFAULT_BED_ROW_FORMAT: `(chrColumn=2, startColumn=6, endColumn=7, thickStartColumn=9, thickEndColumn=10, geneColumn=12)` — 0-based for refGene.txt. Rust CLI defaults: `-c 1 -S 2 -E 3 -g 4` — 1-based for simple BED. Rust subtracts 1 internally → functional defaults `(0, 1, 2, 3)`. OK for parity tests using standard BED files.

### Trap 14: `threads` Defaults and Detection
Java: `Math.max(readThreadsCount(cmd), 1)`. If `-th` without value → `Runtime.getRuntime().availableProcessors()`. If not given → 0 → clamped to 1. Rust: `-t` defaults to `1` explicitly. No auto-detection. Functionally equivalent for default case.

### Trap 15: DISCPAIRQUAL Type Mismatch
Java: `public static final int DISCPAIRQUAL = 35`. Rust: `pub(crate) const DISCPAIRQUAL: f64 = 35.0`. If used in integer comparisons in Java, Rust float comparison may differ for edge cases. Verify usage sites.

### Trap 16: `EXTENSION` Constant Missing from Rust
Java: `public static final int EXTENSION = 5000` — used in ReferenceResource. Rust: Not found in `conf.rs` — may be defined elsewhere or handled differently.

### Trap 17: `validationStringency` Not Mapped
Java: defaults to `LENIENT`, configurable via `-VS`. Affects BAM reading error handling. Rust uses `rust-htslib` with its own defaults. Parity impact unlikely but malformed BAM handling may differ.

### Trap 18: Sample Name Regex Patterns
Java: `SAMPLE_PATTERN` = `/([^\\/\\._]+?)_[^\\/]*.bam/` (non-greedy `+?`).
Rust: `r"([^/\._]+)\.sorted[^/]*\.bam$"` — **DIFFERENT PATTERN** (uses `\.sorted` instead of `_`).
This produces different sample names from the same BAM filename. Usually overridden by `-N`.

### Trap 19: GlobalReadOnlyScope Rust Missing `sample`/`samplem`
Java stores both sample names in GlobalReadOnlyScope. Rust has no sample name fields there. Sample names passed through the pipeline differently. Architectural divergence, not a parity issue unless code reads `instance().sample`.

### Trap 20: `mapq` Field Type — Java `double`, Rust `f64`
Both are 64-bit floating point. Default both `0.0`. Parity safe.
