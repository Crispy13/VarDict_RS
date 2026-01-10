# SimpleMode Deep Dive - Detailed Pipeline Analysis

## Complete Pipeline Walkthrough

### Entry Point: VarDictLauncher.start()

1. **Configuration Parsing** (CmdParser)
   - Command-line arguments → Configuration object
   - Sets: BAM files, BED regions, reference, thresholds

2. **Mode Selection** (VarDictLauncher)
   - If single sample → SimpleMode
   - If paired sample → SomaticMode
   - If amplicon → AmpliconMode
   - etc.

3. **SimpleMode Creation** (SimpleMode constructor)
   ```java
   public SimpleMode(List<List<Region>> segments, ReferenceResource referenceResource) {
       super(segments, referenceResource);
       printHeader();  // Print VCF header
   }
   ```

### Pipeline Execution

#### Path A: Sequential Execution (notParallel)
```java
@Override
public void notParallel() {
    VariantPrinter variantPrinter = VariantPrinter.createPrinter(
        instance().printerTypeOut
    );

    for (List<Region> list : segments) {
        for (Region region : list) {
            processBamInPipeline(region, variantPrinter);
        }
    }
}
```

**Execution:** Regions processed one-by-one in main thread

#### Path B: Parallel Execution (createParallelMode)
```java
protected AbstractParallelMode createParallelMode() {
    return new AbstractParallelMode() {
        @Override
        void produceTasks() throws InterruptedException {
            for (List<Region> list : segments) {
                for (Region region : list) {
                    // Submit each region to thread pool
                    Future<OutputStream> submit = executor.submit(
                        new VardictWorker(region)
                    );
                    toPrint.put(submit);
                }
            }
            toPrint.put(LAST_SIGNAL_FUTURE);
        }
    };
}

private class VardictWorker implements Callable<OutputStream> {
    @Override
    public OutputStream call() {
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        PrintStream out = new PrintStream(baos);
        VariantPrinter variantPrinter = VariantPrinter.createPrinter(
            instance().printerTypeOut
        );
        variantPrinter.setOut(out);
        processBamInPipeline(region, variantPrinter);
        out.close();
        return baos;
    }
}
```

**Execution:** Each region in separate thread, output ordered in queue

---

## Core Pipeline: processBamInPipeline()

This is THE function that does all the work. Let me break it down:

```java
private void processBamInPipeline(Region region, VariantPrinter out) {
    Reference reference = tryToGetReference(region);
    
    // === STEP 1: Create Initial Scope ===
    Scope<InitialData> initialScope = new Scope<>(
        instance().conf.bam.getBam1(),  // BAM file
        region,                          // Region (chr, start, end)
        reference,                       // Reference sequence
        referenceResource,               // Resource mgmt
        0,                              // maxReadLength (will be updated)
        new HashSet<>(),                // splice sites
        out,                            // output printer
        new InitialData()               // empty data
    );

    // === STEP 2: Parse BAM file, filter reads ===
    SAMFileParser parser = new SAMFileParser();
    Scope<InitialData> parsedScope = parser.accept(initialScope);
    // Result: BAM records parsed into InitialData
    
    // === STEP 3: Preprocess SAM records ===
    RecordPreprocessor preprocessor = new RecordPreprocessor();
    Scope<InitialData> preprocessedScope = preprocessor.accept(parsedScope);
    // Result: Reads filtered, bad quality removed
    
    // === STEP 4: Parse CIGAR strings ===
    CigarParser cigarParser = new CigarParser();
    Scope<VariationData> variationScope = cigarParser.accept(preprocessedScope);
    // Result: CIGAR → Variation objects (SNPs, indels found in aligned reads)
    
    // === STEP 5: Realign soft-clipped reads ===
    VariationRealigner realigner = new VariationRealigner();
    Scope<RealignedVariationData> realignedScope = realigner.accept(variationScope);
    // Result: Additional variations from soft-clipped sequence realignment
    
    // === STEP 6: Build Variant objects ===
    ToVarsBuilder varsBuilder = new ToVarsBuilder();
    Scope<AlignedVarsData> alignedVarsScope = varsBuilder.accept(realignedScope);
    // Result: Variant objects with statistics (freq, quality, bias, etc.)
    
    // === STEP 7: Post-process and output ===
    SimplePostProcessModule postProcessor = 
        new SimplePostProcessModule(out);
    postProcessor.accept(alignedVarsScope);
    // Result: Filtered variants printed to output
}
```

---

## STEP-BY-STEP DATA TRANSFORMATIONS

### STEP 1: SAMFileParser - BAM → Reads

**Input:** `Scope<InitialData>` with BAM file path

**Process:**
1. Open BAM file using SAM/BAM reader (HTSJDK)
2. Query for reads overlapping region
3. For each SAM record:
   - Check if properly paired (if paired-end)
   - Check read length, quality flags
   - Extract:
     - Query name (read ID)
     - Sequence
     - Base qualities
     - Mapping quality
     - CIGAR string
     - Position
     - Flags (reverse complement, duplicate, etc.)

**Output:** `Scope<InitialData>` with `reads` list populated

**Key Classes:**
- `SAMFileParser` (module)
- `SAMRecord` (from HTSJDK)
- `InitialData` (stores: `reads`, `maxReadLength`, `refCoverage`)

---

### STEP 2: RecordPreprocessor - Filter Reads

**Input:** Parsed SAM records

**Process:**
1. Check if read passes QC:
   - Not a duplicate
   - Not QC failed
   - Properly mapped (if pair)
   - Adequate mapping quality (≥ minMappingQuality)
2. Calculate position coverage for reference
3. Update max read length

**Output:** Filtered reads + coverage map

**Key Filters:**
- `SAMFlag.READ_UNMAPPED` → skip
- `SAMFlag.READ_FAILS_VENDOR_QUALITY_CHECK` → skip
- `mappingQuality < minMappingQuality` → skip
- `isSecondary || isSupplementary` → skip

---

### STEP 3: CigarParser - Detect Variations in Aligned Reads

**Input:** Filtered SAM records with CIGAR strings

**CIGAR String Example:**
```
ACGTACGTACGT (query)
10M2I3M1D5M   (CIGAR)
├─ 10M: alignment match (may contain SNPs)
├─ 2I: insertion of 2 bases
├─ 3M: 3 more matches
├─ 1D: deletion of 1 base
└─ 5M: 5 more matches
```

**Process:**
1. Parse CIGAR operations
2. For each operation:
   - **M (match):** Compare query to reference base
     - Match → coverage for reference
     - Mismatch → SNP variation
   - **I (insertion):** Extract inserted sequence
     - Create insertion variation at position
   - **D (deletion):** Count deleted bases
     - Create deletion variation
   - **S (soft clip):** Store for later realignment
   - **N (skip):** Skip region (for spliced reads)

**Output:** `Scope<VariationData>`

**Key Classes:**
- `CigarParser` (module)
- `Variation` (low-level representation)
- `VariationMap` (map of position → variations)

**Variation Types Detected:**
1. SNPs (single base mismatch)
2. Insertions (extra bases in query)
3. Deletions (missing bases in query)
4. Soft clips (clipped sequences)

---

### STEP 4: VariationRealigner - Rescue Hidden Indels

**This is a KEY VarDict feature!**

**Input:** Variations found above + soft-clipped sequences

**Problem:** Large indels often manifest as soft clipping:
```
Reference: ACGTAAACGT
Query:     ACGTNN      (soft-clipped, 4 bases soft-clipped)

Actually: ACGT-4ACGT (4-base deletion, not visible in CIGAR without realignment)
```

**Process:**
1. Extract soft-clipped sequences (from S operations)
2. Align clipped sequence locally to reference
3. If alignment explains clipping better:
   - Create new variation (hidden indel)
   - Mark as realigned

**Output:** `Scope<RealignedVariationData>` (now includes rescued indels)

**Key Improvement:** Detects indels that would be missed otherwise

---

### STEP 5: ToVarsBuilder - Statistics & Variant Creation ⭐⭐⭐

**Input:** All variations (SNPs, indels, realigned)

**This is the MOST COMPLEX step!**

**Process:**

#### 5a. Group Variations by Position
```java
Map<Integer, Variation> variations → Map<Integer, List<Variation>>
```

#### 5b. For Each Position, Calculate Statistics

**For each variant at position:**

1. **Count strand occurrence:**
   - How many forward strand reads?
   - How many reverse strand reads?
   - Store in: `varsCountOnForward`, `varsCountOnReverse`

2. **Calculate frequency:**
   ```
   frequency = (forward_count + reverse_count) / total_coverage
   ```
   - Stored in: `frequency`

3. **Position in read analysis:**
   - For each read with this variant:
     - Where in the read does this variant occur? (start? middle? end?)
     - Average position
     - Variance/standard deviation
   - Flag if found at 2+ different positions: `isAtLeastAt2Positions`
   - Stored in: `meanPosition`, `isAtLeastAt2Positions`

4. **Base quality analysis:**
   - Quality (Q-score) of base at variant position
   - Average quality
   - Variance in quality
   - Flag if 2+ different quality levels: `hasAtLeast2DiffQualities`
   - Stored in: `meanQuality`, `hasAtLeast2DiffQualities`

5. **Mapping quality analysis:**
   - MAPQ (mapping quality) of reads with variant
   - Average mapping quality
   - Stored in: `meanMappingQuality`

6. **Strand bias detection:**
   - If all variants on one strand → strong bias
   - Calculate using fisher test (optional)
   - Stored in: `strandBiasFlag` (0/1/2)

7. **High-quality vs Low-quality split:**
   - Split reads by base quality (e.g., ≥Q20 vs <Q20)
   - Calculate frequency for high-quality only
   - Ratio of high to low
   - Stored in: `highQualityReadsFrequency`, `highQualityToLowQualityRatio`

8. **Shift3 (3' positioning):**
   - For deletions, preferred 3' representation
   - How many bases can we shift 3' while maintaining del?
   - Stored in: `shift3`

#### 5c. Create Variant Objects
```java
Variation → Variant (with all statistics)
```

**Output:** `Scope<AlignedVarsData>`
- `alignedVariants`: Map<Integer, Vars> (position → variants at position)
- `Vars` contains:
  - `variants`: List<Variant>
  - `referenceVariant`: Variant (for pileup)
  - `sv`: Structural variant flags

---

### STEP 6: SimplePostProcessModule - Quality Filtering & Output

**Input:** `Scope<AlignedVarsData>` with all variants

**Process:**

```java
@Override
public void accept(Scope<AlignedVarsData> mapScope) {
    for (Map.Entry<Integer, Vars> ent : mapScope.data.alignedVariants.entrySet()) {
        int position = ent.getKey();
        Vars variantsOnPosition = ent.getValue();
        
        // === Filter out-of-region calls ===
        if (position < mapScope.region.start || position > mapScope.region.end) {
            continue;  // Skip if outside region (unless SV)
        }
        
        // === Handle empty variants (pileup mode) ===
        if (variantsOnPosition.variants.isEmpty()) {
            if (!conf.doPileup) continue;  // Skip if not in pileup mode
            // Output reference call
        }
        
        // === Filter and process variants ===
        List<Variant> toOutput = new ArrayList<>();
        for (Variant vref : variantsOnPosition.variants) {
            // Skip if ref contains N (unknown base)
            if (vref.refallele.contains("N")) continue;
            
            // Skip reference calls unless pileup mode
            if (vref.refallele.equals(vref.varallele)) {
                if (!conf.doPileup) continue;
            }
            
            // === KEY QUALITY FILTER ===
            if (!vref.isGoodVar(
                    variantsOnPosition.referenceVariant, 
                    vref.vartype, 
                    mapScope.splice)) {
                if (!conf.doPileup) continue;
            }
            
            // === Handle complex variants ===
            if ("Complex".equals(vref.vartype)) {
                vref.adjComplex();  // Adjust coordinates
            }
            
            // === CRISPR filtering ===
            if (instance().conf.crisprCuttingSite == 0) {
                vref.crispr = 0;
            }
            
            toOutput.add(vref);
        }
        
        // === Create and print output ===
        for (Variant vref : toOutput) {
            SimpleOutputVariant outputVariant = new SimpleOutputVariant(
                vref, 
                mapScope.region, 
                variantsOnPosition.sv, 
                position
            );
            variantPrinter.print(outputVariant);
        }
    }
}
```

**Key Filter: `isGoodVar()`**
```
A variant is "good" if:
1. Not marked as BIAS (unless strand bias threshold relaxed)
2. Position in read is sensible (not too close to ends) ✓ isAtLeastAt2Positions
3. Quality variance acceptable ✓ hasAtLeast2DiffQualities
4. MAPQ acceptable
5. For indels: Not in obvious MSI region (too much repeat)
6. Frequency meets threshold
```

**Output:** VCF/text format with fields:
- Position, Reference, Variant
- Depth, Count, Frequency
- Quality metrics (base qual, mapq, strand bias, etc.)
- Variant type
- CRISPR flag (if applicable)

---

## Data Structure Evolution

```
InitialData
├─ reads: List<SAMRecord>
├─ maxReadLength: int
└─ refCoverage: Map<Int, Int>
     ↓
VariationData
├─ reads: List<SAMRecord>
├─ variations: Map<Int, List<Variation>>
├─ refCoverage: Map<Int, Int>
└─ maxReadLength: int
     ↓
RealignedVariationData
├─ insertionVariants: Map<Int, Map<String, Variation>>
├─ nonInsertionVariants: Map<Int, Map<String, Variation>>
├─ refCoverage: Map<Int, Int>
├─ duprate: Double
└─ maxReadLength: int
     ↓
AlignedVarsData
├─ alignedVariants: Map<Int, Vars>  [KEY DATA STRUCTURE]
└─ refCoverage: Map<Int, Int>

Where Vars contains:
├─ variants: List<Variant>
├─ referenceVariant: Variant
└─ sv: Structural variant flags
```

---

## Key Algorithms

### Variant.isGoodVar() - Quality Filter

```java
public boolean isGoodVar(Variant refVar, String varType, Set<String> splice) {
    // 1. Check strand bias
    if (strandBiasFlag != "0") {
        if (stricter_strand_bias_threshold) return false;
    }
    
    // 2. Check position in read distribution
    if (!isAtLeastAt2Positions) {
        // Variant found at only 1 position → might be error
        return false;
    }
    
    // 3. Check quality variance
    if (!hasAtLeast2DiffQualities) {
        // Variant found with only 1 quality → might be error
        return false;
    }
    
    // 4. Check MSI (microsatellite) region
    if (isInMSIRegion()) {
        // In MSI region → relax thresholds
    }
    
    // 5. For indels, check if it's a valid representation
    if (varType contains "indel") {
        // Check for proper anchoring, length limits, etc.
    }
    
    return true;
}
```

### Soft-Clip Realignment Algorithm

```java
VariationRealigner.realignClips() {
    for (SAMRecord record : records) {
        if (has_soft_clips) {
            softClipSeq = extract_soft_clip_sequence(cigar);
            
            // Try to align clipped bases to reference
            localAlignment = smithWaterman(softClipSeq, reference);
            
            if (alignment_explains_clipping) {
                // Found hidden indel
                createVariation(hidden_indel);
            }
        }
    }
}
```

---

## Configuration Parameters (Important for Simple Mode)

```java
Configuration conf = instance().conf;

// Thresholds
conf.minFreq          // Min allele frequency (0.02 default)
conf.minQuality       // Min base quality (0 default)
conf.minMappingQuality // Min MAPQ (0 default)
conf.minAltCount      // Min variant count (2 default)

// Modes
conf.doPileup         // Output ref calls (false default)
conf.doStrandBias     // Filter strand bias (false default)

// Other
conf.bam              // BAM file(s)
conf.reference        // Reference FASTA
conf.crisprCuttingSite // CRISPR position (0 = disabled)
conf.region           // Region to analyze
conf.printerTypeOut   // Output format (VCF, text, etc.)
```

---

## Summary for Rust Port

### Must-Have Components:
1. ✅ BED region parser
2. ✅ Reference FASTA reader
3. Configuration/CLI parser
4. **BAM reader** (rust-htslib)
5. **SAM record parser**
6. **CIGAR parser** ⭐
7. **Variation tracking** (HashMap-based)
8. **Soft-clip realigner** ⭐
9. **Statistics calculator** ⭐
10. **Variant struct** (Rust struct)
11. **Quality filter** (isGoodVar)
12. **Output formatter** (VCF/text)

### Most Complex Parts:
1. **ToVarsBuilder logic** - Statistics calculation
2. **VariationRealigner** - Smith-Waterman alignment
3. **CIGAR parsing** - Correctly handle all operations
4. **Statistics** - Position, quality, mapq distributions

### Performance Opportunities:
1. Parallel region processing (rayon)
2. Efficient variation map (HashMap vs TreeMap)
3. Batch statistics calculation
4. String interning for alleles

