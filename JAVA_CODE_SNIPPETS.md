# VarDict Java - Key Code Snippets & Reference

## SimpleMode.java - The Pipeline Orchestrator

### Header Output
```java
@Override
public void printHeader() {
    if (instance().conf.printHeader) {
        String header = join("\t",
                "Sample", "Gene", "Chr", "Start", "End", "Ref", "Alt", 
                "Depth", "AltDepth", "RefFwdReads", "RefRevReads", 
                "AltFwdReads", "AltRevReads", "Genotype", "AF", "Bias", 
                "PMean", "PStd", "QMean", "QStd", "MQ", "Sig_Noise", 
                "HiAF", "ExtraAF", "shift3", "MSI", "MSI_NT", "NM",
                "HiCnt", "HiCov", "5pFlankSeq", "3pFlankSeq", "Seg", 
                "VarType", "Duprate", "SV_info");
        if (instance().conf.crisprCuttingSite != 0) {
            header = join("\t", header, "CRISPR");
        }
        System.out.println(header);
    }
}
```

**Output Columns Explained:**
- `Sample` - Sample name
- `Gene` - Gene name (if available)
- `Chr` - Chromosome
- `Start` / `End` - Variant position(s)
- `Ref` / `Alt` - Reference and variant alleles
- `Depth` - Total coverage at position
- `AltDepth` - Variant coverage
- `RefFwdReads` / `RefRevReads` - Reference strand counts
- `AltFwdReads` / `AltRevReads` - Variant strand counts
- `Genotype` - Predicted genotype (HOM, HET)
- `AF` - Allele frequency
- `Bias` - Strand bias (0/1/2)
- `PMean` / `PStd` - Mean/std of position in read
- `QMean` / `QStd` - Mean/std of base quality
- `MQ` - Mean mapping quality
- `Sig_Noise` - Signal-to-noise ratio
- `HiAF` - Frequency in high-quality reads only
- `ExtraAF` - Extra frequency from realigned reads
- `shift3` - 3' shift for deletions
- `MSI` / `MSI_NT` - Microsatellite instability info
- `NM` - Mismatch count in reads
- `HiCnt` - High-quality read count
- `HiCov` - High-quality coverage
- `5pFlankSeq` - 20bp upstream
- `3pFlankSeq` - 20bp downstream
- `Seg` - Region identifier
- `VarType` - Variant type (SNP, Insertion, Deletion, Complex)
- `Duprate` - Duplication rate
- `SV_info` - Structural variant info
- `CRISPR` - CRISPR cut site distance (if enabled)

### Pipeline Creation
```java
private void processBamInPipeline(Region region, VariantPrinter out) {
    Reference reference = tryToGetReference(region);
    
    Scope<InitialData> initialScope = new Scope<>(
        instance().conf.bam.getBam1(),
        region,
        reference,
        referenceResource, 
        0, 
        new HashSet<>(),
        out, 
        new InitialData()
    );

    // Create pipeline with modules chained together
    CompletableFuture<Scope<AlignedVarsData>> pipeline = 
        pipeline(initialScope, new DirectThreadExecutor());
    
    CompletableFuture<Void> simpleProcessOutput = pipeline
            .thenAccept(new SimplePostProcessModule(out))
            .exceptionally(ex -> {
                stopVardictWithException(region, ex);
                throw new RuntimeException(ex);
            });
    
    simpleProcessOutput.join();
}
```

---

## Variant.java - Core Data Structure

### Field Definitions

```java
public class Variant {
    // === IDENTITY ===
    public String descriptionString;      // SNP='A', Insertion='+ACGT', Deletion='-5'
    public String refallele;              // Reference allele
    public String varallele;              // Variant allele
    public String vartype;                // SNP, Insertion, Deletion, Complex
    
    // === POSITION ===
    public int startPosition;
    public int endPosition;
    
    // === COVERAGE ===
    public int positionCoverage;          // Total coverage
    public int varsCountOnForward;        // Variant on forward strand
    public int varsCountOnReverse;        // Variant on reverse strand
    public int refForwardCoverage;        // Reference on forward
    public int refReverseCoverage;        // Reference on reverse
    public int totalPosCoverage;          // Total at position
    
    // === FREQUENCY ===
    public double frequency;              // AF = varsCount / positionCoverage
    public double highQualityReadsFrequency;  // Freq in high-quality reads only
    public double extraFrequency;         // Extra freq from realigned reads
    
    // === POSITION IN READ ===
    public double meanPosition;           // Where in read does variant occur?
    public boolean isAtLeastAt2Positions; // Found at 2+ different positions?
    
    // === BASE QUALITY ===
    public double meanQuality;            // Mean base quality of variant bases
    public boolean hasAtLeast2DiffQualities; // 2+ different quality values?
    
    // === MAPPING QUALITY ===
    public double meanMappingQuality;     // Mean MAPQ of reads with variant
    
    // === STRAND BIAS ===
    public String strandBiasFlag;         // "0" (none), "1" (weak), "2" (strong)
    public double highQualityToLowQualityRatio;
    
    // === MICROSATELLITE ===
    public double msi;                    // MSI score
    public int msint;                     // MSI unit length
    
    // === CONTEXT ===
    public String leftseq;                // 20 bases upstream
    public String rightseq;               // 20 bases downstream
    
    // === INDEL SPECIFICS ===
    public int shift3;                    // 3' shift for deletion representation
    
    // === READ QUALITY ===
    public double numberOfMismatches;     // Avg mismatches in reads with variant
    public int hicnt;                     // Count of high-quality reads
    public int hicov;                     // Coverage by high-quality reads
    
    // === OTHER ===
    public double duprate;                // Duplication rate
    public String genotype;               // HOM, HET, etc.
    public int crispr;                    // CRISPR distance
    public String DEBUG;                  // Debug info
}
```

### Key Methods

```java
/**
 * Determines if variant type (SNP, insertion, deletion, complex)
 */
public String varType() {
    if (descriptionString.length() == 1) {
        return "SNP";  // Single letter = SNP
    } else if (descriptionString.startsWith("+")) {
        return "Insertion";
    } else if (descriptionString.startsWith("-")) {
        return "Deletion";
    } else if (descriptionString.contains("#") || descriptionString.contains("&")) {
        return "Complex";  // Indel with extra bases
    }
    return "Unknown";
}

/**
 * Check if variant passes quality filters
 */
public boolean isGoodVar(Variant refVar, String varType, Set<String> splice) {
    // 1. Strand bias check
    if ("2".equals(strandBiasFlag)) {
        // Strong strand bias - reject unless in specific regions
        if (!isInSpecialRegion(splice)) {
            return false;
        }
    }
    
    // 2. Position variance check
    if (!isAtLeastAt2Positions) {
        // Variant found at only 1 position in all reads
        return false;
    }
    
    // 3. Quality variance check
    if (!hasAtLeast2DiffQualities) {
        // All reads have same quality - might be error
        return false;
    }
    
    // 4. MSI region handling
    if (isMSIRegion()) {
        // Relaxed thresholds in MSI
    }
    
    return true;  // Passed all checks
}

/**
 * For complex variants, adjust representation
 */
public void adjComplex() {
    // Combine adjacent SNPs and indels into single complex variant
    // Recalculate coordinates
}

/**
 * Check if is noise (low quality, low count)
 */
public boolean isNoise() {
    final double qual = this.meanQuality;
    if (((qual < 4.5d || (qual < 12 && !this.hasAtLeast2DiffQualities)) 
            && this.positionCoverage <= 3)) {
        return true;
    }
    return false;
}
```

---

## Variant.java - Quality Analysis Fields

```java
/**
 * These fields are POPULATED by ToVarsBuilder during statistics calculation
 */

// Position in read distribution
double meanPosition;              // Average position of variant in reads
boolean isAtLeastAt2Positions;   // Variant found at positions like 10, 20, 30? (2+ positions)

// Quality score distribution  
double meanQuality;              // Average base quality
boolean hasAtLeast2DiffQualities; // Quality is 20 in some reads, 25 in others? (2+ qualities)

// Mapping quality
double meanMappingQuality;       // Average MAPQ of reads with variant

// Strand bias metrics
String strandBiasFlag;           // "0" = balanced, "1" = slight bias, "2" = strong bias
double highQualityToLowQualityRatio; // HQ(Q≥20) count / LQ(<Q20) count

// High-quality analysis
int hicnt;                       // Number of high-quality reads (Q≥20)
int hicov;                       // Coverage by high-quality reads
double highQualityReadsFrequency; // Frequency in high-quality reads only

// Extra frequency from realignment
double extraFrequency;           // Additional frequency detected in realigned reads

// Microsatellite analysis
double msi;                      // Microsatellite instability (>1 = instable)
int msint;                       // MSI unit length in bp
```

### isGoodVar() - Quality Threshold Logic

```
A variant is GOOD if:
1. Strand bias is acceptable (not strong without cause)
2. isAtLeastAt2Positions = true (found at 2+ positions in reads)
3. hasAtLeast2DiffQualities = true (found with 2+ quality levels)
4. Not in obvious MSI region (or MSI thresholds met)
5. Frequency meets minimum threshold
6. No unexpected strand bias patterns
```

**Why these checks?**
- **Multiple positions:** Real variants spread across read; sequencing errors cluster
- **Multiple qualities:** Real variants have quality variation; errors are uniform
- **Strand balance:** Errors often one-stranded; real variants appear both
- **MSI handling:** Repeat regions have higher error rate; use relaxed thresholds

---

## ToVarsBuilder.java - Variant Creation

### Statistics Calculation Algorithm

```java
/**
 * For each position and each variant at that position:
 * 1. Count occurrences
 * 2. Calculate distribution metrics
 * 3. Create Variant object
 */
public Scope<AlignedVarsData> accept(Scope<RealignedVariationData> scope) {
    initFromScope(scope);
    
    // Group variations by position
    Map<Integer, Vars> alignedVariants = new TreeMap<>();
    
    for (Integer position : variationMap.keySet()) {
        Vars varsAtPosition = new Vars();
        List<Variant> variants = new ArrayList<>();
        
        for (Variation variation : variationMap.get(position)) {
            // === CREATE VARIANT OBJECT ===
            Variant variant = new Variant();
            
            // === BASIC INFO ===
            variant.startPosition = variation.position;
            variant.endPosition = variation.endPosition;
            variant.descriptionString = variation.description;
            
            // === ALLELES ===
            variant.refallele = getRefAllele(position, variation);
            variant.varallele = getVarAllele(variation);
            variant.vartype = variant.varType();
            
            // === COVERAGE ===
            variant.positionCoverage = refCoverage.getOrDefault(position, 0);
            variant.varsCountOnForward = countFromForwardStrand(variation);
            variant.varsCountOnReverse = countFromReverseStrand(variation);
            
            // === FREQUENCY ===
            variant.frequency = (double)(variant.varsCountOnForward + 
                                        variant.varsCountOnReverse) 
                               / variant.positionCoverage;
            
            // === POSITION IN READ ANALYSIS ===
            List<Integer> positions = getPositionsInReads(variation);
            variant.meanPosition = positions.stream()
                                           .mapToInt(Integer::intValue)
                                           .average()
                                           .orElse(0);
            variant.isAtLeastAt2Positions = (positions.stream()
                                                      .distinct()
                                                      .count() >= 2);
            
            // === QUALITY ANALYSIS ===
            List<Integer> qualities = getQualitiesOfVariantBases(variation);
            variant.meanQuality = qualities.stream()
                                          .mapToInt(Integer::intValue)
                                          .average()
                                          .orElse(0);
            variant.hasAtLeast2DiffQualities = (qualities.stream()
                                                        .distinct()
                                                        .count() >= 2);
            
            // === MAPPING QUALITY ===
            List<Integer> mapqs = getMappingQualities(variation);
            variant.meanMappingQuality = mapqs.stream()
                                             .mapToInt(Integer::intValue)
                                             .average()
                                             .orElse(0);
            
            // === STRAND BIAS ===
            variant.strandBiasFlag = calculateStrandBias(variant);
            
            // === HIGH-QUALITY READS ===
            variant.hicnt = countHighQualityReads(variation);
            variant.hicov = countHighQualityCoverage(position);
            variant.highQualityReadsFrequency = 
                (double) variant.hicnt / variant.hicov;
            
            // === FLANKING SEQUENCE ===
            variant.leftseq = getFlankingSeq(position - 20, position - 1);
            variant.rightseq = getFlankingSeq(position + 1, position + 20);
            
            // === MICROSATELLITE ===
            variant.msi = calculateMSI(position);
            variant.msint = calculateMSIUnitLength(position);
            
            variants.add(variant);
        }
        
        varsAtPosition.variants = variants;
        varsAtPosition.referenceVariant = createReferenceVariant(position);
        alignedVariants.put(position, varsAtPosition);
    }
    
    // Return wrapped in scope
    return new Scope<>(scope, new AlignedVarsData(alignedVariants));
}
```

### Key Calculation Examples

#### Strand Bias Calculation
```java
private String calculateStrandBias(Variant v) {
    int fwd = v.varsCountOnForward;
    int rev = v.varsCountOnReverse;
    int total = fwd + rev;
    
    if (total == 0) return "0";
    
    // Check if heavily skewed to one strand
    double ratio = Math.max(fwd, rev) / (double) total;
    
    if (ratio > 0.9) {
        // >90% on one strand - strong bias
        return "2";
    } else if (ratio > 0.75) {
        // 75-90% - weak bias
        return "1";
    } else {
        // Balanced
        return "0";
    }
}
```

#### Microsatellite Detection
```java
private double calculateMSI(int position) {
    // Look at 30bp window around position
    String window = referenceSequence.substring(position - 15, position + 15);
    
    // Count repeating patterns
    // ACACAC = high MSI (repeat length 2)
    // AAAA = high MSI (repeat length 1)
    // ACGTACGTACGT = high MSI (repeat length 4)
    
    int repeatLength = findRepeatLength(window);
    int repeatCount = countRepeats(window, repeatLength);
    
    return (double) repeatCount / repeatLength;
}
```

#### High-Quality Read Filtering
```java
private int countHighQualityReads(Variation v) {
    int count = 0;
    for (SAMRecord record : v.reads) {
        int baseQuality = record.getBaseQualities()[v.positionInRead];
        if (baseQuality >= QUALITY_THRESHOLD) {  // Usually 20
            count++;
        }
    }
    return count;
}
```

---

## SimplePostProcessModule.java - Output Filtering

### Accept Method - The Final Filter

```java
@Override
public void accept(Scope<AlignedVarsData> mapScope) {
    int lastPosition = 0;
    
    // Iterate through all positions with variants
    for (Map.Entry<Integer, Vars> ent : mapScope.data.alignedVariants.entrySet()) {
        try {
            int position = ent.getKey();
            lastPosition = position;
            Vars variantsOnPosition = ent.getValue();
            
            Configuration conf = instance().conf;
            
            List<Variant> vrefs = new ArrayList<>();
            
            // === REGION BOUNDARY CHECK ===
            if (variantsOnPosition.sv.isEmpty()) {
                if (position < mapScope.region.start || 
                    position > mapScope.region.end) {
                    continue;  // Skip out-of-bounds
                }
            }
            
            // === EMPTY VARIANTS (PILEUP MODE) ===
            if (variantsOnPosition.variants.isEmpty()) {
                if (!conf.doPileup) {
                    continue;  // Skip if not in pileup mode
                }
                Variant vref = variantsOnPosition.referenceVariant;
                if (vref == null) {
                    SimpleOutputVariant outputVariant = 
                        new SimpleOutputVariant(vref, mapScope.region, 
                                             variantsOnPosition.sv, position);
                    variantPrinter.print(outputVariant);
                    continue;
                }
                vref.vartype = "";
                vrefs.add(vref);
            } else {
                // === FILTER ACTUAL VARIANTS ===
                List<Variant> vvar = variantsOnPosition.variants;
                for (Variant vref : vvar) {
                    // Skip if ref contains N
                    if (vref.refallele.contains("N")) {
                        continue;
                    }
                    
                    // Handle ref calls (pileup mode)
                    if (vref.refallele.equals(vref.varallele)) {
                        if (!conf.doPileup) {
                            continue;
                        }
                    }
                    
                    // === KEY QUALITY FILTER ===
                    if (!vref.isGoodVar(variantsOnPosition.referenceVariant, 
                                       vref.vartype, mapScope.splice)) {
                        if (!conf.doPileup) {
                            continue;
                        }
                    }
                    
                    // Determine variant type
                    vref.vartype = vref.varType();
                    
                    // Adjust complex variants
                    if ("Complex".equals(vref.vartype)) {
                        vref.adjComplex();
                    }
                    
                    // CRISPR filtering
                    if (instance().conf.crisprCuttingSite == 0) {
                        vref.crispr = 0;
                    }
                    
                    vrefs.add(vref);
                }
            }
            
            // === OUTPUT ===
            for (int vi = 0; vi < vrefs.size(); vi++) {
                Variant vref = vrefs.get(vi);
                SimpleOutputVariant outputVariant = 
                    new SimpleOutputVariant(vref, mapScope.region, 
                                          variantsOnPosition.sv, position);
                variantPrinter.print(outputVariant);
            }
            
        } catch (Exception e) {
            printExceptionAndContinue(e, region, lastPosition);
        }
    }
}
```

---

## Configuration.java - Key Settings

```java
public class Configuration {
    // Files
    public BamFiles bam;              // Input BAM file(s)
    public String reference;          // Reference FASTA
    public String outputFile;         // Output file path
    
    // Region selection
    public Region region;             // Single region
    public List<Region> regions;      // Multiple regions
    public String bedFile;            // BED file path
    
    // Quality thresholds
    public double minFreq = 0.02;     // Min allele frequency (2%)
    public int minQuality = 0;        // Min base quality
    public int minMappingQuality = 0; // Min MAPQ
    public int minAltCount = 2;       // Min variant count (2 reads)
    
    // Analysis modes
    public boolean doPileup = false;       // Output reference calls
    public boolean doStrandBias = false;   // Filter by strand bias
    public int crisprCuttingSite = 0;      // CRISPR cut position (0=disabled)
    
    // Output
    public PrinterTypeOut printerTypeOut = PrinterTypeOut.TEXT;
    public boolean printHeader = true;
    
    // Execution
    public int threads = 1;               // Number of parallel threads
    public boolean useParallel = false;   // Enable parallel processing
}
```

---

## Data Flow Summary

```
Input Files:
  - BAM (indexed)
  - BED (regions)
  - Reference FASTA (indexed)

↓

Configuration (settings, thresholds)

↓ SimpleMode.processBamInPipeline() ↓

Module Chain:
  SAMFileParser 
    → InitialData (reads, coverage)
  ↓
  RecordPreprocessor 
    → InitialData (filtered reads)
  ↓
  CigarParser 
    → VariationData (variations at each position)
  ↓
  VariationRealigner 
    → RealignedVariationData (rescued indels added)
  ↓
  ToVarsBuilder ⭐⭐⭐
    → AlignedVarsData (Variant objects with stats)
  ↓
  SimplePostProcessModule ⭐
    → Filters variants (isGoodVar())
  ↓
  VariantPrinter
    → SimpleOutputVariant → Output (VCF/text format)
```

---

## Key Performance Insights

1. **Bottleneck:** BAM file I/O (sequential access required)
2. **Parallelizable:** Region processing (independent)
3. **Hot Path:** ToVarsBuilder statistics calculation
4. **Memory:** Variation maps (scale with coverage depth)
5. **CPU:** Sorting/grouping variations, quality calculations

## Rust Port Implementation Strategy

1. **Use rust-htslib** for BAM/SAM reading
2. **HashMap<i32, Vec<Variation>>** for variation grouping
3. **Rayon** for parallel region processing
4. **Custom struct** for Variant (replicate Java fields)
5. **SmallVec** for position/quality vectors (avoid heap)
6. **String interning** for alleles (reduce duplicates)

