# ToVarsBuilder.java - Rust Porting Preparation Document

**Module:** `com.astrazeneca.vardict.modules.ToVarsBuilder`  
**Purpose:** Convert low-level Variation objects into high-level Variant objects with complete statistics  
**Complexity:** ⭐⭐⭐ HIGH (Core statistics calculation engine)  
**Line Count:** ~1059 lines (Java)  
**Priority:** CRITICAL - This is the heart of VarDict's variant calling logic

---

## 1. CLASS STRUCTURE

### Main Class: `ToVarsBuilder`

#### Fields (State)
```java
// === From Scope (inherited context) ===
private Configuration conf;                    // Global configuration
private ReferenceResource referenceResource;   // Reference genome accessor
private Reference reference;                   // Current region reference
private Region region;                         // Current genomic region
private String bamFile;                        // Input BAM file path

// === Variation Data (from previous steps) ===
private Map<Integer, Variation> variationMap;  // Position → Variation(s)
private Map<Integer, Integer> refCoverage;     // Position → total coverage
private Map<Integer, Integer> insertionCoverage; // Position → insertion coverage

// === Output Data (created here) ===
private Map<Integer, Vars> alignedVariants;    // Position → Vars (contains Variants)

// === Statistics Tracking ===
// For each variant at each position:
private Map<String, List<Integer>> positionsInReads;    // Track where variant appears in reads
private Map<String, List<Integer>> qualitiesOfBases;    // Track base qualities
private Map<String, List<Integer>> mappingQualities;    // Track MAPQ values
private Map<String, Integer> forwardStrandCounts;       // Count forward strand reads
private Map<String, Integer> reverseStrandCounts;       // Count reverse strand reads
```

**Key Points:**
- ToVarsBuilder is **stateful** - it accumulates data as it processes variations
- Works with **maps keyed by position** (Integer) for efficient lookups
- Maintains **separate tracking** for different statistical dimensions
- Heavily uses the **Scope pattern** to pass context through pipeline

---

## 2. MAIN METHOD SIGNATURES

### Public API

#### Primary Entry Point
```java
@Override
public Scope<AlignedVarsData> accept(Scope<RealignedVariationData> scope)
```
- **Input:** `Scope<RealignedVariationData>` containing all variations from realigner
- **Output:** `Scope<AlignedVarsData>` containing Variant objects with statistics
- **Purpose:** Main processing method - orchestrates entire statistics calculation
- **Flow:** Extract variations → Group by position → Calculate stats → Create Variants

---

## 3. GENERAL FLOW OF `accept()` METHOD

### High-Level Algorithm (Simple Mode Path)

```java
public Scope<AlignedVarsData> accept(Scope<RealignedVariationData> scope) {
    // === STEP 1: Initialize from scope ===
    initFromScope(scope);  
    // Extracts: conf, reference, region, variationMap, refCoverage
    
    // === STEP 2: Initialize output map ===
    Map<Integer, Vars> alignedVariants = new TreeMap<>();
    
    // === STEP 3: Process each position ===
    for (Integer position : variationMap.keySet()) {
        
        // === STEP 3a: Get all variations at this position ===
        Collection<Variation> variationsAtPos = variationMap.get(position);
        
        // === STEP 3b: Create Vars container ===
        Vars varsAtPosition = new Vars();
        List<Variant> variants = new ArrayList<>();
        
        // === STEP 3c: Process each variation at position ===
        for (Variation variation : variationsAtPos) {
            
            // === STEP 3c-i: Create Variant object ===
            Variant variant = new Variant();
            
            // === STEP 3c-ii: Set basic identity ===
            variant.startPosition = variation.position;
            variant.endPosition = variation.endPosition;
            variant.descriptionString = variation.description;
            
            // === STEP 3c-iii: Determine alleles ===
            variant.refallele = getRefAllele(position, variation);
            variant.varallele = getVarAllele(variation);
            variant.vartype = variant.varType();  // "SNP", "Insertion", "Deletion", "Complex"
            
            // === STEP 3c-iv: Calculate coverage & counts ===
            variant.positionCoverage = refCoverage.get(position);
            variant.varsCountOnForward = countForwardStrand(variation);
            variant.varsCountOnReverse = countReverseStrand(variation);
            
            // === STEP 3c-v: Calculate frequency ===
            int totalVarCount = variant.varsCountOnForward + variant.varsCountOnReverse;
            variant.frequency = (double) totalVarCount / variant.positionCoverage;
            
            // === STEP 3c-vi: POSITION IN READ STATISTICS ===
            List<Integer> positions = collectPositionsInReads(variation);
            variant.meanPosition = calculateMean(positions);
            variant.positionStdDev = calculateStdDev(positions);
            variant.isAtLeastAt2Positions = hasAtLeast2Distinct(positions);
            
            // === STEP 3c-vii: BASE QUALITY STATISTICS ===
            List<Integer> qualities = collectBaseQualities(variation);
            variant.meanQuality = calculateMean(qualities);
            variant.qualityStdDev = calculateStdDev(qualities);
            variant.hasAtLeast2DiffQualities = hasAtLeast2Distinct(qualities);
            
            // === STEP 3c-viii: MAPPING QUALITY STATISTICS ===
            List<Integer> mapqs = collectMappingQualities(variation);
            variant.meanMappingQuality = calculateMean(mapqs);
            variant.mappingQualityStdDev = calculateStdDev(mapqs);
            
            // === STEP 3c-ix: STRAND BIAS CALCULATION ===
            variant.strandBiasFlag = calculateStrandBias(variant);
            
            // === STEP 3c-x: HIGH-QUALITY READ FILTERING ===
            variant.hicnt = countHighQualityReads(variation, QUAL_THRESHOLD);
            variant.hicov = getHighQualityCoverage(position, QUAL_THRESHOLD);
            variant.highQualityReadsFrequency = (double) variant.hicnt / variant.hicov;
            variant.highQualityToLowQualityRatio = 
                calculateQualityRatio(variation, QUAL_THRESHOLD);
            
            // === STEP 3c-xi: EXTRA FREQUENCY (from realignment) ===
            variant.extraFrequency = calculateExtraFrequency(variation);
            
            // === STEP 3c-xii: MICROSATELLITE INSTABILITY ===
            MicrosatelliteInfo msiInfo = detectMicrosatellite(position, reference);
            variant.msi = msiInfo.instability;
            variant.msint = msiInfo.unitLength;
            
            // === STEP 3c-xiii: SHIFT3 (for deletions) ===
            if (variant.vartype.equals("Deletion")) {
                variant.shift3 = calculateShift3(position, variant, reference);
            }
            
            // === STEP 3c-xiv: FLANKING SEQUENCES ===
            variant.leftseq = reference.substring(position - 20, position);
            variant.rightseq = reference.substring(position + 1, position + 21);
            
            // === STEP 3c-xv: GENOTYPE PREDICTION ===
            variant.genotype = predictGenotype(variant.frequency);
            
            // === STEP 3c-xvi: NOISE & SIGNAL-TO-NOISE ===
            variant.isNoise = variant.isNoise();
            variant.signalToNoise = calculateSNR(variant);
            
            variants.add(variant);
        }
        
        // === STEP 3d: Create reference variant (for pileup mode) ===
        varsAtPosition.referenceVariant = createReferenceVariant(position);
        
        // === STEP 3e: Add variants to Vars ===
        varsAtPosition.variants = variants;
        
        // === STEP 3f: Store in output map ===
        alignedVariants.put(position, varsAtPosition);
    }
    
    // === STEP 4: Return wrapped in scope ===
    AlignedVarsData outputData = new AlignedVarsData(alignedVariants);
    return new Scope<>(scope, outputData);
}
```

---

## 4. IMPORTANT HELPER METHODS

### 4a. Statistical Calculation Helpers

#### `collectPositionsInReads(Variation v)`
```java
private List<Integer> collectPositionsInReads(Variation v) {
    // For each read containing this variation:
    //   - Find where in the read (1-based) the variant starts
    //   - Add to list
    // Returns: [15, 23, 45, 12, 67, ...] (positions in reads)
}
```

#### `collectBaseQualities(Variation v)`
```java
private List<Integer> collectBaseQualities(Variation v) {
    // For each read containing this variation:
    //   - Get the base quality at variant position
    //   - Convert from Phred to integer (typically 0-40)
    // Returns: [25, 30, 28, 35, 20, ...]
}
```

#### `collectMappingQualities(Variation v)`
```java
private List<Integer> collectMappingQualities(Variation v) {
    // For each read containing this variation:
    //   - Get the MAPQ (mapping quality)
    // Returns: [60, 60, 59, 60, 55, ...]
}
```

#### `hasAtLeast2Distinct(List<Integer> values)`
```java
private boolean hasAtLeast2Distinct(List<Integer> values) {
    // Check if list contains at least 2 different values
    // Example: [20, 20, 25, 20] → true (has 20 and 25)
    // Example: [20, 20, 20, 20] → false (only has 20)
    return values.stream().distinct().count() >= 2;
}
```

### 4b. Strand Bias Calculation

#### `calculateStrandBias(Variant v)`
```java
private String calculateStrandBias(Variant v) {
    int forward = v.varsCountOnForward;
    int reverse = v.varsCountOnReverse;
    int total = forward + reverse;
    
    if (total == 0) return "0";
    
    double ratio = Math.max(forward, reverse) / (double) total;
    
    // Strong bias: >90% on one strand
    if (ratio > 0.9) return "2";
    
    // Weak bias: 75-90% on one strand
    else if (ratio > 0.75) return "1";
    
    // Balanced
    else return "0";
}
```

**Alternate (Fisher Exact Test):**
- Some variants use Fisher's exact test for strand bias
- Compare variant strand distribution vs reference strand distribution
- P-value < threshold → bias flagged

### 4c. Microsatellite Detection

#### `detectMicrosatellite(int position, Reference ref)`
```java
private MicrosatelliteInfo detectMicrosatellite(int position, Reference ref) {
    // === Extract window around position ===
    int windowSize = 30;
    String window = ref.substring(position - windowSize/2, 
                                  position + windowSize/2);
    
    // === Detect repeating unit ===
    // Try units of length 1, 2, 3, 4...
    for (int unitLength = 1; unitLength <= 6; unitLength++) {
        String unit = window.substring(0, unitLength);
        
        // Count consecutive repeats
        int repeatCount = countConsecutiveRepeats(window, unit);
        
        if (repeatCount >= 3) {
            // Found repeat!
            double instability = (double) repeatCount / unitLength;
            return new MicrosatelliteInfo(instability, unitLength);
        }
    }
    
    return new MicrosatelliteInfo(0, 0);
}
```

**Examples:**
- `ACACACACAC` → unit="AC", length=2, repeats=5, MSI=2.5
- `AAAAAA` → unit="A", length=1, repeats=6, MSI=6.0
- `ACGT` → no repeat → MSI=0

### 4d. Shift3 Calculation (Deletions)

#### `calculateShift3(int position, Variant v, Reference ref)`
```java
private int calculateShift3(int position, Variant v, Reference ref) {
    // For deletions, calculate how far right we can shift
    // while maintaining same deletion
    // Example: "ACACAC" with del "-2" can shift 4bp right
    
    int shift = 0;
    int delLength = v.varallele.length() - 1;  // "-5" → 5
    
    String deletedSeq = ref.substring(position, position + delLength);
    
    // Try shifting right
    for (int i = 1; i <= 100; i++) {
        String nextSeq = ref.substring(position + i, position + i + delLength);
        if (nextSeq.equals(deletedSeq)) {
            shift = i;
        } else {
            break;
        }
    }
    
    return shift;
}
```

### 4e. High-Quality Read Filtering

#### `countHighQualityReads(Variation v, int qualThreshold)`
```java
private int countHighQualityReads(Variation v, int qualThreshold) {
    int count = 0;
    for (SAMRecord read : v.reads) {
        int baseQuality = getBaseQualityAt(read, v.positionInRead);
        if (baseQuality >= qualThreshold) {
            count++;
        }
    }
    return count;
}
```

### 4f. Allele Determination

#### `getRefAllele(int position, Variation v)`
```java
private String getRefAllele(int position, Variation v) {
    if (v.isSNP()) {
        return reference.substring(position, position + 1);
    } else if (v.isInsertion()) {
        return reference.substring(position, position + 1);  // Base before insertion
    } else if (v.isDeletion()) {
        int delLength = parseDeleteionLength(v.description);
        return reference.substring(position, position + delLength);
    } else {
        // Complex
        return parseComplexRef(v, reference);
    }
}
```

#### `getVarAllele(Variation v)`
```java
private String getVarAllele(Variation v) {
    if (v.isSNP()) {
        return v.description;  // "A", "T", "C", "G"
    } else if (v.isInsertion()) {
        return v.description.substring(1);  // "+ACGT" → "ACGT"
    } else if (v.isDeletion()) {
        return "";  // Or base before deletion
    } else {
        // Complex
        return parseComplexVar(v);
    }
}
```

---

## 5. DEPENDENCIES ON OTHER CLASSES

### Input Dependencies

#### `Scope<RealignedVariationData>`
```java
class RealignedVariationData {
    Map<Integer, Variation> variationMap;      // All variations found
    Map<Integer, Integer> refCoverage;         // Coverage at each position
    Map<Integer, Integer> insertionCoverage;   // Insertion-specific coverage
    // ... other fields from previous modules
}
```

#### `Variation` (Low-level variation)
```java
class Variation {
    int position;                    // Genomic position
    int endPosition;                 // End position (for indels)
    String description;              // "A" or "+ACGT" or "-5"
    int count;                       // Number of reads with this variation
    List<SAMRecord> reads;           // Reads containing this variation
    Map<Integer, Integer> posInReads; // Position in each read
    Map<Integer, Integer> qualities;  // Quality at each read
    // ... many other tracking fields
}
```

### Output Dependencies

#### `Scope<AlignedVarsData>`
```java
class AlignedVarsData {
    Map<Integer, Vars> alignedVariants;  // Position → Vars
}
```

#### `Vars` (Collection of variants at position)
```java
class Vars {
    List<Variant> variants;          // All variants at this position
    Variant referenceVariant;        // Reference "variant" (for pileup)
    StructuralVariant sv;            // SV flags (if any)
}
```

#### `Variant` (High-level variant with statistics)
```java
class Variant {
    // === IDENTITY === (from Variation)
    String descriptionString;        // "A" or "+ACGT" or "-5"
    String refallele;               // Reference allele
    String varallele;               // Variant allele
    String vartype;                 // "SNP", "Insertion", "Deletion", "Complex"
    int startPosition;
    int endPosition;
    
    // === COVERAGE === (calculated here)
    int positionCoverage;           // Total coverage
    int varsCountOnForward;         // Forward strand count
    int varsCountOnReverse;         // Reverse strand count
    double frequency;               // Allele frequency
    
    // === POSITION STATISTICS === (calculated here)
    double meanPosition;            // Average position in read
    double positionStdDev;          // Std dev of position
    boolean isAtLeastAt2Positions;  // Found at 2+ positions?
    
    // === QUALITY STATISTICS === (calculated here)
    double meanQuality;             // Average base quality
    double qualityStdDev;           // Std dev of quality
    boolean hasAtLeast2DiffQualities; // 2+ quality levels?
    
    // === MAPPING QUALITY === (calculated here)
    double meanMappingQuality;      // Average MAPQ
    
    // === STRAND BIAS === (calculated here)
    String strandBiasFlag;          // "0", "1", or "2"
    
    // === HIGH-QUALITY FILTERING === (calculated here)
    int hicnt;                      // High-quality read count
    int hicov;                      // High-quality coverage
    double highQualityReadsFrequency;
    double highQualityToLowQualityRatio;
    
    // === REALIGNMENT === (calculated here)
    double extraFrequency;          // Extra freq from realignment
    
    // === MICROSATELLITE === (calculated here)
    double msi;                     // MSI instability score
    int msint;                      // MSI unit length
    
    // === OTHER === (calculated here)
    int shift3;                     // 3' shift (deletions)
    String leftseq;                 // 20bp upstream
    String rightseq;                // 20bp downstream
    String genotype;                // "HOM", "HET", etc.
    boolean isNoise;                // Low quality flag
    double signalToNoise;           // SNR
    
    // === METHODS ===
    String varType();               // Determine variant type
    boolean isGoodVar(...);         // Quality filter
    void adjComplex();              // Adjust complex variants
    boolean isNoise();              // Check if noise
}
```

### Other Dependencies

#### `Configuration`
```java
class Configuration {
    double freqThreshold;           // Minimum frequency (e.g., 0.01)
    int qualThreshold;              // Quality threshold (e.g., 20)
    double strandBiasThreshold;     // Strand bias threshold
    boolean doPileup;               // Pileup mode flag
    // ... many other settings
}
```

#### `Reference` & `ReferenceResource`
```java
class Reference {
    String sequence;                // Reference sequence for region
    String substring(int start, int end);
}

class ReferenceResource {
    Reference getReference(Region region);
}
```

#### `Region`
```java
class Region {
    String chr;                     // Chromosome
    int start;                      // Start position
    int end;                        // End position
    String gene;                    // Gene name (optional)
}
```

---

## 6. APPROXIMATE LINE COUNT & COMPLEXITY

### Line Count Breakdown (Java)
```
Total: ~1059 lines

Breakdown:
- Class fields & initialization: ~50 lines
- accept() main method: ~150 lines
- Helper methods (statistics): ~300 lines
- Allele determination: ~100 lines
- Strand bias calculation: ~80 lines
- Microsatellite detection: ~120 lines
- High-quality filtering: ~100 lines
- Shift3 calculation: ~60 lines
- Reference variant creation: ~50 lines
- Utility methods: ~49 lines
```

### Complexity Analysis

#### Algorithmic Complexity
- **Time Complexity:** O(N × M × R)
  - N = number of genomic positions
  - M = number of variations per position (usually small, <10)
  - R = number of reads per variation (usually 10-1000)
- **Space Complexity:** O(N × V)
  - N = number of positions
  - V = number of variants (output)

#### Implementation Complexity: ⭐⭐⭐ HIGH

**Why Complex:**

1. **Statistical Calculations:**
   - Multiple distributions to track (position, quality, mapq)
   - Mean, standard deviation, distinct value counting
   - Need numerical stability for std dev calculations

2. **Strand Bias:**
   - Multiple approaches (ratio-based vs Fisher exact test)
   - Edge cases for low coverage

3. **Microsatellite Detection:**
   - Pattern matching for repeats
   - Multiple unit lengths to check
   - Window-based analysis

4. **High-Quality Filtering:**
   - Separate tracking for high vs low quality
   - Frequency recalculation
   - Ratio calculations

5. **Allele Determination:**
   - Different logic for SNPs, insertions, deletions, complex
   - Reference coordinate handling
   - Edge cases for complex variants

6. **Reference Integration:**
   - Need to extract flanking sequences
   - Coordinate arithmetic
   - Handle out-of-bounds cases

7. **Data Structure Juggling:**
   - Convert Variation → Variant
   - Multiple maps to manage
   - Preserve all tracking information

---

## 7. SIMPLE MODE SPECIFIC CONSIDERATIONS

### What's Different in Simple Mode?

1. **No Somatic Filtering:**
   - Simple mode doesn't compare tumor vs normal
   - All variants treated independently

2. **No Paired Sample Logic:**
   - Frequency calculation is straightforward
   - No germline vs somatic distinction

3. **Simpler Genotype Prediction:**
   ```
   if (frequency > 0.8) → "HOM" (homozygous)
   else if (frequency > 0.2) → "HET" (heterozygous)
   else → "REF" (too low)
   ```

4. **No Sample Comparison:**
   - No "tumor-only" vs "normal-only" variants
   - Simpler output format

### Simple Mode Flow Simplifications

```
Simple Mode:
  Variations → Calculate stats → Create Variants → Filter → Output

Somatic Mode (NOT implementing yet):
  Variations (tumor + normal) → Calculate stats for both 
  → Compare → Classify (germline/somatic) → Filter → Output
```

---

## 8. RUST PORT STRATEGY

### Recommended Rust Structure

```rust
// src/modules/statistics.rs (ToVarsBuilder equivalent)

use std::collections::{HashMap, BTreeMap};
use crate::variants::Variant;
use crate::data::{Variation, Vars};
use crate::scopedata::Scope;

pub struct ToVarsBuilder {
    // Configuration
    conf: Arc<Configuration>,
    
    // Reference access
    reference_resource: Arc<ReferenceResource>,
}

impl ToVarsBuilder {
    pub fn new(conf: Arc<Configuration>, 
               reference_resource: Arc<ReferenceResource>) -> Self {
        Self { conf, reference_resource }
    }
    
    pub fn accept(&self, scope: Scope<RealignedVariationData>) 
        -> Result<Scope<AlignedVarsData>> {
        
        // Main processing
        let mut aligned_variants = BTreeMap::new();
        
        for (position, variations) in scope.data.variation_map {
            let vars_at_pos = self.process_position(
                position, 
                variations,
                &scope
            )?;
            
            aligned_variants.insert(position, vars_at_pos);
        }
        
        Ok(Scope::new(scope, AlignedVarsData { aligned_variants }))
    }
    
    fn process_position(&self, 
                       position: i32,
                       variations: Vec<Variation>,
                       scope: &Scope<RealignedVariationData>) 
        -> Result<Vars> {
        
        let mut variants = Vec::new();
        
        for variation in variations {
            let variant = self.create_variant(
                position, 
                variation,
                scope
            )?;
            
            variants.push(variant);
        }
        
        let ref_variant = self.create_reference_variant(position, scope)?;
        
        Ok(Vars {
            variants,
            reference_variant: Some(ref_variant),
            sv: None,
        })
    }
    
    fn create_variant(&self,
                     position: i32,
                     variation: Variation,
                     scope: &Scope<RealignedVariationData>)
        -> Result<Variant> {
        
        let mut variant = Variant::default();
        
        // Set identity
        variant.description_string = variation.description.clone();
        variant.start_position = variation.position;
        variant.end_position = variation.end_position;
        
        // Determine alleles
        variant.refallele = self.get_ref_allele(position, &variation, scope)?;
        variant.varallele = self.get_var_allele(&variation)?;
        variant.vartype = variant.var_type();
        
        // Calculate coverage
        let cov = scope.data.ref_coverage.get(&position).unwrap_or(&0);
        variant.position_coverage = *cov;
        variant.vars_count_on_forward = self.count_forward_strand(&variation);
        variant.vars_count_on_reverse = self.count_reverse_strand(&variation);
        
        // Calculate frequency
        let total_var = variant.vars_count_on_forward + variant.vars_count_on_reverse;
        variant.frequency = total_var as f64 / variant.position_coverage as f64;
        
        // Position statistics
        let positions = self.collect_positions_in_reads(&variation);
        variant.mean_position = calculate_mean(&positions);
        variant.position_std_dev = calculate_std_dev(&positions);
        variant.is_at_least_at_2_positions = has_at_least_2_distinct(&positions);
        
        // Quality statistics
        let qualities = self.collect_base_qualities(&variation);
        variant.mean_quality = calculate_mean(&qualities);
        variant.quality_std_dev = calculate_std_dev(&qualities);
        variant.has_at_least_2_diff_qualities = has_at_least_2_distinct(&qualities);
        
        // Mapping quality
        let mapqs = self.collect_mapping_qualities(&variation);
        variant.mean_mapping_quality = calculate_mean(&mapqs);
        
        // Strand bias
        variant.strand_bias_flag = calculate_strand_bias(&variant);
        
        // High-quality filtering
        variant.hicnt = self.count_high_quality_reads(&variation, self.conf.qual_threshold);
        variant.hicov = scope.data.high_quality_coverage.get(&position).unwrap_or(&0);
        variant.high_quality_reads_frequency = 
            variant.hicnt as f64 / variant.hicov as f64;
        
        // Microsatellite
        let msi_info = detect_microsatellite(position, &scope.reference);
        variant.msi = msi_info.instability;
        variant.msint = msi_info.unit_length;
        
        // Shift3 (deletions only)
        if variant.vartype == "Deletion" {
            variant.shift3 = calculate_shift3(position, &variant, &scope.reference);
        }
        
        // Flanking sequences
        variant.leftseq = scope.reference.substring(position - 20, position)?;
        variant.rightseq = scope.reference.substring(position + 1, position + 21)?;
        
        // Genotype
        variant.genotype = predict_genotype(variant.frequency);
        
        Ok(variant)
    }
}

// === Helper Functions ===

fn calculate_mean(values: &[i32]) -> f64 {
    if values.is_empty() { return 0.0; }
    values.iter().sum::<i32>() as f64 / values.len() as f64
}

fn calculate_std_dev(values: &[i32]) -> f64 {
    if values.len() < 2 { return 0.0; }
    let mean = calculate_mean(values);
    let variance = values.iter()
        .map(|&x| {
            let diff = x as f64 - mean;
            diff * diff
        })
        .sum::<f64>() / values.len() as f64;
    variance.sqrt()
}

fn has_at_least_2_distinct(values: &[i32]) -> bool {
    let mut seen = std::collections::HashSet::new();
    for &val in values {
        seen.insert(val);
        if seen.len() >= 2 { return true; }
    }
    false
}

fn calculate_strand_bias(variant: &Variant) -> String {
    let forward = variant.vars_count_on_forward;
    let reverse = variant.vars_count_on_reverse;
    let total = forward + reverse;
    
    if total == 0 { return "0".to_string(); }
    
    let ratio = forward.max(reverse) as f64 / total as f64;
    
    if ratio > 0.9 { "2".to_string() }      // Strong bias
    else if ratio > 0.75 { "1".to_string() } // Weak bias
    else { "0".to_string() }                // Balanced
}

fn predict_genotype(frequency: f64) -> String {
    if frequency > 0.8 { "HOM".to_string() }
    else if frequency > 0.2 { "HET".to_string() }
    else { "REF".to_string() }
}

fn calculate_shift3(position: i32, variant: &Variant, reference: &Reference) -> i32 {
    // Implementation for deletion shift calculation
    // ... (see section 4d)
    0
}

fn detect_microsatellite(position: i32, reference: &Reference) -> MicrosatelliteInfo {
    // Implementation for MSI detection
    // ... (see section 4c)
    MicrosatelliteInfo { instability: 0.0, unit_length: 0 }
}

struct MicrosatelliteInfo {
    instability: f64,
    unit_length: i32,
}
```

### Key Rust Considerations

1. **Error Handling:**
   - Use `Result<T, E>` throughout
   - Define custom error types for statistics calculations
   - Handle out-of-bounds reference accesses

2. **Memory Management:**
   - Use `Arc<Configuration>` for shared config
   - Consider `Vec` vs `HashMap` for tracking lists
   - Efficient string handling (avoid unnecessary clones)

3. **Numerical Stability:**
   - Use f64 for all statistical calculations
   - Watch for division by zero
   - Consider using a statistics crate (e.g., `statrs`)

4. **Iterator Patterns:**
   - Use iterators for mean/std dev calculations
   - Consider `rayon` for parallel processing (later optimization)

5. **Testing Strategy:**
   - Unit test each statistical function independently
   - Property tests for mean/std dev
   - Integration test with known VarDict output

---

## 9. IMPLEMENTATION CHECKLIST

### Phase 1: Basic Structure
- [ ] Define `ToVarsBuilder` struct
- [ ] Implement `accept()` method skeleton
- [ ] Set up iteration over positions and variations

### Phase 2: Basic Statistics
- [ ] Implement coverage counting
- [ ] Implement frequency calculation
- [ ] Implement position-in-read statistics (mean, std dev, distinct)
- [ ] Implement base quality statistics
- [ ] Implement mapping quality statistics

### Phase 3: Advanced Statistics
- [ ] Implement strand bias calculation
- [ ] Implement high-quality filtering
- [ ] Implement microsatellite detection
- [ ] Implement shift3 calculation (deletions)

### Phase 4: Allele Determination
- [ ] Implement ref/var allele extraction for SNPs
- [ ] Implement ref/var allele extraction for insertions
- [ ] Implement ref/var allele extraction for deletions
- [ ] Implement ref/var allele extraction for complex variants

### Phase 5: Output Creation
- [ ] Create Variant objects with all fields
- [ ] Create Vars containers
- [ ] Create reference variants (pileup mode)
- [ ] Wrap in Scope for pipeline

### Phase 6: Testing
- [ ] Unit tests for statistical functions
- [ ] Unit tests for allele determination
- [ ] Integration test with real VarDict data
- [ ] Compare output with Java VarDict

### Phase 7: Optimization (Later)
- [ ] Profile hotspots
- [ ] Consider parallel processing
- [ ] Optimize memory allocations

---

## 10. EDGE CASES TO HANDLE

1. **Zero Coverage:**
   - Division by zero when calculating frequency
   - Handle gracefully (frequency = 0.0)

2. **Single Read:**
   - Cannot calculate std dev with n=1
   - Set std dev = 0.0

3. **No Distinct Values:**
   - `isAtLeastAt2Positions` = false
   - `hasAtLeast2DiffQualities` = false

4. **Out of Bounds Reference:**
   - Flanking sequences near chromosome edges
   - Truncate or pad with 'N'

5. **Complex Variants:**
   - Mixed SNP + indel
   - May require special allele parsing

6. **High MSI Regions:**
   - Very long repeats (>50bp)
   - Cap MSI score to prevent overflow

7. **Strand Bias Edge Cases:**
   - All reads on one strand (legitimate for some contexts)
   - Very low coverage (<3 reads)

8. **Deletions at Boundaries:**
   - Deletion extends past region end
   - Shift3 calculation hits end of reference

---

## 11. DEPENDENCIES WITH OTHER MODULES

### Upstream (Inputs)
```
CigarParser → Variations (SNPs, indels)
     ↓
VariationRealigner → Additional variations from soft-clips
     ↓
ToVarsBuilder ← YOU ARE HERE
```

**Requirements from Upstream:**
- `Variation` objects must have:
  - Position information
  - Description string
  - List of reads containing variation
  - Position in each read
  - Quality at each read
- Coverage maps must be complete

### Downstream (Outputs)
```
ToVarsBuilder → Variant objects with statistics
     ↓
SimplePostProcessModule → Quality filtering
     ↓
VariantPrinter → Output formatting
```

**Requirements for Downstream:**
- `Variant` objects must have ALL statistics calculated
- `isGoodVar()` must be implementable (depends on statistics)
- Alleles must be properly formatted for VCF output

---

## 12. TESTING STRATEGY

### Unit Tests

```rust
#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_calculate_mean() {
        assert_eq!(calculate_mean(&[10, 20, 30]), 20.0);
        assert_eq!(calculate_mean(&[]), 0.0);
    }
    
    #[test]
    fn test_calculate_std_dev() {
        let values = vec![10, 20, 30];
        let std_dev = calculate_std_dev(&values);
        assert!((std_dev - 8.165).abs() < 0.01);
    }
    
    #[test]
    fn test_has_at_least_2_distinct() {
        assert!(has_at_least_2_distinct(&[10, 10, 20, 10]));
        assert!(!has_at_least_2_distinct(&[10, 10, 10]));
    }
    
    #[test]
    fn test_strand_bias() {
        let mut variant = Variant::default();
        variant.vars_count_on_forward = 9;
        variant.vars_count_on_reverse = 1;
        assert_eq!(calculate_strand_bias(&variant), "2");  // Strong bias
        
        variant.vars_count_on_forward = 5;
        variant.vars_count_on_reverse = 5;
        assert_eq!(calculate_strand_bias(&variant), "0");  // Balanced
    }
}
```

### Integration Tests

```rust
#[test]
fn test_tovarbuilder_integration() {
    // Create test data
    let variations = create_test_variations();
    let scope = create_test_scope(variations);
    
    // Run ToVarsBuilder
    let builder = ToVarsBuilder::new(config, ref_resource);
    let output = builder.accept(scope).unwrap();
    
    // Verify output
    assert_eq!(output.data.aligned_variants.len(), 10);
    
    let variant = &output.data.aligned_variants[&100].variants[0];
    assert!((variant.frequency - 0.25).abs() < 0.01);
    assert_eq!(variant.vartype, "SNP");
    assert!(variant.is_at_least_at_2_positions);
}
```

### Comparison Tests

```bash
# Run Java VarDict on test BAM
java -jar VarDict.jar -b test.bam -c 1 -S 2 -E 3 -g 4 test.bed > java_output.txt

# Run Rust VarDict on same BAM
cargo run -- -b test.bam -c 1 -S 2 -E 3 -g 4 test.bed > rust_output.txt

# Compare outputs (allow for minor float differences)
diff -u java_output.txt rust_output.txt
```

---

## 13. PERFORMANCE CONSIDERATIONS

### Expected Performance

**Java VarDict:** ~1000 regions/second (single-threaded)

**Rust Target:** Similar or better (with optimizations)

### Optimization Opportunities (Later)

1. **Parallel Processing:**
   - Process positions in parallel (use `rayon`)
   - Each position is independent

2. **Memory Pools:**
   - Reuse `Vec` allocations for tracking lists
   - Avoid repeated allocations

3. **SIMD:**
   - Consider SIMD for mean/std dev calculations (if bottleneck)

4. **Lazy Calculation:**
   - Only calculate statistics that will be used
   - Skip expensive calculations for low-frequency variants

---

## 14. SUMMARY

**ToVarsBuilder** is the **core statistics engine** of VarDict. It converts low-level `Variation` objects (simple position + description) into high-level `Variant` objects with complete statistical analysis.

**Key Responsibilities:**
1. ✅ Calculate allele frequencies
2. ✅ Analyze position distribution in reads
3. ✅ Analyze base quality distribution
4. ✅ Analyze mapping quality
5. ✅ Detect strand bias
6. ✅ Filter by high-quality reads
7. ✅ Detect microsatellite instability
8. ✅ Calculate shift3 for deletions
9. ✅ Extract flanking sequences
10. ✅ Predict genotypes

**Complexity:** HIGH (⭐⭐⭐)
- ~1059 lines in Java
- Multiple statistical calculations
- Multiple data transformations
- Critical for downstream filtering

**Simple Mode:** Simpler than somatic mode (no paired sample logic)

**Rust Port:** Straightforward translation with attention to:
- Numerical stability
- Error handling
- Iterator patterns
- Memory efficiency

---

**Ready to implement!** Start with basic statistics (mean, frequency) and build up to advanced features (MSI, shift3).

**Created:** January 10, 2026  
**Status:** Ready for Rust implementation ✓
