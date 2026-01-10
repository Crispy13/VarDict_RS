# ToVarsBuilder - Preparation for Rust Port

**Module:** `ToVarsBuilder.java` → `to_vars_builder.rs`  
**Complexity:** ⭐⭐⭐ HIGH (1059 lines in Java, ~600-700 for Simple Mode only)  
**Priority:** CRITICAL - Core statistics and variant creation module  
**Status:** Ready for implementation  
**Mode Focus:** **SIMPLE MODE ONLY** - No Somatic, Amplicon, or paired-sample logic  

---

## Overview

ToVarsBuilder is the **most complex module** in the VarDict pipeline. It converts low-level `Variation` objects (from CigarParser/VariationRealigner) into high-level `Variant` objects with complete statistical analysis.

**IMPORTANT:** We only implement **Simple Mode**, which means:
- ✅ Single-sample variant calling
- ✅ Basic statistics (position, quality, strand bias)
- ✅ Simple ratio-based strand bias (NO Fisher's exact test)
- ❌ NO somatic filtering
- ❌ NO paired-sample comparison
- ❌ NO tumor-normal logic
- ❌ NO complex genotype prediction

### Responsibilities

1. **Group variations by genomic position**
2. **Calculate 12+ statistics** for each variant (Simple Mode subset):
   - Allele frequencies (forward/reverse strands) ✅
   - Position-in-read distribution ✅
   - Base quality distribution ✅
   - Mapping quality metrics ✅
   - **Simple ratio-based strand bias** ✅ (NOT Fisher's exact)
   - High-quality read filtering ✅ (basic Q20 only)
   - Microsatellite instability (MSI) detection ✅
   - Shift3 calculation for deletions ✅
3. **Create Variant objects** with all fields populated
4. **Group variants** into `Vars` collections per position

**NOT IMPLEMENTED (Somatic/Amplicon only):**
- ❌ Tumor-normal filtering
- ❌ Fisher's exact test for strand bias
- ❌ Complex genotype prediction
- ❌ Complex high-quality ratio logic
- ❌ Paired-sample comparison

---

## Input/Output

### Input
**Type:** `Scope<RealignedVariationData>`

**Contains:**
- `nonInsertionVariants`: Map<position, Map<variantKey, Variation>>
- `insertionVariants`: Map<position, Map<variantKey, Variation>>  
- `refCoverage`: Map<position, coverage_count>
- Reference sequence
- Configuration settings

### Output
**Type:** `Scope<AlignedVarsData>`

**Contains:**
- `alignedVariants`: Map<position, Vars>
  - Each `Vars` contains:
    - `variants`: Vec<Variant> - All variants at this position
    - `referenceVariant`: Option<Variant> - Reference call (pileup mode)
    - `sv`: StructuralVariantFlags - SV information

---

## Main Algorithm Flow

```rust
fn accept(&mut self, scope: Scope<RealignedVariationData>) -> Scope<AlignedVarsData> {
    // 1. Extract all variations from maps
    let all_variations = combine_variations(
        scope.non_insertion_variants,
        scope.insertion_variants
    );
    
    // 2. Group by position
    let grouped: HashMap<i64, Vec<Variation>> = group_by_position(all_variations);
    
    // 3. For each position, calculate statistics
    let mut aligned_variants = HashMap::new();
    for (position, variations) in grouped {
        // 3a. Group variations by variant key (same type at same position)
        let variant_groups = group_by_variant_key(variations);
        
        // 3b. Calculate statistics for each variant
        let mut variants_at_pos = Vec::new();
        for (key, var_list) in variant_groups {
            let variant = calculate_variant_statistics(
                var_list,
                position,
                scope.ref_coverage.get(&position),
                &scope.reference,
                &scope.config
            );
            variants_at_pos.push(variant);
        }
        
        // 3c. Create Vars object
        let vars = Vars {
            variants: variants_at_pos,
            reference_variant: if config.do_pileup { Some(create_ref_variant()) } else { None },
            sv: StructuralVariantFlags::default(),
        };
        
        aligned_variants.insert(position, vars);
    }
    
    // 4. Return new scope
    Scope {
        data: AlignedVarsData { aligned_variants },
        ..scope
    }
}
```

---

## Statistics Calculation (16 Metrics)

For each variant, calculate:

### 1. **Strand Counts**
```rust
struct StrandCounts {
    forward: usize,  // varsCountOnForward
    reverse: usize,  // varsCountOnReverse
}
```
Count how many reads with variant are on forward vs reverse strand.

### 2. **Allele Frequency**
```rust
frequency = (forward_count + reverse_count) as f64 / total_coverage as f64
```

### 3. **Position-in-Read Analysis**
```rust
struct PositionStats {
    mean_position: f64,           // Average position in read where variant appears
    position_std_dev: f64,         // Standard deviation
    is_at_least_2_positions: bool, // Found at 2+ distinct positions
}
```
- Extract position from each read
- Calculate mean and std deviation
- Check if found at 2+ different positions (quality flag)

### 4. **Base Quality Analysis**
```rust
struct QualityStats {
    mean_quality: f64,              // Average base quality (Q-score)
    quality_std_dev: f64,            // Standard deviation
    has_at_least_2_diff_qualities: bool, // Found at 2+ quality levels
}
```
- Extract base quality from each read
- Calculate mean and std deviation
- Check diversity (prevents single-quality artifacts)

### 5. **Mapping Quality**
```rust
mean_mapping_quality: f64  // Average MAPQ of reads with variant
```

### 6. **Strand Bias Detection**
```rust
enum StrandBiasFlag {
    NoBias = 0,
    WeakBias = 1,
    StrongBias = 2,
}
```
**Simple Mode Uses Simple Ratio-Based Method (NO statistical tests):**
- Calculate ratio of forward/reverse counts
- If all variants on one strand → strong bias
- If >80% on one strand → weak bias
- Otherwise → no bias

```rust
fn check_strand_bias(forward: usize, reverse: usize) -> StrandBiasFlag {
    if forward == 0 || reverse == 0 {
        StrandBiasFlag::StrongBias  // All on one strand
    } else {
        let ratio = forward.max(reverse) as f64 / forward.min(reverse) as f64;
        if ratio > 10.0 {
            StrandBiasFlag::StrongBias   // >90% on one strand
        } else if ratio > 4.0 {
            StrandBiasFlag::WeakBias     // >80% on one strand
        } else {
            StrandBiasFlag::NoBias       // Balanced
        }
    }
}
```

### 7. **High-Quality Read Filtering** (Simple Mode - BASIC)
```rust
struct HighQualityStats {
    high_quality_count: usize,           // Reads with Q ≥ 20 (threshold)
    high_quality_frequency: f64,          // Frequency in high-Q reads only
}
```

**Simple Mode:** Basic Q20 thresholding only (no somatic logic)

### 8. **MSI (Microsatellite Instability) Detection**
```rust
struct MsiInfo {
    msi: f64,      // MSI score (length of repeat)
    msint: f64,    // MSI interval
}
```
Detect if variant is in repetitive sequence region.

### 9. **Shift3 Calculation**
```rust
shift3: i32  // How many bases can deletion shift 3' (right-normalization)
```
For deletions, calculate how far it can be shifted right while maintaining same deletion.

### 10. **Reference/Variant Alleles**
```rust
struct AlleleInfo {
    refallele: String,   // Reference allele sequence
    varallele: String,   // Variant allele sequence  
    vartype: VarType,    // SNP, Insertion, Deletion, Complex
}
```

### 11. **Coverage Metrics**
```rust
struct CoverageInfo {
    position_coverage: usize,     // Total depth at position
    forward_coverage: usize,       // Reads on forward strand
    reverse_coverage: usize,       // Reads on reverse strand
}
```

### 12. **Description String**
```rust
description_string: String
```
Format depends on variant type:
- **SNP:** Single letter ("A", "T", "G", "C")
- **Insertion:** "+ACGT"
- **Deletion:** "-5" (number of bases)
- **Complex:** "ACGT#-3" (insertion + deletion)

### 13. **Genotype Prediction** (Simple Mode - SIMPLIFIED)
```rust
genotype: String  // "0/0", "0/1", "1/1"
```
**Simple Mode Logic (NO complex heuristics):**
- frequency = 0 → "0/0" (homozygous reference)
- 0 < frequency < 0.5 → "0/1" (heterozygous)
- frequency >= 0.5 → "1/1" (homozygous alternate)

**NOT USED:** NO tumor-normal filtering, NO somatic logic, NO complex prediction

### 14. **Left/Right Sequence Context**
```rust
leftseq: String,   // 20bp upstream context
rightseq: String,  // 20bp downstream context
```

### 15. **NM (Edit Distance)**
```rust
nm: f64  // Average edit distance from reference
```

### 16. **Start/End Positions**
```rust
start_position: i64,
end_position: i64,
```

---

## Key Helper Functions

### 1. `group_by_position`
Groups variations by genomic position.

### 2. `group_by_variant_key`  
Groups variations at same position by variant identity.

### 3. `calculate_mean_and_std`
```rust
fn calculate_mean_and_std(values: &[f64]) -> (f64, f64) {
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    let variance = values.iter()
        .map(|v| (v - mean).powi(2))
        .sum::<f64>() / values.len() as f64;
    (mean, variance.sqrt())
}
```

### 4. `check_strand_bias`
```rust
fn check_strand_bias(forward: usize, reverse: usize) -> StrandBiasFlag {
    if forward == 0 || reverse == 0 {
        StrandBiasFlag::StrongBias
    } else {
        let ratio = forward.max(reverse) as f64 / forward.min(reverse) as f64;
        if ratio > 10.0 {
            StrandBiasFlag::StrongBias
        } else if ratio > 3.0 {
            StrandBiasFlag::WeakBias
        } else {
            StrandBiasFlag::NoBias
        }
    }
}
```

### 5. `is_at_least_2_positions`
```rust
fn is_at_least_2_positions(positions: &[i32]) -> bool {
    let unique: HashSet<i32> = positions.iter().copied().collect();
    unique.len() >= 2
}
```

### 6. `detect_msi`
```rust
fn detect_msi(sequence: &str, position: usize) -> (f64, f64) {
    // Find longest tandem repeat around position
    // Return (repeat_length, repeat_count)
}
```

### 7. `calculate_shift3`
```rust
fn calculate_shift3(
    ref_seq: &[u8],
    position: usize,
    deletion_len: usize
) -> i32 {
    // For deletions, how many bases can we shift right?
    let mut shift = 0;
    let mut pos = position;
    while pos + deletion_len < ref_seq.len() {
        if ref_seq[pos] == ref_seq[pos + deletion_len] {
            shift += 1;
            pos += 1;
        } else {
            break;
        }
    }
    shift
}
```

### 8. `create_description_string`
```rust
fn create_description_string(var_type: &VarType, sequence: &str) -> String {
    match var_type {
        VarType::SNV(base) => base.to_string(),
        VarType::Insertion(seq) => format!("+{}", seq),
        VarType::Deletion(len) => format!("-{}", len),
        VarType::Complex { ins, del } => format!("{}#-{}", ins, del),
    }
}
```

---

## Data Structures Needed

### Rust Equivalents

```rust
// Core variant structure
pub struct Variant {
    // Identity
    pub description_string: String,
    pub refallele: String,
    pub varallele: String,
    pub vartype: VarType,
    
    // Position
    pub start_position: i64,
    pub end_position: i64,
    
    // Counts
    pub vars_count_on_forward: usize,
    pub vars_count_on_reverse: usize,
    pub position_coverage: usize,
    
    // Frequencies
    pub frequency: f64,
    pub high_quality_reads_frequency: f64,
    
    // Quality metrics
    pub mean_position: f64,
    pub mean_quality: f64,
    pub mean_mapping_quality: f64,
    pub high_quality_to_low_quality_ratio: f64,
    
    // Flags
    pub strand_bias_flag: StrandBiasFlag,
    pub is_at_least_at_2_positions: bool,
    pub has_at_least_2_diff_qualities: bool,
    
    // Context
    pub leftseq: String,
    pub rightseq: String,
    
    // Special features
    pub msi: f64,
    pub msint: f64,
    pub shift3: i32,
    pub nm: f64,
    
    // Genotype
    pub genotype: String,
    
    // Additional (for filtering)
    pub pstd: bool,  // Position standard deviation flag
    pub qstd: bool,  // Quality standard deviation flag
}

pub enum VarType {
    SNV(char),
    Insertion(String),
    Deletion(usize),
    Complex { insertion: String, deletion: usize },
}

pub enum StrandBiasFlag {
    NoBias = 0,
    WeakBias = 1,
    StrongBias = 2,
}

pub struct Vars {
    pub variants: Vec<Variant>,
    pub reference_variant: Option<Variant>,
    pub sv: StructuralVariantFlags,
}

#[derive(Default)]
pub struct StructuralVariantFlags {
    pub splits: usize,
    pub pairs: usize,
    pub clusters: usize,
}
```

---

## Implementation Strategy

### Phase 1: Basic Structure (Day 1)
- [ ] Create `to_vars_builder.rs` file
- [ ] Define ToVarsBuilder struct
- [ ] Define Variant struct
- [ ] Define Vars struct
- [ ] Define VarType enum
- [ ] Implement basic `accept()` method skeleton

### Phase 2: Statistics Functions (Day 2-3)
- [ ] Implement `calculate_mean_and_std`
- [ ] Implement position-in-read analysis
- [ ] Implement base quality analysis
- [ ] Implement mapping quality calculation
- [ ] Implement strand bias detection
- [ ] Implement high-quality filtering
- [ ] Add unit tests for each

### Phase 3: Variant Creation (Day 4)
- [ ] Implement `create_description_string`
- [ ] Implement allele determination
- [ ] Implement MSI detection
- [ ] Implement shift3 calculation
- [ ] Implement context sequence extraction
- [ ] Create complete Variant objects

### Phase 4: Integration (Day 5)
- [ ] Integrate with existing pipeline
- [ ] Test with real BAM data
- [ ] Compare output with Java VarDict
- [ ] Fix any discrepancies

### Phase 5: Optimization (Day 6+)
- [ ] Optimize hash maps
- [ ] Reduce allocations
- [ ] Parallel processing for positions
- [ ] Profile and optimize hot paths

---

## Edge Cases to Handle

1. **Zero coverage at position** - Skip or handle gracefully
2. **Single read with variant** - May fail "at 2+ positions" check
3. **All reads on one strand** - Strong strand bias
4. **Variant at read ends only** - Quality flag should catch
5. **MSI regions** - Adjust frequency thresholds
6. **Complex indels** - Proper description string formatting
7. **Reference sequence boundaries** - Don't go out of bounds for context
8. **Empty variation lists** - Skip position

---

## Testing Strategy

### Unit Tests
```rust
#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_mean_and_std() {
        let values = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let (mean, std) = calculate_mean_and_std(&values);
        assert_eq!(mean, 3.0);
        assert!((std - 1.414).abs() < 0.01);
    }
    
    #[test]
    fn test_strand_bias_detection() {
        assert_eq!(check_strand_bias(10, 0), StrandBiasFlag::StrongBias);
        assert_eq!(check_strand_bias(10, 10), StrandBiasFlag::NoBias);
        assert_eq!(check_strand_bias(100, 5), StrandBiasFlag::StrongBias);
    }
    
    #[test]
    fn test_is_at_least_2_positions() {
        assert!(is_at_least_2_positions(&[10, 20, 30]));
        assert!(!is_at_least_2_positions(&[10, 10, 10]));
    }
    
    #[test]
    fn test_description_string_snv() {
        let desc = create_description_string(&VarType::SNV('A'), "");
        assert_eq!(desc, "A");
    }
    
    #[test]
    fn test_description_string_insertion() {
        let desc = create_description_string(&VarType::Insertion("ACGT".to_string()), "");
        assert_eq!(desc, "+ACGT");
    }
}
```

### Integration Tests
```rust
#[test]
fn test_simple_snv() {
    // Create mock variations with SNV
    // Run ToVarsBuilder
    // Check output Variant has correct statistics
}

#[test]
fn test_strand_bias_flagging() {
    // Create variations all on forward strand
    // Check StrandBiasFlag is StrongBias
}
```

---

## Dependencies

### Crate-Internal
- `variants::Variant` - Output structure
- `variants::Vars` - Collection structure
- `data::Variation` - Input structure
- `conf::Configuration` - Settings
- `data::Reference` - Reference sequence
- `data::Region` - Genomic region

### External
- `std::collections::HashMap` - Grouping operations
- `std::collections::HashSet` - Uniqueness checks

---

## Performance Considerations

### Time Complexity
**O(N × M × R)** where:
- N = number of positions with variants (~thousands)
- M = variants per position (~1-5 typically)
- R = reads per variant (~10-100 typically)

**Estimated:** ~10-50ms per region for typical data

### Space Complexity
**O(N × M)** for storing all variants

### Optimizations
1. **Reuse allocations** - Pre-allocate Vec capacity
2. **Parallel processing** - Process positions independently
3. **Lazy statistics** - Only calculate when needed
4. **Cache reference sequence** - Don't re-fetch for every position

---

## Comparison with Java

### Java (1059 lines)
- Uses streams and functional programming
- Some code duplication for somatic vs simple mode
- Complex nesting of loops

### Rust Target (~800-1000 lines)
- More explicit control flow
- Better type safety for variant types
- Potential for zero-copy optimizations
- Simpler without somatic mode

---

## Next Steps After Completion

1. **Complete ToVarsBuilder** ← You are here
2. **Post-processing module** (SimplePostProcessModule)
3. **VCF output formatter**
4. **End-to-end integration test**
5. **Benchmarking vs Java version**

---

**Estimated Effort:** 3-5 days for Simple Mode (No somatic/amplicon complexity)  
**Priority:** CRITICAL PATH - blocks output generation  
**Complexity:** ⭐⭐ MEDIUM (simplified without somatic logic)  
**Lines of Code:** ~600-700 in Rust (vs 1059 in Java with all modes)
