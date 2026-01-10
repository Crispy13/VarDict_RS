# VarDict SimpleMode - Quick Reference Card

## 🚀 TL;DR - The Essentials

### What is VarDict?
Variant caller that finds SNPs and indels from BAM files. 10x faster Java version of original Perl implementation.

### Simple Mode Pipeline (5 steps)
```
Parse BAM → Filter reads → Parse CIGAR (SNPs/indels) 
         → Realign soft clips → Calculate stats 
         → Filter quality → Output VCF
```

### Most Important Classes
1. **Variant.java** - Data structure (30+ fields)
2. **ToVarsBuilder.java** - Statistics calculator (1059 lines)
3. **VariationRealigner.java** - Soft-clip realignment
4. **SimplePostProcessModule.java** - Quality filter

---

## 📋 Variant Structure (Most Important!)

```rust
pub struct Variant {
    // Identity
    pub description_string: String,  // "A", "+ACGT", "-5"
    pub ref_allele: String,          // Reference
    pub var_allele: String,          // Variant
    pub var_type: String,            // SNP, Insertion, Deletion, Complex
    
    // Position
    pub start_position: i32,
    pub end_position: i32,
    
    // Coverage
    pub position_coverage: i32,
    pub vars_count_on_forward: i32,
    pub vars_count_on_reverse: i32,
    
    // Frequency
    pub frequency: f64,              // AF
    pub high_quality_reads_frequency: f64,
    pub extra_frequency: f64,
    
    // Quality metrics ⭐⭐⭐
    pub mean_position: f64,          // Position in read
    pub is_at_least_at_2_positions: bool,  // KEY: 2+ positions
    pub mean_quality: f64,           // Base quality
    pub has_at_least_2_diff_qualities: bool,  // KEY: 2+ qualities
    pub mean_mapping_quality: f64,   // MAPQ
    
    // Strand bias
    pub strand_bias_flag: String,    // "0", "1", "2"
    pub high_quality_to_low_quality_ratio: f64,
    
    // Microsatellite
    pub msi: f64,                    // Instability score
    pub msint: i32,                  // Unit length
    
    // Context
    pub left_seq: String,            // 20bp upstream
    pub right_seq: String,           // 20bp downstream
    
    // Indel-specific
    pub shift3: i32,                 // 3' shift
    
    // Other
    pub hicnt: i32,                  // High-quality read count
    pub hicov: i32,                  // High-quality coverage
    pub duprate: f64,                // Duplication rate
    pub crispr: i32,                 // CRISPR distance
}

impl Variant {
    pub fn var_type(&self) -> String {
        // "SNP" if single letter
        // "Insertion" if starts with "+"
        // "Deletion" if starts with "-"
        // "Complex" if has "#" or "&"
    }
    
    pub fn is_good_var(&self) -> bool {
        // Check: 2+ positions, 2+ qualities, strand bias
        // This is the KEY quality filter!
    }
}
```

---

## 🔄 Pipeline Data Flow

```
Input: BAM file, Region(chr, start, end), Reference

Step 1: SAMFileParser
  Input:  Region + BAM
  Output: List<SAMRecord> (reads)
  
Step 2: RecordPreprocessor  
  Input:  SAMRecords
  Output: Filtered reads + coverage map
  
Step 3: CigarParser
  Input:  SAMRecords with CIGAR strings
  Output: Variations (SNPs, indels) grouped by position
  
Step 4: VariationRealigner ⭐
  Input:  Variations + soft-clipped sequences
  Output: Additional variations from hidden indels
  
Step 5: ToVarsBuilder ⭐⭐⭐
  Input:  All variations
  Output: Variant objects with:
          - Allele frequency
          - Position in read distribution
          - Quality distribution
          - Strand bias
          - Mapping quality
          - MSI detection
  
Step 6: SimplePostProcessModule
  Input:  Variant objects
  Output: Filtered variants (isGoodVar = true)
  
Step 7: VariantPrinter
  Output: VCF/text format
```

---

## 🎯 Key Algorithms (Implement in This Order)

### 1. CIGAR Parser (Medium)
Parse CIGAR string like "10M2I3M1D5M"
- M = alignment match (can be SNP)
- I = insertion
- D = deletion  
- S = soft clip (save for later)

**Output:** For each position, list of variations

### 2. Statistics Calculator (Hard) ⭐⭐⭐
For each variant, calculate:
```
List<Integer> positions_in_reads = []
List<Integer> qualities_in_reads = []

for each read with variant:
    positions_in_reads.add(position_in_read)
    qualities_in_reads.add(base_quality)

variant.mean_position = average(positions_in_reads)
variant.is_at_least_at_2_positions = unique_count(positions_in_reads) >= 2
variant.mean_quality = average(qualities_in_reads)
variant.has_at_least_2_diff_qualities = unique_count(qualities_in_reads) >= 2

variant.frequency = variant_count / total_coverage
variant.strand_bias_flag = calculate_bias(forward_count, reverse_count)
```

### 3. Quality Filter (Medium)
```
isGoodVar():
  if strand_bias_flag == "2" and not_in_special_region:
    return false
  if !is_at_least_at_2_positions:
    return false
  if !has_at_least_2_diff_qualities:
    return false
  return true
```

### 4. Soft-Clip Realignment (Hard) ⭐⭐
```
for each read with soft clip:
    extract_soft_clipped_sequence()
    local_alignment = smith_waterman(clipped_seq, reference)
    if alignment_explains_clipping:
        create_hidden_indel_variation()
```

---

## 📊 Output Format

### Header (Tab-separated)
```
Sample Gene Chr Start End Ref Alt Depth AltDepth RefFwdReads RefRevReads 
AltFwdReads AltRevReads Genotype AF Bias PMean PStd QMean QStd MQ 
Sig_Noise HiAF ExtraAF shift3 MSI MSI_NT NM HiCnt HiCov 5pFlankSeq 
3pFlankSeq Seg VarType Duprate SV_info [CRISPR]
```

### Example Row
```
sample  EGFR  chr7  55086707  55086707  G  A  150  45  20  25  22  23  
HET  0.30  0  45.5  12.3  28.4  3.2  55  SNP_NOISE  0.32  0.01  0  
1.2  2  0  0.3  40  120  ACGTACGTAC...  ACGTACGTAC...  chr7:55000000-55200000  
SNP  0.02  .
```

---

## 🔑 Critical Implementation Details

### Variant Description String Format
- `"A"` = SNP (single letter)
- `"+ACG"` = Insertion of ACG
- `"-5"` = Deletion of 5 bases
- `"AC#-3"` = Insertion of AC followed by deletion of 3 bases (complex)
- `"A^G"` = SNP followed by insertion
- `"A^2"` = SNP followed by deletion of 2 bases

### isGoodVar() - The Gate Keeper
This method is **THE** quality filter. Must implement exactly:
1. Check strand bias (too strong = reject)
2. Check `is_at_least_at_2_positions` (false = reject)
3. Check `has_at_least_2_diff_qualities` (false = reject)
4. Check MSI region (relax thresholds if true)
5. Check frequency (must be >= min_freq)

### Statistics Keys
- **is_at_least_at_2_positions:** Variant found at position 15, 20, 25 in different reads = TRUE. Found only at position 10 in all reads = FALSE.
- **has_at_least_2_diff_qualities:** Variant found with Q20 in some reads, Q25 in others = TRUE. Always Q20 = FALSE.
- **strand_bias_flag:** "0" = balanced, "1" = slight bias, "2" = strong bias

---

## 🚦 Configuration Parameters

```
min_freq: f64              // Min allele frequency (default: 0.02)
min_quality: i32           // Min base quality (default: 0)
min_mapping_quality: i32   // Min MAPQ (default: 0)
min_alt_count: i32         // Min variant count (default: 2)

do_pileup: bool            // Output reference calls (default: false)
do_strand_bias: bool       // Strict strand bias filter (default: false)
crispr_cutting_site: i32   // CRISPR cut position (default: 0 = disabled)

print_header: bool         // Output header (default: true)
```

---

## 📁 Rust Port File Structure (Suggested)

```
src/
├── main.rs                      # Entry point
├── lib.rs                       # Library exports
├── config.rs                    # Configuration
├── pipeline.rs                  # Main pipeline
├── variant.rs                   # Variant struct ⭐
├── parser/
│   ├── mod.rs
│   ├── bam.rs                   # BAM parsing
│   ├── cigar.rs                 # CIGAR parsing ⭐
│   └── record.rs                # SAM record filtering
├── modules/
│   ├── mod.rs
│   ├── statistics.rs            # ToVarsBuilder equivalent ⭐⭐⭐
│   ├── realigner.rs             # VariationRealigner ⭐
│   └── filter.rs                # SimplePostProcessModule
├── data/
│   ├── mod.rs
│   ├── variation.rs             # Low-level variation
│   └── scope.rs                 # Pipeline scope
└── output/
    ├── mod.rs
    └── formatter.rs             # VCF formatting
```

---

## ✅ Testing Strategy

### Unit Tests
1. CIGAR parser (test all operations: M, I, D, S, N)
2. Statistics calculator (verify position/quality calculations)
3. Soft-clip realigner (test hidden indel detection)
4. Quality filter (verify acceptance/rejection)

### Integration Tests
1. Run on small BAM file
2. Compare output with Java VarDict (position by position)
3. Verify statistics match exactly
4. Test parallel region processing

---

## 📚 Documentation Cross-Reference

- **Overview:** EXAMINATION_SUMMARY.md
- **Architecture:** VARDICTJAVA_OVERVIEW.md
- **Pipeline steps:** SIMPLEMODE_PIPELINE.md
- **Code examples:** JAVA_CODE_SNIPPETS.md
- **This guide:** README_DOCUMENTATION.md

---

## 🎓 Start Implementing With

1. **CIGAR parser** - Good starting point, well-defined
2. **Variant struct** - Essential, straightforward mapping
3. **Statistics calculator** - Core logic, most complex
4. **Quality filter** - Validation gate
5. **Soft-clip realigner** - Advanced algorithm
6. **Full pipeline** - Integration

---

## 💾 Java Source Reference

**Most important files to study:**
```
VarDictJava/src/main/java/com/astrazeneca/vardict/
├── modes/SimpleMode.java                    # Entry point
├── variations/Variant.java                  # Data struct
├── modules/ToVarsBuilder.java               # Statistics ⭐⭐⭐
├── postprocessmodules/SimplePostProcessModule.java  # Filter
└── modules/VariationRealigner.java          # Realignment
```

---

**Print this and keep it handy while coding!**

Created: January 10, 2026 | Status: Ready for implementation ✓

