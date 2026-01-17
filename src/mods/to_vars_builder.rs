//! ToVarsBuilder - Convert low-level variations to high-level Variant objects with statistics
//!
//! This module implements the statistics calculation and variant creation logic for the VarDict pipeline.
//! It converts Variation objects (from CigarParser/VariationRealigner) into complete Variant objects
//! with allele frequencies, quality metrics, strand bias, and other statistical properties.
//!
//! **Mode:** Simple Mode only (no somatic, amplicon, or paired-sample logic)
//!
//! **Key Operations:**
//! 1. Group variations by genomic position
//! 2. Calculate 12+ statistics per variant
//! 3. Create Variant objects with all fields populated
//! 4. Group variants into Vars collections per position

use std::collections::HashMap;

use crate::mods::simple_variant_caller::SimpleVarKey;

/// Variant type enumeration
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VarType {
    /// Single nucleotide polymorphism
    SNV(char),
    /// Insertion of sequence
    Insertion(String),
    /// Deletion of N bases
    Deletion(usize),
    /// Complex indel (insertion + deletion)
    Complex { insertion: String, deletion: usize },
}

/// Strand bias value (matches Java VarDict format):
/// - 0: Can't assess (one strand only at low counts ≤12)
/// - 1: Has strand bias (high counts but fails balance test)  
/// - 2: No strand bias (balanced, good)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StrandBiasValue {
    CantAssess = 0,   // Low count, only one strand has reads
    HasBias = 1,      // High count but imbalanced  
    NoBias = 2,       // Balanced, no strand bias
}

impl StrandBiasValue {
    pub fn as_int(&self) -> i32 {
        match self {
            StrandBiasValue::CantAssess => 0,
            StrandBiasValue::HasBias => 1,
            StrandBiasValue::NoBias => 2,
        }
    }
}

/// Combined strand bias flag for ref and var (format: "refBias;varBias")
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StrandBiasFlag {
    pub ref_bias: StrandBiasValue,
    pub var_bias: StrandBiasValue,
}

impl Default for StrandBiasFlag {
    fn default() -> Self {
        StrandBiasFlag {
            ref_bias: StrandBiasValue::CantAssess,
            var_bias: StrandBiasValue::CantAssess,
        }
    }
}

impl StrandBiasFlag {
    /// Create from two bias values
    pub fn new(ref_bias: StrandBiasValue, var_bias: StrandBiasValue) -> Self {
        StrandBiasFlag { ref_bias, var_bias }
    }
    
    /// Format as Java-compatible string "refBias;varBias"
    pub fn to_string(&self) -> String {
        format!("{};{}", self.ref_bias.as_int(), self.var_bias.as_int())
    }
    
    /// Check if this is the "2;1" pattern (ref good, var has bias)
    pub fn is_ref_good_var_biased(&self) -> bool {
        self.ref_bias == StrandBiasValue::NoBias && self.var_bias == StrandBiasValue::HasBias
    }
}

/// A single variant with complete statistical information
#[derive(Debug, Clone)]
pub struct Variant {
    // === Identity ===
    pub description_string: String,
    pub refallele: String,
    pub varallele: String,
    pub vartype: VarType,

    // === Position ===
    pub start_position: i64,
    pub end_position: i64,

    // === Counts ===
    pub vars_count_on_forward: usize,
    pub vars_count_on_reverse: usize,
    pub position_coverage: usize,

    // === Frequencies ===
    pub frequency: f64,
    pub high_quality_reads_frequency: f64,

    // === Quality Metrics ===
    pub mean_position: f64,
    pub mean_quality: f64,
    pub mean_mapping_quality: f64,

    // === Flags ===
    pub strand_bias_flag: StrandBiasFlag,
    pub is_at_least_at_2_positions: bool,
    pub has_at_least_2_diff_qualities: bool,

    // === Context ===
    pub leftseq: String,
    pub rightseq: String,

    // === Special Features ===
    pub msi: f64,      // Microsatellite instability score
    pub msint: f64,    // Microsatellite interval
    pub shift3: i32,   // 3' shift for deletions
    pub nm: f64,       // Edit distance

    // === Quality counts ===
    pub high_qual_read_cnt: usize,  // Number of high-quality reads
    pub low_qual_read_cnt: usize,   // Number of low-quality reads
    pub hicov: usize,               // Position coverage by high-quality reads
    
    // === Reference counts (for non-reference variants at same position) ===
    pub ref_forward_count: usize,   // Forward reference reads at this position
    pub ref_reverse_count: usize,   // Reverse reference reads at this position
    
    // === Genotype ===
    pub genotype: String,
}

impl Variant {
    /// Create a new variant with default values
    pub fn new() -> Self {
        Variant {
            description_string: String::new(),
            refallele: String::new(),
            varallele: String::new(),
            vartype: VarType::SNV('N'),
            start_position: 0,
            end_position: 0,
            vars_count_on_forward: 0,
            vars_count_on_reverse: 0,
            position_coverage: 0,
            frequency: 0.0,
            high_quality_reads_frequency: 0.0,
            mean_position: 0.0,
            mean_quality: 0.0,
            mean_mapping_quality: 0.0,
            strand_bias_flag: StrandBiasFlag::default(),
            is_at_least_at_2_positions: false,
            has_at_least_2_diff_qualities: false,
            leftseq: String::new(),
            rightseq: String::new(),
            msi: 0.0,
            msint: 0.0,
            shift3: 0,
            nm: 0.0,
            high_qual_read_cnt: 0,
            low_qual_read_cnt: 0,
            hicov: 0,
            ref_forward_count: 0,
            ref_reverse_count: 0,
            genotype: "0/0".to_string(),
        }
    }

    /// Check if this variant passes quality filters (for testing)
    pub fn is_good_var(&self) -> bool {
        // Simple mode: just check for the "2;1" bias pattern
        // Don't filter on strand bias alone, use frequency too
        !self.strand_bias_flag.is_ref_good_var_biased() || self.frequency > 0.2
    }
}

impl Default for Variant {
    fn default() -> Self {
        Self::new()
    }
}

/// Collection of variants at a single genomic position
#[derive(Debug, Clone, Default)]
pub struct Vars {
    pub variants: Vec<Variant>,
    pub reference_variant: Option<Variant>,
    pub sv_flags: StructuralVariantFlags,
}

/// Structural variant flags
#[derive(Debug, Clone, Default)]
pub struct StructuralVariantFlags {
    pub splits: usize,
    pub pairs: usize,
    pub clusters: usize,
}

/// Input variation data (from CigarParser/VariationRealigner)
#[derive(Debug, Clone)]
pub struct VariationData {
    /// Position in the read (0-based)
    pub position_in_read: u32,
    /// Base quality at this position
    pub quality: u8,
    /// Mapping quality of the read
    pub mapping_quality: u8,
    /// Is this on the reverse strand?
    pub is_reverse: bool,
    /// Read ID (for deduplication)
    pub read_id: String,
}

/// ToVarsBuilder - Main statistics calculation engine
pub struct ToVarsBuilder {
    // Configuration
    min_frequency: f64,
    quality_threshold: u8,
    mapq_threshold: u8,
}

impl ToVarsBuilder {
    /// Create a new ToVarsBuilder with default settings
    pub fn new() -> Self {
        ToVarsBuilder {
            min_frequency: 0.02,          // 2% minimum
            quality_threshold: 20,         // Q20 minimum
            mapq_threshold: 20,            // MAPQ 20 minimum
        }
    }

    /// Set minimum allele frequency threshold
    pub fn with_min_frequency(mut self, freq: f64) -> Self {
        self.min_frequency = freq;
        self
    }

    /// Set base quality threshold
    pub fn with_quality_threshold(mut self, qual: f64) -> Self {
        self.quality_threshold = qual as u8;
        self
    }

    /// Set mapping quality threshold
    pub fn with_mapq_threshold(mut self, mapq: u8) -> Self {
        self.mapq_threshold = mapq;
        self
    }

    /// Calculate statistics for a group of variations at the same position
    /// representing the same variant
    fn calculate_variant_statistics(
        &self,
        variations: &[VariationData],
        position: i64,
        total_coverage: usize,
        var_type: VarType,
    ) -> Variant {
        let mut variant = Variant::new();
        variant.start_position = position;
        variant.end_position = position;
        variant.vartype = var_type.clone();
        variant.description_string = create_description_string(&var_type);
        variant.position_coverage = total_coverage;

        if variations.is_empty() {
            return variant;
        }

        // === Calculate strand counts ===
        let mut forward_count = 0;
        let mut reverse_count = 0;
        let mut positions = Vec::new();
        let mut qualities = Vec::new();
        let mut mapqs = Vec::new();

        for var in variations {
            if var.is_reverse {
                reverse_count += 1;
            } else {
                forward_count += 1;
            }
            positions.push(var.position_in_read as f64);
            qualities.push(var.quality as f64);
            mapqs.push(var.mapping_quality as f64);
        }

        variant.vars_count_on_forward = forward_count;
        variant.vars_count_on_reverse = reverse_count;

        // === Calculate frequency ===
        let total_var_count = forward_count + reverse_count;
        variant.frequency = if total_coverage > 0 {
            total_var_count as f64 / total_coverage as f64
        } else {
            0.0
        };

        // === Position in read analysis ===
        if !positions.is_empty() {
            let (mean_pos, _std_pos) = calculate_mean_and_std(&positions);
            variant.mean_position = mean_pos;
            variant.is_at_least_at_2_positions = has_at_least_2_distinct_f64(&positions);
        }

        // === Base quality analysis ===
        if !qualities.is_empty() {
            let (mean_qual, _std_qual) = calculate_mean_and_std(&qualities);
            variant.mean_quality = mean_qual;
            
            // Check if qualities are distinct
            let qual_ints: Vec<u8> = qualities.iter().map(|q| *q as u8).collect();
            variant.has_at_least_2_diff_qualities = has_at_least_2_distinct(&qual_ints);
        }

        // === Mapping quality analysis ===
        if !mapqs.is_empty() {
            let (mean_mapq, _std_mapq) = calculate_mean_and_std(&mapqs);
            variant.mean_mapping_quality = mean_mapq;
        }

        // === Strand bias calculation ===
        // Calculate variant strand bias from variant read counts
        let var_bias = check_strand_bias(forward_count, reverse_count);
        // Reference bias will be set later when we have reference counts
        // For now, set just the variant bias (ref defaults to CantAssess)
        variant.strand_bias_flag = StrandBiasFlag {
            ref_bias: StrandBiasValue::CantAssess,
            var_bias,
        };

        // === High-quality read filtering (Simple Mode) ===
        let high_quality_count = variations
            .iter()
            .filter(|v| v.quality >= self.quality_threshold)
            .count();
        
        if high_quality_count > 0 && total_coverage > 0 {
            variant.high_quality_reads_frequency = high_quality_count as f64 / total_coverage as f64;
        }

        // === Reference allele and variant allele ===
        variant.refallele = "N".to_string(); // Placeholder - would be determined from context
        variant.varallele = "N".to_string(); // Placeholder

        // === Genotype prediction (Simple Mode) ===
        // No anchor base available in this context
        variant.genotype = determine_genotype(&variant.refallele, &variant.varallele, variant.frequency, None);

        variant
    }

    /// Group variations by position and variant, then calculate statistics
    /// 
    /// Accepts SimpleVarKey enum for performance-optimized variant representation
    pub fn build_variants(
        &self,
        all_variations: Vec<(i64, SimpleVarKey, VariationData)>,
        total_coverage_by_pos: &HashMap<i64, usize>,
    ) -> HashMap<i64, Vars> {
        // Group by position
        let mut by_position: HashMap<i64, Vec<(SimpleVarKey, VariationData)>> = HashMap::new();
        for (pos, var_key, var_data) in all_variations {
            by_position
                .entry(pos)
                .or_insert_with(Vec::new)
                .push((var_key, var_data));
        }

        // Process each position
        let mut result = HashMap::new();
        for (position, variants_at_pos) in by_position {
            let coverage = *total_coverage_by_pos.get(&position).unwrap_or(&0);

            // Group by variant description at this position
            let mut by_variant: HashMap<SimpleVarKey, Vec<VariationData>> = HashMap::new();
            for (var_key, var_data) in variants_at_pos {
                by_variant
                    .entry(var_key)
                    .or_insert_with(Vec::new)
                    .push(var_data);
            }

            // Calculate statistics for each variant
            let mut variant_list = Vec::new();
            for (var_key, var_data_list) in by_variant {
                // Determine variant type from SimpleVarKey
                let var_type = var_key_to_var_type(&var_key);

                let mut variant = self.calculate_variant_statistics(
                    &var_data_list,
                    position,
                    coverage,
                    var_type,
                );

                // Set ref and alt alleles from SimpleVarKey
                variant.refallele = var_key.ref_allele();
                variant.varallele = var_key.alt_allele();
                variant.description_string = var_key.to_key_string();

                // Calculate end position based on reference allele length
                if variant.refallele.len() > 1 {
                    variant.end_position = position + variant.refallele.len() as i64 - 1;
                }

                variant_list.push(variant);
            }

            let vars = Vars {
                variants: variant_list,
                reference_variant: None,
                sv_flags: StructuralVariantFlags::default(),
            };

            result.insert(position, vars);
        }

        result
    }
}

impl Default for ToVarsBuilder {
    fn default() -> Self {
        Self::new()
    }
}

// ============================================================================
// Helper Functions for Statistics Calculation
// ============================================================================

/// Calculate mean and standard deviation
pub fn calculate_mean_and_std(values: &[f64]) -> (f64, f64) {
    if values.is_empty() {
        return (0.0, 0.0);
    }

    let mean = values.iter().sum::<f64>() / values.len() as f64;
    let variance = values
        .iter()
        .map(|v| (v - mean).powi(2))
        .sum::<f64>()
        / values.len() as f64;

    (mean, variance.sqrt())
}

/// Check if a list contains at least 2 distinct values
pub fn has_at_least_2_distinct<T: std::cmp::Eq + std::hash::Hash>(values: &[T]) -> bool {
    let mut seen = std::collections::HashSet::new();
    for v in values {
        seen.insert(v);
        if seen.len() >= 2 {
            return true;
        }
    }
    false
}

/// Check if a list of floats contains at least 2 distinct values (with tolerance for floating point)
pub fn has_at_least_2_distinct_f64(values: &[f64]) -> bool {
    if values.len() < 2 {
        return false;
    }
    
    // Check if all values are approximately the same (within epsilon)
    let first = values[0];
    let epsilon = 1e-9;
    
    for v in &values[1..] {
        if (v - first).abs() > epsilon {
            return true; // Found a different value
        }
    }
    
    false // All values are essentially the same
}

/// Calculate strand bias value based on forward/reverse counts
/// Matches Java VarDict's strandBias function exactly:
/// - For total ≤12: returns 2 if both strands have reads, 0 if only one strand
/// - For total >12: returns 2 if balanced (both ≥5% AND both ≥minBiasReads), else 1
pub fn check_strand_bias(forward: usize, reverse: usize) -> StrandBiasValue {
    let total = forward + reverse;
    let bias_threshold = 0.05;    // Java: instance().conf.bias = 0.05
    let min_bias_reads = 2usize;  // Java: instance().conf.minBiasReads = 2
    
    // using p=0.01, because prop.test(1,12) = 0.01
    if total <= 12 {
        // For low counts, just check if both strands have any reads
        if forward > 0 && reverse > 0 {
            StrandBiasValue::NoBias  // 2: both strands represented
        } else {
            StrandBiasValue::CantAssess  // 0: only one strand
        }
    } else {
        // For higher counts, check if balanced
        let fwd_ratio = forward as f64 / total as f64;
        let rev_ratio = reverse as f64 / total as f64;
        
        if fwd_ratio >= bias_threshold 
            && rev_ratio >= bias_threshold
            && forward >= min_bias_reads 
            && reverse >= min_bias_reads 
        {
            StrandBiasValue::NoBias  // 2: balanced, no bias
        } else {
            StrandBiasValue::HasBias  // 1: strand bias present
        }
    }
}

/// Determine genotype from frequency (Simple Mode)
/// Determine genotype string from reference and variant alleles
/// Format: "REF/ALT" using actual allele sequences
/// For reference calls (ref==alt): "REF/REF"
/// For heterozygous: "REF/ALT"  
/// For homozygous alternate: "ALT/ALT"
/// 
/// Java-compatible format:
/// - Reference calls: "REF/REF"
/// - Pure insertions: "REF/+N" where N is the inserted length
/// - Deletions: "-N<seq>/-N<seq>" format
/// - Same-length complex: "ALT+REF[1:]/ALT" (extended genotype1)
/// - Complex with length diff: Java-compatible extended format with anchor
pub fn determine_genotype(refallele: &str, varallele: &str, frequency: f64, anchor_base: Option<char>) -> String {
    if refallele == varallele {
        // Reference call: G→G becomes "G/G"
        format!("{}/{}", refallele, refallele)
    } else {
        // Determine if this is an insertion, deletion, or substitution
        let ref_len = refallele.len();
        let var_len = varallele.len();
        
        if var_len > ref_len && varallele.starts_with(refallele) {
            // Pure insertion: alt starts with ref and is longer
            // Format: "REF/+N" where N is the insertion length
            let ins_len = var_len - ref_len;
            format!("{}/+{}", refallele, ins_len)
        } else if ref_len > var_len && refallele.starts_with(varallele) {
            // Pure deletion: ref starts with alt and is longer
            // Format: "-N<full_seq>/-N<full_seq>" for homozygous
            let del_len = ref_len - var_len;
            let del_notation = format!("-{}{}", del_len, varallele);
            format!("{}/{}", del_notation, del_notation)
        } else if ref_len > var_len {
            // Complex deletion: ref longer than alt but not a prefix match
            // Format: "-N<alt_seq><ref_from_del_len>/-N<alt_seq>"
            let del_len = ref_len - var_len;
            // genotype1 = -N + alt + ref[del_len:] (skip the deleted bases at start)
            // genotype2 = -N + alt
            let genotype1 = format!("-{}{}{}", del_len, varallele, &refallele[del_len..]);
            let genotype2 = format!("-{}{}", del_len, varallele);
            format!("{}/{}", genotype1, genotype2)
        } else if ref_len == var_len && ref_len > 1 {
            // Same-length complex variant (MNV or substitution)
            // Java extends genotype1 with ref[1:] (reference bases after first char)
            // genotype1 = varallele + refallele[1:]
            // genotype2 = varallele
            let genotype1 = format!("{}{}", varallele, &refallele[1..]);
            let genotype2 = varallele.to_string();
            format!("{}/{}", genotype1, genotype2)
        } else if var_len > ref_len {
            // Complex insertion: alt longer than ref, not a pure insertion
            // Java uses +N notation for genotype2
            // genotype1 = anchor_base + refallele (extended context)
            // genotype2 = +<varallele_len + 1> (Java counts description length - 1)
            let genotype1 = if let Some(anchor) = anchor_base {
                format!("{}{}", anchor, refallele)
            } else {
                refallele.to_string()
            };
            // Java's +N comes from description string length - 1
            // For complex insertion with matched sequence, the description is longer
            let genotype2 = format!("+{}", var_len + 1);
            format!("{}/{}", genotype1, genotype2)
        } else {
            // Single base substitution  
            format!("{}/{}", refallele, varallele)
        }
    }
}

/// Create description string based on variant type
pub fn create_description_string(var_type: &VarType) -> String {
    match var_type {
        VarType::SNV(base) => base.to_string(),
        VarType::Insertion(seq) => format!("+{}", seq),
        VarType::Deletion(len) => format!("-{}", len),
        VarType::Complex {
            insertion,
            deletion,
        } => format!("{}#{}", insertion, if *deletion > 0 {
            format!("-{}", deletion)
        } else {
            String::new()
        }),
    }
}

/// Calculate shift3 (3' shift allowance) for deletions
pub fn calculate_shift3(
    reference: &[u8],
    position: usize,
    deletion_length: usize,
) -> i32 {
    if position + deletion_length >= reference.len() {
        return 0;
    }

    let mut shift = 0;
    let mut pos = position;

    while pos + deletion_length < reference.len() {
        if reference[pos] == reference[pos + deletion_length] {
            shift += 1;
            pos += 1;
        } else {
            break;
        }
    }

    shift
}

/// Determine the variant type from the variant key
/// Variant key format examples:
/// - "A>T" -> SNV
/// - "+ATG" -> Insertion
/// - "-ATG" -> Deletion
/// - "AT>C" or "AT>GC" -> Complex/MNV
fn infer_variant_type(variant_key: &str) -> VarType {
    if variant_key.starts_with('+') {
        // Insertion: extract sequence after the '+'
        let seq = variant_key[1..].to_string();
        VarType::Insertion(seq)
    } else if variant_key.starts_with('-') {
        // Deletion: extract count after the '-'
        if let Ok(count) = variant_key[1..].parse::<usize>() {
            VarType::Deletion(count)
        } else {
            // Fallback if we can't parse the count
            let seq = variant_key[1..].to_string();
            VarType::Deletion(seq.len())
        }
    } else if let Some(pos) = variant_key.find('>') {
        let ref_part = &variant_key[..pos];
        let alt_part = &variant_key[pos + 1..];

        // Check if it's a simple SNV (single base ref and alt)
        if ref_part.len() == 1 && alt_part.len() == 1 {
            if let Some(base) = alt_part.chars().next() {
                VarType::SNV(base)
            } else {
                VarType::SNV('N')
            }
        } else {
            // Multi-base substitution or complex variant
            // For now, treat as Complex with empty insertion and deletion
            VarType::Complex {
                insertion: alt_part.to_string(),
                deletion: ref_part.len(),
            }
        }
    } else {
        // Fallback for unexpected formats
        VarType::SNV('N')
    }
}

/// Convert SimpleVarKey enum to VarType for statistics calculation
fn var_key_to_var_type(var_key: &SimpleVarKey) -> VarType {
    match var_key {
        SimpleVarKey::SNV { alt_base, .. } => VarType::SNV(char::from(*alt_base)),
        SimpleVarKey::Ins { seq } => VarType::Insertion(String::from_utf8_lossy(seq).to_string()),
        SimpleVarKey::Del { len, .. } => VarType::Deletion(*len as usize),
        SimpleVarKey::Complex { ref_seq, alt_seq } => VarType::Complex {
            insertion: String::from_utf8_lossy(alt_seq).to_string(),
            deletion: ref_seq.len(),
        },
    }
}

/// Validate and normalize reference allele by replacing IUPAC ambiguity codes
/// with their first possible base.
/// 
/// IUPAC ambiguity codes:
/// - M (A or C) -> A
/// - R (A or G) -> A
/// - W (A or T) -> A
/// - S (C or G) -> C
/// - Y (C or T) -> C
/// - K (G or T) -> G
/// - V (A or C or G) -> A
/// - H (A or C or T) -> A
/// - D (A or G or T) -> A
/// - B (C or G or T) -> C
/// 
/// From Java ToVarsBuilderTest.java
pub fn validate_ref_allele(allele: &str) -> String {
    allele
        .chars()
        .map(|c| match c {
            'A' | 'C' | 'G' | 'T' | 'N' => c,
            'M' | 'R' | 'W' | 'V' | 'H' | 'D' => 'A', // First base is A
            'S' | 'Y' | 'B' => 'C',                   // First base is C
            'K' => 'G',                               // First base is G
            _ => c, // Keep unknown chars as-is
        })
        .collect()
}

// ============================================================================
// Tests
// ============================================================================

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
    fn test_mean_and_std_empty() {
        let values: Vec<f64> = vec![];
        let (mean, std) = calculate_mean_and_std(&values);
        assert_eq!(mean, 0.0);
        assert_eq!(std, 0.0);
    }

    #[test]
    fn test_mean_and_std_single() {
        let values = vec![5.0];
        let (mean, std) = calculate_mean_and_std(&values);
        assert_eq!(mean, 5.0);
        assert_eq!(std, 0.0);
    }

    #[test]
    fn test_has_at_least_2_distinct_yes() {
        let values = vec![10, 20, 30];
        assert!(has_at_least_2_distinct(&values));
    }

    #[test]
    fn test_has_at_least_2_distinct_no() {
        let values = vec![10, 10, 10];
        assert!(!has_at_least_2_distinct(&values));
    }

    #[test]
    fn test_has_at_least_2_distinct_exactly_2() {
        let values = vec![10, 20];
        assert!(has_at_least_2_distinct(&values));
    }

    #[test]
    fn test_strand_bias_no_bias() {
        // 10 forward, 10 reverse - total 20 > 12, both are >= 5% and >= 2 = NoBias (2)
        assert_eq!(
            check_strand_bias(10, 10),
            StrandBiasValue::NoBias
        );
    }

    #[test]
    fn test_strand_bias_low_count_both_strands() {
        // Total <= 12, but both strands have reads = NoBias (2)
        assert_eq!(
            check_strand_bias(5, 3),
            StrandBiasValue::NoBias
        );
    }

    #[test]
    fn test_strand_bias_low_count_one_strand() {
        // Total <= 12, only one strand has reads = CantAssess (0)
        assert_eq!(
            check_strand_bias(10, 0),
            StrandBiasValue::CantAssess
        );
    }

    #[test]
    fn test_strand_bias_high_count_imbalanced() {
        // Total > 12, but 5/100 = 5% is exactly at threshold, and both have >= 2 reads
        // Forward: 100/105 = 95.2% >= 5%, Reverse: 5/105 = 4.76% < 5% = HasBias
        assert_eq!(
            check_strand_bias(100, 5),
            StrandBiasValue::HasBias
        );
    }

    #[test]
    fn test_strand_bias_balanced_high_count() {
        // Total > 12, both strands >= 5% and >= 2 reads = NoBias
        assert_eq!(
            check_strand_bias(15, 5),  // 15/20=75% and 5/20=25%, both >= 5%
            StrandBiasValue::NoBias
        );
    }

    #[test]
    fn test_determine_genotype_ref_call() {
        // Reference call: G→G should be "G/G"
        assert_eq!(determine_genotype("G", "G", 1.0, None), "G/G");
    }

    #[test]
    fn test_determine_genotype_het() {
        // SNV: single base substitution uses "ref/alt" format (matching Java)
        // Frequency is not used for SNVs
        assert_eq!(determine_genotype("C", "T", 0.3, None), "C/T");
    }

    #[test]
    fn test_determine_genotype_hom_alt() {
        // SNV: still uses "ref/alt" format
        assert_eq!(determine_genotype("C", "T", 0.8, None), "C/T");
    }

    #[test]
    fn test_create_description_snv() {
        assert_eq!(
            create_description_string(&VarType::SNV('A')),
            "A"
        );
    }

    #[test]
    fn test_create_description_insertion() {
        assert_eq!(
            create_description_string(&VarType::Insertion("ACGT".to_string())),
            "+ACGT"
        );
    }

    #[test]
    fn test_create_description_deletion() {
        assert_eq!(
            create_description_string(&VarType::Deletion(5)),
            "-5"
        );
    }

    #[test]
    fn test_create_description_complex() {
        assert_eq!(
            create_description_string(&VarType::Complex {
                insertion: "AC".to_string(),
                deletion: 3
            }),
            "AC#-3"
        );
    }

    #[test]
    fn test_calculate_shift3() {
        let ref_seq = b"AAAA";
        let shift = calculate_shift3(ref_seq, 0, 2);
        assert_eq!(shift, 2); // Can shift right by 2 bases
    }

    #[test]
    fn test_calculate_shift3_no_shift() {
        let ref_seq = b"ATCG";
        let shift = calculate_shift3(ref_seq, 0, 1);
        assert_eq!(shift, 0); // Cannot shift
    }

    #[test]
    fn test_variant_new() {
        let v = Variant::new();
        assert_eq!(v.frequency, 0.0);
        assert_eq!(v.strand_bias_flag, StrandBiasFlag::default());
        assert_eq!(v.genotype, "0/0");
    }

    #[test]
    fn test_variant_default() {
        let v = Variant::default();
        assert_eq!(v.frequency, 0.0);
    }

    #[test]
    fn test_to_vars_builder_new() {
        let builder = ToVarsBuilder::new();
        assert_eq!(builder.min_frequency, 0.02);
        assert_eq!(builder.quality_threshold, 20);
        assert_eq!(builder.mapq_threshold, 20);
    }

    #[test]
    fn test_to_vars_builder_with_min_frequency() {
        let builder = ToVarsBuilder::new()
            .with_min_frequency(0.05);
        assert_eq!(builder.min_frequency, 0.05);
    }

    // === Phase 2 Tests ===

    #[test]
    fn test_infer_variant_type_snv() {
        let var_type = infer_variant_type("A>T");
        assert_eq!(var_type, VarType::SNV('T'));
    }

    #[test]
    fn test_infer_variant_type_insertion() {
        let var_type = infer_variant_type("+ATG");
        match var_type {
            VarType::Insertion(seq) => assert_eq!(seq, "ATG"),
            _ => panic!("Expected insertion"),
        }
    }

    #[test]
    fn test_infer_variant_type_deletion() {
        let var_type = infer_variant_type("-3");
        match var_type {
            VarType::Deletion(len) => assert_eq!(len, 3),
            _ => panic!("Expected deletion"),
        }
    }

    #[test]
    fn test_infer_variant_type_complex() {
        let var_type = infer_variant_type("AT>GC");
        match var_type {
            VarType::Complex {
                insertion,
                deletion,
            } => {
                assert_eq!(insertion, "GC");
                assert_eq!(deletion, 2);
            }
            _ => panic!("Expected complex variant"),
        }
    }

    #[test]
    fn test_has_at_least_2_distinct_f64_yes() {
        let values = vec![1.0, 2.0, 3.0];
        assert!(has_at_least_2_distinct_f64(&values));
    }

    #[test]
    fn test_has_at_least_2_distinct_f64_no() {
        let values = vec![1.0, 1.0, 1.0];
        assert!(!has_at_least_2_distinct_f64(&values));
    }

    #[test]
    fn test_has_at_least_2_distinct_f64_empty() {
        let values: Vec<f64> = vec![];
        assert!(!has_at_least_2_distinct_f64(&values));
    }

    #[test]
    fn test_has_at_least_2_distinct_f64_single() {
        let values = vec![1.0];
        assert!(!has_at_least_2_distinct_f64(&values));
    }

    #[test]
    fn test_calculate_variant_statistics_simple() {
        let builder = ToVarsBuilder::new();
        
        // Create test data: 2 forward strand, 2 reverse strand reads
        let variations = vec![
            VariationData {
                position_in_read: 10,
                quality: 30,
                mapping_quality: 60,
                is_reverse: false,
                read_id: "read1".to_string(),
            },
            VariationData {
                position_in_read: 15,
                quality: 25,
                mapping_quality: 60,
                is_reverse: false,
                read_id: "read2".to_string(),
            },
            VariationData {
                position_in_read: 20,
                quality: 28,
                mapping_quality: 50,
                is_reverse: true,
                read_id: "read3".to_string(),
            },
            VariationData {
                position_in_read: 25,
                quality: 22,
                mapping_quality: 50,
                is_reverse: true,
                read_id: "read4".to_string(),
            },
        ];

        let variant = builder.calculate_variant_statistics(
            &variations,
            1000,    // position
            100,     // coverage
            VarType::SNV('A'),
        );

        // Verify counts
        assert_eq!(variant.vars_count_on_forward, 2);
        assert_eq!(variant.vars_count_on_reverse, 2);
        assert_eq!(variant.position_coverage, 100);

        // Verify frequency: 4 variants out of 100 coverage = 4%
        assert!((variant.frequency - 0.04).abs() < 0.001);

        // Genotype is set to default "N/N" in calculate_variant_statistics
        // The actual genotype is determined later when full context is available
        assert_eq!(variant.genotype, "N/N");

        // Verify strand bias: 2 forward, 2 reverse, total=4 which is <=12
        // Both strands have reads, so should be NoBias (2)
        assert_eq!(variant.strand_bias_flag.var_bias, StrandBiasValue::NoBias);

        // Verify position mean: (10 + 15 + 20 + 25) / 4 = 17.5
        assert!((variant.mean_position - 17.5).abs() < 0.1);

        // Verify quality mean: (30 + 25 + 28 + 22) / 4 = 26.25
        assert!((variant.mean_quality - 26.25).abs() < 0.1);

        // Verify mapping quality mean: (60 + 60 + 50 + 50) / 4 = 55.0
        assert!((variant.mean_mapping_quality - 55.0).abs() < 0.1);
    }
    #[test]
    fn test_calculate_variant_statistics_biased() {
        let builder = ToVarsBuilder::new();
        
        // Create test data: 11 forward, 1 reverse (strongly biased, ratio = 11.0 > 10.0)
        let mut variations = Vec::new();
        for i in 0..11 {
            variations.push(VariationData {
                position_in_read: 10 + i,
                quality: 30,
                mapping_quality: 60,
                is_reverse: false,
                read_id: format!("read_f{}", i),
            });
        }
        variations.push(VariationData {
            position_in_read: 20,
            quality: 20,
            mapping_quality: 40,
            is_reverse: true,
            read_id: "read_r1".to_string(),
        });

        let variant = builder.calculate_variant_statistics(
            &variations,
            1000,
            100,
            VarType::SNV('T'),
        );

        // Verify counts
        assert_eq!(variant.vars_count_on_forward, 11);
        assert_eq!(variant.vars_count_on_reverse, 1);

        // Verify frequency: 12 variants out of 100 = 12%
        assert!((variant.frequency - 0.12).abs() < 0.001);

        // Verify strand bias: 11 forward, 1 reverse, total=12 which is <=12
        // Both strands have reads, so should be NoBias (2) per Java logic
        assert_eq!(variant.strand_bias_flag.var_bias, StrandBiasValue::NoBias);

        // Genotype is set to default "N/N" in calculate_variant_statistics
        // The actual genotype is determined later when full context is available
        assert_eq!(variant.genotype, "N/N");
    }

    #[test]
    fn test_calculate_variant_statistics_empty() {
        let builder = ToVarsBuilder::new();
        
        // Empty variations
        let variations: Vec<VariationData> = vec![];

        let variant = builder.calculate_variant_statistics(
            &variations,
            1000,
            100,
            VarType::SNV('A'),
        );

        // Verify zero counts
        assert_eq!(variant.vars_count_on_forward, 0);
        assert_eq!(variant.vars_count_on_reverse, 0);
        assert_eq!(variant.frequency, 0.0);
        assert_eq!(variant.genotype, "0/0");
    }

    #[test]
    fn test_build_variants_grouping() {
        let builder = ToVarsBuilder::new();
        
        // Create variations as vector of tuples (position, SimpleVarKey, variation_data)
        let all_variations = vec![
            // Position 1000: SNV A>T (2 reads)
            (
                1000i64,
                SimpleVarKey::snv(b'A', b'T'),
                VariationData {
                    position_in_read: 10,
                    quality: 30,
                    mapping_quality: 60,
                    is_reverse: false,
                    read_id: "read1".to_string(),
                },
            ),
            (
                1000i64,
                SimpleVarKey::snv(b'A', b'T'),
                VariationData {
                    position_in_read: 15,
                    quality: 25,
                    mapping_quality: 60,
                    is_reverse: true,
                    read_id: "read2".to_string(),
                },
            ),
            // Position 1000: Insertion +AG (1 read)
            (
                1000i64,
                SimpleVarKey::insertion(b"AG"),
                VariationData {
                    position_in_read: 20,
                    quality: 28,
                    mapping_quality: 50,
                    is_reverse: false,
                    read_id: "read3".to_string(),
                },
            ),
            // Position 1010: SNV C>G (1 read)
            (
                1010i64,
                SimpleVarKey::snv(b'C', b'G'),
                VariationData {
                    position_in_read: 22,
                    quality: 32,
                    mapping_quality: 55,
                    is_reverse: true,
                    read_id: "read4".to_string(),
                },
            ),
        ];

        // Create coverage map
        let mut coverage_map = HashMap::new();
        coverage_map.insert(1000i64, 100);
        coverage_map.insert(1010i64, 100);

        let result = builder.build_variants(all_variations, &coverage_map);

        // Verify we have 2 positions
        assert_eq!(result.len(), 2);

        // Verify position 1000 has 2 variants
        let vars_at_1000 = result.get(&1000).unwrap();
        assert_eq!(vars_at_1000.variants.len(), 2);

        // Verify position 1010 has 1 variant
        let vars_at_1010 = result.get(&1010).unwrap();
        assert_eq!(vars_at_1010.variants.len(), 1);

        // Verify variant types are correct at position 1000
        let var_types: Vec<_> = vars_at_1000.variants.iter()
            .map(|v| &v.vartype)
            .collect();
        assert!(var_types.contains(&&VarType::SNV('T')));
        assert!(var_types.iter().any(|vt| matches!(vt, VarType::Insertion(_))));
    }

    // ========================================================================
    // Tests ported from Java ToVarsBuilderTest.java
    // ========================================================================

    /// Port of Java testValidateRefAllele()
    /// Tests IUPAC ambiguity code resolution
    #[test]
    fn test_validate_ref_allele_single_bases() {
        // From Java: List<String> alleles = Arrays.asList("A", "C", "G", "T", "N", "M", "R", "W", "S", "Y", "K", "V", "H", "D", "B");
        // Expected:  List<String> expected = Arrays.asList("A", "C", "G", "T", "N", "A", "A", "A", "C", "C", "G", "A", "A", "A", "C");
        let test_cases = vec![
            ("A", "A"),
            ("C", "C"),
            ("G", "G"),
            ("T", "T"),
            ("N", "N"),
            ("M", "A"),  // M (A or C) -> A
            ("R", "A"),  // R (A or G) -> A
            ("W", "A"),  // W (A or T) -> A
            ("S", "C"),  // S (C or G) -> C
            ("Y", "C"),  // Y (C or T) -> C
            ("K", "G"),  // K (G or T) -> G
            ("V", "A"),  // V (A or C or G) -> A
            ("H", "A"),  // H (A or C or T) -> A
            ("D", "A"),  // D (A or G or T) -> A
            ("B", "C"),  // B (C or G or T) -> C
        ];

        for (input, expected) in test_cases {
            assert_eq!(
                validate_ref_allele(input),
                expected,
                "Failed for input: {}",
                input
            );
        }
    }

    /// Port of Java testValidateRefAllele() - complex alleles part
    #[test]
    fn test_validate_ref_allele_complex() {
        // From Java: List<String> alleles_complex = Arrays.asList("ANYCGT", "MRACT", "CCGKBG");
        // Expected:  List<String> expected_complex = Arrays.asList("ANCCGT", "AAACT", "CCGGCG");
        let test_cases = vec![
            ("ANYCGT", "ANCCGT"),  // Y->C
            ("MRACT", "AAACT"),    // M->A, R->A
            ("CCGKBG", "CCGGCG"),  // K->G, B->C
        ];

        for (input, expected) in test_cases {
            assert_eq!(
                validate_ref_allele(input),
                expected,
                "Failed for input: {}",
                input
            );
        }
    }

    /// Port of Java createInsertion/createVariant test
    /// Tests the statistics calculation with specific input values
    /// 
    /// Java test input:
    ///   varsCount = 4, varsCountOnForward = 3, varsCountOnReverse = 5
    ///   meanPosition = 9, meanQuality = 10.5, meanMappingQuality = 31
    ///   numberOfMismatches = 8, highQualityReadsCount = 44, lowQualityReadsCount = 35
    ///
    /// Java test expected output:
    ///   frequency = 0.4 (4/10 coverage)
    ///   meanPosition = 2.2 (9/4 = 2.25, but Java might round)
    ///   meanQuality = 2.6 (10.5/4 = 2.625)
    ///   meanMappingQuality = 7.8 (31/4 = 7.75)
    ///   strandBiasFlag = "2" (strong bias: ratio 5/3 = 1.67, but with 0 on one side = strong)
    #[test]
    fn test_statistics_calculation_java_values() {
        // Note: The Java test uses pre-calculated sums (e.g., meanPosition=9 is sum, divide by varsCount=4 to get 2.25)
        // In our implementation, we calculate means directly from individual reads
        // So we create 4 reads with properties that sum to the Java values:
        //   - positions sum to 9 (e.g., 1, 2, 3, 3)
        //   - qualities sum to 10.5 (e.g., 2.0, 2.5, 3.0, 3.0)
        //   - mapqs sum to 31 (e.g., 7, 8, 8, 8)
        
        let builder = ToVarsBuilder::new();
        
        // Create 4 variations: 3 forward, 1 reverse (we can't have 5 reverse with only 4 total)
        // Note: Java test has varsCount=4 but varsCountOnForward=3 + varsCountOnReverse=5 = 8
        // This is inconsistent in the Java test. We'll use 4 total = 3 forward + 1 reverse
        let variations = vec![
            VariationData {
                position_in_read: 1,
                quality: 20,  // Will contribute to mean
                mapping_quality: 7,
                is_reverse: false,
                read_id: "read1".to_string(),
            },
            VariationData {
                position_in_read: 2,
                quality: 25,
                mapping_quality: 8,
                is_reverse: false,
                read_id: "read2".to_string(),
            },
            VariationData {
                position_in_read: 3,
                quality: 30,
                mapping_quality: 8,
                is_reverse: false,
                read_id: "read3".to_string(),
            },
            VariationData {
                position_in_read: 3,
                quality: 30,
                mapping_quality: 8,
                is_reverse: true,
                read_id: "read4".to_string(),
            },
        ];

        let variant = builder.calculate_variant_statistics(
            &variations,
            1234567,  // Same position as Java test
            10,       // Coverage = 10 (so 4/10 = 0.4 frequency)
            VarType::SNV('T'),
        );

        // Check frequency: 4/10 = 0.4
        assert!((variant.frequency - 0.4).abs() < 0.001, 
                "frequency: expected 0.4, got {}", variant.frequency);
        
        // Check counts
        assert_eq!(variant.vars_count_on_forward, 3);
        assert_eq!(variant.vars_count_on_reverse, 1);
        
        // Check mean position: (1+2+3+3)/4 = 2.25
        assert!((variant.mean_position - 2.25).abs() < 0.01,
                "mean_position: expected 2.25, got {}", variant.mean_position);
        
        // Check mean quality: (20+25+30+30)/4 = 26.25
        assert!((variant.mean_quality - 26.25).abs() < 0.01,
                "mean_quality: expected 26.25, got {}", variant.mean_quality);
        
        // Check mean mapping quality: (7+8+8+8)/4 = 7.75
        assert!((variant.mean_mapping_quality - 7.75).abs() < 0.01,
                "mean_mapping_quality: expected 7.75, got {}", variant.mean_mapping_quality);
        
        // Check strand bias: 3 forward, 1 reverse, total=4 which is <=12
        // Both strands have reads, so should be NoBias (2) per Java logic
        assert_eq!(variant.strand_bias_flag.var_bias, StrandBiasValue::NoBias);
    }
}
