//! Output Variant Formatting - Convert Variant objects to tab-delimited output
//!
//! This module implements the output formatting for VarDict variants.
//! Simple Mode produces 36 columns (no Fisher's exact test columns).
//!
//! **Output Format:** Tab-delimited text (one variant per line)
//!
//! **Columns (36 for Simple Mode):**
//! 1. Sample name
//! 2. Gene name
//! 3. Chromosome
//! 4. Start position
//! 5. End position
//! 6. Reference allele
//! 7. Variant allele
//! 8. Total coverage at position
//! 9. Variant coverage
//! 10. Reference forward count
//! 11. Reference reverse count
//! 12. Variant forward count
//! 13. Variant reverse count
//! 14. Genotype (0/0, 0/1, 1/1)
//! 15. Variant frequency
//! 16. Strand bias (e.g., "2;1")
//! 17. Mean position in read
//! 18. Position standard deviation flag (0 or 1)
//! 19. Mean base quality
//! 20. Quality standard deviation flag (0 or 1)
//! 21. Mean mapping quality
//! 22. Quality ratio (high/low quality)
//! 23. High-quality frequency
//! 24. Extra frequency
//! 25. 3' shift
//! 26. MSI score
//! 27. MSI interval
//! 28. Number of mismatches
//! 29. High-quality count
//! 30. High-quality coverage
//! 31. Left sequence context
//! 32. Right sequence context
//! 33. Region (chr:start-end)
//! 34. Variant type
//! 35. Duplicate rate
//! 36. Structural variant info

use crate::mods::to_vars_builder::{Variant, VarType, StrandBiasFlag, var_type_string};

/// Region information for output
#[derive(Debug, Clone, Default)]
pub struct Region {
    pub chr: String,
    pub start: i64,
    pub end: i64,
    pub gene: String,
}

impl Region {
    pub fn new(chr: &str, start: i64, end: i64, gene: &str) -> Self {
        Region {
            chr: chr.to_string(),
            start,
            end,
            gene: gene.to_string(),
        }
    }

    /// Format as "chr:start-end"
    pub fn to_region_string(&self) -> String {
        format!("{}:{}-{}", normalize_chr_for_output(&self.chr), self.start, self.end)
    }
}

/// Normalize chromosome name for output to match VarDict Java (drop leading "chr")
fn normalize_chr_for_output(chr: &str) -> String {
    chr.strip_prefix("chr").unwrap_or(chr).to_string()
}

/// Simple Output Variant - 36 column format for Simple Mode
#[derive(Debug, Clone)]
pub struct SimpleOutputVariant {
    // Identity
    pub sample: String,
    pub gene: String,
    pub chr: String,
    pub start_position: i64,
    pub end_position: i64,
    pub ref_allele: String,
    pub var_allele: String,

    // Coverage
    pub total_coverage: usize,
    pub variant_coverage: usize,
    pub reference_forward_count: usize,
    pub reference_reverse_count: usize,
    pub variant_forward_count: usize,
    pub variant_reverse_count: usize,

    // Genotype & Frequency
    pub genotype: String,
    pub frequency: f64,
    pub bias: String,  // "flag;flag" format

    // Position metrics
    pub pmean: f64,    // Mean position in read
    pub pstd: i32,     // Position std flag (0 or 1)

    // Quality metrics
    pub qual: f64,     // Mean base quality
    pub qstd: i32,     // Quality std flag (0 or 1)
    pub mapq: f64,     // Mean mapping quality
    pub qratio: f64,   // High/low quality ratio
    pub hifreq: f64,   // High-quality frequency
    pub extrafreq: f64,// Extra frequency

    // Special metrics
    pub shift3: i32,
    pub msi: f64,
    pub msint: f64,
    pub nm: f64,       // Number of mismatches
    pub hicnt: usize,  // High-quality count
    pub hicov: usize,  // High-quality coverage

    // Context
    pub left_sequence: String,
    pub right_sequence: String,
    pub region: String,
    pub var_type: String,
    pub duprate: f64,
    pub sv: String,
}

impl SimpleOutputVariant {
    /// Create a SimpleOutputVariant from a Variant and Region
    pub fn from_variant(variant: &Variant, region: &Region, sample: &str, sv: &str) -> Self {
        let var_type_str = var_type_string(&variant.refallele, &variant.varallele);
        
        // Detect reference call (ref == alt)
        let is_ref_call = variant.refallele == variant.varallele;
        
        // For reference calls: bias should be "0;0" (no bias for ref or var)
        // For variants: bias is "ref_bias;var_bias" format
        let bias = if is_ref_call {
            "0;0".to_string()
        } else {
            format_strand_bias(variant.strand_bias_flag)
        };
        
        // For reference calls, counts go to ref_fwd/ref_rev, not var_fwd/var_rev
        let (variant_coverage, ref_fwd, ref_rev, var_fwd, var_rev, frequency) = if is_ref_call {
            (
                0,  // variant_coverage = 0 for ref calls
                variant.vars_count_on_forward,  // ref counts FROM vars counts
                variant.vars_count_on_reverse,
                0,  // var counts = 0
                0,
                0.0,  // frequency = 0 for ref calls
            )
        } else {
            (
                variant.vars_count_on_forward + variant.vars_count_on_reverse,  // Total variant count
                variant.ref_forward_count,  // Reference forward counts from same position
                variant.ref_reverse_count,  // Reference reverse counts from same position
                variant.vars_count_on_forward,
                variant.vars_count_on_reverse,
                variant.frequency,
            )
        };
        
        // For reference calls, vartype should be empty
        let final_var_type = if is_ref_call { String::new() } else { var_type_str };

        let chr = normalize_chr_for_output(&region.chr);
        SimpleOutputVariant {
            sample: sample.to_string(),
            gene: region.gene.clone(),
            chr: chr.clone(),
            start_position: variant.start_position,
            end_position: variant.end_position,
            ref_allele: variant.refallele.clone(),
            var_allele: variant.varallele.clone(),

            total_coverage: variant.position_coverage,
            variant_coverage,
            reference_forward_count: ref_fwd,
            reference_reverse_count: ref_rev,
            variant_forward_count: var_fwd,
            variant_reverse_count: var_rev,

            genotype: variant.genotype.clone(),
            frequency,
            bias,

            pmean: variant.mean_position,
            pstd: if variant.is_at_least_at_2_positions { 1 } else { 0 },
            qual: variant.mean_quality,
            qstd: if variant.has_at_least_2_diff_qualities { 1 } else { 0 },
            mapq: variant.mean_mapping_quality,
            // qratio: high_qual_read_cnt / low_qual_read_cnt (handle divide by zero)
            // For ref calls with no reads, qratio should be 0; otherwise calculate normally
            qratio: if variant.low_qual_read_cnt > 0 {
                variant.high_qual_read_cnt as f64 / variant.low_qual_read_cnt as f64
            } else if variant.high_qual_read_cnt > 0 {
                // All reads are high quality - use high_qual_read_cnt / 0.5 as per Java
                variant.high_qual_read_cnt as f64 * 2.0
            } else {
                0.0
            },
            hifreq: variant.high_quality_reads_frequency,
            extrafreq: variant.extra_frequency,

            shift3: variant.shift3,
            msi: variant.msi,
            msint: variant.msint,
            nm: variant.nm,
            hicnt: variant.high_qual_read_cnt,
            hicov: variant.hicov,

            left_sequence: if variant.leftseq.is_empty() { "0".to_string() } else { variant.leftseq.clone() },
            right_sequence: if variant.rightseq.is_empty() { "0".to_string() } else { variant.rightseq.clone() },
            region: format!("{}:{}-{}", chr, region.start, region.end),
            var_type: final_var_type,
            duprate: variant.duprate,
            sv: if sv.is_empty() { "0".to_string() } else { sv.to_string() },
        }
    }

    /// Create an empty variant (for positions with no variants)
    pub fn empty(position: i64, region: &Region, sample: &str) -> Self {
        let chr = normalize_chr_for_output(&region.chr);
        SimpleOutputVariant {
            sample: sample.to_string(),
            gene: region.gene.clone(),
            chr: chr.clone(),
            start_position: position,
            end_position: position,
            ref_allele: String::new(),
            var_allele: String::new(),

            total_coverage: 0,
            variant_coverage: 0,
            reference_forward_count: 0,
            reference_reverse_count: 0,
            variant_forward_count: 0,
            variant_reverse_count: 0,

            genotype: "0".to_string(),
            frequency: 0.0,
            bias: "0;0".to_string(),

            pmean: 0.0,
            pstd: 0,
            qual: 0.0,
            qstd: 0,
            mapq: 0.0,
            qratio: 0.0,
            hifreq: 0.0,
            extrafreq: 0.0,

            shift3: 0,
            msi: 0.0,
            msint: 0.0,
            nm: 0.0,
            hicnt: 0,
            hicov: 0,

            left_sequence: "0".to_string(),
            right_sequence: "0".to_string(),
            region: format!("{}:{}-{}", chr, region.start, region.end),
            var_type: String::new(),
            duprate: 0.0,
            sv: "0".to_string(),
        }
    }

    /// Format as 36-column tab-delimited string (Simple Mode without Fisher)
    pub fn to_string_36_columns(&self) -> String {
        let parts: Vec<String> = vec![
            // 1-7: Identity
            self.sample.clone(),
            self.gene.clone(),
            self.chr.clone(),
            self.start_position.to_string(),
            self.end_position.to_string(),
            self.ref_allele.clone(),
            self.var_allele.clone(),

            // 8-13: Coverage
            self.total_coverage.to_string(),
            self.variant_coverage.to_string(),
            self.reference_forward_count.to_string(),
            self.reference_reverse_count.to_string(),
            self.variant_forward_count.to_string(),
            self.variant_reverse_count.to_string(),

            // 14-16: Genotype, frequency, bias
            self.genotype.clone(),
            format_f64(self.frequency, 4),
            self.bias.clone(),

            // 17-18: Position metrics
            format_f64(self.pmean, 1),
            self.pstd.to_string(),

            // 19-20: Quality metrics
            format_f64(self.qual, 1),
            self.qstd.to_string(),

            // 21-24: More quality metrics
            format_f64(self.mapq, 1),
            format_f64(self.qratio, 3),
            format_f64(self.hifreq, 4),
            format_f64(self.extrafreq, 4),

            // 25-28: Special metrics
            self.shift3.to_string(),
            format_f64(self.msi, 3),
            format_f64(self.msint, 0),
            format_f64(self.nm, 1),

            // 29-30: High-quality counts
            self.hicnt.to_string(),
            self.hicov.to_string(),

            // 31-36: Context and type
            self.left_sequence.clone(),
            self.right_sequence.clone(),
            self.region.clone(),
            self.var_type.clone(),
            format_f64(self.duprate, 1),
            self.sv.clone(),
        ];

        parts.join("\t")
    }
}

impl std::fmt::Display for SimpleOutputVariant {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_string_36_columns())
    }
}

// ============================================================================
// Helper Functions
// ============================================================================

/// Format a float with specified decimal places, returning "0" for zero values
fn format_f64(value: f64, decimals: usize) -> String {
    if value == 0.0 {
        "0".to_string()
    } else {
        format!("{:.1$}", value, decimals)
    }
}

/// Format strand bias flag as "refBias;varBias" string
fn format_strand_bias(flag: StrandBiasFlag) -> String {
    flag.to_string()
}

/// Format VarType to string representation
fn format_var_type(var_type: &VarType) -> String {
    match var_type {
        VarType::SNV(_) => "SNV".to_string(),
        VarType::Insertion(_) => "Insertion".to_string(),
        VarType::Deletion(_) => "Deletion".to_string(),
        VarType::Complex { .. } => "Complex".to_string(),
    }
}

/// Get column headers for 36-column format
pub fn get_column_headers() -> Vec<&'static str> {
    vec![
        "Sample", "Gene", "Chr", "Start", "End", "Ref", "Alt",
        "Depth", "AltDepth", "RefFwd", "RefRev", "AltFwd", "AltRev",
        "Genotype", "AF", "Bias",
        "PMean", "PStd",
        "Qual", "QStd",
        "MQ", "QRatio", "HiFreq", "ExtraFreq",
        "Shift3", "MSI", "MSILen", "NM",
        "HiCnt", "HiCov",
        "LeftSeq", "RightSeq", "Region", "VarType", "DupRate", "SV"
    ]
}

/// Get column headers as tab-delimited string
pub fn get_header_line() -> String {
    get_column_headers().join("\t")
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_region_to_string() {
        let region = Region::new("chr1", 1000, 2000, "GENE1");
        assert_eq!(region.to_region_string(), "1:1000-2000");
    }

    #[test]
    fn test_format_f64_zero() {
        assert_eq!(format_f64(0.0, 4), "0");
    }

    #[test]
    fn test_format_f64_nonzero() {
        assert_eq!(format_f64(0.1234, 4), "0.1234");
        assert_eq!(format_f64(1.5, 1), "1.5");
        assert_eq!(format_f64(1.234, 2), "1.23");
    }

    #[test]
    fn test_format_strand_bias() {
        use crate::mods::to_vars_builder::StrandBiasValue;
        // Test with new struct format
        assert_eq!(format_strand_bias(StrandBiasFlag::new(StrandBiasValue::CantAssess, StrandBiasValue::CantAssess)), "0;0");
        assert_eq!(format_strand_bias(StrandBiasFlag::new(StrandBiasValue::NoBias, StrandBiasValue::NoBias)), "2;2");
        assert_eq!(format_strand_bias(StrandBiasFlag::new(StrandBiasValue::NoBias, StrandBiasValue::HasBias)), "2;1");
    }

    #[test]
    fn test_format_var_type() {
        assert_eq!(format_var_type(&VarType::SNV('A')), "SNV");
        assert_eq!(format_var_type(&VarType::Insertion("ATG".to_string())), "Insertion");
        assert_eq!(format_var_type(&VarType::Deletion(3)), "Deletion");
        assert_eq!(format_var_type(&VarType::Complex { insertion: "A".to_string(), deletion: 2 }), "Complex");
    }

    #[test]
    fn test_empty_variant() {
        let region = Region::new("chr1", 1000, 2000, "GENE1");
        let output = SimpleOutputVariant::empty(1500, &region, "sample1");

        assert_eq!(output.sample, "sample1");
        assert_eq!(output.start_position, 1500);
        assert_eq!(output.genotype, "0");
        assert_eq!(output.bias, "0;0");
        assert_eq!(output.sv, "0");
    }

    #[test]
    fn test_from_variant() {
        use crate::mods::to_vars_builder::StrandBiasValue;
        let variant = Variant {
            description_string: "A>T".to_string(),
            refallele: "A".to_string(),
            varallele: "T".to_string(),
            vartype: VarType::SNV('T'),
            start_position: 1000,
            end_position: 1000,
            vars_count_on_forward: 5,
            vars_count_on_reverse: 5,
            position_coverage: 100,
            frequency: 0.10,
            high_quality_reads_frequency: 0.08,
            extra_frequency: 0.0,
            mean_position: 25.0,
            mean_quality: 30.0,
            mean_mapping_quality: 60.0,
            strand_bias_flag: StrandBiasFlag::new(StrandBiasValue::NoBias, StrandBiasValue::NoBias),
            is_at_least_at_2_positions: true,
            has_at_least_2_diff_qualities: true,
            leftseq: "ACGT".to_string(),
            rightseq: "TGCA".to_string(),
            msi: 0.0,
            msint: 0.0,
            shift3: 0,
            nm: 1.0,
            high_qual_read_cnt: 10,
            low_qual_read_cnt: 0,
            hicov: 90,
            ref_forward_count: 0,
            ref_reverse_count: 0,
            genotype: "0/1".to_string(),
            duprate: 0.0,
        };

        let region = Region::new("chr1", 900, 1100, "BRCA1");
        let output = SimpleOutputVariant::from_variant(&variant, &region, "sample1", "");

        assert_eq!(output.sample, "sample1");
        assert_eq!(output.gene, "BRCA1");
        assert_eq!(output.chr, "1");
        assert_eq!(output.start_position, 1000);
        assert_eq!(output.ref_allele, "A");
        assert_eq!(output.var_allele, "T");
        assert_eq!(output.variant_forward_count, 5);
        assert_eq!(output.variant_reverse_count, 5);
        assert_eq!(output.genotype, "0/1");
        assert!((output.frequency - 0.10).abs() < 0.001);
        assert_eq!(output.bias, "2;2");  // NoBias for both ref and var
        assert_eq!(output.pstd, 1);
        assert_eq!(output.qstd, 1);
        assert_eq!(output.left_sequence, "ACGT");
        assert_eq!(output.right_sequence, "TGCA");
        assert_eq!(output.var_type, "SNV");
    }

    #[test]
    fn test_to_string_36_columns() {
        let region = Region::new("chr1", 1000, 2000, "GENE1");
        let mut output = SimpleOutputVariant::empty(1500, &region, "sample1");
        output.ref_allele = "A".to_string();
        output.var_allele = "T".to_string();
        output.total_coverage = 100;
        output.variant_coverage = 10;
        output.variant_forward_count = 5;
        output.variant_reverse_count = 5;
        output.genotype = "0/1".to_string();
        output.frequency = 0.10;
        output.var_type = "SNV".to_string();

        let line = output.to_string_36_columns();
        let fields: Vec<&str> = line.split('\t').collect();

        // Should have exactly 36 columns
        assert_eq!(fields.len(), 36);

        // Check key fields
        assert_eq!(fields[0], "sample1");          // Sample
        assert_eq!(fields[1], "GENE1");            // Gene
        assert_eq!(fields[2], "1");                // Chr
        assert_eq!(fields[3], "1500");             // Start
        assert_eq!(fields[5], "A");                // Ref
        assert_eq!(fields[6], "T");                // Alt
        assert_eq!(fields[7], "100");              // Total coverage
        assert_eq!(fields[13], "0/1");             // Genotype
        assert_eq!(fields[14], "0.1000");          // Frequency
        assert_eq!(fields[33], "SNV");             // VarType
    }

    #[test]
    fn test_get_column_headers() {
        let headers = get_column_headers();
        assert_eq!(headers.len(), 36);
        assert_eq!(headers[0], "Sample");
        assert_eq!(headers[13], "Genotype");
        assert_eq!(headers[33], "VarType");
    }

    #[test]
    fn test_get_header_line() {
        let header = get_header_line();
        let fields: Vec<&str> = header.split('\t').collect();
        assert_eq!(fields.len(), 36);
    }
}
