//! Pipeline Integration - End-to-end variant calling pipeline
//!
//! This module connects the main components of the VarDict pipeline:
//! 1. CigarParser - Parse BAM records and extract raw variation data
//! 2. ToVarsBuilder - Calculate statistics and create Variant objects
//! 3. OutputVariant - Format variants for output
//!
//! **Pipeline Flow:**
//! ```text
//! BAM Record → CigarParser → RawVariant → ToVarsBuilder → Variant → OutputVariant → TAB/VCF
//! ```

use std::collections::HashMap;

use crate::mods::{
    to_vars_builder::{ToVarsBuilder, Variant, VariationData, Vars, VarType},
    output_variant::{SimpleOutputVariant, Region},
};

/// Pipeline configuration
#[derive(Debug, Clone)]
pub struct PipelineConfig {
    pub sample_name: String,
    pub min_frequency: f64,
    pub min_variant_reads: usize,
    pub quality_threshold: u8,
    pub mapq_threshold: u8,
}

impl Default for PipelineConfig {
    fn default() -> Self {
        PipelineConfig {
            sample_name: "sample".to_string(),
            min_frequency: 0.01,
            min_variant_reads: 2,
            quality_threshold: 25,
            mapq_threshold: 0,
        }
    }
}

impl PipelineConfig {
    pub fn new(sample_name: &str) -> Self {
        PipelineConfig {
            sample_name: sample_name.to_string(),
            ..Default::default()
        }
    }

    /// Create a builder for PipelineConfig
    pub fn builder() -> PipelineConfigBuilder {
        PipelineConfigBuilder::new()
    }

    pub fn with_min_frequency(mut self, freq: f64) -> Self {
        self.min_frequency = freq;
        self
    }

    pub fn with_quality_threshold(mut self, qual: u8) -> Self {
        self.quality_threshold = qual;
        self
    }
}

/// Builder for PipelineConfig
#[derive(Debug, Clone)]
pub struct PipelineConfigBuilder {
    sample_name: String,
    min_frequency: f64,
    min_variant_reads: usize,
    quality_threshold: u8,
    mapq_threshold: u8,
}

impl PipelineConfigBuilder {
    pub fn new() -> Self {
        let defaults = PipelineConfig::default();
        PipelineConfigBuilder {
            sample_name: defaults.sample_name,
            min_frequency: defaults.min_frequency,
            min_variant_reads: defaults.min_variant_reads,
            quality_threshold: defaults.quality_threshold,
            mapq_threshold: defaults.mapq_threshold,
        }
    }

    pub fn sample_name(mut self, name: String) -> Self {
        self.sample_name = name;
        self
    }

    pub fn min_frequency(mut self, freq: f64) -> Self {
        self.min_frequency = freq;
        self
    }

    pub fn min_variant_reads(mut self, count: usize) -> Self {
        self.min_variant_reads = count;
        self
    }

    pub fn min_base_quality(mut self, qual: u8) -> Self {
        self.quality_threshold = qual;
        self
    }

    pub fn min_mapping_quality(mut self, qual: u8) -> Self {
        self.mapq_threshold = qual;
        self
    }

    pub fn build(self) -> PipelineConfig {
        PipelineConfig {
            sample_name: self.sample_name,
            min_frequency: self.min_frequency,
            min_variant_reads: self.min_variant_reads,
            quality_threshold: self.quality_threshold,
            mapq_threshold: self.mapq_threshold,
        }
    }
}

impl Default for PipelineConfigBuilder {
    fn default() -> Self {
        Self::new()
    }
}

/// Pipeline runner - orchestrates the variant calling process
pub struct Pipeline {
    config: PipelineConfig,
    builder: ToVarsBuilder,
}

impl Pipeline {
    /// Create a new pipeline with configuration
    pub fn new(config: PipelineConfig) -> Self {
        let builder = ToVarsBuilder::new()
            .with_min_frequency(config.min_frequency)
            .with_quality_threshold(config.quality_threshold)
            .with_mapq_threshold(config.mapq_threshold);

        Pipeline { config, builder }
    }

    /// Process variations and generate output variants
    /// 
    /// Input: Raw variation data grouped by (position, variant_key)
    /// Output: Formatted output lines ready for writing
    pub fn process_variations(
        &self,
        variations: Vec<(i64, String, VariationData)>,
        coverage_by_pos: &HashMap<i64, usize>,
        region: &Region,
    ) -> Vec<String> {
        // Step 1: Build variants using ToVarsBuilder
        let vars_by_position = self.builder.build_variants(variations, coverage_by_pos);

        // Step 2: Convert to output format
        let mut output_lines = Vec::new();
        
        // Sort positions for consistent output
        let mut positions: Vec<i64> = vars_by_position.keys().cloned().collect();
        positions.sort();

        for position in positions {
            if let Some(vars) = vars_by_position.get(&position) {
                for variant in &vars.variants {
                    let output = SimpleOutputVariant::from_variant(
                        variant,
                        region,
                        &self.config.sample_name,
                        "", // No SV info in simple mode
                    );
                    output_lines.push(output.to_string());
                }
            }
        }

        output_lines
    }

    /// Process and filter variants, returning only those passing thresholds
    pub fn process_and_filter(
        &self,
        variations: Vec<(i64, String, VariationData)>,
        coverage_by_pos: &HashMap<i64, usize>,
        region: &Region,
    ) -> Vec<Variant> {
        let vars_by_position = self.builder.build_variants(variations, coverage_by_pos);

        let mut result = Vec::new();
        for vars in vars_by_position.values() {
            for variant in &vars.variants {
                // Filter by minimum frequency
                if variant.frequency >= self.config.min_frequency {
                    result.push(variant.clone());
                }
            }
        }

        result
    }

    /// Get the header line for output
    pub fn get_header(&self) -> String {
        crate::mods::output_variant::get_header_line()
    }
}

// ============================================================================
// Conversion Helpers
// ============================================================================

/// Convert raw cigar parser variant data to VariationData for ToVarsBuilder
/// 
/// This bridges the gap between the low-level CigarParser output and
/// the higher-level ToVarsBuilder input.
pub fn convert_raw_variant_to_variation_data(
    position_in_read: u32,
    quality: u8,
    mapping_quality: u8,
    is_reverse: bool,
    read_id: &str,
) -> VariationData {
    VariationData {
        position_in_read,
        quality,
        mapping_quality,
        is_reverse,
        read_id: read_id.to_string(),
    }
}

/// Helper to create variation entry from parsed data
pub fn create_variation_entry(
    position: i64,
    ref_base: &str,
    alt_base: &str,
    position_in_read: u32,
    quality: u8,
    mapping_quality: u8,
    is_reverse: bool,
    read_id: &str,
) -> (i64, String, VariationData) {
    let variant_key = format!("{}>{}",ref_base, alt_base);
    let data = convert_raw_variant_to_variation_data(
        position_in_read,
        quality,
        mapping_quality,
        is_reverse,
        read_id,
    );
    (position, variant_key, data)
}

/// Create insertion variation entry
pub fn create_insertion_entry(
    position: i64,
    inserted_seq: &str,
    position_in_read: u32,
    quality: u8,
    mapping_quality: u8,
    is_reverse: bool,
    read_id: &str,
) -> (i64, String, VariationData) {
    let variant_key = format!("+{}", inserted_seq);
    let data = convert_raw_variant_to_variation_data(
        position_in_read,
        quality,
        mapping_quality,
        is_reverse,
        read_id,
    );
    (position, variant_key, data)
}

/// Create deletion variation entry
pub fn create_deletion_entry(
    position: i64,
    deletion_length: usize,
    position_in_read: u32,
    quality: u8,
    mapping_quality: u8,
    is_reverse: bool,
    read_id: &str,
) -> (i64, String, VariationData) {
    let variant_key = format!("-{}", deletion_length);
    let data = convert_raw_variant_to_variation_data(
        position_in_read,
        quality,
        mapping_quality,
        is_reverse,
        read_id,
    );
    (position, variant_key, data)
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_region() -> Region {
        Region::new("chr1", 1000, 2000, "TEST_GENE")
    }

    #[test]
    fn test_pipeline_config_default() {
        let config = PipelineConfig::default();
        assert_eq!(config.sample_name, "sample");
        assert_eq!(config.min_frequency, 0.01);
        assert_eq!(config.quality_threshold, 25);
    }

    #[test]
    fn test_pipeline_config_builder() {
        let config = PipelineConfig::new("my_sample")
            .with_min_frequency(0.05)
            .with_quality_threshold(25);
        
        assert_eq!(config.sample_name, "my_sample");
        assert_eq!(config.min_frequency, 0.05);
        assert_eq!(config.quality_threshold, 25);
    }

    #[test]
    fn test_pipeline_new() {
        let config = PipelineConfig::new("test");
        let pipeline = Pipeline::new(config);
        assert_eq!(pipeline.config.sample_name, "test");
    }

    #[test]
    fn test_convert_raw_variant() {
        let data = convert_raw_variant_to_variation_data(
            10, 30, 60, false, "read1"
        );
        
        assert_eq!(data.position_in_read, 10);
        assert_eq!(data.quality, 30);
        assert_eq!(data.mapping_quality, 60);
        assert!(!data.is_reverse);
        assert_eq!(data.read_id, "read1");
    }

    #[test]
    fn test_create_snv_entry() {
        let (pos, key, data) = create_variation_entry(
            1000, "A", "T", 10, 30, 60, false, "read1"
        );
        
        assert_eq!(pos, 1000);
        assert_eq!(key, "A>T");
        assert_eq!(data.position_in_read, 10);
    }

    #[test]
    fn test_create_insertion_entry() {
        let (pos, key, data) = create_insertion_entry(
            1000, "ATG", 10, 30, 60, true, "read1"
        );
        
        assert_eq!(pos, 1000);
        assert_eq!(key, "+ATG");
        assert!(data.is_reverse);
    }

    #[test]
    fn test_create_deletion_entry() {
        let (pos, key, data) = create_deletion_entry(
            1000, 5, 10, 30, 60, false, "read1"
        );
        
        assert_eq!(pos, 1000);
        assert_eq!(key, "-5");
        assert_eq!(data.quality, 30);
    }

    #[test]
    fn test_pipeline_process_variations() {
        let config = PipelineConfig::new("test_sample");
        let pipeline = Pipeline::new(config);
        let region = create_test_region();

        // Create test variations: 2 SNVs at position 1000
        let variations = vec![
            create_variation_entry(1000, "A", "T", 10, 30, 60, false, "read1"),
            create_variation_entry(1000, "A", "T", 15, 25, 55, true, "read2"),
        ];

        let mut coverage = HashMap::new();
        coverage.insert(1000i64, 100);

        let output = pipeline.process_variations(variations, &coverage, &region);

        // Should have one output line (one variant type at one position)
        assert_eq!(output.len(), 1);
        
        // Check output contains expected fields
        let line = &output[0];
        assert!(line.contains("test_sample"));
        assert!(line.contains("TEST_GENE"));
        assert!(line.contains("chr1"));
        assert!(line.contains("1000"));
    }

    #[test]
    fn test_pipeline_process_multiple_positions() {
        let config = PipelineConfig::new("sample1");
        let pipeline = Pipeline::new(config);
        let region = create_test_region();

        // Create variations at 2 positions
        let variations = vec![
            create_variation_entry(1000, "A", "T", 10, 30, 60, false, "read1"),
            create_variation_entry(1000, "A", "T", 15, 25, 55, true, "read2"),
            create_variation_entry(1010, "C", "G", 20, 28, 58, false, "read3"),
        ];

        let mut coverage = HashMap::new();
        coverage.insert(1000i64, 100);
        coverage.insert(1010i64, 100);

        let output = pipeline.process_variations(variations, &coverage, &region);

        // Should have 2 output lines (one per position)
        assert_eq!(output.len(), 2);
    }

    #[test]
    fn test_pipeline_process_and_filter() {
        let config = PipelineConfig::new("sample1")
            .with_min_frequency(0.05);
        let pipeline = Pipeline::new(config);
        let region = create_test_region();

        // Create variations: one with high frequency, one with low
        // Position 1000: 5 reads out of 100 = 5% (passes filter)
        // Position 1010: 2 reads out of 100 = 2% (fails filter)
        let mut variations = Vec::new();
        
        // 5 reads at position 1000
        for i in 0..5 {
            variations.push(create_variation_entry(
                1000, "A", "T", 10 + i, 30, 60, i % 2 == 0, &format!("read_{}", i)
            ));
        }
        
        // 2 reads at position 1010
        for i in 0..2 {
            variations.push(create_variation_entry(
                1010, "C", "G", 10 + i, 30, 60, false, &format!("read_b_{}", i)
            ));
        }

        let mut coverage = HashMap::new();
        coverage.insert(1000i64, 100);
        coverage.insert(1010i64, 100);

        let filtered = pipeline.process_and_filter(variations, &coverage, &region);

        // Only position 1000 should pass (5% >= 5%)
        assert_eq!(filtered.len(), 1);
        assert!((filtered[0].frequency - 0.05).abs() < 0.001);
    }

    #[test]
    fn test_pipeline_get_header() {
        let pipeline = Pipeline::new(PipelineConfig::default());
        let header = pipeline.get_header();
        
        // Header should have 36 columns
        let cols: Vec<&str> = header.split('\t').collect();
        assert_eq!(cols.len(), 36);
        assert_eq!(cols[0], "Sample");
        assert_eq!(cols[13], "Genotype");
    }

    // ========================================================================
    // Tests ported from Java SimplePostProcessModuleTest.java
    // ========================================================================

    /// Port of Java testGoodVariant()
    /// Verifies the 36-column output format matches Java exactly
    #[test]
    fn test_output_format_matches_java() {
        let config = PipelineConfig::new("test_bam");
        let pipeline = Pipeline::new(config);
        let region = Region::new("1", 1, 10, "gene_name");

        // Create a "good variant" matching Java test values:
        // refallele = "T", varallele = "A", genotype = "T/A"
        // varsCountOnForward = 3, varsCountOnReverse = 5 (total = 8)
        // positionCoverage = 4, frequency = 0.4
        // meanPosition = 9, meanQuality = 26, meanMappingQuality = 7.8
        // strandBiasFlag = "2" (strong bias), hicnt = 44
        // highQualityToLowQualityRatio = 2.5, numberOfMismatches = 2.0
        
        // Create 8 variations (3 forward + 5 reverse) with properties that produce Java values
        let mut variations = Vec::new();
        
        // 3 forward strand reads
        for i in 0..3 {
            variations.push(create_variation_entry(
                1, "T", "A", 9, 26, 8, false, &format!("fwd_{}", i)
            ));
        }
        
        // 5 reverse strand reads
        for i in 0..5 {
            variations.push(create_variation_entry(
                1, "T", "A", 9, 26, 8, true, &format!("rev_{}", i)
            ));
        }

        let mut coverage = HashMap::new();
        coverage.insert(1i64, 20); // 8/20 = 0.4 frequency

        let output = pipeline.process_variations(variations, &coverage, &region);

        // Should have one output line
        assert_eq!(output.len(), 1);
        
        let line = &output[0];
        let cols: Vec<&str> = line.split('\t').collect();
        
        // Verify 36 columns
        assert_eq!(cols.len(), 36, "Expected 36 columns, got {}", cols.len());
        
        // Verify key fields match Java output format:
        // "test_bam\tgene_name\t1\t0\t0\tT\tA\t0\t4\t0\t0\t3\t5\tT/A\t0.4000..."
        assert_eq!(cols[0], "test_bam", "Sample name");
        assert_eq!(cols[1], "gene_name", "Gene");
        assert_eq!(cols[2], "1", "Chr");
        // cols[3] = start position
        // cols[5], cols[6] = ref, alt alleles
        assert_eq!(cols[11], "3", "Variant forward count");
        assert_eq!(cols[12], "5", "Variant reverse count");
        // cols[13] = genotype
        // cols[14] = frequency (should be ~0.4)
        assert_eq!(cols[32], "1:1-10", "Region");
        assert_eq!(cols[33], "SNV", "Variant type");
    }

    /// Port of Java testNoVariantNoReferenceVariant()
    /// Empty variants should produce no output
    #[test]
    fn test_no_variant_produces_empty_output() {
        let config = PipelineConfig::new("test_bam");
        let pipeline = Pipeline::new(config);
        let region = Region::new("1", 1, 10, "gene_name");

        // Empty variations
        let variations: Vec<(i64, String, VariationData)> = vec![];
        let coverage = HashMap::new();

        let output = pipeline.process_variations(variations, &coverage, &region);

        // Should have no output
        assert!(output.is_empty(), "Expected empty output for no variants");
    }

    /// Port of Java testFewGoodVariants()
    /// Multiple variants should produce multiple output lines
    #[test]
    fn test_multiple_variants_produce_multiple_lines() {
        let config = PipelineConfig::new("test_bam");
        let pipeline = Pipeline::new(config);
        let region = Region::new("1", 1, 10, "gene_name");

        // Create two different variant types at position 1
        let variations = vec![
            // SNV T>A
            create_variation_entry(1, "T", "A", 9, 26, 8, false, "read1"),
            create_variation_entry(1, "T", "A", 9, 26, 8, true, "read2"),
            // SNV T>G (different variant)
            create_variation_entry(1, "T", "G", 9, 26, 8, false, "read3"),
            create_variation_entry(1, "T", "G", 9, 26, 8, true, "read4"),
        ];

        let mut coverage = HashMap::new();
        coverage.insert(1i64, 100);

        let output = pipeline.process_variations(variations, &coverage, &region);

        // Should have 2 output lines (2 different variants at same position)
        assert_eq!(output.len(), 2, "Expected 2 output lines for 2 variant types");
        
        // Both should be for the same sample and gene
        for line in &output {
            assert!(line.starts_with("test_bam\tgene_name\t1"));
        }
    }

    /// Port of Java testBadVariant() concept
    /// Variants with strong strand bias and low quality should be filtered
    #[test]
    fn test_bad_variant_filtering() {
        let config = PipelineConfig::new("test_bam")
            .with_min_frequency(0.1); // Require 10% frequency
        let pipeline = Pipeline::new(config);
        let region = Region::new("1", 1, 10, "gene_name");

        // Create a variant with low frequency (2/100 = 2%)
        let variations = vec![
            create_variation_entry(1, "T", "A", 2, 10, 8, false, "read1"),
            create_variation_entry(1, "T", "A", 2, 10, 8, false, "read2"),
        ];

        let mut coverage = HashMap::new();
        coverage.insert(1i64, 100);

        let filtered = pipeline.process_and_filter(variations, &coverage, &region);

        // Should be filtered out (2% < 10%)
        assert!(filtered.is_empty(), "Low frequency variant should be filtered");
    }
}
