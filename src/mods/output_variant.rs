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

use crate::mods::to_vars_builder::{StrandBiasFlag, VarType, Variant, var_type_string};
use crate::scopedata::global_read_only_scope::{INSTANCE, instance};
use crate::utils::round_half_even;
use statrs::distribution::{Discrete, DiscreteCDF, Hypergeometric};

/// Region information for output
#[derive(Debug, Clone, Default)]
pub struct Region {
    pub chr: String,
    pub start: i64,
    pub end: i64,
    pub display_start: i64,
    pub gene: String,
}

impl Region {
    pub fn new(chr: &str, start: i64, end: i64, gene: &str) -> Self {
        Region {
            chr: chr.to_string(),
            start,
            end,
            display_start: start,
            gene: gene.to_string(),
        }
    }

    /// Format as "chr:start-end"
    pub fn to_region_string(&self) -> String {
        format!("{}:{}-{}", self.chr, self.display_start, self.end)
    }
}

/// Preserve chromosome name for output (Java uses region.chr directly)
fn normalize_chr_for_output(chr: &str) -> String {
    chr.to_string()
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
    pub bias: String, // "flag;flag" format

    // Position metrics
    pub pmean: f64, // Mean position in read
    pub pstd: i32,  // Position std flag (0 or 1)

    // Quality metrics
    pub qual: f64,      // Mean base quality
    pub qstd: i32,      // Quality std flag (0 or 1)
    pub mapq: f64,      // Mean mapping quality
    pub qratio: f64,    // High/low quality ratio
    pub hifreq: f64,    // High-quality frequency
    pub extrafreq: f64, // Extra frequency

    // Special metrics
    pub shift3: i32,
    pub msi: f64,
    pub msint: f64,
    pub nm: f64,      // Number of mismatches
    pub hicnt: usize, // High-quality count
    pub hicov: usize, // High-quality coverage

    // Context
    pub left_sequence: String,
    pub right_sequence: String,
    pub region: String,
    pub var_type: String,
    pub duprate: f64,
    pub crispr: i32,
    pub sv: String,
    pub debug: String,
}

impl SimpleOutputVariant {
    /// Create a SimpleOutputVariant from a Variant and Region
    pub fn from_variant(variant: &Variant, region: &Region, sample: &str, sv: &str) -> Self {
        let var_type_str = format_output_var_type(variant);

        // Detect reference call (ref == alt)
        let is_ref_call = variant.refallele == variant.varallele;

        // Bias is "ref_bias;var_bias" format for all calls (Java uses variant.strandBiasFlag)
        let bias = format_strand_bias(variant.strand_bias_flag);

        // For reference calls, counts go to ref_fwd/ref_rev, not var_fwd/var_rev
        let (variant_coverage, ref_fwd, ref_rev, var_fwd, var_rev, frequency) = if is_ref_call {
            (
                0, // variant_coverage = 0 for ref calls
                variant.ref_forward_count,
                variant.ref_reverse_count,
                0, // var counts = 0
                0,
                0.0, // frequency = 0 for ref calls
            )
        } else {
            (
                variant.position_coverage,
                variant.ref_forward_count, // Reference forward counts from same position
                variant.ref_reverse_count, // Reference reverse counts from same position
                variant.vars_count_on_forward,
                variant.vars_count_on_reverse,
                variant.frequency,
            )
        };

        // For reference calls, vartype should be empty
        let final_var_type = if is_ref_call {
            String::new()
        } else {
            var_type_str
        };

        let chr = normalize_chr_for_output(&region.chr);
        SimpleOutputVariant {
            sample: sample.to_string(),
            gene: region.gene.clone(),
            chr: chr.clone(),
            start_position: variant.start_position,
            end_position: variant.end_position,
            ref_allele: variant.refallele.clone(),
            var_allele: variant.varallele.clone(),

            total_coverage: variant.total_pos_coverage,
            variant_coverage,
            reference_forward_count: ref_fwd,
            reference_reverse_count: ref_rev,
            variant_forward_count: var_fwd,
            variant_reverse_count: var_rev,

            genotype: variant.genotype.clone(),
            frequency,
            bias,

            pmean: variant.mean_position,
            pstd: if variant.is_at_least_at_2_positions {
                1
            } else {
                0
            },
            qual: variant.mean_quality,
            qstd: if variant.has_at_least_2_diff_qualities {
                1
            } else {
                0
            },
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
            nm: if variant.nm > 0.0 { variant.nm } else { 0.0 },
            hicnt: variant.high_qual_read_cnt,
            hicov: variant.hicov,

            left_sequence: if variant.leftseq.is_empty() {
                "0".to_string()
            } else {
                variant.leftseq.clone()
            },
            right_sequence: if variant.rightseq.is_empty() {
                "0".to_string()
            } else {
                variant.rightseq.clone()
            },
            region: format!("{}:{}-{}", chr, region.display_start, region.end),
            var_type: final_var_type,
            duprate: variant.duprate,
            crispr: variant.crispr,
            sv: if sv.is_empty() {
                "0".to_string()
            } else {
                sv.to_string()
            },
            debug: if INSTANCE
                .get()
                .map(|scope| scope.conf.debug)
                .unwrap_or(false)
            {
                format_variant_debug_content(variant)
            } else {
                String::new()
            },
        }
    }

    /// Create an empty position-0 sentinel row.
    ///
    /// Java creates these with a default-constructed Variant (non-null), so all fields use
    /// Variant defaults: genotype=null→"0", strandBiasFlag="0", leftseq=null→"0", rightseq=null→"0".
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
            bias: "0".to_string(),

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
            region: format!("{}:{}-{}", chr, region.display_start, region.end),
            var_type: String::new(),
            duprate: 0.0,
            crispr: 0,
            sv: "0".to_string(),
            debug: String::new(),
        }
    }

    /// Create an empty variant for a zero-coverage real position (null-variant path).
    ///
    /// Java creates these via `new SimpleOutputVariant(null, region, sv, position)` when
    /// referenceVariant is null and variants is empty. Uses SimpleOutputVariant field defaults:
    /// genotype="" (empty), bias="0;0", leftseq="" (empty), rightseq="" (empty).
    pub fn empty_null_variant(position: i64, region: &Region, sample: &str) -> Self {
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

            genotype: String::new(),
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

            left_sequence: String::new(),
            right_sequence: String::new(),
            region: format!("{}:{}-{}", chr, region.display_start, region.end),
            var_type: String::new(),
            duprate: 0.0,
            crispr: 0,
            sv: "0".to_string(),
            debug: String::new(),
        }
    }

    /// Create an empty variant and preserve Java-compatible SV column value
    pub fn empty_with_sv(position: i64, region: &Region, sample: &str, sv: &str) -> Self {
        let mut out = Self::empty(position, region, sample);
        if !sv.is_empty() {
            out.sv = sv.to_string();
        }
        out
    }

    /// Create a null-variant empty row with SV field
    pub fn empty_null_variant_with_sv(
        position: i64,
        region: &Region,
        sample: &str,
        sv: &str,
    ) -> Self {
        let mut out = Self::empty_null_variant(position, region, sample);
        if !sv.is_empty() {
            out.sv = sv.to_string();
        }
        out
    }

    /// Format as 36-column tab-delimited string (Simple Mode without Fisher)
    pub fn to_string_36_columns(&self) -> String {
        let parts: Vec<String> = vec![
            self.sample.clone(),
            self.gene.clone(),
            self.chr.clone(),
            self.start_position.to_string(),
            self.end_position.to_string(),
            self.ref_allele.clone(),
            self.var_allele.clone(),
            self.total_coverage.to_string(),
            self.variant_coverage.to_string(),
            self.reference_forward_count.to_string(),
            self.reference_reverse_count.to_string(),
            self.variant_forward_count.to_string(),
            self.variant_reverse_count.to_string(),
            self.genotype.clone(),
            format_f64(self.frequency, 4),
            self.bias.clone(),
            format_f64(self.pmean, 1),
            self.pstd.to_string(),
            format_f64(self.qual, 1),
            self.qstd.to_string(),
            format_f64(self.mapq, 1),
            format_f64(self.qratio, 3),
            format_f64(self.hifreq, 4),
            format_f64(self.extrafreq, 4),
            self.shift3.to_string(),
            format_f64(self.msi, 3),
            format_f64(self.msint, 0),
            format_f64(self.nm, 1),
            self.hicnt.to_string(),
            self.hicov.to_string(),
            self.left_sequence.clone(),
            self.right_sequence.clone(),
            self.region.clone(),
            self.var_type.clone(),
            format_f64(self.duprate, 1),
            self.sv.clone(),
        ];

        parts.join("\t")
    }

    fn to_string_38_columns(&self) -> String {
        let fisher = FisherExact::new(
            self.reference_forward_count,
            self.reference_reverse_count,
            self.variant_forward_count,
            self.variant_reverse_count,
        );
        let pvalue = fisher.p_value();
        let oddratio = fisher.odd_ratio();

        let hifreq = if self.hifreq == 0.0 {
            "0".to_string()
        } else {
            format!("{:.4}", self.hifreq)
        };

        let nm = if self.nm > 0.0 { self.nm } else { 0.0 };
        let nmf = if nm == 0.0 {
            "0".to_string()
        } else {
            format!("{:.1}", nm)
        };

        let parts: Vec<String> = vec![
            self.sample.clone(),
            self.gene.clone(),
            self.chr.clone(),
            self.start_position.to_string(),
            self.end_position.to_string(),
            self.ref_allele.clone(),
            self.var_allele.clone(),
            self.total_coverage.to_string(),
            self.variant_coverage.to_string(),
            self.reference_forward_count.to_string(),
            self.reference_reverse_count.to_string(),
            self.variant_forward_count.to_string(),
            self.variant_reverse_count.to_string(),
            self.genotype.clone(),
            format_rounded_value_to_print("0.0000", self.frequency),
            self.bias.clone(),
            format_rounded_value_to_print("0.0", self.pmean),
            self.pstd.to_string(),
            format_rounded_value_to_print("0.0", self.qual),
            self.qstd.to_string(),
            format_rounded_value_to_print("0.00000", pvalue),
            oddratio,
            format_rounded_value_to_print("0.0", self.mapq),
            format_rounded_value_to_print("0.000", self.qratio),
            hifreq,
            format_rounded_value_to_print("0.0000", self.extrafreq),
            self.shift3.to_string(),
            format_rounded_value_to_print("0.000", self.msi),
            format_f64(self.msint, 0),
            nmf,
            self.hicnt.to_string(),
            self.hicov.to_string(),
            self.left_sequence.clone(),
            self.right_sequence.clone(),
            self.region.clone(),
            self.var_type.clone(),
            format_rounded_value_to_print("0.00", self.duprate),
            self.sv.clone(),
        ];

        parts.join("\t")
    }

    fn to_string_with_flags_and_debug(
        &self,
        fisher_enabled: bool,
        crispr_enabled: bool,
        debug_enabled: bool,
    ) -> String {
        let mut output_variant = if fisher_enabled {
            self.to_string_38_columns()
        } else {
            self.to_string_36_columns()
        };

        if crispr_enabled {
            output_variant.push('\t');
            output_variant.push_str(&self.crispr.to_string());
        }

        if debug_enabled {
            output_variant.push('\t');
            output_variant.push_str(&self.debug);
        }

        output_variant
    }

    fn to_string_with_flags(&self, fisher_enabled: bool, crispr_enabled: bool) -> String {
        let debug_enabled = INSTANCE
            .get()
            .map(|scope| scope.conf.debug)
            .unwrap_or(false);
        self.to_string_with_flags_and_debug(fisher_enabled, crispr_enabled, debug_enabled)
    }
}

impl std::fmt::Display for SimpleOutputVariant {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let fisher_enabled = INSTANCE
            .get()
            .map(|scope| scope.conf.fisher)
            .unwrap_or(false);
        let crispr_enabled = INSTANCE
            .get()
            .map(|scope| scope.conf.crispr_cutting_site != 0)
            .unwrap_or(false);
        write!(
            f,
            "{}",
            self.to_string_with_flags(fisher_enabled, crispr_enabled)
        )
    }
}

#[derive(Debug, Clone)]
pub struct AmpliconOutputVariant {
    pub sample: String,
    pub gene: String,
    pub chr: String,
    pub start_position: i64,
    pub end_position: i64,
    pub ref_allele: String,
    pub var_allele: String,
    pub total_coverage: usize,
    pub variant_coverage: usize,
    pub reference_forward_count: usize,
    pub reference_reverse_count: usize,
    pub variant_forward_count: usize,
    pub variant_reverse_count: usize,
    pub genotype: String,
    pub frequency: f64,
    pub bias: String,
    pub pmean: f64,
    pub pstd: i32,
    pub qual: f64,
    pub qstd: i32,
    pub mapq: f64,
    pub qratio: f64,
    pub hifreq: f64,
    pub extrafreq: f64,
    pub shift3: i32,
    pub msi: f64,
    pub msint: f64,
    pub nm: f64,
    pub hicnt: usize,
    pub hicov: usize,
    pub left_sequence: String,
    pub right_sequence: String,
    pub region: String,
    pub var_type: String,
    pub good_variants_count: usize,
    pub total_variants_count: usize,
    pub no_coverage: usize,
    pub amplicon_flag: i32,
    pub debug: String,
}

impl AmpliconOutputVariant {
    #[allow(clippy::too_many_arguments)]
    pub fn from_variant(
        variant: Option<&Variant>,
        region: &Region,
        good_variants: &[(Variant, String)],
        bad_variants: &[(Option<Variant>, String)],
        debug_prefix: Option<&str>,
        position: i64,
        good_variants_count: usize,
        no_coverage: usize,
        amplicon_bias_flag: bool,
        sample: &str,
    ) -> Self {
        let chr = normalize_chr_for_output(&region.chr);
        let output_region = if let Some((_, reg)) = good_variants.first() {
            reg.clone()
        } else {
            format!("{}:{}-{}", region.chr, position, position)
        };

        match variant {
            Some(v) => AmpliconOutputVariant {
                sample: sample.to_string(),
                gene: region.gene.clone(),
                chr,
                start_position: v.start_position,
                end_position: v.end_position,
                ref_allele: v.refallele.clone(),
                var_allele: v.varallele.clone(),
                total_coverage: v.total_pos_coverage,
                variant_coverage: v.position_coverage,
                reference_forward_count: v.ref_forward_count,
                reference_reverse_count: v.ref_reverse_count,
                variant_forward_count: v.vars_count_on_forward,
                variant_reverse_count: v.vars_count_on_reverse,
                genotype: if v.genotype.is_empty() {
                    "0".to_string()
                } else {
                    v.genotype.clone()
                },
                frequency: v.frequency,
                bias: v.strand_bias_flag.to_string(),
                pmean: v.mean_position,
                pstd: if v.is_at_least_at_2_positions { 1 } else { 0 },
                qual: v.mean_quality,
                qstd: if v.has_at_least_2_diff_qualities {
                    1
                } else {
                    0
                },
                mapq: v.mean_mapping_quality,
                qratio: if v.low_qual_read_cnt > 0 {
                    v.high_qual_read_cnt as f64 / v.low_qual_read_cnt as f64
                } else if v.high_qual_read_cnt > 0 {
                    v.high_qual_read_cnt as f64 * 2.0
                } else {
                    0.0
                },
                hifreq: v.high_quality_reads_frequency,
                extrafreq: v.extra_frequency,
                shift3: v.shift3,
                msi: v.msi,
                msint: v.msint,
                nm: if v.nm > 0.0 { v.nm } else { 0.0 },
                hicnt: v.high_qual_read_cnt,
                hicov: v.hicov,
                left_sequence: if v.leftseq.is_empty() {
                    "0".to_string()
                } else {
                    v.leftseq.clone()
                },
                right_sequence: if v.rightseq.is_empty() {
                    "0".to_string()
                } else {
                    v.rightseq.clone()
                },
                region: output_region,
                var_type: var_type_string(&v.refallele, &v.varallele),
                good_variants_count,
                total_variants_count: good_variants_count + bad_variants.len(),
                no_coverage,
                amplicon_flag: if amplicon_bias_flag { 1 } else { 0 },
                debug: if instance().conf.debug {
                    build_amplicon_debug(v, good_variants, bad_variants, debug_prefix)
                } else {
                    String::new()
                },
            },
            None => AmpliconOutputVariant {
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
                genotype: String::new(),
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
                left_sequence: String::new(),
                right_sequence: String::new(),
                region: format!("{}:{}-{}", chr, position, position),
                var_type: String::new(),
                good_variants_count,
                total_variants_count: good_variants_count + bad_variants.len(),
                no_coverage,
                amplicon_flag: if amplicon_bias_flag { 1 } else { 0 },
                debug: String::new(),
            },
        }
    }

    pub fn to_string_38_columns(&self) -> String {
        let parts: Vec<String> = vec![
            self.sample.clone(),
            self.gene.clone(),
            self.chr.clone(),
            self.start_position.to_string(),
            self.end_position.to_string(),
            self.ref_allele.clone(),
            self.var_allele.clone(),
            self.total_coverage.to_string(),
            self.variant_coverage.to_string(),
            self.reference_forward_count.to_string(),
            self.reference_reverse_count.to_string(),
            self.variant_forward_count.to_string(),
            self.variant_reverse_count.to_string(),
            self.genotype.clone(),
            format_f64(self.frequency, 4),
            self.bias.clone(),
            format_f64(self.pmean, 1),
            self.pstd.to_string(),
            format_f64(self.qual, 1),
            self.qstd.to_string(),
            format_f64(self.mapq, 1),
            format_f64(self.qratio, 3),
            format_f64(self.hifreq, 4),
            format_f64(self.extrafreq, 4),
            self.shift3.to_string(),
            format_f64(self.msi, 3),
            format_f64(self.msint, 0),
            format_f64(self.nm, 1),
            self.hicnt.to_string(),
            self.hicov.to_string(),
            self.left_sequence.clone(),
            self.right_sequence.clone(),
            self.region.clone(),
            self.var_type.clone(),
            self.good_variants_count.to_string(),
            self.total_variants_count.to_string(),
            self.no_coverage.to_string(),
            self.amplicon_flag.to_string(),
        ];

        parts.join("\t")
    }

    fn to_string_40_columns(&self) -> String {
        let fisher = FisherExact::new(
            self.reference_forward_count,
            self.reference_reverse_count,
            self.variant_forward_count,
            self.variant_reverse_count,
        );
        let pvalue = fisher.p_value();
        let oddratio = fisher.odd_ratio();

        let hifreq = if self.hifreq == 0.0 {
            "0".to_string()
        } else {
            format!("{:.4}", self.hifreq)
        };

        let nm = if self.nm > 0.0 { self.nm } else { 0.0 };
        let nmf = if nm == 0.0 {
            "0".to_string()
        } else {
            format!("{:.1}", nm)
        };

        let parts: Vec<String> = vec![
            self.sample.clone(),
            self.gene.clone(),
            self.chr.clone(),
            self.start_position.to_string(),
            self.end_position.to_string(),
            self.ref_allele.clone(),
            self.var_allele.clone(),
            self.total_coverage.to_string(),
            self.variant_coverage.to_string(),
            self.reference_forward_count.to_string(),
            self.reference_reverse_count.to_string(),
            self.variant_forward_count.to_string(),
            self.variant_reverse_count.to_string(),
            self.genotype.clone(),
            format_rounded_value_to_print("0.0000", self.frequency),
            self.bias.clone(),
            format_rounded_value_to_print("0.0", self.pmean),
            self.pstd.to_string(),
            format_rounded_value_to_print("0.0", self.qual),
            self.qstd.to_string(),
            format_rounded_value_to_print("0.00000", pvalue),
            oddratio,
            format_rounded_value_to_print("0.0", self.mapq),
            format_rounded_value_to_print("0.000", self.qratio),
            hifreq,
            format_rounded_value_to_print("0.0000", self.extrafreq),
            self.shift3.to_string(),
            format_rounded_value_to_print("0.000", self.msi),
            format_f64(self.msint, 0),
            nmf,
            self.hicnt.to_string(),
            self.hicov.to_string(),
            self.left_sequence.clone(),
            self.right_sequence.clone(),
            self.region.clone(),
            self.var_type.clone(),
            self.good_variants_count.to_string(),
            self.total_variants_count.to_string(),
            self.no_coverage.to_string(),
            self.amplicon_flag.to_string(),
        ];

        parts.join("\t")
    }

    fn to_string_with_flags(&self, fisher_enabled: bool, debug_enabled: bool) -> String {
        let output_variant = if fisher_enabled {
            self.to_string_40_columns()
        } else {
            self.to_string_38_columns()
        };

        if debug_enabled {
            format!("{}\t{}", output_variant, self.debug)
        } else {
            output_variant
        }
    }
}

impl std::fmt::Display for AmpliconOutputVariant {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let fisher_enabled = INSTANCE
            .get()
            .map(|scope| scope.conf.fisher)
            .unwrap_or(false);
        let debug_enabled = INSTANCE
            .get()
            .map(|scope| scope.conf.debug)
            .unwrap_or(false);
        write!(
            f,
            "{}",
            self.to_string_with_flags(fisher_enabled, debug_enabled)
        )
    }
}

#[derive(Debug, Clone)]
pub struct SomaticOutputVariant {
    pub sample: String,
    pub gene: String,
    pub chr: String,
    pub start_position: i64,
    pub end_position: i64,
    pub ref_allele: String,
    pub var_allele: String,

    pub var1_total_coverage: usize,
    pub var1_variant_coverage: usize,
    pub var1_ref_forward_coverage: usize,
    pub var1_ref_reverse_coverage: usize,
    pub var1_variant_forward_count: usize,
    pub var1_variant_reverse_count: usize,
    pub var1_genotype: String,
    pub var1_frequency: f64,
    pub var1_strand_bias_flag: String,
    pub var1_mean_position: f64,
    pub var1_is_at_least_at_2_positions: i32,
    pub var1_mean_quality: f64,
    pub var1_has_at_least_2_diff_qualities: i32,
    pub var1_mean_mapping_quality: f64,
    pub var1_high_quality_to_low_quality_ratio: f64,
    pub var1_high_quality_reads_frequency: f64,
    pub var1_extra_frequency: f64,
    pub var1_nm: f64,
    pub var1_duprate: f64,
    pub var1_sv: String,

    pub var2_total_coverage: usize,
    pub var2_variant_coverage: usize,
    pub var2_ref_forward_coverage: usize,
    pub var2_ref_reverse_coverage: usize,
    pub var2_variant_forward_count: usize,
    pub var2_variant_reverse_count: usize,
    pub var2_genotype: String,
    pub var2_frequency: f64,
    pub var2_strand_bias_flag: String,
    pub var2_mean_position: f64,
    pub var2_is_at_least_at_2_positions: i32,
    pub var2_mean_quality: f64,
    pub var2_has_at_least_2_diff_qualities: i32,
    pub var2_mean_mapping_quality: f64,
    pub var2_high_quality_to_low_quality_ratio: f64,
    pub var2_high_quality_reads_frequency: f64,
    pub var2_extra_frequency: f64,
    pub var2_nm: f64,
    pub var2_duprate: f64,
    pub var2_sv: String,

    pub shift3: i32,
    pub msi: f64,
    pub msint: f64,
    pub left_sequence: String,
    pub right_sequence: String,
    pub region: String,
    pub var_label: String,
    pub var_type: String,
    pub debug: String,
}

impl SomaticOutputVariant {
    #[allow(clippy::too_many_arguments)]
    pub fn from_variants(
        begin_variant: Option<&Variant>,
        end_variant: Option<&Variant>,
        tumor_variant: Option<&Variant>,
        normal_variant: Option<&Variant>,
        region: &Region,
        sv1: &str,
        sv2: &str,
        var_label: &str,
        sample: &str,
    ) -> Self {
        let mut output = SomaticOutputVariant {
            sample: sample.to_string(),
            gene: region.gene.clone(),
            chr: normalize_chr_for_output(&region.chr),
            start_position: 0,
            end_position: 0,
            ref_allele: String::new(),
            var_allele: String::new(),

            var1_total_coverage: 0,
            var1_variant_coverage: 0,
            var1_ref_forward_coverage: 0,
            var1_ref_reverse_coverage: 0,
            var1_variant_forward_count: 0,
            var1_variant_reverse_count: 0,
            var1_genotype: "0".to_string(),
            var1_frequency: 0.0,
            var1_strand_bias_flag: "0".to_string(),
            var1_mean_position: 0.0,
            var1_is_at_least_at_2_positions: 0,
            var1_mean_quality: 0.0,
            var1_has_at_least_2_diff_qualities: 0,
            var1_mean_mapping_quality: 0.0,
            var1_high_quality_to_low_quality_ratio: 0.0,
            var1_high_quality_reads_frequency: 0.0,
            var1_extra_frequency: 0.0,
            var1_nm: 0.0,
            var1_duprate: 0.0,
            var1_sv: if sv1.is_empty() {
                "0".to_string()
            } else {
                sv1.to_string()
            },

            var2_total_coverage: 0,
            var2_variant_coverage: 0,
            var2_ref_forward_coverage: 0,
            var2_ref_reverse_coverage: 0,
            var2_variant_forward_count: 0,
            var2_variant_reverse_count: 0,
            var2_genotype: "0".to_string(),
            var2_frequency: 0.0,
            var2_strand_bias_flag: "0".to_string(),
            var2_mean_position: 0.0,
            var2_is_at_least_at_2_positions: 0,
            var2_mean_quality: 0.0,
            var2_has_at_least_2_diff_qualities: 0,
            var2_mean_mapping_quality: 0.0,
            var2_high_quality_to_low_quality_ratio: 0.0,
            var2_high_quality_reads_frequency: 0.0,
            var2_extra_frequency: 0.0,
            var2_nm: 0.0,
            var2_duprate: 0.0,
            var2_sv: if sv2.is_empty() {
                "0".to_string()
            } else {
                sv2.to_string()
            },

            shift3: 0,
            msi: 0.0,
            msint: 0.0,
            left_sequence: String::new(),
            right_sequence: String::new(),
            region: format!("{}:{}-{}", region.chr, region.display_start, region.end),
            var_label: var_label.to_string(),
            var_type: String::new(),
            debug: String::new(),
        };

        if let Some(begin_variant) = begin_variant {
            output.start_position = begin_variant.start_position;
            output.end_position = begin_variant.end_position;
            output.ref_allele = begin_variant.refallele.clone();
            output.var_allele = begin_variant.varallele.clone();
            output.var_type = var_type_string(&begin_variant.refallele, &begin_variant.varallele);
        }

        if let Some(end_variant) = end_variant {
            output.shift3 = end_variant.shift3;
            output.msi = end_variant.msi;
            output.msint = end_variant.msint;
            output.left_sequence = if end_variant.leftseq.is_empty() {
                "0".to_string()
            } else {
                end_variant.leftseq.clone()
            };
            output.right_sequence = if end_variant.rightseq.is_empty() {
                "0".to_string()
            } else {
                end_variant.rightseq.clone()
            };
        }

        if let Some(tumor_variant) = tumor_variant {
            output.var1_total_coverage = tumor_variant.total_pos_coverage;
            output.var1_variant_coverage = tumor_variant.position_coverage;
            output.var1_ref_forward_coverage = tumor_variant.ref_forward_count;
            output.var1_ref_reverse_coverage = tumor_variant.ref_reverse_count;
            output.var1_variant_forward_count = tumor_variant.vars_count_on_forward;
            output.var1_variant_reverse_count = tumor_variant.vars_count_on_reverse;
            output.var1_genotype = format_somatic_genotype(tumor_variant);
            output.var1_frequency = tumor_variant.frequency;
            output.var1_strand_bias_flag = format_somatic_strand_bias(tumor_variant);
            output.var1_mean_position = tumor_variant.mean_position;
            output.var1_is_at_least_at_2_positions = if tumor_variant.is_at_least_at_2_positions {
                1
            } else {
                0
            };
            output.var1_mean_quality = tumor_variant.mean_quality;
            output.var1_has_at_least_2_diff_qualities =
                if tumor_variant.has_at_least_2_diff_qualities {
                    1
                } else {
                    0
                };
            output.var1_mean_mapping_quality = tumor_variant.mean_mapping_quality;
            output.var1_high_quality_to_low_quality_ratio = qratio_from_counts(
                tumor_variant.high_qual_read_cnt,
                tumor_variant.low_qual_read_cnt,
            );
            output.var1_high_quality_reads_frequency = tumor_variant.high_quality_reads_frequency;
            output.var1_extra_frequency = tumor_variant.extra_frequency;
            output.var1_nm = tumor_variant.nm;
            output.var1_duprate = tumor_variant.duprate;
        }

        if let Some(normal_variant) = normal_variant {
            output.var2_total_coverage = normal_variant.total_pos_coverage;
            output.var2_variant_coverage = normal_variant.position_coverage;
            output.var2_ref_forward_coverage = normal_variant.ref_forward_count;
            output.var2_ref_reverse_coverage = normal_variant.ref_reverse_count;
            output.var2_variant_forward_count = normal_variant.vars_count_on_forward;
            output.var2_variant_reverse_count = normal_variant.vars_count_on_reverse;
            output.var2_genotype = format_somatic_genotype(normal_variant);
            output.var2_frequency = normal_variant.frequency;
            output.var2_strand_bias_flag = format_somatic_strand_bias(normal_variant);
            output.var2_mean_position = normal_variant.mean_position;
            output.var2_is_at_least_at_2_positions = if normal_variant.is_at_least_at_2_positions {
                1
            } else {
                0
            };
            output.var2_mean_quality = normal_variant.mean_quality;
            output.var2_has_at_least_2_diff_qualities =
                if normal_variant.has_at_least_2_diff_qualities {
                    1
                } else {
                    0
                };
            output.var2_mean_mapping_quality = normal_variant.mean_mapping_quality;
            output.var2_high_quality_to_low_quality_ratio = qratio_from_counts(
                normal_variant.high_qual_read_cnt,
                normal_variant.low_qual_read_cnt,
            );
            output.var2_high_quality_reads_frequency = normal_variant.high_quality_reads_frequency;
            output.var2_extra_frequency = normal_variant.extra_frequency;
            output.var2_nm = normal_variant.nm;
            output.var2_duprate = normal_variant.duprate;
        }

        output
    }

    pub fn to_string_55_columns(&self) -> String {
        let parts: Vec<String> = vec![
            self.sample.clone(),
            self.gene.clone(),
            self.chr.clone(),
            self.start_position.to_string(),
            self.end_position.to_string(),
            self.ref_allele.clone(),
            self.var_allele.clone(),
            self.var1_total_coverage.to_string(),
            self.var1_variant_coverage.to_string(),
            self.var1_ref_forward_coverage.to_string(),
            self.var1_ref_reverse_coverage.to_string(),
            self.var1_variant_forward_count.to_string(),
            self.var1_variant_reverse_count.to_string(),
            self.var1_genotype.clone(),
            format_f64(self.var1_frequency, 4),
            self.var1_strand_bias_flag.clone(),
            format_f64(self.var1_mean_position, 1),
            self.var1_is_at_least_at_2_positions.to_string(),
            format_f64(self.var1_mean_quality, 1),
            self.var1_has_at_least_2_diff_qualities.to_string(),
            format_f64(self.var1_mean_mapping_quality, 1),
            format_f64(self.var1_high_quality_to_low_quality_ratio, 3),
            format_f64(self.var1_high_quality_reads_frequency, 4),
            format_f64(self.var1_extra_frequency, 4),
            format_f64(
                if self.var1_nm > 0.0 {
                    self.var1_nm
                } else {
                    0.0
                },
                1,
            ),
            self.var2_total_coverage.to_string(),
            self.var2_variant_coverage.to_string(),
            self.var2_ref_forward_coverage.to_string(),
            self.var2_ref_reverse_coverage.to_string(),
            self.var2_variant_forward_count.to_string(),
            self.var2_variant_reverse_count.to_string(),
            self.var2_genotype.clone(),
            format_f64(self.var2_frequency, 4),
            self.var2_strand_bias_flag.clone(),
            format_f64(self.var2_mean_position, 1),
            self.var2_is_at_least_at_2_positions.to_string(),
            format_f64(self.var2_mean_quality, 1),
            self.var2_has_at_least_2_diff_qualities.to_string(),
            format_f64(self.var2_mean_mapping_quality, 1),
            format_f64(self.var2_high_quality_to_low_quality_ratio, 3),
            format_f64(self.var2_high_quality_reads_frequency, 4),
            format_f64(self.var2_extra_frequency, 4),
            format_f64(
                if self.var2_nm > 0.0 {
                    self.var2_nm
                } else {
                    0.0
                },
                1,
            ),
            self.shift3.to_string(),
            format_f64(self.msi, 3),
            format_f64(self.msint, 0),
            self.left_sequence.clone(),
            self.right_sequence.clone(),
            self.region.clone(),
            self.var_label.clone(),
            self.var_type.clone(),
            format_f64(self.var1_duprate, 1),
            self.var1_sv.clone(),
            format_f64(self.var2_duprate, 1),
            self.var2_sv.clone(),
        ];

        parts.join("\t")
    }

    pub fn to_string_61_columns(&self) -> String {
        let var1_nm = if self.var1_nm > 0.0 {
            self.var1_nm
        } else {
            0.0
        };
        let var2_nm = if self.var2_nm > 0.0 {
            self.var2_nm
        } else {
            0.0
        };
        let msi_f = if self.msi == 0.0 {
            "0".to_string()
        } else {
            format!("{:.3}", self.msi)
        };

        let fisher1 = FisherExact::new(
            self.var1_ref_forward_coverage,
            self.var1_ref_reverse_coverage,
            self.var1_variant_forward_count,
            self.var1_variant_reverse_count,
        );
        let pvalue1 = fisher1.p_value();
        let oddratio1 = fisher1.odd_ratio();

        let fisher2 = FisherExact::new(
            self.var2_ref_forward_coverage,
            self.var2_ref_reverse_coverage,
            self.var2_variant_forward_count,
            self.var2_variant_reverse_count,
        );
        let pvalue2 = fisher2.p_value();
        let oddratio2 = fisher2.odd_ratio();

        let tref = self
            .var1_total_coverage
            .saturating_sub(self.var1_variant_coverage);
        let rref = self
            .var2_total_coverage
            .saturating_sub(self.var2_variant_coverage);

        let fisher = FisherExact::new(
            self.var1_variant_coverage,
            tref,
            self.var2_variant_coverage,
            rref,
        );
        let pvalue = fisher.p_value_less().min(fisher.p_value_greater());
        let oddratio = fisher.odd_ratio();

        let parts: Vec<String> = vec![
            self.sample.clone(),
            self.gene.clone(),
            self.chr.clone(),
            self.start_position.to_string(),
            self.end_position.to_string(),
            self.ref_allele.clone(),
            self.var_allele.clone(),
            self.var1_total_coverage.to_string(),
            self.var1_variant_coverage.to_string(),
            self.var1_ref_forward_coverage.to_string(),
            self.var1_ref_reverse_coverage.to_string(),
            self.var1_variant_forward_count.to_string(),
            self.var1_variant_reverse_count.to_string(),
            self.var1_genotype.clone(),
            format_rounded_value_to_print("0.0000", self.var1_frequency),
            self.var1_strand_bias_flag.clone(),
            format_rounded_value_to_print("0.0", self.var1_mean_position),
            self.var1_is_at_least_at_2_positions.to_string(),
            format_rounded_value_to_print("0.0", self.var1_mean_quality),
            self.var1_has_at_least_2_diff_qualities.to_string(),
            format_rounded_value_to_print("0.0", self.var1_mean_mapping_quality),
            format_rounded_value_to_print("0.000", self.var1_high_quality_to_low_quality_ratio),
            format_rounded_value_to_print("0.0000", self.var1_high_quality_reads_frequency),
            format_rounded_value_to_print("0.0000", self.var1_extra_frequency),
            format_rounded_value_to_print("0.0", var1_nm),
            format_rounded_value_to_print("0.00000", pvalue1),
            oddratio1,
            self.var2_total_coverage.to_string(),
            self.var2_variant_coverage.to_string(),
            self.var2_ref_forward_coverage.to_string(),
            self.var2_ref_reverse_coverage.to_string(),
            self.var2_variant_forward_count.to_string(),
            self.var2_variant_reverse_count.to_string(),
            self.var2_genotype.clone(),
            format_rounded_value_to_print("0.0000", self.var2_frequency),
            self.var2_strand_bias_flag.clone(),
            format_rounded_value_to_print("0.0", self.var2_mean_position),
            self.var2_is_at_least_at_2_positions.to_string(),
            format_rounded_value_to_print("0.0", self.var2_mean_quality),
            self.var2_has_at_least_2_diff_qualities.to_string(),
            format_rounded_value_to_print("0.0", self.var2_mean_mapping_quality),
            format_rounded_value_to_print("0.000", self.var2_high_quality_to_low_quality_ratio),
            format_rounded_value_to_print("0.0000", self.var2_high_quality_reads_frequency),
            format_rounded_value_to_print("0.0000", self.var2_extra_frequency),
            format_rounded_value_to_print("0.0", var2_nm),
            format_rounded_value_to_print("0.00000", pvalue2),
            oddratio2,
            self.shift3.to_string(),
            msi_f,
            format_f64(self.msint, 0),
            self.left_sequence.clone(),
            self.right_sequence.clone(),
            self.region.clone(),
            self.var_label.clone(),
            self.var_type.clone(),
            format_rounded_value_to_print("0.00", self.var1_duprate),
            self.var1_sv.clone(),
            format_rounded_value_to_print("0.00", self.var2_duprate),
            self.var2_sv.clone(),
            format_rounded_value_to_print("0.00000", pvalue),
            oddratio,
        ];

        parts.join("\t")
    }
}

impl std::fmt::Display for SomaticOutputVariant {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let output_variant = if instance().conf.fisher {
            self.to_string_61_columns()
        } else {
            self.to_string_55_columns()
        };
        if instance().conf.debug {
            write!(f, "{}\t{}", output_variant, self.debug)
        } else {
            write!(f, "{}", output_variant)
        }
    }
}

fn build_amplicon_debug(
    variant: &Variant,
    good_variants: &[(Variant, String)],
    bad_variants: &[(Option<Variant>, String)],
    debug_prefix: Option<&str>,
) -> String {
    let mut debug = debug_prefix
        .map(|value| value.to_string())
        .unwrap_or_else(|| format_variant_debug_content(variant));

    for (index, (good_variant, region)) in good_variants.iter().enumerate() {
        debug.push('\t');
        debug.push_str(&format!(
            "Good{} {} {}",
            index,
            format_debug_amp_variant(good_variant),
            region
        ));
    }

    for (index, (bad_variant, region)) in bad_variants.iter().enumerate() {
        debug.push('\t');
        if let Some(bad_variant) = bad_variant {
            debug.push_str(&format!(
                "Bad{} {} {}",
                index,
                format_debug_amp_variant(bad_variant),
                region
            ));
        } else {
            debug.push_str(&format!("Bad{} {}", index, region));
        }
    }

    debug
}

pub(crate) fn format_variant_debug_content(variant: &Variant) -> String {
    let mut key = variant.description_string.clone();
    if key.starts_with('+') {
        key = format!("I{}", key);
    }

    let total_variant_reads = variant.vars_count_on_forward + variant.vars_count_on_reverse;
    let pstd = if variant.is_at_least_at_2_positions {
        1
    } else {
        0
    };
    let qstd = if variant.has_at_least_2_diff_qualities {
        1
    } else {
        0
    };

    format!(
        "{}:{}:F-{}:R-{}:{}:{}:{}:{}:{}:{}:{}:{}:{}",
        key,
        total_variant_reads,
        variant.vars_count_on_forward,
        variant.vars_count_on_reverse,
        format_fixed_f64(variant.frequency, 4),
        variant.strand_bias_flag.var_bias.as_int(),
        format_fixed_f64(variant.mean_position, 1),
        pstd,
        format_fixed_f64(variant.mean_quality, 1),
        qstd,
        format_fixed_f64(variant.high_quality_reads_frequency, 4),
        format_fixed_f64(variant.mean_mapping_quality, 1),
        format_fixed_f64(debug_qratio(variant), 3),
    )
}

fn format_debug_amp_variant(variant: &Variant) -> String {
    let genotype = if variant.genotype.is_empty() {
        "0".to_string()
    } else {
        variant.genotype.clone()
    };
    let frequency = if variant.frequency == 0.0 {
        "0".to_string()
    } else {
        format_fixed_f64(variant.frequency, 4)
    };
    let pstd = if variant.is_at_least_at_2_positions {
        1
    } else {
        0
    };
    let qstd = if variant.has_at_least_2_diff_qualities {
        1
    } else {
        0
    };
    let hifreq = if variant.high_quality_reads_frequency == 0.0 {
        "0".to_string()
    } else {
        format_fixed_f64(variant.high_quality_reads_frequency, 4)
    };
    let extrafreq = if variant.extra_frequency == 0.0 {
        "0".to_string()
    } else {
        format_fixed_f64(variant.extra_frequency, 4)
    };

    format!(
        "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}",
        variant.total_pos_coverage,
        variant.position_coverage,
        variant.ref_forward_count,
        variant.ref_reverse_count,
        variant.vars_count_on_forward,
        variant.vars_count_on_reverse,
        genotype,
        frequency,
        variant.strand_bias_flag.to_string(),
        format_fixed_f64(variant.mean_position, 1),
        pstd,
        format_fixed_f64(variant.mean_quality, 1),
        qstd,
        format_fixed_f64(variant.mean_mapping_quality, 1),
        format_fixed_f64(debug_qratio(variant), 3),
        hifreq,
        extrafreq,
    )
}

fn debug_qratio(variant: &Variant) -> f64 {
    qratio_from_counts(variant.high_qual_read_cnt, variant.low_qual_read_cnt)
}

fn qratio_from_counts(high_qual_read_cnt: usize, low_qual_read_cnt: usize) -> f64 {
    if low_qual_read_cnt > 0 {
        high_qual_read_cnt as f64 / low_qual_read_cnt as f64
    } else if high_qual_read_cnt > 0 {
        high_qual_read_cnt as f64 * 2.0
    } else {
        0.0
    }
}

fn is_somatic_placeholder_variant(variant: &Variant) -> bool {
    variant.position_coverage == 0
        && variant.vars_count_on_forward == 0
        && variant.vars_count_on_reverse == 0
        && variant.mean_position == 0.0
        && variant.mean_quality == 0.0
        && variant.mean_mapping_quality == 0.0
        && variant.high_quality_reads_frequency == 0.0
        && variant.extra_frequency == 0.0
        && variant.nm == 0.0
        && !variant.is_at_least_at_2_positions
        && !variant.has_at_least_2_diff_qualities
}

fn format_somatic_genotype(variant: &Variant) -> String {
    if variant.genotype.is_empty()
        || (variant.genotype == "0/0" && is_somatic_placeholder_variant(variant))
    {
        "0".to_string()
    } else {
        variant.genotype.clone()
    }
}

fn format_somatic_strand_bias(variant: &Variant) -> String {
    if is_somatic_placeholder_variant(variant)
        && variant.strand_bias_flag == StrandBiasFlag::default()
    {
        "0".to_string()
    } else {
        variant.strand_bias_flag.to_string()
    }
}

fn format_rounded_value_to_print(pattern: &str, value: f64) -> String {
    if value == value.round() {
        format!("{:.0}", value)
    } else {
        let decimals = pattern
            .split('.')
            .nth(1)
            .map(|part| part.len())
            .unwrap_or(0);
        let rounded = round_half_even(pattern, value);
        let mut text = format!("{rounded:.decimals$}");
        while text.contains('.') && text.ends_with('0') {
            text.pop();
        }
        if text.ends_with('.') {
            text.pop();
        }
        text
    }
}

#[derive(Debug, Clone)]
struct FisherExact {
    m: i32,
    n: i32,
    k: i32,
    x: i32,
    lo: i32,
    hi: i32,
    support: Vec<i32>,
    logdc: Vec<f64>,
    p_value_less: f64,
    p_value_greater: f64,
    p_value_two_sided: f64,
}

impl FisherExact {
    fn new(ref_fwd: usize, ref_rev: usize, alt_fwd: usize, alt_rev: usize) -> Self {
        let m = (ref_fwd + ref_rev) as i32;
        let n = (alt_fwd + alt_rev) as i32;
        let k = (ref_fwd + alt_fwd) as i32;
        let x = ref_fwd as i32;
        let lo = (k - n).max(0);
        let hi = k.min(m);
        let support = (lo..=hi).collect::<Vec<_>>();

        let mut fisher = Self {
            m,
            n,
            k,
            x,
            lo,
            hi,
            support,
            logdc: Vec::new(),
            p_value_less: 0.0,
            p_value_greater: 0.0,
            p_value_two_sided: 0.0,
        };

        fisher.logdc = fisher.logdc_dhyper();
        fisher.calculate_pvalue();
        fisher
    }

    fn odd_ratio(&self) -> String {
        let odd_ratio = self.mle(self.x as f64);
        self.format_odd_ratio_value(odd_ratio)
    }

    /// Ported from: `com.astrazeneca.vardict.data.fishertest.FisherExact.getOddRatio()`
    /// Java source: `FisherExact.java:L136-L144`
    fn format_odd_ratio_value(&self, odd_ratio: f64) -> String {
        if odd_ratio.is_infinite() {
            "Inf".to_string()
        } else if odd_ratio == odd_ratio.round() {
            format!("{}", odd_ratio as i64)
        } else {
            Self::format_java_double(self.round_as_r(odd_ratio))
        }
    }

    fn format_java_double(value: f64) -> String {
        if value.is_finite() && value == value.round() {
            format!("{value:.1}")
        } else {
            value.to_string()
        }
    }

    fn p_value(&self) -> f64 {
        self.round_as_r(self.p_value_two_sided)
    }

    fn p_value_greater(&self) -> f64 {
        self.round_as_r(self.p_value_greater)
    }

    fn p_value_less(&self) -> f64 {
        self.round_as_r(self.p_value_less)
    }

    fn calculate_pvalue(&mut self) {
        self.p_value_less = self.pnhyper(self.x, false);
        self.p_value_greater = self.pnhyper(self.x, true);

        let rel_err = 1.0 + 1e-7;
        let d = self.dnhyper(1.0);
        let x_index = (self.x - self.lo) as usize;
        let threshold = d.get(x_index).copied().unwrap_or(0.0) * rel_err;
        let mut sum = 0.0;
        for value in d {
            if value <= threshold {
                sum += value;
            }
        }
        self.p_value_two_sided = sum;
    }

    fn pnhyper(&self, q: i32, upper_tail: bool) -> f64 {
        if self.m + self.n == 0 {
            return 1.0;
        }

        let distribution =
            match Hypergeometric::new((self.m + self.n) as u64, self.m as u64, self.k as u64) {
                Ok(distribution) => distribution,
                Err(_) => return 0.0,
            };

        if upper_tail {
            if q <= 0 {
                1.0
            } else {
                1.0 - distribution.cdf((q - 1) as u64)
            }
        } else if q < 0 {
            0.0
        } else {
            distribution.cdf(q as u64)
        }
    }

    fn logdc_dhyper(&self) -> Vec<f64> {
        let mut values = Vec::with_capacity(self.support.len());

        for element in &self.support {
            if self.m + self.n == 0 {
                values.push(0.0);
                continue;
            }

            let distribution =
                Hypergeometric::new((self.m + self.n) as u64, self.m as u64, self.k as u64);

            let value = match distribution {
                Ok(distribution) => {
                    let value = distribution.ln_pmf(*element as u64);
                    // Java: if (Double.isNaN(value)) value = 0; — only NaN is replaced,
                    // -Infinity is kept (exp(-inf)=0 is correct in downstream dnhyper)
                    if value.is_nan() { 0.0 } else { value }
                }
                Err(_) => 0.0,
            };
            values.push(round_half_even("0.0000000", value));
        }

        values
    }

    fn dnhyper(&self, ncp: f64) -> Vec<f64> {
        let mut result = Vec::with_capacity(self.support.len());
        for (index, support_value) in self.support.iter().enumerate() {
            result.push(self.logdc[index] + ncp.ln() * *support_value as f64);
        }

        let max_value = result.iter().copied().fold(f64::NEG_INFINITY, f64::max);

        let exponent = result
            .iter()
            .map(|value| (*value - max_value).exp())
            .collect::<Vec<_>>();
        let sum: f64 = exponent.iter().sum();
        exponent
            .iter()
            .map(|value| *value / sum)
            .collect::<Vec<_>>()
    }

    fn mnhyper(&self, ncp: f64) -> f64 {
        if ncp == 0.0 {
            return self.lo as f64;
        }
        if ncp.is_infinite() {
            return self.hi as f64;
        }

        let dnhyper = self.dnhyper(ncp);
        self.support
            .iter()
            .zip(dnhyper.iter())
            .map(|(support_value, probability)| *support_value as f64 * *probability)
            .sum()
    }

    fn mle(&self, x: f64) -> f64 {
        let eps = f64::EPSILON;
        if (x - self.lo as f64).abs() < f64::EPSILON {
            return 0.0;
        }
        if (x - self.hi as f64).abs() < f64::EPSILON {
            return f64::INFINITY;
        }

        let mu = self.mnhyper(1.0);
        if mu > x {
            zeroin_c(0.0, 1.0, |t| self.mnhyper(t) - x, eps.powf(0.25))
        } else if mu < x {
            1.0 / zeroin_c(eps, 1.0, |t| self.mnhyper(1.0 / t) - x, eps.powf(0.25))
        } else {
            1.0
        }
    }

    fn round_as_r(&self, value: f64) -> f64 {
        let mut rounded = round_half_even("0", value * 1e5);
        rounded /= 1e5;
        if rounded == 0.0 {
            0.0
        } else if rounded == 1.0 {
            1.0
        } else {
            rounded
        }
    }
}

fn zeroin_c<F>(ax: f64, bx: f64, f: F, tol: f64) -> f64
where
    F: Fn(f64) -> f64,
{
    let mut a = ax;
    let mut b = bx;
    let mut fa = f(a);
    let mut fb = f(b);
    let mut c = a;
    let mut fc = fa;
    let epsilon = f64::EPSILON;

    loop {
        let prev_step = b - a;

        if fc.abs() < fb.abs() {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }

        let tol_act = 2.0 * epsilon * b.abs() + tol / 2.0;
        let mut new_step = (c - b) / 2.0;

        if new_step.abs() <= tol_act || fb == 0.0 {
            return b;
        }

        if prev_step.abs() >= tol_act && fa.abs() > fb.abs() {
            let cb = c - b;
            let (mut p, mut q) = if a == c {
                let t1 = fb / fa;
                (cb * t1, 1.0 - t1)
            } else {
                let q = fa / fc;
                let t1 = fb / fc;
                let t2 = fb / fa;
                (
                    t2 * (cb * q * (q - t1) - (b - a) * (t1 - 1.0)),
                    (q - 1.0) * (t1 - 1.0) * (t2 - 1.0),
                )
            };

            if p > 0.0 {
                q = -q;
            } else {
                p = -p;
            }

            if p < (0.75 * cb * q - (tol_act * q).abs() / 2.0) && p < (prev_step * q / 2.0).abs() {
                new_step = p / q;
            }
        }

        if new_step.abs() < tol_act {
            new_step = if new_step > 0.0 { tol_act } else { -tol_act };
        }

        a = b;
        fa = fb;
        b += new_step;
        fb = f(b);

        if (fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0) {
            c = a;
            fc = fa;
        }
    }
}

fn format_fixed_f64(value: f64, decimals: usize) -> String {
    format!("{:.1$}", value, decimals)
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

fn format_output_var_type(variant: &Variant) -> String {
    if variant.refallele == variant.varallele {
        return String::new();
    }

    if variant.varallele.starts_with('<')
        && variant.varallele.ends_with('>')
        && variant.varallele.len() >= 5
    {
        return var_type_string(&variant.refallele, &variant.varallele);
    }

    format_var_type(&variant.vartype)
}

/// Get column headers for 36-column format
pub fn get_column_headers() -> Vec<&'static str> {
    vec![
        "Sample",
        "Gene",
        "Chr",
        "Start",
        "End",
        "Ref",
        "Alt",
        "Depth",
        "AltDepth",
        "RefFwdReads",
        "RefRevReads",
        "AltFwdReads",
        "AltRevReads",
        "Genotype",
        "AF",
        "Bias",
        "PMean",
        "PStd",
        "QMean",
        "QStd",
        "MQ",
        "Sig_Noise",
        "HiAF",
        "ExtraAF",
        "shift3",
        "MSI",
        "MSI_NT",
        "NM",
        "HiCnt",
        "HiCov",
        "5pFlankSeq",
        "3pFlankSeq",
        "Seg",
        "VarType",
        "Duprate",
        "SV_info",
    ]
}

/// Get simple-mode column headers with optional CRISPR suffix (Java parity)
pub fn get_simple_column_headers(crispr_enabled: bool) -> Vec<&'static str> {
    let mut headers = get_column_headers();
    if crispr_enabled {
        headers.push("CRISPR");
    }
    headers
}

/// Get column headers as tab-delimited string
pub fn get_header_line() -> String {
    get_column_headers().join("\t")
}

/// Get simple-mode header line with optional CRISPR suffix (Java parity)
pub fn get_simple_header_line(crispr_enabled: bool) -> String {
    get_simple_column_headers(crispr_enabled).join("\t")
}

pub fn get_amplicon_column_headers() -> Vec<&'static str> {
    vec![
        "Sample",
        "Gene",
        "Chr",
        "Start",
        "End",
        "Ref",
        "Alt",
        "Depth",
        "AltDepth",
        "RefFwdReads",
        "RefRevReads",
        "AltFwdReads",
        "AltRevReads",
        "Genotype",
        "AF",
        "Bias",
        "PMean",
        "PStd",
        "QMean",
        "QStd",
        "MQ",
        "Sig_Noise",
        "HiAF",
        "ExtraAF",
        "shift3",
        "MSI",
        "MSI_NT",
        "NM",
        "HiCnt",
        "HiCov",
        "5pFlankSeq",
        "3pFlankSeq",
        "Seg",
        "VarType",
        "GoodVarCount",
        "TotalVarCount",
        "Nocov",
        "Ampflag",
    ]
}

pub fn get_amplicon_header_line() -> String {
    get_amplicon_column_headers().join("\t")
}

pub fn get_somatic_column_headers() -> Vec<&'static str> {
    vec![
        "Sample",
        "Gene",
        "Chr",
        "Start",
        "End",
        "Ref",
        "Alt",
        "Depth",
        "AltDepth",
        "RefFwdReads",
        "RefRevReads",
        "AltFwdReads",
        "AltRevReads",
        "Genotype",
        "AF",
        "Bias",
        "PMean",
        "PStd",
        "QMean",
        "QStd",
        "MQ",
        "Sig_Noise",
        "HiAF",
        "ExtraAF",
        "NM",
        "Depth",
        "AltDepth",
        "RefFwdReads",
        "RefRevReads",
        "AltFwdReads",
        "AltRevReads",
        "Genotype",
        "AF",
        "Bias",
        "PMean",
        "PStd",
        "QMean",
        "QStd",
        "MQ",
        "Sig_Noise",
        "HiAF",
        "ExtraAF",
        "NM",
        "shift3",
        "MSI",
        "MSI_NT",
        "5pFlankSeq",
        "3pFlankSeq",
        "Seg",
        "VarLabel",
        "VarType",
        "Duprate1",
        "SV_info1",
        "Duprate2",
        "SV_info2",
    ]
}

pub fn get_somatic_header_line() -> String {
    get_somatic_column_headers().join("\t")
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mods::to_vars_builder::VarType;
    use crate::scopedata::global_read_only_scope::GlobalReadOnlyScope;

    #[test]
    fn test_region_to_string() {
        let region = Region::new("chr1", 1000, 2000, "GENE1");
        assert_eq!(region.to_region_string(), "chr1:1000-2000");
    }

    #[test]
    fn test_region_to_string_uses_display_start() {
        let mut region = Region::new("20", 1, 1_000_150, "20");
        region.display_start = -149;

        assert_eq!(region.to_region_string(), "20:-149-1000150");
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
        assert_eq!(
            format_strand_bias(StrandBiasFlag::new(
                StrandBiasValue::CantAssess,
                StrandBiasValue::CantAssess
            )),
            "0;0"
        );
        assert_eq!(
            format_strand_bias(StrandBiasFlag::new(
                StrandBiasValue::NoBias,
                StrandBiasValue::NoBias
            )),
            "2;2"
        );
        assert_eq!(
            format_strand_bias(StrandBiasFlag::new(
                StrandBiasValue::NoBias,
                StrandBiasValue::HasBias
            )),
            "2;1"
        );
    }

    #[test]
    fn test_format_var_type() {
        assert_eq!(format_var_type(&VarType::SNV('A')), "SNV");
        assert_eq!(
            format_var_type(&VarType::Insertion("ATG".to_string())),
            "Insertion"
        );
        assert_eq!(format_var_type(&VarType::Deletion(3)), "Deletion");
        assert_eq!(
            format_var_type(&VarType::Complex {
                insertion: "A".to_string(),
                deletion: 2
            }),
            "Complex"
        );
    }

    #[test]
    fn test_simple_output_variant_from_variant_uses_symbolic_alleles_for_var_type() {
        let region = Region::new("1", 1_500_001, 2_000_000, "testbed");
        let variant = Variant {
            description_string: "<dup52673>".to_string(),
            refallele: "C".to_string(),
            varallele: "<DUP>".to_string(),
            vartype: VarType::Complex {
                insertion: String::new(),
                deletion: 0,
            },
            start_position: 1_628_406,
            end_position: 1_681_078,
            vars_count_on_forward: 0,
            vars_count_on_reverse: 3,
            position_coverage: 3,
            total_pos_coverage: 8,
            frequency: 0.375,
            threshold_frequency: 0.375,
            high_quality_reads_frequency: 0.4286,
            extra_frequency: 0.375,
            mean_position: 34.3,
            mean_quality: 37.0,
            mean_mapping_quality: 15.0,
            strand_bias_flag: StrandBiasFlag::default(),
            is_at_least_at_2_positions: true,
            has_at_least_2_diff_qualities: false,
            leftseq: "CTCTGAGTGTGTGGTGCCTG".to_string(),
            rightseq: "TGTGTGTGTGTGAATCTACG".to_string(),
            msi: 0.0,
            msint: 0.0,
            shift3: 0,
            nm: 2.0,
            high_qual_read_cnt: 3,
            low_qual_read_cnt: 0,
            hicov: 7,
            ref_forward_count: 5,
            ref_reverse_count: 3,
            genotype: "C/+52673".to_string(),
            duprate: 0.0,
            crispr: 0,
        };

        let output = SimpleOutputVariant::from_variant(&variant, &region, "sample", "3-0-0");

        assert_eq!(output.var_type, "DUP");
    }

    #[test]
    fn test_simple_output_variant_from_variant_preserves_complex_vartype() {
        let region = Region::new("1", 11_000_001, 11_500_000, "testbed");
        let variant = Variant {
            description_string: "-2TC".to_string(),
            refallele: "CTT".to_string(),
            varallele: "C".to_string(),
            vartype: VarType::Complex {
                insertion: String::new(),
                deletion: 2,
            },
            start_position: 11_464_646,
            end_position: 11_464_648,
            vars_count_on_forward: 11,
            vars_count_on_reverse: 6,
            position_coverage: 17,
            total_pos_coverage: 20,
            frequency: 0.85,
            threshold_frequency: 0.85,
            high_quality_reads_frequency: 0.8095,
            extra_frequency: 0.0,
            mean_position: 19.9,
            mean_quality: 31.6,
            mean_mapping_quality: 31.3,
            strand_bias_flag: StrandBiasFlag::default(),
            is_at_least_at_2_positions: true,
            has_at_least_2_diff_qualities: true,
            leftseq: "GAAAATATAGATCTGTAAAT".to_string(),
            rightseq: "TTTATAAAAATACATTTAAA".to_string(),
            msi: 4.0,
            msint: 1.0,
            shift3: 0,
            nm: 1.2,
            high_qual_read_cnt: 17,
            low_qual_read_cnt: 4,
            hicov: 21,
            ref_forward_count: 3,
            ref_reverse_count: 0,
            genotype: "TTT/-2TC".to_string(),
            duprate: 0.0769,
            crispr: 0,
        };

        let output = SimpleOutputVariant::from_variant(&variant, &region, "sample", "0");

        assert_eq!(output.var_type, "Complex");
    }

    #[test]
    fn test_empty_variant() {
        let region = Region::new("chr1", 1000, 2000, "GENE1");
        let output = SimpleOutputVariant::empty(1500, &region, "sample1");

        assert_eq!(output.sample, "sample1");
        assert_eq!(output.start_position, 1500);
        assert_eq!(output.genotype, "0");
        assert_eq!(output.bias, "0");
        assert_eq!(output.left_sequence, "0");
        assert_eq!(output.right_sequence, "0");
        assert_eq!(output.sv, "0");
        assert_eq!(output.debug, "");
    }

    #[test]
    fn test_empty_position_zero_row_matches_java_simple_placeholder_parity() {
        let region = Region::new("20", 1, 1_000_000, "20");
        let output = SimpleOutputVariant::empty_with_sv(0, &region, "NA12878", "");

        assert_eq!(
            output.to_string_36_columns(),
            "NA12878\t20\t20\t0\t0\t\t\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t20:1-1000000\t\t0\t0"
        );
    }

    #[test]
    fn test_empty_null_variant_row_matches_java_zero_coverage_parity() {
        let region = Region::new("20", 1, 1_000_000, "20");
        let output =
            SimpleOutputVariant::empty_null_variant_with_sv(122276, &region, "NA12878", "");

        // Java null-variant path: genotype="", bias="0;0", left_seq="", right_seq=""
        assert_eq!(output.genotype, "");
        assert_eq!(output.bias, "0;0");
        assert_eq!(output.left_sequence, "");
        assert_eq!(output.right_sequence, "");
        assert_eq!(
            output.to_string_36_columns(),
            "NA12878\t20\t20\t122276\t122276\t\t\t0\t0\t0\t0\t0\t0\t\t0\t0;0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t\t\t20:1-1000000\t\t0\t0"
        );
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
            position_coverage: 10,
            total_pos_coverage: 100,
            frequency: 0.10,
            threshold_frequency: 0.10,
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
            crispr: 0,
        };

        let region = Region::new("chr1", 900, 1100, "BRCA1");
        let output = SimpleOutputVariant::from_variant(&variant, &region, "sample1", "");

        assert_eq!(output.sample, "sample1");
        assert_eq!(output.gene, "BRCA1");
        assert_eq!(output.chr, "chr1");
        assert_eq!(output.start_position, 1000);
        assert_eq!(output.ref_allele, "A");
        assert_eq!(output.var_allele, "T");
        assert_eq!(output.variant_forward_count, 5);
        assert_eq!(output.variant_reverse_count, 5);
        assert_eq!(output.genotype, "0/1");
        assert!((output.frequency - 0.10).abs() < 0.001);
        assert_eq!(output.bias, "2;2"); // NoBias for both ref and var
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
        assert_eq!(fields[0], "sample1"); // Sample
        assert_eq!(fields[1], "GENE1"); // Gene
        assert_eq!(fields[2], "chr1"); // Chr
        assert_eq!(fields[3], "1500"); // Start
        assert_eq!(fields[5], "A"); // Ref
        assert_eq!(fields[6], "T"); // Alt
        assert_eq!(fields[7], "100"); // Total coverage
        assert_eq!(fields[13], "0/1"); // Genotype
        assert_eq!(fields[14], "0.1000"); // Frequency
        assert_eq!(fields[33], "SNV"); // VarType
    }

    #[test]
    fn test_get_column_headers() {
        let headers = get_column_headers();
        assert_eq!(headers.len(), 36);
        assert_eq!(headers[0], "Sample");
        assert_eq!(headers[9], "RefFwdReads");
        assert_eq!(headers[18], "QMean");
        assert_eq!(headers[21], "Sig_Noise");
        assert_eq!(headers[24], "shift3");
        assert_eq!(headers[30], "5pFlankSeq");
        assert_eq!(headers[13], "Genotype");
        assert_eq!(headers[33], "VarType");
        assert_eq!(headers[34], "Duprate");
        assert_eq!(headers[35], "SV_info");
    }

    #[test]
    fn test_get_header_line() {
        let header = get_header_line();
        let fields: Vec<&str> = header.split('\t').collect();
        assert_eq!(fields.len(), 36);
    }

    #[test]
    fn test_get_simple_header_line_crispr_suffix_java_parity() {
        let header = get_simple_header_line(true);
        let fields: Vec<&str> = header.split('\t').collect();
        assert_eq!(fields.len(), 37);
        assert_eq!(fields[36], "CRISPR");
    }

    #[test]
    fn test_simple_output_variant_to_string_38_columns_fisher_parity() {
        let region = Region::new("chr1", 1000, 2000, "GENE1");
        let mut output = SimpleOutputVariant::empty(1500, &region, "sample1");
        output.ref_allele = "A".to_string();
        output.var_allele = "T".to_string();
        output.total_coverage = 100;
        output.variant_coverage = 10;
        output.reference_forward_count = 20;
        output.reference_reverse_count = 15;
        output.variant_forward_count = 6;
        output.variant_reverse_count = 4;
        output.genotype = "0/1".to_string();
        output.frequency = 0.10;
        output.var_type = "SNV".to_string();

        let line = output.to_string_38_columns();
        let fields: Vec<&str> = line.split('\t').collect();

        assert_eq!(fields.len(), 38);
        assert_eq!(fields[0], "sample1");
        assert_eq!(fields[14], "0.1");
        assert!(fields[20].parse::<f64>().is_ok());
        assert!(fields[21] == "Inf" || fields[21].parse::<f64>().is_ok());
        assert_eq!(fields[37], "0");
    }

    #[test]
    fn test_simple_output_variant_crispr_column_behavior_parity() {
        let region = Region::new("chr1", 1000, 2000, "GENE1");
        let mut output = SimpleOutputVariant::empty(1500, &region, "sample1");
        output.crispr = 7;

        let line_no_fisher = output.to_string_with_flags(false, true);
        let fields_no_fisher: Vec<&str> = line_no_fisher.split('\t').collect();
        assert_eq!(fields_no_fisher.len(), 37);
        assert_eq!(fields_no_fisher[36], "7");

        let line_with_fisher = output.to_string_with_flags(true, true);
        let fields_with_fisher: Vec<&str> = line_with_fisher.split('\t').collect();
        assert_eq!(fields_with_fisher.len(), 39);
        assert_eq!(fields_with_fisher[38], "7");
    }

    #[test]
    fn test_simple_output_variant_debug_column_behavior_parity() {
        let region = Region::new("chr1", 1000, 2000, "GENE1");
        let mut output = SimpleOutputVariant::empty(1500, &region, "sample1");
        output.debug = "T:5:F-2:R-3:1.0000:2:21.8:1:42.8:1:1.0000:53.8:10.000".to_string();

        let line = output.to_string_with_flags_and_debug(false, false, true);
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), 37);
        assert_eq!(fields[36], output.debug);
    }

    #[test]
    fn test_amplicon_output_variant_to_string_38_columns() {
        let _ = INSTANCE.get_or_init(GlobalReadOnlyScope::default);
        let region = Region::new("chr1", 1000, 2000, "GENE1");
        let variant = Variant {
            description_string: "A>T".to_string(),
            refallele: "A".to_string(),
            varallele: "T".to_string(),
            vartype: VarType::SNV('T'),
            start_position: 1500,
            end_position: 1500,
            vars_count_on_forward: 5,
            vars_count_on_reverse: 5,
            position_coverage: 10,
            total_pos_coverage: 100,
            frequency: 0.10,
            threshold_frequency: 0.10,
            high_quality_reads_frequency: 0.08,
            extra_frequency: 0.0,
            mean_position: 25.0,
            mean_quality: 30.0,
            mean_mapping_quality: 60.0,
            strand_bias_flag: StrandBiasFlag::default(),
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
            ref_forward_count: 20,
            ref_reverse_count: 15,
            genotype: "0/1".to_string(),
            duprate: 0.0,
            crispr: 0,
        };
        let output = AmpliconOutputVariant::from_variant(
            Some(&variant),
            &region,
            &[],
            &[(None, String::new()), (None, String::new())],
            None,
            1500,
            1,
            0,
            false,
            "sample1",
        );

        let line = output.to_string_38_columns();
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), 38);
        assert_eq!(fields[0], "sample1");
        assert_eq!(fields[1], "GENE1");
        assert_eq!(fields[2], "chr1");
        assert_eq!(fields[34], "1");
        assert_eq!(fields[35], "3");
        assert_eq!(fields[37], "0");
    }

    #[test]
    fn test_get_amplicon_header_line() {
        let header = get_amplicon_header_line();
        let fields: Vec<&str> = header.split('\t').collect();
        assert_eq!(fields.len(), 38);
        assert_eq!(fields[34], "GoodVarCount");
        assert_eq!(fields[37], "Ampflag");
    }

    #[test]
    fn test_amplicon_output_variant_to_string_40_columns_fisher_parity() {
        let _ = INSTANCE.get_or_init(GlobalReadOnlyScope::default);
        let region = Region::new("chr1", 1000, 2000, "GENE1");
        let variant = Variant {
            description_string: "A>T".to_string(),
            refallele: "A".to_string(),
            varallele: "T".to_string(),
            vartype: VarType::SNV('T'),
            start_position: 1500,
            end_position: 1500,
            vars_count_on_forward: 5,
            vars_count_on_reverse: 5,
            position_coverage: 10,
            total_pos_coverage: 100,
            frequency: 0.10,
            threshold_frequency: 0.10,
            high_quality_reads_frequency: 0.08,
            extra_frequency: 0.0,
            mean_position: 25.0,
            mean_quality: 30.0,
            mean_mapping_quality: 60.0,
            strand_bias_flag: StrandBiasFlag::default(),
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
            ref_forward_count: 20,
            ref_reverse_count: 15,
            genotype: "0/1".to_string(),
            duprate: 0.0,
            crispr: 0,
        };
        let output = AmpliconOutputVariant::from_variant(
            Some(&variant),
            &region,
            &[],
            &[(None, String::new()), (None, String::new())],
            None,
            1500,
            1,
            0,
            false,
            "sample1",
        );

        let line = output.to_string_40_columns();
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), 40);
        assert_eq!(fields[0], "sample1");
        assert!(fields[20].parse::<f64>().is_ok());
        assert!(fields[21] == "Inf" || fields[21].parse::<f64>().is_ok());
        assert_eq!(fields[36], "1");
        assert_eq!(fields[39], "0");
    }

    #[test]
    fn test_amplicon_output_variant_tail_columns_contract() {
        let region = Region::new("chr1", 1000, 2000, "GENE1");
        let output = AmpliconOutputVariant::from_variant(
            None,
            &region,
            &[],
            &[],
            None,
            1500,
            2,
            3,
            true,
            "sample1",
        );

        let line = output.to_string_38_columns();
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), 38);
        assert_eq!(fields[34], "2");
        assert_eq!(fields[35], "2");
        assert_eq!(fields[36], "3");
        assert_eq!(fields[37], "1");
    }

    #[test]
    fn test_get_somatic_header_line() {
        let header = get_somatic_header_line();
        let fields: Vec<&str> = header.split('\t').collect();
        assert_eq!(fields.len(), 55);
        assert_eq!(fields[49], "VarLabel");
        assert_eq!(fields[54], "SV_info2");
    }

    #[test]
    fn test_somatic_output_variant_to_string_55_columns() {
        let region = Region::new("chr1", 900, 1100, "GENE1");

        let begin_variant = Variant {
            description_string: "A>T".to_string(),
            refallele: "A".to_string(),
            varallele: "T".to_string(),
            vartype: VarType::SNV('T'),
            start_position: 1000,
            end_position: 1000,
            vars_count_on_forward: 3,
            vars_count_on_reverse: 2,
            position_coverage: 5,
            total_pos_coverage: 20,
            frequency: 0.25,
            threshold_frequency: 0.25,
            high_quality_reads_frequency: 0.20,
            extra_frequency: 0.0,
            mean_position: 30.0,
            mean_quality: 35.0,
            mean_mapping_quality: 60.0,
            strand_bias_flag: StrandBiasFlag::default(),
            is_at_least_at_2_positions: true,
            has_at_least_2_diff_qualities: true,
            leftseq: "ACGT".to_string(),
            rightseq: "TGCA".to_string(),
            msi: 2.0,
            msint: 1.0,
            shift3: 1,
            nm: 1.0,
            high_qual_read_cnt: 5,
            low_qual_read_cnt: 1,
            hicov: 20,
            ref_forward_count: 10,
            ref_reverse_count: 5,
            genotype: "0/1".to_string(),
            duprate: 0.1,
            crispr: 0,
        };

        let end_variant = begin_variant.clone();

        let mut normal_variant = begin_variant.clone();
        normal_variant.frequency = 0.10;
        normal_variant.position_coverage = 2;
        normal_variant.vars_count_on_forward = 1;
        normal_variant.vars_count_on_reverse = 1;
        normal_variant.duprate = 0.0;

        let output = SomaticOutputVariant::from_variants(
            Some(&begin_variant),
            Some(&end_variant),
            Some(&begin_variant),
            Some(&normal_variant),
            &region,
            "sv1",
            "",
            "StrongSomatic",
            "Tumor|Normal",
        );

        let line = output.to_string_55_columns();
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), 55);
        assert_eq!(fields[0], "Tumor|Normal");
        assert_eq!(fields[6], "T");
        assert_eq!(fields[14], "0.2500");
        assert_eq!(fields[32], "0.1000");
        assert_eq!(fields[49], "StrongSomatic");
        assert_eq!(fields[50], "SNV");
        assert_eq!(fields[52], "sv1");
        assert_eq!(fields[54], "0");
    }

    #[test]
    fn test_somatic_output_variant_to_string_61_columns_fisher_parity() {
        let region = Region::new("chr1", 900, 1100, "GENE1");

        let begin_variant = Variant {
            description_string: "A>T".to_string(),
            refallele: "A".to_string(),
            varallele: "T".to_string(),
            vartype: VarType::SNV('T'),
            start_position: 1000,
            end_position: 1000,
            vars_count_on_forward: 3,
            vars_count_on_reverse: 2,
            position_coverage: 5,
            total_pos_coverage: 20,
            frequency: 0.25,
            threshold_frequency: 0.25,
            high_quality_reads_frequency: 0.20,
            extra_frequency: 0.0,
            mean_position: 30.0,
            mean_quality: 35.0,
            mean_mapping_quality: 60.0,
            strand_bias_flag: StrandBiasFlag::default(),
            is_at_least_at_2_positions: true,
            has_at_least_2_diff_qualities: true,
            leftseq: "ACGT".to_string(),
            rightseq: "TGCA".to_string(),
            msi: 2.0,
            msint: 1.0,
            shift3: 1,
            nm: 1.0,
            high_qual_read_cnt: 5,
            low_qual_read_cnt: 1,
            hicov: 20,
            ref_forward_count: 10,
            ref_reverse_count: 5,
            genotype: "0/1".to_string(),
            duprate: 0.1,
            crispr: 0,
        };

        let end_variant = begin_variant.clone();
        let mut normal_variant = begin_variant.clone();
        normal_variant.frequency = 0.10;
        normal_variant.position_coverage = 2;
        normal_variant.vars_count_on_forward = 1;
        normal_variant.vars_count_on_reverse = 1;
        normal_variant.duprate = 0.0;

        let output = SomaticOutputVariant::from_variants(
            Some(&begin_variant),
            Some(&end_variant),
            Some(&begin_variant),
            Some(&normal_variant),
            &region,
            "sv1",
            "",
            "StrongSomatic",
            "Tumor|Normal",
        );

        let line = output.to_string_61_columns();
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), 61);
        assert_eq!(fields[0], "Tumor|Normal");
        assert_eq!(fields[6], "T");
        assert_eq!(fields[14], "0.25");
        assert!(fields[25].parse::<f64>().is_ok());
        assert!(fields[26] == "Inf" || fields[26].parse::<f64>().is_ok());
        assert_eq!(fields[55], "0.1");
        assert_eq!(fields[56], "sv1");
        assert!(fields[59].parse::<f64>().is_ok());
        assert!(fields[60] == "Inf" || fields[60].parse::<f64>().is_ok());
    }

    #[test]
    fn test_fisher_exact_logdc_java_parity() {
        let fisher = FisherExact::new(11, 12, 1, 2);
        let expected = vec![-2.4696392, -1.0345547, -0.8675006, -1.9661129];
        assert_eq!(fisher.logdc.len(), expected.len());

        for (actual, expected_value) in fisher.logdc.iter().zip(expected.iter()) {
            assert!((actual - expected_value).abs() < 1e-7);
        }
    }

    #[test]
    fn test_fisher_exact_counts_java_parity() {
        let cases = vec![
            (
                121usize, 55usize, 18usize, 23usize, 0.00378, 0.99908, 0.00287, 2.79657,
            ),
            (121usize, 5usize, 18usize, 23usize, 0.0, 1.0, 0.0, 29.86184),
            (
                37usize, 76usize, 1usize, 1usize, 1.0, 0.55362, 0.89275, 0.49015,
            ),
            (0usize, 0usize, 1usize, 0usize, 1.0, 1.0, 1.0, 0.0),
            (1usize, 0usize, 1usize, 0usize, 1.0, 1.0, 1.0, 0.0),
            (0usize, 0usize, 0usize, 0usize, 1.0, 1.0, 1.0, 0.0),
            (0usize, 1usize, 1usize, 0usize, 1.0, 0.5, 1.0, 0.0),
            (1usize, 1usize, 1usize, 0usize, 1.0, 0.66667, 1.0, 0.0),
            (1usize, 1usize, 1usize, 1usize, 1.0, 0.83333, 0.83333, 1.0),
            (
                10usize, 10usize, 10usize, 1usize, 0.04722, 0.02599, 0.99802, 0.10703,
            ),
            (
                10usize, 10usize, 10usize, 0usize, 0.01099, 0.00615, 1.0, 0.0,
            ),
            (69usize, 1usize, 74usize, 95usize, 0.0, 1.0, 0.0, 87.68597),
            (
                41usize, 86usize, 1usize, 1usize, 0.54688, 0.54687, 0.89571, 0.47973,
            ),
            (
                130usize, 189usize, 1usize, 0usize, 0.40937, 0.40937, 1.0, 0.0,
            ),
            (
                83usize, 40usize, 1usize, 2usize, 0.25746, 0.96473, 0.25746, 4.09908,
            ),
            (
                74usize, 117usize, 1usize, 0usize, 0.39062, 0.39063, 1.0, 0.0,
            ),
            (60usize, 62usize, 2usize, 0usize, 0.49593, 0.24797, 1.0, 0.0),
            (
                43usize, 83usize, 1usize, 1usize, 1.0, 0.57111, 0.88361, 0.52091,
            ),
            (
                78usize, 40usize, 1usize, 1usize, 1.0, 0.88515, 0.56849, 1.93844,
            ),
        ];

        for (ref_fwd, ref_rev, alt_fwd, alt_rev, pvalue, p_less, p_greater, odd_ratio) in cases {
            let fisher = FisherExact::new(ref_fwd, ref_rev, alt_fwd, alt_rev);
            let actual_pvalue = fisher.p_value();
            let actual_p_less = fisher.p_value_less();
            let actual_p_greater = fisher.p_value_greater();

            assert!(
                (actual_pvalue - pvalue).abs() <= 2e-5,
                "p_value mismatch for ({},{},{},{}): actual={} expected={}",
                ref_fwd,
                ref_rev,
                alt_fwd,
                alt_rev,
                actual_pvalue,
                pvalue
            );
            assert!(
                (actual_p_less - p_less).abs() <= 2e-5,
                "p_value_less mismatch for ({},{},{},{}): actual={} expected={}",
                ref_fwd,
                ref_rev,
                alt_fwd,
                alt_rev,
                actual_p_less,
                p_less
            );
            assert!(
                (actual_p_greater - p_greater).abs() <= 2e-5,
                "p_value_greater mismatch for ({},{},{},{}): actual={} expected={}",
                ref_fwd,
                ref_rev,
                alt_fwd,
                alt_rev,
                actual_p_greater,
                p_greater
            );

            let parsed_odd_ratio = fisher
                .odd_ratio()
                .parse::<f64>()
                .expect("odd ratio should be numeric for Java parity cases");
            assert_eq!(round_half_even("0.00000", parsed_odd_ratio), odd_ratio);
        }
    }

    #[test]
    fn test_fisher_exact_odd_ratio_string_java_parity() {
        let fisher = FisherExact::new(0, 0, 0, 0);

        assert_eq!(fisher.format_odd_ratio_value(0.0), "0");
        assert_eq!(fisher.format_odd_ratio_value(1.0), "1");
        assert_eq!(fisher.format_odd_ratio_value(6.0), "6");
        assert_eq!(fisher.format_odd_ratio_value(5.999_999_999_999_999), "6.0");
        assert_eq!(fisher.format_odd_ratio_value(f64::INFINITY), "Inf");
        assert_eq!(fisher.format_odd_ratio_value(2.449_494_999_999_999_7), "2.44949");
        assert_eq!(fisher.format_odd_ratio_value(0.047_619_999_999_999_996), "0.04762");
        assert_eq!(fisher.format_odd_ratio_value(0.5), "0.5");
    }
}
