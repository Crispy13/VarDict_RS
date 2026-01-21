//! StructuralVariantsProcessor - Equivalent to Java StructuralVariantsProcessor
//!
//! This module processes realigned variation data to:
//! 1. Find structural variants (DEL, INV, DUP) when SV is enabled
//! 2. Adjust SNV counts from soft-clipped reads (always runs)
//!
//! Pipeline position:
//! ```text
//! CigarParser → VariantRealigner → StructuralVariantsProcessor → ToVarsBuilder
//! ```

use std::collections::HashMap;

use crate::scopedata::global_read_only_scope::instance;
use crate::variants::variants::{VarDesc, Variant, SoftClip};

/// Input data for StructuralVariantsProcessor (from VariantRealigner)
#[derive(Default)]
pub struct RealignedVariationData {
    /// Non-insertion variants by position
    pub non_insertion_variants: HashMap<i64, HashMap<VarDesc, Variant>>,
    /// Insertion variants by position
    pub insertion_variants: HashMap<i64, HashMap<VarDesc, Variant>>,
    /// 5' end soft clips by position
    pub soft_clips_5end: HashMap<i64, SoftClip>,
    /// 3' end soft clips by position  
    pub soft_clips_3end: HashMap<i64, SoftClip>,
    /// Reference coverage by position
    pub ref_coverage: HashMap<i64, usize>,
    /// Maximum read length seen
    pub max_read_length: usize,
    /// Duplication rate
    pub duprate: f64,
}

/// Output data from StructuralVariantsProcessor (same structure, possibly modified)
pub type ProcessedVariationData = RealignedVariationData;

/// StructuralVariantsProcessor - Processes variation data for SVs and adjusts SNV counts
pub struct StructuralVariantsProcessor {
    /// Reference sequence
    reference_seq: Vec<u8>,
    /// Reference start position (1-based)
    ref_start: i64,
}

impl StructuralVariantsProcessor {
    /// Create a new StructuralVariantsProcessor
    pub fn new(reference_seq: Vec<u8>, ref_start: i64) -> Self {
        StructuralVariantsProcessor {
            reference_seq,
            ref_start,
        }
    }

    /// Process realigned variation data
    /// 
    /// This is equivalent to Java StructuralVariantsProcessor.process()
    pub fn process(&self, mut data: RealignedVariationData) -> ProcessedVariationData {
        // If SV is enabled, find structural variants
        if !instance().conf.disable_sv {
            self.find_all_svs(&mut data);
        }

        // Always adjust SNV counts from soft clips
        self.adj_snv(&mut data);

        data
    }

    /// Find all structural variants (DEL, INV, DUP)
    /// 
    /// Called when SV detection is enabled
    fn find_all_svs(&self, _data: &mut RealignedVariationData) {
        // SV detection is complex - for now, if someone tries to use it, panic
        unimplemented!(
            "Structural variant detection is not yet implemented in Rust port. \
             Set disable_sv=true in configuration to use Simple Mode."
        );
    }

    /// Adjust SNV counts from short soft-clipped reads
    /// 
    /// This always runs, even when SV detection is disabled.
    /// It looks at short soft clips (≤5 bp) and if the first base matches
    /// a known SNV at the adjacent position, adds the soft clip counts to that SNV.
    fn adj_snv(&self, data: &mut RealignedVariationData) {
        // Process 5' end soft clips
        self.adj_snv_5end(data);
        
        // Process 3' end soft clips
        self.adj_snv_3end(data);
    }

    /// Adjust SNV counts from 5' soft clips
    fn adj_snv_5end(&self, data: &mut RealignedVariationData) {
        // Collect positions to process (avoid borrowing issues)
        let positions: Vec<i64> = data.soft_clips_5end.keys().cloned().collect();

        for position in positions {
            let sclip = match data.soft_clips_5end.get_mut(&position) {
                Some(sc) => sc,
                None => continue,
            };

            // Skip if already used
            if sclip.used() {
                continue;
            }

            // Find consensus sequence from soft clip
            let seq = self.find_conseq(sclip);
            
            // Only process short soft clips (≤5 bp)
            if seq.len() > 5 {
                continue;
            }

            if seq.is_empty() {
                continue;
            }

            // Get first base
            let bp = seq[0];

            // Check if there's a matching SNV at the previous position
            let prev_pos = position - 1;
            
            // Create the variant description key for the first base
            let var_key = VarDesc::snv_key(bp);

            // Check if this base exists as a variant at the previous position
            if let Some(var_map) = data.non_insertion_variants.get_mut(&prev_pos) {
                if var_map.contains_key(&var_key) {
                    // Additional check: if seq length > 1, verify second base matches reference
                    // Java: if reference base is missing, treat as mismatch and skip
                    if seq.len() > 1 {
                        match self.get_ref_base(position - 2) {
                            Some(ref_base) => {
                                if ref_base != seq[1] {
                                    continue;
                                }
                            }
                            None => {
                                continue;
                            }
                        }
                    }

                    // Get the soft clip variant counts
                    let sclip = data.soft_clips_5end.get(&position).unwrap();
                    let vars_count = sclip.var.alt_depth;
                    let high_qual_cnt = sclip.var.high_qual_read_cnt;
                    let low_qual_cnt = sclip.var.low_qual_read_cnt;
                    let mean_pos = sclip.var.mean_pos;
                    let mean_qual = sclip.var.mean_qual;
                    let mean_mapq = sclip.var.mean_mapq;
                    let nm = sclip.var.nm;
                    let fwd_cnt = sclip.var.alt_depth_fwd;
                    let rev_cnt = sclip.var.alt_depth_rev;

                    // Adjust the variant counts
                    if let Some(variant) = var_map.get_mut(&var_key) {
                        adj_cnt(
                            variant,
                            vars_count,
                            high_qual_cnt,
                            low_qual_cnt,
                            mean_pos,
                            mean_qual,
                            mean_mapq,
                            nm,
                            fwd_cnt,
                            rev_cnt,
                        );
                    }

                    // Increment reference coverage
                    *data.ref_coverage.entry(prev_pos).or_insert(0) += vars_count;
                }
            }
        }
    }

    /// Adjust SNV counts from 3' soft clips
    fn adj_snv_3end(&self, data: &mut RealignedVariationData) {
        // Collect positions to process (avoid borrowing issues)
        let positions: Vec<i64> = data.soft_clips_3end.keys().cloned().collect();

        for position in positions {
            let sclip = match data.soft_clips_3end.get_mut(&position) {
                Some(sc) => sc,
                None => continue,
            };

            // Skip if already used
            if sclip.used() {
                continue;
            }

            // Find consensus sequence from soft clip
            let seq = self.find_conseq(sclip);
            
            // Only process short soft clips (≤5 bp)
            if seq.len() > 5 {
                continue;
            }

            if seq.is_empty() {
                continue;
            }

            // Get first base
            let bp = seq[0];

            // Create the variant description key for the first base
            let var_key = VarDesc::snv_key(bp);

            // Check if this base exists as a variant at this position
            if let Some(var_map) = data.non_insertion_variants.get_mut(&position) {
                if var_map.contains_key(&var_key) {
                    // Additional check: if seq length > 1, verify second base matches reference
                    // Java: if reference base is missing, treat as mismatch and skip
                    if seq.len() > 1 {
                        match self.get_ref_base(position + 1) {
                            Some(ref_base) => {
                                if ref_base != seq[1] {
                                    continue;
                                }
                            }
                            None => {
                                continue;
                            }
                        }
                    }

                    // Get the soft clip variant counts
                    let sclip = data.soft_clips_3end.get(&position).unwrap();
                    let vars_count = sclip.var.alt_depth;
                    let high_qual_cnt = sclip.var.high_qual_read_cnt;
                    let low_qual_cnt = sclip.var.low_qual_read_cnt;
                    let mean_pos = sclip.var.mean_pos;
                    let mean_qual = sclip.var.mean_qual;
                    let mean_mapq = sclip.var.mean_mapq;
                    let nm = sclip.var.nm;
                    let fwd_cnt = sclip.var.alt_depth_fwd;
                    let rev_cnt = sclip.var.alt_depth_rev;

                    // Adjust the variant counts
                    if let Some(variant) = var_map.get_mut(&var_key) {
                        adj_cnt(
                            variant,
                            vars_count,
                            high_qual_cnt,
                            low_qual_cnt,
                            mean_pos,
                            mean_qual,
                            mean_mapq,
                            nm,
                            fwd_cnt,
                            rev_cnt,
                        );
                    }

                    // Increment reference coverage
                    *data.ref_coverage.entry(position).or_insert(0) += vars_count;
                }
            }
        }
    }

    /// Find consensus sequence from soft clip data
    /// 
    /// Simplified version of Java findconseq()
    fn find_conseq(&self, sclip: &mut SoftClip) -> Vec<u8> {
        crate::variants::var_utils::find_conseq(sclip, 0)
    }

    /// Get reference base at a position (1-based)
    fn get_ref_base(&self, pos: i64) -> Option<u8> {
        if pos < self.ref_start {
            return None;
        }
        let idx = (pos - self.ref_start) as usize;
        self.reference_seq.get(idx).copied()
    }
}

/// Adjust variant counts by adding values from soft clip
/// 
/// Equivalent to Java VariationUtils.adjCnt()
fn adj_cnt(
    variant: &mut Variant,
    vars_count: usize,
    high_qual_cnt: usize,
    low_qual_cnt: usize,
    mean_pos: f64,
    mean_qual: f64,
    mean_mapq: f64,
    nm: f64,
    fwd_cnt: usize,
    rev_cnt: usize,
) {
    variant.alt_depth += vars_count;
    variant.extra_cnt += vars_count;
    variant.high_qual_read_cnt += high_qual_cnt;
    variant.low_qual_read_cnt += low_qual_cnt;
    variant.mean_pos += mean_pos;
    variant.mean_qual += mean_qual;
    variant.mean_mapq += mean_mapq;
    variant.nm += nm;
    variant.alt_depth_fwd += fwd_cnt;
    variant.alt_depth_rev += rev_cnt;
    variant.pstd = true;
    variant.qstd = true;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_processor_creation() {
        let processor = StructuralVariantsProcessor::new(
            b"ACGTACGT".to_vec(),
            100,
        );
        assert_eq!(processor.ref_start, 100);
    }

    #[test]
    fn test_get_ref_base() {
        let processor = StructuralVariantsProcessor::new(
            b"ACGTACGT".to_vec(),
            100,
        );
        
        assert_eq!(processor.get_ref_base(100), Some(b'A'));
        assert_eq!(processor.get_ref_base(101), Some(b'C'));
        assert_eq!(processor.get_ref_base(107), Some(b'T'));
        assert_eq!(processor.get_ref_base(108), None); // Out of bounds
        assert_eq!(processor.get_ref_base(99), None);  // Before start
    }

    #[test]
    fn test_adj_cnt() {
        let mut variant = Variant::default();
        
        adj_cnt(
            &mut variant,
            10,  // vars_count
            8,   // high_qual_cnt
            2,   // low_qual_cnt
            50.0, // mean_pos
            30.0, // mean_qual
            60.0, // mean_mapq
            1.0,  // nm
            6,    // fwd_cnt
            4,    // rev_cnt
        );

        assert_eq!(variant.alt_depth, 10);
        assert_eq!(variant.high_qual_read_cnt, 8);
        assert_eq!(variant.low_qual_read_cnt, 2);
        assert_eq!(variant.mean_pos, 50.0);
        assert_eq!(variant.alt_depth_fwd, 6);
        assert_eq!(variant.alt_depth_rev, 4);
        assert!(variant.pstd);
        assert!(variant.qstd);
    }

    #[test]
    fn test_empty_data_processing() {
        // This test only verifies the adj_snv logic with empty data
        // The process() method requires GlobalReadOnlyScope initialization
        let processor = StructuralVariantsProcessor::new(
            b"ACGT".to_vec(),
            1,
        );
        
        let mut data = RealignedVariationData::default();
        
        // Directly test adj_snv with empty data
        processor.adj_snv(&mut data);
        
        // Should return empty data without errors
        assert!(data.non_insertion_variants.is_empty());
        assert!(data.insertion_variants.is_empty());
    }
}
