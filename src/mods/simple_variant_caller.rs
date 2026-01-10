//! Simple Variant Caller - Direct BAM to variations conversion for Simple Mode
//!
//! This module provides a simplified variant calling pipeline that:
//! 1. Reads BAM records from a region
//! 2. Parses CIGAR strings to extract variations
//! 3. Collects coverage and variation statistics
//!
//! This is a simplified version suitable for Simple Mode only.

use std::collections::HashMap;

use rust_htslib::bam::{Record, record::Cigar};

use crate::mods::to_vars_builder::VariationData;

/// Result of processing a single read
#[derive(Debug, Default)]
pub struct ReadVariations {
    /// Variations found: (position, variant_key, data)
    pub variations: Vec<(i64, String, VariationData)>,
    /// Positions covered by this read: (position, count=1)
    pub coverage_positions: Vec<i64>,
}

/// Simple variant caller that processes BAM records
pub struct SimpleVariantCaller {
    /// Minimum base quality to consider a variant
    min_base_quality: u8,
    /// Minimum mapping quality to process a read
    min_mapping_quality: u8,
    /// Reference sequence for the current region
    reference: Vec<u8>,
    /// Start position of the reference (1-based)
    ref_start: i64,
}

impl SimpleVariantCaller {
    /// Create a new SimpleVariantCaller
    pub fn new(min_base_quality: u8, min_mapping_quality: u8) -> Self {
        SimpleVariantCaller {
            min_base_quality,
            min_mapping_quality,
            reference: Vec::new(),
            ref_start: 0,
        }
    }

    /// Set the reference sequence for the current region
    pub fn set_reference(&mut self, sequence: Vec<u8>, start: i64) {
        self.reference = sequence;
        self.ref_start = start;
    }

    /// Get the reference base at a genomic position (1-based)
    fn get_ref_base(&self, pos: i64) -> Option<u8> {
        if pos < self.ref_start {
            return None;
        }
        let idx = (pos - self.ref_start) as usize;
        self.reference.get(idx).copied()
    }

    /// Process a single BAM record and extract variations
    pub fn process_record(&self, record: &Record) -> ReadVariations {
        let mut result = ReadVariations::default();

        // Skip unmapped reads
        if record.is_unmapped() {
            return result;
        }

        // Skip reads below mapping quality threshold
        if record.mapq() < self.min_mapping_quality {
            return result;
        }

        // Get read data
        let seq = record.seq();
        let qual = record.qual();
        let is_reverse = record.is_reverse();
        let mapq = record.mapq();
        let read_name = String::from_utf8_lossy(record.qname()).to_string();

        // Parse CIGAR and extract variations
        let cigar = record.cigar();
        
        // Position in reference (0-based from BAM, convert to 1-based)
        let mut ref_pos = record.pos() + 1;
        // Position in read sequence (0-based)
        let mut read_pos: usize = 0;

        for cigar_op in cigar.iter() {
            match cigar_op {
                Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                    // Process each base in the match
                    for i in 0..(*len as usize) {
                        let current_ref_pos = ref_pos + i as i64;
                        let current_read_pos = read_pos + i;

                        // Add to coverage
                        result.coverage_positions.push(current_ref_pos);

                        // Check for mismatch
                        if let Some(ref_base) = self.get_ref_base(current_ref_pos) {
                            // Seq type uses indexing, check bounds first
                            if current_read_pos >= seq.len() {
                                continue;
                            }
                            let read_base = seq[current_read_pos];
                            let base_qual = qual.get(current_read_pos).copied().unwrap_or(0);

                            // Skip low quality bases
                            if base_qual < self.min_base_quality {
                                continue;
                            }

                            // Convert bases to uppercase for comparison
                            let ref_base_upper = ref_base.to_ascii_uppercase();
                            let read_base_upper = read_base.to_ascii_uppercase();

                            if ref_base_upper != read_base_upper && read_base_upper != b'N' {
                                // Found a mismatch (SNV)
                                let variant_key = format!(
                                    "{}>{}",
                                    char::from(ref_base_upper),
                                    char::from(read_base_upper)
                                );

                                let data = VariationData {
                                    position_in_read: current_read_pos as u32,
                                    quality: base_qual,
                                    mapping_quality: mapq,
                                    is_reverse,
                                    read_id: read_name.clone(),
                                };

                                result.variations.push((current_ref_pos, variant_key, data));
                            }
                        }
                    }

                    ref_pos += *len as i64;
                    read_pos += *len as usize;
                }

                Cigar::Ins(len) => {
                    // Insertion: bases in read not in reference
                    let ins_start = read_pos;
                    let ins_end = read_pos + *len as usize;

                    // Get inserted sequence
                    let mut ins_seq = String::new();
                    let mut min_qual = u8::MAX;
                    for i in ins_start..ins_end {
                        if i >= seq.len() {
                            break;
                        }
                        let base = seq[i];
                        ins_seq.push(char::from(base.to_ascii_uppercase()));
                        if let Some(&q) = qual.get(i) {
                            min_qual = min_qual.min(q);
                        }
                    }

                    if min_qual >= self.min_base_quality {
                        let variant_key = format!("+{}", ins_seq);

                        let data = VariationData {
                            position_in_read: ins_start as u32,
                            quality: min_qual,
                            mapping_quality: mapq,
                            is_reverse,
                            read_id: read_name.clone(),
                        };

                        // Insertion is reported at the position before it
                        result.variations.push((ref_pos - 1, variant_key, data));
                    }

                    read_pos += *len as usize;
                    // ref_pos doesn't change for insertion
                }

                Cigar::Del(len) => {
                    // Deletion: bases in reference not in read
                    // Get the deleted sequence from reference
                    let mut del_seq = String::new();
                    for i in 0..(*len as i64) {
                        if let Some(base) = self.get_ref_base(ref_pos + i) {
                            del_seq.push(char::from(base.to_ascii_uppercase()));
                        }
                    }

                    // Use the quality of the flanking base
                    let flank_qual = if read_pos > 0 {
                        qual.get(read_pos - 1).copied().unwrap_or(30)
                    } else {
                        qual.get(read_pos).copied().unwrap_or(30)
                    };

                    if flank_qual >= self.min_base_quality {
                        let variant_key = format!("-{}", del_seq);

                        let data = VariationData {
                            position_in_read: read_pos as u32,
                            quality: flank_qual,
                            mapping_quality: mapq,
                            is_reverse,
                            read_id: read_name.clone(),
                        };

                        result.variations.push((ref_pos, variant_key, data));
                    }

                    ref_pos += *len as i64;
                    // read_pos doesn't change for deletion
                }

                Cigar::SoftClip(len) => {
                    // Soft clipped bases are in the read but not aligned
                    read_pos += *len as usize;
                }

                Cigar::HardClip(_) => {
                    // Hard clipped bases are not in the read sequence
                }

                Cigar::RefSkip(len) => {
                    // Reference skip (e.g., intron in RNA-seq)
                    ref_pos += *len as i64;
                }

                Cigar::Pad(_) => {
                    // Padding - ignore
                }
            }
        }

        result
    }

    /// Process multiple records and aggregate results
    pub fn process_records<'a, I>(
        &self,
        records: I,
    ) -> (Vec<(i64, String, VariationData)>, HashMap<i64, usize>)
    where
        I: Iterator<Item = &'a Record>,
    {
        let mut all_variations = Vec::new();
        let mut coverage: HashMap<i64, usize> = HashMap::new();

        for record in records {
            let read_result = self.process_record(record);

            // Accumulate variations
            all_variations.extend(read_result.variations);

            // Accumulate coverage
            for pos in read_result.coverage_positions {
                *coverage.entry(pos).or_insert(0) += 1;
            }
        }

        (all_variations, coverage)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_variant_caller_creation() {
        let caller = SimpleVariantCaller::new(25, 20);
        assert_eq!(caller.min_base_quality, 25);
        assert_eq!(caller.min_mapping_quality, 20);
    }

    #[test]
    fn test_get_ref_base() {
        let mut caller = SimpleVariantCaller::new(25, 20);
        caller.set_reference(b"ACGTACGT".to_vec(), 100);

        assert_eq!(caller.get_ref_base(100), Some(b'A'));
        assert_eq!(caller.get_ref_base(101), Some(b'C'));
        assert_eq!(caller.get_ref_base(107), Some(b'T'));
        assert_eq!(caller.get_ref_base(99), None);  // Before start
        assert_eq!(caller.get_ref_base(108), None); // After end
    }
}
