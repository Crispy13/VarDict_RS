use std::collections::HashMap;

use anyhow::{Error, anyhow};
use crackle_kit::tracing::{Level, event};

use crate::{
    data::patterns::{BEGIN_MINUS_NUMBER, UP_NUMBER_END, get_cap_group},
    prelude::SmallVecBytes,
    utils::SliceExt,
    variants::{var_utils::get_variants_from_map, variants::Variant},
};

/// Result of finding 3'/5' end matches between two sequences
#[derive(Debug, Clone)]
pub struct Match35 {
    pub matched_5_end: usize,
    pub matched_3_end: usize,
    pub max_matched_length: usize,
}

pub struct VariantRealigner {
    non_insertion_vars: HashMap<i64, HashMap<String, Variant>>,
    ref_coverage: HashMap<i64, usize>,
}

impl VariantRealigner {
    /// Check if a sequence has low complexity (>75% of one base or <3 different bases)
    pub fn is_low_complex_seq(seq: &str) -> bool {
        let len = seq.len();
        if len == 0 {
            return true;
        }

        let a = Self::count(seq, 'A');
        if a as f64 / len as f64 > 0.75 {
            return true;
        }

        let t = Self::count(seq, 'T');
        if t as f64 / len as f64 > 0.75 {
            return true;
        }

        let g = Self::count(seq, 'G');
        if g as f64 / len as f64 > 0.75 {
            return true;
        }

        let c = Self::count(seq, 'C');
        if c as f64 / len as f64 > 0.75 {
            return true;
        }

        // Count how many different bases are present
        let mut nt_cnt = 0;
        if a > 0 {
            nt_cnt += 1;
        }
        if t > 0 {
            nt_cnt += 1;
        }
        if g > 0 {
            nt_cnt += 1;
        }
        if c > 0 {
            nt_cnt += 1;
        }

        nt_cnt < 3
    }

    /// Count occurrences of a character in a string
    fn count(seq: &str, ch: char) -> usize {
        seq.chars().filter(|&c| c == ch).count()
    }

    /// Find if two sequences match with no more than MM mismatches
    /// and total mismatches no more than 15% of length
    pub fn is_match(seq1: &str, seq2: &str, dir: i32) -> bool {
        Self::is_match_with_threshold(seq1, seq2, dir, 3)
    }

    /// Find if two sequences match with specified MM threshold
    pub fn is_match_with_threshold(seq1: &str, seq2: &str, dir: i32, mm_threshold: usize) -> bool {
        // Remove special characters from seq2
        let seq2_clean: String = seq2.chars().filter(|c| *c != '#' && *c != '^').collect();

        let mut mm = 0;
        let min_len = seq1.len().min(seq2_clean.len());

        for n in 0..min_len {
            let seq1_ch = seq1.chars().nth(n).unwrap_or('N');
            let seq2_idx = if dir == 1 {
                n
            } else {
                // For reverse direction, index from the end
                seq2_clean.len() - 1 - n
            };
            let seq2_ch = seq2_clean.chars().nth(seq2_idx).unwrap_or('N');

            if seq1_ch != seq2_ch {
                mm += 1;
            }
        }

        (mm <= mm_threshold) && (mm as f64 / seq1.len() as f64) < 0.15
    }

    /// Find 3'/5' end matches between two sequences
    /// Returns Match35 with (matched_5_end, matched_3_end, max_matched_length)
    pub fn find_35_match(seq5: &str, seq3: &str) -> Match35 {
        const LONG_MISMATCH: usize = 2;

        let mut max_matched_length = 0;
        let mut b3 = 0;
        let mut b5 = 0;

        // Convert to char vectors for efficient indexing
        let seq5_chars: Vec<char> = seq5.chars().collect();
        let seq3_chars: Vec<char> = seq3.chars().collect();

        // Ensure sequences are long enough for comparison
        if seq5_chars.len() <= 8 || seq3_chars.len() <= 8 {
            return Match35 {
                matched_5_end: b5,
                matched_3_end: b3,
                max_matched_length,
            };
        }

        // Java: for (int i = 0; i < seq5.length() - 8; i++)
        for i in 0..(seq5_chars.len() - 8) {
            // Java: for (int j = 1; j < seq3.length() - 8; j++)
            for j in 1..(seq3_chars.len() - 8) {
                let mut num_mismatch = 0;
                let mut total_length = 0;

                // Java: while (totalLength + j <= seq3.length() && i + totalLength <= seq5.length())
                while total_length + j <= seq3_chars.len() && i + total_length <= seq5_chars.len() {
                    // Java: substr(seq3, -j - totalLength, 1)
                    // which is seq3[len - j - totalLength]
                    let seq3_pos = seq3_chars.len() - j - total_length;
                    let seq5_pos = i + total_length;

                    // Java's substr returns empty string for out-of-bounds, which doesn't equal any char
                    // So out-of-bounds is treated as a mismatch
                    let is_mismatch = if seq3_pos >= seq3_chars.len() || seq5_pos >= seq5_chars.len() {
                        true // Out of bounds = mismatch (like Java's empty string comparison)
                    } else {
                        seq3_chars[seq3_pos] != seq5_chars[seq5_pos]
                    };

                    // Java: if (!substr(seq3, -j - totalLength, 1).equals(substr(seq5, i + totalLength, 1)))
                    if is_mismatch {
                        num_mismatch += 1;
                    }

                    // Java: if (numberOfMismatch > longMismatch)
                    if num_mismatch > LONG_MISMATCH {
                        break;
                    }

                    total_length += 1;
                }

                // Java: if (totalLength - numberOfMismatch > maxMatchedLength
                //           && totalLength - numberOfMismatch > 8
                //           && numberOfMismatch / (double) totalLength < 0.1d
                //           && (totalLength + j >= seq3.length() || i + totalLength >= seq5.length()))
                let matched_quality = total_length.saturating_sub(num_mismatch);
                let is_end_match = (total_length + j >= seq3_chars.len()) || (i + total_length >= seq5_chars.len());
                let error_rate = if total_length > 0 {
                    num_mismatch as f64 / total_length as f64
                } else {
                    0.0
                };

                if matched_quality > max_matched_length
                    && matched_quality > 8
                    && error_rate < 0.1
                    && is_end_match
                {
                    max_matched_length = matched_quality;
                    b3 = j;
                    b5 = i;
                    // Java returns immediately when first match is found
                    return Match35 {
                        matched_5_end: b5,
                        matched_3_end: b3,
                        max_matched_length,
                    };
                }
            }
        }

        Match35 {
            matched_5_end: b5,
            matched_3_end: b3,
            max_matched_length,
        }
    }

    fn realign_del(
        &mut self,
        bam_parameters: &[&str],
        pos_to_del_count: &HashMap<i64, HashMap<String, usize>>,
    ) -> Result<(), Error> {
        let _bams: &[&str] = bam_parameters;

        let sorted_pos_to_del = fill_and_sort_tmp(pos_to_del_count);

        for tpl in sorted_pos_to_del {
            let p = tpl.pos;
            let vn = tpl.desc_string;
            let _del_cnt = tpl.count;

            event!(
                Level::INFO,
                "  Realigndel for: {p} {vn} {_del_cnt} cov: {}",
                self.ref_coverage.get(&p).copied().unwrap_or(0)
            );

            let mut del_len = 0i64;

            if let Some(cap) = BEGIN_MINUS_NUMBER.captures(&vn) {
                if let Ok(len_str) = get_cap_group!(cap, 1) {
                    del_len = len_str.as_str().parse::<i64>().unwrap_or(0);
                }
            }

            if let Some(cap) = UP_NUMBER_END.captures(&vn) {
                if let Ok(len_str) = get_cap_group!(cap, 1) {
                    del_len += len_str.as_str().parse::<i64>().unwrap_or(0);
                }
            }

            let _extrains = String::new();
            let _extra = String::new();
            let _inv5 = String::new();
            let _inv3 = String::new();

            // TODO: Complete realign_del logic for simple mode
        }

        Ok(())
    }
}

/// Fill and sort temporary tuple structure from hash tables of insertion or deletions
fn fill_and_sort_tmp(
    changes: &HashMap<i64, HashMap<String, usize>>,
) -> Vec<SortPositionDescription> {
    let mut var_info: Vec<SortPositionDescription> = changes
        .iter()
        .flat_map(|(&pos, v)| {
            v.iter()
                .map(move |(desc, &cnt)| SortPositionDescription::new(pos, desc.clone(), cnt))
        })
        .collect();

    var_info.sort_by(|o1, o2| {
        o2.count
            .cmp(&o1.count)
            .then(o1.pos.cmp(&o2.pos))
            .then(o2.desc_string.cmp(&o1.desc_string))
    });

    var_info
}

pub(crate) struct SortPositionDescription {
    pos: i64,
    desc_string: String,
    count: usize,
}

impl SortPositionDescription {
    fn new(pos: i64, desc_string: String, count: usize) -> Self {
        Self {
            pos,
            desc_string,
            count,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_low_complex_seq() {
        // Test cases from Java VariationRealignerTest
        assert!(VariantRealigner::is_low_complex_seq("AAAAAAAAA"));
        assert!(VariantRealigner::is_low_complex_seq("ATATATATATAT"));
        assert!(VariantRealigner::is_low_complex_seq("CCCCCCCCGA"));
        assert!(!VariantRealigner::is_low_complex_seq("ACGTACGTACGT"));
        assert!(!VariantRealigner::is_low_complex_seq("CCGTAACGGGGT"));
    }

    #[test]
    fn test_is_match() {
        // Test cases from Java VariationRealignerTest
        assert!(VariantRealigner::is_match("AAAAAAAAA", "AAAAAAAAA", 1));
        assert!(VariantRealigner::is_match("AAAAAAAAA", "AAAAAAAAA", -1));
        assert!(VariantRealigner::is_match("ACGTACGTACGTACGT", "AAGTACTTACGTACGT", 1));
        assert!(VariantRealigner::is_match("ACGTACGTACGTACGT", "AAGTACTTACGT", 1));

        assert!(!VariantRealigner::is_match("ACGTACGTACGT", "AAGTACTTACGT", 1));
        assert!(!VariantRealigner::is_match("ACGTCAGCAT", "ACGACTGACT", 1));
    }

    #[test]
    fn test_find_35_match() {
        // Test cases from Java VariationRealignerTest
        let result1 = VariantRealigner::find_35_match("ACGTACGTACGTACGTACGTACGT", "TGCATGCATGCATGCATGCATGCA");
        assert_eq!(result1.matched_5_end, 0);
        assert_eq!(result1.matched_3_end, 1);
        assert_eq!(result1.max_matched_length, 24);

        let result2 = VariantRealigner::find_35_match("ACGTACGTACGTACGTACGTACGT", "TGCATGCATGCCTGCATGCATGCA");
        assert_eq!(result2.matched_5_end, 0);
        assert_eq!(result2.matched_3_end, 1);
        assert_eq!(result2.max_matched_length, 23);

        let result3 = VariantRealigner::find_35_match("ACGTACGTACGTACGTACGTACGT", "TGCATGCATGCATGGGTGCATGCA");
        assert_eq!(result3.matched_5_end, 0);
        assert_eq!(result3.matched_3_end, 1);
        assert_eq!(result3.max_matched_length, 22);

        let result4 = VariantRealigner::find_35_match("ACGTACGTACGTACGTACGTACGT", "TGCATGCATGCATGGGTGCAAAAA");
        assert_eq!(result4.matched_5_end, 0);
        assert_eq!(result4.matched_3_end, 13);
        assert_eq!(result4.max_matched_length, 12);

        let result5 = VariantRealigner::find_35_match("ACGTACGTACGTACGTACGTACGT", "AAAATGCATGCATGGGTGCAAAAA");
        assert_eq!(result5.matched_5_end, 14);
        assert_eq!(result5.matched_3_end, 11);
        assert_eq!(result5.max_matched_length, 10);
    }
}
