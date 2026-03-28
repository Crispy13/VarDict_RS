use std::collections::HashMap;

use crackle_kit::nuc_base_map::NucBaseMap;

use crate::{
    conf::Configuration,
    prelude::{InnerMap, LibDefaultHasher},
    scopedata::global_read_only_scope::instance,
    variants::variants::{SoftClip, VarDesc, Variant},
};

pub(crate) fn get_variants_from_map<'a>(
    var_map: &'a mut HashMap<i64, InnerMap<VarDesc, Variant>, LibDefaultHasher>,
    start: i64,
    var_desc: &VarDesc,
) -> &'a mut Variant {
    let pos_map = var_map.entry(start).or_insert_with(Default::default);

    get_variant_from_pos_map(pos_map, var_desc)
}

pub(crate) fn get_variant_from_pos_map<'a>(
    pos_map: &'a mut InnerMap<VarDesc, Variant>,
    var_desc: &VarDesc,
) -> &'a mut Variant {
    let variant_ptr = if let Some(variant) = pos_map.get_mut(var_desc) {
        variant as *mut Variant
    } else {
        pos_map.entry(var_desc.clone()).or_default() as *mut Variant
    };

    // SAFETY: `variant_ptr` always points to an entry stored inside `pos_map`.
    // We do not mutate `pos_map` again after taking the pointer, and the returned
    // reference lifetime is still tied to the mutable borrow of `pos_map`.
    unsafe { &mut *variant_ptr }
}

/// Get `Variant` from `SoftClip.seq` field
pub(crate) fn get_variation_from_seq(
    softclip: &mut SoftClip,
    idx: usize,
    base: u8,
) -> &mut Variant {
    softclip
        .seq
        .entry(idx)
        .or_insert_with(|| NucBaseMap::default())
        .get_checked_or_insert_with(base, || Variant::default())
        .unwrap()
}

/// Find consensus sequence from soft clip data (ported from Java VariationUtils.findconseq)
pub(crate) fn find_conseq(softclip: &mut SoftClip, dir: i32) -> Vec<u8> {
    if softclip.consensus_seq_is_set() {
        return softclip.consensus_seq().to_vec();
    }

    let mut total = 0usize;
    let mut matched = 0usize;
    let mut seq: Vec<u8> = Vec::new();
    let mut flag = false;

    let nt_size = softclip.nt.len();

    for (pos, base_counts) in &softclip.nt {
        let mut max_count = 0usize;
        let mut max_quality = 0.0f64;
        let mut chosen_base: Option<u8> = None;
        let mut total_count = 0usize;

        let idx = *pos as usize;
        // Java iterates TreeMap<Character, Integer> entries in natural key order.
        // For nucleotide keys this is: A, C, G, N, T.
        for &base in [b'A', b'C', b'G', b'N', b'T'].iter() {
            let Some(&count) = base_counts.get(base) else {
                continue;
            };
            total_count += count;

            let mut choose = count > max_count;
            if !choose {
                if let Some(seq_map) = softclip.seq.get(&idx) {
                    if let Some(var) = seq_map.get(base) {
                        if var.mean_qual > max_quality {
                            choose = true;
                        }
                    }
                }
            }

            if choose {
                max_count = count;
                chosen_base = Some(base);
                if let Some(seq_map) = softclip.seq.get(&idx) {
                    if let Some(var) = seq_map.get(base) {
                        max_quality = var.mean_qual;
                    }
                }
            }
        }

        if *pos == 3
            && nt_size >= 6
            && softclip.var.alt_depth > 0
            && (total_count as f64 / softclip.var.alt_depth as f64) < 0.2
            && total_count <= 2
        {
            break;
        }

        if total_count > 0 {
            if (total_count - max_count > 2 || max_count <= total_count - max_count)
                && (max_count as f64 / total_count as f64) < 0.8
            {
                if flag {
                    break;
                }
                flag = true;
            }

            total += total_count;
            matched += max_count;
            if let Some(base) = chosen_base {
                seq.push(base);
            }
        }
    }

    let mut seq_out = Vec::new();
    let seq_len = seq.len();
    if total > 0
        && (matched as f64 / total as f64) > 0.9
        && (seq_len as f64 / 1.5) > (nt_size as f64 - seq_len as f64)
        && ((seq_len as f64 / nt_size as f64) > 0.8
            || nt_size.saturating_sub(seq_len) < 10
            || seq_len > 25)
    {
        seq_out = seq;
    }

    if !seq_out.is_empty() && seq_out.len() > Configuration::SEED_2 as usize {
        if has_poly_a7(&seq_out) || has_poly_t7(&seq_out) || is_low_complex_seq(&seq_out) {
            softclip.mark_used();
        }
    }

    if !seq_out.is_empty() && seq_out.len() >= Configuration::ADSEED as usize {
        let seed_len = Configuration::ADSEED as usize;
        let seed = &seq_out[..seed_len];
        let seed_str = String::from_utf8_lossy(seed).to_string();

        if dir == 3 {
            if instance().adaptor_forward.contains_key(&seed_str) {
                seq_out.clear();
            }
        } else if dir == 5 {
            let rev_seed: Vec<u8> = seed.iter().rev().copied().collect();
            let rev_seed_str = String::from_utf8_lossy(&rev_seed).to_string();
            if instance().adaptor_reverse.contains_key(&rev_seed_str) {
                seq_out.clear();
            }
        }
    }

    softclip.set_consensus_seq(seq_out.clone());
    seq_out
}

#[inline]
pub(crate) fn is_has_and_not_equals(b: u8, contig_ref_seq: &[u8], index: usize) -> bool {
    match contig_ref_seq.get(index) {
        Some(&rb) => b != rb,
        None => false,
    }
}

#[inline]
pub(crate) fn is_has_and_equals(b: u8, contig_ref_seq: &[u8], index: usize) -> bool {
    match contig_ref_seq.get(index) {
        Some(&rb) => b == rb,
        None => false,
    }
}

#[inline]
pub(crate) fn is_has_and_equals_ref_and_seq_base(
    contig_ref_seq: &[u8],
    index1: usize,
    seq: &[u8],
    index2: usize,
) -> bool {
    match contig_ref_seq.get(index1) {
        Some(&rb) => rb == *seq.get(index2).unwrap(),
        None => false,
    }
}

#[inline]
pub(crate) fn is_has_and_not_equals_ref_and_seq_base(
    contig_ref_seq: &[u8],
    index1: usize,
    seq: &[u8],
    index2: usize,
) -> bool {
    match contig_ref_seq.get(index1) {
        Some(&rb) => rb != *seq.get(index2).unwrap(),
        None => false,
    }
}

#[inline]
/// Check the bases of the two indices are the same
pub(crate) fn is_has_and_equals_two_index(
    index1: usize,
    contig_ref_seq: &[u8],
    index2: usize,
) -> bool {
    contig_ref_seq
        .get(index1)
        .and_then(|b1| contig_ref_seq.get(index2).and_then(|b2| Some(b1 == b2)))
        .unwrap_or(false)
}

fn has_poly_a7(seq: &[u8]) -> bool {
    seq.len() >= 8
        && seq
            .get(1..8)
            .map(|s| s.iter().all(|&b| b == b'A'))
            .unwrap_or(false)
}

fn has_poly_t7(seq: &[u8]) -> bool {
    seq.len() >= 8
        && seq
            .get(1..8)
            .map(|s| s.iter().all(|&b| b == b'T'))
            .unwrap_or(false)
}

fn is_low_complex_seq(seq: &[u8]) -> bool {
    let len = seq.len();
    if len == 0 {
        return true;
    }

    let a = count_base(seq, b'A');
    let t = count_base(seq, b'T');
    let g = count_base(seq, b'G');
    let c = count_base(seq, b'C');

    let mut ntcnt = 0;
    if a > 0 {
        ntcnt += 1;
    }
    if t > 0 {
        ntcnt += 1;
    }
    if g > 0 {
        ntcnt += 1;
    }
    if c > 0 {
        ntcnt += 1;
    }

    if a as f64 / len as f64 > 0.75 {
        return true;
    }
    if t as f64 / len as f64 > 0.75 {
        return true;
    }
    if g as f64 / len as f64 > 0.75 {
        return true;
    }
    if c as f64 / len as f64 > 0.75 {
        return true;
    }

    ntcnt < 3
}

fn count_base(seq: &[u8], base: u8) -> usize {
    seq.iter().filter(|&&b| b == base).count()
}

pub(crate) struct HomoPolymerChecker {
    first_base: Option<u8>,
    has_multiple_bases: bool,
}

impl HomoPolymerChecker {
    pub(crate) fn new() -> Self {
        Self {
            first_base: None,
            has_multiple_bases: false,
        }
    }

    pub(crate) fn record_base(&mut self, b: u8) {
        let normalized_base = b.to_ascii_uppercase();

        match self.first_base {
            Some(first_base) if first_base != normalized_base => {
                self.has_multiple_bases = true;
            }
            Some(_) => {}
            None => {
                self.first_base = Some(normalized_base);
            }
        }
    }

    pub(crate) fn is_homopolymer(&self) -> bool {
        !self.has_multiple_bases
    }
}

#[cfg(test)]
mod tests {
    use super::HomoPolymerChecker;

    #[test]
    fn homopolymer_checker_normalizes_case() {
        let mut checker = HomoPolymerChecker::new();
        checker.record_base(b'a');
        checker.record_base(b'A');

        assert!(checker.is_homopolymer());
    }

    #[test]
    fn homopolymer_checker_detects_mixed_standard_bases() {
        let mut checker = HomoPolymerChecker::new();
        checker.record_base(b'A');
        checker.record_base(b'C');

        assert!(!checker.is_homopolymer());
    }

    #[test]
    fn homopolymer_checker_accepts_ambiguous_iupac_bases() {
        let mut checker = HomoPolymerChecker::new();
        checker.record_base(b'Y');
        checker.record_base(b'Y');

        assert!(checker.is_homopolymer());
    }

    #[test]
    fn homopolymer_checker_marks_ambiguous_and_standard_mix_as_non_homopolymer() {
        let mut checker = HomoPolymerChecker::new();
        checker.record_base(b'Y');
        checker.record_base(b'C');

        assert!(!checker.is_homopolymer());
    }

    #[test]
    fn homopolymer_checker_without_bases_is_non_expanding() {
        let checker = HomoPolymerChecker::new();

        assert!(checker.is_homopolymer());
    }
}
