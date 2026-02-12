use std::collections::HashMap;

use crackle_kit::nuc_base_map::NucBaseMap;

use crate::{
    conf::Configuration,
    scopedata::global_read_only_scope::instance,
    utils::SliceExt,
    variants::variants::{SoftClip, VarDesc, Variant},
};

pub(crate) fn get_variants_from_map<'a>(
    var_map: &'a mut HashMap<i64, HashMap<VarDesc, Variant>>,
    start: i64,
    var_desc: &VarDesc,
) -> &'a mut Variant {
    let pos_map = var_map
        .entry(start)
        .or_insert_with(|| HashMap::with_capacity(1));

    if pos_map.contains_key(var_desc) {
        pos_map.get_mut(var_desc).unwrap()
    } else {
        pos_map.entry(var_desc.clone()).or_default()
    }

    // // Get a raw pointer to avoid the borrow checker blocking the 'None' branch
    // let var_ptr: *const Variant = match pos_map.get(desc_string) {
    //     // Fast Path: Found it. 1 Hash. 0 Allocations.
    //     Some(var) => var,

    //     // Slow Path: Not found.
    //     // We use 'entry' here to insert AND get the reference in one go.
    //     // 2 Hashes total (1 check above + 1 insert here).
    //     None => pos_map.entry(desc_string.to_string()).or_default(),
    // };

    // // SAFETY:
    // // 1. We know 'var_ptr' points to valid memory inside 'pos_map'.
    // // 2. We do not mutate 'pos_map' again after obtaining this pointer.
    // // 3. The returned reference lifetime 'a is tied to the map, preventing
    // //    the caller from invalidating the pointer while holding the reference.
    // unsafe { &*var_ptr }
}

/// Get `Variant` from `SoftClip.seq` field
pub(crate) fn get_variation_from_seq(softclip: &mut SoftClip, idx: usize, base: u8) -> &mut Variant {
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
        for &base in [b'A', b'T', b'G', b'C', b'N'].iter() {
            let Some(&count) = base_counts.get(base) else { continue; };
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
    seq.len() >= 8 && seq.get(1..8).map(|s| s.iter().all(|&b| b == b'A')).unwrap_or(false)
}

fn has_poly_t7(seq: &[u8]) -> bool {
    seq.len() >= 8 && seq.get(1..8).map(|s| s.iter().all(|&b| b == b'T')).unwrap_or(false)
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
    if a > 0 { ntcnt += 1; }
    if t > 0 { ntcnt += 1; }
    if g > 0 { ntcnt += 1; }
    if c > 0 { ntcnt += 1; }

    if a as f64 / len as f64 > 0.75 { return true; }
    if t as f64 / len as f64 > 0.75 { return true; }
    if g as f64 / len as f64 > 0.75 { return true; }
    if c as f64 / len as f64 > 0.75 { return true; }

    ntcnt < 3
}

fn count_base(seq: &[u8], base: u8) -> usize {
    seq.iter().filter(|&&b| b == base).count()
}

pub(crate) struct HomoPolymerChecker(u8);

impl HomoPolymerChecker {
    pub(crate) fn new() -> Self {
        Self(0)
    }

    pub(crate) fn record_base(&mut self, b: u8) {
        let bit = match b.to_ascii_lowercase() {
            b'a' => 0b1,
            b'c' => 0b10,
            b'g' => 0b100,
            b't' => 0b1000,
            b'n' => 0b10000,
            _ => panic!("Invalid base: {}", b as char),
        };

        self.0 |= bit
    }

    pub(crate) fn is_homopolymer(&self) -> bool {
        if self.0 == 0 {
            panic!("No base is recorded.")
        }
        self.0.count_ones() == 1
    }
}
