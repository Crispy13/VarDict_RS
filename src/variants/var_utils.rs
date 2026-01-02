use std::collections::HashMap;

use crate::{utils::SliceExt, variants::variants::{VarDesc, Variant}};

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
