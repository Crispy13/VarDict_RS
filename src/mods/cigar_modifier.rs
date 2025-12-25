use std::{borrow::Cow, collections::HashMap};

use rust_htslib::bam::{Record, record::CigarStringView};

use crate::data::{reference::Reference, region::Region};

pub struct CigarModifier<'a> {
    align_start_pos: usize,
    cigar_str: Cow<'a, CigarStringView>,
    original_cigar: &'a CigarStringView,
    query_sequence: &'a [u8],
    query_quality: &'a [u8],
    ref_data: &'a Reference,
    indel: u32,
    max_read_length: usize,
    region: &'a Region,
}

impl<'a> CigarModifier<'a> {
    pub(crate) fn new(
        align_start_pos: usize,
        cigar_str: &'a CigarStringView,
        query_sequence: &'a [u8],
        query_quality: &'a [u8],
        ref_data: &'a Reference,
        indel: u32,
        max_read_length: usize,
        region: &'a Region,
    ) -> Self {
        Self {
            align_start_pos,
            cigar_str: Cow::Borrowed(cigar_str),
            original_cigar: cigar_str,
            query_sequence,
            query_quality,
            ref_data,
            indel,
            max_read_length,
            region,
        }
    }
}
