use std::collections::{HashMap, VecDeque};

use rust_htslib::bam::record::Cigar;

use crate::variants::variants::Variant;

pub(crate) mod patterns;
pub mod reference;
pub mod region;

pub(crate) type VariantMap = HashMap<String, Variant>;

pub(crate) struct ModifiedCigar<'a> {
    pub(crate) align_start_pos: i64,
    pub(crate) cigar: VecDeque<Cigar>,
    pub(crate) query_seq: &'a [u8],
    pub(crate) query_qual: &'a [u8],
}

impl<'a> ModifiedCigar<'a> {
    pub(crate) fn new(
        align_start_pos: i64,
        cigar: VecDeque<Cigar>,
        query_seq: &'a [u8],
        query_qual: &'a [u8],
    ) -> Self {
        Self {
            align_start_pos,
            cigar,
            query_seq,
            query_qual,
        }
    }
}
