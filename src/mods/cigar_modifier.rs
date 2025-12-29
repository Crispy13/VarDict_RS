use std::{
    borrow::Cow,
    collections::{HashMap, VecDeque},
};

use anyhow::{Error, anyhow};
use crackle_kit::{
    data::bases::rev_comp::RevComplementor,
    tracing::{Level, event},
};
use rust_htslib::bam::{
    Record,
    record::{Cigar, CigarString, CigarStringView},
};

use crate::{
    conf::Configuration,
    data::{reference::Reference, region::Region},
    scopedata::global_read_only_scope::{INSTANCE, instance},
    utils::{BytesExt, SliceExt, SliceExt2},
};

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
    rev_complementor: &'a mut RevComplementor,
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
        rev_complementor: &'a mut RevComplementor,
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
            rev_complementor,
        }
    }

    pub fn modify_cigar(&mut self) -> Result<(), Error> {
        let mut flag = true;

        let mut cigar_vec = VecDeque::from_iter(self.cigar_str.0.iter().copied());
        let mut cigar_pos = self.cigar_str.pos() as u32;

        // if CIGAR starts with deletion cut it off
        if let Some(Cigar::Del(l)) = cigar_vec.front() {
            cigar_pos += *l;
            cigar_vec.pop_front().unwrap();
        }

        // if CIGAR ends with deletion cut it off
        if let Some(Cigar::Del(_l)) = cigar_vec.back() {
            cigar_vec.pop_back().unwrap();
        }

        // replace insertion at the beginning and end with soft clipping
        if let Some(Cigar::Ins(l)) = cigar_vec.front() {
            *cigar_vec.front_mut().unwrap() = Cigar::SoftClip(*l);
        }

        if let Some(Cigar::Ins(l)) = cigar_vec.back() {
            *cigar_vec.back_mut().unwrap() = Cigar::SoftClip(*l);
        }

        // if the length Soft clip of 5' end greater or equal to 10:
        if let Some(Cigar::SoftClip(l)) = cigar_vec
            .front()
            .copied()
            .and_then(|c| if c.len() >= 10 { Some(c) } else { None })
        {
            let l = l as usize;
            if !instance().conf.chimeric_filter && l >= Configuration::SEED_2 as usize {
                let softclip_seq = self.query_sequence.get_or_err(0..l)?;
                let seq = self.rev_complementor.reverse_complement(softclip_seq);
                let rc_seed = seq.get_or_err(0..Configuration::SEED_2 as usize)?;

                if let Some(poss) = self.ref_data.seed.get(rc_seed) {
                    if poss.len() == 1
                        && ((cigar_pos - poss.get(0).copied().unwrap() as u32) as usize)
                            < 2 * self.max_read_length
                    {
                        cigar_vec.pop_front().unwrap();
                        self.query_sequence = self.query_sequence.get_or_err(0..l)?;
                        self.query_quality = self.query_quality.get_or_err(0..l)?;

                        event!(
                            Level::INFO,
                            "{} at 5' is a chimeric at {} by SEED {}",
                            self.query_sequence.try_as_str()?,
                            cigar_pos,
                            Configuration::SEED_2,
                        )
                    }
                }
            }
        } else if let Some(Cigar::SoftClip(l)) = cigar_vec
            .back()
            .copied()
            .and_then(|c| if c.len() >= 10 { Some(c) } else { None })
        {
            let l = l as usize;
            if !instance().conf.chimeric_filter && l >= Configuration::SEED_2 as usize {
                let softclip_seq = self
                    .query_sequence
                    .get_with_int(-(l as i32)..)?;

                let seq = self.rev_complementor.reverse_complement(softclip_seq);
                let rc_seed = seq.get_with_int(-(Configuration::SEED_2)..)?;

            }
        }

        todo!()
    }
}
