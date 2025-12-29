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
    variants::var_utils::is_has_and_not_equals,
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
                        && ((cigar_pos as i32 - poss.get(0).copied().unwrap() as i32).abs()
                            as usize)
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
                let softclip_seq = self.query_sequence.get_with_int(-(l as i32)..)?;

                let seq = self.rev_complementor.reverse_complement(softclip_seq);
                let rc_seed = seq.get_with_int(-(Configuration::SEED_2)..)?;

                if let Some(poss) = self.ref_data.seed.get(rc_seed) {
                    if poss.len() == 1
                        && ((cigar_pos as i32 - poss.get(0).copied().unwrap() as i32).abs()
                            as usize)
                            < 2 * self.max_read_length
                    {
                        cigar_vec.pop_back().unwrap();
                        self.query_sequence = self.query_sequence.get_with_int(0..-(l as i32))?;
                        self.query_quality = self.query_quality.get_with_int(0..-(l as i32))?;

                        event!(
                            Level::INFO,
                            "{} at 3' is a chimeric at {} by SEED {}",
                            self.query_sequence.try_as_str()?,
                            cigar_pos,
                            Configuration::SEED_2,
                        )
                    }
                }
            }
        }

        while flag && self.indel > 0 {
            flag = false;

            // check read starts with softclip then insertion or deletion.
            match (cigar_vec.get(0), cigar_vec.get(1)) {
                (Some(&Cigar::SoftClip(sl)), Some(c2 @ (&Cigar::Ins(idl) | &Cigar::Del(idl)))) => {
                    let tslen = sl + if matches!(c2, Cigar::Ins(_)) { idl } else { 0 };
                    cigar_pos += if matches!(c2, Cigar::Del(_)) { idl } else { 0 };

                    cigar_vec.pop_front().unwrap();
                    *cigar_vec.front_mut().unwrap() = Cigar::SoftClip(tslen);

                    flag = true;
                }
                _ => {}
            }

            let mut cigar_vec_iter = cigar_vec.iter().rev();

            match (cigar_vec_iter.next(), cigar_vec_iter.next()) {
                (Some(&Cigar::SoftClip(sl)), Some(c2 @ (&Cigar::Ins(idl) | &Cigar::Del(idl)))) => {
                    let tslen = sl + if matches!(c2, Cigar::Ins(_)) { idl } else { 0 };

                    cigar_vec.pop_back().unwrap();
                    *cigar_vec.back_mut().unwrap() = Cigar::SoftClip(tslen);

                    flag = true;
                }
                _ => {}
            }

            match (cigar_vec.get(0), cigar_vec.get(1), cigar_vec.get(2)) {
                (
                    Some(&Cigar::SoftClip(sl)),
                    Some(&Cigar::Match(ml)),
                    Some(c3 @ (&Cigar::Ins(idl) | &Cigar::Del(idl))),
                ) => {
                    if ml <= 10 {
                        let tslen = sl + ml + if matches!(c3, Cigar::Ins(_)) { idl } else { 0 };
                        cigar_pos += ml + if matches!(c3, Cigar::Del(_)) { idl } else { 0 };

                        cigar_vec.drain(..2);
                        *cigar_vec.front_mut().unwrap() = Cigar::SoftClip(tslen);

                        flag = true;
                    }
                }
                _ => {}
            }

            let mut cigar_vec_iter = cigar_vec.iter().rev();
            match (
                cigar_vec_iter.next(),
                cigar_vec_iter.next(),
                cigar_vec_iter.next(),
            ) {
                (
                    Some(&Cigar::SoftClip(sl)),
                    Some(&Cigar::Match(ml)),
                    Some(c3 @ (&Cigar::Ins(idl) | &Cigar::Del(idl))),
                ) => {
                    if ml <= 10 {
                        let tslen = sl + ml + if matches!(c3, Cigar::Ins(_)) { idl } else { 0 };
                        cigar_pos += ml + if matches!(c3, Cigar::Del(_)) { idl } else { 0 };

                        cigar_vec.drain(..2);
                        *cigar_vec.back_mut().unwrap() = Cigar::SoftClip(tslen);

                        flag = true;
                    }
                }
                _ => {}
            }

            match (cigar_vec.get(0), cigar_vec.get(1), cigar_vec.get(2)) {
                (
                    Some(&Cigar::Match(ml1)),
                    Some(c_id @ (&Cigar::Ins(idl) | &Cigar::Del(idl))),
                    Some(&Cigar::Match(mut ml2)),
                ) => {
                    let mut tslen = ml1
                        + if matches!(c_id, Cigar::Ins(_)) {
                            idl
                        } else {
                            0
                        };
                    cigar_pos += ml1
                        + if matches!(c_id, Cigar::Del(_)) {
                            idl
                        } else {
                            0
                        };

                    let mut tn = 0;
                    while tn < ml2
                        && is_has_and_not_equals(
                            self.query_sequence
                                .get_or_err((tslen + tn) as usize)
                                .copied()?,
                            &self.ref_data.ref_seq,
                            (cigar_pos + tn) as usize,
                        )
                    {
                        tn += 1;
                    }

                    tslen += tn;
                    ml2 -= tn;
                    cigar_pos += tn;

                    cigar_vec.pop_front().unwrap();
                    *cigar_vec.get_mut(0).unwrap() = Cigar::SoftClip(tslen);
                    *cigar_vec.get_mut(1).unwrap() = Cigar::Match(ml2);

                    flag = true;
                }
                _ => {}
            }

            let mut cigar_vec_iter = cigar_vec.iter().rev();
            match (cigar_vec_iter.next(), cigar_vec_iter.next()) {
                (Some(&Cigar::Match(ml)), Some(c_id @ (&Cigar::Ins(idl) | &Cigar::Del(idl)))) => {
                    let tslen = ml
                        + if matches!(c_id, Cigar::Ins(_)) {
                            idl
                        } else {
                            0
                        };

                    *cigar_vec.back_mut().unwrap() = Cigar::SoftClip(tslen);

                    flag = true;
                }
                _ => {}
            }

            match cigar_vec {
                _ => {}
            }
        }

        todo!()
    }
}

// Returns the index 'i' where the pattern starts
fn find_m_dmi_mdm(cigar: &VecDeque<Cigar>) -> Option<(usize, [u32; 7])> {
    // We need at least 7 elements
    if cigar.len() < 7 {
        return None;
    }

    // Loop through valid start positions
    for i in 0..cigar.len() - 6 {
        // Use pattern matching on references
        match (
            &cigar[i],
            &cigar[i + 1],
            &cigar[i + 2],
            &cigar[i + 3],
            &cigar[i + 4],
            &cigar[i + 5],
            &cigar[i + 6],
        ) {
            (
                &Cigar::Match(i1), // 1
                &Cigar::Del(i2),   // 2
                &Cigar::Match(i3), // 3
                &Cigar::Ins(i4),   // 4
                &Cigar::Match(i5), // 5
                &Cigar::Del(i6),   // 6
                &Cigar::Match(i7), // 7
            ) => return Some((i, [i1, i2, i3, i4, i5, i6, i7])), // FOUND IT!
            _ => continue,
        }
    }
    None
}
