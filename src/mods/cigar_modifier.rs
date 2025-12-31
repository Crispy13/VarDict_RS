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
    variants::var_utils::{is_has_and_equals, is_has_and_not_equals},
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

            if let Some(si_and_c_lens) = find_m_dmi_mdm(&cigar_vec) {
                flag =
                    self.two_dels_ins_to_complex(cigar_pos, &mut cigar_vec, si_and_c_lens, flag)?;
            } else if let Some(si_and_c_lens) = find_m_dm_dm_dm(&cigar_vec) {
                flag = self.three_deletions(cigar_pos, &mut cigar_vec, si_and_c_lens, flag)?;
            } else if let Some(si_and_cigars) = find_m_id_m_id_m_id_m(&cigar_vec) {
                flag = self.three_indels(cigar_pos, &mut cigar_vec, si_and_cigars, flag)?;
            }

            if let Some(si_and_cigars) = find_d_m_di_i(&cigar_vec) {}
        }

        todo!()
    }

    fn two_dels_ins_to_complex(
        &self,
        cigar_pos: u32,
        cigar_vd: &mut VecDeque<Cigar>,
        si_and_c_lens: (usize, [u32; 7]), // Cigars: M D M I M D M
        mut flag: bool,
    ) -> Result<bool, Error> {
        // length of both matched sequences and insertion
        let (si, c_lens) = si_and_c_lens;
        let mut tslen = (c_lens[2] + c_lens[3] + c_lens[4]) as i32;

        // length of deletions and internal matched sequences
        let mut dlen = (c_lens[1] + c_lens[2] + c_lens[4] + c_lens[5]) as i32;

        // length of internal matched sequences
        let mid = (c_lens[2] + c_lens[4]) as i32;

        // offset of first deletion in the reference sequence
        let mut refoff = (cigar_pos + c_lens[0]) as i32;

        // offset of first deletion in the read
        let mut rdoff = c_lens[0] as i32;

        // offset of first deletion in the read corrected by possibly matching bases
        let mut rdoff_corrected = c_lens[0] as i32;

        let mut rm = c_lens[6] as i32;

        if si > 0 {
            // If the complex is not at start of CIGAR string
            for cigar in (0..si).map(|i| cigar_vd[i]) {
                match cigar {
                    Cigar::Match(l) => {
                        refoff += l as i32;
                        rdoff += l as i32;
                    }
                    Cigar::RefSkip(l) | Cigar::Del(l) => {
                        refoff += l as i32;
                    }
                    Cigar::SoftClip(l) | Cigar::Ins(l) => {
                        rdoff += l as i32;
                    }
                    _ => {}
                }
            }
        }

        //number of bases after refoff/rdoff that match in reference and read
        let mut rn = 0;
        while rdoff + rn < self.query_sequence.len() as i32
            && is_has_and_equals(
                self.query_sequence
                    .get_or_err((rdoff + rn) as usize)
                    .copied()?,
                &self.ref_data.ref_seq,
                (refoff + rn) as usize,
            )
        {
            rn += 1;
        }

        rdoff_corrected += rn;
        dlen -= rn;
        tslen -= rn;

        if mid <= 15 {
            //If length of internal matched sequences is no more than 15, replace M-D-M-I-M-D complex with M-D-I
            // cigar_vd.clear();

            // 1. Overwrite the first slot with the new extended Match
            cigar_vd[si] = Cigar::Match(rdoff_corrected as u32);

            // 2. Overwrite the next slots based on the logic
            if tslen <= 0 {
                // Result: M, D, M (3 items total)
                dlen -= tslen;
                rm += tslen;

                cigar_vd[si + 1] = Cigar::Del(dlen as u32);
                cigar_vd[si + 2] = Cigar::Match(rm as u32);

                // We used 3 slots. We originally had 7.
                // We need to remove the remaining 4 items (7 - 3 = 4).
                // Remove indices from si+3 up to si+7
                cigar_vd.drain(si + 3..si + 7);
            } else {
                // Result: M, D, I, M (4 items total)
                cigar_vd[si + 1] = Cigar::Del(dlen as u32);
                cigar_vd[si + 2] = Cigar::Ins(tslen as u32);
                cigar_vd[si + 3] = Cigar::Match(rm as u32);

                // We used 4 slots. We originally had 7.
                // We need to remove the remaining 3 items (7 - 4 = 3).
                cigar_vd.drain(si + 4..si + 7);
            }

            flag = true;
        }

        Ok(flag)
    }

    fn three_deletions(
        &self,
        cigar_pos: u32,
        cigar_vd: &mut VecDeque<Cigar>,
        si_and_c_lens: (usize, [u32; 7]), // Cigars: M D M D M D M
        mut flag: bool,
    ) -> Result<bool, Error> {
        let (si, c_lens) = si_and_c_lens;

        //length of both matched sequences and insertion
        let mut tslen = (c_lens[2] + c_lens[4]) as i32;

        //length of deletions and internal matched sequences
        let mut dlen = (c_lens[1] + c_lens[2] + c_lens[3] + c_lens[4] + c_lens[5]) as i32;

        //length of internal matched sequences
        let mid = (c_lens[2] + c_lens[4]) as i32;

        //offset of first deletion in the reference sequence
        let mut refoff = (cigar_pos + c_lens[0]) as i32;

        //offset of first deletion in the read
        let mut rdoff = (c_lens[0]) as i32;

        //offset of first deletion in the read corrected by possibly matching bases
        let mut rdoff_corrected = (c_lens[0]) as i32;

        let mut rm = c_lens[6] as i32;

        if si > 0 {
            // If the complex is not at start of CIGAR string
            for cigar in (0..si).map(|i| cigar_vd[i]) {
                match cigar {
                    Cigar::Match(l) => {
                        refoff += l as i32;
                        rdoff += l as i32;
                    }
                    Cigar::RefSkip(l) | Cigar::Del(l) => {
                        refoff += l as i32;
                    }
                    Cigar::SoftClip(l) | Cigar::Ins(l) => {
                        rdoff += l as i32;
                    }
                    _ => {}
                }
            }
        }

        //number of bases after refoff/rdoff that match in reference and read
        let mut rn = 0;
        while rdoff + rn < self.query_sequence.len() as i32
            && is_has_and_equals(
                self.query_sequence
                    .get_or_err((rdoff + rn) as usize)
                    .copied()?,
                &self.ref_data.ref_seq,
                (refoff + rn) as usize,
            )
        {
            rn += 1;
        }

        rdoff_corrected += rn;
        dlen -= rn;
        tslen -= rn;

        if mid <= 15 {
            //If length of internal matched sequences is no more than 15, replace M-D-M-I-M-D complex with M-D-I

            // 1. Overwrite the first slot with the new extended Match
            cigar_vd[si] = Cigar::Match(rdoff_corrected as u32);

            // 2. Overwrite the next slots based on the logic
            if tslen <= 0 {
                // Result: M, D, M (3 items total)
                dlen -= tslen;
                rm += tslen;

                cigar_vd[si + 1] = Cigar::Del(dlen as u32);
                cigar_vd[si + 2] = Cigar::Match(rm as u32);

                // We used 3 slots. We originally had 7.
                // We need to remove the remaining 4 items (7 - 3 = 4).
                // Remove indices from si+3 up to si+7
                cigar_vd.drain(si + 3..si + 7);
            } else {
                // Result: M, D, I, M (4 items total)
                cigar_vd[si + 1] = Cigar::Del(dlen as u32);
                cigar_vd[si + 2] = Cigar::Ins(tslen as u32);
                cigar_vd[si + 3] = Cigar::Match(rm as u32);

                // We used 4 slots. We originally had 7.
                // We need to remove the remaining 3 items (7 - 4 = 3).
                cigar_vd.drain(si + 4..si + 7);
            }

            flag = true;
        }

        Ok(flag)
    }

    fn three_indels(
        &self,
        cigar_pos: u32,
        cigar_vd: &mut VecDeque<Cigar>,
        si_and_cigars: (usize, [Cigar; 7]), // Cigars: M D M D M D M
        mut flag: bool,
    ) -> Result<bool, Error> {
        let (si, cigars) = si_and_cigars;

        let mut tslen = (cigars[2].len() + cigars[4].len()) as i32;
        let mut dlen = (cigars[2].len() + cigars[4].len()) as i32;

        match cigars[1] {
            Cigar::Ins(l) => {
                tslen += l as i32;
            }
            Cigar::Del(l) => {
                dlen += l as i32;
            }
            _ => {}
        }
        match cigars[3] {
            Cigar::Ins(l) => {
                tslen += l as i32;
            }
            Cigar::Del(l) => {
                dlen += l as i32;
            }
            _ => {}
        }
        match cigars[5] {
            Cigar::Ins(l) => {
                tslen += l as i32;
            }
            Cigar::Del(l) => {
                dlen += l as i32;
            }
            _ => {}
        }

        let mid = (cigars[2].len() + cigars[4].len()) as i32;

        let mut refoff = (cigar_pos + cigars[0].len()) as i32;
        let mut rdoff = (cigars[0].len()) as i32;
        let mut rdoff_corrected = cigars[0].len() as i32;
        let mut rm = cigars[6].len() as i32;

        if si > 0 {
            // If the complex is not at start of CIGAR string
            for cigar in (0..si).map(|i| cigar_vd[i]) {
                match cigar {
                    Cigar::Match(l) => {
                        refoff += l as i32;
                        rdoff += l as i32;
                    }
                    Cigar::RefSkip(l) | Cigar::Del(l) => {
                        refoff += l as i32;
                    }
                    Cigar::SoftClip(l) | Cigar::Ins(l) => {
                        rdoff += l as i32;
                    }
                    _ => {}
                }
            }
        }

        let mut rn = 0;
        while rdoff + rn < self.query_sequence.len() as i32
            && is_has_and_equals(
                self.query_sequence
                    .get_or_err((rdoff + rn) as usize)
                    .copied()?,
                &self.ref_data.ref_seq,
                (refoff + rn) as usize,
            )
        {
            rn += 1;
        }

        rdoff_corrected += rn;
        dlen -= rn;
        tslen -= rn;

        if mid <= 15 {
            let mut drain_off = 0;
            // 1. Overwrite the first slot with the new extended Match
            cigar_vd[si] = Cigar::Match(rdoff_corrected as u32);

            // 2. Overwrite the next slots based on the logic
            if tslen <= 0 {
                // Result: M, D, M (3 items total)
                dlen -= tslen;
                rm += tslen;

                if dlen == 0 {
                    rdoff_corrected = rdoff_corrected + rm;
                    cigar_vd[si] = Cigar::Match(rdoff_corrected as u32);

                    drain_off = 1;
                } else if dlen < 0 {
                    tslen = -dlen;
                    rm += dlen;

                    if rm < 0 {
                        rdoff_corrected = rdoff_corrected + rm;
                        cigar_vd[si] = Cigar::Match(rdoff_corrected as u32);
                        cigar_vd[si + 1] = Cigar::Ins(tslen as u32);

                        drain_off = 2;
                    } else {
                        cigar_vd[si + 1] = Cigar::Ins(tslen as u32);
                        cigar_vd[si + 2] = Cigar::Match(rm as u32);

                        drain_off = 3;
                    }
                } else {
                    cigar_vd[si + 1] = Cigar::Del(dlen as u32);
                    cigar_vd[si + 2] = Cigar::Match(rm as u32);

                    drain_off = 3;
                }
            } else {
                if dlen == 0 {
                    cigar_vd[si + 1] = Cigar::Ins(tslen as u32);
                    cigar_vd[si + 2] = Cigar::Match(rm as u32);

                    drain_off = 3;
                } else if dlen < 0 {
                    rm += dlen;
                    cigar_vd[si + 1] = Cigar::Ins(tslen as u32);
                    cigar_vd[si + 2] = Cigar::Match(rm as u32);

                    drain_off = 3;
                } else {
                    cigar_vd[si + 1] = Cigar::Del(dlen as u32);
                    cigar_vd[si + 2] = Cigar::Ins(tslen as u32);
                    cigar_vd[si + 3] = Cigar::Match(rm as u32);

                    drain_off = 4;
                }
            }

            cigar_vd.drain(si + drain_off..si + 7);

            flag = true;
        }

        Ok(flag)
    }

    fn combine_to_close_to_correct(
        cigar_vd: &mut VecDeque<Cigar>,
        si_and_cigars: (usize, [Cigar; 3], Option<Cigar>),
        mut flag: bool,
    ) -> Result<bool, Error> {
        let (si, cigars, opt_i) = si_and_cigars;
        let g1 = cigars[0].len() as i32; // D
        let g2 = cigars[1].len() as i32; // M
        let g3 = cigars[2].len() as i32; // DI

        if g2 <= 15 {
            //length of both deletions and matched sequence
            let mut dlen = g1 + g2;

            //matched sequence length
            let mut ilen = g2;

            let mut drain_end_off = 4;
            match opt_i {
                Some(Cigar::Ins(l)) => {
                    ilen += g3;
                }
                Some(Cigar::Del(l)) => {
                    dlen += g3;
                    ilen += l as i32;
                }
                Some(c) => {
                    unreachable!()
                }
                None => {
                    drain_end_off = 3;
                }
            }

            cigar_vd[si] = Cigar::Del(dlen as u32);
            cigar_vd[si + 1] = Cigar::Ins(ilen as u32);

            cigar_vd.drain(si..si + drain_end_off);

            flag = true;
        }

        Ok(flag)
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

// Returns the index 'i' where the pattern starts
fn find_m_id_m_id_m_id_m(cigar: &VecDeque<Cigar>) -> Option<(usize, [Cigar; 7])> {
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
                &c1 @ Cigar::Match(i1),                  // 1
                &c2 @ (Cigar::Ins(i2) | Cigar::Del(i2)), // 2
                &c3 @ Cigar::Match(i3),                  // 3
                &c4 @ (Cigar::Ins(i4) | Cigar::Del(i4)), // 4
                &c5 @ Cigar::Match(i5),                  // 5
                &c6 @ (Cigar::Ins(i6) | Cigar::Del(i6)), // 4
                &c7 @ Cigar::Match(i7),                  // 7
            ) => return Some((i, [c1, c2, c3, c4, c5, c6, c7])), // FOUND IT!
            _ => continue,
        }
    }
    None
}

fn find_m_dm_dm_dm(cigar: &VecDeque<Cigar>) -> Option<(usize, [u32; 7])> {
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
                &Cigar::Del(i4),   // 4
                &Cigar::Match(i5), // 5
                &Cigar::Del(i6),   // 4
                &Cigar::Match(i7), // 7
            ) => return Some((i, [i1, i2, i3, i4, i5, i6, i7])), // FOUND IT!
            _ => continue,
        }
    }
    None
}

fn find_d_m_di_i(cigar: &VecDeque<Cigar>) -> Option<(usize, [Cigar; 3], Option<Cigar>)> {
    // We need at least 7 elements
    if cigar.len() < 3 {
        return None;
    }

    // Loop through valid start positions
    for i in 0..cigar.len() - 2 {
        // Use pattern matching on references
        match (&cigar[i], &cigar[i + 1], &cigar[i + 2]) {
            (
                &c1 @ Cigar::Del(_),                   // 1
                &c2 @ Cigar::Match(_),                 // 2
                &c3 @ (Cigar::Ins(_) | Cigar::Del(_)), // 3
            ) => {
                return Some((
                    i,
                    [c1, c2, c3],
                    cigar.get(i + 3).copied().and_then(|c| {
                        if matches!(c, Cigar::Ins(_)) {
                            Some(c)
                        } else {
                            None
                        }
                    }),
                ));
            }
            _ => continue,
        }
    }

    None
}
#[cfg(test)]
mod tests {
    use super::*;
}
