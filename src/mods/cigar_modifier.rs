use std::{
    borrow::Cow,
    collections::{HashMap, HashSet, VecDeque},
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
    data::{ModifiedCigar, reference::Reference, region::Region},
    scopedata::global_read_only_scope::{INSTANCE, instance},
    utils::{BytesExt, SliceExt, SliceExt2},
    variants::var_utils::{
        HomoPolymerChecker, is_has_and_equals, is_has_and_equals_ref_and_seq_base,
        is_has_and_not_equals, is_has_and_not_equals_ref_and_seq_base,
    },
};

pub struct CigarModifier<'a, 'b> {
    pos: i64,
    cigar_str: Cow<'a, CigarStringView>,
    original_cigar: &'a CigarStringView,
    query_sequence: &'b [u8],
    query_quality: &'b [u8],
    ref_data: &'a Reference,
    indel: u32,
    max_read_length: usize,
    region: &'a Region,
    rev_complementor: &'a mut RevComplementor,
}

impl<'a, 'b> CigarModifier<'a, 'b> {
    pub(crate) fn new(
        pos: i64,
        cigar_str: &'a CigarStringView,
        query_sequence: &'b [u8],
        query_quality: &'b [u8],
        ref_data: &'a Reference,
        indel: u32,
        max_read_length: usize,
        region: &'a Region,
        rev_complementor: &'a mut RevComplementor,
    ) -> Self {
        Self {
            pos,
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

    fn contig_ref_seq(&self) -> &Vec<u8> {
        &self.ref_data.ref_seq
    }

    pub fn modify_cigar(&mut self) -> Result<ModifiedCigar<'b>, Error> {
        let mut flag = true;

        let mut cigar_vec = VecDeque::from_iter(self.cigar_str.0.iter().copied());
        let mut ref_start_pos = self.cigar_str.pos() as u32;

        // if CIGAR starts with deletion cut it off
        if let Some(Cigar::Del(l)) = cigar_vec.front() {
            ref_start_pos += *l;
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
                        && ((ref_start_pos as i32 - poss.get(0).copied().unwrap() as i32).abs()
                            as usize)
                            < 2 * self.max_read_length
                    {
                        cigar_vec.pop_front().unwrap();
                        self.query_sequence = self.query_sequence.get_or_err(l..)?;
                        self.query_quality = self.query_quality.get_or_err(l..)?;

                        event!(
                            Level::INFO,
                            "{} at 5' is a chimeric at {} by SEED {}",
                            self.query_sequence.try_as_str()?,
                            ref_start_pos,
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
                        && ((ref_start_pos as i32 - poss.get(0).copied().unwrap() as i32).abs()
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
                            ref_start_pos,
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
                    ref_start_pos += if matches!(c2, Cigar::Del(_)) { idl } else { 0 };

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
                        ref_start_pos += ml + if matches!(c3, Cigar::Del(_)) { idl } else { 0 };

                        cigar_vec.drain(..2);
                        *cigar_vec.front_mut().unwrap() = Cigar::SoftClip(tslen);

                        flag = true;
                    }
                }
                _ => {}
            }

            // NUMBER_IorD_NUMBER_M_NUMBER_S_END
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
                        ref_start_pos += ml + if matches!(c3, Cigar::Del(_)) { idl } else { 0 };

                        cigar_vec.drain(cigar_vec.len() - 2..);
                        *cigar_vec.back_mut().unwrap() = Cigar::SoftClip(tslen);

                        flag = true;
                    }
                }
                _ => {}
            }

            // BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M
            match (cigar_vec.get(0), cigar_vec.get(1), cigar_vec.get(2)) {
                (
                    Some(&Cigar::Match(ml1)),
                    Some(c_id @ (&Cigar::Ins(idl) | &Cigar::Del(idl))),
                    Some(&Cigar::Match(mut ml2)),
                ) if ml1 < 10 => {
                    let mut tslen = ml1
                        + if matches!(c_id, Cigar::Ins(_)) {
                            idl
                        } else {
                            0
                        };
                    ref_start_pos += ml1
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
                            (ref_start_pos + tn) as usize,
                        )
                    {
                        tn += 1;
                    }

                    tslen += tn;
                    ml2 -= tn;
                    ref_start_pos += tn;

                    cigar_vec.pop_front().unwrap();
                    *cigar_vec.get_mut(0).unwrap() = Cigar::SoftClip(tslen);
                    *cigar_vec.get_mut(1).unwrap() = Cigar::Match(ml2);

                    flag = true;
                }
                _ => {}
            }

            // NUMBER_IorD_DIGIT_M_END
            let mut cigar_vec_iter = cigar_vec.iter().rev();
            match (cigar_vec_iter.next(), cigar_vec_iter.next()) {
                (Some(&Cigar::Match(ml)), Some(c_id @ (&Cigar::Ins(idl) | &Cigar::Del(idl))))
                    if ml < 10 =>
                {
                    let tslen = ml
                        + if matches!(c_id, Cigar::Ins(_)) {
                            idl
                        } else {
                            0
                        };

                    cigar_vec.pop_back().unwrap();
                    *cigar_vec.back_mut().unwrap() = Cigar::SoftClip(tslen);

                    flag = true;
                }
                _ => {}
            }

            if let Some(si_and_c_lens) = find_m_dmi_mdm(&cigar_vec) {
                flag = self.two_dels_ins_to_complex(
                    ref_start_pos,
                    &mut cigar_vec,
                    si_and_c_lens,
                    flag,
                )?;
            } else if let Some(si_and_c_lens) = find_m_dm_dm_dm(&cigar_vec) {
                flag = self.three_deletions(ref_start_pos, &mut cigar_vec, si_and_c_lens, flag)?;
            } else if let Some(si_and_cigars) = find_m_id_m_id_m_id_m(&cigar_vec) {
                flag = self.three_indels(ref_start_pos, &mut cigar_vec, si_and_cigars, flag)?;
            }

            if let Some(si_and_cigars) = find_d_m_di_i(&cigar_vec) {
                flag = self.combine_to_close_to_correct(&mut cigar_vec, si_and_cigars, flag)?;
            }

            if let Some(si_and_cigars) = find_d_i_m_id_i(&cigar_vec) {
                flag = self.combine_to_close_to_one(&mut cigar_vec, si_and_cigars, flag)?;
            }

            if let Some(si_and_cigars) = find_d_d(&cigar_vec) {
                let (si, cigars) = si_and_cigars;
                let dlen = (cigars[0].len() + cigars[1].len()) as i32;
                cigar_vec[si] = Cigar::Del(dlen as u32);
                cigar_vec.remove(si + 1).unwrap();

                flag = true;
            }

            if let Some(si_and_cigars) = find_i_i(&cigar_vec) {
                let (si, cigars) = si_and_cigars;
                let ilen = (cigars[0].len() + cigars[1].len()) as i32;
                cigar_vec[si] = Cigar::Ins(ilen as u32);
                cigar_vec.remove(si + 1).unwrap();

                flag = true;
            }
        }

        let mut cigar_iter_rev = cigar_vec.iter().rev();
        match (cigar_iter_rev.next(), cigar_iter_rev.next()) {
            (Some(&Cigar::SoftClip(sl)), Some(&Cigar::Match(ml))) => {
                self.capture_mis_softly3_ms(ref_start_pos, &mut cigar_vec, sl, ml)?;
            }
            (Some(&Cigar::Match(ml)), _) => {
                self.capture_mis_softly3_mismatches(ref_start_pos, &mut cigar_vec, ml)?;
            }
            _ => {}
        }

        match (cigar_vec.get(0), cigar_vec.get(1)) {
            (Some(&Cigar::SoftClip(sl)), Some(&Cigar::Match(ml))) => {
                self.combine_dig_s_dig_m(&mut ref_start_pos, &mut cigar_vec, sl, ml)?;
            }
            (Some(&Cigar::Match(ml)), _) => {
                self.combine_begin_dig_m(&mut ref_start_pos, &mut cigar_vec, ml)?;
            }
            _ => {}
        }

        let mc = ModifiedCigar::new(
            ref_start_pos as i64,
            cigar_vec,
            self.query_sequence,
            self.query_quality,
        );

        Ok(mc)
    }

    fn capture_mis_softly3_ms(
        &self,
        align_start_pos: u32,
        cigar_vd: &mut VecDeque<Cigar>,
        sl: u32,
        ml: u32,
    ) -> Result<(), Error> {
        // capture_mis_softly_ms
        let mut mch = ml as i32;
        let mut soft = sl as i32;

        // offset of soft-clipped sequence in the reference string (position + length of
        // matched)
        let mut refoff = (align_start_pos + ml) as i32;

        // offset of soft-clipped sequence in the read
        let mut rdoff = ml as i32;

        for &cigar in cigar_vd.iter().take(cigar_vd.len().saturating_sub(2)) {
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

        // number of bases after refoff/rdoff that match in reference and read sequences
        let mut rn = 0;

        let mut homopolymer_checker = HomoPolymerChecker::new();

        while rn < soft
            && is_has_and_equals_ref_and_seq_base(
                &self.ref_data.ref_seq,
                (refoff + rn) as usize,
                self.query_sequence,
                (rdoff + rn) as usize,
            )
            && self.query_quality.get_or_err((rdoff + rn) as usize)? - 33
                > Configuration::LOW_QUAL as u8
        {
            rn += 1;
        }

        if rn > 0 {
            mch += rn;
            soft -= rn;

            if soft > 0 {
                *cigar_vd.get_mut(cigar_vd.len() - 1).unwrap() = Cigar::SoftClip(soft as u32);
                *cigar_vd.get_mut(cigar_vd.len() - 2).unwrap() = Cigar::Match(mch as u32);
            } else {
                cigar_vd.pop_back().unwrap();
                *cigar_vd.get_mut(cigar_vd.len() - 1).unwrap() = Cigar::Match(mch as u32);
            }
            rn = 0;
        }

        if soft > 0 {
            while rn + 1 < soft
                && is_has_and_equals_ref_and_seq_base(
                    &self.ref_data.ref_seq,
                    (refoff + rn + 1) as usize,
                    self.query_sequence,
                    (rdoff + rn + 1) as usize,
                )
                && self.query_quality.get_or_err((rdoff + rn + 1) as usize)? - 33
                    > Configuration::LOW_QUAL as u8
            {
                rn += 1;
                homopolymer_checker.record_base(
                    self.ref_data
                        .ref_seq
                        .get_or_err((refoff + rn + 1) as usize)
                        .copied()?,
                );
            }

            if rn > 4 && !homopolymer_checker.is_homopolymer() {
                mch += rn + 1;
                soft -= rn + 1;

                let cigar_vd_len = cigar_vd.len();
                if soft > 0 {
                    cigar_vd[cigar_vd_len - 1] = Cigar::SoftClip(soft as u32);
                    cigar_vd[cigar_vd_len - 2] = Cigar::Match(mch as u32);
                } else {
                    cigar_vd.pop_back().unwrap();
                    cigar_vd[cigar_vd_len - 1] = Cigar::Match(mch as u32);
                }
            }

            if rn == 0 {
                let mut rrn = 0;
                let mut rmch = 0;

                while rrn < mch && rn < mch {
                    if self
                        .ref_data
                        .ref_seq
                        .get((refoff - rrn - 1) as usize)
                        .is_none()
                    {
                        break;
                    }

                    if rrn < rdoff
                        && is_has_and_not_equals_ref_and_seq_base(
                            &self.ref_data.ref_seq,
                            (refoff - rrn - 1) as usize,
                            self.query_sequence,
                            (rdoff - rrn - 1) as usize,
                        )
                    {
                        rn = rrn + 1;
                        rmch = 0;
                    } else if rrn < rdoff
                        && is_has_and_equals_ref_and_seq_base(
                            &self.ref_data.ref_seq,
                            (refoff - rrn - 1) as usize,
                            self.query_sequence,
                            (rdoff - rrn - 1) as usize,
                        )
                    {
                        rmch += 1;
                    }

                    rrn += 1;

                    // Stop at three consecure matches
                    if rmch >= 3 {
                        break;
                    }
                }

                if rn > 0 && rn < mch {
                    soft += rn;
                    mch -= rn;
                    *cigar_vd.get_mut(cigar_vd.len() - 1).unwrap() = Cigar::SoftClip(soft as u32);
                    *cigar_vd.get_mut(cigar_vd.len() - 1).unwrap() = Cigar::Match(mch as u32);
                }
            }
        }

        Ok(())
    }

    fn capture_mis_softly3_mismatches(
        &self,
        align_start_pos: u32,
        cigar_vd: &mut VecDeque<Cigar>,
        ml: u32,
    ) -> Result<(), Error> {
        // capture_mis_softly3_mismatches
        let mut mch = ml as i32;
        let mut refoff = (align_start_pos + ml) as i32;
        let mut rdoff = mch;

        // If the complex is not at start of CIGAR string
        for &cigar in cigar_vd.iter().take(cigar_vd.len().saturating_sub(1)) {
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

        let mut rn = 0;
        let mut rrn = 0;
        let mut rmch = 0;

        while rrn < mch && rn < mch {
            if self
                .contig_ref_seq()
                .get((refoff - rrn - 1) as usize)
                .is_none()
            {
                break;
            }

            if rrn < rdoff
                && is_has_and_not_equals_ref_and_seq_base(
                    self.contig_ref_seq(),
                    (refoff - rrn - 1) as usize,
                    self.query_sequence,
                    (rdoff - rrn - 1) as usize,
                )
            {
                rn = rrn + 1;
                rmch = 0;
            } else if rrn < rdoff
                && is_has_and_equals_ref_and_seq_base(
                    self.contig_ref_seq(),
                    (refoff - rrn - 1) as usize,
                    self.query_sequence,
                    (rdoff - rrn - 1) as usize,
                )
            {
                rmch += 1;
            }

            rrn += 1;

            // Stop at three consecure matches
            if rmch >= 3 {
                break;
            }
        }

        mch -= rn;
        if rn > 0 && rn <= 3 {
            *cigar_vd.get_mut(cigar_vd.len() - 1).unwrap() = Cigar::Match(mch as u32);
            cigar_vd.push_back(Cigar::SoftClip(rn as u32));
        }

        Ok(())
    }

    fn two_dels_ins_to_complex(
        &self,
        align_start_pos: u32,
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
        let mut refoff = (align_start_pos + c_lens[0]) as i32;

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
        align_start_pos: u32,
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
        let mut refoff = (align_start_pos + c_lens[0]) as i32;

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
        align_start_pos: u32,
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

        let mut refoff = (align_start_pos + cigars[0].len()) as i32;
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
        &self,
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

            let mut drain_end_off = 3;
            match cigars[2] {
                Cigar::Ins(l) => {
                    ilen += g3;
                }
                Cigar::Del(l) => {
                    dlen += g3;

                    if let Some(c) = opt_i {
                        ilen += c.len() as i32;
                        drain_end_off += 1;
                    }
                }
                oth => {
                    unreachable!("{oth:?}")
                }
            }

            cigar_vd[si] = Cigar::Del(dlen as u32);
            cigar_vd[si + 1] = Cigar::Ins(ilen as u32);

            cigar_vd.drain(si + 2..si + drain_end_off);

            flag = true;
        }

        Ok(flag)
    }

    fn combine_to_close_to_one(
        &self,
        cigar_vd: &mut VecDeque<Cigar>,
        si_and_cigars: (usize, [Cigar; 3], Option<Cigar>),
        mut flag: bool,
    ) -> Result<bool, Error> {
        let (si, cigars, opt_i) = si_and_cigars;

        let g2 = cigars[0].len() as i32;
        let g3 = cigars[1].len() as i32;
        let g4 = cigars[2].len() as i32;

        if g3 <= 15 {
            let mut drain_end_off = 3;

            //length of matched sequence and deletion
            let mut dlen = g3;
            //length of first insertion and matched sequence
            let mut ilen = g2 + g3;

            match cigars[2] {
                Cigar::Ins(_) => {
                    ilen += g4;
                }
                Cigar::Del(_) => {
                    dlen += g4;

                    //last insertion string
                    if let Some(c) = opt_i {
                        ilen += c.len() as i32;
                        drain_end_off += 1;
                    }
                }
                _ => {}
            }

            cigar_vd[si] = Cigar::Del(dlen as u32);
            cigar_vd[si + 1] = Cigar::Ins(ilen as u32);

            cigar_vd.drain(si + 2..si + drain_end_off);

            flag = true;
        }

        Ok(flag)
    }

    fn combine_dig_s_dig_m(
        &self,
        cigar_pos: &mut u32,
        cigar_vd: &mut VecDeque<Cigar>,
        sl: u32,
        ml: u32,
    ) -> Result<(), Error> {
        //length of matched sequence
        let mut mch = ml as i32;
        //length of soft-clipping
        let mut soft = sl as i32;

        //number of bases before matched sequence that match in reference and read sequences
        let mut rn = 0;
        let mut homop = HomoPolymerChecker::new();

        while rn < soft
            && is_has_and_equals_ref_and_seq_base(
                self.contig_ref_seq(),
                ((*cigar_pos as i32) - rn - 1) as usize,
                self.query_sequence,
                (soft - rn - 1) as usize,
            )
            && self.query_quality.get_or_err((soft - rn - 1) as usize)? - 33
                > Configuration::LOW_QUAL as u8
        {
            rn += 1;
        }

        if rn > 0 {
            mch += rn;
            soft -= rn;

            if soft > 0 {
                cigar_vd[0] = Cigar::SoftClip(soft as u32);
                cigar_vd[1] = Cigar::Match(mch as u32);
            } else {
                cigar_vd.pop_front().unwrap();
                cigar_vd[0] = Cigar::Match(mch as u32);
            }

            *cigar_pos -= rn as u32;
            rn = 0;
        }

        if soft > 0 {
            while rn + 1 < soft
                && is_has_and_equals_ref_and_seq_base(
                    self.contig_ref_seq(),
                    ((*cigar_pos as i32) - rn - 2) as usize,
                    self.query_sequence,
                    (soft - rn - 2) as usize,
                )
                && self.query_quality.get_or_err((soft - rn - 2) as usize)? - 33
                    > Configuration::LOW_QUAL as u8
            {
                rn += 1;
                homop.record_base(
                    self.contig_ref_seq()
                        .get_or_err(((*cigar_pos as i32) - rn - 2) as usize)
                        .copied()?,
                );
            }

            if (rn > 4 && !homop.is_homopolymer())
                || is_has_and_equals_ref_and_seq_base(
                    self.contig_ref_seq(),
                    ((*cigar_pos as i32) - 1) as usize,
                    self.query_sequence,
                    (soft - 1) as usize,
                )
            {
                mch += rn + 1;
                soft -= rn + 1;

                if soft > 0 {
                    cigar_vd[0] = Cigar::SoftClip(soft as u32);
                    cigar_vd[1] = Cigar::Match(mch as u32);
                } else {
                    cigar_vd.pop_back().unwrap();
                    cigar_vd[0] = Cigar::Match(mch as u32);
                }

                *cigar_pos -= rn as u32 + 1;
            }

            if rn == 0 {
                let mut rrn = 0;
                let mut rmch = 0;

                while rrn < mch && rn < mch {
                    if self
                        .contig_ref_seq()
                        .get((*cigar_pos as i32 + rrn) as usize)
                        .is_none()
                    {
                        break;
                    }

                    if is_has_and_not_equals_ref_and_seq_base(
                        self.contig_ref_seq(),
                        (*cigar_pos as i32 + rrn) as usize,
                        self.query_sequence,
                        (soft + rrn) as usize,
                    ) {
                        rn = rrn + 1;
                        rmch = 0;
                    } else if is_has_and_equals_ref_and_seq_base(
                        self.contig_ref_seq(),
                        (*cigar_pos as i32 + rrn) as usize,
                        self.query_sequence,
                        (soft + rrn) as usize,
                    ) {
                        rmch += 1;
                    }

                    rrn += 1;

                    if rmch >= 3 {
                        break;
                    }
                }

                if rn > 0 && rn < mch {
                    soft += rn;
                    mch -= rn;

                    cigar_vd[0] = Cigar::SoftClip(soft as u32);
                    cigar_vd[1] = Cigar::Match(mch as u32);

                    *cigar_pos += rn as u32;
                }
            }
        }

        Ok(())
    }

    fn combine_begin_dig_m(
        &self,
        cigar_pos: &mut u32,
        cigar_vd: &mut VecDeque<Cigar>,
        ml: u32,
    ) -> Result<(), Error> {
        let mut mch = ml as i32;
        let mut rn = 0;
        let mut rrn = 0;
        let mut rmch = 0;

        while rrn < mch && rn < mch {
            if self
                .contig_ref_seq()
                .get((*cigar_pos as i32 + rrn) as usize)
                .is_none()
            {
                break;
            }

            if is_has_and_not_equals_ref_and_seq_base(
                self.contig_ref_seq(),
                (*cigar_pos as i32 + rrn) as usize,
                self.query_sequence,
                (rrn) as usize,
            ) {
                rn = rrn + 1;
                rmch = 0;
            } else if is_has_and_equals_ref_and_seq_base(
                self.contig_ref_seq(),
                (*cigar_pos as i32 + rrn) as usize,
                self.query_sequence,
                (rrn) as usize,
            ) {
                rmch += 1;
            }

            rrn += 1;

            if rmch >= 3 {
                break;
            }
        }

        if rn > 0 && rn <= 3 {
            mch -= rn;
            cigar_vd[0] = Cigar::Match(mch as u32);
            cigar_vd.push_front(Cigar::SoftClip(rn as u32));
            *cigar_pos += rn as u32;
        }

        Ok(())
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

fn find_d_i_m_id_i(cigar: &VecDeque<Cigar>) -> Option<(usize, [Cigar; 3], Option<Cigar>)> {
    // We need at least 4 elements
    if cigar.len() < 4 {
        return None;
    }

    // Loop through valid start positions
    for i in 0..cigar.len() - 3 {
        let c1 = cigar[i];
        match c1 {
            Cigar::Del(_) | Cigar::HardClip(_) => {
                continue;
            }
            _ => {}
        }

        // Use pattern matching on references
        match (&cigar[i + 1], &cigar[i + 2], &cigar[i + 3]) {
            (
                &c2 @ Cigar::Ins(_),                   // 1
                &c3 @ Cigar::Match(_),                 // 2
                &c4 @ (Cigar::Ins(_) | Cigar::Del(_)), // 3
            ) => {
                return Some((
                    i + 1,
                    [c2, c3, c4],
                    cigar.get(i + 4).copied().and_then(|c| {
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

fn find_d_d(cigar: &VecDeque<Cigar>) -> Option<(usize, [Cigar; 2])> {
    // We need at least 2 elements
    if cigar.len() < 2 {
        return None;
    }

    // Loop through valid start positions
    for i in 0..cigar.len() - 1 {
        // Use pattern matching on references
        match (&cigar[i], &cigar[i + 1]) {
            (
                &c1 @ Cigar::Del(_), // 1
                &c2 @ Cigar::Del(_), // 2
            ) => return Some((i, [c1, c2])),
            _ => continue,
        }
    }
    None
}

fn find_i_i(cigar: &VecDeque<Cigar>) -> Option<(usize, [Cigar; 2])> {
    // We need at least 2 elements
    if cigar.len() < 2 {
        return None;
    }

    // Loop through valid start positions
    for i in 0..cigar.len() - 1 {
        // Use pattern matching on references
        match (&cigar[i], &cigar[i + 1]) {
            (
                &c1 @ Cigar::Ins(_), // 1
                &c2 @ Cigar::Ins(_), // 2
            ) => return Some((i, [c1, c2])),
            _ => continue,
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use crate::{
        mods::cigar_parser::CigarParser,
        scopedata::global_read_only_scope::{GlobalReadOnlyScope, INSTANCE},
    };
    use rust_htslib::bam::record::{Cigar, CigarString, CigarStringView};

    use super::*;

    /// Helper function to parse a CIGAR string like "15M1I1M1I2M1I27M" into CigarString
    fn parse_cigar_string(cigar_str: &str) -> CigarString {
        let mut cigars = Vec::new();
        let mut num = String::new();
        
        for c in cigar_str.chars() {
            if c.is_ascii_digit() {
                num.push(c);
            } else {
                let len: u32 = num.parse().expect("Invalid CIGAR length");
                num.clear();
                
                let cigar = match c {
                    'M' => Cigar::Match(len),
                    'I' => Cigar::Ins(len),
                    'D' => Cigar::Del(len),
                    'N' => Cigar::RefSkip(len),
                    'S' => Cigar::SoftClip(len),
                    'H' => Cigar::HardClip(len),
                    'P' => Cigar::Pad(len),
                    '=' => Cigar::Equal(len),
                    'X' => Cigar::Diff(len),
                    _ => panic!("Unknown CIGAR operator: {}", c),
                };
                cigars.push(cigar);
            }
        }
        
        CigarString(cigars)
    }

    /// Helper to convert CigarStringView back to a string
    fn cigar_to_string(cigar: &CigarStringView) -> String {
        cigar.iter().map(|c| match c {
            Cigar::Match(l) => format!("{}M", l),
            Cigar::Ins(l) => format!("{}I", l),
            Cigar::Del(l) => format!("{}D", l),
            Cigar::RefSkip(l) => format!("{}N", l),
            Cigar::SoftClip(l) => format!("{}S", l),
            Cigar::HardClip(l) => format!("{}H", l),
            Cigar::Pad(l) => format!("{}P", l),
            Cigar::Equal(l) => format!("{}=", l),
            Cigar::Diff(l) => format!("{}X", l),
        }).collect()
    }

    #[test]
    fn test_parse_cigar_string_helper() {
        let cigar = parse_cigar_string("15M1I1M1I2M1I27M");
        assert_eq!(cigar.0.len(), 7);
        assert!(matches!(cigar.0[0], Cigar::Match(15)));
        assert!(matches!(cigar.0[1], Cigar::Ins(1)));
        assert!(matches!(cigar.0[2], Cigar::Match(1)));
        assert!(matches!(cigar.0[3], Cigar::Ins(1)));
        assert!(matches!(cigar.0[4], Cigar::Match(2)));
        assert!(matches!(cigar.0[5], Cigar::Ins(1)));
        assert!(matches!(cigar.0[6], Cigar::Match(27)));
    }

    // Note: find_offset test requires complex setup of CigarParser state
    // which is tightly coupled to the parsing loop. The Java test creates
    // a standalone CigarParser and calls findOffset directly, but in Rust
    // the function is an internal method that relies on self.contig_ref_seq().
    // This test is marked as ignored until the architecture allows easier testing.
    #[test]
    #[ignore = "requires refactoring to make find_offset testable in isolation"]
    #[should_panic(expected = "Requires architecture changes to test find_offset in isolation")]
    fn find_offset() {
        let conf = Configuration {
            goodq: 23.0,
            vext: 3,
            ..Default::default()
        };

        INSTANCE.get_or_init(|| GlobalReadOnlyScope {
            conf,
            ..Default::default()
        });

        let _ref_pos = 1;
        let _read_pos = 2;
        let _cigar_len = 3;
        let _query_sequence = "ACGTACGT";
        let _query_quality = "<<<<<<<<";
        let _ref_seq = "AA";

        // Java test expects: Offset(2, "GT", "<<", 2)
        // To implement: need to refactor find_offset to accept reference directly
        // or create a proper CigarParser with initialized reference state.
        todo!("Requires architecture changes to test find_offset in isolation")
    }

    #[test]
    fn test_cigar_modifier_mapped_read_no_change() {
        use crate::data::reference::Reference;
        use crate::data::region::Region;
        use crackle_kit::data::bases::rev_comp::RevComplementor;
        use rust_htslib::bam::{Read, Record, Reader};

        let conf = Configuration {
            chimeric_filter: true,
            ..Default::default()
        };

        let _ = INSTANCE.set(GlobalReadOnlyScope {
            conf,
            ..Default::default()
        });

        let bam_path = "/home/eck/workspace/vardict_rs/test_data/test_168714.bam";
        let mut reader = Reader::from_path(bam_path).expect("Failed to open test BAM");

        let mut target: Option<Record> = None;
        for result in reader.records() {
            let record = result.expect("Failed to read BAM record");
            let qname = std::str::from_utf8(record.qname()).unwrap_or("");
            if qname == "SRR098401.7003120" && !record.is_unmapped() {
                target = Some(record);
                break;
            }
        }

        let record = target.expect("Mapped SRR098401.7003120 read not found");

        let ref_seq = vec![b'N'; 500];
        let reference = Reference::new_with_start(ref_seq, 1);
        let region = Region::new("20".to_string(), 168600, 168800, "test_region".to_string());

        let mut rev_complementor = RevComplementor::new();
        let query_seq_owned: Vec<u8> = record.seq().into_decoded_base_iter().collect();
        let query_qual_owned: Vec<u8> = record.qual().to_vec();

        let cigar_view = record.cigar();
        let mut modifier = CigarModifier::new(
            record.pos(),
            &cigar_view,
            query_seq_owned.as_slice(),
            query_qual_owned.as_slice(),
            &reference,
            0,
            query_seq_owned.len(),
            &region,
            &mut rev_complementor,
        );

        let modified = modifier.modify_cigar().expect("modify_cigar failed");
        let modified_cigar = CigarString(modified.cigar.into_iter().collect()).into_view(0);

        assert_eq!(modified.align_start_pos, record.pos());
        assert_eq!(cigar_to_string(&modified_cigar), cigar_to_string(&cigar_view));
        assert_eq!(modified.query_seq, query_seq_owned.as_slice());
        assert_eq!(modified.query_qual, query_qual_owned.as_slice());
    }
}
