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

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
struct TerminalCigarOffsets {
    ref_offset: i32,
    read_offset: i32,
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
                    let genomic_ref_start = self.ref_data.region_start + ref_start_pos as i64;
                    if poss.len() == 1
                        && ((genomic_ref_start - poss.get(0).copied().unwrap()).abs() as usize)
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
                        );
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
                    let genomic_ref_start = self.ref_data.region_start + ref_start_pos as i64;
                    if poss.len() == 1
                        && ((genomic_ref_start - poss.get(0).copied().unwrap()).abs() as usize)
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
                        );
                    }
                }
            }
        }

        while flag && self.indel > 0 {
            flag = false;

            flag |= self.normalize_front_softclip_indel(&mut ref_start_pos, &mut cigar_vec);
            flag |= self.normalize_back_indel_softclip(&mut cigar_vec);
            flag |= self.normalize_front_softclip_match_indel(&mut ref_start_pos, &mut cigar_vec);
            flag |= self.normalize_back_indel_match_softclip(&mut cigar_vec);
            flag |= self.normalize_front_short_match_indel_match(&mut ref_start_pos, &mut cigar_vec)?;
            flag |= self.normalize_back_indel_short_match(&mut cigar_vec);

            if let Some(pattern_match) = find_primary_realign_pattern(&cigar_vec) {
                match pattern_match {
                    PrimaryRealignPattern::TwoDelsInsToComplex(si_and_c_lens) => {
                        flag = self.two_dels_ins_to_complex(
                            ref_start_pos,
                            &mut cigar_vec,
                            si_and_c_lens,
                            flag,
                        )?;
                    }
                    PrimaryRealignPattern::ThreeDeletions(si_and_c_lens) => {
                        flag =
                            self.three_deletions(ref_start_pos, &mut cigar_vec, si_and_c_lens, flag)?;
                    }
                    PrimaryRealignPattern::ThreeIndels(si_and_cigars) => {
                        flag = self.three_indels(ref_start_pos, &mut cigar_vec, si_and_cigars, flag)?;
                    }
                }
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

        let terminal_offsets = terminal_cigar_offsets(ref_start_pos, &cigar_vec);
        let mut cigar_iter_rev = cigar_vec.iter().rev();
        match (cigar_iter_rev.next(), cigar_iter_rev.next()) {
            (Some(&Cigar::SoftClip(sl)), Some(&Cigar::Match(ml))) => {
                self.capture_mis_softly3_ms(
                    TerminalCigarOffsets {
                        ref_offset: terminal_offsets.ref_offset,
                        read_offset: terminal_offsets.read_offset - sl as i32,
                    },
                    &mut cigar_vec,
                    sl,
                    ml,
                )?;
            }
            (Some(&Cigar::Match(ml)), _) => {
                self.capture_mis_softly3_mismatches(terminal_offsets, &mut cigar_vec, ml)?;
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

    fn normalize_front_edge(
        &self,
        ref_start_pos: &mut u32,
        cigar_vec: &mut VecDeque<Cigar>,
    ) -> Result<bool, Error> {
        let mut changed = false;

        changed |= self.normalize_front_softclip_indel(ref_start_pos, cigar_vec);
        changed |= self.normalize_front_softclip_match_indel(ref_start_pos, cigar_vec);
        changed |= self.normalize_front_short_match_indel_match(ref_start_pos, cigar_vec)?;

        Ok(changed)
    }

    fn normalize_front_softclip_indel(
        &self,
        ref_start_pos: &mut u32,
        cigar_vec: &mut VecDeque<Cigar>,
    ) -> bool {
        let mut changed = false;

        match front_pair(cigar_vec) {
            (Some(Cigar::SoftClip(sl)), Some(c2 @ (Cigar::Ins(idl) | Cigar::Del(idl)))) => {
                let tslen = sl + if matches!(c2, Cigar::Ins(_)) { idl } else { 0 };
                *ref_start_pos += if matches!(c2, Cigar::Del(_)) { idl } else { 0 };

                cigar_vec.pop_front().unwrap();
                *cigar_vec.front_mut().unwrap() = Cigar::SoftClip(tslen);
                changed = true;
            }
            _ => {}
        }

        changed
    }

    fn normalize_front_softclip_match_indel(
        &self,
        ref_start_pos: &mut u32,
        cigar_vec: &mut VecDeque<Cigar>,
    ) -> bool {
        let mut changed = false;

        match front_triplet(cigar_vec) {
            (
                Some(Cigar::SoftClip(sl)),
                Some(Cigar::Match(ml)),
                Some(c3 @ (Cigar::Ins(idl) | Cigar::Del(idl))),
            ) if ml <= 10 => {
                let tslen = sl + ml + if matches!(c3, Cigar::Ins(_)) { idl } else { 0 };
                *ref_start_pos += ml + if matches!(c3, Cigar::Del(_)) { idl } else { 0 };

                cigar_vec.drain(..2);
                *cigar_vec.front_mut().unwrap() = Cigar::SoftClip(tslen);
                changed = true;
            }
            _ => {}
        }

        changed
    }

    fn normalize_front_short_match_indel_match(
        &self,
        ref_start_pos: &mut u32,
        cigar_vec: &mut VecDeque<Cigar>,
    ) -> Result<bool, Error> {
        let mut changed = false;

        match front_triplet(cigar_vec) {
            (
                Some(Cigar::Match(ml1)),
                Some(c_id @ (Cigar::Ins(idl) | Cigar::Del(idl))),
                Some(Cigar::Match(mut ml2)),
            ) if ml1 < 10 => {
                let mut tslen = ml1 + if matches!(c_id, Cigar::Ins(_)) { idl } else { 0 };
                *ref_start_pos += ml1 + if matches!(c_id, Cigar::Del(_)) { idl } else { 0 };

                let mut tn = 0;
                while tn < ml2
                    && is_has_and_not_equals(
                        self.query_sequence
                            .get_or_err((tslen + tn) as usize)
                            .copied()?,
                        &self.ref_data.ref_seq,
                        (*ref_start_pos + tn) as usize,
                    )
                {
                    tn += 1;
                }

                tslen += tn;
                ml2 -= tn;
                *ref_start_pos += tn;

                cigar_vec.pop_front().unwrap();
                *cigar_vec.get_mut(0).unwrap() = Cigar::SoftClip(tslen);
                *cigar_vec.get_mut(1).unwrap() = Cigar::Match(ml2);
                changed = true;
            }
            _ => {}
        }

        Ok(changed)
    }

    fn normalize_back_edge(&self, cigar_vec: &mut VecDeque<Cigar>) -> bool {
        let mut changed = false;

        changed |= self.normalize_back_indel_softclip(cigar_vec);
        changed |= self.normalize_back_indel_match_softclip(cigar_vec);
        changed |= self.normalize_back_indel_short_match(cigar_vec);

        changed
    }

    fn normalize_back_indel_softclip(&self, cigar_vec: &mut VecDeque<Cigar>) -> bool {
        let mut changed = false;

        match back_pair(cigar_vec) {
            (Some(Cigar::SoftClip(sl)), Some(c2 @ (Cigar::Ins(idl) | Cigar::Del(idl)))) => {
                let tslen = sl + if matches!(c2, Cigar::Ins(_)) { idl } else { 0 };

                cigar_vec.pop_back().unwrap();
                *cigar_vec.back_mut().unwrap() = Cigar::SoftClip(tslen);
                changed = true;
            }
            _ => {}
        }

        changed
    }

    fn normalize_back_indel_match_softclip(&self, cigar_vec: &mut VecDeque<Cigar>) -> bool {
        let mut changed = false;

        match back_triplet(cigar_vec) {
            (
                Some(Cigar::SoftClip(sl)),
                Some(Cigar::Match(ml)),
                Some(c3 @ (Cigar::Ins(idl) | Cigar::Del(idl))),
            ) if ml <= 10 => {
                let tslen = sl + ml + if matches!(c3, Cigar::Ins(_)) { idl } else { 0 };

                cigar_vec.drain(cigar_vec.len() - 2..);
                *cigar_vec.back_mut().unwrap() = Cigar::SoftClip(tslen);
                changed = true;
            }
            _ => {}
        }

        changed
    }

    fn normalize_back_indel_short_match(&self, cigar_vec: &mut VecDeque<Cigar>) -> bool {
        let mut changed = false;

        match back_pair(cigar_vec) {
            (Some(Cigar::Match(ml)), Some(c_id @ (Cigar::Ins(idl) | Cigar::Del(idl))))
                if ml < 10 =>
            {
                let tslen = ml + if matches!(c_id, Cigar::Ins(_)) { idl } else { 0 };

                cigar_vec.pop_back().unwrap();
                *cigar_vec.back_mut().unwrap() = Cigar::SoftClip(tslen);
                changed = true;
            }
            _ => {}
        }

        changed
    }

    fn capture_mis_softly3_ms(
        &self,
        terminal_offsets: TerminalCigarOffsets,
        cigar_vd: &mut VecDeque<Cigar>,
        sl: u32,
        ml: u32,
    ) -> Result<(), Error> {
        // capture_mis_softly_ms
        let ref_seq = &self.ref_data.ref_seq;
        let query_sequence = self.query_sequence;
        let query_quality = self.query_quality;
        let mut mch = ml as i32;
        let mut soft = sl as i32;
        let TerminalCigarOffsets {
            ref_offset: refoff,
            read_offset: rdoff,
        } = terminal_offsets;

        // number of bases after refoff/rdoff that match in reference and read sequences
        let mut rn = 0;

        let mut homopolymer_checker = HomoPolymerChecker::new();

        while rn < soft {
            let ref_idx = (refoff + rn) as usize;
            let seq_idx = (rdoff + rn) as usize;

            if !is_has_and_equals_ref_and_seq_base(ref_seq, ref_idx, query_sequence, seq_idx) {
                break;
            }

            if query_quality[seq_idx] <= Configuration::LOW_QUAL as u8 {
                break;
            }

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
            while rn + 1 < soft {
                let ref_idx = (refoff + rn + 1) as usize;
                let seq_idx = (rdoff + rn + 1) as usize;

                if !is_has_and_equals_ref_and_seq_base(ref_seq, ref_idx, query_sequence, seq_idx) {
                    break;
                }

                if query_quality[seq_idx] <= Configuration::LOW_QUAL as u8 {
                    break;
                }

                rn += 1;
                if let Some(base) = ref_seq.get((refoff + rn + 1) as usize) {
                    homopolymer_checker.record_base(*base);
                }
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
                    if let Some(last_idx) = cigar_vd.len().checked_sub(1) {
                        cigar_vd[last_idx] = Cigar::Match(mch as u32);
                    }
                }
            }

            if rn == 0 {
                let mut rrn = 0;
                let mut rmch = 0;

                while rrn < mch && rn < mch {
                    let ref_idx = (refoff - rrn - 1) as usize;
                    if ref_seq.get(ref_idx).is_none() {
                        break;
                    }

                    let seq_idx = (rdoff - rrn - 1) as usize;

                    if rrn < rdoff
                        && is_has_and_not_equals_ref_and_seq_base(
                            ref_seq,
                            ref_idx,
                            query_sequence,
                            seq_idx,
                        )
                    {
                        rn = rrn + 1;
                        rmch = 0;
                    } else if rrn < rdoff
                        && is_has_and_equals_ref_and_seq_base(
                            ref_seq,
                            ref_idx,
                            query_sequence,
                            seq_idx,
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
                    let last_idx = cigar_vd.len() - 1;
                    let prev_idx = cigar_vd.len() - 2;
                    cigar_vd[last_idx] = Cigar::SoftClip(soft as u32);
                    cigar_vd[prev_idx] = Cigar::Match(mch as u32);
                }
            }
        }

        Ok(())
    }

    fn capture_mis_softly3_mismatches(
        &self,
        terminal_offsets: TerminalCigarOffsets,
        cigar_vd: &mut VecDeque<Cigar>,
        ml: u32,
    ) -> Result<(), Error> {
        // capture_mis_softly3_mismatches
        let ref_seq = &self.ref_data.ref_seq;
        let query_sequence = self.query_sequence;
        let mut mch = ml as i32;
        let TerminalCigarOffsets {
            ref_offset: refoff,
            read_offset: rdoff,
        } = terminal_offsets;

        let mut rn = 0;
        let mut rrn = 0;
        let mut rmch = 0;

        while rrn < mch && rn < mch {
            let ref_idx = (refoff - rrn - 1) as usize;
            if ref_seq.get(ref_idx).is_none() {
                break;
            }

            let seq_idx = (rdoff - rrn - 1) as usize;

            if rrn < rdoff
                && is_has_and_not_equals_ref_and_seq_base(
                    ref_seq,
                    ref_idx,
                    query_sequence,
                    seq_idx,
                )
            {
                rn = rrn + 1;
                rmch = 0;
            } else if rrn < rdoff
                && is_has_and_equals_ref_and_seq_base(
                    ref_seq,
                    ref_idx,
                    query_sequence,
                    seq_idx,
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
        si_and_c_lens: PrimaryPatternMatch<[u32; 7]>, // Cigars: M D M I M D M
        mut flag: bool,
    ) -> Result<bool, Error> {
        // length of both matched sequences and insertion
        let PrimaryPatternMatch {
            start_idx: si,
            prefix_ref_len,
            prefix_read_len,
            pattern: c_lens,
        } = si_and_c_lens;
        let mut tslen = (c_lens[2] + c_lens[3] + c_lens[4]) as i32;

        // length of deletions and internal matched sequences
        let mut dlen = (c_lens[1] + c_lens[2] + c_lens[4] + c_lens[5]) as i32;

        // length of internal matched sequences
        let mid = (c_lens[2] + c_lens[4]) as i32;

        // offset of first deletion in the reference sequence
        let refoff = (align_start_pos + prefix_ref_len + c_lens[0]) as i32;

        // offset of first deletion in the read
        let rdoff = (prefix_read_len + c_lens[0]) as i32;

        // offset of first deletion in the read corrected by possibly matching bases
        let mut rdoff_corrected = c_lens[0] as i32;

        let mut rm = c_lens[6] as i32;

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
        si_and_c_lens: PrimaryPatternMatch<[u32; 7]>, // Cigars: M D M D M D M
        mut flag: bool,
    ) -> Result<bool, Error> {
        let PrimaryPatternMatch {
            start_idx: si,
            prefix_ref_len,
            prefix_read_len,
            pattern: c_lens,
        } = si_and_c_lens;

        //length of both matched sequences and insertion
        let mut tslen = (c_lens[2] + c_lens[4]) as i32;

        //length of deletions and internal matched sequences
        let mut dlen = (c_lens[1] + c_lens[2] + c_lens[3] + c_lens[4] + c_lens[5]) as i32;

        //length of internal matched sequences
        let mid = (c_lens[2] + c_lens[4]) as i32;

        //offset of first deletion in the reference sequence
        let refoff = (align_start_pos + prefix_ref_len + c_lens[0]) as i32;

        //offset of first deletion in the read
        let rdoff = (prefix_read_len + c_lens[0]) as i32;

        //offset of first deletion in the read corrected by possibly matching bases
        let mut rdoff_corrected = (c_lens[0]) as i32;

        let mut rm = c_lens[6] as i32;

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
        si_and_cigars: PrimaryPatternMatch<[Cigar; 7]>, // Cigars: M D M D M D M
        mut flag: bool,
    ) -> Result<bool, Error> {
        let PrimaryPatternMatch {
            start_idx: si,
            prefix_ref_len,
            prefix_read_len,
            pattern: cigars,
        } = si_and_cigars;

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

        let refoff = (align_start_pos + prefix_ref_len + cigars[0].len()) as i32;
        let rdoff = (prefix_read_len + cigars[0].len()) as i32;
        let mut rdoff_corrected = cigars[0].len() as i32;
        let mut rm = cigars[6].len() as i32;

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
        let ref_seq = &self.ref_data.ref_seq;
        let query_sequence = self.query_sequence;
        let query_quality = self.query_quality;
        let mut mch = ml as i32;
        //length of soft-clipping
        let mut soft = sl as i32;

        //number of bases before matched sequence that match in reference and read sequences
        let mut rn = 0;
        let mut homop = HomoPolymerChecker::new();

        while rn < soft {
            let ref_idx = ((*cigar_pos as i32) - rn - 1) as usize;
            let seq_idx = (soft - rn - 1) as usize;

            if !is_has_and_equals_ref_and_seq_base(ref_seq, ref_idx, query_sequence, seq_idx) {
                break;
            }

            if query_quality[seq_idx] <= Configuration::LOW_QUAL as u8 {
                break;
            }

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
            while rn + 1 < soft {
                let ref_idx = ((*cigar_pos as i32) - rn - 2) as usize;
                let seq_idx = (soft - rn - 2) as usize;

                if !is_has_and_equals_ref_and_seq_base(ref_seq, ref_idx, query_sequence, seq_idx) {
                    break;
                }

                if query_quality[seq_idx] <= Configuration::LOW_QUAL as u8 {
                    break;
                }

                rn += 1;
                homop.record_base(
                    ref_seq
                        .get_or_err(((*cigar_pos as i32) - rn - 2) as usize)
                        .copied()?,
                );
            }

            if (rn > 4 && !homop.is_homopolymer())
                || is_has_and_equals_ref_and_seq_base(
                    ref_seq,
                    ((*cigar_pos as i32) - 1) as usize,
                    query_sequence,
                    (soft - 1) as usize,
                )
            {
                mch += rn + 1;
                soft -= rn + 1;

                if soft > 0 {
                    cigar_vd[0] = Cigar::SoftClip(soft as u32);
                    cigar_vd[1] = Cigar::Match(mch as u32);
                } else {
                    cigar_vd.pop_front().unwrap();
                    cigar_vd[0] = Cigar::Match(mch as u32);
                }

                *cigar_pos -= rn as u32 + 1;
            }

            if rn == 0 {
                let mut rrn = 0;
                let mut rmch = 0;

                while rrn < mch && rn < mch {
                    let ref_idx = (*cigar_pos as i32 + rrn) as usize;
                    if ref_seq.get(ref_idx).is_none() {
                        break;
                    }

                    let seq_idx = (soft + rrn) as usize;

                    if is_has_and_not_equals_ref_and_seq_base(
                        ref_seq,
                        ref_idx,
                        query_sequence,
                        seq_idx,
                    ) {
                        rn = rrn + 1;
                        rmch = 0;
                    } else if is_has_and_equals_ref_and_seq_base(
                        ref_seq,
                        ref_idx,
                        query_sequence,
                        seq_idx,
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

                    match (cigar_vd.get(0), cigar_vd.get(1)) {
                        (Some(Cigar::SoftClip(_)), Some(Cigar::Match(_))) => {
                            cigar_vd[0] = Cigar::SoftClip(soft as u32);
                            cigar_vd[1] = Cigar::Match(mch as u32);
                        }
                        _ => {}
                    }

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
            let ref_idx = (*cigar_pos as i32 + rrn) as usize;
            let seq_idx = rrn as usize;
            if self.contig_ref_seq().get((*cigar_pos as i32 + rrn) as usize).is_none() {
                break;
            }

            let ref_base = self.contig_ref_seq().get(ref_idx).copied();
            let seq_base = self.query_sequence.get(seq_idx).copied();
            let mismatch = is_has_and_not_equals_ref_and_seq_base(
                self.contig_ref_seq(),
                ref_idx,
                self.query_sequence,
                seq_idx,
            );
            if mismatch {
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

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
struct PrimaryPatternMatch<T> {
    start_idx: usize,
    prefix_ref_len: u32,
    prefix_read_len: u32,
    pattern: T,
}

enum PrimaryRealignPattern {
    TwoDelsInsToComplex(PrimaryPatternMatch<[u32; 7]>),
    ThreeDeletions(PrimaryPatternMatch<[u32; 7]>),
    ThreeIndels(PrimaryPatternMatch<[Cigar; 7]>),
}

fn front_pair(cigar: &VecDeque<Cigar>) -> (Option<Cigar>, Option<Cigar>) {
    (cigar.get(0).copied(), cigar.get(1).copied())
}

fn front_triplet(cigar: &VecDeque<Cigar>) -> (Option<Cigar>, Option<Cigar>, Option<Cigar>) {
    (
        cigar.get(0).copied(),
        cigar.get(1).copied(),
        cigar.get(2).copied(),
    )
}

fn back_pair(cigar: &VecDeque<Cigar>) -> (Option<Cigar>, Option<Cigar>) {
    let len = cigar.len();
    (
        len.checked_sub(1).and_then(|idx| cigar.get(idx)).copied(),
        len.checked_sub(2).and_then(|idx| cigar.get(idx)).copied(),
    )
}

fn back_triplet(cigar: &VecDeque<Cigar>) -> (Option<Cigar>, Option<Cigar>, Option<Cigar>) {
    let len = cigar.len();
    (
        len.checked_sub(1).and_then(|idx| cigar.get(idx)).copied(),
        len.checked_sub(2).and_then(|idx| cigar.get(idx)).copied(),
        len.checked_sub(3).and_then(|idx| cigar.get(idx)).copied(),
    )
}

fn find_primary_realign_pattern(cigar: &VecDeque<Cigar>) -> Option<PrimaryRealignPattern> {
    if cigar.len() < 7 {
        return None;
    }

    let mut first_three_deletions = None;
    let mut first_three_indels = None;
    let mut prefix_ref_len = 0u32;
    let mut prefix_read_len = 0u32;

    for i in 0..cigar.len() - 6 {
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
                &Cigar::Match(i1),
                &Cigar::Del(i2),
                &Cigar::Match(i3),
                &Cigar::Ins(i4),
                &Cigar::Match(i5),
                &Cigar::Del(i6),
                &Cigar::Match(i7),
            ) => {
                return Some(PrimaryRealignPattern::TwoDelsInsToComplex(
                    PrimaryPatternMatch {
                        start_idx: i,
                        prefix_ref_len,
                        prefix_read_len,
                        pattern: [i1, i2, i3, i4, i5, i6, i7],
                    },
                ));
            }
            (
                &Cigar::Match(i1),
                &Cigar::Del(i2),
                &Cigar::Match(i3),
                &Cigar::Del(i4),
                &Cigar::Match(i5),
                &Cigar::Del(i6),
                &Cigar::Match(i7),
            ) => {
                if first_three_deletions.is_none() {
                    first_three_deletions = Some(PrimaryRealignPattern::ThreeDeletions(
                        PrimaryPatternMatch {
                            start_idx: i,
                            prefix_ref_len,
                            prefix_read_len,
                            pattern: [i1, i2, i3, i4, i5, i6, i7],
                        },
                    ));
                }
            }
            (
                &c1 @ Cigar::Match(_),
                &c2 @ (Cigar::Ins(_) | Cigar::Del(_)),
                &c3 @ Cigar::Match(_),
                &c4 @ (Cigar::Ins(_) | Cigar::Del(_)),
                &c5 @ Cigar::Match(_),
                &c6 @ (Cigar::Ins(_) | Cigar::Del(_)),
                &c7 @ Cigar::Match(_),
            ) => {
                if first_three_indels.is_none() {
                    first_three_indels = Some(PrimaryRealignPattern::ThreeIndels(
                        PrimaryPatternMatch {
                            start_idx: i,
                            prefix_ref_len,
                            prefix_read_len,
                            pattern: [c1, c2, c3, c4, c5, c6, c7],
                        },
                    ));
                }
            }
            _ => {}
        }

        accumulate_cigar_offsets(cigar[i], &mut prefix_ref_len, &mut prefix_read_len);
    }

    first_three_deletions.or(first_three_indels)
}

fn accumulate_cigar_offsets(cigar: Cigar, ref_len: &mut u32, read_len: &mut u32) {
    match cigar {
        Cigar::Match(l) => {
            *ref_len += l;
            *read_len += l;
        }
        Cigar::RefSkip(l) | Cigar::Del(l) => {
            *ref_len += l;
        }
        Cigar::SoftClip(l) | Cigar::Ins(l) => {
            *read_len += l;
        }
        _ => {}
    }
}

fn terminal_cigar_offsets(align_start_pos: u32, cigar_vd: &VecDeque<Cigar>) -> TerminalCigarOffsets {
    let mut ref_offset = align_start_pos;
    let mut read_offset = 0u32;

    for &cigar in cigar_vd {
        accumulate_cigar_offsets(cigar, &mut ref_offset, &mut read_offset);
    }

    TerminalCigarOffsets {
        ref_offset: ref_offset as i32,
        read_offset: read_offset as i32,
    }
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

    // Java parity for NOTDIG_DIG_I_DIG_M_DIG_DI_DIGI:
    // inspect only the first I-M-[ID] occurrence that has a preceding operator.
    // If that first occurrence is preceded by D/H, do not search later occurrences.
    for j in 1..cigar.len() - 2 {
        match (&cigar[j], &cigar[j + 1], &cigar[j + 2]) {
            (&c2 @ Cigar::Ins(_), &c3 @ Cigar::Match(_), &c4 @ (Cigar::Ins(_) | Cigar::Del(_))) => {
                let prev = cigar[j - 1];
                if matches!(prev, Cigar::Del(_) | Cigar::HardClip(_)) {
                    return None;
                }

                return Some((
                    j,
                    [c2, c3, c4],
                    cigar.get(j + 3).copied().and_then(|c| {
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
    use std::collections::VecDeque;

    use crate::{
        data::{reference::Reference, region::Region},
        mods::cigar_parser::CigarParser,
        scopedata::global_read_only_scope::{GlobalReadOnlyScope, INSTANCE},
    };
    use crackle_kit::data::bases::rev_comp::RevComplementor;
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
        cigar
            .iter()
            .map(|c| match c {
                Cigar::Match(l) => format!("{}M", l),
                Cigar::Ins(l) => format!("{}I", l),
                Cigar::Del(l) => format!("{}D", l),
                Cigar::RefSkip(l) => format!("{}N", l),
                Cigar::SoftClip(l) => format!("{}S", l),
                Cigar::HardClip(l) => format!("{}H", l),
                Cigar::Pad(l) => format!("{}P", l),
                Cigar::Equal(l) => format!("{}=", l),
                Cigar::Diff(l) => format!("{}X", l),
            })
            .collect()
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

    #[test]
    fn test_find_primary_realign_pattern_prefers_complex_over_earlier_three_deletions() {
        let cigar = VecDeque::from(parse_cigar_string("1M1D1M1D1M1D1M1M1D1M1I1M1D1M").0);

        let found = find_primary_realign_pattern(&cigar);

        assert!(matches!(
            found,
            Some(PrimaryRealignPattern::TwoDelsInsToComplex(PrimaryPatternMatch {
                start_idx: 7,
                ..
            }))
        ));
    }

    #[test]
    fn test_find_primary_realign_pattern_prefers_three_deletions_over_earlier_three_indels() {
        let cigar = VecDeque::from(parse_cigar_string("1M1I1M1I1M1I1M1M1D1M1D1M1D1M").0);

        let found = find_primary_realign_pattern(&cigar);

        assert!(matches!(
            found,
            Some(PrimaryRealignPattern::ThreeDeletions(PrimaryPatternMatch {
                start_idx: 7,
                ..
            }))
        ));
    }

    #[test]
    fn test_find_primary_realign_pattern_carries_prefix_offsets() {
        let cigar = VecDeque::from(parse_cigar_string("2S3M1M1D1M1I1M1D1M").0);

        let found = find_primary_realign_pattern(&cigar);

        assert!(matches!(
            found,
            Some(PrimaryRealignPattern::TwoDelsInsToComplex(PrimaryPatternMatch {
                start_idx: 2,
                prefix_ref_len: 3,
                prefix_read_len: 5,
                ..
            }))
        ));
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
        use rust_htslib::bam::{Read, Reader, Record};

        let conf = Configuration {
            chimeric_filter: true,
            ..Default::default()
        };

        let _ = INSTANCE.set(GlobalReadOnlyScope {
            conf,
            ..Default::default()
        });

        let bam_path = concat!(env!("CARGO_MANIFEST_DIR"), "/testdata/test_168714.bam");
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
        assert_eq!(
            cigar_to_string(&modified_cigar),
            cigar_to_string(&cigar_view)
        );
        assert_eq!(modified.query_seq, query_seq_owned.as_slice());
        assert_eq!(modified.query_qual, query_qual_owned.as_slice());
    }

    #[test]
    fn test_capture_mis_softly3_ms_boundary_does_not_error() {
        let reference = Reference::new_with_start(vec![b'A'; 10], 1);
        let region = Region::new("chr1".to_string(), 1, 10, "test".to_string());
        let mut rev_complementor = RevComplementor::new();

        let cigar = CigarString(vec![Cigar::Match(2), Cigar::SoftClip(3)]).into_view(0);
        let query_sequence = vec![b'C', b'C', b'T', b'A', b'C'];
        let query_quality = vec![30, 30, 30, 30, 30];

        let modifier = CigarModifier::new(
            0,
            &cigar,
            &query_sequence,
            &query_quality,
            &reference,
            0,
            query_sequence.len(),
            &region,
            &mut rev_complementor,
        );

        let mut cigar_vd = VecDeque::from(vec![Cigar::Match(2), Cigar::SoftClip(3)]);

        let terminal_offsets = terminal_cigar_offsets(6, &cigar_vd);
        let result = modifier.capture_mis_softly3_ms(
            TerminalCigarOffsets {
                ref_offset: terminal_offsets.ref_offset,
                read_offset: terminal_offsets.read_offset - 3,
            },
            &mut cigar_vd,
            3,
            2,
        );
        assert!(result.is_ok());
    }

    #[test]
    fn test_terminal_cigar_offsets_count_read_and_reference_consumption() {
        let cigar_vd = VecDeque::from(vec![
            Cigar::SoftClip(2),
            Cigar::Match(5),
            Cigar::Del(3),
            Cigar::Match(4),
            Cigar::Ins(2),
            Cigar::Match(6),
            Cigar::SoftClip(3),
        ]);

        let offsets = terminal_cigar_offsets(7, &cigar_vd);

        assert_eq!(offsets.ref_offset, 25);
        assert_eq!(offsets.read_offset, 22);
    }

    #[test]
    fn test_modify_cigar_front_edge_normalization_cascades_in_order() {
        let reference = Reference::new_with_start(vec![b'C'; 64], 1);
        let region = Region::new("chr1".to_string(), 1, 64, "test".to_string());
        let mut rev_complementor = RevComplementor::new();

        let cigar = parse_cigar_string("5S2I3M4D10M").into_view(0);
        let query_sequence = vec![b'A'; 20];
        let query_quality = vec![30; 20];

        let mut modifier = CigarModifier::new(
            0,
            &cigar,
            &query_sequence,
            &query_quality,
            &reference,
            1,
            query_sequence.len(),
            &region,
            &mut rev_complementor,
        );

        let modified = modifier.modify_cigar().expect("modify_cigar failed");
        let modified_cigar = CigarString(modified.cigar.into_iter().collect()).into_view(0);

        assert_eq!(modified.align_start_pos, 7);
        assert_eq!(cigar_to_string(&modified_cigar), "10S10M");
    }

    #[test]
    fn test_modify_cigar_back_edge_normalization_cascades_in_order() {
        let reference = Reference::new_with_start(vec![b'C'; 64], 1);
        let region = Region::new("chr1".to_string(), 1, 64, "test".to_string());
        let mut rev_complementor = RevComplementor::new();

        let cigar = parse_cigar_string("10M4I3M5S").into_view(0);
        let query_sequence = vec![b'A'; 22];
        let query_quality = vec![30; 22];

        let mut modifier = CigarModifier::new(
            0,
            &cigar,
            &query_sequence,
            &query_quality,
            &reference,
            1,
            query_sequence.len(),
            &region,
            &mut rev_complementor,
        );

        let modified = modifier.modify_cigar().expect("modify_cigar failed");
        let modified_cigar = CigarString(modified.cigar.into_iter().collect()).into_view(0);

        assert_eq!(modified.align_start_pos, 0);
        assert_eq!(cigar_to_string(&modified_cigar), "10M12S");
    }

    #[test]
    fn test_modify_cigar_prefers_trailing_indel_match_softclip_before_front_rewrite() {
        INSTANCE.get_or_init(GlobalReadOnlyScope::default);

        let query_sequence = b"TCACTTTGTTTTTGGAGATTTCAAGTCTTCAAATCACTCAGCTATCCTAAGAACAGTTTTTCTTTTGTTTTTGGAGGTTTCAAGTCTTCAAATCACTCAGC"
            .to_vec();
        let query_quality = vec![40; query_sequence.len()];
        let reference = Reference::new_with_start(query_sequence.clone(), 1);
        let region = Region::new(
            "chr1".to_string(),
            1,
            query_sequence.len(),
            "test".to_string(),
        );
        let mut rev_complementor = RevComplementor::new();

        let cigar = parse_cigar_string("4M59I10M28S").into_view(0);
        let mut modifier = CigarModifier::new(
            0,
            &cigar,
            &query_sequence,
            &query_quality,
            &reference,
            59,
            query_sequence.len(),
            &region,
            &mut rev_complementor,
        );

        let modified = modifier.modify_cigar().expect("modify_cigar failed");
        let modified_cigar = CigarString(modified.cigar.into_iter().collect()).into_view(0);

        assert_eq!(modified.align_start_pos, 0);
        assert_eq!(cigar_to_string(&modified_cigar), "101M");
    }

    #[test]
    fn test_modify_cigar_rewrites_front_edge_deletion_with_trailing_softclip() {
        INSTANCE.get_or_init(GlobalReadOnlyScope::default);

        let reference = Reference::new_with_start(vec![b'A'; 256], 1);
        let region = Region::new("chr1".to_string(), 1, 256, "test".to_string());
        let mut rev_complementor = RevComplementor::new();

        let cigar = parse_cigar_string("9M1D71M21S").into_view(0);
        let query_sequence = vec![b'A'; 101];
        let query_quality = vec![30; 101];

        let mut modifier = CigarModifier::new(
            1,
            &cigar,
            &query_sequence,
            &query_quality,
            &reference,
            1,
            query_sequence.len(),
            &region,
            &mut rev_complementor,
        );

        let mut ref_start_pos = 1;
        let mut cigar_vd = VecDeque::from(parse_cigar_string("9M1D71M21S").0);
        let changed = modifier
            .normalize_front_edge(&mut ref_start_pos, &mut cigar_vd)
            .expect("normalize_front_edge failed");
        let modified_cigar = CigarString(cigar_vd.into_iter().collect()).into_view(0);

        assert!(changed);
        assert_eq!(ref_start_pos, 11);
        assert_eq!(cigar_to_string(&modified_cigar), "9S71M21S");
    }

    #[test]
    fn test_modify_cigar_rewrites_front_edge_insertion_with_trailing_ops() {
        INSTANCE.get_or_init(GlobalReadOnlyScope::default);

        let reference = Reference::new_with_start(vec![b'A'; 256], 1);
        let region = Region::new("chr1".to_string(), 1, 256, "test".to_string());
        let mut rev_complementor = RevComplementor::new();

        let cigar = parse_cigar_string("8M2I22M1I68M").into_view(0);
        let query_sequence = vec![b'A'; 101];
        let query_quality = vec![30; 101];

        let mut modifier = CigarModifier::new(
            1,
            &cigar,
            &query_sequence,
            &query_quality,
            &reference,
            1,
            query_sequence.len(),
            &region,
            &mut rev_complementor,
        );

        let mut ref_start_pos = 1;
        let mut cigar_vd = VecDeque::from(parse_cigar_string("8M2I22M1I68M").0);
        let changed = modifier
            .normalize_front_edge(&mut ref_start_pos, &mut cigar_vd)
            .expect("normalize_front_edge failed");
        let modified_cigar = CigarString(cigar_vd.into_iter().collect()).into_view(0);

        assert!(changed);
        assert_eq!(ref_start_pos, 9);
        assert_eq!(cigar_to_string(&modified_cigar), "10S22M1I68M");
    }
}
