use std::{collections::VecDeque, sync::Arc};

use anyhow::{Error, anyhow};
use crackle_kit::{
    data::bases::rev_comp::RevComplementor,
    tracing::{Level, event},
};
use rust_htslib::bam::{
    Record, Writer,
    record::{Cigar, CigarString, CigarStringView},
};

use crate::{
    data::{reference::Reference, region::Region},
    mods::cigar_modifier::CigarModifier,
    scopedata::global_read_only_scope::{GlobalReadOnlyScope, instance},
    utils::aligner::Aligner,
};

pub struct CigarParser {
    query_seq_buf: Vec<u8>,
    query_qual_buf: Vec<u8>,
    aligner: Aligner,
    instance: Arc<GlobalReadOnlyScope>,
    reference: Reference,
    max_read_len: usize,
    region: Region,
    rev_complementor: RevComplementor,
    discordant_count: usize,
}

impl CigarParser {
    fn parse_cigar(&mut self, record: &mut Record) -> Result<(), Error> {
        self.query_seq_buf
            .extend(record.seq().into_decoded_base_iter());
        self.query_qual_buf.extend(record.qual());

        let mut query_seq = self.query_seq_buf.as_slice();
        let mut query_qual = self.query_qual_buf.as_slice();
        let mapping_quality = record.mapq();

        record.cache_cigar_if_empty();
        let mut cigar = record.cigar();

        let ins_del_len = get_ins_del_len(&cigar);

        let tot_nm = match record.aux_option(self.aligner.nm_tag())? {
            Some(rust_htslib::bam::record::Aux::I32(nm)) => nm - ins_del_len as i32,
            Some(oth) => Err(anyhow!("Got non i32 type for NM tag: {:?}", oth))?,
            None => {
                if !cigar.is_empty() {
                    event!(Level::WARN, "No NM tag for mismatches: {:?}", record);
                }

                if record.is_unmapped() || cigar.is_empty() {
                    return Ok(());
                }

                0
            }
        };

        let is_mate_on_the_same_contig = record.tid() == record.mtid();
        let nm = tot_nm;
        let direction = record.is_reverse();

        if self.instance.amplicon_based_calling {
            todo!()
        }

        let mut align_start_pos = 0;
        let mut read_pos_including_softclip = 0;
        let mut read_pos_excluding_softclip = 0;

        if self.instance.conf.perform_local_realignment {
            // Modify the CIGAR for potential mis-alignment for indels at the end of reads to softclipping and let VarDict's
            // algorithm to figure out indels

            let mut cigar_modifier = CigarModifier::new(
                record.query_alignment_start(),
                &cigar,
                query_seq,
                query_qual,
                &self.reference,
                ins_del_len,
                self.max_read_len,
                &self.region,
                &mut self.rev_complementor,
            );

            let mc = cigar_modifier.modify_cigar()?;

            align_start_pos = mc.align_start_pos;
            cigar = CigarString(mc.cigar.into_iter().collect::<Vec<_>>())
                .into_view(mc.align_start_pos as i64);

            query_qual = mc.query_qual;
            query_seq = mc.query_seq;
        } else {
            align_start_pos = record.query_alignment_start();
        }

        self.clean_up_cigar(record);

        let offset = 0;

        //determine discordant reads
        if record.tid() != record.mtid() {
            self.discordant_count += 1;
        }

        //Ignore reads that are softclipped at both ends and both greater than 10 bp
        match cigar.0.as_slice() {
            [Cigar::SoftClip(sl1), .., Cigar::SoftClip(sl2)] if *sl1 >= 10 && *sl2 >= 10 => {
                return Ok(());
            }
            _ => {}
        }

        // Only match and insertion counts toward read length
        // For total length, including soft-clipped bases
        let read_match_ins_len = get_match_insertion_length(&cigar);

        if instance().conf.min_match != 0 && read_match_ins_len < instance().conf.min_match as usize
        {
            return Ok(());
        }

        // The total length, including soft-clipped bases
        let read_len_including_softclips = get_soft_clipped_length(&cigar);
        self.max_read_len = read_len_including_softclips.max(self.max_read_len);

        // If supplementary alignment is present
        if instance().conf.sam_filter != 0 && record.is_supplementary() {
            return Ok(()); // Ignore the supplementary for now so that it won't skew the coverage
        }

        // Skip sites that are not in region of interest in CRISPR mode
        if self.skip_sites_out_region_of_interest(cigar.as_slice(), align_start_pos) {
            return Ok(());
        }

        todo!()
    }

    fn clean_up_cigar(&self, record: &mut Record) {
        record.cache_cigar_if_empty();
        let cigar_vec = record.cigar_cached().unwrap().0.as_slice();

        if cigar_vec.len() > 0 {
            let mut no_matches_yet = true;

            // find start index of leading part (non matches)
            let leading_part_idx = cigar_vec
                .iter()
                .rposition(|c| c.consumes_read_bases() && c.consumes_reference_bases())
                .unwrap_or(0);

            let new_cigar_vec = cigar_vec
                .into_iter()
                .enumerate()
                .filter_map(|(i, &cigar)| {
                    if no_matches_yet {
                        match cigar {
                            Cigar::Ins(l) => Some(Cigar::SoftClip(l)),
                            Cigar::HardClip(_) => None, // remove hard clip
                            cigar
                                if cigar.consumes_read_bases()
                                    && cigar.consumes_reference_bases() =>
                            {
                                no_matches_yet = false;
                                Some(cigar)
                            }
                            oth => Some(oth),
                        }
                    } else {
                        if i <= leading_part_idx {
                            Some(cigar)
                        } else {
                            match cigar {
                                Cigar::Ins(l) => Some(Cigar::SoftClip(l)),
                                Cigar::HardClip(_) => None, // remove hard clip
                                oth => Some(oth),
                            }
                        }
                    }
                })
                .collect::<Vec<_>>();

            record.set_cigar(Some(&CigarString(new_cigar_vec)));
        }
    }

    fn skip_sites_out_region_of_interest(&self, cigars: &[Cigar], align_start_pos: usize) -> bool {
        let cut_site = instance().conf.crispr_cutting_site;
        let filter_bp = instance().conf.crispr_cutting_site;

        if cut_site != 0 {
            //The total aligned length, excluding soft-clipped bases and insertions
            let rlen3 = cigars
                .iter()
                .map(|cigar| match cigar {
                    Cigar::Match(l) | Cigar::Del(l) | Cigar::Equal(l) | Cigar::Diff(l) => *l,
                    _ => 0,
                })
                .sum::<u32>();

            if filter_bp != 0 {
                return !(cut_site - align_start_pos as i32 > filter_bp
                    && align_start_pos as i32 + rlen3 as i32 - cut_site > filter_bp);
            }
        }

        false
    }
}

macro_rules! get_cached_cigar_or_make {
    ($record:expr) => {{
        if $record.cigar_cached().is_none() {
            $record.cache_cigar();
        }

        $record.cigar_cached().unwrap()
    }};
}
use get_cached_cigar_or_make;

#[inline]
fn get_ins_del_len(cigar: &CigarStringView) -> u32 {
    cigar
        .iter()
        .filter_map(|c| match c {
            rust_htslib::bam::record::Cigar::Ins(l) | rust_htslib::bam::record::Cigar::Del(l) => {
                Some(*l)
            }
            _ => None,
        })
        .sum::<u32>()
}

#[inline]
fn get_match_insertion_length(cigar: &CigarStringView) -> usize {
    cigar
        .iter()
        .map(|c| match c {
            Cigar::Match(l) | Cigar::Ins(l) => *l as usize,
            _ => 0,
        })
        .sum::<usize>()
}

#[inline]
fn get_soft_clipped_length(cigar: &CigarStringView) -> usize {
    cigar
        .iter()
        .map(|c| match c {
            Cigar::Match(l) | Cigar::Ins(l) | Cigar::SoftClip(l) => *l as usize,
            _ => 0,
        })
        .sum::<usize>()
}
