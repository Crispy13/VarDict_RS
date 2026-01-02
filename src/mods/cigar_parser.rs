use bio_types::genome::AbstractInterval;
use std::{
    collections::{HashMap, HashSet, VecDeque},
    sync::Arc,
};

use anyhow::{Error, anyhow};
use crackle_kit::{
    data::bases::rev_comp::RevComplementor,
    tracing::{Level, event},
};
use rust_htslib::bam::{
    Record, Writer,
    ext::BamRecordExtensions,
    record::{self, Cigar, CigarString, CigarStringView},
};

use crate::{
    data::{reference::Reference, region::Region},
    mods::cigar_modifier::CigarModifier,
    scopedata::global_read_only_scope::{GlobalReadOnlyScope, instance},
    utils::aligner::Aligner,
};

type SplicingKey = (i64, i64);
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
    splice_count: HashMap<SplicingKey, Vec<usize>>,
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

        let mut pos = 0;
        let mut read_pos_including_softclip = 0;
        let mut read_pos_excluding_softclip = 0;

        if self.instance.conf.perform_local_realignment {
            // Modify the CIGAR for potential mis-alignment for indels at the end of reads to softclipping and let VarDict's
            // algorithm to figure out indels

            let mut cigar_modifier = CigarModifier::new(
                record.pos(),
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

            pos = mc.align_start_pos;
            cigar = CigarString(mc.cigar.into_iter().collect::<Vec<_>>())
                .into_view(mc.align_start_pos as i64);

            query_qual = mc.query_qual;
            query_seq = mc.query_seq;
        } else {
            pos = record.pos();
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
        if self.skip_sites_out_region_of_interest(cigar.as_slice(), pos) {
            return Ok(());
        }

        // true if mate is in forward forection
        let mate_direction = if !record.is_mate_reverse() {
            false
        } else {
            true
        };

        if record.is_paired() && record.is_mate_unmapped() {
            // TODO
        } else if record.mapq() > 10 && !instance().conf.disable_sv {
            // Consider high mapping quality mates only
            todo!()
        }

        let mpos = record.mpos();
        let adj_pos = pos;
        let offset = 0;

        'process_cigar: {
            //Loop over CIGAR records
            for c in cigar.iter().copied() {
                if self.skip_overlapping_reads(record, adj_pos, pos, direction, mpos) {
                    break;
                }

                //Letter from CIGAR
                match c {
                    Cigar::RefSkip(l) => {
                        self.process_not_matched(&mut adj_pos, &mut offset, l);
                    }
                    Cigar::SoftClip(l) => {}
                    _ => {}
                }
            }
        }

        todo!()
    }

    /// Process CIGAR soft-clipped part. Will ignore large soft-clips and create Variations for mis-softclipping reads
    /// due to alignment
    fn process_soft_clip(
        contig: &str,
        record: &Record,
        query_sequence: &[u8],
        mapq: u8,
        contig_ref_seq: &[u8],
        query_quality: &[u8],
        nm: usize,
        is_reverse: bool,
        pos:i64,
        total_length_including_soft_clipped: usize,
        ci: usize,
    ) {
        //First record in CIGAR
        if ci == 0 {
            // 5' soft clipped
            // Ignore large soft clip due to chimeric reads in library construction
            
        }
    }

    /// N in CIGAR - skipped region from reference
    /// Skip the region and add string start-end to %SPLICE
    fn process_not_matched(&mut self, adj_pos: &mut i64, offset: &mut i32, cigar_len: u32) {
        let key = (*adj_pos - 1, *adj_pos + cigar_len as i64 - 1);

        self.splice_count
            .entry(key)
            .and_modify(|v| v[0] += 1)
            .or_insert_with(|| vec![1]);

        *adj_pos += cigar_len as i64;
        *offset = 0;
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

    /// Skip overlapping reads to avoid double counts for coverage (alignment or second in pair flag is used)
    fn skip_overlapping_reads(
        &self,
        record: &Record,
        adj_pos: i64,
        pos: i64,
        direction: bool, // true is reverse direction
        mate_pos: i64,
    ) -> bool {
        if instance().conf.unique_mode_alignment_enabled
            && is_paired_and_same_chromosome(record)
            && !direction
            && adj_pos as i64 >= mate_pos
        {
            return true;
        }

        if instance().conf.unique_mode_second_in_pair_enabled
            && record.is_last_in_template()
            && is_paired_and_same_chromosome(record)
            && are_reads_overlap(record, adj_pos, pos, mate_pos)
        {
            return true;
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

#[inline]
fn is_paired_and_same_chromosome(record: &Record) -> bool {
    record.is_paired() && record.mtid() == record.tid()
}

/// Check if reads are overlapping. Two cases are considered: if position after mate start and if before.
fn are_reads_overlap(record: &Record, adj_pos: i64, pos: i64, mate_pos: i64) -> bool {
    if pos >= mate_pos {
        adj_pos as i64 >= mate_pos
            && adj_pos as i64
                <= mate_pos
                    + record
                        .cigar_cached()
                        .unwrap()
                        .iter()
                        .map(|c| if c.consumes_read_bases() { c.len() } else { 0 })
                        .sum::<u32>() as i64
    } else {
        adj_pos as i64 >= mate_pos && record.mpos() <= record.reference_end()
    }
}
