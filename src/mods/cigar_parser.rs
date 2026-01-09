use bio_types::genome::AbstractInterval;
use std::{
    collections::{BTreeMap, HashMap, HashSet, VecDeque},
    ops::AddAssign,
    sync::Arc,
};

use anyhow::{Error, anyhow};
use crackle_kit::{
    data::bases::{comp::complement_base, rev_comp::RevComplementor},
    nuc_base_map::NucBaseMap,
    tracing::{Level, event},
};
use rust_htslib::bam::{
    Record, Writer,
    ext::BamRecordExtensions,
    record::{self, Cigar, CigarString, CigarStringView},
};

use crate::{
    conf::Configuration,
    data::{
        patterns::{SA_CIGAR_D_S_3CLIP, SA_CIGAR_D_S_5CLIP},
        reference::Reference,
        region::Region,
    },
    mods::cigar_modifier::CigarModifier,
    scopedata::global_read_only_scope::{GlobalReadOnlyScope, instance},
    utils::{BytesExt, SliceExt, SliceExt2, aligner::Aligner},
    variants::{
        var_utils::{get_variants_from_map, get_variation_from_seq, is_has_and_equals},
        variants::{SoftClip, VarDesc, Variant},
    },
};

type SplicingKey = (i64, i64);
pub struct CigarParser {
    query_seq_buf: Option<Vec<u8>>,
    query_qual_buf: Option<Vec<u8>>,
    aligner: Aligner,
    instance: Arc<GlobalReadOnlyScope>,
    reference: Reference,
    max_read_len: usize,
    region: Region,
    rev_complementor: RevComplementor,
    discordant_count: usize,
    splice_count: HashMap<SplicingKey, Vec<usize>>,

    /// keep track the read position (offset), including softclipped
    read_pos_including_softclip: usize,

    /// keep track the position (offset) in the alignment, excluding softclipped
    read_pos_excluding_softclip: usize,

    /// ref start position of the current read
    start: i64,
    offset: usize,

    cigar: CigarStringView,
    cigar_len: u32,

    non_insertion_vars: HashMap<i64, HashMap<VarDesc, Variant>>,

    ref_coverage: HashMap<i64, usize>,

    soft_clips5_end: HashMap<i64, SoftClip>,
    soft_clips3_end: HashMap<i64, SoftClip>,
}

impl Default for CigarParser {
    fn default() -> Self {
        Self {
            query_seq_buf: Default::default(),
            query_qual_buf: Default::default(),
            aligner: Default::default(),
            instance: Default::default(),
            reference: Default::default(),
            max_read_len: Default::default(),
            region: Default::default(),
            discordant_count: Default::default(),
            splice_count: Default::default(),
            read_pos_including_softclip: Default::default(),
            read_pos_excluding_softclip: Default::default(),
            start: Default::default(),
            offset: Default::default(),
            cigar_len: Default::default(),
            non_insertion_vars: Default::default(),
            ref_coverage: Default::default(),
            soft_clips5_end: Default::default(),
            soft_clips3_end: Default::default(),
            rev_complementor: RevComplementor::new(),
        }
    }
}

impl CigarParser {
    fn parse_cigar(&mut self, record: &mut Record) -> Result<(), Error> {
        let mut query_seq_buf = self.query_seq_buf.take().unwrap();
        let mut query_qual_buf = self.query_qual_buf.take().unwrap();

        query_seq_buf.extend(record.seq().into_decoded_base_iter());
        query_qual_buf.extend(record.qual());

        let mut query_seq = query_seq_buf.as_slice();
        let mut query_qual = query_qual_buf.as_slice();

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
        let is_reverse = record.is_reverse();

        if self.instance.amplicon_based_calling {
            todo!()
        }

        let mut pos = 0;
        self.read_pos_including_softclip = 0;
        self.read_pos_excluding_softclip = 0;

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

        let mut offset = 0;

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
        if self.skip_sites_out_region_of_interest(cigar.as_slice()) {
            return Ok(());
        }

        // true if mate is in forward forection
        let mate_is_reverse = record.is_mate_reverse();

        if record.is_paired() && record.is_mate_unmapped() {
            // TODO
        } else if record.mapq() > 10 && !instance().conf.disable_sv {
            // Consider high mapping quality mates only
            todo!()
        }

        let mpos = record.mpos();
        let adj_pos = pos;
        self.offset = 0;

        'process_cigar: {
            //Loop over CIGAR records
            for (ci, c) in cigar.iter().copied().enumerate() {
                if self.skip_overlapping_reads(record, adj_pos, pos, is_reverse, mpos) {
                    break;
                }

                self.cigar_len = c.len();
                //Letter from CIGAR
                match c {
                    Cigar::RefSkip(l) => {
                        self.process_not_matched(l);
                    }
                    Cigar::SoftClip(_) => self.process_soft_clip(
                        // &self.region.chrom,
                        record,
                        query_seq,
                        mapping_quality,
                        query_qual,
                        nm as usize,
                        is_reverse,
                        pos,
                        read_len_including_softclips,
                        ci,
                        // &mut self.cigar_len,
                        // self.max_read_len,
                        &cigar,
                    )?,
                    Cigar::HardClip(_) => {
                        // hard clipping - skip
                        offset = 0;
                    }
                    Cigar::Ins(_) => {
                        offset = 0;
                    }
                    _ => {}
                }
            }
        }

        // return buf resource to self.
        query_qual_buf.clear();
        query_seq_buf.clear();
        self.query_qual_buf.insert(query_qual_buf);
        self.query_seq_buf.insert(query_seq_buf);

        todo!()
    }

    /// Process CIGAR soft-clipped part. Will ignore large soft-clips and create Variations for mis-softclipping reads
    /// due to alignment
    fn process_soft_clip(
        &mut self,
        record: &Record,
        query_sequence: &[u8],
        mapq: u8,
        query_quality: &[u8],
        num_mismatch: usize,
        is_reverse: bool,
        pos: i64,
        total_length_including_soft_clipped: usize,
        ci: usize,
        cigar: &CigarStringView,
    ) -> Result<(), Error> {
        let mut cigar_len = &mut self.cigar_len;
        let contig = self.region.chrom.as_str();
        let max_read_len = self.max_read_len;

        //First record in CIGAR
        if ci == 0 {
            // 5' soft clipped
            // Ignore large soft clip due to chimeric reads in library construction
            if !instance().conf.chimeric_filter {
                if *cigar_len >= 20
                    && let Some(sa_tag_val) = record.aux_option(b"SA")?
                {
                    if is_read_chimeric_with_sa(
                        record,
                        pos,
                        sa_tag_val.try_get_str()?,
                        is_reverse,
                        true,
                        max_read_len,
                        cigar,
                    )? {
                        self.read_pos_including_softclip += *cigar_len as usize;
                        self.offset = 0;

                        // Had to reset the start due to softclipping adjustment
                        self.start = 0;

                        return Ok(());
                    }
                    // trying to detect chimeric reads even when there's no supplementary
                    // alignment from aligner
                } else if *cigar_len >= Configuration::SEED_1 as u32 {
                    let ref_seed_map = &self.reference.seed;
                    let rev_comp_seq = record
                        .seq()
                        .into_decoded_base_iter()
                        .take(*cigar_len as usize)
                        .rev()
                        .map(complement_base)
                        // .take(Configuration::SEED_1 as usize)
                        .collect::<Vec<_>>();

                    if let Some(poss) =
                        ref_seed_map.get(rev_comp_seq.get_or_err(
                            0..(rev_comp_seq.len().min(Configuration::SEED_1 as usize)),
                        )?)
                    {
                        if poss.len() == 1
                            && self.start - poss.get(0).unwrap() < 2 * self.max_read_len as i64
                        {
                            self.read_pos_including_softclip += *cigar_len as usize;
                            self.offset = 0;
                            // Had to reset the start due to softclipping adjustment
                            self.start = pos;

                            event!(
                                Level::INFO,
                                "{} at 5' is a chimeric at {} by SEED {}",
                                rev_comp_seq.as_slice().try_as_str()?,
                                self.start,
                                Configuration::SEED_1,
                            );

                            return Ok(());
                        }
                    }
                }
            }
            // Align softclipped but matched sequences due to mis-softclipping
            /*
                Conditions:
                1). segment length > 1
                2). start between 0 and chromosome length,
                3). reference genome is known at this position
                4). reference and read bases match
                5). read quality is more than 10
            */
            while *cigar_len >= 1
                && self.start > 1
                && self.start - 1 <= *instance().chr_lens.get(contig).unwrap() as i64
                && is_has_and_equals(
                    query_sequence
                        .get_or_err(*cigar_len as usize - 1)
                        .copied()?,
                    &self.reference.ref_seq,
                    (self.start - 1) as usize,
                )
                && query_quality.get_or_err(*cigar_len as usize - 1)? - 33 > 10
            {
                //create variant if it is not present
                let ref_b = self
                    .reference
                    .ref_seq
                    .get_or_err(self.start as usize - 1)
                    .copied()?;

                let var: &mut Variant = get_variants_from_map(
                    &mut self.non_insertion_vars,
                    self.start,
                    &VarDesc::SNV { ref_base: ref_b },
                );
                //add count
                add_cnt(
                    var,
                    is_reverse,
                    *cigar_len as usize,
                    query_quality.get_or_err(*cigar_len as usize - 1).copied()? - 33,
                    mapq,
                    num_mismatch,
                );
                //increase coverage
                inc_cnt(&mut self.ref_coverage, self.start - 1, 1);

                self.start -= 1;
                *cigar_len -= 1;
            }

            if *cigar_len > 0 {
                //If there remains a soft-clipped sequence at the beginning (not everything was matched)
                let mut read_qual_sum = 0_usize;
                let mut num_high_qual_base = 0;
                let mut num_low_qual_base = 0;

                // Loop over remaining soft-clipped sequence
                for si in (0..*cigar_len as usize).rev() {
                    // Stop if unknown base (N - any of ATGC) is found
                    if query_sequence.get_or_err(si).copied()? == b'N' {
                        break;
                    }

                    // base quality
                    let bq = query_quality.get_or_err(si)? - 33;
                    if bq <= 12 {
                        num_low_qual_base += 1;
                        break; //Stop if a low-quality base is found
                    }

                    read_qual_sum += bq as usize;
                    num_high_qual_base += 1;
                }

                {
                    let cigar_len = self.cigar_len;
                    self.sclip5_high_quality_processing(
                        query_sequence,
                        mapq,
                        query_quality,
                        num_mismatch,
                        is_reverse,
                        read_qual_sum,
                        num_high_qual_base,
                        num_low_qual_base,
                        cigar_len,
                    )?;
                }

                cigar_len = &mut self.cigar_len;
            }

            *cigar_len = cigar.get(ci).unwrap().len();
        } else if ci == cigar.len() - 1 {
            // 3' soft clip
            // Ignore large soft clip due to chimeric reads in library construction
            if !instance().conf.chimeric_filter {
                if *cigar_len >= 20
                    && let Some(sa_tag_val) = record.aux_option(b"SA")?
                {
                    if is_read_chimeric_with_sa(
                        record,
                        pos,
                        sa_tag_val.try_get_str()?,
                        is_reverse,
                        false,
                        max_read_len,
                        cigar,
                    )? {
                        self.read_pos_including_softclip += *cigar_len as usize;
                        self.offset = 0;

                        // Had to reset the start due to softclipping adjustment
                        self.start = pos;

                        return Ok(());
                    }
                } else if *cigar_len >= Configuration::SEED_1 as u32 {
                    let ref_seed_map = &self.reference.seed;
                    let rev_comp_seq = record
                        .seq()
                        .into_decoded_base_iter()
                        .rev()
                        .take(*cigar_len as usize)
                        .map(complement_base)
                        // .take(Configuration::SEED_1 as usize)
                        .collect::<Vec<_>>();

                    // TODO: Start from here.
                    if let Some(poss) = ref_seed_map
                        .get(rev_comp_seq.get_with_int(-(Configuration::SEED_1 as i32)..)?)
                    {
                        if poss.len() == 1
                            && self.start - poss.get(0).unwrap() < 2 * self.max_read_len as i64
                        {
                            self.read_pos_including_softclip += *cigar_len as usize;
                            self.offset = 0;
                            // Had to reset the start due to softclipping adjustment
                            self.start = pos;

                            event!(
                                Level::INFO,
                                "{} at 3' is a chimeric at {} by SEED {}",
                                rev_comp_seq.as_slice().try_as_str()?,
                                self.start,
                                Configuration::SEED_1,
                            );
                            return Ok(());
                        }
                    }
                }
            }

            /*
            Conditions:
            1). read position is less than sequence length
            2). reference base is defined for start
            3). reference base at start matches read base at n
            4). read quality is more than 10
             */

            while self.read_pos_including_softclip < query_sequence.len()
                && is_has_and_equals(
                    query_sequence
                        .get_or_err(self.read_pos_including_softclip)
                        .copied()?,
                    &self.reference.ref_seq,
                    (self.start) as usize,
                )
                && query_quality.get_or_err(self.read_pos_including_softclip)? - 33 > 10
            {
                //Initialize entry if not present
                let ref_b = self
                    .reference
                    .ref_seq
                    .get_or_err(self.start as usize - 1)
                    .copied()?;

                let var: &mut Variant = get_variants_from_map(
                    &mut self.non_insertion_vars,
                    self.start,
                    &VarDesc::SNV { ref_base: ref_b },
                );
                //add count
                add_cnt(
                    var,
                    is_reverse,
                    total_length_including_soft_clipped - self.read_pos_excluding_softclip,
                    query_quality
                        .get_or_err(self.read_pos_including_softclip)
                        .copied()?
                        - 33,
                    mapq,
                    num_mismatch,
                );
                // Add coverage
                inc_cnt(&mut self.ref_coverage, self.start, 1);
                self.read_pos_including_softclip += 1;
                self.read_pos_excluding_softclip += 1;
                self.start += 1;
                *cigar_len -= 1;
            }

            // If there remains a soft-clipped sequence at the end (not everything was
            // matched)
            if query_sequence.len() - self.read_pos_including_softclip > 0 {
                let mut read_qual_sum = 0;
                let mut num_high_qual_base = 0;
                let mut num_low_qual_base = 0;
                for si in 0..*cigar_len {
                    // Loop over remaining soft-clipped sequence
                    // Stop if unknown base (N - any of ATGC) is found

                    // At this point, self.read_pos_including_softclip is start offset of soft clip. why?
                    if query_sequence
                        .get_or_err(self.read_pos_including_softclip + si as usize)
                        .copied()?
                        == b'N'
                    {
                        break;
                    }

                    let base_quality = query_quality
                        .get_or_err(self.read_pos_including_softclip + si as usize)?
                        - 33;

                    if base_quality <= 12 {
                        num_low_qual_base += 1;
                    }
                    // Stop if a low-quality base is found
                    if num_low_qual_base > 1 {
                        break;
                    }

                    read_qual_sum += base_quality as usize;
                    num_high_qual_base += 1;
                }

                {
                    let cigar_len = self.cigar_len;
                    self.sclip3_high_quality_processing(
                        query_sequence,
                        mapq,
                        query_quality,
                        num_mismatch,
                        is_reverse,
                        read_qual_sum,
                        num_high_qual_base,
                        num_low_qual_base,
                        cigar_len,
                    )?;
                }

                cigar_len = &mut self.cigar_len;
            }
        }

        // Move read position by m (length of segment in CIGAR)
        self.read_pos_including_softclip += *cigar_len as usize;
        self.offset = 0;
        self.start = pos; // Had to reset the start due to softclipping adjustment

        Ok(())
    }

    /// Process CIGAR deletions part. Will ignore indels next to introns and create
    /// Variations for deletions
    fn process_deletion(
        &mut self,
        query_sequence: &[u8],
        mapq: u8,
        contig_ref_seq: &[u8],
        query_quality: &[u8],
        nm: usize,
        is_reverse: bool,
        read_len_including_match_ins: usize,
        ci: usize,
    ) -> Result<usize, Error> {
        // Ignore deletions right after introns at exon edge in RNA-seq
        if skip_indel_next_to_intron(&self.cigar, ci)? {
            self.read_pos_excluding_softclip += self.cigar_len as usize;

            return Ok(ci);
        }

        // $s description string of deleted segment
        

        todo!()
    }

    /// N in CIGAR - skipped region from reference
    /// Skip the region and add string start-end to %SPLICE
    fn process_not_matched(&mut self, cigar_len: u32) {
        let key = (self.start - 1, self.start + cigar_len as i64 - 1);

        self.splice_count
            .entry(key)
            .and_modify(|v| v[0] += 1)
            .or_insert_with(|| vec![1]);

        self.start += cigar_len as i64;
        self.offset = 0;
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

    fn skip_sites_out_region_of_interest(&self, cigars: &[Cigar]) -> bool {
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
                return !(cut_site - self.start as i32 > filter_bp
                    && self.start as i32 + rlen3 as i32 - cut_site > filter_bp);
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

    fn contig_ref_seq(&self) -> &Vec<u8> {
        &self.reference.ref_seq
    }

    /// Process soft clip on 5' if it is has high quality reads
    fn sclip5_high_quality_processing(
        &mut self,
        query_sequence: &[u8],
        mapq: u8,
        query_quality: &[u8],
        num_mismatch: usize,
        is_reverse: bool,
        read_qual_sum: usize,
        num_high_qual_base: usize,
        num_low_qual_base: usize,
        cigar_len: u32,
    ) -> Result<(), Error> {
        // If we have at least 1 high-quality soft-clipped base of region of interest
        if num_high_qual_base >= 1
            && num_high_qual_base > num_low_qual_base
            && self.start as usize >= self.region.start
            && self.start as usize <= self.region.end
        {
            //add record to $sclip5
            let sclip = self
                .soft_clips5_end
                .entry(self.start)
                .or_insert_with(|| SoftClip::default());

            for si in (0..cigar_len).rev() {
                if !(cigar_len - si <= num_high_qual_base as u32) {
                    break;
                }

                let b = query_sequence.get_or_err(si as usize).copied()?;
                let idx = cigar_len - 1 - si; // distannce from start of match.
                let cnts = sclip
                    .nt
                    .entry(idx as i64)
                    .or_insert_with(|| NucBaseMap::default());

                // increase count of current base.
                *cnts.get_mut(b).unwrap() += 1;

                let seq_var = get_variation_from_seq(sclip, idx as usize, b);
                add_cnt(
                    seq_var,
                    is_reverse,
                    si as usize - (cigar_len as usize - num_high_qual_base),
                    *query_quality.get_or_err(si as usize)? - 33,
                    mapq,
                    num_mismatch,
                );
            }

            add_cnt(
                &mut sclip.var,
                is_reverse,
                cigar_len as usize,
                (read_qual_sum / num_high_qual_base) as u8,
                mapq,
                num_mismatch,
            );
        }

        Ok(())
    }

    /// Process soft clip on 3' if it is has high quality reads
    fn sclip3_high_quality_processing(
        &mut self,
        query_sequence: &[u8],
        mapq: u8,
        query_quality: &[u8],
        num_mismatch: usize,
        is_reverse: bool,
        read_qual_sum: usize,
        num_high_qual_base: usize,
        num_low_qual_base: usize,
        cigar_len: u32,
    ) -> Result<(), Error> {
        // If we have at least 1 high-quality soft-clipped base of region of interest
        if num_high_qual_base >= 1
            && num_high_qual_base > num_low_qual_base
            && self.start as usize >= self.region.start
            && self.start as usize <= self.region.end
        {
            //add record to $sclip5
            let sclip = self
                .soft_clips3_end
                .entry(self.start)
                .or_insert_with(|| SoftClip::default());

            for si in (0..num_high_qual_base) {
                let b = query_sequence
                    .get_or_err(self.read_pos_including_softclip + si as usize)
                    .copied()?;
                let idx = si; // distannce from start of match.
                let cnts = sclip
                    .nt
                    .entry(idx as i64)
                    .or_insert_with(|| NucBaseMap::default());

                // increase count of current base.
                *cnts.get_mut(b).unwrap() += 1;

                let seq_var = get_variation_from_seq(sclip, idx as usize, b);
                add_cnt(
                    seq_var,
                    is_reverse,
                    num_high_qual_base - si,
                    *query_quality.get_or_err(self.read_pos_including_softclip + si)? - 33,
                    mapq,
                    num_mismatch,
                );
            }

            add_cnt(
                &mut sclip.var,
                is_reverse,
                cigar_len as usize,
                (read_qual_sum / num_high_qual_base) as u8,
                mapq,
                num_mismatch,
            );
        }

        Ok(())
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

/// Check if read is chimeric and contains SA tag
fn is_read_chimeric_with_sa(
    record: &Record,
    pos: i64,
    sa_tag_val: &str,
    is_reverse: bool,
    is_5side: bool,
    max_read_len: usize,
    cigar: &CigarStringView,
) -> Result<bool, Error> {
    let mut sa_tag_split = sa_tag_val.split(",");

    let mut get_next = || {
        sa_tag_split
            .next()
            .ok_or_else(|| anyhow!("Failed to parse SA tag: {}", sa_tag_val))
    };

    let sa_chr = get_next()?;
    let sa_pos = get_next()?.parse::<i64>()?;
    let sa_dir = get_next()?;
    let sa_cigar = get_next()?;
    let sa_dir_is_forward = sa_dir == "+";

    let cap_opt = if is_5side {
        SA_CIGAR_D_S_5CLIP.captures(sa_cigar)
    } else {
        SA_CIGAR_D_S_3CLIP.captures(sa_cigar)
    };

    let is_chimeric_with_sa = (is_reverse && sa_dir_is_forward)
        || (!is_reverse && !sa_dir_is_forward)
            && sa_chr == record.contig() // TODO: Shouldn't two chromosomes are different?
            && (sa_pos - pos).abs() < 2 * max_read_len as i64
            && cap_opt.is_some();

    if is_chimeric_with_sa {
        event!(
            Level::INFO,
            "{} {} {} {} {} is ignored as chimeric with SA: {},{},{}",
            record.qname().try_as_str()?,
            record.contig(),
            pos,
            record.mapq(),
            cigar,
            sa_pos,
            sa_dir,
            sa_cigar,
        )
    }

    Ok(is_chimeric_with_sa)
}

/// Increment variant counters.
fn add_cnt(var: &mut Variant, is_reverse: bool, read_pos: usize, bq: u8, mapq: u8, nm: usize) {
    var.alt_depth += 1;
    var.inc_dir(is_reverse);
    var.mean_pos += read_pos as f64;
    var.mean_qual += bq as f64;
    var.mean_mapq += mapq as f64;
    var.nm += nm as f64;

    if bq as f64 >= instance().conf.goodq {
        var.high_qual_read_cnt += 1;
    } else {
        var.low_qual_read_cnt += 1;
    }
}

/// Increase count for given key
#[inline]
fn inc_cnt(coverage_map: &mut HashMap<i64, usize>, pos: i64, depth: usize) {
    coverage_map
        .entry(pos)
        .and_modify(|v| v.add_assign(depth))
        .or_insert_with(|| depth);
}

/// Skip the insertions and deletions that are right after or before introns
/// (they indicate of aligner problem)
fn skip_indel_next_to_intron(cigar: &CigarStringView, ci: usize) -> Result<bool, Error> {
    // if the current cigar is not the last item.
    // and the next cigar is ref skip
    // or
    // the current cigar is not the first one and the previous one is refskip
    if (ci > 0 && matches!(cigar.get_or_err(ci - 1)?, Cigar::RefSkip(_)))
        || (ci < cigar.len() - 1 && matches!(cigar.get_or_err(ci + 1)?, Cigar::RefSkip(_)))
    {
        Ok(true)
    } else {
        Ok(false)
    }
}
