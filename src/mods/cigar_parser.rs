use bio_types::genome::AbstractInterval;
use smallvec::SmallVec;
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
    prelude::SmallVecBytes,
    scopedata::global_read_only_scope::{GlobalReadOnlyScope, instance},
    utils::{BytesExt, SliceExt, SliceExt2, aligner::Aligner},
    variants::{
        var_utils::{get_variants_from_map, get_variation_from_seq, is_has_and_equals, is_has_and_not_equals},
        variants::{InsOrDelLen, SoftClip, VarDesc, Variant},
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
    insertion_vars: HashMap<i64, HashMap<VarDesc, Variant>>,

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
            insertion_vars: Default::default(),
            ref_coverage: Default::default(),
            soft_clips5_end: Default::default(),
            soft_clips3_end: Default::default(),
            rev_complementor: RevComplementor::new(),
            cigar: CigarString(vec![]).into_view(0),
        }
    }
}

impl CigarParser {
    /// Create a new CigarParser for processing a region
    pub fn new(
        region: Region,
        reference: Reference,
        instance: Arc<GlobalReadOnlyScope>,
    ) -> Self {
        Self {
            query_seq_buf: Some(Vec::with_capacity(512)),
            query_qual_buf: Some(Vec::with_capacity(512)),
            aligner: Aligner::default(),
            instance,
            reference,
            max_read_len: 0,
            region,
            discordant_count: 0,
            splice_count: HashMap::new(),
            read_pos_including_softclip: 0,
            read_pos_excluding_softclip: 0,
            start: 0,
            offset: 0,
            cigar_len: 0,
            non_insertion_vars: HashMap::new(),
            insertion_vars: HashMap::new(),
            ref_coverage: HashMap::new(),
            soft_clips5_end: HashMap::new(),
            soft_clips3_end: HashMap::new(),
            rev_complementor: RevComplementor::new(),
            cigar: CigarString(vec![]).into_view(0),
        }
    }

    /// Process multiple BAM records and accumulate variant data
    ///
    /// This is the main entry point for pipeline usage.
    pub fn process_records<'a, I>(&mut self, records: I) -> Result<(), Error>
    where
        I: Iterator<Item = &'a mut Record>,
    {
        for record in records {
            self.parse_cigar(record)?;
        }
        Ok(())
    }

    /// Get the collected non-insertion variants
    pub fn get_non_insertion_vars(&self) -> &HashMap<i64, HashMap<VarDesc, Variant>> {
        &self.non_insertion_vars
    }

    /// Take ownership of the collected non-insertion variants
    pub fn take_non_insertion_vars(&mut self) -> HashMap<i64, HashMap<VarDesc, Variant>> {
        std::mem::take(&mut self.non_insertion_vars)
    }

    /// Get the collected insertion variants
    pub fn get_insertion_vars(&self) -> &HashMap<i64, HashMap<VarDesc, Variant>> {
        &self.insertion_vars
    }

    /// Take ownership of the collected insertion variants
    pub fn take_insertion_vars(&mut self) -> HashMap<i64, HashMap<VarDesc, Variant>> {
        std::mem::take(&mut self.insertion_vars)
    }

    /// Get the reference coverage map
    pub fn get_ref_coverage(&self) -> &HashMap<i64, usize> {
        &self.ref_coverage
    }

    /// Take ownership of the reference coverage map
    pub fn take_ref_coverage(&mut self) -> HashMap<i64, usize> {
        std::mem::take(&mut self.ref_coverage)
    }

    /// Get the 5' soft clips
    pub fn get_soft_clips_5end(&self) -> &HashMap<i64, SoftClip> {
        &self.soft_clips5_end
    }

    /// Take ownership of the 5' soft clips
    pub fn take_soft_clips_5end(&mut self) -> HashMap<i64, SoftClip> {
        std::mem::take(&mut self.soft_clips5_end)
    }

    /// Get the 3' soft clips
    pub fn get_soft_clips_3end(&self) -> &HashMap<i64, SoftClip> {
        &self.soft_clips3_end
    }

    /// Take ownership of the 3' soft clips
    pub fn take_soft_clips_3end(&mut self) -> HashMap<i64, SoftClip> {
        std::mem::take(&mut self.soft_clips3_end)
    }

    /// Get the maximum read length encountered
    pub fn get_max_read_len(&self) -> usize {
        self.max_read_len
    }

    /// Get the discordant read count
    pub fn get_discordant_count(&self) -> usize {
        self.discordant_count
    }

    fn parse_cigar(&mut self, record: &mut Record) -> Result<(), Error> {
        event!(Level::DEBUG, "[parse_cigar] Starting for record at pos {}", record.pos());
        
        let mut query_seq_buf = self.query_seq_buf.take().unwrap();
        let mut query_qual_buf = self.query_qual_buf.take().unwrap();

        query_seq_buf.extend(record.seq().into_decoded_base_iter());
        query_qual_buf.extend(record.qual());

        let mut query_seq = query_seq_buf.as_slice();
        let mut query_qual = query_qual_buf.as_slice();

        let mapping_quality = record.mapq();

        record.cache_cigar_if_empty();
        let mut cigar = record.cigar();
        
        event!(Level::DEBUG, "[parse_cigar] CIGAR: {:?}", cigar);

        let ins_del_len = get_ins_del_len(&cigar);

        let tot_nm = match record.aux_option(self.aligner.nm_tag())? {
            Some(rust_htslib::bam::record::Aux::I32(nm)) => nm - ins_del_len as i32,
            Some(oth) => Err(anyhow!("Got non i32 type for NM tag: {:?}", oth))?,
            None => {
                if !cigar.is_empty() {
                    event!(Level::WARN, "No NM tag for mismatches: {:?}", record);
                }

                if record.is_unmapped() || cigar.is_empty() {
                    event!(Level::DEBUG, "[parse_cigar] Early return: unmapped or empty cigar");
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

            // Convert to 1-based position (BAM is 0-based, Java VarDict uses 1-based)
            pos = mc.align_start_pos + 1;
            cigar = CigarString(mc.cigar.into_iter().collect::<Vec<_>>())
                .into_view(mc.align_start_pos as i64 + 1);

            query_qual = mc.query_qual;
            query_seq = mc.query_seq;
        } else {
            // Convert to 1-based position (BAM is 0-based, Java VarDict uses 1-based)
            pos = record.pos() + 1;
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
        self.start = pos; // Initialize start position from record position
        
        // Store cigar in self so helper functions can access it
        self.cigar = cigar.clone();

        'process_cigar: {
            //Loop over CIGAR records
            let mut ci = 0;
            while ci < cigar.len() {
                let c = cigar.get(ci).copied().unwrap();
                if self.skip_overlapping_reads(record, adj_pos, pos, is_reverse, mpos) {
                    break;
                }

                self.cigar_len = c.len();
                //Letter from CIGAR
                event!(Level::DEBUG, "[parse_cigar] Processing CIGAR op: {:?} at ci={}", c, ci);
                
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
                        // Before processing insertion, emit a ref call at the position before the insertion
                        // Java emits a reference call at the anchor position (position - 1 where insertion occurs)
                        let anchor_pos = self.start - 1;
                        event!(Level::DEBUG, "[Ins] Adding ref call at anchor position {} (self.start={})", anchor_pos, self.start);
                        if anchor_pos >= self.region.start as i64 && anchor_pos <= self.region.end as i64 {
                            if let Some(ref_base) = self.reference.get(anchor_pos) {
                                event!(Level::DEBUG, "[Ins] ref_base at {} = {}", anchor_pos, ref_base as char);
                                let var_desc = VarDesc::SNV { ref_base };
                                let var = get_variants_from_map(&mut self.non_insertion_vars, anchor_pos, &var_desc);
                                // For insertion anchor: read_pos is 0 (first position), tp will be 1
                                // BAM quality is already Phred (not +33 adjusted), so use directly
                                let bq = query_qual.first().copied().unwrap_or(0);
                                // Use bq=0 for high quality tracking - insertion anchors don't count as high quality
                                // The base quality is still recorded for mean_qual calculation
                                add_cnt_anchor(var, is_reverse, 0, bq, mapping_quality, nm as usize, Some(read_match_ins_len));
                                inc_cnt(&mut self.ref_coverage, anchor_pos, 1);
                            }
                        }
                        ci = self.process_insertion(
                            query_seq,
                            mapping_quality,
                            query_qual,
                            nm as usize,
                            is_reverse,
                            pos,
                            read_match_ins_len,
                            ci,
                            &cigar,
                        )?;
                    }
                    Cigar::Del(_) => {
                        offset = 0;
                        ci = self.process_deletion(
                            query_seq,
                            mapping_quality,
                            query_qual,
                            nm as usize,
                            is_reverse,
                            read_match_ins_len,
                            ci,
                        )?;
                    }
                    _ => {
                        // Match case: now dealing with matching part
                        let mut nmoff = 0;
                        let mut moffset = 0;
                        
                        // Loop over bases of CIGAR segment
                        let mut i = offset;
                        while i < self.cigar_len as usize && self.read_pos_including_softclip < query_seq.len() {
                            // Flag to trim reads at opt_t_bases from start or end
                            let trim = self.is_trim_at_opt_t_bases(
                                is_reverse,
                                read_len_including_softclips,
                            );
                            
                            // variation string. Initialize to first base of the read sequence
                            let ch1 = *query_seq.get_or_err(self.read_pos_including_softclip)?;
                            let mut s = SmallVecBytes::new();
                            s.push(ch1);
                            let mut start_with_deletion = false;
                            
                            // Skip if base is unknown
                            if ch1 == b'N' {
                                self.start += 1;
                                self.read_pos_including_softclip += 1;
                                self.read_pos_excluding_softclip += 1;
                                i += 1;
                                continue;
                            }

                            // Sum of qualities for bases
                            let mut q = (*query_qual.get_or_err(self.read_pos_including_softclip)?) as f64;
                            // Number of bases for quality calculation
                            let mut qbases = 1;
                            // Number of bases in insertion for quality calculation
                            let mut qibases = 0;
                            // For more than one nucleotide mismatch
                            let mut ss = SmallVecBytes::new();

                            // Multi-base mismatch detection loop
                            let mut loop_iterations = 0;
                            
                            while (self.start + 1) >= self.region.start as i64
                                && (self.start + 1) <= self.region.end as i64
                                && (i + 1) < self.cigar_len as usize
                                && self.read_pos_including_softclip + 1 < query_seq.len()  // Bounds check for read
                                && q >= instance().conf.goodq
                                && self.reference.has_and_not_equals(
                                    self.start,
                                    *query_seq.get_or_err(self.read_pos_including_softclip)?,
                                )
                                && !self.reference.has_and_equals(self.start, b'N')
                            {
                                event!(Level::DEBUG, "[MNV loop] iter={} pos={} read_base={} ref_base={:?} ss_len={}", 
                                    loop_iterations, self.start, 
                                    *query_seq.get(self.read_pos_including_softclip).unwrap_or(&b'?') as char,
                                    self.reference.get(self.start).map(|b| b as char),
                                    ss.len());
                                loop_iterations += 1;
                                // Require higher quality for MNV
                                // BAM stores Phred quality directly (no +33 ASCII offset)
                                if *query_qual.get_or_err(self.read_pos_including_softclip + 1)?
                                    < instance().conf.goodq as u8 + 5
                                {
                                    break;
                                }
                                
                                // Break if base is unknown in the read
                                let nuc = *query_seq.get_or_err(self.read_pos_including_softclip + 1)?;
                                if nuc == b'N' {
                                    break;
                                }
                                
                                if self.reference.has_and_equals(self.start + 1, b'N') {
                                    break;
                                }

                                // Check if base doesn't match reference
                                let ref_base = match self.reference.get(self.start + 1) {
                                    Some(b) => b,
                                    None => break,
                                };
                                if nuc != ref_base {
                                    // Consecutive mismatch
                                    ss.push(nuc);
                                    q += (*query_qual.get_or_err(self.read_pos_including_softclip + 1)?) as f64;
                                    qbases += 1;
                                    self.read_pos_including_softclip += 1;
                                    self.read_pos_excluding_softclip += 1;
                                    i += 1;
                                    self.start += 1;
                                    nmoff += 1;
                                } else {
                                    // Look ahead for more mismatches
                                    let mut ssn = 0;
                                    for ssi in 1..=(instance().conf.vext as usize) {
                                        if i + 1 + ssi >= self.cigar_len as usize {
                                            break;
                                        }
                                        if self.read_pos_including_softclip + 1 + ssi < query_seq.len()
                                            && self.reference.has_and_not_equals(
                                                self.start + 1 + ssi as i64,
                                                *query_seq.get_or_err(self.read_pos_including_softclip + 1 + ssi)?,
                                            )
                                        {
                                            ssn = ssi + 1;
                                            break;
                                        }
                                    }
                                    
                                    if ssn == 0 {
                                        break;
                                    }
                                    
                                    // Require higher quality for MNV
                                    // BAM stores Phred quality directly (no +33 ASCII offset)
                                    if *query_qual.get_or_err(self.read_pos_including_softclip + ssn)?
                                        < instance().conf.goodq as u8 + 5
                                    {
                                        break;
                                    }
                                    
                                    for ssi in 1..=ssn {
                                        ss.push(*query_seq.get_or_err(self.read_pos_including_softclip + ssi)?);
                                        q += (*query_qual.get_or_err(self.read_pos_including_softclip + ssi)?) as f64;
                                        qbases += 1;
                                    }
                                    self.read_pos_including_softclip += ssn;
                                    self.read_pos_excluding_softclip += ssn;
                                    i += ssn;
                                    self.start += ssn as i64;
                                }
                            }

                            // If multi-base mismatch is found, append it to s
                            if !ss.is_empty() {
                                s.push(b'&');
                                s.extend_from_slice(&ss);
                            }
                            
                            let mut ddlen = 0;

                            // Check for adjacent deletions
                            if self.is_closer_then_vext_and_good_base(
                                query_seq,
                                query_qual,
                                ci,
                                i,
                                &ss,
                                Cigar::Del(0),
                            )? {
                                // Extend through end of current segment
                                while i + 1 < self.cigar_len as usize {
                                    s.push(*query_seq.get_or_err(self.read_pos_including_softclip + 1)?);
                                    q += (*query_qual.get_or_err(self.read_pos_including_softclip + 1)?) as f64;
                                    qbases += 1;

                                    i += 1;
                                    self.read_pos_including_softclip += 1;
                                    self.read_pos_excluding_softclip += 1;
                                    self.start += 1;
                                }

                                // Reformat s with deletion notation
                                let next_del_len = self.cigar.get(ci + 1).map(|c| c.len()).unwrap_or(0);
                                let mut new_s = SmallVecBytes::new();
                                new_s.extend_from_slice(b"-");
                                let del_str = next_del_len.to_string();
                                new_s.extend_from_slice(del_str.as_bytes());
                                new_s.push(b'&');
                                
                                // Remove '&' if present and rebuild
                                if let Some(amp_pos) = s.iter().position(|&b| b == b'&') {
                                    new_s.extend_from_slice(&s[0..amp_pos]);
                                    new_s.extend_from_slice(&s[amp_pos + 1..]);
                                } else {
                                    new_s.extend_from_slice(&s);
                                }
                                s = new_s;
                                
                                start_with_deletion = true;
                                ddlen = next_del_len as usize;
                                ci += 1;

                                // Check for insertion two segments ahead
                                if self.is_two_insertions_ahead(ci) {
                                    let ins_len = self.cigar.get(ci + 1).map(|c| c.len()).unwrap_or(0) as usize;
                                    s.push(b'^');
                                    s.extend_from_slice(
                                        query_seq.get_or_err(
                                            (self.read_pos_including_softclip + 1)
                                                ..(self.read_pos_including_softclip + 1 + ins_len),
                                        )?,
                                    );

                                    for _ in 1..=ins_len {
                                        q += (*query_qual.get_or_err(self.read_pos_including_softclip + 1)?) as f64;
                                        qibases += 1;
                                    }

                                    self.read_pos_including_softclip += ins_len;
                                    self.read_pos_excluding_softclip += ins_len;
                                    ci += 1;
                                }

                                if self.is_next_after_num_matched(ci, 1)? {
                                    if let Some((toffset, tnmoff, tseq, tqual)) = 
                                        self.find_offset(
                                            (self.start + ddlen as i64 + 1) as usize,
                                            self.read_pos_including_softclip + 1,
                                            self.cigar.get(ci + 1).map(|c| c.len()).unwrap_or(0) as usize,
                                            query_seq,
                                            query_qual,
                                        )? 
                                    {
                                        moffset = toffset;
                                        nmoff += tnmoff;
                                        s.push(b'&');
                                        s.extend_from_slice(&tseq);
                                        for &bq in &tqual {
                                            // BAM stores Phred quality directly
                                            q += bq as f64;
                                            qibases += 1;
                                        }
                                    }
                                }
                            } else if self.is_closer_then_vext_and_good_base(
                                query_seq,
                                query_qual,
                                ci,
                                i,
                                &ss,
                                Cigar::Ins(0),
                            )? {
                                while i + 1 < self.cigar_len as usize {
                                    s.push(*query_seq.get_or_err(self.read_pos_including_softclip + 1)?);
                                    q += (*query_qual.get_or_err(self.read_pos_including_softclip + 1)?) as f64;
                                    qbases += 1;

                                    i += 1;
                                    self.read_pos_including_softclip += 1;
                                    self.read_pos_excluding_softclip += 1;
                                    self.start += 1;
                                }
                                
                                let next_len = self.cigar.get(ci + 1).map(|c| c.len()).unwrap_or(0) as usize;
                                
                                // Remove '&' from s
                                let mut new_s = SmallVecBytes::new();
                                if let Some(amp_pos) = s.iter().position(|&b| b == b'&') {
                                    new_s.extend_from_slice(&s[0..amp_pos]);
                                    new_s.extend_from_slice(&s[amp_pos + 1..]);
                                } else {
                                    new_s.extend_from_slice(&s);
                                }
                                
                                let mut insertion_seq = SmallVecBytes::new();
                                insertion_seq.extend_from_slice(
                                    query_seq.get_or_err(
                                        (self.read_pos_including_softclip + 1)
                                            ..(self.read_pos_including_softclip + 1 + next_len),
                                    )?,
                                );
                                new_s.extend_from_slice(&insertion_seq);
                                
                                let mut final_s = SmallVecBytes::new();
                                final_s.extend_from_slice(&insertion_seq[0..next_len.min(new_s.len())]);
                                final_s.push(b'&');
                                if next_len < new_s.len() {
                                    final_s.extend_from_slice(&new_s[next_len..]);
                                }
                                final_s.insert(0, b'+');
                                s = final_s;

                                for _ in 1..=next_len {
                                    q += (*query_qual.get_or_err(self.read_pos_including_softclip + 1)?) as f64;
                                    qibases += 1;
                                }
                                
                                self.read_pos_including_softclip += next_len;
                                self.read_pos_excluding_softclip += next_len;
                                ci += 1;
                                qibases -= 1;
                                qbases += 1;
                            }

                            if !trim {
                                let pos = self.start - qbases as i64 + 1;
                                if pos >= self.region.start as i64
                                    && pos <= self.region.end as i64
                                    && !s.iter().any(|&b| b == b'N')
                                {
                                    event!(Level::DEBUG, "[parse_cigar] Calling add_variation_for_matching_part: pos={}, s={:?}", 
                                        pos, String::from_utf8_lossy(&s));
                                    
                                    self.add_variation_for_matching_part(
                                        mapping_quality,
                                        nm as usize,
                                        is_reverse,
                                        read_match_ins_len,
                                        self.read_pos_excluding_softclip,  // Use current position (after MNV loop, before final +1)
                                        nmoff,
                                        &s,
                                        start_with_deletion,
                                        q,
                                        qbases,
                                        qibases,
                                        ddlen,
                                        pos,
                                    );
                                } else {
                                    event!(Level::DEBUG, "[parse_cigar] Skipping variation: pos={} (region: {}-{}), has_N={}, trim={}", 
                                        pos, self.region.start, self.region.end, 
                                        s.iter().any(|&b| b == b'N'), trim);
                                }
                            }

                            // If variation starts with deletion
                            if start_with_deletion {
                                self.start += ddlen as i64;
                            }

                            // Shift reference position by 1 if CIGAR segment is not insertion
                            if !matches!(c, Cigar::Ins(_)) {
                                self.start += 1;
                            }
                            
                            // Shift read position by 1 if CIGAR segment is not deletion
                            if !matches!(c, Cigar::Del(_)) {
                                self.read_pos_including_softclip += 1;
                                self.read_pos_excluding_softclip += 1;
                            }

                            // Check for overlapping reads
                            if self.skip_overlapping_reads(record, adj_pos, pos, is_reverse, mpos) {
                                break 'process_cigar;
                            }

                            i += 1;
                        }

                        if moffset != 0 {
                            offset = moffset;
                            self.read_pos_including_softclip += moffset;
                            self.start += moffset as i64;
                            self.read_pos_excluding_softclip += moffset;
                        }
                    }
                }

                if self.start > self.region.end as i64 {
                    break;
                }
                
                ci += 1;
            }
        }

        // return buf resource to self.
        query_qual_buf.clear();
        query_seq_buf.clear();
        let _ = self.query_qual_buf.insert(query_qual_buf);
        let _ = self.query_seq_buf.insert(query_seq_buf);

        Ok(())
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
                && *query_quality.get_or_err(*cigar_len as usize - 1)? > 10
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
                // BAM stores Phred quality directly
                add_cnt(
                    var,
                    is_reverse,
                    *cigar_len as usize,
                    query_quality.get_or_err(*cigar_len as usize - 1).copied()?,
                    mapq,
                    num_mismatch,
                    None,
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

                    // base quality - BAM stores Phred quality directly
                    let bq = *query_quality.get_or_err(si)?;
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
                && *query_quality.get_or_err(self.read_pos_including_softclip)? > 10
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
                //add count - BAM stores Phred quality directly
                add_cnt(
                    var,
                    is_reverse,
                    total_length_including_soft_clipped - self.read_pos_excluding_softclip,
                    query_quality
                        .get_or_err(self.read_pos_including_softclip)
                        .copied()?,
                    mapq,
                    num_mismatch,
                    None,
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

                    // BAM stores Phred quality directly
                    let base_quality = *query_quality
                        .get_or_err(self.read_pos_including_softclip + si as usize)?;

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
        query_seq: &[u8],
        mapq: u8,
        query_qual: &[u8],
        nm: usize,
        is_reverse: bool,
        read_len_including_match_ins: usize,
        mut ci: usize,
    ) -> Result<usize, Error> {
        // Ignore deletions right after introns at exon edge in RNA-seq
        if skip_indel_next_to_intron(&self.cigar, ci)? {
            self.read_pos_excluding_softclip += self.cigar_len as usize;

            return Ok(ci);
        }

        let mut var_desc = VarDesc::Del {
            len: self.cigar_len,
            match_seq: SmallVec::new(),
            ins_or_del_len: InsOrDelLen::None,
            mismatch_seq: SmallVec::new(),
        };

        let mut seq_to_append_if_next_matched = SmallVecBytes::new();
        let quality_of_last_segment_before_del = query_qual
            .get_or_err(self.read_pos_including_softclip - 1)
            .copied()?;

        let mut qual_seg = SmallVecBytes::new();

        // For multiple indels within $VEXT bp
        // offset for reference position if next segment is matched
        let mut ref_offset = 0;
        // offset for read position if next segment is matched
        let mut read_offset = 0;
        let mut nmoff = 0;

        /*
         * Condition:
         * 1). CIGAR string has next entry
         * 2). length of next CIGAR segment is less than conf.vext
         * 3). next segment is matched
         * 4). CIGAR string has one more entry after next one
         * 5). this entry is insertion or deletion
         */

        if is_followed_by_match_and_indel(&self.cigar, ci) {
            // `is_followed_by_match_and_indel` proves that there are at least 2 cigar elems next to `ci`.
            let n_cigar = self.cigar.get(ci + 1).unwrap();
            let nn_cigar = self.cigar.get(ci + 2).unwrap();

            let mlen = n_cigar.len() as usize;
            let indel_len = nn_cigar.len() as usize;
            let begin = self.read_pos_including_softclip;

            append_segments(
                query_seq,
                query_qual,
                &self.cigar,
                ci,
                &mut var_desc,
                &mut qual_seg,
                begin,
                mlen,
                indel_len,
                false,
            );

            // add length of next segment to both read and reference offsets
            // add length of next-next segment to reference position (for insertion) or to
            // read position(for deletion)
            ref_offset += mlen
                + if matches!(nn_cigar, Cigar::Del(_)) {
                    indel_len
                } else {
                    0
                };
            read_offset += mlen
                + if matches!(nn_cigar, Cigar::Ins(_)) {
                    indel_len
                } else {
                    0
                };

            if is_next_after_num_matched(&self.cigar, ci, 3) {
                let mut vsn = 0;
                let tn = self.read_pos_including_softclip + read_offset;
                let ts = self.start as usize + ref_offset + self.cigar_len as usize;

                let nnn_cigar_len = self.cigar.get(ci + 3).unwrap().len();

                let mut vi = 0;
                while vsn <= instance().conf.vext as usize && vi < nnn_cigar_len as usize {
                    let b = *query_seq.get_or_err(tn + vi)?;
                    // BAM stores Phred quality directly
                    let q = *query_qual.get_or_err(tn + vi)?;

                    if b == b'N' {
                        break;
                    }

                    if (q as f64) < instance().conf.goodq {
                        break;
                    }

                    if self.reference.has_and_equals((ts + vi) as i64, b'N') {
                        break;
                    }

                    match self.reference.get((ts + vi) as i64) {
                        Some(ref_b) => {
                            if b != ref_b {
                                self.offset = vi + 1;
                                nmoff += 1;
                                vsn = 0;
                            } else {
                                vsn += 1;
                            }
                        }
                        None => {}
                    };

                    vi += 1;
                }

                if self.offset != 0 {
                    seq_to_append_if_next_matched
                        .extend_from_slice(query_seq.get_or_err(tn..(tn + self.offset))?);
                    qual_seg.extend_from_slice(query_qual.get_or_err(tn..(tn + self.offset))?);
                }
            }
            // skip next 2 CIGAR segments
            ci += 2;
        } else if is_next_ins(&self.cigar, ci) {
            /*
             * Condition:
             * 1). CIGAR string has next entry
             * 2). next CIGAR segment is an insertion
             */
            let ins_len = self.cigar.get(ci + 1).unwrap().len() as usize;

            // Append '^' + next segment (inserted) to var_desc
            if let VarDesc::Del { ins_or_del_len, .. } = &mut var_desc {
                *ins_or_del_len = InsOrDelLen::InsSeq(
                    query_seq
                        .get_or_err(
                            self.read_pos_including_softclip
                                ..self.read_pos_including_softclip + ins_len,
                        )?
                        .into(),
                );
            }

            // Append next segment to quality string
            qual_seg.extend_from_slice(
                query_qual.get_or_err(
                    self.read_pos_including_softclip..self.read_pos_including_softclip + ins_len,
                )?,
            );

            // Shift read position by length of insertion
            read_offset += ins_len;

            if is_next_after_num_matched(&self.cigar, ci, 2) {
                let mlen = self.cigar.get(ci + 2).unwrap().len() as usize;
                let mut vsn = 0;
                let tn = self.read_pos_including_softclip + read_offset;
                let ts = self.start as usize + self.cigar_len as usize;

                let mut vi = 0;
                while vsn <= instance().conf.vext as usize && vi < mlen {
                    let seq_ch = *query_seq.get_or_err(tn + vi)?;
                    if seq_ch == b'N' {
                        break;
                    }
                    if ((*query_qual.get_or_err(tn + vi)?) as f64)
                        < instance().conf.goodq
                    {
                        break;
                    }
                    match self.reference.get((ts + vi) as i64) {
                        Some(ref_ch) => {
                            if ref_ch == b'N' {
                                break;
                            }
                            if seq_ch != ref_ch {
                                self.offset = vi + 1;
                                nmoff += 1;
                                vsn = 0;
                            } else {
                                vsn += 1;
                            }
                        }
                        None => {}
                    }
                    vi += 1;
                }

                if self.offset != 0 {
                    seq_to_append_if_next_matched
                        .extend_from_slice(query_seq.get_or_err(tn..tn + self.offset)?);
                    qual_seg.extend_from_slice(query_qual.get_or_err(tn..tn + self.offset)?);
                }
            }
            // skip next CIGAR segment
            ci += 1;
        } else if is_next_matched(&self.cigar, ci) {
            /*
             * Condition:
             * 1). CIGAR string has next entry
             * 2). next CIGAR segment is matched
             */
            let mlen = self.cigar.get(ci + 1).unwrap().len() as usize;
            let mut vsn = 0;

            // Loop over next CIGAR segment (no more than conf.vext bases ahead)
            let mut vi = 0;
            while vsn <= instance().conf.vext as usize && vi < mlen {
                let seq_ch = *query_seq.get_or_err(self.read_pos_including_softclip + vi)?;
                // If base is unknown, exit loop
                if seq_ch == b'N' {
                    break;
                }
                // If base quality is less than GOODQ, exit loop
                if ((*query_qual.get_or_err(self.read_pos_including_softclip + vi)?) as f64)
                    < instance().conf.goodq
                {
                    break;
                }
                // If reference sequence has base at this position
                let ref_pos = self.start + self.cigar_len as i64 + vi as i64;
                match self.reference.get(ref_pos) {
                    Some(ref_ch) => {
                        if ref_ch == b'N' {
                            break;
                        }
                        if seq_ch != ref_ch {
                            self.offset = vi + 1;
                            nmoff += 1;
                            vsn = 0;
                        } else {
                            vsn += 1;
                        }
                    }
                    None => {}
                }
                vi += 1;
            }

            // If next CIGAR segment has good matching base
            if self.offset != 0 {
                // Append first offset bases of next segment to ss and q
                seq_to_append_if_next_matched.extend_from_slice(
                    query_seq.get_or_err(
                        self.read_pos_including_softclip
                            ..self.read_pos_including_softclip + self.offset,
                    )?,
                );
                qual_seg.extend_from_slice(
                    query_qual.get_or_err(
                        self.read_pos_including_softclip
                            ..self.read_pos_including_softclip + self.offset,
                    )?,
                );
            }
        }

        // Append '&' and seq_to_append_if_next_matched to var_desc if offset > 0
        if self.offset > 0 {
            if let VarDesc::Del { mismatch_seq, .. } = &mut var_desc {
                mismatch_seq.extend_from_slice(&seq_to_append_if_next_matched);
            }
        }

        // Quality of first matched base after deletion
        // Append best of q1 and q2
        if self.read_pos_including_softclip + self.offset >= query_qual.len() {
            qual_seg.push(quality_of_last_segment_before_del);
        } else {
            let quality_after_del =
                *query_qual.get_or_err(self.read_pos_including_softclip + self.offset)?;
            qual_seg.push(quality_of_last_segment_before_del.max(quality_after_del));
        }

        // If reference position is inside region of interest
        if self.start >= self.region.start as i64 && self.start <= self.region.end as i64 {
            self.add_variation_for_deletion(
                mapq,
                nm,
                is_reverse,
                read_len_including_match_ins,
                &var_desc,
                &qual_seg,
                nmoff,
            );
        }

        // Adjust reference position by cigar_len + offset + ref_offset
        self.start += self.cigar_len as i64 + self.offset as i64 + ref_offset as i64;

        // Adjust read position by offset + read_offset
        self.read_pos_including_softclip += self.offset + read_offset;
        self.read_pos_excluding_softclip += self.offset + read_offset;

        Ok(ci)
    }

    /// Add variation record for a deletion
    fn add_variation_for_deletion(
        &mut self,
        mapq: u8,
        nm: usize,
        is_reverse: bool,
        read_len_including_match_ins: usize,
        var_desc: &VarDesc,
        qual_seg: &[u8],
        nmoff: usize,
    ) {
        // Get or create variation structure for this deletion
        let var = get_variants_from_map(&mut self.non_insertion_vars, self.start, var_desc);

        // Increment direction count
        var.inc_dir(is_reverse);

        // Increase variant count
        var.alt_depth += 1;

        // Minimum of positions from start of read and end of read
        // Java: tp = n < rlen1 - n ? n + 1 : rlen1 - n
        let from_start = self.read_pos_excluding_softclip;
        let from_end = read_len_including_match_ins.saturating_sub(self.read_pos_excluding_softclip);
        let tp = if from_start < from_end {
            from_start + 1
        } else {
            from_end
        };

        // Average quality of bases in quality segment - BAM stores Phred quality directly
        let tmpq: f64 = if qual_seg.is_empty() {
            0.0
        } else {
            qual_seg.iter().map(|&q| q as f64).sum::<f64>() / qual_seg.len() as f64
        };

        // pstd: true if variant is covered by reads with different positions
        if !var.pstd && var.pp != 0 && tp != var.pp {
            var.pstd = true;
        }

        // qstd: true if variant is covered by reads with different qualities
        if !var.qstd && var.pq != 0.0 && (tmpq - var.pq).abs() > f64::EPSILON {
            var.qstd = true;
        }

        var.mean_pos += tp as f64;
        var.mean_qual += tmpq;
        var.mean_mapq += mapq as f64;
        var.pp = tp;
        var.pq = tmpq;
        var.nm += nm.saturating_sub(nmoff) as f64;

        if tmpq >= instance().conf.goodq {
            var.high_qual_read_cnt += 1;
        } else {
            var.low_qual_read_cnt += 1;
        }

        // Increase coverage count for reference bases missing from the read
        for i in 0..self.cigar_len as i64 {
            inc_cnt(&mut self.ref_coverage, self.start + i, 1);
        }
    }

    /// Check if position should be trimmed at opt_T bases
    fn is_trim_at_opt_t_bases(&self, is_reverse: bool, total_length_including_softclipped: usize) -> bool {
        // Note: trim_bases_after not available in Rust config yet
        false
    }

    /// Check if current position is closer than vext and has good base
    /// Checks for adjacent deletions or insertions
    fn is_closer_then_vext_and_good_base(
        &self,
        query_seq: &[u8],
        query_qual: &[u8],
        ci: usize,
        i: usize,
        ss: &[u8],
        cigar_type: Cigar,
    ) -> Result<bool, Error> {
        // Condition:
        // 1). index i is no farther than vext from end of CIGAR segment
        // 2). CIGAR string contains next entry
        // 3). next CIGAR entry is a deletion or insertion (based on cigar_type)
        // 4). reference sequence contains a base at start
        // 5). either multi-nucleotide mismatch is found or read base doesn't match reference base
        // 6). read base has good quality

        let is_del = matches!(cigar_type, Cigar::Del(_));
        let is_ins = matches!(cigar_type, Cigar::Ins(_));

        if (self.cigar_len as usize) - i <= instance().conf.vext as usize
            && self.cigar.len() > ci + 1
        {
            let next_cigar = self.cigar.get(ci + 1);
            let has_correct_type = if is_del {
                next_cigar.map(|c| matches!(c, Cigar::Del(_))).unwrap_or(false)
            } else if is_ins {
                next_cigar.map(|c| matches!(c, Cigar::Ins(_))).unwrap_or(false)
            } else {
                false
            };

            if has_correct_type && self.reference.get(self.start).is_some() {
                // Check if there's a mismatch or multi-base mismatch
                let has_mismatch = !ss.is_empty()
                    || self.reference.has_and_not_equals(
                        self.start,
                        *query_seq.get_or_err(self.read_pos_including_softclip)?,
                    );

                // BAM stores Phred quality directly
                let has_good_quality = *query_qual.get_or_err(self.read_pos_including_softclip)?
                    >= instance().conf.goodq as u8;

                return Ok(has_mismatch && has_good_quality);
            }
        }
        Ok(false)
    }

    /// Find mismatches in next segment with good quality
    fn find_offset(
        &self,
        ts: usize,
        tn: usize,
        seg_len: usize,
        query_seq: &[u8],
        query_qual: &[u8],
    ) -> Result<Option<(usize, usize, SmallVecBytes, SmallVecBytes)>, Error> {
        let mut offset = 0;
        let mut nmoff = 0;
        let mut seq_to_append = SmallVecBytes::new();
        let mut qual_to_append = SmallVecBytes::new();

        let mut vsn = 0;
        let mut vi = 0;

        while vsn <= instance().conf.vext as usize && vi < seg_len {
            let b = *query_seq.get_or_err(tn + vi)?;
            // BAM stores Phred quality directly
            let q = *query_qual.get_or_err(tn + vi)?;

            if b == b'N' {
                break;
            }

            if (q as f64) < instance().conf.goodq {
                break;
            }

            match self.reference.get((ts + vi) as i64) {
                Some(ref_b) => {
                    if b != ref_b {
                        offset = vi + 1;
                        nmoff += 1;
                        vsn = 0;
                    } else {
                        vsn += 1;
                    }
                }
                None => break,
            };

            vi += 1;
        }

        if offset != 0 {
            seq_to_append.extend_from_slice(query_seq.get_or_err(tn..(tn + offset))?);
            qual_to_append.extend_from_slice(query_qual.get_or_err(tn..(tn + offset))?);
            Ok(Some((offset, nmoff, seq_to_append, qual_to_append)))
        } else {
            Ok(None)
        }
    }

    /// Check if there are two insertions ahead
    fn is_two_insertions_ahead(&self, ci: usize) -> bool {
        self.cigar.len() > ci + 1 && matches!(self.cigar.get(ci + 1), Some(Cigar::Ins(_)))
    }

    /// Check if there's a matched segment after num positions
    fn is_next_after_num_matched(&self, ci: usize, num: usize) -> Result<bool, Error> {
        Ok(self.cigar.len() > ci + num
            && matches!(self.cigar.get(ci + num), Some(Cigar::Match(_))))
    }

    /// Process insertion in CIGAR
    /// 
    /// Java equivalent: processInsertion() in CigarParser.java
    /// 
    /// Key logic:
    /// 1. Get the inserted bases from the read
    /// 2. If next CIGAR is Match, look for mismatches using vext algorithm
    /// 3. Combine insertion + following mismatches into complex variant
    /// 4. Store variant at position `start - 1` (position before insertion)
    fn process_insertion(
        &mut self,
        query_seq: &[u8],
        mapq: u8,
        query_qual: &[u8],
        nm: usize,
        is_reverse: bool,
        _pos: i64,
        read_len_including_match_ins: usize,
        mut ci: usize,
        cigar: &CigarStringView,
    ) -> Result<usize, Error> {
        let ins_len = self.cigar_len as usize;
        
        event!(Level::DEBUG, "[process_insertion] Called: ci={} ins_len={} start={} read_pos={} cigar_len={}", 
            ci, ins_len, self.start, self.read_pos_including_softclip, cigar.len());
        
        for i in 0..cigar.len() {
            event!(Level::DEBUG, "[process_insertion] cigar[{}] = {:?}", i, cigar.get(i));
        }
        
        // Get the inserted sequence from the read
        let ins_seq = query_seq.get_or_err(
            self.read_pos_including_softclip..self.read_pos_including_softclip + ins_len
        )?;
        
        // Get quality of inserted segment
        let ins_qual = query_qual.get_or_err(
            self.read_pos_including_softclip..self.read_pos_including_softclip + ins_len
        )?;
        
        // Build the variant description string
        let mut desc: SmallVec<[u8; 64]> = SmallVec::new();
        desc.extend_from_slice(ins_seq);
        
        // Build quality string
        let mut qual_seg: Vec<u8> = ins_qual.to_vec();
        
        // Sequence to append if next segment has mismatches
        let mut ss: Vec<u8> = Vec::new();
        self.offset = 0;
        let mut nmoff = 0;
        
        // Check if next CIGAR segment is Match - look for mismatches to combine
        // Note: We always check next match for insertions, regardless of perform_local_realignment
        let next_is_match = cigar.get(ci + 1)
            .map_or(false, |c| matches!(c, Cigar::Match(_)));
        
        event!(Level::DEBUG, "[process_insertion] next_is_match={} cigar[ci+1]={:?}", 
            next_is_match, cigar.get(ci + 1));
        
        if next_is_match {
            let mlen = cigar.get(ci + 1).unwrap().len() as usize;
            let mut vsn = 0;
            
            // Position in read after insertion
            let tn = self.read_pos_including_softclip + ins_len;
            // Reference position (insertion doesn't consume reference)
            let ts = self.start;
            
            event!(Level::DEBUG, "[process_insertion] Looking for mismatches: ins_len={} mlen={} tn={} ts={}", 
                ins_len, mlen, tn, ts);
            
            // Loop over next CIGAR segment looking for mismatches (Java: findOffset-like logic)
            let mut vi = 0;
            while vsn <= instance().conf.vext as usize && vi < mlen {
                let seq_ch = *query_seq.get_or_err(tn + vi)?;
                
                // If base is unknown, exit loop
                if seq_ch == b'N' {
                    break;
                }
                
                // If base quality is less than GOODQ, exit loop - BAM stores Phred quality directly
                let qual = *query_qual.get_or_err(tn + vi)?;
                if (qual as f64) < instance().conf.goodq {
                    event!(Level::DEBUG, "[process_insertion] vi={} qual={} < goodq={}, breaking", vi, qual, instance().conf.goodq);
                    break;
                }
                
                // Check reference base
                match self.reference.get(ts + vi as i64) {
                    Some(ref_ch) => {
                        event!(Level::DEBUG, "[process_insertion] vi={} read={} ref={} at pos {}", 
                            vi, seq_ch as char, ref_ch as char, ts + vi as i64);
                        
                        if ref_ch == b'N' {
                            break;
                        }
                        if seq_ch != ref_ch {
                            // Mismatch found
                            self.offset = vi + 1;
                            nmoff += 1;
                            vsn = 0;
                        } else {
                            vsn += 1;
                        }
                    }
                    None => {
                        event!(Level::DEBUG, "[process_insertion] vi={} ref not found at pos {}", vi, ts + vi as i64);
                    }
                }
                vi += 1;
            }
            
            event!(Level::DEBUG, "[process_insertion] After loop: offset={} ss_len={} mlen={}", self.offset, ss.len(), mlen);
            
            // If mismatches found in following Match segment, append them
            if self.offset > 0 {
                ss.extend_from_slice(query_seq.get_or_err(tn..tn + self.offset)?);
                qual_seg.extend_from_slice(query_qual.get_or_err(tn..tn + self.offset)?);
                
                // Increase coverage for positions corresponding to first offset bases
                for osi in 0..self.offset {
                    inc_cnt(&mut self.ref_coverage, self.start + osi as i64, 1);
                }
                
                // If we consumed the entire next Match segment, skip it
                if self.offset == mlen {
                    ci += 1;
                    event!(Level::DEBUG, "[process_insertion] Consumed entire Match({}) segment, advancing ci to {}", mlen, ci);
                }
            }
        }
        
        // Skip if contains N
        let full_desc: Vec<u8> = if self.offset > 0 && !ss.is_empty() {
            // Combine insertion + following mismatches
            let mut combined = desc.to_vec();
            combined.extend_from_slice(&ss);
            combined
        } else {
            desc.to_vec()
        };
        
        if full_desc.iter().any(|&b| b == b'N') {
            self.read_pos_including_softclip += ins_len;
            self.read_pos_excluding_softclip += ins_len;
            return Ok(ci);
        }
        
        // Determine variant position and type based on whether mismatches follow
        // If mismatches found in following Match, create Complex at Match position
        // Otherwise, create simple Insertion at position before
        let (variant_pos, var_desc) = if self.offset > 0 && !ss.is_empty() {
            // Complex variant: insertion + following match bases
            // Position is the start of the Match segment (self.start)
            let var_pos = self.start;
            
            // Build reference sequence (the bases the Match covers)
            // Reference length = offset (number of bases consumed from Match)
            let mut ref_seq: SmallVec<[u8; 32]> = SmallVec::new();
            for i in 0..self.offset {
                if let Some(ref_base) = self.reference.get(self.start + i as i64) {
                    ref_seq.push(ref_base);
                }
            }
            
            // Alt sequence = inserted bases + match bases
            let alt_seq: SmallVec<[u8; 32]> = full_desc.iter().copied().collect();
            
            (var_pos, VarDesc::Complex { ref_seq, alt_seq })
        } else {
            // Simple insertion at position before (start - 1)
            let insertion_pos = self.start - 1;
            (insertion_pos, VarDesc::Ins { seq: full_desc.into() })
        };
        
        // Check if variant position is within region of interest
        if variant_pos >= self.region.start as i64 && variant_pos <= self.region.end as i64 {
            // Average quality - BAM already stores Phred quality directly (no need for -33)
            let avg_qual = if !qual_seg.is_empty() {
                let sum: u32 = qual_seg.iter().map(|&q| q as u32).sum();
                (sum as f64 / qual_seg.len() as f64) as u8
            } else {
                0
            };
            
            // Store the variant
            // Complex variants go to non_insertion_vars, simple insertions to insertion_vars
            let var = match &var_desc {
                VarDesc::Complex { .. } => {
                    get_variants_from_map(&mut self.non_insertion_vars, variant_pos, &var_desc)
                }
                _ => {
                    get_variants_from_map(&mut self.insertion_vars, variant_pos, &var_desc)
                }
            };
            
            // Use actual read position for mean_pos calculation
            // For insertions/complex variants, this is where the variant starts in the read
            let read_pos = self.read_pos_including_softclip;
            add_cnt(var, is_reverse, read_pos, avg_qual, mapq, nm.saturating_sub(nmoff), Some(read_len_including_match_ins));
            
            // Note: ref_coverage for Complex variants is already updated in the offset loop above (line ~1765)
            // Do NOT increment here - that would cause double-counting
        }
        
        // Advance read position by insertion length + any offset from following mismatches
        self.read_pos_including_softclip += ins_len + self.offset;
        self.read_pos_excluding_softclip += ins_len + self.offset;
        
        // Advance reference position by offset (mismatches consumed from following Match)
        self.start += self.offset as i64;
        
        Ok(ci)
    }

    /// Add variation record for a matched segment (SNV or MNV)
    fn add_variation_for_matching_part(
        &mut self,
        mapq: u8,
        nm: usize,
        is_reverse: bool,
        read_len_including_match_ins: usize,
        read_pos: usize,  // Current read position (0-based)
        nmoff: usize,
        s: &[u8],
        start_with_deletion: bool,
        q: f64,
        qbases: usize,
        qibases: usize,
        ddlen: usize,
        pos: i64,
    ) {
        // Build VarDesc from s string
        // s format examples: "A", "A&TGC" (MNV), "+ATC" (insertion), "-2&AT" (deletion+match), etc.
        
        if s.is_empty() || s.iter().all(|&b| b == b'N') {
            return;
        }

        // Average quality
        let avg_qual = if qbases + qibases > 0 {
            (q / (qbases + qibases) as f64) as u8
        } else {
            q as u8
        };
        
        // Check if this is a complex variant (contains & for MNV)
        if let Some(amp_pos) = s.iter().position(|&b| b == b'&') {
            // MNV format: first_base&rest_of_mnv
            // e.g., "A&TGC" means 4-base MNV: A, T, G, C at consecutive positions
            let first_base = s[0];
            let rest = &s[amp_pos + 1..]; // TGC
            
            // Build the full alt sequence
            let mut alt_seq: SmallVec<[u8; 32]> = SmallVec::new();
            alt_seq.push(first_base);
            alt_seq.extend_from_slice(rest);
            
            let mnv_len = alt_seq.len();
            
            // Build reference sequence for this MNV using coordinate-translated access
            let mut ref_seq: SmallVec<[u8; 32]> = SmallVec::new();
            for offset in 0..mnv_len {
                if let Some(ref_base) = self.reference.get(pos + offset as i64) {
                    ref_seq.push(ref_base);
                } else {
                    // Reference doesn't cover this position
                    event!(Level::DEBUG, "[add_variation_for_matching_part] MNV at pos {} extends beyond reference", pos);
                    return;
                }
            }
            
            // Create Complex/MNV VarDesc
            let var_desc = VarDesc::Complex {
                ref_seq,
                alt_seq,
            };
            
            // Get or create variant and add count - use read_pos for tp calculation
            let var = get_variants_from_map(&mut self.non_insertion_vars, pos, &var_desc);
            add_cnt(var, is_reverse, read_pos, avg_qual, mapq, nm, Some(read_len_including_match_ins));
            
            // Update reference coverage for all positions in MNV
            for offset in 0..mnv_len {
                inc_cnt(&mut self.ref_coverage, pos + offset as i64, 1);
            }
        } else if s.len() == 1 {
            // Simple base: single base (could be match or SNV)
            let alt_base = s[0];
            
            // Get reference base at this position using coordinate-translated access
            let ref_base = match self.reference.get(pos) {
                Some(b) => b,
                None => {
                    event!(Level::DEBUG, "[add_variation_for_matching_part] SNV at pos {} not in reference", pos);
                    return;
                }
            };
            
            // Create VarDesc keyed by alt_base (the read base)
            // This stores both matches (alt == ref) and mismatches (alt != ref)
            let var_desc = VarDesc::SNV { ref_base: alt_base };
            
            // Get or create variant and add count - use read_pos for tp calculation
            let var = get_variants_from_map(&mut self.non_insertion_vars, pos, &var_desc);
            add_cnt(var, is_reverse, read_pos, avg_qual, mapq, nm, Some(read_len_including_match_ins));
            
            // Update reference coverage
            inc_cnt(&mut self.ref_coverage, pos, 1);
        } else if s.starts_with(b"+") {
            // Insertion: +ATC format
            let ins_seq: SmallVec<[u8; 32]> = s[1..].iter().copied().collect();
            
            let var_desc = VarDesc::Ins { seq: ins_seq };
            
            // Store insertions in non_insertion_vars (they're distinguished by VarDesc type)
            // use read_pos for tp calculation
            let var = get_variants_from_map(&mut self.non_insertion_vars, pos, &var_desc);
            add_cnt(var, is_reverse, read_pos, avg_qual, mapq, nm, Some(read_len_including_match_ins));
            
            inc_cnt(&mut self.ref_coverage, pos, 1);
        } else if s.starts_with(b"-") {
            // Deletion: -2 or -2&AT format (deletion with optional following sequence)
            // This is handled elsewhere in deletion processing, but add fallback
            // For now, skip complex deletion handling
        } else {
            // Unknown format - skip
        }
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
    
    /// Get a reference base at a genomic position
    /// Converts genomic position to slice index using region.start()
    fn ref_base_at(&self, genomic_pos: i64) -> Option<u8> {
        if genomic_pos < self.region.start() as i64 {
            return None;
        }
        let idx = (genomic_pos - self.region.start() as i64) as usize;
        self.reference.ref_seq.get(idx).copied()
    }
    
    /// Check if a reference position has a base that equals the given base
    fn ref_has_and_equals(&self, genomic_pos: i64, base: u8) -> bool {
        self.ref_base_at(genomic_pos).map_or(false, |b| b == base)
    }
    
    /// Check if a reference position has a base that does NOT equal the given base
    fn ref_has_and_not_equals(&self, genomic_pos: i64, base: u8) -> bool {
        self.ref_base_at(genomic_pos).map_or(false, |b| b != base)
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
                // BAM stores Phred quality directly
                add_cnt(
                    seq_var,
                    is_reverse,
                    si as usize - (cigar_len as usize - num_high_qual_base),
                    *query_quality.get_or_err(si as usize)?,
                    mapq,
                    num_mismatch,
                    None,
                );
            }

            add_cnt(
                &mut sclip.var,
                is_reverse,
                cigar_len as usize,
                (read_qual_sum / num_high_qual_base) as u8,
                mapq,
                num_mismatch,
                None,
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
                // BAM stores Phred quality directly
                add_cnt(
                    seq_var,
                    is_reverse,
                    num_high_qual_base - si,
                    *query_quality.get_or_err(self.read_pos_including_softclip + si)?,
                    mapq,
                    num_mismatch,
                    None,
                );
            }

            add_cnt(
                &mut sclip.var,
                is_reverse,
                cigar_len as usize,
                (read_qual_sum / num_high_qual_base) as u8,
                mapq,
                num_mismatch,
                None,
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
/// read_pos: position in read (excluding soft clips)  
/// bq: base quality
/// read_len: optional total read length for calculating tp correctly
fn add_cnt(var: &mut Variant, is_reverse: bool, read_pos: usize, bq: u8, mapq: u8, nm: usize, read_len: Option<usize>) {
    var.alt_depth += 1;
    var.inc_dir(is_reverse);
    
    // Calculate tp: minimum distance from either end of read (like Java)
    // Java: tp = n < rlen1 - n ? n + 1 : rlen1 - n
    let tp = if let Some(rlen) = read_len {
        let from_start = read_pos;
        let from_end = rlen.saturating_sub(read_pos);
        if from_start < from_end {
            from_start + 1
        } else {
            from_end
        }
    } else {
        // Fallback if read_len not provided - use read_pos + 1
        // This is not ideal but maintains backward compatibility
        read_pos + 1
    };
    
    var.mean_pos += tp as f64;
    var.mean_qual += bq as f64;
    var.mean_mapq += mapq as f64;
    var.nm += nm as f64;

    if bq as f64 >= instance().conf.goodq {
        var.high_qual_read_cnt += 1;
    } else {
        var.low_qual_read_cnt += 1;
    }
}

/// Increment variant counters for insertion/deletion anchor positions.
/// These don't count toward high quality reads.
/// read_pos: position in read (excluding soft clips)  
/// bq: base quality
/// read_len: optional total read length for calculating tp correctly
fn add_cnt_anchor(var: &mut Variant, is_reverse: bool, read_pos: usize, bq: u8, mapq: u8, nm: usize, read_len: Option<usize>) {
    var.alt_depth += 1;
    var.inc_dir(is_reverse);
    
    // Calculate tp: minimum distance from either end of read (like Java)
    let tp = if let Some(rlen) = read_len {
        let from_start = read_pos;
        let from_end = rlen.saturating_sub(read_pos);
        if from_start < from_end {
            from_start + 1
        } else {
            from_end
        }
    } else {
        read_pos + 1
    };
    
    var.mean_pos += tp as f64;
    var.mean_qual += bq as f64;
    var.mean_mapq += mapq as f64;
    var.nm += nm as f64;
    // Note: Don't increment high_qual_read_cnt or low_qual_read_cnt for anchors
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
    // Handle empty cigar case
    if cigar.len() == 0 {
        return Ok(false);
    }
    
    // if the current cigar is not the last item.
    // and the next cigar is ref skip
    // or
    // the current cigar is not the first one and the previous one is refskip
    let has_prev_refskip = ci > 0 && cigar.get(ci - 1).map_or(false, |c| matches!(c, Cigar::RefSkip(_)));
    let has_next_refskip = ci + 1 < cigar.len() && cigar.get(ci + 1).map_or(false, |c| matches!(c, Cigar::RefSkip(_)));
    
    Ok(has_prev_refskip || has_next_refskip)
}

fn is_followed_by_match_and_indel(cigar: &CigarStringView, ci: usize) -> bool {
    if !instance().conf.perform_local_realignment || ci + 2 >= cigar.len() {
        return false;
    }

    let n_cigar = cigar.get(ci + 1).unwrap();
    let nn_cigar = cigar.get(ci + 2).unwrap();

    matches!(n_cigar, &Cigar::Match(l) if l <= instance().conf.vext as u32)
        && matches!(nn_cigar, Cigar::Ins(_) | Cigar::Del(_))
        && {
            match cigar.get(ci + 3) {
                Some(Cigar::Ins(_) | Cigar::Del(_)) => false,
                Some(_) => true,
                None => true,
            }
        }
}

/// Append sequence for deletion or insertion cases to create description string
/// and quality string.
fn append_segments(
    query_seq: &[u8],
    query_qual: &[u8],
    cigar: &CigarStringView,
    ci: usize,
    var_desc: &mut VarDesc,
    qual_seg: &mut SmallVecBytes,
    begin: usize,
    mlen: usize,
    indel_len: usize,
    is_ins: bool,
) -> Result<(), Error> {
    let VarDesc::Del {
        len: del_len,
        match_seq,
        ins_or_del_len,
        mismatch_seq: _,
    } = var_desc
    else {
        panic!("Not deletion description: {:?}", var_desc)
    };

    // begin is n + m for insertion and n for deletion
    // append to s '#' and part of read sequence corresponding to next CIGAR segment
    // (matched one)
    match_seq.extend_from_slice(query_seq.get_or_err(begin..(begin + mlen))?);
    // append quality string of next matched segment from read
    qual_seg.extend_from_slice(query_qual.get_or_err(begin..(begin + mlen))?);

    // if an insertion is two segments ahead, append '^' + part of sequence
    // corresponding
    // to next-next segment otherwise (deletion) append '^' + length of a next-next
    // segment
    if matches!(cigar.get(ci + 2).unwrap(), Cigar::Ins(_)) {
        *ins_or_del_len = InsOrDelLen::InsSeq(
            query_seq
                .get_or_err(begin + mlen..begin + mlen + indel_len)?
                .into(),
        )
    } else {
        *ins_or_del_len = InsOrDelLen::DelLen(indel_len);
    }

    // if an insertion is two segments ahead, append part of quality string sequence
    // corresponding to next-next segment otherwise (deletion)
    // append first quality score of next segment or return empty string
    // Quality handling
    if matches!(cigar.get(ci + 2).unwrap(), Cigar::Ins(_)) {
        // ci+2 is Insertion → always append slice
        qual_seg.extend_from_slice(query_qual.get_or_err(begin + mlen..begin + mlen + indel_len)?);
    } else {
        // ci+2 is Deletion
        if is_ins {
            // called from insertion context → single char
            qual_seg.push(query_qual.get_or_err(begin + mlen).copied()?);
        }
        // else: called from deletion context → nothing
    }

    Ok(())
}

#[inline]
fn is_next_after_num_matched(cigar: &CigarStringView, ci: usize, number: usize) -> bool {
    cigar
        .get(ci + number)
        .map_or(false, |c| matches!(c, Cigar::Match(_)))
}

#[inline]
fn is_next_ins(cigar: &CigarStringView, ci: usize) -> bool {
    if !instance().conf.perform_local_realignment {
        return false
    }

    cigar.get(ci+1).map_or(false, |c| {
        matches!(c, Cigar::Ins(_))
    })
}

#[inline]
fn is_next_matched(cigar: &CigarStringView, ci: usize) -> bool {
    if !instance().conf.perform_local_realignment {
        return false
    }

    cigar.get(ci + 1).map_or(false, |c| matches!(c, Cigar::Match(_)))
}

/// Check if a character is a valid DNA base (A, T, G, or C)
#[inline]
fn is_atgc(c: u8) -> bool {
    matches!(c, b'A' | b'T' | b'G' | b'C')
}

/// Check if sequence matches pattern: ^[ATGC]&[ATGC]+$
/// First character must be ATGC, second must be '&', rest must be ATGC
/// Used to detect complex variant descriptions like "A&ACGT"
pub fn is_begin_atgc_amp_atgcs_end(sequence: &[u8]) -> bool {
    if sequence.len() > 2 {
        let first_char = sequence[0];
        let second_char = sequence[1];
        if second_char == b'&' && is_atgc(first_char) {
            for &c in &sequence[2..] {
                if !is_atgc(c) {
                    return false;
                }
            }
            return true;
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::CigarString;

    /// Helper to create a CigarStringView from CIGAR elements
    fn make_cigar(elements: Vec<Cigar>) -> CigarStringView {
        CigarString(elements).into_view(0)
    }

    // Tests ported from CigarModifierTest.java

    /// Test getInsertionDeletionLength - calculates sum of I + D lengths
    /// Java: Cigar with 1M2S4I8D16N32H expects 12 (4+8)
    #[test]
    fn test_get_insertion_deletion_length() {
        let cigar = make_cigar(vec![
            Cigar::Match(1),
            Cigar::SoftClip(2),
            Cigar::Ins(4),
            Cigar::Del(8),
            Cigar::RefSkip(16),
            Cigar::HardClip(32),
        ]);

        assert_eq!(get_ins_del_len(&cigar), 12); // 4 + 8
    }

    /// Test getMatchInsertionLength - calculates sum of M + I lengths
    /// Java: Cigar with 1M2S4I8D16N32H expects 5 (1+4)
    #[test]
    fn test_get_match_insertion_length() {
        let cigar = make_cigar(vec![
            Cigar::Match(1),
            Cigar::SoftClip(2),
            Cigar::Ins(4),
            Cigar::Del(8),
            Cigar::RefSkip(16),
            Cigar::HardClip(32),
        ]);

        assert_eq!(get_match_insertion_length(&cigar), 5); // 1 + 4
    }

    /// Test getSoftClippedLength - calculates sum of M + I + S lengths
    /// Java: Cigar with 1M2S4I8D16N32H expects 7 (1+2+4)
    #[test]
    fn test_get_soft_clipped_length() {
        let cigar = make_cigar(vec![
            Cigar::Match(1),
            Cigar::SoftClip(2),
            Cigar::Ins(4),
            Cigar::Del(8),
            Cigar::RefSkip(16),
            Cigar::HardClip(32),
        ]);

        assert_eq!(get_soft_clipped_length(&cigar), 7); // 1 + 2 + 4
    }

    /// Test getCigarOperator - retrieves operator at each position
    #[test]
    fn test_get_cigar_operator() {
        let cigar = make_cigar(vec![
            Cigar::Match(1),
            Cigar::SoftClip(2),
            Cigar::Ins(4),
            Cigar::Del(8),
            Cigar::RefSkip(16),
            Cigar::HardClip(32),
        ]);

        // Check each element matches expected type
        assert!(matches!(cigar.get(0), Some(Cigar::Match(1))));
        assert!(matches!(cigar.get(1), Some(Cigar::SoftClip(2))));
        assert!(matches!(cigar.get(2), Some(Cigar::Ins(4))));
        assert!(matches!(cigar.get(3), Some(Cigar::Del(8))));
        assert!(matches!(cigar.get(4), Some(Cigar::RefSkip(16))));
        assert!(matches!(cigar.get(5), Some(Cigar::HardClip(32))));
    }

    // Tests for is_begin_atgc_amp_atgcs_end - ported from CigarModifierTest.java

    /// Test "A&ACGT" -> true (first char ATGC, second is &, rest is ATGC)
    #[test]
    fn test_is_begin_atgc_amp_atgcs_end_valid() {
        assert!(is_begin_atgc_amp_atgcs_end(b"A&ACGT"));
    }

    /// Test "A&" -> false (no chars after &)
    #[test]
    fn test_is_begin_atgc_amp_atgcs_end_no_suffix() {
        assert!(!is_begin_atgc_amp_atgcs_end(b"A&"));
    }

    /// Test "A&ASGT" -> false (S is not ATGC)
    #[test]
    fn test_is_begin_atgc_amp_atgcs_end_invalid_char() {
        assert!(!is_begin_atgc_amp_atgcs_end(b"A&ASGT"));
    }
}