use bio_types::genome::AbstractInterval;
use smallvec::SmallVec;
use std::{
    collections::{BTreeMap, HashMap, HashSet, VecDeque},
    env,
    ops::AddAssign,
    panic::{AssertUnwindSafe, catch_unwind},
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
        variants::{InsOrDelLen, Mate, SoftClip, VarDesc, Variant},
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
    splice_count_insert_index: HashMap<SplicingKey, usize>,
    next_splice_count_insert_index: usize,

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
    non_insertion_vars_insert_index: HashMap<i64, usize>,
    next_non_insertion_vars_insert_index: usize,
    insertion_vars: HashMap<i64, HashMap<VarDesc, Variant>>,

    /// Track MNPs (multi-nucleotide polymorphisms) by position and description string
    mnp: HashMap<i64, HashMap<String, usize>>,

    /// Track insertion counts by position and description string (Java: positionToInsertionCount)
    position_to_insertion_count: HashMap<i64, HashMap<String, usize>>,

    /// Track deletion counts by position and description string (Java: positionToDeletionCount)
    position_to_deletions_count: HashMap<i64, HashMap<String, usize>>,

    ref_coverage: HashMap<i64, usize>,

    soft_clips5_end: HashMap<i64, SoftClip>,
    soft_clips3_end: HashMap<i64, SoftClip>,
    svdelfend: i64,
    svdelrend: i64,
    svdupfend: i64,
    svduprend: i64,
    svinvfend5: i64,
    svinvrend5: i64,
    svinvfend3: i64,
    svinvrend3: i64,
    svfdel: Vec<SoftClip>,
    svrdel: Vec<SoftClip>,
    svfdup: Vec<SoftClip>,
    svrdup: Vec<SoftClip>,
    svfinv5: Vec<SoftClip>,
    svrinv5: Vec<SoftClip>,
    svfinv3: Vec<SoftClip>,
    svrinv3: Vec<SoftClip>,
    current_qname: Option<String>,
    last_modified_pos: Option<i64>,
    last_modified_cigar: Option<String>,
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
            splice_count_insert_index: Default::default(),
            next_splice_count_insert_index: 0,
            read_pos_including_softclip: Default::default(),
            read_pos_excluding_softclip: Default::default(),
            start: Default::default(),
            offset: Default::default(),
            cigar_len: Default::default(),
            non_insertion_vars: Default::default(),
            non_insertion_vars_insert_index: Default::default(),
            next_non_insertion_vars_insert_index: 0,
            insertion_vars: Default::default(),
            mnp: Default::default(),
            position_to_insertion_count: Default::default(),
            position_to_deletions_count: Default::default(),
            ref_coverage: Default::default(),
            soft_clips5_end: Default::default(),
            soft_clips3_end: Default::default(),
            svdelfend: 0,
            svdelrend: 0,
            svdupfend: 0,
            svduprend: 0,
            svinvfend5: 0,
            svinvrend5: 0,
            svinvfend3: 0,
            svinvrend3: 0,
            svfdel: Default::default(),
            svrdel: Default::default(),
            svfdup: Default::default(),
            svrdup: Default::default(),
            svfinv5: Default::default(),
            svrinv5: Default::default(),
            svfinv3: Default::default(),
            svrinv3: Default::default(),
            rev_complementor: RevComplementor::new(),
            cigar: CigarString(vec![]).into_view(0),
            current_qname: None,
            last_modified_pos: None,
            last_modified_cigar: None,
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
            splice_count_insert_index: HashMap::new(),
            next_splice_count_insert_index: 0,
            read_pos_including_softclip: 0,
            read_pos_excluding_softclip: 0,
            start: 0,
            offset: 0,
            cigar_len: 0,
            non_insertion_vars: HashMap::new(),
            non_insertion_vars_insert_index: HashMap::new(),
            next_non_insertion_vars_insert_index: 0,
            insertion_vars: HashMap::new(),
            mnp: HashMap::new(),
            position_to_insertion_count: HashMap::new(),
            position_to_deletions_count: HashMap::new(),
            ref_coverage: HashMap::new(),
            soft_clips5_end: HashMap::new(),
            soft_clips3_end: HashMap::new(),
            svdelfend: 0,
            svdelrend: 0,
            svdupfend: 0,
            svduprend: 0,
            svinvfend5: 0,
            svinvrend5: 0,
            svinvfend3: 0,
            svinvrend3: 0,
            svfdel: Vec::new(),
            svrdel: Vec::new(),
            svfdup: Vec::new(),
            svrdup: Vec::new(),
            svfinv5: Vec::new(),
            svrinv5: Vec::new(),
            svfinv3: Vec::new(),
            svrinv3: Vec::new(),
            rev_complementor: RevComplementor::new(),
            cigar: CigarString(vec![]).into_view(0),
            current_qname: None,
            last_modified_pos: None,
            last_modified_cigar: None,
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
            self.process_record(record)?;
        }
        Ok(())
    }

    /// Process a single BAM record (streaming)
    pub fn process_record(&mut self, record: &mut Record) -> Result<(), Error> {
        let record_name = String::from_utf8_lossy(record.qname()).to_string();
        match catch_unwind(AssertUnwindSafe(|| self.parse_cigar(record))) {
            Ok(Ok(())) => Ok(()),
            Ok(Err(error)) => crate::utils::print_exception_and_continue(
                &error,
                "record",
                &record_name,
                Some(&self.region),
                &self.instance.conf,
            ),
            Err(payload) => {
                let panic_message = if let Some(message) = payload.downcast_ref::<&str>() {
                    (*message).to_string()
                } else if let Some(message) = payload.downcast_ref::<String>() {
                    message.clone()
                } else {
                    "unknown panic payload".to_string()
                };
                let error = anyhow!(
                    "panic while processing record '{}': {}",
                    record_name,
                    panic_message
                );
                crate::utils::print_exception_and_continue(
                    &error,
                    "record",
                    &record_name,
                    Some(&self.region),
                    &self.instance.conf,
                )
            }
        }
    }

    /// Get the collected non-insertion variants
    pub fn get_non_insertion_vars(&self) -> &HashMap<i64, HashMap<VarDesc, Variant>> {
        &self.non_insertion_vars
    }

    /// Take ownership of the collected non-insertion variants
    pub fn take_non_insertion_vars(&mut self) -> HashMap<i64, HashMap<VarDesc, Variant>> {
        std::mem::take(&mut self.non_insertion_vars)
    }

    pub fn take_non_insertion_vars_insert_index(&mut self) -> HashMap<i64, usize> {
        std::mem::take(&mut self.non_insertion_vars_insert_index)
    }

    fn get_non_insertion_variant(&mut self, pos: i64, var_desc: &VarDesc) -> &mut Variant {
        if !self.non_insertion_vars.contains_key(&pos) {
            self.non_insertion_vars_insert_index
                .insert(pos, self.next_non_insertion_vars_insert_index);
            self.next_non_insertion_vars_insert_index += 1;
        }

        get_variants_from_map(&mut self.non_insertion_vars, pos, var_desc)
    }

    /// Get the collected insertion variants
    pub fn get_insertion_vars(&self) -> &HashMap<i64, HashMap<VarDesc, Variant>> {
        &self.insertion_vars
    }

    /// Take ownership of the collected insertion variants
    pub fn take_insertion_vars(&mut self) -> HashMap<i64, HashMap<VarDesc, Variant>> {
        std::mem::take(&mut self.insertion_vars)
    }

    pub fn take_mnp(&mut self) -> HashMap<i64, HashMap<String, usize>> {
        std::mem::take(&mut self.mnp)
    }

    pub fn take_position_to_insertion_count(&mut self) -> HashMap<i64, HashMap<String, usize>> {
        std::mem::take(&mut self.position_to_insertion_count)
    }

    pub fn take_position_to_deletions_count(&mut self) -> HashMap<i64, HashMap<String, usize>> {
        std::mem::take(&mut self.position_to_deletions_count)
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

    /// Take forward deletion discordant clusters
    pub fn take_svfdel(&mut self) -> Vec<SoftClip> {
        std::mem::take(&mut self.svfdel)
    }

    /// Take reverse deletion discordant clusters
    pub fn take_svrdel(&mut self) -> Vec<SoftClip> {
        std::mem::take(&mut self.svrdel)
    }

    pub fn take_svfdup(&mut self) -> Vec<SoftClip> {
        std::mem::take(&mut self.svfdup)
    }

    pub fn take_svrdup(&mut self) -> Vec<SoftClip> {
        std::mem::take(&mut self.svrdup)
    }

    pub fn take_svfinv5(&mut self) -> Vec<SoftClip> {
        std::mem::take(&mut self.svfinv5)
    }

    pub fn take_svrinv5(&mut self) -> Vec<SoftClip> {
        std::mem::take(&mut self.svrinv5)
    }

    pub fn take_svfinv3(&mut self) -> Vec<SoftClip> {
        std::mem::take(&mut self.svfinv3)
    }

    pub fn take_svrinv3(&mut self) -> Vec<SoftClip> {
        std::mem::take(&mut self.svrinv3)
    }

    /// Take ownership of splice counts
    pub fn take_splice_count(&mut self) -> HashMap<SplicingKey, Vec<usize>> {
        std::mem::take(&mut self.splice_count)
    }

    pub fn take_splice_count_insert_index(&mut self) -> HashMap<SplicingKey, usize> {
        std::mem::take(&mut self.splice_count_insert_index)
    }

    fn parse_cigar_with_amp_case(
        &self,
        record: &Record,
        cigar: &CigarStringView,
        is_mate_reference_name_equal: bool,
    ) -> bool {
        let (distance_to_amplicon, overlap_fraction) = self
            .instance
            .amplicon_based_calling
            .as_deref()
            .and_then(|value| {
                let split: Vec<&str> = value.split(':').collect();
                let distance = split.first()?.parse::<i64>().ok()?;
                let overlap = split.get(1)?.parse::<f64>().ok()?;
                Some((distance, overlap))
            })
            .unwrap_or((10, 0.95));

        let read_len_match_del = get_aligned_length(cigar);
        let mut seg_start = record.pos() + 1;
        let mut seg_end = seg_start + read_len_match_del - 1;

        if matches!(cigar.iter().next(), Some(Cigar::SoftClip(_))) {
            let ts1 = seg_start.max(self.region.start as i64);
            let te1 = seg_end.min(self.region.end as i64);
            let overlap = ((ts1 - te1).abs() as f64) / ((seg_end - seg_start) as f64);
            if !(overlap > overlap_fraction) {
                return true;
            }
        } else if matches!(cigar.iter().last(), Some(Cigar::SoftClip(_))) {
            let ts1 = seg_start.max(self.region.start as i64);
            let te1 = seg_end.min(self.region.end as i64);
            let overlap = ((te1 - ts1).abs() as f64) / ((seg_end - seg_start) as f64);
            if !(overlap > overlap_fraction) {
                return true;
            }
        } else {
            if is_mate_reference_name_equal && record.insert_size() != 0 {
                if record.insert_size() > 0 {
                    seg_end = seg_start + record.insert_size() - 1;
                } else {
                    seg_start = record.mpos() + 1;
                    seg_end = record.mpos() + 1 - record.insert_size() - 1;
                }
            }

            let ts1 = seg_start.max(self.region.start as i64);
            let te1 = seg_end.min(self.region.end as i64);
            let overlap = (((ts1 - te1) as f64) / ((seg_end - seg_start) as f64)).abs();

            if ((seg_start - self.region.start as i64).abs() > distance_to_amplicon
                || (seg_end - self.region.end as i64).abs() > distance_to_amplicon)
                || overlap <= overlap_fraction
            {
                return true;
            }
        }

        false
    }

    fn add_discordant_cluster(
        sdref: &mut SoftClip,
        start: i64,
        end: i64,
        mate_start: i64,
        mate_end: i64,
        dir: i64,
        read_len: i64,
        mlen: i32,
        softp: i64,
        pmean: f64,
        qmean: f64,
        mapq: f64,
        nm: f64,
        goodq: f64,
    ) {
        sdref.var.alt_depth += 1;
        sdref.var.extra_cnt += 1;
        sdref.var.mean_pos += pmean;
        sdref.var.mean_qual += qmean;
        sdref.var.mean_mapq += mapq;
        sdref.var.nm += nm;
        if dir == 1 {
            sdref.var.alt_depth_fwd += 1;
        } else {
            sdref.var.alt_depth_rev += 1;
        }
        if qmean >= goodq {
            sdref.var.high_qual_read_cnt += 1;
        } else {
            sdref.var.low_qual_read_cnt += 1;
        }
        if sdref.start == 0 || sdref.start >= start {
            sdref.start = start;
        }
        if sdref.end == 0 || sdref.end <= end {
            sdref.end = end;
        }
        sdref.mates.push(Mate {
            mate_start,
            mate_end,
            mate_len: mlen,
            start,
            end,
            mean_pos: pmean,
            mean_qual: qmean,
            mean_mapq: mapq,
            nm,
        });
        if sdref.mstart == 0 || sdref.mstart >= mate_start {
            sdref.mstart = mate_start;
        }
        if sdref.mend == 0 || sdref.mend <= mate_end {
            sdref.mend = mate_start + read_len;
        }

        if softp != 0 {
            if dir == 1 {
                if (softp - sdref.end).abs() < 10 {
                    *sdref.soft.entry(softp).or_insert(0) += 1;
                }
            } else if (softp - sdref.start).abs() < 10 {
                *sdref.soft.entry(softp).or_insert(0) += 1;
            }
        }
    }

    fn add_discordant_count(svref: &mut SoftClip) {
        svref.disc += 1;
    }

    fn aux_as_i64(aux: &record::Aux<'_>) -> Option<i64> {
        match aux {
            record::Aux::I8(v) => Some(*v as i64),
            record::Aux::U8(v) => Some(*v as i64),
            record::Aux::I16(v) => Some(*v as i64),
            record::Aux::U16(v) => Some(*v as i64),
            record::Aux::I32(v) => Some(*v as i64),
            record::Aux::U32(v) => Some(*v as i64),
            _ => None,
        }
    }

    fn mc_has_softclip_both_ends(mc: &str) -> bool {
        let mut seen_softclip = false;
        let bytes = mc.as_bytes();
        for idx in 1..bytes.len() {
            if bytes[idx] == b'S' && bytes[idx - 1].is_ascii_digit() {
                if seen_softclip {
                    return true;
                }
                seen_softclip = true;
            }
        }
        false
    }

    fn prepare_sv_deletion_structures_for_analysis(
        &mut self,
        record: &Record,
        query_qual: &[u8],
        number_of_mismatches: i32,
        read_direction: bool,
        mate_direction: bool,
        start: i64,
        total_length_including_softclip: usize,
        cigar: &CigarStringView,
    ) {
        if record.tid() != record.mtid() {
            return;
        }

        let min_map_base = Configuration::MINMAPBASE;
        if query_qual.len() <= min_map_base {
            return;
        }

        let mate_start = record.mpos() + 1;
        let mend = mate_start + total_length_including_softclip as i64;
        let aligned_len = get_aligned_length_mnd(cigar);
        let end = start + aligned_len;

        let mut soft5 = 0i64;
        if let Some(Cigar::SoftClip(len)) = cigar.iter().next() {
            let tt = *len as usize;
            if tt > 0 && tt <= query_qual.len() && query_qual[tt - 1] as f64 > instance().conf.goodq {
                soft5 = start;
            }
        }

        let mut soft3 = 0i64;
        if let Some(Cigar::SoftClip(len)) = cigar.iter().last() {
            let tt = *len as usize;
            if tt > 0 && tt <= query_qual.len() {
                let qi = query_qual.len().saturating_sub(tt);
                if qi < query_qual.len() && query_qual[qi] as f64 > instance().conf.goodq {
                    soft3 = end;
                }
            }
        }

        let read_dir_num = if read_direction { -1i64 } else { 1i64 };
        let mate_dir_num = if mate_direction { 1i64 } else { -1i64 };
        let mlen = record.insert_size();

        if let Ok(Some(mc_aux)) = record.aux_option(b"MC") {
            if let Ok(mc_tag) = mc_aux.try_get_str() {
                if Self::mc_has_softclip_both_ends(mc_tag) {
                    return;
                }
            }
        }

        if let Ok(Some(mq_aux)) = record.aux_option(b"MQ") {
            if let Some(mq) = Self::aux_as_i64(&mq_aux) {
                if mq < 15 {
                    return;
                }
            }
        }

        let min_cluster_dist = Configuration::MINSVCDIST * self.max_read_len.max(1) as f64;
        let min_d = 75i64;

        let qmean = query_qual[min_map_base] as f64;
        let pmean = self.max_read_len.max(1) as f64 / 2.0;
        if read_dir_num * mate_dir_num == -1 && mlen * read_dir_num > 0 {
            let span = if mate_start > start {
                mend - start
            } else {
                end - mate_start
            };
            if span.abs()
                <= (instance().conf.inssize + instance().conf.insstdamt * instance().conf.insstd) as i64
            {
                return;
            }

            if read_dir_num == 1 {
                if self.svfdel.is_empty() || (start - self.svdelfend) as f64 > min_cluster_dist {
                    self.svfdel.push(SoftClip::default());
                }
                if let Some(last) = self.svfdel.last_mut() {
                    Self::add_discordant_cluster(
                        last,
                        start,
                        end,
                        mate_start,
                        mend,
                        read_dir_num,
                        total_length_including_softclip as i64,
                        span as i32,
                        soft3,
                        pmean,
                        qmean,
                        record.mapq() as f64,
                        number_of_mismatches as f64,
                        instance().conf.goodq,
                    );
                }
                self.svdelfend = end;
            } else {
                if self.svrdel.is_empty() || (start - self.svdelrend) as f64 > min_cluster_dist {
                    self.svrdel.push(SoftClip::default());
                }
                if let Some(last) = self.svrdel.last_mut() {
                    Self::add_discordant_cluster(
                        last,
                        start,
                        end,
                        mate_start,
                        mend,
                        read_dir_num,
                        total_length_including_softclip as i64,
                        span as i32,
                        soft5,
                        pmean,
                        qmean,
                        record.mapq() as f64,
                        number_of_mismatches as f64,
                        instance().conf.goodq,
                    );
                }
                self.svdelrend = end;
            }

            if !self.svfdel.is_empty() && ((start - self.svdelfend).abs() as f64) <= min_cluster_dist {
                if let Some(last) = self.svfdel.last_mut() {
                    Self::add_discordant_count(last);
                }
            }
            if !self.svrdel.is_empty() && ((start - self.svdelrend).abs() as f64) <= min_cluster_dist {
                if let Some(last) = self.svrdel.last_mut() {
                    Self::add_discordant_count(last);
                }
            }
            if !self.svfdup.is_empty() && (start - self.svdupfend).abs() <= min_d {
                if let Some(last) = self.svfdup.last_mut() {
                    Self::add_discordant_count(last);
                }
            }
            if !self.svrdup.is_empty() && (start - self.svduprend).abs() <= min_d {
                if let Some(last) = self.svrdup.last_mut() {
                    Self::add_discordant_count(last);
                }
            }
            if !self.svfinv5.is_empty() && (start - self.svinvfend5).abs() <= min_d {
                if let Some(last) = self.svfinv5.last_mut() {
                    Self::add_discordant_count(last);
                }
            }
            if !self.svrinv5.is_empty() && (start - self.svinvrend5).abs() <= min_d {
                if let Some(last) = self.svrinv5.last_mut() {
                    Self::add_discordant_count(last);
                }
            }
            if !self.svfinv3.is_empty() && (start - self.svinvfend3).abs() <= min_d {
                if let Some(last) = self.svfinv3.last_mut() {
                    Self::add_discordant_count(last);
                }
            }
            if !self.svrinv3.is_empty() && (start - self.svinvrend3).abs() <= min_d {
                if let Some(last) = self.svrinv3.last_mut() {
                    Self::add_discordant_count(last);
                }
            }
            return;
        }

        if read_dir_num * mate_dir_num == -1 && mlen * read_dir_num < 0 {
            if read_dir_num == 1 {
                if self.svfdup.is_empty() || (start - self.svdupfend) as f64 > min_cluster_dist {
                    self.svfdup.push(SoftClip::default());
                }
                if let Some(last) = self.svfdup.last_mut() {
                    Self::add_discordant_cluster(
                        last,
                        start,
                        end,
                        mate_start,
                        mend,
                        read_dir_num,
                        total_length_including_softclip as i64,
                        mlen as i32,
                        soft3,
                        pmean,
                        qmean,
                        record.mapq() as f64,
                        number_of_mismatches as f64,
                        instance().conf.goodq,
                    );
                }
                self.svdupfend = end;
            } else {
                if self.svrdup.is_empty() || (start - self.svduprend) as f64 > min_cluster_dist {
                    self.svrdup.push(SoftClip::default());
                }
                if let Some(last) = self.svrdup.last_mut() {
                    Self::add_discordant_cluster(
                        last,
                        start,
                        end,
                        mate_start,
                        mend,
                        read_dir_num,
                        total_length_including_softclip as i64,
                        mlen as i32,
                        soft5,
                        pmean,
                        qmean,
                        record.mapq() as f64,
                        number_of_mismatches as f64,
                        instance().conf.goodq,
                    );
                }
                self.svduprend = end;
            }

            if !self.svfdup.is_empty() && ((start - self.svdupfend).abs() as f64) <= min_cluster_dist {
                if let Some(last) = self.svfdup.last_mut() {
                    Self::add_discordant_count(last);
                }
            }
            if !self.svrdup.is_empty() && ((start - self.svduprend).abs() as f64) <= min_cluster_dist {
                if let Some(last) = self.svrdup.last_mut() {
                    Self::add_discordant_count(last);
                }
            }
            if !self.svfdel.is_empty() && (start - self.svdelfend).abs() <= min_d {
                if let Some(last) = self.svfdel.last_mut() {
                    Self::add_discordant_count(last);
                }
            }
            if !self.svrdel.is_empty() && (start - self.svdelrend).abs() <= min_d {
                if let Some(last) = self.svrdel.last_mut() {
                    Self::add_discordant_count(last);
                }
            }
            if !self.svfinv5.is_empty() && (start - self.svinvfend5).abs() <= min_d {
                if let Some(last) = self.svfinv5.last_mut() {
                    Self::add_discordant_count(last);
                }
            }
            if !self.svrinv5.is_empty() && (start - self.svinvrend5).abs() <= min_d {
                if let Some(last) = self.svrinv5.last_mut() {
                    Self::add_discordant_count(last);
                }
            }
            if !self.svfinv3.is_empty() && (start - self.svinvfend3).abs() <= min_d {
                if let Some(last) = self.svfinv3.last_mut() {
                    Self::add_discordant_count(last);
                }
            }
            if !self.svrinv3.is_empty() && (start - self.svinvrend3).abs() <= min_d {
                if let Some(last) = self.svrinv3.last_mut() {
                    Self::add_discordant_count(last);
                }
            }
        } else if read_dir_num * mate_dir_num == 1 {
            let max_read_len_i64 = self.max_read_len.max(1) as i64;
            if read_dir_num == 1 && mlen != 0 {
                if mlen < -3 * max_read_len_i64 {
                    if self.svfinv3.is_empty() || (start - self.svinvfend3) as f64 > min_cluster_dist {
                        self.svfinv3.push(SoftClip::default());
                    }
                    if let Some(last) = self.svfinv3.last_mut() {
                        Self::add_discordant_cluster(
                            last,
                            start,
                            end,
                            mate_start,
                            mend,
                            read_dir_num,
                            total_length_including_softclip as i64,
                            mlen as i32,
                            soft3,
                            pmean,
                            qmean,
                            record.mapq() as f64,
                            number_of_mismatches as f64,
                            instance().conf.goodq,
                        );
                        Self::add_discordant_count(last);
                    }
                    self.svinvfend3 = end;
                } else if mlen > 3 * max_read_len_i64 {
                    if self.svfinv5.is_empty() || (start - self.svinvfend5) as f64 > min_cluster_dist {
                        self.svfinv5.push(SoftClip::default());
                    }
                    if let Some(last) = self.svfinv5.last_mut() {
                        Self::add_discordant_cluster(
                            last,
                            start,
                            end,
                            mate_start,
                            mend,
                            read_dir_num,
                            total_length_including_softclip as i64,
                            mlen as i32,
                            soft3,
                            pmean,
                            qmean,
                            record.mapq() as f64,
                            number_of_mismatches as f64,
                            instance().conf.goodq,
                        );
                        Self::add_discordant_count(last);
                    }
                    self.svinvfend5 = end;
                }
            } else if mlen != 0 {
                if mlen < -3 * max_read_len_i64 {
                    if self.svrinv3.is_empty() || (start - self.svinvrend3) as f64 > min_cluster_dist {
                        self.svrinv3.push(SoftClip::default());
                    }
                    if let Some(last) = self.svrinv3.last_mut() {
                        Self::add_discordant_cluster(
                            last,
                            start,
                            end,
                            mate_start,
                            mend,
                            read_dir_num,
                            total_length_including_softclip as i64,
                            mlen as i32,
                            soft5,
                            pmean,
                            qmean,
                            record.mapq() as f64,
                            number_of_mismatches as f64,
                            instance().conf.goodq,
                        );
                        Self::add_discordant_count(last);
                    }
                    self.svinvrend3 = end;
                } else if mlen > 3 * max_read_len_i64 {
                    if self.svrinv5.is_empty() || (start - self.svinvrend5) as f64 > min_cluster_dist {
                        self.svrinv5.push(SoftClip::default());
                    }
                    if let Some(last) = self.svrinv5.last_mut() {
                        Self::add_discordant_cluster(
                            last,
                            start,
                            end,
                            mate_start,
                            mend,
                            read_dir_num,
                            total_length_including_softclip as i64,
                            mlen as i32,
                            soft5,
                            pmean,
                            qmean,
                            record.mapq() as f64,
                            number_of_mismatches as f64,
                            instance().conf.goodq,
                        );
                        Self::add_discordant_count(last);
                    }
                    self.svinvrend5 = end;
                }
            }

            if mlen != 0 {
                if !self.svfdel.is_empty() && (start - self.svdelfend) <= min_d {
                    if let Some(last) = self.svfdel.last_mut() {
                        Self::add_discordant_count(last);
                    }
                }
                if !self.svrdel.is_empty() && (start - self.svdelrend) <= min_d {
                    if let Some(last) = self.svrdel.last_mut() {
                        Self::add_discordant_count(last);
                    }
                }
                if !self.svfdup.is_empty() && (start - self.svdupfend) <= min_d {
                    if let Some(last) = self.svfdup.last_mut() {
                        Self::add_discordant_count(last);
                    }
                }
                if !self.svrdup.is_empty() && (start - self.svduprend) <= min_d {
                    if let Some(last) = self.svrdup.last_mut() {
                        Self::add_discordant_count(last);
                    }
                }
            }
        }
    }


    fn parse_cigar(&mut self, record: &mut Record) -> Result<(), Error> {
        event!(Level::DEBUG, "Starting for record at pos {}", record.pos());
        self.current_qname = Some(String::from_utf8_lossy(record.qname()).to_string());
        let trace_target = env::var("VARDICT_TRACE_QNAME")
            .ok()
            .map(|needle| {
                needle == "*" || record.qname() == needle.as_bytes()
            })
            .unwrap_or(record.qname() == b"SRR098401.96368837");
        if trace_target {
            event!(
                Level::WARN,
                "[trace_target] qname={} parse_enter pos={} cigar={}",
                String::from_utf8_lossy(record.qname()),
                record.pos() + 1,
                record.cigar().to_string()
            );
        }
        
        // Build query sequence and quality as owned vectors
        let mut query_seq_owned: Vec<u8> = record.seq().into_decoded_base_iter().collect();
        let mut query_qual_owned: Vec<u8> = record.qual().to_vec();

        let is_reverse = record.is_reverse();

        let mut query_seq = query_seq_owned.as_slice();
        let mut query_qual = query_qual_owned.as_slice();

        let mapping_quality = record.mapq();

        record.cache_cigar_if_empty();
        let mut cigar = record.cigar();
        
        event!(Level::DEBUG, "CIGAR: {:?}", cigar);
        if record.qname() == b"read_1" {
            event!(
                Level::DEBUG,
                "[debug read_1] pos={} cigar_before={:?}",
                record.pos(),
                cigar
            );
        }

        let ins_del_len = get_ins_del_len(&cigar);

        let tot_nm = match record.aux_option(self.aligner.nm_tag())? {
            Some(rust_htslib::bam::record::Aux::I32(nm)) => nm - ins_del_len as i32,
            Some(rust_htslib::bam::record::Aux::I8(nm)) => nm as i32 - ins_del_len as i32,
            Some(rust_htslib::bam::record::Aux::I16(nm)) => nm as i32 - ins_del_len as i32,
            Some(rust_htslib::bam::record::Aux::U8(nm)) => nm as i32 - ins_del_len as i32,
            Some(rust_htslib::bam::record::Aux::U16(nm)) => nm as i32 - ins_del_len as i32,
            Some(rust_htslib::bam::record::Aux::U32(nm)) => nm as i32 - ins_del_len as i32,
            Some(oth) => Err(anyhow!("Got unexpected type for NM tag: {:?}", oth))?,
            None => {
                if !cigar.is_empty() {
                    event!(Level::DEBUG, "No NM tag for mismatches; continuing with NM=0");
                }

                let has_alignment = !cigar.is_empty() && record.pos() >= 0;
                if record.is_unmapped() && !has_alignment {
                    return Ok(());
                }
                if cigar.is_empty() {
                    return Ok(());
                }

                0
            }
        };

        let nm = tot_nm;

        if trace_target {
            event!(
                Level::INFO,
                "[trace_target] qname={} initial_pos={} cigar={} nm={} ins_del_len={}",
                String::from_utf8_lossy(record.qname()),
                record.pos() + 1,
                cigar,
                nm,
                ins_del_len
            );
        }

        if nm > instance().conf.mismatch {
            if trace_target {
                event!(Level::INFO, "[trace_target] return: nm > mismatch");
            }
            return Ok(());
        }

        if self.instance.amplicon_based_calling.is_some()
            && self.parse_cigar_with_amp_case(record, &cigar, record.tid() == record.mtid())
        {
            if trace_target {
                event!(Level::INFO, "[trace_target] return: amplicon gate filtered read");
            }
            return Ok(());
        }

        let mut pos = 0;
        self.read_pos_including_softclip = 0;
        self.read_pos_excluding_softclip = 0;

        if self.instance.conf.perform_local_realignment {
            let region_offset = self.reference.region_start - 1;
            let local_pos = record.pos() - region_offset;
            if local_pos < 0 {
                // Reference slice doesn't cover the read start; skip local realignment for this read.
                pos = record.pos() + 1;
                self.last_modified_pos = Some(pos);
                self.last_modified_cigar = Some(cigar.to_string());
            } else {
                let mut local_cigar = CigarString(record.cigar().iter().copied().collect::<Vec<_>>())
                    .into_view(local_pos);
            // Modify the CIGAR for potential mis-alignment for indels at the end of reads to softclipping and let VarDict's
            // algorithm to figure out indels

                let mut cigar_modifier = CigarModifier::new(
                    local_pos,
                    &local_cigar,
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
                pos = mc.align_start_pos + self.reference.region_start;
                cigar = CigarString(mc.cigar.into_iter().collect::<Vec<_>>())
                    .into_view(pos);

                self.last_modified_pos = Some(pos);
                self.last_modified_cigar = Some(cigar.to_string());

            event!(Level::DEBUG, "Modified CIGAR: {:?}", cigar);
            if record.qname() == b"read_1" {
                event!(
                    Level::DEBUG,
                    "[debug read_1] cigar_after={:?} align_start={}",
                    cigar,
                    mc.align_start_pos
                );
            }

                query_qual = mc.query_qual;
                query_seq = mc.query_seq;
            }

        } else {
            // Convert to 1-based position (BAM is 0-based, Java VarDict uses 1-based)
            pos = record.pos() + 1;
            self.last_modified_pos = Some(pos);
            self.last_modified_cigar = Some(cigar.to_string());
        }

        self.clean_up_cigar(record);

        self.start = pos;
        self.offset = 0;

        //determine discordant reads
        if record.tid() != record.mtid() {
            self.discordant_count += 1;
        }

        //Ignore reads that are softclipped at both ends and both greater than 10 bp
        match cigar.0.as_slice() {
            [Cigar::SoftClip(sl1), .., Cigar::SoftClip(sl2)] if *sl1 >= 10 && *sl2 >= 10 => {
                if trace_target {
                    event!(
                        Level::INFO,
                        "[trace_target] return: both-end softclip check sl1={} sl2={} cigar={}",
                        sl1,
                        sl2,
                        cigar
                    );
                }
                return Ok(());
            }
            _ => {}
        }

        // Only match and insertion counts toward read length
        // For total length, including soft-clipped bases
        let read_match_ins_len = get_match_insertion_length(&cigar);

        if instance().conf.min_match != 0 && read_match_ins_len < instance().conf.min_match as usize
        {
            if trace_target {
                event!(
                    Level::INFO,
                    "[trace_target] return: min_match read_match_ins_len={} min_match={}",
                    read_match_ins_len,
                    instance().conf.min_match
                );
            }
            return Ok(());
        }

        // The total length, including soft-clipped bases
        let read_len_including_softclips = get_soft_clipped_length(&cigar);
        self.max_read_len = read_len_including_softclips.max(self.max_read_len);

        // If supplementary alignment is present
        if instance().conf.sam_filter != 0 {
            const SUPPLEMENTARY_ALIGNMENT: u16 = 0x800;
            if (record.flags() & SUPPLEMENTARY_ALIGNMENT) != 0 {
                if trace_target {
                    event!(Level::INFO, "[trace_target] return: supplementary alignment");
                }
                return Ok(());
            }
        }

        // Skip sites that are not in region of interest in CRISPR mode
        if self.skip_sites_out_region_of_interest(cigar.0.as_slice()) {
            if trace_target {
                event!(Level::INFO, "[trace_target] return: skip_sites_out_region_of_interest");
            }
            return Ok(());
        }

        if !instance().conf.disable_sv {
            if !(record.is_paired() && record.is_mate_unmapped()) && mapping_quality > 10 {
                self.prepare_sv_deletion_structures_for_analysis(
                    record,
                    query_qual,
                    nm,
                    is_reverse,
                    !record.is_mate_reverse(),
                    self.start,
                    read_len_including_softclips,
                    &cigar,
                );
            }
        }

        if trace_target {
            event!(
                Level::INFO,
                "[trace_target] entering process loop start={} read_match_ins_len={} soft_len={} modified_cigar={}",
                self.start,
                read_match_ins_len,
                read_len_including_softclips,
                self.last_modified_cigar.as_deref().unwrap_or("-")
            );
        }

        let alignment_start = self.start;
        let mpos = if record.mpos() >= 0 { record.mpos() + 1 } else { 0 };

        self.cigar = cigar;
        let cigar_view = self.cigar.clone();

        'process_cigar: {
            let mut ci = 0;
            while ci < self.cigar.len() {
                if self.skip_overlapping_reads(record, self.start, alignment_start, is_reverse, mpos)
                {
                    break 'process_cigar;
                }

                let c = *self.cigar.get(ci).unwrap();
                self.cigar_len = c.len();

                match c {
                    Cigar::RefSkip(_) => {
                        self.process_not_matched(self.cigar_len);
                        ci += 1;
                        continue;
                    }
                    Cigar::SoftClip(_) => {
                        self.process_soft_clip(
                            record,
                            query_seq,
                            mapping_quality,
                            query_qual,
                            nm,
                            is_reverse,
                            pos,
                            read_len_including_softclips,
                            ci,
                            &cigar_view,
                        )?;
                        ci += 1;
                        continue;
                    }
                    Cigar::HardClip(_) => {
                        self.offset = 0;
                        ci += 1;
                        continue;
                    }
                    Cigar::Ins(_) => {
                        self.offset = 0;
                        ci = self.process_insertion(
                            query_seq,
                            mapping_quality,
                            query_qual,
                            nm,
                            is_reverse,
                            pos,
                            read_match_ins_len,
                            ci,
                            &cigar_view,
                        )?;
                        ci += 1;
                        continue;
                    }
                    Cigar::Del(_) => {
                        self.offset = 0;
                        ci = self.process_deletion(
                            query_seq,
                            mapping_quality,
                            query_qual,
                            nm,
                            is_reverse,
                            read_match_ins_len,
                            ci,
                        )?;
                        ci += 1;
                        continue;
                    }
                    _ => {}
                }

                let mut nmoff = 0;
                let mut moffset = 0;

                let mut i = self.offset;
                while i < self.cigar_len as usize {
                    let trim =
                        self.is_trim_at_opt_t_bases(is_reverse, read_len_including_softclips);

                    let ch1 = *query_seq.get_or_err(self.read_pos_including_softclip)?;
                    let mut s = SmallVecBytes::new();
                    s.push(ch1);
                    let mut start_with_deletion = false;

                    if ch1 == b'N' {
                        if instance().conf.include_n_in_total_depth {
                            inc_cnt(&mut self.ref_coverage, self.start, 1);
                        }
                        self.start += 1;
                        self.read_pos_including_softclip += 1;
                        self.read_pos_excluding_softclip += 1;
                        i += 1;
                        continue;
                    }

                    let mut q = *query_qual.get_or_err(self.read_pos_including_softclip)? as f64;
                    let mut qbases = 1;
                    let mut qibases = 0;
                    let mut ss = SmallVecBytes::new();

                    loop {
                        if (self.start + 1) < self.region.start as i64
                            || (self.start + 1) > self.region.end as i64
                            || (i + 1) >= self.cigar_len as usize
                            || q < instance().conf.goodq
                        {
                            break;
                        }

                        let current_base =
                            *query_seq.get_or_err(self.read_pos_including_softclip)?;
                        if !self.reference.has_and_not_equals(self.start, current_base) {
                            break;
                        }
                        if self.reference.has_and_equals(self.start, b'N') {
                            break;
                        }

                        let next_qual =
                            *query_qual.get_or_err(self.read_pos_including_softclip + 1)?;
                        if (next_qual as f64) < instance().conf.goodq + 5.0 {
                            break;
                        }

                        let nuc = *query_seq.get_or_err(self.read_pos_including_softclip + 1)?;
                        if nuc == b'N' {
                            break;
                        }
                        if self.reference.has_and_equals(self.start + 1, b'N') {
                            break;
                        }

                        match self.reference.get(self.start + 1) {
                            Some(ref_base) => {
                                if nuc != ref_base {
                                    ss.push(nuc);
                                    q += next_qual as f64;
                                    qbases += 1;
                                    self.read_pos_including_softclip += 1;
                                    self.read_pos_excluding_softclip += 1;
                                    i += 1;
                                    self.start += 1;
                                    nmoff += 1;
                                } else {
                                    let mut ssn = 0;
                                    for ssi in 1..=instance().conf.vext as usize {
                                        if i + 1 + ssi >= self.cigar_len as usize {
                                            break;
                                        }
                                        if self.read_pos_including_softclip + 1 + ssi
                                            < query_seq.len()
                                            && self.reference.has_and_not_equals(
                                                self.start + 1 + ssi as i64,
                                                *query_seq.get_or_err(
                                                    self.read_pos_including_softclip + 1 + ssi,
                                                )?,
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
                                    if (*query_qual.get_or_err(self.read_pos_including_softclip + ssn)? as f64)
                                        < instance().conf.goodq + 5.0
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
                            None => {
                                break;
                            }
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

                                    for qi in 1..=ins_len {
                                        let idx = self.read_pos_including_softclip + 1 + qi;
                                        q += (*query_qual.get_or_err(idx)?) as f64;
                                        qibases += 1;
                                    }

                                    self.read_pos_including_softclip += ins_len;
                                    self.read_pos_excluding_softclip += ins_len;
                                    ci += 1;
                                }

                                if self.is_next_after_num_matched(ci, 1)? {
                                    let next_len = self
                                        .cigar
                                        .get(ci + 1)
                                        .map(|c| c.len())
                                        .unwrap_or(0) as usize;
                                    if let Some((toffset, tnmoff, tseq, tqual)) =
                                        self.find_offset(
                                            (self.start + ddlen as i64 + 1) as usize,
                                            self.read_pos_including_softclip + 1,
                                            next_len,
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
                                
                                // Java: remove first '&', append insertion, then insert '&' at next_len and prefix '+'
                                let mut merged = SmallVecBytes::new();
                                if let Some(amp_pos) = s.iter().position(|&b| b == b'&') {
                                    merged.extend_from_slice(&s[0..amp_pos]);
                                    merged.extend_from_slice(&s[amp_pos + 1..]);
                                } else {
                                    merged.extend_from_slice(&s);
                                }

                                let mut insertion_seq = SmallVecBytes::new();
                                insertion_seq.extend_from_slice(
                                    query_seq.get_or_err(
                                        (self.read_pos_including_softclip + 1)
                                            ..(self.read_pos_including_softclip + 1 + next_len),
                                    )?,
                                );
                                merged.extend_from_slice(&insertion_seq);

                                let split = next_len.min(merged.len());
                                let mut final_s = SmallVecBytes::new();
                                final_s.push(b'+');
                                final_s.extend_from_slice(&merged[..split]);
                                final_s.push(b'&');
                                if split < merged.len() {
                                    final_s.extend_from_slice(&merged[split..]);
                                }
                                s = final_s;

                                for qi in 1..=next_len {
                                    let idx = self.read_pos_including_softclip + 1 + qi;
                                    q += (*query_qual.get_or_err(idx)?) as f64;
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
                                    self.add_variation_for_matching_part(
                                        mapping_quality,
                                        nm,
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
                            if self.skip_overlapping_reads(record, self.start, alignment_start, is_reverse, mpos) {
                                break 'process_cigar;
                            }

                            i += 1;
                        }

                        if moffset != 0 {
                            self.offset = moffset;
                            self.read_pos_including_softclip += moffset;
                            self.start += moffset as i64;
                            self.read_pos_excluding_softclip += moffset;
                        }

                if self.start > self.region.end as i64 {
                    break;
                }
                
                ci += 1;
            }
        }

        // Buffers are now managed by parse_cigar wrapper, no need to restore them here

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
        num_mismatch: i32,
        is_reverse: bool,
        pos: i64,
        total_length_including_soft_clipped: usize,
        ci: usize,
        cigar: &CigarStringView,
    ) -> Result<(), Error> {
        let mut cigar_len = self.cigar_len;
        let contig = self.region.chrom.clone();
        let max_read_len = self.max_read_len;

        //First record in CIGAR
        if ci == 0 {
            // 5' soft clipped
            // Ignore large soft clip due to chimeric reads in library construction
            if !instance().conf.chimeric_filter {
                if cigar_len >= 20
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
                        self.read_pos_including_softclip += cigar_len as usize;
                        self.offset = 0;

                        // Had to reset the start due to softclipping adjustment
                        self.start = pos;

                        self.cigar_len = cigar_len;
                        return Ok(());
                    }
                    // trying to detect chimeric reads even when there's no supplementary
                    // alignment from aligner
                } else if cigar_len >= Configuration::SEED_1 as u32 {
                    let ref_seed_map = &self.reference.seed;
                    let rev_comp_seq = record
                        .seq()
                        .into_decoded_base_iter()
                        .take(cigar_len as usize)
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
                            && (self.start - *poss.get(0).unwrap()).abs()
                                < 2 * self.max_read_len as i64
                        {
                            self.read_pos_including_softclip += cigar_len as usize;
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

                            self.cigar_len = cigar_len;
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
            let chr_len = *instance()
                .chr_lens
                .get(contig.as_str())
                .unwrap_or(&0) as i64;
            while cigar_len >= 1
                && self.start > 1
                && self.start - 1 <= chr_len
                && self.reference.has_and_equals(
                    self.start - 1,
                    query_sequence
                        .get_or_err(cigar_len as usize - 1)
                        .copied()?,
                )
                    && *query_quality.get_or_err(cigar_len as usize - 1)? > 10
            {
                //create variant if it is not present
                let ref_b = self
                    .reference
                    .get(self.start - 1)
                    .ok_or_else(|| anyhow!("Reference base missing at {}", self.start - 1))?;

                let read_pos = cigar_len as usize;
                let base_qual = query_quality.get_or_err(cigar_len as usize - 1).copied()? as f64;
                {
                    let var =
                        self.get_non_insertion_variant(self.start - 1, &VarDesc::SNV { ref_base: ref_b });
                    //add count
                    // BAM stores Phred quality directly
                    add_cnt(var, is_reverse, read_pos, base_qual, mapq, num_mismatch, None);
                }
                //increase coverage
                inc_cnt(&mut self.ref_coverage, self.start - 1, 1);

                self.start -= 1;
                cigar_len -= 1;
            }

            if cigar_len > 0 {
                //If there remains a soft-clipped sequence at the beginning (not everything was matched)
                let mut read_qual_sum = 0_usize;
                let mut num_high_qual_base = 0;
                let mut num_low_qual_base = 0;

                // Loop over remaining soft-clipped sequence
                for si in (0..cigar_len as usize).rev() {
                    // Stop if unknown base (N - any of ATGC) is found
                    if query_sequence.get_or_err(si).copied()? == b'N' {
                        break;
                    }

                    // base quality - BAM stores Phred quality directly
                    let bq = *query_quality.get_or_err(si)?;
                    if bq <= 12 {
                        num_low_qual_base += 1;
                    }
                    // Stop if more than one low-quality base is found
                    if num_low_qual_base > 1 {
                        break;
                    }

                    read_qual_sum += bq as usize;
                    num_high_qual_base += 1;
                }

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

            cigar_len = cigar.get(ci).unwrap().len();
        } else if ci == cigar.len() - 1 {
            // 3' soft clip
            // Ignore large soft clip due to chimeric reads in library construction
            if !instance().conf.chimeric_filter {
                if cigar_len >= 20
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
                        self.read_pos_including_softclip += cigar_len as usize;
                        self.offset = 0;

                        // Had to reset the start due to softclipping adjustment
                        self.start = pos;

                        self.cigar_len = cigar_len;
                        return Ok(());
                    }
                } else if cigar_len >= Configuration::SEED_1 as u32 {
                    let ref_seed_map = &self.reference.seed;
                    let rev_comp_seq = record
                        .seq()
                        .into_decoded_base_iter()
                        .rev()
                        .take(cigar_len as usize)
                        .map(complement_base)
                        // .take(Configuration::SEED_1 as usize)
                        .collect::<Vec<_>>();

                    if let Some(poss) = ref_seed_map
                        .get(rev_comp_seq.get_with_int(-(Configuration::SEED_1 as i32)..)?)
                    {
                        if poss.len() == 1
                            && (self.start - *poss.get(0).unwrap()).abs()
                                < 2 * self.max_read_len as i64
                        {
                            self.read_pos_including_softclip += cigar_len as usize;
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
                            self.cigar_len = cigar_len;
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
                && self.reference.has_and_equals(
                    self.start,
                    query_sequence
                        .get_or_err(self.read_pos_including_softclip)
                        .copied()?,
                )
                && *query_quality.get_or_err(self.read_pos_including_softclip)? > 10
            {
                //Initialize entry if not present
                let ref_b = self
                    .reference
                    .get(self.start)
                    .ok_or_else(|| anyhow!("Reference base missing at {}", self.start))?;

                let read_pos = total_length_including_soft_clipped - self.read_pos_excluding_softclip;
                let base_qual = query_quality
                    .get_or_err(self.read_pos_including_softclip)
                    .copied()? as f64;
                {
                    let var =
                        self.get_non_insertion_variant(self.start, &VarDesc::SNV { ref_base: ref_b });
                    //add count - BAM stores Phred quality directly
                    add_cnt(var, is_reverse, read_pos, base_qual, mapq, num_mismatch, None);
                }
                // Add coverage
                inc_cnt(&mut self.ref_coverage, self.start, 1);
                self.read_pos_including_softclip += 1;
                self.read_pos_excluding_softclip += 1;
                self.start += 1;
                cigar_len -= 1;
            }

            // If there remains a soft-clipped sequence at the end (not everything was
            // matched) - Java checks this BEFORE advancing read_pos_including_softclip
            if query_sequence.len() - self.read_pos_including_softclip > 0 {
                let mut read_qual_sum = 0;
                let mut num_high_qual_base = 0;
                let mut num_low_qual_base = 0;
                for si in 0..cigar_len {
                    // Loop over remaining soft-clipped sequence
                    // Stop if unknown base (N - any of ATGC) is found

                    // At this point, self.read_pos_including_softclip is start offset of soft clip
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

                // Note: self.start here points to the first soft-clipped position (from while loop)
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
        }

        // Move read position by m (length of segment in CIGAR)
        self.read_pos_including_softclip += cigar_len as usize;
        self.offset = 0;
        self.start = pos; // Had to reset the start due to softclipping adjustment

        self.cigar_len = cigar_len;

        Ok(())
    }

    /// Process CIGAR deletions part. Will ignore indels next to introns and create
    /// Variations for deletions
    fn process_deletion(
        &mut self,
        query_seq: &[u8],
        mapq: u8,
        query_qual: &[u8],
        nm: i32,
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
            )?;

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
        // Anchor position is the deleted segment start (Java uses start)
        let anchor_pos = self.start;

        // Only record the variant if the anchor falls inside the region of interest
        if anchor_pos >= self.region.start as i64 && anchor_pos <= self.region.end as i64 {
            self.add_variation_for_deletion(
                mapq,
                nm,
                is_reverse,
                read_len_including_match_ins,
                &var_desc,
                &qual_seg,
                nmoff,
                anchor_pos,
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
        nm: i32,
        is_reverse: bool,
        read_len_including_match_ins: usize,
        var_desc: &VarDesc,
        qual_seg: &[u8],
        nmoff: usize,
        anchor_pos: i64,
    ) {
        let desc_string = var_desc.to_string();
        Self::increment_position_count(
            &mut self.position_to_deletions_count,
            anchor_pos,
            &desc_string,
        );

        // Minimum of positions from start of read and end of read
        // Java: tp = n < rlen1 - n ? n + 1 : rlen1 - n
        let from_start = self.read_pos_excluding_softclip;
        let from_end = read_len_including_match_ins.saturating_sub(from_start);
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

        {
            // Get or create variation structure for this deletion using the anchor position
            let var = self.get_non_insertion_variant(anchor_pos, var_desc);

            // Increment direction count
            var.inc_dir(is_reverse);

            // Increase variant count
            var.alt_depth += 1;

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
            var.nm += (nm - nmoff as i32) as f64;

            if tmpq >= instance().conf.goodq {
                var.high_qual_read_cnt += 1;
            } else {
                var.low_qual_read_cnt += 1;
            }
        }

        // Increase coverage count for reference bases missing from the read
        for i in 0..self.cigar_len as i64 {
            inc_cnt(&mut self.ref_coverage, self.start + i, 1);
        }
    }

    /// Check if position should be trimmed at opt_T bases
    fn is_trim_at_opt_t_bases(&self, is_reverse: bool, total_length_including_softclipped: usize) -> bool {
        if instance().conf.trim_bases_after != 0 {
            let trim_bases = instance().conf.trim_bases_after as usize;
            if !is_reverse {
                return self.read_pos_including_softclip > trim_bases;
            } else {
                return total_length_including_softclipped
                    .saturating_sub(self.read_pos_including_softclip)
                    > trim_bases;
            }
        }
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
        if !instance().conf.perform_local_realignment {
            return Ok(false);
        }
        // Do not adjust complex if we have hard-clips after insertion/deletion
        if self.cigar.len() > ci + 2 {
            if matches!(self.cigar.get(ci + 2), Some(Cigar::HardClip(_))) {
                return Ok(false);
            }
        }
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
                let has_good_quality = (*query_qual.get_or_err(self.read_pos_including_softclip)? as f64)
                    >= instance().conf.goodq;

                return Ok(has_mismatch && has_good_quality);
            }
        }
        Ok(false)
    }

    /// Find mismatches in next segment with good quality
    fn find_offset(
        &mut self,
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
            for osi in 0..offset {
                inc_cnt(&mut self.ref_coverage, (ts + osi) as i64, 1);
            }
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
        nm: i32,
        is_reverse: bool,
        pos: i64,
        read_len_including_match_ins: usize,
        mut ci: usize,
        cigar: &CigarStringView,
    ) -> Result<usize, Error> {
        let trace_ins = env::var("VARDICT_TRACE_INS_QNAME")
            .ok()
            .map(|needle| {
                needle == "*"
                    || self
                        .current_qname
                        .as_deref()
                        .map_or(false, |name| name == needle)
            })
            .unwrap_or(false);

        let ins_len = self.cigar_len as usize;

        if trace_ins {
            event!(
                Level::WARN,
                "[trace_ins] enter qname={} ci={} start={} read_pos_incl={} read_pos_excl={} ins_len={} cigar={}",
                self.current_qname.as_deref().unwrap_or("-"),
                ci,
                self.start,
                self.read_pos_including_softclip,
                self.read_pos_excluding_softclip,
                ins_len,
                self.last_modified_cigar.as_deref().unwrap_or("-")
            );
        }

        // Ignore insertions right after introns at exon edge in RNA-seq
        if skip_indel_next_to_intron(&self.cigar, ci)? {
            self.read_pos_including_softclip += ins_len;
            return Ok(ci);
        }

        // Inserted segment of read sequence
        let mut desc = SmallVecBytes::new();
        desc.extend_from_slice(query_seq.get_or_err(
            self.read_pos_including_softclip..self.read_pos_including_softclip + ins_len,
        )?);

        // Quality of this segment
        let mut qual_seg = SmallVecBytes::new();
        qual_seg.extend_from_slice(query_qual.get_or_err(
            self.read_pos_including_softclip..self.read_pos_including_softclip + ins_len,
        )?);

        // Sequence to be appended if next segment is matched
        let mut ss = SmallVecBytes::new();

        // For multiple indels within VEXT bp
        let mut multoffs = 0usize;
        let mut multoffp = 0usize;
        let mut nmoff = 0usize;
        self.offset = 0;

        if is_followed_by_match_and_indel(&self.cigar, ci) {
            let mlen = cigar.get(ci + 1).unwrap().len() as usize;
            let indel_len = cigar.get(ci + 2).unwrap().len() as usize;
            let begin = self.read_pos_including_softclip + ins_len;

            // Append matched segment after insertion
            desc.push(b'#');
            desc.extend_from_slice(query_seq.get_or_err(begin..begin + mlen)?);
            qual_seg.extend_from_slice(query_qual.get_or_err(begin..begin + mlen)?);

            // Append '^' + next-next segment (sequence for insertion, length for deletion)
            desc.push(b'^');
            if matches!(cigar.get(ci + 2), Some(Cigar::Ins(_))) {
                desc.extend_from_slice(query_seq.get_or_err(begin + mlen..begin + mlen + indel_len)?);
                qual_seg.extend_from_slice(query_qual.get_or_err(begin + mlen..begin + mlen + indel_len)?);
            } else {
                let len_str = indel_len.to_string();
                desc.extend_from_slice(len_str.as_bytes());
                qual_seg.push(*query_qual.get_or_err(begin + mlen)?);
            }

            // add length of next segment to both multoffs and multoffp
            multoffs += mlen
                + if matches!(cigar.get(ci + 2), Some(Cigar::Del(_))) {
                    indel_len
                } else {
                    0
                };
            multoffp += mlen
                + if matches!(cigar.get(ci + 2), Some(Cigar::Ins(_))) {
                    indel_len
                } else {
                    0
                };

            if is_next_after_num_matched(cigar, ci, 3) {
                let seg_len = cigar.get(ci + 3).unwrap().len() as usize;
                if let Some((offset, tnm, seq, qual)) = self.find_offset(
                    self.start as usize + multoffs,
                    self.read_pos_including_softclip + ins_len + multoffp,
                    seg_len,
                    query_seq,
                    query_qual,
                )? {
                    self.offset = offset;
                    nmoff += tnm;
                    ss = seq;
                    qual_seg.extend_from_slice(&qual);
                }
            }
            ci += 2;
        } else if is_next_matched(cigar, ci) {
            let mlen = cigar.get(ci + 1).unwrap().len() as usize;
            let mut vsn = 0;
            let mut vi = 0;
            while vsn <= instance().conf.vext as usize && vi < mlen {
                let read_idx = self.read_pos_including_softclip + ins_len + vi;
                let seq_ch = *query_seq.get_or_err(read_idx)?;
                if seq_ch == b'N' {
                    break;
                }
                let qual = *query_qual.get_or_err(read_idx)? as f64;
                if qual < instance().conf.goodq {
                    break;
                }
                match self.reference.get(self.start + vi as i64) {
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
                let start_idx = self.read_pos_including_softclip + ins_len;
                ss.extend_from_slice(query_seq.get_or_err(start_idx..start_idx + self.offset)?);
                qual_seg.extend_from_slice(query_qual.get_or_err(start_idx..start_idx + self.offset)?);
                for osi in 0..self.offset {
                    inc_cnt(&mut self.ref_coverage, self.start + osi as i64, 1);
                }
            }
        }

        if self.offset > 0 {
            desc.push(b'&');
            desc.extend_from_slice(&ss);
        }

        // If start-1 is in region of interest and no N
        let mut insertion_pos = self.start - 1;
        if insertion_pos >= self.region.start as i64
            && insertion_pos <= self.region.end as i64
            && !desc.iter().any(|&b| b == b'N')
        {
            if is_begin_atgc_end(&desc) {
                let (adj_pos, adj_seq) = adj_ins_pos(insertion_pos, &desc, &self.reference);
                let idx = (self.read_pos_including_softclip as i64 - 1)
                    - (self.start - 1 - adj_pos);
                if idx > 0 {
                    insertion_pos = adj_pos;
                    desc = adj_seq;
                }
            }

            let desc_string = format!("+{}", String::from_utf8_lossy(desc.as_slice()));
            if trace_ins {
                event!(
                    Level::WARN,
                    "[trace_ins] emit qname={} insertion_pos={} desc={} offset={} multoffs={} multoffp={} nmoff={}",
                    self.current_qname.as_deref().unwrap_or("-"),
                    insertion_pos,
                    desc_string,
                    self.offset,
                    multoffs,
                    multoffp,
                    nmoff
                );
            }
            Self::increment_position_count(
                &mut self.position_to_insertion_count,
                insertion_pos,
                &desc_string,
            );

            let var_desc = VarDesc::Ins { seq: desc.clone() };
            let var = get_variants_from_map(&mut self.insertion_vars, insertion_pos, &var_desc);

            let tmpq = if qual_seg.is_empty() {
                0.0
            } else {
                qual_seg.iter().map(|&q| q as f64).sum::<f64>() / qual_seg.len() as f64
            };

            add_cnt(
                var,
                is_reverse,
                self.read_pos_excluding_softclip,
                tmpq,
                mapq,
                nm - nmoff as i32,
                Some(read_len_including_match_ins),
            );

            // Adjust the reference count for insertion reads
            let index_in_query = (self.read_pos_including_softclip as i64 - 1)
                - (self.start - 1 - insertion_pos);
            if insertion_pos > pos && index_in_query >= 0 {
                let idx = index_in_query as usize;
                if idx < query_seq.len() {
                    let read_base = query_seq[idx].to_ascii_uppercase();
                    if self.reference.has_and_equals(insertion_pos, read_base) {
                        if let Some(vars_at) = self.non_insertion_vars.get_mut(&insertion_pos) {
                            let key = VarDesc::snv_key(read_base);
                            if let Some(ref_var) = vars_at.get_mut(&key) {
                                let base_qual = query_qual[idx] as f64;
                                sub_cnt(
                                    ref_var,
                                    is_reverse,
                                    self.read_pos_excluding_softclip,
                                    base_qual,
                                    mapq,
                                    nm - nmoff as i32,
                                    Some(read_len_including_match_ins),
                                );
                            }
                        }
                    }
                }
            }

            // Adjust count if the insertion is at the edge so that the AF won't > 1
            if ci == 1
                && matches!(
                    cigar.get(0).copied(),
                    Some(Cigar::SoftClip(_) | Cigar::HardClip(_))
                )
            {
                if let Some(ref_base) = self.reference.get(insertion_pos) {
                    let ref_var = get_variants_from_map(
                        &mut self.non_insertion_vars,
                        insertion_pos,
                        &VarDesc::snv_key(ref_base),
                    );
                    ref_var.inc_dir(is_reverse);
                    ref_var.alt_depth += 1;

                    let from_start = self.read_pos_excluding_softclip;
                    let from_end = read_len_including_match_ins
                        .saturating_sub(self.read_pos_excluding_softclip);
                    let tp = if from_start < from_end {
                        from_start + 1
                    } else {
                        from_end
                    };

                    if !ref_var.pstd && ref_var.pp != 0 && tp != ref_var.pp {
                        ref_var.pstd = true;
                    }
                    if !ref_var.qstd && ref_var.pq != 0.0 && (tmpq - ref_var.pq).abs() > f64::EPSILON {
                        ref_var.qstd = true;
                    }
                    ref_var.mean_pos += tp as f64;
                    ref_var.mean_qual += tmpq;
                    ref_var.mean_mapq += mapq as f64;
                    ref_var.pp = tp;
                    ref_var.pq = tmpq;
                    ref_var.nm += (nm - nmoff as i32) as f64;
                    inc_cnt(&mut self.ref_coverage, insertion_pos, 1);
                }
            }
        }

        if trace_ins {
            event!(
                Level::WARN,
                "[trace_ins] exit qname={} start={} read_pos_incl={} read_pos_excl={} ci={}",
                self.current_qname.as_deref().unwrap_or("-"),
                self.start,
                self.read_pos_including_softclip,
                self.read_pos_excluding_softclip,
                ci
            );
        }

        // adjust read position by m (CIGAR segment length) + offset + multoffp
        self.read_pos_including_softclip += ins_len + self.offset + multoffp;
        self.read_pos_excluding_softclip += ins_len + self.offset + multoffp;
        // adjust reference position by offset + multoffs
        self.start += self.offset as i64 + multoffs as i64;

        Ok(ci)
    }

    /// Add variation record for a matched segment (SNV or MNV)
    fn add_variation_for_matching_part(
        &mut self,
        mapq: u8,
        nm: i32,
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
        let debug_pos = env::var("VARDICT_DEBUG_POS")
            .ok()
            .and_then(|v| v.parse::<i64>().ok());
        let trace_merge_ins = env::var("VARDICT_TRACE_MERGED_INS")
            .ok()
            .map(|v| v == "1" || v.eq_ignore_ascii_case("true"))
            .unwrap_or(false);

        // Build VarDesc from s string
        // s format examples: "A", "A&TGC" (MNV), "+ATC" (insertion), "-2&AT" (deletion+match), etc.
        
        if s.is_empty() || s.iter().all(|&b| b == b'N') {
            return;
        }

        // Average quality
        let avg_qual = if qbases + qibases > 0 {
            q / (qbases + qibases) as f64
        } else {
            q
        };
        
        let nm_adjusted = nm - nmoff as i32;
        let mut did_add_variant = false;
        let s_str = String::from_utf8_lossy(s);

        // Check if this is an insertion first
        if s.starts_with(b"+") {
            // Insertion: +ATC format
            let ins_seq: SmallVec<[u8; 32]> = s[1..].iter().copied().collect();

            if trace_merge_ins {
                event!(
                    Level::WARN,
                    "[trace_merged_ins] qname={} pos={} s={} read_pos={} start={} read_pos_incl={} read_pos_excl={}",
                    self.current_qname.as_deref().unwrap_or("-"),
                    pos,
                    s_str,
                    read_pos,
                    self.start,
                    self.read_pos_including_softclip,
                    self.read_pos_excluding_softclip
                );
            }

            let var_desc = VarDesc::Ins { seq: ins_seq };

            Self::increment_position_count(
                &mut self.position_to_insertion_count,
                pos,
                &s_str,
            );

            // Store insertions in insertion_vars (Java: addVariationForMatchingPart uses insertionVariants)
            // use read_pos for tp calculation
            let var = get_variants_from_map(&mut self.insertion_vars, pos, &var_desc);
            add_cnt(
                var,
                is_reverse,
                read_pos,
                avg_qual,
                mapq,
                nm_adjusted,
                Some(read_len_including_match_ins),
            );
            did_add_variant = true;
        } else if s.starts_with(b"-") {
            // Deletion from matching part: preserve raw description string
            let var_desc = VarDesc::Raw { desc: s.to_vec().into() };

            let var = self.get_non_insertion_variant(pos, &var_desc);
            add_cnt(
                var,
                is_reverse,
                read_pos,
                avg_qual,
                mapq,
                nm_adjusted,
                Some(read_len_including_match_ins),
            );
            did_add_variant = true;
        } else if s.iter().any(|&b| b == b'&') {
            // Complex/raw variant (MNV with '&')
            let var_desc = VarDesc::Raw { desc: s.to_vec().into() };

            let var = self.get_non_insertion_variant(pos, &var_desc);
            add_cnt(
                var,
                is_reverse,
                read_pos,
                avg_qual,
                mapq,
                nm_adjusted,
                Some(read_len_including_match_ins),
            );
            did_add_variant = true;
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

            if pos == 76962 {
                event!(
                    Level::DEBUG,
                    "[SNV-TRACE] pos={} qname={} alt={} ref={} read_pos={} avg_qual={:.2} mapq={} nm={} is_reverse={}",
                    pos,
                    self.current_qname.as_deref().unwrap_or("-"),
                    alt_base as char,
                    ref_base as char,
                    read_pos,
                    avg_qual,
                    mapq,
                    nm_adjusted,
                    is_reverse
                );
            }
            
            // Get or create variant and add count - use read_pos for tp calculation
            let var = self.get_non_insertion_variant(pos, &var_desc);
            add_cnt(
                var,
                is_reverse,
                read_pos,
                avg_qual,
                mapq,
                nm_adjusted,
                Some(read_len_including_match_ins),
            );
            did_add_variant = true;
        } else {
            // Unknown format - skip
        }

        if did_add_variant {
            let shift = if s.starts_with(b"+") && s.iter().any(|&b| b == b'&') {
                1
            } else {
                0
            };

            let covered_bases = qbases.saturating_sub(shift);
            for qi in 1..=covered_bases {
                inc_cnt(&mut self.ref_coverage, self.start - qi as i64 + 1, 1);
            }

            if start_with_deletion {
                Self::increment_position_count(
                    &mut self.position_to_deletions_count,
                    pos,
                    &s_str,
                );
                for qi in 1..ddlen {
                    inc_cnt(&mut self.ref_coverage, self.start + qi as i64, 1);
                }
            }
        }

        if debug_pos == Some(pos) && s.starts_with(b"C&") {
            event!(
                Level::DEBUG,
                "[cigar_debug_pos] qname={} pos={} s={} read_pos_excl={} read_pos_incl={} start={} cigar={}",
                self.current_qname.as_deref().unwrap_or("-"),
                pos,
                String::from_utf8_lossy(s),
                read_pos,
                self.read_pos_including_softclip,
                self.start,
                self.last_modified_cigar.as_deref().unwrap_or("-")
            );
        }

        if is_begin_atgc_amp_atgcs_end(s) {
            let desc = String::from_utf8_lossy(s).to_string();
            let pos_map = self.mnp.entry(pos).or_insert_with(HashMap::new);
            *pos_map.entry(desc).or_insert(0) += 1;
        }
    }

    fn increment_position_count(
        map: &mut HashMap<i64, HashMap<String, usize>>,
        pos: i64,
        key: &str,
    ) {
        let pos_map = map.entry(pos).or_insert_with(HashMap::new);
        *pos_map.entry(key.to_string()).or_insert(0) += 1;
    }

    /// N in CIGAR - skipped region from reference
    /// Skip the region and add string start-end to %SPLICE
    fn process_not_matched(&mut self, cigar_len: u32) {
        let key = (self.start - 1, self.start + cigar_len as i64 - 1);

        if !self.splice_count.contains_key(&key) {
            self.splice_count_insert_index
                .insert(key, self.next_splice_count_insert_index);
            self.next_splice_count_insert_index += 1;
        }

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
        let filter_bp = instance().conf.crispr_filtering_bp;

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
        current_start: i64,
        pos: i64,
        direction: bool, // true is reverse direction
        mate_pos: i64,
    ) -> bool {
        if instance().conf.unique_mode_alignment_enabled
            && is_paired_and_same_chromosome(record)
            && !direction
            && current_start >= mate_pos
        {
            return true;
        }

        let is_second_in_pair = (record.flags() & 0x80) != 0;
        if instance().conf.unique_mode_second_in_pair_enabled
            && is_second_in_pair
            && is_paired_and_same_chromosome(record)
            && are_reads_overlap(record, current_start, pos, mate_pos)
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
        num_mismatch: i32,
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

                let b = query_sequence.get_or_err(si as usize).copied()?.to_ascii_uppercase();
                let idx = cigar_len - 1 - si; // distance from start of match.
                let cnts = sclip
                    .nt
                    .entry(idx as i64)
                    .or_insert_with(|| NucBaseMap::default());

                // increase count of current base (skip if not A/T/C/G/N)
                // Use get_or_insert_with to initialize the value if not present
                if let Some(cnt) = cnts.get_or_insert_with(b, || 0) {
                    *cnt += 1;
                }

                let seq_var = get_variation_from_seq(sclip, idx as usize, b);
                // BAM stores Phred quality directly
                add_cnt(
                    seq_var,
                    is_reverse,
                    si as usize - (cigar_len as usize - num_high_qual_base),
                    *query_quality.get_or_err(si as usize)? as f64,
                    mapq,
                    num_mismatch,
                    None,
                );
            }

            add_cnt(
                &mut sclip.var,
                is_reverse,
                cigar_len as usize,
                read_qual_sum as f64 / num_high_qual_base as f64,
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
        num_mismatch: i32,
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
                    .copied()?
                    .to_ascii_uppercase();
                let idx = si; // distance from start of match.
                let cnts = sclip
                    .nt
                    .entry(idx as i64)
                    .or_insert_with(|| NucBaseMap::default());

                // increase count of current base (skip if not A/T/C/G/N)
                // Use get_or_insert_with to initialize the value if not present
                if let Some(cnt) = cnts.get_or_insert_with(b, || 0) {
                    *cnt += 1;
                }

                let seq_var = get_variation_from_seq(sclip, idx as usize, b);
                // BAM stores Phred quality directly
                add_cnt(
                    seq_var,
                    is_reverse,
                    num_high_qual_base - si,
                    *query_quality.get_or_err(self.read_pos_including_softclip + si)? as f64,
                    mapq,
                    num_mismatch,
                    None,
                );
            }

            add_cnt(
                &mut sclip.var,
                is_reverse,
                cigar_len as usize,
                read_qual_sum as f64 / num_high_qual_base as f64,
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
fn get_aligned_length(cigar: &CigarStringView) -> i64 {
    cigar
        .iter()
        .map(|c| match c {
            Cigar::Match(l) | Cigar::Del(l) => *l as i64,
            _ => 0,
        })
        .sum::<i64>()
}

#[inline]
fn get_aligned_length_mnd(cigar: &CigarStringView) -> i64 {
    cigar
        .iter()
        .map(|c| match c {
            Cigar::Match(l) | Cigar::Del(l) | Cigar::RefSkip(l) => *l as i64,
            _ => 0,
        })
        .sum::<i64>()
}

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
fn are_reads_overlap(record: &Record, current_start: i64, pos: i64, mate_pos: i64) -> bool {
    if pos >= mate_pos {
        let ref_len = record
            .cigar()
            .iter()
            .map(|c| if c.consumes_reference_bases() { c.len() } else { 0 })
            .sum::<u32>() as i64;
        current_start >= mate_pos && current_start <= mate_pos + ref_len - 1
    } else {
        current_start >= mate_pos && record.mpos() + 1 <= record.reference_end()
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

    let is_chimeric_with_sa = ((is_reverse && sa_dir_is_forward)
        || (!is_reverse && !sa_dir_is_forward))
        && sa_chr == record.contig()
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
fn add_cnt(var: &mut Variant, is_reverse: bool, read_pos: usize, bq: f64, mapq: u8, nm: i32, read_len: Option<usize>) {
    var.alt_depth += 1;
    var.inc_dir(is_reverse);

    if let Some(rlen) = read_len {
        // Java addVariationForMatchingPart/addVariationForDeletion logic
        let from_start = read_pos;
        let from_end = rlen.saturating_sub(read_pos);
        let tp = if from_start < from_end {
            from_start + 1
        } else {
            from_end
        };

        let tmpq = bq;

        if !var.pstd && var.pp != 0 && tp != var.pp {
            var.pstd = true;
        }

        if !var.qstd && var.pq != 0.0 && (tmpq - var.pq).abs() > f64::EPSILON {
            var.qstd = true;
        }

        var.mean_pos += tp as f64;
        var.mean_qual += tmpq;
        var.mean_mapq += mapq as f64;
        var.pp = tp;
        var.pq = tmpq;
        var.nm += nm as f64;

        if tmpq >= instance().conf.goodq {
            var.high_qual_read_cnt += 1;
        } else {
            var.low_qual_read_cnt += 1;
        }
    } else {
        // Java addCnt logic (no pstd/qstd updates)
        var.mean_pos += read_pos as f64;
        var.mean_qual += bq;
        var.mean_mapq += mapq as f64;
        var.nm += nm as f64;

        if bq >= instance().conf.goodq {
            var.high_qual_read_cnt += 1;
        } else {
            var.low_qual_read_cnt += 1;
        }
    }
}

/// Increment variant counters without adjusting high/low quality counts.
/// Used for edge insertion adjustment to match Java's ref-call metrics.
fn add_cnt_no_qual(var: &mut Variant, is_reverse: bool, read_pos: usize, bq: f64, mapq: u8, nm: i32, read_len: Option<usize>) {
    var.alt_depth += 1;
    var.inc_dir(is_reverse);

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

    if !var.pstd && var.pp != 0 && tp != var.pp {
        var.pstd = true;
    }

    let tmpq = bq;
    if !var.qstd && var.pq != 0.0 && (tmpq - var.pq).abs() > f64::EPSILON {
        var.qstd = true;
    }

    var.mean_pos += tp as f64;
    var.mean_qual += tmpq;
    var.mean_mapq += mapq as f64;
    var.pp = tp;
    var.pq = tmpq;
    var.nm += nm as f64;
}

/// Decrement variant counters (Java subCnt equivalent).
/// read_pos: position in read (excluding soft clips)
/// bq: base quality
/// read_len: optional total read length for calculating tp correctly
fn sub_cnt(var: &mut Variant, is_reverse: bool, read_pos: usize, bq: f64, mapq: u8, nm: i32, read_len: Option<usize>) {
    if var.alt_depth > 0 {
        var.alt_depth = var.alt_depth.saturating_sub(1);
    }

    if is_reverse {
        var.alt_depth_rev = var.alt_depth_rev.saturating_sub(1);
    } else {
        var.alt_depth_fwd = var.alt_depth_fwd.saturating_sub(1);
    }

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

    var.mean_pos -= tp as f64;
    var.mean_qual -= bq;
    var.mean_mapq -= mapq as f64;
    var.nm -= nm as f64;

    if bq >= instance().conf.goodq {
        var.high_qual_read_cnt = var.high_qual_read_cnt.saturating_sub(1);
    } else {
        var.low_qual_read_cnt = var.low_qual_read_cnt.saturating_sub(1);
    }
}

/// Increment variant counters for insertion/deletion anchor positions.
/// These don't count toward high quality reads.
/// read_pos: position in read (excluding soft clips)  
/// bq: base quality
/// read_len: optional total read length for calculating tp correctly
fn add_cnt_anchor(var: &mut Variant, is_reverse: bool, read_pos: usize, bq: f64, mapq: u8, nm: i32, read_len: Option<usize>) {
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
    
    // pstd: true if variant is covered by reads with different positions
    if !var.pstd && var.pp != 0 && tp != var.pp {
        var.pstd = true;
    }
    
    // qstd: true if variant is covered by reads with different qualities
    let tmpq = bq;
    if !var.qstd && var.pq != 0.0 && (tmpq - var.pq).abs() > f64::EPSILON {
        var.qstd = true;
    }
    
    var.mean_pos += tp as f64;
    var.mean_qual += tmpq;
    var.mean_mapq += mapq as f64;
    var.pp = tp;
    var.pq = tmpq;
    var.nm += nm as f64;
    // Note: Don't increment high_qual_read_cnt or low_qual_read_cnt for anchors
}

fn format_i64_usize_map(map: &HashMap<i64, usize>) -> String {
    let mut keys: Vec<i64> = map.keys().copied().collect();
    keys.sort_unstable();
    let mut out = String::from("{");
    for (i, k) in keys.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        let v = map.get(k).copied().unwrap_or(0);
        out.push_str(&format!("{}:{}", k, v));
    }
    out.push('}');
    out
}

fn format_var_desc(desc: &VarDesc) -> String {
    match desc {
        VarDesc::SNV { ref_base } => format!("SNV({})", *ref_base as char),
        VarDesc::Ins { seq } => format!("INS({})", String::from_utf8_lossy(seq.as_slice())),
        VarDesc::Del { len, match_seq, ins_or_del_len, mismatch_seq } => {
            let match_s = String::from_utf8_lossy(match_seq.as_slice());
            let mismatch_s = String::from_utf8_lossy(mismatch_seq.as_slice());
            let iod = match ins_or_del_len {
                InsOrDelLen::None => "None".to_string(),
                InsOrDelLen::InsSeq(s) => format!("InsSeq({})", String::from_utf8_lossy(s.as_slice())),
                InsOrDelLen::DelLen(l) => format!("DelLen({})", l),
            };
            format!("DEL(len={},match={},ins_or_del_len={},mismatch={})", len, match_s, iod, mismatch_s)
        }
        VarDesc::Complex { ref_seq, alt_seq } => format!(
            "COMPLEX(ref={},alt={})",
            String::from_utf8_lossy(ref_seq.as_slice()),
            String::from_utf8_lossy(alt_seq.as_slice())
        ),
        VarDesc::Raw { desc } => format!("RAW({})", String::from_utf8_lossy(desc.as_slice())),
    }
}

fn format_variant(var: &Variant) -> String {
    format!(
        "{{alt_depth:{},alt_depth_fwd:{},alt_depth_rev:{},extra_cnt:{},mean_pos:{},mean_qual:{},mean_mapq:{},nm:{},low_qual_read_cnt:{},high_qual_read_cnt:{},pstd:{},qstd:{},pp:{},pq:{}}}",
        var.alt_depth,
        var.alt_depth_fwd,
        var.alt_depth_rev,
        var.extra_cnt,
        var.mean_pos,
        var.mean_qual,
        var.mean_mapq,
        var.nm,
        var.low_qual_read_cnt,
        var.high_qual_read_cnt,
        var.pstd,
        var.qstd,
        var.pp,
        var.pq
    )
}

fn format_variation_pos_map(map: &HashMap<i64, HashMap<VarDesc, Variant>>) -> String {
    let mut pos_keys: Vec<i64> = map.keys().copied().collect();
    pos_keys.sort_unstable();
    let mut out = String::from("{");
    for (pi, pos) in pos_keys.iter().enumerate() {
        if pi > 0 {
            out.push(',');
        }
        let mut entries: Vec<(String, &Variant)> = map
            .get(pos)
            .map(|inner| inner.iter().map(|(k, v)| (format_var_desc(k), v)).collect())
            .unwrap_or_default();
        entries.sort_by(|a, b| a.0.cmp(&b.0));
        out.push_str(&format!("{}:{{", pos));
        for (i, (k, v)) in entries.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            out.push_str(&format!("\"{}\":{}", k, format_variant(v)));
        }
        out.push('}');
    }
    out.push('}');
    out
}


fn format_splice_count(map: &HashMap<(i64, i64), Vec<usize>>) -> String {
    let mut keys: Vec<(i64, i64)> = map.keys().copied().collect();
    keys.sort_unstable();
    let mut out = String::from("{");
    for (i, k) in keys.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        let vals = map.get(k).map(|v| v.clone()).unwrap_or_default();
        out.push_str(&format!("\"{}-{}\":{:?}", k.0, k.1, vals));
    }
    out.push('}');
    out
}

/// Increase count for given key
#[inline]
fn inc_cnt(coverage_map: &mut HashMap<i64, usize>, pos: i64, depth: usize) {
    coverage_map
        .entry(pos)
        .and_modify(|v| v.add_assign(depth))
        .or_insert_with(|| depth);
}

fn cleanup_cigar_view(cigar: &CigarStringView, pos: i64) -> CigarStringView {
    let mut elems: Vec<Cigar> = cigar.iter().copied().collect();

    // Leading elements.
    let mut idx = 0;
    let mut no_matches_yet = true;
    while idx < elems.len() && no_matches_yet {
        match elems[idx] {
            Cigar::Ins(len) => {
                elems[idx] = Cigar::SoftClip(len);
            }
            Cigar::HardClip(_) => {
                elems.remove(idx);
                continue;
            }
            ref op if consumes_read_and_ref(op) => {
                no_matches_yet = false;
            }
            _ => {}
        }
        idx += 1;
    }

    // Trailing elements.
    let mut idx = elems.len();
    let mut no_matches_yet = true;
    while idx > 0 && no_matches_yet {
        idx -= 1;
        match elems[idx] {
            Cigar::Ins(len) => {
                elems[idx] = Cigar::SoftClip(len);
            }
            Cigar::HardClip(_) => {
                elems.remove(idx);
                continue;
            }
            ref op if consumes_read_and_ref(op) => {
                no_matches_yet = false;
            }
            _ => {}
        }
    }

    CigarString(elems).into_view(pos)
}

fn consumes_read_and_ref(op: &Cigar) -> bool {
    matches!(op, Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_))
}

/// Java parity: getCigarOperator(Cigar, ci)
/// Treat insertions at the first/last CIGAR element as soft-clipping.
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
    if !instance().conf.perform_local_realignment {
        return false;
    }
    if ci + 3 >= cigar.len() {
        return false;
    }

    let n_cigar = cigar.get(ci + 1).unwrap();
    let nn_cigar = cigar.get(ci + 2).unwrap();

    matches!(n_cigar, &Cigar::Match(l) if l <= instance().conf.vext as u32)
        && matches!(nn_cigar, Cigar::Ins(_) | Cigar::Del(_))
        && !matches!(cigar.get(ci + 3), Some(Cigar::Ins(_) | Cigar::Del(_)))
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
        return Err(anyhow!("Not deletion description: {:?}", var_desc));
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
    let next_next_cigar = cigar
        .get(ci + 2)
        .ok_or_else(|| anyhow!("Missing ci+2 CIGAR segment at index {}", ci + 2))?;

    if matches!(next_next_cigar, Cigar::Ins(_)) {
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
    if matches!(next_next_cigar, Cigar::Ins(_)) {
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

/// Check sequence matches pattern: ^[ATGC]+$
fn is_begin_atgc_end(sequence: &[u8]) -> bool {
    if sequence.is_empty() {
        return false;
    }
    sequence.iter().all(|&c| is_atgc(c))
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

/// Adjust insertion position and sequence (Java: VariationRealigner.adjInsPos)
fn adj_ins_pos(mut bi: i64, ins: &[u8], reference: &Reference) -> (i64, SmallVec<[u8; 32]>) {
    let len = ins.len();
    if len == 0 {
        return (bi, SmallVec::new());
    }

    let mut n = 1usize;
    let mut adjusted = ins.to_vec();

    loop {
        let ref_base = reference.get(bi);
        let ins_base = adjusted.get(len - n).copied();
        let is_eq = match (ref_base, ins_base) {
            (Some(r), Some(b)) => r == b,
            (None, _) => false,
            (_, None) => false,
        };
        if !is_eq {
            break;
        }
        n += 1;
        if n > len {
            n = 1;
        }
        bi -= 1;
    }

    if n > 1 {
        adjusted.rotate_right(n - 1);
    }

    (bi, adjusted.into())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::CigarString;

    /// Helper to create a CigarStringView from CIGAR elements
    fn make_cigar(elements: Vec<Cigar>) -> CigarStringView {
        CigarString(elements).into_view(0)
    }

    #[test]
    fn test_append_segments_returns_error_for_non_deletion_var_desc() {
        let query_seq = b"ACGT";
        let query_qual = b"!!!!";
        let cigar = make_cigar(vec![Cigar::Del(1), Cigar::Match(1), Cigar::Ins(1)]);
        let mut var_desc = VarDesc::Ins {
            seq: SmallVec::new(),
        };
        let mut qual_seg = SmallVecBytes::new();

        let result = append_segments(
            query_seq,
            query_qual,
            &cigar,
            0,
            &mut var_desc,
            &mut qual_seg,
            0,
            1,
            1,
            false,
        );

        assert!(result.is_err());
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

    /// Test getAlignedLength (Java parity) - calculates sum of M + D lengths
    #[test]
    fn test_get_aligned_length() {
        let cigar = make_cigar(vec![
            Cigar::Match(1),
            Cigar::SoftClip(2),
            Cigar::Ins(4),
            Cigar::Del(8),
            Cigar::RefSkip(16),
            Cigar::HardClip(32),
        ]);

        assert_eq!(get_aligned_length(&cigar), 9); // 1 + 8
    }

    /// Test getAlignedLengthMND (Java ALIGNED_LENGTH_MND parity) - calculates sum of M + N + D lengths
    #[test]
    fn test_get_aligned_length_mnd() {
        let cigar = make_cigar(vec![
            Cigar::Match(1),
            Cigar::SoftClip(2),
            Cigar::Ins(4),
            Cigar::Del(8),
            Cigar::RefSkip(16),
            Cigar::HardClip(32),
        ]);

        assert_eq!(get_aligned_length_mnd(&cigar), 25); // 1 + 8 + 16
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

    #[test]
    fn test_cigar_parser_mapped_read_ref_loaded_no_indels() {
        use crate::data::reference::FastaReader;
        use crate::data::region::Region;
        use crate::scopedata::global_read_only_scope::{GlobalReadOnlyScope, INSTANCE, instance};
        use rust_htslib::bam::{Read, Reader};
        use std::sync::Arc;

        let conf = Configuration::default();

        let _ = INSTANCE.set(GlobalReadOnlyScope {
            conf,
            ..Default::default()
        });

        let bam_path = "/home/eck/workspace/vardict_rs/test_data/test_168714.bam";
        let fasta_path = "/home/eck/workspace/vardict_rs/VarDictJava/tests/integration/reference/hs37d5.fa";

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
        let alignment_start = record.pos() + 1;
        let read_len = record.seq_len() as i64;

        let region = Region::new("20".to_string(), 168600, 168800, "test_region".to_string());
        let mut ref_start = region.start.saturating_sub(1200);
        if ref_start == 0 {
            ref_start = 1;
        }
        let ref_end = region.end + 1200;

        let fasta = FastaReader::open(fasta_path).expect("Failed to open reference FASTA");
        let reference = fasta
            .get_reference(region.chr(), ref_start, ref_end)
            .expect("Failed to fetch reference sequence");

        let instance = Arc::new(instance().clone());
        let mut parser = CigarParser::new(region, reference, instance);

        let mut records = vec![record];
        parser
            .process_records(records.iter_mut())
            .expect("process_records failed");

        println!("=== Rust CigarParser non_insertion_vars (mapped read) ===");
        println!(
            "Alignment start: {}, read length: {}",
            alignment_start,
            read_len
        );
        println!(
            "Non-insertion positions: {}",
            parser.get_non_insertion_vars().len()
        );
        for (pos, vars) in parser.get_non_insertion_vars() {
            if vars.is_empty() {
                continue;
            }
            println!("Position: {}", pos);
            for (desc, var) in vars {
                println!("  {:?} -> {:?}", desc, var);
            }
        }

        assert!(!parser.reference.ref_seq.is_empty());
        assert!(!parser.get_non_insertion_vars().is_empty());
        assert!(parser.get_insertion_vars().is_empty());
        assert!(parser.get_soft_clips_5end().is_empty());
        assert!(parser.get_soft_clips_3end().is_empty());
    }

    #[test]
    fn test_cigar_parser_non_insertion_variants_first_bed_region() {
        use crate::conf::Configuration;
        use crate::data::bam_reader::BamReader;
        use crate::data::reference::FastaReader;
        use crate::data::region::Region;
        use crate::mods::vardict_pipeline::VarDictPipeline;
        use crate::scopedata::global_read_only_scope::{GlobalReadOnlyScope, INSTANCE, instance};
        use crate::variants::variants::VarDesc;
        use std::collections::HashMap;
        use std::fs::File;
        use std::io::{BufRead, BufReader, Write};
        use std::sync::Arc;

        let mut conf = Configuration::default();
        conf.disable_sv = true;
        let mut chr_lens = HashMap::new();
        chr_lens.insert("20".to_string(), 63_025_520);

        let _ = INSTANCE.set(GlobalReadOnlyScope {
            conf,
            chr_lens,
            ..Default::default()
        });

        let bam_path = "/home/eck/workspace/vardict_rs/VarDictJava/tests/integration/input/NA12878.chrom20.ILLUMINA.bwa.CEU.exome.20121211.bam";
        let bed_path = "/home/eck/workspace/vardict_rs/VarDictJava/tests/integration/input/20120518.consensus.annotation.bed.chr20";
        let fasta_path = "/home/eck/workspace/vardict_rs/VarDictJava/tests/integration/reference/hs37d5.fa";

        let file = File::open(bed_path).expect("Failed to open BED file");
        let mut reader = BufReader::new(file);
        let mut line = String::new();
        let mut region_opt = None;

        while reader.read_line(&mut line).unwrap_or(0) > 0 {
            let trimmed = line.trim();
            if !trimmed.is_empty()
                && !trimmed.starts_with('#')
                && !trimmed.starts_with("track")
                && !trimmed.starts_with("browser")
            {
                let fields: Vec<&str> = trimmed.split('\t').collect();
                let mut chr = fields.get(0).unwrap().to_string();
                if chr.starts_with("chr") {
                    chr = chr.trim_start_matches("chr").to_string();
                }
                let mut start: usize = fields.get(1).unwrap().parse().unwrap();
                let mut end: usize = fields.get(2).unwrap().parse().unwrap();
                if start < end {
                    start += 1;
                }
                if start == 0 {
                    start = 1;
                }
                if end < start {
                    std::mem::swap(&mut start, &mut end);
                }
                let gene = fields.get(3).unwrap_or(&"").to_string();
                region_opt = Some(Region::new(chr, start, end, gene));
                break;
            }
            line.clear();
        }

        let region = region_opt.expect("No BED regions found");

        let mut ref_start = region.start.saturating_sub(1200);
        if ref_start == 0 {
            ref_start = 1;
        }
        let ref_end = region.end + 1200;

        let fasta = FastaReader::open(fasta_path).expect("Failed to open reference FASTA");
        let reference = fasta
            .get_reference(region.chr(), ref_start, ref_end)
            .expect("Failed to fetch reference sequence");

        let instance = Arc::new(instance().clone());
        let pipeline = VarDictPipeline::new("test");
        let sam_filter = 0x504u32;
        let mut bam_reader = BamReader::open(bam_path).expect("Failed to open BAM");
        let (mut records, _lines) = pipeline
            .collect_filtered_records(&region, &mut bam_reader, sam_filter)
            .expect("Failed to collect filtered reads");

        let mut parser = CigarParser::new(region.clone(), reference.clone(), instance.clone());
        parser
            .process_records(records.iter_mut())
            .expect("process_records failed");

        let non_insertion = parser.get_non_insertion_vars();

        // Collect read names contributing to position 68352 (desc "T") for comparison with Java
        let target_pos = 68352i64;
        let target_desc = VarDesc::snv_key(b'T');
        let mut per_read_parser = CigarParser::new(region.clone(), reference.clone(), instance.clone());
        let mut contributing_reads: Vec<String> = Vec::new();
        let mut extra_debug: Vec<String> = Vec::new();
        let extra_targets: std::collections::HashSet<&str> = [
            "SRR098401.91595700",
            "SRR098401.84181448",
            "SRR098401.3139480",
            "SRR098401.19274949",
            "SRR098401.58730609",
            "SRR098401.60668499",
            "SRR098401.1491532",
        ]
        .into_iter()
        .collect();

        for record in records.iter() {
            let before = per_read_parser
                .get_non_insertion_vars()
                .get(&target_pos)
                .and_then(|m| m.get(&target_desc))
                .map(|v| v.alt_depth)
                .unwrap_or(0);

            let mut single = vec![record.clone()];
            per_read_parser
                .process_records(single.iter_mut())
                .expect("per-read process_records failed");

            let after = per_read_parser
                .get_non_insertion_vars()
                .get(&target_pos)
                .and_then(|m| m.get(&target_desc))
                .map(|v| v.alt_depth)
                .unwrap_or(0);

            let qname = String::from_utf8_lossy(record.qname()).to_string();
            if after > before {
                let qname = String::from_utf8_lossy(record.qname()).to_string();
                contributing_reads.push(qname);
            }

            if extra_targets.contains(qname.as_str()) {
                extra_debug.push(format!(
                    "{}\tflag={}\tpos={}\tmapq={}\tcigar={}\tcontributed={}",
                    qname,
                    record.flags(),
                    record.pos() + 1,
                    record.mapq(),
                    record.cigar().to_string(),
                    after > before
                ));
            }
        }

        let reads_out_path = "/home/eck/workspace/vardict_rs/tmp_compare/rust.cigarparser.non_insertion.pos68352.reads.txt";
        std::fs::create_dir_all("/home/eck/workspace/vardict_rs/tmp_compare")
            .expect("Failed to create tmp_compare");
        let mut reads_writer = std::io::BufWriter::new(
            std::fs::File::create(reads_out_path)
                .expect("Failed to create rust read dump"),
        );
        for name in &contributing_reads {
            writeln!(reads_writer, "{}", name).expect("Failed to write read name");
        }

        let extra_debug_path = "/home/eck/workspace/vardict_rs/tmp_compare/rust.cigarparser.non_insertion.pos68352.extra_reads.detail.txt";
        let mut extra_writer = std::io::BufWriter::new(
            std::fs::File::create(extra_debug_path)
                .expect("Failed to create rust extra read debug"),
        );
        for line in &extra_debug {
            writeln!(extra_writer, "{}", line).expect("Failed to write extra debug line");
        }

        let out_dir = "/home/eck/workspace/vardict_rs/tmp_compare";
        let out_path = "/home/eck/workspace/vardict_rs/tmp_compare/rust.cigarparser.non_insertion.txt";
        std::fs::create_dir_all(out_dir).expect("Failed to create tmp_compare");
        let out_file = std::fs::File::create(out_path).expect("Failed to create rust cigarparser dump");
        let mut writer = std::io::BufWriter::new(out_file);

        for (pos, vars) in non_insertion.iter() {
            if vars.is_empty() {
                continue;
            }
            writeln!(writer, "POSITION\t{}", pos).expect("Failed to write position");
            for (desc, var) in vars.iter() {
                writeln!(
                    writer,
                    "VAR\t{}\t{}\t{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{}",
                    desc.to_key_string(),
                    var.alt_depth,
                    var.alt_depth_fwd,
                    var.alt_depth_rev,
                    var.mean_pos,
                    var.mean_qual,
                    var.mean_mapq,
                    var.nm,
                    var.low_qual_read_cnt,
                    var.high_qual_read_cnt,
                    var.pstd,
                    var.qstd,
                    var.pp,
                    var.pq,
                    0
                )
                .expect("Failed to write variant line");
            }
        }

        assert_variant(
            non_insertion,
            68352,
            b'T',
            127,
            84,
            43,
            2231.000,
            4185.000,
            7558.000,
            10.000,
            3,
            124,
            true,
            true,
            1,
            31.000,
        );

        assert_variant(
            non_insertion,
            68353,
            b'G',
            124,
            81,
            43,
            2261.000,
            4639.000,
            7378.000,
            9.000,
            2,
            122,
            true,
            true,
            1,
            34.000,
        );

        assert_variant(
            non_insertion,
            68359,
            b'C',
            129,
            78,
            51,
            2422.000,
            4505.000,
            7635.000,
            15.000,
            7,
            122,
            true,
            true,
            1,
            34.000,
        );

        assert_variant(
            non_insertion,
            68367,
            b'T',
            135,
            75,
            60,
            2557.000,
            4655.000,
            7964.000,
            17.000,
            6,
            129,
            true,
            true,
            1,
            30.000,
        );
    }

    fn assert_variant(
        non_insertion: &std::collections::HashMap<i64, std::collections::HashMap<crate::variants::variants::VarDesc, crate::variants::variants::Variant>>,
        pos: i64,
        ref_base: u8,
        alt_depth: usize,
        alt_depth_fwd: usize,
        alt_depth_rev: usize,
        mean_pos: f64,
        mean_qual: f64,
        mean_mapq: f64,
        nm: f64,
        low_qual_read_cnt: usize,
        high_qual_read_cnt: usize,
        pstd: bool,
        qstd: bool,
        pp: usize,
        pq: f64,
    ) {
        let vars_at = non_insertion
            .get(&pos)
            .unwrap_or_else(|| panic!("Missing position in non_insertion_vars: {}", pos));
        let key = VarDesc::snv_key(ref_base);
        let var = vars_at
            .get(&key)
            .unwrap_or_else(|| panic!("Missing SNV {:?} at position {}", ref_base as char, pos));

        assert_eq!(var.alt_depth, alt_depth, "alt_depth mismatch at {}", pos);
        assert_eq!(var.alt_depth_fwd, alt_depth_fwd, "alt_depth_fwd mismatch at {}", pos);
        assert_eq!(var.alt_depth_rev, alt_depth_rev, "alt_depth_rev mismatch at {}", pos);
        assert!((var.mean_pos - mean_pos).abs() < 0.001, "mean_pos mismatch at {}", pos);
        assert!((var.mean_qual - mean_qual).abs() < 0.001, "mean_qual mismatch at {}", pos);
        assert!((var.mean_mapq - mean_mapq).abs() < 0.001, "mean_mapq mismatch at {}", pos);
        assert!((var.nm - nm).abs() < 0.001, "nm mismatch at {}", pos);
        assert_eq!(var.low_qual_read_cnt, low_qual_read_cnt, "low_qual_read_cnt mismatch at {}", pos);
        assert_eq!(var.high_qual_read_cnt, high_qual_read_cnt, "high_qual_read_cnt mismatch at {}", pos);
        assert_eq!(var.pstd, pstd, "pstd mismatch at {}", pos);
        assert_eq!(var.qstd, qstd, "qstd mismatch at {}", pos);
        assert_eq!(var.pp, pp, "pp mismatch at {}", pos);
        assert!((var.pq - pq).abs() < 0.001, "pq mismatch at {}", pos);
    }
}