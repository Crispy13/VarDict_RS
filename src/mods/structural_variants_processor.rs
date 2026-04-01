//! StructuralVariantsProcessor - Equivalent to Java StructuralVariantsProcessor
//!
//! This module processes realigned variation data to:
//! 1. Find structural variants (DEL, INV, DUP) when SV is enabled
//! 2. Adjust SNV counts from soft-clipped reads (always runs)
//!
//! Pipeline position:
//! ```text
//! CigarParser → VariantRealigner → StructuralVariantsProcessor → ToVarsBuilder
//! ```

use crackle_kit::{
    data::bases::rev_comp::RevComplementor,
    tracing::{Level, event},
};
use indexmap::IndexMap;
use std::collections::{HashMap, HashSet};
use std::env;
use std::fs::{self, OpenOptions};
use std::io::Write;
use std::sync::Arc;

use crate::conf::Configuration;
use crate::data::reference::{Reference, ReferenceSeedMap};
use crate::data::region::Region;
use crate::data::shared_reference::SharedReferenceHandle;
use crate::mods::variant_realigner::VariantRealigner;
use crate::prelude::{InnerMap, LibDefaultHasher};
use crate::scopedata::global_read_only_scope::instance;
use crate::variants::variants::{SoftClip, StructuralVariantCounts, VarDesc, Variant};

/// Input data for StructuralVariantsProcessor (from VariantRealigner)
#[derive(Default)]
pub struct RealignedVariationData {
    /// Non-insertion variants by position
    pub non_insertion_variants: HashMap<i64, InnerMap<VarDesc, Variant>, LibDefaultHasher>,
    /// Java VariationMap.sv equivalent counts by position.
    pub sv_counts: HashMap<i64, StructuralVariantCounts, LibDefaultHasher>,
    /// Insertion variants by position
    pub insertion_variants: HashMap<i64, InnerMap<VarDesc, Variant>, LibDefaultHasher>,
    /// 5' end soft clips by position
    pub soft_clips_5end: HashMap<i64, SoftClip, LibDefaultHasher>,
    /// 3' end soft clips by position  
    pub soft_clips_3end: HashMap<i64, SoftClip, LibDefaultHasher>,
    /// Reference coverage by position
    pub ref_coverage: HashMap<i64, usize, LibDefaultHasher>,
    /// Splice junction positions carried forward for RNA-seq-specific SV filtering.
    pub splice: HashSet<String>,
    /// Maximum read length seen
    pub max_read_length: usize,
    /// Duplication rate
    pub duprate: f64,
    /// Forward deletion soft-clip clusters
    pub svfdel: Vec<SoftClip>,
    /// Reverse deletion soft-clip clusters
    pub svrdel: Vec<SoftClip>,
    /// Forward duplication discordant clusters
    pub svfdup: Vec<SoftClip>,
    /// Reverse duplication discordant clusters
    pub svrdup: Vec<SoftClip>,
    /// Forward inversion 5' discordant clusters
    pub svfinv5: Vec<SoftClip>,
    /// Reverse inversion 5' discordant clusters
    pub svrinv5: Vec<SoftClip>,
    /// Forward inversion 3' discordant clusters
    pub svfinv3: Vec<SoftClip>,
    /// Reverse inversion 3' discordant clusters
    pub svrinv3: Vec<SoftClip>,
    /// Forward inter-chromosomal fusion clusters needed for Java SOFTP2SV parity.
    pub svffus: HashMap<i32, Vec<SoftClip>, LibDefaultHasher>,
    /// Reverse inter-chromosomal fusion clusters needed for Java SOFTP2SV parity.
    pub svrfus: HashMap<i32, Vec<SoftClip>, LibDefaultHasher>,
    /// Java-compatible SOFTP2SV top-entry used-state by soft clip position.
    pub softp2sv_first_used: HashMap<i64, bool, LibDefaultHasher>,
}

/// Output data from StructuralVariantsProcessor (same structure, possibly modified)
pub type ProcessedVariationData = RealignedVariationData;

#[derive(Clone, Copy)]
struct ProcessMemorySnapshot {
    vmrss_kb: usize,
    vmswap_kb: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct ReferenceWindow {
    start: i64,
    end: i64,
}

fn current_process_memory_snapshot() -> Option<ProcessMemorySnapshot> {
    let status = fs::read_to_string("/proc/self/status").ok()?;
    let mut vmrss_kb = None;
    let mut vmswap_kb = None;

    for line in status.lines() {
        if let Some(value) = line.strip_prefix("VmRSS:") {
            vmrss_kb = value.split_whitespace().next()?.parse::<usize>().ok();
        } else if let Some(value) = line.strip_prefix("VmSwap:") {
            vmswap_kb = value.split_whitespace().next()?.parse::<usize>().ok();
        }

        if vmrss_kb.is_some() && vmswap_kb.is_some() {
            break;
        }
    }

    Some(ProcessMemorySnapshot {
        vmrss_kb: vmrss_kb?,
        vmswap_kb: vmswap_kb.unwrap_or(0),
    })
}

fn append_rss_stage_log_if_enabled(region: Option<&Region>, stage: &str) {
    let Some(region) = region else {
        return;
    };

    let Some(path) = env::var_os("VARDICT_RSS_MEMORY_LOG") else {
        return;
    };
    if path.is_empty() {
        return;
    }

    let Some(snapshot) = current_process_memory_snapshot() else {
        return;
    };

    let path = std::path::PathBuf::from(path);
    let should_write_header = fs::metadata(&path)
        .map(|meta| meta.len() == 0)
        .unwrap_or(true);
    let mut file = match OpenOptions::new().create(true).append(true).open(&path) {
        Ok(file) => file,
        Err(_) => return,
    };

    if should_write_header {
        let _ = writeln!(
            file,
            "elapsed_ms\tregion_chr\tregion_start\tregion_end\tvmrss_kb\tvmswap_kb\tstage"
        );
    }

    let _ = writeln!(
        file,
        "0\t{}\t{}\t{}\t{}\t{}\t{}",
        region.chr(),
        region.start(),
        region.end(),
        snapshot.vmrss_kb,
        snapshot.vmswap_kb,
        stage,
    );
}

/// StructuralVariantsProcessor - Processes variation data for SVs and adjusts SNV counts
pub struct StructuralVariantsProcessor {
    /// Reference sequence
    reference_seq: Arc<Vec<u8>>,
    /// Reference seed map
    reference_seed: Arc<ReferenceSeedMap>,
    /// Reference start position (1-based)
    ref_start: i64,
    /// Original region reference sequence to return to downstream stages.
    original_reference_seq: Arc<Vec<u8>>,
    /// Original region reference seed map to return to downstream stages.
    original_reference_seed: Arc<ReferenceSeedMap>,
    /// Original region reference start position to return to downstream stages.
    original_ref_start: i64,
    /// Chromosome name for on-demand reference extension
    chromosome: Option<String>,
    /// BAM paths used by realignment parity paths
    bam_paths: Vec<String>,
    /// Shared reference handle for on-demand reference extension
    shared_reference: Option<SharedReferenceHandle>,
    /// Historical remote windows loaded earlier in SV processing.
    historical_reference_windows: Vec<ReferenceWindow>,
    /// DEL variants whose trailing rightseq span was inside the active remote SV window.
    historical_del_rightseq_variants:
        HashMap<i64, HashSet<String, LibDefaultHasher>, LibDefaultHasher>,
    /// Tracks reference regions explicitly loaded by SV methods, mirroring Java's
    /// Reference.loadedRegions list.
    loaded_regions: Vec<(i64, i64)>,
    /// Remote windows whose BAM-derived coverage has already been merged into the live data.
    loaded_reference_coverage_windows: Vec<ReferenceWindow>,
}

impl StructuralVariantsProcessor {
    /// Create a new StructuralVariantsProcessor
    pub fn new(reference_seq: Vec<u8>, reference_seed: ReferenceSeedMap, ref_start: i64) -> Self {
        let reference_seq = Arc::new(reference_seq);
        let reference_seed = Arc::new(reference_seed);
        StructuralVariantsProcessor {
            reference_seq: Arc::clone(&reference_seq),
            reference_seed: Arc::clone(&reference_seed),
            ref_start,
            original_reference_seq: reference_seq,
            original_reference_seed: reference_seed,
            original_ref_start: ref_start,
            chromosome: None,
            bam_paths: Vec::new(),
            shared_reference: None,
            historical_reference_windows: Vec::new(),
            historical_del_rightseq_variants: Default::default(),
            loaded_regions: Vec::new(),
            loaded_reference_coverage_windows: Vec::new(),
        }
    }

    /// Create a StructuralVariantsProcessor with optional on-demand reference extension context.
    pub fn new_with_context(
        reference_seq: Arc<Vec<u8>>,
        reference_seed: Arc<ReferenceSeedMap>,
        ref_start: i64,
        chromosome: Option<String>,
        bam_paths: Vec<String>,
        shared_reference: Option<SharedReferenceHandle>,
    ) -> Self {
        StructuralVariantsProcessor {
            reference_seq: Arc::clone(&reference_seq),
            reference_seed: Arc::clone(&reference_seed),
            ref_start,
            original_reference_seq: reference_seq,
            original_reference_seed: reference_seed,
            original_ref_start: ref_start,
            chromosome,
            bam_paths,
            shared_reference,
            historical_reference_windows: Vec::new(),
            historical_del_rightseq_variants: Default::default(),
            loaded_regions: Vec::new(),
            loaded_reference_coverage_windows: Vec::new(),
        }
    }

    /// Process realigned variation data
    ///
    /// This is equivalent to Java StructuralVariantsProcessor.process()
    pub fn process(&mut self, mut data: RealignedVariationData) -> ProcessedVariationData {
        self.process_internal(&mut data, None);
        data
    }

    pub fn process_with_region(
        &mut self,
        mut data: RealignedVariationData,
        region: &Region,
    ) -> ProcessedVariationData {
        self.process_internal(&mut data, Some(region));
        data
    }

    fn process_internal(&mut self, data: &mut RealignedVariationData, region: Option<&Region>) {
        // If SV is enabled, find structural variants
        if !instance().conf.disable_sv {
            self.find_all_svs(data, region);
        }

        // adj_snv() must validate nearby reference bases across the whole original region.
        // The SV discovery passes can temporarily narrow the working reference window for
        // targeted lookups, so restore the original region reference before soft-clip rescue.
        self.restore_original_reference_window();

        // Always adjust SNV counts from soft clips
        self.adj_snv(data);
        append_rss_stage_log_if_enabled(region, "sv_adj_snv_complete");
    }

    pub fn into_reference(self) -> Reference {
        Reference {
            ref_seq: self.original_reference_seq,
            seed: self.original_reference_seed,
            region_start: self.original_ref_start,
        }
    }

    pub fn historical_reference_windows(&self) -> Vec<(i64, i64)> {
        self.historical_reference_windows
            .iter()
            .map(|window| (window.start, window.end))
            .collect()
    }

    pub fn historical_del_rightseq_variants(
        &self,
    ) -> HashMap<i64, HashSet<String, LibDefaultHasher>, LibDefaultHasher> {
        self.historical_del_rightseq_variants.clone()
    }

    /// Find all structural variants (DEL, INV, DUP)
    ///
    /// Called when SV detection is enabled
    fn find_all_svs(&mut self, data: &mut RealignedVariationData, region: Option<&Region>) {
        let nnte = data.max_read_length as i64;
        let orig_end = self.original_ref_start + self.original_reference_seq.len() as i64 - 1;
        self.loaded_regions
            .push((self.original_ref_start - 2 * nnte, orig_end + 2 * nnte));

        self.find_del(data);
        append_rss_stage_log_if_enabled(region, "sv_find_del_complete");
        self.find_inv(data, region);
        append_rss_stage_log_if_enabled(region, "sv_find_inv_complete");
        self.restore_original_reference_window_preserving_history();
        self.find_svs_del_candidates(data, region);
        append_rss_stage_log_if_enabled(region, "sv_find_svs_del_candidates_complete");
        self.find_del_disc(data);
        append_rss_stage_log_if_enabled(region, "sv_find_del_disc_complete");
        self.find_inv_disc(data, region);
        append_rss_stage_log_if_enabled(region, "sv_find_inv_disc_complete");
        self.find_dup_disc(data);
        self.record_del_rightseq_from_loaded_regions(data);
        append_rss_stage_log_if_enabled(region, "sv_find_dup_disc_complete");
    }

    fn find_del(&mut self, data: &mut RealignedVariationData) {
        let minr = instance().conf.minr;

        for idx in 0..data.svfdel.len() {
            let (used, vars_count, end, mstart, mend, mean_qual, mean_pos, mean_mapq, nm, soft_map) = {
                let del = &data.svfdel[idx];
                (
                    del.used(),
                    del.var.alt_depth,
                    del.end,
                    del.mstart,
                    del.mend,
                    del.var.mean_qual,
                    del.var.mean_pos,
                    del.var.mean_mapq,
                    del.var.nm,
                    del.soft.clone(),
                )
            };

            if used || vars_count < minr {
                continue;
            }

            let region_was_loaded = self.is_region_loaded(mstart, mend);
            let span_preloaded = self.is_span_loaded(mstart, mend);
            self.ensure_reference_span(mstart - 300, mend + 300);
            if !region_was_loaded {
                let nnte = data.max_read_length as i64;
                self.loaded_regions
                    .push((mstart - nnte - 300, mend + nnte + 300));
            }

            let softp = Self::select_primary_soft_pos(&soft_map).unwrap_or(0);
            event!(
                Level::DEBUG,
                phase = "find_del_forward_cluster",
                idx,
                softp,
                vars_count,
                end,
                mstart,
            );
            if softp != 0 {
                let seq = {
                    let Some(scv) = data.soft_clips_3end.get_mut(&softp) else {
                        continue;
                    };
                    if scv.used() {
                        continue;
                    }
                    self.find_conseq(scv)
                };

                if seq.is_empty() || seq.len() < Configuration::SEED_2 as usize {
                    continue;
                }

                if !span_preloaded {
                    self.load_uncovered_reference_coverage(data, mstart - 200, mend + 200);
                }

                let m = self.find_match(&seq, softp, 1, Configuration::SEED_1 as usize, 3);
                if m.base_position == 0 {
                    event!(Level::DEBUG, phase = "find_del_forward_nomatch", idx, softp,);
                    continue;
                }

                if !(m.base_position - softp > 30
                    && Self::is_overlap(
                        softp,
                        m.base_position,
                        end,
                        mstart,
                        data.max_read_length as i64,
                    ))
                {
                    event!(
                        Level::DEBUG,
                        phase = "find_del_forward_overlap_fail",
                        idx,
                        softp,
                        bp = m.base_position,
                    );
                    continue;
                }

                let mut bp = m.base_position - 1;
                let dellen = bp - softp + 1;
                if dellen <= 0 {
                    continue;
                }
                let mut p5 = softp;
                while self.get_ref_base(bp).is_some()
                    && self.get_ref_base(p5 - 1).is_some()
                    && self.get_ref_base(bp) == self.get_ref_base(p5 - 1)
                {
                    bp -= 1;
                    if bp != 0 {
                        p5 -= 1;
                    }
                }

                let del_key = format!("-{}", dellen);
                let vref =
                    Self::get_or_create_variation(&mut data.non_insertion_variants, p5, &del_key);
                vref.alt_depth = 0;

                let split_count = data
                    .soft_clips_3end
                    .get(&softp)
                    .map(|s| s.var.alt_depth)
                    .unwrap_or(0);
                Self::add_sv_counts(
                    &mut data.non_insertion_variants,
                    &mut data.sv_counts,
                    p5,
                    vars_count,
                    split_count,
                    1,
                );

                let current_cov = data.ref_coverage.get(&p5).copied().unwrap_or(0);
                if current_cov <= vars_count {
                    data.ref_coverage.insert(p5, vars_count);
                }
                if let Some(bp_cov) = data.ref_coverage.get(&bp).copied() {
                    let p5_cov = data.ref_coverage.get(&p5).copied().unwrap_or(0);
                    if p5_cov < bp_cov {
                        data.ref_coverage.insert(p5, bp_cov);
                    }
                }

                if let Some(scv) = data.soft_clips_3end.get(&softp) {
                    let scv_var = scv.var.clone();
                    let variation = Self::get_or_create_variation(
                        &mut data.non_insertion_variants,
                        p5,
                        &del_key,
                    );
                    adj_cnt_from_variant(variation, &scv_var);

                    if let Some(ref_base) = self.get_ref_base(p5) {
                        let ref_key = (ref_base as char).to_string();
                        if let Some(pos_map) = data.non_insertion_variants.get_mut(&p5) {
                            if let Some(reference_var) =
                                Self::get_variation_mut_by_key_string(pos_map, &ref_key)
                            {
                                sub_cnt_from_variant(reference_var, &scv_var);
                            }
                        }
                    }
                }

                let mut tv = Variant::default();
                tv.alt_depth = vars_count;
                tv.high_qual_read_cnt = vars_count;
                tv.alt_depth_fwd = vars_count / 2;
                tv.alt_depth_rev = vars_count - vars_count / 2;
                tv.mean_qual = mean_qual * vars_count as f64 / vars_count as f64;
                tv.mean_pos = mean_pos * vars_count as f64 / vars_count as f64;
                tv.mean_mapq = mean_mapq * vars_count as f64 / vars_count as f64;
                tv.nm = nm * vars_count as f64 / vars_count as f64;

                let variation =
                    Self::get_or_create_variation(&mut data.non_insertion_variants, p5, &del_key);
                adj_cnt_from_variant(variation, &tv);

                if let Some(del) = data.svfdel.get_mut(idx) {
                    del.mark_used();
                }
                event!(
                    Level::DEBUG,
                    phase = "find_del_forward_emit",
                    idx,
                    p5,
                    bp,
                    dellen,
                    vars_count,
                    split_count,
                );
                Self::mark_sv(p5, bp, &mut data.svrdel, data.max_read_length as i64);
            } else {
                let mut candidate_positions: Vec<i64> =
                    data.soft_clips_3end.keys().copied().collect();
                candidate_positions.sort_unstable();

                for candidate in candidate_positions {
                    if !(candidate >= end - 3 && candidate - end < 3 * data.max_read_length as i64)
                    {
                        continue;
                    }

                    let seq = {
                        let Some(scv) = data.soft_clips_3end.get_mut(&candidate) else {
                            continue;
                        };
                        if scv.used() {
                            continue;
                        }
                        self.find_conseq(scv)
                    };
                    if seq.is_empty() || seq.len() < Configuration::SEED_2 as usize {
                        continue;
                    }

                    let mut m =
                        self.find_match(&seq, candidate, 1, Configuration::SEED_1 as usize, 3);
                    if m.base_position == 0 {
                        m = self.find_match(&seq, candidate, 1, Configuration::SEED_2 as usize, 0);
                    }
                    if m.base_position == 0 {
                        event!(
                            Level::DEBUG,
                            phase = "find_del_forward_nosoftp_nomatch",
                            idx,
                            candidate,
                        );
                        continue;
                    }

                    if !(m.base_position - candidate > 30
                        && Self::is_overlap(
                            candidate,
                            m.base_position,
                            end,
                            mstart,
                            data.max_read_length as i64,
                        ))
                    {
                        event!(
                            Level::DEBUG,
                            phase = "find_del_forward_nosoftp_overlap_fail",
                            idx,
                            candidate,
                            bp = m.base_position,
                        );
                        continue;
                    }

                    let bp = m.base_position - 1;
                    let dellen = bp - candidate + 1;
                    if dellen <= 0 {
                        continue;
                    }

                    let del_key = format!("-{}", dellen);
                    let vref = Self::get_or_create_variation(
                        &mut data.non_insertion_variants,
                        candidate,
                        &del_key,
                    );
                    vref.alt_depth = 0;

                    let split_count = data
                        .soft_clips_3end
                        .get(&candidate)
                        .map(|s| s.var.alt_depth)
                        .unwrap_or(0);
                    Self::add_sv_counts(
                        &mut data.non_insertion_variants,
                        &mut data.sv_counts,
                        candidate,
                        vars_count,
                        split_count,
                        1,
                    );

                    let current_cov = data.ref_coverage.get(&candidate).copied().unwrap_or(0);
                    if current_cov <= vars_count {
                        data.ref_coverage.insert(candidate, vars_count);
                    }
                    if let Some(bp_cov) = data.ref_coverage.get(&bp).copied() {
                        let p5_cov = data.ref_coverage.get(&candidate).copied().unwrap_or(0);
                        if p5_cov < bp_cov {
                            data.ref_coverage.insert(candidate, bp_cov);
                        }
                    }

                    if let Some(scv) = data.soft_clips_3end.get(&candidate) {
                        let variation = Self::get_or_create_variation(
                            &mut data.non_insertion_variants,
                            candidate,
                            &del_key,
                        );
                        adj_cnt_from_variant(variation, &scv.var);
                    }

                    let mut tv = Variant::default();
                    tv.alt_depth = vars_count;
                    tv.high_qual_read_cnt = vars_count;
                    tv.alt_depth_fwd = vars_count / 2;
                    tv.alt_depth_rev = vars_count - vars_count / 2;
                    tv.mean_qual = mean_qual * vars_count as f64 / vars_count as f64;
                    tv.mean_pos = mean_pos * vars_count as f64 / vars_count as f64;
                    tv.mean_mapq = mean_mapq * vars_count as f64 / vars_count as f64;
                    tv.nm = nm * vars_count as f64 / vars_count as f64;

                    let variation = Self::get_or_create_variation(
                        &mut data.non_insertion_variants,
                        candidate,
                        &del_key,
                    );
                    adj_cnt_from_variant(variation, &tv);

                    if let Some(del) = data.svfdel.get_mut(idx) {
                        del.mark_used();
                    }
                    event!(
                        Level::DEBUG,
                        phase = "find_del_forward_nosoftp_emit",
                        idx,
                        candidate,
                        bp,
                        dellen,
                        vars_count,
                        split_count,
                    );
                    Self::mark_sv(candidate, bp, &mut data.svrdel, data.max_read_length as i64);
                    break;
                }
            }
        }

        for idx in 0..data.svrdel.len() {
            let (
                used,
                vars_count,
                mend,
                start,
                mstart,
                mean_qual,
                mean_pos,
                mean_mapq,
                nm,
                soft_map,
            ) = {
                let del = &data.svrdel[idx];
                (
                    del.used(),
                    del.var.alt_depth,
                    del.mend,
                    del.start,
                    del.mstart,
                    del.var.mean_qual,
                    del.var.mean_pos,
                    del.var.mean_mapq,
                    del.var.nm,
                    del.soft.clone(),
                )
            };

            if used || vars_count < minr {
                continue;
            }

            let region_was_loaded = self.is_region_loaded(mstart, mend);
            let span_preloaded = self.is_span_loaded(mstart, mend);
            self.ensure_reference_span(mstart - 300, mend + 300);
            if !region_was_loaded {
                self.loaded_regions.push((mstart - 300, mend + 300));
            }

            let softp = Self::select_primary_soft_pos(&soft_map).unwrap_or(0);
            event!(
                Level::DEBUG,
                phase = "find_del_reverse_cluster",
                idx,
                softp,
                vars_count,
                start,
                mend,
            );
            if softp != 0 {
                let seq = {
                    let Some(scv) = data.soft_clips_5end.get_mut(&softp) else {
                        continue;
                    };
                    if scv.used() {
                        continue;
                    }
                    self.find_conseq(scv)
                };

                if seq.is_empty() || seq.len() < Configuration::SEED_2 as usize {
                    continue;
                }

                if !span_preloaded {
                    self.load_uncovered_reference_coverage(data, mstart - 200, mend + 200);
                }

                let mut m = self.find_match(&seq, softp, -1, Configuration::SEED_1 as usize, 3);
                if m.base_position == 0 {
                    m = self.find_match(&seq, softp, -1, Configuration::SEED_2 as usize, 0);
                }
                if m.base_position == 0 {
                    event!(Level::DEBUG, phase = "find_del_reverse_nomatch", idx, softp,);
                    continue;
                }

                if !(softp - m.base_position > 30
                    && Self::is_overlap(
                        m.base_position,
                        softp,
                        mend,
                        start,
                        data.max_read_length as i64,
                    ))
                {
                    event!(
                        Level::DEBUG,
                        phase = "find_del_reverse_overlap_fail",
                        idx,
                        softp,
                        bp = m.base_position,
                    );
                    continue;
                }

                let bp = m.base_position + 1;
                let p3 = softp - 1;
                let dellen = p3 - bp + 1;
                if dellen <= 0 {
                    continue;
                }

                let del_key = format!("-{}", dellen);
                let vref =
                    Self::get_or_create_variation(&mut data.non_insertion_variants, bp, &del_key);
                vref.alt_depth = 0;

                let split_count = data
                    .soft_clips_5end
                    .get(&softp)
                    .map(|s| s.var.alt_depth)
                    .unwrap_or(0);
                Self::add_sv_counts(
                    &mut data.non_insertion_variants,
                    &mut data.sv_counts,
                    bp,
                    vars_count,
                    split_count,
                    1,
                );

                if let Some(scv) = data.soft_clips_5end.get(&softp) {
                    let variation = Self::get_or_create_variation(
                        &mut data.non_insertion_variants,
                        bp,
                        &del_key,
                    );
                    adj_cnt_from_variant(variation, &scv.var);
                }

                let current_cov = data.ref_coverage.get(&bp).copied().unwrap_or(0);
                if current_cov <= vars_count {
                    data.ref_coverage.insert(bp, vars_count);
                }
                if let Some(p3_cov) = data.ref_coverage.get(&p3).copied() {
                    let bp_cov = data.ref_coverage.get(&bp).copied().unwrap_or(0);
                    if p3_cov > bp_cov {
                        data.ref_coverage.insert(bp, p3_cov);
                    }
                }

                let mut tv = Variant::default();
                tv.alt_depth = vars_count;
                tv.high_qual_read_cnt = vars_count;
                tv.alt_depth_fwd = vars_count / 2;
                tv.alt_depth_rev = vars_count - vars_count / 2;
                tv.mean_qual = mean_qual * vars_count as f64 / vars_count as f64;
                tv.mean_pos = mean_pos * vars_count as f64 / vars_count as f64;
                tv.mean_mapq = mean_mapq * vars_count as f64 / vars_count as f64;
                tv.nm = nm * vars_count as f64 / vars_count as f64;

                let variation =
                    Self::get_or_create_variation(&mut data.non_insertion_variants, bp, &del_key);
                adj_cnt_from_variant(variation, &tv);

                if let Some(del) = data.svrdel.get_mut(idx) {
                    del.mark_used();
                }
                event!(
                    Level::DEBUG,
                    phase = "find_del_reverse_emit",
                    idx,
                    bp,
                    p3,
                    dellen,
                    vars_count,
                    split_count,
                );
                Self::mark_sv(bp, p3, &mut data.svfdel, data.max_read_length as i64);
            } else {
                let mut candidate_positions: Vec<i64> =
                    data.soft_clips_5end.keys().copied().collect();
                candidate_positions.sort_unstable();

                for candidate in candidate_positions {
                    if !(candidate <= start + 3
                        && start - candidate < 3 * data.max_read_length as i64)
                    {
                        continue;
                    }

                    let seq = {
                        let Some(scv) = data.soft_clips_5end.get_mut(&candidate) else {
                            continue;
                        };
                        if scv.used() {
                            continue;
                        }
                        self.find_conseq(scv)
                    };

                    if seq.is_empty() || seq.len() < Configuration::SEED_2 as usize {
                        continue;
                    }

                    let mut m =
                        self.find_match(&seq, candidate, -1, Configuration::SEED_1 as usize, 3);
                    if m.base_position == 0 {
                        m = self.find_match(&seq, candidate, -1, Configuration::SEED_2 as usize, 0);
                    }
                    if m.base_position == 0 {
                        event!(
                            Level::DEBUG,
                            phase = "find_del_reverse_nosoftp_nomatch",
                            idx,
                            candidate,
                        );
                        continue;
                    }

                    if !(candidate - m.base_position > 30
                        && Self::is_overlap(
                            m.base_position,
                            candidate,
                            mend,
                            start,
                            data.max_read_length as i64,
                        ))
                    {
                        event!(
                            Level::DEBUG,
                            phase = "find_del_reverse_nosoftp_overlap_fail",
                            idx,
                            candidate,
                            bp = m.base_position,
                        );
                        continue;
                    }

                    let bp = m.base_position + 1;
                    let p3 = candidate - 1;
                    let dellen = p3 - bp + 1;
                    if dellen <= 0 {
                        continue;
                    }

                    let del_key = format!("-{}", dellen);
                    let vref = Self::get_or_create_variation(
                        &mut data.non_insertion_variants,
                        bp,
                        &del_key,
                    );
                    vref.alt_depth = 0;

                    let split_count = data
                        .soft_clips_5end
                        .get(&candidate)
                        .map(|s| s.var.alt_depth)
                        .unwrap_or(0);
                    Self::add_sv_counts(
                        &mut data.non_insertion_variants,
                        &mut data.sv_counts,
                        bp,
                        vars_count,
                        split_count,
                        1,
                    );

                    if let Some(scv) = data.soft_clips_5end.get(&candidate) {
                        let variation = Self::get_or_create_variation(
                            &mut data.non_insertion_variants,
                            bp,
                            &del_key,
                        );
                        adj_cnt_from_variant(variation, &scv.var);
                    }

                    if !data.ref_coverage.contains_key(&bp) {
                        data.ref_coverage.insert(bp, vars_count);
                    }
                    if let Some(p3_cov) = data.ref_coverage.get(&p3).copied() {
                        let bp_cov = data.ref_coverage.get(&bp).copied().unwrap_or(0);
                        if p3_cov > bp_cov {
                            data.ref_coverage.insert(bp, p3_cov);
                        }
                    }
                    Self::inc_ref_coverage(&mut data.ref_coverage, bp, split_count);

                    let mut tv = Variant::default();
                    tv.alt_depth = vars_count;
                    tv.high_qual_read_cnt = vars_count;
                    tv.alt_depth_fwd = vars_count / 2;
                    tv.alt_depth_rev = vars_count - vars_count / 2;
                    tv.mean_qual = mean_qual * vars_count as f64 / vars_count as f64;
                    tv.mean_pos = mean_pos * vars_count as f64 / vars_count as f64;
                    tv.mean_mapq = mean_mapq * vars_count as f64 / vars_count as f64;
                    tv.nm = nm * vars_count as f64 / vars_count as f64;

                    let variation = Self::get_or_create_variation(
                        &mut data.non_insertion_variants,
                        bp,
                        &del_key,
                    );
                    adj_cnt_from_variant(variation, &tv);

                    if let Some(del) = data.svrdel.get_mut(idx) {
                        del.mark_used();
                    }
                    event!(
                        Level::DEBUG,
                        phase = "find_del_reverse_nosoftp_emit",
                        idx,
                        candidate,
                        bp,
                        p3,
                        dellen,
                        vars_count,
                        split_count,
                    );
                    Self::mark_sv(bp, p3, &mut data.svfdel, data.max_read_length as i64);
                    break;
                }
            }
        }
    }

    fn find_inv(&mut self, data: &mut RealignedVariationData, region: Option<&Region>) {
        self.find_inv_sub(
            data,
            InversionClusterKind::Forward5,
            1,
            InversionSide::End5,
            region,
        );
        self.find_inv_sub(
            data,
            InversionClusterKind::Reverse5,
            -1,
            InversionSide::End5,
            region,
        );
        self.find_inv_sub(
            data,
            InversionClusterKind::Forward3,
            1,
            InversionSide::End3,
            region,
        );
        self.find_inv_sub(
            data,
            InversionClusterKind::Reverse3,
            -1,
            InversionSide::End3,
            region,
        );
    }

    fn find_inv_sub(
        &mut self,
        data: &mut RealignedVariationData,
        kind: InversionClusterKind,
        dir: i64,
        side: InversionSide,
        region: Option<&Region>,
    ) {
        let min_cluster_dist = (Configuration::MINSVCDIST * data.max_read_length as f64) as i64;
        let diag = std::env::var("VARDICT_DIAG_INV_SUB").is_ok();

        for idx in 0..Self::inv_cluster_len(data, kind) {
            let Some(inv) = Self::inv_cluster_snapshot(data, kind, idx) else {
                continue;
            };

            if diag {
                eprintln!(
                    "[DIAG find_inv_sub] idx={} kind={:?} dir={} side={:?} used={} vars_count={} start={} end={} mstart={} mend={} mlen={} primary_softp={:?}",
                    idx, kind, dir, side, inv.used, inv.vars_count, inv.start, inv.end, inv.mstart, inv.mend, inv.mlen, inv.primary_softp
                );
            }

            if inv.used || inv.vars_count < instance().conf.minr {
                if diag {
                    eprintln!("[DIAG find_inv_sub]   SKIP: used={} vars_count={} minr={}", inv.used, inv.vars_count, instance().conf.minr);
                }
                continue;
            }

            let region_was_loaded = self.is_region_loaded(inv.mstart, inv.mend);
            let inv_span_preloaded = self.is_span_loaded(inv.mstart, inv.mend);
            append_rss_stage_log_if_enabled(region, "sv_find_inv_before_ensure_reference_span");
            self.ensure_reference_span(inv.mstart - 500, inv.mend + 500);
            append_rss_stage_log_if_enabled(region, "sv_find_inv_after_ensure_reference_span");
            if !region_was_loaded {
                self.loaded_regions.push((inv.mstart - 500, inv.mend + 500));
            }

            if !inv_span_preloaded {
                append_rss_stage_log_if_enabled(
                    region,
                    "sv_find_inv_before_load_partial_ref_coverage",
                );
                self.load_uncovered_reference_coverage(data, inv.mstart - 200, inv.mend + 200);
                append_rss_stage_log_if_enabled(
                    region,
                    "sv_find_inv_after_load_partial_ref_coverage",
                );
            }

            let mut softp = inv.primary_softp.unwrap_or(0);
            let mut bp = 0i64;
            let mut extra = Vec::new();
            let mut scv_var: Option<Variant> = None;
            let mut source_softp = 0i64;

            if diag {
                eprintln!("[DIAG find_inv_sub]   softp={}", softp);
            }

            if softp != 0 {
                if dir == 1 {
                    let Some(scv) = data.soft_clips_3end.get_mut(&softp) else {
                        if diag { eprintln!("[DIAG find_inv_sub]   softp={} NOT in soft_clips_3end -> continue", softp); }
                        continue;
                    };
                    if scv.used() {
                        if diag { eprintln!("[DIAG find_inv_sub]   softp={} scv.used=true -> continue", softp); }
                        continue;
                    }
                    source_softp = softp;
                    let scv_seq = self.ensure_conseq(scv);
                    if diag { eprintln!("[DIAG find_inv_sub]   softp={} conseq_len={} seq={}", softp, scv_seq.len(), String::from_utf8_lossy(&scv_seq[..scv_seq.len().min(30)])); }
                    if scv_seq.is_empty() {
                        if diag { eprintln!("[DIAG find_inv_sub]   softp={} empty conseq -> continue", softp); }
                        continue;
                    }
                    let mut m =
                        self.find_match_rev(scv_seq, softp, dir, Configuration::SEED_1 as usize, 3);
                    if diag { eprintln!("[DIAG find_inv_sub]   find_match_rev(SEED1) bp={}", m.base_position); }
                    if m.base_position == 0 {
                        m = self.find_match_rev(
                            scv_seq,
                            softp,
                            dir,
                            Configuration::SEED_2 as usize,
                            0,
                        );
                        if diag { eprintln!("[DIAG find_inv_sub]   find_match_rev(SEED2) bp={}", m.base_position); }
                    }
                    if m.base_position == 0 {
                        if diag { eprintln!("[DIAG find_inv_sub]   bp=0 after both seeds -> continue"); }
                        continue;
                    }
                    bp = m.base_position;
                    extra = m.matched_sequence;
                    scv_var = Some(scv.var.clone());
                } else {
                    let Some(scv) = data.soft_clips_5end.get_mut(&softp) else {
                        if diag { eprintln!("[DIAG find_inv_sub]   softp={} NOT in soft_clips_5end -> continue", softp); }
                        continue;
                    };
                    if scv.used() {
                        if diag { eprintln!("[DIAG find_inv_sub]   softp={} scv.used=true (5end) -> continue", softp); }
                        continue;
                    }
                    source_softp = softp;
                    let scv_seq = self.ensure_conseq(scv);
                    if diag { eprintln!("[DIAG find_inv_sub]   softp={} 5end conseq_len={}", softp, scv_seq.len()); }
                    if scv_seq.is_empty() {
                        if diag { eprintln!("[DIAG find_inv_sub]   softp={} empty conseq (5end) -> continue", softp); }
                        continue;
                    }
                    let mut m =
                        self.find_match_rev(scv_seq, softp, dir, Configuration::SEED_1 as usize, 3);
                    if diag { eprintln!("[DIAG find_inv_sub]   5end find_match_rev(SEED1) bp={}", m.base_position); }
                    if m.base_position == 0 {
                        m = self.find_match_rev(
                            scv_seq,
                            softp,
                            dir,
                            Configuration::SEED_2 as usize,
                            0,
                        );
                        if diag { eprintln!("[DIAG find_inv_sub]   5end find_match_rev(SEED2) bp={}", m.base_position); }
                    }
                    if m.base_position == 0 {
                        if diag { eprintln!("[DIAG find_inv_sub]   5end bp=0 -> continue"); }
                        continue;
                    }
                    bp = m.base_position;
                    extra = m.matched_sequence;
                    scv_var = Some(scv.var.clone());
                }
            } else {
                let sp = if dir == 1 { inv.end } else { inv.start };
                for i in 1..=2 * data.max_read_length {
                    let cp = sp + i as i64 * dir;

                    if dir == 1 {
                        let Some(scv) = data.soft_clips_3end.get_mut(&cp) else {
                            continue;
                        };
                        if scv.used() {
                            if diag { eprintln!("[DIAG find_inv_sub]   search cp={} used=true(3end) -> skip", cp); }
                            continue;
                        }
                        source_softp = cp;
                        let scv_seq = self.ensure_conseq(scv);
                        if diag { eprintln!("[DIAG find_inv_sub]   search cp={} conseq_len={} seq={}", cp, scv_seq.len(), String::from_utf8_lossy(&scv_seq[..scv_seq.len().min(30)])); }
                        if scv_seq.is_empty() {
                            continue;
                        }
                        let mut m = self.find_match_rev(
                            scv_seq,
                            cp,
                            dir,
                            Configuration::SEED_1 as usize,
                            3,
                        );
                        if diag { eprintln!("[DIAG find_inv_sub]   search cp={} find_match_rev(SEED1) bp={}", cp, m.base_position); }
                        // NOTE: Preserving Java's behavior where bp/extra are always
                        // overwritten before the bp==0 check. A later failed match must
                        // reset bp to 0 so the post-loop `if bp == 0` guard fires.
                        // Java: StructuralVariantsProcessor.java:L779-L786
                        bp = m.base_position;
                        extra = m.matched_sequence;
                        if bp == 0 {
                            m = self.find_match_rev(
                                scv_seq,
                                cp,
                                dir,
                                Configuration::SEED_2 as usize,
                                0,
                            );
                            if diag { eprintln!("[DIAG find_inv_sub]   search cp={} find_match_rev(SEED2) bp={}", cp, m.base_position); }
                            bp = m.base_position;
                            extra = m.matched_sequence;
                        }
                        if bp == 0 {
                            continue;
                        }
                        scv_var = Some(scv.var.clone());
                    } else {
                        let Some(scv) = data.soft_clips_5end.get_mut(&cp) else {
                            continue;
                        };
                        if scv.used() {
                            continue;
                        }
                        source_softp = cp;
                        let scv_seq = self.ensure_conseq(scv);
                        if scv_seq.is_empty() {
                            continue;
                        }
                        let mut m = self.find_match_rev(
                            scv_seq,
                            cp,
                            dir,
                            Configuration::SEED_1 as usize,
                            3,
                        );
                        // NOTE: Preserving Java's behavior where bp/extra are always
                        // overwritten before the bp==0 check. See StructuralVariantsProcessor.java:L779-L786
                        bp = m.base_position;
                        extra = m.matched_sequence;
                        if bp == 0 {
                            m = self.find_match_rev(
                                scv_seq,
                                cp,
                                dir,
                                Configuration::SEED_2 as usize,
                                0,
                            );
                            bp = m.base_position;
                            extra = m.matched_sequence;
                        }
                        if bp == 0 {
                            continue;
                        }
                        scv_var = Some(scv.var.clone());
                    }

                    softp = cp;

                    if (dir == 1 && (bp - inv.mend).abs() < min_cluster_dist)
                        || (dir == -1 && (bp - inv.mstart).abs() < min_cluster_dist)
                    {
                        break;
                    }
                }
                if bp == 0 {
                    continue;
                }
            }

            let Some(scv_var) = scv_var else {
                continue;
            };

            if side == InversionSide::End5 {
                if dir == -1 {
                    bp -= 1;
                }
            } else if dir == 1 {
                bp += 1;
                if bp != 0 {
                    softp -= 1;
                }
            } else {
                softp -= 1;
            }

            if side == InversionSide::End3 {
                std::mem::swap(&mut bp, &mut softp);
            }

            if (dir == -1 && side == InversionSide::End5)
                || (dir == 1 && side == InversionSide::End3)
            {
                while self.get_ref_base(softp).is_some()
                    && self.get_ref_base(bp).is_some()
                    && self.get_ref_base(softp).map(Self::complement_base_u8)
                        == self.get_ref_base(bp)
                {
                    softp += 1;
                    if softp != 0 {
                        bp -= 1;
                    }
                }
            }

            while self.get_ref_base(softp - 1).is_some()
                && self.get_ref_base(bp + 1).is_some()
                && self.get_ref_base(softp - 1).map(Self::complement_base_u8)
                    == self.get_ref_base(bp + 1)
            {
                softp -= 1;
                if softp != 0 {
                    bp += 1;
                }
            }

            let mlen_abs = inv.mlen.abs() as f64;
            let ratio = if mlen_abs == 0.0 {
                f64::INFINITY
            } else {
                (bp - softp) as f64 / mlen_abs
            };

            if diag {
                eprintln!(
                    "[DIAG find_inv_sub]   FINAL CHECK: bp={} softp={} diff={} mlen={} ratio={:.4} check=(bp>softp={} diff>150={} ratio<1.5={})",
                    bp, softp, bp - softp, inv.mlen, ratio,
                    bp > softp, bp - softp > 150, ratio < 1.5
                );
            }

            if !(bp > softp && bp - softp > 150 && ratio < 1.5) {
                if diag { eprintln!("[DIAG find_inv_sub]   FAILED final check -> continue"); }
                continue;
            }

            if diag {
                eprintln!("[DIAG find_inv_sub]   SUCCESS: creating INV at softp={} bp={}", softp, bp);
            }

            let len = bp - softp + 1;
            let flank = Configuration::SVFLANK as i64;

            let mut rc = RevComplementor::new();
            let ins = if len - 2 * flank <= 0 {
                let seq_rc = rc.reverse_complement(&self.join_ref(softp, bp));
                String::from_utf8_lossy(&seq_rc).to_string()
            } else {
                let ins5 = rc
                    .reverse_complement(&self.join_ref(bp - flank + 1, bp))
                    .to_vec();
                let ins3 = rc
                    .reverse_complement(&self.join_ref(softp, softp + flank - 1))
                    .to_vec();
                format!(
                    "{}<inv{}>{}",
                    String::from_utf8_lossy(&ins5),
                    len - 2 * flank,
                    String::from_utf8_lossy(&ins3),
                )
            };

            let mut ins_final = ins;
            if dir == 1 && !extra.is_empty() {
                let rc_extra = rc.reverse_complement(&extra).to_vec();
                ins_final = format!("{}{}", String::from_utf8_lossy(&rc_extra), ins_final);
            } else if dir == -1 && !extra.is_empty() {
                ins_final = format!("{}{}", ins_final, String::from_utf8_lossy(&extra));
            }

            let gt = format!("-{}^{}", len, ins_final);
            Self::add_sv_counts(
                &mut data.non_insertion_variants,
                &mut data.sv_counts,
                softp,
                inv.vars_count,
                scv_var.alt_depth,
                1,
            );

            let vref = Self::get_or_create_variation(&mut data.non_insertion_variants, softp, &gt);
            vref.pstd = true;
            vref.qstd = true;

            adj_cnt_from_variant(vref, &scv_var);

            if dir == -1 {
                if let Some(ref_base) = self.get_ref_base(softp) {
                    let ref_key = (ref_base as char).to_string();
                    if let Some(pos_map) = data.non_insertion_variants.get_mut(&softp) {
                        if let Some(reference_var) =
                            Self::get_variation_mut_by_key_string(pos_map, &ref_key)
                        {
                            sub_cnt_from_variant(reference_var, &scv_var);
                        }
                    }
                }
            }

            let cov = data
                .ref_coverage
                .get(&(softp - 1))
                .copied()
                .unwrap_or(inv.vars_count);
            data.ref_coverage.insert(softp, cov);

            Self::mark_inv_cluster_used(data, kind, idx);

            if dir == 1 {
                if let Some(scv) = data.soft_clips_3end.get_mut(&source_softp) {
                    scv.mark_used();
                }
            } else if let Some(scv) = data.soft_clips_5end.get_mut(&source_softp) {
                scv.mark_used();
            }

            let mut dels5: HashMap<i64, InnerMap<String, usize>, LibDefaultHasher> =
                Default::default();
            let mut del_map: InnerMap<String, usize> = Default::default();
            del_map.insert(gt.clone(), inv.vars_count);
            dels5.insert(softp, del_map);

            let realigner = VariantRealigner::new_with_context(
                Arc::clone(&self.reference_seq),
                Arc::clone(&self.reference_seed),
                self.ref_start,
                self.chromosome.clone(),
                self.bam_paths.clone(),
            )
            .with_reference_fallback(
                Arc::clone(&self.original_reference_seq),
                self.original_ref_start,
            );
            append_rss_stage_log_if_enabled(region, "sv_find_inv_before_process_deletions");
            realigner.process_deletions(data, &dels5);
            append_rss_stage_log_if_enabled(region, "sv_find_inv_after_process_deletions");
            return;
        }
    }

    /// Ported from: `com.astrazeneca.vardict.modules.StructuralVariantsProcessor.findsv()`
    /// Java source: `StructuralVariantsProcessor.java:L899-L1196`
    fn find_svs_del_candidates(
        &mut self,
        data: &mut RealignedVariationData,
        region: Option<&Region>,
    ) {
        let minr = instance().conf.minr;
        let region_bounds = region
            .map(|current_region| (current_region.start() as i64, current_region.end() as i64));

        let mut tmp5: Vec<SortPositionSoftClip> = data
            .soft_clips_5end
            .iter()
            .filter(|(_, sclip)| !sclip.used())
            .filter(|(position, _)| {
                region_bounds.map_or(true, |(start, end)| {
                    **position >= start && **position <= end
                })
            })
            .map(|(position, sclip)| SortPositionSoftClip {
                position: *position,
                count: sclip.var.alt_depth,
            })
            .collect();
        tmp5.sort_by(|a, b| b.count.cmp(&a.count));

        for tuple5 in tmp5 {
            let p5 = tuple5.position;
            let cnt5 = tuple5.count;
            if cnt5 < minr {
                break;
            }

            let used_val = data
                .soft_clips_5end
                .get(&p5)
                .map(|s| s.used())
                .unwrap_or(true);
            let softp2sv_val = Self::is_softp2sv_first_used(data, p5);

            if Self::should_trace_findsv_candidate(p5, None) {
                event!(
                    Level::TRACE,
                    phase = "findsv_5_softclip",
                    p5,
                    cnt5,
                    used = used_val,
                    softp2sv_used = softp2sv_val,
                );
            }

            if used_val {
                continue;
            }
            if softp2sv_val {
                continue;
            }

            let seq = {
                let Some(sc5v) = data.soft_clips_5end.get_mut(&p5) else {
                    continue;
                };
                self.find_conseq(sc5v)
            };
            if Self::should_trace_findsv_candidate(p5, None) {
                event!(
                    Level::TRACE,
                    phase = "findsv_5_conseq",
                    p5,
                    cnt5,
                    seq_len = seq.len(),
                    seq = %String::from_utf8_lossy(&seq),
                );
            }
            if seq.is_empty() || seq.len() < Configuration::SEED_2 as usize {
                continue;
            }

            let m = self.find_match(&seq, p5, -1, Configuration::SEED_1 as usize, 3);
            let mut bp = m.base_position;
            if Self::should_trace_findsv_candidate(p5, Some(bp)) {
                event!(
                    Level::TRACE,
                    phase = "findsv_5_direct_match",
                    p5,
                    cnt5,
                    bp,
                );
            }
            event!(Level::DEBUG, phase = "findsv_5_candidate", p5, cnt5, bp,);
            if bp != 0 {
                if bp < p5 {
                    let pairs_data = Self::check_pairs(
                        bp,
                        p5,
                        &mut data.svfdel,
                        &mut data.svrdel,
                        data.max_read_length as i64,
                    );
                    if pairs_data.pairs == 0 {
                        event!(Level::DEBUG, phase = "findsv_5_pairs_zero", p5, bp,);
                        continue;
                    }

                    let p5_adj = p5 - 1;
                    let bp_adj = bp + 1;
                    let dellen = p5_adj - bp_adj + 1;
                    if dellen <= 0 {
                        continue;
                    }

                    let del_key = format!("-{}", dellen);
                    let vref = Self::get_or_create_variation(
                        &mut data.non_insertion_variants,
                        bp_adj,
                        &del_key,
                    );
                    vref.alt_depth = 0;

                    Self::add_sv_counts(
                        &mut data.non_insertion_variants,
                        &mut data.sv_counts,
                        bp_adj,
                        pairs_data.pairs,
                        cnt5,
                        1,
                    );

                    if !data.ref_coverage.contains_key(&bp_adj) {
                        data.ref_coverage.insert(bp_adj, pairs_data.pairs + cnt5);
                    }
                    if let Some(cov_p5) = data.ref_coverage.get(&(p5_adj + 1)).copied() {
                        let cov_bp = data.ref_coverage.get(&bp_adj).copied().unwrap_or(0);
                        if cov_bp < cov_p5 {
                            data.ref_coverage.insert(bp_adj, cov_p5);
                        }
                    }

                    if let Some(sc5v) = data.soft_clips_5end.get(&p5) {
                        let variation = Self::get_or_create_variation(
                            &mut data.non_insertion_variants,
                            bp_adj,
                            &del_key,
                        );
                        adj_cnt_from_variant(variation, &sc5v.var);
                    }

                    let mut tmp = Variant::default();
                    tmp.alt_depth = pairs_data.pairs;
                    tmp.high_qual_read_cnt = pairs_data.pairs;
                    tmp.alt_depth_fwd = pairs_data.pairs / 2;
                    tmp.alt_depth_rev = pairs_data.pairs - pairs_data.pairs / 2;
                    tmp.mean_pos = pairs_data.pmean;
                    tmp.mean_qual = pairs_data.qmean;
                    tmp.mean_mapq = pairs_data.q_mean;
                    tmp.nm = pairs_data.nm;

                    let variation = Self::get_or_create_variation(
                        &mut data.non_insertion_variants,
                        bp_adj,
                        &del_key,
                    );
                    adj_cnt_from_variant(variation, &tmp);
                } else {
                    // candidate duplication
                }
            } else {
                // Java: StructuralVariantsProcessor.java ~L978 — single findMatchRev with SEED_1/MM=3
                let m_rev = self.find_match_rev(&seq, p5, -1, Configuration::SEED_1 as usize, 3);
                bp = m_rev.base_position;
                let extra = m_rev.matched_sequence;
                if Self::should_trace_findsv_candidate(p5, Some(bp)) {
                    let java_len_check = bp > p5 && bp - p5 > 150;
                    event!(
                        Level::TRACE,
                        phase = "findsv_5_reverse_match",
                        p5,
                        cnt5,
                        bp,
                        extra_len = extra.len(),
                        extra = %String::from_utf8_lossy(&extra),
                        java_len_check,
                    );
                }
                if bp == 0 {
                    continue;
                }
                let flank_ok = (bp - p5).abs() > Configuration::SVFLANK as i64;
                if Self::should_trace_findsv_candidate(p5, Some(bp)) {
                    event!(
                        Level::TRACE,
                        phase = "findsv_5_length_gate",
                        p5,
                        bp,
                        flank_ok,
                        flank_gap = (bp - p5).abs(),
                        java_len_check = bp > p5 && bp - p5 > 150,
                    );
                }
                if !flank_ok {
                    continue;
                }

                let sc5_var = data.soft_clips_5end.get(&p5).map(|sc| sc.var.clone());

                let mut p5_inv = p5;
                if bp <= p5_inv {
                    std::mem::swap(&mut bp, &mut p5_inv);
                }
                bp -= 1;

                while self.get_ref_base(bp + 1).is_some()
                    && self.get_ref_base(p5_inv - 1).is_some()
                    && self.get_ref_base(bp + 1).map(Self::complement_base_u8)
                        == self.get_ref_base(p5_inv - 1)
                {
                    p5_inv -= 1;
                    if p5_inv != 0 {
                        bp += 1;
                    }
                }

                let flank = Configuration::SVFLANK as i64;
                let mut rc = RevComplementor::new();
                let ins5 = rc
                    .reverse_complement(&self.join_ref(bp - flank + 1, bp))
                    .to_vec();
                let ins3 = rc
                    .reverse_complement(&self.join_ref(p5_inv, p5_inv + flank - 1))
                    .to_vec();
                let inv_len = bp - p5_inv + 1;
                let mid = inv_len - ins5.len() as i64 - ins3.len() as i64;

                let mut vn = format!(
                    "-{}^{}<inv{}>{}{}",
                    inv_len,
                    String::from_utf8_lossy(&ins5),
                    mid,
                    String::from_utf8_lossy(&ins3),
                    String::from_utf8_lossy(&extra),
                );
                if mid <= 0 {
                    let tins = rc.reverse_complement(&self.join_ref(p5_inv, bp)).to_vec();
                    vn = format!(
                        "-{}^{}{}",
                        inv_len,
                        String::from_utf8_lossy(&tins),
                        String::from_utf8_lossy(&extra),
                    );
                }
                let _ =
                    Self::get_or_create_variation(&mut data.non_insertion_variants, p5_inv, &vn);
                Self::add_sv_counts(
                    &mut data.non_insertion_variants,
                    &mut data.sv_counts,
                    p5_inv,
                    0,
                    cnt5,
                    0,
                );

                if let Some(sc5_var) = sc5_var {
                    let variation = Self::get_or_create_variation(
                        &mut data.non_insertion_variants,
                        p5_inv,
                        &vn,
                    );
                    adj_cnt_from_variant(variation, &sc5_var);
                }

                if let Some(sc5v) = data.soft_clips_5end.get_mut(&p5) {
                    sc5v.mark_used();
                }

                Self::inc_ref_coverage(&mut data.ref_coverage, p5_inv, cnt5);
                if let Some(bp_cov) = data.ref_coverage.get(&bp).copied() {
                    let p5_cov = data.ref_coverage.get(&p5_inv).copied().unwrap_or(0);
                    if p5_cov < bp_cov {
                        data.ref_coverage.insert(p5_inv, bp_cov);
                    }
                }
            }
        }

        let mut tmp3: Vec<SortPositionSoftClip> = data
            .soft_clips_3end
            .iter()
            .filter(|(_, sclip)| !sclip.used())
            .filter(|(position, _)| {
                region_bounds.map_or(true, |(start, end)| {
                    **position >= start && **position <= end
                })
            })
            .map(|(position, sclip)| SortPositionSoftClip {
                position: *position,
                count: sclip.var.alt_depth,
            })
            .collect();
        tmp3.sort_by(|a, b| b.count.cmp(&a.count));

        for tuple3 in tmp3 {
            let mut p3 = tuple3.position;
            let cnt3 = tuple3.count;
            if cnt3 < minr {
                break;
            }

            let used_val = data
                .soft_clips_3end
                .get(&p3)
                .map(|s| s.used())
                .unwrap_or(true);
            let softp2sv_val = Self::is_softp2sv_first_used(data, p3);

            if Self::should_trace_findsv_candidate(p3, None) {
                event!(
                    Level::TRACE,
                    phase = "findsv_3_softclip",
                    p3,
                    cnt3,
                    used = used_val,
                    softp2sv_used = softp2sv_val,
                );
            }

            if used_val {
                continue;
            }
            if softp2sv_val {
                continue;
            }

            let seq = {
                let Some(sc3v) = data.soft_clips_3end.get_mut(&p3) else {
                    continue;
                };
                self.find_conseq(sc3v)
            };
            if Self::should_trace_findsv_candidate(p3, None) {
                event!(
                    Level::TRACE,
                    phase = "findsv_3_conseq",
                    p3,
                    cnt3,
                    seq_len = seq.len(),
                    seq = %String::from_utf8_lossy(&seq),
                );
            }
            if seq.is_empty() || seq.len() < Configuration::SEED_2 as usize {
                continue;
            }

            let m = self.find_match(&seq, p3, 1, Configuration::SEED_1 as usize, 3);
            let mut bp = m.base_position;
            if Self::should_trace_findsv_candidate(p3, Some(bp)) {
                event!(
                    Level::TRACE,
                    phase = "findsv_3_direct_match",
                    p3,
                    cnt3,
                    bp,
                );
            }
            event!(Level::DEBUG, phase = "findsv_3_candidate", p3, cnt3, bp,);
            if bp != 0 {
                if bp > p3 {
                    let pairs_data = Self::check_pairs(
                        p3,
                        bp,
                        &mut data.svfdel,
                        &mut data.svrdel,
                        data.max_read_length as i64,
                    );
                    if pairs_data.pairs == 0 {
                        event!(Level::DEBUG, phase = "findsv_3_pairs_zero", p3, bp,);
                        continue;
                    }

                    let dellen = bp - p3;
                    bp -= 1;

                    while self.get_ref_base(bp).is_some()
                        && self.get_ref_base(p3 - 1).is_some()
                        && self.get_ref_base(bp) == self.get_ref_base(p3 - 1)
                    {
                        bp -= 1;
                        if bp != 0 {
                            p3 -= 1;
                        }
                    }

                    if dellen <= 0 {
                        continue;
                    }

                    let del_key = format!("-{}", dellen);
                    let vref = Self::get_or_create_variation(
                        &mut data.non_insertion_variants,
                        p3,
                        &del_key,
                    );
                    vref.alt_depth = 0;

                    Self::add_sv_counts(
                        &mut data.non_insertion_variants,
                        &mut data.sv_counts,
                        p3,
                        pairs_data.pairs,
                        cnt3,
                        1,
                    );

                    if !data.ref_coverage.contains_key(&p3) {
                        data.ref_coverage.insert(p3, pairs_data.pairs + cnt3);
                    }
                    if let Some(cov_p3) = data.ref_coverage.get(&p3).copied() {
                        if let Some(cov_bp) = data.ref_coverage.get(&bp).copied() {
                            if cov_bp < cov_p3 {
                                data.ref_coverage.insert(bp, cov_p3);
                            }
                        }
                    }

                    if let Some(sc3v) = data.soft_clips_3end.get(&tuple3.position) {
                        let variation = Self::get_or_create_variation(
                            &mut data.non_insertion_variants,
                            p3,
                            &del_key,
                        );
                        adj_cnt_from_variant(variation, &sc3v.var);
                    }

                    let mut tmp = Variant::default();
                    tmp.alt_depth = pairs_data.pairs;
                    tmp.high_qual_read_cnt = pairs_data.pairs;
                    tmp.alt_depth_fwd = pairs_data.pairs / 2;
                    tmp.alt_depth_rev = pairs_data.pairs - pairs_data.pairs / 2;
                    tmp.mean_pos = pairs_data.pmean;
                    tmp.mean_qual = pairs_data.qmean;
                    tmp.mean_mapq = pairs_data.q_mean;
                    tmp.nm = pairs_data.nm;

                    let variation = Self::get_or_create_variation(
                        &mut data.non_insertion_variants,
                        p3,
                        &del_key,
                    );
                    adj_cnt_from_variant(variation, &tmp);
                } else {
                    // candidate duplication
                }
            } else {
                // Java: StructuralVariantsProcessor.java ~L1114 — single findMatchRev with SEED_1/MM=3
                let m_rev = self.find_match_rev(&seq, p3, 1, Configuration::SEED_1 as usize, 3);
                bp = m_rev.base_position;
                let extra = m_rev.matched_sequence;
                if Self::should_trace_findsv_candidate(p3, Some(bp)) {
                    let java_len_check = bp > p3 && bp - p3 > 150;
                    event!(
                        Level::TRACE,
                        phase = "findsv_3_reverse_match",
                        p3,
                        cnt3,
                        bp,
                        extra_len = extra.len(),
                        extra = %String::from_utf8_lossy(&extra),
                        java_len_check,
                    );
                }
                if bp == 0 {
                    continue;
                }
                let flank_ok = (bp - p3).abs() > Configuration::SVFLANK as i64;
                if Self::should_trace_findsv_candidate(p3, Some(bp)) {
                    event!(
                        Level::TRACE,
                        phase = "findsv_3_length_gate",
                        p3,
                        bp,
                        flank_ok,
                        flank_gap = (bp - p3).abs(),
                        java_len_check = bp > p3 && bp - p3 > 150,
                    );
                }
                if !flank_ok {
                    continue;
                }

                let sc3_var = data.soft_clips_3end.get(&p3).map(|sc| sc.var.clone());

                if bp < p3 {
                    std::mem::swap(&mut bp, &mut p3);
                    p3 += 1;
                    bp -= 1;
                }

                while self.get_ref_base(bp + 1).is_some()
                    && self.get_ref_base(p3 - 1).is_some()
                    && self.get_ref_base(bp + 1).map(Self::complement_base_u8)
                        == self.get_ref_base(p3 - 1)
                {
                    p3 -= 1;
                    if p3 != 0 {
                        bp += 1;
                    }
                }

                let flank = Configuration::SVFLANK as i64;
                let mut rc = RevComplementor::new();
                let ins5 = rc
                    .reverse_complement(&self.join_ref(bp - flank + 1, bp))
                    .to_vec();
                let ins3 = rc
                    .reverse_complement(&self.join_ref(p3, p3 + flank - 1))
                    .to_vec();
                let inv_len = bp - p3 + 1;
                let mid = inv_len - 2 * flank;

                let mut vn = format!(
                    "-{}^{}{}<inv{}>{}",
                    inv_len,
                    String::from_utf8_lossy(&extra),
                    String::from_utf8_lossy(&ins5),
                    mid,
                    String::from_utf8_lossy(&ins3),
                );
                if mid <= 0 {
                    let tins = rc.reverse_complement(&self.join_ref(p3, bp)).to_vec();
                    vn = format!(
                        "-{}^{}{}",
                        inv_len,
                        String::from_utf8_lossy(&extra),
                        String::from_utf8_lossy(&tins),
                    );
                }

                let _ = Self::get_or_create_variation(&mut data.non_insertion_variants, p3, &vn);
                Self::add_sv_counts(
                    &mut data.non_insertion_variants,
                    &mut data.sv_counts,
                    p3,
                    0,
                    cnt3,
                    0,
                );

                if let Some(sc3_var) = sc3_var {
                    let variation =
                        Self::get_or_create_variation(&mut data.non_insertion_variants, p3, &vn);
                    adj_cnt_from_variant(variation, &sc3_var);
                }

                if let Some(sc3v) = data.soft_clips_3end.get_mut(&p3) {
                    sc3v.mark_used();
                }

                Self::inc_ref_coverage(&mut data.ref_coverage, p3, cnt3);
                if let Some(bp_cov) = data.ref_coverage.get(&bp).copied() {
                    let p3_cov = data.ref_coverage.get(&p3).copied().unwrap_or(0);
                    if p3_cov < bp_cov {
                        data.ref_coverage.insert(p3, bp_cov);
                    }
                }
            }
        }
    }

    /// Java parity helper for `SOFTP2SV{p}->[0].used` checks used by `findsv()`.
    ///
    /// Java builds `SOFTP2SV` by soft-position, sorts each bucket by `varsCount` descending,
    /// then skips candidate processing when the first item is already used.
    /// Rust computes the same check on demand from SV cluster vectors.
    fn is_softp2sv_first_used(data: &RealignedVariationData, softp: i64) -> bool {
        if let Some(is_used) = data.softp2sv_first_used.get(&softp) {
            return *is_used;
        }

        let mut best: Option<(usize, bool)> = None;

        let mut consider = |sv: &SoftClip| {
            if sv.softp as i64 != softp {
                return;
            }
            let candidate = (sv.var.alt_depth, sv.used());
            match best {
                None => best = Some(candidate),
                Some((current_depth, _)) => {
                    if candidate.0 > current_depth {
                        best = Some(candidate);
                    }
                }
            }
        };

        for sv in &data.svfinv3 {
            consider(sv);
        }
        for sv in &data.svrinv3 {
            consider(sv);
        }
        for sv in &data.svfinv5 {
            consider(sv);
        }
        for sv in &data.svrinv5 {
            consider(sv);
        }
        for sv in &data.svfdel {
            consider(sv);
        }
        for sv in &data.svrdel {
            consider(sv);
        }
        for sv in &data.svfdup {
            consider(sv);
        }
        for sv in &data.svrdup {
            consider(sv);
        }

        best.map(|(_, used)| used).unwrap_or(false)
    }

    fn should_trace_findsv_candidate(softp: i64, bp: Option<i64>) -> bool {
        let Ok(raw_target) = env::var("VARDICT_TRACE_FINDSV_POS") else {
            return false;
        };
        let Ok(target) = raw_target.parse::<i64>() else {
            return false;
        };
        let radius = env::var("VARDICT_TRACE_FINDSV_RADIUS")
            .ok()
            .and_then(|value| value.parse::<i64>().ok())
            .unwrap_or(100);

        if (softp - target).abs() <= radius {
            return true;
        }

        bp.is_some_and(|candidate_bp| (candidate_bp - target).abs() <= radius)
    }

    /// Ported from: `com.astrazeneca.vardict.modules.StructuralVariantsProcessor.findDELdisc()`
    /// Java source: `StructuralVariantsProcessor.java:L1198-L1338`
    fn find_del_disc(&mut self, data: &mut RealignedVariationData) {
        let min_dist = 8 * data.max_read_length as i64;
        let minr = instance().conf.minr;

        for idx in 0..data.svfdel.len() {
            let (
                used,
                vars_count,
                end,
                mstart,
                mean_mapq,
                mean_qual,
                mean_pos,
                nm,
                softp,
                del_mlen,
            ) = {
                let del = &data.svfdel[idx];
                (
                    del.used(),
                    del.var.alt_depth,
                    del.end,
                    del.mstart,
                    del.var.mean_mapq,
                    del.var.mean_qual,
                    del.var.mean_pos,
                    del.var.nm,
                    del.softp,
                    del.mlen,
                )
            };

            if used || vars_count < minr + 5 {
                continue;
            }
            if !data.splice.is_empty() && i64::from(del_mlen).abs() < 250_000 {
                continue;
            }
            if mstart <= end + min_dist {
                continue;
            }
            if vars_count == 0 || mean_mapq / vars_count as f64 <= Configuration::DISCPAIRQUAL {
                continue;
            }

            let mlen = mstart - end - data.max_read_length as i64 / (vars_count + 1) as i64;
            if !(mlen > 0 && mlen > min_dist) {
                continue;
            }

            let mut bp = end + (data.max_read_length as i64 / (vars_count + 1) as i64) / 2;
            if softp != 0 {
                bp = softp as i64;
            }

            self.ensure_reference_span(bp - 150, bp + 150);
            let ext = mlen.min(1000);
            self.loaded_regions.push((bp - 150 - ext, bp + 150 + ext));

            let del_key = format!("-{}", mlen);
            let vref =
                Self::get_or_create_variation(&mut data.non_insertion_variants, bp, &del_key);
            vref.alt_depth = 0;

            let splits = data
                .soft_clips_3end
                .get(&(end + 1))
                .map(|s| s.var.alt_depth)
                .unwrap_or(0)
                + data
                    .soft_clips_5end
                    .get(&mstart)
                    .map(|s| s.var.alt_depth)
                    .unwrap_or(0);
            Self::add_sv_counts(
                &mut data.non_insertion_variants,
                &mut data.sv_counts,
                bp,
                vars_count,
                splits,
                1,
            );

            let mut tv = Variant::default();
            tv.alt_depth = 2 * vars_count;
            tv.high_qual_read_cnt = 2 * vars_count;
            tv.alt_depth_fwd = vars_count;
            tv.alt_depth_rev = vars_count;
            tv.mean_qual = 2.0 * mean_qual;
            tv.mean_pos = 2.0 * mean_pos;
            tv.mean_mapq = 2.0 * mean_mapq;
            tv.nm = 2.0 * nm;
            let variation =
                Self::get_or_create_variation(&mut data.non_insertion_variants, bp, &del_key);
            adj_cnt_from_variant(variation, &tv);

            if !data.ref_coverage.contains_key(&bp) {
                data.ref_coverage.insert(bp, 2 * vars_count);
            }

            if let Some(del) = data.svfdel.get_mut(idx) {
                del.mark_used();
            }
            Self::mark_sv(end, mstart, &mut data.svrdel, data.max_read_length as i64);
        }

        for idx in 0..data.svrdel.len() {
            let (
                used,
                vars_count,
                start,
                del_mstart,
                del_mend,
                mean_mapq,
                mean_qual,
                mean_pos,
                nm,
                softp,
                del_mlen,
            ) = {
                let del = &data.svrdel[idx];
                (
                    del.used(),
                    del.var.alt_depth,
                    del.start,
                    del.mstart,
                    del.mend,
                    del.var.mean_mapq,
                    del.var.mean_qual,
                    del.var.mean_pos,
                    del.var.nm,
                    del.softp,
                    del.mlen,
                )
            };

            if used || vars_count < minr + 5 {
                continue;
            }
            if !data.splice.is_empty() && i64::from(del_mlen).abs() < 250_000 {
                continue;
            }
            if start <= del_mend + min_dist {
                continue;
            }
            if vars_count == 0 || mean_mapq / vars_count as f64 <= Configuration::DISCPAIRQUAL {
                continue;
            }

            let mlen = start - del_mend - data.max_read_length as i64 / (vars_count + 1) as i64;
            if !(mlen > 0 && mlen > min_dist) {
                continue;
            }

            let bp = del_mend + (data.max_read_length as i64 / (vars_count + 1) as i64) / 2;

            self.ensure_reference_span(bp - 150, bp + 150);
            let ext = mlen.min(1000);
            self.loaded_regions.push((bp - 150 - ext, bp + 150 + ext));

            let del_key = format!("-{}", mlen);
            let vref =
                Self::get_or_create_variation(&mut data.non_insertion_variants, bp, &del_key);
            vref.alt_depth = 0;

            let splits = data
                .soft_clips_3end
                .get(&(del_mend + 1))
                .map(|s| s.var.alt_depth)
                .unwrap_or(0)
                + data
                    .soft_clips_5end
                    .get(&start)
                    .map(|s| s.var.alt_depth)
                    .unwrap_or(0);
            Self::add_sv_counts(
                &mut data.non_insertion_variants,
                &mut data.sv_counts,
                bp,
                vars_count,
                splits,
                1,
            );

            if softp != 0 {
                if let Some(sc) = data.soft_clips_5end.get_mut(&(softp as i64)) {
                    sc.mark_used();
                }
            }

            let mut tv = Variant::default();
            tv.alt_depth = 2 * vars_count;
            tv.high_qual_read_cnt = 2 * vars_count;
            tv.alt_depth_fwd = vars_count;
            tv.alt_depth_rev = vars_count;
            tv.mean_qual = 2.0 * mean_qual;
            tv.mean_pos = 2.0 * mean_pos;
            tv.mean_mapq = 2.0 * mean_mapq;
            tv.nm = 2.0 * nm;
            let variation =
                Self::get_or_create_variation(&mut data.non_insertion_variants, bp, &del_key);
            adj_cnt_from_variant(variation, &tv);

            if !data.ref_coverage.contains_key(&bp) {
                data.ref_coverage.insert(bp, 2 * vars_count);
            }
            if let Some(start_cov) = data.ref_coverage.get(&start).copied() {
                let bp_cov = data.ref_coverage.get(&bp).copied().unwrap_or(0);
                if bp_cov < start_cov {
                    data.ref_coverage.insert(bp, start_cov);
                }
            }

            if let Some(del) = data.svrdel.get_mut(idx) {
                del.mark_used();
            }
            self.ensure_reference_span(del_mstart - 100, del_mend + 100);
            self.loaded_regions.push((del_mstart - 300, del_mend + 300));
            Self::mark_sv(
                del_mend,
                start,
                &mut data.svfdel,
                data.max_read_length as i64,
            );
        }
    }

    fn find_inv_disc(&mut self, data: &mut RealignedVariationData, region: Option<&Region>) {
        let minr = instance().conf.minr;
        let mut rev_complementor = RevComplementor::new();
        let diag = std::env::var("VARDICT_DIAG_INV_SUB").is_ok();

        append_rss_stage_log_if_enabled(region, "sv_find_inv_disc_enter");

        for f_idx in 0..data.svfinv5.len() {
            let (f_used, cnt, me, ms, end, start, nm, pmean, qmean, q_mean) = {
                let invf5 = &data.svfinv5[f_idx];
                (
                    invf5.used(),
                    invf5.var.alt_depth,
                    invf5.mend,
                    invf5.mstart,
                    invf5.end,
                    invf5.start,
                    invf5.var.nm,
                    invf5.var.mean_pos,
                    invf5.var.mean_qual,
                    invf5.var.mean_mapq,
                )
            };

            if diag {
                eprintln!("[DIAG find_inv_disc] f_idx={} f_used={} cnt={} end={} me={} ms={} start={}", f_idx, f_used, cnt, end, me, ms, start);
            }

            if f_used || cnt == 0 {
                continue;
            }
            if q_mean / cnt as f64 <= Configuration::DISCPAIRQUAL {
                continue;
            }

            for r_idx in 0..data.svrinv5.len() {
                let (r_used, rcnt, rstart, rms, rnm, rpmean, rqmean, r_q_mean) = {
                    let invr5 = &data.svrinv5[r_idx];
                    (
                        invr5.used(),
                        invr5.var.alt_depth,
                        invr5.start,
                        invr5.mstart,
                        invr5.var.nm,
                        invr5.var.mean_pos,
                        invr5.var.mean_qual,
                        invr5.var.mean_mapq,
                    )
                };

                if diag {
                    eprintln!("[DIAG find_inv_disc]   r_idx={} r_used={} rcnt={} rstart={} rms={}", r_idx, r_used, rcnt, rstart, rms);
                }

                if r_used || rcnt == 0 {
                    continue;
                }
                if r_q_mean / rcnt as f64 <= Configuration::DISCPAIRQUAL {
                    continue;
                }
                if cnt + rcnt <= minr + 5 {
                    continue;
                }
                if !Self::is_overlap(end, me, rstart, rms, data.max_read_length as i64) {
                    continue;
                }

                let bp = ((end + rstart) / 2).abs();
                let pe = ((me + rms) / 2).abs();
                if pe < bp {
                    continue;
                }
                append_rss_stage_log_if_enabled(region, "sv_find_inv_disc_5_before_ref_fetch");

                let len = pe - bp + 1;
                if len <= 0 {
                    continue;
                }

                let flank = Configuration::SVFLANK as i64;
                let ins = if len - 2 * flank <= 0 {
                    let rc = rev_complementor.reverse_complement(&self.join_ref(bp, pe));
                    String::from_utf8_lossy(&rc).to_string()
                } else {
                    let ins5 = rev_complementor
                        .reverse_complement(&self.join_ref(bp, bp + flank - 1))
                        .to_vec();
                    let ins3 = rev_complementor
                        .reverse_complement(&self.join_ref(pe - flank + 1, pe))
                        .to_vec();
                    format!(
                        "{}<inv{}>{}",
                        String::from_utf8_lossy(&ins3),
                        len - 2 * flank,
                        String::from_utf8_lossy(&ins5),
                    )
                };
                append_rss_stage_log_if_enabled(region, "sv_find_inv_disc_5_after_ref_fetch");

                let inv_key = format!("-{}^{}", len, ins);

                if diag {
                    eprintln!("[DIAG find_inv_disc]   CREATING INV5: f_idx={} r_idx={} bp={} pe={} len={} cnt={} rcnt={}", f_idx, r_idx, bp, pe, len, cnt, rcnt);
                }

                let vref =
                    Self::get_or_create_variation(&mut data.non_insertion_variants, bp, &inv_key);
                vref.pstd = true;
                vref.qstd = true;

                let mut tmp = Variant::default();
                tmp.alt_depth = cnt + rcnt;
                tmp.high_qual_read_cnt = cnt + rcnt;
                tmp.alt_depth_fwd = cnt;
                tmp.alt_depth_rev = rcnt;
                tmp.mean_qual = qmean + rqmean;
                tmp.mean_pos = pmean + rpmean;
                tmp.mean_mapq = q_mean + r_q_mean;
                tmp.nm = nm + rnm;
                adj_cnt_from_variant(vref, &tmp);

                let splits = data
                    .soft_clips_5end
                    .get(&start)
                    .map(|s| s.var.alt_depth)
                    .unwrap_or(0)
                    + data
                        .soft_clips_5end
                        .get(&ms)
                        .map(|s| s.var.alt_depth)
                        .unwrap_or(0);
                Self::add_sv_counts(
                    &mut data.non_insertion_variants,
                    &mut data.sv_counts,
                    bp,
                    cnt,
                    splits,
                    1,
                );

                if !data.ref_coverage.contains_key(&bp) {
                    data.ref_coverage.insert(bp, 2 * cnt);
                }

                if let Some(invf5) = data.svfinv5.get_mut(f_idx) {
                    invf5.mark_used();
                }
                if let Some(invr5) = data.svrinv5.get_mut(r_idx) {
                    invr5.mark_used();
                }
                event!(
                    Level::DEBUG,
                    phase = "find_inv_disc_emit_5",
                    bp,
                    pe,
                    cnt,
                    rcnt,
                    inv_key = %inv_key,
                );
                append_rss_stage_log_if_enabled(region, "sv_find_inv_disc_5_emit");
                Self::mark_sv(bp, pe, &mut data.svfinv3, data.max_read_length as i64);
                Self::mark_sv(bp, pe, &mut data.svrinv3, data.max_read_length as i64);
            }
        }

        for f_idx in 0..data.svfinv3.len() {
            let (f_used, cnt, me, end, nm, pmean, qmean, q_mean) = {
                let invf3 = &data.svfinv3[f_idx];
                (
                    invf3.used(),
                    invf3.var.alt_depth,
                    invf3.mend,
                    invf3.end,
                    invf3.var.nm,
                    invf3.var.mean_pos,
                    invf3.var.mean_qual,
                    invf3.var.mean_mapq,
                )
            };

            if diag {
                eprintln!("[DIAG find_inv_disc] INV3 f_idx={} f_used={} cnt={} end={} me={}", f_idx, f_used, cnt, end, me);
            }

            if f_used || cnt == 0 {
                continue;
            }

            for r_idx in 0..data.svrinv3.len() {
                let (r_used, rcnt, rstart, rms, rnm, rpmean, rqmean, r_q_mean) = {
                    let invr3 = &data.svrinv3[r_idx];
                    (
                        invr3.used(),
                        invr3.var.alt_depth,
                        invr3.start,
                        invr3.mstart,
                        invr3.var.nm,
                        invr3.var.mean_pos,
                        invr3.var.mean_qual,
                        invr3.var.mean_mapq,
                    )
                };

                if r_used || rcnt == 0 {
                    continue;
                }
                if r_q_mean / rcnt as f64 <= Configuration::DISCPAIRQUAL {
                    continue;
                }
                if cnt + rcnt <= minr + 5 {
                    continue;
                }
                if !Self::is_overlap(me, end, rms, rstart, data.max_read_length as i64) {
                    continue;
                }

                let pe = ((end + rstart) / 2).abs();
                let bp = ((me + rms) / 2).abs();
                if pe < bp {
                    continue;
                }
                append_rss_stage_log_if_enabled(region, "sv_find_inv_disc_3_before_ref_fetch");

                let len = pe - bp + 1;
                if len <= 0 {
                    continue;
                }

                let flank = Configuration::SVFLANK as i64;
                let ins = if len - 2 * flank <= 0 {
                    let rc = rev_complementor.reverse_complement(&self.join_ref(bp, pe));
                    String::from_utf8_lossy(&rc).to_string()
                } else {
                    let ins5 = rev_complementor
                        .reverse_complement(&self.join_ref(bp, bp + flank - 1))
                        .to_vec();
                    let ins3 = rev_complementor
                        .reverse_complement(&self.join_ref(pe - flank + 1, pe))
                        .to_vec();
                    format!(
                        "{}<inv{}>{}",
                        String::from_utf8_lossy(&ins3),
                        len - 2 * flank,
                        String::from_utf8_lossy(&ins5),
                    )
                };
                append_rss_stage_log_if_enabled(region, "sv_find_inv_disc_3_after_ref_fetch");

                let inv_key = format!("-{}^{}", len, ins);
                let vref =
                    Self::get_or_create_variation(&mut data.non_insertion_variants, bp, &inv_key);
                vref.pstd = true;
                vref.qstd = true;

                let mut tmp = Variant::default();
                tmp.alt_depth = cnt + rcnt;
                tmp.high_qual_read_cnt = cnt + rcnt;
                tmp.alt_depth_fwd = cnt;
                tmp.alt_depth_rev = rcnt;
                tmp.mean_qual = qmean + rqmean;
                tmp.mean_pos = pmean + rpmean;
                tmp.mean_mapq = q_mean + r_q_mean;
                tmp.nm = nm + rnm;
                adj_cnt_from_variant(vref, &tmp);

                let splits = data
                    .soft_clips_3end
                    .get(&(end + 1))
                    .map(|s| s.var.alt_depth)
                    .unwrap_or(0)
                    + data
                        .soft_clips_3end
                        .get(&(me + 1))
                        .map(|s| s.var.alt_depth)
                        .unwrap_or(0);
                Self::add_sv_counts(
                    &mut data.non_insertion_variants,
                    &mut data.sv_counts,
                    bp,
                    cnt,
                    splits,
                    1,
                );

                if !data.ref_coverage.contains_key(&bp) {
                    data.ref_coverage.insert(bp, 2 * cnt);
                }

                if let Some(invf3) = data.svfinv3.get_mut(f_idx) {
                    invf3.mark_used();
                }
                if let Some(invr3) = data.svrinv3.get_mut(r_idx) {
                    invr3.mark_used();
                }
                event!(
                    Level::DEBUG,
                    phase = "find_inv_disc_emit_3",
                    bp,
                    pe,
                    cnt,
                    rcnt,
                    inv_key = %inv_key,
                );
                append_rss_stage_log_if_enabled(region, "sv_find_inv_disc_3_emit");
                Self::mark_sv(bp, pe, &mut data.svfinv5, data.max_read_length as i64);
                Self::mark_sv(bp, pe, &mut data.svrinv5, data.max_read_length as i64);
            }
        }
    }

    fn find_dup_disc(&mut self, data: &mut RealignedVariationData) {
        let minr = instance().conf.minr;

        for idx in 0..data.svfdup.len() {
            let (used, ms, me, cnt, mut end, _start, pmean, qmean, q_mean, nm, softp, soft_map) = {
                let dup = &data.svfdup[idx];
                (
                    dup.used(),
                    dup.mstart,
                    dup.mend,
                    dup.var.alt_depth,
                    dup.end,
                    dup.start,
                    dup.var.mean_pos,
                    dup.var.mean_qual,
                    dup.var.mean_mapq,
                    dup.var.nm,
                    dup.softp,
                    dup.soft.clone(),
                )
            };

            if used || cnt < minr + 5 || cnt == 0 {
                continue;
            }
            if q_mean / cnt as f64 <= Configuration::DISCPAIRQUAL {
                continue;
            }

            let read_len_adj = data.max_read_length as i64 / cnt as i64;
            let mut mlen = end - ms + read_len_adj;
            let mut bp = ms - read_len_adj / 2;
            let mut pe = end;

            self.ensure_reference_span(bp - 150, bp + 150);
            self.load_uncovered_reference_coverage(data, ms - 200, me + 200);

            let mut cntf = cnt;
            let mut cntr = cnt;
            let mut qmeanf = qmean;
            let mut qmeanr = qmean;
            let mut qmean_f_mapq = q_mean;
            let mut qmean_r_mapq = q_mean;
            let mut pmeanf = pmean;
            let mut pmeanr = pmean;
            let mut nmf = nm;
            let mut nmr = nm;

            if !soft_map.is_empty() {
                let Some(mut pe_soft) = Self::select_primary_soft_pos(&soft_map) else {
                    continue;
                };
                pe = pe_soft;

                let seq = {
                    let Some(sc3) = data.soft_clips_3end.get_mut(&pe) else {
                        continue;
                    };
                    if sc3.used() {
                        continue;
                    }
                    cntf = sc3.var.alt_depth;
                    qmeanf = sc3.var.mean_qual;
                    qmean_f_mapq = sc3.var.mean_mapq;
                    pmeanf = sc3.var.mean_pos;
                    nmf = sc3.var.nm;
                    self.find_conseq(sc3)
                };

                let mut tbp = self
                    .find_match(&seq, bp, 1, Configuration::SEED_1 as usize, 3)
                    .base_position;

                if tbp != 0 && tbp < pe {
                    if let Some(sc3) = data.soft_clips_3end.get_mut(&pe) {
                        sc3.mark_used();
                    }

                    while self.get_ref_base(pe_soft - 1).is_some()
                        && self.get_ref_base(tbp - 1).is_some()
                        && self.get_ref_base(pe_soft - 1) == self.get_ref_base(tbp - 1)
                    {
                        tbp -= 1;
                        if tbp != 0 {
                            pe_soft -= 1;
                        }
                    }

                    mlen = pe_soft - tbp;
                    bp = tbp;
                    pe = pe_soft - 1;
                    end = pe;

                    if let Some(sc5) = data.soft_clips_5end.get(&bp) {
                        cntr = sc5.var.alt_depth;
                        qmeanr = sc5.var.mean_qual;
                        qmean_r_mapq = sc5.var.mean_mapq;
                        pmeanr = sc5.var.mean_pos;
                        nmr = sc5.var.nm;
                    }
                }
            }

            self.ensure_reference_span(bp - 150, pe + 150);

            let mut ins = self.join_ref(bp, bp + Configuration::SVFLANK as i64 - 1);
            let dup_len = mlen - 2 * Configuration::SVFLANK as i64;
            ins.extend_from_slice(format!("<dup{}>", dup_len).as_bytes());
            ins.extend_from_slice(&self.join_ref(pe - Configuration::SVFLANK as i64 + 1, pe));

            let splits = if softp != 0 {
                data.soft_clips_3end
                    .get(&(softp as i64))
                    .map(|s| s.var.alt_depth)
                    .unwrap_or(0)
            } else {
                0
            };
            Self::add_sv_counts(
                &mut data.non_insertion_variants,
                &mut data.sv_counts,
                bp,
                cnt,
                splits,
                1,
            );

            let tcnt = cntr + cntf;
            let mut tmp = Variant::default();
            tmp.alt_depth = tcnt;
            tmp.high_qual_read_cnt = tcnt;
            tmp.alt_depth_fwd = cntf;
            tmp.alt_depth_rev = cntr;
            tmp.mean_qual = qmeanf + qmeanr;
            tmp.mean_pos = pmeanf + pmeanr;
            tmp.mean_mapq = qmean_f_mapq + qmean_r_mapq;
            tmp.nm = nmf + nmr;

            let iref_key = VarDesc::Ins {
                seq: ins.iter().copied().collect(),
            };
            let iref = data
                .insertion_variants
                .entry(bp)
                .or_default()
                .entry(iref_key)
                .or_default();
            iref.alt_depth = 0;
            adj_cnt_from_variant(iref, &tmp);

            if let Some(dup) = data.svfdup.get_mut(idx) {
                dup.mark_used();
            }

            if !data.ref_coverage.contains_key(&bp) {
                data.ref_coverage.insert(bp, tcnt);
            }
            if let Some(end_cov) = data.ref_coverage.get(&end).copied() {
                let bp_cov = data.ref_coverage.get(&bp).copied().unwrap_or(0);
                if bp_cov < end_cov {
                    data.ref_coverage.insert(bp, end_cov);
                }
            }

            let (clusters, _) =
                Self::mark_dup_sv(bp, pe, &mut data.svrdup, data.max_read_length as i64);
            if clusters != 0 {
                Self::add_sv_counts(
                    &mut data.non_insertion_variants,
                    &mut data.sv_counts,
                    bp,
                    0,
                    0,
                    clusters,
                );
            }
        }

        for idx in 0..data.svrdup.len() {
            let (used, ms, me, cnt, _end, start, pmean, qmean, q_mean, nm, soft_map) = {
                let dup = &data.svrdup[idx];
                (
                    dup.used(),
                    dup.mstart,
                    dup.mend,
                    dup.var.alt_depth,
                    dup.end,
                    dup.start,
                    dup.var.mean_pos,
                    dup.var.mean_qual,
                    dup.var.mean_mapq,
                    dup.var.nm,
                    dup.soft.clone(),
                )
            };

            if used || cnt < minr + 5 || cnt == 0 {
                continue;
            }
            if q_mean / cnt as f64 <= Configuration::DISCPAIRQUAL {
                continue;
            }

            let read_len_adj = data.max_read_length as i64 / cnt as i64;
            let mut mlen = me - start + read_len_adj;
            let mut bp = start - read_len_adj / 2;
            let mut pe = mlen + bp - 1;
            let mut tpe = pe;

            self.ensure_reference_span(pe - 150, pe + 150);
            self.load_uncovered_reference_coverage(data, ms - 200, me + 200);

            let mut cntf = cnt;
            let mut cntr = cnt;
            let mut qmeanf = qmean;
            let mut qmeanr = qmean;
            let mut qmean_f_mapq = q_mean;
            let mut qmean_r_mapq = q_mean;
            let mut pmeanf = pmean;
            let mut pmeanr = pmean;
            let mut nmf = nm;
            let mut nmr = nm;

            if !soft_map.is_empty() {
                let Some(bp_soft) = Self::select_primary_soft_pos(&soft_map) else {
                    continue;
                };
                bp = bp_soft;

                let seq = {
                    let Some(sc5) = data.soft_clips_5end.get_mut(&bp) else {
                        continue;
                    };
                    if sc5.used() {
                        continue;
                    }
                    cntr = sc5.var.alt_depth;
                    qmeanr = sc5.var.mean_qual;
                    qmean_r_mapq = sc5.var.mean_mapq;
                    pmeanr = sc5.var.mean_pos;
                    nmr = sc5.var.nm;
                    self.find_conseq(sc5)
                };

                let tbp = self
                    .find_match(&seq, pe, -1, Configuration::SEED_1 as usize, 3)
                    .base_position;
                if tbp != 0 && tbp > bp {
                    if let Some(sc5) = data.soft_clips_5end.get_mut(&bp) {
                        sc5.mark_used();
                    }

                    pe = tbp;
                    mlen = pe - bp + 1;
                    tpe = pe + 1;
                    while self.get_ref_base(tpe).is_some()
                        && self.get_ref_base(bp + (tpe - pe - 1)).is_some()
                        && self.get_ref_base(tpe) == self.get_ref_base(bp + (tpe - pe - 1))
                    {
                        tpe += 1;
                    }

                    if let Some(sc3) = data.soft_clips_3end.get(&tpe) {
                        cntf = sc3.var.alt_depth;
                        qmeanf = sc3.var.mean_qual;
                        qmean_f_mapq = sc3.var.mean_mapq;
                        pmeanf = sc3.var.mean_pos;
                        nmf = sc3.var.nm;
                    }
                }
            }

            self.ensure_reference_span(bp - 150, pe + 150);

            let mut ins = self.join_ref(bp, bp + Configuration::SVFLANK as i64 - 1);
            let dup_len = mlen - 2 * Configuration::SVFLANK as i64;
            ins.extend_from_slice(format!("<dup{}>", dup_len).as_bytes());
            ins.extend_from_slice(&self.join_ref(pe - Configuration::SVFLANK as i64 + 1, pe));

            let mut splits = data
                .soft_clips_5end
                .get(&bp)
                .map(|s| s.var.alt_depth)
                .unwrap_or(0);
            splits += data
                .soft_clips_3end
                .get(&tpe)
                .map(|s| s.var.alt_depth)
                .unwrap_or(0);
            Self::add_sv_counts(
                &mut data.non_insertion_variants,
                &mut data.sv_counts,
                bp,
                cnt,
                splits,
                1,
            );

            let tcnt = cntr + cntf;
            let mut tmp = Variant::default();
            tmp.alt_depth = tcnt;
            tmp.high_qual_read_cnt = tcnt;
            tmp.alt_depth_fwd = cntf;
            tmp.alt_depth_rev = cntr;
            tmp.mean_qual = qmeanf + qmeanr;
            tmp.mean_pos = pmeanf + pmeanr;
            tmp.mean_mapq = qmean_f_mapq + qmean_r_mapq;
            tmp.nm = nmf + nmr;

            let iref_key = VarDesc::Ins {
                seq: ins.iter().copied().collect(),
            };
            let iref = data
                .insertion_variants
                .entry(bp)
                .or_default()
                .entry(iref_key)
                .or_default();
            iref.alt_depth = 0;
            adj_cnt_from_variant(iref, &tmp);

            if let Some(dup) = data.svrdup.get_mut(idx) {
                dup.mark_used();
            }

            if !data.ref_coverage.contains_key(&bp) {
                data.ref_coverage.insert(bp, tcnt);
            }
            if let Some(me_cov) = data.ref_coverage.get(&me).copied() {
                let bp_cov = data.ref_coverage.get(&bp).copied().unwrap_or(0);
                if bp_cov < me_cov {
                    data.ref_coverage.insert(bp, me_cov);
                }
            }

            let (clusters, _) =
                Self::mark_dup_sv(bp, pe, &mut data.svfdup, data.max_read_length as i64);
            if clusters != 0 {
                Self::add_sv_counts(
                    &mut data.non_insertion_variants,
                    &mut data.sv_counts,
                    bp,
                    0,
                    0,
                    clusters,
                );
            }
        }
    }

    fn inv_cluster_len(data: &RealignedVariationData, kind: InversionClusterKind) -> usize {
        match kind {
            InversionClusterKind::Forward5 => data.svfinv5.len(),
            InversionClusterKind::Reverse5 => data.svrinv5.len(),
            InversionClusterKind::Forward3 => data.svfinv3.len(),
            InversionClusterKind::Reverse3 => data.svrinv3.len(),
        }
    }

    fn inv_cluster_snapshot(
        data: &RealignedVariationData,
        kind: InversionClusterKind,
        idx: usize,
    ) -> Option<InversionClusterSnapshot> {
        let cluster = match kind {
            InversionClusterKind::Forward5 => data.svfinv5.get(idx),
            InversionClusterKind::Reverse5 => data.svrinv5.get(idx),
            InversionClusterKind::Forward3 => data.svfinv3.get(idx),
            InversionClusterKind::Reverse3 => data.svrinv3.get(idx),
        }?;

        Some(InversionClusterSnapshot {
            used: cluster.used(),
            vars_count: cluster.var.alt_depth,
            start: cluster.start,
            end: cluster.end,
            mstart: cluster.mstart,
            mend: cluster.mend,
            mlen: cluster.mlen,
            primary_softp: Self::select_primary_soft_pos(&cluster.soft),
        })
    }

    fn mark_inv_cluster_used(
        data: &mut RealignedVariationData,
        kind: InversionClusterKind,
        idx: usize,
    ) {
        let cluster = match kind {
            InversionClusterKind::Forward5 => data.svfinv5.get_mut(idx),
            InversionClusterKind::Reverse5 => data.svrinv5.get_mut(idx),
            InversionClusterKind::Forward3 => data.svfinv3.get_mut(idx),
            InversionClusterKind::Reverse3 => data.svrinv3.get_mut(idx),
        };
        if let Some(cluster) = cluster {
            cluster.mark_used();
        }
    }

    fn complement_base_u8(base: u8) -> u8 {
        match base.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'N',
        }
    }

    fn find_match_rev(
        &self,
        seq: &[u8],
        position: i64,
        dir: i64,
        seed_len: usize,
        mm: usize,
    ) -> MatchResult {
        self.find_match_rev_internal(seq, position, dir, seed_len, mm, true, false, false)
    }

    fn find_match_rev_findsv(&self, seq: &[u8], position: i64, dir: i64) -> MatchResult {
        for (seed_len, mm, include_shared_reference_fallback) in [
            (Configuration::SEED_1 as usize, 3usize, false),
            (Configuration::SEED_2 as usize, 0usize, false),
            (Configuration::SEED_1 as usize, 3usize, true),
            (Configuration::SEED_2 as usize, 0usize, true),
        ] {
            let result = self.find_match_rev_internal(
                seq,
                position,
                dir,
                seed_len,
                mm,
                true,
                include_shared_reference_fallback,
                include_shared_reference_fallback,
            );
            if result.base_position != 0 {
                return result;
            }
        }

        MatchResult::default()
    }

    fn find_match_rev_internal(
        &self,
        seq: &[u8],
        _position: i64,
        dir: i64,
        seed_len: usize,
        mm: usize,
        include_historical_windows: bool,
        include_shared_reference_fallback: bool,
        force_shared_reference_fallback: bool,
    ) -> MatchResult {
        let diag = std::env::var("VARDICT_DIAG_INV_SUB").is_ok();
        let mut seq_work = seq.to_vec();
        if dir == 1 {
            seq_work.reverse();
        }
        seq_work = seq_work
            .iter()
            .map(|b| Self::complement_base_u8(*b))
            .collect();

        if seq_work.len() < seed_len {
            return MatchResult::default();
        }

        let mut persistent_extra: Vec<u8> = Vec::new();
        let mut diag_seeds_checked = 0u32;
        let mut diag_seeds_skip_multi = 0u32;
        let mut diag_seeds_skip_zero = 0u32;

        for i in (0..=seq_work.len() - seed_len).rev() {
            let seed = &seq_work[i..i + seed_len];
            let seeds = if force_shared_reference_fallback {
                self.seed_positions_with_findsv_scope(
                    seed,
                    include_historical_windows,
                    include_shared_reference_fallback,
                )
            } else {
                self.seed_positions_with_scope(
                    seed,
                    include_historical_windows,
                    include_shared_reference_fallback,
                )
            };
            diag_seeds_checked += 1;
            if seeds.len() != 1 {
                if seeds.is_empty() { diag_seeds_skip_zero += 1; } else { diag_seeds_skip_multi += 1; }
                continue;
            }

            let first_seed = seeds[0];
            let mut bp = if dir == 1 {
                first_seed + seq_work.len() as i64 - i as i64 - 1
            } else {
                first_seed - i as i64
            };

            let initial_match = self.is_match_ref(&seq_work, bp, -dir, mm);
            if diag {
                eprintln!(
                    "[DIAG find_match_rev] SEED HIT: i={} seed={} first_seed={} bp={} match={} pos={} dir={} checked={} skip_zero={} skip_multi={} hist={}",
                    i, String::from_utf8_lossy(seed), first_seed, bp, initial_match, _position, dir,
                    diag_seeds_checked, diag_seeds_skip_zero, diag_seeds_skip_multi,
                    include_historical_windows
                );
            }
            if initial_match {
                return MatchResult {
                    base_position: bp,
                    matched_sequence: persistent_extra,
                };
            }

            let mut window_start = 0usize;
            let mut window_end = seq_work.len();
            let mut eqcnt = 0usize;
            for j in 1..=15 {
                bp -= dir;
                if dir == -1 {
                    window_start += 1;
                } else {
                    window_end = window_end.saturating_sub(1);
                }

                if window_start >= window_end {
                    break;
                }

                let sseq = &seq_work[window_start..window_end];

                let extra = if dir == -1 {
                    let ch0 = sseq[0];
                    if self.is_has_and_not_equals(bp, ch0) {
                        continue;
                    }
                    eqcnt += 1;
                    let ch1 = sseq.get(1).copied().unwrap_or(b'N');
                    if sseq.len() < 2 || self.is_has_and_not_equals(bp + 1, ch1) {
                        continue;
                    }
                    &seq_work[..j]
                } else {
                    let ch_last = sseq[sseq.len() - 1];
                    if self.is_has_and_not_equals(bp, ch_last) {
                        continue;
                    }
                    eqcnt += 1;
                    if sseq.len() < 2 {
                        continue;
                    }
                    let ch_prev = sseq[sseq.len() - 2];
                    if self.is_has_and_not_equals(bp - 1, ch_prev) {
                        continue;
                    }
                    &seq_work[seq_work.len() - j..]
                };

                persistent_extra = extra.to_vec();

                if eqcnt >= 3 && (eqcnt as f64 / j as f64) > 0.5 {
                    break;
                }

                let retry_match = self.is_match_ref(&sseq, bp, -dir, 1);
                if retry_match {
                    return MatchResult {
                        base_position: bp,
                        matched_sequence: extra.to_vec(),
                    };
                }
            }
        }

        MatchResult::default()
    }

    fn inc_ref_coverage(
        ref_coverage: &mut HashMap<i64, usize, LibDefaultHasher>,
        pos: i64,
        cnt: usize,
    ) {
        let entry = ref_coverage.entry(pos).or_insert(0);
        *entry += cnt;
    }

    fn check_pairs(
        start: i64,
        end: i64,
        svfdel: &mut [SoftClip],
        svrdel: &mut [SoftClip],
        max_read_length: i64,
    ) -> PairsData {
        let mut out = PairsData::default();

        for cluster in svfdel.iter_mut().chain(svrdel.iter_mut()) {
            if cluster.used() {
                continue;
            }

            let mut s = (cluster.start + cluster.end) / 2;
            let mut e = (cluster.mstart + cluster.mend) / 2;
            if s > e {
                std::mem::swap(&mut s, &mut e);
            }

            if !Self::is_overlap(start, end, s, e, max_read_length) {
                continue;
            }

            event!(
                Level::DEBUG,
                phase = "check_pairs_overlap",
                query_start = start,
                query_end = end,
                cluster_start = cluster.start,
                cluster_end = cluster.end,
                cluster_mstart = cluster.mstart,
                cluster_mend = cluster.mend,
                cluster_cnt = cluster.var.alt_depth,
            );

            if cluster.var.alt_depth > out.pairs {
                out.pairs = cluster.var.alt_depth;
                out.pmean = cluster.var.mean_pos;
                out.qmean = cluster.var.mean_qual;
                out.q_mean = cluster.var.mean_mapq;
                out.nm = cluster.var.nm;
                event!(
                    Level::DEBUG,
                    phase = "check_pairs_pick",
                    picked_pairs = out.pairs,
                    query_start = start,
                    query_end = end,
                );
            }
            cluster.mark_used();
        }

        out
    }

    fn select_primary_soft_pos(soft: &IndexMap<i64, usize>) -> Option<i64> {
        let mut best: Option<(i64, usize)> = None;
        for (pos, count) in soft.iter() {
            match best {
                None => best = Some((*pos, *count)),
                Some((_, best_count)) if *count > best_count => best = Some((*pos, *count)),
                _ => {}
            }
        }
        best.map(|(pos, _)| pos)
    }

    fn get_or_create_variation<'a>(
        map: &'a mut HashMap<i64, InnerMap<VarDesc, Variant>, LibDefaultHasher>,
        pos: i64,
        key_str: &str,
    ) -> &'a mut Variant {
        let pos_map = map.entry(pos).or_default();
        let key = pos_map
            .keys()
            .find(|k| k.to_key_string() == key_str)
            .cloned()
            .unwrap_or_else(|| VarDesc::Raw {
                desc: key_str.as_bytes().to_vec().into(),
            });
        pos_map.entry(key).or_default()
    }

    fn get_variation_mut_by_key_string<'a>(
        pos_map: &'a mut InnerMap<VarDesc, Variant>,
        key_str: &str,
    ) -> Option<&'a mut Variant> {
        let key = pos_map
            .keys()
            .find(|k| k.to_key_string() == key_str)
            .cloned()?;
        pos_map.get_mut(&key)
    }

    fn add_sv_counts(
        map: &mut HashMap<i64, InnerMap<VarDesc, Variant>, LibDefaultHasher>,
        sv_counts: &mut HashMap<i64, StructuralVariantCounts, LibDefaultHasher>,
        pos: i64,
        pairs: usize,
        splits: usize,
        clusters: usize,
    ) {
        let pos_map = map.entry(pos).or_default();
        let key = pos_map
            .keys()
            .find(|k| matches!(k, VarDesc::Raw { desc } if desc.as_slice() == b"SV"))
            .cloned()
            .unwrap_or_else(|| VarDesc::Raw {
                desc: b"SV".to_vec().into(),
            });
        pos_map.entry(key).or_default();
        let sv = sv_counts.entry(pos).or_default();
        sv.pairs += pairs;
        sv.splits += splits;
        sv.clusters += clusters;
    }

    fn mark_sv(start: i64, end: i64, clusters: &mut [SoftClip], rlen: i64) {
        for cluster in clusters.iter_mut() {
            let (start2, end2) = if cluster.start < cluster.mstart {
                (cluster.end, cluster.mstart)
            } else {
                (cluster.mend, cluster.start)
            };
            if Self::is_overlap(start, end, start2, end2, rlen) {
                cluster.mark_used();
            }
        }
    }

    fn mark_dup_sv(start: i64, end: i64, clusters: &mut [SoftClip], rlen: i64) -> (usize, usize) {
        let mut cnt = 0usize;
        let mut pairs = 0usize;

        for cluster in clusters.iter_mut() {
            let (start2, end2) = if cluster.start < cluster.mstart {
                (cluster.start, cluster.mend)
            } else {
                (cluster.mstart, cluster.end)
            };
            if Self::is_overlap(start, end, start2, end2, rlen) {
                cluster.mark_used();
                cnt += 1;
                pairs += cluster.var.alt_depth;
            }
        }

        (cnt, pairs)
    }

    fn join_ref(&self, start: i64, end: i64) -> Vec<u8> {
        if end < start {
            return Vec::new();
        }

        let requested_start = start.max(1);
        let requested_end = end.max(requested_start);

        if !self.reference_seq.is_empty() {
            let current_end = self.ref_start + self.reference_seq.len() as i64 - 1;
            if requested_start >= self.ref_start && requested_end <= current_end {
                let start_idx = (requested_start - self.ref_start) as usize;
                let end_idx = (requested_end - self.ref_start + 1) as usize;
                return self.reference_seq[start_idx..end_idx].to_vec();
            }
        }

        if let (Some(shared_reference), Some(chromosome)) =
            (self.shared_reference.as_ref(), self.chromosome.as_deref())
        {
            if let Some(sequence) = shared_reference
                .get_subseq(chromosome, requested_start as usize, requested_end as usize)
                .map(|seq| seq.to_vec())
            {
                return sequence;
            }
        }

        let mut out = Vec::with_capacity((requested_end - requested_start + 1) as usize);
        for pos in requested_start..=requested_end {
            let Some(base) = self.get_ref_base(pos) else {
                break;
            };
            out.push(base);
        }
        out
    }

    fn is_overlap(start1: i64, end1: i64, start2: i64, end2: i64, rlen: i64) -> bool {
        if start1 >= end2 || start2 >= end1 {
            return false;
        }

        let mut positions = [start1, end1, start2, end2];
        positions.sort_unstable();
        let ins = positions[2] - positions[1];

        if end1 != start1
            && end2 != start2
            && (ins as f64 / (end1 - start1) as f64) > 0.75
            && (ins as f64 / (end2 - start2) as f64) > 0.75
        {
            return true;
        }

        positions[1] - positions[0] + positions[3] - positions[2] < 3 * rlen
    }

    fn find_match(
        &self,
        seq: &[u8],
        position: i64,
        dir: i64,
        seed_len: usize,
        mm: usize,
    ) -> MatchResult {
        self.find_match_internal(seq, position, dir, seed_len, mm, true, false)
    }

    fn find_match_internal(
        &self,
        seq: &[u8],
        _position: i64,
        dir: i64,
        seed_len: usize,
        mm: usize,
        include_historical_windows: bool,
        include_shared_reference_fallback: bool,
    ) -> MatchResult {
        let mut seq_work = seq.to_vec();
        if dir == -1 {
            seq_work.reverse();
        }

        let mut persistent_extra: Vec<u8> = Vec::new();

        if seq_work.len() < seed_len {
            return MatchResult::default();
        }

        for i in (0..=seq_work.len() - seed_len).rev() {
            let seed = &seq_work[i..i + seed_len];
            let seeds = self.seed_positions_with_scope(
                seed,
                include_historical_windows,
                include_shared_reference_fallback,
            );
            if seeds.len() != 1 {
                continue;
            }

            let first_seed = seeds[0];
            let mut bp = if dir == 1 {
                first_seed - i as i64
            } else {
                first_seed + seq_work.len() as i64 - i as i64 - 1
            };

            if self.is_match_ref(&seq_work, bp, dir, mm) {
                persistent_extra.clear();
                let mut mm_idx: i64 = if dir == -1 { -1 } else { 0 };
                loop {
                    let Some(ch) = Self::char_at(&seq_work, mm_idx) else {
                        break;
                    };
                    if self.is_has_and_not_equals(bp, ch) {
                        persistent_extra.push(ch);
                        bp += dir;
                        mm_idx += dir;
                    } else {
                        break;
                    }
                }
                if !persistent_extra.is_empty() && dir == -1 {
                    persistent_extra.reverse();
                }
                return MatchResult {
                    base_position: bp,
                    matched_sequence: persistent_extra,
                };
            }

            let mut sseq = seq_work.clone();
            let mut eqcnt = 0usize;
            for ii in 1..=15 {
                bp += dir;
                sseq = if dir == 1 {
                    Self::substr_bytes(&sseq, 1, None)
                } else {
                    Self::substr_bytes(&sseq, 0, Some(-1))
                };
                if sseq.is_empty() {
                    break;
                }

                let extra = if dir == 1 {
                    let Some(ch0) = sseq.first().copied() else {
                        continue;
                    };
                    if self.is_has_and_not_equals(bp, ch0) {
                        continue;
                    }
                    eqcnt += 1;
                    let ch1 = sseq.get(1).copied().unwrap_or(b'N');
                    if sseq.len() < 2 || self.is_has_and_not_equals(bp + 1, ch1) {
                        continue;
                    }
                    Self::substr_bytes(&seq_work, 0, Some(ii as i64))
                } else {
                    let Some(ch_last) = Self::char_at(&sseq, -1) else {
                        continue;
                    };
                    if self.is_has_and_not_equals(bp, ch_last) {
                        continue;
                    }
                    eqcnt += 1;
                    let Some(ch_prev) = Self::char_at(&sseq, -2) else {
                        continue;
                    };
                    if self.is_has_and_not_equals(bp - 1, ch_prev) {
                        continue;
                    }
                    Self::substr_bytes(&seq_work, -(ii as i64), None)
                };

                persistent_extra = extra.clone();

                if eqcnt >= 3 && (eqcnt as f64 / ii as f64) > 0.5 {
                    break;
                }

                if self.is_match_ref(&sseq, bp, dir, 1) {
                    return MatchResult {
                        base_position: bp,
                        matched_sequence: extra,
                    };
                }
            }
        }

        MatchResult::default()
    }

    fn is_has_and_not_equals(&self, pos: i64, ch: u8) -> bool {
        self.get_ref_base(pos).map(|b| b != ch).unwrap_or(false)
    }

    /// Ported from: `com.astrazeneca.vardict.modules.VariationRealigner.ismatchref()`
    /// Java source: `VariationRealigner.java:L2914-L2931`
    fn is_match_ref(&self, seq: &[u8], position: i64, dir: i64, mm: usize) -> bool {
        if seq.is_empty() {
            return false;
        }
        let mut mismatches = 0usize;
        for n in 0..seq.len() {
            let ref_pos = position + dir * n as i64;
            let ref_base = match self.get_ref_base(ref_pos) {
                Some(b) => b,
                None => return false,
            };
            let seq_base = if dir == 1 {
                seq[n]
            } else {
                Self::char_at(seq, -(n as i64 + 1)).unwrap_or(b'N')
            };
            if seq_base != ref_base {
                mismatches += 1;
            }
        }
        mismatches <= mm && (mismatches as f64 / seq.len() as f64) < 0.15
    }

    fn char_at(seq: &[u8], idx: i64) -> Option<u8> {
        if seq.is_empty() {
            return None;
        }
        let len = seq.len() as i64;
        let pos = if idx < 0 { len + idx } else { idx };
        if pos < 0 || pos >= len {
            None
        } else {
            Some(seq[pos as usize])
        }
    }

    fn substr_bytes(seq: &[u8], begin: i64, len: Option<i64>) -> Vec<u8> {
        let seq_len = seq.len() as i64;
        let mut b = begin;
        if b < 0 {
            b = seq_len + b;
        }
        if b > seq_len {
            return Vec::new();
        }

        match len {
            None => {
                if b < 0 {
                    b = 0;
                }
                let b_usize = b as usize;
                seq.get(b_usize..).unwrap_or(&[]).to_vec()
            }
            Some(l) if l > 0 => {
                if b < 0 {
                    return Vec::new();
                }
                let b_usize = b as usize;
                let end = (b + l).min(seq_len).max(b);
                seq.get(b_usize..end as usize).unwrap_or(&[]).to_vec()
            }
            Some(0) => Vec::new(),
            Some(l) => {
                if b < 0 {
                    return Vec::new();
                }
                let b_usize = b as usize;
                let end = seq_len + l;
                if end < b {
                    return Vec::new();
                }
                seq.get(b_usize..end as usize).unwrap_or(&[]).to_vec()
            }
        }
    }

    /// Adjust SNV counts from short soft-clipped reads
    ///
    /// This always runs, even when SV detection is disabled.
    /// It looks at short soft clips (≤5 bp) and if the first base matches
    /// a known SNV at the adjacent position, adds the soft clip counts to that SNV.
    fn adj_snv(&self, data: &mut RealignedVariationData) {
        // Process 5' end soft clips
        self.adj_snv_5end(data);

        // Process 3' end soft clips
        self.adj_snv_3end(data);
    }

    /// Adjust SNV counts from 5' soft clips
    fn adj_snv_5end(&self, data: &mut RealignedVariationData) {
        // Collect positions to process (avoid borrowing issues)
        let positions: Vec<i64> = data.soft_clips_5end.keys().cloned().collect();

        for position in positions {
            let sclip = match data.soft_clips_5end.get_mut(&position) {
                Some(sc) => sc,
                None => continue,
            };

            // Skip if already used
            if sclip.used() {
                continue;
            }

            // Find consensus sequence from soft clip
            let seq = self.find_conseq(sclip);

            // Only process short soft clips (≤5 bp)
            if seq.len() > 5 {
                continue;
            }

            if seq.is_empty() {
                continue;
            }

            // Get first base
            let bp = seq[0];

            let vars_count = sclip.var.alt_depth;
            let high_qual_cnt = sclip.var.high_qual_read_cnt;
            let low_qual_cnt = sclip.var.low_qual_read_cnt;
            let mean_pos = sclip.var.mean_pos;
            let mean_qual = sclip.var.mean_qual;
            let mean_mapq = sclip.var.mean_mapq;
            let nm = sclip.var.nm;
            let fwd_cnt = sclip.var.alt_depth_fwd;
            let rev_cnt = sclip.var.alt_depth_rev;

            // Check if there's a matching SNV at the previous position
            let prev_pos = position - 1;

            // Create the variant description key for the first base
            let var_key = VarDesc::snv_key(bp);

            // Check if this base exists as a variant at the previous position
            if let Some(var_map) = data.non_insertion_variants.get_mut(&prev_pos) {
                if var_map.contains_key(&var_key) {
                    // Additional check: if seq length > 1, verify second base matches reference
                    // Java: if reference base is missing, treat as mismatch and skip
                    if seq.len() > 1 {
                        match self.get_ref_base(position - 2) {
                            Some(ref_base) => {
                                if ref_base != seq[1] {
                                    continue;
                                }
                            }
                            None => {
                                continue;
                            }
                        }
                    }

                    // Adjust the variant counts
                    if let Some(variant) = var_map.get_mut(&var_key) {
                        adj_cnt(
                            variant,
                            vars_count,
                            high_qual_cnt,
                            low_qual_cnt,
                            mean_pos,
                            mean_qual,
                            mean_mapq,
                            nm,
                            fwd_cnt,
                            rev_cnt,
                        );
                    }

                    // Increment reference coverage
                    *data.ref_coverage.entry(prev_pos).or_insert(0) += vars_count;
                }
            }
        }
    }

    /// Adjust SNV counts from 3' soft clips
    fn adj_snv_3end(&self, data: &mut RealignedVariationData) {
        // Collect positions to process (avoid borrowing issues)
        let positions: Vec<i64> = data.soft_clips_3end.keys().cloned().collect();

        for position in positions {
            let sclip = match data.soft_clips_3end.get_mut(&position) {
                Some(sc) => sc,
                None => continue,
            };

            // Skip if already used
            if sclip.used() {
                continue;
            }

            // Find consensus sequence from soft clip
            let seq = self.find_conseq(sclip);

            // Only process short soft clips (≤5 bp)
            if seq.len() > 5 {
                continue;
            }

            if seq.is_empty() {
                continue;
            }

            // Get first base
            let bp = seq[0];

            let vars_count = sclip.var.alt_depth;
            let high_qual_cnt = sclip.var.high_qual_read_cnt;
            let low_qual_cnt = sclip.var.low_qual_read_cnt;
            let mean_pos = sclip.var.mean_pos;
            let mean_qual = sclip.var.mean_qual;
            let mean_mapq = sclip.var.mean_mapq;
            let nm = sclip.var.nm;
            let fwd_cnt = sclip.var.alt_depth_fwd;
            let rev_cnt = sclip.var.alt_depth_rev;

            // Create the variant description key for the first base
            let var_key = VarDesc::snv_key(bp);

            // Check if this base exists as a variant at this position
            if let Some(var_map) = data.non_insertion_variants.get_mut(&position) {
                if var_map.contains_key(&var_key) {
                    // Additional check: if seq length > 1, verify second base matches reference
                    // Java: if reference base is missing, treat as mismatch and skip
                    if seq.len() > 1 {
                        match self.get_ref_base(position + 1) {
                            Some(ref_base) => {
                                if ref_base != seq[1] {
                                    continue;
                                }
                            }
                            None => {
                                continue;
                            }
                        }
                    }

                    // Adjust the variant counts
                    if let Some(variant) = var_map.get_mut(&var_key) {
                        adj_cnt(
                            variant,
                            vars_count,
                            high_qual_cnt,
                            low_qual_cnt,
                            mean_pos,
                            mean_qual,
                            mean_mapq,
                            nm,
                            fwd_cnt,
                            rev_cnt,
                        );
                    }

                    // Increment reference coverage
                    *data.ref_coverage.entry(position).or_insert(0) += vars_count;
                }
            }
        }
    }

    /// Find consensus sequence from soft clip data
    ///
    /// Simplified version of Java findconseq()
    fn find_conseq(&self, sclip: &mut SoftClip) -> Vec<u8> {
        crate::variants::var_utils::find_conseq(sclip, 0)
    }

    fn ensure_conseq<'a>(&self, sclip: &'a mut SoftClip) -> &'a [u8] {
        if !sclip.consensus_seq_is_set() {
            let _ = crate::variants::var_utils::find_conseq(sclip, 0);
        }
        sclip.consensus_seq()
    }

    /// Get reference base at a position (1-based)
    fn get_ref_base(&self, pos: i64) -> Option<u8> {
        Self::get_ref_base_from_window(&self.reference_seq, self.ref_start, pos).or_else(|| {
            if self.original_ref_start == self.ref_start
                && Arc::ptr_eq(&self.original_reference_seq, &self.reference_seq)
            {
                self.shared_reference_base(pos)
            } else {
                Self::get_ref_base_from_window(
                    &self.original_reference_seq,
                    self.original_ref_start,
                    pos,
                )
                .or_else(|| self.shared_reference_base(pos))
            }
        })
    }

    fn get_ref_base_from_window(reference_seq: &[u8], ref_start: i64, pos: i64) -> Option<u8> {
        if pos < ref_start {
            return None;
        }
        let idx = (pos - ref_start) as usize;
        reference_seq.get(idx).copied()
    }

    fn seed_positions(&self, seed: &[u8], include_historical_windows: bool) -> Vec<i64> {
        self.seed_positions_with_scope(seed, include_historical_windows, true)
    }

    fn seed_positions_with_scope(
        &self,
        seed: &[u8],
        include_historical_windows: bool,
        include_shared_reference_fallback: bool,
    ) -> Vec<i64> {
        let mut positions = Vec::new();
        let reference_window_shifted = !(self.original_ref_start == self.ref_start
            && Arc::ptr_eq(&self.original_reference_seed, &self.reference_seed));
        let has_remote_reference_context =
            reference_window_shifted || !self.historical_reference_windows.is_empty();

        if let Some(current) = self.reference_seed.get(seed) {
            positions.extend(current.iter().copied());
        }

        if reference_window_shifted {
            if let Some(original) = self.original_reference_seed.get(seed) {
                for pos in original {
                    if !positions.contains(pos) {
                        positions.push(*pos);
                    }
                }
            }
        }

        if include_historical_windows {
            self.extend_seed_positions_from_historical_windows(seed, &mut positions);

            // The plain candidate-creation reverse matcher should stay aligned
            // with Java's mutable REF hash: current plus previously loaded
            // windows are visible, but not an unrestricted chromosome-wide
            // seed scan. Reserve the shared-reference fallback for the wider
            // historical matcher used by inversion rescue paths.
            if include_shared_reference_fallback
                && has_remote_reference_context
                && positions.len() <= 1
            {
                self.extend_seed_positions_from_shared_reference(seed, &mut positions);
            }
        }

        positions
    }

    fn seed_positions_with_findsv_scope(
        &self,
        seed: &[u8],
        include_historical_windows: bool,
        include_shared_reference_fallback: bool,
    ) -> Vec<i64> {
        let mut positions = self.seed_positions_with_scope(seed, include_historical_windows, false);

        if include_historical_windows && include_shared_reference_fallback && positions.len() <= 1 {
            self.extend_seed_positions_from_shared_reference(seed, &mut positions);
        }

        positions
    }

    fn ensure_reference_span(&mut self, start: i64, end: i64) {
        if start > end {
            return;
        }

        let original_end = if self.original_reference_seq.is_empty() {
            self.original_ref_start - 1
        } else {
            self.original_ref_start + self.original_reference_seq.len() as i64 - 1
        };

        let requested_start = start.max(1);
        let requested_end = end.max(requested_start);

        let current_start = self.ref_start;
        let current_end = if self.reference_seq.is_empty() {
            self.ref_start - 1
        } else {
            self.ref_start + self.reference_seq.len() as i64 - 1
        };

        if !self.reference_seq.is_empty()
            && requested_start >= current_start
            && requested_end <= current_end
        {
            return;
        }

        // If the original reference covers the requested span, restore it instead of
        // loading from shared_reference. This avoids unnecessary remote lookups.
        if !self.original_reference_seq.is_empty()
            && requested_start >= self.original_ref_start
            && requested_end <= original_end
        {
            self.restore_original_reference_window_preserving_history();
            return;
        }

        // Keep prior remote windows only as coordinates for later seed lookups while still
        // replacing the active contiguous window to avoid the old memory blow-up.
        self.remember_current_reference_window();

        let Some(chromosome) = self.chromosome.as_deref() else {
            return;
        };
        let Some(shared_reference) = self.shared_reference.as_ref() else {
            return;
        };

        let new_start = requested_start;
        let new_end = requested_end;

        let Some(sequence) = shared_reference
            .get_subseq(chromosome, new_start as usize, new_end as usize)
            .map(|seq| seq.to_vec())
        else {
            return;
        };

        let chr_len = instance().chr_lens.get(chromosome).copied();
        let mut reference = Reference::new_with_start(sequence, new_start);
        reference.build_seed_map(new_end, chr_len);

        self.reference_seq = reference.ref_seq;
        self.reference_seed = reference.seed;
        self.ref_start = new_start;
    }

    fn restore_original_reference_window(&mut self) {
        self.reference_seq = Arc::clone(&self.original_reference_seq);
        self.reference_seed = Arc::clone(&self.original_reference_seed);
        self.ref_start = self.original_ref_start;
    }

    fn restore_original_reference_window_preserving_history(&mut self) {
        self.remember_current_reference_window();
        self.restore_original_reference_window();
    }

    fn is_span_loaded(&self, start: i64, end: i64) -> bool {
        if start > end {
            return false;
        }

        let mut windows = self.historical_reference_windows.clone();

        if !self.original_reference_seq.is_empty() {
            windows.push(ReferenceWindow {
                start: self.original_ref_start,
                end: self.original_ref_start + self.original_reference_seq.len() as i64 - 1,
            });
        }

        if !self.reference_seq.is_empty() {
            windows.push(ReferenceWindow {
                start: self.ref_start,
                end: self.ref_start + self.reference_seq.len() as i64 - 1,
            });
        }

        windows
            .iter()
            .any(|window| start >= window.start && end <= window.end)
    }

    /// Returns true when [start, end] is fully within the original reference
    /// or a previously tracked loaded region. Does NOT check the current
    /// ephemeral reference window - only stable regions.
    fn is_region_loaded(&self, start: i64, end: i64) -> bool {
        if start > end {
            return false;
        }
        if !self.original_reference_seq.is_empty() {
            let original_end =
                self.original_ref_start + self.original_reference_seq.len() as i64 - 1;
            if start >= self.original_ref_start && end <= original_end {
                return true;
            }
        }
        self.loaded_regions
            .iter()
            .any(|&(loaded_start, loaded_end)| start >= loaded_start && end <= loaded_end)
    }

    fn is_reference_coverage_loaded(&self, start: i64, end: i64) -> bool {
        if start > end {
            return false;
        }

        if !self.original_reference_seq.is_empty() {
            let original_window = ReferenceWindow {
                start: self.original_ref_start,
                end: self.original_ref_start + self.original_reference_seq.len() as i64 - 1,
            };
            if start >= original_window.start && end <= original_window.end {
                return true;
            }
        }

        self.loaded_reference_coverage_windows
            .iter()
            .any(|window| start >= window.start && end <= window.end)
    }

    fn uncovered_reference_coverage_spans(&self, start: i64, end: i64) -> Vec<ReferenceWindow> {
        if start > end {
            return Vec::new();
        }

        if self.is_reference_coverage_loaded(start, end) {
            return Vec::new();
        }

        let mut covered_windows = self.loaded_reference_coverage_windows.clone();
        if !self.original_reference_seq.is_empty() {
            covered_windows.push(ReferenceWindow {
                start: self.original_ref_start,
                end: self.original_ref_start + self.original_reference_seq.len() as i64 - 1,
            });
        }

        if covered_windows.is_empty() {
            return vec![ReferenceWindow { start, end }];
        }

        covered_windows.sort_by_key(|window| window.start);

        let mut merged_covered: Vec<ReferenceWindow> = Vec::with_capacity(covered_windows.len());
        for window in covered_windows {
            if let Some(last) = merged_covered.last_mut() {
                if window.start <= last.end + 1 {
                    last.end = last.end.max(window.end);
                    continue;
                }
            }
            merged_covered.push(window);
        }

        let mut cursor = start;
        let mut uncovered = Vec::new();
        for window in merged_covered {
            if window.end < cursor {
                continue;
            }
            if window.start > end {
                break;
            }
            if window.start > cursor {
                uncovered.push(ReferenceWindow {
                    start: cursor,
                    end: (window.start - 1).min(end),
                });
            }
            cursor = cursor.max(window.end.saturating_add(1));
            if cursor > end {
                break;
            }
        }

        if cursor <= end {
            uncovered.push(ReferenceWindow { start: cursor, end });
        }

        uncovered
    }

    fn load_uncovered_reference_coverage(
        &mut self,
        data: &mut RealignedVariationData,
        start: i64,
        end: i64,
    ) {
        for window in self.uncovered_reference_coverage_spans(start, end) {
            let realigner = VariantRealigner::new_with_context(
                Arc::clone(&self.reference_seq),
                Arc::clone(&self.reference_seed),
                self.ref_start,
                self.chromosome.clone(),
                self.bam_paths.clone(),
            );
            realigner.load_partial_ref_coverage(data, window.start, window.end);
            self.remember_reference_coverage_window(window.start, window.end);
        }
    }

    fn shared_reference_base(&self, pos: i64) -> Option<u8> {
        let chromosome = self.chromosome.as_deref()?;
        let shared_reference = self.shared_reference.as_ref()?;
        let pos = usize::try_from(pos).ok()?;
        shared_reference.get_base(chromosome, pos)
    }

    fn remember_current_reference_window(&mut self) {
        if self.reference_seq.is_empty() {
            return;
        }

        if self.ref_start == self.original_ref_start
            && Arc::ptr_eq(&self.reference_seq, &self.original_reference_seq)
        {
            return;
        }

        let end = self.ref_start + self.reference_seq.len() as i64 - 1;
        self.remember_reference_window(self.ref_start, end);
    }

    fn remember_reference_window(&mut self, start: i64, end: i64) {
        if start > end {
            return;
        }

        let window = ReferenceWindow { start, end };
        if self.historical_reference_windows.contains(&window) {
            return;
        }

        self.historical_reference_windows.push(window);
    }

    fn record_del_rightseq_from_loaded_regions(&mut self, data: &RealignedVariationData) {
        let mut candidates: Vec<(i64, String)> = Vec::new();
        for (&position, variations) in &data.non_insertion_variants {
            for key in variations.keys() {
                let del_key = key.to_key_string();
                if del_key.starts_with('-') {
                    candidates.push((position, del_key));
                }
            }
        }
        for (position, del_key) in candidates {
            if del_key.starts_with('-') {
                let numeric_part = del_key[1..]
                    .split(|c: char| c == '#' || c == '^' || c == '&')
                    .next()
                    .unwrap_or("");
                if let Some(del_len) = numeric_part.parse::<i64>().ok() {
                    let right_start = position + del_len;
                    let right_end = right_start + 19;
                    if self.is_region_loaded(right_start, right_end) {
                        self.historical_del_rightseq_variants
                            .entry(position)
                            .or_default()
                            .insert(del_key);
                    }
                }
            }
        }
    }

    fn remember_reference_coverage_window(&mut self, start: i64, end: i64) {
        if start > end {
            return;
        }

        self.loaded_reference_coverage_windows
            .push(ReferenceWindow { start, end });
        self.loaded_reference_coverage_windows
            .sort_by_key(|window| window.start);

        let mut merged_windows: Vec<ReferenceWindow> =
            Vec::with_capacity(self.loaded_reference_coverage_windows.len());
        for window in self.loaded_reference_coverage_windows.drain(..) {
            if let Some(last) = merged_windows.last_mut() {
                if window.start <= last.end + 1 {
                    last.end = last.end.max(window.end);
                    continue;
                }
            }
            merged_windows.push(window);
        }

        self.loaded_reference_coverage_windows = merged_windows;
    }

    fn extend_seed_positions_from_historical_windows(&self, seed: &[u8], positions: &mut Vec<i64>) {
        if seed.is_empty() {
            return;
        }

        let Some(shared_reference) = self.shared_reference.as_ref() else {
            return;
        };
        let Some(chromosome) = self.chromosome.as_deref() else {
            return;
        };

        for window in &self.historical_reference_windows {
            let Some(sequence) =
                shared_reference.get_subseq(chromosome, window.start as usize, window.end as usize)
            else {
                continue;
            };

            if sequence.len() < seed.len() {
                continue;
            }

            for offset in 0..=sequence.len() - seed.len() {
                if &sequence[offset..offset + seed.len()] != seed {
                    continue;
                }

                let pos = window.start + offset as i64;
                if !positions.contains(&pos) {
                    positions.push(pos);
                }
            }
        }
    }

    /// Java's StructuralVariantsProcessor only looks up seeds from the bounded
    /// `REF.seed` HashMap (the loaded reference window).  It never does a
    /// chromosome-wide linear scan.  The previous Rust implementation scanned
    /// the full chromosome (~249 MB for chr1) for every seed lookup that missed
    /// the local window maps — O(chromosome_len × seed_len) per call —
    /// causing 28× slowdown vs Java in pileup mode.
    ///
    /// Disabled to match Java behaviour.  If a seed is not present in the
    /// current or historical window seed maps, it is simply "not found," the
    /// same result Java produces.
    fn extend_seed_positions_from_shared_reference(&self, _seed: &[u8], _positions: &mut Vec<i64>) {
        // Intentionally a no-op — see doc comment above.
    }
}

#[derive(Default)]
struct MatchResult {
    base_position: i64,
    matched_sequence: Vec<u8>,
}

struct SortPositionSoftClip {
    position: i64,
    count: usize,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum InversionSide {
    End5,
    End3,
}

#[derive(Clone, Copy, Debug)]
enum InversionClusterKind {
    Forward5,
    Reverse5,
    Forward3,
    Reverse3,
}

struct InversionClusterSnapshot {
    used: bool,
    vars_count: usize,
    start: i64,
    end: i64,
    mstart: i64,
    mend: i64,
    mlen: i32,
    primary_softp: Option<i64>,
}

#[derive(Default)]
struct PairsData {
    pairs: usize,
    pmean: f64,
    qmean: f64,
    q_mean: f64,
    nm: f64,
}

/// Adjust variant counts by adding values from soft clip
///
/// Equivalent to Java VariationUtils.adjCnt()
fn adj_cnt(
    variant: &mut Variant,
    vars_count: usize,
    high_qual_cnt: usize,
    low_qual_cnt: usize,
    mean_pos: f64,
    mean_qual: f64,
    mean_mapq: f64,
    nm: f64,
    fwd_cnt: usize,
    rev_cnt: usize,
) {
    variant.alt_depth += vars_count;
    variant.extra_cnt += vars_count;
    variant.high_qual_read_cnt += high_qual_cnt;
    variant.low_qual_read_cnt += low_qual_cnt;
    variant.mean_pos += mean_pos;
    variant.mean_qual += mean_qual;
    variant.mean_mapq += mean_mapq;
    variant.nm += nm;
    variant.alt_depth_fwd += fwd_cnt;
    variant.alt_depth_rev += rev_cnt;
    variant.pstd = true;
    variant.qstd = true;
}

fn adj_cnt_from_variant(dest: &mut Variant, src: &Variant) {
    adj_cnt(
        dest,
        src.alt_depth,
        src.high_qual_read_cnt,
        src.low_qual_read_cnt,
        src.mean_pos,
        src.mean_qual,
        src.mean_mapq,
        src.nm,
        src.alt_depth_fwd,
        src.alt_depth_rev,
    );
}

fn sub_cnt_from_variant(dest: &mut Variant, src: &Variant) {
    dest.alt_depth = dest.alt_depth.saturating_sub(src.alt_depth);
    dest.high_qual_read_cnt = dest
        .high_qual_read_cnt
        .saturating_sub(src.high_qual_read_cnt);
    dest.low_qual_read_cnt = dest.low_qual_read_cnt.saturating_sub(src.low_qual_read_cnt);
    dest.mean_pos = (dest.mean_pos - src.mean_pos).max(0.0);
    dest.mean_qual = (dest.mean_qual - src.mean_qual).max(0.0);
    dest.mean_mapq = (dest.mean_mapq - src.mean_mapq).max(0.0);
    dest.nm = (dest.nm - src.nm).max(0.0);
    dest.alt_depth_fwd = dest.alt_depth_fwd.saturating_sub(src.alt_depth_fwd);
    dest.alt_depth_rev = dest.alt_depth_rev.saturating_sub(src.alt_depth_rev);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::shared_reference::{ChromosomeData, SharedReference};
    use crate::scopedata::global_read_only_scope::{GlobalReadOnlyScope, INSTANCE};
    use std::fs;

    #[test]
    fn test_processor_creation() {
        let processor =
            StructuralVariantsProcessor::new(b"ACGTACGT".to_vec(), Default::default(), 100);
        assert_eq!(processor.ref_start, 100);
    }

    #[test]
    fn test_get_ref_base() {
        let processor =
            StructuralVariantsProcessor::new(b"ACGTACGT".to_vec(), Default::default(), 100);

        assert_eq!(processor.get_ref_base(100), Some(b'A'));
        assert_eq!(processor.get_ref_base(101), Some(b'C'));
        assert_eq!(processor.get_ref_base(107), Some(b'T'));
        assert_eq!(processor.get_ref_base(108), None); // Out of bounds
        assert_eq!(processor.get_ref_base(99), None); // Before start
    }

    #[test]
    fn test_into_reference_preserves_original_region_reference() {
        let mut processor =
            StructuralVariantsProcessor::new(b"ACGTACGT".to_vec(), Default::default(), 100);

        processor.reference_seq = Arc::new(b"GT".to_vec());
        processor.reference_seed = Arc::new(Default::default());
        processor.ref_start = 102;

        let reference = processor.into_reference();

        assert_eq!(reference.region_start, 100);
        assert_eq!(*reference.ref_seq, b"ACGTACGT".to_vec());
    }

    #[test]
    fn test_restore_original_reference_window_preserving_history_keeps_remote_window_history() {
        let mut processor =
            StructuralVariantsProcessor::new(b"ACGTACGT".to_vec(), Default::default(), 100);

        processor.reference_seq = Arc::new(b"TTGGAACC".to_vec());
        processor.reference_seed = Arc::new(Default::default());
        processor.ref_start = 500;

        processor.restore_original_reference_window_preserving_history();

        assert_eq!(
            processor.historical_reference_windows,
            vec![ReferenceWindow {
                start: 500,
                end: 507,
            }]
        );
        assert_eq!(processor.ref_start, 100);
        assert_eq!(processor.reference_seq.as_ref(), b"ACGTACGT");
    }

    #[test]
    fn test_restore_original_reference_window_drops_unrecorded_remote_window_history() {
        let mut processor =
            StructuralVariantsProcessor::new(b"ACGTACGT".to_vec(), Default::default(), 100);

        processor.reference_seq = Arc::new(b"TTGGAACC".to_vec());
        processor.reference_seed = Arc::new(Default::default());
        processor.ref_start = 500;

        processor.restore_original_reference_window();

        assert!(processor.historical_reference_windows.is_empty());
        assert_eq!(processor.ref_start, 100);
        assert_eq!(processor.reference_seq.as_ref(), b"ACGTACGT");
    }

    #[test]
    fn test_ensure_reference_span_restores_original_window_for_in_region_requests() {
        let mut processor =
            StructuralVariantsProcessor::new(b"ACGTACGT".to_vec(), Default::default(), 100);

        processor.reference_seq = Arc::new(b"GT".to_vec());
        processor.reference_seed = Arc::new(Default::default());
        processor.ref_start = 102;

        processor.ensure_reference_span(103, 105);

        assert_eq!(processor.ref_start, 100);
        assert_eq!(processor.reference_seq.as_ref(), b"ACGTACGT");
    }

    #[test]
    fn test_seed_positions_search_historical_windows_via_shared_reference() {
        let chrom = "testchr".to_string();
        let full_sequence = b"AAAAAAAAAAAAAAAAAAAAACGTACGTACGAACGTACGATTTTTTTTTTTTTTTTTTTT";

        let mut original_reference = Reference::new_with_start(full_sequence[..20].to_vec(), 1);
        original_reference.build_seed_map(20, Some(full_sequence.len()));

        let mut current_reference = Reference::new_with_start(full_sequence[40..].to_vec(), 41);
        current_reference.build_seed_map(full_sequence.len() as i64, Some(full_sequence.len()));

        let mut chromosomes: HashMap<String, ChromosomeData, LibDefaultHasher> = Default::default();
        chromosomes.insert(
            chrom.clone(),
            ChromosomeData {
                sequence: Arc::new(full_sequence.to_vec()),
                length: full_sequence.len(),
            },
        );

        let shared_reference = Arc::new(SharedReference {
            chromosomes,
            chromosome_names: vec![chrom.clone()],
            total_size: full_sequence.len(),
        });

        let mut processor = StructuralVariantsProcessor::new_with_context(
            Arc::clone(&original_reference.ref_seq),
            Arc::clone(&original_reference.seed),
            1,
            Some(chrom),
            Vec::new(),
            Some(shared_reference),
        );

        processor.reference_seq = Arc::clone(&current_reference.ref_seq);
        processor.reference_seed = Arc::clone(&current_reference.seed);
        processor.ref_start = 41;
        processor
            .historical_reference_windows
            .push(ReferenceWindow { start: 21, end: 40 });

        let seed = b"ACGTACGTACGA";
        let positions = processor.seed_positions(seed, true);

        assert_eq!(positions, vec![21]);
        assert_eq!(processor.get_ref_base(25), Some(b'A'));
    }

    #[test]
    fn test_seed_positions_search_full_shared_reference_when_not_loaded() {
        let chrom = "testchr".to_string();
        let full_sequence = b"AAAAAAAAAAAAAAAAAAAAACGTACGTACGAACGTACGATTTTTTTTTTTTTTTTTTTT";

        let mut original_reference = Reference::new_with_start(full_sequence[..20].to_vec(), 1);
        original_reference.build_seed_map(20, Some(full_sequence.len()));

        let mut current_reference = Reference::new_with_start(full_sequence[40..].to_vec(), 41);
        current_reference.build_seed_map(full_sequence.len() as i64, Some(full_sequence.len()));

        let mut chromosomes: HashMap<String, ChromosomeData, LibDefaultHasher> = Default::default();
        chromosomes.insert(
            chrom.clone(),
            ChromosomeData {
                sequence: Arc::new(full_sequence.to_vec()),
                length: full_sequence.len(),
            },
        );

        let shared_reference = Arc::new(SharedReference {
            chromosomes,
            chromosome_names: vec![chrom.clone()],
            total_size: full_sequence.len(),
        });

        let mut processor = StructuralVariantsProcessor::new_with_context(
            Arc::clone(&original_reference.ref_seq),
            Arc::clone(&original_reference.seed),
            1,
            Some(chrom),
            Vec::new(),
            Some(shared_reference),
        );

        processor.reference_seq = Arc::clone(&current_reference.ref_seq);
        processor.reference_seed = Arc::clone(&current_reference.seed);
        processor.ref_start = 41;

        let seed = b"ACGTACGTACGA";
        let positions = processor.seed_positions(seed, true);

        // Java never scans the full chromosome for seeds — only bounded window
        // HashMap lookups.  Seed at position 21 is outside all loaded windows,
        // so it remains unfound (empty).
        assert_eq!(positions, vec![]);
    }

    #[test]
    fn test_seed_positions_keep_local_window_uniqueness() {
        let chrom = "testchr".to_string();
        let full_sequence = b"TTTTTTTTTTACGTACGTACGACCCCCCCCCCCCCCCCCCCCACGTACGTACGAGGGGGGGGGG";

        let mut original_reference = Reference::new_with_start(full_sequence[..35].to_vec(), 1);
        original_reference.build_seed_map(35, Some(full_sequence.len()));

        let mut chromosomes: HashMap<String, ChromosomeData, LibDefaultHasher> = Default::default();
        chromosomes.insert(
            chrom.clone(),
            ChromosomeData {
                sequence: Arc::new(full_sequence.to_vec()),
                length: full_sequence.len(),
            },
        );

        let shared_reference = Arc::new(SharedReference {
            chromosomes,
            chromosome_names: vec![chrom.clone()],
            total_size: full_sequence.len(),
        });

        let processor = StructuralVariantsProcessor::new_with_context(
            Arc::clone(&original_reference.ref_seq),
            Arc::clone(&original_reference.seed),
            1,
            Some(chrom),
            Vec::new(),
            Some(shared_reference),
        );

        let seed = b"ACGTACGTACGA";
        let positions = processor.seed_positions(seed, true);

        assert_eq!(positions, vec![11]);
    }

    #[test]
    fn test_seed_positions_search_full_shared_reference_after_restore_with_history() {
        let chrom = "testchr".to_string();
        let full_sequence = b"AAAAAAAAAAAAAAAAAAAAACGTACGTACGAACGTACGATTTTTTTTTTTTTTTTTTTT";

        let mut restored_reference = Reference::new_with_start(full_sequence[40..].to_vec(), 41);
        restored_reference.build_seed_map(full_sequence.len() as i64, Some(full_sequence.len()));

        let mut chromosomes: HashMap<String, ChromosomeData, LibDefaultHasher> = Default::default();
        chromosomes.insert(
            chrom.clone(),
            ChromosomeData {
                sequence: Arc::new(full_sequence.to_vec()),
                length: full_sequence.len(),
            },
        );

        let shared_reference = Arc::new(SharedReference {
            chromosomes,
            chromosome_names: vec![chrom.clone()],
            total_size: full_sequence.len(),
        });

        let mut processor = StructuralVariantsProcessor::new_with_context(
            Arc::clone(&restored_reference.ref_seq),
            Arc::clone(&restored_reference.seed),
            41,
            Some(chrom),
            Vec::new(),
            Some(shared_reference),
        );
        processor
            .historical_reference_windows
            .push(ReferenceWindow { start: 1, end: 20 });

        let seed = b"ACGTACGTACGA";
        let positions = processor.seed_positions(seed, true);

        // Java never scans the full chromosome for seeds.  Seed at position 21
        // is outside the historical window [1,20] and the current window [41+],
        // so it remains unfound.
        assert_eq!(positions, vec![]);
    }

    #[test]
    fn test_seed_positions_can_skip_shared_reference_fallback_for_plain_reverse_match() {
        let chrom = "testchr".to_string();
        let full_sequence = b"AAAAAAAAAAAAAAAAAAAAACGTACGTACGATTTTTTTTTTGGGGGGGGGGACGTACGTACGA";

        let mut original_reference = Reference::new_with_start(full_sequence[..20].to_vec(), 1);
        original_reference.build_seed_map(20, Some(full_sequence.len()));

        let mut current_reference = Reference::new_with_start(full_sequence[40..].to_vec(), 41);
        current_reference.build_seed_map(full_sequence.len() as i64, Some(full_sequence.len()));

        let mut chromosomes: HashMap<String, ChromosomeData, LibDefaultHasher> = Default::default();
        chromosomes.insert(
            chrom.clone(),
            ChromosomeData {
                sequence: Arc::new(full_sequence.to_vec()),
                length: full_sequence.len(),
            },
        );

        let shared_reference = Arc::new(SharedReference {
            chromosomes,
            chromosome_names: vec![chrom.clone()],
            total_size: full_sequence.len(),
        });

        let mut processor = StructuralVariantsProcessor::new_with_context(
            Arc::clone(&original_reference.ref_seq),
            Arc::clone(&original_reference.seed),
            1,
            Some(chrom),
            Vec::new(),
            Some(shared_reference),
        );

        processor.reference_seq = Arc::clone(&current_reference.ref_seq);
        processor.reference_seed = Arc::clone(&current_reference.seed);
        processor.ref_start = 41;
        processor
            .historical_reference_windows
            .push(ReferenceWindow { start: 21, end: 40 });

        let seed = b"ACGTACGTACGA";

        assert_eq!(
            processor.seed_positions_with_scope(seed, true, false),
            vec![21]
        );
        // With the shared-reference fallback disabled (matching Java), only
        // seeds found in loaded windows are returned — position 53 from the
        // full-chromosome scan is no longer produced.
        assert_eq!(
            processor.seed_positions_with_scope(seed, true, true),
            vec![21]
        );
    }

    #[test]
    fn test_find_match_searches_historical_windows() {
        let chrom = "testchr".to_string();
        let full_sequence = b"AAAAAAAAAAAAAAAAAAAAACGTACGTACGAACGTACGATTTTTTTTTTTTTTTTTTTT";

        let mut original_reference = Reference::new_with_start(full_sequence[..20].to_vec(), 1);
        original_reference.build_seed_map(20, Some(full_sequence.len()));

        let mut current_reference = Reference::new_with_start(full_sequence[40..].to_vec(), 41);
        current_reference.build_seed_map(full_sequence.len() as i64, Some(full_sequence.len()));

        let mut chromosomes: HashMap<String, ChromosomeData, LibDefaultHasher> = Default::default();
        chromosomes.insert(
            chrom.clone(),
            ChromosomeData {
                sequence: Arc::new(full_sequence.to_vec()),
                length: full_sequence.len(),
            },
        );

        let shared_reference = Arc::new(SharedReference {
            chromosomes,
            chromosome_names: vec![chrom.clone()],
            total_size: full_sequence.len(),
        });

        let mut processor = StructuralVariantsProcessor::new_with_context(
            Arc::clone(&original_reference.ref_seq),
            Arc::clone(&original_reference.seed),
            1,
            Some(chrom),
            Vec::new(),
            Some(shared_reference),
        );

        processor.reference_seq = Arc::clone(&current_reference.ref_seq);
        processor.reference_seed = Arc::clone(&current_reference.seed);
        processor.ref_start = 41;
        processor
            .historical_reference_windows
            .push(ReferenceWindow { start: 21, end: 40 });

        let match_result = processor.find_match(b"ACGTACGTACGA", 50, 1, 12, 0);

        assert_eq!(match_result.base_position, 21);
        assert!(match_result.matched_sequence.is_empty());
    }

    #[test]
    fn test_find_del_forward_nosoftp_uses_java_like_match_scope() {
        let _ = INSTANCE.get_or_init(GlobalReadOnlyScope::default);

        let chrom = "testchr".to_string();
        let seq = b"ACGTACGTACGA";
        let mut full_sequence = vec![b'A'; 1000];
        full_sequence[349..349 + seq.len()].copy_from_slice(seq);
        full_sequence[799..799 + seq.len()].copy_from_slice(seq);

        let mut original_reference = Reference::new_with_start(full_sequence[..100].to_vec(), 1);
        original_reference.build_seed_map(100, Some(full_sequence.len()));

        let mut chromosomes: HashMap<String, ChromosomeData, LibDefaultHasher> = Default::default();
        chromosomes.insert(
            chrom.clone(),
            ChromosomeData {
                sequence: Arc::new(full_sequence),
                length: 1000,
            },
        );

        let shared_reference = Arc::new(SharedReference {
            chromosomes,
            chromosome_names: vec![chrom.clone()],
            total_size: 1000,
        });

        let mut processor = StructuralVariantsProcessor::new_with_context(
            Arc::clone(&original_reference.ref_seq),
            Arc::clone(&original_reference.seed),
            1,
            Some(chrom),
            Vec::new(),
            Some(shared_reference),
        );

        let mut data = RealignedVariationData {
            max_read_length: 100,
            ..Default::default()
        };

        let mut del = SoftClip::default();
        del.end = 300;
        del.mstart = 340;
        del.mend = 360;
        del.var.alt_depth = 10;
        data.svfdel.push(del);

        let mut scv = SoftClip::default();
        scv.var.alt_depth = 1;
        scv.var.alt_depth_rev = 1;
        scv.set_consensus_seq(seq.to_vec());
        data.soft_clips_3end.insert(300, scv);

        processor.find_del(&mut data);

        let emitted = data
            .non_insertion_variants
            .get(&300)
            .and_then(|pos_map| pos_map.iter().find(|(key, _)| key.to_string() == "-50"))
            .map(|(_, variant)| variant.alt_depth);

        assert_eq!(emitted, Some(11));
        assert_eq!(data.sv_counts.get(&300).map(|sv| sv.pairs), Some(10));
        assert_eq!(data.sv_counts.get(&300).map(|sv| sv.splits), Some(1));
        assert!(data.svfdel[0].used());
    }

    #[test]
    fn test_find_match_rev_findsv_can_use_shared_reference_fallback() {
        let chrom = "testchr".to_string();
        let target = b"ACGTACGTACGA";
        let mut soft_seq = target.to_vec();
        soft_seq.reverse();
        for base in &mut soft_seq {
            *base = StructuralVariantsProcessor::complement_base_u8(*base);
        }

        let mut full_sequence = vec![b'A'; 120];
        full_sequence[52..52 + target.len()].copy_from_slice(target);

        let mut original_reference = Reference::new_with_start(full_sequence[..20].to_vec(), 1);
        original_reference.build_seed_map(20, Some(full_sequence.len()));

        let mut current_reference = Reference::new_with_start(full_sequence[40..80].to_vec(), 41);
        current_reference.build_seed_map(full_sequence.len() as i64, Some(full_sequence.len()));

        let mut chromosomes: HashMap<String, ChromosomeData, LibDefaultHasher> = Default::default();
        chromosomes.insert(
            chrom.clone(),
            ChromosomeData {
                sequence: Arc::new(full_sequence),
                length: 120,
            },
        );

        let shared_reference = Arc::new(SharedReference {
            chromosomes,
            chromosome_names: vec![chrom.clone()],
            total_size: 120,
        });

        let mut processor = StructuralVariantsProcessor::new_with_context(
            Arc::clone(&original_reference.ref_seq),
            Arc::clone(&original_reference.seed),
            1,
            Some(chrom),
            Vec::new(),
            Some(shared_reference),
        );

        processor.reference_seq = Arc::clone(&current_reference.ref_seq);
        processor.reference_seed = Arc::clone(&current_reference.seed);
        processor.ref_start = 41;

        assert_eq!(
            processor
                .find_match_rev(&soft_seq, 45, 1, Configuration::SEED_1 as usize, 3)
                .base_position,
            0
        );
        assert_eq!(processor.find_match_rev_findsv(&soft_seq, 45, 1).base_position, 64);
    }

    #[test]
    fn test_reference_coverage_tracking_is_separate_from_historical_reference_windows() {
        let mut processor =
            StructuralVariantsProcessor::new(b"ACGTACGT".to_vec(), Default::default(), 100);

        processor
            .historical_reference_windows
            .push(ReferenceWindow { start: 21, end: 40 });

        assert!(processor.is_span_loaded(21, 40));
        assert!(!processor.is_reference_coverage_loaded(21, 40));

        processor.remember_reference_coverage_window(21, 40);

        assert!(processor.is_reference_coverage_loaded(21, 40));
    }

    #[test]
    fn test_reference_coverage_windows_merge_overlaps() {
        let mut processor =
            StructuralVariantsProcessor::new(b"ACGTACGT".to_vec(), Default::default(), 100);

        processor.remember_reference_coverage_window(21, 40);
        processor.remember_reference_coverage_window(35, 55);

        assert_eq!(
            processor.loaded_reference_coverage_windows,
            vec![ReferenceWindow { start: 21, end: 55 }]
        );
        assert!(processor.is_reference_coverage_loaded(30, 50));
    }

    #[test]
    fn test_is_region_loaded_checks_original_and_loaded_regions() {
        let mut processor =
            StructuralVariantsProcessor::new(b"ACGTACGT".to_vec(), Default::default(), 100);
        processor.loaded_regions.push((210, 229));

        assert!(processor.is_region_loaded(100, 107));
        assert!(processor.is_region_loaded(210, 229));
        assert!(!processor.is_region_loaded(208, 229));
    }

    #[test]
    fn test_record_del_rightseq_from_loaded_regions_marks_loaded_region_del() {
        let mut processor =
            StructuralVariantsProcessor::new(b"ACGTACGT".to_vec(), Default::default(), 100);
        processor.loaded_regions.push((210, 229));

        let mut data = RealignedVariationData::default();
        let variation = StructuralVariantsProcessor::get_or_create_variation(
            &mut data.non_insertion_variants,
            200,
            "-10",
        );
        variation.alt_depth = 1;

        processor.record_del_rightseq_from_loaded_regions(&data);

        assert!(
            processor
                .historical_del_rightseq_variants
                .get(&200)
                .is_some_and(|descriptions| descriptions.contains("-10"))
        );
    }

    #[test]
    fn test_record_del_rightseq_from_loaded_regions_ignores_unloaded_del() {
        let mut processor =
            StructuralVariantsProcessor::new(b"ACGTACGT".to_vec(), Default::default(), 100);

        let mut data = RealignedVariationData::default();
        let variation = StructuralVariantsProcessor::get_or_create_variation(
            &mut data.non_insertion_variants,
            200,
            "-10",
        );
        variation.alt_depth = 1;

        processor.record_del_rightseq_from_loaded_regions(&data);

        assert!(
            processor
                .historical_del_rightseq_variants
                .get(&200)
                .is_none()
        );
    }

    #[test]
    fn test_record_del_rightseq_from_loaded_regions_marks_complex_del_key() {
        let mut processor =
            StructuralVariantsProcessor::new(b"ACGTACGT".to_vec(), Default::default(), 100);
        processor.loaded_regions.push((1925, 1944));

        let mut data = RealignedVariationData::default();
        let variation = StructuralVariantsProcessor::get_or_create_variation(
            &mut data.non_insertion_variants,
            200,
            "-1725#TTGTGAAGATATTT^-1713&AGGCCTATTTAGG",
        );
        variation.alt_depth = 1;

        processor.record_del_rightseq_from_loaded_regions(&data);

        assert!(
            processor
                .historical_del_rightseq_variants
                .get(&200)
                .is_some_and(|descriptions| {
                    descriptions.contains("-1725#TTGTGAAGATATTT^-1713&AGGCCTATTTAGG")
                })
        );
    }

    #[test]
    fn test_uncovered_reference_coverage_spans_skip_loaded_segments() {
        let mut processor =
            StructuralVariantsProcessor::new(b"ACGTACGT".to_vec(), Default::default(), 100);

        processor.remember_reference_coverage_window(21, 40);
        processor.remember_reference_coverage_window(60, 80);

        assert_eq!(
            processor.uncovered_reference_coverage_spans(10, 90),
            vec![
                ReferenceWindow { start: 10, end: 20 },
                ReferenceWindow { start: 41, end: 59 },
                ReferenceWindow { start: 81, end: 90 },
            ]
        );
    }

    #[test]
    fn test_adj_cnt() {
        let mut variant = Variant::default();

        adj_cnt(
            &mut variant,
            10,   // vars_count
            8,    // high_qual_cnt
            2,    // low_qual_cnt
            50.0, // mean_pos
            30.0, // mean_qual
            60.0, // mean_mapq
            1.0,  // nm
            6,    // fwd_cnt
            4,    // rev_cnt
        );

        assert_eq!(variant.alt_depth, 10);
        assert_eq!(variant.high_qual_read_cnt, 8);
        assert_eq!(variant.low_qual_read_cnt, 2);
        assert_eq!(variant.mean_pos, 50.0);
        assert_eq!(variant.alt_depth_fwd, 6);
        assert_eq!(variant.alt_depth_rev, 4);
        assert!(variant.pstd);
        assert!(variant.qstd);
    }

    #[test]
    fn test_empty_data_processing() {
        // This test only verifies the adj_snv logic with empty data
        // The process() method requires GlobalReadOnlyScope initialization
        let processor = StructuralVariantsProcessor::new(b"ACGT".to_vec(), Default::default(), 1);

        let mut data = RealignedVariationData::default();

        // Directly test adj_snv with empty data
        processor.adj_snv(&mut data);

        // Should return empty data without errors
        assert!(data.non_insertion_variants.is_empty());
        assert!(data.insertion_variants.is_empty());
    }
}
