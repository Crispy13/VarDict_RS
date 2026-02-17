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

use std::collections::HashMap;
use crackle_kit::tracing::{event, Level};

use crate::conf::Configuration;
use crate::scopedata::global_read_only_scope::instance;
use crate::variants::variants::{VarDesc, Variant, SoftClip};

/// Input data for StructuralVariantsProcessor (from VariantRealigner)
#[derive(Default)]
pub struct RealignedVariationData {
    /// Non-insertion variants by position
    pub non_insertion_variants: HashMap<i64, HashMap<VarDesc, Variant>>,
    /// Insertion variants by position
    pub insertion_variants: HashMap<i64, HashMap<VarDesc, Variant>>,
    /// 5' end soft clips by position
    pub soft_clips_5end: HashMap<i64, SoftClip>,
    /// 3' end soft clips by position  
    pub soft_clips_3end: HashMap<i64, SoftClip>,
    /// Reference coverage by position
    pub ref_coverage: HashMap<i64, usize>,
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
}

/// Output data from StructuralVariantsProcessor (same structure, possibly modified)
pub type ProcessedVariationData = RealignedVariationData;

/// StructuralVariantsProcessor - Processes variation data for SVs and adjusts SNV counts
pub struct StructuralVariantsProcessor {
    /// Reference sequence
    reference_seq: Vec<u8>,
    /// Reference seed map
    reference_seed: HashMap<Vec<u8>, Vec<i64>>,
    /// Reference start position (1-based)
    ref_start: i64,
}

impl StructuralVariantsProcessor {
    /// Create a new StructuralVariantsProcessor
    pub fn new(
        reference_seq: Vec<u8>,
        reference_seed: HashMap<Vec<u8>, Vec<i64>>,
        ref_start: i64,
    ) -> Self {
        StructuralVariantsProcessor {
            reference_seq,
            reference_seed,
            ref_start,
        }
    }

    /// Process realigned variation data
    /// 
    /// This is equivalent to Java StructuralVariantsProcessor.process()
    pub fn process(&self, mut data: RealignedVariationData) -> ProcessedVariationData {
        // If SV is enabled, find structural variants
        if !instance().conf.disable_sv {
            self.find_all_svs(&mut data);
        }

        // Always adjust SNV counts from soft clips
        self.adj_snv(&mut data);

        data
    }

    /// Find all structural variants (DEL, INV, DUP)
    /// 
    /// Called when SV detection is enabled
    fn find_all_svs(&self, data: &mut RealignedVariationData) {
        self.find_del(data);
        self.find_svs_del_candidates(data);
        self.find_del_disc(data);
        self.find_dup_disc(data);
    }

    fn find_del(&self, data: &mut RealignedVariationData) {
        let minr = instance().conf.minr;

        for idx in 0..data.svfdel.len() {
            let (used, vars_count, end, mstart, mean_qual, mean_pos, mean_mapq, nm, soft_map) = {
                let del = &data.svfdel[idx];
                (
                    del.used(),
                    del.var.alt_depth,
                    del.end,
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

                let m = self.find_match(&seq, softp, 1, Configuration::SEED_1 as usize, 3);
                if m.base_position == 0 {
                    event!(
                        Level::DEBUG,
                        phase = "find_del_forward_nomatch",
                        idx,
                        softp,
                    );
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
                    let variation =
                        Self::get_or_create_variation(&mut data.non_insertion_variants, p5, &del_key);
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
                let mut candidate_positions: Vec<i64> = data.soft_clips_3end.keys().copied().collect();
                candidate_positions.sort_unstable();

                for candidate in candidate_positions {
                    if !(candidate >= end - 3
                        && candidate - end < 3 * data.max_read_length as i64)
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

                    let mut m = self.find_match(&seq, candidate, 1, Configuration::SEED_1 as usize, 3);
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
                    let vref =
                        Self::get_or_create_variation(&mut data.non_insertion_variants, candidate, &del_key);
                    vref.alt_depth = 0;

                    let split_count = data
                        .soft_clips_3end
                        .get(&candidate)
                        .map(|s| s.var.alt_depth)
                        .unwrap_or(0);
                    Self::add_sv_counts(
                        &mut data.non_insertion_variants,
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
            let (used, vars_count, mend, start, mean_qual, mean_pos, mean_mapq, nm, soft_map) = {
                let del = &data.svrdel[idx];
                (
                    del.used(),
                    del.var.alt_depth,
                    del.mend,
                    del.start,
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

                let mut m =
                    self.find_match(&seq, softp, -1, Configuration::SEED_1 as usize, 3);
                if m.base_position == 0 {
                    m = self.find_match(&seq, softp, -1, Configuration::SEED_2 as usize, 0);
                }
                if m.base_position == 0 {
                    event!(
                        Level::DEBUG,
                        phase = "find_del_reverse_nomatch",
                        idx,
                        softp,
                    );
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
                    bp,
                    vars_count,
                    split_count,
                    1,
                );

                if let Some(scv) = data.soft_clips_5end.get(&softp) {
                    let variation =
                        Self::get_or_create_variation(&mut data.non_insertion_variants, bp, &del_key);
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
                let mut candidate_positions: Vec<i64> = data.soft_clips_5end.keys().copied().collect();
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
                    let vref =
                        Self::get_or_create_variation(&mut data.non_insertion_variants, bp, &del_key);
                    vref.alt_depth = 0;

                    let split_count = data
                        .soft_clips_5end
                        .get(&candidate)
                        .map(|s| s.var.alt_depth)
                        .unwrap_or(0);
                    Self::add_sv_counts(
                        &mut data.non_insertion_variants,
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

    fn find_svs_del_candidates(&self, data: &mut RealignedVariationData) {
        let minr = instance().conf.minr;

        let mut tmp5: Vec<SortPositionSoftClip> = data
            .soft_clips_5end
            .iter()
            .filter(|(_, sclip)| !sclip.used())
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

            if data
                .soft_clips_5end
                .get(&p5)
                .map(|s| s.used())
                .unwrap_or(true)
            {
                continue;
            }

            let seq = {
                let Some(sc5v) = data.soft_clips_5end.get_mut(&p5) else {
                    continue;
                };
                self.find_conseq(sc5v)
            };
            if seq.is_empty() || seq.len() < Configuration::SEED_2 as usize {
                continue;
            }

            let m = self.find_match(&seq, p5, -1, Configuration::SEED_1 as usize, 3);
            let bp = m.base_position;
            event!(
                Level::DEBUG,
                phase = "findsv_5_candidate",
                p5,
                cnt5,
                bp,
            );
            if bp != 0 && bp < p5 {
                let pairs_data = Self::check_pairs(
                    bp,
                    p5,
                    &mut data.svfdel,
                    &mut data.svrdel,
                    data.max_read_length as i64,
                );
                if pairs_data.pairs == 0 {
                    event!(
                        Level::DEBUG,
                        phase = "findsv_5_pairs_zero",
                        p5,
                        bp,
                    );
                    continue;
                }

                let p5_adj = p5 - 1;
                let bp_adj = bp + 1;
                let dellen = p5_adj - bp_adj + 1;
                if dellen <= 0 {
                    continue;
                }

                let del_key = format!("-{}", dellen);
                let vref =
                    Self::get_or_create_variation(&mut data.non_insertion_variants, bp_adj, &del_key);
                vref.alt_depth = 0;

                Self::add_sv_counts(
                    &mut data.non_insertion_variants,
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
                    let variation =
                        Self::get_or_create_variation(&mut data.non_insertion_variants, bp_adj, &del_key);
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

                let variation =
                    Self::get_or_create_variation(&mut data.non_insertion_variants, bp_adj, &del_key);
                adj_cnt_from_variant(variation, &tmp);
            }
        }

        let mut tmp3: Vec<SortPositionSoftClip> = data
            .soft_clips_3end
            .iter()
            .filter(|(_, sclip)| !sclip.used())
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

            if data
                .soft_clips_3end
                .get(&p3)
                .map(|s| s.used())
                .unwrap_or(true)
            {
                continue;
            }

            let seq = {
                let Some(sc3v) = data.soft_clips_3end.get_mut(&p3) else {
                    continue;
                };
                self.find_conseq(sc3v)
            };
            if seq.is_empty() || seq.len() < Configuration::SEED_2 as usize {
                continue;
            }

            let m = self.find_match(&seq, p3, 1, Configuration::SEED_1 as usize, 3);
            let mut bp = m.base_position;
            event!(
                Level::DEBUG,
                phase = "findsv_3_candidate",
                p3,
                cnt3,
                bp,
            );
            if bp != 0 && bp > p3 {
                let pairs_data = Self::check_pairs(
                    p3,
                    bp,
                    &mut data.svfdel,
                    &mut data.svrdel,
                    data.max_read_length as i64,
                );
                if pairs_data.pairs == 0 {
                    event!(
                        Level::DEBUG,
                        phase = "findsv_3_pairs_zero",
                        p3,
                        bp,
                    );
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
                let vref =
                    Self::get_or_create_variation(&mut data.non_insertion_variants, p3, &del_key);
                vref.alt_depth = 0;

                Self::add_sv_counts(
                    &mut data.non_insertion_variants,
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
                    let variation =
                        Self::get_or_create_variation(&mut data.non_insertion_variants, p3, &del_key);
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

                let variation =
                    Self::get_or_create_variation(&mut data.non_insertion_variants, p3, &del_key);
                adj_cnt_from_variant(variation, &tmp);
            }
        }
    }

    fn find_del_disc(&self, data: &mut RealignedVariationData) {
        let min_dist = 8 * data.max_read_length as i64;
        let minr = instance().conf.minr;

        for idx in 0..data.svfdel.len() {
            let (used, vars_count, end, mstart, mean_mapq, mean_qual, mean_pos, nm, softp) = {
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
                )
            };

            if used || vars_count < minr + 5 {
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

            let del_key = format!("-{}", mlen);
            let vref = Self::get_or_create_variation(&mut data.non_insertion_variants, bp, &del_key);
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
            Self::add_sv_counts(&mut data.non_insertion_variants, bp, vars_count, splits, 1);

            let mut tv = Variant::default();
            tv.alt_depth = 2 * vars_count;
            tv.high_qual_read_cnt = 2 * vars_count;
            tv.alt_depth_fwd = vars_count;
            tv.alt_depth_rev = vars_count;
            tv.mean_qual = 2.0 * mean_qual;
            tv.mean_pos = 2.0 * mean_pos;
            tv.mean_mapq = 2.0 * mean_mapq;
            tv.nm = 2.0 * nm;
            let variation = Self::get_or_create_variation(&mut data.non_insertion_variants, bp, &del_key);
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
            let (used, vars_count, start, mend, mean_mapq, mean_qual, mean_pos, nm, softp) = {
                let del = &data.svrdel[idx];
                (
                    del.used(),
                    del.var.alt_depth,
                    del.start,
                    del.mend,
                    del.var.mean_mapq,
                    del.var.mean_qual,
                    del.var.mean_pos,
                    del.var.nm,
                    del.softp,
                )
            };

            if used || vars_count < minr + 5 {
                continue;
            }
            if start <= mend + min_dist {
                continue;
            }
            if vars_count == 0 || mean_mapq / vars_count as f64 <= Configuration::DISCPAIRQUAL {
                continue;
            }

            let mlen = start - mend - data.max_read_length as i64 / (vars_count + 1) as i64;
            if !(mlen > 0 && mlen > min_dist) {
                continue;
            }

            let bp = mend + (data.max_read_length as i64 / (vars_count + 1) as i64) / 2;

            let del_key = format!("-{}", mlen);
            let vref = Self::get_or_create_variation(&mut data.non_insertion_variants, bp, &del_key);
            vref.alt_depth = 0;

            let splits = data
                .soft_clips_3end
                .get(&(mend + 1))
                .map(|s| s.var.alt_depth)
                .unwrap_or(0)
                + data
                    .soft_clips_5end
                    .get(&start)
                    .map(|s| s.var.alt_depth)
                    .unwrap_or(0);
            Self::add_sv_counts(&mut data.non_insertion_variants, bp, vars_count, splits, 1);

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
            let variation = Self::get_or_create_variation(&mut data.non_insertion_variants, bp, &del_key);
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
            Self::mark_sv(mend, start, &mut data.svfdel, data.max_read_length as i64);
        }
    }

    fn find_dup_disc(&self, data: &mut RealignedVariationData) {
        let minr = instance().conf.minr;

        for idx in 0..data.svfdup.len() {
            let (
                used,
                ms,
                _me,
                cnt,
                mut end,
                _start,
                pmean,
                qmean,
                q_mean,
                nm,
                softp,
                soft_map,
            ) = {
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

            let mut ins = self.join_ref(bp, bp + Configuration::SVFLANK as i64 - 1);
            let dup_len = mlen - 2 * Configuration::SVFLANK as i64;
            ins.extend_from_slice(format!("<dup{}>", dup_len).as_bytes());
            ins.extend_from_slice(&self.join_ref(
                pe - Configuration::SVFLANK as i64 + 1,
                pe,
            ));

            let splits = if softp != 0 {
                data.soft_clips_3end
                    .get(&(softp as i64))
                    .map(|s| s.var.alt_depth)
                    .unwrap_or(0)
            } else {
                0
            };
            Self::add_sv_counts(&mut data.non_insertion_variants, bp, cnt, splits, 1);

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
                seq: ins.iter().map(|b| b.to_ascii_uppercase()).collect(),
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

            let (clusters, _) = Self::mark_dup_sv(bp, pe, &mut data.svrdup, data.max_read_length as i64);
            if clusters != 0 {
                Self::add_sv_counts(&mut data.non_insertion_variants, bp, 0, 0, clusters);
            }
        }

        for idx in 0..data.svrdup.len() {
            let (
                used,
                _ms,
                me,
                cnt,
                _end,
                start,
                pmean,
                qmean,
                q_mean,
                nm,
                soft_map,
            ) = {
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
                        && self.get_ref_base(tpe)
                            == self.get_ref_base(bp + (tpe - pe - 1))
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

            let mut ins = self.join_ref(bp, bp + Configuration::SVFLANK as i64 - 1);
            let dup_len = mlen - 2 * Configuration::SVFLANK as i64;
            ins.extend_from_slice(format!("<dup{}>", dup_len).as_bytes());
            ins.extend_from_slice(&self.join_ref(
                pe - Configuration::SVFLANK as i64 + 1,
                pe,
            ));

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
            Self::add_sv_counts(&mut data.non_insertion_variants, bp, cnt, splits, 1);

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
                seq: ins.iter().map(|b| b.to_ascii_uppercase()).collect(),
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

            let (clusters, _) = Self::mark_dup_sv(bp, pe, &mut data.svfdup, data.max_read_length as i64);
            if clusters != 0 {
                Self::add_sv_counts(&mut data.non_insertion_variants, bp, 0, 0, clusters);
            }
        }
    }

    fn inc_ref_coverage(ref_coverage: &mut HashMap<i64, usize>, pos: i64, cnt: usize) {
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

    fn select_primary_soft_pos(soft: &HashMap<i64, usize>) -> Option<i64> {
        soft.iter().max_by_key(|(_, count)| *count).map(|(pos, _)| *pos)
    }

    fn get_or_create_variation<'a>(
        map: &'a mut HashMap<i64, HashMap<VarDesc, Variant>>,
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

    fn add_sv_counts(
        map: &mut HashMap<i64, HashMap<VarDesc, Variant>>,
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
        let sv = pos_map.entry(key).or_default();
        sv.alt_depth += pairs;
        sv.high_qual_read_cnt += splits;
        sv.low_qual_read_cnt += clusters;
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

        let mut out = Vec::with_capacity((end - start + 1) as usize);
        for pos in start..=end {
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
        _position: i64,
        dir: i64,
        seed_len: usize,
        mm: usize,
    ) -> MatchResult {
        let mut seq_work = seq.to_vec();
        if dir == -1 {
            seq_work.reverse();
        }

        if seq_work.len() < seed_len {
            return MatchResult::default();
        }

        for i in (0..=seq_work.len() - seed_len).rev() {
            let seed = &seq_work[i..i + seed_len];
            let Some(seeds) = self.reference_seed.get(seed) else {
                continue;
            };
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
                let mut extra: Vec<u8> = Vec::new();
                let mut mm_idx: i64 = if dir == -1 { -1 } else { 0 };
                loop {
                    let Some(ch) = Self::char_at(&seq_work, mm_idx) else {
                        break;
                    };
                    if self.is_has_and_not_equals(bp, ch) {
                        extra.push(ch);
                        bp += dir;
                        mm_idx += dir;
                    } else {
                        break;
                    }
                }
                if !extra.is_empty() && dir == -1 {
                    extra.reverse();
                }
                return MatchResult {
                    base_position: bp,
                    matched_sequence: extra,
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

                    // Get the soft clip variant counts
                    let sclip = data.soft_clips_5end.get(&position).unwrap();
                    let vars_count = sclip.var.alt_depth;
                    let high_qual_cnt = sclip.var.high_qual_read_cnt;
                    let low_qual_cnt = sclip.var.low_qual_read_cnt;
                    let mean_pos = sclip.var.mean_pos;
                    let mean_qual = sclip.var.mean_qual;
                    let mean_mapq = sclip.var.mean_mapq;
                    let nm = sclip.var.nm;
                    let fwd_cnt = sclip.var.alt_depth_fwd;
                    let rev_cnt = sclip.var.alt_depth_rev;

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

                    // Get the soft clip variant counts
                    let sclip = data.soft_clips_3end.get(&position).unwrap();
                    let vars_count = sclip.var.alt_depth;
                    let high_qual_cnt = sclip.var.high_qual_read_cnt;
                    let low_qual_cnt = sclip.var.low_qual_read_cnt;
                    let mean_pos = sclip.var.mean_pos;
                    let mean_qual = sclip.var.mean_qual;
                    let mean_mapq = sclip.var.mean_mapq;
                    let nm = sclip.var.nm;
                    let fwd_cnt = sclip.var.alt_depth_fwd;
                    let rev_cnt = sclip.var.alt_depth_rev;

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

    /// Get reference base at a position (1-based)
    fn get_ref_base(&self, pos: i64) -> Option<u8> {
        if pos < self.ref_start {
            return None;
        }
        let idx = (pos - self.ref_start) as usize;
        self.reference_seq.get(idx).copied()
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_processor_creation() {
        let processor = StructuralVariantsProcessor::new(
            b"ACGTACGT".to_vec(),
            HashMap::new(),
            100,
        );
        assert_eq!(processor.ref_start, 100);
    }

    #[test]
    fn test_get_ref_base() {
        let processor = StructuralVariantsProcessor::new(
            b"ACGTACGT".to_vec(),
            HashMap::new(),
            100,
        );
        
        assert_eq!(processor.get_ref_base(100), Some(b'A'));
        assert_eq!(processor.get_ref_base(101), Some(b'C'));
        assert_eq!(processor.get_ref_base(107), Some(b'T'));
        assert_eq!(processor.get_ref_base(108), None); // Out of bounds
        assert_eq!(processor.get_ref_base(99), None);  // Before start
    }

    #[test]
    fn test_adj_cnt() {
        let mut variant = Variant::default();
        
        adj_cnt(
            &mut variant,
            10,  // vars_count
            8,   // high_qual_cnt
            2,   // low_qual_cnt
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
        let processor = StructuralVariantsProcessor::new(
            b"ACGT".to_vec(),
            HashMap::new(),
            1,
        );
        
        let mut data = RealignedVariationData::default();
        
        // Directly test adj_snv with empty data
        processor.adj_snv(&mut data);
        
        // Should return empty data without errors
        assert!(data.non_insertion_variants.is_empty());
        assert!(data.insertion_variants.is_empty());
    }
}
