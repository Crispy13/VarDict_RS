use rust_htslib::bam::{Record, record::Cigar};
use std::collections::HashMap;
use std::sync::Arc;

use crate::{
    conf::Configuration,
    data::bam_reader::BamReader,
    data::patterns::{
        AMP_ATGC, ATGSs_AMP_ATGSs_END, BEGIN_PLUS_ATGC, CARET_ATGC_END, CARET_ATGNC, DUP_NUM_ATGC,
        HASH_ATGC, UP_NUMBER_END,
    },
    data::reference::{Reference, ReferenceSeedMap},
    data::region::Region,
    mods::structural_variants_processor::RealignedVariationData,
    mods::vardict_pipeline::VarDictPipeline,
    prelude::{LibDefaultHasher, SmallVecBytes},
    utils::vec_map::{Entry as VecEntry, VecMap},
    variants::{
        var_utils::{find_conseq, get_variants_from_map},
        variants::{InsOrDelLen, SoftClip, VarDesc, Variant},
    },
};

type VariantMap = VecMap<VarDesc, Variant>;
type VariantMapByPos = HashMap<i64, VariantMap, LibDefaultHasher>;
type CountMap = VecMap<String, usize>;
type CountMapByPos = HashMap<i64, CountMap, LibDefaultHasher>;

/// Result of finding 3'/5' end matches between two sequences
#[derive(Debug, Clone)]
pub struct Match35 {
    pub matched_5_end: usize,
    pub matched_3_end: usize,
    pub max_matched_length: usize,
}

#[derive(Debug, Clone)]
struct Match {
    base_position: i64,
    matched_sequence: Vec<u8>,
}

#[derive(Debug, Clone)]
struct BaseInsertion {
    base_insert: i64,
    insertion_sequence: Vec<u8>,
    base_insert2: i64,
}

#[derive(Debug, Clone)]
struct Mismatch {
    mismatch_sequence: String,
    mismatch_position: i64,
    end: u8, // 3 or 5
}

#[derive(Default, Debug, Clone)]
struct MismatchResult {
    mismatches: Vec<Mismatch>,
    scp: Vec<i64>,
    nm: usize,
    misp: i64,
    misnt: Option<u8>,
}

#[derive(Default, Debug, Clone)]
struct SvCluster {
    mate_start: i64,
    mate_end: i64,
    cnt: usize,
    mate_len: i32,
    start: i64,
    end: i64,
    mean_pos: f64,
    mean_qual: f64,
    mean_mapq: f64,
    nm: f64,
}

fn find_conseq_transient(softclip: &mut SoftClip, dir: i32) -> Vec<u8> {
    find_conseq(softclip, dir)
}

fn cache_transient_conseq(softclip: &mut SoftClip, seq: &[u8]) {
    if !softclip.consensus_seq_is_set() {
        softclip.set_consensus_seq(seq.to_vec());
    }
}

fn normalize_insertion_key_bytes(bytes: &[u8]) -> SmallVecBytes {
    let mut normalized = SmallVecBytes::with_capacity(bytes.len());
    let mut inside_marker = false;
    for &byte in bytes {
        match byte {
            b'<' => {
                inside_marker = true;
                normalized.push(byte);
            }
            b'>' => {
                inside_marker = false;
                normalized.push(byte);
            }
            _ if inside_marker => normalized.push(byte),
            _ => normalized.push(byte.to_ascii_uppercase()),
        }
    }
    normalized
}

pub struct VariantRealigner {
    reference_seq: Arc<Vec<u8>>,
    reference_seed: Arc<ReferenceSeedMap>,
    ref_start: i64,
    original_reference_seq: Arc<Vec<u8>>,
    original_ref_start: i64,
    chromosome: Option<String>,
    bam_paths: Vec<String>,
}

impl VariantRealigner {
    pub fn new(reference_seq: Vec<u8>, reference_seed: ReferenceSeedMap, ref_start: i64) -> Self {
        Self::new_with_context(
            Arc::new(reference_seq),
            Arc::new(reference_seed),
            ref_start,
            None,
            Vec::new(),
        )
    }

    pub fn new_with_context(
        reference_seq: Arc<Vec<u8>>,
        reference_seed: Arc<ReferenceSeedMap>,
        ref_start: i64,
        chromosome: Option<String>,
        bam_paths: Vec<String>,
    ) -> Self {
        Self {
            original_reference_seq: Arc::clone(&reference_seq),
            original_ref_start: ref_start,
            reference_seq,
            reference_seed,
            ref_start,
            chromosome,
            bam_paths,
        }
    }

    pub fn with_reference_fallback(
        mut self,
        original_reference_seq: Arc<Vec<u8>>,
        original_ref_start: i64,
    ) -> Self {
        self.original_reference_seq = original_reference_seq;
        self.original_ref_start = original_ref_start;
        self
    }

    pub fn process_deletions(
        &self,
        data: &mut RealignedVariationData,
        position_to_deletions_count: &CountMapByPos,
    ) {
        self.process_deletions_with_passing_check(data, position_to_deletions_count, true);
    }

    fn process_deletions_with_passing_check(
        &self,
        data: &mut RealignedVariationData,
        position_to_deletions_count: &CountMapByPos,
        enable_passing_reads_check: bool,
    ) {
        let positions: Vec<i64> = data.non_insertion_variants.keys().copied().collect();
        for position in positions {
            if let Some(pos_map) = data.non_insertion_variants.get_mut(&position) {
                Self::merge_duplicate_deletion_keys(pos_map);
            }
        }

        let mut del_keys: Vec<(i64, VarDesc, String, usize)> = Vec::new();

        if !position_to_deletions_count.is_empty() {
            for (pos, desc_map) in position_to_deletions_count {
                for (desc_str, count) in desc_map {
                    if *count == 0 {
                        continue;
                    }

                    let desc_key = data
                        .non_insertion_variants
                        .get(pos)
                        .and_then(|var_map| Self::find_desc_by_key_string(var_map, desc_str))
                        .cloned()
                        .unwrap_or_else(|| VarDesc::Raw {
                            desc: desc_str.as_bytes().to_vec().into(),
                        });

                    del_keys.push((*pos, desc_key, desc_str.clone(), *count));
                }
            }
        } else {
            for (pos, var_map) in data.non_insertion_variants.iter() {
                for (desc, variant) in var_map {
                    let desc_str = desc.to_key_string();
                    if desc_str.starts_with('-') {
                        del_keys.push((*pos, desc.clone(), desc_str, variant.alt_depth));
                    }
                }
            }
        }

        del_keys.sort_by(|a, b| {
            b.3.cmp(&a.3)
                .then_with(|| a.0.cmp(&b.0))
                .then_with(|| b.2.cmp(&a.2))
        });

        if !position_to_deletions_count.is_empty() {
            for (pos, desc, _, _) in &del_keys {
                data.non_insertion_variants
                    .entry(*pos)
                    .or_default()
                    .entry(desc.clone())
                    .or_default();
            }
        }

        for (pos, desc, _, count) in &del_keys {
            self.realign_deletion_mismatches(*pos, desc, *count, data, enable_passing_reads_check);
        }

        // Java realigndel post-pass:
        // for i = tmp.size()-1; i > 0; i-- {
        //   if vn =~ /^(-\d+)&[ATGC]+$/ and vars(vn) < vars($1) then merge vn into $1
        // }
        for idx in (1..del_keys.len()).rev() {
            let (pos, desc, desc_str, _) = &del_keys[idx];
            let Some(base_del_desc) = Self::minus_amp_base_desc(desc_str) else {
                continue;
            };

            let Some(pos_map) = data.non_insertion_variants.get_mut(pos) else {
                continue;
            };

            let Some(base_key) = pos_map
                .keys()
                .find(|key| key.to_key_string() == base_del_desc)
                .cloned()
            else {
                continue;
            };

            let Some(vref_alt_depth) = pos_map.get(desc).map(|v| v.alt_depth) else {
                continue;
            };
            let Some(tref_alt_depth) = pos_map.get(&base_key).map(|v| v.alt_depth) else {
                continue;
            };

            if vref_alt_depth < tref_alt_depth {
                if let Some(vref) = pos_map.remove(desc) {
                    if let Some(tref) = pos_map.get_mut(&base_key) {
                        adj_cnt(tref, &vref);
                    } else {
                        pos_map.insert(desc.clone(), vref);
                    }
                }
            }
        }
    }

    fn find_desc_by_key_string<'a>(
        var_map: &'a VariantMap,
        key_string: &str,
    ) -> Option<&'a VarDesc> {
        var_map
            .keys()
            .find(|desc| desc.to_key_string() == key_string)
    }

    pub fn filter_all_sv_structures(&self, data: &mut RealignedVariationData) {
        Self::filter_sv(&mut data.svfinv3, data.max_read_length);
        Self::filter_sv(&mut data.svrinv3, data.max_read_length);
        Self::filter_sv(&mut data.svfinv5, data.max_read_length);
        Self::filter_sv(&mut data.svrinv5, data.max_read_length);
        Self::filter_sv(&mut data.svfdel, data.max_read_length);
        Self::filter_sv(&mut data.svrdel, data.max_read_length);
        Self::filter_sv(&mut data.svfdup, data.max_read_length);
        Self::filter_sv(&mut data.svrdup, data.max_read_length);
        for sv_list in data.svffus.values_mut() {
            Self::filter_sv(sv_list, data.max_read_length);
        }
        for sv_list in data.svrfus.values_mut() {
            Self::filter_sv(sv_list, data.max_read_length);
        }

        data.softp2sv_first_used.clear();
        Self::collect_softp2sv_first_used(&mut data.softp2sv_first_used, &data.svfinv3);
        Self::collect_softp2sv_first_used(&mut data.softp2sv_first_used, &data.svrinv3);
        Self::collect_softp2sv_first_used(&mut data.softp2sv_first_used, &data.svfinv5);
        Self::collect_softp2sv_first_used(&mut data.softp2sv_first_used, &data.svrinv5);
        Self::collect_softp2sv_first_used(&mut data.softp2sv_first_used, &data.svfdel);
        Self::collect_softp2sv_first_used(&mut data.softp2sv_first_used, &data.svrdel);
        Self::collect_softp2sv_first_used(&mut data.softp2sv_first_used, &data.svfdup);
        Self::collect_softp2sv_first_used(&mut data.softp2sv_first_used, &data.svrdup);
        for sv_list in data.svffus.values() {
            Self::collect_softp2sv_first_used(&mut data.softp2sv_first_used, sv_list);
        }
        for sv_list in data.svrfus.values() {
            Self::collect_softp2sv_first_used(&mut data.softp2sv_first_used, sv_list);
        }
    }

    fn collect_softp2sv_first_used(
        dest: &mut HashMap<i64, bool, LibDefaultHasher>,
        sv_list: &[SoftClip],
    ) {
        let mut best_depth: HashMap<i64, usize, LibDefaultHasher> = Default::default();

        for sv in sv_list {
            let softp = sv.softp as i64;
            if softp == 0 {
                continue;
            }

            let depth = sv.var.alt_depth;
            match best_depth.get(&softp).copied() {
                None => {
                    best_depth.insert(softp, depth);
                    dest.insert(softp, sv.used());
                }
                Some(current_depth) if depth > current_depth => {
                    best_depth.insert(softp, depth);
                    dest.insert(softp, sv.used());
                }
                _ => {}
            }
        }
    }

    fn filter_sv(sv_list: &mut [SoftClip], max_read_length: usize) {
        for sv in sv_list.iter_mut() {
            let cluster = Self::check_cluster(&mut sv.mates, max_read_length);

            if cluster.mate_start != 0 {
                sv.mstart = cluster.mate_start;
                sv.mend = cluster.mate_end;
                sv.var.alt_depth = cluster.cnt;
                sv.mlen = cluster.mate_len;
                sv.start = cluster.start;
                sv.end = cluster.end;
                sv.var.mean_pos = cluster.mean_pos;
                sv.var.mean_qual = cluster.mean_qual;
                sv.var.mean_mapq = cluster.mean_mapq;
                sv.var.nm = cluster.nm;
            } else {
                sv.mark_used();
            }

            if sv.disc != 0 {
                let ratio = sv.var.alt_depth as f64 / sv.disc as f64;
                if ratio < 0.5 && !(ratio >= 0.35 && sv.var.alt_depth >= 5) {
                    sv.mark_used();
                }
            }

            let mut soft: Vec<(i64, usize)> = sv.soft.iter().map(|(p, c)| (*p, *c)).collect();
            soft.sort_by(|a, b| b.1.cmp(&a.1));
            sv.softp = soft.first().map(|(p, _)| *p as i32).unwrap_or(0);
        }
    }

    fn check_cluster(mates: &mut [crate::variants::variants::Mate], read_len: usize) -> SvCluster {
        if mates.is_empty() {
            return SvCluster::default();
        }

        mates.sort_by(|a, b| a.mate_start.cmp(&b.mate_start));

        let first = &mates[0];
        let mut clusters = vec![SvCluster {
            mate_start: first.mate_start,
            mate_end: first.mate_end,
            cnt: 0,
            mate_len: 0,
            start: first.start,
            end: first.end,
            mean_pos: 0.0,
            mean_qual: 0.0,
            mean_mapq: 0.0,
            nm: 0.0,
        }];

        let mut cur = 0usize;
        let min_sv_dist = (Configuration::MINSVCDIST * read_len as f64) as i64;

        for mate in mates.iter() {
            if mate.mate_start - clusters[cur].mate_end > min_sv_dist {
                cur += 1;
                clusters.push(SvCluster {
                    mate_start: mate.mate_start,
                    mate_end: mate.mate_end,
                    cnt: 0,
                    mate_len: 0,
                    start: mate.start,
                    end: mate.end,
                    mean_pos: 0.0,
                    mean_qual: 0.0,
                    mean_mapq: 0.0,
                    nm: 0.0,
                });
            }

            let current = &mut clusters[cur];
            current.cnt += 1;
            current.mate_len += mate.mate_len;

            if mate.mate_end > current.mate_end {
                current.mate_end = mate.mate_end;
            }
            if mate.start < current.start {
                current.start = mate.start;
            }
            if mate.end > current.end {
                current.end = mate.end;
            }
            current.mean_pos += mate.mean_pos;
            current.mean_qual += mate.mean_qual;
            current.mean_mapq += mate.mean_mapq;
            current.nm += mate.nm;
        }

        clusters.sort_by(|a, b| b.cnt.cmp(&a.cnt));
        let first_cluster = &clusters[0];

        if first_cluster.cnt as f64 / mates.len() as f64 >= 0.60 {
            SvCluster {
                mate_start: first_cluster.mate_start,
                mate_end: first_cluster.mate_end,
                cnt: first_cluster.cnt,
                mate_len: first_cluster.mate_len / first_cluster.cnt as i32,
                start: first_cluster.start,
                end: first_cluster.end,
                mean_pos: first_cluster.mean_pos,
                mean_qual: first_cluster.mean_qual,
                mean_mapq: first_cluster.mean_mapq,
                nm: first_cluster.nm,
            }
        } else {
            SvCluster::default()
        }
    }

    fn minus_amp_base_desc(desc: &str) -> Option<String> {
        let (base, suffix) = desc.split_once('&')?;
        if !base.starts_with('-') || base.len() <= 1 {
            return None;
        }
        if !base[1..].as_bytes().iter().all(|b| b.is_ascii_digit()) {
            return None;
        }
        if suffix.is_empty()
            || !suffix
                .as_bytes()
                .iter()
                .all(|b| matches!(*b, b'A' | b'T' | b'G' | b'C'))
        {
            return None;
        }
        Some(base.to_string())
    }

    fn merge_duplicate_deletion_keys(pos_map: &mut VariantMap) {
        let mut grouped: HashMap<String, Vec<VarDesc>, LibDefaultHasher> = Default::default();
        for key in pos_map.keys() {
            let key_str = key.to_key_string();
            if key_str.starts_with('-') {
                grouped.entry(key_str).or_default().push(key.clone());
            }
        }

        for keys in grouped.into_values() {
            if keys.len() <= 1 {
                continue;
            }

            let mut primary = keys[0].clone();
            let mut primary_is_raw = matches!(primary, VarDesc::Raw { .. });
            let mut primary_depth = pos_map.get(&primary).map(|v| v.alt_depth).unwrap_or(0);

            for key in &keys[1..] {
                let key_is_raw = matches!(key, VarDesc::Raw { .. });
                let key_depth = pos_map.get(key).map(|v| v.alt_depth).unwrap_or(0);
                if (primary_is_raw && !key_is_raw)
                    || (primary_is_raw == key_is_raw && key_depth > primary_depth)
                {
                    primary = key.clone();
                    primary_is_raw = key_is_raw;
                    primary_depth = key_depth;
                }
            }

            let mut merged_variant: Option<Variant> = None;
            for key in keys {
                if let Some(variant) = pos_map.remove(&key) {
                    if let Some(current) = merged_variant.as_mut() {
                        let previous_current_pstd = current.pstd;
                        let previous_current_pp = current.pp;
                        let incoming_pstd = variant.pstd;
                        let incoming_pp = variant.pp;

                        adj_cnt(current, &variant);
                        current.extra_cnt = current
                            .extra_cnt
                            .saturating_sub(variant.alt_depth)
                            .saturating_add(variant.extra_cnt);
                        current.pstd = previous_current_pstd
                            || incoming_pstd
                            || (previous_current_pp > 0
                                && incoming_pp > 0
                                && previous_current_pp != incoming_pp);
                    } else {
                        merged_variant = Some(variant);
                    }
                }
            }

            if let Some(variant) = merged_variant {
                pos_map.insert(primary, variant);
            }
        }
    }

    fn merge_suffix_shift_insertion_keys(
        pos_map: &mut VariantMap,
        position_to_insertion_count: Option<&CountMap>,
    ) {
        let mut keys: Vec<(VarDesc, Vec<u8>)> = pos_map
            .keys()
            .filter_map(|key| {
                if let VarDesc::Ins { seq } = key {
                    Some((key.clone(), seq.as_slice().to_vec()))
                } else {
                    None
                }
            })
            .collect();

        if keys.len() <= 1 {
            return;
        }

        keys.sort_by(|a, b| b.1.len().cmp(&a.1.len()));

        let mut merges: Vec<(VarDesc, VarDesc)> = Vec::new();

        for i in 0..keys.len() {
            let (long_key, long_seq) = (&keys[i].0, &keys[i].1);
            let long_depth = pos_map.get(long_key).map(|v| v.alt_depth).unwrap_or(0);
            if long_depth < 3 {
                continue;
            }

            for j in (i + 1)..keys.len() {
                let (short_key, short_seq) = (&keys[j].0, &keys[j].1);
                let short_depth = pos_map.get(short_key).map(|v| v.alt_depth).unwrap_or(0);
                let short_supported = position_to_insertion_count.is_some_and(|counts| {
                    counts.contains_key(&format!("+{}", String::from_utf8_lossy(short_seq)))
                });
                let long_supported = position_to_insertion_count.is_some_and(|counts| {
                    counts.contains_key(&format!("+{}", String::from_utf8_lossy(long_seq)))
                });

                if short_depth != 1 {
                    continue;
                }
                if short_supported && long_supported {
                    continue;
                }
                if long_seq.len() != short_seq.len() + 1 {
                    continue;
                }
                if long_depth < short_depth.saturating_mul(3) {
                    continue;
                }
                if long_seq[1..] != short_seq[..] {
                    continue;
                }
                merges.push((short_key.clone(), long_key.clone()));
                break;
            }
        }

        for (from_key, to_key) in merges {
            if from_key == to_key {
                continue;
            }

            let Some(src) = pos_map.remove(&from_key) else {
                continue;
            };

            if let Some(dest) = pos_map.get_mut(&to_key) {
                adj_cnt(dest, &src);
                dest.extra_cnt = dest
                    .extra_cnt
                    .saturating_sub(src.alt_depth)
                    .saturating_add(src.extra_cnt);
            } else {
                pos_map.insert(from_key, src);
            }
        }
    }

    fn has_sv_marker(variation_map: &VariantMap) -> bool {
        variation_map
            .keys()
            .any(|desc| matches!(desc, VarDesc::Raw { desc } if desc.as_slice() == b"SV"))
    }

    fn ensure_sv_marker(data: &mut RealignedVariationData, position: i64) {
        let variation_map = data.non_insertion_variants.entry(position).or_default();
        variation_map
            .entry(VarDesc::Raw {
                desc: b"SV".to_vec().into(),
            })
            .or_default();
        data.sv_counts.entry(position).or_default();
    }

    fn add_sv_split_count(data: &mut RealignedVariationData, position: i64, splits: usize) {
        Self::ensure_sv_marker(data, position);
        if splits == 0 {
            return;
        }
        if let Some(sv) = data.sv_counts.get_mut(&position) {
            sv.splits += splits;
        }
    }

    fn add_sv_counts(
        data: &mut RealignedVariationData,
        position: i64,
        pairs: usize,
        splits: usize,
        clusters: usize,
    ) {
        Self::ensure_sv_marker(data, position);
        if let Some(sv) = data.sv_counts.get_mut(&position) {
            sv.pairs += pairs;
            sv.splits += splits;
            sv.clusters += clusters;
        }
    }

    fn adjust_sv_splits(data: &mut RealignedVariationData, position: i64, delta: isize) {
        if delta == 0 {
            return;
        }

        let Some(sv) = data.sv_counts.get_mut(&position) else {
            return;
        };

        if delta > 0 {
            sv.splits += delta as usize;
        } else {
            sv.splits = sv.splits.saturating_sub((-delta) as usize);
        }
    }

    fn remove_sv_marker(data: &mut RealignedVariationData, position: i64) {
        let key = VarDesc::Raw {
            desc: b"SV".to_vec().into(),
        };
        data.sv_counts.remove(&position);
        if let Some(map) = data.non_insertion_variants.get_mut(&position) {
            map.remove(&key);
            if map.is_empty() {
                data.non_insertion_variants.remove(&position);
            }
        }
    }

    fn move_sv_marker(data: &mut RealignedVariationData, from: i64, to: i64) {
        if from == to {
            return;
        }

        let key = VarDesc::Raw {
            desc: b"SV".to_vec().into(),
        };

        let from_sv = data
            .non_insertion_variants
            .get_mut(&from)
            .and_then(|map| map.remove(&key));

        if let Some(map) = data.non_insertion_variants.get_mut(&from) {
            if map.is_empty() {
                data.non_insertion_variants.remove(&from);
            }
        }

        let from_counts = data.sv_counts.remove(&from);

        let Some(from_sv) = from_sv else {
            if let Some(from_counts) = from_counts {
                data.sv_counts.insert(to, from_counts);
            }
            return;
        };

        let to_map = data.non_insertion_variants.entry(to).or_default();
        to_map.insert(key, from_sv);
        if let Some(from_counts) = from_counts {
            data.sv_counts.insert(to, from_counts);
        }
    }

    fn is_overlap(start1: i64, end1: i64, start2: i64, end2: i64, read_len: i64) -> bool {
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

        positions[1] - positions[0] + positions[3] - positions[2] < 3 * read_len
    }

    fn mark_sv_clusters(
        start: i64,
        end: i64,
        clusters: &mut [crate::variants::variants::SoftClip],
        read_len: i64,
    ) -> (usize, usize, usize) {
        let mut cov = 0usize;
        let mut cnt = 0usize;
        let mut pairs = 0usize;

        for cluster in clusters.iter_mut() {
            let (start2, end2) = if cluster.start < cluster.mstart {
                (cluster.end, cluster.mstart)
            } else {
                (cluster.mend, cluster.start)
            };

            if Self::is_overlap(start, end, start2, end2, read_len) {
                cluster.mark_used();
                cnt += 1;
                pairs += cluster.var.alt_depth;
                if cluster.end != cluster.start {
                    cov += ((cluster.var.alt_depth as i64 * read_len)
                        / (cluster.end - cluster.start)) as usize
                        + 1;
                }
            }
        }

        (cov, cnt, pairs)
    }

    fn mark_dup_sv(
        start: i64,
        end: i64,
        clusters: &mut [crate::variants::variants::SoftClip],
        read_len: i64,
    ) -> (usize, usize) {
        let mut cnt = 0usize;
        let mut pairs = 0usize;

        for cluster in clusters.iter_mut() {
            let (start2, end2) = if cluster.start < cluster.mstart {
                (cluster.start, cluster.mend)
            } else {
                (cluster.mstart, cluster.end)
            };

            if Self::is_overlap(start, end, start2, end2, read_len) {
                cluster.mark_used();
                cnt += 1;
                pairs += cluster.var.alt_depth;
            }
        }

        (cnt, pairs)
    }

    fn rm_cnt(dest: &mut Variant, src: &Variant) {
        dest.alt_depth = dest.alt_depth.saturating_sub(src.alt_depth);
        dest.alt_depth_fwd = dest.alt_depth_fwd.saturating_sub(src.alt_depth_fwd);
        dest.alt_depth_rev = dest.alt_depth_rev.saturating_sub(src.alt_depth_rev);
        dest.extra_cnt = dest.extra_cnt.saturating_sub(src.extra_cnt);
        dest.mean_pos = (dest.mean_pos - src.mean_pos).max(0.0);
        dest.mean_qual = (dest.mean_qual - src.mean_qual).max(0.0);
        dest.mean_mapq = (dest.mean_mapq - src.mean_mapq).max(0.0);
        dest.low_qual_read_cnt = dest.low_qual_read_cnt.saturating_sub(src.low_qual_read_cnt);
        dest.high_qual_read_cnt = dest
            .high_qual_read_cnt
            .saturating_sub(src.high_qual_read_cnt);
    }

    fn adjust_sv_split_count(data: &mut RealignedVariationData, position: i64, delta: isize) {
        if delta == 0 {
            return;
        }

        let Some(sv) = data.sv_counts.get_mut(&position) else {
            return;
        };

        if delta > 0 {
            sv.splits += delta as usize;
        } else {
            sv.splits = sv.splits.saturating_sub((-delta) as usize);
        }
    }

    fn has_deletion_like_non_insertion_key(data: &RealignedVariationData, position: i64) -> bool {
        data.non_insertion_variants
            .get(&position)
            .map(|variation_map| {
                variation_map
                    .keys()
                    .any(|key| key.to_key_string().starts_with('-'))
            })
            .unwrap_or(false)
    }

    fn is_single_base_ampersand_deletion_key(key: &VarDesc) -> bool {
        matches!(key, VarDesc::Raw { desc } if desc.as_slice().starts_with(b"-1&"))
    }

    fn updated_insertion_count_after_realign(
        before: Option<&VariantMap>,
        after: Option<&VariantMap>,
        original_key: &VarDesc,
    ) -> Option<usize> {
        let before_map = before?;
        let after_map = after?;

        if before_map.len() != 1 || !before_map.contains_key(original_key) {
            return None;
        }

        let mut positive_candidate: Option<usize> = None;

        for (key, variant) in after_map {
            let before_count = before_map.get(key).map(|v| v.alt_depth).unwrap_or(0);
            if variant.alt_depth > before_count {
                if positive_candidate.is_some() {
                    return None;
                }
                positive_candidate = Some(variant.alt_depth);
            }
        }

        positive_candidate
    }

    fn resolved_insertion_key_after_realign(
        before: Option<&VariantMap>,
        after: Option<&VariantMap>,
        original_key: &VarDesc,
    ) -> Option<VarDesc> {
        let after_map = after?;

        if after_map.contains_key(original_key) {
            return Some(original_key.clone());
        }

        let mut positive_candidate: Option<VarDesc> = None;
        for (key, variant) in after_map {
            let before_count = before
                .and_then(|map| map.get(key))
                .map(|v| v.alt_depth)
                .unwrap_or(0);
            let delta = variant.alt_depth as isize - before_count as isize;
            if delta > 0 {
                if positive_candidate.is_some() {
                    return None;
                }
                positive_candidate = Some(key.clone());
            }
        }

        positive_candidate
    }

    fn emit_realigner_insertion_diag(
        &self,
        _data: &RealignedVariationData,
        _position: i64,
        _source: &str,
    ) {
        let _ = self;
    }

    fn emit_refcov_inc_diag(
        _data: &RealignedVariationData,
        _position: i64,
        _amount: usize,
        _source: &str,
    ) {
    }

    fn emit_refcov_inc_diag_with_before(
        _before: usize,
        _position: i64,
        _amount: usize,
        _source: &str,
    ) {
    }

    pub fn load_partial_ref_coverage(
        &self,
        data: &mut RealignedVariationData,
        start: i64,
        end: i64,
    ) {
        // Preserve Java boundary behavior for on-demand coverage loads: a fully
        // left-of-chromosome window must stay empty, not clamp to position 1
        // and replay already-counted reference evidence.
        if self.bam_paths.is_empty() || start > end || end < 1 {
            return;
        }

        let Some(chromosome) = self.chromosome.as_deref() else {
            return;
        };

        let start = start.max(1);
        let end = end.max(start);
        let region = Region::new(
            chromosome.to_string(),
            start as usize,
            end as usize,
            String::new(),
        );

        let mut reference =
            Reference::new_with_start(self.reference_seq.as_ref().clone(), self.ref_start);
        let ref_end = reference.region_start + reference.ref_seq.len() as i64 - 1;
        let chr_len = crate::scopedata::global_read_only_scope::instance()
            .chr_lens
            .get(chromosome)
            .copied();
        reference.build_seed_map(ref_end, chr_len);

        let pipeline = VarDictPipeline::new("realigner_partial");
        let Ok(extra) = pipeline.run_partial_cigar_for_bams(&region, &reference, &self.bam_paths)
        else {
            return;
        };

        Self::merge_variant_maps(&mut data.non_insertion_variants, extra.non_insertion_vars);
        Self::merge_variant_maps(&mut data.insertion_variants, extra.insertion_vars);

        Self::merge_soft_clips(&mut data.soft_clips_5end, extra.soft_clips_5end);
        Self::merge_soft_clips(&mut data.soft_clips_3end, extra.soft_clips_3end);

        for (position, coverage) in extra.ref_coverage {
            Self::emit_refcov_inc_diag(data, position, coverage, "realigner_merge");
            *data.ref_coverage.entry(position).or_insert(0) += coverage;
        }
    }

    fn extract_inv_flanks(desc_str: &str) -> Option<(String, String)> {
        let caret_idx = desc_str.find('^')?;
        let tail = &desc_str[caret_idx + 1..];
        let tail_lower = tail.to_ascii_lowercase();
        let inv_tag_idx = tail_lower.find("<inv")?;
        let inv_tail = &tail[inv_tag_idx..];
        let close_rel = inv_tail.find('>')?;

        let left = &tail[..inv_tag_idx];
        let right = &inv_tail[close_rel + 1..];
        if left.is_empty() || right.is_empty() {
            return None;
        }

        let is_atgnc = |s: &str| {
            s.as_bytes()
                .iter()
                .all(|base| matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'N'))
        };

        if !is_atgnc(left) || !is_atgnc(right) {
            return None;
        }

        Some((left.to_string(), right.to_string()))
    }

    fn merge_variant_maps(
        dest: &mut VariantMapByPos,
        src: VariantMapByPos,
    ) {
        for (position, src_map) in src {
            let dest_map = dest.entry(position).or_default();
            for (desc, src_var) in src_map {
                match dest_map.entry(desc) {
                    VecEntry::Occupied(mut occupied) => adj_cnt(occupied.get_mut(), &src_var),
                    VecEntry::Vacant(vacant) => {
                        vacant.insert(src_var);
                    }
                }
            }
        }
    }

    fn merge_soft_clips(
        dest: &mut HashMap<i64, SoftClip, LibDefaultHasher>,
        src: HashMap<i64, SoftClip, LibDefaultHasher>,
    ) {
        for (position, src_sc) in src {
            if let Some(dest_sc) = dest.get_mut(&position) {
                Self::merge_soft_clip_entry(dest_sc, &src_sc);
            } else {
                dest.insert(position, src_sc);
            }
        }
    }

    fn merge_soft_clip_entry(dest: &mut SoftClip, src: &SoftClip) {
        adj_cnt(&mut dest.var, &src.var);

        for (offset, src_base_map) in &src.nt {
            let dest_base_map = dest.nt.entry(*offset).or_default();
            for base in [b'A', b'C', b'G', b'T', b'N'] {
                if let Some(src_count) = src_base_map.get(base) {
                    if let Some(dest_count) = dest_base_map.get_or_insert_with(base, || 0) {
                        *dest_count += *src_count;
                    }
                }
            }
        }

        for (offset, src_base_map) in &src.seq {
            let dest_base_map = dest.seq.entry(*offset).or_default();
            for base in [b'A', b'C', b'G', b'T', b'N'] {
                if let Some(src_var) = src_base_map.get(base) {
                    if let Some(dest_var) = dest_base_map.get_or_insert_with(base, Variant::default)
                    {
                        adj_cnt(dest_var, src_var);
                    }
                }
            }
        }

        if !dest.consensus_seq_is_set() && src.consensus_seq_is_set() {
            dest.set_consensus_seq(src.consensus_seq().to_vec());
        }

        if src.used() {
            dest.mark_used();
        }

        if dest.start == 0 || (src.start != 0 && src.start < dest.start) {
            dest.start = src.start;
        }
        if dest.end == 0 || src.end > dest.end {
            dest.end = src.end;
        }
        if dest.mstart == 0 || (src.mstart != 0 && src.mstart < dest.mstart) {
            dest.mstart = src.mstart;
        }
        if dest.mend == 0 || src.mend > dest.mend {
            dest.mend = src.mend;
        }
        if dest.mlen == 0 {
            dest.mlen = src.mlen;
        }
        dest.disc += src.disc;
        if dest.softp == 0 {
            dest.softp = src.softp;
        }

        for (soft_pos, soft_cnt) in &src.soft {
            *dest.soft.entry(*soft_pos).or_insert(0) += *soft_cnt;
        }
        if !src.mates.is_empty() {
            dest.mates.extend(src.mates.clone());
        }
    }

    pub fn process_insertions(
        &self,
        data: &mut RealignedVariationData,
        position_to_insertion_count: &CountMapByPos,
    ) {
        if position_to_insertion_count.is_empty() {
            return;
        }

        let positions: Vec<i64> = data.insertion_variants.keys().copied().collect();
        for position in positions {
            if let Some(pos_map) = data.insertion_variants.get_mut(&position) {
                if let Some(insertion_count_at_position) =
                    position_to_insertion_count.get(&position)
                {
                    Self::merge_suffix_shift_insertion_keys(
                        pos_map,
                        Some(insertion_count_at_position),
                    );
                }
            }
        }

        let mut tmp: Vec<(i64, String, usize)> = Vec::new();
        for (pos, desc_map) in position_to_insertion_count {
            for (desc, count) in desc_map {
                tmp.push((*pos, desc.clone(), *count));
            }
        }
        tmp.sort_by(|a, b| {
            b.2.cmp(&a.2)
                .then_with(|| a.0.cmp(&b.0))
                .then_with(|| b.1.cmp(&a.1))
        });

        for (position, vn, insertion_count) in tmp.iter().cloned() {
            if insertion_count == 0 {
                continue;
            }
            let Some(insert_match) = BEGIN_PLUS_ATGC.captures(&vn) else {
                continue;
            };
            let insert = insert_match.get(1).map(|m| m.as_str()).unwrap_or("");
            if insert.is_empty() {
                continue;
            }

            let mut ins3 = String::new();
            let mut inslen = insert.len();
            if let Some(cap) = DUP_NUM_ATGC.captures(&vn) {
                let dup_len: usize = cap
                    .get(1)
                    .map(|m| m.as_str())
                    .unwrap_or("0")
                    .parse()
                    .unwrap_or(0);
                ins3 = cap.get(2).map(|m| m.as_str()).unwrap_or("").to_string();
                inslen += dup_len + ins3.len();
            }

            let extra = AMP_ATGC
                .captures(&vn)
                .and_then(|c| c.get(1))
                .map(|m| m.as_str().to_string())
                .unwrap_or_default();

            let compm = HASH_ATGC
                .captures(&vn)
                .and_then(|c| c.get(1))
                .map(|m| m.as_str().to_string())
                .unwrap_or_default();

            let _newins = CARET_ATGC_END
                .captures(&vn)
                .and_then(|c| c.get(1))
                .map(|m| m.as_str().to_string())
                .unwrap_or_default();

            let newdel: i64 = UP_NUMBER_END
                .captures(&vn)
                .and_then(|c| c.get(1))
                .map(|m| m.as_str().parse::<i64>().unwrap_or(0))
                .unwrap_or(0);

            let mut tn = vn.clone();
            if tn.starts_with('+') {
                tn.remove(0);
            }
            if let Some(idx) = tn.find('&') {
                tn.replace_range(idx..=idx, "");
            }
            if let Some(idx) = tn.find('#') {
                tn.replace_range(idx..=idx, "");
            }
            if let Some(mat) = UP_NUMBER_END.find(&tn) {
                tn.replace_range(mat.start()..mat.end(), "");
            }
            if let Some(idx) = tn.find('^') {
                tn.replace_range(idx..=idx, "");
            }

            let wustart = (position - 150).max(1);
            let mut wupseq = self.get_ref_range(wustart, position);
            wupseq.extend_from_slice(tn.as_bytes());

            let ref_end = self.ref_start + self.reference_seq.len() as i64 - 1;
            let mut sanend = position + vn.len() as i64 + 100;
            if sanend > ref_end {
                sanend = ref_end;
            }

            let (sanpseq, p3) = if !ins3.is_empty() {
                let mut sanp = Vec::new();
                let p3 =
                    position + inslen as i64 - ins3.len() as i64 + Configuration::SVFLANK as i64;
                if ins3.len() > Configuration::SVFLANK as usize {
                    sanp.extend_from_slice(&Self::substr_bytes(
                        ins3.as_bytes(),
                        Configuration::SVFLANK as i64 - ins3.len() as i64,
                        None,
                    ));
                }
                sanp.extend_from_slice(&self.get_ref_range(position + 1, position + 101));
                (sanp, p3 + 1)
            } else {
                let mut sanp = tn.as_bytes().to_vec();
                let start = position + extra.len() as i64 + 1 + compm.len() as i64 + newdel;
                sanp.extend_from_slice(&self.get_ref_range(start, sanend));
                (sanp, position + 1)
            };

            let r3 = self.find_mm3(p3, &sanpseq, data);
            let mm5_pos = position + extra.len() as i64 + compm.len() as i64 + newdel;
            let r5 = self.find_mm5(mm5_pos, &wupseq, data);

            let mut all_mm = Vec::new();
            all_mm.extend(r3.mismatches.iter().cloned());
            all_mm.extend(r5.mismatches.iter().cloned());

            let vn_key = VarDesc::Ins {
                seq: normalize_insertion_key_bytes(vn[1..].as_bytes()),
            };
            let mut active_vn_key = vn_key.clone();
            let _ = crate::variants::var_utils::get_variants_from_map(
                &mut data.insertion_variants,
                position,
                &vn_key,
            );

            for mm in all_mm {
                let mm_bytes: Vec<u8> = mm
                    .mismatch_sequence
                    .as_bytes()
                    .iter()
                    .map(|b| b.to_ascii_uppercase())
                    .collect();
                if mm_bytes.is_empty() {
                    continue;
                }

                let key = if mm_bytes.len() == 1 {
                    VarDesc::SNV {
                        ref_base: mm_bytes[0],
                    }
                } else {
                    let mut raw = Vec::with_capacity(mm_bytes.len() + 1);
                    raw.push(mm_bytes[0]);
                    raw.push(b'&');
                    raw.extend_from_slice(&mm_bytes[1..]);
                    VarDesc::Raw { desc: raw.into() }
                };

                let tv_snapshot = data
                    .non_insertion_variants
                    .get(&mm.mismatch_position)
                    .and_then(|m| m.get(&key))
                    .cloned();
                let Some(tv) = tv_snapshot else {
                    continue;
                };

                if tv.alt_depth == 0 {
                    continue;
                }
                let mean_qual = tv.mean_qual / tv.alt_depth as f64;
                if mean_qual
                    < crate::scopedata::global_read_only_scope::instance()
                        .conf
                        .goodq
                {
                    continue;
                }
                let nm_threshold = if mm.end == 3 { r3.nm } else { r5.nm };
                let mean_pos = tv.mean_pos / tv.alt_depth as f64;
                if mean_pos > nm_threshold as f64 + 4.0 {
                    continue;
                }
                if tv.alt_depth >= insertion_count + insert.len()
                    || (tv.alt_depth as f64 / insertion_count as f64) >= 8.0
                {
                    continue;
                }

                if mm.mismatch_position > position && mm.end == 5 {
                    Self::emit_refcov_inc_diag(
                        data,
                        position,
                        tv.alt_depth,
                        "realigner_process_insertions_mm_end5",
                    );
                    *data.ref_coverage.entry(position).or_insert(0) += tv.alt_depth;
                }

                let tv_owned = {
                    let tv_owned = data
                        .non_insertion_variants
                        .get_mut(&mm.mismatch_position)
                        .and_then(|m| m.remove(&key));
                    if let Some(map) = data.non_insertion_variants.get_mut(&mm.mismatch_position) {
                        if map.is_empty() {
                            data.non_insertion_variants.remove(&mm.mismatch_position);
                        }
                    }
                    tv_owned
                };
                let Some(tv_owned) = tv_owned else {
                    continue;
                };

                let ref_var = if mm.mismatch_position > position && mm.end == 3 {
                    self.get_ref_base(position)
                        .map(|ref_base| VarDesc::SNV { ref_base })
                        .and_then(|ref_key| {
                            data.non_insertion_variants
                                .get_mut(&position)
                                .and_then(|m| m.remove(&ref_key))
                                .map(|var| (ref_key, var))
                        })
                } else {
                    None
                };

                if let Some(vref) = data
                    .insertion_variants
                    .get_mut(&position)
                    .and_then(|m| m.get_mut(&active_vn_key))
                {
                    if let Some((ref_key, mut ref_var)) = ref_var {
                        adj_cnt_with_ref(vref, &tv_owned, Some(&mut ref_var));
                        if let Some(pos_map) = data.non_insertion_variants.get_mut(&position) {
                            pos_map.insert(ref_key, ref_var);
                        }
                    } else {
                        adj_cnt(vref, &tv_owned);
                    }
                }
            }

            if r3.misp != 0 && r3.mismatches.len() == 1 {
                if let Some(base) = r3.misnt {
                    let key = VarDesc::SNV { ref_base: base };
                    if let Some(map) = data.non_insertion_variants.get_mut(&r3.misp) {
                        if let Some(var) = map.get(&key) {
                            if var.alt_depth < insertion_count {
                                map.remove(&key);
                            }
                        }
                    }
                }
            }

            if r5.misp != 0 && r5.mismatches.len() == 1 {
                if let Some(base) = r5.misnt {
                    let key = VarDesc::SNV { ref_base: base };
                    if let Some(map) = data.non_insertion_variants.get_mut(&r5.misp) {
                        if let Some(var) = map.get(&key) {
                            if var.alt_depth < insertion_count {
                                map.remove(&key);
                            }
                        }
                    }
                }
            }

            for sc5pp in r5.scp.iter().copied() {
                let refcov_before = if position == 6970385 {
                    data.ref_coverage.get(&6970385_i64).copied().unwrap_or(0)
                } else {
                    0
                };
                if let Some(tv) = data.soft_clips_5end.get_mut(&sc5pp) {
                    if tv.used() {
                        continue;
                    }
                    let seq = find_conseq_transient(tv, 0);
                    let matched = !seq.is_empty() && Self::is_match_bytes(&seq, &wupseq, -1);
                    if matched {
                        let added_alt_depth = tv.var.alt_depth;
                        cache_transient_conseq(tv, &seq);
                        if sc5pp > position {
                            Self::emit_refcov_inc_diag_with_before(
                                refcov_before,
                                position,
                                tv.var.alt_depth,
                                "realigner_process_insertions_sc5",
                            );
                            *data.ref_coverage.entry(position).or_insert(0) += tv.var.alt_depth;
                        }
                        if let Some(vref) = data
                            .insertion_variants
                            .get_mut(&position)
                            .and_then(|m| m.get_mut(&active_vn_key))
                        {
                            adj_cnt(vref, &tv.var);
                        }
                        tv.mark_used();
                    }
                }
            }

            for sc3pp in r3.scp.iter().copied() {
                let refcov_before = if position == 6970385 {
                    data.ref_coverage.get(&6970385_i64).copied().unwrap_or(0)
                } else {
                    0
                };
                if let Some(tv) = data.soft_clips_3end.get_mut(&sc3pp) {
                    if tv.used() {
                        continue;
                    }
                    let seq = find_conseq_transient(tv, 0);
                    let mseq = if !ins3.is_empty() {
                        sanpseq.clone()
                    } else {
                        let offset = sc3pp - position - 1;
                        Self::substr_bytes(&sanpseq, offset, None)
                    };
                    let matched = !seq.is_empty() && Self::is_match_bytes(&seq, &mseq, 1);
                    if matched {
                        let added_alt_depth = tv.var.alt_depth;
                        cache_transient_conseq(tv, &seq);
                        let mean_pos = if tv.var.alt_depth > 0 {
                            tv.var.mean_pos / tv.var.alt_depth as f64
                        } else {
                            0.0
                        };
                        if sc3pp <= position || insert.len() as f64 > mean_pos {
                            Self::emit_refcov_inc_diag_with_before(
                                refcov_before,
                                position,
                                tv.var.alt_depth,
                                "realigner_process_insertions_sc3",
                            );
                            *data.ref_coverage.entry(position).or_insert(0) += tv.var.alt_depth;
                        }

                        let use_ref_var = sc3pp > position && !(insert.len() as f64 > mean_pos);

                        let ref_var = if use_ref_var {
                            self.get_ref_base(position)
                                .map(|ref_base| VarDesc::SNV { ref_base })
                                .and_then(|ref_key| {
                                    data.non_insertion_variants
                                        .get_mut(&position)
                                        .and_then(|m| m.remove(&ref_key))
                                        .map(|var| (ref_key, var))
                                })
                        } else {
                            None
                        };

                        if let Some(vref) = data
                            .insertion_variants
                            .get_mut(&position)
                            .and_then(|m| m.get_mut(&active_vn_key))
                        {
                            if let Some((ref_key, mut ref_var)) = ref_var {
                                adj_cnt_with_ref(vref, &tv.var, Some(&mut ref_var));
                                if let Some(pos_map) =
                                    data.non_insertion_variants.get_mut(&position)
                                {
                                    pos_map.insert(ref_key, ref_var);
                                }
                            } else {
                                adj_cnt(vref, &tv.var);
                            }
                        }
                        tv.mark_used();
                        if insert.len() + 1 == vn.len()
                            && insert.len() > data.max_read_length
                            && sc3pp >= position + 1 + insert.len() as i64
                        {
                            let mut flag = 0;
                            let offset = ((sc3pp - position - 1) % insert.len() as i64) as usize;
                            let mut tvn_bytes = vn.as_bytes().to_vec();
                            for (seqi, &b) in seq.iter().enumerate() {
                                if seqi + offset >= insert.len() {
                                    break;
                                }
                                if b != insert.as_bytes()[seqi + offset] {
                                    flag += 1;
                                    let shift = seqi + offset + 1;
                                    if shift < tvn_bytes.len() {
                                        tvn_bytes[shift] = b;
                                    }
                                }
                            }
                            if flag > 0 {
                                if let Some(map) = data.insertion_variants.get_mut(&position) {
                                    if let Some(variation) = map.remove(&active_vn_key) {
                                        let new_key = VarDesc::Ins {
                                            seq: normalize_insertion_key_bytes(&tvn_bytes[1..]),
                                        };
                                        map.insert(new_key.clone(), variation);
                                        active_vn_key = new_key;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if !r3.scp.is_empty() && !r5.scp.is_empty() {
                let first3 = r3.scp[0];
                let first5 = r5.scp[0];
                if first3 > first5 + 3
                    && ((first3 - first5) as f64) < data.max_read_length as f64 * 0.75
                {
                    if let Some(ref_base) = self.get_ref_base(position) {
                        if let Some(map) = data.non_insertion_variants.get_mut(&position) {
                            if let Some(ref_var) = map.get_mut(&VarDesc::SNV { ref_base }) {
                                Self::adj_ref_factor(
                                    ref_var,
                                    (first3 - first5 - 1) as f64 / data.max_read_length as f64,
                                );
                            }
                        }
                    }
                    if let Some(vref) = data
                        .insertion_variants
                        .get_mut(&position)
                        .and_then(|m| m.get_mut(&active_vn_key))
                    {
                        Self::adj_ref_factor(
                            vref,
                            -((first3 - first5 - 1) as f64 / data.max_read_length as f64),
                        );
                    }
                }
            }
        }

        for idx in (1..tmp.len()).rev() {
            let (position, vn, _) = tmp[idx].clone();
            let Some(map) = data.insertion_variants.get_mut(&position) else {
                continue;
            };
            let vref_key = if vn.starts_with('+') {
                VarDesc::Ins {
                    seq: normalize_insertion_key_bytes(vn[1..].as_bytes()),
                }
            } else {
                continue;
            };
            let vref_alt = map.get(&vref_key).map(|v| v.alt_depth);
            if let Some(cap) = ATGSs_AMP_ATGSs_END.captures(&vn) {
                let tn = cap.get(1).map(|m| m.as_str()).unwrap_or("");
                if !tn.is_empty() {
                    let tref_key = VarDesc::Ins {
                        seq: normalize_insertion_key_bytes(tn[1..].as_bytes()),
                    };
                    let tref_alt = map.get(&tref_key).map(|v| v.alt_depth);
                    if let (Some(vref_alt), Some(tref_alt)) = (vref_alt, tref_alt) {
                        if vref_alt < tref_alt {
                            let vref = map.remove(&vref_key);
                            if let Some(vref) = vref {
                                if let Some(tref) = map.get_mut(&tref_key) {
                                    let ref_var = self
                                        .get_ref_base(position)
                                        .map(|ref_base| VarDesc::SNV { ref_base })
                                        .and_then(|ref_key| {
                                            data.non_insertion_variants
                                                .get_mut(&position)
                                                .and_then(|m| m.remove(&ref_key))
                                                .map(|var| (ref_key, var))
                                        });

                                    if let Some((ref_key, mut ref_var)) = ref_var {
                                        adj_cnt_with_ref(tref, &vref, Some(&mut ref_var));
                                        if let Some(pos_map) =
                                            data.non_insertion_variants.get_mut(&position)
                                        {
                                            pos_map.insert(ref_key, ref_var);
                                        }
                                    } else {
                                        adj_cnt(tref, &vref);
                                    }
                                } else {
                                    map.insert(vref_key, vref);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    pub fn realign_long_insertions_30(&self, data: &mut RealignedVariationData) {
        let mut tmp5: Vec<(i64, usize)> = data
            .soft_clips_5end
            .iter()
            .map(|(pos, sc)| (*pos, sc.var.alt_depth))
            .collect();
        tmp5.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

        let mut tmp3: Vec<(i64, usize)> = data
            .soft_clips_3end
            .iter()
            .map(|(pos, sc)| (*pos, sc.var.alt_depth))
            .collect();
        tmp3.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

        for (p5, cnt5) in tmp5 {
            if data
                .soft_clips_5end
                .get(&p5)
                .map(|sc| sc.used())
                .unwrap_or(false)
            {
                continue;
            }

            for (p3, cnt3) in tmp3.iter().copied() {
                if data
                    .soft_clips_5end
                    .get(&p5)
                    .map(|sc| sc.used())
                    .unwrap_or(false)
                {
                    break;
                }
                if data
                    .soft_clips_3end
                    .get(&p3)
                    .map(|sc| sc.used())
                    .unwrap_or(false)
                {
                    continue;
                }
                if (p5 - p3) as f64 > data.max_read_length as f64 * 2.5 {
                    continue;
                }
                if p3 - p5 > data.max_read_length as i64 - 10 {
                    continue;
                }

                let seq5 = {
                    let Some(sc5v) = data.soft_clips_5end.get_mut(&p5) else {
                        continue;
                    };
                    find_conseq(sc5v, 5)
                };
                let seq3 = {
                    let Some(sc3v) = data.soft_clips_3end.get_mut(&p3) else {
                        continue;
                    };
                    find_conseq(sc3v, 3)
                };

                if seq5.len() <= 10 || seq3.len() <= 10 {
                    continue;
                }
                if cnt3 == 0 {
                    continue;
                }
                let ratio = cnt5 as f64 / cnt3 as f64;
                if ratio < 0.08 || ratio > 12.0 {
                    continue;
                }

                let match35 = Self::find_35_match(
                    &String::from_utf8_lossy(&seq5),
                    &String::from_utf8_lossy(&seq3),
                );
                let score = match35.max_matched_length;
                if score == 0 {
                    continue;
                }
                let smscore = score / 2;

                let mut ins = if match35.matched_3_end + smscore > 1 {
                    Self::substr_bytes(
                        &seq3,
                        0,
                        Some(-(match35.matched_3_end as i64 + smscore as i64) + 1),
                    )
                } else {
                    seq3.clone()
                };

                if match35.matched_5_end + smscore > 0 {
                    let mut part = Self::substr_bytes(
                        &seq5,
                        0,
                        Some((match35.matched_5_end + smscore) as i64),
                    );
                    part.reverse();
                    ins.extend_from_slice(&part);
                }

                if Self::is_low_complex_seq(&String::from_utf8_lossy(&ins)) {
                    continue;
                }

                let mut ins_desc = ins;
                let bi: i64;

                if p5 > p3 {
                    if seq3.len() > ins_desc.len() {
                        let tail = &seq3[ins_desc.len()..];
                        let ref_seq =
                            self.join_ref(p5, p5 + seq3.len() as i64 - ins_desc.len() as i64 + 2);
                        if !Self::is_match_bytes(tail, &ref_seq, 1) {
                            continue;
                        }
                    }
                    if seq5.len() > ins_desc.len() {
                        let tail = &seq5[ins_desc.len()..];
                        let start = p3 - (seq5.len() as i64 - ins_desc.len() as i64) - 2;
                        let ref_seq = self.join_ref(start, p3 - 1);
                        if !Self::is_match_bytes(tail, &ref_seq, -1) {
                            continue;
                        }
                    }

                    let tmp = self.join_ref(p3, p5 - 1);
                    if tmp.len() > ins_desc.len() {
                        let mut desc = (p3 - p5).to_string().into_bytes();
                        desc.push(b'^');
                        desc.extend_from_slice(&ins_desc);
                        ins_desc = desc;
                        bi = p3;
                    } else if tmp.len() < ins_desc.len() {
                        let mut desc = Self::substr_bytes(
                            &ins_desc,
                            0,
                            Some((ins_desc.len() - tmp.len()) as i64),
                        );
                        desc.push(b'&');
                        let tail = Self::substr_bytes(&ins_desc, (p3 - p5) as i64, None);
                        desc.extend_from_slice(&tail);
                        let mut with_plus = Vec::new();
                        with_plus.push(b'+');
                        with_plus.extend_from_slice(&desc);
                        ins_desc = with_plus;
                        bi = p3 - 1;
                    } else {
                        let mut desc = Vec::new();
                        desc.push(b'-');
                        desc.extend_from_slice(ins_desc.len().to_string().as_bytes());
                        desc.push(b'^');
                        desc.extend_from_slice(&ins_desc);
                        ins_desc = desc;
                        bi = p3;
                    }
                } else {
                    if seq3.len() > ins_desc.len() {
                        let tail = &seq3[ins_desc.len()..];
                        let ref_seq =
                            self.join_ref(p5, p5 + seq3.len() as i64 - ins_desc.len() as i64 + 2);
                        if !Self::is_match_bytes(tail, &ref_seq, 1) {
                            continue;
                        }
                    }
                    if seq5.len() > ins_desc.len() {
                        let tail = &seq5[ins_desc.len()..];
                        let start = p3 - (seq5.len() as i64 - ins_desc.len() as i64) - 2;
                        let ref_seq = self.join_ref(start, p3 - 1);
                        if !Self::is_match_bytes(tail, &ref_seq, -1) {
                            continue;
                        }
                    }

                    if ins_desc.len() <= (p3 - p5) as usize {
                        let mut rpt = 2i64;
                        let mut tnr = 3i64;
                        while (((p3 - p5 + ins_desc.len() as i64) as f64 / tnr as f64)
                            / ins_desc.len() as f64)
                            > 1.0
                        {
                            if (p3 - p5 + ins_desc.len() as i64) % tnr == 0 {
                                rpt += 1;
                            }
                            tnr += 1;
                        }
                        let to = p5 as f64 + (p3 - p5 + ins_desc.len() as i64) as f64 / rpt as f64
                            - ins_desc.len() as f64;
                        let tmp = self.join_ref_float(p5, to);
                        let mut desc = Vec::new();
                        desc.push(b'+');
                        desc.extend_from_slice(&tmp);
                        desc.extend_from_slice(&ins_desc);
                        ins_desc = desc;
                    } else {
                        let tmp = self.join_ref(p5, p3 - 1);
                        if (ins_desc.len() as i64 - tmp.len() as i64) % 2 == 0 {
                            let tex = (ins_desc.len() - tmp.len()) / 2;
                            let left = Self::substr_bytes(&ins_desc, 0, Some(tex as i64));
                            let right = Self::substr_bytes(&ins_desc, tex as i64, None);
                            let mut tmp_plus_left = tmp.clone();
                            tmp_plus_left.extend_from_slice(&left);
                            if tmp_plus_left == right {
                                let mut desc = Vec::new();
                                desc.push(b'+');
                                desc.extend_from_slice(&right);
                                ins_desc = desc;
                            } else {
                                let mut desc = Vec::new();
                                desc.push(b'+');
                                desc.extend_from_slice(&tmp);
                                desc.extend_from_slice(&ins_desc);
                                ins_desc = desc;
                            }
                        } else {
                            let mut desc = Vec::new();
                            desc.push(b'+');
                            desc.extend_from_slice(&tmp);
                            desc.extend_from_slice(&ins_desc);
                            ins_desc = desc;
                        }
                    }
                    bi = p5 - 1;
                }

                let sc5_var = data
                    .soft_clips_5end
                    .get(&p5)
                    .map(|sc| sc.var.clone())
                    .unwrap_or_default();
                let sc3_var = data
                    .soft_clips_3end
                    .get(&p3)
                    .map(|sc| sc.var.clone())
                    .unwrap_or_default();

                if let Some(sc3v) = data.soft_clips_3end.get_mut(&p3) {
                    sc3v.mark_used();
                }
                if let Some(sc5v) = data.soft_clips_5end.get_mut(&p5) {
                    sc5v.mark_used();
                }

                Self::emit_refcov_inc_diag(
                    data,
                    bi,
                    sc5_var.alt_depth,
                    "realigner_softclip_bridge",
                );
                *data.ref_coverage.entry(bi).or_insert(0) += sc5_var.alt_depth;

                let ins_starts_with_plus = ins_desc.first() == Some(&b'+');
                let ins_starts_with_minus = ins_desc.first() == Some(&b'-');

                if ins_starts_with_plus {
                    let key = VarDesc::Ins {
                        seq: normalize_insertion_key_bytes(&ins_desc[1..]),
                    };
                    {
                        let vref = get_variants_from_map(&mut data.insertion_variants, bi, &key);
                        vref.pstd = true;
                        vref.qstd = true;

                        let ref_key = self
                            .get_ref_base(bi)
                            .map(|ref_base| VarDesc::SNV { ref_base });
                        let mut ref_var = ref_key.as_ref().and_then(|k| {
                            data.non_insertion_variants
                                .get_mut(&bi)
                                .and_then(|m| m.remove(k))
                        });
                        if let Some(ref_var_mut) = ref_var.as_mut() {
                            adj_cnt_with_ref(vref, &sc3_var, Some(ref_var_mut));
                        } else {
                            adj_cnt(vref, &sc3_var);
                        }
                        adj_cnt(vref, &sc5_var);

                        if let (Some(ref_key), Some(ref_var)) = (ref_key, ref_var) {
                            if let Some(pos_map) = data.non_insertion_variants.get_mut(&bi) {
                                pos_map.insert(ref_key, ref_var);
                            }
                        }
                    }

                    let gap_len = p3 - p5;
                    if !self.bam_paths.is_empty()
                        && gap_len >= 5
                        && gap_len < data.max_read_length.saturating_sub(10) as i64
                    {
                        if let Some(ref_base) = self.get_ref_base(bi) {
                            let ref_key = VarDesc::SNV { ref_base };
                            let h_snapshot = data
                                .non_insertion_variants
                                .get(&bi)
                                .and_then(|m| m.get(&ref_key))
                                .cloned();
                            let vref_snapshot = data
                                .insertion_variants
                                .get(&bi)
                                .and_then(|m| m.get(&key))
                                .cloned();

                            if let (Some(h), Some(vref_now)) = (h_snapshot, vref_snapshot) {
                                if h.alt_depth != 0
                                    && self.no_passing_reads(p5, p3)
                                    && vref_now.alt_depth > 2 * h.alt_depth
                                {
                                    if let Some(pos_map) = data.non_insertion_variants.get_mut(&bi)
                                    {
                                        if let Some(mut h_work) = pos_map.remove(&ref_key) {
                                            let h_src = h_work.clone();
                                            if let Some(vref_mut) = data
                                                .insertion_variants
                                                .get_mut(&bi)
                                                .and_then(|m| m.get_mut(&key))
                                            {
                                                adj_cnt_with_ref(
                                                    vref_mut,
                                                    &h_src,
                                                    Some(&mut h_work),
                                                );
                                            }
                                            pos_map.insert(ref_key, h_work);
                                        }
                                    }
                                }
                            }
                        }
                    }

                    let current_ins_count = data
                        .insertion_variants
                        .get(&bi)
                        .and_then(|m| m.get(&key))
                        .map(|v| v.alt_depth)
                        .unwrap_or(0);

                    let mut tins: CountMapByPos = Default::default();
                    let mut map: CountMap = Default::default();
                    map.insert(
                        String::from_utf8_lossy(&ins_desc).to_string(),
                        current_ins_count,
                    );
                    tins.insert(bi, map);
                    self.process_insertions(data, &tins);
                    self.emit_realigner_insertion_diag(data, bi, "softclip_bridge");
                } else if ins_starts_with_minus {
                    let vref_key =
                        self.del_desc_to_key(&ins_desc)
                            .unwrap_or_else(|| VarDesc::Raw {
                                desc: ins_desc.clone().into(),
                            });
                    let ref_key = self
                        .get_ref_base(bi)
                        .map(|ref_base| VarDesc::SNV { ref_base });
                    let mut ref_var = ref_key.as_ref().and_then(|k| {
                        data.non_insertion_variants
                            .get_mut(&bi)
                            .and_then(|m| m.remove(k))
                    });
                    {
                        let vref =
                            get_variants_from_map(&mut data.non_insertion_variants, bi, &vref_key);
                        vref.pstd = true;
                        vref.qstd = true;
                        if let Some(ref_var_mut) = ref_var.as_mut() {
                            adj_cnt_with_ref(vref, &sc3_var, Some(ref_var_mut));
                        } else {
                            adj_cnt(vref, &sc3_var);
                        }
                        adj_cnt(vref, &sc5_var);
                    }
                    if let (Some(ref_key), Some(ref_var)) = (ref_key, ref_var) {
                        if let Some(pos_map) = data.non_insertion_variants.get_mut(&bi) {
                            pos_map.insert(ref_key, ref_var);
                        }
                    }
                    self.realign_deletion_for_desc(bi, &vref_key, data);
                } else {
                    let vref_key = VarDesc::Raw {
                        desc: ins_desc.clone().into(),
                    };
                    let vref =
                        get_variants_from_map(&mut data.non_insertion_variants, bi, &vref_key);
                    vref.pstd = true;
                    vref.qstd = true;
                    adj_cnt(vref, &sc3_var);
                    adj_cnt(vref, &sc5_var);
                }

                break;
            }
        }
    }

    /// Java parity: realignlgdel (5' and 3' soft-clip paths).
    pub fn realign_large_deletions(&self, data: &mut RealignedVariationData) {
        let conf = &crate::scopedata::global_read_only_scope::instance().conf;
        let longmm = 3usize;
        let indel_size = 50i64;
        let base_region_start = data.ref_coverage.keys().min().copied().unwrap_or(1);
        let base_region_end = data
            .ref_coverage
            .keys()
            .max()
            .copied()
            .unwrap_or(base_region_start);

        let mut tmp5: Vec<(i64, usize)> = data
            .soft_clips_5end
            .iter()
            .map(|(pos, sc)| (*pos, sc.var.alt_depth))
            .collect();
        tmp5.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

        // Java initializes svcov once before each lgdel loop and reuses the
        // last computed value when a later candidate resolves via findbp.
        let mut svcov = 0usize;

        for (position, cnt) in tmp5 {
            if cnt < conf.minr {
                break;
            }

            if data
                .soft_clips_5end
                .get(&position)
                .map(|sc| sc.used())
                .unwrap_or(false)
            {
                continue;
            }

            let seq = {
                let Some(sc5v) = data.soft_clips_5end.get_mut(&position) else {
                    continue;
                };
                find_conseq(sc5v, 5)
            };
            if seq.is_empty() || seq.len() < 7 {
                continue;
            }

            let mut bp = self.find_bp(&seq, position - 5, -1);
            let mut matched_extra: Vec<u8> = Vec::new();

            if bp == 0 {
                if Self::is_low_complex_seq(&String::from_utf8_lossy(&seq)) {
                    continue;
                }

                let matched =
                    self.find_match(&seq, position, -1, Configuration::SEED_1 as usize, 1);
                bp = matched.base_position;
                matched_extra = matched.matched_sequence;
                if !(bp != 0
                    && position - bp > 15
                    && position - bp < Configuration::SVMAXLEN as i64)
                {
                    continue;
                }

                bp += 1;

                let (cov_f, clusters_f, pairs_f) = Self::mark_sv_clusters(
                    bp,
                    position,
                    &mut data.svfdel,
                    data.max_read_length as i64,
                );
                let (cov_r, clusters_r, pairs_r) = Self::mark_sv_clusters(
                    bp,
                    position,
                    &mut data.svrdel,
                    data.max_read_length as i64,
                );
                svcov = cov_f + cov_r;
                let clusters = clusters_f + clusters_r;
                let pairs = pairs_f + pairs_r;

                if svcov == 0 && cnt <= conf.minr {
                    continue;
                }
                Self::add_sv_counts(data, bp, pairs, cnt, clusters);

                if bp < base_region_start {
                    let tts = bp - data.max_read_length as i64;
                    let mut tte = bp + data.max_read_length as i64;
                    if bp + data.max_read_length as i64 >= base_region_start {
                        tte = base_region_start - 1;
                    }
                    self.load_partial_ref_coverage(data, tts, tte);
                }
            }

            let mut dellen = position - bp;
            let mut extra = Vec::<u8>::new();
            let gt = if matched_extra.is_empty() {
                let mut en = 0usize;
                while en < seq.len()
                    && self
                        .get_ref_base(bp - en as i64 - 1)
                        .map(|base| seq[en] != base)
                        .unwrap_or(false)
                {
                    extra.push(seq[en]);
                    en += 1;
                }

                if !extra.is_empty() {
                    extra.reverse();
                    let gt = format!("-{}&{}", dellen, String::from_utf8_lossy(&extra));
                    bp -= extra.len() as i64;
                    gt
                } else {
                    format!("-{}", dellen)
                }
            } else {
                dellen -= matched_extra.len() as i64;
                if dellen == 0 {
                    format!(
                        "-{}^{}",
                        matched_extra.len(),
                        String::from_utf8_lossy(&matched_extra)
                    )
                } else {
                    format!("-{}&{}", dellen, String::from_utf8_lossy(&matched_extra))
                }
            };

            if dellen < 0 {
                continue;
            }

            let mut n = 0i64;
            if extra.is_empty() && matched_extra.is_empty() {
                while self.get_ref_base(bp + n).is_some()
                    && self.get_ref_base(bp + dellen + n).is_some()
                    && self.get_ref_base(bp + n) == self.get_ref_base(bp + dellen + n)
                {
                    n += 1;
                }
            }

            let mut sc3p = bp + n;
            let mut mcnt = 0usize;
            let mut mismatch_len = 0usize;
            while mcnt <= longmm
                && self.get_ref_base(bp + n).is_some()
                && self.get_ref_base(bp + dellen + n).is_some()
                && self.get_ref_base(bp + n) != self.get_ref_base(bp + dellen + n)
            {
                n += 1;
                mcnt += 1;
                mismatch_len += 1;
            }

            if mismatch_len == 1 {
                let mut nm = 0usize;
                while self.get_ref_base(bp + n).is_some()
                    && self.get_ref_base(bp + dellen + n).is_some()
                    && self.get_ref_base(bp + n) == self.get_ref_base(bp + dellen + n)
                {
                    n += 1;
                    nm += 1;
                }
                if nm >= 3 && !data.soft_clips_3end.contains_key(&sc3p) {
                    sc3p = bp + n;
                }
            }

            if data
                .non_insertion_variants
                .get(&bp)
                .map(Self::has_sv_marker)
                .unwrap_or(false)
                && !data.soft_clips_3end.contains_key(&sc3p)
            {
                if svcov == 0 && cnt <= conf.minr {
                    Self::remove_sv_marker(data, bp);
                    continue;
                }
            }

            let key = VarDesc::Raw {
                desc: gt.as_bytes().to_vec().into(),
            };

            let sc5_contribution = data
                .soft_clips_5end
                .get(&position)
                .map(|sc| sc.var.clone())
                .unwrap_or_default();

            {
                let tv = get_variants_from_map(&mut data.non_insertion_variants, bp, &key);
                tv.qstd = true;
                tv.pstd = true;
                adj_cnt(tv, &sc5_contribution);
            }

            if let Some(sc5v) = data.soft_clips_5end.get_mut(&position) {
                sc5v.mark_used();
            }

            if !data.ref_coverage.contains_key(&bp) {
                if let Some(cov_p) = data.ref_coverage.get(&position).copied() {
                    data.ref_coverage.insert(bp, cov_p);
                }
            }

            if dellen > 0 && dellen < indel_size {
                for tp in bp..(bp + dellen) {
                    Self::emit_refcov_inc_diag(
                        data,
                        tp,
                        sc5_contribution.alt_depth,
                        "realigner_lgdel_5_sc5_span",
                    );
                    *data.ref_coverage.entry(tp).or_insert(0) += sc5_contribution.alt_depth;
                }
            }

            let use_sc3 = data
                .soft_clips_3end
                .get(&sc3p)
                .map(|sc| !sc.used())
                .unwrap_or(false);
            if use_sc3 {
                let sc3_contribution = data
                    .soft_clips_3end
                    .get(&sc3p)
                    .map(|sc| sc.var.clone())
                    .unwrap_or_default();

                if sc3p > bp {
                    let ref_key = self
                        .get_ref_base(bp)
                        .map(|ref_base| VarDesc::SNV { ref_base });
                    let mut ref_var = ref_key.as_ref().and_then(|k| {
                        data.non_insertion_variants
                            .get_mut(&bp)
                            .and_then(|m| m.remove(k))
                    });

                    {
                        let tv = get_variants_from_map(&mut data.non_insertion_variants, bp, &key);
                        if let Some(ref_var_mut) = ref_var.as_mut() {
                            adj_cnt_with_ref(tv, &sc3_contribution, Some(ref_var_mut));
                        } else {
                            adj_cnt(tv, &sc3_contribution);
                        }
                    }

                    if let (Some(ref_key), Some(ref_var)) = (ref_key, ref_var) {
                        if let Some(pos_map) = data.non_insertion_variants.get_mut(&bp) {
                            pos_map.insert(ref_key, ref_var);
                        }
                    }
                } else {
                    let tv = get_variants_from_map(&mut data.non_insertion_variants, bp, &key);
                    adj_cnt(tv, &sc3_contribution);
                }

                if sc3p == bp && dellen > 0 && dellen < indel_size {
                    for tp in bp..(bp + dellen) {
                        Self::emit_refcov_inc_diag(
                            data,
                            tp,
                            sc3_contribution.alt_depth,
                            "realigner_lgdel_5_sc3_span",
                        );
                        *data.ref_coverage.entry(tp).or_insert(0) += sc3_contribution.alt_depth;
                    }
                }

                if sc3p > bp && dellen > 0 {
                    for ip in (bp + 1)..sc3p {
                        let Some(ref_base) = self.get_ref_base(dellen + ip) else {
                            continue;
                        };
                        let ref_key = VarDesc::SNV { ref_base };
                        let mut remove_pos = false;
                        if let Some(map) = data.non_insertion_variants.get_mut(&ip) {
                            if let Some(vv) = map.get_mut(&ref_key) {
                                Self::rm_cnt(vv, &sc3_contribution);
                                if vv.alt_depth == 0 {
                                    map.remove(&ref_key);
                                }
                            }
                            if map.is_empty() {
                                remove_pos = true;
                            }
                        }
                        if remove_pos {
                            data.non_insertion_variants.remove(&ip);
                        }
                    }
                }

                if let Some(sc3v) = data.soft_clips_3end.get_mut(&sc3p) {
                    sc3v.mark_used();
                }
            }

            let tv_before_realign = data
                .non_insertion_variants
                .get(&bp)
                .and_then(|m| m.get(&key))
                .map(|v| v.alt_depth)
                .unwrap_or(0);

            let mut dels: CountMapByPos = Default::default();
            let mut inner: CountMap = Default::default();
            inner.insert(gt.clone(), tv_before_realign);
            dels.insert(bp, inner);
            self.process_deletions_with_passing_check(data, &dels, true);

            let tv_after_realign = data
                .non_insertion_variants
                .get(&bp)
                .and_then(|m| m.get(&key))
                .map(|v| v.alt_depth)
                .unwrap_or(tv_before_realign);
            let sv_delta = tv_after_realign as isize - tv_before_realign as isize;
            Self::adjust_sv_split_count(data, bp, sv_delta);

            if let Some(tv) = data
                .non_insertion_variants
                .get_mut(&bp)
                .and_then(|m| m.get_mut(&key))
            {
                if svcov > tv.alt_depth && tv.alt_depth > 0 {
                    Self::add_var_factor(tv, (svcov - tv.alt_depth) as f64 / tv.alt_depth as f64);
                }
            }
        }

        let mut tmp3: Vec<(i64, usize)> = data
            .soft_clips_3end
            .iter()
            .map(|(pos, sc)| (*pos, sc.var.alt_depth))
            .collect();
        tmp3.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

        let mut svcov = 0usize;

        for (position, cnt) in tmp3 {
            if cnt < conf.minr {
                break;
            }

            if data
                .soft_clips_3end
                .get(&position)
                .map(|sc| sc.used())
                .unwrap_or(false)
            {
                continue;
            }

            let seq = {
                let Some(sc3v) = data.soft_clips_3end.get_mut(&position) else {
                    continue;
                };
                find_conseq(sc3v, 3)
            };
            if seq.is_empty() || seq.len() < 7 {
                continue;
            }

            let mut breakpoint = self.find_bp(&seq, position + 5, 1);
            let mut matched_extra: Vec<u8> = Vec::new();

            if breakpoint == 0 {
                if Self::is_low_complex_seq(&String::from_utf8_lossy(&seq)) {
                    continue;
                }

                let matched = self.find_match(&seq, position, 1, Configuration::SEED_1 as usize, 1);
                breakpoint = matched.base_position;
                matched_extra = matched.matched_sequence;
                if !(breakpoint != 0
                    && breakpoint - position > 15
                    && position - breakpoint < Configuration::SVMAXLEN as i64)
                {
                    continue;
                }

                let (cov_f, clusters_f, pairs_f) = Self::mark_sv_clusters(
                    position,
                    breakpoint,
                    &mut data.svfdel,
                    data.max_read_length as i64,
                );
                let (cov_r, clusters_r, pairs_r) = Self::mark_sv_clusters(
                    position,
                    breakpoint,
                    &mut data.svrdel,
                    data.max_read_length as i64,
                );
                svcov = cov_f + cov_r;
                let clusters = clusters_f + clusters_r;
                let pairs = pairs_f + pairs_r;

                if svcov == 0 && cnt <= conf.minr {
                    continue;
                }
                Self::add_sv_counts(data, position, pairs, cnt, clusters);

                if breakpoint > base_region_end {
                    let mut tts = breakpoint - data.max_read_length as i64;
                    let tte = breakpoint + data.max_read_length as i64;
                    if breakpoint - data.max_read_length as i64 <= base_region_end {
                        tts = base_region_end + 1;
                    }
                    self.load_partial_ref_coverage(data, tts, tte);
                }
            }

            let mut deleted_len = breakpoint - position;
            let mut extra = Vec::<u8>::new();

            if !matched_extra.is_empty() {
                deleted_len -= matched_extra.len() as i64;
            } else {
                let mut en = 0usize;
                while en < seq.len()
                    && self
                        .get_ref_base(breakpoint + en as i64)
                        .map(|base| base != seq[en])
                        .unwrap_or(false)
                {
                    extra.push(seq[en]);
                    en += 1;
                }
            }

            if deleted_len <= 0 {
                continue;
            }

            let mut gt = format!("-{}", deleted_len);
            let mut anchor_pos = position;
            let mut sc5_pos = breakpoint;

            if !extra.is_empty() {
                gt = format!("-{}&{}", deleted_len, String::from_utf8_lossy(&extra));
                sc5_pos += extra.len() as i64;
            } else if !matched_extra.is_empty() {
                gt = format!(
                    "-{}&{}",
                    deleted_len,
                    String::from_utf8_lossy(&matched_extra)
                );
            } else {
                while self.get_ref_base(anchor_pos - 1).is_some()
                    && self.get_ref_base(anchor_pos + deleted_len - 1).is_some()
                    && self.get_ref_base(anchor_pos - 1)
                        == self.get_ref_base(anchor_pos + deleted_len - 1)
                {
                    anchor_pos -= 1;
                    if anchor_pos != 0 {
                        sc5_pos -= 1;
                    }
                }

                if anchor_pos != position {
                    Self::move_sv_marker(data, position, anchor_pos);
                }
            }

            if data
                .non_insertion_variants
                .get(&anchor_pos)
                .map(Self::has_sv_marker)
                .unwrap_or(false)
                && !data.soft_clips_5end.contains_key(&sc5_pos)
            {
                if svcov == 0 && cnt <= conf.minr {
                    Self::remove_sv_marker(data, anchor_pos);
                    continue;
                }
            }

            if deleted_len >= indel_size && !data.ref_coverage.contains_key(&anchor_pos) {
                if let Some(cov_prev) = data.ref_coverage.get(&(position - 1)).copied() {
                    data.ref_coverage.insert(anchor_pos, cov_prev);
                } else {
                    data.ref_coverage.insert(anchor_pos, cnt);
                }
            }

            if deleted_len < indel_size {
                let span = deleted_len + extra.len() as i64 + matched_extra.len() as i64;
                for tp in anchor_pos..(anchor_pos + span) {
                    Self::emit_refcov_inc_diag(
                        data,
                        tp,
                        cnt,
                        "realigner_lgdel_3_span",
                    );
                    *data.ref_coverage.entry(tp).or_insert(0) += cnt;
                }
            }

            let mut sc3_contribution = match data.soft_clips_3end.get(&position) {
                Some(sc3v) => sc3v.var.clone(),
                None => continue,
            };
            sc3_contribution.mean_pos += deleted_len as f64 * sc3_contribution.alt_depth as f64;

            let key = VarDesc::Raw {
                desc: gt.as_bytes().to_vec().into(),
            };

            {
                let tv = get_variants_from_map(&mut data.non_insertion_variants, anchor_pos, &key);
                tv.qstd = true;
                tv.pstd = true;
                adj_cnt(tv, &sc3_contribution);
            }

            if let Some(sc3v) = data.soft_clips_3end.get_mut(&position) {
                sc3v.var.mean_pos += deleted_len as f64 * sc3v.var.alt_depth as f64;
                sc3v.mark_used();
            }

            let tv_before_realign = data
                .non_insertion_variants
                .get(&anchor_pos)
                .and_then(|m| m.get(&key))
                .map(|v| v.alt_depth)
                .unwrap_or(0);

            let mut dels: CountMapByPos = Default::default();
            let mut inner: CountMap = Default::default();
            inner.insert(gt.clone(), tv_before_realign);
            dels.insert(anchor_pos, inner);
            self.process_deletions_with_passing_check(data, &dels, true);

            let tv_after_realign = data
                .non_insertion_variants
                .get(&anchor_pos)
                .and_then(|m| m.get(&key))
                .map(|v| v.alt_depth)
                .unwrap_or(tv_before_realign);
            let sv_delta = tv_after_realign as isize - tv_before_realign as isize;
            Self::adjust_sv_split_count(data, anchor_pos, sv_delta);

            if let Some(tv) = data
                .non_insertion_variants
                .get_mut(&anchor_pos)
                .and_then(|m| m.get_mut(&key))
            {
                if svcov > tv.alt_depth && tv.alt_depth > 0 {
                    Self::add_var_factor(tv, (svcov - tv.alt_depth) as f64 / tv.alt_depth as f64);
                }
            }
        }
    }

    pub fn realign_long_insertions(&self, data: &mut RealignedVariationData) {
        let conf = &crate::scopedata::global_read_only_scope::instance().conf;
        let base_region_start = data.ref_coverage.keys().min().copied().unwrap_or(1);
        let base_region_end = data
            .ref_coverage
            .keys()
            .max()
            .copied()
            .unwrap_or(base_region_start);

        let mut tmp: Vec<(i64, usize)> = data
            .soft_clips_5end
            .iter()
            .map(|(pos, sc)| (*pos, sc.var.alt_depth))
            .collect();
        tmp.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

        for (p, cnt) in tmp {
            if cnt < conf.minr {
                break;
            }
            if data
                .soft_clips_5end
                .get(&p)
                .map(|sc| sc.used())
                .unwrap_or(false)
            {
                continue;
            }

            let seq = {
                let Some(sc5v) = data.soft_clips_5end.get_mut(&p) else {
                    continue;
                };
                find_conseq_transient(sc5v, 0)
            };
            if seq.is_empty() || seq.len() < 12 {
                continue;
            }

            let tpl = self.find_bi(&seq, p, -1);
            let mut bi = tpl.base_insert;
            let mut ins = tpl.insertion_sequence;

            if bi == 0 {
                if Self::is_low_complex_seq(&String::from_utf8_lossy(&seq)) {
                    continue;
                }
                let match_res = self.find_match(&seq, p, -1, Configuration::SEED_1 as usize, 1);
                bi = match_res.base_position;
                let extra = match_res.matched_sequence;
                if !(bi != 0 && bi - p > 15 && bi - p < Configuration::SVMAXLEN as i64) {
                    continue;
                }

                if bi > base_region_end {
                    let mut tts = bi - data.max_read_length as i64;
                    let tte = bi + data.max_read_length as i64;
                    if bi - data.max_read_length as i64 <= base_region_end {
                        tts = base_region_end + 1;
                    }
                    self.load_partial_ref_coverage(data, tts, tte);
                }

                if bi - p > conf.sv_min_len as i64 + 2 * Configuration::SVFLANK as i64 {
                    ins = self.join_ref(p, p + Configuration::SVFLANK as i64 - 1);
                    let dup_len = bi - p - 2 * Configuration::SVFLANK as i64 + 1;
                    ins.extend_from_slice(format!("<dup{}>", dup_len).as_bytes());
                    let tail = self.join_ref_for_5_lgins(
                        bi - Configuration::SVFLANK as i64 + 1,
                        bi,
                        &seq,
                        &extra,
                    );
                    ins.extend_from_slice(&tail);
                } else {
                    ins = self.join_ref_for_5_lgins(p, bi, &seq, &extra);
                }
                ins.extend_from_slice(&extra);

                let ref_cov_bi = data.ref_coverage.get(&bi).copied();
                let ref_cov_p1 = data.ref_coverage.get(&(p - 1)).copied();
                if ref_cov_p1.is_none()
                    || (ref_cov_bi.is_some()
                        && ref_cov_p1.is_some()
                        && ref_cov_p1.unwrap() < ref_cov_bi.unwrap())
                {
                    if let Some(val) = ref_cov_bi {
                        data.ref_coverage.insert(p - 1, val);
                    } else {
                        data.ref_coverage.insert(p - 1, cnt);
                    }
                } else if cnt > ref_cov_p1.unwrap_or(0) {
                    Self::emit_refcov_inc_diag(
                        data,
                        p - 1,
                        cnt,
                        "realigner_lgins_5_prev_anchor",
                    );
                    *data.ref_coverage.entry(p - 1).or_insert(0) += cnt;
                }

                let (clusters_f, pairs_f) =
                    Self::mark_dup_sv(p, bi, &mut data.svfdup, data.max_read_length as i64);
                let (clusters_r, pairs_r) =
                    Self::mark_dup_sv(p, bi, &mut data.svrdup, data.max_read_length as i64);
                let clusters = clusters_f + clusters_r;
                let pairs = pairs_f + pairs_r;

                bi = p - 1;
                Self::add_sv_counts(data, bi, pairs, cnt, clusters);
            }

            let sc5_var = data
                .soft_clips_5end
                .get(&p)
                .map(|sc| sc.var.clone())
                .unwrap_or_default();
            let sc5_seq = data
                .soft_clips_5end
                .get(&p)
                .map(|sc| sc.seq.clone())
                .unwrap_or_default();

            let iref_key = VarDesc::Ins {
                seq: ins.iter().copied().collect(),
            };
            let iref = get_variants_from_map(&mut data.insertion_variants, bi, &iref_key);
            iref.pstd = true;
            iref.qstd = true;
            adj_cnt(iref, &sc5_var);
            let original_ins_count = iref.alt_depth;

            if let Some(variation_map) = data.non_insertion_variants.get(&bi) {
                if !Self::has_sv_marker(variation_map) {
                    Self::emit_refcov_inc_diag(
                        data,
                        bi,
                        sc5_var.alt_depth,
                        "realigner_lgins_5_anchor",
                    );
                    *data.ref_coverage.entry(bi).or_insert(0) += sc5_var.alt_depth;
                }
            }

            let mut len = ins.len();
            if ins.iter().any(|b| *b == b'&') {
                len = len.saturating_sub(1);
            }
            let seq_len = sc5_seq.keys().next_back().map(|k| *k + 1).unwrap_or(0);

            for ii in (len + 1)..seq_len {
                let pii = bi - ii as i64 + len as i64;
                let Some(map) = sc5_seq.get(&ii) else {
                    continue;
                };
                for base in [b'A', b'C', b'G', b'T', b'N'] {
                    let Some(tv) = map.get(base) else {
                        continue;
                    };
                    let key = VarDesc::SNV { ref_base: base };
                    let tvr = get_variants_from_map(&mut data.non_insertion_variants, pii, &key);
                    adj_cnt(tvr, tv);
                    tvr.pstd = true;
                    tvr.qstd = true;
                    Self::emit_refcov_inc_diag(
                        data,
                        pii,
                        tv.alt_depth,
                        "realigner_lgins_5_tail_snv",
                    );
                    *data.ref_coverage.entry(pii).or_insert(0) += tv.alt_depth;
                }
            }

            if bi + len as i64 != 0 {
                if let Some(sc5v) = data.soft_clips_5end.get_mut(&p) {
                    cache_transient_conseq(sc5v, &seq);
                    sc5v.mark_used();
                }
            }

            let mut tins: CountMapByPos = Default::default();
            let mut map: CountMap = Default::default();
            map.insert(
                format!("+{}", String::from_utf8_lossy(&ins)),
                original_ins_count,
            );
            tins.insert(bi, map);
            let before_insertions = data.insertion_variants.get(&bi).cloned();
            self.process_insertions(data, &tins);
            self.emit_realigner_insertion_diag(data, bi, "realign_large_deletions_5end");

            let updated_ins_count = Self::resolved_insertion_key_after_realign(
                before_insertions.as_ref(),
                data.insertion_variants.get(&bi),
                &iref_key,
            )
            .and_then(|resolved_key| {
                data.insertion_variants
                    .get(&bi)
                    .and_then(|m| m.get(&resolved_key))
                    .map(|v| v.alt_depth)
            })
            .unwrap_or(original_ins_count);
            let sv_delta = updated_ins_count as isize - original_ins_count as isize;
            Self::adjust_sv_split_count(data, bi, sv_delta);

            let rpflag = ins
                .iter()
                .enumerate()
                .all(|(idx, base)| self.get_ref_base(bi + 1 + idx as i64) == Some(*base));
            if rpflag
                && !self.bam_paths.is_empty()
                && ins.len() >= 5
                && ins.len() < data.max_read_length.saturating_sub(10)
            {
                let resolved_key = Self::resolved_insertion_key_after_realign(
                    before_insertions.as_ref(),
                    data.insertion_variants.get(&bi),
                    &iref_key,
                );
                let ref_key = self
                    .get_ref_base(bi)
                    .map(|ref_base| VarDesc::SNV { ref_base });

                if let (Some(resolved_key), Some(ref_key)) = (resolved_key, ref_key) {
                    let mref_snapshot = data
                        .non_insertion_variants
                        .get(&bi)
                        .and_then(|m| m.get(&ref_key))
                        .cloned();
                    let kref_count = data
                        .insertion_variants
                        .get(&bi)
                        .and_then(|m| m.get(&resolved_key))
                        .map(|v| v.alt_depth)
                        .unwrap_or(0);

                    if let Some(mref) = mref_snapshot {
                        let no_passing = self.no_passing_reads(bi, bi + ins.len() as i64);
                        if mref.alt_depth != 0 && no_passing && kref_count > 2 * mref.alt_depth {
                            if let Some(mut ref_var) = data
                                .non_insertion_variants
                                .get_mut(&bi)
                                .and_then(|m| m.remove(&ref_key))
                            {
                                let ref_src = ref_var.clone();
                                if let Some(kref) = data
                                    .insertion_variants
                                    .get_mut(&bi)
                                    .and_then(|m| m.get_mut(&resolved_key))
                                {
                                    adj_cnt_with_ref(kref, &ref_src, Some(&mut ref_var));
                                }
                                if let Some(pos_map) = data.non_insertion_variants.get_mut(&bi) {
                                    pos_map.insert(ref_key, ref_var);
                                }
                            }
                        }
                    }
                }
            }
        }

        let mut tmp3: Vec<(i64, usize)> = data
            .soft_clips_3end
            .iter()
            .map(|(pos, sc)| (*pos, sc.var.alt_depth))
            .collect();
        tmp3.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

        for (position, cnt) in tmp3 {
            if cnt < conf.minr {
                break;
            }
            if data
                .soft_clips_3end
                .get(&position)
                .map(|sc| sc.used())
                .unwrap_or(false)
            {
                continue;
            }

            let seq = {
                let Some(sc3v) = data.soft_clips_3end.get_mut(&position) else {
                    continue;
                };
                find_conseq_transient(sc3v, 0)
            };
            if seq.is_empty() || seq.len() < 12 {
                continue;
            }

            let mut p = position;
            let tpl = self.find_bi(&seq, p, 1);
            let mut bi = tpl.base_insert;
            let mut ins = tpl.insertion_sequence;

            if bi == 0 {
                if Self::is_low_complex_seq(&String::from_utf8_lossy(&seq)) {
                    continue;
                }
                let match_res = self.find_match(&seq, p, 1, Configuration::SEED_1 as usize, 1);
                bi = match_res.base_position;
                let extra = match_res.matched_sequence;
                if !(bi != 0 && p - bi > 15 && p - bi < Configuration::SVMAXLEN as i64) {
                    continue;
                }

                if bi < base_region_start {
                    let tts = bi - data.max_read_length as i64;
                    let mut tte = bi + data.max_read_length as i64;
                    if bi + data.max_read_length as i64 >= base_region_start {
                        tte = base_region_start - 1;
                    }
                    self.load_partial_ref_coverage(data, tts, tte);
                }

                let mut shift5 = 0i64;
                while self.get_ref_base(p - 1).is_some()
                    && self.get_ref_base(bi - 1).is_some()
                    && self.get_ref_base(p - 1) == self.get_ref_base(bi - 1)
                {
                    p -= 1;
                    bi -= 1;
                    shift5 += 1;
                }

                if p - bi > conf.sv_min_len as i64 + 2 * Configuration::SVFLANK as i64 {
                    ins = self.join_ref_for_3_lgins(
                        bi,
                        bi + Configuration::SVFLANK as i64 - 1,
                        shift5,
                        &seq,
                        &extra,
                    );
                    let dup_len = p - bi - 2 * Configuration::SVFLANK as i64;
                    ins.extend_from_slice(format!("<dup{}>", dup_len).as_bytes());
                    let tail = self.join_ref(p - Configuration::SVFLANK as i64, p - 1);
                    ins.extend_from_slice(&tail);
                } else {
                    ins = self.join_ref_for_3_lgins(bi, p - 1, shift5, &seq, &extra);
                }
                ins.extend_from_slice(&extra);

                let (clusters_f, pairs_f) =
                    Self::mark_dup_sv(bi, p - 1, &mut data.svfdup, data.max_read_length as i64);
                let (clusters_r, pairs_r) =
                    Self::mark_dup_sv(bi, p - 1, &mut data.svrdup, data.max_read_length as i64);
                let clusters = clusters_f + clusters_r;
                let pairs = pairs_f + pairs_r;

                bi -= 1;
                Self::add_sv_counts(data, bi, pairs, cnt, clusters);

                let ref_cov_p = data.ref_coverage.get(&p).copied();
                let ref_cov_bi = data.ref_coverage.get(&bi).copied();
                if ref_cov_bi.is_none()
                    || (ref_cov_p.is_some()
                        && ref_cov_bi.is_some()
                        && ref_cov_bi.unwrap() < ref_cov_p.unwrap())
                {
                    if let Some(val) = ref_cov_p {
                        data.ref_coverage.insert(bi, val);
                    } else {
                        data.ref_coverage.insert(bi, cnt);
                    }
                } else if cnt > ref_cov_bi.unwrap_or(0) {
                    Self::emit_refcov_inc_diag(
                        data,
                        bi,
                        cnt,
                        "realigner_lgins_3_anchor",
                    );
                    *data.ref_coverage.entry(bi).or_insert(0) += cnt;
                }
            }

            let sc3_var = data
                .soft_clips_3end
                .get(&position)
                .map(|sc| sc.var.clone())
                .unwrap_or_default();
            let sc3_seq = data
                .soft_clips_3end
                .get(&position)
                .map(|sc| sc.seq.clone())
                .unwrap_or_default();

            let iref_key = VarDesc::Ins {
                seq: ins.iter().copied().collect(),
            };
            let iref = get_variants_from_map(&mut data.insertion_variants, bi, &iref_key);
            iref.pstd = true;
            iref.qstd = true;

            let mean_pos = if cnt > 0 {
                sc3_var.mean_pos / cnt as f64
            } else {
                0.0
            };
            let mut ref_key = self
                .get_ref_base(bi)
                .map(|ref_base| VarDesc::SNV { ref_base });
            if (p - bi) as f64 > mean_pos {
                ref_key = None;
            }

            let mut ref_var = ref_key.as_ref().and_then(|k| {
                data.non_insertion_variants
                    .get_mut(&bi)
                    .and_then(|m| m.remove(k))
            });
            if let Some(ref_var_mut) = ref_var.as_mut() {
                adj_cnt_with_ref(iref, &sc3_var, Some(ref_var_mut));
            } else {
                adj_cnt(iref, &sc3_var);
            }

            if let (Some(ref_key), Some(ref_var)) = (ref_key, ref_var) {
                if let Some(pos_map) = data.non_insertion_variants.get_mut(&bi) {
                    pos_map.insert(ref_key, ref_var);
                }
            }

            let original_ins_count = iref.alt_depth;

            let mut len = ins.len();
            if ins.iter().any(|b| *b == b'&') {
                len = len.saturating_sub(1);
            }
            let seq_len = sc3_seq.keys().next_back().map(|k| *k + 1).unwrap_or(0);

            for ii in len..seq_len {
                let pii = p + ii as i64 - len as i64;
                let Some(map) = sc3_seq.get(&ii) else {
                    continue;
                };
                for base in [b'A', b'C', b'G', b'T', b'N'] {
                    let Some(tv) = map.get(base) else {
                        continue;
                    };
                    let key = VarDesc::SNV { ref_base: base };
                    Self::emit_refcov_inc_diag(
                        data,
                        pii,
                        tv.alt_depth,
                        "realigner_lgins_3_tail_snv",
                    );
                    let vref = get_variants_from_map(&mut data.non_insertion_variants, pii, &key);
                    adj_cnt(vref, tv);
                    vref.pstd = true;
                    vref.qstd = true;
                    *data.ref_coverage.entry(pii).or_insert(0) += tv.alt_depth;
                }
            }

            if let Some(sc3v) = data.soft_clips_3end.get_mut(&position) {
                sc3v.mark_used();
            }

            let mut tins: CountMapByPos = Default::default();
            let mut map: CountMap = Default::default();
            map.insert(
                format!("+{}", String::from_utf8_lossy(&ins)),
                original_ins_count,
            );
            tins.insert(bi, map);
            let before_insertions = data.insertion_variants.get(&bi).cloned();
            self.process_insertions(data, &tins);
            self.emit_realigner_insertion_diag(data, bi, "realign_large_deletions_3end");

            let updated_ins_count = Self::resolved_insertion_key_after_realign(
                before_insertions.as_ref(),
                data.insertion_variants.get(&bi),
                &iref_key,
            )
            .and_then(|resolved_key| {
                data.insertion_variants
                    .get(&bi)
                    .and_then(|m| m.get(&resolved_key))
                    .map(|v| v.alt_depth)
            })
            .unwrap_or(original_ins_count);
            let sv_delta = updated_ins_count as isize - original_ins_count as isize;
            Self::adjust_sv_split_count(data, bi, sv_delta);

            let rpflag = ins
                .iter()
                .enumerate()
                .all(|(idx, base)| self.get_ref_base(bi + 1 + idx as i64) == Some(*base));
            if rpflag
                && !self.bam_paths.is_empty()
                && ins.len() >= 5
                && ins.len() < data.max_read_length.saturating_sub(10)
            {
                let resolved_key = Self::resolved_insertion_key_after_realign(
                    before_insertions.as_ref(),
                    data.insertion_variants.get(&bi),
                    &iref_key,
                );
                let ref_key = self
                    .get_ref_base(bi)
                    .map(|ref_base| VarDesc::SNV { ref_base });

                if let (Some(resolved_key), Some(ref_key)) = (resolved_key, ref_key) {
                    let mref_snapshot = data
                        .non_insertion_variants
                        .get(&bi)
                        .and_then(|m| m.get(&ref_key))
                        .cloned();
                    let kref_count = data
                        .insertion_variants
                        .get(&bi)
                        .and_then(|m| m.get(&resolved_key))
                        .map(|v| v.alt_depth)
                        .unwrap_or(0);

                    if let Some(mref) = mref_snapshot {
                        let no_passing = self.no_passing_reads(bi, bi + ins.len() as i64);
                        if mref.alt_depth != 0 && no_passing && kref_count > 2 * mref.alt_depth {
                            if let Some(mut ref_var) = data
                                .non_insertion_variants
                                .get_mut(&bi)
                                .and_then(|m| m.remove(&ref_key))
                            {
                                let ref_src = ref_var.clone();
                                if let Some(kref) = data
                                    .insertion_variants
                                    .get_mut(&bi)
                                    .and_then(|m| m.get_mut(&resolved_key))
                                {
                                    adj_cnt_with_ref(kref, &ref_src, Some(&mut ref_var));
                                }
                                if let Some(pos_map) = data.non_insertion_variants.get_mut(&bi) {
                                    pos_map.insert(ref_key, ref_var);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /// Adjust MNPs (multi-nucleotide polymorphisms) when there are breakpoints within MNP
    pub fn adjust_mnp(
        &self,
        data: &mut RealignedVariationData,
        mnp: &CountMapByPos,
    ) {
        let mut tmp: Vec<(i64, String, usize)> = Vec::new();
        for (pos, desc_map) in mnp {
            for (desc, count) in desc_map {
                tmp.push((*pos, desc.clone(), *count));
            }
        }
        tmp.sort_by(|a, b| {
            b.2.cmp(&a.2)
                .then_with(|| a.0.cmp(&b.0))
                .then_with(|| b.1.cmp(&a.1))
        });

        for (position, vn, _) in tmp {
            let vref_key = VarDesc::Raw {
                desc: vn.as_bytes().to_vec().into(),
            };
            if data
                .non_insertion_variants
                .get(&position)
                .and_then(|m| m.get(&vref_key))
                .is_none()
            {
                continue;
            }

            let mnt = vn.replacen('&', "", 1);
            if mnt.len() < 2 {
                continue;
            }

            for i in 0..(mnt.len() - 1) {
                let left = &mnt[..i + 1];
                let left_desc = if left.len() > 1 {
                    format!("{}&{}", &left[..1], &left[1..])
                } else {
                    left.to_string()
                };
                let right = &mnt[i + 1..];
                let right_desc = if right.len() > 1 {
                    format!("{}&{}", &right[..1], &right[1..])
                } else {
                    right.to_string()
                };

                let left_key = if left_desc.len() == 1 {
                    VarDesc::SNV {
                        ref_base: left_desc.as_bytes()[0],
                    }
                } else {
                    VarDesc::Raw {
                        desc: left_desc.as_bytes().to_vec().into(),
                    }
                };

                if let Some(tref) = data
                    .non_insertion_variants
                    .get(&position)
                    .and_then(|m| m.get(&left_key))
                    .cloned()
                {
                    let current_vref_cnt = data
                        .non_insertion_variants
                        .get(&position)
                        .and_then(|m| m.get(&vref_key))
                        .map(|v| v.alt_depth)
                        .unwrap_or(0);
                    if tref.alt_depth > 0
                        && tref.alt_depth < current_vref_cnt
                        && tref.mean_pos / tref.alt_depth as f64 <= (i + 1) as f64
                    {
                        if let Some(vars_on_pos) = data.non_insertion_variants.get_mut(&position) {
                            if let Some(vref) = vars_on_pos.get_mut(&vref_key) {
                                adj_cnt(vref, &tref);
                            }
                            vars_on_pos.remove(&left_key);
                        }
                    }
                }

                let right_pos = position + i as i64 + 1;
                let right_key = if right_desc.len() == 1 {
                    VarDesc::SNV {
                        ref_base: right_desc.as_bytes()[0],
                    }
                } else {
                    VarDesc::Raw {
                        desc: right_desc.as_bytes().to_vec().into(),
                    }
                };

                if let Some(tref) = data
                    .non_insertion_variants
                    .get(&right_pos)
                    .and_then(|m| m.get(&right_key))
                    .cloned()
                {
                    let current_vref_cnt = data
                        .non_insertion_variants
                        .get(&position)
                        .and_then(|m| m.get(&vref_key))
                        .map(|v| v.alt_depth)
                        .unwrap_or(0);
                    if tref.alt_depth < current_vref_cnt {
                        if let Some(vars_on_pos) = data.non_insertion_variants.get_mut(&position) {
                            if let Some(vref) = vars_on_pos.get_mut(&vref_key) {
                                adj_cnt(vref, &tref);
                            }
                        }
                        Self::emit_refcov_inc_diag(
                            data,
                            position,
                            tref.alt_depth,
                            "realigner_adjust_mnp_right",
                        );
                        *data.ref_coverage.entry(position).or_insert(0) += tref.alt_depth;

                        if let Some(vars_right) = data.non_insertion_variants.get_mut(&right_pos) {
                            vars_right.remove(&right_key);
                        }
                    }
                }
            }

            let refcov_before_sc3 = if position == 6970385 {
                data.ref_coverage.get(&6970385_i64).copied().unwrap_or(0)
            } else {
                0
            };
            if let Some(sc3v) = data.soft_clips_3end.get_mut(&position) {
                if !sc3v.used() {
                    let seq = find_conseq_transient(sc3v, 0);
                    if seq.starts_with(mnt.as_bytes()) {
                        let tail = &seq[mnt.len()..];
                        if tail.is_empty()
                            || self.is_match_ref(tail, position + mnt.len() as i64, 1, 3)
                        {
                            cache_transient_conseq(sc3v, &seq);
                            if let Some(vars_on_pos) =
                                data.non_insertion_variants.get_mut(&position)
                            {
                                if let Some(vref) = vars_on_pos.get_mut(&vref_key) {
                                    adj_cnt(vref, &sc3v.var);
                                }
                            }
                            Self::emit_refcov_inc_diag_with_before(
                                refcov_before_sc3,
                                position,
                                sc3v.var.alt_depth,
                                "realigner_adjust_mnp_sc3",
                            );
                            *data.ref_coverage.entry(position).or_insert(0) += sc3v.var.alt_depth;
                            sc3v.mark_used();
                        }
                    }
                }
            }

            let pos_5end = position + mnt.len() as i64;
            let refcov_before_sc5 = if position == 6970385 {
                data.ref_coverage.get(&6970385_i64).copied().unwrap_or(0)
            } else {
                0
            };
            if let Some(sc5v) = data.soft_clips_5end.get_mut(&pos_5end) {
                if !sc5v.used() {
                    let seq = find_conseq_transient(sc5v, 0);
                    if !seq.is_empty() && seq.len() >= mnt.len() {
                        let mut reversed_seq = seq.clone();
                        reversed_seq.reverse();
                        if reversed_seq.ends_with(mnt.as_bytes()) {
                            let prefix = &reversed_seq[..reversed_seq.len() - mnt.len()];
                            if prefix.is_empty() || self.is_match_ref(prefix, position - 1, -1, 3) {
                                cache_transient_conseq(sc5v, &seq);
                                if let Some(vars_on_pos) =
                                    data.non_insertion_variants.get_mut(&position)
                                {
                                    if let Some(vref) = vars_on_pos.get_mut(&vref_key) {
                                        adj_cnt(vref, &sc5v.var);
                                    }
                                }
                                Self::emit_refcov_inc_diag_with_before(
                                    refcov_before_sc5,
                                    position,
                                    sc5v.var.alt_depth,
                                    "realigner_adjust_mnp_sc5",
                                );
                                *data.ref_coverage.entry(position).or_insert(0) +=
                                    sc5v.var.alt_depth;
                                sc5v.mark_used();
                            }
                        }
                    }
                }
            }
        }
    }

    fn realign_deletion_mismatches(
        &self,
        pos: i64,
        desc: &VarDesc,
        fixed_dcnt: usize,
        data: &mut RealignedVariationData,
        enable_passing_reads_check: bool,
    ) {
        let desc_str = desc.to_key_string();
        let desc_bytes = desc_str.as_bytes();
        if desc_bytes.first() != Some(&b'-') {
            return;
        }

        let mut idx = 1usize;
        while idx < desc_bytes.len() && desc_bytes[idx].is_ascii_digit() {
            idx += 1;
        }
        if idx == 1 {
            return;
        }

        let mut dellen = std::str::from_utf8(&desc_bytes[1..idx])
            .ok()
            .and_then(|value| value.parse::<i64>().ok())
            .unwrap_or(0);
        if let Some(cap) = UP_NUMBER_END.captures(&desc_str) {
            if let Some(extra) = cap
                .get(1)
                .and_then(|group| group.as_str().parse::<i64>().ok())
            {
                dellen += extra;
            }
        }

        let mut extra_seq: String = desc_str[idx..]
            .chars()
            .filter(|ch| !matches!(*ch, '^' | '&' | '#'))
            .collect();
        let mut extrains_len = CARET_ATGNC
            .captures(&desc_str)
            .and_then(|cap| cap.get(1))
            .map(|group| group.as_str().len() as i64)
            .unwrap_or(0);

        let inv_flanks = Self::extract_inv_flanks(&desc_str);
        if inv_flanks.is_some() {
            extra_seq.clear();
            extrains_len = 0;
        }

        let mut wupseq = self.get_ref_range(pos - 200, pos - 1);
        wupseq.extend_from_slice(extra_seq.as_bytes());
        if let Some((_, inv3)) = inv_flanks.as_ref() {
            wupseq = inv3.as_bytes().to_vec();
        }

        let san_start = pos + dellen + extra_seq.len() as i64 - extrains_len;
        let ref_end = self.ref_start + self.reference_seq.len() as i64 - 1;
        let mut san_end = pos + 200;
        if san_end > ref_end {
            san_end = ref_end;
        }
        let mut sanpseq = extra_seq.as_bytes().to_vec();
        sanpseq.extend_from_slice(&self.get_ref_range(san_start, san_end));
        if let Some((inv5, _)) = inv_flanks.as_ref() {
            sanpseq = inv5.as_bytes().to_vec();
        }

        let r3 = self.find_mm3(pos, &sanpseq, data);
        let r5 = self.find_mm5(
            pos + dellen + extra_seq.len() as i64 - extrains_len - 1,
            &wupseq,
            data,
        );

        if data
            .non_insertion_variants
            .get(&pos)
            .and_then(|m| m.get(desc))
            .is_none()
        {
            return;
        }

        let dcnt = fixed_dcnt;
        if dcnt == 0 {
            return;
        }

        let mut all_mm = Vec::new();
        all_mm.extend(r3.mismatches.iter().cloned());
        all_mm.extend(r5.mismatches.iter().cloned());

        for mm in all_mm {
            let mm_bytes: Vec<u8> = mm
                .mismatch_sequence
                .as_bytes()
                .iter()
                .map(|b| b.to_ascii_uppercase())
                .collect();
            if mm_bytes.is_empty() {
                continue;
            }

            let key = if mm_bytes.len() == 1 {
                VarDesc::SNV {
                    ref_base: mm_bytes[0],
                }
            } else {
                let mut raw = Vec::with_capacity(mm_bytes.len() + 1);
                raw.push(mm_bytes[0]);
                raw.push(b'&');
                raw.extend_from_slice(&mm_bytes[1..]);
                VarDesc::Raw { desc: raw.into() }
            };
            let tv_snapshot = data
                .non_insertion_variants
                .get(&mm.mismatch_position)
                .and_then(|m| m.get(&key))
                .cloned();
            let Some(tv) = tv_snapshot else {
                continue;
            };

            if tv.alt_depth == 0 {
                continue;
            }

            let mean_qual = tv.mean_qual / tv.alt_depth as f64;
            if mean_qual
                < crate::scopedata::global_read_only_scope::instance()
                    .conf
                    .goodq
            {
                continue;
            }

            let nm_threshold = if mm.end == 3 { r3.nm } else { r5.nm };
            let mean_pos = tv.mean_pos / tv.alt_depth as f64;
            if mean_pos > nm_threshold as f64 + 4.0 {
                continue;
            }

            if tv.alt_depth >= dcnt + dellen as usize || (tv.alt_depth as f64 / dcnt as f64) >= 8.0
            {
                continue;
            }

            if mm.mismatch_position > pos && mm.end == 5 {
                let f = if mean_pos > 0.0 {
                    ((mm.mismatch_position - pos) as f64) / mean_pos
                } else {
                    1.0
                };
                let f = f.clamp(0.0, 1.0);
                let add = (tv.alt_depth as f64 * f) as usize;
                Self::emit_refcov_inc_diag(data, pos, add, "realigner_del_mm_fraction");
                *data.ref_coverage.entry(pos).or_insert(0) += add;

                if let Some(ref_base) = self.get_ref_base(pos) {
                    if let Some(pos_map) = data.non_insertion_variants.get_mut(&pos) {
                        let ref_key = VarDesc::SNV { ref_base };
                        if let Some(ref_var) = pos_map.get_mut(&ref_key) {
                            adj_ref_cnt(&tv, ref_var, dellen);
                        }
                    }
                }
            }

            let tv_owned = data
                .non_insertion_variants
                .get_mut(&mm.mismatch_position)
                .and_then(|m| m.remove(&key));
            let Some(tv_owned) = tv_owned else {
                continue;
            };

            if let Some(map) = data.non_insertion_variants.get_mut(&mm.mismatch_position) {
                if map.is_empty() {
                    data.non_insertion_variants.remove(&mm.mismatch_position);
                }
            }

            let ref_key = if mm.mismatch_position > pos && mm.end == 3 {
                self.get_ref_base(pos)
                    .map(|ref_base| VarDesc::SNV { ref_base })
            } else {
                None
            };

            if let Some(pos_map) = data.non_insertion_variants.get_mut(&pos) {
                let mut ref_var = ref_key.as_ref().and_then(|k| pos_map.remove(k));
                if let Some(vref) = pos_map.get_mut(desc) {
                    if let Some(ref_var_mut) = ref_var.as_mut() {
                        adj_cnt_with_ref(vref, &tv_owned, Some(ref_var_mut));
                    } else {
                        adj_cnt(vref, &tv_owned);
                    }
                }
                if let Some(ref_key) = ref_key {
                    if let Some(ref_var) = ref_var {
                        pos_map.insert(ref_key, ref_var);
                    }
                }
            }
        }

        if r3.misp != 0 && r3.mismatches.len() == 1 {
            if let Some(base) = r3.misnt {
                let key = VarDesc::SNV { ref_base: base };
                if let Some(map) = data.non_insertion_variants.get_mut(&r3.misp) {
                    if let Some(var) = map.get(&key) {
                        if var.alt_depth < dcnt {
                            map.remove(&key);
                        }
                    }
                }
            }
        }

        if r5.misp != 0 && r5.mismatches.len() == 1 {
            if let Some(base) = r5.misnt {
                let key = VarDesc::SNV { ref_base: base };
                if let Some(map) = data.non_insertion_variants.get_mut(&r5.misp) {
                    if let Some(var) = map.get(&key) {
                        if var.alt_depth < dcnt {
                            map.remove(&key);
                        }
                    }
                }
            }
        }

        for sc5pp in r5.scp {
            let refcov_before = if pos == 6970385 {
                data.ref_coverage.get(&6970385_i64).copied().unwrap_or(0)
            } else {
                0
            };
            if let Some(tv) = data.soft_clips_5end.get_mut(&sc5pp) {
                if tv.used() {
                    continue;
                }
                if dcnt <= 2 && tv.var.alt_depth / dcnt > 5 {
                    continue;
                }
                let seq = find_conseq_transient(tv, 0);
                if seq.is_empty() {
                    continue;
                }
                let is_match = Self::is_match_bytes(&seq, &wupseq, -1);
                if is_match {
                    cache_transient_conseq(tv, &seq);
                    if sc5pp > pos {
                        Self::emit_refcov_inc_diag_with_before(
                            refcov_before,
                            pos,
                            tv.var.alt_depth,
                            "realigner_del_sc5",
                        );
                        *data.ref_coverage.entry(pos).or_insert(0) += tv.var.alt_depth;
                    }
                    if let Some(vref) = data
                        .non_insertion_variants
                        .get_mut(&pos)
                        .and_then(|m| m.get_mut(desc))
                    {
                        adj_cnt(vref, &tv.var);
                    }
                    tv.mark_used();
                }
            }
        }

        for sc3pp in r3.scp {
            let refcov_before = if pos == 6970385 {
                data.ref_coverage.get(&6970385_i64).copied().unwrap_or(0)
            } else {
                0
            };
            if let Some(tv) = data.soft_clips_3end.get_mut(&sc3pp) {
                if tv.used() {
                    continue;
                }
                if dcnt <= 2 && tv.var.alt_depth / dcnt > 5 {
                    continue;
                }
                let seq = find_conseq_transient(tv, 0);
                if seq.is_empty() {
                    continue;
                }
                let mseq = Self::substr_bytes(&sanpseq, sc3pp - pos, None);
                let is_match = Self::is_match_bytes(&seq, &mseq, 1);
                if is_match {
                    cache_transient_conseq(tv, &seq);
                    if sc3pp <= pos {
                        Self::emit_refcov_inc_diag_with_before(
                            refcov_before,
                            pos,
                            tv.var.alt_depth,
                            "realigner_del_sc3",
                        );
                        *data.ref_coverage.entry(pos).or_insert(0) += tv.var.alt_depth;
                    }

                    let ref_key = if sc3pp > pos {
                        self.get_ref_base(pos)
                            .map(|ref_base| VarDesc::SNV { ref_base })
                    } else {
                        None
                    };

                    if let Some(pos_map) = data.non_insertion_variants.get_mut(&pos) {
                        let mut ref_var = ref_key.as_ref().and_then(|k| pos_map.remove(k));
                        if let Some(vref) = pos_map.get_mut(desc) {
                            if let Some(ref_var_mut) = ref_var.as_mut() {
                                adj_cnt_with_ref(vref, &tv.var, Some(ref_var_mut));
                            } else {
                                adj_cnt(vref, &tv.var);
                            }
                        }
                        if let Some(ref_key) = ref_key {
                            if let Some(ref_var) = ref_var {
                                pos_map.insert(ref_key, ref_var);
                            }
                        }
                    }
                    tv.mark_used();
                }
            }
        }

        // Java realigndel noPassingReads block:
        // if (pe - p >= 5 && pe - p < maxReadLength - 10
        //     && h != null && h.varsCount != 0
        //     && noPassingReads(chr, p, pe, bams)
        //     && vref.varsCount > 2 * h.varsCount * (1 - (pe - p) / maxReadLength)) {
        //     adjCnt(vref, h, h);
        // }
        let pe = pos + dellen + extra_seq.len() as i64 - extrains_len;
        let gap_len = pe - pos;
        if enable_passing_reads_check
            && !self.bam_paths.is_empty()
            && gap_len >= 5
            && gap_len < data.max_read_length.saturating_sub(10) as i64
        {
            if let Some(ref_base) = self.get_ref_base(pos) {
                let ref_key = VarDesc::SNV { ref_base };
                if let Some(pos_map) = data.non_insertion_variants.get_mut(&pos) {
                    let h_snapshot = pos_map.get(&ref_key).cloned();
                    let vref_snapshot = pos_map.get(desc).cloned();

                    if let (Some(h), Some(vref_now)) = (h_snapshot, vref_snapshot) {
                        if h.alt_depth != 0
                            && self.no_passing_reads(pos, pe)
                            && vref_now.alt_depth as f64
                                > 2.0
                                    * h.alt_depth as f64
                                    * (1.0 - gap_len as f64 / data.max_read_length as f64)
                        {
                            if let Some(mut h_work) = pos_map.remove(&ref_key) {
                                let h_src = h_work.clone();
                                if let Some(vref_mut) = pos_map.get_mut(desc) {
                                    adj_cnt_with_ref(vref_mut, &h_src, Some(&mut h_work));
                                }
                                pos_map.insert(ref_key, h_work);
                            }
                        }
                    }
                }
            }
        }
    }

    fn realign_with_softclips_5end(
        &self,
        pos: i64,
        desc: &VarDesc,
        data: &mut RealignedVariationData,
        wupseq: &[u8],
    ) {
        let positions: Vec<i64> = data.soft_clips_5end.keys().cloned().collect();
        for sc_pos in positions {
            let refcov_before = if pos == 6970385 {
                data.ref_coverage.get(&6970385_i64).copied().unwrap_or(0)
            } else {
                0
            };
            let Some(sclip) = data.soft_clips_5end.get_mut(&sc_pos) else {
                continue;
            };
            if sclip.used() {
                continue;
            }

            let seq = find_conseq_transient(sclip, 0);
            if seq.is_empty() {
                continue;
            }

            if !Self::is_match_bytes(&seq, wupseq, -1) {
                continue;
            }

            cache_transient_conseq(sclip, &seq);
            if let Some(var_map) = data.non_insertion_variants.get_mut(&pos) {
                if let Some(variant) = var_map.get_mut(desc) {
                    adj_cnt(variant, &sclip.var);
                }
            }
            if sc_pos > pos {
                Self::emit_refcov_inc_diag_with_before(
                    refcov_before,
                    pos,
                    sclip.var.alt_depth,
                    "realigner_softclips_5end",
                );
                *data.ref_coverage.entry(pos).or_insert(0) += sclip.var.alt_depth;
            }
            sclip.mark_used();
        }
    }

    fn realign_with_softclips_3end(
        &self,
        pos: i64,
        desc: &VarDesc,
        data: &mut RealignedVariationData,
        sanpseq: &[u8],
    ) {
        let positions: Vec<i64> = data.soft_clips_3end.keys().cloned().collect();
        for sc_pos in positions {
            let refcov_before = if pos == 6970385 {
                data.ref_coverage.get(&6970385_i64).copied().unwrap_or(0)
            } else {
                0
            };
            let Some(sclip) = data.soft_clips_3end.get_mut(&sc_pos) else {
                continue;
            };
            if sclip.used() {
                continue;
            }

            let seq = find_conseq_transient(sclip, 0);
            if seq.is_empty() {
                continue;
            }

            if !Self::is_match_bytes(&seq, sanpseq, 1) {
                continue;
            }

            cache_transient_conseq(sclip, &seq);
            if sc_pos <= pos {
                Self::emit_refcov_inc_diag_with_before(
                    refcov_before,
                    pos,
                    sclip.var.alt_depth,
                    "realigner_softclips_3end",
                );
                *data.ref_coverage.entry(pos).or_insert(0) += sclip.var.alt_depth;
            }
            if let Some(var_map) = data.non_insertion_variants.get_mut(&pos) {
                if let Some(variant) = var_map.get_mut(desc) {
                    adj_cnt(variant, &sclip.var);
                }
            }
            sclip.mark_used();
        }
    }

    fn is_match_bytes(seq1: &[u8], seq2: &[u8], dir: i32) -> bool {
        let s1 = String::from_utf8_lossy(seq1);
        let s2 = String::from_utf8_lossy(seq2);
        Self::is_match(&s1, &s2, dir)
    }

    fn is_has_and_equals(&self, pos: i64, ch: u8) -> bool {
        self.get_ref_base(pos).map(|b| b == ch).unwrap_or(false)
    }

    fn is_has_and_not_equals(&self, pos: i64, ch: u8) -> bool {
        self.get_ref_base(pos).map(|b| b != ch).unwrap_or(false)
    }

    fn find_mm5(
        &self,
        position: i64,
        wupseq: &[u8],
        data: &mut RealignedVariationData,
    ) -> MismatchResult {
        let seq: Vec<u8> = wupseq
            .iter()
            .copied()
            .filter(|b| *b != b'#' && *b != b'^')
            .collect();

        let longmm = 3usize;
        let mut mismatches = Vec::new();
        let mut n = 0usize;
        let mut mn = 0usize;
        let mut mcnt = 0usize;
        let mut buf: Vec<u8> = Vec::new();
        let mut sc5p: Vec<i64> = Vec::new();

        while let Some(ch) = char_at_neg(&seq, n as isize) {
            if self.is_has_and_not_equals(position - n as i64, ch) && mcnt < longmm {
                buf.insert(0, ch);
                mismatches.push(Mismatch {
                    mismatch_sequence: String::from_utf8_lossy(&buf).to_string(),
                    mismatch_position: position - n as i64,
                    end: 5,
                });
                n += 1;
                mcnt += 1;
            } else {
                break;
            }
        }

        sc5p.push(position + 1);
        let mut misp = 0i64;
        let mut misnt = None;

        if buf.len() == 1 {
            while let Some(ch) = char_at_neg(&seq, n as isize) {
                if self.is_has_and_equals(position - n as i64, ch) {
                    n += 1;
                    if n != 0 {
                        mn += 1;
                    }
                } else {
                    break;
                }
            }

            if mn > 1 {
                let mut n2 = 0usize;
                while -1isize - n as isize - 1 - n2 as isize >= 0 {
                    let Some(ch) = char_at_neg(&seq, (n + 1 + n2) as isize) else {
                        break;
                    };
                    if !self.is_has_and_equals(position - n as i64 - 1 - n2 as i64, ch) {
                        break;
                    }
                    n2 += 1;
                }

                if n2 > 2 {
                    sc5p.push(position - n as i64 - n2 as i64);
                    misp = position - n as i64;
                    misnt = char_at_neg(&seq, n as isize);
                    let mark_pos = position - n as i64 - n2 as i64;
                    if let Some(sc) = data.soft_clips_5end.get_mut(&mark_pos) {
                        sc.mark_used();
                    }
                    mn += n2;
                } else {
                    sc5p.push(position - n as i64);
                    let mark_pos = position - n as i64;
                    if let Some(sc) = data.soft_clips_5end.get_mut(&mark_pos) {
                        sc.mark_used();
                    }
                }
            }
        }

        MismatchResult {
            mismatches,
            scp: sc5p,
            nm: mn,
            misp,
            misnt,
        }
    }

    fn find_mm3(
        &self,
        p: i64,
        sanpseq: &[u8],
        data: &mut RealignedVariationData,
    ) -> MismatchResult {
        let seq: Vec<u8> = sanpseq
            .iter()
            .copied()
            .filter(|b| *b != b'#' && *b != b'^')
            .collect();
        let longmm = 3usize;
        let mut mismatches = Vec::new();
        let mut n = 0usize;
        let mut mn = 0usize;
        let mut mcnt = 0usize;
        let mut sc3p: Vec<i64> = Vec::new();
        let mut buf: Vec<u8> = Vec::new();

        while n < seq.len() && self.get_ref_base(p + n as i64) == Some(seq[n]) {
            n += 1;
        }

        sc3p.push(p + n as i64);
        let tbp = p + n as i64;

        while mcnt <= longmm && n < seq.len() {
            let ref_base = self.get_ref_base(p + n as i64);
            if matches!(ref_base, Some(base) if base == seq[n]) {
                break;
            }
            buf.push(seq[n]);
            mismatches.push(Mismatch {
                mismatch_sequence: String::from_utf8_lossy(&buf).to_string(),
                mismatch_position: tbp,
                end: 3,
            });
            n += 1;
            mcnt += 1;
        }

        let mut misp = 0i64;
        let mut misnt = None;

        if buf.len() == 1 {
            while n < seq.len() && self.is_has_and_equals(p + n as i64, seq[n]) {
                n += 1;
                if n != 0 {
                    mn += 1;
                }
            }

            if mn > 1 {
                let mut n2 = 0usize;
                while n + n2 + 1 < seq.len()
                    && self.is_has_and_equals(p + n as i64 + 1 + n2 as i64, seq[n + n2 + 1])
                {
                    n2 += 1;
                }
                if n2 > 2 && n + n2 + 1 < seq.len() {
                    sc3p.push(p + n as i64 + n2 as i64);
                    misp = p + n as i64;
                    misnt = Some(seq[n]);
                    if let Some(sc) = data.soft_clips_3end.get_mut(&(p + n as i64 + n2 as i64)) {
                        sc.mark_used();
                    }
                    mn += n2;
                } else {
                    sc3p.push(p + n as i64);
                    if let Some(sc) = data.soft_clips_3end.get_mut(&(p + n as i64)) {
                        sc.mark_used();
                    }
                }
            }
        }

        MismatchResult {
            mismatches,
            scp: sc3p,
            nm: mn,
            misp,
            misnt,
        }
    }

    fn join_ref(&self, from: i64, to: i64) -> Vec<u8> {
        self.get_ref_range(from, to)
    }

    fn join_ref_float(&self, from: i64, to: f64) -> Vec<u8> {
        let mut out = Vec::new();
        let mut i = from;
        while (i as f64) < to {
            if let Some(base) = self.get_ref_base(i) {
                out.push(base);
            }
            i += 1;
        }
        out
    }

    fn join_ref_for_5_lgins(&self, from: i64, to: i64, seq: &[u8], extra: &[u8]) -> Vec<u8> {
        let mut out = Vec::new();
        let seq_len = seq.len();
        let extra_len = extra.len();
        let usable_len = seq_len.saturating_sub(extra_len) as i64;

        for i in from..=to {
            if to - i < usable_len {
                let idx = (to - i) as usize + extra_len;
                if let Some(base) = seq.get(idx) {
                    out.push(*base);
                }
            } else if let Some(base) = self.get_ref_base(i) {
                out.push(base);
            }
        }

        out
    }

    fn join_ref_for_3_lgins(
        &self,
        from: i64,
        to: i64,
        shift5: i64,
        seq: &[u8],
        extra: &[u8],
    ) -> Vec<u8> {
        let mut out = Vec::new();
        let seq_len = seq.len();
        let extra_len = extra.len();
        let usable_len = seq_len.saturating_sub(extra_len) as i64;

        for i in from..=to {
            let rel = i - from;
            if rel >= shift5 && rel - shift5 < usable_len {
                let idx = (rel - shift5) as usize + extra_len;
                if let Some(base) = seq.get(idx) {
                    out.push(*base);
                }
            } else if let Some(base) = self.get_ref_base(i) {
                out.push(base);
            }
        }

        out
    }

    fn adj_ins_pos(&self, mut bi: i64, ins: &[u8]) -> BaseInsertion {
        let mut n = 1usize;
        let len = ins.len();
        let mut ins_seq = ins.to_vec();

        while let Some(ref_base) = self.get_ref_base(bi) {
            if len == 0 {
                break;
            }
            let idx = len.saturating_sub(n);
            if idx >= len || ref_base != ins_seq[idx] {
                break;
            }
            n += 1;
            if n > len {
                n = 1;
            }
            bi -= 1;
        }

        if n > 1 {
            let tail = Self::substr_bytes(&ins_seq, 1 - n as i64, None);
            let head = Self::substr_bytes(&ins_seq, 0, Some(1 - n as i64));
            let mut rotated = tail;
            rotated.extend_from_slice(&head);
            ins_seq = rotated;
        }

        BaseInsertion {
            base_insert: bi,
            insertion_sequence: ins_seq,
            base_insert2: bi,
        }
    }

    fn find_bi(&self, seq: &[u8], position: i64, dir: i64) -> BaseInsertion {
        let maxmm = 3usize;
        let dir_ext = if dir == -1 { 1 } else { 0 };
        let mut score = 0i64;
        let mut bi = 0i64;
        let mut ins: Vec<u8> = Vec::new();
        let mut bi2 = 0i64;

        let ref_end = self.ref_start + self.reference_seq.len() as i64 - 1;
        let chr_len = self
            .chromosome
            .as_ref()
            .and_then(|chrom| {
                crate::scopedata::global_read_only_scope::instance()
                    .chr_lens
                    .get(chrom)
                    .copied()
            })
            .map(|len| len as i64)
            .unwrap_or(ref_end);

        for n in 6..seq.len() {
            if position + 6 >= chr_len {
                break;
            }
            let mut mm = 0usize;
            let mut i = 0usize;
            let mut m: std::collections::HashSet<u8> = std::collections::HashSet::new();

            while i + n < seq.len() {
                let ref_pos = position + dir * i as i64 - dir_ext;
                if ref_pos < 1 || ref_pos > chr_len {
                    break;
                }
                let ref_base = self.get_ref_base(ref_pos);
                if ref_base.map_or(true, |base| seq[i + n] != base) {
                    mm += 1;
                } else {
                    m.insert(seq[i + n]);
                }
                if mm > maxmm {
                    break;
                }
                i += 1;
            }

            let mnt = m.len();
            if mnt < 2 {
                continue;
            }
            let mm_rate_ok = i > 0 && (mm as f64 / i as f64) < 0.15;
            let end_match = i + n >= seq.len().saturating_sub(1);
            let long_match = i + n == seq.len();
            if (mnt >= 3 && end_match && i >= 8 && mm_rate_ok)
                || (mnt >= 2 && mm == 0 && long_match && n >= 20 && i >= 8)
            {
                let mut insert = Self::substr_bytes(seq, 0, Some(n as i64));
                let mut extra: Vec<u8> = Vec::new();
                let mut ept = 0usize;
                while n + ept + 1 < seq.len() {
                    let pos1 = position + ept as i64 * dir - dir_ext;
                    let pos2 = position + (ept + 1) as i64 * dir - dir_ext;
                    let first_matches = self
                        .get_ref_base(pos1)
                        .map(|base| seq[n + ept] == base)
                        .unwrap_or(false);
                    let second_matches = self
                        .get_ref_base(pos2)
                        .map(|base| seq[n + ept + 1] == base)
                        .unwrap_or(false);
                    if first_matches && second_matches {
                        break;
                    }
                    extra.push(seq[n + ept]);
                    ept += 1;
                }

                if dir == -1 {
                    insert.extend_from_slice(&extra);
                    insert.reverse();
                    if !extra.is_empty() {
                        let pos = insert.len().saturating_sub(extra.len());
                        insert.insert(pos, b'&');
                    }
                    if mm == 0 && long_match {
                        bi = position - 1 - extra.len() as i64;
                        ins = insert;
                        bi2 = position - 1;
                        if extra.is_empty() {
                            let tpl = self.adj_ins_pos(bi, &ins);
                            bi = tpl.base_insert;
                            ins = tpl.insertion_sequence;
                            bi2 = tpl.base_insert2;
                        }
                        return BaseInsertion {
                            base_insert: bi,
                            insertion_sequence: ins,
                            base_insert2: bi2,
                        };
                    } else if (i as i64 - mm as i64) > score {
                        bi = position - 1 - extra.len() as i64;
                        ins = insert;
                        bi2 = position - 1;
                        score = i as i64 - mm as i64;
                    }
                } else {
                    let mut s = -1i64;
                    if !extra.is_empty() {
                        insert.push(b'&');
                        insert.extend_from_slice(&extra);
                    } else {
                        while s >= -(n as i64) {
                            let seq_ch = Self::char_at(&insert, s);
                            let ref_ch = self.get_ref_base(position + s);
                            if seq_ch.is_some()
                                && ref_ch.is_some()
                                && seq_ch.unwrap() == ref_ch.unwrap()
                            {
                                s -= 1;
                            } else {
                                break;
                            }
                        }
                        if s < -1 {
                            let tins = Self::substr_bytes(&insert, s + 1, Some(1 - s));
                            let truncate_at = (insert.len() as i64 + s + 1).max(0) as usize;
                            insert.truncate(truncate_at);
                            let mut rotated = tins;
                            rotated.extend_from_slice(&insert);
                            insert = rotated;
                        }
                    }

                    if mm == 0 && long_match {
                        bi = position + s;
                        ins = insert;
                        bi2 = position + s + extra.len() as i64;
                        if extra.is_empty() {
                            let tpl = self.adj_ins_pos(bi, &ins);
                            bi = tpl.base_insert;
                            ins = tpl.insertion_sequence;
                            bi2 = tpl.base_insert2;
                        }
                        return BaseInsertion {
                            base_insert: bi,
                            insertion_sequence: ins,
                            base_insert2: bi2,
                        };
                    } else if (i as i64 - mm as i64) > score {
                        bi = position + s;
                        ins = insert;
                        bi2 = position + s + extra.len() as i64;
                        score = i as i64 - mm as i64;
                    }
                }
            }
        }

        if bi2 == bi && !ins.is_empty() && bi != 0 {
            let tpl = self.adj_ins_pos(bi, &ins);
            bi = tpl.base_insert;
            ins = tpl.insertion_sequence;
        }

        BaseInsertion {
            base_insert: bi,
            insertion_sequence: ins,
            base_insert2: bi2,
        }
    }

    fn find_bp(&self, sequence: &[u8], start_position: i64, direction: i64) -> i64 {
        let max_mm = 3i64;
        let mut bp = 0i64;
        let mut score = 0i64;
        let chr_len = self
            .chromosome
            .as_ref()
            .and_then(|chrom| {
                crate::scopedata::global_read_only_scope::instance()
                    .chr_lens
                    .get(chrom)
                    .copied()
            })
            .map(|len| len as i64)
            .unwrap_or(0);

        if chr_len <= 0 {
            return 0;
        }

        for n in 0..50i64 {
            let mut mm = 0i64;
            let mut i = 0usize;
            let mut matched_bases: std::collections::HashSet<u8> = std::collections::HashSet::new();

            while i < sequence.len() {
                let ref_pos = start_position + direction * n + direction * i as i64;
                if ref_pos < 1 || ref_pos > chr_len {
                    break;
                }

                let seq_base = sequence[i];
                if self
                    .get_ref_base(ref_pos)
                    .map(|base| base == seq_base)
                    .unwrap_or(false)
                {
                    matched_bases.insert(seq_base);
                } else {
                    mm += 1;
                }

                if mm > max_mm - n / 100 {
                    break;
                }

                i += 1;
            }

            if matched_bases.len() < 3 {
                continue;
            }

            if i == 0 {
                continue;
            }

            if mm <= max_mm - n / 100
                && i >= sequence.len().saturating_sub(2)
                && i >= (8 + n / 10) as usize
                && (mm as f64 / i as f64) < 0.12
            {
                let lbp =
                    start_position + direction * n - if direction < 0 { direction } else { 0 };
                if mm == 0 && i == sequence.len() {
                    return lbp;
                } else if (i as i64 - mm) > score {
                    bp = lbp;
                    score = i as i64 - mm;
                }
            }
        }

        bp
    }

    fn find_match(
        &self,
        seq: &[u8],
        _position: i64,
        dir: i64,
        seed_len: usize,
        mm: usize,
    ) -> Match {
        let mut seq_work = seq.to_vec();
        if dir == -1 {
            seq_work.reverse();
        }

        let mut persistent_extra: Vec<u8> = Vec::new();

        if seq_work.len() < seed_len {
            return Match {
                base_position: 0,
                matched_sequence: Vec::new(),
            };
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
                return Match {
                    base_position: bp,
                    matched_sequence: persistent_extra,
                };
            } else {
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
                        let Some(ch0) = sseq.get(0).copied() else {
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
                        return Match {
                            base_position: bp,
                            matched_sequence: extra,
                        };
                    }
                }
            }
        }

        Match {
            base_position: 0,
            matched_sequence: Vec::new(),
        }
    }

    fn del_desc_to_key(&self, desc: &[u8]) -> Option<VarDesc> {
        if desc.first() != Some(&b'-') {
            return None;
        }
        let mut idx = 1usize;
        while idx < desc.len() && desc[idx].is_ascii_digit() {
            idx += 1;
        }
        if idx == 1 {
            return None;
        }
        let len = std::str::from_utf8(&desc[1..idx])
            .ok()?
            .parse::<u32>()
            .ok()?;
        let mut ins_or_del_len = InsOrDelLen::None;
        if idx < desc.len() && desc[idx] == b'^' {
            let seq_bytes = &desc[(idx + 1)..];
            if !seq_bytes.is_empty() {
                let seq: SmallVecBytes = seq_bytes.iter().map(|b| b.to_ascii_uppercase()).collect();
                ins_or_del_len = InsOrDelLen::InsSeq(seq);
            }
        }

        Some(VarDesc::Del {
            len,
            match_seq: SmallVecBytes::new(),
            ins_or_del_len,
            mismatch_seq: SmallVecBytes::new(),
        })
    }

    fn realign_deletion_for_desc(
        &self,
        pos: i64,
        desc: &VarDesc,
        data: &mut RealignedVariationData,
    ) {
        let dcnt = data
            .non_insertion_variants
            .get(&pos)
            .and_then(|m| m.get(desc))
            .map(|v| v.alt_depth)
            .unwrap_or(0);
        self.realign_deletion_mismatches(pos, desc, dcnt, data, false);
    }

    fn get_ref_range(&self, start: i64, end: i64) -> Vec<u8> {
        if start > end {
            return Vec::new();
        }
        let mut seq = Vec::new();
        for pos in start..=end {
            if let Some(base) = self.get_ref_base(pos) {
                seq.push(base);
            }
        }
        seq
    }

    fn no_passing_reads(&self, start: i64, end: i64) -> bool {
        let Some(chr) = self.chromosome.as_deref() else {
            return false;
        };
        if start <= 0 || end <= start || self.bam_paths.is_empty() {
            return false;
        }

        let mut cnt = 0usize;
        let mut midcnt = 0usize;
        let dlen = (end - start) as u32;

        for bam in &self.bam_paths {
            let mut reader = match BamReader::open(bam) {
                Ok(reader) => reader,
                Err(_) => return false,
            };

            if reader.fetch(chr, start as usize, end as usize).is_err() {
                return false;
            }

            let mut record = Record::new();
            loop {
                match reader.read(&mut record) {
                    Ok(true) => {}
                    Ok(false) => break,
                    Err(_) => return false,
                }

                if Self::record_contains_java_deletion_token(&record, dlen) {
                    continue;
                }

                let read_start = record.pos() + 1;
                let read_end = read_start
                    + Self::aligned_length_excluding_softclips_and_insertions(&record) as i64;

                if read_end > end + 2 && read_start < start - 2 {
                    cnt += 1;
                }

                if read_start < start - 2 && read_end > start && read_end < end {
                    midcnt += 1;
                }
            }
        }

        cnt == 0 && midcnt + 1 > 0
    }

    fn record_contains_java_deletion_token(record: &Record, dlen: u32) -> bool {
        let token = format!("{}D", dlen);
        record.cigar().to_string().contains(&token)
    }

    fn aligned_length_excluding_softclips_and_insertions(record: &Record) -> u32 {
        let mut len = 0u32;
        for op in record.cigar().iter() {
            match *op {
                Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) | Cigar::Del(l) => {
                    len += l;
                }
                _ => {}
            }
        }
        len
    }

    fn get_ref_base(&self, pos: i64) -> Option<u8> {
        Self::get_ref_base_from_window(&self.reference_seq, self.ref_start, pos).or_else(|| {
            if self.original_ref_start == self.ref_start
                && Arc::ptr_eq(&self.original_reference_seq, &self.reference_seq)
            {
                None
            } else {
                Self::get_ref_base_from_window(
                    &self.original_reference_seq,
                    self.original_ref_start,
                    pos,
                )
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
                char_at_neg(seq, n as isize).unwrap_or(b'N')
            };
            if seq_base != ref_base {
                mismatches += 1;
            }
        }
        mismatches <= mm && (mismatches as f64 / seq.len() as f64) < 0.15
    }
    /// Check if a sequence has low complexity (>75% of one base or <3 different bases)
    pub fn is_low_complex_seq(seq: &str) -> bool {
        let len = seq.len();
        if len == 0 {
            return true;
        }

        let a = Self::count(seq, 'A');
        if a as f64 / len as f64 > 0.75 {
            return true;
        }

        let t = Self::count(seq, 'T');
        if t as f64 / len as f64 > 0.75 {
            return true;
        }

        let g = Self::count(seq, 'G');
        if g as f64 / len as f64 > 0.75 {
            return true;
        }

        let c = Self::count(seq, 'C');
        if c as f64 / len as f64 > 0.75 {
            return true;
        }

        // Count how many different bases are present
        let mut nt_cnt = 0;
        if a > 0 {
            nt_cnt += 1;
        }
        if t > 0 {
            nt_cnt += 1;
        }
        if g > 0 {
            nt_cnt += 1;
        }
        if c > 0 {
            nt_cnt += 1;
        }

        nt_cnt < 3
    }

    /// Count occurrences of a character in a string
    fn count(seq: &str, ch: char) -> usize {
        seq.chars().filter(|&c| c == ch).count()
    }

    /// Find if two sequences match with no more than MM mismatches
    /// and total mismatches no more than 15% of length
    pub fn is_match(seq1: &str, seq2: &str, dir: i32) -> bool {
        Self::is_match_with_threshold(seq1, seq2, dir, 3)
    }

    /// Find if two sequences match with specified MM threshold
    pub fn is_match_with_threshold(seq1: &str, seq2: &str, dir: i32, mm_threshold: usize) -> bool {
        // Remove special characters from seq2
        let seq2_clean: String = seq2.chars().filter(|c| *c != '#' && *c != '^').collect();

        let mut mm = 0;
        let min_len = seq1.len().min(seq2_clean.len());

        for n in 0..min_len {
            let seq1_ch = seq1.chars().nth(n).unwrap_or('N');
            let seq2_idx = if dir == 1 {
                n
            } else {
                // For reverse direction, index from the end
                seq2_clean.len() - 1 - n
            };
            let seq2_ch = seq2_clean.chars().nth(seq2_idx).unwrap_or('N');

            if seq1_ch != seq2_ch {
                mm += 1;
            }
        }

        (mm <= mm_threshold) && (mm as f64 / seq1.len() as f64) < 0.15
    }

    /// Find 3'/5' end matches between two sequences
    /// Returns Match35 with (matched_5_end, matched_3_end, max_matched_length)
    pub fn find_35_match(seq5: &str, seq3: &str) -> Match35 {
        const LONG_MISMATCH: usize = 2;

        let mut max_matched_length = 0;
        let mut b3 = 0;
        let mut b5 = 0;

        // Convert to char vectors for efficient indexing
        let seq5_chars: Vec<char> = seq5.chars().collect();
        let seq3_chars: Vec<char> = seq3.chars().collect();

        // Ensure sequences are long enough for comparison
        if seq5_chars.len() <= 8 || seq3_chars.len() <= 8 {
            return Match35 {
                matched_5_end: b5,
                matched_3_end: b3,
                max_matched_length,
            };
        }

        // Java: for (int i = 0; i < seq5.length() - 8; i++)
        for i in 0..(seq5_chars.len() - 8) {
            // Java: for (int j = 1; j < seq3.length() - 8; j++)
            for j in 1..(seq3_chars.len() - 8) {
                let mut num_mismatch = 0;
                let mut total_length = 0;

                // Java: while (totalLength + j <= seq3.length() && i + totalLength <= seq5.length())
                while total_length + j <= seq3_chars.len() && i + total_length <= seq5_chars.len() {
                    // Java: substr(seq3, -j - totalLength, 1)
                    // which is seq3[len - j - totalLength]
                    let seq3_pos = seq3_chars.len() - j - total_length;
                    let seq5_pos = i + total_length;

                    // Java's substr returns empty string for out-of-bounds, which doesn't equal any char
                    // So out-of-bounds is treated as a mismatch
                    let is_mismatch =
                        if seq3_pos >= seq3_chars.len() || seq5_pos >= seq5_chars.len() {
                            true // Out of bounds = mismatch (like Java's empty string comparison)
                        } else {
                            seq3_chars[seq3_pos] != seq5_chars[seq5_pos]
                        };

                    // Java: if (!substr(seq3, -j - totalLength, 1).equals(substr(seq5, i + totalLength, 1)))
                    if is_mismatch {
                        num_mismatch += 1;
                    }

                    // Java: if (numberOfMismatch > longMismatch)
                    if num_mismatch > LONG_MISMATCH {
                        break;
                    }

                    total_length += 1;
                }

                // Java: if (totalLength - numberOfMismatch > maxMatchedLength
                //           && totalLength - numberOfMismatch > 8
                //           && numberOfMismatch / (double) totalLength < 0.1d
                //           && (totalLength + j >= seq3.length() || i + totalLength >= seq5.length()))
                let matched_quality = total_length.saturating_sub(num_mismatch);
                let is_end_match = (total_length + j >= seq3_chars.len())
                    || (i + total_length >= seq5_chars.len());
                let error_rate = if total_length > 0 {
                    num_mismatch as f64 / total_length as f64
                } else {
                    0.0
                };

                if matched_quality > max_matched_length
                    && matched_quality > 8
                    && error_rate < 0.1
                    && is_end_match
                {
                    max_matched_length = matched_quality;
                    b3 = j;
                    b5 = i;
                    // Java returns immediately when first match is found
                    return Match35 {
                        matched_5_end: b5,
                        matched_3_end: b3,
                        max_matched_length,
                    };
                }
            }
        }

        Match35 {
            matched_5_end: b5,
            matched_3_end: b3,
            max_matched_length,
        }
    }

    fn adj_ref_factor(ref_var: &mut Variant, factor_f: f64) {
        let mut factor = factor_f;
        if factor > 1.0 {
            factor = 1.0;
        }
        if factor < -1.0 {
            return;
        }

        let old_cnt = ref_var.alt_depth as i64;
        let new_cnt = old_cnt - (factor * old_cnt as f64) as i64;
        ref_var.alt_depth = new_cnt.max(0) as usize;
        ref_var.high_qual_read_cnt = (ref_var.high_qual_read_cnt as i64
            - (factor * ref_var.high_qual_read_cnt as f64) as i64)
            .max(0) as usize;
        ref_var.low_qual_read_cnt = (ref_var.low_qual_read_cnt as i64
            - (factor * ref_var.low_qual_read_cnt as f64) as i64)
            .max(0) as usize;

        let mut factor_cnt = if old_cnt != 0 {
            (ref_var.alt_depth as f64 - old_cnt as f64).abs() / old_cnt as f64
        } else {
            1.0
        };
        if factor < 0.0 {
            factor_cnt = -factor_cnt;
        }

        ref_var.mean_pos -= ref_var.mean_pos * factor_cnt;
        ref_var.mean_qual -= ref_var.mean_qual * factor_cnt;
        ref_var.mean_mapq -= ref_var.mean_mapq * factor_cnt;
        ref_var.nm -= factor * ref_var.nm;
        ref_var.alt_depth_fwd = (ref_var.alt_depth_fwd as i64
            - (factor * ref_var.alt_depth_fwd as f64) as i64)
            .max(0) as usize;
        ref_var.alt_depth_rev = (ref_var.alt_depth_rev as i64
            - (factor * ref_var.alt_depth_rev as f64) as i64)
            .max(0) as usize;

        correct_cnt(ref_var);
    }

    fn add_var_factor(var: &mut Variant, factor_f: f64) {
        if factor_f < -1.0 {
            return;
        }

        let old_alt_depth = var.alt_depth as i64;
        var.alt_depth = (old_alt_depth + (factor_f * old_alt_depth as f64) as i64).max(0) as usize;
        var.high_qual_read_cnt = (var.high_qual_read_cnt as i64
            + (factor_f * var.high_qual_read_cnt as f64) as i64)
            .max(0) as usize;
        var.low_qual_read_cnt = (var.low_qual_read_cnt as i64
            + (factor_f * var.low_qual_read_cnt as f64) as i64)
            .max(0) as usize;
        var.mean_pos += factor_f * var.mean_pos;
        var.mean_qual += factor_f * var.mean_qual;
        var.mean_mapq += factor_f * var.mean_mapq;
        var.nm += factor_f * var.nm;
        var.alt_depth_fwd = (var.alt_depth_fwd as i64
            + (factor_f * var.alt_depth_fwd as f64) as i64)
            .max(0) as usize;
        var.alt_depth_rev = (var.alt_depth_rev as i64
            + (factor_f * var.alt_depth_rev as f64) as i64)
            .max(0) as usize;
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
}

fn adj_cnt_with_ref(dest: &mut Variant, src: &Variant, reference: Option<&mut Variant>) {
    adj_cnt(dest, src);

    let Some(reference) = reference else {
        return;
    };

    reference.alt_depth = reference.alt_depth.saturating_sub(src.alt_depth);
    reference.high_qual_read_cnt = reference
        .high_qual_read_cnt
        .saturating_sub(src.high_qual_read_cnt);
    reference.low_qual_read_cnt = reference
        .low_qual_read_cnt
        .saturating_sub(src.low_qual_read_cnt);
    reference.mean_pos -= src.mean_pos;
    reference.mean_qual -= src.mean_qual;
    reference.mean_mapq -= src.mean_mapq;
    reference.nm -= src.nm;
    reference.alt_depth_fwd = reference.alt_depth_fwd.saturating_sub(src.alt_depth_fwd);
    reference.alt_depth_rev = reference.alt_depth_rev.saturating_sub(src.alt_depth_rev);
    correct_cnt(reference);
}

fn adj_ref_cnt(tv: &Variant, reference: &mut Variant, len: i64) {
    if tv.alt_depth == 0 {
        return;
    }

    let mean_pos = tv.mean_pos / tv.alt_depth as f64;
    let mut f = if mean_pos != 0.0 {
        (mean_pos - len as f64 + 1.0) / mean_pos
    } else {
        0.0
    };

    if f < 0.0 {
        return;
    }
    if f > 1.0 {
        f = 1.0;
    }

    let delta_alt = (f * tv.alt_depth as f64) as usize;
    let delta_high = (f * tv.high_qual_read_cnt as f64) as usize;
    let delta_low = (f * tv.low_qual_read_cnt as f64) as usize;
    let delta_fwd = (f * tv.alt_depth_fwd as f64) as usize;
    let delta_rev = (f * tv.alt_depth_rev as f64) as usize;

    reference.alt_depth = reference.alt_depth.saturating_sub(delta_alt);
    reference.high_qual_read_cnt = reference.high_qual_read_cnt.saturating_sub(delta_high);
    reference.low_qual_read_cnt = reference.low_qual_read_cnt.saturating_sub(delta_low);
    reference.mean_pos -= f * tv.mean_pos;
    reference.mean_qual -= f * tv.mean_qual;
    reference.mean_mapq -= f * tv.mean_mapq;
    reference.nm -= f * tv.nm;
    reference.alt_depth_fwd = reference.alt_depth_fwd.saturating_sub(delta_fwd);
    reference.alt_depth_rev = reference.alt_depth_rev.saturating_sub(delta_rev);
    correct_cnt(reference);
}

fn adj_cnt(dest: &mut Variant, src: &Variant) {
    dest.alt_depth += src.alt_depth;
    dest.extra_cnt += src.alt_depth;
    dest.high_qual_read_cnt += src.high_qual_read_cnt;
    dest.low_qual_read_cnt += src.low_qual_read_cnt;
    dest.mean_pos += src.mean_pos;
    dest.mean_qual += src.mean_qual;
    dest.mean_mapq += src.mean_mapq;
    dest.nm += src.nm;
    dest.alt_depth_fwd += src.alt_depth_fwd;
    dest.alt_depth_rev += src.alt_depth_rev;
    dest.pstd = true;
    dest.qstd = true;
}

fn correct_cnt(var: &mut Variant) {
    if var.mean_pos < 0.0 {
        var.mean_pos = 0.0;
    }
    if var.mean_qual < 0.0 {
        var.mean_qual = 0.0;
    }
    if var.mean_mapq < 0.0 {
        var.mean_mapq = 0.0;
    }
}

fn char_at_neg(seq: &[u8], offset_from_end: isize) -> Option<u8> {
    if offset_from_end < 0 {
        return None;
    }
    let len = seq.len() as isize;
    let pos = len - 1 - offset_from_end;
    if pos < 0 || pos >= len {
        None
    } else {
        Some(seq[pos as usize])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::reference::FastaReader;
    use crate::data::region::Region;
    use crate::mods::cigar_parser::CigarParser;
    use crate::scopedata::global_read_only_scope::{GlobalReadOnlyScope, INSTANCE, instance_arc};
    use rust_htslib::bam::{Read, Reader};
    use std::sync::Arc;

    #[test]
    fn test_is_low_complex_seq() {
        // Test cases from Java VariationRealignerTest
        assert!(VariantRealigner::is_low_complex_seq("AAAAAAAAA"));
        assert!(VariantRealigner::is_low_complex_seq("ATATATATATAT"));
        assert!(VariantRealigner::is_low_complex_seq("CCCCCCCCGA"));
        assert!(!VariantRealigner::is_low_complex_seq("ACGTACGTACGT"));
        assert!(!VariantRealigner::is_low_complex_seq("CCGTAACGGGGT"));
    }

    #[test]
    fn test_is_match() {
        // Test cases from Java VariationRealignerTest
        assert!(VariantRealigner::is_match("AAAAAAAAA", "AAAAAAAAA", 1));
        assert!(VariantRealigner::is_match("AAAAAAAAA", "AAAAAAAAA", -1));
        assert!(VariantRealigner::is_match(
            "ACGTACGTACGTACGT",
            "AAGTACTTACGTACGT",
            1
        ));
        assert!(VariantRealigner::is_match(
            "ACGTACGTACGTACGT",
            "AAGTACTTACGT",
            1
        ));

        assert!(!VariantRealigner::is_match(
            "ACGTACGTACGT",
            "AAGTACTTACGT",
            1
        ));
        assert!(!VariantRealigner::is_match("ACGTCAGCAT", "ACGACTGACT", 1));
    }

    #[test]
    fn test_find_35_match() {
        // Test cases from Java VariationRealignerTest
        let result1 =
            VariantRealigner::find_35_match("ACGTACGTACGTACGTACGTACGT", "TGCATGCATGCATGCATGCATGCA");
        assert_eq!(result1.matched_5_end, 0);
        assert_eq!(result1.matched_3_end, 1);
        assert_eq!(result1.max_matched_length, 24);

        let result2 =
            VariantRealigner::find_35_match("ACGTACGTACGTACGTACGTACGT", "TGCATGCATGCCTGCATGCATGCA");
        assert_eq!(result2.matched_5_end, 0);
        assert_eq!(result2.matched_3_end, 1);
        assert_eq!(result2.max_matched_length, 23);

        let result3 =
            VariantRealigner::find_35_match("ACGTACGTACGTACGTACGTACGT", "TGCATGCATGCATGGGTGCATGCA");
        assert_eq!(result3.matched_5_end, 0);
        assert_eq!(result3.matched_3_end, 1);
        assert_eq!(result3.max_matched_length, 22);

        let result4 =
            VariantRealigner::find_35_match("ACGTACGTACGTACGTACGTACGT", "TGCATGCATGCATGGGTGCAAAAA");
        assert_eq!(result4.matched_5_end, 0);
        assert_eq!(result4.matched_3_end, 13);
        assert_eq!(result4.max_matched_length, 12);

        let result5 =
            VariantRealigner::find_35_match("ACGTACGTACGTACGTACGTACGT", "AAAATGCATGCATGGGTGCAAAAA");
        assert_eq!(result5.matched_5_end, 14);
        assert_eq!(result5.matched_3_end, 11);
        assert_eq!(result5.max_matched_length, 10);
    }

    #[test]
    fn test_merge_variant_maps_preserves_new_variant_fields() {
        let mut dest: VariantMapByPos = Default::default();
        let mut src: VariantMapByPos = Default::default();

        let key = VarDesc::snv_key(b'A');
        let mut source_variant = Variant::default();
        source_variant.alt_depth = 10;
        source_variant.extra_cnt = 0;
        source_variant.high_qual_read_cnt = 9;
        source_variant.low_qual_read_cnt = 1;
        source_variant.mean_pos = 314.0;
        source_variant.mean_qual = 393.0;
        source_variant.mean_mapq = 46.0;
        source_variant.nm = 3.0;
        source_variant.alt_depth_fwd = 9;
        source_variant.alt_depth_rev = 1;
        source_variant.pstd = true;
        source_variant.qstd = true;
        source_variant.pp = 25;
        source_variant.pq = 42.0;

        src.entry(2550194)
            .or_default()
            .insert(key.clone(), source_variant.clone());

        VariantRealigner::merge_variant_maps(&mut dest, src);

        let merged = dest
            .get(&2550194)
            .and_then(|vars| vars.get(&key))
            .expect("merged variant");
        assert_eq!(merged.alt_depth, source_variant.alt_depth);
        assert_eq!(merged.extra_cnt, source_variant.extra_cnt);
        assert_eq!(merged.high_qual_read_cnt, source_variant.high_qual_read_cnt);
        assert_eq!(merged.low_qual_read_cnt, source_variant.low_qual_read_cnt);
        assert_eq!(merged.mean_pos, source_variant.mean_pos);
        assert_eq!(merged.mean_qual, source_variant.mean_qual);
        assert_eq!(merged.mean_mapq, source_variant.mean_mapq);
        assert_eq!(merged.nm, source_variant.nm);
        assert_eq!(merged.alt_depth_fwd, source_variant.alt_depth_fwd);
        assert_eq!(merged.alt_depth_rev, source_variant.alt_depth_rev);
        assert_eq!(merged.pstd, source_variant.pstd);
        assert_eq!(merged.qstd, source_variant.qstd);
        assert_eq!(merged.pp, source_variant.pp);
        assert_eq!(merged.pq, source_variant.pq);
    }

    #[test]
    fn test_get_ref_base_falls_back_to_original_window() {
        let realigner = VariantRealigner::new_with_context(
            Arc::new(b"TTTT".to_vec()),
            Arc::new(Default::default()),
            200,
            None,
            Vec::new(),
        )
        .with_reference_fallback(Arc::new(b"ACGT".to_vec()), 100);

        assert_eq!(realigner.get_ref_base(100), Some(b'A'));
        assert_eq!(realigner.get_ref_base(103), Some(b'T'));
        assert_eq!(realigner.get_ref_base(200), Some(b'T'));
        assert_eq!(realigner.get_ref_base(204), None);
    }

    #[test]
    fn test_load_partial_ref_coverage_skips_windows_left_of_chromosome_start() {
        let realigner = VariantRealigner::new_with_context(
            Arc::new(b"ACGT".to_vec()),
            Arc::new(Default::default()),
            1,
            Some("MT".to_string()),
            vec!["unused.bam".to_string()],
        );
        let mut data = RealignedVariationData::default();
        let key = VarDesc::snv_key(b'G');

        let mut existing = Variant::default();
        existing.alt_depth = 350;
        existing.alt_depth_fwd = 347;
        existing.alt_depth_rev = 3;
        data.non_insertion_variants
            .entry(1)
            .or_default()
            .insert(key.clone(), existing);
        data.ref_coverage.insert(1, 350);

        realigner.load_partial_ref_coverage(&mut data, -199, 0);

        let variant = data
            .non_insertion_variants
            .get(&1)
            .and_then(|vars| vars.get(&key))
            .expect("existing breakpoint reference variant");
        assert_eq!(variant.alt_depth, 350);
        assert_eq!(variant.alt_depth_fwd, 347);
        assert_eq!(variant.alt_depth_rev, 3);
        assert_eq!(data.ref_coverage.get(&1), Some(&350));
    }

    #[test]
    fn test_move_sv_marker_overwrites_destination_counts() {
        let mut data = RealignedVariationData::default();
        let key = VarDesc::Raw {
            desc: b"SV".to_vec().into(),
        };

        data.non_insertion_variants
            .entry(100)
            .or_default()
            .insert(key.clone(), Variant::default());
        data.non_insertion_variants
            .entry(99)
            .or_default()
            .insert(key.clone(), Variant::default());
        data.sv_counts.insert(
            100,
            crate::variants::variants::StructuralVariantCounts {
                pairs: 6,
                splits: 4,
                clusters: 2,
            },
        );
        data.sv_counts.insert(
            99,
            crate::variants::variants::StructuralVariantCounts {
                pairs: 1,
                splits: 1,
                clusters: 1,
            },
        );

        VariantRealigner::move_sv_marker(&mut data, 100, 99);

        assert!(
            data.non_insertion_variants
                .get(&99)
                .and_then(|map| map.get(&key))
                .is_some()
        );
        let moved = data.sv_counts.get(&99).expect("destination SV counts");
        assert_eq!(moved.pairs, 6);
        assert_eq!(moved.splits, 4);
        assert_eq!(moved.clusters, 2);
        assert!(
            data.non_insertion_variants
                .get(&100)
                .and_then(|map| map.get(&key))
                .is_none()
        );
        assert!(!data.sv_counts.contains_key(&100));
    }

    #[test]
    fn test_merge_suffix_shift_insertion_keys_keeps_independently_supported_short_suffix_key() {
        let mut pos_map: VariantMap = VecMap::default();
        pos_map.insert(
            VarDesc::Ins {
                seq: b"AA".to_vec().into(),
            },
            Variant {
                alt_depth: 6,
                alt_depth_fwd: 2,
                alt_depth_rev: 4,
                high_qual_read_cnt: 6,
                mean_pos: 132.0,
                mean_qual: 203.0,
                mean_mapq: 236.0,
                nm: 1.0,
                ..Default::default()
            },
        );
        pos_map.insert(
            VarDesc::Ins {
                seq: b"A".to_vec().into(),
            },
            Variant {
                alt_depth: 1,
                alt_depth_fwd: 1,
                high_qual_read_cnt: 1,
                mean_pos: 45.0,
                mean_qual: 42.0,
                mean_mapq: 29.0,
                ..Default::default()
            },
        );

        let mut insertion_counts: CountMap = VecMap::default();
        insertion_counts.insert("+AA".to_string(), 6);
        insertion_counts.insert("+A".to_string(), 1);

        VariantRealigner::merge_suffix_shift_insertion_keys(&mut pos_map, Some(&insertion_counts));

        let short = pos_map
            .get(&VarDesc::Ins {
                seq: b"A".to_vec().into(),
            })
            .expect("short insertion key should remain when independently supported");
        assert_eq!(short.alt_depth, 1);
        assert_eq!(short.alt_depth_fwd, 1);
        assert_eq!(short.high_qual_read_cnt, 1);

        let long = pos_map
            .get(&VarDesc::Ins {
                seq: b"AA".to_vec().into(),
            })
            .expect("long insertion key missing");
        assert_eq!(long.alt_depth, 6);
        assert_eq!(long.alt_depth_fwd, 2);
        assert_eq!(long.alt_depth_rev, 4);
        assert_eq!(long.high_qual_read_cnt, 6);
        assert_eq!(long.extra_cnt, 0);
        assert_eq!(long.mean_pos, 132.0);
        assert_eq!(long.mean_qual, 203.0);
        assert_eq!(long.mean_mapq, 236.0);
    }

    #[test]
    fn test_merge_suffix_shift_insertion_keys_merges_unsupported_short_suffix_key() {
        let mut pos_map: VariantMap = VecMap::default();
        pos_map.insert(
            VarDesc::Ins {
                seq: b"AA".to_vec().into(),
            },
            Variant {
                alt_depth: 6,
                alt_depth_fwd: 2,
                alt_depth_rev: 4,
                high_qual_read_cnt: 6,
                mean_pos: 132.0,
                mean_qual: 203.0,
                mean_mapq: 236.0,
                nm: 1.0,
                ..Default::default()
            },
        );
        pos_map.insert(
            VarDesc::Ins {
                seq: b"A".to_vec().into(),
            },
            Variant {
                alt_depth: 1,
                alt_depth_fwd: 1,
                high_qual_read_cnt: 1,
                mean_pos: 45.0,
                mean_qual: 42.0,
                mean_mapq: 29.0,
                ..Default::default()
            },
        );

        let mut insertion_counts: CountMap = VecMap::default();
        insertion_counts.insert("+AA".to_string(), 6);

        VariantRealigner::merge_suffix_shift_insertion_keys(&mut pos_map, Some(&insertion_counts));

        assert!(!pos_map.contains_key(&VarDesc::Ins {
            seq: b"A".to_vec().into(),
        }));

        let merged = pos_map
            .get(&VarDesc::Ins {
                seq: b"AA".to_vec().into(),
            })
            .expect("merged insertion key missing");
        assert_eq!(merged.alt_depth, 7);
        assert_eq!(merged.alt_depth_fwd, 3);
        assert_eq!(merged.alt_depth_rev, 4);
        assert_eq!(merged.high_qual_read_cnt, 7);
        assert_eq!(merged.extra_cnt, 0);
        assert_eq!(merged.mean_pos, 177.0);
        assert_eq!(merged.mean_qual, 245.0);
        assert_eq!(merged.mean_mapq, 265.0);
    }

    #[test]
    fn test_variant_realigner_mapped_read_no_indels() {
        let conf = crate::conf::Configuration {
            perform_local_realignment: false,
            disable_sv: true,
            ..Default::default()
        };

        let _ = INSTANCE.set(GlobalReadOnlyScope {
            conf,
            ..Default::default()
        });

        let bam_path = concat!(env!("CARGO_MANIFEST_DIR"), "/testdata/test_168714.bam");
        let fasta_path = concat!(env!("CARGO_MANIFEST_DIR"), "/testdata/hs37d5.fa");

        let mut reader = Reader::from_path(bam_path).expect("Failed to open test BAM");
        let mut target: Option<rust_htslib::bam::Record> = None;
        for result in reader.records() {
            let record = result.expect("Failed to read BAM record");
            let qname = std::str::from_utf8(record.qname()).unwrap_or("");
            if qname == "SRR098401.7003120" && !record.is_unmapped() {
                target = Some(record);
                break;
            }
        }

        let record = target.expect("Mapped SRR098401.7003120 read not found");

        let region = Region::new("20".to_string(), 168600, 168800, "test_region".to_string());
        let mut ref_start = region.start.saturating_sub(1200);
        if ref_start == 0 {
            ref_start = 1;
        }
        let ref_end = region.end + 1200;

        let mut fasta = FastaReader::open(fasta_path).expect("Failed to open reference FASTA");
        let reference = fasta
            .get_reference(region.chr(), ref_start, ref_end)
            .expect("Failed to fetch reference sequence");

        let scope_instance = instance_arc();
        let mut parser = CigarParser::new(region.clone(), reference.clone(), scope_instance);

        let mut records = vec![record];
        parser
            .process_records(records.iter_mut())
            .expect("process_records failed");

        let position_to_deletions_count = parser.take_position_to_deletions_count();

        let mut sv_input = crate::mods::structural_variants_processor::RealignedVariationData {
            non_insertion_variants: parser.take_non_insertion_vars(),
            sv_counts: Default::default(),
            insertion_variants: parser.take_insertion_vars(),
            soft_clips_5end: parser.take_soft_clips_5end(),
            soft_clips_3end: parser.take_soft_clips_3end(),
            ref_coverage: parser.take_ref_coverage(),
            splice: Default::default(),
            max_read_length: parser.get_max_read_len(),
            duprate: 0.0,
            svfdel: parser.take_svfdel(),
            svrdel: parser.take_svrdel(),
            svfdup: parser.take_svfdup(),
            svrdup: parser.take_svrdup(),
            svfinv5: parser.take_svfinv5(),
            svrinv5: parser.take_svrinv5(),
            svfinv3: parser.take_svfinv3(),
            svrinv3: parser.take_svrinv3(),
            svffus: parser.take_svffus(),
            svrfus: parser.take_svrfus(),
            softp2sv_first_used: Default::default(),
        };

        if instance_arc().conf.perform_local_realignment {
            let realigner = VariantRealigner::new_with_context(
                Arc::clone(&reference.ref_seq),
                Arc::clone(&reference.seed),
                reference.region_start,
                None,
                Vec::new(),
            );
            realigner.process_deletions(&mut sv_input, &position_to_deletions_count);
        }

        assert!(!reference.ref_seq.is_empty());
        assert!(!sv_input.non_insertion_variants.is_empty());
        assert!(sv_input.insertion_variants.is_empty());
        assert!(sv_input.soft_clips_5end.is_empty());
        assert!(sv_input.soft_clips_3end.is_empty());
    }
}
