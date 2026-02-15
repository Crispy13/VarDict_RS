use std::collections::HashMap;
use rust_htslib::bam::{record::Cigar, Record};

use crate::{
    conf::Configuration,
    data::bam_reader::BamReader,
    data::patterns::{
        AMP_ATGC, ATGSs_AMP_ATGSs_END, BEGIN_PLUS_ATGC, CARET_ATGC_END, DUP_NUM_ATGC, HASH_ATGC,
        UP_NUMBER_END,
    },
    mods::structural_variants_processor::RealignedVariationData,
    prelude::SmallVecBytes,
    variants::{
        var_utils::{find_conseq, get_variants_from_map},
        variants::{InsOrDelLen, VarDesc, Variant},
    },
};

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

pub struct VariantRealigner {
    reference_seq: Vec<u8>,
    reference_seed: HashMap<Vec<u8>, Vec<i64>>,
    ref_start: i64,
    chromosome: Option<String>,
    bam_paths: Vec<String>,
}

impl VariantRealigner {
    pub fn new(reference_seq: Vec<u8>, reference_seed: HashMap<Vec<u8>, Vec<i64>>, ref_start: i64) -> Self {
        Self::new_with_context(reference_seq, reference_seed, ref_start, None, Vec::new())
    }

    pub fn new_with_context(
        reference_seq: Vec<u8>,
        reference_seed: HashMap<Vec<u8>, Vec<i64>>,
        ref_start: i64,
        chromosome: Option<String>,
        bam_paths: Vec<String>,
    ) -> Self {
        Self {
            reference_seq,
            reference_seed,
            ref_start,
            chromosome,
            bam_paths,
        }
    }

    pub fn process_deletions(
        &self,
        data: &mut RealignedVariationData,
    ) {
        let mut del_keys: Vec<(i64, VarDesc, String)> = Vec::new();
        for (pos, var_map) in data.non_insertion_variants.iter() {
            for desc in var_map.keys() {
                if matches!(desc, VarDesc::Del { .. }) {
                    del_keys.push((*pos, desc.clone(), desc.to_key_string()));
                }
            }
        }

        del_keys.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.2.cmp(&b.2)));

        for (pos, desc, _) in &del_keys {
            self.realign_deletion_mismatches(*pos, desc, data);
        }

        // Java realigndel post-pass:
        // for i = tmp.size()-1; i > 0; i-- {
        //   if vn =~ /^(-\d+)&[ATGC]+$/ and vars(vn) < vars($1) then merge vn into $1
        // }
        for idx in (1..del_keys.len()).rev() {
            let (pos, desc, desc_str) = &del_keys[idx];
            let Some(base_del_desc) = Self::minus_amp_base_desc(desc_str) else {
                continue;
            };
            let Some(base_key) = self.del_desc_to_key(base_del_desc.as_bytes()) else {
                continue;
            };

            let Some(pos_map) = data.non_insertion_variants.get_mut(pos) else {
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

    fn has_sv_marker(variation_map: &HashMap<VarDesc, Variant>) -> bool {
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
    }

    fn add_sv_split_count(data: &mut RealignedVariationData, position: i64, splits: usize) {
        Self::ensure_sv_marker(data, position);
        if splits == 0 {
            return;
        }
        if let Some(variation_map) = data.non_insertion_variants.get_mut(&position) {
            let key = VarDesc::Raw {
                desc: b"SV".to_vec().into(),
            };
            let sv = variation_map.entry(key).or_default();
            sv.alt_depth += splits;
        }
    }

    pub fn process_insertions(
        &self,
        data: &mut RealignedVariationData,
        position_to_insertion_count: &std::collections::HashMap<i64, std::collections::HashMap<String, usize>>,
    ) {
        if position_to_insertion_count.is_empty() {
            return;
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
            let Some(insert_match) = BEGIN_PLUS_ATGC.captures(&vn) else { continue; };
            let insert = insert_match.get(1).map(|m| m.as_str()).unwrap_or("");
            if insert.is_empty() {
                continue;
            }

            let mut ins3 = String::new();
            let mut inslen = insert.len();
            if let Some(cap) = DUP_NUM_ATGC.captures(&vn) {
                let dup_len: usize = cap.get(1).map(|m| m.as_str()).unwrap_or("0").parse().unwrap_or(0);
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
                let p3 = position + inslen as i64 - ins3.len() as i64 + Configuration::SVFLANK as i64;
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
                seq: vn[1..].as_bytes().iter().map(|b| b.to_ascii_uppercase()).collect(),
            };
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
                    VarDesc::SNV { ref_base: mm_bytes[0] }
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
                let Some(tv) = tv_snapshot else { continue; };

                if tv.alt_depth == 0 {
                    continue;
                }
                let mean_qual = tv.mean_qual / tv.alt_depth as f64;
                if mean_qual < crate::scopedata::global_read_only_scope::instance().conf.goodq {
                    continue;
                }
                let nm_threshold = if mm.end == 3 { r3.nm } else { r5.nm };
                let mean_pos = tv.mean_pos / tv.alt_depth as f64;
                if mean_pos > nm_threshold as f64 + 4.0 {
                    continue;
                }
                if tv.alt_depth >= insertion_count + insert.len() || (tv.alt_depth as f64 / insertion_count as f64) >= 8.0 {
                    continue;
                }

                if mm.mismatch_position > position && mm.end == 5 {
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
                let Some(tv_owned) = tv_owned else { continue; };

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
                    .and_then(|m| m.get_mut(&vn_key))
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
                        if map.is_empty() {
                            data.non_insertion_variants.remove(&r3.misp);
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
                        if map.is_empty() {
                            data.non_insertion_variants.remove(&r5.misp);
                        }
                    }
                }
            }

            for sc5pp in r5.scp.iter().copied() {
                if let Some(tv) = data.soft_clips_5end.get_mut(&sc5pp) {
                    if tv.used() {
                        continue;
                    }
                    let seq = find_conseq(tv, 0);
                    if !seq.is_empty() && Self::is_match_bytes(&seq, &wupseq, -1) {
                        if sc5pp > position {
                            *data.ref_coverage.entry(position).or_insert(0) += tv.var.alt_depth;
                        }
                        if let Some(vref) = data
                            .insertion_variants
                            .get_mut(&position)
                            .and_then(|m| m.get_mut(&vn_key))
                        {
                            adj_cnt(vref, &tv.var);
                        }
                        tv.mark_used();
                    }
                }
            }

            for sc3pp in r3.scp.iter().copied() {
                if let Some(tv) = data.soft_clips_3end.get_mut(&sc3pp) {
                    if tv.used() {
                        continue;
                    }
                    let seq = find_conseq(tv, 0);
                    let mseq = if !ins3.is_empty() {
                        sanpseq.clone()
                    } else {
                        let offset = sc3pp - position - 1;
                        Self::substr_bytes(&sanpseq, offset, None)
                    };
                    if !seq.is_empty() && Self::is_match_bytes(&seq, &mseq, 1) {
                        let mean_pos = if tv.var.alt_depth > 0 {
                            tv.var.mean_pos / tv.var.alt_depth as f64
                        } else {
                            0.0
                        };
                        if sc3pp <= position || insert.len() as f64 > mean_pos {
                            *data.ref_coverage.entry(position).or_insert(0) += tv.var.alt_depth;
                        }

                        let use_ref_var = sc3pp > position
                            && !(insert.len() as f64 > mean_pos);

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
                            .and_then(|m| m.get_mut(&vn_key))
                        {
                            if let Some((ref_key, mut ref_var)) = ref_var {
                                adj_cnt_with_ref(vref, &tv.var, Some(&mut ref_var));
                                if let Some(pos_map) = data.non_insertion_variants.get_mut(&position) {
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
                                    let old_key = VarDesc::Ins {
                                        seq: insert.as_bytes().iter().map(|b| b.to_ascii_uppercase()).collect(),
                                    };
                                    if let Some(variation) = map.remove(&old_key) {
                                        let new_key = VarDesc::Ins { seq: tvn_bytes[1..].iter().map(|b| b.to_ascii_uppercase()).collect() };
                                        map.insert(new_key, variation);
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
                if first3 > first5 + 3 && first3 - first5 < (data.max_read_length as f64 * 0.75) as i64 {
                    if let Some(ref_base) = self.get_ref_base(position) {
                        if let Some(map) = data.non_insertion_variants.get_mut(&position) {
                            if let Some(ref_var) = map.get_mut(&VarDesc::SNV { ref_base }) {
                                Self::adj_ref_factor(ref_var, (first3 - first5 - 1) as f64 / data.max_read_length as f64);
                            }
                        }
                    }
                    if let Some(vref) = data
                        .insertion_variants
                        .get_mut(&position)
                        .and_then(|m| m.get_mut(&vn_key))
                    {
                        Self::adj_ref_factor(vref, -((first3 - first5 - 1) as f64 / data.max_read_length as f64));
                    }
                }
            }
        }

        for idx in (1..tmp.len()).rev() {
            let (position, vn, _) = tmp[idx].clone();
            let Some(map) = data.insertion_variants.get_mut(&position) else { continue; };
            let vref_key = if vn.starts_with('+') {
                VarDesc::Ins { seq: vn[1..].as_bytes().iter().map(|b| b.to_ascii_uppercase()).collect() }
            } else {
                continue;
            };
            let vref_alt = map.get(&vref_key).map(|v| v.alt_depth);
            if let Some(cap) = ATGSs_AMP_ATGSs_END.captures(&vn) {
                let tn = cap.get(1).map(|m| m.as_str()).unwrap_or("");
                if !tn.is_empty() {
                    let tref_key = VarDesc::Ins { seq: tn[1..].as_bytes().iter().map(|b| b.to_ascii_uppercase()).collect() };
                    let tref_alt = map.get(&tref_key).map(|v| v.alt_depth);
                    if let (Some(vref_alt), Some(tref_alt)) = (vref_alt, tref_alt) {
                        if vref_alt < tref_alt {
                            let vref = map.remove(&vref_key);
                            if let Some(vref) = vref {
                                if let Some(tref) = map.get_mut(&tref_key) {
                                    let ref_var = self.get_ref_base(position)
                                        .map(|ref_base| VarDesc::SNV { ref_base })
                                        .and_then(|ref_key| {
                                            data.non_insertion_variants
                                                .get_mut(&position)
                                                .and_then(|m| m.remove(&ref_key))
                                                .map(|var| (ref_key, var))
                                        });

                                    if let Some((ref_key, mut ref_var)) = ref_var {
                                        adj_cnt_with_ref(tref, &vref, Some(&mut ref_var));
                                        if let Some(pos_map) = data.non_insertion_variants.get_mut(&position) {
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
                    let Some(sc5v) = data.soft_clips_5end.get_mut(&p5) else { continue; };
                    find_conseq(sc5v, 5)
                };
                let seq3 = {
                    let Some(sc3v) = data.soft_clips_3end.get_mut(&p3) else { continue; };
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
                        let ref_seq = self.join_ref(
                            p5,
                            p5 + seq3.len() as i64 - ins_desc.len() as i64 + 2,
                        );
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
                        let ref_seq = self.join_ref(
                            p5,
                            p5 + seq3.len() as i64 - ins_desc.len() as i64 + 2,
                        );
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
                        let to = p5 as f64
                            + (p3 - p5 + ins_desc.len() as i64) as f64 / rpt as f64
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

                *data.ref_coverage.entry(bi).or_insert(0) += sc5_var.alt_depth;

                let ins_starts_with_plus = ins_desc.first() == Some(&b'+');
                let ins_starts_with_minus = ins_desc.first() == Some(&b'-');

                if ins_starts_with_plus {
                    let key = VarDesc::Ins {
                        seq: ins_desc[1..]
                            .iter()
                            .map(|b| b.to_ascii_uppercase())
                            .collect(),
                    };
                    let vref = get_variants_from_map(&mut data.insertion_variants, bi, &key);
                    vref.pstd = true;
                    vref.qstd = true;

                    let ref_key = self.get_ref_base(bi).map(|ref_base| VarDesc::SNV { ref_base });
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

                    let mut tins = HashMap::new();
                    let mut map = HashMap::new();
                    map.insert(
                        String::from_utf8_lossy(&ins_desc).to_string(),
                        vref.alt_depth,
                    );
                    tins.insert(bi, map);
                    self.process_insertions(data, &tins);
                } else if ins_starts_with_minus {
                    let vref_key = self
                        .del_desc_to_key(&ins_desc)
                        .unwrap_or_else(|| VarDesc::Raw { desc: ins_desc.clone().into() });
                    if let Some(pos_map) = data.non_insertion_variants.get_mut(&bi) {
                        let ref_key =
                            self.get_ref_base(bi).map(|ref_base| VarDesc::SNV { ref_base });
                        let mut ref_var = ref_key.as_ref().and_then(|k| pos_map.remove(k));
                        let vref = pos_map.entry(vref_key.clone()).or_default();
                        vref.pstd = true;
                        vref.qstd = true;
                        if let Some(ref_var_mut) = ref_var.as_mut() {
                            adj_cnt_with_ref(vref, &sc3_var, Some(ref_var_mut));
                        } else {
                            adj_cnt(vref, &sc3_var);
                        }
                        adj_cnt(vref, &sc5_var);
                        if let (Some(ref_key), Some(ref_var)) = (ref_key, ref_var) {
                            pos_map.insert(ref_key, ref_var);
                        }
                    }
                    self.realign_deletion_for_desc(bi, &vref_key, data);
                } else {
                    let vref_key = VarDesc::Raw {
                        desc: ins_desc.clone().into(),
                    };
                    let vref = get_variants_from_map(&mut data.non_insertion_variants, bi, &vref_key);
                    vref.pstd = true;
                    vref.qstd = true;
                    adj_cnt(vref, &sc3_var);
                    adj_cnt(vref, &sc5_var);
                }

                break;
            }
        }
    }

    pub fn realign_long_insertions(&self, data: &mut RealignedVariationData) {
        let conf = &crate::scopedata::global_read_only_scope::instance().conf;

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
                let Some(sc5v) = data.soft_clips_5end.get_mut(&p) else { continue; };
                find_conseq(sc5v, 0)
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
                if !(bi != 0
                    && bi - p > 15
                    && bi - p < Configuration::SVMAXLEN as i64)
                {
                    continue;
                }

                if bi - p
                    > conf.sv_min_len as i64 + 2 * Configuration::SVFLANK as i64
                {
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
                    *data.ref_coverage.entry(p - 1).or_insert(0) += cnt;
                }

                bi = p - 1;
                Self::add_sv_split_count(data, bi, cnt);
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
                seq: ins.iter().map(|b| b.to_ascii_uppercase()).collect(),
            };
            let iref = get_variants_from_map(&mut data.insertion_variants, bi, &iref_key);
            iref.pstd = true;
            iref.qstd = true;
            adj_cnt(iref, &sc5_var);

            if let Some(variation_map) = data.non_insertion_variants.get(&bi) {
                if !Self::has_sv_marker(variation_map) {
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
                let Some(map) = sc5_seq.get(&ii) else { continue; };
                for base in [b'A', b'C', b'G', b'T', b'N'] {
                    let Some(tv) = map.get(base) else { continue; };
                    let key = VarDesc::SNV { ref_base: base };
                    let tvr = get_variants_from_map(&mut data.non_insertion_variants, pii, &key);
                    adj_cnt(tvr, tv);
                    tvr.pstd = true;
                    tvr.qstd = true;
                    *data.ref_coverage.entry(pii).or_insert(0) += tv.alt_depth;
                }
            }

            if bi + len as i64 != 0 {
                if let Some(sc5v) = data.soft_clips_5end.get_mut(&p) {
                    sc5v.mark_used();
                }
            }

            let mut tins = HashMap::new();
            let mut map = HashMap::new();
            map.insert(
                format!("+{}", String::from_utf8_lossy(&ins)),
                iref.alt_depth,
            );
            tins.insert(bi, map);
            self.process_insertions(data, &tins);

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
                let Some(sc3v) = data.soft_clips_3end.get_mut(&position) else { continue; };
                find_conseq(sc3v, 0)
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
                if !(bi != 0
                    && p - bi > 15
                    && p - bi < Configuration::SVMAXLEN as i64)
                {
                    continue;
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

                if p - bi
                    > conf.sv_min_len as i64 + 2 * Configuration::SVFLANK as i64
                {
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

                bi -= 1;
                Self::add_sv_split_count(data, bi, cnt);

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
                seq: ins.iter().map(|b| b.to_ascii_uppercase()).collect(),
            };
            let iref = get_variants_from_map(&mut data.insertion_variants, bi, &iref_key);
            iref.pstd = true;
            iref.qstd = true;

            let mean_pos = if cnt > 0 {
                sc3_var.mean_pos / cnt as f64
            } else {
                0.0
            };
            let mut ref_key = self.get_ref_base(bi).map(|ref_base| VarDesc::SNV { ref_base });
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

            let mut len = ins.len();
            if ins.iter().any(|b| *b == b'&') {
                len = len.saturating_sub(1);
            }
            let seq_len = sc3_seq.keys().next_back().map(|k| *k + 1).unwrap_or(0);

            for ii in len..seq_len {
                let pii = p + ii as i64 - len as i64;
                let Some(map) = sc3_seq.get(&ii) else { continue; };
                for base in [b'A', b'C', b'G', b'T', b'N'] {
                    let Some(tv) = map.get(base) else { continue; };
                    let key = VarDesc::SNV { ref_base: base };
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

            let mut tins = HashMap::new();
            let mut map = HashMap::new();
            map.insert(
                format!("+{}", String::from_utf8_lossy(&ins)),
                iref.alt_depth,
            );
            tins.insert(bi, map);
            self.process_insertions(data, &tins);

        }
    }

    /// Adjust MNPs (multi-nucleotide polymorphisms) when there are breakpoints within MNP
    pub fn adjust_mnp(
        &self,
        data: &mut RealignedVariationData,
        mnp: &std::collections::HashMap<i64, std::collections::HashMap<String, usize>>,
    ) {
        let mut tmp: Vec<(i64, String)> = Vec::new();
        for (pos, desc_map) in mnp {
            for desc in desc_map.keys() {
                tmp.push((*pos, desc.clone()));
            }
        }
        tmp.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));

        for (position, vn) in tmp {
            let vref_key = VarDesc::Raw { desc: vn.as_bytes().to_vec().into() };
            let vref_cnt = match data
                .non_insertion_variants
                .get(&position)
                .and_then(|m| m.get(&vref_key))
            {
                Some(v) => v.alt_depth,
                None => continue,
            };

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
                        if tref.alt_depth > 0
                            && tref.alt_depth < vref_cnt
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
                        if tref.alt_depth < vref_cnt {
                            if let Some(vars_on_pos) = data.non_insertion_variants.get_mut(&position) {
                                if let Some(vref) = vars_on_pos.get_mut(&vref_key) {
                                    adj_cnt(vref, &tref);
                                }
                            }
                            *data.ref_coverage.entry(position).or_insert(0) += tref.alt_depth;

                            if let Some(vars_right) = data.non_insertion_variants.get_mut(&right_pos) {
                                vars_right.remove(&right_key);
                                if vars_right.is_empty() {
                                    data.non_insertion_variants.remove(&right_pos);
                                }
                            }
                        }
                    }
            }

            if let Some(sc3v) = data.soft_clips_3end.get_mut(&position) {
                if !sc3v.used() {
                    let seq = find_conseq(sc3v, 0);
                    if seq.starts_with(mnt.as_bytes()) {
                        let tail = &seq[mnt.len()..];
                        if tail.is_empty()
                            || self.is_match_ref(tail, position + mnt.len() as i64, 1, 3)
                        {
                            if let Some(vars_on_pos) = data.non_insertion_variants.get_mut(&position) {
                                if let Some(vref) = vars_on_pos.get_mut(&vref_key) {
                                    adj_cnt(vref, &sc3v.var);
                                }
                            }
                            *data.ref_coverage.entry(position).or_insert(0) += sc3v.var.alt_depth;
                            sc3v.mark_used();
                        }
                    }
                }
            }

            let pos_5end = position + mnt.len() as i64;
            if let Some(sc5v) = data.soft_clips_5end.get_mut(&pos_5end) {
                if !sc5v.used() {
                    let mut seq = find_conseq(sc5v, 0);
                    if !seq.is_empty() && seq.len() >= mnt.len() {
                        seq.reverse();
                        if seq.ends_with(mnt.as_bytes()) {
                            let prefix = &seq[..seq.len() - mnt.len()];
                            if prefix.is_empty() || self.is_match_ref(prefix, position - 1, -1, 3) {
                                if let Some(vars_on_pos) = data.non_insertion_variants.get_mut(&position) {
                                    if let Some(vref) = vars_on_pos.get_mut(&vref_key) {
                                        adj_cnt(vref, &sc5v.var);
                                    }
                                }
                                *data.ref_coverage.entry(position).or_insert(0) += sc5v.var.alt_depth;
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
        data: &mut RealignedVariationData,
    ) {
        let (dellen, extra_seq, extrains_len) = match desc {
            VarDesc::Del {
                len,
                match_seq,
                ins_or_del_len,
                mismatch_seq,
                ..
            } => {
                let mut total = *len as i64;
                if let InsOrDelLen::DelLen(extra) = ins_or_del_len {
                    total += *extra as i64;
                }

                let mut extra_bytes = Vec::new();
                if !match_seq.is_empty() {
                    extra_bytes.extend_from_slice(match_seq.as_slice());
                }
                match ins_or_del_len {
                    InsOrDelLen::InsSeq(seq) => {
                        extra_bytes.extend_from_slice(seq.as_slice());
                    }
                    InsOrDelLen::DelLen(del_len) => {
                        extra_bytes.extend_from_slice(del_len.to_string().as_bytes());
                    }
                    InsOrDelLen::None => {}
                }
                if !mismatch_seq.is_empty() {
                    extra_bytes.extend_from_slice(mismatch_seq.as_slice());
                }

                let extra = String::from_utf8_lossy(&extra_bytes).to_string();
                let extrains_len = match ins_or_del_len {
                    InsOrDelLen::InsSeq(seq) => seq.len() as i64,
                    _ => 0,
                };
                (total, extra, extrains_len)
            }
            _ => return,
        };

        let mut wupseq = self.get_ref_range(pos - 200, pos - 1);
        wupseq.extend_from_slice(extra_seq.as_bytes());

        let san_start = pos + dellen + extra_seq.len() as i64 - extrains_len;
        let mut sanpseq = extra_seq.as_bytes().to_vec();
        sanpseq.extend_from_slice(&self.get_ref_range(san_start, san_start + 200));

        let r3 = self.find_mm3(pos, &sanpseq, data);
        let r5 = self.find_mm5(pos + dellen + extra_seq.len() as i64 - extrains_len - 1, &wupseq, data);

        let dcnt = match data
            .non_insertion_variants
            .get(&pos)
            .and_then(|m| m.get(desc))
        {
            Some(vref) => vref.alt_depth,
            None => return,
        };
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
            let Some(tv) = tv_snapshot else { continue; };

            if tv.alt_depth == 0 {
                continue;
            }

            let mean_qual = tv.mean_qual / tv.alt_depth as f64;
            if mean_qual < crate::scopedata::global_read_only_scope::instance().conf.goodq {
                continue;
            }

            let nm_threshold = if mm.end == 3 { r3.nm } else { r5.nm };
            let mean_pos = tv.mean_pos / tv.alt_depth as f64;
            if mean_pos > nm_threshold as f64 + 4.0 {
                continue;
            }

            if tv.alt_depth >= dcnt + dellen as usize || (tv.alt_depth as f64 / dcnt as f64) >= 8.0 {
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
            let Some(tv_owned) = tv_owned else { continue; };

            if let Some(map) = data.non_insertion_variants.get_mut(&mm.mismatch_position) {
                if map.is_empty() {
                    data.non_insertion_variants.remove(&mm.mismatch_position);
                }
            }

            let ref_key = if mm.mismatch_position > pos && mm.end == 3 {
                self.get_ref_base(pos).map(|ref_base| VarDesc::SNV { ref_base })
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
                    if map.is_empty() {
                        data.non_insertion_variants.remove(&r3.misp);
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
                    if map.is_empty() {
                        data.non_insertion_variants.remove(&r5.misp);
                    }
                }
            }
        }

        for sc5pp in r5.scp {
            if let Some(tv) = data.soft_clips_5end.get_mut(&sc5pp) {
                if tv.used() {
                    continue;
                }
                if dcnt <= 2 && tv.var.alt_depth / dcnt > 5 {
                    continue;
                }
                let seq = find_conseq(tv, 0);
                if seq.is_empty() {
                    continue;
                }
                if Self::is_match_bytes(&seq, &wupseq, -1) {
                    if sc5pp > pos {
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
            if let Some(tv) = data.soft_clips_3end.get_mut(&sc3pp) {
                if tv.used() {
                    continue;
                }
                if dcnt <= 2 && tv.var.alt_depth / dcnt > 5 {
                    continue;
                }
                let seq = find_conseq(tv, 0);
                if seq.is_empty() {
                    continue;
                }
                let offset = if sc3pp > pos {
                    (sc3pp - pos) as usize
                } else {
                    0
                };
                let mseq = if offset <= sanpseq.len() {
                    &sanpseq[offset..]
                } else {
                    &[]
                };
                if Self::is_match_bytes(&seq, mseq, 1) {
                    if sc3pp <= pos {
                        *data.ref_coverage.entry(pos).or_insert(0) += tv.var.alt_depth;
                    }

                    let ref_key = if sc3pp > pos {
                        self.get_ref_base(pos).map(|ref_base| VarDesc::SNV { ref_base })
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
        if !self.bam_paths.is_empty()
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
            let Some(sclip) = data.soft_clips_5end.get_mut(&sc_pos) else { continue; };
            if sclip.used() {
                continue;
            }

            let seq = find_conseq(sclip, 0);
            if seq.is_empty() {
                continue;
            }

            if !Self::is_match_bytes(&seq, wupseq, -1) {
                continue;
            }

            if let Some(var_map) = data.non_insertion_variants.get_mut(&pos) {
                if let Some(variant) = var_map.get_mut(desc) {
                    adj_cnt(variant, &sclip.var);
                }
            }
            if sc_pos > pos {
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
            let Some(sclip) = data.soft_clips_3end.get_mut(&sc_pos) else { continue; };
            if sclip.used() {
                continue;
            }

            let seq = find_conseq(sclip, 0);
            if seq.is_empty() {
                continue;
            }

            if !Self::is_match_bytes(&seq, sanpseq, 1) {
                continue;
            }

            if sc_pos <= pos {
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
                while let Some(ch) = char_at_neg(&seq, (n + 1 + n2) as isize) {
                    if self.is_has_and_equals(position - n as i64 - 1 - n2 as i64, ch) {
                        n2 += 1;
                    } else {
                        break;
                    }
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
            if ref_base.is_none() || ref_base.unwrap() == seq[n] {
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
                while n + n2 + 1 < seq.len() && self.is_has_and_equals(p + n as i64 + 1 + n2 as i64, seq[n + n2 + 1]) {
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

    fn join_ref_for_5_lgins(
        &self,
        from: i64,
        to: i64,
        seq: &[u8],
        extra: &[u8],
    ) -> Vec<u8> {
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
                            if seq_ch.is_some() && ref_ch.is_some() && seq_ch.unwrap() == ref_ch.unwrap() {
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

    fn find_match(&self, seq: &[u8], _position: i64, dir: i64, seed_len: usize, mm: usize) -> Match {
        let mut seq_work = seq.to_vec();
        if dir == -1 {
            seq_work.reverse();
        }

        if seq_work.len() < seed_len {
            return Match {
                base_position: 0,
                matched_sequence: Vec::new(),
            };
        }

        for i in (0..=seq_work.len() - seed_len).rev() {
            let seed = &seq_work[i..i + seed_len];
            let Some(seeds) = self.reference_seed.get(seed) else { continue; };
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
                    let Some(ch) = Self::char_at(&seq_work, mm_idx) else { break; };
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
                return Match {
                    base_position: bp,
                    matched_sequence: extra,
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
                        let Some(ch0) = sseq.get(0).copied() else { continue; };
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
                        let Some(ch_last) = Self::char_at(&sseq, -1) else { continue; };
                        if self.is_has_and_not_equals(bp, ch_last) {
                            continue;
                        }
                        eqcnt += 1;
                        let Some(ch_prev) = Self::char_at(&sseq, -2) else { continue; };
                        if self.is_has_and_not_equals(bp - 1, ch_prev) {
                            continue;
                        }
                        Self::substr_bytes(&seq_work, -(ii as i64), None)
                    };

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
        let len = std::str::from_utf8(&desc[1..idx]).ok()?.parse::<u32>().ok()?;
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
        self.realign_deletion_mismatches(pos, desc, data);
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

                if Self::record_contains_exact_deletion_len(&record, dlen) {
                    continue;
                }

                let read_start = record.pos() + 1;
                let read_end = read_start + Self::aligned_length_excluding_softclips_and_insertions(&record) as i64;

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

    fn record_contains_exact_deletion_len(record: &Record, dlen: u32) -> bool {
        record
            .cigar()
            .iter()
            .any(|op| matches!(*op, Cigar::Del(len) if len == dlen))
    }

    fn aligned_length_excluding_softclips_and_insertions(record: &Record) -> u32 {
        let mut len = 0u32;
        for op in record.cigar().iter() {
            match *op {
                Cigar::Match(l)
                | Cigar::Equal(l)
                | Cigar::Diff(l)
                | Cigar::Del(l)
                | Cigar::RefSkip(l) => {
                    len += l;
                }
                _ => {}
            }
        }
        len
    }

    fn get_ref_base(&self, pos: i64) -> Option<u8> {
        if pos < self.ref_start {
            return None;
        }
        let idx = (pos - self.ref_start) as usize;
        self.reference_seq.get(idx).copied()
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
                    let is_mismatch = if seq3_pos >= seq3_chars.len() || seq5_pos >= seq5_chars.len() {
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
                let is_end_match = (total_length + j >= seq3_chars.len()) || (i + total_length >= seq5_chars.len());
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
        ref_var.high_qual_read_cnt = (ref_var.high_qual_read_cnt as i64 - (factor * ref_var.high_qual_read_cnt as f64) as i64).max(0) as usize;
        ref_var.low_qual_read_cnt = (ref_var.low_qual_read_cnt as i64 - (factor * ref_var.low_qual_read_cnt as f64) as i64).max(0) as usize;

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
        ref_var.alt_depth_fwd = (ref_var.alt_depth_fwd as i64 - (factor * ref_var.alt_depth_fwd as f64) as i64).max(0) as usize;
        ref_var.alt_depth_rev = (ref_var.alt_depth_rev as i64 - (factor * ref_var.alt_depth_rev as f64) as i64).max(0) as usize;

        correct_cnt(ref_var);
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

    let Some(reference) = reference else { return; };

    reference.alt_depth = reference.alt_depth.saturating_sub(src.alt_depth);
    reference.high_qual_read_cnt = reference.high_qual_read_cnt.saturating_sub(src.high_qual_read_cnt);
    reference.low_qual_read_cnt = reference.low_qual_read_cnt.saturating_sub(src.low_qual_read_cnt);
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
    if var.nm < 0.0 {
        var.nm = 0.0;
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
    use crate::scopedata::global_read_only_scope::{GlobalReadOnlyScope, INSTANCE, instance};
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
        assert!(VariantRealigner::is_match("ACGTACGTACGTACGT", "AAGTACTTACGTACGT", 1));
        assert!(VariantRealigner::is_match("ACGTACGTACGTACGT", "AAGTACTTACGT", 1));

        assert!(!VariantRealigner::is_match("ACGTACGTACGT", "AAGTACTTACGT", 1));
        assert!(!VariantRealigner::is_match("ACGTCAGCAT", "ACGACTGACT", 1));
    }

    #[test]
    fn test_find_35_match() {
        // Test cases from Java VariationRealignerTest
        let result1 = VariantRealigner::find_35_match("ACGTACGTACGTACGTACGTACGT", "TGCATGCATGCATGCATGCATGCA");
        assert_eq!(result1.matched_5_end, 0);
        assert_eq!(result1.matched_3_end, 1);
        assert_eq!(result1.max_matched_length, 24);

        let result2 = VariantRealigner::find_35_match("ACGTACGTACGTACGTACGTACGT", "TGCATGCATGCCTGCATGCATGCA");
        assert_eq!(result2.matched_5_end, 0);
        assert_eq!(result2.matched_3_end, 1);
        assert_eq!(result2.max_matched_length, 23);

        let result3 = VariantRealigner::find_35_match("ACGTACGTACGTACGTACGTACGT", "TGCATGCATGCATGGGTGCATGCA");
        assert_eq!(result3.matched_5_end, 0);
        assert_eq!(result3.matched_3_end, 1);
        assert_eq!(result3.max_matched_length, 22);

        let result4 = VariantRealigner::find_35_match("ACGTACGTACGTACGTACGTACGT", "TGCATGCATGCATGGGTGCAAAAA");
        assert_eq!(result4.matched_5_end, 0);
        assert_eq!(result4.matched_3_end, 13);
        assert_eq!(result4.max_matched_length, 12);

        let result5 = VariantRealigner::find_35_match("ACGTACGTACGTACGTACGTACGT", "AAAATGCATGCATGGGTGCAAAAA");
        assert_eq!(result5.matched_5_end, 14);
        assert_eq!(result5.matched_3_end, 11);
        assert_eq!(result5.max_matched_length, 10);
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

        let bam_path = "/home/eck/workspace/vardict_rs/test_data/test_168714.bam";
        let fasta_path = "/home/eck/workspace/vardict_rs/VarDictJava/tests/integration/reference/hs37d5.fa";

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

        let fasta = FastaReader::open(fasta_path).expect("Failed to open reference FASTA");
        let reference = fasta
            .get_reference(region.chr(), ref_start, ref_end)
            .expect("Failed to fetch reference sequence");

        let scope_instance = Arc::new(instance().clone());
        let mut parser = CigarParser::new(region.clone(), reference.clone(), scope_instance);

        let mut records = vec![record];
        parser
            .process_records(records.iter_mut())
            .expect("process_records failed");

        let mut sv_input = crate::mods::structural_variants_processor::RealignedVariationData {
            non_insertion_variants: parser.take_non_insertion_vars(),
            insertion_variants: parser.take_insertion_vars(),
            soft_clips_5end: parser.take_soft_clips_5end(),
            soft_clips_3end: parser.take_soft_clips_3end(),
            ref_coverage: parser.take_ref_coverage(),
            max_read_length: parser.get_max_read_len(),
            duprate: 0.0,
        };

        if instance().conf.perform_local_realignment {
            let realigner = VariantRealigner::new(
                reference.ref_seq.clone(),
                reference.seed.clone(),
                reference.region_start,
            );
            realigner.process_deletions(&mut sv_input);
        }

        assert!(!reference.ref_seq.is_empty());
        assert!(!sv_input.non_insertion_variants.is_empty());
        assert!(sv_input.insertion_variants.is_empty());
        assert!(sv_input.soft_clips_5end.is_empty());
        assert!(sv_input.soft_clips_3end.is_empty());
    }
}
