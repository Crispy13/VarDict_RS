use crate::{
    conf::Configuration,
    data::patterns::{
        AMP_ATGC, ATGSs_AMP_ATGSs_END, BEGIN_PLUS_ATGC, CARET_ATGC_END, DUP_NUM_ATGC, HASH_ATGC,
        UP_NUMBER_END,
    },
    mods::structural_variants_processor::RealignedVariationData,
    variants::{
        var_utils::find_conseq,
        variants::{InsOrDelLen, SoftClip, VarDesc, Variant},
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
    ref_start: i64,
}

impl VariantRealigner {
    pub fn new(reference_seq: Vec<u8>, ref_start: i64) -> Self {
        Self {
            reference_seq,
            ref_start,
        }
    }

    pub fn process_deletions(
        &self,
        data: &mut RealignedVariationData,
    ) {
        // Collect deletion keys to avoid borrow checker issues while mutating maps
        let mut del_keys: Vec<(i64, VarDesc)> = Vec::new();
        for (pos, var_map) in data.non_insertion_variants.iter() {
            for desc in var_map.keys() {
                if matches!(desc, VarDesc::Del { .. }) {
                    del_keys.push((*pos, desc.clone()));
                }
            }
        }

        for (pos, desc) in del_keys {
            let dellen = match &desc {
                VarDesc::Del { len, ins_or_del_len, .. } => {
                    let mut total = *len as i64;
                    if let crate::variants::variants::InsOrDelLen::DelLen(extra) = ins_or_del_len {
                        total += *extra as i64;
                    }
                    total
                }
                _ => 0,
            };

            let wupseq = self.get_ref_range(pos - 200, pos - 1);
            let sanpseq = self.get_ref_range(pos + dellen, pos + dellen + 200);

            self.realign_with_softclips_5end(pos, &desc, data, &wupseq);
            self.realign_with_softclips_3end(pos, &desc, data, &sanpseq);

            // Run mismatch-based realignment similar to Java realigndel
            self.realign_deletion_mismatches(pos, &desc, data);
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
        tmp.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));

        for (position, vn, insertion_count) in tmp.iter().cloned() {
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

            let insert_key = VarDesc::Ins {
                seq: insert.as_bytes().iter().map(|b| b.to_ascii_uppercase()).collect(),
            };
            let vref = crate::variants::var_utils::get_variants_from_map(
                &mut data.insertion_variants,
                position,
                &insert_key,
            ) as *mut Variant;

            for mm in all_mm {
                let mut mm_bytes: Vec<u8> = mm
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

                let ref_key = if mm.mismatch_position > position && mm.end == 3 {
                    self.get_ref_base(position).map(|ref_base| VarDesc::SNV { ref_base })
                } else {
                    None
                };

                unsafe {
                    if let Some(ref_key) = ref_key {
                        if let Some(pos_map) = data.non_insertion_variants.get_mut(&position) {
                            let mut ref_var = pos_map.remove(&ref_key);
                            if let Some(vref) = vref.as_mut() {
                                if let Some(ref_var_mut) = ref_var.as_mut() {
                                    adj_cnt_with_ref(vref, &tv_owned, Some(ref_var_mut));
                                } else {
                                    adj_cnt(vref, &tv_owned);
                                }
                            }
                            if let Some(ref_var) = ref_var {
                                pos_map.insert(ref_key, ref_var);
                            }
                        }
                    } else if let Some(vref) = vref.as_mut() {
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
                        unsafe {
                            if let Some(vref) = vref.as_mut() {
                                adj_cnt(vref, &tv.var);
                            }
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
                        if sc3pp <= position || insert.len() as f64 > (tv.var.mean_pos / tv.var.alt_depth as f64) {
                            *data.ref_coverage.entry(position).or_insert(0) += tv.var.alt_depth;
                        }

                        let mut lref = None;
                        if sc3pp > position {
                            if let Some(ref_base) = self.get_ref_base(position) {
                                lref = data
                                    .non_insertion_variants
                                    .get_mut(&position)
                                    .and_then(|m| m.get_mut(&VarDesc::SNV { ref_base }));
                            }
                        }
                        if insert.len() as f64 > (tv.var.mean_pos / tv.var.alt_depth as f64) {
                            lref = None;
                        }

                        unsafe {
                            if let Some(vref) = vref.as_mut() {
                                if let Some(ref_var) = lref.as_deref_mut() {
                                    adj_cnt_with_ref(vref, &tv.var, Some(ref_var));
                                } else {
                                    adj_cnt(vref, &tv.var);
                                }
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
                    unsafe {
                        if let Some(vref) = vref.as_mut() {
                            Self::adj_ref_factor(vref, -((first3 - first5 - 1) as f64 / data.max_read_length as f64));
                        }
                    }
                }
            }
        }

        let mut tmp2: Vec<(i64, String)> = Vec::new();
        for (pos, desc_map) in position_to_insertion_count {
            for desc in desc_map.keys() {
                tmp2.push((*pos, desc.clone()));
            }
        }
        tmp2.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));

        for (position, vn) in tmp2.into_iter().rev() {
            let Some(map) = data.insertion_variants.get_mut(&position) else { continue; };
            let vref_key = if vn.starts_with('+') {
                VarDesc::Ins { seq: vn[1..].as_bytes().iter().map(|b| b.to_ascii_uppercase()).collect() }
            } else {
                continue;
            };
            let Some(vref) = map.get(&vref_key) else { continue; };
            if let Some(cap) = ATGSs_AMP_ATGSs_END.captures(&vn) {
                let tn = cap.get(1).map(|m| m.as_str()).unwrap_or("");
                if !tn.is_empty() {
                    let tref_key = VarDesc::Ins { seq: tn[1..].as_bytes().iter().map(|b| b.to_ascii_uppercase()).collect() };
                    if let Some(tref) = map.get(&tref_key) {
                        if vref.alt_depth < tref.alt_depth {
                            let ref_key = self.get_ref_base(position).map(|ref_base| VarDesc::SNV { ref_base });
                            if let Some(ref_key) = ref_key {
                                let mut ref_var = data
                                    .non_insertion_variants
                                    .get_mut(&position)
                                    .and_then(|m| m.remove(&ref_key));
                                if let Some(ref_var_mut) = ref_var.as_mut() {
                                    let mut tref_mut = tref.clone();
                                    let mut vref_mut = vref.clone();
                                    adj_cnt_with_ref(&mut tref_mut, &vref_mut, Some(ref_var_mut));
                                    map.insert(tref_key.clone(), tref_mut);
                                }
                                if let Some(ref_var) = ref_var {
                                    data.non_insertion_variants.get_mut(&position).unwrap().insert(ref_key, ref_var);
                                }
                            } else {
                                let mut tref_mut = tref.clone();
                                let vref_mut = vref.clone();
                                adj_cnt(&mut tref_mut, &vref_mut);
                                map.insert(tref_key.clone(), tref_mut);
                            }
                            map.remove(&vref_key);
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
                ins_or_del_len,
                mismatch_seq,
                ..
            } => {
                let mut total = *len as i64;
                if let InsOrDelLen::DelLen(extra) = ins_or_del_len {
                    total += *extra as i64;
                }
                let extra = String::from_utf8_lossy(mismatch_seq).to_string();
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
                    if let Some(sc) = data.soft_clips_5end.get_mut(&(position - n as i64 - n2 as i64)) {
                        sc.mark_used();
                    }
                    mn += n2;
                } else {
                    sc5p.push(position - n as i64);
                    if let Some(sc) = data.soft_clips_5end.get_mut(&(position - n as i64)) {
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

    fn substr_bytes(seq: &[u8], begin: i64, len: Option<i64>) -> Vec<u8> {
        let seq_len = seq.len() as i64;
        let mut b = begin;
        if b < 0 {
            b = seq_len + b;
        }
        if b < 0 || b > seq_len {
            return Vec::new();
        }
        let b_usize = b as usize;

        match len {
            None => seq.get(b_usize..).unwrap_or(&[]).to_vec(),
            Some(l) if l > 0 => {
                let end = (b + l).min(seq_len).max(b);
                seq.get(b_usize..end as usize).unwrap_or(&[]).to_vec()
            }
            Some(0) => Vec::new(),
            Some(l) => {
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
            let realigner = VariantRealigner::new(reference.ref_seq.clone(), reference.region_start);
            realigner.process_deletions(&mut sv_input);
        }

        assert!(!reference.ref_seq.is_empty());
        assert!(!sv_input.non_insertion_variants.is_empty());
        assert!(sv_input.insertion_variants.is_empty());
        assert!(sv_input.soft_clips_5end.is_empty());
        assert!(sv_input.soft_clips_3end.is_empty());
    }
}
