//! VarDict Pipeline - Real VarDict Java-equivalent pipeline for Simple Mode
//!
//! This module implements the actual VarDict pipeline flow:
//! ```text
//! BAM Record → CigarParser → VariantRealigner → StructuralVariantsProcessor → ToVarsBuilder → SimplePostProcessor → Output
//! ```
//!
//! This is the proper VarDict Simple Mode pipeline, replacing the simplified
//! `simple_variant_caller.rs` approach.

use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::collections::HashSet;
use std::env;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::Arc;

use anyhow::{Result, Error};
use crackle_kit::tracing::{Level, event};
use rust_htslib::bam::Record;

use crate::data::reference::Reference;
use crate::data::region::Region;
use crate::data::shared_reference::SharedReferenceHandle;
use crate::data::bam_reader::BamReader;
use crate::mods::cigar_parser::CigarParser;
use crate::mods::output_variant::{AmpliconOutputVariant, SimpleOutputVariant, SomaticOutputVariant, Region as OutputRegion};
use crate::mods::structural_variants_processor::{StructuralVariantsProcessor, RealignedVariationData};
use crate::mods::variant_realigner::VariantRealigner;
use crate::scopedata::global_read_only_scope::instance;
use crate::mods::to_vars_builder::{
    ToVarsBuilder, Variant, VariationData, Vars, VarType, determine_genotype, var_type_string,
    check_strand_bias, StrandBiasFlag,
};
use crate::mods::simple_variant_caller::SimpleVarKey;
use crate::scopedata::global_read_only_scope::GlobalReadOnlyScope;
use crate::utils::round_half_even;
use crate::variants::variants::{VarDesc, Variant as RawVariant, SoftClip};
use rand::Rng;

/// Data produced by CigarParser - mirrors Java VariationData
#[derive(Default)]
pub struct CigarParserOutput {
    /// Non-insertion variants by position
    pub non_insertion_vars: HashMap<i64, HashMap<VarDesc, RawVariant>>,
    /// Insertion order of non-insertion variant positions
    pub non_insertion_vars_insert_index: HashMap<i64, usize>,
    /// Insertion variants by position (key is position before insertion)
    pub insertion_vars: HashMap<i64, HashMap<VarDesc, RawVariant>>,
    /// 5' end soft clips by position
    pub soft_clips_5end: HashMap<i64, SoftClip>,
    /// 3' end soft clips by position
    pub soft_clips_3end: HashMap<i64, SoftClip>,
    /// Reference coverage by position
    pub ref_coverage: HashMap<i64, usize>,
    /// MNP map (position -> description -> count)
    pub mnp: HashMap<i64, HashMap<String, usize>>,
    /// Insertion counts by position and description (Java: positionToInsertionCount)
    pub position_to_insertion_count: HashMap<i64, HashMap<String, usize>>,
    /// Deletion counts by position and description (Java: positionToDeletionCount)
    pub position_to_deletions_count: HashMap<i64, HashMap<String, usize>>,
    /// Maximum read length seen
    pub max_read_len: usize,
    /// Discordant read count
    pub discordant_count: usize,
    /// Splice positions ("start-end")
    pub splice: HashSet<String>,
    /// Splice counts by intron key ("start-end")
    pub splice_count: HashMap<String, usize>,
    /// Duplication rate
    pub duprate: f64,
    /// Total reads seen by preprocessor
    pub total_reads: usize,
    /// Duplicate reads filtered by preprocessor
    pub duplicate_reads: usize,
}

impl CigarParserOutput {
    pub fn write_jsonl_snapshot_if_enabled(&self, region: &Region) -> Result<()> {
        let path = match env::var("VARDICT_CIGAR_PARSER_JSONL") {
            Ok(val) => val.trim().to_string(),
            Err(_) => String::new(),
        };
        if path.is_empty() {
            return Ok(());
        }

        self.write_jsonl_snapshot(&path, region)
    }

    fn write_jsonl_snapshot(&self, path: &str, region: &Region) -> Result<()> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        let meta = format!(
            "{{\"region\":\"{}\",\"maxReadLength\":{},\"duprate\":\"{}\"}}",
            json_escape(&region.to_region_string()),
            self.max_read_len,
            fmt_f64(self.duprate),
        );
        write_json_line(&mut writer, "META", 0, "-", &meta)?;

        write_variant_map(&mut writer, "NONINS", &self.non_insertion_vars)?;
        write_variant_map(&mut writer, "INS", &self.insertion_vars)?;
        write_ref_cov(&mut writer, &self.ref_coverage)?;
        write_count_map(&mut writer, "MNP", &self.mnp)?;
        write_count_map(&mut writer, "INSCOUNT", &self.position_to_insertion_count)?;
        write_count_map(&mut writer, "DELCOUNT", &self.position_to_deletions_count)?;
        write_soft_clips(&mut writer, "SCLIP5", &self.soft_clips_5end, false)?;
        write_soft_clips(&mut writer, "SCLIP3", &self.soft_clips_3end, false)?;
        write_splice(&mut writer, &self.splice)?;
        write_splice_count(&mut writer, &self.splice_count)?;

        writer.flush()?;
        Ok(())
    }
}

fn write_realigned_jsonl_snapshot_if_enabled(
    data: &RealignedVariationData,
    region: &Region,
) -> Result<()> {
    let path = match env::var("VARDICT_VARIANT_REALIGNER_JSONL") {
        Ok(val) => val.trim().to_string(),
        Err(_) => String::new(),
    };
    if path.is_empty() {
        return Ok(());
    }

    write_realigned_jsonl_snapshot(data, region, &path)
}

fn write_realigned_jsonl_snapshot(
    data: &RealignedVariationData,
    region: &Region,
    path: &str,
) -> Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    let meta = format!(
        "{{\"region\":\"{}\",\"maxReadLength\":{},\"duprate\":\"{}\"}}",
        json_escape(&region.to_region_string()),
        data.max_read_length,
        fmt_f64(data.duprate),
    );
    write_json_line(&mut writer, "META", 0, "-", &meta)?;

    write_variant_map(&mut writer, "NONINS", &data.non_insertion_variants)?;
    write_variant_map(&mut writer, "INS", &data.insertion_variants)?;
    write_ref_cov(&mut writer, &data.ref_coverage)?;
    write_soft_clips(&mut writer, "SCLIP5", &data.soft_clips_5end, true)?;
    write_soft_clips(&mut writer, "SCLIP3", &data.soft_clips_3end, true)?;

    writer.flush()?;
    Ok(())
}

fn write_structural_variants_jsonl_snapshot_if_enabled(
    data: &RealignedVariationData,
    region: &Region,
) -> Result<()> {
    let path = match env::var("VARDICT_STRUCTURAL_VARIANTS_JSONL") {
        Ok(val) => val.trim().to_string(),
        Err(_) => String::new(),
    };
    if path.is_empty() {
        return Ok(());
    }

    write_structural_variants_jsonl_snapshot(data, region, &path)
}

fn write_structural_variants_jsonl_snapshot(
    data: &RealignedVariationData,
    region: &Region,
    path: &str,
) -> Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    let meta = format!(
        "{{\"region\":\"{}\",\"maxReadLength\":{},\"duprate\":\"{}\"}}",
        json_escape(&region.to_region_string()),
        data.max_read_length,
        fmt_f64(data.duprate),
    );
    write_json_line(&mut writer, "META", 0, "-", &meta)?;

    write_variant_map(&mut writer, "NONINS", &data.non_insertion_variants)?;
    write_variant_map(&mut writer, "INS", &data.insertion_variants)?;
    write_ref_cov(&mut writer, &data.ref_coverage)?;
    write_soft_clips(&mut writer, "SCLIP5", &data.soft_clips_5end, false)?;
    write_soft_clips(&mut writer, "SCLIP3", &data.soft_clips_3end, false)?;

    writer.flush()?;
    Ok(())
}

fn write_tovars_jsonl_snapshot_if_enabled(
    data: &AlignedVarsData,
    region: &Region,
    max_read_len: usize,
    duprate: f64,
) -> Result<()> {
    let path = match env::var("VARDICT_TO_VARS_JSONL") {
        Ok(val) => val.trim().to_string(),
        Err(_) => String::new(),
    };
    if path.is_empty() {
        return Ok(());
    }

    write_tovars_jsonl_snapshot(data, region, max_read_len, duprate, &path)
}

fn write_tovars_jsonl_snapshot(
    data: &AlignedVarsData,
    region: &Region,
    max_read_len: usize,
    duprate: f64,
    path: &str,
) -> Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    let meta = format!(
        "{{\"region\":\"{}\",\"maxReadLength\":{},\"duprate\":\"{}\"}}",
        json_escape(&region.to_region_string()),
        max_read_len,
        fmt_f64(duprate),
    );
    write_json_line(&mut writer, "META", 0, "-", &meta)?;

    write_tovars_variants(&mut writer, &data.aligned_variants)?;

    writer.flush()?;
    Ok(())
}

struct RecordPreprocessorJsonlEntry {
    pos: i64,
    key: String,
    data: String,
}

fn write_record_preprocessor_jsonl_snapshot(
    path: &str,
    region: &Region,
    total_reads: usize,
    duplicate_reads: usize,
    entries: &[RecordPreprocessorJsonlEntry],
) -> Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    let meta = format!(
        "{{\"region\":\"{}\",\"totalReads\":{},\"duplicateReads\":{}}}",
        json_escape(&region.to_region_string()),
        total_reads,
        duplicate_reads,
    );
    write_json_line(&mut writer, "META", 0, "-", &meta)?;

    for entry in entries {
        write_json_line(&mut writer, "RECORD", entry.pos, &entry.key, &entry.data)?;
    }

    writer.flush()?;
    Ok(())
}

fn write_json_line<W: Write>(
    writer: &mut W,
    line_type: &str,
    pos: i64,
    key: &str,
    data: &str,
) -> Result<()> {
    writeln!(
        writer,
        "{{\"type\":\"{}\",\"pos\":{},\"key\":\"{}\",\"data\":{}}}",
        line_type,
        pos,
        json_escape(key),
        data
    )?;
    Ok(())
}

fn write_variant_map<W: Write>(
    writer: &mut W,
    line_type: &str,
    map: &HashMap<i64, HashMap<VarDesc, RawVariant>>,
) -> Result<()> {
    let mut positions: Vec<i64> = map.keys().copied().collect();
    positions.sort_unstable();

    for pos in positions {
        if let Some(vars) = map.get(&pos) {
            let mut keys: Vec<&VarDesc> = vars.keys().collect();
            keys.sort_by(|a, b| a.to_key_string().cmp(&b.to_key_string()));
            for key in keys {
                let key_str = key.to_key_string();
                let var = vars.get(key).expect("variant missing for key");
                let data = format!("{{\"variant\":{}}}", variant_json(var));
                write_json_line(writer, line_type, pos, &key_str, &data)?;
            }
        }
    }

    Ok(())
}

fn write_tovars_variants<W: Write>(
    writer: &mut W,
    map: &HashMap<i64, Vars>,
) -> Result<()> {
    let mut positions: Vec<i64> = map.keys().copied().collect();
    positions.sort_unstable();

    for pos in positions {
        let vars = map.get(&pos).expect("vars missing for position");
        let mut variants: Vec<&Variant> = vars.variants.iter().collect();
        variants.sort_by(|a, b| a.description_string.cmp(&b.description_string));
        for variant in variants {
            let data = format!("{{\"variant\":{}}}", tovars_variant_json(variant));
            write_json_line(writer, "VAR", pos, &variant.description_string, &data)?;
        }
        if let Some(ref_variant) = &vars.reference_variant {
            let data = format!("{{\"variant\":{}}}", tovars_variant_json(ref_variant));
            write_json_line(writer, "REF", pos, &ref_variant.description_string, &data)?;
        }
    }

    Ok(())
}

fn write_ref_cov<W: Write>(writer: &mut W, map: &HashMap<i64, usize>) -> Result<()> {
    let mut positions: Vec<i64> = map.keys().copied().collect();
    positions.sort_unstable();
    for pos in positions {
        let count = map.get(&pos).copied().unwrap_or(0);
        let data = format!("{{\"count\":{}}}", count);
        write_json_line(writer, "REFCOV", pos, "-", &data)?;
    }
    Ok(())
}

fn write_count_map<W: Write>(
    writer: &mut W,
    line_type: &str,
    map: &HashMap<i64, HashMap<String, usize>>,
) -> Result<()> {
    let mut positions: Vec<i64> = map.keys().copied().collect();
    positions.sort_unstable();
    for pos in positions {
        if let Some(inner) = map.get(&pos) {
            let mut keys: Vec<&String> = inner.keys().collect();
            keys.sort();
            for key in keys {
                let count = inner.get(key).copied().unwrap_or(0);
                let data = format!("{{\"count\":{}}}", count);
                write_json_line(writer, line_type, pos, key, &data)?;
            }
        }
    }
    Ok(())
}

fn write_soft_clips<W: Write>(
    writer: &mut W,
    line_type: &str,
    map: &HashMap<i64, SoftClip>,
    compute_consensus_if_unset: bool,
) -> Result<()> {
    let mut positions: Vec<i64> = map.keys().copied().collect();
    positions.sort_unstable();
    for pos in positions {
        if let Some(sc) = map.get(&pos) {
            let data = soft_clip_json(sc, compute_consensus_if_unset);
            write_json_line(writer, line_type, pos, "-", &data)?;
        }
    }
    Ok(())
}

fn write_splice<W: Write>(writer: &mut W, splice: &HashSet<String>) -> Result<()> {
    let mut keys: Vec<&String> = splice.iter().collect();
    keys.sort();
    for key in keys {
        let pos = parse_splice_pos(key);
        let data = String::from("{}");
        write_json_line(writer, "SPLICE", pos, key, &data)?;
    }
    Ok(())
}

fn write_splice_count<W: Write>(
    writer: &mut W,
    splice_count: &HashMap<String, usize>,
) -> Result<()> {
    let mut keys: Vec<&String> = splice_count.keys().collect();
    keys.sort();
    for key in keys {
        let count = splice_count.get(key).copied().unwrap_or(0);
        let pos = parse_splice_pos(key);
        let data = format!("{{\"count\":{}}}", count);
        write_json_line(writer, "SPLICECOUNT", pos, key, &data)?;
    }
    Ok(())
}

fn parse_splice_pos(key: &str) -> i64 {
    key.split_once('-')
        .and_then(|(start, _)| start.parse::<i64>().ok())
        .unwrap_or(0)
}

fn variant_json(v: &RawVariant) -> String {
    format!(
        "{{\"varsCount\":{},\"varsCountOnForward\":{},\"varsCountOnReverse\":{},\"extracnt\":{},\"meanPosition\":\"{}\",\"meanQuality\":\"{}\",\"meanMappingQuality\":\"{}\",\"numberOfMismatches\":\"{}\",\"lowQualityReadsCount\":{},\"highQualityReadsCount\":{},\"pstd\":{},\"qstd\":{},\"pp\":{},\"pq\":\"{}\"}}",
        v.alt_depth,
        v.alt_depth_fwd,
        v.alt_depth_rev,
        v.extra_cnt,
        fmt_f64(v.mean_pos),
        fmt_f64(v.mean_qual),
        fmt_f64(v.mean_mapq),
        fmt_f64(v.nm),
        v.low_qual_read_cnt,
        v.high_qual_read_cnt,
        v.pstd,
        v.qstd,
        v.pp,
        fmt_f64(v.pq),
    )
}

fn tovars_variant_json(v: &Variant) -> String {
    let msint = v.msint.round() as i64;
    let var_type = var_type_string(&v.refallele, &v.varallele);
    format!(
        "{{\"descriptionString\":\"{}\",\"positionCoverage\":{},\"varsCountOnForward\":{},\"varsCountOnReverse\":{},\"strandBiasFlag\":\"{}\",\"frequency\":\"{}\",\"meanPosition\":\"{}\",\"pstd\":{},\"meanQuality\":\"{}\",\"qstd\":{},\"meanMappingQuality\":\"{}\",\"highQualityReadsFrequency\":\"{}\",\"extraFrequency\":\"{}\",\"shift3\":{},\"msi\":\"{}\",\"msint\":{},\"numberOfMismatches\":\"{}\",\"hicnt\":{},\"hicov\":{},\"leftseq\":\"{}\",\"rightseq\":\"{}\",\"startPosition\":{},\"endPosition\":{},\"refReverseCoverage\":{},\"refForwardCoverage\":{},\"totalPosCoverage\":{},\"duprate\":\"{}\",\"genotype\":\"{}\",\"varallele\":\"{}\",\"refallele\":\"{}\",\"varType\":\"{}\",\"crispr\":{}}}",
        json_escape(&v.description_string),
        v.position_coverage,
        v.vars_count_on_forward,
        v.vars_count_on_reverse,
        json_escape(&v.strand_bias_flag.to_string()),
        fmt_f64_with("0.0000", v.frequency),
        fmt_f64_with("0.0", v.mean_position),
        v.is_at_least_at_2_positions,
        fmt_f64_with("0.0", v.mean_quality),
        v.has_at_least_2_diff_qualities,
        fmt_f64_with("0.0", v.mean_mapping_quality),
        fmt_f64_with("0.0000", v.high_quality_reads_frequency),
        fmt_f64_with("0.0000", v.extra_frequency),
        v.shift3,
        fmt_f64_with("0.000", v.msi),
        msint,
        fmt_f64_with("0.0", v.nm),
        v.high_qual_read_cnt,
        v.hicov,
        json_escape(&v.leftseq),
        json_escape(&v.rightseq),
        v.start_position,
        v.end_position,
        v.ref_reverse_count,
        v.ref_forward_count,
        v.total_pos_coverage,
        fmt_f64_with("0.000", v.duprate),
        json_escape(&v.genotype),
        json_escape(&v.varallele),
        json_escape(&v.refallele),
        json_escape(&var_type),
        v.crispr,
    )
}

fn soft_clip_json(sc: &SoftClip, compute_consensus_if_unset: bool) -> String {
    let consensus_bytes = if sc.consensus_seq_is_set() {
        sc.consensus_seq().to_vec()
    } else if compute_consensus_if_unset {
        let mut soft_clip_for_snapshot = SoftClip::default();
        soft_clip_for_snapshot.var = sc.var.clone();
        soft_clip_for_snapshot.nt = sc.nt.clone();
        soft_clip_for_snapshot.seq = sc.seq.clone();
        crate::variants::var_utils::find_conseq(&mut soft_clip_for_snapshot, 0)
    } else {
        Vec::new()
    };
    let consensus = String::from_utf8_lossy(&consensus_bytes);
    let nt = soft_clip_nt_json(&sc.nt);
    let seq = soft_clip_seq_json(&sc.seq);
    format!(
        "{{\"variant\":{},\"nt\":{},\"seq\":{},\"sequence\":\"{}\",\"used\":{}}}",
        variant_json(&sc.var),
        nt,
        seq,
        json_escape(&consensus),
        sc.used()
    )
}

fn soft_clip_nt_json(map: &std::collections::BTreeMap<i64, crackle_kit::nuc_base_map::NucBaseMap<usize>>) -> String {
    let mut out = String::from("[");
    let mut first = true;
    for (offset, base_map) in map.iter() {
        for base in [b'A', b'C', b'G', b'T', b'N'] {
            if let Some(val) = base_map.get(base) {
                if !first {
                    out.push(',');
                }
                first = false;
                out.push_str(&format!(
                    "{{\"offset\":{},\"base\":\"{}\",\"count\":{}}}",
                    offset,
                    base as char,
                    val
                ));
            }
        }
    }
    out.push(']');
    out
}

fn soft_clip_seq_json(map: &std::collections::BTreeMap<usize, crackle_kit::nuc_base_map::NucBaseMap<RawVariant>>) -> String {
    let mut out = String::from("[");
    let mut first = true;
    for (offset, base_map) in map.iter() {
        for base in [b'A', b'C', b'G', b'T', b'N'] {
            if let Some(val) = base_map.get(base) {
                if !first {
                    out.push(',');
                }
                first = false;
                out.push_str(&format!(
                    "{{\"offset\":{},\"base\":\"{}\",\"variant\":{}}}",
                    offset,
                    base as char,
                    variant_json(val)
                ));
            }
        }
    }
    out.push(']');
    out
}

fn fmt_f64(value: f64) -> String {
    fmt_f64_with("0.000", value)
}

fn fmt_f64_with(pattern: &str, value: f64) -> String {
    let rounded = round_half_even(pattern, value);
    let decimals = pattern
        .split('.')
        .nth(1)
        .map(|s| s.len())
        .unwrap_or(0);
    if decimals == 0 {
        return format!("{:.0}", rounded);
    }
    format!("{:.*}", decimals, rounded)
}

fn json_escape(value: &str) -> String {
    let mut out = String::with_capacity(value.len());
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            _ => out.push(ch),
        }
    }
    out
}

/// Data after realignment - mirrors Java RealignedVariationData
#[derive(Debug, Clone, Default)]
pub struct RealignedOutput {
    /// Non-insertion variants (may be modified by realigner)
    pub non_insertion_vars: HashMap<i64, HashMap<VarDesc, RawVariant>>,
    /// Insertion order of non-insertion variant positions
    pub non_insertion_vars_insert_index: HashMap<i64, usize>,
    /// Insertion variants
    pub insertion_vars: HashMap<i64, HashMap<VarDesc, RawVariant>>,
    /// Reference coverage by position
    pub ref_coverage: HashMap<i64, usize>,
    /// Duplication rate
    pub duprate: f64,
    /// Maximum read length
    pub max_read_len: usize,
    /// Splice positions ("start-end")
    pub splice: HashSet<String>,
}

/// Final aligned variants data - mirrors Java AlignedVarsData
#[derive(Debug, Clone, Default)]
pub struct AlignedVarsData {
    /// Variants by position
    pub aligned_variants: HashMap<i64, Vars>,
    /// Insertion order of positions into aligned_variants
    pub aligned_variants_order: Vec<i64>,
    /// Reference coverage by position
    pub ref_coverage: HashMap<i64, usize>,
}

#[derive(Debug, Clone)]
pub struct RegionAlignedVarsOutput {
    pub aligned_vars: AlignedVarsData,
    pub splice: HashSet<String>,
    pub max_read_length: usize,
}

#[derive(Debug, Clone, Default)]
pub struct SomaticCombineLookupResult {
    pub combined_variant: Option<Variant>,
    pub max_read_length: usize,
}

type SomaticCombineLookup =
    dyn Fn(&str, i64, &str, usize) -> SomaticCombineLookupResult;

fn java_hashmap_capacity(size: usize) -> usize {
    let mut capacity = 16usize;
    if size == 0 {
        return capacity;
    }

    let mut threshold = capacity - (capacity >> 2);
    while size > threshold {
        capacity <<= 1;
        threshold = capacity - (capacity >> 2);
    }

    capacity
}

fn java_hashmap_bucket_index(key: i64, capacity: usize) -> usize {
    let key = key as i32 as u32;
    let hash = key ^ (key >> 16);
    (hash as usize) & (capacity - 1)
}

fn java_hashmap_iteration_order<I>(
    keys: I,
    size: usize,
    insertion_index: Option<&HashMap<i64, usize>>,
) -> Vec<i64>
where
    I: Iterator<Item = i64>,
{
    let capacity = java_hashmap_capacity(size);
    let mut entries: Vec<(i64, usize, usize)> = Vec::new();

    for key in keys {
        let bucket = java_hashmap_bucket_index(key, capacity);
        let order = insertion_index
            .and_then(|map| map.get(&key).copied())
            .unwrap_or(usize::MAX);
        entries.push((key, bucket, order));
    }

    entries.sort_by(|a, b| {
        a.1.cmp(&b.1)
            .then_with(|| a.2.cmp(&b.2))
            .then_with(|| a.0.cmp(&b.0))
    });

    entries.into_iter().map(|(key, _, _)| key).collect()
}

/// Main VarDict pipeline configuration
pub struct VarDictPipeline {
    pub sample_name: String,
    pub min_frequency: f64,
    pub min_base_quality: f64,
    pub min_mapping_quality: u8,
    pub do_pileup: bool,
}

struct RecordPreprocessorState {
    total_reads: usize,
    duplicate_reads: usize,
    duplicates: HashSet<String>,
    first_matching_position: i64,
}

impl RecordPreprocessorState {
    fn new() -> Self {
        Self {
            total_reads: 0,
            duplicate_reads: 0,
            duplicates: HashSet::new(),
            first_matching_position: -1,
        }
    }
}

impl VarDictPipeline {
    /// Create a new pipeline with default settings
    pub fn new(sample_name: &str) -> Self {
        Self {
            sample_name: sample_name.to_string(),
            min_frequency: 0.01,
            min_base_quality: 0.0,
            min_mapping_quality: 0,
            do_pileup: false,
        }
    }

    /// Set minimum variant frequency
    pub fn with_min_frequency(mut self, min_frequency: f64) -> Self {
        self.min_frequency = min_frequency;
        self
    }

    /// Set minimum base quality
    pub fn with_min_base_quality(mut self, min_base_quality: f64) -> Self {
        self.min_base_quality = min_base_quality;
        self
    }

    /// Set minimum mapping quality
    pub fn with_min_mapping_quality(mut self, min_mapping_quality: u8) -> Self {
        self.min_mapping_quality = min_mapping_quality;
        self
    }

    /// Enable pileup mode
    pub fn with_pileup(mut self, enable: bool) -> Self {
        self.do_pileup = enable;
        self
    }

    /// Process a region using BamReader and SharedReference
    ///
    /// This is the main entry point for parallel pipeline usage.
    /// It fetches reads from the BAM file and processes them through the pipeline.
    pub fn process_region_from_bam(
        &self,
        region: &Region,
        shared_reference: &SharedReferenceHandle,
        bam_reader: &mut BamReader,
        instance: Arc<GlobalReadOnlyScope>,
    ) -> Result<Vec<String>> {
        let region_output = self.process_region_to_aligned_vars_from_bam(
            region,
            shared_reference,
            bam_reader,
            Arc::clone(&instance),
        )?;

        let start_post = std::time::Instant::now();
        let output_lines = self.run_simple_post_processor(
            region_output.aligned_vars,
            region,
            &region_output.splice,
        )?;
        let elapsed_post = start_post.elapsed();

        event!(Level::INFO, "[TIMING] PostProcessor: {:.3}s - {} lines",
            elapsed_post.as_secs_f64(),
            output_lines.len());

        Ok(output_lines)
    }

    pub fn process_region_to_aligned_vars_from_bam(
        &self,
        region: &Region,
        shared_reference: &SharedReferenceHandle,
        bam_reader: &mut BamReader,
        instance: Arc<GlobalReadOnlyScope>,
    ) -> Result<RegionAlignedVarsOutput> {
        let bam_paths = instance.bam_paths.clone();
        self.process_region_to_aligned_vars_from_bam_with_paths(
            region,
            shared_reference,
            bam_reader,
            instance,
            &bam_paths,
        )
    }

    pub fn process_region_to_aligned_vars_from_bam_with_paths(
        &self,
        region: &Region,
        shared_reference: &SharedReferenceHandle,
        bam_reader: &mut BamReader,
        instance: Arc<GlobalReadOnlyScope>,
        bam_paths: &[String],
    ) -> Result<RegionAlignedVarsOutput> {
        let extend = (instance.conf.number_nucleotide_to_extend + instance.conf.reference_extension)
            .max(0) as usize;
        let extended_start = if region.start() > extend {
            region.start() - extend
        } else {
            1
        };
        let mut extended_end = region.end() + extend;
        if let Some(&chr_len) = instance.chr_lens.get(region.chr()) {
            if extended_end > chr_len {
                extended_end = chr_len;
            }
        }
        
        // Get reference sequence for this region from shared reference
        let ref_seq = match shared_reference.get_subseq(region.chr(), extended_start, extended_end) {
            Some(seq) => seq.to_vec(),
            None => {
                return Err(anyhow::anyhow!(
                    "Failed to get reference for {}:{}-{}",
                    region.chr(), extended_start, extended_end
                ));
            }
        };
        let mut reference = Reference::new_with_start(ref_seq, extended_start as i64);
        let chr_len = instance.chr_lens.get(region.chr()).copied();
        reference.build_seed_map(extended_end as i64, chr_len);

        // Get SAM filter from instance configuration
        let sam_filter = instance.conf.sam_filter;

        let cigar_output = self.run_cigar_parser_from_bam(
            region,
            &reference,
            Arc::clone(&instance),
            bam_reader,
            sam_filter,
        )?;

        cigar_output.write_jsonl_snapshot_if_enabled(region)?;

        self.process_region_to_aligned_vars_from_cigar_output(
            cigar_output,
            region,
            &reference,
            bam_paths,
        )
    }

    pub fn process_region_to_aligned_vars_from_bam_paths(
        &self,
        region: &Region,
        shared_reference: &SharedReferenceHandle,
        bam_paths: &[String],
        instance: Arc<GlobalReadOnlyScope>,
    ) -> Result<RegionAlignedVarsOutput> {
        let extend = (instance.conf.number_nucleotide_to_extend + instance.conf.reference_extension)
            .max(0) as usize;
        let extended_start = if region.start() > extend {
            region.start() - extend
        } else {
            1
        };
        let mut extended_end = region.end() + extend;
        if let Some(&chr_len) = instance.chr_lens.get(region.chr()) {
            if extended_end > chr_len {
                extended_end = chr_len;
            }
        }

        let ref_seq = match shared_reference.get_subseq(region.chr(), extended_start, extended_end) {
            Some(seq) => seq.to_vec(),
            None => {
                return Err(anyhow::anyhow!(
                    "Failed to get reference for {}:{}-{}",
                    region.chr(), extended_start, extended_end
                ));
            }
        };
        let mut reference = Reference::new_with_start(ref_seq, extended_start as i64);
        let chr_len = instance.chr_lens.get(region.chr()).copied();
        reference.build_seed_map(extended_end as i64, chr_len);

        let sam_filter = instance.conf.sam_filter;
        let cigar_output = self.run_cigar_parser_from_bam_paths(
            region,
            &reference,
            Arc::clone(&instance),
            bam_paths,
            sam_filter,
        )?;

        self.process_region_to_aligned_vars_from_cigar_output(
            cigar_output,
            region,
            &reference,
            bam_paths,
        )
    }

    #[deprecated]
    pub(crate) fn collect_filtered_records(
        &self,
        region: &Region,
        bam_reader: &mut BamReader,
        sam_filter: u32,
    ) -> Result<(Vec<Record>, Vec<String>)> {
        let mut records = Vec::new();
        let mut lines = Vec::new();
        let mut record = Record::new();
        let mut preprocess_state = RecordPreprocessorState::new();
        let header_view = Arc::new(rust_htslib::bam::HeaderView::from_header(bam_reader.header()));

        bam_reader.fetch(region.chr(), region.start(), region.end())?;

        while bam_reader.read(&mut record).unwrap_or(false) {
            let mate_ref_name = Self::mate_reference_name(&record, bam_reader);
            if !self.passes_preprocess(&record, sam_filter, &mut preprocess_state, &mate_ref_name) {
                continue;
            }

            let qname = std::str::from_utf8(record.qname()).unwrap_or("");
            let alignment_start = record.pos() + 1;
            let mate_start = if record.mpos() >= 0 { record.mpos() + 1 } else { 0 };
            let cigar = record.cigar().to_string();
            let seq_bytes = record.seq().as_bytes();
            let seq_str = String::from_utf8_lossy(&seq_bytes).to_string();
            let qual_str: String = record
                .qual()
                .iter()
                .map(|q| (*q as u8 + 33) as char)
                .collect();
            let mapq = record.mapq();

            lines.push(format!(
                "FILTERED_READ\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                qname,
                record.flags(),
                alignment_start,
                mate_start,
                cigar,
                seq_str,
                qual_str,
                mapq
            ));

            let mut cloned = record.clone();
            cloned.set_header(header_view.clone());
            records.push(cloned);
        }

        Ok((records, lines))
    }

    fn passes_preprocess(
        &self,
        record: &Record,
        sam_filter: u32,
        state: &mut RecordPreprocessorState,
        mate_ref_name: &str,
    ) -> bool {
        // Java preprocessRecord: downsampling
        if let Some(downsample) = instance().conf.downsampling {
            if rand::random::<f64>() <= downsample {
                return false;
            }
        }

        // 1. Java SamView.read(): Skip records that match the filter flags
        if sam_filter != 0 {
            if (record.flags() & (sam_filter as u16)) != 0 {
                return false;
            }
        }

        // 2. Java preprocessRecord line 117: Ignore low mapping quality reads
        let min_mapq = instance().conf.mapping_quality.unwrap_or(self.min_mapping_quality);
        if min_mapq > 0 && record.mapq() < min_mapq {
            return false;
        }

        // 3. Java preprocessRecord line 122: Skip not primary alignment reads
        const SECONDARY_ALIGNMENT: u16 = 0x100;
        if (record.flags() & SECONDARY_ALIGNMENT) != 0 && sam_filter != 0 {
            return false;
        }

        // 4. Java preprocessRecord line 124: Skip reads where sequence is not stored in read
        let seq = record.seq();
        if seq.len() == 0 || (seq.len() == 1 && seq.as_bytes()[0] == b'*') {
            return false;
        }

        state.total_reads += 1;

        // Java preprocessRecord: duplicate removal (-t)
        if instance().conf.remove_duplicated_reads {
            let alignment_start = record.pos() + 1;
            if alignment_start != state.first_matching_position {
                state.duplicates.clear();
            }

            let mate_start = if record.mpos() >= 0 { record.mpos() + 1 } else { 0 };

            if mate_start < 10 {
                let dup_key = format!("{}-{}-{}", alignment_start, mate_ref_name, mate_start);
                if state.duplicates.contains(&dup_key) {
                    state.duplicate_reads += 1;
                    return false;
                }
                state.duplicates.insert(dup_key);
                state.first_matching_position = alignment_start;
            } else if record.is_paired() && record.is_mate_unmapped() {
                let dup_key = format!("{}-{}", alignment_start, record.cigar().to_string());
                if state.duplicates.contains(&dup_key) {
                    state.duplicate_reads += 1;
                    return false;
                }
                state.duplicates.insert(dup_key);
                state.first_matching_position = alignment_start;
            }
        }

        true
    }

    fn mate_reference_name(record: &Record, bam_reader: &BamReader) -> String {
        if record.mtid() < 0 {
            return "*".to_string();
        }

        if record.mtid() == record.tid() {
            return "=".to_string();
        }

        let mtid = record.mtid() as usize;
        bam_reader
            .target_names()
            .get(mtid)
            .cloned()
            .unwrap_or_else(|| "*".to_string())
    }

    /// Process a batch of records for a region
    ///
    /// This is the main entry point - processes BAM records and returns output lines.
    pub fn process_region<I>(
        &self,
        records: I,
        region: &Region,
        reference: &Reference,
        instance: Arc<GlobalReadOnlyScope>,
    ) -> Result<Vec<String>>
    where
        I: Iterator<Item = Record>,
    {
        let bam_paths = instance.bam_paths.clone();

        let mut working_reference = reference.clone();
        let region_end_for_seed = working_reference.region_start + working_reference.ref_seq.len() as i64 - 1;
        let chr_len = instance.chr_lens.get(region.chr()).copied();
        working_reference.build_seed_map(region_end_for_seed, chr_len);

        // Step 1: Parse CIGAR strings (CigarParser)
        let start_cigar = std::time::Instant::now();
        let cigar_output = self.run_cigar_parser(records, region, &working_reference, instance)?;
        let elapsed_cigar = start_cigar.elapsed();
        // Step 2: Write JSONL snapshot if enabled
        cigar_output.write_jsonl_snapshot_if_enabled(region)?;


        event!(Level::INFO, "[TIMING] CigarParser: {:.3}s - {} non_insertion_vars, {} ref_coverage positions",
            elapsed_cigar.as_secs_f64(),
            cigar_output.non_insertion_vars.len(),
            cigar_output.ref_coverage.len());

        self.process_region_from_cigar_output(
            cigar_output,
            region,
            &working_reference,
            &bam_paths,
        )
    }

    fn process_region_from_cigar_output(
        &self,
        cigar_output: CigarParserOutput,
        region: &Region,
        reference: &Reference,
        bam_paths: &[String],
    ) -> Result<Vec<String>> {
        let region_output = self.process_region_to_aligned_vars_from_cigar_output(
            cigar_output,
            region,
            reference,
            bam_paths,
        )?;

        let start_post = std::time::Instant::now();
        let output_lines = self.run_simple_post_processor(
            region_output.aligned_vars,
            region,
            &region_output.splice,
        )?;
        let elapsed_post = start_post.elapsed();

        event!(Level::INFO, "[TIMING] PostProcessor: {:.3}s - {} lines",
            elapsed_post.as_secs_f64(),
            output_lines.len());

        Ok(output_lines)
    }

    fn process_region_to_aligned_vars_from_cigar_output(
        &self,
        cigar_output: CigarParserOutput,
        region: &Region,
        reference: &Reference,
        bam_paths: &[String],
    ) -> Result<RegionAlignedVarsOutput> {
        let start_realign = std::time::Instant::now();
        let realigned_output =
            self.run_variant_realigner_and_sv_processor(cigar_output, region, reference, bam_paths)?;
        let elapsed_realign = start_realign.elapsed();

        event!(Level::INFO, "[TIMING] VariantRealigner+SVProcessor: {:.3}s - {} non_insertion_vars, {} ref_coverage",
            elapsed_realign.as_secs_f64(),
            realigned_output.non_insertion_vars.len(),
            realigned_output.ref_coverage.len());

        let start_tovars = std::time::Instant::now();
        let splice = realigned_output.splice.clone();
        let max_read_length = realigned_output.max_read_len;
        let aligned_vars = self.run_to_vars_builder(realigned_output, reference, region)?;
        let elapsed_tovars = start_tovars.elapsed();

        event!(Level::INFO, "[TIMING] ToVarsBuilder: {:.3}s - {} variants",
            elapsed_tovars.as_secs_f64(),
            aligned_vars.aligned_variants.len());

        Ok(RegionAlignedVarsOutput {
            aligned_vars,
            splice,
            max_read_length,
        })
    }

    pub fn run_amplicon_post_processor(
        &self,
        group_region: &Region,
        vars_per_amplicon: &[HashMap<i64, Vars>],
        amplicon_regions: &[Region],
        splice: &HashSet<String>,
    ) -> Vec<String> {
        let mut output_lines = Vec::new();
        let output_region = OutputRegion {
            chr: group_region.chr().to_string(),
            start: group_region.start() as i64,
            end: group_region.end() as i64,
            gene: group_region.gene().to_string(),
        };

        let mut amplicons_on_positions: std::collections::BTreeMap<i64, Vec<(usize, &Region)>> =
            std::collections::BTreeMap::new();
        for (amplicon_number, amp_region) in amplicon_regions.iter().enumerate() {
            for position in amp_region.insert_start()..=amp_region.insert_end() {
                amplicons_on_positions
                    .entry(position as i64)
                    .or_default()
                    .push((amplicon_number, amp_region));
            }
        }

        for (position, amplicon_regions_at_pos) in amplicons_on_positions {
            let mut gvs: Vec<(Variant, String)> = Vec::new();
            let mut ref_variants: Vec<Variant> = Vec::new();
            let mut vref_list: Vec<Variant> = Vec::new();
            let mut goodmap: HashSet<String> = HashSet::new();
            let mut vcovs: Vec<usize> = Vec::new();
            let mut good_variants_on_amp: std::collections::BTreeMap<usize, Vec<Variant>> =
                std::collections::BTreeMap::new();
            let mut maxcov = 0usize;

            for (amplicon_number, amp_region) in amplicon_regions_at_pos.iter().copied() {
                let vars_at_amplicon = vars_per_amplicon
                    .get(amplicon_number)
                    .and_then(|vars| vars.get(&position));
                let variants_on_amplicon = vars_at_amplicon.map(|vars| &vars.variants);
                let ref_amplicon = vars_at_amplicon.and_then(|vars| vars.reference_variant.as_ref());

                if let Some(variants) = variants_on_amplicon {
                    if !variants.is_empty() {
                        let mut good_vars = Vec::new();
                        for variant in variants {
                            vcovs.push(variant.total_pos_coverage);
                            if variant.total_pos_coverage > maxcov {
                                maxcov = variant.total_pos_coverage;
                            }
                            if self.is_good_var(variant, ref_amplicon, splice) {
                                gvs.push((
                                    variant.clone(),
                                    format!(
                                        "{}:{}-{}",
                                        amp_region.chr(),
                                        amp_region.start(),
                                        amp_region.end()
                                    ),
                                ));
                                good_vars.push(variant.clone());
                                goodmap.insert(format!(
                                    "{}-{}-{}",
                                    amplicon_number,
                                    variant.refallele,
                                    variant.varallele
                                ));
                            }
                        }
                        if !good_vars.is_empty() {
                            good_variants_on_amp.insert(amplicon_number, good_vars);
                        }
                    } else if let Some(reference_variant) = ref_amplicon {
                        vcovs.push(reference_variant.total_pos_coverage);
                    } else {
                        vcovs.push(0);
                    }
                } else {
                    vcovs.push(0);
                }

                if let Some(reference_variant) = ref_amplicon {
                    ref_variants.push(reference_variant.clone());
                }
            }

            let nocov = vcovs
                .iter()
                .filter(|coverage| (**coverage as f64) < (maxcov as f64 / 50.0))
                .count();

            gvs.sort_by(|left, right| {
                right
                    .0
                    .frequency
                    .partial_cmp(&left.0.frequency)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            ref_variants.sort_by(|left, right| right.total_pos_coverage.cmp(&left.total_pos_coverage));

            if gvs.is_empty() {
                if self.do_pileup {
                    if let Some(reference_variant) = ref_variants.first() {
                        vref_list.push(reference_variant.clone());
                    } else {
                        output_lines.push(
                            AmpliconOutputVariant::from_variant(
                                None,
                                &output_region,
                                &[],
                                &[],
                                None,
                                position,
                                0,
                                nocov,
                                false,
                                &self.sample_name,
                            )
                            .to_string(),
                        );
                        continue;
                    }
                } else {
                    continue;
                }
            } else {
                self.fill_vref_list(&gvs, &mut vref_list);
            }

            let mut flag = self.is_amp_bias_flag(&good_variants_on_amp);
            let mut good_variants = gvs.clone();

            for mut vref in vref_list {
                if flag {
                    let top_description = &gvs[0].0.description_string;
                    let mut gcnt: Vec<(Variant, String)> = Vec::new();
                    for (amplicon_number, amp_region) in amplicon_regions_at_pos.iter().copied() {
                        if let Some(vars_at_amplicon) = vars_per_amplicon
                            .get(amplicon_number)
                            .and_then(|vars| vars.get(&position))
                        {
                            if let Some(variant) = vars_at_amplicon
                                .variants
                                .iter()
                                .find(|variant| variant.description_string == *top_description)
                            {
                                if self.is_good_var(
                                    variant,
                                    vars_at_amplicon.reference_variant.as_ref(),
                                    splice,
                                ) {
                                    gcnt.push((
                                        variant.clone(),
                                        format!(
                                            "{}:{}-{}",
                                            amp_region.chr(),
                                            amp_region.start(),
                                            amp_region.end()
                                        ),
                                    ));
                                }
                            }
                        }
                    }
                    if gcnt.len() == gvs.len() {
                        flag = false;
                    }
                    gcnt.sort_by(|left, right| {
                        right
                            .0
                            .frequency
                            .partial_cmp(&left.0.frequency)
                            .unwrap_or(std::cmp::Ordering::Equal)
                    });
                    good_variants = gcnt;
                }

                let initial_gvscnt = self.count_variant_on_amplicons(&vref, &good_variants_on_amp);
                let mut current_gvscnt = initial_gvscnt;
                let mut bad_variants: Vec<(Option<Variant>, String)> = Vec::new();

                if initial_gvscnt != amplicon_regions_at_pos.len() || flag {
                    for (amplicon_number, amp_region) in amplicon_regions_at_pos.iter().copied() {
                        if goodmap.contains(&format!(
                            "{}-{}-{}",
                            amplicon_number,
                            vref.refallele,
                            vref.varallele
                        )) {
                            continue;
                        }
                        if self.do_pileup && vref.refallele == vref.varallele {
                            continue;
                        }
                        if vref.start_position >= amp_region.insert_start() as i64
                            && vref.end_position <= amp_region.insert_end() as i64
                        {
                            let region_string = format!(
                                "{}:{}-{}",
                                amp_region.chr(),
                                amp_region.start(),
                                amp_region.end()
                            );
                            let bad_variant = vars_per_amplicon
                                .get(amplicon_number)
                                .and_then(|vars| vars.get(&position))
                                .and_then(|vars_at_amplicon| {
                                    vars_at_amplicon
                                        .variants
                                        .first()
                                        .cloned()
                                        .or_else(|| vars_at_amplicon.reference_variant.clone())
                                });
                            bad_variants.push((bad_variant, region_string));
                        } else if (vref.start_position < amp_region.insert_end() as i64
                            && (amp_region.insert_end() as i64) < vref.end_position)
                            || (vref.start_position < amp_region.insert_start() as i64
                                && (amp_region.insert_start() as i64) < vref.end_position)
                        {
                            if current_gvscnt > 1 {
                                current_gvscnt -= 1;
                            }
                        }
                    }
                }

                if flag && current_gvscnt < initial_gvscnt {
                    flag = false;
                }

                if var_type_string(&vref.refallele, &vref.varallele) == "Complex" {
                    vref.adj_complex();
                }

                let debug_prefix = if instance().conf.debug {
                    amplicon_regions_at_pos.iter().copied().find_map(|(amplicon_number, _)| {
                        let vars_at_amplicon = vars_per_amplicon
                            .get(amplicon_number)
                            .and_then(|vars| vars.get(&position))?;

                        let matches_current = vars_at_amplicon
                            .variants
                            .iter()
                            .any(|variant| {
                                variant.refallele == vref.refallele
                                    && variant.varallele == vref.varallele
                            })
                            || vars_at_amplicon
                                .reference_variant
                                .as_ref()
                                .map_or(false, |variant| {
                                    variant.refallele == vref.refallele
                                        && variant.varallele == vref.varallele
                                });

                        if matches_current {
                            Some(build_amplicon_debug_prefix(vars_at_amplicon))
                        } else {
                            None
                        }
                    })
                } else {
                    None
                };

                output_lines.push(
                    AmpliconOutputVariant::from_variant(
                        Some(&vref),
                        &output_region,
                        &good_variants,
                        &bad_variants,
                        debug_prefix.as_deref(),
                        position,
                        current_gvscnt,
                        nocov,
                        flag,
                        &self.sample_name,
                    )
                    .to_string(),
                );
            }
        }

        output_lines
    }

    fn count_variant_on_amplicons(
        &self,
        variant: &Variant,
        good_variants_on_amp: &std::collections::BTreeMap<usize, Vec<Variant>>,
    ) -> usize {
        good_variants_on_amp
            .values()
            .flat_map(|variants| variants.iter())
            .filter(|candidate| {
                candidate.refallele == variant.refallele && candidate.varallele == variant.varallele
            })
            .count()
    }

    fn fill_vref_list(&self, gvs: &[(Variant, String)], vref_list: &mut Vec<Variant>) {
        for (good_variant, _) in gvs {
            let already_added = vref_list.iter().any(|existing| {
                existing.varallele == good_variant.varallele
                    && existing.refallele == good_variant.refallele
            });
            if !already_added {
                vref_list.push(good_variant.clone());
            }
        }
    }

    fn is_amp_bias_flag(
        &self,
        good_variants_on_amp: &std::collections::BTreeMap<usize, Vec<Variant>>,
    ) -> bool {
        if good_variants_on_amp.is_empty() {
            return false;
        }

        let amplicon_list: Vec<usize> = good_variants_on_amp.keys().copied().collect();
        for index in 0..amplicon_list.len().saturating_sub(1) {
            let current_amplicon = amplicon_list[index];
            let next_amplicon = amplicon_list[index + 1];
            let Some(current_variants) = good_variants_on_amp.get(&current_amplicon) else {
                return true;
            };
            let Some(next_variants) = good_variants_on_amp.get(&next_amplicon) else {
                return true;
            };

            if current_variants.len() != next_variants.len() {
                return true;
            }

            let mut current_sorted = current_variants.clone();
            let mut next_sorted = next_variants.clone();
            current_sorted.sort_by(|left, right| right.total_pos_coverage.cmp(&left.total_pos_coverage));
            next_sorted.sort_by(|left, right| right.total_pos_coverage.cmp(&left.total_pos_coverage));

            for position in 0..current_sorted.len() {
                if current_sorted[position].description_string != next_sorted[position].description_string {
                    return true;
                }
            }
        }

        false
    }

    /// Step 1: Run CigarParser on records
    fn run_cigar_parser<I>(
        &self,
        records: I,
        region: &Region,
        reference: &Reference,
        instance: Arc<GlobalReadOnlyScope>,
    ) -> Result<CigarParserOutput>
    where
        I: Iterator<Item = Record>,
    {
        // Create CigarParser with region and reference
        let mut cigar_parser = CigarParser::new(
            region.clone(),
            reference.clone(),
            instance,
        );

        // Process each record
        // Note: We need mutable records for parse_cigar
        let mut records_vec: Vec<Record> = records.collect();
        for record in records_vec.iter_mut() {
            // CigarParser.parse_cigar expects &mut Record
            // The process_records method handles the iteration internally
        }
        
        // Process all records
        cigar_parser.process_records(records_vec.iter_mut())?;

        Ok(self.build_cigar_output(&mut cigar_parser, 0, 0))
    }

    fn run_cigar_parser_from_bam(
        &self,
        region: &Region,
        reference: &Reference,
        instance: Arc<GlobalReadOnlyScope>,
        bam_reader: &mut BamReader,
        sam_filter: u32,
    ) -> Result<CigarParserOutput> {
        let mut cigar_parser = CigarParser::new(
            region.clone(),
            reference.clone(),
            instance,
        );

        let mut preprocess_state = RecordPreprocessorState::new();
        let jsonl_path = env::var("VARDICT_RECORD_PREPROCESSOR_JSONL")
            .ok()
            .map(|val| val.trim().to_string())
            .filter(|val| !val.is_empty());
        let mut jsonl_entries = jsonl_path.as_ref().map(|_| Vec::new());

        bam_reader.fetch(region.chr(), region.start(), region.end())?;

        let mut record = Record::new();
        while bam_reader.read(&mut record).unwrap_or(false) {
            let mate_ref_name = Self::mate_reference_name(&record, bam_reader);
            let passes_sam_filter = sam_filter == 0 || (record.flags() & (sam_filter as u16)) == 0;
            let passed = self.passes_preprocess(
                &record,
                sam_filter,
                &mut preprocess_state,
                &mate_ref_name,
            );
            if let Some(entries) = jsonl_entries.as_mut() {
                if !passes_sam_filter {
                    continue;
                }
                let qname = String::from_utf8_lossy(record.qname()).to_string();
                let alignment_start = record.pos() + 1;
                let mate_start = if record.mpos() >= 0 { record.mpos() + 1 } else { 0 };
                let cigar = record.cigar().to_string();
                let seq_bytes = record.seq().as_bytes();
                let seq_str = String::from_utf8_lossy(&seq_bytes).to_string();
                let qual_str: String = record
                    .qual()
                    .iter()
                    .map(|q| (*q as u8 + 33) as char)
                    .collect();
                let data = format!(
                    "{{\"passed\":{},\"flag\":{},\"pos\":{},\"mpos\":{},\"mapq\":{},\"cigar\":\"{}\",\"mateRef\":\"{}\",\"sequence\":\"{}\",\"quality\":\"{}\",\"totalReads\":{},\"duplicateReads\":{}}}",
                    passed,
                    record.flags(),
                    alignment_start,
                    mate_start,
                    record.mapq(),
                    json_escape(&cigar),
                    json_escape(&mate_ref_name),
                    json_escape(&seq_str),
                    json_escape(&qual_str),
                    preprocess_state.total_reads,
                    preprocess_state.duplicate_reads,
                );
                entries.push(RecordPreprocessorJsonlEntry {
                    pos: alignment_start,
                    key: qname,
                    data,
                });
            }
            if !passed {
                continue;
            }

            cigar_parser.process_record(&mut record)?;
        }

        if let (Some(path), Some(entries)) = (jsonl_path.as_ref(), jsonl_entries.as_ref()) {
            write_record_preprocessor_jsonl_snapshot(
                path,
                region,
                preprocess_state.total_reads,
                preprocess_state.duplicate_reads,
                entries,
            )?;
        }

        Ok(self.build_cigar_output(
            &mut cigar_parser,
            preprocess_state.total_reads,
            preprocess_state.duplicate_reads,
        ))
    }

    fn run_cigar_parser_from_bam_paths(
        &self,
        region: &Region,
        reference: &Reference,
        instance: Arc<GlobalReadOnlyScope>,
        bam_paths: &[String],
        sam_filter: u32,
    ) -> Result<CigarParserOutput> {
        let mut cigar_parser = CigarParser::new(
            region.clone(),
            reference.clone(),
            instance,
        );

        let mut preprocess_state = RecordPreprocessorState::new();
        for bam_path in bam_paths {
            let mut bam_reader = BamReader::open(bam_path).map_err(|error| {
                anyhow::anyhow!(
                    "Failed to open BAM for combined processing {}: {}",
                    bam_path,
                    error
                )
            })?;

            bam_reader.fetch(region.chr(), region.start(), region.end())?;

            let mut record = Record::new();
            while bam_reader.read(&mut record).unwrap_or(false) {
                let mate_ref_name = Self::mate_reference_name(&record, &bam_reader);
                if !self.passes_preprocess(
                    &record,
                    sam_filter,
                    &mut preprocess_state,
                    &mate_ref_name,
                ) {
                    continue;
                }
                cigar_parser.process_record(&mut record)?;
            }
        }

        Ok(self.build_cigar_output(
            &mut cigar_parser,
            preprocess_state.total_reads,
            preprocess_state.duplicate_reads,
        ))
    }

    pub(crate) fn run_partial_cigar_for_bams(
        &self,
        region: &Region,
        reference: &Reference,
        bam_paths: &[String],
    ) -> Result<CigarParserOutput> {
        let scope = Arc::new(instance().clone());
        let sam_filter = instance().conf.sam_filter;
        self.run_cigar_parser_from_bam_paths(region, reference, scope, bam_paths, sam_filter)
    }


    fn build_cigar_output(
        &self,
        cigar_parser: &mut CigarParser,
        total_reads: usize,
        duplicate_reads: usize,
    ) -> CigarParserOutput {
        let splice_count_raw = cigar_parser.take_splice_count();
        let mut splice: HashSet<String> = HashSet::new();
        let mut splice_count: HashMap<String, usize> = HashMap::new();
        for ((start, end), counts) in splice_count_raw.into_iter() {
            let key = format!("{}-{}", start, end);
            splice.insert(key.clone());
            let count = counts.get(0).copied().unwrap_or(0);
            splice_count.insert(key, count);
        }

        let duprate = if instance().conf.remove_duplicated_reads && total_reads != 0 {
            (duplicate_reads as f64 / total_reads as f64 * 1000.0).round() / 1000.0
        } else {
            0.0
        };

        CigarParserOutput {
            non_insertion_vars: cigar_parser.take_non_insertion_vars(),
            non_insertion_vars_insert_index: cigar_parser.take_non_insertion_vars_insert_index(),
            insertion_vars: cigar_parser.take_insertion_vars(),
            soft_clips_5end: cigar_parser.take_soft_clips_5end(),
            soft_clips_3end: cigar_parser.take_soft_clips_3end(),
            ref_coverage: cigar_parser.take_ref_coverage(),
            mnp: cigar_parser.take_mnp(),
            position_to_insertion_count: cigar_parser.take_position_to_insertion_count(),
            position_to_deletions_count: cigar_parser.take_position_to_deletions_count(),
            max_read_len: cigar_parser.get_max_read_len(),
            discordant_count: cigar_parser.get_discordant_count(),
            splice,
            splice_count,
            duprate,
            total_reads,
            duplicate_reads,
        }
    }

    /// Step 2: Run VariantRealigner and StructuralVariantsProcessor
    fn run_variant_realigner_and_sv_processor(
        &self,
        input: CigarParserOutput,
        region: &Region,
        reference: &Reference,
        bam_paths: &[String],
    ) -> Result<RealignedOutput> {
        // TODO: Integrate actual VariantRealigner for soft clip realignment
        // For now, we pass through to StructuralVariantsProcessor

        let CigarParserOutput {
            non_insertion_vars,
            non_insertion_vars_insert_index,
            insertion_vars,
            soft_clips_5end,
            soft_clips_3end,
            ref_coverage,
            mnp,
            position_to_insertion_count,
            position_to_deletions_count,
            max_read_len,
            duprate,
            splice,
            ..
        } = input;

        // Convert CigarParserOutput to RealignedVariationData for SV processor
        let mut sv_input = RealignedVariationData {
            non_insertion_variants: non_insertion_vars,
            insertion_variants: insertion_vars,
            soft_clips_5end,
            soft_clips_3end,
            ref_coverage,
            max_read_length: max_read_len,
            duprate,
        };

        // Perform minimal deletion realignment using soft clips when enabled
        // Re-enable realigner to match Java behavior
        let realigner = VariantRealigner::new_with_context(
            reference.ref_seq.clone(),
            reference.seed.clone(),
            reference.region_start,
            Some(region.chr().to_string()),
            bam_paths.to_vec(),
        );
        realigner.adjust_mnp(&mut sv_input, &mnp);

        if instance().conf.perform_local_realignment {
            realigner.process_deletions(&mut sv_input, &position_to_deletions_count);
            realigner.process_insertions(&mut sv_input, &position_to_insertion_count);
            realigner.realign_long_insertions_30(&mut sv_input);
            realigner.realign_long_insertions(&mut sv_input);
        }

        write_realigned_jsonl_snapshot_if_enabled(&sv_input, region)?;
        
        // Run StructuralVariantsProcessor (adjSNV always runs, SV detection is unimplemented)
        let sv_processor = StructuralVariantsProcessor::new(
            reference.ref_seq.clone(),
            reference.region_start,
        );
        let processed = sv_processor.process(sv_input);

        write_structural_variants_jsonl_snapshot_if_enabled(&processed, region)?;
        
        // Convert back to RealignedOutput
        Ok(RealignedOutput {
            non_insertion_vars: processed.non_insertion_variants,
            non_insertion_vars_insert_index,
            insertion_vars: processed.insertion_variants,
            ref_coverage: processed.ref_coverage,
            duprate: processed.duprate,
            max_read_len: processed.max_read_length,
            splice,
        })
    }

    /// Step 3: Run ToVarsBuilder to calculate statistics
    fn run_to_vars_builder(
        &self,
        input: RealignedOutput,
        reference: &Reference,
        region: &Region,
    ) -> Result<AlignedVarsData> {
        let mut aligned_variants: HashMap<i64, Vars> = HashMap::new();
        let mut aligned_variants_order: Vec<i64> = Vec::new();
        let RealignedOutput {
            non_insertion_vars,
            non_insertion_vars_insert_index,
            insertion_vars,
            ref_coverage,
            duprate,
            max_read_len,
            ..
        } = input;
        let mut non_insertion_vars = non_insertion_vars;
        let mut insertion_vars = insertion_vars;

        let debug_pos = env::var("VARDICT_DEBUG_POS")
            .ok()
            .and_then(|value| value.parse::<i64>().ok());

        if let Some(pos) = debug_pos {
            event!(
                Level::DEBUG,
                debug_pos = pos,
                nonins_has_pos = non_insertion_vars.contains_key(&pos),
                ins_has_pos = insertion_vars.contains_key(&pos),
                refcov_has_pos = ref_coverage.contains_key(&pos),
                nonins_len = non_insertion_vars.len(),
                ins_len = insertion_vars.len(),
                "to_vars_builder: debug position presence before iteration"
            );
        }

        let mut position_keys: Vec<i64> = non_insertion_vars.keys().copied().collect();
        let mut seen_positions: HashSet<i64> = position_keys.iter().copied().collect();
        for pos in insertion_vars.keys().copied() {
            if seen_positions.insert(pos) {
                position_keys.push(pos);
            }
        }

        let positions = java_hashmap_iteration_order(
            position_keys.into_iter(),
            seen_positions.len(),
            Some(&non_insertion_vars_insert_index),
        );

        for position in positions {
            let trace_this_pos = debug_pos == Some(position);
            let vars_at_pos = non_insertion_vars
                .get(&position)
                .cloned()
                .unwrap_or_default();

            if trace_this_pos {
                event!(
                    Level::DEBUG,
                    position,
                    nonins_count = vars_at_pos.len(),
                    has_insertion = insertion_vars.contains_key(&position),
                    has_refcov = ref_coverage.contains_key(&position),
                    "to_vars_builder: position encountered"
                );
            }

            if vars_at_pos.is_empty() && !insertion_vars.contains_key(&position) {
                if trace_this_pos {
                    event!(
                        Level::DEBUG,
                        position,
                        "to_vars_builder: skipped because no non-insertion and no insertion variants"
                    );
                }
                continue;
            }

            if position < region.start() as i64 || position > region.end() as i64 {
                if trace_this_pos {
                    event!(
                        Level::DEBUG,
                        position,
                        region_start = region.start(),
                        region_end = region.end(),
                        "to_vars_builder: skipped because position outside region"
                    );
                }
                continue;
            }

            if !ref_coverage.contains_key(&position) {
                if trace_this_pos {
                    event!(
                        Level::DEBUG,
                        position,
                        "to_vars_builder: skipped because reference coverage missing"
                    );
                }
                continue;
            }

            if self.is_same_variation_on_ref(
                position,
                &vars_at_pos,
                insertion_vars.get(&position),
                reference,
            ) {
                if trace_this_pos {
                    event!(
                        Level::DEBUG,
                        position,
                        "to_vars_builder: skipped because only reference variation"
                    );
                }
                continue;
            }

            let mut total_pos_coverage = match ref_coverage.get(&position) {
                Some(coverage) if *coverage > 0 => *coverage,
                _ => continue,
            };

            if trace_this_pos {
                event!(
                    Level::DEBUG,
                    position,
                    ref_cov_at_pos = total_pos_coverage,
                    ref_cov_next = ref_coverage.get(&(position + 1)).copied().unwrap_or(0),
                    "to_vars_builder: reference coverage snapshot"
                );
            }

            let hicov = self.calc_hicov(insertion_vars.get(&position), &vars_at_pos);

            let mut var_list: Vec<Variant> = Vec::new();
            let mut debug_lines: Vec<String> = Vec::new();

            let mut keys: Vec<VarDesc> = vars_at_pos.keys().cloned().collect();
            keys.sort_by(|a, b| a.to_key_string().cmp(&b.to_key_string()));

            let sv_string = self.create_variant_records(
                position,
                &vars_at_pos,
                total_pos_coverage,
                &mut var_list,
                &mut debug_lines,
                &keys,
                hicov,
                duprate,
            );

            total_pos_coverage = self.create_insertion_records(
                position,
                total_pos_coverage,
                insertion_vars.get(&position),
                &mut non_insertion_vars,
                &ref_coverage,
                reference,
                &mut var_list,
                &mut debug_lines,
                hicov,
                duprate,
            );

            if trace_this_pos {
                let pre_sort_summary = var_list
                    .iter()
                    .map(|variant| {
                        format!(
                            "{}|pcov={}|tot={}|fwd={}|rev={}|freq={:.4}",
                            variant.description_string,
                            variant.position_coverage,
                            variant.total_pos_coverage,
                            variant.vars_count_on_forward,
                            variant.vars_count_on_reverse,
                            variant.frequency
                        )
                    })
                    .collect::<Vec<_>>()
                    .join(",");
                event!(
                    Level::DEBUG,
                    position,
                    variant_count = var_list.len(),
                    total_pos_coverage,
                    pre_sort_summary = %pre_sort_summary,
                    "to_vars_builder: variants created before sorting"
                );
            }

            self.sort_variants(&mut var_list);

            let maxfreq = self.collect_vars_at_position(
                &mut aligned_variants,
                &mut aligned_variants_order,
                position,
                reference,
                &var_list,
            );

            if let Some(sv) = sv_string {
                if let Some(vars_entry) = aligned_variants.get_mut(&position) {
                    vars_entry.sv = sv;
                }
            }

            if !self.do_pileup
                && maxfreq <= instance().conf.freq
                && instance().amplicon_based_calling.is_none()
            {
                if instance().bam_paths.len() < 2 {
                    if trace_this_pos {
                        event!(
                            Level::DEBUG,
                            position,
                            maxfreq,
                            min_freq = instance().conf.freq,
                            "to_vars_builder: removing position due to maxfreq threshold"
                        );
                    }
                    aligned_variants.remove(&position);
                    continue;
                }
            }

            if let Some(variations_at_pos) = aligned_variants.get_mut(&position) {
                self.collect_reference_variants(
                    position,
                    total_pos_coverage,
                    variations_at_pos,
                    &ref_coverage,
                    &mut non_insertion_vars,
                    reference,
                    region,
                    &mut debug_lines,
                    duprate,
                );
                if trace_this_pos {
                    let post_ref_summary = variations_at_pos
                        .variants
                        .iter()
                        .map(|variant| {
                            format!(
                                "{}|{}>{}|pcov={}|tot={}|fwd={}|rev={}|freq={:.4}",
                                variant.description_string,
                                variant.refallele,
                                variant.varallele,
                                variant.position_coverage,
                                variant.total_pos_coverage,
                                variant.vars_count_on_forward,
                                variant.vars_count_on_reverse,
                                variant.frequency
                            )
                        })
                        .collect::<Vec<_>>()
                        .join(",");
                    let reference_variant_cov = variations_at_pos
                        .reference_variant
                        .as_ref()
                        .map(|variant| variant.total_pos_coverage)
                        .unwrap_or(0);
                    event!(
                        Level::DEBUG,
                        position,
                        variant_count_after_ref = variations_at_pos.variants.len(),
                        has_reference_variant = variations_at_pos.reference_variant.is_some(),
                        reference_variant_cov,
                        post_ref_summary = %post_ref_summary,
                        "to_vars_builder: position retained after reference collection"
                    );
                }
            }
        }

        let aligned_data = AlignedVarsData {
            aligned_variants,
            aligned_variants_order,
            ref_coverage,
        };

        write_tovars_jsonl_snapshot_if_enabled(&aligned_data, region, max_read_len, duprate)?;

        Ok(aligned_data)
    }

    fn is_same_variation_on_ref(
        &self,
        position: i64,
        vars_at_pos: &HashMap<VarDesc, RawVariant>,
        insertion_vars: Option<&HashMap<VarDesc, RawVariant>>,
        reference: &Reference,
    ) -> bool {
        let mut keys = HashSet::new();
        for desc in vars_at_pos.keys() {
            keys.insert(desc.to_key_string());
        }

        if insertion_vars.is_some() {
            keys.insert("I".to_string());
        }

        if keys.len() == 1 {
            if let Some(ref_base) = reference.get(position).map(|b| (b as char).to_string()) {
                if keys.contains(&ref_base)
                    && !self.do_pileup
                    && instance().bam_paths.len() < 2
                    && instance().amplicon_based_calling.is_none()
                {
                    return true;
                }
            }
        }

        false
    }

    fn calc_hicov(
        &self,
        _insertion_vars: Option<&HashMap<VarDesc, RawVariant>>,
        non_insertion_vars: &HashMap<VarDesc, RawVariant>,
    ) -> usize {
        let mut hicov = 0usize;
        for (desc, raw_var) in non_insertion_vars {
            let is_sv = matches!(desc, VarDesc::Raw { desc } if desc.as_slice() == b"SV");
            let is_insertion_like = matches!(desc, VarDesc::Ins { .. })
                || matches!(desc, VarDesc::Raw { desc } if desc.as_slice().starts_with(b"+"));
            if is_sv || is_insertion_like {
                continue;
            }
            hicov += raw_var.high_qual_read_cnt;
        }
        hicov
    }

    fn create_variant_records(
        &self,
        position: i64,
        vars_at_pos: &HashMap<VarDesc, RawVariant>,
        total_pos_coverage: usize,
        var_list: &mut Vec<Variant>,
        _debug_lines: &mut Vec<String>,
        keys: &[VarDesc],
        hicov: usize,
        duprate: f64,
    ) -> Option<String> {
        use crate::mods::to_vars_builder::{check_strand_bias, StrandBiasFlag, StrandBiasValue};

        let mut sv_string: Option<String> = None;

        for desc in keys {
            if matches!(desc, VarDesc::Raw { desc } if desc.as_slice() == b"SV") {
                if let Some(sv_var) = vars_at_pos.get(desc) {
                    sv_string = Some(format!("{}-0-0", sv_var.alt_depth));
                }
                continue;
            }

            let Some(raw_var) = vars_at_pos.get(desc) else {
                continue;
            };
            let fwd = raw_var.alt_depth_fwd;
            let rev = raw_var.alt_depth_rev;
            let total_count = if raw_var.alt_depth > 0 {
                raw_var.alt_depth
            } else {
                fwd + rev
            };
            if total_count == 0 {
                continue;
            }
            let bias = check_strand_bias(fwd, rev);

            let base_quality = round_half_even("0.0", raw_var.mean_qual / total_count as f64);
            let mapping_quality = round_half_even("0.0", raw_var.mean_mapq / total_count as f64);
            let hicnt = raw_var.high_qual_read_cnt;
            let locnt = raw_var.low_qual_read_cnt;

            let mut ttcov = total_pos_coverage;
            if total_count > total_pos_coverage
                && raw_var.extra_cnt > 0
                && total_count - total_pos_coverage < raw_var.extra_cnt
            {
                ttcov = total_count;
            }

            let mut variant = Variant::new();
            variant.description_string = desc.to_key_string();
            variant.position_coverage = total_count;
            variant.vars_count_on_forward = fwd;
            variant.vars_count_on_reverse = rev;
            variant.strand_bias_flag = StrandBiasFlag::new(StrandBiasValue::CantAssess, bias);
            variant.frequency = round_half_even("0.0000", total_count as f64 / ttcov as f64);
            variant.mean_position = round_half_even("0.0", raw_var.mean_pos / total_count as f64);
            variant.is_at_least_at_2_positions = raw_var.pstd;
            variant.mean_quality = base_quality;
            variant.has_at_least_2_diff_qualities = raw_var.qstd;
            variant.mean_mapping_quality = mapping_quality;
            variant.high_quality_reads_frequency = if hicov > 0 {
                round_half_even("0.0000", hicnt as f64 / hicov as f64)
            } else {
                0.0
            };
            variant.extra_frequency = if raw_var.extra_cnt > 0 {
                round_half_even("0.0000", raw_var.extra_cnt as f64 / ttcov as f64)
            } else {
                0.0
            };
            variant.shift3 = 0;
            variant.msi = 0.0;
            variant.nm = round_half_even("0.0", raw_var.nm / total_count as f64);
            variant.high_qual_read_cnt = hicnt;
            variant.low_qual_read_cnt = locnt;
            variant.hicov = hicov;
            variant.duprate = duprate;
            variant.start_position = position;
            variant.end_position = position;

            var_list.push(variant);
        }

        sv_string
    }

    fn create_insertion_records(
        &self,
        position: i64,
        mut total_pos_coverage: usize,
        insertion_vars: Option<&HashMap<VarDesc, RawVariant>>,
        non_insertion_vars: &mut HashMap<i64, HashMap<VarDesc, RawVariant>>,
        ref_coverage: &HashMap<i64, usize>,
        reference: &Reference,
        var_list: &mut Vec<Variant>,
        _debug_lines: &mut Vec<String>,
        hicov: usize,
        duprate: f64,
    ) -> usize {
        use crate::mods::to_vars_builder::{check_strand_bias, StrandBiasFlag, StrandBiasValue};

        let Some(insertion_variations) = insertion_vars else {
            return total_pos_coverage;
        };

        let mut running_hicov = hicov;

        let mut keys: Vec<VarDesc> = insertion_variations.keys().cloned().collect();
        keys.sort_by(|a, b| a.to_key_string().cmp(&b.to_key_string()));

        for desc in keys {
            let desc_str = desc.to_key_string();
            if desc_str.contains('&') {
                let coverage_position = if let Some(without_plus) = desc_str.strip_prefix('+') {
                    if let Some((prefix, _)) = without_plus.split_once('&') {
                        position + prefix.len() as i64
                    } else {
                        position + 1
                    }
                } else {
                    position + 1
                };
                if let Some(&coverage) = ref_coverage.get(&coverage_position) {
                    total_pos_coverage = coverage;
                } else if coverage_position != position + 1 {
                    if let Some(&fallback_coverage) = ref_coverage.get(&(position + 1)) {
                        total_pos_coverage = fallback_coverage;
                    }
                }
            }

            let Some(cnt) = insertion_variations.get(&desc) else {
                continue;
            };

            let fwd = cnt.alt_depth_fwd;
            let rev = cnt.alt_depth_rev;
            let bias = check_strand_bias(fwd, rev);
            let total_count = if cnt.alt_depth > 0 { cnt.alt_depth } else { fwd + rev };
            if total_count == 0 {
                continue;
            }

            let vqual = round_half_even("0.0", cnt.mean_qual / total_count as f64);
            let mq = round_half_even("0.0", cnt.mean_mapq / total_count as f64);
            let hicnt = cnt.high_qual_read_cnt;
            let locnt = cnt.low_qual_read_cnt;
            if running_hicov < hicnt {
                running_hicov = hicnt;
            }

            let mut ttcov = total_pos_coverage;
            if total_count > total_pos_coverage
                && cnt.extra_cnt != 0
                && total_count - total_pos_coverage < cnt.extra_cnt
            {
                ttcov = total_count;
            }

            if ttcov < total_count {
                ttcov = total_count;
                if let Some(&next_cov) = ref_coverage.get(&(position + 1)) {
                    if ttcov < next_cov.saturating_sub(total_count) {
                        ttcov = next_cov;
                        if let Some(next_map) = non_insertion_vars.get_mut(&(position + 1)) {
                            if let Some(ref_base) = reference.get(position + 1) {
                                let key = VarDesc::SNV { ref_base };
                                if let Some(next_var) = next_map.get_mut(&key) {
                                    next_var.alt_depth_fwd =
                                        next_var.alt_depth_fwd.saturating_sub(fwd);
                                    next_var.alt_depth_rev =
                                        next_var.alt_depth_rev.saturating_sub(rev);
                                }
                            }
                        }
                    }
                }
                total_pos_coverage = ttcov;
            }

            let mut variant = Variant::new();
            variant.description_string = desc_str;
            variant.position_coverage = total_count;
            variant.vars_count_on_forward = fwd;
            variant.vars_count_on_reverse = rev;
            variant.strand_bias_flag = StrandBiasFlag::new(StrandBiasValue::CantAssess, bias);
            variant.frequency = round_half_even("0.0000", total_count as f64 / ttcov as f64);
            variant.mean_position = round_half_even("0.0", cnt.mean_pos / total_count as f64);
            variant.is_at_least_at_2_positions = cnt.pstd;
            variant.mean_quality = vqual;
            variant.has_at_least_2_diff_qualities = cnt.qstd;
            variant.mean_mapping_quality = mq;
            variant.high_quality_reads_frequency = if running_hicov > 0 {
                round_half_even("0.0000", hicnt as f64 / running_hicov as f64)
            } else {
                0.0
            };
            variant.extra_frequency = if cnt.extra_cnt != 0 {
                round_half_even("0.0000", cnt.extra_cnt as f64 / ttcov as f64)
            } else {
                0.0
            };
            variant.shift3 = 0;
            variant.msi = 0.0;
            variant.nm = round_half_even("0.0", cnt.nm / total_count as f64);
            variant.high_qual_read_cnt = hicnt;
            variant.low_qual_read_cnt = locnt;
            variant.hicov = running_hicov;
            variant.duprate = duprate;
            variant.start_position = position;
            variant.end_position = position;

            var_list.push(variant);
        }

        total_pos_coverage
    }

    fn sort_variants(&self, variants: &mut [Variant]) {
        variants.sort_by(|a, b| {
            let a_score = a.mean_quality * a.position_coverage as f64;
            let b_score = b.mean_quality * b.position_coverage as f64;
            match b_score.partial_cmp(&a_score).unwrap_or(std::cmp::Ordering::Equal) {
                std::cmp::Ordering::Equal => a.description_string.cmp(&b.description_string),
                other => other,
            }
        });
    }

    fn collect_vars_at_position(
        &self,
        aligned_variants: &mut HashMap<i64, Vars>,
        aligned_variants_order: &mut Vec<i64>,
        position: i64,
        reference: &Reference,
        variants: &[Variant],
    ) -> f64 {
        let mut maxfreq = 0.0;
        let ref_base = reference.get(position).map(|b| (b as char).to_string());
        let entry = match aligned_variants.entry(position) {
            Entry::Occupied(entry) => entry.into_mut(),
            Entry::Vacant(entry) => {
                aligned_variants_order.push(position);
                entry.insert(Vars::default())
            }
        };

        for variant in variants {
            if let Some(ref_base) = &ref_base {
                if variant.description_string == *ref_base {
                    entry.reference_variant = Some(variant.clone());
                    continue;
                }
            }

            if variant.frequency > maxfreq {
                maxfreq = variant.frequency;
            }
            entry.variants.push(variant.clone());
        }

        maxfreq
    }

    fn collect_reference_variants(
        &self,
        position: i64,
        mut total_pos_coverage: usize,
        variations_at_pos: &mut Vars,
        ref_coverage: &HashMap<i64, usize>,
        non_insertion_vars: &mut HashMap<i64, HashMap<VarDesc, RawVariant>>,
        reference: &Reference,
        region: &Region,
        _debug_lines: &mut Vec<String>,
        duprate: f64,
    ) {
        use crate::data::patterns::{
            AMP_ATGC, ANY_SV, BEGIN_DIGITS, BEGIN_MINUS_NUMBER, BEGIN_MINUS_NUMBER_CARET,
            CARET_ATGNC, DUP_NUM, HASH_GROUP_CARET_GROUP, INV_NUM, SOME_SV_NUMBERS,
        };
        use crate::mods::to_vars_builder::{StrandBiasFlag, StrandBiasValue};

        let mut reference_forward_coverage = 0usize;
        let mut reference_reverse_coverage = 0usize;

        let mut genotype1 = if let Some(ref_var) = &variations_at_pos.reference_variant {
            if ref_var.frequency >= instance().conf.freq {
                ref_var.description_string.clone()
            } else if !variations_at_pos.variants.is_empty() {
                variations_at_pos.variants[0].description_string.clone()
            } else {
                ref_var.description_string.clone()
            }
        } else if !variations_at_pos.variants.is_empty() {
            variations_at_pos.variants[0].description_string.clone()
        } else {
            String::new()
        };

        if let Some(ref_var) = &variations_at_pos.reference_variant {
            reference_forward_coverage = ref_var.vars_count_on_forward;
            reference_reverse_coverage = ref_var.vars_count_on_reverse;
        }

        if genotype1.starts_with('+') {
            if let Some(caps) = DUP_NUM.captures(&genotype1) {
                if let Ok(dup_len) = caps.get(1).map(|m| m.as_str()).unwrap_or("0").parse::<i32>() {
                    genotype1 = format!("+{}", crate::conf::Configuration::SVFLANK + dup_len);
                }
            } else if genotype1.len() > 1 {
                genotype1 = format!("+{}", genotype1.len() - 1);
            }
        }

        if let Some(&ref_cov_at_pos) = ref_coverage.get(&position) {
            if total_pos_coverage > ref_cov_at_pos {
                if let Some(next_map) = non_insertion_vars.get(&(position + 1)) {
                    if let Some(ref_base) = reference.get(position + 1) {
                        let key = VarDesc::SNV { ref_base };
                        if let Some(tpref) = next_map.get(&key) {
                            reference_forward_coverage = tpref.alt_depth_fwd;
                            reference_reverse_coverage = tpref.alt_depth_rev;
                        }
                    }
                }
            }
        }

        let mut positions_for_changed_ref_variant: Vec<i64> = Vec::new();

        if !variations_at_pos.variants.is_empty() {
            for vref in variations_at_pos.variants.iter_mut() {
                let mut genotype1current = genotype1.clone();
                let mut genotype2 = vref.description_string.clone();

                if genotype2.starts_with('+') {
                    genotype2 = format!("+{}", genotype2.len().saturating_sub(1));
                }

                let description_string = vref.description_string.clone();
                let mut deletion_length = 0usize;
                if let Some(caps) = BEGIN_MINUS_NUMBER.captures(&description_string) {
                    if let Ok(val) = caps.get(1).map(|m| m.as_str()).unwrap_or("0").parse::<usize>() {
                        deletion_length = val;
                    }
                }

                let mut end_position = position;
                if description_string.starts_with('-') && deletion_length > 0 {
                    end_position = position + deletion_length as i64 - 1;
                }

                let mut refallele = String::new();
                let mut varallele = String::new();
                let mut shift3 = 0i32;
                let mut msi = 0.0;
                let mut msint = 0.0;
                let mut start_position = position;

                if description_string.starts_with('+') {
                    if !description_string.contains('&')
                        && !description_string.contains('#')
                        && !description_string.to_ascii_lowercase().contains("<dup")
                    {
                        let (msi_val, shift_val, msint_val) = self.proceed_vref_is_insertion(
                            position,
                            &description_string,
                            reference,
                            region,
                        );
                        msi = msi_val;
                        shift3 = shift_val;
                        msint = msint_val;
                    }

                    if instance().conf.move_indels_to_3 {
                        start_position += shift3 as i64;
                        end_position += shift3 as i64;
                    }

                    refallele = reference
                        .get(position)
                        .map(|b| (b as char).to_string())
                        .unwrap_or_default();
                    varallele = format!("{}{}", refallele, description_string.trim_start_matches('+'));

                    if varallele.len() > instance().conf.sv_min_len {
                        end_position += varallele.len() as i64;
                        varallele = "<DUP>".to_string();
                    }

                    if let Some(caps) = DUP_NUM.captures(&varallele) {
                        if let Ok(dup_count) = caps.get(1).map(|m| m.as_str()).unwrap_or("0").parse::<i32>() {
                            end_position = start_position + (2 * crate::conf::Configuration::SVFLANK + dup_count) as i64 - 1;
                            genotype2 = format!("+{}", 2 * crate::conf::Configuration::SVFLANK + dup_count);
                            varallele = "<DUP>".to_string();
                        }
                    }
                } else if description_string.starts_with('-') {
                    let matcher_inv = INV_NUM.captures(&description_string);
                    let matcher_start_minus = BEGIN_MINUS_NUMBER_CARET.is_match(&description_string);

                    if deletion_length < instance().conf.sv_min_len {
                        if deletion_length > 0 {
                            let prefix = format!("-{}", deletion_length);
                            varallele = description_string.replacen(&prefix, "", 1);
                        } else {
                            varallele = description_string.clone();
                        }

                        let (msi_val, shift_val, msint_val) =
                            self.proceed_vref_is_deletion(position, deletion_length, reference);
                        msi = msi_val;
                        shift3 = shift_val;
                        msint = msint_val;

                        if matcher_inv.is_some() {
                            varallele = "<INV>".to_string();
                            genotype2 = format!("<INV{}>", deletion_length);
                        }
                    } else if matcher_start_minus {
                        varallele = "<INV>".to_string();
                        genotype2 = format!("<INV{}>", deletion_length);
                    } else {
                        varallele = "<DEL>".to_string();
                    }

                    if !description_string.contains('&')
                        && !description_string.contains('#')
                        && !description_string.contains('^')
                    {
                        if instance().conf.move_indels_to_3 {
                            start_position += shift3 as i64;
                        }
                        if varallele != "<DEL>" {
                            if let Some(base) = reference.get(position - 1) {
                                varallele = (base as char).to_string();
                            }
                        }
                        if let Some(base) = reference.get(position - 1) {
                            refallele.push(base as char);
                        }
                        start_position -= 1;
                    }

                    if SOME_SV_NUMBERS.is_match(&description_string) {
                        refallele = reference
                            .get(position)
                            .map(|b| (b as char).to_string())
                            .unwrap_or_default();
                    } else if deletion_length < instance().conf.sv_min_len {
                        refallele.push_str(&self.get_reference_range(
                            reference,
                            position,
                            position + deletion_length as i64 - 1,
                        ));
                    }
                } else {
                    let (msi_val, msint_val, shift_val) =
                        self.detect_microsatellite_snp(reference, position);
                    msi = msi_val;
                    msint = msint_val;
                    shift3 = shift_val;

                    refallele = reference
                        .get(position)
                        .map(|b| (b as char).to_string())
                        .unwrap_or_default();
                    varallele = description_string.clone();
                }

                if let Some(caps) = AMP_ATGC.captures(&description_string) {
                    let extra = caps.get(1).map(|m| m.as_str()).unwrap_or("");
                    if !extra.is_empty() {
                        varallele = varallele.replacen('&', "", 1);
                        let tch = self.get_reference_range(
                            reference,
                            end_position + 1,
                            end_position + extra.len() as i64,
                        );
                        refallele.push_str(&tch);
                        genotype1current.push_str(&tch);
                        end_position += extra.len() as i64;

                        let varallele_match = varallele.clone();
                        if let Some(caps2) = AMP_ATGC.captures(&varallele_match) {
                            let vextra = caps2.get(1).map(|m| m.as_str()).unwrap_or("");
                            if !vextra.is_empty() {
                                varallele = varallele.replacen('&', "", 1);
                                let tch2 = self.get_reference_range(
                                    reference,
                                    end_position + 1,
                                    end_position + vextra.len() as i64,
                                );
                                refallele.push_str(&tch2);
                                genotype1current.push_str(&tch2);
                                end_position += vextra.len() as i64;
                            }
                        }

                        if description_string.starts_with('+') {
                            if !refallele.is_empty() {
                                refallele = refallele[1..].to_string();
                            }
                            if !varallele.is_empty() {
                                varallele = varallele[1..].to_string();
                            }
                            start_position += 1;
                        }

                        if varallele == "<DEL>" && !refallele.is_empty() {
                            refallele = reference
                                .get(start_position)
                                .map(|b| (b as char).to_string())
                                .unwrap_or_default();
                            if let Some(&coverage) = ref_coverage.get(&(start_position - 1)) {
                                total_pos_coverage = coverage;
                            }
                            if vref.position_coverage > total_pos_coverage {
                                total_pos_coverage = vref.position_coverage;
                            }
                            vref.frequency = vref.position_coverage as f64 / total_pos_coverage as f64;
                        }
                    }
                }

                if let Some(caps) = HASH_GROUP_CARET_GROUP.captures(&description_string) {
                    let matched_seq = caps.get(1).map(|m| m.as_str()).unwrap_or("");
                    let tail = caps.get(2).map(|m| m.as_str()).unwrap_or("");

                    end_position += matched_seq.len() as i64;
                    refallele.push_str(&self.get_reference_range(
                        reference,
                        end_position - matched_seq.len() as i64 + 1,
                        end_position,
                    ));

                    if let Some(digits) = BEGIN_DIGITS.captures(tail) {
                        if let Ok(deletion) = digits.get(1).map(|m| m.as_str()).unwrap_or("0").parse::<i64>() {
                            refallele.push_str(&self.get_reference_range(
                                reference,
                                end_position + 1,
                                end_position + deletion,
                            ));
                            end_position += deletion;
                        }
                    }

                    varallele = varallele.replacen('#', "", 1);
                    varallele = remove_caret_and_digits(&varallele);
                    genotype1current = genotype1current.replace('#', "m").replace('^', "i");
                    genotype2 = genotype2.replace('#', "m").replace('^', "i");
                }

                if CARET_ATGNC.is_match(&description_string) {
                    varallele = varallele.replacen('^', "", 1);
                    genotype1current = genotype1current.replace('^', "i");
                    genotype2 = genotype2.replace('^', "i");
                }

                let cut_site = instance().conf.crispr_cutting_site as i64;
                if cut_site != 0 && refallele.len() > 1 && varallele.len() > 1 {
                    let mut n = 0usize;
                    while refallele.len() > n + 1
                        && varallele.len() > n + 1
                        && refallele.as_bytes()[n] == varallele.as_bytes()[n]
                    {
                        n += 1;
                    }
                    if n != 0 {
                        start_position += n as i64;
                        refallele = refallele[n..].to_string();
                        varallele = varallele[n..].to_string();
                    }
                }

                if cut_site != 0
                    && refallele.len() != varallele.len()
                    && refallele.chars().next() == varallele.chars().next()
                {
                    if start_position != cut_site && end_position != cut_site {
                        let mut n = 0i64;
                        let dis = (cut_site - start_position).abs().min((cut_site - end_position).abs());
                        if start_position < cut_site {
                            while start_position + n < cut_site
                                && n < shift3 as i64
                                && end_position + n != cut_site
                            {
                                n += 1;
                            }
                            if (start_position + n - cut_site).abs() > dis
                                && (end_position + n - cut_site).abs() > dis
                            {
                                n = 0;
                            }
                        }
                        if end_position < cut_site && n == 0 {
                            if (end_position - cut_site).abs() <= (start_position - cut_site).abs() {
                                while end_position + n < cut_site && n < shift3 as i64 {
                                    n += 1;
                                }
                            }
                        }
                        if n > 0 {
                            start_position += n;
                            end_position += n;
                            refallele.clear();
                            for pos in start_position..=end_position {
                                if let Some(base) = reference.get(pos) {
                                    refallele.push(base as char);
                                }
                            }
                            let mut tva = String::new();
                            if refallele.len() < varallele.len() {
                                tva = varallele[1..].to_string();
                                if tva.len() > 1 {
                                    let ttn = (n as usize) % tva.len();
                                    if ttn != 0 {
                                        tva = format!("{}{}", &tva[ttn..], &tva[..ttn]);
                                    }
                                }
                            }
                            if let Some(base) = reference.get(start_position) {
                                varallele = format!("{}{}", base as char, tva);
                            }
                            vref.crispr = n as i32;
                        }
                    }
                }

                vref.leftseq = self.get_reference_range(
                    reference,
                    (start_position - 20).max(1),
                    start_position - 1,
                );

                let chr_len = instance()
                    .chr_lens
                    .get(region.chr())
                    .copied()
                    .unwrap_or(0) as i64;
                let fallback_len = reference.region_start + reference.ref_seq.len() as i64 - 1;
                let chr_len = if chr_len > 0 { chr_len } else { fallback_len };
                let right_end = (end_position + 20).min(chr_len);
                vref.rightseq = self.get_reference_range(reference, end_position + 1, right_end);

                let mut genotype = format!("{}/{}", genotype1current, genotype2)
                    .replace('&', "")
                    .replace('#', "")
                    .replace('^', "i");

                vref.extra_frequency = round_half_even("0.0000", vref.extra_frequency);
                vref.frequency = round_half_even("0.0000", vref.frequency);
                vref.high_quality_reads_frequency =
                    round_half_even("0.0000", vref.high_quality_reads_frequency);
                vref.msi = round_half_even("0.000", msi);
                vref.msint = msint;
                vref.shift3 = shift3;
                vref.start_position = start_position;
                vref.end_position = end_position;
                vref.refallele = self.validate_refallele(&refallele);
                vref.varallele = varallele;
                vref.genotype = genotype;
                vref.total_pos_coverage = total_pos_coverage;
                vref.ref_forward_count = reference_forward_coverage;
                vref.ref_reverse_count = reference_reverse_coverage;

                let ref_bias = if let Some(ref_var) = &variations_at_pos.reference_variant {
                    ref_var.strand_bias_flag.var_bias
                } else {
                    StrandBiasValue::CantAssess
                };
                vref.strand_bias_flag = StrandBiasFlag::new(ref_bias, vref.strand_bias_flag.var_bias);

                if start_position != position && self.do_pileup {
                    positions_for_changed_ref_variant.push(position);
                }
            }

            if instance().conf.disable_sv {
                variations_at_pos
                    .variants
                    .retain(|vref| !ANY_SV.is_match(&vref.varallele));
            }
        } else if let Some(ref_var) = variations_at_pos.reference_variant.as_mut() {
            self.update_ref_variant(
                position,
                total_pos_coverage,
                ref_var,
                reference,
                reference_forward_coverage,
                reference_reverse_coverage,
                duprate,
            );
        } else {
            variations_at_pos.reference_variant = Some(Variant::new());
        }

        if let Some(ref_var) = variations_at_pos.reference_variant.as_mut() {
            if self.do_pileup
                && (positions_for_changed_ref_variant.contains(&position)
                    || instance().amplicon_based_calling.is_some())
            {
                self.update_ref_variant(
                    position,
                    total_pos_coverage,
                    ref_var,
                    reference,
                    reference_forward_coverage,
                    reference_reverse_coverage,
                    duprate,
                );
            }
        }
    }

    fn update_ref_variant(
        &self,
        position: i64,
        total_pos_coverage: usize,
        vref: &mut Variant,
        reference: &Reference,
        reference_forward_coverage: usize,
        reference_reverse_coverage: usize,
        duprate: f64,
    ) {
        use crate::mods::to_vars_builder::{check_strand_bias, StrandBiasFlag, StrandBiasValue};

        vref.total_pos_coverage = total_pos_coverage;
        vref.position_coverage = 0;
        vref.frequency = 0.0;
        vref.ref_forward_count = reference_forward_coverage;
        vref.ref_reverse_count = reference_reverse_coverage;
        vref.vars_count_on_forward = 0;
        vref.vars_count_on_reverse = 0;
        vref.msi = 0.0;
        vref.msint = 0.0;
        vref.shift3 = 0;
        vref.start_position = position;
        vref.end_position = position;
        vref.high_quality_reads_frequency =
            round_half_even("0.0000", vref.high_quality_reads_frequency);

        let reference_base = reference
            .get(position)
            .map(|b| (b as char).to_string())
            .unwrap_or_default();

        vref.refallele = self.validate_refallele(&reference_base);
        vref.varallele = self.validate_refallele(&reference_base);
        vref.genotype = format!("{}/{}", reference_base, reference_base);
        vref.leftseq.clear();
        vref.rightseq.clear();
        vref.duprate = duprate;
        vref.crispr = 0;
        vref.strand_bias_flag = StrandBiasFlag::new(
            check_strand_bias(reference_forward_coverage, reference_reverse_coverage),
            StrandBiasValue::CantAssess,
        );
    }

    fn validate_refallele(&self, refallele: &str) -> String {
        let mut out = refallele.to_string();
        let replacements = [
            ('M', 'A'),
            ('R', 'A'),
            ('W', 'A'),
            ('S', 'C'),
            ('Y', 'C'),
            ('K', 'G'),
            ('V', 'A'),
            ('H', 'A'),
            ('D', 'A'),
            ('B', 'C'),
        ];

        for (from, to) in replacements {
            if out.contains(from) {
                out = out.replacen(from, &to.to_string(), 1);
            }
        }
        out
    }

    fn proceed_vref_is_deletion(
        &self,
        position: i64,
        del_len: usize,
        reference: &Reference,
    ) -> (f64, i32, f64) {
        let (msi, msint, shift3) = self.detect_microsatellite(reference, position, del_len);
        (msi, shift3, msint)
    }

    fn proceed_vref_is_insertion(
        &self,
        position: i64,
        desc: &str,
        reference: &Reference,
        region: &Region,
    ) -> (f64, i32, f64) {
        let tseq1 = desc.trim_start_matches('+');
        let leftseq = self.get_reference_range(reference, (position - 50).max(1), position);
        let chr_len = instance()
            .chr_lens
            .get(region.chr())
            .copied()
            .unwrap_or(0) as i64;
        let fallback_len = reference.region_start + reference.ref_seq.len() as i64 - 1;
        let chr_len = if chr_len > 0 { chr_len } else { fallback_len };
        let tseq2 = self.get_reference_range(reference, position + 1, (position + 70).min(chr_len));

        let (mut msi, mut msint, shift3) = self.find_msi(tseq1, &tseq2, Some(&leftseq));

        let (tmsi, tmsint, _) = self.find_msi(&leftseq, &tseq2, None);
        if msi < tmsi {
            msi = tmsi;
            msint = tmsint;
        }

        if !tseq1.is_empty() && msi <= (shift3 as f64) / (tseq1.len() as f64) {
            msi = (shift3 as f64) / (tseq1.len() as f64);
        }

        (msi, shift3 as i32, msint)
    }

    /// Build Vars struct from raw variant data at a position
    fn build_vars_at_position(
        &self,
        position: i64,
        var_map: HashMap<VarDesc, RawVariant>,
        ref_coverage: &HashMap<i64, usize>,
        reference: &Reference,
        ref_counts_by_pos: &mut HashMap<i64, (usize, usize)>,
        hicov_by_pos: &HashMap<i64, usize>,
        duprate: f64,
    ) -> Vars {
        use crate::mods::to_vars_builder::{check_strand_bias, StrandBiasFlag, StrandBiasValue};
        
        let mut total_coverage = ref_coverage.get(&position).copied().unwrap_or(0);
        let actual_ref_base = reference.get(position).unwrap_or(b'N');
        let position_hicov = hicov_by_pos.get(&position).copied().unwrap_or(0);

        let has_amp_insertion = var_map.iter().any(|(desc, _)| match desc {
            VarDesc::Ins { seq } => seq.iter().any(|&b| b == b'&'),
            VarDesc::Raw { desc } => desc.starts_with(b"+") && desc.iter().any(|&b| b == b'&'),
            _ => false,
        });

        if has_amp_insertion {
            if let Some(&coverage) = ref_coverage.get(&(position + 1)) {
                total_coverage = coverage;
            }
        }
        
        // First, identify the reference variant and get its forward/reverse counts
        // Reference variant is an SNV where read_base == actual_ref_base
        let mut ref_fwd_count = 0usize;
        let mut ref_rev_count = 0usize;
        for (desc, raw_var) in &var_map {
            if let VarDesc::SNV { ref_base: read_base } = desc {
                if *read_base == actual_ref_base {
                    ref_fwd_count = raw_var.alt_depth_fwd;
                    ref_rev_count = raw_var.alt_depth_rev;
                    break;
                }
            }
        }
        
        // If reference counts were not found in this var_map (e.g., insertion-only map),
        // fall back to the non-insertion reference counts at the same position.
        if ref_fwd_count == 0 && ref_rev_count == 0 {
            if let Some(&(fwd, rev)) = ref_counts_by_pos.get(&position) {
                ref_fwd_count = fwd;
                ref_rev_count = rev;
            }
        }

        // Java logic: if total coverage exceeds ref coverage at this position and
        // there is a reference variant at position+1, use its strand counts.
        let ref_cov_at_pos = ref_coverage.get(&position).copied().unwrap_or(0);
        if total_coverage > ref_cov_at_pos {
            if let Some(&(fwd, rev)) = ref_counts_by_pos.get(&(position + 1)) {
                ref_fwd_count = fwd;
                ref_rev_count = rev;
            }
        }
        
        // Calculate reference strand bias (used as first part of "refBias;varBias" flag)
        let ref_strand_bias = check_strand_bias(ref_fwd_count, ref_rev_count);
        
        let mut variants = Vec::new();
        let mut reference_variant_opt = None;

        for (desc, raw_var) in var_map {
            let total_count = raw_var.alt_depth_fwd + raw_var.alt_depth_rev;
            if total_count == 0 {
                continue;
            }
            let mut ttcov = total_coverage;

            if total_count > total_coverage
                && raw_var.extra_cnt > 0
                && total_count - total_coverage < raw_var.extra_cnt
            {
                ttcov = total_count;
            }

            let is_insertion = matches!(desc, VarDesc::Ins { .. })
                || matches!(&desc, VarDesc::Raw { desc } if desc.starts_with(b"+"))
                || matches!(&desc, VarDesc::Complex { ref_seq, alt_seq } if alt_seq.len() > ref_seq.len());

            if is_insertion && ttcov < total_count {
                ttcov = total_count;
                if let Some(&next_cov) = ref_coverage.get(&(position + 1)) {
                    if next_cov > total_count && ttcov < next_cov - total_count {
                        ttcov = next_cov;
                        if let Some((ref_fwd, ref_rev)) = ref_counts_by_pos.get_mut(&(position + 1)) {
                            *ref_fwd = ref_fwd.saturating_sub(raw_var.alt_depth_fwd);
                            *ref_rev = ref_rev.saturating_sub(raw_var.alt_depth_rev);
                        }
                    }
                }
                total_coverage = ttcov;
            }

            let extra_frequency = if raw_var.extra_cnt > 0 && ttcov > 0 {
                raw_var.extra_cnt as f64 / ttcov as f64
            } else {
                0.0
            };

            let mut variant = self.convert_raw_variant(
                &desc,
                &raw_var,
                position,
                ttcov,
                reference,
                position_hicov,
                extra_frequency,
                duprate,
            );
            
            // For non-reference variants, set the reference forward/reverse counts
            // and update strand bias to include ref bias
            if variant.refallele != variant.varallele {
                variant.ref_forward_count = ref_fwd_count;
                variant.ref_reverse_count = ref_rev_count;
                
                // Update strand bias flag to include reference bias as first part
                // Java: vref.strandBiasFlag = referenceVariant.strandBiasFlag + ";" + vref.strandBiasFlag
                variant.strand_bias_flag = StrandBiasFlag::new(
                    ref_strand_bias,
                    variant.strand_bias_flag.var_bias,
                );
                variants.push(variant);
            } else {
                // This is a reference call - set bias as "ref_bias;0" (Java uses ref bias)
                let ref_bias = check_strand_bias(raw_var.alt_depth_fwd, raw_var.alt_depth_rev);
                variant.strand_bias_flag = StrandBiasFlag::new(
                    ref_bias,
                    StrandBiasValue::CantAssess,
                );
                reference_variant_opt = Some(variant);
            }
        }

        // Sort variants to match Java ordering (meanQuality * variantCount, then descriptionString)
        variants.sort_by(|a, b| {
            let a_count = a.vars_count_on_forward + a.vars_count_on_reverse;
            let b_count = b.vars_count_on_forward + b.vars_count_on_reverse;
            let a_score = a.mean_quality * a_count as f64;
            let b_score = b.mean_quality * b_count as f64;
            match b_score.partial_cmp(&a_score).unwrap_or(std::cmp::Ordering::Equal) {
                std::cmp::Ordering::Equal => a.description_string.cmp(&b.description_string),
                other => other,
            }
        });


        Vars {
            variants,
            reference_variant: reference_variant_opt,
            sv: String::new(),
            sv_flags: Default::default(),
        }
    }

    /// Apply Java ToVarsBuilder genotype selection logic for non-reference variants.
    fn apply_java_genotypes(
        &self,
        variants: &mut [Variant],
        reference_variant: Option<&Variant>,
        position: i64,
        reference: &Reference,
    ) {
        use crate::conf::Configuration;
        use crate::data::patterns::{AMP_ATGC, BEGIN_MINUS_NUMBER, DUP_NUM};

        if variants.is_empty() && reference_variant.is_none() {
            return;
        }

        let conf = &instance().conf;

        let mut genotype1 = if let Some(ref_var) = reference_variant {
            if ref_var.frequency >= conf.freq {
                ref_var.description_string.clone()
            } else if !variants.is_empty() {
                variants[0].description_string.clone()
            } else {
                ref_var.description_string.clone()
            }
        } else if !variants.is_empty() {
            variants[0].description_string.clone()
        } else {
            String::new()
        };

        if genotype1.starts_with('+') {
            if let Some(cap) = DUP_NUM.captures(&genotype1) {
                if let Ok(dup_len) = cap.get(1).map(|m| m.as_str()).unwrap_or("0").parse::<i32>() {
                    genotype1 = format!("+{}", Configuration::SVFLANK + dup_len);
                }
            } else if genotype1.len() > 1 {
                genotype1 = format!("+{}", genotype1.len() - 1);
            }
        }

        for variant in variants.iter_mut() {
            let mut genotype1current = genotype1.clone();
            let mut genotype2 = variant.description_string.clone();

            let mut end_position = position;
            if let Some(caps) = BEGIN_MINUS_NUMBER.captures(&variant.description_string) {
                if let Some(m) = caps.get(1) {
                    if let Ok(del_len) = m.as_str().parse::<i64>() {
                        if del_len > 0 {
                            end_position = position + del_len - 1;
                        }
                    }
                }
            }

            if let Some(caps) = AMP_ATGC.captures(&variant.description_string) {
                if let Some(extra_match) = caps.get(1) {
                    let extra_len = extra_match.as_str().len() as i64;
                    if extra_len > 0 {
                        let extra_seq = self.get_reference_range(
                            reference,
                            end_position + 1,
                            end_position + extra_len,
                        );
                        genotype1current.push_str(&extra_seq);
                        end_position += extra_len;
                    }
                }
            }

            if genotype2.starts_with('+') {
                if genotype2.len() > 1 {
                    genotype2 = format!("+{}", genotype2.len() - 1);
                }
            }

            if let Some(cap) = DUP_NUM.captures(&genotype2) {
                if let Ok(dup_len) = cap.get(1).map(|m| m.as_str()).unwrap_or("0").parse::<i32>() {
                    genotype2 = format!("+{}", 2 * Configuration::SVFLANK + dup_len);
                }
            }

            let genotype = format!("{}/{}", genotype1current, genotype2)
                .replace('&', "")
                .replace('#', "")
                .replace('^', "i");

            variant.genotype = genotype;
        }
    }

    /// Convert RawVariant (from CigarParser) to Variant (for output)
    fn convert_raw_variant(
        &self,
        desc: &VarDesc,
        raw: &RawVariant,
        position: i64,
        total_coverage: usize,
        reference: &Reference,
        position_hicov: usize,
        extra_frequency: f64,
        duprate: f64,
    ) -> Variant {
        if (168600..=168720).contains(&position) {
            event!(
                Level::DEBUG,
                "convert_raw_variant: pos={}, desc={}, desc_type={}",
                position,
                desc.to_key_string(),
                desc.variant_type()
            );
        }
        let mut position = position;
        let total_count = raw.alt_depth_fwd + raw.alt_depth_rev;
        let frequency = if total_coverage > 0 {
            round_half_even("0.0000", total_count as f64 / total_coverage as f64)
        } else {
            0.0
        };

        let mut hicov = position_hicov;
        if hicov < raw.high_qual_read_cnt {
            hicov = raw.high_qual_read_cnt;
        }

        // Look up reference base for this position
        let actual_ref_base = reference.get(position).unwrap_or(b'N');
        
        // Determine variant type and alleles based on VarDesc
        let (mut var_type, mut refallele, mut varallele) = match desc {
            VarDesc::SNV { ref_base: read_base } => {
                // read_base is actually the observed read base (alt)
                // Look up actual reference from Reference struct
                let ref_char = actual_ref_base as char;
                let var_char = *read_base as char;
                (
                    VarType::SNV(var_char),
                    ref_char.to_string(),
                    var_char.to_string(),
                )
            }
            VarDesc::Ins { seq } => {
                let desc_str = format!("+{}", String::from_utf8_lossy(seq));
                let (inferred, ref_str, var_str, start_position) =
                    self.convert_desc_string_to_alleles(&desc_str, position, reference);
                position = start_position;
                (inferred, ref_str, var_str)
            }
            VarDesc::Del { .. } => {
                let desc_str = desc.to_key_string();
                let (inferred, ref_str, var_str, start_position) =
                    self.convert_desc_string_to_alleles(&desc_str, position, reference);
                position = start_position;
                (inferred, ref_str, var_str)
            }
            VarDesc::Complex { alt_seq, .. } => {
                let desc_str = if alt_seq.len() > 1 {
                    let mut s = String::new();
                    s.push(alt_seq[0] as char);
                    s.push('&');
                    s.push_str(&String::from_utf8_lossy(&alt_seq[1..]));
                    s
                } else {
                    String::from_utf8_lossy(alt_seq).to_string()
                };
                let (inferred, ref_str, var_str, start_position) =
                    self.convert_desc_string_to_alleles(&desc_str, position, reference);
                position = start_position;
                (inferred, ref_str, var_str)
            }
            VarDesc::Raw { desc } => {
                let desc_str = String::from_utf8_lossy(desc).to_string();
                let (inferred, ref_str, var_str, start_position) =
                    self.convert_desc_string_to_alleles(&desc_str, position, reference);
                position = start_position;
                (inferred, ref_str, var_str)
            }
        };

        // NOTE: Java VarDict does NOT normalize deletion anchors, so we skip that step
        // to match Java output exactly. Normalization would shift positions and trim shared bases.

        // Calculate end position based on reference allele length
        let end_position = if refallele.len() > 1 {
            position + refallele.len() as i64 - 1
        } else {
            position
        };

        // Get anchor base (base at position-1) for complex insertion genotype formatting
        let anchor_base = reference.get(position - 1).map(|b| b as char);

        // Calculate genotype before moving refallele and varallele
        let genotype = determine_genotype(&refallele, &varallele, frequency, anchor_base);
        
        // Calculate means by dividing sums by count
        // RawVariant stores sums, we need actual means
        let mean_position = if total_count > 0 {
            round_half_even("0.0", raw.mean_pos / total_count as f64)
        } else {
            0.0
        };
        
        let mean_quality = if total_count > 0 {
            round_half_even("0.0", raw.mean_qual / total_count as f64)
        } else {
            0.0
        };
        
        let mean_mapping_quality = if total_count > 0 {
            round_half_even("0.0", raw.mean_mapq / total_count as f64)
        } else {
            0.0
        };
        
        let nm_mean = if total_count > 0 {
            round_half_even("0.0", raw.nm / total_count as f64)
        } else {
            0.0
        };

        // Only populate leftseq/rightseq for non-reference variants
        // Reference calls (ref==alt) get empty strings, which output as "0"
        let is_ref_call = refallele == varallele;
        let (leftseq, rightseq) = if is_ref_call {
            (String::new(), String::new())
        } else {
            (
                self.get_flanking_sequence(reference, position, 20, true),
                self.get_flanking_sequence(reference, end_position, 20, false),
            )
        };

        // Compute MSI for deletions and same-length substitutions only
        // Don't compute for insertions (alt longer than ref) or ref calls
        // Extract deletion length from description like Java does
        let is_insertion = varallele.len() > refallele.len();
        let (msi, msint, shift3) = if is_ref_call || is_insertion {
            (0.0, 0.0, 0)
        } else {
            // Get deletion length from VarDesc, similar to how Java extracts from description string
            match desc {
                VarDesc::Del { len, .. } => {
                    // Pure deletion or deletion with mismatches - use deletion MSI calculation
                    let del_len = *len as usize;
                    if del_len > 0 {
                        self.detect_microsatellite(reference, position, del_len)
                    } else {
                        (0.0, 0.0, 0)
                    }
                }
                VarDesc::Complex { ref_seq, alt_seq } => {
                    // Complex variant - check if it's effectively a deletion (ref > alt)
                    // or a substitution (same length)
                    if ref_seq.len() > alt_seq.len() {
                        // Net deletion: use deletion MSI calculation with ref_seq length
                        self.detect_microsatellite(reference, position, ref_seq.len())
                    } else {
                        // Same length or net insertion: use SNP/MNP MSI calculation
                        // Java: tseq1 = joinRef(ref, position - 30, position + 1)
                        //       tseq2 = joinRef(ref, position + 2, position + 70)
                        self.detect_microsatellite_snp(reference, position)
                    }
                }
                VarDesc::SNV { .. } => {
                    // SNV: use SNP/MNP MSI calculation
                    self.detect_microsatellite_snp(reference, position)
                }
                VarDesc::Raw { desc } => {
                    if desc.starts_with(b"-") {
                        let tail = String::from_utf8_lossy(&desc[1..]);
                        if let Some(del_len) = parse_leading_digits(&tail) {
                            if del_len > 0 {
                                self.detect_microsatellite(reference, position, del_len)
                            } else {
                                (0.0, 0.0, 0)
                            }
                        } else {
                            (0.0, 0.0, 0)
                        }
                    } else {
                        self.detect_microsatellite_snp(reference, position)
                    }
                }
                _ => (0.0, 0.0, 0),
            }
        };

        Variant {
            description_string: desc.to_key_string(),
            refallele,
            varallele,
            vartype: var_type,
            start_position: position,
            end_position,
            vars_count_on_forward: raw.alt_depth_fwd,
            vars_count_on_reverse: raw.alt_depth_rev,
            position_coverage: total_count,
            total_pos_coverage: total_coverage,
            frequency,
            high_quality_reads_frequency: if total_coverage > 0 {
                if hicov > 0 {
                    round_half_even("0.0000", raw.high_qual_read_cnt as f64 / hicov as f64)
                } else {
                    0.0
                }
            } else {
                0.0
            },
            extra_frequency: if extra_frequency > 0.0 {
                round_half_even("0.0000", extra_frequency)
            } else {
                0.0
            },
            mean_position,
            mean_quality,
            mean_mapping_quality,
            strand_bias_flag: {
                use crate::mods::to_vars_builder::{check_strand_bias, StrandBiasFlag, StrandBiasValue};
                // Calculate var bias from variant counts
                let var_bias = check_strand_bias(raw.alt_depth_fwd, raw.alt_depth_rev);
                // Reference bias will be set later when we have ref counts
                StrandBiasFlag::new(StrandBiasValue::CantAssess, var_bias)
            },
            is_at_least_at_2_positions: raw.pstd,
            has_at_least_2_diff_qualities: raw.qstd,
            leftseq,
            rightseq,
            msi,
            msint,
            shift3,
            nm: nm_mean,
            high_qual_read_cnt: raw.high_qual_read_cnt,
            low_qual_read_cnt: raw.low_qual_read_cnt,
            hicov,
            ref_forward_count: 0,  // Will be set by build_vars_at_position for non-ref variants
            ref_reverse_count: 0,
            genotype,
            duprate,
            crispr: 0,
        }
    }

    /// Convert Java-style description string to alleles and adjusted start position
    fn convert_desc_string_to_alleles(
        &self,
        desc_str: &str,
        position: i64,
        reference: &Reference,
    ) -> (VarType, String, String, i64) {
        let mut start_position = position;
        let mut end_position = position;
        let mut ref_str = String::new();
        let mut var_str = String::new();

        if desc_str.starts_with('+') {
            // Insertion
            let ins_str = desc_str.trim_start_matches('+');
            ref_str = reference
                .get(position)
                .map(|b| (b as char).to_string())
                .unwrap_or_default();
            var_str = format!("{}{}", ref_str, ins_str);
            let mut had_amp = false;

            if let Some(extra) = extract_amp_seq(desc_str) {
                had_amp = true;
                var_str = var_str.replacen("&", "", 1);
                let extra_len = extra.len() as i64;
                ref_str.push_str(&self.get_reference_range(
                    reference,
                    end_position + 1,
                    end_position + extra_len,
                ));
                end_position += extra_len;

                while let Some(vextra) = extract_amp_seq(&var_str) {
                    var_str = var_str.replacen("&", "", 1);
                    let vextra_len = vextra.len() as i64;
                    ref_str.push_str(&self.get_reference_range(
                        reference,
                        end_position + 1,
                        end_position + vextra_len,
                    ));
                    end_position += vextra_len;
                }

                if !ref_str.is_empty() && !var_str.is_empty() {
                    ref_str = ref_str[1..].to_string();
                    var_str = var_str[1..].to_string();
                    start_position += 1;
                }
            }

            if let Some((matched_seq, tail)) = extract_hash_caret(desc_str) {
                let matched_len = matched_seq.len() as i64;
                end_position += matched_len;
                ref_str.push_str(&self.get_reference_range(
                    reference,
                    end_position - matched_len + 1,
                    end_position,
                ));

                if let Some(del_len) = parse_leading_digits(&tail) {
                    let del_len_i64 = del_len as i64;
                    ref_str.push_str(&self.get_reference_range(
                        reference,
                        end_position + 1,
                        end_position + del_len_i64,
                    ));
                    end_position += del_len_i64;
                }

                var_str = var_str.replacen('#', "", 1);
                var_str = remove_caret_and_digits(&var_str);
            }

            if desc_str.contains('^') {
                var_str = var_str.replacen('^', "", 1);
            }

            let inferred = infer_var_type_from_alleles(&ref_str, &var_str);
            return (inferred, ref_str, var_str, start_position);
        }

        if desc_str.starts_with('-') {
            let del_len = parse_leading_digits(&desc_str[1..]).unwrap_or(0);
            if del_len > 0 {
                end_position = position + del_len as i64 - 1;
            }

            // Remove leading -<digits>
            if del_len > 0 {
                let prefix = format!("-{}", del_len);
                var_str = desc_str.replacen(&prefix, "", 1);
            } else {
                var_str = desc_str.to_string();
            }

            let has_suffix = desc_str.contains('&') || desc_str.contains('#') || desc_str.contains('^');
            if !has_suffix {
                // Simple deletion: anchor base at position-1
                let anchor_pos = position - 1;
                if anchor_pos >= 1 {
                    if let Some(anchor_base) = reference.get(anchor_pos) {
                        ref_str.push(anchor_base as char);
                    }
                }
                for i in 0..del_len as i64 {
                    if let Some(b) = reference.get(position + i) {
                        ref_str.push(b as char);
                    }
                }
                var_str = if !ref_str.is_empty() {
                    ref_str[0..1].to_string()
                } else {
                    String::new()
                };
                start_position = anchor_pos;
            } else {
                // Deletion with suffix: build ref from deleted bases
                for i in 0..del_len as i64 {
                    if let Some(b) = reference.get(position + i) {
                        ref_str.push(b as char);
                    }
                }
            }

            if let Some(extra) = extract_amp_seq(desc_str) {
                var_str = var_str.replacen("&", "", 1);
                let extra_len = extra.len() as i64;
                ref_str.push_str(&self.get_reference_range(
                    reference,
                    end_position + 1,
                    end_position + extra_len,
                ));
                end_position += extra_len;

                while let Some(vextra) = extract_amp_seq(&var_str) {
                    var_str = var_str.replacen("&", "", 1);
                    let vextra_len = vextra.len() as i64;
                    ref_str.push_str(&self.get_reference_range(
                        reference,
                        end_position + 1,
                        end_position + vextra_len,
                    ));
                    end_position += vextra_len;
                }
            }

            if let Some((matched_seq, tail)) = extract_hash_caret(desc_str) {
                let matched_len = matched_seq.len() as i64;
                end_position += matched_len;
                ref_str.push_str(&self.get_reference_range(
                    reference,
                    end_position - matched_len + 1,
                    end_position,
                ));

                if let Some(extra_del_len) = parse_leading_digits(&tail) {
                    let del_len_i64 = extra_del_len as i64;
                    ref_str.push_str(&self.get_reference_range(
                        reference,
                        end_position + 1,
                        end_position + del_len_i64,
                    ));
                    end_position += del_len_i64;
                }

                var_str = var_str.replacen('#', "", 1);
                var_str = remove_caret_and_digits(&var_str);
            }

            if desc_str.contains('^') {
                var_str = var_str.replacen('^', "", 1);
            }

            return (VarType::Deletion(del_len), ref_str, var_str, start_position);
        }

        // SNP/MNP or other substitution
        ref_str = reference
            .get(position)
            .map(|b| (b as char).to_string())
            .unwrap_or_default();
        var_str = desc_str.to_string();

        if let Some(extra) = extract_amp_seq(desc_str) {
            var_str = var_str.replacen("&", "", 1);
            let extra_len = extra.len() as i64;
            ref_str.push_str(&self.get_reference_range(
                reference,
                end_position + 1,
                end_position + extra_len,
            ));
            end_position += extra_len;

            while let Some(vextra) = extract_amp_seq(&var_str) {
                var_str = var_str.replacen("&", "", 1);
                let vextra_len = vextra.len() as i64;
                ref_str.push_str(&self.get_reference_range(
                    reference,
                    end_position + 1,
                    end_position + vextra_len,
                ));
                end_position += vextra_len;
            }
        }

        if let Some((matched_seq, tail)) = extract_hash_caret(desc_str) {
            let matched_len = matched_seq.len() as i64;
            end_position += matched_len;
            ref_str.push_str(&self.get_reference_range(
                reference,
                end_position - matched_len + 1,
                end_position,
            ));

            if let Some(extra_del_len) = parse_leading_digits(&tail) {
                let del_len_i64 = extra_del_len as i64;
                ref_str.push_str(&self.get_reference_range(
                    reference,
                    end_position + 1,
                    end_position + del_len_i64,
                ));
                end_position += del_len_i64;
            }

            var_str = var_str.replacen('#', "", 1);
            var_str = remove_caret_and_digits(&var_str);
        }

        if desc_str.contains('^') {
            var_str = var_str.replacen('^', "", 1);
        }

        let inferred = infer_var_type_from_alleles(&ref_str, &var_str);
        (inferred, ref_str, var_str, start_position)
    }
    
    /// Get 20bp flanking sequence from reference
    /// For left: get bases before the position
    /// For right: get bases after the position
    fn get_flanking_sequence(&self, reference: &Reference, position: i64, length: usize, is_left: bool) -> String {
        let mut seq = String::new();
        
        if is_left {
            // Get bases BEFORE the position (position-20 to position-1)
            for i in (1..=length as i64).rev() {
                let pos = position - i;
                if let Some(base) = reference.get(pos) {
                    seq.push(base as char);
                }
            }
        } else {
            // Get bases AFTER the position (position+1 to position+20)
            for i in 1..=length as i64 {
                let pos = position + i;
                if let Some(base) = reference.get(pos) {
                    seq.push(base as char);
                }
            }
        }
        
        seq
    }

    /// Detect microsatellite instability for a variant
    /// Java-compatible implementation matching findMSI algorithm
    /// For deletions: tseq1 = deleted bases, tseq2 = sequence after
    /// Returns (msi, msint, shift3) where:
    /// - msi: number of repeats (instability score)
    /// - msint: unit length (1 for homopolymer, 2 for dinucleotide, etc.)
    fn detect_microsatellite(&self, reference: &Reference, position: i64, del_len: usize) -> (f64, f64, i32) {
        // Get left sequence (70 bases before position)
        let leftseq = self.get_reference_range(reference, position - 70, position - 1);
        
        // Get tseq (from position to position + deletion_len + 70)
        // For deletions, tseq1 is the deleted portion, tseq2 is what follows
        let tseq = self.get_reference_range(reference, position, position + (del_len as i64) - 1 + 70);
        
        if tseq.len() < del_len {
            return (0.0, 0.0, 0);
        }
        
        let tseq1 = &tseq[..del_len.min(tseq.len())];
        let tseq2 = if del_len < tseq.len() { &tseq[del_len..] } else { "" };
        
        // First call: findMSI(tseq1, tseq2, leftseq)
        let (mut msi, msint1, shift3) = self.find_msi(tseq1, tseq2, Some(&leftseq));
        let mut msint = msint1;
        
        // Second call: findMSI(leftseq, tseq2, None) 
        let (tmsi, tmsint, _) = self.find_msi(&leftseq, tseq2, None);
        
        // Java: if (msi < tmsi) { msi = tmsi; msint = tmsint; }
        if msi < tmsi {
            msi = tmsi;
            msint = tmsint;
            // Don't change shift3 - Java keeps the original shift3
        }
        
        // Java: if (msi <= shift3 / (double) dellen) { msi = shift3 / (double) dellen; }
        if del_len > 0 && msi <= (shift3 as f64) / (del_len as f64) {
            msi = (shift3 as f64) / (del_len as f64);
        }
        
        (msi, msint, shift3 as i32)
    }
    
    /// Detect microsatellite instability for SNP/MNP variants
    /// Java-compatible implementation for variants that don't start with + or -
    /// Java: tseq1 = joinRef(ref, position - 30, position + 1)
    ///       tseq2 = joinRef(ref, position + 2, position + 70)
    fn detect_microsatellite_snp(&self, reference: &Reference, position: i64) -> (f64, f64, i32) {
        // tseq1 = reference from (position - 30) to (position + 1)
        let tseq1 = self.get_reference_range(reference, (position - 30).max(1), position + 1);
        
        // tseq2 = reference from (position + 2) to (position + 70)
        let tseq2 = self.get_reference_range(reference, position + 2, position + 70);
        
        // Call findMSI with no left sequence
        let (msi, msint, shift3) = self.find_msi(&tseq1, &tseq2, None);
        
        (msi, msint, shift3 as i32)
    }
    
    /// Get a range of bases from reference as a string
    fn get_reference_range(&self, reference: &Reference, start: i64, end: i64) -> String {
        let mut seq = String::new();
        for pos in start..=end {
            if let Some(base) = reference.get(pos) {
                seq.push(base as char);
            }
        }
        seq
    }

    /// Java-compatible findMSI implementation
    /// Finds microsatellite repeats by checking:
    /// 1. Repeats at the END of tseq1 (optionally with left prepended)
    /// 2. Repeats at the START of tseq2
    /// Returns (msi, msint, shift3) where shift3 is count of matching bases at start of tseq
    fn find_msi(&self, tseq1: &str, tseq2: &str, left: Option<&str>) -> (f64, f64, usize) {
        if tseq1.is_empty() {
            return (0.0, 0.0, 0);
        }
        
        let mut max_msi = 0.0;
        let mut max_msint = 0.0;
        
        // Try unit lengths from 1 to 6
        for nmsi in 1..=6.min(tseq1.len()) {
            // Get the last nmsi bases of tseq1 as the repeat unit
            let msint = &tseq1[tseq1.len() - nmsi..];
            
            // Count repeats at end of tseq1 (and optionally left+tseq1)
            let search_str = if let Some(l) = left {
                format!("{}{}", l, tseq1)
            } else {
                tseq1.to_string()
            };
            
            let end_repeats = self.count_trailing_repeats(&search_str, msint);
            
            // Count repeats at start of tseq2
            let start_repeats = self.count_leading_repeats(tseq2, msint);
            
            let cur_msi = end_repeats + start_repeats;
            
            if cur_msi > max_msi {
                max_msi = cur_msi;
                max_msint = nmsi as f64;
            }
        }
        
        // Calculate shift3: count how many chars at start of tseq match tseq2
        // Java: while (shift3 < tseq2.length() && tseq.charAt(shift3) == tseq2.charAt(shift3)) { shift3++; }
        let tseq = format!("{}{}", tseq1, tseq2);
        let tseq_bytes = tseq.as_bytes();
        let tseq2_bytes = tseq2.as_bytes();
        let mut shift3 = 0usize;
        while shift3 < tseq2_bytes.len() && shift3 < tseq_bytes.len() && tseq_bytes[shift3] == tseq2_bytes[shift3] {
            shift3 += 1;
        }
        
        (max_msi, max_msint, shift3)
    }
    
    /// Count how many times a unit repeats at the end of a string
    fn count_trailing_repeats(&self, s: &str, unit: &str) -> f64 {
        if unit.is_empty() || s.len() < unit.len() {
            return 0.0;
        }
        
        let mut count = 0;
        let mut pos = s.len();
        
        while pos >= unit.len() {
            let start = pos - unit.len();
            if &s[start..pos] == unit {
                count += 1;
                pos = start;
            } else {
                break;
            }
        }
        
        count as f64
    }
    
    /// Count how many times a unit repeats at the start of a string
    fn count_leading_repeats(&self, s: &str, unit: &str) -> f64 {
        if unit.is_empty() || s.len() < unit.len() {
            return 0.0;
        }
        
        let mut count = 0;
        let mut pos = 0;
        
        while pos + unit.len() <= s.len() {
            if &s[pos..pos + unit.len()] == unit {
                count += 1;
                pos += unit.len();
            } else {
                break;
            }
        }
        
        count as f64
    }

            pub fn run_somatic_post_processor(
                &self,
                normal_data: AlignedVarsData,
                tumor_data: AlignedVarsData,
                region: &Region,
                splice: &HashSet<String>,
            ) -> Vec<String> {
                self.run_somatic_post_processor_with_combine_lookup(
                    normal_data,
                    tumor_data,
                    region,
                    splice,
                    0,
                    None,
                )
            }

            pub fn run_somatic_post_processor_with_combine_lookup(
                &self,
                normal_data: AlignedVarsData,
                tumor_data: AlignedVarsData,
                region: &Region,
                splice: &HashSet<String>,
                initial_max_read_length: usize,
                combine_lookup: Option<&SomaticCombineLookup>,
            ) -> Vec<String> {
                const STRONG_SOMATIC: &str = "StrongSomatic";
                const SAMPLE_SPECIFIC: &str = "SampleSpecific";
                const DELETION: &str = "Deletion";

                let output_region = OutputRegion {
                    chr: region.chr().to_string(),
                    start: region.start() as i64,
                    end: region.end() as i64,
                    gene: region.gene().to_string(),
                };

                let mut output_lines = Vec::new();
                let mut max_read_length = initial_max_read_length;

                let mut all_positions: std::collections::BTreeSet<i64> =
                    tumor_data.aligned_variants.keys().copied().collect();
                all_positions.extend(normal_data.aligned_variants.keys().copied());

                for position in all_positions {
                    if position < region.start() as i64 || position > region.end() as i64 {
                        continue;
                    }

                    let tumor_vars = tumor_data.aligned_variants.get(&position);
                    let normal_vars = normal_data.aligned_variants.get(&position);

                    match (tumor_vars, normal_vars) {
                        (None, None) => {}
                        (None, Some(variants)) => {
                            self.calling_for_one_sample(
                                variants,
                                true,
                                DELETION,
                                &output_region,
                                splice,
                                &mut output_lines,
                            );
                        }
                        (Some(variants), None) => {
                            self.calling_for_one_sample(
                                variants,
                                false,
                                SAMPLE_SPECIFIC,
                                &output_region,
                                splice,
                                &mut output_lines,
                            );
                        }
                        (Some(tumor), Some(normal)) => {
                            self.calling_for_both_samples(
                                position,
                                tumor,
                                normal,
                                &output_region,
                                splice,
                                STRONG_SOMATIC,
                                &mut max_read_length,
                                combine_lookup,
                                &mut output_lines,
                            );
                        }
                    }
                }

                output_lines
            }

            fn calling_for_one_sample(
                &self,
                variants: &Vars,
                is_first_cover: bool,
                var_label: &str,
                region: &OutputRegion,
                splice: &HashSet<String>,
                output_lines: &mut Vec<String>,
            ) {
                if variants.variants.is_empty() {
                    return;
                }

                for variant in &variants.variants {
                    if variant.refallele == variant.varallele {
                        continue;
                    }

                    if !self.is_good_var(variant, variants.reference_variant.as_ref(), splice) {
                        continue;
                    }

                    let mut variant = variant.clone();
                    if var_type_string(&variant.refallele, &variant.varallele) == "Complex" {
                        variant.adj_complex();
                    }

                    let output = if is_first_cover {
                        SomaticOutputVariant::from_variants(
                            Some(&variant),
                            Some(&variant),
                            None,
                            Some(&variant),
                            region,
                            "",
                            &variants.sv,
                            var_label,
                            &self.sample_name,
                        )
                    } else {
                        SomaticOutputVariant::from_variants(
                            Some(&variant),
                            Some(&variant),
                            Some(&variant),
                            None,
                            region,
                            &variants.sv,
                            "",
                            var_label,
                            &self.sample_name,
                        )
                    };
                    output_lines.push(output.to_string());
                }
            }

            fn calling_for_both_samples(
                &self,
                position: i64,
                tumor_vars: &Vars,
                normal_vars: &Vars,
                region: &OutputRegion,
                splice: &HashSet<String>,
                strong_somatic_label: &str,
                max_read_length: &mut usize,
                combine_lookup: Option<&SomaticCombineLookup>,
                output_lines: &mut Vec<String>,
            ) {
                if tumor_vars.variants.is_empty() && normal_vars.variants.is_empty() {
                    return;
                }

                if !tumor_vars.variants.is_empty() {
                    self.print_variations_from_first_sample(
                        position,
                        tumor_vars,
                        normal_vars,
                        region,
                        splice,
                        strong_somatic_label,
                        max_read_length,
                        combine_lookup,
                        output_lines,
                    );
                } else if !normal_vars.variants.is_empty() {
                    self.print_variations_from_second_sample(
                        position,
                        tumor_vars,
                        normal_vars,
                        region,
                        splice,
                        max_read_length,
                        combine_lookup,
                        output_lines,
                    );
                }
            }

            fn print_variations_from_first_sample(
                &self,
                position: i64,
                tumor_vars: &Vars,
                normal_vars: &Vars,
                region: &OutputRegion,
                splice: &HashSet<String>,
                strong_somatic_label: &str,
                max_read_length: &mut usize,
                combine_lookup: Option<&SomaticCombineLookup>,
                output_lines: &mut Vec<String>,
            ) {
                const LIKELY_LOH: &str = "LikelyLOH";
                const GERMLINE: &str = "Germline";
                const STRONG_LOH: &str = "StrongLOH";
                const FALSE_VALUE: &str = "FALSE";
                let debug_pos = env::var("VARDICT_DEBUG_POS")
                    .ok()
                    .and_then(|value| value.parse::<i64>().ok());
                let trace_this_pos = debug_pos == Some(position);

                let mut number_of_processed_variation = 0usize;
                while number_of_processed_variation < tumor_vars.variants.len()
                    && self.is_good_var(
                        &tumor_vars.variants[number_of_processed_variation],
                        tumor_vars.reference_variant.as_ref(),
                        splice,
                    )
                {
                    let mut tumor_variant = tumor_vars.variants[number_of_processed_variation].clone();
                    if tumor_variant.refallele == tumor_variant.varallele {
                        number_of_processed_variation += 1;
                        continue;
                    }

                    let description_string = tumor_variant.description_string.clone();
                    if var_type_string(&tumor_variant.refallele, &tumor_variant.varallele) == "Complex" {
                        tumor_variant.adj_complex();
                    }

                    if let Some(mut normal_variant) =
                        Self::find_variant_by_description(normal_vars, &description_string).cloned()
                    {
                        let var_label = self.determinate_somatic_type(
                            normal_vars,
                            &tumor_variant,
                            &mut normal_variant,
                            splice,
                        );
                        if trace_this_pos {
                            event!(
                                Level::DEBUG,
                                position,
                                description = %description_string,
                                branch = "first_sample_direct_match",
                                var_label = %var_label,
                                tumor_pcov = tumor_variant.position_coverage,
                                tumor_tot = tumor_variant.total_pos_coverage,
                                tumor_fwd = tumor_variant.vars_count_on_forward,
                                tumor_rev = tumor_variant.vars_count_on_reverse,
                                normal_pcov = normal_variant.position_coverage,
                                normal_tot = normal_variant.total_pos_coverage,
                                normal_fwd = normal_variant.vars_count_on_forward,
                                normal_rev = normal_variant.vars_count_on_reverse,
                                "somatic_output_selection"
                            );
                        }
                        let output = SomaticOutputVariant::from_variants(
                            Some(&tumor_variant),
                            Some(&normal_variant),
                            Some(&tumor_variant),
                            Some(&normal_variant),
                            region,
                            &tumor_vars.sv,
                            &normal_vars.sv,
                            &var_label,
                            &self.sample_name,
                        );
                        output_lines.push(output.to_string());
                    } else {
                        let mut normal_variant_for_combine = Variant::default();
                        normal_variant_for_combine.description_string = description_string.clone();

                        let mut var_label = strong_somatic_label.to_string();
                        if Self::should_run_combine_analysis(&tumor_variant) {
                            if let Some(lookup) = combine_lookup {
                                let combine_type = self.combine_analysis_with_lookup(
                                    &tumor_variant,
                                    &mut normal_variant_for_combine,
                                    region.chr.as_str(),
                                    position,
                                    &description_string,
                                    splice,
                                    max_read_length,
                                    lookup,
                                );
                                if combine_type == FALSE_VALUE {
                                    number_of_processed_variation += 1;
                                    continue;
                                }
                                if !combine_type.is_empty() {
                                    var_label = combine_type;
                                }
                            }
                        }

                        let normal_variant_for_print =
                            if let Some(first_normal_variant) = normal_vars.variants.first() {
                                let mut placeholder = Variant::default();
                                placeholder.total_pos_coverage = first_normal_variant.total_pos_coverage;
                                placeholder.ref_forward_count = first_normal_variant.ref_forward_count;
                                placeholder.ref_reverse_count = first_normal_variant.ref_reverse_count;
                                Some(placeholder)
                            } else {
                                normal_vars.reference_variant.clone()
                            };

                        let output = if var_label == strong_somatic_label {
                            SomaticOutputVariant::from_variants(
                                Some(&tumor_variant),
                                Some(&tumor_variant),
                                Some(&tumor_variant),
                                normal_variant_for_print.as_ref(),
                                region,
                                &tumor_vars.sv,
                                &normal_vars.sv,
                                strong_somatic_label,
                                &self.sample_name,
                            )
                        } else {
                            if trace_this_pos {
                                event!(
                                    Level::DEBUG,
                                    position,
                                    description = %description_string,
                                    branch = "first_sample_combine_lookup",
                                    var_label = %var_label,
                                    tumor_pcov = tumor_variant.position_coverage,
                                    tumor_tot = tumor_variant.total_pos_coverage,
                                    tumor_fwd = tumor_variant.vars_count_on_forward,
                                    tumor_rev = tumor_variant.vars_count_on_reverse,
                                    normal_pcov = normal_variant_for_combine.position_coverage,
                                    normal_tot = normal_variant_for_combine.total_pos_coverage,
                                    normal_fwd = normal_variant_for_combine.vars_count_on_forward,
                                    normal_rev = normal_variant_for_combine.vars_count_on_reverse,
                                    "somatic_output_selection"
                                );
                            }
                            SomaticOutputVariant::from_variants(
                                Some(&tumor_variant),
                                Some(&tumor_variant),
                                Some(&tumor_variant),
                                Some(&normal_variant_for_combine),
                                region,
                                &tumor_vars.sv,
                                &normal_vars.sv,
                                &var_label,
                                &self.sample_name,
                            )
                        };
                        output_lines.push(output.to_string());
                    }

                    number_of_processed_variation += 1;
                }

                if number_of_processed_variation == 0 {
                    if normal_vars.variants.is_empty() {
                        return;
                    }

                    for normal_variant in &normal_vars.variants {
                        if !self.is_good_var(normal_variant, normal_vars.reference_variant.as_ref(), splice) {
                            continue;
                        }

                        let mut normal_variant = normal_variant.clone();
                        let description_string = normal_variant.description_string.clone();

                        if let Some(mut tumor_variant) =
                            Self::find_variant_by_description(tumor_vars, &description_string).cloned()
                        {
                            if tumor_variant.refallele == tumor_variant.varallele {
                                continue;
                            }

                            let var_label = if tumor_variant.frequency < instance().conf.lofreq {
                                LIKELY_LOH
                            } else {
                                GERMLINE
                            };

                            if var_type_string(&normal_variant.refallele, &normal_variant.varallele) == "Complex" {
                                tumor_variant.adj_complex();
                            }

                            if trace_this_pos {
                                event!(
                                    Level::DEBUG,
                                    position,
                                    description = %description_string,
                                    branch = "second_sample_match_when_first_not_processed",
                                    var_label,
                                    tumor_pcov = tumor_variant.position_coverage,
                                    tumor_tot = tumor_variant.total_pos_coverage,
                                    tumor_fwd = tumor_variant.vars_count_on_forward,
                                    tumor_rev = tumor_variant.vars_count_on_reverse,
                                    normal_pcov = normal_variant.position_coverage,
                                    normal_tot = normal_variant.total_pos_coverage,
                                    normal_fwd = normal_variant.vars_count_on_forward,
                                    normal_rev = normal_variant.vars_count_on_reverse,
                                    "somatic_output_selection"
                                );
                            }

                            let output = SomaticOutputVariant::from_variants(
                                Some(&tumor_variant),
                                Some(&normal_variant),
                                Some(&tumor_variant),
                                Some(&normal_variant),
                                region,
                                &tumor_vars.sv,
                                &normal_vars.sv,
                                var_label,
                                &self.sample_name,
                            );
                            output_lines.push(output.to_string());
                        } else {
                            if normal_variant.refallele == normal_variant.varallele {
                                continue;
                            }

                            let first_tumor_variant = tumor_vars.variants.first();
                            let total_coverage = first_tumor_variant
                                .map(|variant| variant.total_pos_coverage)
                                .unwrap_or(0);

                            let tumor_reference_variant = tumor_vars.reference_variant.as_ref();
                            let reference_forward = tumor_reference_variant
                                .map(|variant| variant.vars_count_on_forward)
                                .unwrap_or(0);
                            let reference_reverse = tumor_reference_variant
                                .map(|variant| variant.vars_count_on_reverse)
                                .unwrap_or(0);

                            let genotype = if let Some(variant) = first_tumor_variant {
                                variant.genotype.clone()
                            } else if let Some(reference_variant) = tumor_reference_variant {
                                format!(
                                    "{}/{}",
                                    reference_variant.description_string,
                                    reference_variant.description_string
                                )
                            } else {
                                "N/N".to_string()
                            };

                            if var_type_string(&normal_variant.refallele, &normal_variant.varallele) == "Complex" {
                                normal_variant.adj_complex();
                            }

                            let mut tumor_variant_for_print = Variant::default();
                            tumor_variant_for_print.total_pos_coverage = total_coverage;
                            tumor_variant_for_print.ref_forward_count = reference_forward;
                            tumor_variant_for_print.ref_reverse_count = reference_reverse;
                            tumor_variant_for_print.genotype = genotype;

                            let output = SomaticOutputVariant::from_variants(
                                Some(&normal_variant),
                                Some(&normal_variant),
                                Some(&tumor_variant_for_print),
                                Some(&normal_variant),
                                region,
                                "",
                                &normal_vars.sv,
                                STRONG_LOH,
                                &self.sample_name,
                            );
                            output_lines.push(output.to_string());
                        }
                    }
                }
            }

            fn print_variations_from_second_sample(
                &self,
                position: i64,
                tumor_vars: &Vars,
                normal_vars: &Vars,
                region: &OutputRegion,
                splice: &HashSet<String>,
                max_read_length: &mut usize,
                combine_lookup: Option<&SomaticCombineLookup>,
                output_lines: &mut Vec<String>,
            ) {
                const STRONG_LOH: &str = "StrongLOH";
                const FALSE_VALUE: &str = "FALSE";

                for normal_variant in &normal_vars.variants {
                    if normal_variant.refallele == normal_variant.varallele {
                        continue;
                    }

                    if !self.is_good_var(normal_variant, normal_vars.reference_variant.as_ref(), splice) {
                        continue;
                    }

                    let mut normal_variant = normal_variant.clone();
                    let description_string = normal_variant.description_string.clone();
                    let mut tumor_variant_for_combine =
                        Self::find_variant_by_description(tumor_vars, &description_string)
                            .cloned()
                            .unwrap_or_else(|| {
                                let mut variant = Variant::default();
                                variant.description_string = description_string;
                                variant
                            });
                    tumor_variant_for_combine.position_coverage = 0;

                    let mut var_label = STRONG_LOH.to_string();
                    let mut combine_type = String::new();
                    if Self::should_run_combine_analysis(&normal_variant) {
                        if let Some(lookup) = combine_lookup {
                            combine_type = self.combine_analysis_with_lookup(
                                &normal_variant,
                                &mut tumor_variant_for_combine,
                                region.chr.as_str(),
                                position,
                                &normal_variant.description_string,
                                splice,
                                max_read_length,
                                lookup,
                            );
                            if combine_type == FALSE_VALUE {
                                continue;
                            }
                        }
                    }

                    let tumor_variant_for_print = if !combine_type.is_empty() {
                        var_label = combine_type;
                        Some(tumor_variant_for_combine)
                    } else {
                        tumor_vars.reference_variant.clone()
                    };

                    if var_type_string(&normal_variant.refallele, &normal_variant.varallele) == "Complex" {
                        normal_variant.adj_complex();
                    }

                    let output = SomaticOutputVariant::from_variants(
                        Some(&normal_variant),
                        Some(&normal_variant),
                        tumor_variant_for_print.as_ref(),
                        Some(&normal_variant),
                        region,
                        "",
                        &normal_vars.sv,
                        &var_label,
                        &self.sample_name,
                    );
                    output_lines.push(output.to_string());
                }
            }

            fn should_run_combine_analysis(variant: &Variant) -> bool {
                let var_type = var_type_string(&variant.refallele, &variant.varallele);
                if var_type == "SNV" {
                    return false;
                }

                let description = variant.description_string.as_str();
                let has_minus_num_num = Self::contains_minus_num_num(description);

                (description.len() > 10 || has_minus_num_num)
                    && variant.position_coverage < instance().conf.minr + 3
                    && !description.contains('<')
            }

            fn contains_minus_num_num(description: &str) -> bool {
                let bytes = description.as_bytes();
                if bytes.len() < 3 {
                    return false;
                }

                for idx in 0..=bytes.len() - 3 {
                    if bytes[idx] == b'-'
                        && bytes[idx + 1].is_ascii_digit()
                        && bytes[idx + 2].is_ascii_digit()
                    {
                        return true;
                    }
                }

                false
            }

            fn combine_analysis_with_lookup(
                &self,
                variant1: &Variant,
                variant2: &mut Variant,
                chr_name: &str,
                position: i64,
                description_string: &str,
                _splice: &HashSet<String>,
                max_read_length: &mut usize,
                combine_lookup: &SomaticCombineLookup,
            ) -> String {
                const FALSE_VALUE: &str = "FALSE";
                const GERMLINE: &str = "Germline";

                if variant1.end_position.saturating_sub(variant1.start_position) as usize
                    > instance().conf.sv_min_len
                {
                    return String::new();
                }

                let lookup_result =
                    combine_lookup(chr_name, position, description_string, *max_read_length);
                if lookup_result.max_read_length > 0 {
                    *max_read_length = lookup_result.max_read_length;
                }

                let Some(vref) = lookup_result.combined_variant else {
                    return FALSE_VALUE.to_string();
                };

                if vref.position_coverage.saturating_sub(variant1.position_coverage)
                    >= instance().conf.minr
                {
                    variant2.total_pos_coverage =
                        vref.total_pos_coverage.saturating_sub(variant1.total_pos_coverage);
                    variant2.position_coverage =
                        vref.position_coverage.saturating_sub(variant1.position_coverage);
                    variant2.ref_forward_count =
                        vref.ref_forward_count.saturating_sub(variant1.ref_forward_count);
                    variant2.ref_reverse_count =
                        vref.ref_reverse_count.saturating_sub(variant1.ref_reverse_count);
                    variant2.vars_count_on_forward =
                        vref.vars_count_on_forward.saturating_sub(variant1.vars_count_on_forward);
                    variant2.vars_count_on_reverse =
                        vref.vars_count_on_reverse.saturating_sub(variant1.vars_count_on_reverse);

                    if variant2.position_coverage != 0 {
                        let vref_cov = vref.position_coverage as f64;
                        let variant1_cov = variant1.position_coverage as f64;
                        let variant2_cov = variant2.position_coverage as f64;

                        variant2.mean_position =
                            (vref.mean_position * vref_cov - variant1.mean_position * variant1_cov)
                                / variant2_cov;
                        variant2.mean_quality =
                            (vref.mean_quality * vref_cov - variant1.mean_quality * variant1_cov)
                                / variant2_cov;
                        variant2.mean_mapping_quality =
                            (vref.mean_mapping_quality * vref_cov
                                - variant1.mean_mapping_quality * variant1_cov)
                                / variant2_cov;
                        variant2.high_quality_reads_frequency =
                            (vref.high_quality_reads_frequency * vref_cov
                                - variant1.high_quality_reads_frequency * variant1_cov)
                                / variant2_cov;
                        variant2.extra_frequency =
                            (vref.extra_frequency * vref_cov - variant1.extra_frequency * variant1_cov)
                                / variant2_cov;
                        variant2.nm =
                            (vref.nm * vref_cov - variant1.nm * variant1_cov) / variant2_cov;
                    } else {
                        variant2.mean_position = 0.0;
                        variant2.mean_quality = 0.0;
                        variant2.mean_mapping_quality = 0.0;
                        variant2.high_quality_reads_frequency = 0.0;
                        variant2.extra_frequency = 0.0;
                        variant2.nm = 0.0;
                    }

                    variant2.is_at_least_at_2_positions = true;
                    variant2.has_at_least_2_diff_qualities = true;

                    if variant2.total_pos_coverage == 0 {
                        return FALSE_VALUE.to_string();
                    }

                    variant2.frequency =
                        variant2.position_coverage as f64 / variant2.total_pos_coverage as f64;
                    variant2.high_qual_read_cnt = variant1.high_qual_read_cnt;
                    variant2.low_qual_read_cnt = variant1.low_qual_read_cnt;
                    variant2.genotype = vref.genotype.clone();
                    variant2.description_string = description_string.to_string();

                    variant2.strand_bias_flag = StrandBiasFlag::new(
                        check_strand_bias(variant2.ref_forward_count, variant2.ref_reverse_count),
                        check_strand_bias(
                            variant2.vars_count_on_forward,
                            variant2.vars_count_on_reverse,
                        ),
                    );

                    return GERMLINE.to_string();
                }

                if vref.position_coverage + 2 < variant1.position_coverage {
                    return FALSE_VALUE.to_string();
                }

                String::new()
            }

            fn determinate_somatic_type(
                &self,
                variants: &Vars,
                standard_variant: &Variant,
                variant_to_compare: &mut Variant,
                splice: &HashSet<String>,
            ) -> String {
                const STRONG_SOMATIC: &str = "StrongSomatic";
                const LIKELY_LOH: &str = "LikelyLOH";
                const GERMLINE: &str = "Germline";
                const LIKELY_SOMATIC: &str = "LikelySomatic";
                const AF_DIFF: &str = "AFDiff";

                let standard_var_type = var_type_string(&standard_variant.refallele, &standard_variant.varallele);

                let mut var_label = if self.is_good_var_with_type(
                    variant_to_compare,
                    variants.reference_variant.as_ref(),
                    splice,
                    Some(&standard_var_type),
                ) {
                    if standard_variant.frequency > (1.0 - instance().conf.lofreq)
                        && variant_to_compare.frequency < 0.8
                        && variant_to_compare.frequency > 0.2
                    {
                        LIKELY_LOH.to_string()
                    } else if variant_to_compare.frequency < instance().conf.lofreq
                        || variant_to_compare.position_coverage <= 1
                    {
                        LIKELY_SOMATIC.to_string()
                    } else {
                        GERMLINE.to_string()
                    }
                } else if variant_to_compare.frequency < instance().conf.lofreq
                    || variant_to_compare.position_coverage <= 1
                {
                    LIKELY_SOMATIC.to_string()
                } else {
                    AF_DIFF.to_string()
                };

                if self.is_noise(variant_to_compare) && standard_var_type == "SNV" {
                    var_label = STRONG_SOMATIC.to_string();
                }

                var_label
            }

            fn is_noise(&self, variant: &mut Variant) -> bool {
                let mean_quality = variant.mean_quality;
                let quality_is_low = (mean_quality < 4.5
                    || (mean_quality < 12.0 && !variant.has_at_least_2_diff_qualities))
                    && variant.position_coverage <= 3;
                let low_freq_with_low_quality = mean_quality < instance().conf.goodq
                    && variant.frequency < 2.0 * instance().conf.lofreq
                    && variant.position_coverage <= 1;

                if quality_is_low || low_freq_with_low_quality {
                    variant.total_pos_coverage = variant
                        .total_pos_coverage
                        .saturating_sub(variant.position_coverage);
                    variant.position_coverage = 0;
                    variant.vars_count_on_forward = 0;
                    variant.vars_count_on_reverse = 0;
                    variant.frequency = 0.0;
                    variant.high_quality_reads_frequency = 0.0;
                    return true;
                }

                false
            }

            fn find_variant_by_description<'a>(vars: &'a Vars, description: &str) -> Option<&'a Variant> {
                vars.variants
                    .iter()
                    .find(|variant| variant.description_string == description)
                    .or_else(|| {
                        vars.reference_variant
                            .as_ref()
                            .filter(|variant| variant.description_string == description)
                    })
            }

    /// Step 4: Run SimplePostProcessor to filter and format output
    fn run_simple_post_processor(
        &self,
        data: AlignedVarsData,
        region: &Region,
        splice: &HashSet<String>,
    ) -> Result<Vec<String>> {
        let mut output_lines = Vec::new();

        // Output uses the same coordinate base as the parsed regions
        let output_region = OutputRegion {
            chr: region.chr().to_string(),
            start: region.start() as i64,
            end: region.end() as i64,
            gene: region.gene().to_string(),
        };

        let mut aligned_order_index: HashMap<i64, usize> = HashMap::new();
        for (idx, pos) in data.aligned_variants_order.iter().enumerate() {
            aligned_order_index.insert(*pos, idx);
        }

        let ordered_positions = java_hashmap_iteration_order(
            data.aligned_variants.keys().copied(),
            data.aligned_variants.len(),
            Some(&aligned_order_index),
        );

        for position in ordered_positions {
            let Some(vars) = data.aligned_variants.get(&position) else {
                continue;
            };
            event!(Level::DEBUG, "[PostProcessor] Processing position {}: {} variants", position, vars.variants.len());

            // Skip positions outside region only when SV marker is absent (Java parity)
            if (position < region.start() as i64 || position > region.end() as i64)
                && vars.sv.is_empty()
            {
                event!(Level::DEBUG, "[PostProcessor] Skipping position {} - outside region {}-{}", position, region.start(), region.end());
                continue;
            }

            // Skip empty variants unless pileup mode
            if vars.variants.is_empty() {
                event!(Level::DEBUG, "[PostProcessor] Position {} has 0 variants", position);
                if !self.do_pileup {
                    continue;
                }
                // In pileup mode, output reference (or empty if none)
                if let Some(ref ref_var) = vars.reference_variant {
                    let output = SimpleOutputVariant::from_variant(
                        ref_var,
                        &output_region,
                        &self.sample_name,
                        &vars.sv,
                    );
                    output_lines.push(output.to_string());
                } else {
                    let output = SimpleOutputVariant::empty_with_sv(
                        position,
                        &output_region,
                        &self.sample_name,
                        &vars.sv,
                    );
                    output_lines.push(output.to_string());
                }
                continue;
            }

            for variant in &vars.variants {
                event!(Level::DEBUG, "[PostProcessor] Variant: pos={} ref={} alt={} freq={:.3} good={} type={:?} hicnt={} meanpos={:.1} meanq={:.1} fwd={} rev={}", 
                    variant.start_position, variant.refallele, variant.varallele, variant.frequency,
                    self.is_good_var(variant, vars.reference_variant.as_ref(), splice), variant.vartype,
                    variant.high_qual_read_cnt, variant.mean_position, variant.mean_quality,
                    variant.vars_count_on_forward, variant.vars_count_on_reverse);

                // Skip if ref contains N
                if variant.refallele.contains('N') {
                    event!(Level::DEBUG, "[PostProcessor] Skipping - ref contains N");
                    continue;
                }

                // Skip reference calls unless pileup mode
                if variant.refallele == variant.varallele {
                    event!(Level::DEBUG, "[PostProcessor] Skipping - ref call (ref==alt)");
                    if !self.do_pileup {
                        continue;
                    }
                }

                // If variant start position shifted (pileup + single variant), output reference
                if variant.start_position != position && self.do_pileup && vars.variants.len() == 1 {
                    if let Some(ref ref_var) = vars.reference_variant {
                        let output = SimpleOutputVariant::from_variant(
                            ref_var,
                            &output_region,
                            &self.sample_name,
                            &vars.sv,
                        );
                        output_lines.push(output.to_string());
                    } else {
                        let output = SimpleOutputVariant::empty_with_sv(
                            position,
                            &output_region,
                            &self.sample_name,
                            &vars.sv,
                        );
                        output_lines.push(output.to_string());
                    }
                }

                let var_type = var_type_string(&variant.refallele, &variant.varallele);

                // Apply quality filter (isGoodVar equivalent)
                if !self.is_good_var(variant, vars.reference_variant.as_ref(), splice) {
                    event!(Level::DEBUG, "[PostProcessor] Skipping - failed isGoodVar filter");
                    if !self.do_pileup {
                        continue;
                    }
                }

                let mut variant = variant.clone();
                if var_type == "Complex" {
                    variant.adj_complex();
                }

                event!(Level::DEBUG, "[PostProcessor] Adding variant to output");

                // Generate output
                let output = SimpleOutputVariant::from_variant(
                    &variant,
                    &output_region,
                    &self.sample_name,
                    &vars.sv,
                );
                output_lines.push(output.to_string());
            }
            }

        Ok(output_lines)
    }

    /// Quality filter - equivalent to Java Variant.isGoodVar()
    /// 
    /// Java checks: frequency >= conf.freq, hicnt >= conf.minr, 
    /// meanPosition >= conf.readPosFilter, meanQuality >= conf.goodq,
    /// highQualityToLowQualityRatio >= conf.qratio
    fn is_good_var(
        &self,
        variant: &Variant,
        ref_variant: Option<&Variant>,
        splice: &HashSet<String>,
    ) -> bool {
        self.is_good_var_with_type(variant, ref_variant, splice, None)
    }

    fn is_good_var_with_type(
        &self,
        variant: &Variant,
        ref_variant: Option<&Variant>,
        splice: &HashSet<String>,
        forced_type: Option<&str>,
    ) -> bool {
        if variant.refallele.is_empty() {
            return false;
        }

        let var_type = forced_type
            .map(|value| value.to_string())
            .unwrap_or_else(|| var_type_string(&variant.refallele, &variant.varallele));

        if variant.frequency < instance().conf.freq
            || variant.high_qual_read_cnt < instance().conf.minr
            || variant.mean_position < instance().conf.read_pos_filter
            || variant.mean_quality < instance().conf.goodq
        {
            event!(Level::DEBUG, 
                "[is_good_var] FAILED: freq={} (need>={}), hicnt={} (need>={}), meanpos={} (need>={}), meanq={} (need>={})",
                variant.frequency, instance().conf.freq,
                variant.high_qual_read_cnt, instance().conf.minr,
                variant.mean_position, instance().conf.read_pos_filter,
                variant.mean_quality, instance().conf.goodq
            );
            return false;
        }

        if let Some(ref_var) = ref_variant {
            if ref_var.high_qual_read_cnt > instance().conf.minr && variant.frequency < 0.25 {
                let d = variant.mean_mapping_quality
                    + variant.refallele.len() as f64
                    + variant.varallele.len() as f64;
                let f = (1.0 + d) / (ref_var.mean_mapping_quality + 1.0);
                if (d - 2.0 < 5.0 && ref_var.mean_mapping_quality > 20.0) || f < 0.25 {
                    return false;
                }
            }
        }

        if var_type == "Deletion" {
            let splice_key = format!("{}-{}", variant.start_position, variant.end_position);
            if splice.contains(&splice_key) {
                return false;
            }
        }

        let hl_ratio = if variant.low_qual_read_cnt > 0 {
            variant.high_qual_read_cnt as f64 / variant.low_qual_read_cnt as f64
        } else if variant.high_qual_read_cnt > 0 {
            variant.high_qual_read_cnt as f64 * 2.0
        } else {
            0.0
        };

        if hl_ratio < instance().conf.qratio {
            return false;
        }

        if variant.frequency > 0.30 {
            return true;
        }

        if variant.mean_mapping_quality < instance().conf.mapq as f64 {
            return false;
        }

        if variant.msi >= 15.0
            && variant.frequency <= instance().conf.monomer_msi_frequency
            && variant.msint == 1.0
        {
            return false;
        }

        if variant.msi >= 12.0
            && variant.frequency <= instance().conf.non_monomer_msi_frequency
            && variant.msint > 1.0
        {
            return false;
        }

        if variant.strand_bias_flag.is_ref_good_var_biased() && variant.frequency < 0.20 {
            if var_type.is_empty()
                || var_type == "SNV"
                || (variant.refallele.len() < 3 && variant.varallele.len() < 3)
            {
                return false;
            }
        }

        true
    }
}

fn build_amplicon_debug_prefix(vars_at_amplicon: &Vars) -> String {
    let mut entries: Vec<Variant> = Vec::new();
    if let Some(reference_variant) = vars_at_amplicon.reference_variant.as_ref() {
        entries.push(reference_variant.clone());
    }
    entries.extend(vars_at_amplicon.variants.iter().cloned());
    entries.sort_by(|left, right| {
        let left_is_insertion = left.description_string.starts_with('+');
        let right_is_insertion = right.description_string.starts_with('+');
        match left_is_insertion.cmp(&right_is_insertion) {
            std::cmp::Ordering::Equal => left.description_string.cmp(&right.description_string),
            other => other,
        }
    });

    entries
        .iter()
        .map(format_amplicon_debug_prefix_variant)
        .collect::<Vec<_>>()
        .join(" & ")
}

fn format_amplicon_debug_prefix_variant(variant: &Variant) -> String {
    let mut key = variant.description_string.clone();
    if key.starts_with('+') {
        key = format!("I{}", key);
    }
    let pstd = if variant.is_at_least_at_2_positions { 1 } else { 0 };
    let qstd = if variant.has_at_least_2_diff_qualities { 1 } else { 0 };

    format!(
        "{}:{}:F-{}:R-{}:{}:{}:{}:{}:{}:{}:{}:{}:{}",
        key,
        variant.vars_count_on_forward + variant.vars_count_on_reverse,
        variant.vars_count_on_forward,
        variant.vars_count_on_reverse,
        format!("{:.4}", variant.frequency),
        variant.strand_bias_flag.var_bias.as_int(),
        format!("{:.1}", variant.mean_position),
        pstd,
        format!("{:.1}", variant.mean_quality),
        qstd,
        format!("{:.4}", variant.high_quality_reads_frequency),
        format!("{:.1}", variant.mean_mapping_quality),
        format!("{:.3}", amplicon_debug_qratio(variant)),
    )
}

fn amplicon_debug_qratio(variant: &Variant) -> f64 {
    if variant.low_qual_read_cnt > 0 {
        variant.high_qual_read_cnt as f64 / variant.low_qual_read_cnt as f64
    } else if variant.high_qual_read_cnt > 0 {
        variant.high_qual_read_cnt as f64 * 2.0
    } else {
        0.0
    }
}

fn extract_amp_seq(s: &str) -> Option<String> {
    let idx = s.find('&')?;
    let tail = &s[idx + 1..];
    let seq: String = tail
        .chars()
        .take_while(|c| matches!(c, 'A' | 'T' | 'G' | 'C'))
        .collect();
    if seq.is_empty() {
        None
    } else {
        Some(seq)
    }
}

fn extract_hash_caret(s: &str) -> Option<(String, String)> {
    let hash_idx = s.find('#')?;
    let tail = &s[hash_idx + 1..];
    let caret_idx = tail.find('^')?;
    let matched = &tail[..caret_idx];
    let tail_after = &tail[caret_idx + 1..];
    if matched.is_empty() || tail_after.is_empty() {
        None
    } else {
        Some((matched.to_string(), tail_after.to_string()))
    }
}

fn parse_leading_digits(s: &str) -> Option<usize> {
    let digits: String = s.chars().take_while(|c| c.is_ascii_digit()).collect();
    if digits.is_empty() {
        None
    } else {
        digits.parse::<usize>().ok()
    }
}

fn remove_caret_and_digits(s: &str) -> String {
    if let Some(idx) = s.find('^') {
        let bytes = s.as_bytes();
        let mut end = idx + 1;
        while end < bytes.len() && bytes[end].is_ascii_digit() {
            end += 1;
        }
        let mut out = String::new();
        out.push_str(&s[..idx]);
        out.push_str(&s[end..]);
        out
    } else {
        s.to_string()
    }
}

fn infer_var_type_from_alleles(refallele: &str, varallele: &str) -> VarType {
    if refallele.len() == 1 && varallele.len() == 1 {
        return VarType::SNV(varallele.chars().next().unwrap_or('N'));
    }

    if varallele.len() > refallele.len() && varallele.starts_with(refallele) {
        return VarType::Insertion(varallele[refallele.len()..].to_string());
    }

    if refallele.len() > varallele.len() && refallele.starts_with(varallele) {
        return VarType::Deletion(refallele.len().saturating_sub(varallele.len()));
    }

    VarType::Complex {
        insertion: varallele.to_string(),
        deletion: refallele.len(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ensure_test_scope_initialized() {
        use crate::scopedata::global_read_only_scope::{GlobalReadOnlyScope, INSTANCE};
        let _ = INSTANCE.get_or_init(GlobalReadOnlyScope::default);
    }

    fn make_amplicon_region(gene: &str, insert_start: usize, insert_end: usize) -> Region {
        Region::new_with_insert(
            "chr1".to_string(),
            100,
            110,
            gene.to_string(),
            insert_start,
            insert_end,
        )
    }

    fn make_test_variant(
        ref_allele: &str,
        var_allele: &str,
        description: &str,
        start_position: i64,
        end_position: i64,
        total_coverage: usize,
        variant_coverage: usize,
        frequency: f64,
    ) -> Variant {
        let mut variant = Variant::default();
        variant.description_string = description.to_string();
        variant.refallele = ref_allele.to_string();
        variant.varallele = var_allele.to_string();
        variant.start_position = start_position;
        variant.end_position = end_position;
        variant.total_pos_coverage = total_coverage;
        variant.position_coverage = variant_coverage;
        variant.vars_count_on_forward = variant_coverage;
        variant.frequency = frequency;
        variant.high_qual_read_cnt = 200;
        variant.low_qual_read_cnt = 1;
        variant.mean_position = 100.0;
        variant.mean_quality = 100.0;
        variant.mean_mapping_quality = 100.0;
        variant.is_at_least_at_2_positions = true;
        variant.has_at_least_2_diff_qualities = true;
        variant.hicov = total_coverage;
        variant.high_quality_reads_frequency = frequency;
        variant.genotype = format!("{}/{}", ref_allele, var_allele);
        variant
    }

    fn make_vars_at_position(
        position: i64,
        variants: Vec<Variant>,
        reference_variant: Option<Variant>,
    ) -> HashMap<i64, Vars> {
        let mut vars = HashMap::new();
        vars.insert(
            position,
            Vars {
                variants,
                reference_variant,
                ..Vars::default()
            },
        );
        vars
    }

    fn parse_amplicon_output_fields(line: &str) -> Vec<String> {
        let fields: Vec<String> = line.split('\t').map(str::to_string).collect();
        assert_eq!(fields.len(), 38);
        fields
    }

    fn parse_somatic_output_fields(line: &str) -> Vec<String> {
        let fields: Vec<String> = line.split('\t').map(str::to_string).collect();
        assert_eq!(fields.len(), 55);
        fields
    }

    #[test]
    fn test_pipeline_creation() {
        let pipeline = VarDictPipeline::new("test_sample")
            .with_min_frequency(0.05)
            .with_min_base_quality(20.0)
            .with_min_mapping_quality(10);

        assert_eq!(pipeline.sample_name, "test_sample");
        assert!((pipeline.min_frequency - 0.05).abs() < 0.001);
        assert!((pipeline.min_base_quality - 20.0).abs() < 0.001);
        assert_eq!(pipeline.min_mapping_quality, 10);
    }

    #[test]
    fn test_create_variant_java_values() {
        use crate::mods::to_vars_builder::StrandBiasValue;
        use crate::prelude::SmallVecBytes;
        use crate::variants::variants::Variant as RawVariant;

        let pipeline = VarDictPipeline::new("test");
        let position = 1_234_567i64;

        let mut raw = RawVariant::default();
        raw.alt_depth = 4;
        raw.alt_depth_fwd = 3;
        raw.alt_depth_rev = 5;
        raw.mean_pos = 9.0;
        raw.mean_qual = 10.5;
        raw.mean_mapq = 31.0;
        raw.nm = 8.0;
        raw.high_qual_read_cnt = 44;
        raw.low_qual_read_cnt = 35;

        let desc = VarDesc::Raw {
            desc: SmallVecBytes::from_slice(b"T"),
        };
        let mut vars_at_pos: HashMap<VarDesc, RawVariant> = HashMap::new();
        vars_at_pos.insert(desc.clone(), raw);

        let mut var_list = Vec::new();
        let mut debug_lines = Vec::new();
        let keys = vec![desc];

        let _ = pipeline.create_variant_records(
            position,
            &vars_at_pos,
            10,
            &mut var_list,
            &mut debug_lines,
            &keys,
            0,
            0.0,
        );

        assert_eq!(var_list.len(), 1);
        let v = &var_list[0];
        assert_eq!(v.description_string, "T");
        assert_eq!(v.position_coverage, 4);
        assert_eq!(v.vars_count_on_forward, 3);
        assert_eq!(v.vars_count_on_reverse, 5);
        assert_eq!(v.strand_bias_flag.var_bias, StrandBiasValue::NoBias);
        assert!((v.frequency - 0.4).abs() < 0.0001);
        assert!((v.mean_position - 2.2).abs() < 0.001);
        assert!((v.mean_quality - 2.6).abs() < 0.001);
        assert!((v.mean_mapping_quality - 7.8).abs() < 0.001);
        assert!((v.nm - 2.0).abs() < 0.001);
        assert_eq!(v.high_qual_read_cnt, 44);
        assert_eq!(v.low_qual_read_cnt, 35);
        assert_eq!(v.hicov, 0);
    }

    #[test]
    fn test_create_insertion_java_values() {
        use crate::mods::to_vars_builder::StrandBiasValue;
        use crate::prelude::SmallVecBytes;
        use crate::variants::variants::Variant as RawVariant;

        let pipeline = VarDictPipeline::new("test");
        let position = 1_234_567i64;

        let mut raw = RawVariant::default();
        raw.alt_depth = 4;
        raw.alt_depth_fwd = 3;
        raw.alt_depth_rev = 5;
        raw.mean_pos = 9.0;
        raw.mean_qual = 10.5;
        raw.mean_mapq = 31.0;
        raw.nm = 8.0;
        raw.high_qual_read_cnt = 44;
        raw.low_qual_read_cnt = 35;

        let desc = VarDesc::Raw {
            desc: SmallVecBytes::from_slice(b"T"),
        };
        let mut insertion_variations: HashMap<VarDesc, RawVariant> = HashMap::new();
        insertion_variations.insert(desc.clone(), raw);
        let mut insertion_vars: HashMap<i64, HashMap<VarDesc, RawVariant>> = HashMap::new();
        insertion_vars.insert(position, insertion_variations);

        let mut non_insertion_vars: HashMap<i64, HashMap<VarDesc, RawVariant>> = HashMap::new();
        let ref_coverage: HashMap<i64, usize> = HashMap::new();
        let reference = Reference::from_seq_with_start(b"A", 1);

        let mut var_list = Vec::new();
        let mut debug_lines = Vec::new();

        let updated_tcov = pipeline.create_insertion_records(
            position,
            10,
            insertion_vars.get(&position),
            &mut non_insertion_vars,
            &ref_coverage,
            &reference,
            &mut var_list,
            &mut debug_lines,
            0,
            0.0,
        );

        assert_eq!(updated_tcov, 10);
        assert_eq!(var_list.len(), 1);
        let v = &var_list[0];
        assert_eq!(v.description_string, "T");
        assert_eq!(v.position_coverage, 4);
        assert_eq!(v.vars_count_on_forward, 3);
        assert_eq!(v.vars_count_on_reverse, 5);
        assert_eq!(v.strand_bias_flag.var_bias, StrandBiasValue::NoBias);
        assert!((v.frequency - 0.4).abs() < 0.0001);
        assert!((v.mean_position - 2.2).abs() < 0.001);
        assert!((v.mean_quality - 2.6).abs() < 0.001);
        assert!((v.mean_mapping_quality - 7.8).abs() < 0.001);
        assert!((v.nm - 2.0).abs() < 0.001);
        assert_eq!(v.high_qual_read_cnt, 44);
        assert_eq!(v.low_qual_read_cnt, 35);
        assert_eq!(v.hicov, 44);
        assert!((v.high_quality_reads_frequency - 1.0).abs() < 0.001);
        assert!((v.extra_frequency - 0.0).abs() < 0.0001);
    }

    #[test]
    fn test_strand_bias_using_check_strand_bias() {
        use crate::mods::to_vars_builder::{check_strand_bias, StrandBiasValue};

        // Low count (<=12) - both strands have reads = NoBias (2)
        assert_eq!(check_strand_bias(5, 5), StrandBiasValue::NoBias);

        // Low count - only one strand = CantAssess (0)
        assert_eq!(check_strand_bias(10, 0), StrandBiasValue::CantAssess);

        // High count - balanced = NoBias (2)
        assert_eq!(check_strand_bias(50, 50), StrandBiasValue::NoBias);

        // High count - imbalanced (95:5, 5/100=5%, 95/100=95%, but 5/100=0.05 fails >= 5% threshold)
        // Wait, 5/100 = 0.05 = 5%, so exactly at threshold. Let's check 4 vs 96:
        assert_eq!(check_strand_bias(96, 4), StrandBiasValue::HasBias);  // 4/100 = 4% < 5%
    }

    #[test]
    fn test_genotype_determination() {
        // SNV: single base substitution uses "ref/alt" format (matching Java)
        assert_eq!(determine_genotype("C", "T", 0.9, None), "C/T");
        assert_eq!(determine_genotype("C", "T", 0.5, None), "C/T");
        assert_eq!(determine_genotype("C", "T", 0.3, None), "C/T");  // Frequency is not used for SNVs
        assert_eq!(determine_genotype("G", "G", 1.0, None), "G/G");  // Ref call
    }

    #[test]
    fn test_is_good_var_strand_bias() {
        use crate::mods::to_vars_builder::{StrandBiasFlag, StrandBiasValue};
        use std::collections::HashSet;
        
        let pipeline = VarDictPipeline::new("test");
        
        // Good variant - no bias on both ref and var (2;2)
        let mut good_var = Variant::default();
        good_var.strand_bias_flag = StrandBiasFlag::new(StrandBiasValue::NoBias, StrandBiasValue::NoBias);
        good_var.is_at_least_at_2_positions = true;
        good_var.has_at_least_2_diff_qualities = true;
        good_var.frequency = 0.3;
        good_var.high_qual_read_cnt = 5;
        good_var.mean_position = 10.0;
        good_var.mean_quality = 30.0;
        good_var.refallele = "A".to_string();
        good_var.varallele = "G".to_string();
        assert!(pipeline.is_good_var(&good_var, None, &HashSet::new()));

        // Bad variant - ref good (2), var has bias (1), low frequency (2;1 pattern)
        let mut bad_var = Variant::default();
        bad_var.strand_bias_flag = StrandBiasFlag::new(StrandBiasValue::NoBias, StrandBiasValue::HasBias);
        bad_var.is_at_least_at_2_positions = true;
        bad_var.has_at_least_2_diff_qualities = true;
        bad_var.frequency = 0.05; // Low frequency + "2;1" pattern = bad
        bad_var.high_qual_read_cnt = 5;
        bad_var.mean_position = 10.0;
        bad_var.mean_quality = 30.0;
        bad_var.refallele = "A".to_string();
        bad_var.varallele = "G".to_string();
        assert!(!pipeline.is_good_var(&bad_var, None, &HashSet::new()));

        // High frequency can overcome strand bias (2;1 but freq > 0.20)
        let mut high_freq_bias = Variant::default();
        high_freq_bias.strand_bias_flag = StrandBiasFlag::new(StrandBiasValue::NoBias, StrandBiasValue::HasBias);
        high_freq_bias.is_at_least_at_2_positions = true;
        high_freq_bias.has_at_least_2_diff_qualities = true;
        high_freq_bias.frequency = 0.5; // High enough to pass despite 2;1 pattern
        high_freq_bias.high_qual_read_cnt = 5;
        high_freq_bias.mean_position = 10.0;
        high_freq_bias.mean_quality = 30.0;
        high_freq_bias.refallele = "A".to_string();
        high_freq_bias.varallele = "G".to_string();
        assert!(pipeline.is_good_var(&high_freq_bias, None, &HashSet::new()));
    }

    #[test]
    fn test_record_preprocessor_unmapped_with_alignment_passes() {
        use crate::data::bam_reader::BamReader;

        let bam_path = "/home/eck/workspace/vardict_rs/test_data/test_168714.bam";
        let mut bam_reader = BamReader::open(bam_path).expect("Failed to open test BAM");

        bam_reader
            .fetch("20", 168600, 168800)
            .expect("Failed to fetch test region");

        let sam_filter = 0x504u32;
        let min_mapq = 0u8;

        let mut record = rust_htslib::bam::Record::new();
        let mut unmapped_with_alignment = 0usize;
        let mut mapped = 0usize;
        let mut filtered = 0usize;

        while bam_reader.read(&mut record).unwrap_or(false) {
            let qname = std::str::from_utf8(record.qname()).unwrap_or("");
            if !qname.contains("SRR098401.7003120") {
                continue;
            }

            if sam_filter != 0 {
                if (record.flags() & (sam_filter as u16)) != 0 {
                    filtered += 1;
                    continue;
                }
            }

            if min_mapq > 0 && record.mapq() < min_mapq {
                filtered += 1;
                continue;
            }

            const SECONDARY_ALIGNMENT: u16 = 0x100;
            if (record.flags() & SECONDARY_ALIGNMENT) != 0 && sam_filter != 0 {
                filtered += 1;
                continue;
            }

            let seq = record.seq();
            if seq.len() == 0 || (seq.len() == 1 && seq.as_bytes()[0] == b'*') {
                filtered += 1;
                continue;
            }

            if record.is_unmapped() && !record.cigar().is_empty() && record.pos() >= 0 {
                unmapped_with_alignment += 1;
            } else if !record.is_unmapped() {
                mapped += 1;
            }
        }

        assert_eq!(filtered, 1, "Unmapped read should be filtered by samfilter");
        assert_eq!(unmapped_with_alignment, 0, "Unmapped read should not be counted");
        assert_eq!(mapped, 1, "Expected 1 properly mapped read");
    }

    #[test]
    fn test_record_preprocessor_dump_all_bed_regions() {
        use std::fs;
        use std::io::{BufRead, BufReader, Write};
        use crate::data::bam_reader::BamReader;

        let bam_path = "/home/eck/workspace/vardict_rs/VarDictJava/tests/integration/input/NA12878.chrom20.ILLUMINA.bwa.CEU.exome.20121211.bam";
        let bed_path = "/home/eck/workspace/vardict_rs/VarDictJava/tests/integration/input/20120518.consensus.annotation.bed.chr20";

        let file = fs::File::open(bed_path).expect("Failed to open BED file");
        let mut reader = BufReader::new(file);
        let mut line = String::new();

        let mut bam_reader = BamReader::open(bam_path).expect("Failed to open BAM");
        let pipeline = VarDictPipeline::new("test");
        let sam_filter = 0x504u32;

        let out_dir = "/home/eck/workspace/vardict_rs/tmp_compare";
        let out_path = "/home/eck/workspace/vardict_rs/tmp_compare/rust.preproc.all.txt";
        fs::create_dir_all(out_dir).expect("Failed to create tmp_compare");
        let out_file = fs::File::create(out_path).expect("Failed to create rust preproc dump");
        let mut writer = std::io::BufWriter::new(out_file);
        let mut total_lines = 0usize;

        while reader.read_line(&mut line).unwrap_or(0) > 0 {
            let trimmed = line.trim();
            if trimmed.is_empty()
                || trimmed.starts_with('#')
                || trimmed.starts_with("track")
                || trimmed.starts_with("browser")
            {
                line.clear();
                continue;
            }

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

            let region = Region::new(chr.clone(), start, end, gene.clone());
            writeln!(writer, "REGION\t{}\t{}\t{}\t{}", chr, start, end, gene)
                .expect("Failed to write region line");
            total_lines += 1;

            let (_records, lines) = pipeline
                .collect_filtered_records(&region, &mut bam_reader, sam_filter)
                .expect("Failed to collect filtered reads");
            for line in lines {
                writeln!(writer, "{}", line).expect("Failed to write read line");
                total_lines += 1;
            }

            line.clear();
        }

        assert!(total_lines > 0, "Expected at least one kept read across all bed regions");
    }

    #[test]
    fn test_amplicon_nocov_strict_lt_maxcov_over_50() {
        use std::collections::HashSet;

        ensure_test_scope_initialized();

        let pipeline = VarDictPipeline::new("sample");
        let position = 105i64;
        let amplicon_regions = vec![
            make_amplicon_region("G1", 105, 105),
            make_amplicon_region("G2", 105, 105),
            make_amplicon_region("G3", 105, 105),
        ];
        let group_region = amplicon_regions.last().unwrap().clone();

        let good_variant = make_test_variant("A", "T", "A>T", 105, 105, 100, 20, 1.0);
        let ref_cov_1 = make_test_variant("A", "A", "A", 105, 105, 1, 1, 0.0);
        let ref_cov_2 = make_test_variant("A", "A", "A", 105, 105, 2, 2, 0.0);

        let vars_per_amplicon = vec![
            make_vars_at_position(position, vec![good_variant], None),
            make_vars_at_position(position, Vec::new(), Some(ref_cov_1)),
            make_vars_at_position(position, Vec::new(), Some(ref_cov_2)),
        ];

        let output_lines = pipeline.run_amplicon_post_processor(
            &group_region,
            &vars_per_amplicon,
            &amplicon_regions,
            &HashSet::new(),
        );

        assert_eq!(output_lines.len(), 1);
        let fields = parse_amplicon_output_fields(&output_lines[0]);
        assert_eq!(fields[6], "T");
        assert_eq!(fields[34], "1");
        assert_eq!(fields[35], "3");
        assert_eq!(fields[36], "1");
        assert_eq!(fields[37], "0");
    }

    #[test]
    fn test_amplicon_primer_overlap_decrements_gvscnt_without_ampbias() {
        use std::collections::HashSet;

        ensure_test_scope_initialized();

        let pipeline = VarDictPipeline::new("sample");
        let position = 105i64;
        let amplicon_regions = vec![
            make_amplicon_region("G1", 100, 110),
            make_amplicon_region("G2", 100, 110),
            make_amplicon_region("G3", 105, 105),
        ];
        let group_region = amplicon_regions.last().unwrap().clone();

        let deletion_1 = make_test_variant("ATC", "A", "delX", 104, 106, 40, 18, 1.0);
        let deletion_2 = make_test_variant("ATC", "A", "delX", 104, 106, 38, 16, 1.0);

        let vars_per_amplicon = vec![
            make_vars_at_position(position, vec![deletion_1], None),
            make_vars_at_position(position, vec![deletion_2], None),
            HashMap::new(),
        ];

        let output_lines = pipeline.run_amplicon_post_processor(
            &group_region,
            &vars_per_amplicon,
            &amplicon_regions,
            &HashSet::new(),
        );

        assert_eq!(output_lines.len(), 1);
        let fields = parse_amplicon_output_fields(&output_lines[0]);
        assert_eq!(fields[33], "Deletion");
        assert_eq!(fields[34], "1");
        assert_eq!(fields[35], "1");
        assert_eq!(fields[37], "0");
    }

    #[test]
    fn test_amplicon_flag_resets_when_current_gvscnt_decreases() {
        use std::collections::HashSet;

        ensure_test_scope_initialized();

        let pipeline = VarDictPipeline::new("sample");
        let position = 105i64;
        let amplicon_regions = vec![
            make_amplicon_region("G1", 100, 110),
            make_amplicon_region("G2", 100, 110),
            make_amplicon_region("G3", 105, 105),
        ];
        let group_region = amplicon_regions.last().unwrap().clone();

        let deletion_1 = make_test_variant("ATC", "A", "delX", 104, 106, 40, 18, 1.0);
        let deletion_2 = make_test_variant("ATC", "A", "delX", 104, 106, 38, 16, 1.0);
        let snv_other = make_test_variant("A", "G", "snvY", 105, 105, 35, 14, 1.0);

        let vars_per_amplicon = vec![
            make_vars_at_position(position, vec![deletion_1], None),
            make_vars_at_position(position, vec![deletion_2], None),
            make_vars_at_position(position, vec![snv_other], None),
        ];

        let output_lines = pipeline.run_amplicon_post_processor(
            &group_region,
            &vars_per_amplicon,
            &amplicon_regions,
            &HashSet::new(),
        );

        assert_eq!(output_lines.len(), 2);
        let parsed: Vec<Vec<String>> = output_lines
            .iter()
            .map(|line| parse_amplicon_output_fields(line))
            .collect();
        let deletion_fields = parsed
            .iter()
            .find(|fields| fields[33] == "Deletion")
            .expect("expected deletion output line");

        assert_eq!(deletion_fields[34], "1");
        assert_eq!(deletion_fields[35], "1");
        assert_eq!(deletion_fields[37], "0");
    }

    #[test]
    fn test_amplicon_pileup_empty_gvs_with_reference_uses_top_ref() {
        use std::collections::HashSet;

        ensure_test_scope_initialized();

        let pipeline = VarDictPipeline::new("sample").with_pileup(true);
        let position = 105i64;
        let amplicon_regions = vec![
            make_amplicon_region("G1", 105, 105),
            make_amplicon_region("G2", 105, 105),
        ];
        let group_region = amplicon_regions.last().unwrap().clone();

        let ref_high_cov = make_test_variant("A", "A", "A", 105, 105, 30, 30, 0.0);
        let ref_low_cov = make_test_variant("A", "A", "A", 105, 105, 10, 10, 0.0);

        let vars_per_amplicon = vec![
            make_vars_at_position(position, Vec::new(), Some(ref_high_cov)),
            make_vars_at_position(position, Vec::new(), Some(ref_low_cov)),
        ];

        let output_lines = pipeline.run_amplicon_post_processor(
            &group_region,
            &vars_per_amplicon,
            &amplicon_regions,
            &HashSet::new(),
        );

        assert_eq!(output_lines.len(), 1);
        let fields = parse_amplicon_output_fields(&output_lines[0]);
        assert_eq!(fields[5], "A");
        assert_eq!(fields[6], "A");
        assert_eq!(fields[7], "30");
        assert_eq!(fields[32], "chr1:105-105");
        assert_eq!(fields[34], "0");
        assert_eq!(fields[35], "0");
        assert_eq!(fields[36], "0");
        assert_eq!(fields[37], "0");
    }

    #[test]
    fn test_amplicon_pileup_empty_gvs_without_reference_outputs_empty_variant() {
        use std::collections::HashSet;

        ensure_test_scope_initialized();

        let pipeline = VarDictPipeline::new("sample").with_pileup(true);
        let amplicon_regions = vec![
            make_amplicon_region("G1", 105, 105),
            make_amplicon_region("G2", 105, 105),
        ];
        let group_region = amplicon_regions.last().unwrap().clone();

        let vars_per_amplicon = vec![HashMap::new(), HashMap::new()];

        let output_lines = pipeline.run_amplicon_post_processor(
            &group_region,
            &vars_per_amplicon,
            &amplicon_regions,
            &HashSet::new(),
        );

        assert_eq!(output_lines.len(), 1);
        let fields = parse_amplicon_output_fields(&output_lines[0]);
        assert_eq!(fields[3], "105");
        assert_eq!(fields[4], "105");
        assert_eq!(fields[5], "");
        assert_eq!(fields[6], "");
        assert_eq!(fields[7], "0");
        assert_eq!(fields[8], "0");
        assert_eq!(fields[34], "0");
        assert_eq!(fields[35], "0");
        assert_eq!(fields[36], "0");
        assert_eq!(fields[37], "0");
    }

    #[test]
    fn test_somatic_postprocessor_sample_specific_when_only_tumor_has_variant() {
        ensure_test_scope_initialized();

        let pipeline = VarDictPipeline::new("sample");
        let region = Region::new("chr1".to_string(), 100, 110, "GENE".to_string());
        let position = 105i64;

        let tumor_variant = make_test_variant("A", "T", "A>T", 105, 105, 100, 30, 0.30);
        let tumor_aligned = AlignedVarsData {
            aligned_variants: make_vars_at_position(position, vec![tumor_variant], None),
            ..AlignedVarsData::default()
        };
        let normal_aligned = AlignedVarsData::default();

        let output_lines = pipeline.run_somatic_post_processor(
            normal_aligned,
            tumor_aligned,
            &region,
            &std::collections::HashSet::new(),
        );

        assert_eq!(output_lines.len(), 1);
        let fields = parse_somatic_output_fields(&output_lines[0]);
        assert_eq!(fields[49], "SampleSpecific");
        assert_eq!(fields[2], "chr1");
        assert_eq!(fields[3], "105");
        assert_eq!(fields[6], "T");
        assert_eq!(fields[7], "100");
        assert_eq!(fields[25], "0");
    }

    #[test]
    fn test_somatic_postprocessor_germline_when_variant_present_in_both_samples() {
        ensure_test_scope_initialized();

        let pipeline = VarDictPipeline::new("sample");
        let region = Region::new("chr1".to_string(), 100, 110, "GENE".to_string());
        let position = 105i64;

        let tumor_variant = make_test_variant("A", "T", "A>T", 105, 105, 100, 60, 0.60);
        let normal_variant = make_test_variant("A", "T", "A>T", 105, 105, 100, 55, 0.55);

        let tumor_aligned = AlignedVarsData {
            aligned_variants: make_vars_at_position(position, vec![tumor_variant], None),
            ..AlignedVarsData::default()
        };
        let normal_aligned = AlignedVarsData {
            aligned_variants: make_vars_at_position(position, vec![normal_variant], None),
            ..AlignedVarsData::default()
        };

        let output_lines = pipeline.run_somatic_post_processor(
            normal_aligned,
            tumor_aligned,
            &region,
            &std::collections::HashSet::new(),
        );

        assert_eq!(output_lines.len(), 1);
        let fields = parse_somatic_output_fields(&output_lines[0]);
        assert_eq!(fields[49], "Germline");
        assert_eq!(fields[14], "0.6000");
        assert_eq!(fields[32], "0.5500");
    }

    #[test]
    fn test_somatic_combine_analysis_promotes_to_germline() {
        ensure_test_scope_initialized();

        let pipeline = VarDictPipeline::new("sample");
        let region = Region::new("chr1".to_string(), 100, 130, "GENE".to_string());
        let position = 105i64;

        let tumor_variant = make_test_variant("ATTT", "A", "-12AA", 105, 108, 100, 2, 0.02);
        let normal_ref = make_test_variant("A", "A", "A", 105, 105, 120, 120, 0.0);

        let tumor_aligned = AlignedVarsData {
            aligned_variants: make_vars_at_position(position, vec![tumor_variant.clone()], None),
            ..AlignedVarsData::default()
        };
        let normal_aligned = AlignedVarsData {
            aligned_variants: make_vars_at_position(position, Vec::new(), Some(normal_ref)),
            ..AlignedVarsData::default()
        };

        let mut combined_variant = tumor_variant.clone();
        combined_variant.total_pos_coverage = 150;
        combined_variant.position_coverage = 5;
        combined_variant.ref_forward_count = 60;
        combined_variant.ref_reverse_count = 70;
        combined_variant.vars_count_on_forward = 3;
        combined_variant.vars_count_on_reverse = 2;
        combined_variant.mean_position = 80.0;
        combined_variant.mean_quality = 70.0;
        combined_variant.mean_mapping_quality = 65.0;
        combined_variant.high_quality_reads_frequency = 0.3;
        combined_variant.extra_frequency = 0.1;
        combined_variant.nm = 2.0;
        combined_variant.genotype = "ATTT/A".to_string();

        let lookup = move |_: &str, _: i64, _: &str, max_read_length: usize| {
            SomaticCombineLookupResult {
                combined_variant: Some(combined_variant.clone()),
                max_read_length,
            }
        };

        let output_lines = pipeline.run_somatic_post_processor_with_combine_lookup(
            normal_aligned,
            tumor_aligned,
            &region,
            &std::collections::HashSet::new(),
            150,
            Some(&lookup),
        );

        assert_eq!(output_lines.len(), 1);
        let fields = parse_somatic_output_fields(&output_lines[0]);
        assert_eq!(fields[49], "Germline");
        assert_eq!(fields[25], "50");
        assert_eq!(fields[26], "3");
    }

    #[test]
    fn test_somatic_combine_analysis_false_skips_second_sample_variant() {
        ensure_test_scope_initialized();

        let pipeline = VarDictPipeline::new("sample");
        let region = Region::new("chr1".to_string(), 100, 130, "GENE".to_string());
        let position = 105i64;

        let tumor_ref = make_test_variant("A", "A", "A", 105, 105, 100, 100, 0.0);
        let normal_variant = make_test_variant("ATTT", "A", "-12AA", 105, 108, 120, 2, 0.02);

        let tumor_aligned = AlignedVarsData {
            aligned_variants: make_vars_at_position(position, Vec::new(), Some(tumor_ref)),
            ..AlignedVarsData::default()
        };
        let normal_aligned = AlignedVarsData {
            aligned_variants: make_vars_at_position(position, vec![normal_variant], None),
            ..AlignedVarsData::default()
        };

        let lookup = |_: &str, _: i64, _: &str, max_read_length: usize| SomaticCombineLookupResult {
            combined_variant: None,
            max_read_length,
        };

        let output_lines = pipeline.run_somatic_post_processor_with_combine_lookup(
            normal_aligned,
            tumor_aligned,
            &region,
            &std::collections::HashSet::new(),
            150,
            Some(&lookup),
        );

        assert!(output_lines.is_empty());
    }
}
