//! VarDict Pipeline - Real VarDict Java-equivalent pipeline for Simple Mode
//!
//! This module implements the actual VarDict pipeline flow:
//! ```text
//! BAM Record → CigarParser → VariantRealigner → StructuralVariantsProcessor → ToVarsBuilder → SimplePostProcessor → Output
//! ```
//!
//! This is the proper VarDict Simple Mode pipeline, replacing the simplified
//! `simple_variant_caller.rs` approach.

use std::collections::HashMap;
use std::sync::Arc;

use anyhow::{Result, Error};
use crackle_kit::tracing::{Level, event};
use rust_htslib::bam::Record;

use crate::data::reference::Reference;
use crate::data::region::Region;
use crate::data::shared_reference::SharedReferenceHandle;
use crate::data::bam_reader::BamReader;
use crate::mods::cigar_parser::CigarParser;
use crate::mods::output_variant::{SimpleOutputVariant, Region as OutputRegion};
use crate::mods::structural_variants_processor::{StructuralVariantsProcessor, RealignedVariationData};
use crate::mods::to_vars_builder::{ToVarsBuilder, Variant, VariationData, Vars, VarType, determine_genotype};
use crate::mods::simple_variant_caller::SimpleVarKey;
use crate::scopedata::global_read_only_scope::GlobalReadOnlyScope;
use crate::variants::variants::{VarDesc, Variant as RawVariant, SoftClip};

/// Data produced by CigarParser - mirrors Java VariationData
#[derive(Default)]
pub struct CigarParserOutput {
    /// Non-insertion variants by position
    pub non_insertion_vars: HashMap<i64, HashMap<VarDesc, RawVariant>>,
    /// Insertion variants by position (key is position before insertion)
    pub insertion_vars: HashMap<i64, HashMap<VarDesc, RawVariant>>,
    /// 5' end soft clips by position
    pub soft_clips_5end: HashMap<i64, SoftClip>,
    /// 3' end soft clips by position
    pub soft_clips_3end: HashMap<i64, SoftClip>,
    /// Reference coverage by position
    pub ref_coverage: HashMap<i64, usize>,
    /// Maximum read length seen
    pub max_read_len: usize,
    /// Discordant read count
    pub discordant_count: usize,
}

/// Data after realignment - mirrors Java RealignedVariationData
#[derive(Debug, Clone, Default)]
pub struct RealignedOutput {
    /// Non-insertion variants (may be modified by realigner)
    pub non_insertion_vars: HashMap<i64, HashMap<VarDesc, RawVariant>>,
    /// Insertion variants
    pub insertion_vars: HashMap<i64, HashMap<VarDesc, RawVariant>>,
    /// Reference coverage by position
    pub ref_coverage: HashMap<i64, usize>,
    /// Duplication rate
    pub duprate: f64,
    /// Maximum read length
    pub max_read_len: usize,
}

/// Final aligned variants data - mirrors Java AlignedVarsData
#[derive(Debug, Clone, Default)]
pub struct AlignedVarsData {
    /// Aligned variants by position
    pub aligned_variants: HashMap<i64, Vars>,
    /// Reference coverage
    pub ref_coverage: HashMap<i64, usize>,
}

/// VarDict Simple Mode Pipeline
///
/// Orchestrates the complete variant calling pipeline following Java VarDict flow.
pub struct VarDictPipeline {
    /// Sample name for output
    sample_name: String,
    /// Minimum allele frequency
    min_frequency: f64,
    /// Minimum base quality
    min_base_quality: u8,
    /// Minimum mapping quality
    min_mapping_quality: u8,
    /// Enable pileup mode (output reference positions too)
    do_pileup: bool,
}

impl VarDictPipeline {
    /// Create a new VarDict pipeline
    pub fn new(sample_name: &str) -> Self {
        VarDictPipeline {
            sample_name: sample_name.to_string(),
            min_frequency: 0.01,
            min_base_quality: 25,
            min_mapping_quality: 0,
            do_pileup: false,
        }
    }

    /// Set minimum allele frequency threshold
    pub fn with_min_frequency(mut self, freq: f64) -> Self {
        self.min_frequency = freq;
        self
    }

    /// Set minimum base quality threshold
    pub fn with_min_base_quality(mut self, qual: u8) -> Self {
        self.min_base_quality = qual;
        self
    }

    /// Set minimum mapping quality threshold
    pub fn with_min_mapping_quality(mut self, mapq: u8) -> Self {
        self.min_mapping_quality = mapq;
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
        // Extend reference loading to include flanking sequence for leftseq/rightseq
        // VarDict Java loads extra sequence around the region for flanking context
        let flank_size = 20usize; // 20bp flanking for leftseq/rightseq
        let extended_start = if region.start() > flank_size { 
            region.start() - flank_size 
        } else { 
            1 
        };
        let extended_end = region.end() + flank_size;
        
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
        let reference = Reference::new_with_start(ref_seq, extended_start as i64);

        // Fetch reads for this region
        bam_reader.fetch(region.chr(), region.start(), region.end())?;

        // Collect all records into a vector (needed for mutable iteration)
        let mut records = Vec::new();
        let mut record = Record::new();
        while bam_reader.read(&mut record).unwrap_or(false) {
            records.push(record.clone());
        }

        // Process through the pipeline
        self.process_region(records.into_iter(), region, &reference, instance)
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
        // Step 1: Parse CIGAR strings (CigarParser)
        let cigar_output = self.run_cigar_parser(records, region, reference, instance)?;
        
        event!(Level::DEBUG, "CigarParser output: {} non_insertion_vars, {} ref_coverage positions",
            cigar_output.non_insertion_vars.len(),
            cigar_output.ref_coverage.len());

        // Step 2: Realign soft clips and process structural variants
        let realigned_output = self.run_variant_realigner_and_sv_processor(cigar_output, region, reference)?;
        
        event!(Level::DEBUG, "Realigner output: {} non_insertion_vars, {} ref_coverage",
            realigned_output.non_insertion_vars.len(),
            realigned_output.ref_coverage.len());

        // Step 3: Build variant objects with statistics (ToVarsBuilder)
        let aligned_vars = self.run_to_vars_builder(realigned_output, reference)?;
        
        event!(Level::DEBUG, "ToVarsBuilder output: {} variants",
            aligned_vars.aligned_variants.len());

        // Step 4: Post-process and generate output (SimplePostProcessor)
        let output_lines = self.run_simple_post_processor(aligned_vars, region)?;
        
        event!(Level::DEBUG, "PostProcessor output: {} lines", output_lines.len());

        Ok(output_lines)
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

        // Extract results including soft clips
        Ok(CigarParserOutput {
            non_insertion_vars: cigar_parser.take_non_insertion_vars(),
            insertion_vars: cigar_parser.take_insertion_vars(),
            soft_clips_5end: cigar_parser.take_soft_clips_5end(),
            soft_clips_3end: cigar_parser.take_soft_clips_3end(),
            ref_coverage: cigar_parser.take_ref_coverage(),
            max_read_len: cigar_parser.get_max_read_len(),
            discordant_count: cigar_parser.get_discordant_count(),
        })
    }

    /// Step 2: Run VariantRealigner and StructuralVariantsProcessor
    fn run_variant_realigner_and_sv_processor(
        &self,
        input: CigarParserOutput,
        region: &Region,
        reference: &Reference,
    ) -> Result<RealignedOutput> {
        // TODO: Integrate actual VariantRealigner for soft clip realignment
        // For now, we pass through to StructuralVariantsProcessor
        
        // Convert CigarParserOutput to RealignedVariationData for SV processor
        let sv_input = RealignedVariationData {
            non_insertion_variants: input.non_insertion_vars,
            insertion_variants: input.insertion_vars,
            soft_clips_5end: input.soft_clips_5end,
            soft_clips_3end: input.soft_clips_3end,
            ref_coverage: input.ref_coverage,
            max_read_length: input.max_read_len,
            duprate: 0.0,
        };
        
        // Run StructuralVariantsProcessor (adjSNV always runs, SV detection is unimplemented)
        let sv_processor = StructuralVariantsProcessor::new(
            reference.ref_seq.clone(),
            region.start() as i64,
        );
        let processed = sv_processor.process(sv_input);
        
        // Convert back to RealignedOutput
        Ok(RealignedOutput {
            non_insertion_vars: processed.non_insertion_variants,
            insertion_vars: processed.insertion_variants,
            ref_coverage: processed.ref_coverage,
            duprate: processed.duprate,
            max_read_len: processed.max_read_length,
        })
    }

    /// Step 3: Run ToVarsBuilder to calculate statistics
    fn run_to_vars_builder(
        &self,
        input: RealignedOutput,
        reference: &Reference,
    ) -> Result<AlignedVarsData> {
        // Convert RawVariant data to the format ToVarsBuilder expects
        let mut aligned_variants: HashMap<i64, Vars> = HashMap::new();
        
        // First, collect reference forward/reverse counts from all non-insertion positions
        // This allows us to look up reference counts from adjacent positions when needed
        let mut ref_counts_by_pos: HashMap<i64, (usize, usize)> = HashMap::new();
        
        for (pos, var_map) in &input.non_insertion_vars {
            let actual_ref_base = reference.get(*pos).unwrap_or(b'N');
            for (desc, raw_var) in var_map {
                if let VarDesc::SNV { ref_base: read_base } = desc {
                    if *read_base == actual_ref_base {
                        ref_counts_by_pos.insert(*pos, (raw_var.alt_depth_fwd, raw_var.alt_depth_rev));
                        break;
                    }
                }
            }
        }

        // Process non-insertion variants
        for (pos, var_map) in input.non_insertion_vars {
            let vars = self.build_vars_at_position(pos, var_map, &input.ref_coverage, reference, &ref_counts_by_pos);
            if !vars.variants.is_empty() {
                aligned_variants.insert(pos, vars);
            }
        }

        // Process insertion variants
        for (pos, var_map) in input.insertion_vars {
            let vars = self.build_vars_at_position(pos, var_map, &input.ref_coverage, reference, &ref_counts_by_pos);
            if !vars.variants.is_empty() {
                aligned_variants.entry(pos)
                    .or_insert_with(Vars::default)
                    .variants
                    .extend(vars.variants);
            }
        }

        Ok(AlignedVarsData {
            aligned_variants,
            ref_coverage: input.ref_coverage,
        })
    }

    /// Build Vars struct from raw variant data at a position
    fn build_vars_at_position(
        &self,
        position: i64,
        var_map: HashMap<VarDesc, RawVariant>,
        ref_coverage: &HashMap<i64, usize>,
        reference: &Reference,
        ref_counts_by_pos: &HashMap<i64, (usize, usize)>,
    ) -> Vars {
        let total_coverage = ref_coverage.get(&position).copied().unwrap_or(0);
        let actual_ref_base = reference.get(position).unwrap_or(b'N');
        
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
        
        // Check if we need to look up reference counts from adjacent positions
        // This is only needed for insertion variants where the reference is at the anchor position
        let has_insertion_variant = var_map.iter().any(|(desc, _)| {
            matches!(desc, VarDesc::Ins { .. }) || 
            matches!(desc, VarDesc::Complex { ref_seq, alt_seq, .. } if alt_seq.len() > ref_seq.len())
        });
        
        // If we didn't find reference counts at this position, check adjacent positions
        // This handles cases like complex insertions where the reference is at the anchor position
        // Java does similar fallback logic in ToVarsBuilder.collectReferenceVariants
        // Only apply for insertion variants, not for MNVs or deletions
        if ref_fwd_count == 0 && ref_rev_count == 0 && has_insertion_variant {
            // Try position - 1 (the anchor for insertions)
            if let Some(&(fwd, rev)) = ref_counts_by_pos.get(&(position - 1)) {
                ref_fwd_count = fwd;
                ref_rev_count = rev;
            } else if let Some(&(fwd, rev)) = ref_counts_by_pos.get(&(position + 1)) {
                // Try position + 1 (Java's fallback for some cases)
                ref_fwd_count = fwd;
                ref_rev_count = rev;
            }
        }
        
        let mut variants = Vec::new();

        for (desc, raw_var) in var_map {
            let mut variant = self.convert_raw_variant(&desc, &raw_var, position, total_coverage, reference);
            
            // For non-reference variants, set the reference forward/reverse counts
            if variant.refallele != variant.varallele {
                variant.ref_forward_count = ref_fwd_count;
                variant.ref_reverse_count = ref_rev_count;
            }
            
            variants.push(variant);
        }

        Vars {
            variants,
            reference_variant: None,
            sv_flags: Default::default(),
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
    ) -> Variant {
        let total_count = raw.alt_depth_fwd + raw.alt_depth_rev;
        let frequency = if total_coverage > 0 {
            total_count as f64 / total_coverage as f64
        } else {
            0.0
        };

        // Look up reference base for this position
        let actual_ref_base = reference.get(position).unwrap_or(b'N');
        
        // Determine variant type and alleles based on VarDesc
        let (var_type, refallele, varallele) = match desc {
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
                let ref_char = actual_ref_base as char;
                let ins_str = String::from_utf8_lossy(seq);
                (
                    VarType::Insertion(ins_str.to_string()),
                    ref_char.to_string(),
                    format!("{}{}", ref_char, ins_str),
                )
            }
            VarDesc::Del { len, mismatch_seq, .. } => {
                // Check if there are mismatches following the deletion
                if !mismatch_seq.is_empty() {
                    // Complex variant: deletion + following mismatches
                    // Build reference sequence (deleted bases + following reference bases)
                    let mut ref_seq = Vec::new();
                    for i in 0..(*len as i64) {
                        if let Some(b) = reference.get(position + i) {
                            ref_seq.push(b);
                        }
                    }
                    // Add reference bases for the mismatch positions
                    for i in 0..mismatch_seq.len() {
                        if let Some(b) = reference.get(position + (*len as i64) + i as i64) {
                            ref_seq.push(b);
                        }
                    }
                    
                    // Alt sequence is just the mismatch sequence (deletion removes bases)
                    let alt_seq: Vec<u8> = mismatch_seq.iter().copied().collect();
                    
                    (
                        VarType::Complex {
                            insertion: String::from_utf8_lossy(&alt_seq).to_string(),
                            deletion: ref_seq.len(),
                        },
                        String::from_utf8_lossy(&ref_seq).to_string(),
                        String::from_utf8_lossy(&alt_seq).to_string(),
                    )
                } else {
                    // Simple deletion: ref allele includes deleted bases
                    let mut ref_str = String::new();
                    if let Some(prev_base) = reference.get(position - 1) {
                        ref_str.push(prev_base as char);
                    }
                    for i in 0..(*len as i64) {
                        if let Some(b) = reference.get(position + i) {
                            ref_str.push(b as char);
                        }
                    }
                    let var_str = if ref_str.len() > 0 { ref_str[0..1].to_string() } else { String::new() };
                    (
                        VarType::Deletion(*len as usize),
                        ref_str,
                        var_str,
                    )
                }
            }
            VarDesc::Complex { ref_seq, alt_seq } => {
                (
                    VarType::Complex {
                        insertion: String::from_utf8_lossy(alt_seq).to_string(),
                        deletion: ref_seq.len(),
                    },
                    String::from_utf8_lossy(ref_seq).to_string(),
                    String::from_utf8_lossy(alt_seq).to_string(),
                )
            }
        };

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
            raw.mean_pos / total_count as f64
        } else {
            0.0
        };
        
        let mean_quality = if total_count > 0 {
            raw.mean_qual / total_count as f64
        } else {
            0.0
        };
        
        let mean_mapping_quality = if total_count > 0 {
            raw.mean_mapq / total_count as f64
        } else {
            0.0
        };
        
        let nm_mean = if total_count > 0 {
            raw.nm / total_count as f64
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
        let (msi, msint) = if is_ref_call || is_insertion {
            (0.0, 0.0)
        } else {
            // Get deletion length from VarDesc, similar to how Java extracts from description string
            match desc {
                VarDesc::Del { len, .. } => {
                    // Pure deletion or deletion with mismatches - use deletion MSI calculation
                    let del_len = *len as usize;
                    if del_len > 0 {
                        self.detect_microsatellite(reference, position, del_len)
                    } else {
                        (0.0, 0.0)
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
                _ => (0.0, 0.0),
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
            position_coverage: total_coverage,
            frequency,
            high_quality_reads_frequency: if total_coverage > 0 {
                raw.high_qual_read_cnt as f64 / total_coverage as f64
            } else {
                0.0
            },
            mean_position,
            mean_quality,
            mean_mapping_quality,
            strand_bias_flag: calculate_strand_bias(raw.alt_depth_fwd, raw.alt_depth_rev),
            is_at_least_at_2_positions: raw.pstd,
            has_at_least_2_diff_qualities: raw.qstd,
            leftseq,
            rightseq,
            msi,
            msint,
            shift3: 0,
            nm: nm_mean,
            high_qual_read_cnt: raw.high_qual_read_cnt,
            low_qual_read_cnt: raw.low_qual_read_cnt,
            ref_forward_count: 0,  // Will be set by build_vars_at_position for non-ref variants
            ref_reverse_count: 0,
            genotype,
        }
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
    /// Returns (msi, msint) where:
    /// - msi: number of repeats (instability score)
    /// - msint: unit length (1 for homopolymer, 2 for dinucleotide, etc.)
    fn detect_microsatellite(&self, reference: &Reference, position: i64, del_len: usize) -> (f64, f64) {
        // Get left sequence (70 bases before position)
        let leftseq = self.get_reference_range(reference, position - 70, position - 1);
        
        // Get tseq (from position to position + deletion_len + 70)
        // For deletions, tseq1 is the deleted portion, tseq2 is what follows
        let tseq = self.get_reference_range(reference, position, position + (del_len as i64) - 1 + 70);
        
        if tseq.len() < del_len {
            return (0.0, 0.0);
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
        
        (msi, msint)
    }
    
    /// Detect microsatellite instability for SNP/MNP variants
    /// Java-compatible implementation for variants that don't start with + or -
    /// Java: tseq1 = joinRef(ref, position - 30, position + 1)
    ///       tseq2 = joinRef(ref, position + 2, position + 70)
    fn detect_microsatellite_snp(&self, reference: &Reference, position: i64) -> (f64, f64) {
        // tseq1 = reference from (position - 30) to (position + 1)
        let tseq1 = self.get_reference_range(reference, (position - 30).max(1), position + 1);
        
        // tseq2 = reference from (position + 2) to (position + 70)
        let tseq2 = self.get_reference_range(reference, position + 2, position + 70);
        
        // Call findMSI with no left sequence
        let (msi, msint, _) = self.find_msi(&tseq1, &tseq2, None);
        
        (msi, msint)
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

    /// Step 4: Run SimplePostProcessor to filter and format output
    fn run_simple_post_processor(
        &self,
        data: AlignedVarsData,
        region: &Region,
    ) -> Result<Vec<String>> {
        let mut output_lines = Vec::new();

        // Sort positions for consistent output
        let mut positions: Vec<i64> = data.aligned_variants.keys().copied().collect();
        positions.sort();

        // Convert to 1-based coordinates for output (BED is 0-based, VarDict output is 1-based)
        let output_region = OutputRegion {
            chr: region.chr().to_string(),
            start: region.start() as i64 + 1,  // +1 for 1-based output
            end: region.end() as i64,
            gene: region.gene().to_string(),
        };

        for position in positions {
            if let Some(vars) = data.aligned_variants.get(&position) {
                event!(Level::DEBUG, "[PostProcessor] Processing position {}: {} variants", position, vars.variants.len());
                
                // Skip positions outside region (unless SV)
                if position < region.start() as i64 || position > region.end() as i64 {
                    event!(Level::DEBUG, "[PostProcessor] Skipping position {} - outside region {}-{}", position, region.start(), region.end());
                    continue;
                }

                // Skip empty variants unless pileup mode
                if vars.variants.is_empty() {
                    event!(Level::DEBUG, "[PostProcessor] Position {} has 0 variants", position);
                    if !self.do_pileup {
                        continue;
                    }
                    // In pileup mode, output reference
                    if let Some(ref ref_var) = vars.reference_variant {
                        let output = SimpleOutputVariant::from_variant(
                            ref_var,
                            &output_region,
                            &self.sample_name,
                            "",
                        );
                        output_lines.push(output.to_string());
                    }
                    continue;
                }

                for variant in &vars.variants {
                    event!(Level::DEBUG, "[PostProcessor] Variant: pos={} ref={} alt={} freq={:.3} good={} type={:?}", 
                        variant.start_position, variant.refallele, variant.varallele, variant.frequency,
                        self.is_good_var(variant, vars.reference_variant.as_ref()), variant.vartype);
                    
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

                    // Apply quality filter (isGoodVar equivalent)
                    if !self.is_good_var(variant, vars.reference_variant.as_ref()) {
                        event!(Level::DEBUG, "[PostProcessor] Skipping - failed isGoodVar filter");
                        if !self.do_pileup {
                            continue;
                        }
                    }

                    // Apply frequency filter
                    if variant.frequency < self.min_frequency {
                        event!(Level::DEBUG, "[PostProcessor] Skipping - freq {:.3} < min {:.3}", variant.frequency, self.min_frequency);
                        if !self.do_pileup {
                            continue;
                        }
                    }

                    event!(Level::DEBUG, "[PostProcessor] Adding variant to output");

                    // Generate output
                    let output = SimpleOutputVariant::from_variant(
                        variant,
                        &output_region,
                        &self.sample_name,
                        "",
                    );
                    output_lines.push(output.to_string());
                }
            }
        }

        Ok(output_lines)
    }

    /// Quality filter - equivalent to Java Variant.isGoodVar()
    fn is_good_var(&self, variant: &Variant, _ref_variant: Option<&Variant>) -> bool {
        use crate::mods::to_vars_builder::StrandBiasFlag;

        // Check strand bias
        if variant.strand_bias_flag == StrandBiasFlag::StrongBias {
            // Allow high-frequency variants even with strand bias
            if variant.frequency < 0.1 {
                return false;
            }
        }

        // Check position variance
        if !variant.is_at_least_at_2_positions {
            // Single position variants are suspicious for low frequency
            if variant.frequency < 0.35 {
                return false;
            }
        }

        // Check quality variance
        if !variant.has_at_least_2_diff_qualities {
            // Single quality variants are suspicious for low frequency
            if variant.frequency < 0.35 {
                return false;
            }
        }

        true
    }
}

/// Calculate strand bias flag from forward/reverse counts
fn calculate_strand_bias(forward: usize, reverse: usize) -> crate::mods::to_vars_builder::StrandBiasFlag {
    use crate::mods::to_vars_builder::StrandBiasFlag;

    let total = forward + reverse;
    if total == 0 {
        return StrandBiasFlag::NoBias;
    }

    // For low counts (1 or 2 reads), strand bias is not meaningful
    // Java behavior: don't flag bias for low-count variants
    if total <= 2 {
        return StrandBiasFlag::NoBias;
    }

    let forward_ratio = forward as f64 / total as f64;

    // Strong bias: >90% on one strand
    if forward_ratio > 0.9 || forward_ratio < 0.1 {
        StrandBiasFlag::StrongBias
    } else if forward_ratio > 0.75 || forward_ratio < 0.25 {
        StrandBiasFlag::WeakBias
    } else {
        StrandBiasFlag::NoBias
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pipeline_creation() {
        let pipeline = VarDictPipeline::new("test_sample")
            .with_min_frequency(0.05)
            .with_min_base_quality(20)
            .with_min_mapping_quality(10);

        assert_eq!(pipeline.sample_name, "test_sample");
        assert!((pipeline.min_frequency - 0.05).abs() < 0.001);
        assert_eq!(pipeline.min_base_quality, 20);
        assert_eq!(pipeline.min_mapping_quality, 10);
    }

    #[test]
    fn test_strand_bias_calculation() {
        use crate::mods::to_vars_builder::StrandBiasFlag;

        // No reads
        assert_eq!(calculate_strand_bias(0, 0), StrandBiasFlag::NoBias);

        // Balanced
        assert_eq!(calculate_strand_bias(50, 50), StrandBiasFlag::NoBias);

        // Weak bias
        assert_eq!(calculate_strand_bias(80, 20), StrandBiasFlag::WeakBias);
        assert_eq!(calculate_strand_bias(20, 80), StrandBiasFlag::WeakBias);

        // Strong bias
        assert_eq!(calculate_strand_bias(95, 5), StrandBiasFlag::StrongBias);
        assert_eq!(calculate_strand_bias(5, 95), StrandBiasFlag::StrongBias);
        assert_eq!(calculate_strand_bias(100, 0), StrandBiasFlag::StrongBias);
    }

    #[test]
    fn test_genotype_determination() {
        // SNV: single base substitution uses "varallele/varallele" format (matching Java)
        assert_eq!(determine_genotype("C", "T", 0.9, None), "T/T");
        assert_eq!(determine_genotype("C", "T", 0.5, None), "T/T");
        assert_eq!(determine_genotype("C", "T", 0.3, None), "T/T");  // Frequency is not used for SNVs
        assert_eq!(determine_genotype("G", "G", 1.0, None), "G/G");  // Ref call
    }

    #[test]
    fn test_is_good_var_strand_bias() {
        let pipeline = VarDictPipeline::new("test");
        
        // Good variant - no bias, good position/quality variance
        let mut good_var = Variant::default();
        good_var.strand_bias_flag = crate::mods::to_vars_builder::StrandBiasFlag::NoBias;
        good_var.is_at_least_at_2_positions = true;
        good_var.has_at_least_2_diff_qualities = true;
        good_var.frequency = 0.3;
        assert!(pipeline.is_good_var(&good_var, None));

        // Bad variant - strong bias, low frequency
        let mut bad_var = Variant::default();
        bad_var.strand_bias_flag = crate::mods::to_vars_builder::StrandBiasFlag::StrongBias;
        bad_var.is_at_least_at_2_positions = true;
        bad_var.has_at_least_2_diff_qualities = true;
        bad_var.frequency = 0.05; // Low frequency + strong bias = bad
        assert!(!pipeline.is_good_var(&bad_var, None));

        // High frequency can overcome strand bias
        let mut high_freq_bias = Variant::default();
        high_freq_bias.strand_bias_flag = crate::mods::to_vars_builder::StrandBiasFlag::StrongBias;
        high_freq_bias.is_at_least_at_2_positions = true;
        high_freq_bias.has_at_least_2_diff_qualities = true;
        high_freq_bias.frequency = 0.5; // High enough to pass despite bias
        assert!(pipeline.is_good_var(&high_freq_bias, None));
    }
}
