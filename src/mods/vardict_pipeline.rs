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
use crate::mods::variant_realigner::VariantRealigner;
use crate::scopedata::global_read_only_scope::instance;
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
    /// Variants by position
    pub aligned_variants: HashMap<i64, Vars>,
    /// Reference coverage by position
    pub ref_coverage: HashMap<i64, usize>,
}

/// Main VarDict pipeline configuration
pub struct VarDictPipeline {
    pub sample_name: String,
    pub min_frequency: f64,
    pub min_base_quality: f64,
    pub min_mapping_quality: u8,
    pub do_pileup: bool,
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

        // Get SAM filter from instance configuration
        let sam_filter = instance.conf.sam_filter;

        let records = self.collect_filtered_records(region, bam_reader, sam_filter)?.0;

        // Process through the pipeline
        self.process_region(records.into_iter(), region, &reference, instance)
    }

    pub(crate) fn collect_filtered_records(
        &self,
        region: &Region,
        bam_reader: &mut BamReader,
        sam_filter: u32,
    ) -> Result<(Vec<Record>, Vec<String>)> {
        let mut records = Vec::new();
        let mut lines = Vec::new();
        let mut record = Record::new();

        bam_reader.fetch(region.chr(), region.start(), region.end())?;

        while bam_reader.read(&mut record).unwrap_or(false) {
            // 1. Java SamView.read(): Skip records that match the filter flags
            if sam_filter != 0 {
                if (record.flags() & (sam_filter as u16)) != 0 {
                    continue;
                }
            }

            // 2. Java preprocessRecord line 117: Ignore low mapping quality reads
            if self.min_mapping_quality > 0 && record.mapq() < self.min_mapping_quality {
                continue;
            }

            // 3. Java preprocessRecord line 122: Skip not primary alignment reads
            const SECONDARY_ALIGNMENT: u16 = 0x100;
            if (record.flags() & SECONDARY_ALIGNMENT) != 0 && sam_filter != 0 {
                continue;
            }

            // 4. Java preprocessRecord line 124: Skip reads where sequence is not stored in read
            let seq = record.seq();
            if seq.len() == 0 || (seq.len() == 1 && seq.as_bytes()[0] == b'*') {
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

            records.push(record.clone());
        }

        Ok((records, lines))
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
        let start_total = std::time::Instant::now();
        
        // Step 1: Parse CIGAR strings (CigarParser)
        let start_cigar = std::time::Instant::now();
        let cigar_output = self.run_cigar_parser(records, region, reference, instance)?;
        let elapsed_cigar = start_cigar.elapsed();
        
        event!(Level::INFO, "[TIMING] CigarParser: {:.3}s - {} non_insertion_vars, {} ref_coverage positions",
            elapsed_cigar.as_secs_f64(),
            cigar_output.non_insertion_vars.len(),
            cigar_output.ref_coverage.len());

        // Step 2: Realign soft clips and process structural variants
        let start_realign = std::time::Instant::now();
        let realigned_output = self.run_variant_realigner_and_sv_processor(cigar_output, region, reference)?;
        let elapsed_realign = start_realign.elapsed();
        
        event!(Level::INFO, "[TIMING] VariantRealigner+SVProcessor: {:.3}s - {} non_insertion_vars, {} ref_coverage",
            elapsed_realign.as_secs_f64(),
            realigned_output.non_insertion_vars.len(),
            realigned_output.ref_coverage.len());

        // Step 3: Build variant objects with statistics (ToVarsBuilder)
        let start_tovars = std::time::Instant::now();
        let aligned_vars = self.run_to_vars_builder(realigned_output, reference)?;
        let elapsed_tovars = start_tovars.elapsed();
        
        event!(Level::INFO, "[TIMING] ToVarsBuilder: {:.3}s - {} variants",
            elapsed_tovars.as_secs_f64(),
            aligned_vars.aligned_variants.len());

        // Step 4: Post-process and generate output (SimplePostProcessor)
        let start_post = std::time::Instant::now();
        let output_lines = self.run_simple_post_processor(aligned_vars, region)?;
        let elapsed_post = start_post.elapsed();
        
        event!(Level::INFO, "[TIMING] PostProcessor: {:.3}s - {} lines",
            elapsed_post.as_secs_f64(),
            output_lines.len());

        // let elapsed_total = start_total.elapsed();
        
        // // Log timing summary for regions taking >10ms
        // if elapsed_total.as_millis() > 10 {
        //     event!(Level::TRACE, "[REGION] {}:{}-{} total={:.3}s cigar={:.3}s realign={:.3}s tovars={:.3}s post={:.3}s",
        //         region.chr(), region.start(), region.end(),
        //         elapsed_total.as_secs_f64(),
        //         elapsed_cigar.as_secs_f64(),
        //         elapsed_realign.as_secs_f64(),
        //         elapsed_tovars.as_secs_f64(),
        //         elapsed_post.as_secs_f64());
        // }

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
        let mut sv_input = RealignedVariationData {
            non_insertion_variants: input.non_insertion_vars,
            insertion_variants: input.insertion_vars,
            soft_clips_5end: input.soft_clips_5end,
            soft_clips_3end: input.soft_clips_3end,
            ref_coverage: input.ref_coverage,
            max_read_length: input.max_read_len,
            duprate: 0.0,
        };

        // Perform minimal deletion realignment using soft clips when enabled
        // Re-enable realigner to match Java behavior
        if instance().conf.perform_local_realignment {
            let realigner = VariantRealigner::new(reference.ref_seq.clone(), reference.region_start);
            realigner.process_deletions(&mut sv_input);
        }
        
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
        let mut hicov_by_pos: HashMap<i64, usize> = HashMap::new();
        
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

            let hicov: usize = var_map.values().map(|raw_var| raw_var.high_qual_read_cnt).sum();
            hicov_by_pos.insert(*pos, hicov);
        }

        // Process non-insertion variants
        for (pos, var_map) in input.non_insertion_vars {
            let vars = self.build_vars_at_position(
                pos,
                var_map,
                &input.ref_coverage,
                reference,
                &ref_counts_by_pos,
                &hicov_by_pos,
            );
            if !vars.variants.is_empty() {
                aligned_variants.insert(pos, vars);
            }
        }

        // Process insertion variants
        for (pos, var_map) in input.insertion_vars {
            let vars = self.build_vars_at_position(
                pos,
                var_map,
                &input.ref_coverage,
                reference,
                &ref_counts_by_pos,
                &hicov_by_pos,
            );
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
        hicov_by_pos: &HashMap<i64, usize>,
    ) -> Vars {
        use crate::mods::to_vars_builder::{check_strand_bias, StrandBiasFlag, StrandBiasValue};
        
        let total_coverage = ref_coverage.get(&position).copied().unwrap_or(0);
        let actual_ref_base = reference.get(position).unwrap_or(b'N');
        let position_hicov = hicov_by_pos.get(&position).copied().unwrap_or(0);
        
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
        
        // Calculate reference strand bias (used as first part of "refBias;varBias" flag)
        let ref_strand_bias = check_strand_bias(ref_fwd_count, ref_rev_count);
        
        let mut variants = Vec::new();
        let mut reference_variant_opt = None;

        for (desc, raw_var) in var_map {
            let mut variant = self.convert_raw_variant(
                &desc,
                &raw_var,
                position,
                total_coverage,
                reference,
                position_hicov,
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
            } else {
                // This is a reference call - save it as the reference variant
                reference_variant_opt = Some(variant.clone());
            }
            
            variants.push(variant);
        }

        Vars {
            variants,
            reference_variant: reference_variant_opt,
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
        position_hicov: usize,
    ) -> Variant {
        let mut position = position;
        let total_count = raw.alt_depth_fwd + raw.alt_depth_rev;
        let frequency = if total_coverage > 0 {
            total_count as f64 / total_coverage as f64
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
                    // For complex variants in VarDict raw output format, we don't use anchor base
                    // Position is stored as anchor (position - 1), so actual deletion starts at position + 1
                    let del_start = position + 1;
                    // Update position to be the actual deletion start for output
                    position = del_start;
                    
                    // Build reference sequence (deleted bases + following reference bases that get replaced)
                    let mut ref_seq = Vec::new();
                    for i in 0..(*len as i64) {
                        if let Some(b) = reference.get(del_start + i) {
                            ref_seq.push(b);
                        }
                    }
                    // Add reference bases for the mismatch positions (after the deleted stretch)
                    for i in 0..mismatch_seq.len() {
                        if let Some(b) = reference.get(del_start + (*len as i64) + i as i64) {
                            ref_seq.push(b);
                        }
                    }

                    // Alt sequence is just the mismatch sequence (no anchor)
                    let alt_seq: Vec<u8> = mismatch_seq.to_vec();

                    (
                        VarType::Complex {
                            insertion: String::from_utf8_lossy(&alt_seq).to_string(),
                            deletion: ref_seq.len(),
                        },
                        String::from_utf8_lossy(&ref_seq).to_string(),
                        String::from_utf8_lossy(&alt_seq).to_string(),
                    )
                } else {
                    // Simple deletion: ref allele includes anchor + deleted bases
                    let mut ref_str = String::new();
                    if let Some(anchor_base) = reference.get(position) {
                        ref_str.push(anchor_base as char);
                    }
                    for i in 1..=(*len as i64) {
                        if let Some(b) = reference.get(position + i) {
                            ref_str.push(b as char);
                        }
                    }
                    // Alt allele is just the anchor base
                    let var_str = if !ref_str.is_empty() {
                        ref_str[0..1].to_string()
                    } else {
                        String::new()
                    };
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
            position_coverage: total_coverage,
            frequency,
            high_quality_reads_frequency: if total_coverage > 0 {
                if hicov > 0 {
                    raw.high_qual_read_cnt as f64 / hicov as f64
                } else {
                    0.0
                }
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

        // Output uses the same coordinate base as the parsed regions
        let output_region = OutputRegion {
            chr: region.chr().to_string(),
            start: region.start() as i64 + 1,
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
    /// 
    /// Java checks: frequency >= conf.freq, hicnt >= conf.minr, 
    /// meanPosition >= conf.readPosFilter, meanQuality >= conf.goodq,
    /// highQualityToLowQualityRatio >= conf.qratio
    fn is_good_var(&self, variant: &Variant, ref_variant: Option<&Variant>) -> bool {
        // Check frequency (already checked separately, but kept for completeness)
        if variant.frequency < self.min_frequency {
            return false;
        }

        // Check high-quality read count (minr = 2 by default)
        let min_reads = 2;
        if variant.high_qual_read_cnt < min_reads {
            return false;
        }

        // Check mean position (readPosFilter = 5 by default)
        let read_pos_filter = 5.0;
        if variant.mean_position < read_pos_filter {
            return false;
        }

        // Check mean quality (goodq = 22.5 by default)
        let goodq = 22.5;
        if variant.mean_quality < goodq {
            return false;
        }

        // Check high-quality to low-quality ratio (qratio = 1.5 by default)
        let qratio = 1.5;
        if variant.low_qual_read_cnt > 0 {
            let ratio = variant.high_qual_read_cnt as f64 / variant.low_qual_read_cnt as f64;
            if ratio < qratio {
                return false;
            }
        }

        // Check mapping quality vs reference (Java logic for low-frequency variants)
        if let Some(ref_var) = ref_variant {
            if ref_var.high_qual_read_cnt >= min_reads && variant.frequency < 0.25 {
                let d = variant.mean_mapping_quality + variant.refallele.len() as f64 + variant.varallele.len() as f64;
                let f = (1.0 + d) / (ref_var.mean_mapping_quality + 1.0);
                if (d - 2.0 < 5.0 && ref_var.mean_mapping_quality > 20.0) || f < 0.25 {
                    return false;
                }
            }
        }

        // High frequency variants pass without further checks
        if variant.frequency > 0.30 {
            return true;
        }

        // Check mapping quality threshold for low-frequency variants (mapq = 0 by default)
        // Java: if (meanMappingQuality < instance().conf.mapq) return false;
        // With default mapq=0, this is effectively a no-op
        let mapq_threshold = 0.0;
        if variant.mean_mapping_quality < mapq_threshold {
            return false;
        }

        // MSI (microsatellite instability) filters
        // Java: if (msi >= 15 && frequency <= monomerMsiFrequency(0.005) && msint == 1) return false
        let monomer_msi_freq = 0.005;
        if variant.msi >= 15.0 && variant.frequency <= monomer_msi_freq && variant.msint == 1.0 {
            return false;
        }
        // Java: if (msi >= 12 && frequency <= nonMonomerMsiFrequency(0.002) && msint > 1) return false
        let non_monomer_msi_freq = 0.002;
        if variant.msi >= 12.0 && variant.frequency <= non_monomer_msi_freq && variant.msint > 1.0 {
            return false;
        }

        // Strand bias filter: "2;1" pattern (ref good, var biased) at low frequency for small variants
        // Java: if (strandBiasFlag.equals("2;1") && frequency < 0.20d)
        //         if (type == null || type.equals("SNV") || (refallele.length() < 3 && varallele.length() < 3))
        //           return false
        if variant.strand_bias_flag.is_ref_good_var_biased() && variant.frequency < 0.20 {
            let is_small_variant = match &variant.vartype {
                VarType::SNV(_) => true,
                _ => variant.refallele.len() < 3 && variant.varallele.len() < 3,
            };
            if is_small_variant {
                return false;
            }
        }

        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
        assert!(pipeline.is_good_var(&good_var, None));

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
        assert!(!pipeline.is_good_var(&bad_var, None));

        // High frequency can overcome strand bias (2;1 but freq > 0.20)
        let mut high_freq_bias = Variant::default();
        high_freq_bias.strand_bias_flag = StrandBiasFlag::new(StrandBiasValue::NoBias, StrandBiasValue::HasBias);
        high_freq_bias.is_at_least_at_2_positions = true;
        high_freq_bias.has_at_least_2_diff_qualities = true;
        high_freq_bias.frequency = 0.5; // High enough to pass despite 2;1 pattern
        high_freq_bias.high_qual_read_cnt = 5;
        high_freq_bias.mean_position = 10.0;
        high_freq_bias.mean_quality = 30.0;
        assert!(pipeline.is_good_var(&high_freq_bias, None));
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
}
