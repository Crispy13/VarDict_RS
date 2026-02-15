//! VarDict-rs: A Rust implementation of VarDict variant caller (Simple Mode)
//!
//! Usage: vardict -G <reference.fa> -b <input.bam> [options] <region or BED file>

use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::PathBuf;

use anyhow::{Context, Result, anyhow};
use clap::Parser;

use crackle_kit::tracing::level_filters::LevelFilter;
use crackle_kit::tracing_kit::{setup_logging_stderr_only, setup_logging_stderr_only_verbose};
use vardict_rs::data::bam_reader::BamReader;
use vardict_rs::data::region::Region;
use vardict_rs::mods::pipeline::{Pipeline, PipelineConfig};
use vardict_rs::mods::vardict_pipeline::VarDictPipeline;

const DEFAULT_AMPLICON_PARAMETERS: &str = "10:0.95";

/// VarDict-rs: Variant caller for NGS data (Simple Mode)
#[derive(Parser, Debug)]
#[command(name = "vardict")]
#[command(author = "VarDict-rs Authors")]
#[command(version = "0.1.0")]
#[command(about = "A sensitive variant caller for NGS data", long_about = None)]
struct Args {
    /// Reference FASTA file (must be indexed with .fai)
    #[arg(short = 'G', long = "ref", required = true)]
    reference: PathBuf,

    /// Indexed BAM file
    #[arg(short = 'b', long = "bam", required = true)]
    bam: PathBuf,

    /// Region of interest (chr:start-end) or BED file
    /// If a file path, reads regions from BED format
    /// If chr:start-end format, processes single region
    #[arg(short = 'R', long = "region")]
    region: Option<String>,

    /// BED file with regions to process (positional argument)
    #[arg()]
    bed_file: Option<PathBuf>,

    /// Sample name (overrides automatic detection from BAM filename)
    #[arg(short = 'N', long = "name")]
    sample_name: Option<String>,

    /// Minimum allele frequency threshold (default: 0.01)
    #[arg(short = 'f', long = "freq", default_value = "0.01")]
    min_frequency: f64,

    /// Minimum number of variant reads (default: 2)
    #[arg(short = 'r', long = "minr", default_value = "2")]
    min_variant_reads: usize,

    /// Minimum base quality for a base to be considered (default: 22.5)
    #[arg(short = 'q', long = "qual", default_value = "22.5")]
    min_base_quality: f64,

    /// Minimum mapping quality for a read to be considered (default: 0)
    #[arg(short = 'Q', long = "mapq", default_value = "0")]
    min_mapping_quality: u8,

    /// Number of nucleotides to extend regions (Java: -x)
    #[arg(short = 'x', long = "extend", default_value = "0")]
    number_nucleotide_to_extend: i32,

    /// Reference extension for fetching sequence (Java: -Y)
    #[arg(short = 'Y', long = "reference-extension", default_value = "1200")]
    reference_extension: i32,

    /// Print header line
    #[arg(short = 'H', long = "header")]
    print_header: bool,

    /// Extension of bp to look for mismatches after indel (default: 2)
    #[arg(short = 'X', long = "vext", default_value = "2")]
    vext: i32,

    /// If set, reads with mismatches more than INT will be filtered and ignored (default: 8)
    #[arg(short = 'm', long = "mismatch", default_value = "8")]
    mismatch: i32,

    /// The hexical to filter reads (samtools style). Default: 0x504
    /// Filters: 2nd alignments, unmapped, duplicates. Use 0 to disable.
    #[arg(short = 'F', long = "filter", default_value = "0x504")]
    sam_filter: String,

    /// Downsampling fraction (Java: -Z). When set, randomly drop reads by this fraction.
    #[arg(short = 'Z', long = "downsample")]
    downsampling: Option<f64>,

    /// Turn off structural variant calling
    #[arg(short = 'U', long = "nosv")]
    no_sv: bool,

    /// Amplicon mode parameters (Java: -a), e.g. "10:0.95"
    #[arg(short = 'a', long = "amplicon")]
    amplicon_based_calling: Option<String>,

    /// Debug mode - print additional information
    #[arg(short = 'D', long = "debug")]
    debug: bool,

    /// Output all variants including reference calls (pileup mode)
    #[arg(short = 'p', long = "pileup")]
    pileup: bool,

    /// Indicate whether coordinates are zero-based: 1 for zero-based, 0 for one-based
    #[arg(short = 'z', long = "zero", value_parser = clap::value_parser!(u8).range(0..=1))]
    zero_based: Option<u8>,

    /// Perform local realignment (default: 1). Use 0 to disable.
    #[arg(short = 'k', long = "realign", default_value = "1")]
    local_realignment: u8,

    /// The column for chromosome in BED file (default: 1)
    #[arg(short = 'c', long = "col-chr", default_value = "1")]
    col_chr: usize,

    /// The column for start position in BED file (default: 2)
    #[arg(short = 'S', long = "col-start", default_value = "2")]
    col_start: usize,

    /// The column for end position in BED file (default: 3)
    #[arg(short = 'E', long = "col-end", default_value = "3")]
    col_end: usize,

    /// The column for gene name in BED file (default: 4)
    #[arg(short = 'g', long = "col-gene", default_value = "4")]
    col_gene: usize,

    /// Number of threads for parallel processing (default: 1)
    #[arg(short = 't', long = "threads", default_value = "1")]
    num_threads: usize,

    /// Remove duplicated reads (Java: -t). TODO: option only; logic not yet implemented.
    #[arg(long = "remove-duplicates")]
    remove_duplicates: bool,


    /// Log level
    #[arg(long, default_value_t = LevelFilter::WARN)]
    log_level: LevelFilter,
}

fn main() -> Result<()> {
    let args = Args::parse();

    setup_logging_stderr_only(args.log_level)?;

    // Validate input files exist
    if !args.reference.exists() {
        return Err(anyhow!("Reference file not found: {:?}", args.reference));
    }
    if !args.bam.exists() {
        return Err(anyhow!("BAM file not found: {:?}", args.bam));
    }

    // Check for FASTA index
    let fai_path = args.reference.with_extension("fa.fai");
    let fai_path2 = {
        let mut p = args.reference.clone();
        p.set_file_name(format!("{}.fai", args.reference.file_name().unwrap().to_string_lossy()));
        p
    };
    if !fai_path.exists() && !fai_path2.exists() {
        return Err(anyhow!(
            "Reference index (.fai) not found. Please run: samtools faidx {:?}",
            args.reference
        ));
    }

    // Check for BAM index
    let bai_path = args.bam.with_extension("bam.bai");
    let bai_path2 = {
        let mut p = args.bam.clone();
        p.set_file_name(format!("{}.bai", args.bam.file_name().unwrap().to_string_lossy()));
        p
    };
    if !bai_path.exists() && !bai_path2.exists() {
        return Err(anyhow!(
            "BAM index (.bai) not found. Please run: samtools index {:?}",
            args.bam
        ));
    }

    // Determine sample name - clone bam path first since we'll consume sample_name
    let bam_path = args.bam.clone();
    let sample_name = args.sample_name.clone().unwrap_or_else(|| {
        bam_path
            .file_stem()
            .and_then(|s| s.to_str())
            .map(|s| {
                // Extract sample name before first underscore
                s.split('_').next().unwrap_or(s).to_string()
            })
            .unwrap_or_else(|| "SAMPLE".to_string())
    });

    // Get regions to process
    let region_load = get_regions(&args)?;
    if region_load.regions.is_empty() {
        return Err(anyhow!("No regions specified. Provide -R option or a BED file."));
    }

    // Build pipeline configuration
    let config = PipelineConfig::builder()
        .sample_name(sample_name.clone())
        .min_frequency(args.min_frequency)
        .min_variant_reads(args.min_variant_reads)
        .min_base_quality(args.min_base_quality)
        .min_mapping_quality(args.min_mapping_quality)
        .pileup(args.pileup)
        .build();

    let execution_mode = resolve_execution_mode(&args, &region_load.amplicon_based_calling);

    // Print header if requested
    if args.print_header {
        if execution_mode == ExecutionMode::Amplicon {
            println!("{}", vardict_rs::mods::output_variant::get_amplicon_header_line());
        } else {
            let pipeline = Pipeline::new(config.clone());
            println!("{}", pipeline.get_header());
        }
    }

    // Always use SharedReference (loaded into memory for fast access)
    run_variant_calling(
        &args,
        config,
        region_load.regions,
        region_load.amplicon_based_calling,
        region_load.amplicon_region_groups,
        execution_mode,
    )?;

    Ok(())
}

/// Run variant calling using SharedReference (loaded into memory)
/// 
/// SharedReference is the default for both single and multi-threaded modes.
/// The reference is loaded once and shared across all threads for fast access.
fn run_variant_calling(
    args: &Args,
    config: PipelineConfig,
    regions: Vec<Region>,
    amplicon_based_calling: Option<String>,
    amplicon_region_groups: Option<Vec<Vec<Region>>>,
    execution_mode: ExecutionMode,
) -> Result<()> {
    use vardict_rs::data::shared_reference::load_shared_reference_chroms;
    use vardict_rs::mods::parallel_pipeline::ParallelPipeline;
    use vardict_rs::scopedata::global_read_only_scope::{GlobalReadOnlyScope, INSTANCE};
    use vardict_rs::conf::Configuration;
    use std::sync::Arc;
    use std::time::Instant;

    let start_total = Instant::now();
    let num_threads = args.num_threads.max(1);
    
    if args.debug {
        eprintln!("Loading reference genome into memory...");
    }

    // Get unique chromosomes from regions
    let start_ref_load = Instant::now();
    let chroms: std::collections::HashSet<&str> = regions.iter()
        .map(|r| r.chr())
        .collect();
    let chrom_vec: Vec<&str> = chroms.into_iter().collect();

    // Load only the needed chromosomes for efficiency
    let reference = load_shared_reference_chroms(
        args.reference.to_str().unwrap(),
        &chrom_vec,
    ).context("Failed to load reference genome")?;
    
    let elapsed_ref_load = start_ref_load.elapsed();

    // Initialize GlobalReadOnlyScope (required by VarDictPipeline)
    // Must be done AFTER loading reference to populate chr_lens
    let mut conf = Configuration::default();
    let sam_filter = if let Some(hex) = args.sam_filter.strip_prefix("0x")
        .or_else(|| args.sam_filter.strip_prefix("0X"))
    {
        u32::from_str_radix(hex, 16)
            .context("Failed to parse sam_filter as hex")?
    } else {
        args.sam_filter
            .parse::<u32>()
            .context("Failed to parse sam_filter as decimal")?
    };
    conf.goodq = args.min_base_quality;
    conf.freq = if args.pileup { -1.0 } else { args.min_frequency };
    conf.minr = if args.pileup { 0 } else { args.min_variant_reads };
    conf.vext = args.vext;
    conf.mismatch = args.mismatch;
    conf.mapping_quality = if args.min_mapping_quality > 0 {
        Some(args.min_mapping_quality)
    } else {
        None
    };
    conf.sam_filter = sam_filter;
    conf.downsampling = args.downsampling;
    conf.remove_duplicated_reads = args.remove_duplicates;
    conf.disable_sv = args.no_sv;
    conf.amplicon_based_calling = amplicon_based_calling.clone();
    conf.perform_local_realignment = args.local_realignment == 1;
    conf.number_nucleotide_to_extend = args.number_nucleotide_to_extend;
    conf.reference_extension = args.reference_extension;
    let mut scope = GlobalReadOnlyScope::default();
    scope.amplicon_based_calling = conf.amplicon_based_calling.clone();
    scope.conf = conf;
    scope.chr_lens = reference.get_chromosome_lengths();
    scope.bam_paths = vec![args.bam.to_string_lossy().to_string()];
    let _ = INSTANCE.set(scope);

    let region_batches = select_region_batches_for_execution(
        execution_mode,
        regions,
        amplicon_region_groups,
    );

    if args.debug {
        eprintln!("[TIMING] Reference loading: {:.3}s", elapsed_ref_load.as_secs_f64());
        eprintln!("Loaded {} chromosome(s), {:.2} MB total",
            reference.num_chromosomes(),
            reference.total_size() as f64 / 1_048_576.0);
        eprintln!("Execution mode: {:?}", execution_mode);
        eprintln!("Execution batches: {}", region_batches.len());
        if num_threads > 1 {
            eprintln!(
                "Processing {} regions with {} threads...",
                region_batches.iter().map(Vec::len).sum::<usize>(),
                num_threads
            );
        } else {
            eprintln!(
                "Processing {} regions...",
                region_batches.iter().map(Vec::len).sum::<usize>()
            );
        }
    }

    // Process regions
    let start_processing = Instant::now();
    let bam_path = args.bam.to_str().unwrap().to_string();

    match execution_mode {
        ExecutionMode::Simple => {
            let pipeline = ParallelPipeline::new(reference, config, num_threads);
            let mut results = Vec::new();
            for batch in region_batches {
                let mut batch_results = pipeline.process_regions_vardict(bam_path.clone(), batch);
                results.append(&mut batch_results);
            }

            let elapsed_processing = start_processing.elapsed();

            let start_output = Instant::now();
            let mut stdout = io::stdout().lock();
            for result in results {
                if let Some(error) = result.error {
                    if args.debug {
                        eprintln!("Error processing {}:{}-{}: {}",
                            result.region.chr(), result.region.start(), result.region.end(), error);
                    }
                } else {
                    for line in result.output_lines {
                        writeln!(stdout, "{}", line)?;
                    }
                }
            }
            let elapsed_output = start_output.elapsed();

            let elapsed_total = start_total.elapsed();

            if args.debug {
                eprintln!("[TIMING] Processing all regions: {:.3}s", elapsed_processing.as_secs_f64());
                eprintln!("[TIMING] Output writing: {:.3}s", elapsed_output.as_secs_f64());
                eprintln!("[TIMING] TOTAL execution: {:.3}s", elapsed_total.as_secs_f64());
            }
        }
        ExecutionMode::Amplicon => {
            use std::collections::HashSet;

            let vardict_pipeline = VarDictPipeline::new(&config.sample_name)
                .with_min_frequency(config.min_frequency)
                .with_min_base_quality(config.quality_threshold)
                .with_min_mapping_quality(config.mapq_threshold);

            let global_scope = Arc::new(INSTANCE.get().expect("GlobalReadOnlyScope not initialized").clone());
            let mut stdout = io::stdout().lock();

            for amplicon_group in region_batches {
                if amplicon_group.is_empty() {
                    continue;
                }

                let mut vars_per_amplicon = Vec::with_capacity(amplicon_group.len());
                let mut splice: HashSet<String> = HashSet::new();

                for region in &amplicon_group {
                    let mut bam_reader = BamReader::open(&bam_path)?;
                    let aligned_output = vardict_pipeline.process_region_to_aligned_vars_from_bam(
                        region,
                        &reference,
                        &mut bam_reader,
                        Arc::clone(&global_scope),
                    )?;
                    vars_per_amplicon.push(aligned_output.aligned_vars.aligned_variants);
                    splice.extend(aligned_output.splice.into_iter());
                }

                let current_region = amplicon_group
                    .last()
                    .cloned()
                    .ok_or_else(|| anyhow!("Empty amplicon group encountered"))?;
                let output_lines = vardict_pipeline.run_amplicon_post_processor(
                    &current_region,
                    &vars_per_amplicon,
                    &amplicon_group,
                    &splice,
                );

                for line in output_lines {
                    writeln!(stdout, "{}", line)?;
                }
            }

            let elapsed_processing = start_processing.elapsed();
            let elapsed_total = start_total.elapsed();

            if args.debug {
                eprintln!("[TIMING] Processing all regions: {:.3}s", elapsed_processing.as_secs_f64());
                eprintln!("[TIMING] Output writing: {:.3}s", 0.0f64);
                eprintln!("[TIMING] TOTAL execution: {:.3}s", elapsed_total.as_secs_f64());
            }
        }
    }

    Ok(())
}

struct RegionLoadResult {
    regions: Vec<Region>,
    amplicon_based_calling: Option<String>,
    amplicon_region_groups: Option<Vec<Vec<Region>>>,
}

struct ParsedBedResult {
    regions: Vec<Region>,
    amplicon_based_calling: Option<String>,
    amplicon_region_groups: Option<Vec<Vec<Region>>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ExecutionMode {
    Simple,
    Amplicon,
}

fn resolve_execution_mode(args: &Args, amplicon_based_calling: &Option<String>) -> ExecutionMode {
    if args.region.is_some() {
        ExecutionMode::Simple
    } else if amplicon_based_calling.is_some() {
        ExecutionMode::Amplicon
    } else {
        ExecutionMode::Simple
    }
}

fn select_region_batches_for_execution(
    execution_mode: ExecutionMode,
    regions: Vec<Region>,
    amplicon_region_groups: Option<Vec<Vec<Region>>>,
) -> Vec<Vec<Region>> {
    match execution_mode {
        ExecutionMode::Simple => vec![regions],
        ExecutionMode::Amplicon => {
            if let Some(groups) = amplicon_region_groups {
                if groups.is_empty() {
                    vec![regions]
                } else {
                    groups
                }
            } else {
                vec![regions]
            }
        }
    }
}

/// Parse regions from command line arguments
fn get_regions(args: &Args) -> Result<RegionLoadResult> {
    let mut regions = Vec::new();
    let bam_targets = BamReader::open(&args.bam)
        .context("Failed to open BAM for region normalization")?
        .target_names();

    // Check for -R option first
    if let Some(ref region_str) = args.region {
        let region_path = PathBuf::from(region_str);
        if region_path.exists() {
            let parsed = parse_bed_file(&region_path, args, Some(&bam_targets))?;
            return Ok(RegionLoadResult {
                regions: parsed.regions,
                amplicon_based_calling: parsed.amplicon_based_calling,
                amplicon_region_groups: parsed.amplicon_region_groups,
            });
        }
        let region = parse_region_string(
            region_str,
            args.zero_based.unwrap_or(0) == 1,
            Some(&bam_targets),
        )?;
        regions.push(region);
        return Ok(RegionLoadResult {
            regions,
            amplicon_based_calling: None,
            amplicon_region_groups: None,
        });
    }

    // Check for BED file
    if let Some(ref bed_path) = args.bed_file {
        if !bed_path.exists() {
            return Err(anyhow!("BED file not found: {:?}", bed_path));
        }
        let parsed = parse_bed_file(bed_path, args, Some(&bam_targets))?;
        return Ok(RegionLoadResult {
            regions: parsed.regions,
            amplicon_based_calling: parsed.amplicon_based_calling,
            amplicon_region_groups: parsed.amplicon_region_groups,
        });
    }

    Ok(RegionLoadResult {
        regions,
        amplicon_based_calling: None,
        amplicon_region_groups: None,
    })
}

/// Parse a region string like "chr1:1000-2000" or "chr1:1000"
fn parse_region_string(
    s: &str,
    zero_based: bool,
    bam_targets: Option<&[String]>,
) -> Result<Region> {
    // Format: chr:start-end or chr:start
    let parts: Vec<&str> = s.split(':').collect();
    if parts.len() != 2 {
        return Err(anyhow!("Invalid region format: {}. Expected chr:start-end", s));
    }

    let chr = normalize_region_chrom(parts[0], bam_targets);
    let pos_parts: Vec<&str> = parts[1].split('-').collect();

    let (start, end) = match pos_parts.len() {
        1 => {
            let pos: usize = pos_parts[0].parse()
                .context("Invalid position in region")?;
            (pos, pos)
        }
        2 => {
            let start: usize = pos_parts[0].parse()
                .context("Invalid start position in region")?;
            let end: usize = pos_parts[1].parse()
                .context("Invalid end position in region")?;
            (start, end)
        }
        _ => return Err(anyhow!("Invalid region format: {}", s)),
    };

    // Convert to 1-based if zero-based input
    let (start, end) = if zero_based {
        (start + 1, end)
    } else {
        (start, end)
    };

    // Java uses the chromosome as the gene field for -R regions.
    Ok(Region::new(chr.clone(), start, end, chr))
}

/// Parse a BED file and return regions
fn parse_bed_file(
    path: &PathBuf,
    args: &Args,
    bam_targets: Option<&[String]>,
) -> Result<ParsedBedResult> {
    let file = File::open(path).context("Failed to open BED file")?;
    let reader = BufReader::new(file);
    let mut bed_lines = Vec::new();
    let mut amplicon_parameters = args.amplicon_based_calling.clone();
    let mut zero_based = args.zero_based.map(|value| value == 1);

    for (_line_num, line) in reader.lines().enumerate() {
        let line = line.context("Failed to read BED line")?;
        let line = line.trim();

        // Skip empty lines and comments
        if line.is_empty() || line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
            continue;
        }

        if amplicon_parameters.is_none() {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() == 8 {
                let col6 = fields[6].parse::<i32>();
                let col7 = fields[7].parse::<i32>();
                if col6.is_ok() && col7.is_ok() {
                    let start_region = fields[1].parse::<i32>().map_err(|error| {
                        anyhow!(
                            "Incorrect format of BED file for amplicon mode. It must be 8 columns and 2, 3, 7 and 8 columns must contain region and amplicon starts and ends. {}",
                            error
                        )
                    })?;
                    let end_region = fields[2].parse::<i32>().map_err(|error| {
                        anyhow!(
                            "Incorrect format of BED file for amplicon mode. It must be 8 columns and 2, 3, 7 and 8 columns must contain region and amplicon starts and ends. {}",
                            error
                        )
                    })?;
                    let start_amplicon = col6.unwrap();
                    let end_amplicon = col7.unwrap();
                    if start_amplicon >= start_region && end_amplicon <= end_region {
                        amplicon_parameters = Some(DEFAULT_AMPLICON_PARAMETERS.to_string());
                        if zero_based.is_none() {
                            zero_based = Some(true);
                        }
                    }
                }
            }
        }

        bed_lines.push(line.to_string());
    }

    let use_zero_based = zero_based.unwrap_or(false);
    let (regions, amplicon_region_groups) = if amplicon_parameters.is_some() {
        let region_groups = parse_amplicon_region_groups(&bed_lines, bam_targets, use_zero_based)?;
        let flattened_regions = region_groups
            .iter()
            .flat_map(|group| group.iter().cloned())
            .collect();
        (flattened_regions, Some(region_groups))
    } else {
        (
            parse_standard_regions(&bed_lines, args, bam_targets, use_zero_based)?,
            None,
        )
    };

    Ok(ParsedBedResult {
        regions,
        amplicon_based_calling: amplicon_parameters,
        amplicon_region_groups,
    })
}

fn parse_standard_regions(
    bed_lines: &[String],
    args: &Args,
    bam_targets: Option<&[String]>,
    zero_based: bool,
) -> Result<Vec<Region>> {
    let mut regions = Vec::new();

    for (line_num, line) in bed_lines.iter().enumerate() {
        let fields: Vec<&str> = line.split('\t').collect();

        // Get chromosome (1-indexed column number to 0-indexed)
        let chr_idx = args.col_chr.saturating_sub(1);
        let start_idx = args.col_start.saturating_sub(1);
        let end_idx = args.col_end.saturating_sub(1);
        let gene_idx = args.col_gene.saturating_sub(1);

        if fields.len() <= chr_idx || fields.len() <= start_idx || fields.len() <= end_idx {
            if args.debug {
                eprintln!("Skipping malformed BED line {}: {}", line_num + 1, line);
            }
            continue;
        }

        let chr = normalize_region_chrom(fields[chr_idx], bam_targets);
        let start: usize = fields[start_idx].parse()
            .with_context(|| format!("Invalid start at line {}", line_num + 1))?;
        let end: usize = fields[end_idx].parse()
            .with_context(|| format!("Invalid end at line {}", line_num + 1))?;
        
        let gene = if fields.len() > gene_idx {
            fields[gene_idx].to_string()
        } else {
            String::new()
        };

        // BED is 0-based, half-open; convert to 1-based inclusive
        let (start, end) = if zero_based {
            (start + 1, end) // BED format: 0-based start, end is exclusive
        } else {
            (start, end)
        };

        regions.push(Region::new(chr, start, end, gene));
    }

    Ok(regions)
}

fn parse_amplicon_region_groups(
    bed_lines: &[String],
    bam_targets: Option<&[String]>,
    zero_based: bool,
) -> Result<Vec<Vec<Region>>> {
    let mut chromosome_order = Vec::new();
    let mut regions_by_chrom: std::collections::HashMap<String, Vec<(usize, usize, Region)>> =
        std::collections::HashMap::new();

    for (line_num, line) in bed_lines.iter().enumerate() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 8 {
            return Err(anyhow!(
                "Incorrect format of BED file for amplicon mode at line {}: expected at least 8 columns",
                line_num + 1
            ));
        }

        let chr = normalize_region_chrom(fields[0], bam_targets);
        let mut start: usize = fields[1]
            .parse()
            .with_context(|| format!("Invalid start at line {}", line_num + 1))?;
        let end: usize = fields[2]
            .parse()
            .with_context(|| format!("Invalid end at line {}", line_num + 1))?;
        let gene = fields[3].to_string();
        let mut insert_start: usize = fields[6]
            .parse()
            .with_context(|| format!("Invalid amplicon start at line {}", line_num + 1))?;
        let insert_end: usize = fields[7]
            .parse()
            .with_context(|| format!("Invalid amplicon end at line {}", line_num + 1))?;

        if zero_based && start < end {
            start += 1;
            insert_start += 1;
        }

        if !regions_by_chrom.contains_key(&chr) {
            chromosome_order.push(chr.clone());
        }

        let region = Region::new_with_insert(
            chr,
            start,
            end,
            gene,
            insert_start,
            insert_end,
        );

        let chrom_key = region.chr().to_string();
        regions_by_chrom
            .entry(chrom_key)
            .or_default()
            .push((insert_start, insert_end, region));
    }

    let mut region_groups = Vec::new();
    let mut previous_chr: Option<String> = None;
    let mut previous_end: Option<usize> = None;

    for chrom in chromosome_order {
        if let Some(chr_regions) = regions_by_chrom.get_mut(&chrom) {
            chr_regions.sort_by_key(|(insert_start, _, _)| *insert_start);

            for (insert_start, insert_end, region) in chr_regions.iter() {
                let starts_new_group = match (&previous_chr, previous_end) {
                    (Some(prev_chr), Some(prev_end)) => {
                        region.chr() != prev_chr || *insert_start > prev_end
                    }
                    _ => true,
                };

                if starts_new_group {
                    region_groups.push(Vec::new());
                }

                if let Some(current_group) = region_groups.last_mut() {
                    current_group.push(region.clone());
                }

                previous_chr = Some(region.chr().to_string());
                previous_end = Some(*insert_end);
            }
        }
    }

    Ok(region_groups)
}

fn normalize_region_chrom(chrom: &str, bam_targets: Option<&[String]>) -> String {
    let Some(targets) = bam_targets else {
        return chrom.to_string();
    };

    if targets.iter().any(|name| name == chrom) {
        return chrom.to_string();
    }

    if let Some(stripped) = chrom.strip_prefix("chr") {
        if targets.iter().any(|name| name == stripped) {
            return stripped.to_string();
        }
    }

    let with_chr = format!("chr{}", chrom);
    if targets.iter().any(|name| name == &with_chr) {
        return with_chr;
    }

    chrom.to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_region_string_full() {
        let region = parse_region_string("chr1:1000-2000", false, None).unwrap();
        assert_eq!(region.chr(), "chr1");
        assert_eq!(region.start(), 1000);
        assert_eq!(region.end(), 2000);
    }

    #[test]
    fn test_parse_region_string_single_position() {
        let region = parse_region_string("chr1:1000", false, None).unwrap();
        assert_eq!(region.chr(), "chr1");
        assert_eq!(region.start(), 1000);
        assert_eq!(region.end(), 1000);
    }

    #[test]
    fn test_parse_region_string_zero_based() {
        let region = parse_region_string("chr1:999-2000", true, None).unwrap();
        assert_eq!(region.chr(), "chr1");
        assert_eq!(region.start(), 1000); // 999 + 1
        assert_eq!(region.end(), 2000);
    }

    #[test]
    fn test_parse_region_string_invalid() {
        assert!(parse_region_string("invalid", false, None).is_err());
        assert!(parse_region_string("chr1", false, None).is_err());
    }

    #[test]
    fn test_parse_args_amplicon_based_calling() {
        let args = Args::parse_from([
            "vardict",
            "-G",
            "reference.fa",
            "-b",
            "input.bam",
            "-R",
            "chr1:1-10",
            "-a",
            "10:0.95",
        ]);

        assert_eq!(args.amplicon_based_calling.as_deref(), Some("10:0.95"));
    }

    #[test]
    fn test_parse_args_amplicon_based_calling_absent_by_default() {
        let args = Args::parse_from([
            "vardict",
            "-G",
            "reference.fa",
            "-b",
            "input.bam",
            "-R",
            "chr1:1-10",
        ]);

        assert!(args.amplicon_based_calling.is_none());
    }

    #[test]
    fn test_parse_args_zero_based_option_values() {
        let args_zero = Args::parse_from([
            "vardict",
            "-G",
            "reference.fa",
            "-b",
            "input.bam",
            "-R",
            "chr1:1-10",
            "-z",
            "1",
        ]);
        assert_eq!(args_zero.zero_based, Some(1));

        let args_one_based = Args::parse_from([
            "vardict",
            "-G",
            "reference.fa",
            "-b",
            "input.bam",
            "-R",
            "chr1:1-10",
            "-z",
            "0",
        ]);
        assert_eq!(args_one_based.zero_based, Some(0));
    }

    #[test]
    fn test_parse_bed_file_auto_detects_amplicon_and_zero_based_default() {
        let args = Args::parse_from([
            "vardict",
            "-G",
            "reference.fa",
            "-b",
            "input.bam",
        ]);

        let bed_path = PathBuf::from(format!(
            "tmp/test_amplicon_auto_{}_{}.bed",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all("tmp").unwrap();
        std::fs::write(&bed_path, "chr1\t100\t200\tGENE\t.\t.\t120\t180\n").unwrap();

        let parsed = parse_bed_file(&bed_path, &args, None).unwrap();
        assert_eq!(
            parsed.amplicon_based_calling.as_deref(),
            Some(DEFAULT_AMPLICON_PARAMETERS)
        );
        assert_eq!(parsed.regions.len(), 1);
        assert_eq!(parsed.regions[0].start(), 101);
        assert_eq!(parsed.amplicon_region_groups.as_ref().map(Vec::len), Some(1));

        let _ = std::fs::remove_file(&bed_path);
    }

    #[test]
    fn test_parse_bed_file_amplicon_groups_by_insert_overlap() {
        let args = Args::parse_from([
            "vardict",
            "-G",
            "reference.fa",
            "-b",
            "input.bam",
            "-a",
            "10:0.95",
        ]);

        let bed_path = PathBuf::from(format!(
            "tmp/test_amplicon_groups_{}_{}.bed",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        std::fs::create_dir_all("tmp").unwrap();
        std::fs::write(
            &bed_path,
            concat!(
                "chr1\t100\t200\tG1\t.\t.\t120\t140\n",
                "chr1\t150\t250\tG2\t.\t.\t135\t160\n",
                "chr1\t260\t320\tG3\t.\t.\t200\t210\n"
            ),
        )
        .unwrap();

        let parsed = parse_bed_file(&bed_path, &args, None).unwrap();
        let groups = parsed.amplicon_region_groups.expect("expected amplicon groups");

        assert_eq!(groups.len(), 2);
        assert_eq!(groups[0].len(), 2);
        assert_eq!(groups[1].len(), 1);
        assert_eq!(groups[0][0].gene(), "G1");
        assert_eq!(groups[0][1].gene(), "G2");
        assert_eq!(groups[1][0].gene(), "G3");

        let _ = std::fs::remove_file(&bed_path);
    }

    #[test]
    fn test_get_regions_region_option_ignores_amplicon_setting() {
        let args = Args::parse_from([
            "vardict",
            "-G",
            "reference.fa",
            "-b",
            "test_data/test_168714.bam",
            "-R",
            "chr20:168700-168710",
            "-a",
            "10:0.95",
        ]);

        let loaded = get_regions(&args).unwrap();
        assert_eq!(loaded.amplicon_based_calling, None);
        assert_eq!(loaded.regions.len(), 1);
    }

    #[test]
    fn test_resolve_execution_mode_region_forces_simple() {
        let args = Args::parse_from([
            "vardict",
            "-G",
            "reference.fa",
            "-b",
            "input.bam",
            "-R",
            "chr1:1-10",
            "-a",
            "10:0.95",
        ]);

        let mode = resolve_execution_mode(&args, &Some("10:0.95".to_string()));
        assert_eq!(mode, ExecutionMode::Simple);
    }

    #[test]
    fn test_resolve_execution_mode_amplicon_without_region() {
        let args = Args::parse_from([
            "vardict",
            "-G",
            "reference.fa",
            "-b",
            "input.bam",
        ]);

        let mode = resolve_execution_mode(&args, &Some("10:0.95".to_string()));
        assert_eq!(mode, ExecutionMode::Amplicon);
    }

    #[test]
    fn test_select_region_batches_for_execution_simple_mode_uses_flat_regions() {
        let regions = vec![
            Region::new("chr1".to_string(), 10, 20, "G1".to_string()),
            Region::new("chr1".to_string(), 30, 40, "G2".to_string()),
        ];

        let batches = select_region_batches_for_execution(ExecutionMode::Simple, regions, None);
        assert_eq!(batches.len(), 1);
        assert_eq!(batches[0].len(), 2);
    }

    #[test]
    fn test_select_region_batches_for_execution_amplicon_mode_uses_groups() {
        let regions = vec![
            Region::new("chr1".to_string(), 10, 20, "G1".to_string()),
            Region::new("chr1".to_string(), 30, 40, "G2".to_string()),
            Region::new("chr1".to_string(), 50, 60, "G3".to_string()),
        ];
        let groups = vec![
            vec![Region::new("chr1".to_string(), 10, 20, "G1".to_string())],
            vec![
                Region::new("chr1".to_string(), 30, 40, "G2".to_string()),
                Region::new("chr1".to_string(), 50, 60, "G3".to_string()),
            ],
        ];

        let batches = select_region_batches_for_execution(
            ExecutionMode::Amplicon,
            regions,
            Some(groups),
        );
        assert_eq!(batches.len(), 2);
        assert_eq!(batches[0].len(), 1);
        assert_eq!(batches[1].len(), 2);
    }
}
