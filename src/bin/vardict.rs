//! VarDict-rs: A Rust implementation of VarDict variant caller (Simple Mode)
//!
//! Usage: vardict -G <reference.fa> -b <input.bam> [options] <region or BED file>

#[cfg(feature = "mimalloc-global")]
use mimalloc::MiMalloc;

#[cfg(feature = "mimalloc-global")]
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

#[cfg(feature = "jemalloc-global")]
use tikv_jemallocator::Jemalloc;

#[cfg(feature = "jemalloc-global")]
#[global_allocator]
static GLOBAL_JE: Jemalloc = Jemalloc;

#[cfg(feature = "dhat-heap")]
#[global_allocator]
static ALLOC: dhat::Alloc = dhat::Alloc;

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::sync::LazyLock;

use anyhow::{Context, Result, anyhow};
use clap::Parser;
use regex::Regex;

use crackle_kit::tracing::level_filters::LevelFilter;
use crackle_kit::tracing::{Level, event};
use crackle_kit::tracing_kit::setup_logging_stderr_only;
use vardict_rs::data::bam_reader::BamReader;
use vardict_rs::data::region::Region;
use vardict_rs::mods::pipeline::PipelineConfig;
use vardict_rs::mods::vardict_pipeline::{SomaticCombineLookupResult, VarDictPipeline};

const DEFAULT_AMPLICON_PARAMETERS: &str = "10:0.95";

static SAMPLE_PATTERN: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"([^/\._]+)\.sorted[^/]*\.bam$").expect("valid regex"));
static SAMPLE_PATTERN2: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"([^/]+)[_\.][^/]*bam$").expect("valid regex"));

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

    /// Indexed BAM file, or paired BAMs for somatic mode: tumor.bam|normal.bam
    #[arg(short = 'b', long = "bam", required = true)]
    bam: String,

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

    /// Move indels to 3' end (Java: -3)
    #[arg(short = '3', long = "move-indels-3prime")]
    move_indels_3prime: bool,

    /// The read position filter (Java: -P)
    #[arg(short = 'P', long = "read-pos-filter", default_value = "5.0")]
    read_pos_filter: f64,

    /// The Qratio of good quality reads to bad quality reads (Java: -o)
    #[arg(short = 'o', long = "qratio", default_value = "1.5")]
    qratio: f64,

    /// Minimum match count threshold (Java: -M)
    #[arg(short = 'M', long = "min-match", default_value = "0")]
    min_match: i32,

    /// Trim bases after this position (Java: -T)
    #[arg(short = 'T', long = "trim-bases-after", default_value = "0")]
    trim_bases_after: i32,

    /// Minimum reads per strand to avoid strand bias (Java: -B)
    #[arg(short = 'B', long = "min-bias-reads", default_value = "2")]
    min_bias_reads: usize,

    /// Turn off structural variant calling
    #[arg(short = 'U', long = "nosv")]
    no_sv: bool,

    /// Keep only one read from overlapping pairs based on alignment position (Java: -u)
    #[arg(short = 'u', long = "unique-overlap")]
    unique_mode_alignment: bool,

    /// Keep only one read from overlapping pairs based on second-in-pair flag (Java: -UN)
    #[arg(long = "unique-second-in-pair", visible_alias = "UN")]
    unique_mode_second_in_pair: bool,

    /// Turn on deleting duplicate variants (Java: --deldupvar)
    #[arg(long = "deldupvar")]
    delete_duplicate_variants: bool,

    /// Amplicon mode parameters (Java: -a), e.g. "10:0.95"
    #[arg(short = 'a', long = "amplicon")]
    amplicon_based_calling: Option<String>,

    /// Debug mode - print additional information
    #[arg(short = 'D', long = "debug")]
    debug: bool,

    /// Experimental feature: compute Fisher exact test columns in Java parity mode
    #[arg(long = "fisher")]
    fisher: bool,

    /// Mean insert size (Java: -w)
    #[arg(short = 'w', long = "insert-size", default_value = "300")]
    insert_size: i32,

    /// Insert size standard deviation (Java: -W)
    #[arg(short = 'W', long = "insert-std", default_value = "100")]
    insert_std: i32,

    /// Number of insert size standard deviations for discordant filtering (Java: -A)
    #[arg(short = 'A', long = "insert-std-amt", default_value = "4")]
    insert_std_amt: i32,

    /// Minimum structural variant length (Java: -L)
    #[arg(short = 'L', long = "sv-min-len", default_value = "1000")]
    sv_min_len: usize,

    /// Enable chimeric read filtering
    #[arg(long = "chimeric")]
    chimeric: bool,

    /// CRISPR cutting site position (Java: -J / --crispr)
    #[arg(short = 'J', long = "crispr", default_value = "0")]
    crispr_cutting_site: i32,

    /// CRISPR filtering bp overlap (Java: -j)
    #[arg(short = 'j', default_value = "0")]
    crispr_filtering_bp: i32,

    /// Output all variants including reference calls (pileup mode)
    #[arg(short = 'p', long = "pileup")]
    pileup: bool,

    /// Indicate whether coordinates are zero-based: 1 for zero-based, 0 for one-based
    #[arg(short = 'z', long = "zero", value_parser = clap::value_parser!(u8).range(0..=1))]
    zero_based: Option<u8>,

    /// Count N bases in total depth (Java: -K)
    #[arg(short = 'K', long = "include-n")]
    include_n: bool,

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

    /// Output splicing read counts only (Java: -i)
    #[arg(short = 'i', long = "splice")]
    output_splicing: bool,

    /// Log level
    #[arg(long, default_value_t = LevelFilter::WARN)]
    log_level: LevelFilter,
}

fn normalize_legacy_cli_args(raw_args: Vec<OsString>) -> Vec<OsString> {
    raw_args
        .into_iter()
        .map(|arg| {
            if arg == "-UN" {
                OsString::from("--unique-second-in-pair")
            } else if arg == "-fisher" {
                OsString::from("--fisher")
            } else {
                arg
            }
        })
        .collect()
}

fn parse_args() -> Args {
    let raw_args = std::env::args_os().collect::<Vec<_>>();
    let normalized_args = normalize_legacy_cli_args(raw_args);
    Args::parse_from(normalized_args)
}

fn main() -> Result<()> {
    #[cfg(feature = "dhat-heap")]
    let _profiler = dhat::Profiler::new_heap();

    let args = parse_args();

    setup_logging_stderr_only(args.log_level)?;

    let bam_inputs = parse_bam_inputs(&args.bam)?;

    // Validate input files exist
    if !args.reference.exists() {
        return Err(anyhow!("Reference file not found: {:?}", args.reference));
    }
    validate_bam_with_index(&bam_inputs.primary_bam)?;
    if let Some(ref secondary_bam) = bam_inputs.secondary_bam {
        validate_bam_with_index(secondary_bam)?;
    }

    // Check for FASTA index
    let fai_path = args.reference.with_extension("fa.fai");
    let fai_path2 = {
        let mut p = args.reference.clone();
        let reference_name = args
            .reference
            .file_name()
            .ok_or_else(|| anyhow!("Invalid reference path: {:?}", args.reference))?;
        p.set_file_name(format!("{}.fai", reference_name.to_string_lossy()));
        p
    };
    if !fai_path.exists() && !fai_path2.exists() {
        return Err(anyhow!(
            "Reference index (.fai) not found for {:?}. Please generate the FASTA index with samtools and retry.",
            args.reference
        ));
    }

    // Determine sample name from first BAM unless overridden
    let sample_name = args
        .sample_name
        .clone()
        .unwrap_or_else(|| infer_sample_name_from_bam(&bam_inputs.primary_bam));

    // Get regions to process
    let region_load = get_regions(&args, &bam_inputs.primary_bam)?;
    if region_load.regions.is_empty() {
        return Err(anyhow!(
            "No regions specified. Provide -R option or a BED file."
        ));
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

    let execution_mode = resolve_execution_mode(
        &args,
        &region_load.amplicon_based_calling,
        bam_inputs.secondary_bam.is_some(),
    );

    // Print header if requested
    if args.print_header {
        match execution_mode {
            ExecutionMode::Amplicon => {
                println!(
                    "{}",
                    vardict_rs::mods::output_variant::get_amplicon_header_line()
                );
            }
            ExecutionMode::Somatic => {
                println!(
                    "{}",
                    vardict_rs::mods::output_variant::get_somatic_header_line()
                );
            }
            ExecutionMode::Splicing => {
                println!("Sample\tChr\tIntron\tIntron count");
            }
            ExecutionMode::Simple => {
                println!(
                    "{}",
                    vardict_rs::mods::output_variant::get_simple_header_line(
                        args.crispr_cutting_site != 0,
                    )
                );
            }
        }
    }

    let bam_paths = bam_inputs.to_paths();

    // Always use SharedReference (loaded into memory for fast access)
    run_variant_calling(
        &args,
        config,
        bam_paths,
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
    bam_paths: Vec<PathBuf>,
    regions: Vec<Region>,
    amplicon_based_calling: Option<String>,
    amplicon_region_groups: Option<Vec<Vec<Region>>>,
    execution_mode: ExecutionMode,
) -> Result<()> {
    use std::collections::HashSet;
    use std::sync::Arc;
    use std::time::Instant;
    use vardict_rs::conf::Configuration;
    use vardict_rs::data::shared_reference::load_shared_reference_chroms;
    use vardict_rs::mods::parallel_pipeline::{
        OrderedStreamConsumer, ParallelPipeline, RegionResult,
    };
    use vardict_rs::scopedata::global_read_only_scope::{GlobalReadOnlyScope, INSTANCE};

    if bam_paths.is_empty() {
        return Err(anyhow!("No BAM paths available for execution"));
    }
    let bam_paths_string = bam_paths
        .iter()
        .map(|path| path.to_string_lossy().to_string())
        .collect::<Vec<_>>();

    let start_total = Instant::now();
    let num_threads = args.num_threads.max(1);

    if args.debug {
        event!(Level::INFO, "Loading reference genome into memory...");
    }

    // Get unique chromosomes from regions
    let start_ref_load = Instant::now();
    let chroms: std::collections::HashSet<&str> = regions.iter().map(|r| r.chr()).collect();
    let chrom_vec: Vec<&str> = chroms.into_iter().collect();
    let reference_path = args
        .reference
        .to_str()
        .ok_or_else(|| anyhow!("Reference path is not valid UTF-8: {:?}", args.reference))?;

    // Load only the needed chromosomes for efficiency
    let reference = load_shared_reference_chroms(reference_path, &chrom_vec)
        .context("Failed to load reference genome")?;

    let elapsed_ref_load = start_ref_load.elapsed();

    // Initialize GlobalReadOnlyScope (required by VarDictPipeline)
    // Must be done AFTER loading reference to populate chr_lens
    let mut conf = Configuration::default();
    let sam_filter = if let Some(hex) = args
        .sam_filter
        .strip_prefix("0x")
        .or_else(|| args.sam_filter.strip_prefix("0X"))
    {
        u32::from_str_radix(hex, 16).context("Failed to parse sam_filter as hex")?
    } else {
        args.sam_filter
            .parse::<u32>()
            .context("Failed to parse sam_filter as decimal")?
    };
    conf.goodq = args.min_base_quality;
    conf.freq = if args.pileup {
        -1.0
    } else {
        args.min_frequency
    };
    conf.minr = if args.pileup {
        0
    } else {
        args.min_variant_reads
    };
    conf.vext = args.vext;
    conf.mismatch = args.mismatch;
    conf.min_match = args.min_match;
    conf.trim_bases_after = args.trim_bases_after;
    conf.read_pos_filter = args.read_pos_filter;
    conf.qratio = args.qratio;
    conf.min_bias_reads = args.min_bias_reads;
    conf.mapping_quality = if args.min_mapping_quality > 0 {
        Some(args.min_mapping_quality)
    } else {
        None
    };
    conf.sam_filter = sam_filter;
    conf.downsampling = args.downsampling;
    conf.move_indels_to_3 = args.move_indels_3prime;
    conf.chimeric_filter = args.chimeric;
    conf.remove_duplicated_reads = args.remove_duplicates;
    conf.disable_sv = args.no_sv;
    conf.unique_mode_alignment_enabled = args.unique_mode_alignment;
    conf.unique_mode_second_in_pair_enabled = args.unique_mode_second_in_pair;
    conf.delete_duplicate_variants = args.delete_duplicate_variants;
    conf.debug = args.debug;
    conf.fisher = args.fisher;
    conf.inssize = args.insert_size;
    conf.insstd = args.insert_std;
    conf.insstdamt = args.insert_std_amt;
    conf.sv_min_len = args.sv_min_len;
    conf.crispr_cutting_site = args.crispr_cutting_site;
    conf.crispr_filtering_bp = args.crispr_filtering_bp;
    conf.include_n_in_total_depth = args.include_n;
    conf.amplicon_based_calling = amplicon_based_calling.clone();
    conf.perform_local_realignment = args.local_realignment == 1;
    conf.number_nucleotide_to_extend = args.number_nucleotide_to_extend;
    conf.reference_extension = args.reference_extension;
    let mut scope = GlobalReadOnlyScope::default();
    scope.amplicon_based_calling = conf.amplicon_based_calling.clone();
    scope.conf = conf;
    scope.chr_lens = reference.get_chromosome_lengths();
    scope.bam_paths = bam_paths_string.clone();
    let _ = INSTANCE.set(scope);

    let region_batches = select_region_batches_for_execution(
        execution_mode,
        regions,
        amplicon_region_groups,
        num_threads,
    );

    if args.debug {
        event!(
            Level::INFO,
            "[TIMING] Reference loading: {:.3}s",
            elapsed_ref_load.as_secs_f64()
        );
        event!(
            Level::INFO,
            "Loaded {} chromosome(s), {:.2} MB total",
            reference.num_chromosomes(),
            reference.total_size() as f64 / 1_048_576.0
        );
        event!(Level::INFO, "Execution mode: {:?}", execution_mode);
        event!(Level::INFO, "Execution batches: {}", region_batches.len());
        if num_threads > 1 {
            event!(
                Level::INFO,
                "Processing {} regions with {} threads...",
                region_batches.iter().map(Vec::len).sum::<usize>(),
                num_threads
            );
        } else {
            event!(
                Level::INFO,
                "Processing {} regions...",
                region_batches.iter().map(Vec::len).sum::<usize>()
            );
        }
    }

    // Process regions
    let start_processing = Instant::now();
    let primary_bam_path = bam_paths_string[0].clone();

    match execution_mode {
        ExecutionMode::Simple => {
            let pipeline = ParallelPipeline::new(reference, config, num_threads);
            let all_regions: Vec<Region> = region_batches.into_iter().flatten().collect();

            if num_threads > 1 {
                let (sender, receiver) = crossbeam_channel::bounded::<(usize, RegionResult)>(10);
                let debug = args.debug;

                let consumer_handle =
                    std::thread::spawn(move || OrderedStreamConsumer::new(receiver, debug).run());

                pipeline.process_regions_vardict_streaming(primary_bam_path, all_regions, sender);

                consumer_handle.join().expect("consumer thread panicked")?;
            } else {
                let mut stdout = io::stdout().lock();
                pipeline.process_regions_vardict_direct_write(
                    primary_bam_path,
                    all_regions,
                    &mut stdout,
                    args.debug,
                )?;
            }

            let elapsed_processing = start_processing.elapsed();

            let elapsed_total = start_total.elapsed();

            if args.debug {
                event!(
                    Level::INFO,
                    "[TIMING] Processing all regions: {:.3}s",
                    elapsed_processing.as_secs_f64()
                );
                event!(Level::INFO, "[TIMING] Output writing: {:.3}s", 0.0f64);
                event!(
                    Level::INFO,
                    "[TIMING] TOTAL execution: {:.3}s",
                    elapsed_total.as_secs_f64()
                );
            }
        }
        ExecutionMode::Amplicon => {
            let vardict_pipeline = VarDictPipeline::new(&config.sample_name)
                .with_min_frequency(config.min_frequency)
                .with_min_base_quality(config.quality_threshold)
                .with_min_mapping_quality(config.mapq_threshold);

            let global_scope = Arc::new(
                INSTANCE
                    .get()
                    .expect("GlobalReadOnlyScope not initialized")
                    .clone(),
            );
            let mut stdout = io::stdout().lock();

            for amplicon_group in region_batches {
                if amplicon_group.is_empty() {
                    continue;
                }

                let mut vars_per_amplicon = Vec::with_capacity(amplicon_group.len());
                let mut splice: HashSet<String> = HashSet::new();

                for region in &amplicon_group {
                    let mut bam_reader = BamReader::open(&primary_bam_path)?;
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
                event!(
                    Level::INFO,
                    "[TIMING] Processing all regions: {:.3}s",
                    elapsed_processing.as_secs_f64()
                );
                event!(Level::INFO, "[TIMING] Output writing: {:.3}s", 0.0f64);
                event!(
                    Level::INFO,
                    "[TIMING] TOTAL execution: {:.3}s",
                    elapsed_total.as_secs_f64()
                );
            }
        }
        ExecutionMode::Somatic => {
            if bam_paths_string.len() < 2 {
                return Err(anyhow!(
                    "Somatic execution requires paired BAM input in the format tumor.bam|normal.bam"
                ));
            }

            let tumor_bam_path = bam_paths_string[0].clone();
            let normal_bam_path = bam_paths_string[1].clone();

            let vardict_pipeline = VarDictPipeline::new(&config.sample_name)
                .with_min_frequency(config.min_frequency)
                .with_min_base_quality(config.quality_threshold)
                .with_min_mapping_quality(config.mapq_threshold)
                .with_pileup(config.pileup);

            let global_scope = Arc::new(
                INSTANCE
                    .get()
                    .expect("GlobalReadOnlyScope not initialized")
                    .clone(),
            );
            let mut stdout = io::stdout().lock();

            for region_group in region_batches {
                for region in &region_group {
                    let mut splice: HashSet<String> = HashSet::new();
                    let tumor_bam_paths = vec![tumor_bam_path.clone()];
                    let normal_bam_paths = vec![normal_bam_path.clone()];
                    let combined_bam_paths = vec![tumor_bam_path.clone(), normal_bam_path.clone()];

                    let mut tumor_bam_reader = BamReader::open(&tumor_bam_path)?;
                    let tumor_output = vardict_pipeline
                        .process_region_to_aligned_vars_from_bam_with_paths(
                            region,
                            &reference,
                            &mut tumor_bam_reader,
                            Arc::clone(&global_scope),
                            &tumor_bam_paths,
                        )?;
                    splice.extend(tumor_output.splice.iter().cloned());

                    let mut normal_bam_reader = BamReader::open(&normal_bam_path)?;
                    let normal_output = vardict_pipeline
                        .process_region_to_aligned_vars_from_bam_with_paths(
                            region,
                            &reference,
                            &mut normal_bam_reader,
                            Arc::clone(&global_scope),
                            &normal_bam_paths,
                        )?;
                    splice.extend(normal_output.splice.iter().cloned());

                    let combined_output = vardict_pipeline
                        .process_region_to_aligned_vars_from_bam_paths(
                            region,
                            &reference,
                            &combined_bam_paths,
                            Arc::clone(&global_scope),
                        )?;

                    let initial_max_read_length = tumor_output
                        .max_read_length
                        .max(normal_output.max_read_length);
                    let combined_max_read_length = combined_output.max_read_length;
                    let combined_aligned_variants = combined_output.aligned_vars.aligned_variants;

                    let combine_lookup =
                        move |_chr_name: &str,
                              position: i64,
                              description_string: &str,
                              max_read_length: usize| {
                            let combined_variant =
                                combined_aligned_variants.get(&position).and_then(|vars| {
                                    vars.variants
                                        .iter()
                                        .find(|variant| {
                                            variant.description_string == description_string
                                        })
                                        .cloned()
                                });

                            SomaticCombineLookupResult {
                                combined_variant,
                                max_read_length: combined_max_read_length.max(max_read_length),
                            }
                        };

                    let output_lines = vardict_pipeline
                        .run_somatic_post_processor_with_combine_lookup(
                            normal_output.aligned_vars,
                            tumor_output.aligned_vars,
                            region,
                            &splice,
                            initial_max_read_length,
                            Some(&combine_lookup),
                        );

                    for line in output_lines {
                        writeln!(stdout, "{}", line)?;
                    }
                }
            }

            let elapsed_processing = start_processing.elapsed();
            let elapsed_total = start_total.elapsed();

            if args.debug {
                event!(
                    Level::INFO,
                    "[TIMING] Processing all regions: {:.3}s",
                    elapsed_processing.as_secs_f64()
                );
                event!(Level::INFO, "[TIMING] Output writing: {:.3}s", 0.0f64);
                event!(
                    Level::INFO,
                    "[TIMING] TOTAL execution: {:.3}s",
                    elapsed_total.as_secs_f64()
                );
            }
        }
        ExecutionMode::Splicing => {
            let vardict_pipeline = VarDictPipeline::new(&config.sample_name)
                .with_min_frequency(config.min_frequency)
                .with_min_base_quality(config.quality_threshold)
                .with_min_mapping_quality(config.mapq_threshold)
                .with_pileup(config.pileup);

            let global_scope = Arc::new(
                INSTANCE
                    .get()
                    .expect("GlobalReadOnlyScope not initialized")
                    .clone(),
            );
            let mut stdout = io::stdout().lock();

            for region_group in region_batches {
                for region in &region_group {
                    let mut bam_reader = BamReader::open(&primary_bam_path)?;
                    let output_lines = vardict_pipeline.process_region_splicing_from_bam(
                        region,
                        &reference,
                        &mut bam_reader,
                        Arc::clone(&global_scope),
                    )?;

                    for line in output_lines {
                        writeln!(stdout, "{}", line)?;
                    }
                }
            }

            let elapsed_processing = start_processing.elapsed();
            let elapsed_total = start_total.elapsed();

            if args.debug {
                event!(
                    Level::INFO,
                    "[TIMING] Processing all regions: {:.3}s",
                    elapsed_processing.as_secs_f64()
                );
                event!(Level::INFO, "[TIMING] Output writing: {:.3}s", 0.0f64);
                event!(
                    Level::INFO,
                    "[TIMING] TOTAL execution: {:.3}s",
                    elapsed_total.as_secs_f64()
                );
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
    Somatic,
    Splicing,
}

fn resolve_execution_mode(
    args: &Args,
    amplicon_based_calling: &Option<String>,
    has_paired_bam: bool,
) -> ExecutionMode {
    if args.output_splicing {
        ExecutionMode::Splicing
    } else if args.region.is_none() && amplicon_based_calling.is_some() {
        ExecutionMode::Amplicon
    } else if has_paired_bam {
        ExecutionMode::Somatic
    } else {
        ExecutionMode::Simple
    }
}

fn select_region_batches_for_execution(
    execution_mode: ExecutionMode,
    regions: Vec<Region>,
    amplicon_region_groups: Option<Vec<Vec<Region>>>,
    num_threads: usize,
) -> Vec<Vec<Region>> {
    match execution_mode {
        ExecutionMode::Simple | ExecutionMode::Somatic | ExecutionMode::Splicing => {
            let batch_size = (num_threads * 4).max(8);
            regions
                .chunks(batch_size)
                .map(|chunk| chunk.to_vec())
                .collect()
        }
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

#[derive(Debug, Clone)]
struct BamInputs {
    primary_bam: PathBuf,
    secondary_bam: Option<PathBuf>,
}

impl BamInputs {
    fn to_paths(&self) -> Vec<PathBuf> {
        let mut paths = vec![self.primary_bam.clone()];
        if let Some(ref secondary_bam) = self.secondary_bam {
            paths.push(secondary_bam.clone());
        }
        paths
    }
}

fn parse_bam_inputs(raw_bam: &str) -> Result<BamInputs> {
    let parts = raw_bam
        .split('|')
        .map(str::trim)
        .filter(|part| !part.is_empty())
        .collect::<Vec<_>>();

    if parts.is_empty() {
        return Err(anyhow!(
            "BAM argument is empty. Provide -b <bam> or -b <tumor.bam|normal.bam>"
        ));
    }

    if parts.len() > 2 {
        return Err(anyhow!(
            "Invalid BAM argument: expected one BAM or two BAMs separated by '|', got {} entries",
            parts.len()
        ));
    }

    Ok(BamInputs {
        primary_bam: PathBuf::from(parts[0]),
        secondary_bam: parts.get(1).map(|value| PathBuf::from(*value)),
    })
}

fn validate_bam_with_index(path: &Path) -> Result<()> {
    if !path.exists() {
        return Err(anyhow!("BAM file not found: {:?}", path));
    }

    let bai_path = path.with_extension("bam.bai");
    let bai_path2 = {
        let mut p = path.to_path_buf();
        let file_name = path
            .file_name()
            .ok_or_else(|| anyhow!("Invalid BAM path: {:?}", path))?
            .to_string_lossy()
            .to_string();
        p.set_file_name(format!("{}.bai", file_name));
        p
    };
    if !bai_path.exists() && !bai_path2.exists() {
        return Err(anyhow!(
            "BAM index (.bai) not found. Please run: samtools index {:?}",
            path
        ));
    }

    Ok(())
}

fn infer_sample_name_from_bam(path: &Path) -> String {
    let bam_path = path.to_string_lossy();

    if let Some(captures) = SAMPLE_PATTERN.captures(&bam_path) {
        if let Some(sample) = captures.get(1) {
            return sample.as_str().to_string();
        }
    }

    if let Some(captures) = SAMPLE_PATTERN2.captures(&bam_path) {
        if let Some(sample) = captures.get(1) {
            return sample.as_str().to_string();
        }
    }

    path.file_stem()
        .and_then(|stem| stem.to_str())
        .map(|stem| stem.to_string())
        .unwrap_or_else(|| "SAMPLE".to_string())
}

/// Parse regions from command line arguments
fn get_regions(args: &Args, primary_bam_path: &Path) -> Result<RegionLoadResult> {
    let mut regions = Vec::new();
    let bam_targets = BamReader::open(primary_bam_path)
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
        regions.push(extend_region(region, args.number_nucleotide_to_extend));
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

fn extend_region(region: Region, number_nucleotide_to_extend: i32) -> Region {
    let x = number_nucleotide_to_extend.max(0) as usize;
    if x == 0 {
        return region;
    }

    let display_start = region.start() as i64 - x as i64;
    let clamped_start = region.start().saturating_sub(x);
    let extended_end = region.end().saturating_add(x);

    Region::new_extended(
        region.chr().to_string(),
        clamped_start,
        extended_end,
        region.gene().to_string(),
        display_start,
    )
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
        return Err(anyhow!(
            "Invalid region format: {}. Expected chr:start-end",
            s
        ));
    }

    let chr = normalize_region_chrom(parts[0], bam_targets);
    let pos_parts: Vec<&str> = parts[1].split('-').collect();

    let (start, end) = match pos_parts.len() {
        1 => {
            let pos: usize = pos_parts[0].parse().context("Invalid position in region")?;
            (pos, pos)
        }
        2 => {
            let start: usize = pos_parts[0]
                .parse()
                .context("Invalid start position in region")?;
            let end: usize = pos_parts[1]
                .parse()
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
        if line.is_empty()
            || line.starts_with('#')
            || line.starts_with("track")
            || line.starts_with("browser")
        {
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
                event!(
                    Level::WARN,
                    "Skipping malformed BED line {}: {}",
                    line_num + 1,
                    line
                );
            }
            continue;
        }

        let chr = normalize_region_chrom(fields[chr_idx], bam_targets);
        let start: usize = fields[start_idx]
            .parse()
            .with_context(|| format!("Invalid start at line {}", line_num + 1))?;
        let end: usize = fields[end_idx]
            .parse()
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

        regions.push(extend_region(
            Region::new(chr, start, end, gene),
            args.number_nucleotide_to_extend,
        ));
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

        let region = Region::new_with_insert(chr, start, end, gene, insert_start, insert_end);

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

    fn parse_args_for_test<I, T>(args: I) -> Args
    where
        I: IntoIterator<Item = T>,
        T: Into<OsString>,
    {
        let raw_args = args.into_iter().map(Into::into).collect::<Vec<_>>();
        let normalized_args = normalize_legacy_cli_args(raw_args);
        Args::parse_from(normalized_args)
    }

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
    fn test_infer_sample_name_from_bam_matches_java_secondary_pattern() {
        let sample = infer_sample_name_from_bam(Path::new(
            "/tmp/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam",
        ));

        assert_eq!(
            sample,
            "NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211"
        );
    }

    #[test]
    fn test_infer_sample_name_from_bam_matches_java_sorted_pattern() {
        let sample = infer_sample_name_from_bam(Path::new("/tmp/L861Q.sorted.bam"));

        assert_eq!(sample, "L861Q");
    }

    #[test]
    fn test_infer_sample_name_from_bam_falls_back_to_stem() {
        let sample = infer_sample_name_from_bam(Path::new("/tmp/no_delimiter_bamname.bamx"));

        assert_eq!(sample, "no_delimiter_bamname");
    }

    #[test]
    fn test_parse_region_string_invalid() {
        assert!(parse_region_string("invalid", false, None).is_err());
        assert!(parse_region_string("chr1", false, None).is_err());
    }

    #[test]
    fn test_extend_region_preserves_negative_display_start() {
        let region = Region::new("20".to_string(), 1, 1_000_000, "20".to_string());

        let extended = extend_region(region, 150);

        assert_eq!(extended.start(), 0);
        assert_eq!(extended.display_start(), -149);
        assert_eq!(extended.end(), 1_000_150);
    }

    #[test]
    fn test_parse_standard_regions_applies_x_extension() {
        let args = parse_args_for_test(["vardict", "-G", "ref.fa", "-b", "reads.bam", "-x", "150"]);
        let bed_lines = vec!["20\t0\t1000000\t20".to_string()];

        let regions = parse_standard_regions(&bed_lines, &args, None, true).unwrap();

        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].start(), 0);
        assert_eq!(regions[0].display_start(), -149);
        assert_eq!(regions[0].end(), 1_000_150);
    }

    #[test]
    fn test_parse_args_amplicon_based_calling() {
        let args = parse_args_for_test([
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
        let args = parse_args_for_test([
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
        let args_zero = parse_args_for_test([
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

        let args_one_based = parse_args_for_test([
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
        let args = parse_args_for_test(["vardict", "-G", "reference.fa", "-b", "input.bam"]);

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
        assert_eq!(
            parsed.amplicon_region_groups.as_ref().map(Vec::len),
            Some(1)
        );

        let _ = std::fs::remove_file(&bed_path);
    }

    #[test]
    fn test_parse_bed_file_amplicon_groups_by_insert_overlap() {
        let args = parse_args_for_test([
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
        let groups = parsed
            .amplicon_region_groups
            .expect("expected amplicon groups");

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
        let args = parse_args_for_test([
            "vardict",
            "-G",
            "reference.fa",
            "-b",
            "testdata/test_168714.bam",
            "-R",
            "chr20:168700-168710",
            "-a",
            "10:0.95",
        ]);

        let loaded = get_regions(&args, Path::new("testdata/test_168714.bam")).unwrap();
        assert_eq!(loaded.amplicon_based_calling, None);
        assert_eq!(loaded.regions.len(), 1);
    }

    #[test]
    fn test_resolve_execution_mode_region_forces_simple() {
        let args = parse_args_for_test([
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

        let mode = resolve_execution_mode(&args, &Some("10:0.95".to_string()), false);
        assert_eq!(mode, ExecutionMode::Simple);
    }

    #[test]
    fn test_resolve_execution_mode_amplicon_without_region() {
        let args = parse_args_for_test(["vardict", "-G", "reference.fa", "-b", "input.bam"]);

        let mode = resolve_execution_mode(&args, &Some("10:0.95".to_string()), false);
        assert_eq!(mode, ExecutionMode::Amplicon);
    }

    #[test]
    fn test_resolve_execution_mode_amplicon_without_region_overrides_paired_bam() {
        let args = parse_args_for_test([
            "vardict",
            "-G",
            "reference.fa",
            "-b",
            "tumor.bam|normal.bam",
        ]);

        let mode = resolve_execution_mode(&args, &Some("10:0.95".to_string()), true);
        assert_eq!(mode, ExecutionMode::Amplicon);
    }

    #[test]
    fn test_resolve_execution_mode_paired_bam_forces_somatic() {
        let args = parse_args_for_test([
            "vardict",
            "-G",
            "reference.fa",
            "-b",
            "tumor.bam|normal.bam",
            "-R",
            "chr1:1-10",
            "-a",
            "10:0.95",
        ]);

        let mode = resolve_execution_mode(&args, &Some("10:0.95".to_string()), true);
        assert_eq!(mode, ExecutionMode::Somatic);
    }

    #[test]
    fn test_resolve_execution_mode_splicing_has_highest_priority() {
        let args = parse_args_for_test([
            "vardict",
            "-G",
            "reference.fa",
            "-b",
            "tumor.bam|normal.bam",
            "-R",
            "chr1:1-10",
            "-a",
            "10:0.95",
            "-i",
        ]);

        let mode = resolve_execution_mode(&args, &Some("10:0.95".to_string()), true);
        assert_eq!(mode, ExecutionMode::Splicing);
    }

    #[test]
    fn test_parse_args_unique_mode_alignment_flag() {
        let args = parse_args_for_test([
            "vardict",
            "-G",
            "reference.fa",
            "-b",
            "input.bam",
            "-R",
            "chr1:1-10",
            "-u",
        ]);

        assert!(args.unique_mode_alignment);
        assert!(!args.unique_mode_second_in_pair);
    }

    #[test]
    fn test_parse_args_unique_mode_second_in_pair_legacy_short() {
        let args = parse_args_for_test([
            "vardict",
            "-G",
            "reference.fa",
            "-b",
            "input.bam",
            "-R",
            "chr1:1-10",
            "-UN",
        ]);

        assert!(!args.no_sv);
        assert!(args.unique_mode_second_in_pair);
    }

    #[test]
    fn test_parse_args_fisher_long_option() {
        let args = parse_args_for_test([
            "vardict",
            "-G",
            "reference.fa",
            "-b",
            "input.bam",
            "-R",
            "chr1:1-10",
            "--fisher",
        ]);

        assert!(args.fisher);
    }

    #[test]
    fn test_parse_args_fisher_legacy_single_dash_option() {
        let args = parse_args_for_test([
            "vardict",
            "-G",
            "reference.fa",
            "-b",
            "input.bam",
            "-R",
            "chr1:1-10",
            "-fisher",
        ]);

        assert!(args.fisher);
    }

    #[test]
    fn test_parse_args_crispr_options() {
        let args = parse_args_for_test([
            "vardict",
            "-G",
            "reference.fa",
            "-b",
            "input.bam",
            "-R",
            "chr1:1-10",
            "-J",
            "50454941",
            "-j",
            "25",
        ]);

        assert_eq!(args.crispr_cutting_site, 50_454_941);
        assert_eq!(args.crispr_filtering_bp, 25);
    }

    #[test]
    fn test_parse_bam_inputs_single_bam() {
        let inputs = parse_bam_inputs("sample.bam").unwrap();
        assert_eq!(inputs.primary_bam, PathBuf::from("sample.bam"));
        assert!(inputs.secondary_bam.is_none());
    }

    #[test]
    fn test_parse_bam_inputs_somatic_pair() {
        let inputs = parse_bam_inputs("tumor.bam|normal.bam").unwrap();
        assert_eq!(inputs.primary_bam, PathBuf::from("tumor.bam"));
        assert_eq!(inputs.secondary_bam, Some(PathBuf::from("normal.bam")));
    }

    #[test]
    fn test_parse_bam_inputs_rejects_more_than_two_entries() {
        let parsed = parse_bam_inputs("a.bam|b.bam|c.bam");
        assert!(parsed.is_err());
    }

    #[test]
    fn test_select_region_batches_for_execution_simple_mode_uses_flat_regions() {
        let regions = vec![
            Region::new("chr1".to_string(), 10, 20, "G1".to_string()),
            Region::new("chr1".to_string(), 30, 40, "G2".to_string()),
        ];

        let batches = select_region_batches_for_execution(ExecutionMode::Simple, regions, None, 1);
        assert_eq!(batches.len(), 1);
        assert_eq!(batches[0].len(), 2);
    }

    #[test]
    fn test_select_region_batches_for_execution_simple_mode_chunks_large_inputs() {
        let regions = (0..65)
            .map(|index| {
                Region::new(
                    "chr1".to_string(),
                    1 + index * 10,
                    10 + index * 10,
                    "G".to_string(),
                )
            })
            .collect::<Vec<_>>();

        let batches = select_region_batches_for_execution(ExecutionMode::Simple, regions, None, 1);

        assert_eq!(batches.len(), 9);
        assert_eq!(batches[0].len(), 8);
        assert_eq!(batches[7].len(), 8);
        assert_eq!(batches[8].len(), 1);
        assert_eq!(batches[0][0].start(), 1);
        assert_eq!(batches[8][0].start(), 641);
    }

    #[test]
    fn test_select_region_batches_for_execution_simple_mode_scales_with_threads() {
        let regions = (0..10)
            .map(|index| {
                Region::new(
                    "chr1".to_string(),
                    1 + index * 10,
                    10 + index * 10,
                    "G".to_string(),
                )
            })
            .collect::<Vec<_>>();

        let batches = select_region_batches_for_execution(ExecutionMode::Simple, regions, None, 4);

        assert_eq!(batches.len(), 1);
        assert_eq!(batches[0].len(), 10);
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

        let batches =
            select_region_batches_for_execution(ExecutionMode::Amplicon, regions, Some(groups), 1);
        assert_eq!(batches.len(), 2);
        assert_eq!(batches[0].len(), 1);
        assert_eq!(batches[1].len(), 2);
    }
}
