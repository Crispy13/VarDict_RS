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

    /// Debug mode - print additional information
    #[arg(short = 'D', long = "debug")]
    debug: bool,

    /// Output all variants including reference calls (pileup mode)
    #[arg(short = 'p', long = "pileup")]
    pileup: bool,

    /// Indicate coordinates are zero-based (default: 1 for BED)
    #[arg(short = 'z', long = "zero")]
    zero_based: bool,

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
    let regions = get_regions(&args)?;
    if regions.is_empty() {
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

    // Print header if requested
    if args.print_header {
        let pipeline = Pipeline::new(config.clone());
        println!("{}", pipeline.get_header());
    }

    // Always use SharedReference (loaded into memory for fast access)
    run_variant_calling(&args, config, regions)?;

    Ok(())
}

/// Run variant calling using SharedReference (loaded into memory)
/// 
/// SharedReference is the default for both single and multi-threaded modes.
/// The reference is loaded once and shared across all threads for fast access.
fn run_variant_calling(args: &Args, config: PipelineConfig, regions: Vec<Region>) -> Result<()> {
    use vardict_rs::data::shared_reference::load_shared_reference_chroms;
    use vardict_rs::mods::parallel_pipeline::ParallelPipeline;
    use vardict_rs::scopedata::global_read_only_scope::{GlobalReadOnlyScope, INSTANCE};
    use vardict_rs::conf::Configuration;
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
    conf.perform_local_realignment = args.local_realignment == 1;
    conf.number_nucleotide_to_extend = args.number_nucleotide_to_extend;
    conf.reference_extension = args.reference_extension;
    let mut scope = GlobalReadOnlyScope::default();
    scope.conf = conf;
    scope.chr_lens = reference.get_chromosome_lengths();
    scope.bam_paths = vec![args.bam.to_string_lossy().to_string()];
    let _ = INSTANCE.set(scope);

    if args.debug {
        eprintln!("[TIMING] Reference loading: {:.3}s", elapsed_ref_load.as_secs_f64());
        eprintln!("Loaded {} chromosome(s), {:.2} MB total",
            reference.num_chromosomes(),
            reference.total_size() as f64 / 1_048_576.0);
        if num_threads > 1 {
            eprintln!("Processing {} regions with {} threads...", regions.len(), num_threads);
        } else {
            eprintln!("Processing {} regions...", regions.len());
        }
    }

    // Create parallel pipeline (works for single thread too)
    let pipeline = ParallelPipeline::new(reference, config, num_threads);

    // Process regions
    let start_processing = Instant::now();
    let bam_path = args.bam.to_str().unwrap().to_string();
    let results = pipeline.process_regions_vardict(bam_path, regions);
    let elapsed_processing = start_processing.elapsed();

    // Output results
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

    Ok(())
}

/// Parse regions from command line arguments
fn get_regions(args: &Args) -> Result<Vec<Region>> {
    let mut regions = Vec::new();
    let bam_targets = BamReader::open(&args.bam)
        .context("Failed to open BAM for region normalization")?
        .target_names();

    // Check for -R option first
    if let Some(ref region_str) = args.region {
        let region_path = PathBuf::from(region_str);
        if region_path.exists() {
            regions = parse_bed_file(&region_path, args, Some(&bam_targets))?;
            return Ok(regions);
        }
        let region = parse_region_string(region_str, args.zero_based, Some(&bam_targets))?;
        regions.push(region);
        return Ok(regions);
    }

    // Check for BED file
    if let Some(ref bed_path) = args.bed_file {
        if !bed_path.exists() {
            return Err(anyhow!("BED file not found: {:?}", bed_path));
        }
        regions = parse_bed_file(bed_path, args, Some(&bam_targets))?;
        return Ok(regions);
    }

    Ok(regions)
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
) -> Result<Vec<Region>> {
    let file = File::open(path).context("Failed to open BED file")?;
    let reader = BufReader::new(file);
    let mut regions = Vec::new();

    for (line_num, line) in reader.lines().enumerate() {
        let line = line.context("Failed to read BED line")?;
        let line = line.trim();

        // Skip empty lines and comments
        if line.is_empty() || line.starts_with('#') || line.starts_with("track") || line.starts_with("browser") {
            continue;
        }

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
        let (start, end) = if args.zero_based {
            (start + 1, end) // BED format: 0-based start, end is exclusive
        } else {
            (start, end)
        };

        regions.push(Region::new(chr, start, end, gene));
    }

    Ok(regions)
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
}
