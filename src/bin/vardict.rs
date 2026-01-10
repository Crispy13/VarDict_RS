//! VarDict-rs: A Rust implementation of VarDict variant caller (Simple Mode)
//!
//! Usage: vardict -G <reference.fa> -b <input.bam> [options] <region or BED file>

use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::PathBuf;

use anyhow::{Context, Result, anyhow};
use clap::Parser;

use vardict_rs::data::bam_reader::{BamReader, passes_filter};
use vardict_rs::data::reference::FastaReader;
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

    /// Minimum base quality for a base to be considered (default: 25)
    #[arg(short = 'q', long = "qual", default_value = "25")]
    min_base_quality: u8,

    /// Minimum mapping quality for a read to be considered (default: 0)
    #[arg(short = 'Q', long = "mapq", default_value = "0")]
    min_mapping_quality: u8,

    /// Print header line
    #[arg(short = 'H', long = "header")]
    print_header: bool,

    /// Extension of bp to look for mismatches after indel (default: 2)
    #[arg(short = 'e', long = "vext", default_value = "2")]
    vext: i32,

    /// The hexical to filter reads (samtools style). Default: 0x504
    /// Filters: 2nd alignments, unmapped, duplicates. Use 0 to disable.
    #[arg(short = 'F', long = "filter", default_value = "0x504")]
    sam_filter: String,

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

    /// Perform local realignment (default: true, set to false for Ion/PacBio)
    #[arg(short = 'k', long = "realign", default_value = "true")]
    local_realignment: bool,

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
}

fn main() -> Result<()> {
    let args = Args::parse();

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
        .build();

    // Print header if requested
    if args.print_header {
        let pipeline = Pipeline::new(config.clone());
        println!("{}", pipeline.get_header());
    }

    // Open FASTA reference reader
    let fasta_reader = FastaReader::open(args.reference.to_str().unwrap())
        .context("Failed to open reference FASTA file")?;

    // Open BAM reader
    let mut bam_reader = BamReader::open(args.bam.to_str().unwrap())
        .context("Failed to open BAM file")?;

    // Process each region
    let mut stdout = io::stdout().lock();
    
    for region in &regions {
        if args.debug {
            eprintln!("Processing region: {}:{}-{}", region.chr(), region.start(), region.end());
        }

        match process_region(&args, &config, region, &sample_name, &mut bam_reader, &fasta_reader) {
            Ok(output_lines) => {
                for line in output_lines {
                    writeln!(stdout, "{}", line)?;
                }
            }
            Err(e) => {
                if args.debug {
                    eprintln!("Error processing region {}:{}-{}: {}", 
                        region.chr(), region.start(), region.end(), e);
                }
            }
        }
    }

    Ok(())
}

/// Parse regions from command line arguments
fn get_regions(args: &Args) -> Result<Vec<Region>> {
    let mut regions = Vec::new();

    // Check for -R option first
    if let Some(ref region_str) = args.region {
        let region = parse_region_string(region_str, args.zero_based)?;
        regions.push(region);
        return Ok(regions);
    }

    // Check for BED file
    if let Some(ref bed_path) = args.bed_file {
        if !bed_path.exists() {
            return Err(anyhow!("BED file not found: {:?}", bed_path));
        }
        regions = parse_bed_file(bed_path, args)?;
        return Ok(regions);
    }

    Ok(regions)
}

/// Parse a region string like "chr1:1000-2000" or "chr1:1000"
fn parse_region_string(s: &str, zero_based: bool) -> Result<Region> {
    // Format: chr:start-end or chr:start
    let parts: Vec<&str> = s.split(':').collect();
    if parts.len() != 2 {
        return Err(anyhow!("Invalid region format: {}. Expected chr:start-end", s));
    }

    let chr = parts[0].to_string();
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

    Ok(Region::new(chr, start, end, String::new()))
}

/// Parse a BED file and return regions
fn parse_bed_file(path: &PathBuf, args: &Args) -> Result<Vec<Region>> {
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

        let chr = fields[chr_idx].to_string();
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

/// Process a single region and return output lines
fn process_region(
    args: &Args,
    config: &PipelineConfig,
    region: &Region,
    _sample_name: &str,
    bam_reader: &mut BamReader,
    fasta_reader: &FastaReader,
) -> Result<Vec<String>> {
    use vardict_rs::mods::simple_variant_caller::SimpleVariantCaller;
    use vardict_rs::mods::output_variant::Region as OutputRegion;
    
    if args.debug {
        eprintln!("  Region: {}:{}-{} (gene: {})", 
            region.chr(), region.start(), region.end(), region.gene());
    }

    // 1. Fetch reference sequence for this region
    let ref_seq = fasta_reader.fetch_seq(
        region.chr(),
        region.start() as usize,
        region.end() as usize,
    )?;

    if args.debug {
        eprintln!("  Reference length: {} bp", ref_seq.len());
    }

    // 2. Create variant caller and set reference
    let mut caller = SimpleVariantCaller::new(
        config.quality_threshold,
        config.mapq_threshold,
    );
    caller.set_reference(ref_seq, region.start() as i64);

    // 3. Fetch reads from BAM
    bam_reader.fetch(region.chr(), region.start(), region.end())?;

    // Parse SAM filter
    let sam_filter: u32 = args.sam_filter.parse().unwrap_or(0x504);

    // 4. Read all records into a vector (needed for processing)
    let mut records = Vec::new();
    let mut record = rust_htslib::bam::Record::new();
    while bam_reader.read(&mut record)? {
        // Apply SAM flag filter (skip filtered reads)
        if passes_filter(&record, sam_filter, config.mapq_threshold) {
            records.push(record.clone());
        }
    }

    if args.debug {
        eprintln!("  Reads in region: {}", records.len());
    }

    if records.is_empty() {
        return Ok(vec![]);
    }

    // 5. Process all records
    let (variations, coverage) = caller.process_records(records.iter());

    if args.debug {
        eprintln!("  Variations found: {}", variations.len());
        eprintln!("  Coverage positions: {}", coverage.len());
    }

    if variations.is_empty() {
        return Ok(vec![]);
    }

    // 6. Create pipeline and process
    let pipeline = Pipeline::new(config.clone());
    
    let output_region = OutputRegion::new(
        region.chr(),
        region.start() as i64,
        region.end() as i64,
        region.gene(),
    );

    let output_lines = pipeline.process_variations(variations, &coverage, &output_region);

    Ok(output_lines)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_region_string_full() {
        let region = parse_region_string("chr1:1000-2000", false).unwrap();
        assert_eq!(region.chr(), "chr1");
        assert_eq!(region.start(), 1000);
        assert_eq!(region.end(), 2000);
    }

    #[test]
    fn test_parse_region_string_single_position() {
        let region = parse_region_string("chr1:1000", false).unwrap();
        assert_eq!(region.chr(), "chr1");
        assert_eq!(region.start(), 1000);
        assert_eq!(region.end(), 1000);
    }

    #[test]
    fn test_parse_region_string_zero_based() {
        let region = parse_region_string("chr1:999-2000", true).unwrap();
        assert_eq!(region.chr(), "chr1");
        assert_eq!(region.start(), 1000); // 999 + 1
        assert_eq!(region.end(), 2000);
    }

    #[test]
    fn test_parse_region_string_invalid() {
        assert!(parse_region_string("invalid", false).is_err());
        assert!(parse_region_string("chr1", false).is_err());
    }
}
