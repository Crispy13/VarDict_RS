//! Multi-threaded Pipeline for Parallel Variant Calling
//!
//! This module provides parallel processing of regions using a thread pool
//! with shared reference genome data.
//!
//! Architecture:
//! ```text
//! Main Thread                   Worker Threads
//! ────────────                  ──────────────
//! Load reference ──► Arc<SharedReference> ◄── read only access
//!       │
//!       ▼
//! Spawn thread pool
//!       │
//!       ├──► Worker 1: Process regions [0, 1, 2, ...]
//!       ├──► Worker 2: Process regions [n, n+1, n+2, ...]
//!       ├──► Worker 3: Process regions [m, m+1, m+2, ...]
//!       └──► ...
//!             │
//!             ▼
//!       Collect results ──► Output
//! ```

use std::path::Path;
use std::sync::Arc;

use anyhow::Result;
use rayon::prelude::*;
use rust_htslib::bam::{HeaderView, Record, ext::BamRecordExtensions};

use crate::data::bam_reader::BamReader;
use crate::data::region::Region;
use crate::data::shared_reference::SharedReferenceHandle;
use crate::mods::pipeline::{Pipeline, PipelineConfig};
use crate::mods::simple_variant_caller::SimpleVariantCaller;
use crate::mods::vardict_pipeline::VarDictPipeline;
use crate::scopedata::global_read_only_scope::instance;

const SMALL_BATCH_PREFETCH_MULTIPLIER: usize = 4;
const MAX_PREFETCH_GROUP_SPAN_BP: usize = 2_500_000;

/// Result from processing a single region
#[derive(Debug)]
pub struct RegionResult {
    /// Region that was processed
    pub region: Region,
    /// Output lines (variant calls)
    pub output_lines: Vec<String>,
    /// Any error message
    pub error: Option<String>,
}

/// Multi-threaded pipeline runner
pub struct ParallelPipeline {
    /// Shared reference genome (thread-safe)
    reference: SharedReferenceHandle,
    /// Pipeline configuration
    config: PipelineConfig,
    /// Number of worker threads
    num_threads: usize,
}

impl ParallelPipeline {
    /// Create a new parallel pipeline
    pub fn new(
        reference: SharedReferenceHandle,
        config: PipelineConfig,
        num_threads: usize,
    ) -> Self {
        ParallelPipeline {
            reference,
            config,
            num_threads: num_threads.max(1),
        }
    }

    /// Load reference and create pipeline
    pub fn with_reference<P: AsRef<Path>>(
        reference_path: P,
        config: PipelineConfig,
        num_threads: usize,
    ) -> Result<Self> {
        let reference = crate::data::shared_reference::load_shared_reference(reference_path)?;
        Ok(Self::new(reference, config, num_threads))
    }

    /// Load specific chromosomes and create pipeline
    pub fn with_chromosomes<P: AsRef<Path>>(
        reference_path: P,
        chromosomes: &[&str],
        config: PipelineConfig,
        num_threads: usize,
    ) -> Result<Self> {
        let reference = crate::data::shared_reference::load_shared_reference_chroms(
            reference_path,
            chromosomes,
        )?;
        Ok(Self::new(reference, config, num_threads))
    }

    /// Process regions in parallel
    ///
    /// Arguments:
    /// * `bam_path` - Path to the BAM file (each thread opens its own handle)
    /// * `regions` - List of regions to process
    ///
    /// Returns results for each region
    pub fn process_regions<P: AsRef<Path> + Send + Sync + Clone + 'static>(
        &self,
        bam_path: P,
        regions: Vec<Region>,
    ) -> Vec<RegionResult> {
        if regions.is_empty() {
            return Vec::new();
        }

        // Calculate regions per thread
        let regions_per_thread = (regions.len() + self.num_threads - 1) / self.num_threads;

        let chunk_size = regions_per_thread.max(1);

        // Partition regions into chunks for each thread and keep original chunk start index
        let region_chunks: Vec<(usize, Vec<Region>)> = regions
            .chunks(chunk_size)
            .enumerate()
            .map(|(chunk_index, chunk)| (chunk_index * chunk_size, chunk.to_vec()))
            .collect();

        let mut chunk_results: Vec<(usize, Vec<RegionResult>)> = region_chunks
            .into_par_iter()
            .map(|(chunk_start_index, chunk)| {
                let reference = Arc::clone(&self.reference);
                let config = self.config.clone();
                let bam_path = bam_path.clone();
                (
                    chunk_start_index,
                    process_region_chunk(chunk, reference, config, bam_path),
                )
            })
            .collect();

        // Preserve original BED input order when flattening chunk results
        chunk_results.sort_by_key(|(chunk_start_index, _)| *chunk_start_index);

        chunk_results
            .into_iter()
            .flat_map(|(_, results)| results)
            .collect()
    }

    /// Process regions and return all output lines
    pub fn process_regions_to_output<P: AsRef<Path> + Send + Sync + Clone + 'static>(
        &self,
        bam_path: P,
        regions: Vec<Region>,
    ) -> Vec<String> {
        let results = self.process_regions(bam_path, regions);

        results.into_iter().flat_map(|r| r.output_lines).collect()
    }

    /// Process regions using the real VarDict pipeline (CigarParser → VariantRealigner → ToVarsBuilder)
    ///
    /// This is the preferred method that follows the Java VarDict Simple Mode flow.
    pub fn process_regions_vardict<P: AsRef<Path> + Send + Sync + Clone + 'static>(
        &self,
        bam_path: P,
        regions: Vec<Region>,
    ) -> Vec<RegionResult> {
        if regions.is_empty() {
            return Vec::new();
        }

        if should_use_prefetched_region_groups(regions.len(), self.num_threads) {
            return process_prefetched_region_groups_vardict(
                regions,
                Arc::clone(&self.reference),
                self.config.clone(),
                bam_path,
            );
        }

        // Calculate regions per thread
        let regions_per_thread = (regions.len() + self.num_threads - 1) / self.num_threads;

        let chunk_size = regions_per_thread.max(1);

        // Partition regions into chunks for each thread and keep original chunk start index
        let region_chunks: Vec<(usize, Vec<Region>)> = regions
            .chunks(chunk_size)
            .enumerate()
            .map(|(chunk_index, chunk)| (chunk_index * chunk_size, chunk.to_vec()))
            .collect();

        let mut chunk_results: Vec<(usize, Vec<RegionResult>)> = region_chunks
            .into_par_iter()
            .map(|(chunk_start_index, chunk)| {
                let reference = Arc::clone(&self.reference);
                let config = self.config.clone();
                let bam_path = bam_path.clone();
                (
                    chunk_start_index,
                    process_region_chunk_vardict(chunk, reference, config, bam_path),
                )
            })
            .collect();

        // Preserve original BED input order when flattening chunk results
        chunk_results.sort_by_key(|(chunk_start_index, _)| *chunk_start_index);

        chunk_results
            .into_iter()
            .flat_map(|(_, results)| results)
            .collect()
    }

    /// Process regions using VarDict pipeline and return all output lines
    pub fn process_regions_vardict_to_output<P: AsRef<Path> + Send + Sync + Clone + 'static>(
        &self,
        bam_path: P,
        regions: Vec<Region>,
    ) -> Vec<String> {
        let results = self.process_regions_vardict(bam_path, regions);

        results.into_iter().flat_map(|r| r.output_lines).collect()
    }

    /// Get header line
    pub fn get_header(&self) -> String {
        crate::mods::output_variant::get_header_line()
    }

    /// Get the shared reference
    pub fn reference(&self) -> &SharedReferenceHandle {
        &self.reference
    }
}

fn should_use_prefetched_region_groups(region_count: usize, num_threads: usize) -> bool {
    region_count > 1 && region_count <= num_threads.saturating_mul(SMALL_BATCH_PREFETCH_MULTIPLIER)
}

fn build_prefetch_region_groups(regions: Vec<Region>) -> Vec<Vec<(usize, Region)>> {
    let mut groups = Vec::new();
    let mut current_group: Vec<(usize, Region)> = Vec::new();
    let mut current_group_start = 0usize;
    let mut current_group_end = 0usize;

    for (index, region) in regions.into_iter().enumerate() {
        let starts_new_group = if let Some((_, first_region)) = current_group.first() {
            let next_group_end = current_group_end.max(region.end());
            let next_span = next_group_end.saturating_sub(current_group_start);
            first_region.chr() != region.chr() || next_span > MAX_PREFETCH_GROUP_SPAN_BP
        } else {
            false
        };

        if starts_new_group {
            groups.push(current_group);
            current_group = Vec::new();
        }

        if current_group.is_empty() {
            current_group_start = region.start();
            current_group_end = region.end();
        } else {
            current_group_end = current_group_end.max(region.end());
        }

        current_group.push((index, region));
    }

    if !current_group.is_empty() {
        groups.push(current_group);
    }

    groups
}

fn collect_prefetched_records(
    bam_reader: &mut BamReader,
    chrom: &str,
    start: usize,
    end: usize,
) -> Result<Vec<Record>> {
    let mut records = Vec::new();
    let header_view = Arc::new(HeaderView::from_header(bam_reader.header()));

    bam_reader.fetch(chrom, start, end)?;

    let mut record = Record::new();
    while bam_reader.read(&mut record)? {
        let mut cloned = record.clone();
        cloned.set_header(Arc::clone(&header_view));
        records.push(cloned);
    }

    Ok(records)
}

fn record_overlaps_region(record: &Record, region: &Region) -> bool {
    if record.is_unmapped() || record.pos() < 0 {
        return false;
    }

    let alignment_start = record.pos() + 1;
    let alignment_end = record.reference_end();
    alignment_start <= region.end() as i64 && alignment_end >= region.start() as i64
}

fn partition_prefetched_records_by_region(
    records: Vec<Record>,
    regions: &[(usize, Region)],
) -> Vec<Vec<Record>> {
    let mut buckets = vec![Vec::new(); regions.len()];

    for record in &records {
        for (bucket_index, (_, region)) in regions.iter().enumerate() {
            if record_overlaps_region(record, region) {
                buckets[bucket_index].push(record.clone());
            }
        }
    }

    buckets
}

fn process_prefetched_region_groups_vardict<P: AsRef<Path>>(
    regions: Vec<Region>,
    reference: SharedReferenceHandle,
    config: PipelineConfig,
    bam_path: P,
) -> Vec<RegionResult> {
    use crackle_kit::tracing::{Level, event};

    let start_thread = std::time::Instant::now();
    let mut indexed_results = Vec::new();
    let bam_path_str = match bam_path.as_ref().to_str() {
        Some(path) => path,
        None => {
            return regions
                .into_iter()
                .map(|region| RegionResult {
                    region,
                    output_lines: Vec::new(),
                    error: Some("BAM path is not valid UTF-8".to_string()),
                })
                .collect();
        }
    };

    let global_scope = Arc::new(instance().clone());
    let mut bam_reader = match BamReader::open(bam_path_str) {
        Ok(reader) => reader,
        Err(error) => {
            return regions
                .into_iter()
                .map(|region| RegionResult {
                    region,
                    output_lines: Vec::new(),
                    error: Some(format!("Failed to open BAM: {}", error)),
                })
                .collect();
        }
    };
    let target_names = bam_reader.target_names();

    let region_groups = build_prefetch_region_groups(regions);
    event!(
        Level::INFO,
        "[TIMING] Prefetch path engaged for {} region group(s)",
        region_groups.len()
    );

    for group in region_groups {
        let group_start = group
            .iter()
            .map(|(_, region)| region.start())
            .min()
            .unwrap_or(1);
        let group_end = group
            .iter()
            .map(|(_, region)| region.end())
            .max()
            .unwrap_or(group_start);
        let group_chr = match group.first() {
            Some((_, region)) => region.chr().to_string(),
            None => continue,
        };

        let prefetched_records =
            match collect_prefetched_records(&mut bam_reader, &group_chr, group_start, group_end) {
                Ok(records) => records,
                Err(error) => {
                    indexed_results.extend(group.into_iter().map(|(index, region)| {
                        (
                            index,
                            RegionResult {
                                region,
                                output_lines: Vec::new(),
                                error: Some(format!("Processing error: {}", error)),
                            },
                        )
                    }));
                    continue;
                }
            };

        let record_buckets = partition_prefetched_records_by_region(prefetched_records, &group);
        let group_config = config.clone();
        let group_reference = Arc::clone(&reference);
        let group_scope = Arc::clone(&global_scope);

        let mut group_results: Vec<(usize, RegionResult)> = group
            .into_par_iter()
            .zip(record_buckets.into_par_iter())
            .map(|((index, region), records)| {
                let pipeline = VarDictPipeline::new(&group_config.sample_name)
                    .with_min_frequency(group_config.min_frequency)
                    .with_min_base_quality(group_config.quality_threshold)
                    .with_min_mapping_quality(group_config.mapq_threshold);

                let result = match pipeline.process_region_from_cached_records(
                    &region,
                    &group_reference,
                    records,
                    &target_names,
                    Arc::clone(&group_scope),
                ) {
                    Ok(output_lines) => RegionResult {
                        region,
                        output_lines,
                        error: None,
                    },
                    Err(error) => RegionResult {
                        region,
                        output_lines: Vec::new(),
                        error: Some(format!("Processing error: {}", error)),
                    },
                };

                (index, result)
            })
            .collect();

        indexed_results.append(&mut group_results);
    }

    indexed_results.sort_by_key(|(index, _)| *index);

    let elapsed_thread = start_thread.elapsed();
    event!(
        Level::INFO,
        "[TIMING] Prefetched worker completed in {:.3}s",
        elapsed_thread.as_secs_f64()
    );

    indexed_results
        .into_iter()
        .map(|(_, result)| result)
        .collect()
}

/// Worker function to process a chunk of regions using VarDict pipeline
fn process_region_chunk_vardict<P: AsRef<Path>>(
    regions: Vec<Region>,
    reference: SharedReferenceHandle,
    config: PipelineConfig,
    bam_path: P,
) -> Vec<RegionResult> {
    use crackle_kit::tracing::{Level, event};

    let start_thread = std::time::Instant::now();
    let mut results = Vec::new();
    let bam_path_str = match bam_path.as_ref().to_str() {
        Some(path) => path,
        None => {
            for region in regions {
                results.push(RegionResult {
                    region,
                    output_lines: Vec::new(),
                    error: Some("BAM path is not valid UTF-8".to_string()),
                });
            }
            return results;
        }
    };

    let vardict_pipeline = VarDictPipeline::new(&config.sample_name)
        .with_min_frequency(config.min_frequency)
        .with_min_base_quality(config.quality_threshold)
        .with_min_mapping_quality(config.mapq_threshold);

    let global_scope = Arc::new(instance().clone());

    // Reuse one indexed BAM handle across the worker's region chunk to avoid
    // reopening the BAM and index for every single-region shard.
    let mut bam_reader = match BamReader::open(bam_path_str) {
        Ok(reader) => reader,
        Err(e) => {
            for region in regions {
                results.push(RegionResult {
                    region,
                    output_lines: Vec::new(),
                    error: Some(format!("Failed to open BAM: {}", e)),
                });
            }
            return results;
        }
    };

    for region in regions {
        let result = match vardict_pipeline.process_region_from_bam(
            &region,
            &reference,
            &mut bam_reader,
            Arc::clone(&global_scope),
        ) {
            Ok(output_lines) => RegionResult {
                region: region.clone(),
                output_lines,
                error: None,
            },
            Err(e) => RegionResult {
                region: region.clone(),
                output_lines: Vec::new(),
                error: Some(format!("Processing error: {}", e)),
            },
        };

        results.push(result);
    }

    let elapsed_thread = start_thread.elapsed();
    event!(
        Level::INFO,
        "[TIMING] Worker completed in {:.3}s",
        elapsed_thread.as_secs_f64()
    );

    results
}

/// Worker function to process a chunk of regions
fn process_region_chunk<P: AsRef<Path>>(
    regions: Vec<Region>,
    reference: SharedReferenceHandle,
    config: PipelineConfig,
    bam_path: P,
) -> Vec<RegionResult> {
    let bam_path_str = match bam_path.as_ref().to_str() {
        Some(path) => path,
        None => {
            let mut results = Vec::new();
            for region in regions {
                results.push(RegionResult {
                    region,
                    output_lines: Vec::new(),
                    error: Some("BAM path is not valid UTF-8".to_string()),
                });
            }
            return results;
        }
    };

    // Each thread opens its own BAM reader
    let mut bam_reader = match BamReader::open(bam_path_str) {
        Ok(reader) => reader,
        Err(e) => {
            let mut results = Vec::new();
            for region in regions {
                results.push(RegionResult {
                    region,
                    output_lines: Vec::new(),
                    error: Some(format!("Failed to open BAM: {}", e)),
                });
            }
            return results;
        }
    };

    // Create variant caller for this thread
    let caller = SimpleVariantCaller::new(config.quality_threshold as u8, config.mapq_threshold);

    // Create pipeline for this thread
    let pipeline = Pipeline::new(config.clone());

    let mut results = Vec::new();
    for region in regions {
        let result = process_single_region(
            &region,
            &reference,
            &caller,
            &pipeline,
            &mut bam_reader,
            &config,
        );

        results.push(result);
    }

    results
}

/// Process a single region using VarDictPipeline (Java-equivalent flow)
fn process_single_region(
    region: &Region,
    reference: &SharedReferenceHandle,
    _caller: &SimpleVariantCaller,
    _pipeline: &Pipeline,
    bam_reader: &mut BamReader,
    config: &PipelineConfig,
) -> RegionResult {
    // Create VarDictPipeline with configuration
    let pipeline = VarDictPipeline::new(&config.sample_name)
        .with_min_frequency(config.min_frequency)
        .with_min_base_quality(config.quality_threshold)
        .with_min_mapping_quality(config.mapq_threshold)
        .with_pileup(config.pileup); // Use config's pileup setting

    // Get GlobalReadOnlyScope instance
    let gros = Arc::new(instance().clone());

    // Process region through the real VarDict pipeline
    match pipeline.process_region_from_bam(region, reference, bam_reader, gros) {
        Ok(output_lines) => RegionResult {
            region: region.clone(),
            output_lines,
            error: None,
        },
        Err(e) => RegionResult {
            region: region.clone(),
            output_lines: Vec::new(),
            error: Some(format!("Pipeline error: {}", e)),
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::shared_reference::SharedReference;
    use std::collections::HashMap;

    #[test]
    fn test_parallel_pipeline_creation() {
        // Create mock shared reference
        let mut chromosomes: HashMap<
            String,
            crate::data::shared_reference::ChromosomeData,
            crate::prelude::LibDefaultHasher,
        > = Default::default();
        chromosomes.insert(
            "chr1".to_string(),
            crate::data::shared_reference::ChromosomeData {
                sequence: b"ACGTACGT".repeat(100),
                length: 800,
            },
        );

        let reference = Arc::new(SharedReference {
            chromosomes,
            chromosome_names: vec!["chr1".to_string()],
            total_size: 800,
        });

        let config = PipelineConfig::default();
        let pipeline = ParallelPipeline::new(reference, config, 4);

        assert_eq!(pipeline.num_threads, 4);
    }

    #[test]
    fn test_region_partitioning() {
        // Test that regions are properly partitioned across threads
        let regions: Vec<Region> = (0..10)
            .map(|i| Region::new("chr1".to_string(), i * 1000, (i + 1) * 1000, String::new()))
            .collect();

        let num_threads = 4;
        let regions_per_thread = (regions.len() + num_threads - 1) / num_threads;

        let chunks: Vec<Vec<Region>> = regions
            .chunks(regions_per_thread)
            .map(|c| c.to_vec())
            .collect();

        assert_eq!(chunks.len(), 4);
        assert!(chunks.iter().all(|c| c.len() <= 3));
    }

    #[test]
    fn test_prefetch_region_groups_split_by_span_and_contig() {
        let regions = vec![
            Region::new("chr1".to_string(), 1, 1_000_000, String::new()),
            Region::new("chr1".to_string(), 1_000_001, 2_000_000, String::new()),
            Region::new("chr1".to_string(), 2_000_001, 2_600_000, String::new()),
            Region::new("chr2".to_string(), 1, 10_000, String::new()),
        ];

        let groups = build_prefetch_region_groups(regions);

        assert_eq!(groups.len(), 3);
        assert_eq!(groups[0].len(), 2);
        assert_eq!(groups[1].len(), 1);
        assert_eq!(groups[2].len(), 1);
        assert_eq!(groups[0][0].1.chr(), "chr1");
        assert_eq!(groups[2][0].1.chr(), "chr2");
    }
}
