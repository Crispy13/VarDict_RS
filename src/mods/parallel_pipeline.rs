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

use std::sync::Arc;
use std::thread;
use std::path::Path;

use anyhow::{Context, Result};

use crate::data::shared_reference::{SharedReference, SharedReferenceHandle};
use crate::data::region::Region;
use crate::data::bam_reader::BamReader;
use crate::mods::pipeline::{Pipeline, PipelineConfig};
use crate::mods::simple_variant_caller::SimpleVariantCaller;
use crate::mods::vardict_pipeline::VarDictPipeline;
use crate::mods::output_variant::Region as OutputRegion;
use crate::scopedata::global_read_only_scope::{GlobalReadOnlyScope, instance};

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
            reference_path, chromosomes
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
        use std::sync::mpsc;

        if regions.is_empty() {
            return Vec::new();
        }

        // Calculate regions per thread
        let regions_per_thread = (regions.len() + self.num_threads - 1) / self.num_threads;
        
        // Create channel for results
        let (tx, rx) = mpsc::channel();

        // Partition regions into chunks for each thread
        let region_chunks: Vec<Vec<Region>> = regions
            .chunks(regions_per_thread.max(1))
            .map(|c| c.to_vec())
            .collect();

        // Spawn worker threads
        let mut handles = Vec::new();
        
        for (thread_id, chunk) in region_chunks.into_iter().enumerate() {
            let tx = tx.clone();
            let reference = Arc::clone(&self.reference);
            let config = self.config.clone();
            let bam_path = bam_path.clone();

            let handle = thread::spawn(move || {
                process_region_chunk(
                    thread_id,
                    chunk,
                    reference,
                    config,
                    bam_path,
                    tx,
                )
            });
            
            handles.push(handle);
        }

        // Drop the original sender so rx.iter() will end
        drop(tx);

        // Collect results
        let mut results: Vec<RegionResult> = rx.iter().collect();

        // Wait for all threads to complete
        for handle in handles {
            if let Err(e) = handle.join() {
                eprintln!("Worker thread panicked: {:?}", e);
            }
        }

        // Sort results by region for consistent output order
        results.sort_by(|a, b| {
            (&a.region.chr(), a.region.start(), a.region.end())
                .cmp(&(&b.region.chr(), b.region.start(), b.region.end()))
        });

        results
    }

    /// Process regions and return all output lines
    pub fn process_regions_to_output<P: AsRef<Path> + Send + Sync + Clone + 'static>(
        &self,
        bam_path: P,
        regions: Vec<Region>,
    ) -> Vec<String> {
        let results = self.process_regions(bam_path, regions);
        
        results.into_iter()
            .flat_map(|r| r.output_lines)
            .collect()
    }

    /// Process regions using the real VarDict pipeline (CigarParser → VariantRealigner → ToVarsBuilder)
    /// 
    /// This is the preferred method that follows the Java VarDict Simple Mode flow.
    pub fn process_regions_vardict<P: AsRef<Path> + Send + Sync + Clone + 'static>(
        &self,
        bam_path: P,
        regions: Vec<Region>,
    ) -> Vec<RegionResult> {
        use std::sync::mpsc;

        if regions.is_empty() {
            return Vec::new();
        }

        // Calculate regions per thread
        let regions_per_thread = (regions.len() + self.num_threads - 1) / self.num_threads;
        
        // Create channel for results
        let (tx, rx) = mpsc::channel();

        // Partition regions into chunks for each thread
        let region_chunks: Vec<Vec<Region>> = regions
            .chunks(regions_per_thread.max(1))
            .map(|c| c.to_vec())
            .collect();

        // Spawn worker threads
        let mut handles = Vec::new();
        
        for (thread_id, chunk) in region_chunks.into_iter().enumerate() {
            let tx = tx.clone();
            let reference = Arc::clone(&self.reference);
            let config = self.config.clone();
            let bam_path = bam_path.clone();

            let handle = thread::spawn(move || {
                process_region_chunk_vardict(
                    thread_id,
                    chunk,
                    reference,
                    config,
                    bam_path,
                    tx,
                )
            });
            
            handles.push(handle);
        }

        // Drop the original sender so rx.iter() will end
        drop(tx);

        // Collect results
        let mut results: Vec<RegionResult> = rx.iter().collect();

        // Wait for all threads to complete
        for handle in handles {
            if let Err(e) = handle.join() {
                eprintln!("Worker thread panicked: {:?}", e);
            }
        }

        // Sort results by region for consistent output order
        results.sort_by(|a, b| {
            (&a.region.chr(), a.region.start(), a.region.end())
                .cmp(&(&b.region.chr(), b.region.start(), b.region.end()))
        });

        results
    }

    /// Process regions using VarDict pipeline and return all output lines
    pub fn process_regions_vardict_to_output<P: AsRef<Path> + Send + Sync + Clone + 'static>(
        &self,
        bam_path: P,
        regions: Vec<Region>,
    ) -> Vec<String> {
        let results = self.process_regions_vardict(bam_path, regions);
        
        results.into_iter()
            .flat_map(|r| r.output_lines)
            .collect()
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

/// Worker function to process a chunk of regions using VarDict pipeline
fn process_region_chunk_vardict<P: AsRef<Path>>(
    thread_id: usize,
    regions: Vec<Region>,
    reference: SharedReferenceHandle,
    config: PipelineConfig,
    bam_path: P,
    tx: std::sync::mpsc::Sender<RegionResult>,
) {
    use crackle_kit::tracing::{Level, event};
    
    let start_thread = std::time::Instant::now();
    
    // Each thread opens its own BAM reader
    let mut bam_reader = match BamReader::open(bam_path.as_ref().to_str().unwrap()) {
        Ok(reader) => reader,
        Err(e) => {
            // Send error for all regions
            for region in regions {
                let _ = tx.send(RegionResult {
                    region,
                    output_lines: Vec::new(),
                    error: Some(format!("Failed to open BAM: {}", e)),
                });
            }
            return;
        }
    };

    // Create VarDict pipeline for this thread
    let vardict_pipeline = VarDictPipeline::new(&config.sample_name)
        .with_min_frequency(config.min_frequency)
        .with_min_base_quality(config.quality_threshold)
        .with_min_mapping_quality(config.mapq_threshold);

    // Get or create GlobalReadOnlyScope
    // In production, this would be set up once at startup
    let global_scope = Arc::new(instance().clone());

    // Process each region
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
        
        let _ = tx.send(result);
    }
    
    let elapsed_thread = start_thread.elapsed();
    event!(Level::INFO, "[TIMING] Thread {} completed in {:.3}s", thread_id, elapsed_thread.as_secs_f64());
}


/// Worker function to process a chunk of regions
fn process_region_chunk<P: AsRef<Path>>(
    _thread_id: usize,
    regions: Vec<Region>,
    reference: SharedReferenceHandle,
    config: PipelineConfig,
    bam_path: P,
    tx: std::sync::mpsc::Sender<RegionResult>,
) {
    // Each thread opens its own BAM reader
    let mut bam_reader = match BamReader::open(bam_path.as_ref().to_str().unwrap()) {
        Ok(reader) => reader,
        Err(e) => {
            // Send error for all regions
            for region in regions {
                let _ = tx.send(RegionResult {
                    region,
                    output_lines: Vec::new(),
                    error: Some(format!("Failed to open BAM: {}", e)),
                });
            }
            return;
        }
    };

    // Create variant caller for this thread
    let caller = SimpleVariantCaller::new(
        config.quality_threshold as u8,
        config.mapq_threshold,
    );

    // Create pipeline for this thread
    let pipeline = Pipeline::new(config.clone());

    // Process each region
    for region in regions {
        let result = process_single_region(
            &region,
            &reference,
            &caller,
            &pipeline,
            &mut bam_reader,
            &config,
        );
        
        let _ = tx.send(result);
    }
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
    use std::collections::HashMap;

    #[test]
    fn test_parallel_pipeline_creation() {
        // Create mock shared reference
        let mut chromosomes = HashMap::new();
        chromosomes.insert("chr1".to_string(), crate::data::shared_reference::ChromosomeData {
            sequence: b"ACGTACGT".repeat(100),
            length: 800,
        });

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
}
