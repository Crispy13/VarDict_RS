//! Shared Reference Data for Multi-threaded Processing
//!
//! This module provides a thread-safe, read-only reference genome that can be
//! loaded once and shared across all worker threads.
//!
//! Design rationale:
//! - Load entire reference genome (≈3GB for human) once at startup
//! - Use Arc<SharedReference> to share read-only data across threads
//! - No locks needed since data is immutable after initialization

use std::collections::HashMap;
use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use anyhow::{Context, Result, anyhow};
use bio::io::fasta::{IndexedReader, Sequence};
use crackle_kit::tracing::{Level, event};

use crate::prelude::LibDefaultHasher;

/// Normalize chromosome name to match reference
///
/// Tries the original name first, then with/without "chr" prefix
fn normalize_chrom_name(sequences: &[Sequence], chrom: &str) -> Option<String> {
    // Try exact match first
    if chromosome_len(sequences, chrom).is_some() {
        return Some(chrom.to_string());
    }

    // Try with "chr" prefix stripped
    if let Some(stripped) = chrom.strip_prefix("chr") {
        if chromosome_len(sequences, stripped).is_some() {
            return Some(stripped.to_string());
        }
    }

    // Try with "chr" prefix added
    let with_chr = format!("chr{}", chrom);
    if chromosome_len(sequences, &with_chr).is_some() {
        return Some(with_chr);
    }

    None
}

fn chromosome_len(sequences: &[Sequence], chrom: &str) -> Option<usize> {
    sequences
        .iter()
        .find(|sequence| sequence.name == chrom && sequence.len > 0)
        .map(|sequence| sequence.len as usize)
}

fn read_fetched_sequence(reader: &mut IndexedReader<File>) -> Result<Vec<u8>> {
    let mut sequence = Vec::new();
    reader
        .read(&mut sequence)
        .context("Failed to read fetched FASTA sequence")?;

    for base in &mut sequence {
        *base = base.to_ascii_uppercase();
    }

    Ok(sequence)
}

/// Chromosome sequence data
#[derive(Debug)]
pub struct ChromosomeData {
    /// Sequence bytes (uppercase)
    pub sequence: Arc<Vec<u8>>,
    /// Length of the sequence
    pub length: usize,
}

impl ChromosomeData {
    /// Get a base at a 1-based position
    #[inline]
    pub fn get_base(&self, pos: usize) -> Option<u8> {
        if pos == 0 || pos > self.length {
            return None;
        }
        self.sequence.get(pos - 1).copied()
    }

    /// Get a subsequence (1-based, inclusive)
    pub fn get_subseq(&self, start: usize, end: usize) -> Option<&[u8]> {
        if start == 0 || end < start || end > self.length {
            return None;
        }
        Some(&self.sequence[start - 1..end])
    }
}

/// Shared reference genome for multi-threaded access
///
/// This structure is designed to be wrapped in Arc<> and shared across threads.
/// All data is immutable after construction, so no synchronization is needed.
#[derive(Debug)]
pub struct SharedReference {
    /// Chromosome name -> sequence data
    pub chromosomes: HashMap<String, ChromosomeData, LibDefaultHasher>,
    /// Index of chromosome names for iteration
    pub chromosome_names: Vec<String>,
    /// Total size of all sequences in bytes
    pub total_size: usize,
}

impl SharedReference {
    /// Load a specific chromosome from a FASTA file
    ///
    /// Handles chromosome name normalization (e.g., "chr20" -> "20" or vice versa)
    pub fn load_chromosome<P: AsRef<Path>>(path: P, chrom: &str) -> Result<Self> {
        let fasta_path = path.as_ref();
        let mut reader = IndexedReader::from_file(&fasta_path)
            .with_context(|| format!("Failed to open FASTA: {:?}", path.as_ref()))?;
        let sequences = reader.index.sequences();

        // Normalize chromosome name to match reference
        let ref_chrom = normalize_chrom_name(&sequences, chrom).ok_or_else(|| {
            anyhow!(
                "Chromosome '{}' not found in reference (tried with/without 'chr' prefix)",
                chrom
            )
        })?;

        reader
            .fetch_all(&ref_chrom)
            .with_context(|| format!("Failed to fetch chromosome {}", ref_chrom))?;
        let sequence = read_fetched_sequence(&mut reader)
            .with_context(|| format!("Failed to load chromosome {}", ref_chrom))?;
        let seq_len = sequence.len();
        let sequence = Arc::new(sequence);

        let mut chromosomes: HashMap<String, ChromosomeData, LibDefaultHasher> = Default::default();
        // Store under BOTH the original name and the reference name for lookup flexibility
        chromosomes.insert(
            chrom.to_string(),
            ChromosomeData {
                sequence: Arc::clone(&sequence),
                length: seq_len,
            },
        );
        if ref_chrom != chrom {
            chromosomes.insert(
                ref_chrom.clone(),
                ChromosomeData {
                    sequence,
                    length: seq_len,
                },
            );
        }

        Ok(SharedReference {
            chromosomes,
            chromosome_names: vec![chrom.to_string()],
            total_size: seq_len,
        })
    }

    /// Load multiple chromosomes from a FASTA file
    ///
    /// Handles chromosome name normalization (e.g., "chr20" -> "20" or vice versa)
    pub fn load_chromosomes<P: AsRef<Path>>(path: P, chroms: &[&str]) -> Result<Self> {
        let fasta_path = path.as_ref();
        let mut reader = IndexedReader::from_file(&fasta_path)
            .with_context(|| format!("Failed to open FASTA: {:?}", path.as_ref()))?;
        let sequences = reader.index.sequences();

        let mut chromosomes: HashMap<String, ChromosomeData, LibDefaultHasher> = Default::default();
        let mut chromosome_names = Vec::new();
        let mut total_size = 0usize;

        for chrom in chroms {
            // Normalize chromosome name to match reference
            let ref_chrom = match normalize_chrom_name(&sequences, chrom) {
                Some(name) => name,
                None => {
                    event!(
                        Level::WARN,
                        "Chromosome '{}' not found in reference (tried with/without 'chr' prefix)",
                        chrom
                    );
                    continue;
                }
            };

            reader
                .fetch_all(&ref_chrom)
                .with_context(|| format!("Failed to fetch chromosome {}", ref_chrom))?;
            let sequence = read_fetched_sequence(&mut reader)
                .with_context(|| format!("Failed to load chromosome {}", ref_chrom))?;
            let seq_len = sequence.len();
            total_size += seq_len;
            let sequence = Arc::new(sequence);

            // Store under the original name (from BED file) for lookup
            chromosomes.insert(
                chrom.to_string(),
                ChromosomeData {
                    sequence: Arc::clone(&sequence),
                    length: seq_len,
                },
            );
            // Also store under the reference name if different
            if ref_chrom != *chrom {
                chromosomes.insert(
                    ref_chrom.clone(),
                    ChromosomeData {
                        sequence,
                        length: seq_len,
                    },
                );
            }
            chromosome_names.push(chrom.to_string());
        }

        Ok(SharedReference {
            chromosomes,
            chromosome_names,
            total_size,
        })
    }

    /// Load all chromosomes from a FASTA file
    ///
    /// Note: For human genome (~3GB), this will allocate ~3GB of memory.
    /// This is suitable for modern hardware with 8GB+ RAM.
    pub fn load_all<P: AsRef<Path>>(path: P) -> Result<Self> {
        let fasta_path = path.as_ref();
        let mut reader = IndexedReader::from_file(&fasta_path)
            .with_context(|| format!("Failed to open FASTA: {:?}", path.as_ref()))?;
        let sequences = reader.index.sequences();

        let mut chromosomes: HashMap<String, ChromosomeData, LibDefaultHasher> = Default::default();
        let mut chromosome_names = Vec::with_capacity(sequences.len());
        let mut total_size = 0usize;

        for sequence_info in sequences {
            if sequence_info.len == 0 {
                continue;
            }

            reader
                .fetch_all(&sequence_info.name)
                .with_context(|| format!("Failed to fetch chromosome {}", sequence_info.name))?;
            let sequence = read_fetched_sequence(&mut reader)
                .with_context(|| format!("Failed to load chromosome {}", sequence_info.name))?;
            let seq_len = sequence.len();
            total_size += seq_len;
            let sequence = Arc::new(sequence);

            chromosomes.insert(
                sequence_info.name.clone(),
                ChromosomeData {
                    sequence,
                    length: seq_len,
                },
            );
            chromosome_names.push(sequence_info.name);
        }

        event!(
            Level::INFO,
            "Loaded {} chromosomes, total size: {:.2} GB",
            chromosome_names.len(),
            total_size as f64 / 1_073_741_824.0
        );

        Ok(SharedReference {
            chromosomes,
            chromosome_names,
            total_size,
        })
    }

    /// Get a chromosome's data
    pub fn get_chromosome(&self, name: &str) -> Option<&ChromosomeData> {
        self.chromosomes.get(name)
    }

    /// Get a base at a specific position (1-based)
    #[inline]
    pub fn get_base(&self, chrom: &str, pos: usize) -> Option<u8> {
        self.chromosomes.get(chrom)?.get_base(pos)
    }

    /// Get a subsequence (1-based, inclusive)
    pub fn get_subseq(&self, chrom: &str, start: usize, end: usize) -> Option<&[u8]> {
        self.chromosomes.get(chrom)?.get_subseq(start, end)
    }

    /// Get the list of chromosome names
    pub fn chromosome_names(&self) -> &[String] {
        &self.chromosome_names
    }

    /// Get total size of all sequences
    pub fn total_size(&self) -> usize {
        self.total_size
    }

    /// Get number of chromosomes
    pub fn num_chromosomes(&self) -> usize {
        self.chromosomes.len()
    }

    /// Get chromosome lengths as a HashMap
    pub fn get_chromosome_lengths(&self) -> HashMap<String, usize, LibDefaultHasher> {
        self.chromosomes
            .iter()
            .map(|(name, data)| (name.clone(), data.length))
            .collect()
    }
}

/// Thread-safe reference handle
///
/// This is a cheap clone that shares the underlying data.
pub type SharedReferenceHandle = Arc<SharedReference>;

/// Create a shared reference handle from a FASTA file
pub fn load_shared_reference<P: AsRef<Path>>(path: P) -> Result<SharedReferenceHandle> {
    let reference = SharedReference::load_all(path)?;
    Ok(Arc::new(reference))
}

/// Create a shared reference for specific chromosomes
pub fn load_shared_reference_chroms<P: AsRef<Path>>(
    path: P,
    chroms: &[&str],
) -> Result<SharedReferenceHandle> {
    let reference = SharedReference::load_chromosomes(path, chroms)?;
    Ok(Arc::new(reference))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chromosome_data_get_base() {
        let data = ChromosomeData {
            sequence: Arc::new(b"ACGTACGT".to_vec()),
            length: 8,
        };

        assert_eq!(data.get_base(1), Some(b'A'));
        assert_eq!(data.get_base(4), Some(b'T'));
        assert_eq!(data.get_base(8), Some(b'T'));
        assert_eq!(data.get_base(0), None); // 0 is invalid
        assert_eq!(data.get_base(9), None); // Out of bounds
    }

    #[test]
    fn test_chromosome_data_get_subseq() {
        let data = ChromosomeData {
            sequence: Arc::new(b"ACGTACGT".to_vec()),
            length: 8,
        };

        assert_eq!(data.get_subseq(1, 4), Some(b"ACGT".as_slice()));
        assert_eq!(data.get_subseq(5, 8), Some(b"ACGT".as_slice()));
        assert_eq!(data.get_subseq(0, 4), None); // 0 is invalid
        assert_eq!(data.get_subseq(1, 9), None); // Out of bounds
    }

    #[test]
    fn test_shared_reference_threadsafe() {
        use std::thread;

        // Create a mock shared reference
        let mut chromosomes: HashMap<String, ChromosomeData, LibDefaultHasher> = Default::default();
        chromosomes.insert(
            "chr1".to_string(),
            ChromosomeData {
                sequence: Arc::new(b"ACGTACGTACGT".to_vec()),
                length: 12,
            },
        );
        chromosomes.insert(
            "chr2".to_string(),
            ChromosomeData {
                sequence: Arc::new(b"GGGGCCCCAAAA".to_vec()),
                length: 12,
            },
        );

        let reference = Arc::new(SharedReference {
            chromosomes,
            chromosome_names: vec!["chr1".to_string(), "chr2".to_string()],
            total_size: 24,
        });

        // Spawn threads that read from the shared reference
        let mut handles = Vec::new();
        for i in 0..4 {
            let ref_clone = Arc::clone(&reference);
            let handle = thread::spawn(move || {
                let base1 = ref_clone.get_base("chr1", 1 + i);
                let base2 = ref_clone.get_base("chr2", 1 + i);
                (base1, base2)
            });
            handles.push(handle);
        }

        // Collect results
        let results: Vec<_> = handles.into_iter().map(|h| h.join().unwrap()).collect();

        // Verify results
        assert_eq!(results[0], (Some(b'A'), Some(b'G')));
        assert_eq!(results[1], (Some(b'C'), Some(b'G')));
        assert_eq!(results[2], (Some(b'G'), Some(b'G')));
        assert_eq!(results[3], (Some(b'T'), Some(b'G')));

        // Only one strong reference should remain
        assert_eq!(Arc::strong_count(&reference), 1);
    }
}
