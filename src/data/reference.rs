use std::collections::HashMap;
use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use anyhow::{Context, Result, anyhow};
use bio::io::fasta::IndexedReader;

use crate::conf::Configuration;
use crate::prelude::LibDefaultHasher;

#[derive(Default, Clone, Debug)]
pub struct ReferenceSeedMap {
    seed1: HashMap<[u8; 17], Vec<i64>, LibDefaultHasher>,
    seed2: HashMap<[u8; 12], Vec<i64>, LibDefaultHasher>,
}

impl ReferenceSeedMap {
    pub fn get(&self, key: &[u8]) -> Option<&Vec<i64>> {
        match key.len() {
            17 => self.seed1.get(<&[u8; 17]>::try_from(key).ok()?),
            12 => self.seed2.get(<&[u8; 12]>::try_from(key).ok()?),
            _ => None,
        }
    }

    pub(crate) fn insert_seed1(&mut self, key: [u8; 17], pos: i64) {
        self.seed1.entry(key).or_default().push(pos);
    }

    pub(crate) fn insert_seed2(&mut self, key: [u8; 12], pos: i64) {
        self.seed2.entry(key).or_default().push(pos);
    }

    pub(crate) fn reserve(&mut self, region_len: usize) {
        let seed_capacity = region_len.saturating_mul(60) / 100;
        self.seed1.reserve(seed_capacity);
        self.seed2.reserve(seed_capacity);
    }

    pub(crate) fn clear(&mut self) {
        self.seed1.clear();
        self.seed2.clear();
    }

    pub(crate) fn shrink_to_fit(&mut self) {
        self.seed1.shrink_to_fit();
        self.seed2.shrink_to_fit();
    }

    pub fn is_empty(&self) -> bool {
        self.seed1.is_empty() && self.seed2.is_empty()
    }
}

/// Reference sequence data
///
/// Fields use Arc internally so that `Reference::clone()` is O(1) (refcount bump)
/// rather than deep-copying the sequence and seed map.
#[derive(Default, Clone)]
pub struct Reference {
    pub ref_seq: Arc<Vec<u8>>,
    pub seed: Arc<ReferenceSeedMap>,
    /// Start position of this reference slice in genomic coordinates (1-based)
    pub region_start: i64,
}

impl Reference {
    /// Create a new Reference from a sequence slice
    pub fn from_seq(seq: &[u8]) -> Self {
        Reference {
            ref_seq: Arc::new(seq.to_vec()),
            seed: Default::default(),
            region_start: 0,
        }
    }

    /// Create a new Reference from a sequence slice with region start position
    pub fn from_seq_with_start(seq: &[u8], region_start: i64) -> Self {
        Reference {
            ref_seq: Arc::new(seq.to_vec()),
            seed: Default::default(),
            region_start,
        }
    }

    /// Create a new Reference with owned sequence
    pub fn new(ref_seq: Vec<u8>) -> Self {
        Reference {
            ref_seq: Arc::new(ref_seq),
            seed: Default::default(),
            region_start: 0,
        }
    }

    /// Create a new Reference with owned sequence and region start
    pub fn new_with_start(ref_seq: Vec<u8>, region_start: i64) -> Self {
        Reference {
            ref_seq: Arc::new(ref_seq),
            seed: Default::default(),
            region_start,
        }
    }

    /// Build the reference seed map using SEED_1 and SEED_2 lengths.
    pub fn build_seed_map(&mut self, region_end: i64, chr_len: Option<usize>) {
        let seed_map = Arc::make_mut(&mut self.seed);
        seed_map.clear();

        if self.ref_seq.is_empty() {
            return;
        }

        let seed1 = Configuration::SEED_1 as usize;
        let seed2 = Configuration::SEED_2 as usize;
        let seq_len = self.ref_seq.len();
        seed_map.reserve(seq_len);
        let at_end = chr_len.map_or(false, |len| region_end as usize == len);
        let site_end = if at_end {
            seq_len
        } else {
            seq_len.saturating_sub(seed1)
        };

        for i in 0..site_end {
            if at_end && i > seq_len.saturating_sub(seed1) {
                continue;
            }

            if i + seed1 <= seq_len {
                let key: [u8; 17] = self.ref_seq[i..i + seed1].try_into().unwrap();
                seed_map.insert_seed1(key, self.region_start + i as i64);
            }

            if i + seed2 <= seq_len {
                let key: [u8; 12] = self.ref_seq[i..i + seed2].try_into().unwrap();
                seed_map.insert_seed2(key, self.region_start + i as i64);
            }
        }
    }

    /// Clear the seed map after realignment to release memory early.
    pub fn clear_seed_map(&mut self) {
        if Arc::strong_count(&self.seed) == 1 {
            let seed_map = Arc::make_mut(&mut self.seed);
            seed_map.clear();
            seed_map.shrink_to_fit();
        } else {
            self.seed = Default::default();
        }
    }

    /// Get the base at a genomic position (0-based or 1-based depending on region_start)
    /// Returns None if position is out of bounds
    pub fn get(&self, genomic_pos: i64) -> Option<u8> {
        if genomic_pos < self.region_start {
            return None;
        }
        let idx = (genomic_pos - self.region_start) as usize;
        self.ref_seq.get(idx).copied()
    }

    /// Get the base at a genomic position as i64 (for compatibility)
    /// Returns None if position is out of bounds
    pub fn get_i64(&self, genomic_pos: i64) -> Option<u8> {
        self.get(genomic_pos)
    }

    /// Check if a genomic position contains a specific base
    pub fn has_and_equals(&self, genomic_pos: i64, base: u8) -> bool {
        self.get(genomic_pos).map_or(false, |b| b == base)
    }

    /// Check if a genomic position does NOT contain a specific base
    pub fn has_and_not_equals(&self, genomic_pos: i64, base: u8) -> bool {
        self.get(genomic_pos).map_or(false, |b| b != base)
    }
}

/// FASTA file reader with indexed access
pub struct FastaReader {
    reader: IndexedReader<File>,
}

impl FastaReader {
    /// Open a FASTA file with its index
    ///
    /// The .fai index must exist next to the FASTA file.
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let fasta_path = path.as_ref();
        let reader = IndexedReader::from_file(&fasta_path)
            .with_context(|| format!("Failed to open FASTA: {:?}", path.as_ref()))?;
        Ok(FastaReader { reader })
    }

    /// Fetch a reference sequence for a region
    ///
    /// Arguments:
    /// * `chrom` - Chromosome/contig name
    /// * `start` - Start position (1-based, inclusive)
    /// * `end` - End position (1-based, inclusive)
    ///
    /// Returns the sequence as uppercase bytes
    pub fn fetch_seq(&mut self, chrom: &str, start: usize, end: usize) -> Result<Vec<u8>> {
        if start == 0 {
            return Err(anyhow!("Start position must be 1-based (got 0)"));
        }
        if end < start {
            return Err(anyhow!(
                "End position ({}) cannot be less than start ({})",
                end,
                start
            ));
        }

        // Convert 1-based inclusive to 0-based half-open for IndexedReader.
        let begin_0based = start - 1;
        let stop_0based = end;

        self.reader
            .fetch(chrom, begin_0based as u64, stop_0based as u64)
            .with_context(|| format!("Failed to fetch {}:{}-{}", chrom, start, end))?;
        let mut seq = Vec::new();
        self.reader
            .read(&mut seq)
            .with_context(|| format!("Failed to read {}:{}-{}", chrom, start, end))?;

        for base in &mut seq {
            *base = base.to_ascii_uppercase();
        }

        Ok(seq)
    }

    /// Get the length of a chromosome/contig
    pub fn seq_len(&self, chrom: &str) -> usize {
        self.reader
            .index
            .sequences()
            .into_iter()
            .find(|sequence| sequence.name == chrom)
            .map(|sequence| sequence.len as usize)
            .unwrap_or(0)
    }

    /// Create a Reference struct with sequence for a region
    pub fn get_reference(&mut self, chrom: &str, start: usize, end: usize) -> Result<Reference> {
        let ref_seq = self.fetch_seq(chrom, start, end)?;
        let mut reference = Reference {
            ref_seq: Arc::new(ref_seq),
            seed: Default::default(),
            region_start: start as i64,
        };
        let chr_len = self.seq_len(chrom);
        reference.build_seed_map(end as i64, Some(chr_len));
        Ok(reference)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reference_default() {
        let r = Reference::default();
        assert!(r.ref_seq.is_empty());
        assert!(r.seed.is_empty());
    }

    #[test]
    fn test_build_seed_map_supports_slice_lookup_and_tracks_positions() {
        let mut reference = Reference::new_with_start(b"ACGTACGTACGTACGTACGT".to_vec(), 10);
        reference.build_seed_map(29, Some(29));

        let seed_17 = b"ACGTACGTACGTACGTA";
        let seed_12 = b"ACGTACGTACGT";

        let positions_17 = reference.seed.get(seed_17.as_slice()).expect("17-mer seed");
        let positions_12 = reference.seed.get(seed_12.as_slice()).expect("12-mer seed");

        assert!(positions_17.contains(&10));
        assert!(positions_12.contains(&10));
        assert!(positions_12.len() >= positions_17.len());
    }
}
