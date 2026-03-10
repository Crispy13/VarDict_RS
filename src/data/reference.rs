use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result, anyhow};
use rust_htslib::faidx;

use crate::conf::Configuration;
use crate::prelude::LibDefaultHasher;
pub type ReferenceSeedMap = HashMap<Vec<u8>, Vec<i64>, LibDefaultHasher>;

/// Reference sequence data
#[derive(Default, Clone)]
pub struct Reference {
    pub ref_seq: Vec<u8>,
    pub seed: ReferenceSeedMap,
    /// Start position of this reference slice in genomic coordinates (1-based)
    pub region_start: i64,
}

impl Reference {
    /// Create a new Reference from a sequence slice
    pub fn from_seq(seq: &[u8]) -> Self {
        Reference {
            ref_seq: seq.to_vec(),
            seed: Default::default(),
            region_start: 0,
        }
    }

    /// Create a new Reference from a sequence slice with region start position
    pub fn from_seq_with_start(seq: &[u8], region_start: i64) -> Self {
        Reference {
            ref_seq: seq.to_vec(),
            seed: Default::default(),
            region_start,
        }
    }

    /// Create a new Reference with owned sequence
    pub fn new(ref_seq: Vec<u8>) -> Self {
        Reference {
            ref_seq,
            seed: Default::default(),
            region_start: 0,
        }
    }

    /// Create a new Reference with owned sequence and region start
    pub fn new_with_start(ref_seq: Vec<u8>, region_start: i64) -> Self {
        Reference {
            ref_seq,
            seed: Default::default(),
            region_start,
        }
    }

    /// Build the reference seed map using SEED_1 and SEED_2 lengths.
    pub fn build_seed_map(&mut self, region_end: i64, chr_len: Option<usize>) {
        self.seed.clear();

        if self.ref_seq.is_empty() {
            return;
        }

        let seed1 = Configuration::SEED_1 as usize;
        let seed2 = Configuration::SEED_2 as usize;
        let seq_len = self.ref_seq.len();
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
                let key = self.ref_seq[i..i + seed1].to_vec();
                self.seed
                    .entry(key)
                    .or_insert_with(Vec::new)
                    .push(self.region_start + i as i64);
            }

            if i + seed2 <= seq_len {
                let key = self.ref_seq[i..i + seed2].to_vec();
                self.seed
                    .entry(key)
                    .or_insert_with(Vec::new)
                    .push(self.region_start + i as i64);
            }
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
    reader: faidx::Reader,
}

impl FastaReader {
    /// Open a FASTA file with its index
    ///
    /// The .fai index must exist (run `samtools faidx <fasta>` to create)
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let reader = faidx::Reader::from_path(path.as_ref())
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
    pub fn fetch_seq(&self, chrom: &str, start: usize, end: usize) -> Result<Vec<u8>> {
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

        // Convert 1-based inclusive to 0-based half-open for faidx
        let begin_0based = start - 1;
        let end_0based = end - 1; // faidx uses 0-based inclusive end

        let seq = self
            .reader
            .fetch_seq(chrom, begin_0based, end_0based)
            .with_context(|| format!("Failed to fetch {}:{}-{}", chrom, start, end))?;

        // Convert to uppercase
        Ok(seq.into_iter().map(|b| b.to_ascii_uppercase()).collect())
    }

    /// Get the length of a chromosome/contig
    pub fn seq_len(&self, chrom: &str) -> usize {
        self.reader.fetch_seq_len(chrom) as usize
    }

    /// Create a Reference struct with sequence for a region
    pub fn get_reference(&self, chrom: &str, start: usize, end: usize) -> Result<Reference> {
        let ref_seq = self.fetch_seq(chrom, start, end)?;
        let mut reference = Reference {
            ref_seq,
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
