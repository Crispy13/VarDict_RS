use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result, anyhow};
use rust_htslib::faidx;

/// Reference sequence data
#[derive(Default)]
pub(crate) struct Reference {
    pub(crate) ref_seq: Vec<u8>,
    pub(crate) seed: HashMap<Vec<u8>, Vec<i64>>,
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
            return Err(anyhow!("End position ({}) cannot be less than start ({})", end, start));
        }

        // Convert 1-based inclusive to 0-based half-open for faidx
        let begin_0based = start - 1;
        let end_0based = end - 1; // faidx uses 0-based inclusive end

        let seq = self.reader.fetch_seq(chrom, begin_0based, end_0based)
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
        Ok(Reference {
            ref_seq,
            seed: HashMap::new(), // Seeds computed separately if needed
        })
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
}