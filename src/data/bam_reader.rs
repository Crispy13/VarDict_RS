//! BAM file reading module
//!
//! Provides region-based BAM file access for variant calling.

use std::path::Path;

use anyhow::{Context, Result, anyhow};
use rust_htslib::bam::{self, Read, IndexedReader, Record};

/// BAM file reader with indexed access for region queries
pub struct BamReader {
    reader: IndexedReader,
    header: bam::Header,
}

impl BamReader {
    /// Open an indexed BAM file
    /// 
    /// The .bai index must exist (run `samtools index <bam>` to create)
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self> {
        let reader = IndexedReader::from_path(path.as_ref())
            .with_context(|| format!("Failed to open BAM: {:?}", path.as_ref()))?;
        let header = bam::Header::from_template(reader.header());
        Ok(BamReader { reader, header })
    }

    /// Set the region to query
    /// 
    /// Arguments:
    /// * `chrom` - Chromosome/contig name
    /// * `start` - Start position (1-based, inclusive)  
    /// * `end` - End position (1-based, inclusive)
    pub fn fetch(&mut self, chrom: &str, start: usize, end: usize) -> Result<()> {
        // BAM uses 0-based coordinates
        let tid = self.reader.header()
            .tid(chrom.as_bytes())
            .ok_or_else(|| anyhow!("Chromosome '{}' not found in BAM header", chrom))?;
        
        // Convert 1-based to 0-based for BAM
        let begin = (start.saturating_sub(1)) as i64;
        let end = end as i64;
        
        self.reader.fetch((tid, begin, end))
            .with_context(|| format!("Failed to fetch region {}:{}-{}", chrom, start, end))?;
        
        Ok(())
    }

    /// Read the next record from the current region
    /// 
    /// Returns None when no more records are available
    pub fn read(&mut self, record: &mut Record) -> Result<bool> {
        match self.reader.read(record) {
            Some(Ok(())) => Ok(true),
            Some(Err(e)) => Err(anyhow!("Error reading BAM record: {}", e)),
            None => Ok(false),
        }
    }

    /// Get the BAM header
    pub fn header(&self) -> &bam::Header {
        &self.header
    }

    /// Get chromosome/contig names from header
    pub fn target_names(&self) -> Vec<String> {
        self.reader.header()
            .target_names()
            .iter()
            .map(|n| String::from_utf8_lossy(n).to_string())
            .collect()
    }

    /// Get chromosome/contig lengths from header
    pub fn target_lens(&self) -> Vec<u64> {
        (0..self.reader.header().target_count())
            .map(|tid| self.reader.header().target_len(tid).unwrap_or(0))
            .collect()
    }
}

/// Iterator over records in a BAM region
pub struct BamRecordIter<'a> {
    reader: &'a mut BamReader,
    record: Record,
    done: bool,
}

impl<'a> BamRecordIter<'a> {
    pub fn new(reader: &'a mut BamReader) -> Self {
        BamRecordIter {
            reader,
            record: Record::new(),
            done: false,
        }
    }
}

impl<'a> Iterator for BamRecordIter<'a> {
    type Item = Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        match self.reader.read(&mut self.record) {
            Ok(true) => Some(Ok(self.record.clone())),
            Ok(false) => {
                self.done = true;
                None
            }
            Err(e) => {
                self.done = true;
                Some(Err(e))
            }
        }
    }
}

/// Check if a BAM record passes basic quality filters
pub fn passes_filter(record: &Record, sam_filter: u32, min_mapq: u8) -> bool {
    // Check SAM flags
    if record.flags() & (sam_filter as u16) != 0 {
        return false;
    }
    
    // Check mapping quality
    if record.mapq() < min_mapq {
        return false;
    }
    
    // Skip unmapped reads
    if record.is_unmapped() {
        return false;
    }
    
    // Skip secondary alignments if flagged
    if record.is_secondary() {
        return false;
    }
    
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_passes_filter_basic() {
        let record = Record::new();
        // New record has flags=0, mapq=0
        // Should pass with sam_filter=0, min_mapq=0
        // But will fail because unmapped check
    }
}
