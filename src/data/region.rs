#[derive(Debug, Clone, Default)]
pub struct Region {
    pub(crate) chrom: String,
    pub(crate) start: usize,
    pub(crate) end: usize,
    pub(crate) gene: String,
    pub(crate) ins_start: usize,
    pub(crate) ins_end: usize,
}

impl Region {
    /// Create a new region
    pub fn new(chrom: String, start: usize, end: usize, gene: String) -> Self {
        Self {
            chrom,
            start,
            end,
            gene,
            ins_start: 0,
            ins_end: 0,
        }
    }

    /// Get the chromosome name
    pub fn chr(&self) -> &str {
        &self.chrom
    }

    /// Get the start position (1-based)
    pub fn start(&self) -> usize {
        self.start
    }

    /// Get the end position (1-based, inclusive)
    pub fn end(&self) -> usize {
        self.end
    }

    /// Get the gene name
    pub fn gene(&self) -> &str {
        &self.gene
    }

    /// Get the region length
    pub fn len(&self) -> usize {
        if self.end >= self.start {
            self.end - self.start + 1
        } else {
            0
        }
    }

    /// Check if region is empty
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Format as "chr:start-end"
    pub fn to_region_string(&self) -> String {
        format!("{}:{}-{}", self.chrom, self.start, self.end)
    }
}
