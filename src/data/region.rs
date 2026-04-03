#[derive(Debug, Clone, Default)]
pub struct Region {
    pub(crate) chrom: String,
    pub(crate) start: usize,
    pub(crate) end: usize,
    pub(crate) display_start: i64,
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
            display_start: start as i64,
            gene,
            ins_start: 0,
            ins_end: 0,
        }
    }

    /// Create a new region with distinct display and computational starts.
    pub fn new_extended(
        chrom: String,
        start: usize,
        end: usize,
        gene: String,
        display_start: i64,
    ) -> Self {
        Self {
            chrom,
            start,
            end,
            display_start,
            gene,
            ins_start: 0,
            ins_end: 0,
        }
    }

    /// Create a new region with amplicon insert interval coordinates.
    pub fn new_with_insert(
        chrom: String,
        start: usize,
        end: usize,
        gene: String,
        ins_start: usize,
        ins_end: usize,
    ) -> Self {
        Self {
            chrom,
            start,
            end,
            display_start: start as i64,
            gene,
            ins_start,
            ins_end,
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

    /// Get the display start position, which may be negative for extended regions.
    pub fn display_start(&self) -> i64 {
        self.display_start
    }

    /// Get the gene name
    pub fn gene(&self) -> &str {
        &self.gene
    }

    /// Get amplicon insert start position (1-based)
    pub fn insert_start(&self) -> usize {
        self.ins_start
    }

    /// Get amplicon insert end position (1-based, inclusive)
    pub fn insert_end(&self) -> usize {
        self.ins_end
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
        format!("{}:{}-{}", self.chrom, self.display_start, self.end)
    }
}

#[cfg(test)]
mod tests {
    use super::Region;

    #[test]
    fn new_defaults_display_start_to_start() {
        let region = Region::new("chr1".to_string(), 100, 200, "GENE".to_string());

        assert_eq!(region.display_start(), 100);
        assert_eq!(region.to_region_string(), "chr1:100-200");
    }

    #[test]
    fn new_extended_preserves_negative_display_start() {
        let region = Region::new_extended("20".to_string(), 1, 1_000_150, "20".to_string(), -149);

        assert_eq!(region.start(), 1);
        assert_eq!(region.display_start(), -149);
        assert_eq!(region.to_region_string(), "20:-149-1000150");
    }
}
