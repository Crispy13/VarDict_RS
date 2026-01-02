#[derive(Default)]
pub(crate) struct Variant {
    pub(crate) alt_depth: usize,
    pub(crate) alt_depth_fwd: usize,
    pub(crate) alt_depth_rev: usize,

    /// Sum of variant positions in read
    pub(crate) mean_pos: f64,

    /// Sum of base qualities for variant 
    pub(crate) mean_qual: f64,

    /// Sum of mapping qualities for variant
    pub(crate) mean_mapq: f64,

    /// Sum of number of mismatches for variant
    pub(crate) nm: f64,

    /// Number of low-quality reads with the variant
    pub(crate) low_qual_read_cnt: usize,

    /// Number of high-quality reads with the variant
    pub(crate) high_qual_read_cnt: usize,


}

impl Variant {
    pub(crate) fn is_empty(&self) -> bool {
        todo!()
    }

    pub(crate) fn inc_dir(&mut self, is_reverse:bool) {
        if is_reverse { 
            self.alt_depth_rev += 1;
        } else {
            self.alt_depth_fwd += 1;
        }
    }
}


/// Variant Information
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub(crate) enum VarDesc {
    SNV{ ref_base: u8 },
}

impl VarDesc {
    
}
