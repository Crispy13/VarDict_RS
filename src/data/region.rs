#[derive(Debug, Clone, Default)]
pub struct Region {
    pub(crate) chrom: String,
    pub(crate) start: usize,
    pub(crate) end: usize,
    pub(crate) gene: String,
    pub(crate) ins_start: usize,
    pub(crate) ins_end: usize,
}
