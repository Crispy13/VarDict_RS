#[derive(Debug, Clone)]
pub struct Region {
    chrom: String,
    start: usize,
    end: usize,
    gene: String,
    ins_start: usize,
    ins_end: usize,
}
