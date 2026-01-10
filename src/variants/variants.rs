use std::{borrow::Cow, collections::{BTreeMap, HashMap}};

use crackle_kit::nuc_base_map::NucBaseMap;
use smallvec::SmallVec;

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

    /// Flag: true if variant is covered by reads with different positions
    pub(crate) pstd: bool,

    /// Flag: true if variant is covered by reads with different qualities  
    pub(crate) qstd: bool,

    /// Previous position (for pstd calculation)
    pub(crate) pp: usize,

    /// Previous quality (for qstd calculation)
    pub(crate) pq: f64,
}

impl Variant {
    pub(crate) fn is_empty(&self) -> bool {
        todo!()
    }

    pub(crate) fn inc_dir(&mut self, is_reverse: bool) {
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
    SNV { ref_base: u8 },
    Del {
        /// Length of deletion
        len: u32,
        /// Matched sequence after deletion (from D+M+I/D pattern) - the '#' part
        match_seq: SmallVec<[u8; 32]>,
        /// For D+M+I: insertion sequence; for D+M+D: deletion length - the '^' part
        ins_or_del_len: InsOrDelLen,
        /// Mismatched sequence to append (the '&' part)
        mismatch_seq: SmallVec<[u8; 32]>,
    }
}

impl VarDesc {}

#[derive(Debug, Clone, Hash, PartialEq, Eq, Default)]
pub(crate) enum InsOrDelLen {
    #[default]
    None,
    InsSeq(SmallVec<[u8;32]>),
    DelLen(usize),
}

#[derive(Default)]
pub(crate) struct SoftClip {
    pub(crate) var: Variant,

    /// Map of position(offset) of high quality base in the base sequence (from SAM record)
    /// to base on this position and it's count
    pub(crate) nt: BTreeMap<i64, NucBaseMap<usize>>,

    /// Map of position of high quality base in the base sequence (from SAM record)
    /// to base on this position and it's variation
    pub(crate) seq: BTreeMap<usize, NucBaseMap<Variant>>,

    /// The consensus sequence in soft-clipped reads.
    consensus_seq: Vec<u8>,

    used: bool,

    /// Additional fields for structural variation
    start: i64,
    end: i64,
    mstart: i64,
    mend: i64,
    mlen: i32,
    disc: i32,
    softp: i32,

    /// Map of softclip positions to their counts for SV
    soft: HashMap<i64, usize>,
    mates: Vec<Mate>,
}


struct Mate {
    // TODO: for SV.
}