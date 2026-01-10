use std::collections::{BTreeMap, HashMap};

use crackle_kit::nuc_base_map::NucBaseMap;
use smallvec::SmallVec;

use crate::prelude::SmallVecBytes;

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

/// Variant Description - used as key in variant maps and for tracking complex variants
/// 
/// Uses SmallVec for inline storage of short sequences (most variants are small)
/// This avoids heap allocation for the common case.
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub enum VarDesc {
    /// Single nucleotide variant: ref_base (alt determined from read data)
    /// Used as key for looking up variants at a position
    SNV { 
        ref_base: u8,
    },
    /// Insertion of sequence after position
    Ins {
        /// Inserted sequence
        seq: SmallVec<[u8; 32]>,
    },
    /// Deletion - complex structure for tracking parsing state
    Del {
        /// Length of deletion
        len: u32,
        /// Matched sequence following deletion
        match_seq: SmallVecBytes,
        /// Insertion or deletion length within complex deletion
        ins_or_del_len: InsOrDelLen,
        /// Mismatched sequence following deletion
        mismatch_seq: SmallVecBytes,
    },
    /// Complex variant (MNV or indel combination)
    Complex {
        /// Reference allele
        ref_seq: SmallVec<[u8; 32]>,
        /// Alternative allele
        alt_seq: SmallVec<[u8; 32]>,
    },
}

impl VarDesc {
    /// Create SNV from ref base (for use as HashMap key)
    pub fn snv_key(ref_base: u8) -> Self {
        VarDesc::SNV { 
            ref_base: ref_base.to_ascii_uppercase(),
        }
    }

    /// Create insertion
    pub fn insertion(seq: &[u8]) -> Self {
        VarDesc::Ins {
            seq: seq.iter().map(|b| b.to_ascii_uppercase()).collect(),
        }
    }

    /// Create deletion with tracking state
    pub fn deletion(len: u32) -> Self {
        VarDesc::Del {
            len,
            match_seq: SmallVecBytes::new(),
            ins_or_del_len: InsOrDelLen::None,
            mismatch_seq: SmallVecBytes::new(),
        }
    }

    /// Create complex variant
    pub fn complex(ref_seq: &[u8], alt_seq: &[u8]) -> Self {
        VarDesc::Complex {
            ref_seq: ref_seq.iter().map(|b| b.to_ascii_uppercase()).collect(),
            alt_seq: alt_seq.iter().map(|b| b.to_ascii_uppercase()).collect(),
        }
    }

    /// Get variant type as string for output
    pub fn variant_type(&self) -> &'static str {
        match self {
            VarDesc::SNV { .. } => "SNV",
            VarDesc::Ins { .. } => "Insertion",
            VarDesc::Del { .. } => "Deletion",
            VarDesc::Complex { .. } => "Complex",
        }
    }

    /// Get reference allele string
    pub fn ref_allele(&self) -> String {
        match self {
            VarDesc::SNV { ref_base, .. } => String::from(char::from(*ref_base)),
            VarDesc::Ins { .. } => String::new(), // Insertions have no ref (or context base)
            VarDesc::Del { len, .. } => format!("-{}", len),
            VarDesc::Complex { ref_seq, .. } => String::from_utf8_lossy(ref_seq).to_string(),
        }
    }

    /// Get alternative allele string
    pub fn alt_allele(&self) -> String {
        match self {
            VarDesc::SNV { ref_base, .. } => format!("?{}", char::from(*ref_base)), // Alt not stored
            VarDesc::Ins { seq } => format!("+{}", String::from_utf8_lossy(seq)),
            VarDesc::Del { len, .. } => format!("-{}", len),
            VarDesc::Complex { alt_seq, .. } => String::from_utf8_lossy(alt_seq).to_string(),
        }
    }

    /// Convert to variant key string (for backward compatibility)
    pub fn to_key_string(&self) -> String {
        match self {
            VarDesc::SNV { ref_base } => {
                // Only ref_base is stored, alt is determined from read data
                format!("{}",  char::from(*ref_base))
            }
            VarDesc::Ins { seq } => {
                format!("+{}", String::from_utf8_lossy(seq))
            }
            VarDesc::Del { len, .. } => {
                format!("-{}", len)
            }
            VarDesc::Complex { ref_seq, alt_seq } => {
                format!("{}>{}", 
                    String::from_utf8_lossy(ref_seq),
                    String::from_utf8_lossy(alt_seq))
            }
        }
    }
}

impl Default for VarDesc {
    fn default() -> Self {
        VarDesc::SNV { ref_base: b'N' }
    }
}

impl std::fmt::Display for VarDesc {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_key_string())
    }
}

/// Legacy enum for insertion/deletion length tracking (for complex patterns)
#[derive(Debug, Clone, Hash, PartialEq, Eq, Default)]
pub enum InsOrDelLen {
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