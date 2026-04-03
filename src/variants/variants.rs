use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::fmt::Write;

use crackle_kit::nuc_base_map::NucBaseMap;
use indexmap::IndexMap;
use smallvec::SmallVec;

use crate::prelude::SmallVecBytes;

const MAX_SEGMENTS: usize = 8;

fn format_u32_stack(n: u32, buf: &mut [u8; 10]) -> &[u8] {
    let mut value = n;
    let mut index = buf.len();

    loop {
        index -= 1;
        buf[index] = b'0' + (value % 10) as u8;
        value /= 10;
        if value == 0 {
            break;
        }
    }

    &buf[index..]
}

fn format_usize_stack(n: usize, buf: &mut [u8; 20]) -> &[u8] {
    let mut value = n;
    let mut index = buf.len();

    loop {
        index -= 1;
        buf[index] = b'0' + (value % 10) as u8;
        value /= 10;
        if value == 0 {
            break;
        }
    }

    &buf[index..]
}

#[derive(Default, Debug, Clone)]
pub struct Variant {
    pub alt_depth: usize,
    pub alt_depth_fwd: usize,
    pub alt_depth_rev: usize,

    /// Adjusted count for indels due to realignment (Java: extracnt)
    pub extra_cnt: usize,

    /// Sum of variant positions in read
    pub mean_pos: f64,

    /// Sum of base qualities for variant
    pub mean_qual: f64,

    /// Sum of mapping qualities for variant
    pub mean_mapq: f64,

    /// Sum of number of mismatches for variant
    pub nm: f64,

    /// Number of low-quality reads with the variant
    pub low_qual_read_cnt: usize,

    /// Number of high-quality reads with the variant
    pub high_qual_read_cnt: usize,

    /// Flag: true if variant is covered by reads with different positions
    pub pstd: bool,

    /// Flag: true if variant is covered by reads with different qualities  
    pub qstd: bool,

    /// Previous position (for pstd calculation)
    pub pp: usize,

    /// Previous quality (for qstd calculation)
    pub pq: f64,
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

#[derive(Default, Debug, Clone, Copy, PartialEq, Eq)]
pub struct StructuralVariantCounts {
    pub splits: usize,
    pub pairs: usize,
    pub clusters: usize,
}

/// Variant Description - used as key in variant maps and for tracking complex variants
///
/// Uses SmallVec for inline storage of short sequences (most variants are small)
/// This avoids heap allocation for the common case.
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub enum VarDesc {
    /// Single nucleotide variant: ref_base (alt determined from read data)
    /// Used as key for looking up variants at a position
    SNV { ref_base: u8 },
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
    /// Raw description string (Java-style), e.g. "A&TGC", "-3&AT", etc.
    Raw {
        /// Raw description bytes
        desc: SmallVecBytes,
    },
}

#[derive(Clone, Copy)]
struct KeySegments<'a> {
    segments: [&'a [u8]; MAX_SEGMENTS],
    len: usize,
}

impl<'a> KeySegments<'a> {
    fn new() -> Self {
        Self {
            segments: [&[]; MAX_SEGMENTS],
            len: 0,
        }
    }

    fn push(&mut self, segment: &'a [u8]) {
        if segment.is_empty() {
            return;
        }

        debug_assert!(self.len < MAX_SEGMENTS);
        self.segments[self.len] = segment;
        self.len += 1;
    }

    fn flat_iter(&self) -> impl Iterator<Item = u8> + '_ {
        self.segments[..self.len]
            .iter()
            .flat_map(|segment| segment.iter().copied())
    }
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
            VarDesc::Raw { .. } => "Raw",
        }
    }

    /// Get reference allele string
    pub fn ref_allele(&self) -> String {
        match self {
            VarDesc::SNV { ref_base, .. } => String::from(char::from(*ref_base)),
            VarDesc::Ins { .. } => String::new(), // Insertions have no ref (or context base)
            VarDesc::Del { len, .. } => format!("-{}", len),
            VarDesc::Complex { ref_seq, .. } => String::from_utf8_lossy(ref_seq).to_string(),
            VarDesc::Raw { .. } => String::new(),
        }
    }

    /// Get alternative allele string
    pub fn alt_allele(&self) -> String {
        match self {
            VarDesc::SNV { ref_base, .. } => format!("?{}", char::from(*ref_base)), // Alt not stored
            VarDesc::Ins { seq } => format!("+{}", String::from_utf8_lossy(seq)),
            VarDesc::Del { len, .. } => format!("-{}", len),
            VarDesc::Complex { alt_seq, .. } => String::from_utf8_lossy(alt_seq).to_string(),
            VarDesc::Raw { desc } => String::from_utf8_lossy(desc).to_string(),
        }
    }

    /// Convert to variant key string (for backward compatibility)
    pub fn to_key_string(&self) -> String {
        use std::fmt::Write;
        match self {
            VarDesc::SNV { ref_base } => String::from(*ref_base as char),
            VarDesc::Ins { seq } => {
                let mut s = String::with_capacity(1 + seq.len());
                s.push('+');
                // Genomic sequences are ASCII-safe (project convention)
                s.push_str(std::str::from_utf8(seq).unwrap());
                s
            }
            VarDesc::Del {
                len,
                match_seq,
                ins_or_del_len,
                mismatch_seq,
            } => {
                let mut s = String::with_capacity(4 + match_seq.len() + mismatch_seq.len());
                s.push('-');
                write!(s, "{len}").unwrap();
                if !match_seq.is_empty() {
                    s.push('#');
                    s.push_str(std::str::from_utf8(match_seq).unwrap());
                }
                match ins_or_del_len {
                    InsOrDelLen::InsSeq(seq) => {
                        s.push('^');
                        s.push_str(std::str::from_utf8(seq).unwrap());
                    }
                    InsOrDelLen::DelLen(n) => {
                        s.push('^');
                        write!(s, "{n}").unwrap();
                    }
                    InsOrDelLen::None => {}
                }
                if !mismatch_seq.is_empty() {
                    s.push('&');
                    s.push_str(std::str::from_utf8(mismatch_seq).unwrap());
                }
                s
            }
            VarDesc::Complex { ref_seq, alt_seq } => {
                let mut s = String::with_capacity(ref_seq.len() + 1 + alt_seq.len());
                s.push_str(std::str::from_utf8(ref_seq).unwrap());
                s.push('>');
                s.push_str(std::str::from_utf8(alt_seq).unwrap());
                s
            }
            VarDesc::Raw { desc } => String::from_utf8_lossy(desc).into_owned(),
        }
    }

    fn key_segments<'a>(
        &'a self,
        u32_buf: &'a mut [u8; 10],
        usize_buf: &'a mut [u8; 20],
    ) -> KeySegments<'a> {
        let mut segments = KeySegments::new();

        match self {
            VarDesc::SNV { ref_base } => {
                segments.push(std::slice::from_ref(ref_base));
            }
            VarDesc::Ins { seq } => {
                segments.push(b"+");
                segments.push(seq.as_slice());
            }
            VarDesc::Del {
                len,
                match_seq,
                ins_or_del_len,
                mismatch_seq,
            } => {
                segments.push(b"-");
                segments.push(format_u32_stack(*len, u32_buf));

                if !match_seq.is_empty() {
                    segments.push(b"#");
                    segments.push(match_seq.as_slice());
                }

                match ins_or_del_len {
                    InsOrDelLen::InsSeq(seq) => {
                        segments.push(b"^");
                        segments.push(seq.as_slice());
                    }
                    InsOrDelLen::DelLen(n) => {
                        segments.push(b"^");
                        segments.push(format_usize_stack(*n, usize_buf));
                    }
                    InsOrDelLen::None => {}
                }

                if !mismatch_seq.is_empty() {
                    segments.push(b"&");
                    segments.push(mismatch_seq.as_slice());
                }
            }
            VarDesc::Complex { ref_seq, alt_seq } => {
                segments.push(ref_seq.as_slice());
                segments.push(b">");
                segments.push(alt_seq.as_slice());
            }
            VarDesc::Raw { desc } => {
                segments.push(desc.as_slice());
            }
        }

        segments
    }

    pub fn cmp_as_key_string(&self, other: &Self) -> Ordering {
        let mut self_u32_buf = [0u8; 10];
        let mut self_usize_buf = [0u8; 20];
        let mut other_u32_buf = [0u8; 10];
        let mut other_usize_buf = [0u8; 20];

        let self_segments = self.key_segments(&mut self_u32_buf, &mut self_usize_buf);
        let other_segments = other.key_segments(&mut other_u32_buf, &mut other_usize_buf);

        self_segments.flat_iter().cmp(other_segments.flat_iter())
    }

    /// Returns true iff `self.to_key_string() == String::from(b as char)`.
    #[inline]
    pub fn key_matches_byte(&self, b: u8) -> bool {
        match self {
            VarDesc::SNV { ref_base } => *ref_base == b,
            VarDesc::Raw { desc } => desc.len() == 1 && desc[0] == b,
            _ => false,
        }
    }

    /// Returns true iff the key string would start with '-'.
    /// Zero-allocation equivalent of checking whether the key begins with '-'.
    #[inline]
    pub fn is_deletion(&self) -> bool {
        match self {
            VarDesc::Del { .. } => true,
            VarDesc::Raw { desc } => desc.first() == Some(&b'-'),
            _ => false,
        }
    }

    /// Returns true iff `to_key_string() == s`.
    /// Zero-allocation comparison using key_segments infrastructure.
    pub fn key_equals(&self, s: &str) -> bool {
        let mut u32_buf = [0u8; 10];
        let mut usize_buf = [0u8; 20];
        let segments = self.key_segments(&mut u32_buf, &mut usize_buf);
        let s_bytes = s.as_bytes();
        let total_len: usize = segments.segments[..segments.len]
            .iter()
            .map(|segment| segment.len())
            .sum();

        if total_len != s_bytes.len() {
            return false;
        }

        let mut offset = 0;
        for &segment in &segments.segments[..segments.len] {
            if &s_bytes[offset..offset + segment.len()] != segment {
                return false;
            }
            offset += segment.len();
        }

        true
    }
}

impl Default for VarDesc {
    fn default() -> Self {
        VarDesc::SNV { ref_base: b'N' }
    }
}

impl std::fmt::Display for VarDesc {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            VarDesc::SNV { ref_base } => f.write_char(*ref_base as char),
            VarDesc::Ins { seq } => {
                f.write_char('+')?;
                f.write_str(std::str::from_utf8(seq).unwrap())
            }
            VarDesc::Del {
                len,
                match_seq,
                ins_or_del_len,
                mismatch_seq,
            } => {
                write!(f, "-{len}")?;
                if !match_seq.is_empty() {
                    f.write_char('#')?;
                    f.write_str(std::str::from_utf8(match_seq).unwrap())?;
                }
                match ins_or_del_len {
                    InsOrDelLen::InsSeq(seq) => {
                        f.write_char('^')?;
                        f.write_str(std::str::from_utf8(seq).unwrap())?;
                    }
                    InsOrDelLen::DelLen(n) => write!(f, "^{n}")?,
                    InsOrDelLen::None => {}
                }
                if !mismatch_seq.is_empty() {
                    f.write_char('&')?;
                    f.write_str(std::str::from_utf8(mismatch_seq).unwrap())?;
                }
                Ok(())
            }
            VarDesc::Complex { ref_seq, alt_seq } => {
                f.write_str(std::str::from_utf8(ref_seq).unwrap())?;
                f.write_char('>')?;
                f.write_str(std::str::from_utf8(alt_seq).unwrap())
            }
            VarDesc::Raw { desc } => {
                write!(f, "{}", String::from_utf8_lossy(desc))
            }
        }
    }
}

/// Legacy enum for insertion/deletion length tracking (for complex patterns)
#[derive(Debug, Clone, Hash, PartialEq, Eq, Default)]
pub enum InsOrDelLen {
    #[default]
    None,
    InsSeq(SmallVec<[u8; 32]>),
    DelLen(usize),
}

#[derive(Default, Clone)]
pub(crate) struct SoftClip {
    pub(crate) var: Variant,

    /// Map of position(offset) of high quality base in the base sequence (from SAM record)
    /// to base on this position and it's count
    pub(crate) nt: BTreeMap<i64, NucBaseMap<usize>>,

    /// Map of position of high quality base in the base sequence (from SAM record)
    /// to base on this position and it's variation
    pub(crate) seq: BTreeMap<usize, NucBaseMap<Variant>>,

    /// The consensus sequence in soft-clipped reads.
    ///
    /// Java caches even empty consensus (sequence = "").
    /// `None` means not computed yet; `Some(vec![])` means computed empty.
    consensus_seq: Option<Vec<u8>>,

    used: bool,

    /// Additional fields for structural variation
    pub(crate) start: i64,
    pub(crate) end: i64,
    pub(crate) mstart: i64,
    pub(crate) mend: i64,
    pub(crate) mlen: i32,
    pub(crate) disc: i32,
    pub(crate) softp: i32,

    /// Map of softclip positions to their counts for SV.
    /// Java uses LinkedHashMap here, so insertion order is significant.
    pub(crate) soft: IndexMap<i64, usize>,
    pub(crate) mates: Vec<Mate>,
}

#[derive(Default, Clone)]
pub(crate) struct Mate {
    pub(crate) mate_start: i64,
    pub(crate) mate_end: i64,
    pub(crate) mate_len: i32,
    pub(crate) start: i64,
    pub(crate) end: i64,
    pub(crate) mean_pos: f64,
    pub(crate) mean_qual: f64,
    pub(crate) mean_mapq: f64,
    pub(crate) nm: f64,
}

impl SoftClip {
    /// Get the consensus sequence
    pub fn consensus_seq(&self) -> &[u8] {
        self.consensus_seq.as_deref().unwrap_or(&[])
    }

    /// Returns true when consensus has already been computed,
    /// including an explicitly empty consensus.
    pub fn consensus_seq_is_set(&self) -> bool {
        self.consensus_seq.is_some()
    }

    /// Set the consensus sequence
    pub fn set_consensus_seq(&mut self, seq: Vec<u8>) {
        self.consensus_seq = Some(seq);
    }

    /// Check if this soft clip has been used
    pub fn used(&self) -> bool {
        self.used
    }

    /// Mark this soft clip as used
    pub fn mark_used(&mut self) {
        self.used = true;
    }
}

#[cfg(test)]
mod tests {
    use super::{InsOrDelLen, VarDesc, format_u32_stack, format_usize_stack};
    use crate::prelude::SmallVecBytes;

    fn raw(desc: &[u8]) -> VarDesc {
        VarDesc::Raw {
            desc: SmallVecBytes::from_slice(desc),
        }
    }

    fn del(
        len: u32,
        match_seq: &[u8],
        ins_or_del_len: InsOrDelLen,
        mismatch_seq: &[u8],
    ) -> VarDesc {
        VarDesc::Del {
            len,
            match_seq: SmallVecBytes::from_slice(match_seq),
            ins_or_del_len,
            mismatch_seq: SmallVecBytes::from_slice(mismatch_seq),
        }
    }

    fn ins_seq(seq: &[u8]) -> InsOrDelLen {
        InsOrDelLen::InsSeq(seq.iter().copied().collect())
    }

    #[test]
    fn test_key_matches_byte_parity_with_to_key_string() {
        let cases: Vec<(VarDesc, u8)> = vec![
            (VarDesc::snv_key(b'A'), b'A'),
            (VarDesc::snv_key(b'A'), b'C'),
            (VarDesc::snv_key(b'T'), b'T'),
            (VarDesc::insertion(b"ACG"), b'A'),
            (VarDesc::deletion(3), b'-'),
            (VarDesc::complex(b"AT", b"GC"), b'A'),
            (
                VarDesc::Raw {
                    desc: SmallVecBytes::from_slice(b"A"),
                },
                b'A',
            ),
            (
                VarDesc::Raw {
                    desc: SmallVecBytes::from_slice(b"A"),
                },
                b'C',
            ),
            (
                VarDesc::Raw {
                    desc: SmallVecBytes::from_slice(b"AT"),
                },
                b'A',
            ),
        ];

        for (desc, byte) in &cases {
            let string_match = desc.to_key_string() == String::from(*byte as char);
            let byte_match = desc.key_matches_byte(*byte);
            assert_eq!(
                string_match, byte_match,
                "Mismatch for {:?} vs byte {}: string={}, byte={}",
                desc, *byte as char, string_match, byte_match
            );
        }
    }

    #[test]
    fn test_is_deletion_parity_with_to_key_string_starts_with() {
        let cases: Vec<VarDesc> = vec![
            VarDesc::snv_key(b'A'),
            VarDesc::insertion(b"ACG"),
            VarDesc::deletion(3),
            VarDesc::complex(b"AT", b"GC"),
            raw(b"-3"),
            raw(b"-3#ACG^TT&A"),
            raw(b"A>T"),
            raw(b""),
        ];

        for desc in &cases {
            let key_string = desc.to_key_string();
            assert_eq!(
                key_string.starts_with('-'),
                desc.is_deletion(),
                "Mismatch for {:?}",
                desc
            );
        }
    }

    #[test]
    fn test_key_equals_parity_with_to_key_string_equality() {
        let cases: Vec<(VarDesc, &[&str])> = vec![
            (VarDesc::snv_key(b'A'), &["A", "C", "+A"]),
            (VarDesc::insertion(b"ACG"), &["+ACG", "+AC", "ACG"]),
            (
                del(12, b"GG", InsOrDelLen::DelLen(45), b"TT"),
                &["-12#GG^45&TT", "-12#GG^45", "-12"],
            ),
            (VarDesc::complex(b"AT", b"GC"), &["AT>GC", "AT", "GC"]),
            (raw(b"-3#ACG^TT&A"), &["-3#ACG^TT&A", "-3#ACG", "A>T"]),
            (raw(b"SV"), &["SV", "sv", "A"]),
        ];

        for (desc, candidates) in &cases {
            for candidate in *candidates {
                assert_eq!(
                    desc.to_key_string() == *candidate,
                    desc.key_equals(candidate),
                    "Mismatch for {:?} vs {candidate}",
                    desc
                );
            }
        }
    }

    #[test]
    fn test_format_u32_stack_round_trip() {
        let values = [0_u32, 1, 9, 10, 11, 99, 100, 999, 1_000, u32::MAX];

        for value in values {
            let mut buf = [0u8; 10];
            let formatted = format_u32_stack(value, &mut buf);
            assert_eq!(std::str::from_utf8(formatted).unwrap(), value.to_string());
        }
    }

    #[test]
    fn test_format_usize_stack_round_trip() {
        let values = [
            0_usize,
            1,
            9,
            10,
            11,
            99,
            100,
            999,
            1_000,
            12_345,
            usize::MAX,
        ];

        for value in values {
            let mut buf = [0u8; 20];
            let formatted = format_usize_stack(value, &mut buf);
            assert_eq!(std::str::from_utf8(formatted).unwrap(), value.to_string());
        }
    }

    #[test]
    fn test_cmp_as_key_string_matches_string_ordering() {
        let cases: Vec<VarDesc> = vec![
            VarDesc::snv_key(b'A'),
            VarDesc::snv_key(b'C'),
            VarDesc::snv_key(b'G'),
            VarDesc::snv_key(b'N'),
            VarDesc::snv_key(b'T'),
            VarDesc::insertion(b""),
            VarDesc::insertion(b"A"),
            VarDesc::insertion(b"AC"),
            VarDesc::insertion(b"ACG"),
            VarDesc::insertion(b"TT"),
            VarDesc::insertion(b"NNNN"),
            VarDesc::insertion(b"ACGTACGT"),
            del(0, b"", InsOrDelLen::None, b""),
            del(1, b"", InsOrDelLen::None, b""),
            del(2, b"A", InsOrDelLen::None, b""),
            del(3, b"AC", InsOrDelLen::None, b"G"),
            del(4, b"", ins_seq(b"T"), b""),
            del(5, b"G", ins_seq(b"TT"), b"A"),
            del(6, b"", InsOrDelLen::DelLen(0), b""),
            del(7, b"CC", InsOrDelLen::DelLen(3), b""),
            del(8, b"", InsOrDelLen::DelLen(12), b"T"),
            del(9, b"AA", InsOrDelLen::DelLen(45), b"GG"),
            del(10, b"TT", ins_seq(b""), b""),
            del(11, b"", ins_seq(b"ACG"), b"TT"),
            del(12, b"GG", ins_seq(b"TTA"), b"C"),
            del(99, b"A", InsOrDelLen::None, b"TT"),
            del(123, b"ACG", InsOrDelLen::DelLen(6789), b"TG"),
            del(u32::MAX, b"T", ins_seq(b"G"), b"A"),
            VarDesc::complex(b"", b"A"),
            VarDesc::complex(b"A", b""),
            VarDesc::complex(b"A", b"T"),
            VarDesc::complex(b"AT", b"GC"),
            VarDesc::complex(b"G", b"GA"),
            VarDesc::complex(b"ACGT", b"T"),
            VarDesc::complex(b"NN", b"NNN"),
            VarDesc::complex(b"TTAA", b"CCGG"),
            raw(b""),
            raw(b"A"),
            raw(b"+"),
            raw(b"+ACG"),
            raw(b"-3"),
            raw(b"-3#ACG"),
            raw(b"-3#ACG^TT&A"),
            raw(b"-12^34"),
            raw(b"A>T"),
            raw(b">"),
            raw(b"#"),
            raw(b"&"),
            raw(b"^12"),
            raw(b"Z"),
            raw(b"a"),
            raw(b"TT"),
            raw(b"SV"),
            raw(b"N"),
        ];

        for (i, a) in cases.iter().enumerate() {
            for (j, b) in cases.iter().enumerate() {
                let str_ord = a.to_key_string().cmp(&b.to_key_string());
                let seg_ord = a.cmp_as_key_string(b);
                assert_eq!(
                    str_ord, seg_ord,
                    "Mismatch for cases[{i}] vs cases[{j}]: {:?} vs {:?}",
                    a, b
                );
            }
        }
    }
}
