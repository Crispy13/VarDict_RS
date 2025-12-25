use anyhow::{Error, anyhow};
use crackle_kit::tracing::{Level, event};
use rust_htslib::bam::{Record, Writer, record::CigarStringView};

use crate::utils::aligner::Aligner;

pub struct CigarParser {
    query_seq_buf: Vec<u8>,
    aligner: Aligner,
}

impl CigarParser {
    fn parse_cigar(&mut self, record: &mut Record) -> Result<(), Error> {
        self.query_seq_buf
            .extend(record.seq().into_decoded_base_iter());

        let query_seq = self.query_seq_buf.as_slice();
        let mapping_quality = record.mapq();

        let cigar = get_cached_cigar_or_make!(record);

        let ins_del_len = get_ins_del_len(cigar);

        let tot_nm = match record.aux_option(self.aligner.nm_tag())? {
            Some(rust_htslib::bam::record::Aux::I32(nm)) => nm - ins_del_len as i32,
            Some(oth) => Err(anyhow!("Got non i32 type for NM tag: {:?}", oth))?,
            None => {
                if !cigar.is_empty() {
                    event!(Level::WARN, "No NM tag for mismatches: {:?}", record);
                }

                if record.is_unmapped() || cigar.is_empty() {
                    return Ok(());
                }

                0
            }
        };

        todo!()
    }
}

macro_rules! get_cached_cigar_or_make {
    ($record:expr) => {{
        if $record.cigar_cached().is_none() {
            $record.cache_cigar();
        }

        $record.cigar_cached().unwrap()
    }};
}
use get_cached_cigar_or_make;

#[inline]
fn get_ins_del_len(cigar: &CigarStringView) -> u32 {
    cigar
        .iter()
        .filter_map(|c| match c {
            rust_htslib::bam::record::Cigar::Ins(l) | rust_htslib::bam::record::Cigar::Del(l) => {
                Some(*l)
            }
            _ => None,
        })
        .sum::<u32>()
}
