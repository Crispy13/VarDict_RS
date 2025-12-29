use std::sync::Arc;

use anyhow::{Error, anyhow};
use crackle_kit::{data::bases::rev_comp::RevComplementor, tracing::{Level, event}};
use rust_htslib::bam::{Record, Writer, record::CigarStringView};

use crate::{
    data::{reference::Reference, region::Region}, mods::cigar_modifier::CigarModifier,
    scopedata::global_read_only_scope::GlobalReadOnlyScope, utils::aligner::Aligner,
};

pub struct CigarParser {
    query_seq_buf: Vec<u8>,
    aligner: Aligner,
    instance: Arc<GlobalReadOnlyScope>,
    reference: Reference,
    max_read_len: usize,
    region: Region,
    rev_complementor: RevComplementor,
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

        let query_qual = record.qual();
        let is_mate_on_the_same_contig = record.tid() == record.mtid();
        let nm = tot_nm;
        let direction = record.is_reverse();

        if self.instance.amplicon_based_calling {
            todo!()
        }

        let mut pos = 0;
        let mut read_pos_including_softclip = 0;
        let mut read_pos_excluding_softclip = 0;

        if self.instance.conf.perform_local_realignment {
            // Modify the CIGAR for potential mis-alignment for indels at the end of reads to softclipping and let VarDict's
            // algorithm to figure out indels

            let cigar_modifier = CigarModifier::new(
                record.query_alignment_start(),
                cigar,
                query_seq,
                query_qual,
                &self.reference,
                ins_del_len,
                self.max_read_len,
                &self.region,
                &mut self.rev_complementor,
            );
        }

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
