#[derive(Clone)]
pub struct Configuration {
    pub perform_local_realignment: bool,

    /// For downsampling fraction (Java: -Z)
    pub downsampling: Option<f64>,

    /// Indicate to turn off chimeric reads filtering.  Chimeric reads are artifacts from library construction,
    /// where a read can be split into two segments, each will be aligned within 1-2 read length distance,
    /// but in opposite direction.
    pub chimeric_filter: bool,
    /// Indicate to remove duplicated reads (Java: -t)
    pub remove_duplicated_reads: bool,
    pub min_match: i32,

    pub disable_sv: bool,

    pub unique_mode_alignment_enabled: bool,
    pub unique_mode_second_in_pair_enabled: bool,

    /// The hexical to filter reads.
    pub sam_filter: u32,

    /// Trim bases after this position (Java: trimBasesAfter). 0 disables.
    pub trim_bases_after: i32,

    pub crispr_cutting_site: i32,
    pub crispr_filtering_bp: i32,

    /// Count N bases in total depth (Java: includeNInTotalDepth)
    pub include_n_in_total_depth: bool,
    // pub seed_2: i32,

    /// The phred score for a base to be considered a good call
    pub goodq: f64,

    /// Number of nucleotides to extend regions (Java: -x)
    pub number_nucleotide_to_extend: i32,

    /// Reference extension for fetching sequence (Java: -Y)
    pub reference_extension: i32,

    /// The threshold for allele frequency. If -p it is set to -1.
    pub freq: f64,

    /// Minimum number of variant reads
    pub minr: usize,

    /// The read position filter
    pub read_pos_filter: f64,

    /// The Qratio of (good_quality_reads)/(bad_quality_reads+0.5)
    pub qratio: f64,

    /// Mean mapping quality threshold
    pub mapq: f64,

    /// The variant frequency threshold to determine variant as good in case of monomer MSI
    pub monomer_msi_frequency: f64,

    /// The variant frequency threshold to determine variant as good in case of non-monomer MSI
    pub non_monomer_msi_frequency: f64,

    /// Extension of bp to look for mismatches after insertion or deletion
    pub vext: i32,

    /// If set, reads with mismatches more than INT will be filtered and ignored
    pub mismatch: i32,

    /// If set, reads with mapping quality less than INT will be filtered and ignored (Java: -Q)
    pub mapping_quality: Option<u8>,

    /// Move indels to 3' end (Java: moveIndelsTo3 / -3)
    pub move_indels_to_3: bool,

    /// Minimum length for structural variants (Java: SVMINLEN / -L)
    pub sv_min_len: usize,
}

impl Default for Configuration {
    fn default() -> Self {
        Self {
            perform_local_realignment: Default::default(),
            downsampling: None,
            chimeric_filter: false,
            remove_duplicated_reads: false,
            min_match: 0,
            sam_filter: 0x504,
            trim_bases_after: 0,
            crispr_cutting_site: 0,
            crispr_filtering_bp: 0,
            include_n_in_total_depth: false,
            disable_sv: false,
            unique_mode_alignment_enabled: false,
            unique_mode_second_in_pair_enabled: false,
            goodq: 22.5,
            number_nucleotide_to_extend: 0,
            reference_extension: 1200,
            freq: 0.01,
            minr: 2,
            read_pos_filter: 5.0,
            qratio: 1.5,
            mapq: 0.0,
            monomer_msi_frequency: 0.25,
            non_monomer_msi_frequency: 0.1,
            vext: 2,
            mismatch: 8,
            mapping_quality: None,
            move_indels_to_3: false,
            sv_min_len: 1000,
        }
    }
}

impl Configuration {
    pub(crate) const SEED_1: i32 = 17;
    pub(crate) const SEED_2: i32 = 12;
    pub(crate) const ADSEED: i32 = 6;
    pub(crate) const SVMAXLEN: i32 = 150000;
    pub(crate) const SVFLANK: i32 = 50;

    /// Any base with quality <=10 will be consider low quality in soft-clipped seq and extension will stop.
    pub(crate) const LOW_QUAL: i32 = 10;
}
