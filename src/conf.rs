pub(crate) struct Configuration {
    pub(crate) perform_local_realignment: bool,

    /// Indicate to turn off chimeric reads filtering.  Chimeric reads are artifacts from library construction,
    /// where a read can be split into two segments, each will be aligned within 1-2 read length distance,
    /// but in opposite direction.
    pub(crate) chimeric_filter: bool,
    pub(crate) min_match: i32,

    pub(crate) disable_sv: bool,

    pub(crate) unique_mode_alignment_enabled: bool,
    pub(crate) unique_mode_second_in_pair_enabled: bool,

    /// The hexical to filter reads.
    pub(crate) sam_filter: u32,

    pub(crate) crispr_cutting_site: i32,
    pub(crate) crispr_filtering_bp: i32,
    // pub(crate) seed_2: i32,

    /// The phred score for a base to be considered a good call
    pub(crate) goodq: f64,

    /// Extension of bp to look for mismatches after insertion or deletion
    pub(crate) vext: i32,
}

impl Default for Configuration {
    fn default() -> Self {
        Self {
            perform_local_realignment: Default::default(),
            chimeric_filter: false,
            min_match: 0,
            sam_filter: 0x504,
            crispr_cutting_site: 0,
            crispr_filtering_bp: 0,
            disable_sv: false,
            unique_mode_alignment_enabled: false,
            unique_mode_second_in_pair_enabled: false,
            goodq: 22.5,
        }
    }
}

impl Configuration {
    pub(crate) const SEED_1: i32 = 17;
    pub(crate) const SEED_2: i32 = 12;

    /// Any base with quality <=10 will be consider low quality in soft-clipped seq and extension will stop.
    pub(crate) const LOW_QUAL: i32 = 10;
}
