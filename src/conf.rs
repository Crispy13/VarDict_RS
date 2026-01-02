pub(crate) struct Configuration {
    pub(crate) perform_local_realignment: bool,
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
}

impl Default for Configuration {
    fn default() -> Self {
        Self {
            perform_local_realignment: Default::default(),
            chimeric_filter: Default::default(),
            min_match: 0,
            sam_filter: 0x504,
            crispr_cutting_site: 0,
            crispr_filtering_bp: 0,
            disable_sv: false,
            unique_mode_alignment_enabled: false,
            unique_mode_second_in_pair_enabled: false,
        }
    }
}

impl Configuration {
    pub(crate) const SEED_2: i32 = 12;

    /// Any base with quality <=10 will be consider low quality in soft-clipped seq and extension will stop.
    pub(crate) const LOW_QUAL: i32 = 10;
}
