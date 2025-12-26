pub(crate) struct Configuration {
    pub(crate) perform_local_realignment: bool,
    pub(crate) chimeric_filter: bool,

    // pub(crate) seed_2: i32,
}

impl Configuration {
    pub(crate) const SEED_2:i32 = 12;
}

