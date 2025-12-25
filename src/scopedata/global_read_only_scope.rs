use crate::conf::Configuration;

pub(crate) struct GlobalReadOnlyScope {
    pub(crate) amplicon_based_calling: bool,
    pub(crate) conf: Configuration,
}