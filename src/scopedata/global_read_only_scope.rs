use std::{collections::HashMap, sync::OnceLock};

use crate::conf::Configuration;

pub static INSTANCE: OnceLock<GlobalReadOnlyScope> = OnceLock::new();

/// Get GlobalReadOnlyScope object.
///
/// Panic if it has not been initialized.
pub fn instance() -> &'static GlobalReadOnlyScope {
    INSTANCE.get().unwrap()
}

#[derive(Default, Clone)]
pub struct GlobalReadOnlyScope {
    pub amplicon_based_calling: bool,
    pub chr_lens: HashMap<String, usize>,
    pub bam_paths: Vec<String>,
    pub conf: Configuration,
    pub adaptor_forward: HashMap<String, usize>,
    pub adaptor_reverse: HashMap<String, usize>,
}
