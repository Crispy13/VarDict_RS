use std::{collections::HashMap, sync::OnceLock};

use crate::conf::Configuration;
use crate::prelude::LibDefaultHasher;

pub static INSTANCE: OnceLock<GlobalReadOnlyScope> = OnceLock::new();

/// Get GlobalReadOnlyScope object.
///
/// Panic if it has not been initialized.
pub fn instance() -> &'static GlobalReadOnlyScope {
    INSTANCE.get().unwrap()
}

#[derive(Default, Clone)]
pub struct GlobalReadOnlyScope {
    pub amplicon_based_calling: Option<String>,
    pub chr_lens: HashMap<String, usize, LibDefaultHasher>,
    pub bam_paths: Vec<String>,
    pub conf: Configuration,
    pub adaptor_forward: HashMap<String, usize, LibDefaultHasher>,
    pub adaptor_reverse: HashMap<String, usize, LibDefaultHasher>,
}
