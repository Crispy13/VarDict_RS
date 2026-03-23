use std::{collections::HashMap, sync::{Arc, OnceLock}};

use crate::conf::Configuration;
use crate::prelude::LibDefaultHasher;

pub static INSTANCE: OnceLock<GlobalReadOnlyScope> = OnceLock::new();

static INSTANCE_ARC: OnceLock<Arc<GlobalReadOnlyScope>> = OnceLock::new();

/// Get GlobalReadOnlyScope object.
///
/// Panic if it has not been initialized.
pub fn instance() -> &'static GlobalReadOnlyScope {
    INSTANCE.get().unwrap()
}

/// Get a shared Arc handle to the GlobalReadOnlyScope.
///
/// The Arc is created once and cached; subsequent calls return `Arc::clone`.
pub fn instance_arc() -> Arc<GlobalReadOnlyScope> {
    Arc::clone(INSTANCE_ARC.get_or_init(|| Arc::new(instance().clone())))
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
