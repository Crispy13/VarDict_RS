use std::{collections::HashMap, sync::{LazyLock, OnceLock}};

use crate::conf::Configuration;

pub(crate) static INSTANCE: OnceLock<GlobalReadOnlyScope> = OnceLock::new();

/// Get GlobalReadOnlyScope object.
///
/// Panic if it has not been initialized.
pub(crate) fn instance() -> &'static GlobalReadOnlyScope {
    INSTANCE.get().unwrap()
}

#[derive(Default)]
pub(crate) struct GlobalReadOnlyScope {
    pub(crate) amplicon_based_calling: bool,
    pub(crate) chr_lens: HashMap<String, usize>,
    pub(crate) conf: Configuration,

}
