use std::{collections::HashMap, env, sync::OnceLock};

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
    pub conf: Configuration,
    pub adaptor_forward: HashMap<String, usize>,
    pub adaptor_reverse: HashMap<String, usize>,
    pub debug_dump_steps: bool,
    pub debug_dump_region: Option<DebugDumpRegion>,

}

#[derive(Clone, Debug, Default)]
pub struct DebugDumpRegion {
    pub chr: String,
    pub start: i64,
    pub end: i64,
}

pub fn parse_debug_dump_region_env() -> Option<DebugDumpRegion> {
    let value = env::var("VARDICT_DEBUG_DUMP_REGION").ok()?;
    let value = value.trim();
    if value.is_empty() {
        return None;
    }

    let (chr, range) = value.split_once(':')?;
    let (start_str, end_str) = range.split_once('-')?;
    let start = start_str.replace(',', "").parse::<i64>().ok()?;
    let end = end_str.replace(',', "").parse::<i64>().ok()?;
    let (start, end) = if start <= end { (start, end) } else { (end, start) };

    Some(DebugDumpRegion {
        chr: chr.trim().to_string(),
        start,
        end,
    })
}

impl GlobalReadOnlyScope {
    pub fn should_dump_steps_for(&self, chr: &str, start: i64, end: i64) -> bool {
        if !self.debug_dump_steps {
            return false;
        }

        if let Some(region) = &self.debug_dump_region {
            if chr != region.chr {
                return false;
            }

            let (query_start, query_end) = if start <= end { (start, end) } else { (end, start) };
            let overlaps = query_start <= region.end && query_end >= region.start;
            return overlaps;
        }

        true
    }
}
