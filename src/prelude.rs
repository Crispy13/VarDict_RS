use rustc_hash::FxBuildHasher;
use smallvec::SmallVec;

pub type LibDefaultHasher = FxBuildHasher;

pub type SmallVecBytes = SmallVec<[u8; 32]>;

/// Per-position inner map used for variant and count storage.
///
/// Currently backed by `VecMap` (Vec-based linear scan) which is faster and
/// uses ~2x less memory than HashMap for the typical 1-4 entries per position.
/// Benchmark data: `cargo bench --bench vecmap_vs_hashmap_bench`.
///
/// To switch to HashMap, replace these two aliases:
/// ```ignore
/// pub(crate) type InnerMap<K, V> = std::collections::HashMap<K, V, LibDefaultHasher>;
/// pub(crate) use std::collections::hash_map::Entry as InnerMapEntry;
/// ```
pub(crate) type InnerMap<K, V> = crate::utils::vec_map::VecMap<K, V>;
pub(crate) use crate::utils::vec_map::Entry as InnerMapEntry;
