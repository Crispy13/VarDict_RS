use std::collections::HashMap;
use std::hash::Hash;

use rustc_hash::FxBuildHasher;
use smallvec::SmallVec;

pub type LibDefaultHasher = FxBuildHasher;

pub type SmallVecBytes = SmallVec<[u8; 32]>;

// ---------------------------------------------------------------------------
// InnerMap — per-position map trait + backend selection
// ---------------------------------------------------------------------------
//
// Both VecMap and HashMap implement `InnerMapTrait`, which guarantees they
// expose the same core API. Production code uses the concrete `InnerMap` alias.
//
// **Current backend: VecMap** — Vec-based linear scan, faster and ~2x less
// memory than HashMap for the typical 1–4 entries per position.
// Benchmark: `cargo bench --bench vecmap_vs_hashmap_bench`
//
// To switch to HashMap, change these two lines:
//   pub(crate) type InnerMap<K, V> = HashMap<K, V, LibDefaultHasher>;
//   pub(crate) use std::collections::hash_map::Entry as InnerMapEntry;

pub(crate) type InnerMap<K, V> = crate::utils::vec_map::VecMap<K, V>;
pub(crate) use crate::utils::vec_map::Entry as InnerMapEntry;

/// Shared contract for per-position inner maps.
///
/// Both `VecMap` and `HashMap<K, V, FxBuildHasher>` implement this trait,
/// enforcing API compatibility at compile time so the backend can be swapped
/// by changing the `InnerMap` type alias above.
#[allow(dead_code)]
pub(crate) trait InnerMapTrait<K, V> {
    fn get(&self, key: &K) -> Option<&V>;
    fn get_mut(&mut self, key: &K) -> Option<&mut V>;
    fn insert(&mut self, key: K, value: V) -> Option<V>;
    fn remove(&mut self, key: &K) -> Option<V>;
    fn contains_key(&self, key: &K) -> bool;
    fn len(&self) -> usize;
    fn is_empty(&self) -> bool;
    fn clear(&mut self);
    fn retain<F: FnMut(&K, &mut V) -> bool>(&mut self, f: F);
    fn shrink_to_fit(&mut self);
}

impl<K: PartialEq, V> InnerMapTrait<K, V> for crate::utils::vec_map::VecMap<K, V> {
    fn get(&self, key: &K) -> Option<&V> {
        self.get(key)
    }
    fn get_mut(&mut self, key: &K) -> Option<&mut V> {
        self.get_mut(key)
    }
    fn insert(&mut self, key: K, value: V) -> Option<V> {
        self.insert(key, value)
    }
    fn remove(&mut self, key: &K) -> Option<V> {
        self.remove(key)
    }
    fn contains_key(&self, key: &K) -> bool {
        self.contains_key(key)
    }
    fn len(&self) -> usize {
        self.len()
    }
    fn is_empty(&self) -> bool {
        self.is_empty()
    }
    fn clear(&mut self) {
        self.clear()
    }
    fn retain<F: FnMut(&K, &mut V) -> bool>(&mut self, f: F) {
        self.retain(f)
    }
    fn shrink_to_fit(&mut self) {
        self.shrink_to_fit()
    }
}

impl<K: Eq + Hash, V> InnerMapTrait<K, V> for HashMap<K, V, LibDefaultHasher> {
    fn get(&self, key: &K) -> Option<&V> {
        self.get(key)
    }
    fn get_mut(&mut self, key: &K) -> Option<&mut V> {
        self.get_mut(key)
    }
    fn insert(&mut self, key: K, value: V) -> Option<V> {
        self.insert(key, value)
    }
    fn remove(&mut self, key: &K) -> Option<V> {
        self.remove(key)
    }
    fn contains_key(&self, key: &K) -> bool {
        self.contains_key(key)
    }
    fn len(&self) -> usize {
        self.len()
    }
    fn is_empty(&self) -> bool {
        self.is_empty()
    }
    fn clear(&mut self) {
        self.clear()
    }
    fn retain<F: FnMut(&K, &mut V) -> bool>(&mut self, f: F) {
        self.retain(f)
    }
    fn shrink_to_fit(&mut self) {
        self.shrink_to_fit()
    }
}
