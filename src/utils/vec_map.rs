use std::{fmt, iter::FusedIterator, mem, slice, vec};

/// Compact map backed by a `Vec<(K, V)>` for tiny key sets.
///
/// # DIVERGENCE FROM JAVA
///
/// **This is an optimized data structure, NOT a faithful port of Java logic.**
/// Java VarDictJava uses `HashMap<String, Variation>` (via `VariationMap`) for inner
/// variant maps. This Rust port replaces those with `VecMap` for memory efficiency.
///
/// - Java `HashMap` lazily allocates buckets; Rust `hashbrown::HashMap` pre-allocates
///   ~128–388 bytes per instance even for 1–4 entries. With ~5M inner maps, this caused
///   Rust to use ~2x Java's memory.
/// - VecMap assumes inner maps typically hold 1–4 entries (one per allele at a genomic
///   position). This assumption is based on diploid biology, NOT guaranteed for all inputs.
/// - Edge cases (homopolymer runs, high-coverage amplicon >1000x, noisy/repetitive regions)
///   can produce 5–10+ entries per position.
///
/// # TODO
///
/// - **Needs more tests**: Current coverage is 8 unit tests. Should add stress tests with
///   >10 entries and benchmark VecMap vs HashMap at various sizes.
/// - **Revert to HashMap if benchmarks indicate regression**: If profiling on production
///   workloads shows VecMap scan is slower than HashMap (crossover at ~15–20 entries),
///   revert `RawVarMap` and `CountMap` type aliases in `vardict_pipeline.rs` back to
///   `HashMap<..., LibDefaultHasher>`.
/// - **Consider hybrid approach**: VecMap for ≤N entries, promote to HashMap above N.
#[derive(Clone, Default, PartialEq)]
pub struct VecMap<K, V> {
    entries: Vec<(K, V)>,
}

impl<K, V> VecMap<K, V> {
    /// Creates an empty map.
    pub const fn new() -> Self {
        Self {
            entries: Vec::new(),
        }
    }

    /// Creates an empty map with capacity for at least `capacity` entries.
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            entries: Vec::with_capacity(capacity),
        }
    }

    /// Returns the number of stored entries.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Returns `true` when the map contains no entries.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Removes every entry from the map.
    pub fn clear(&mut self) {
        self.entries.clear();
    }

    /// Retains only the entries for which the predicate returns `true`.
    pub fn retain<F>(&mut self, mut f: F)
    where
        F: FnMut(&K, &mut V) -> bool,
    {
        self.entries.retain_mut(|(key, value)| f(key, value));
    }

    /// Shrinks the backing allocation to fit the current length.
    pub fn shrink_to_fit(&mut self) {
        self.entries.shrink_to_fit();
    }

    /// Returns an iterator over key-value pairs.
    pub fn iter(&self) -> Iter<'_, K, V> {
        Iter {
            inner: self.entries.iter(),
        }
    }

    /// Returns a mutable iterator over key-value pairs.
    pub fn iter_mut(&mut self) -> IterMut<'_, K, V> {
        IterMut {
            inner: self.entries.iter_mut(),
        }
    }

    /// Returns an iterator over keys.
    pub fn keys(&self) -> Keys<'_, K, V> {
        Keys {
            inner: self.entries.iter(),
        }
    }

    /// Returns an iterator over values.
    pub fn values(&self) -> Values<'_, K, V> {
        Values {
            inner: self.entries.iter(),
        }
    }

    /// Returns a mutable iterator over values.
    pub fn values_mut(&mut self) -> ValuesMut<'_, K, V> {
        ValuesMut {
            inner: self.entries.iter_mut(),
        }
    }

    /// Consumes the map and returns its backing iterator.
    pub fn into_iter(self) -> vec::IntoIter<(K, V)> {
        self.entries.into_iter()
    }
}

impl<K: PartialEq, V> VecMap<K, V> {
    fn position(&self, key: &K) -> Option<usize> {
        self.entries
            .iter()
            .position(|(existing_key, _)| existing_key == key)
    }

    /// Returns the value for `key` when present.
    pub fn get(&self, key: &K) -> Option<&V> {
        self.entries
            .iter()
            .find(|(existing_key, _)| existing_key == key)
            .map(|(_, value)| value)
    }

    /// Returns a mutable reference to the value for `key` when present.
    pub fn get_mut(&mut self, key: &K) -> Option<&mut V> {
        self.entries
            .iter_mut()
            .find(|(existing_key, _)| existing_key == key)
            .map(|(_, value)| value)
    }

    /// Inserts a key-value pair, replacing and returning the old value when the key exists.
    pub fn insert(&mut self, key: K, value: V) -> Option<V> {
        if let Some((_, stored_value)) = self
            .entries
            .iter_mut()
            .find(|(existing_key, _)| existing_key == &key)
        {
            return Some(mem::replace(stored_value, value));
        }

        self.entries.push((key, value));
        None
    }

    /// Removes the entry for `key`, returning its value when present.
    pub fn remove(&mut self, key: &K) -> Option<V> {
        self.position(key)
            .map(|index| self.entries.swap_remove(index).1)
    }

    /// Returns `true` when `key` exists in the map.
    pub fn contains_key(&self, key: &K) -> bool {
        self.position(key).is_some()
    }

    /// Returns an entry view for in-place manipulation.
    pub fn entry(&mut self, key: K) -> Entry<'_, K, V> {
        match self.position(&key) {
            Some(index) => Entry::Occupied(OccupiedEntry { map: self, index }),
            None => Entry::Vacant(VacantEntry { map: self, key }),
        }
    }
}

impl<K: fmt::Debug, V: fmt::Debug> fmt::Debug for VecMap<K, V> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_map().entries(self.iter()).finish()
    }
}

impl<K: PartialEq, V> FromIterator<(K, V)> for VecMap<K, V> {
    fn from_iter<T: IntoIterator<Item = (K, V)>>(iter: T) -> Self {
        let mut iterator = iter.into_iter();
        let (lower_bound, _) = iterator.size_hint();
        let mut map = Self::with_capacity(lower_bound);

        for (key, value) in iterator {
            map.insert(key, value);
        }

        map
    }
}

impl<K, V> IntoIterator for VecMap<K, V> {
    type Item = (K, V);
    type IntoIter = vec::IntoIter<(K, V)>;

    fn into_iter(self) -> Self::IntoIter {
        self.entries.into_iter()
    }
}

impl<'a, K, V> IntoIterator for &'a VecMap<K, V> {
    type Item = (&'a K, &'a V);
    type IntoIter = Iter<'a, K, V>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

impl<'a, K, V> IntoIterator for &'a mut VecMap<K, V> {
    type Item = (&'a K, &'a mut V);
    type IntoIter = IterMut<'a, K, V>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter_mut()
    }
}

/// Entry view into a `VecMap`.
pub enum Entry<'a, K, V> {
    Occupied(OccupiedEntry<'a, K, V>),
    Vacant(VacantEntry<'a, K, V>),
}

impl<'a, K, V> Entry<'a, K, V> {
    /// Ensures a value is present by inserting `default` when vacant.
    pub fn or_insert(self, default: V) -> &'a mut V {
        match self {
            Self::Occupied(entry) => entry.into_mut(),
            Self::Vacant(entry) => entry.insert(default),
        }
    }

    /// Ensures a value is present by lazily inserting the produced value when vacant.
    pub fn or_insert_with<F>(self, default: F) -> &'a mut V
    where
        F: FnOnce() -> V,
    {
        match self {
            Self::Occupied(entry) => entry.into_mut(),
            Self::Vacant(entry) => entry.insert(default()),
        }
    }

    /// Applies `f` to an occupied value and returns the entry view unchanged.
    pub fn and_modify<F>(self, f: F) -> Self
    where
        F: FnOnce(&mut V),
    {
        match self {
            Self::Occupied(mut entry) => {
                f(entry.get_mut());
                Self::Occupied(entry)
            }
            Self::Vacant(entry) => Self::Vacant(entry),
        }
    }
}

impl<'a, K, V: Default> Entry<'a, K, V> {
    /// Ensures a default value is present and returns a mutable reference to it.
    pub fn or_default(self) -> &'a mut V {
        self.or_insert_with(V::default)
    }
}

/// Occupied entry view into a `VecMap`.
pub struct OccupiedEntry<'a, K, V> {
    map: &'a mut VecMap<K, V>,
    index: usize,
}

impl<'a, K, V> OccupiedEntry<'a, K, V> {
    /// Returns a shared reference to the stored value.
    pub fn get(&self) -> &V {
        &self.map.entries[self.index].1
    }

    /// Returns a mutable reference to the stored value.
    pub fn get_mut(&mut self) -> &mut V {
        &mut self.map.entries[self.index].1
    }

    /// Converts the entry into a mutable reference to the stored value.
    pub fn into_mut(self) -> &'a mut V {
        &mut self.map.entries[self.index].1
    }

    /// Replaces the stored value and returns the previous one.
    pub fn insert(&mut self, value: V) -> V {
        mem::replace(&mut self.map.entries[self.index].1, value)
    }

    /// Returns a shared reference to the stored key.
    pub fn key(&self) -> &K {
        &self.map.entries[self.index].0
    }
}

/// Vacant entry view into a `VecMap`.
pub struct VacantEntry<'a, K, V> {
    map: &'a mut VecMap<K, V>,
    key: K,
}

impl<'a, K, V> VacantEntry<'a, K, V> {
    /// Inserts `value` and returns a mutable reference to it.
    pub fn insert(self, value: V) -> &'a mut V {
        self.map.entries.push((self.key, value));
        let last_index = self.map.entries.len() - 1;
        &mut self.map.entries[last_index].1
    }

    /// Returns a shared reference to the vacant key.
    pub fn key(&self) -> &K {
        &self.key
    }
}

/// Iterator over `VecMap` keys and values.
pub struct Iter<'a, K, V> {
    inner: slice::Iter<'a, (K, V)>,
}

impl<'a, K, V> Iterator for Iter<'a, K, V> {
    type Item = (&'a K, &'a V);

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next().map(|(key, value)| (key, value))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.inner.size_hint()
    }
}

impl<'a, K, V> DoubleEndedIterator for Iter<'a, K, V> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.inner.next_back().map(|(key, value)| (key, value))
    }
}

impl<'a, K, V> ExactSizeIterator for Iter<'a, K, V> {}
impl<'a, K, V> FusedIterator for Iter<'a, K, V> {}

/// Mutable iterator over `VecMap` keys and values.
pub struct IterMut<'a, K, V> {
    inner: slice::IterMut<'a, (K, V)>,
}

impl<'a, K, V> Iterator for IterMut<'a, K, V> {
    type Item = (&'a K, &'a mut V);

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next().map(|(key, value)| (&*key, value))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.inner.size_hint()
    }
}

impl<'a, K, V> DoubleEndedIterator for IterMut<'a, K, V> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.inner.next_back().map(|(key, value)| (&*key, value))
    }
}

impl<'a, K, V> ExactSizeIterator for IterMut<'a, K, V> {}
impl<'a, K, V> FusedIterator for IterMut<'a, K, V> {}

/// Iterator over `VecMap` keys.
pub struct Keys<'a, K, V> {
    inner: slice::Iter<'a, (K, V)>,
}

impl<'a, K, V> Iterator for Keys<'a, K, V> {
    type Item = &'a K;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next().map(|(key, _)| key)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.inner.size_hint()
    }
}

impl<'a, K, V> DoubleEndedIterator for Keys<'a, K, V> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.inner.next_back().map(|(key, _)| key)
    }
}

impl<'a, K, V> ExactSizeIterator for Keys<'a, K, V> {}
impl<'a, K, V> FusedIterator for Keys<'a, K, V> {}

/// Iterator over `VecMap` values.
pub struct Values<'a, K, V> {
    inner: slice::Iter<'a, (K, V)>,
}

impl<'a, K, V> Iterator for Values<'a, K, V> {
    type Item = &'a V;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next().map(|(_, value)| value)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.inner.size_hint()
    }
}

impl<'a, K, V> DoubleEndedIterator for Values<'a, K, V> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.inner.next_back().map(|(_, value)| value)
    }
}

impl<'a, K, V> ExactSizeIterator for Values<'a, K, V> {}
impl<'a, K, V> FusedIterator for Values<'a, K, V> {}

/// Mutable iterator over `VecMap` values.
pub struct ValuesMut<'a, K, V> {
    inner: slice::IterMut<'a, (K, V)>,
}

impl<'a, K, V> Iterator for ValuesMut<'a, K, V> {
    type Item = &'a mut V;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next().map(|(_, value)| value)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.inner.size_hint()
    }
}

impl<'a, K, V> DoubleEndedIterator for ValuesMut<'a, K, V> {
    fn next_back(&mut self) -> Option<Self::Item> {
        self.inner.next_back().map(|(_, value)| value)
    }
}

impl<'a, K, V> ExactSizeIterator for ValuesMut<'a, K, V> {}
impl<'a, K, V> FusedIterator for ValuesMut<'a, K, V> {}

#[cfg(test)]
mod tests {
    use super::{Entry, VecMap};

    #[test]
    fn basic_insert_get_get_mut_and_remove() {
        let mut map = VecMap::with_capacity(2);
        assert!(map.is_empty());
        assert_eq!(map.insert("a", 1), None);
        assert_eq!(map.insert("b", 2), None);
        assert_eq!(map.len(), 2);
        assert!(map.contains_key(&"a"));
        assert_eq!(map.get(&"a"), Some(&1));
        assert_eq!(map.get(&"missing"), None);

        if let Some(value) = map.get_mut(&"b") {
            *value += 5;
        }

        assert_eq!(map.get(&"b"), Some(&7));
        assert_eq!(map.remove(&"a"), Some(1));
        assert_eq!(map.remove(&"a"), None);
        assert_eq!(map.len(), 1);
    }

    #[test]
    fn entry_api_supports_default_insert_modify_and_replace() {
        let mut map = VecMap::new();

        *map.entry("alpha").or_default() += 1;
        assert_eq!(map.get(&"alpha"), Some(&1));

        *map.entry("beta").or_insert(5) += 2;
        *map.entry("gamma").or_insert_with(|| 9) += 1;

        let entry = map.entry("alpha").and_modify(|value| *value += 4);
        let mut occupied = match entry {
            Entry::Occupied(entry) => entry,
            Entry::Vacant(_) => panic!("entry unexpectedly vacant"),
        };

        assert_eq!(occupied.key(), &"alpha");
        assert_eq!(occupied.get(), &5);
        assert_eq!(occupied.insert(10), 5);
        *occupied.get_mut() += 1;
        assert_eq!(map.get(&"alpha"), Some(&11));

        match map.entry("beta") {
            Entry::Occupied(entry) => {
                *entry.into_mut() += 3;
            }
            Entry::Vacant(_) => panic!("entry unexpectedly vacant"),
        }

        assert_eq!(map.get(&"beta"), Some(&10));
        assert_eq!(map.get(&"gamma"), Some(&10));

        let vacant_key = match map.entry("delta") {
            Entry::Occupied(_) => panic!("entry unexpectedly occupied"),
            Entry::Vacant(entry) => *entry.key(),
        };
        assert_eq!(vacant_key, "delta");
    }

    #[test]
    fn iterators_cover_iter_iter_mut_keys_values_and_into_iter() {
        let mut map = VecMap::from_iter([("a", 1), ("b", 2), ("c", 3)]);

        let pairs: Vec<_> = map.iter().map(|(key, value)| (*key, *value)).collect();
        assert_eq!(pairs, vec![("a", 1), ("b", 2), ("c", 3)]);

        for (_, value) in map.iter_mut() {
            *value += 10;
        }

        for value in map.values_mut() {
            *value += 1;
        }

        assert_eq!(map.keys().copied().collect::<Vec<_>>(), vec!["a", "b", "c"]);
        assert_eq!(map.values().copied().collect::<Vec<_>>(), vec![12, 13, 14]);

        let borrowed_pairs: Vec<_> = (&map)
            .into_iter()
            .map(|(key, value)| (*key, *value))
            .collect();
        assert_eq!(borrowed_pairs, vec![("a", 12), ("b", 13), ("c", 14)]);

        let consumed_pairs: Vec<_> = map.clone().into_iter().collect();
        assert_eq!(consumed_pairs, vec![("a", 12), ("b", 13), ("c", 14)]);

        let mutable_pairs: Vec<_> = (&mut map)
            .into_iter()
            .map(|(key, value)| {
                *value += 1;
                (*key, *value)
            })
            .collect();
        assert_eq!(mutable_pairs, vec![("a", 13), ("b", 14), ("c", 15)]);
    }

    #[test]
    fn keys_support_next_back_min_max_and_copied_collect() {
        let map = VecMap::from_iter([(3, "c"), (1, "a"), (2, "b")]);
        let mut keys = map.keys();

        assert_eq!(keys.next(), Some(&3));
        assert_eq!(keys.next_back(), Some(&2));
        assert_eq!(map.keys().min(), Some(&1));
        assert_eq!(map.keys().max(), Some(&3));
        assert_eq!(map.keys().copied().collect::<Vec<_>>(), vec![3, 1, 2]);
    }

    #[test]
    fn clone_and_default_behave_as_expected() {
        let empty = VecMap::<i32, i32>::default();
        assert!(empty.is_empty());

        let original = VecMap::from_iter([(1, 10), (2, 20)]);
        let cloned = original.clone();
        assert_eq!(cloned, original);
    }

    #[test]
    fn from_iterator_keeps_last_value_for_duplicate_keys() {
        let map = VecMap::from_iter([("a", 1), ("b", 2), ("a", 3)]);

        assert_eq!(map.len(), 2);
        assert_eq!(map.get(&"a"), Some(&3));
        assert_eq!(map.get(&"b"), Some(&2));
    }

    #[test]
    fn clear_retain_and_shrink_to_fit_work() {
        let mut map = VecMap::from_iter([(1, 10), (2, 20), (3, 30), (4, 40)]);

        map.retain(|key, value| {
            *value += 1;
            key % 2 == 0
        });
        map.shrink_to_fit();

        assert_eq!(map.len(), 2);
        assert_eq!(map.get(&2), Some(&21));
        assert_eq!(map.get(&4), Some(&41));
        assert_eq!(map.get(&1), None);

        map.clear();
        assert!(map.is_empty());
    }

    #[test]
    fn duplicate_insert_replaces_without_growing_and_missing_remove_is_none() {
        let mut map = VecMap::new();

        assert_eq!(map.insert("dup", 1), None);
        assert_eq!(map.insert("dup", 2), Some(1));
        assert_eq!(map.len(), 1);
        assert_eq!(map.remove(&"missing"), None);
        assert_eq!(map.get(&"dup"), Some(&2));
    }
}
