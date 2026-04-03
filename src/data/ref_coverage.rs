use std::collections::HashMap;

use crate::prelude::LibDefaultHasher;

const LEFT_PAD: i64 = 500;
const RIGHT_PAD: i64 = 1500;
const SENTINEL: u32 = u32::MAX;

/// Reference coverage indexed by genomic position.
///
/// Most positions fall inside the current region window, so those counts are
/// stored in a dense vector. Rare out-of-window positions fall back to the
/// same HashMap hasher configuration used previously.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RefCoverage {
    dense: Vec<u32>,
    offset: i64,
    overflow: HashMap<i64, usize, LibDefaultHasher>,
}

impl Default for RefCoverage {
    fn default() -> Self {
        Self {
            dense: Vec::new(),
            offset: 0,
            overflow: HashMap::default(),
        }
    }
}

impl RefCoverage {
    /// Create coverage storage for a region plus the padding used by VarDict.
    pub fn new(region_start: usize, region_end: usize) -> Self {
        let region_start = i64::try_from(region_start).expect("region start exceeds i64");
        let region_end = i64::try_from(region_end).expect("region end exceeds i64");
        let offset = region_start - LEFT_PAD;
        let dense_end = region_end.max(region_start) + RIGHT_PAD;
        let dense_len = usize::try_from(dense_end - offset + 1)
            .expect("reference coverage window exceeds usize");

        Self {
            dense: vec![SENTINEL; dense_len],
            offset,
            overflow: HashMap::with_hasher(LibDefaultHasher::default()),
        }
    }

    /// Return the stored coverage at `pos`, or `None` when the position was never written.
    pub fn get(&self, pos: i64) -> Option<usize> {
        if let Some(index) = self.dense_index(pos) {
            let value = self.dense[index];
            return (value != SENTINEL).then_some(value as usize);
        }

        self.overflow.get(&pos).copied()
    }

    /// Set the coverage at `pos`, preserving explicitly stored zero values.
    pub fn set(&mut self, pos: i64, val: usize) {
        if let Some(index) = self.dense_index(pos) {
            self.dense[index] = u32::try_from(val).expect("reference coverage exceeds u32");
            return;
        }

        self.overflow.insert(pos, val);
    }

    /// Increment the coverage at `pos` by `delta`.
    pub fn inc(&mut self, pos: i64, delta: usize) {
        if delta == 0 {
            return;
        }

        if let Some(index) = self.dense_index(pos) {
            let delta = u32::try_from(delta).expect("reference coverage exceeds u32");
            let current = if self.dense[index] == SENTINEL {
                0
            } else {
                self.dense[index]
            };
            self.dense[index] = current
                .checked_add(delta)
                .expect("reference coverage exceeds u32");
            return;
        }

        *self.overflow.entry(pos).or_insert(0) += delta;
    }

    /// Return true when `pos` has any stored coverage value, including zero.
    pub fn contains_key(&self, pos: i64) -> bool {
        if let Some(index) = self.dense_index(pos) {
            return self.dense[index] != SENTINEL;
        }

        self.overflow.contains_key(&pos)
    }

    /// Return the smallest stored position across dense and overflow storage.
    pub fn min_key(&self) -> Option<i64> {
        let dense_min = self
            .dense
            .iter()
            .position(|&value| value != SENTINEL)
            .map(|index| self.offset + i64::try_from(index).expect("dense index exceeds i64"));
        let overflow_min = self.overflow.keys().min().copied();

        match (dense_min, overflow_min) {
            (Some(left), Some(right)) => Some(left.min(right)),
            (Some(left), None) => Some(left),
            (None, Some(right)) => Some(right),
            (None, None) => None,
        }
    }

    /// Return the largest stored position across dense and overflow storage.
    pub fn max_key(&self) -> Option<i64> {
        let dense_max = self
            .dense
            .iter()
            .rposition(|&value| value != SENTINEL)
            .map(|index| self.offset + i64::try_from(index).expect("dense index exceeds i64"));
        let overflow_max = self.overflow.keys().max().copied();

        match (dense_max, overflow_max) {
            (Some(left), Some(right)) => Some(left.max(right)),
            (Some(left), None) => Some(left),
            (None, Some(right)) => Some(right),
            (None, None) => None,
        }
    }

    /// Iterate over all stored positions in sorted genomic order.
    pub fn iter_sorted(&self) -> impl Iterator<Item = (i64, usize)> + '_ {
        let mut items = Vec::with_capacity(self.len());
        for (index, &value) in self.dense.iter().enumerate() {
            if value != SENTINEL {
                items.push((
                    self.offset + i64::try_from(index).expect("dense index exceeds i64"),
                    value as usize,
                ));
            }
        }
        items.extend(self.overflow.iter().map(|(&pos, &value)| (pos, value)));
        items.sort_unstable_by_key(|(pos, _)| *pos);
        items.into_iter()
    }

    /// Return the number of stored positions.
    pub fn len(&self) -> usize {
        self.dense.iter().filter(|&&value| value != SENTINEL).count() + self.overflow.len()
    }

    /// Return true when no positions are stored.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    fn dense_index(&self, pos: i64) -> Option<usize> {
        let delta = pos - self.offset;
        if delta < 0 {
            return None;
        }

        let index = usize::try_from(delta).ok()?;
        (index < self.dense.len()).then_some(index)
    }
}

impl IntoIterator for RefCoverage {
    type Item = (i64, usize);
    type IntoIter = std::vec::IntoIter<(i64, usize)>;

    fn into_iter(self) -> Self::IntoIter {
        let mut items = Vec::with_capacity(self.dense.len() + self.overflow.len());
        let offset = self.offset;

        items.extend(
            self.dense
                .into_iter()
                .enumerate()
                .filter_map(move |(index, value)| {
                    (value != SENTINEL).then_some((
                        offset + i64::try_from(index).expect("dense index exceeds i64"),
                        value as usize,
                    ))
                }),
        );
        items.extend(self.overflow);

        items.into_iter()
    }
}

#[cfg(test)]
mod tests {
    use super::RefCoverage;

    #[test]
    fn dense_get_set_and_inc_round_trip() {
        let mut coverage = RefCoverage::new(1_000, 1_100);

        coverage.set(1_025, 3);
        coverage.inc(1_025, 4);
        coverage.inc(1_026, 2);

        assert_eq!(coverage.get(1_025), Some(7));
        assert_eq!(coverage.get(1_026), Some(2));
        assert_eq!(coverage.get(1_027), None);
    }

    #[test]
    fn overflow_positions_fall_back_to_hashmap_storage() {
        let mut coverage = RefCoverage::new(1_000, 1_100);

        coverage.inc(100, 5);
        coverage.set(3_000, 9);

        assert_eq!(coverage.get(100), Some(5));
        assert_eq!(coverage.get(3_000), Some(9));
    }

    #[test]
    fn zero_values_are_tracked_as_present() {
        let mut coverage = RefCoverage::new(100, 200);

        coverage.set(150, 0);

        assert_eq!(coverage.get(150), Some(0));
        assert!(coverage.contains_key(150));
        assert_eq!(coverage.get(151), None);
        assert!(!coverage.contains_key(151));
    }

    #[test]
    fn overflow_zero_values_are_tracked_as_present() {
        let mut coverage = RefCoverage::new(1_000, 1_100);

        coverage.set(3_000, 0);

        assert_eq!(coverage.get(3_000), Some(0));
        assert!(coverage.contains_key(3_000));
    }

    #[test]
    fn min_key_checks_dense_and_overflow() {
        let mut coverage = RefCoverage::new(1_000, 1_100);

        coverage.set(1_020, 4);
        coverage.set(250, 3);

        assert_eq!(coverage.min_key(), Some(250));

        coverage.set(250, 0);
        assert_eq!(coverage.min_key(), Some(250));
    }

    #[test]
    fn max_key_checks_dense_and_overflow() {
        let mut coverage = RefCoverage::new(1_000, 1_100);

        coverage.set(1_020, 4);
        coverage.set(3_000, 3);

        assert_eq!(coverage.max_key(), Some(3_000));

        coverage.set(3_000, 0);
        assert_eq!(coverage.max_key(), Some(3_000));
    }

    #[test]
    fn iter_sorted_merges_dense_and_overflow_in_position_order() {
        let mut coverage = RefCoverage::new(1_000, 1_100);

        coverage.set(1_005, 2);
        coverage.set(1_003, 7);
        coverage.set(3_000, 5);
        coverage.set(250, 1);

        let observed: Vec<(i64, usize)> = coverage.iter_sorted().collect();
        assert_eq!(observed, vec![(250, 1), (1_003, 7), (1_005, 2), (3_000, 5)]);
    }
}
