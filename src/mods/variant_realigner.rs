use std::collections::HashMap;

pub struct VariantRealigner {}

impl VariantRealigner {
    fn realign_del(&self, pos_to_del_count: &HashMap<i64, HashMap<String, usize>>) {
        
    }
}

/// Fill the temp tuple structure from hash tables of insertion or deletions positions to description string
/// and counts of variation.
fn fill_and_sort_tmp(
    changes: &HashMap<i64, HashMap<String, usize>>,
) -> Vec<SortPositionDescription<'_>> {
    let mut var_info = changes
        .iter()
        .flat_map(|(&pos, v)| {
            v.iter()
                .map(move |(desc, &cnt)| SortPositionDescription::new(pos, desc, cnt))
        })
        .collect::<Vec<_>>();

    var_info.sort_by(|o1, o2| {
        o2.count
            .cmp(&o1.count)
            .then(o1.pos.cmp(&o2.pos))
            .then(o2.desc_string.cmp(&o1.desc_string))
    });

    var_info
}

pub(crate) struct SortPositionDescription<'a> {
    pos: i64,
    desc_string: &'a str,
    count: usize,
}

impl<'a> SortPositionDescription<'a> {
    fn new(pos: i64, desc_string: &'a str, count: usize) -> Self {
        Self {
            pos,
            desc_string,
            count,
        }
    }
}
