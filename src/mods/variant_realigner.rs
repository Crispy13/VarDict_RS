use std::collections::HashMap;

use anyhow::{Error, anyhow};
use crackle_kit::tracing::{Level, event};

use crate::{
    data::patterns::{BEGIN_MINUS_NUMBER, UP_NUMBER_END, get_cap_group},
    variants::{var_utils::get_variants_from_map, variants::Variant},
};

pub struct VariantRealigner {
    non_insertion_vars: HashMap<i64, HashMap<String, Variant>>,
    ref_coverage: HashMap<i64, usize>,
}

impl VariantRealigner {
    fn realign_del(
        &mut self,
        bam_parameters: &[&str],
        pos_to_del_count: &HashMap<i64, HashMap<String, usize>>,
    ) -> Result<(), Error> {
        let bams: &[&str] = todo!();

        let sorted_pos_to_del = fill_and_sort_tmp(pos_to_del_count);

        let mut last_pos = 0;
        for tpl in sorted_pos_to_del {
            let p = tpl.pos;
            last_pos = p;

            let vn = tpl.desc_string;
            let del_cnt = tpl.count;

            event!(
                Level::INFO,
                "  Realigndel for: {p} {vn} {del_cnt} cov: {}",
                self.ref_coverage.get(&p).copied().unwrap_or(0)
            );

            let var = get_variants_from_map(&mut self.non_insertion_vars, p, vn);

            let mut del_len = 0;

            match BEGIN_MINUS_NUMBER.captures(vn) {
                Some(cap) => {
                    del_len = get_cap_group!(cap, 1)?.as_str().parse::<i64>()?;
                }
                None => {}
            }

            match UP_NUMBER_END.captures(vn) {
                Some(cap) => {
                    del_len += get_cap_group!(cap, 1)?.as_str().parse::<i64>()?;
                }
                None => {}
            }
        }

        todo!()
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
