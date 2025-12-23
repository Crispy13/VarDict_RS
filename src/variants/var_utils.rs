use std::collections::HashMap;

use crate::variants::variants::Variant;

pub(crate) fn get_variants_from_map<'a>(
    var_map: &'a mut HashMap<i64, HashMap<String, Variant>>,
    start: i64,
    desc_string: &str,
) -> &'a Variant {
    let pos_map = var_map
        .entry(start)
        .or_insert_with(|| HashMap::with_capacity(1));

    if pos_map.contains_key(desc_string) {
        pos_map.get(desc_string).unwrap()
    } else {
        pos_map.entry(desc_string.to_string()).or_default()
    }

    // // Get a raw pointer to avoid the borrow checker blocking the 'None' branch
    // let var_ptr: *const Variant = match pos_map.get(desc_string) {
    //     // Fast Path: Found it. 1 Hash. 0 Allocations.
    //     Some(var) => var,

    //     // Slow Path: Not found.
    //     // We use 'entry' here to insert AND get the reference in one go.
    //     // 2 Hashes total (1 check above + 1 insert here).
    //     None => pos_map.entry(desc_string.to_string()).or_default(),
    // };

    // // SAFETY:
    // // 1. We know 'var_ptr' points to valid memory inside 'pos_map'.
    // // 2. We do not mutate 'pos_map' again after obtaining this pointer.
    // // 3. The returned reference lifetime 'a is tied to the map, preventing
    // //    the caller from invalidating the pointer while holding the reference.
    // unsafe { &*var_ptr }
}
