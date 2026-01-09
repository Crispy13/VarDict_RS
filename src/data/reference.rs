use std::collections::HashMap;


#[derive(Default)]
pub(crate) struct Reference {
    pub(crate) ref_seq: Vec<u8>,
    pub(crate) seed: HashMap<Vec<u8>,Vec<i64>>,
}