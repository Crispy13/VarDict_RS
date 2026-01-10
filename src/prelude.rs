use std::{collections::HashMap, hash::RandomState};

use smallvec::SmallVec;

pub type LibDefaultHasher = RandomState;

pub type SmallVecBytes = SmallVec<[u8; 32]>;