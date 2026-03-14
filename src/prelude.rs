use rustc_hash::FxBuildHasher;
use smallvec::SmallVec;

pub type LibDefaultHasher = FxBuildHasher;

pub type SmallVecBytes = SmallVec<[u8; 32]>;
