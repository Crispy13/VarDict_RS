use std::collections::HashMap;

use crate::variants::variants::Variant;

pub mod reference;
pub(crate) mod patterns;

pub(crate) type VariantMap = HashMap<String, Variant>;