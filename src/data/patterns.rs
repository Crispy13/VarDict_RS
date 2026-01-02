use std::sync::LazyLock;

use regex::Regex;

macro_rules! get_cap_group {
    ($cap:expr, $idx:expr) => {
        $cap.get($idx).ok_or_else(|| {
            anyhow!(
                "Can't get group {} from {}",
                $idx,
                $cap.get_match().as_str()
            )
        })
    };
}
pub(crate) use get_cap_group;

pub(crate) static BEGIN_MINUS_NUMBER: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"^-(\d+)").unwrap());

pub(crate) static UP_NUMBER_END: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"\^(\d+)$").unwrap());

pub(crate) static SA_CIGAR_D_S_3CLIP: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"\d\dS$").unwrap());

pub(crate) static SA_CIGAR_D_S_5CLIP: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"^\d\d+S$").unwrap());
