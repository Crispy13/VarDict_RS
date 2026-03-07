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

pub(crate) static BEGIN_DIGITS: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"^(\d+)").unwrap());

pub(crate) static BEGIN_MINUS_NUMBER_CARET: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"^-\d+\^").unwrap());

pub(crate) static HASH_GROUP_CARET_GROUP: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"#(.+)\^(.+)").unwrap());

pub(crate) static UP_NUMBER_END: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"\^(\d+)$").unwrap());

pub(crate) static CARET_ATGNC: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"\^([ATGNC]+)").unwrap());

pub(crate) static CARET_ATGC_END: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"\^([ATGC]+)$").unwrap());

pub(crate) static AMP_ATGC: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"&([ATGC]+)").unwrap());

pub(crate) static BEGIN_PLUS_ATGC: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"^\+([ATGC]+)").unwrap());

pub(crate) static HASH_ATGC: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"#([ATGC]+)").unwrap());

pub(crate) static ATGSs_AMP_ATGSs_END: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"(\+[ATGC]+)&[ATGC]+$").unwrap());

pub(crate) static DUP_NUM_ATGC: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"(?i)<dup(\d+)>([atgc]+)$").unwrap());

pub(crate) static DUP_NUM: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"(?i)<dup(\d+)").unwrap());

pub(crate) static INV_NUM: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"(?i)<inv(\d+)").unwrap());

pub(crate) static SOME_SV_NUMBERS: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"<(...)\d+>").unwrap());

pub(crate) static ANY_SV: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"<(...)>").unwrap());

pub(crate) static BEGIN_ATGC_END: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"^[ATGC]+$").unwrap());

pub(crate) static SA_CIGAR_D_S_3CLIP: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"\d\dS$").unwrap());

pub(crate) static SA_CIGAR_D_S_5CLIP: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"^\d\d+S").unwrap());
