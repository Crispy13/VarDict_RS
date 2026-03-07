use anyhow::{Error, anyhow};
use crackle_kit::tracing::{Level, event};
use std::{
    cmp::Ordering,
    fmt::Debug,
    ops::{Range, RangeFrom},
    slice::SliceIndex,
    str::Utf8Error,
    sync::atomic::Ordering as AtomicOrdering,
};

pub mod aligner;

pub fn print_exception_and_continue(
    exception: &Error,
    place: &str,
    place_def: &str,
    region: Option<&crate::data::region::Region>,
    conf: &crate::conf::Configuration,
) -> Result<(), Error> {
    match region {
        Some(region) => {
            event!(
                Level::ERROR,
                "There was Exception while processing {} on {} on region {}. The processing will be continued from the next {}.\n{:#}",
                place,
                place_def,
                region.to_region_string(),
                place,
                exception
            );
        }
        None => {
            event!(
                Level::ERROR,
                "There was Exception while processing {} on {} but region is undefined. The processing will be continued from the next {}.\n{:#}",
                place,
                place_def,
                place,
                exception
            );
        }
    }

    let current_count = conf.exception_counter.fetch_add(1, AtomicOrdering::SeqCst) + 1;

    if current_count > crate::conf::Configuration::MAX_EXCEPTION_COUNT {
        event!(
            Level::ERROR,
            "VarDict-rs fails (there were {} continued exceptions during the run).",
            current_count
        );
        return Err(anyhow!("Too many continued exceptions: {}", current_count));
    }

    Ok(())
}

pub fn round_half_even(pattern: &str, value: f64) -> f64 {
    let decimals = pattern.split('.').nth(1).map(|s| s.len()).unwrap_or(0);
    if decimals == 0 {
        return round_half_even_with_scale(value, 1.0).0;
    }

    let scale = 10_f64.powi(decimals as i32);
    round_half_even_with_scale(value, scale).0 / scale
}

fn round_half_even_with_scale(value: f64, scale: f64) -> (f64, bool) {
    if !value.is_finite() {
        return (value, false);
    }

    if !scale.is_finite() || scale <= 0.0 {
        return (value, false);
    }

    let scale_int = scale.round() as u64;
    if (scale - scale_int as f64).abs() < f64::EPSILON {
        if let Some((floor, frac_cmp_half)) =
            exact_scaled_floor_and_frac_cmp_half(value.abs(), scale_int)
        {
            let rounded = match frac_cmp_half {
                Ordering::Less => floor,
                Ordering::Greater => floor.saturating_add(1),
                Ordering::Equal => {
                    if floor % 2 == 0 {
                        floor
                    } else {
                        floor.saturating_add(1)
                    }
                }
            };

            let sign = if value < 0.0 { -1.0 } else { 1.0 };
            return (rounded as f64 * sign, true);
        }
    }

    let sign = if value < 0.0 { -1.0 } else { 1.0 };
    let abs_value = value.abs();
    let abs_scaled = abs_value * scale;
    let floor = abs_scaled.floor();
    let frac = abs_scaled - floor;
    let tie_tolerance = 1e-12;
    let rounded = if (frac - 0.5).abs() <= tie_tolerance {
        let tie_value = (floor + 0.5) / scale;
        if abs_value < tie_value {
            floor
        } else if abs_value > tie_value {
            floor + 1.0
        } else if (floor as i64) % 2 == 0 {
            floor
        } else {
            floor + 1.0
        }
    } else if frac > 0.5 {
        floor + 1.0
    } else {
        floor
    };

    (rounded * sign, true)
}

fn exact_scaled_floor_and_frac_cmp_half(abs_value: f64, scale: u64) -> Option<(u128, Ordering)> {
    let bits = abs_value.to_bits();
    let exp_bits = ((bits >> 52) & 0x7ff) as i32;
    let mantissa = bits & ((1u64 << 52) - 1);

    let (significand, exponent) = if exp_bits == 0 {
        if mantissa == 0 {
            return Some((0, Ordering::Less));
        }
        (mantissa, 1 - 1023 - 52)
    } else {
        ((1u64 << 52) | mantissa, exp_bits - 1023 - 52)
    };

    let numerator = (significand as u128).checked_mul(scale as u128)?;

    if exponent >= 0 {
        let shift = exponent as u32;
        if shift >= 128 {
            return None;
        }
        let floor = numerator.checked_shl(shift)?;
        return Some((floor, Ordering::Less));
    }

    let shift = (-exponent) as u32;
    if shift >= 128 {
        return Some((0, Ordering::Less));
    }

    let denominator = 1u128 << shift;
    let floor = numerator / denominator;
    let remainder = numerator % denominator;
    let cmp = remainder.saturating_mul(2).cmp(&denominator);
    Some((floor, cmp))
}

fn subbyte(s: &[u8], mut begin: i32, len: i32) -> Option<&[u8]> {
    if begin < 0 {
        begin = s.len() as i32 + begin;
    }

    if len > 0 {
        s.get(begin as usize..((begin + len) as usize).min(s.len()))
    } else if len == 0 {
        return Some(b"");
    } else {
        let end = s.len() as i32 + len;
        if end < begin {
            return Some(b"");
        }

        s.get(begin as usize..end as usize)
    }
}

// Define the trait
pub trait SliceExt<I: SliceIndex<Self> + Clone> {
    fn get_or_err(&self, index: I) -> Result<&I::Output, Error>;
}

// impl<T> SliceExt<Range<usize>> for [T] {
//     fn get_or_err(
//         &self,
//         index: Range<usize>,
//     ) -> Result<&<Range<usize> as SliceIndex<[T]>>::Output, Error> {
//         match self.get(index.clone()) {
//             Some(v) => Ok(v),
//             None => Err(anyhow!(
//                 "Indexing failed: idx:{:?} target_len:{}",
//                 index,
//                 self.len()
//             )),
//         }
//     }
// }

impl<T, I: SliceIndex<Self> + Clone + Debug> SliceExt<I> for [T] {
    fn get_or_err(&self, index: I) -> Result<&<I as SliceIndex<[T]>>::Output, Error> {
        match self.get(index.clone()) {
            Some(v) => Ok(v),
            None => Err(anyhow!(
                "Indexing failed: idx:{:?} target_len:{}",
                index,
                self.len()
            )),
        }
    }
}

pub trait SliceExt2<I> {
    type Output: ?Sized;

    fn get_with_int(&self, index: I) -> Result<&Self::Output, Error>;
}

macro_rules! impl_get_with_int_range {
    ($ty:ty) => {
        impl<T> SliceExt2<Range<$ty>> for [T] {
            type Output = [T];

            fn get_with_int(&self, index: Range<$ty>) -> Result<&Self::Output, Error> {
                let start = index.start;
                let end = index.end;

                let start_us = if start < 0 {
                    (<$ty>::try_from(self.len())? + start) as usize
                } else {
                    start as usize
                };

                let end_us = if end < 0 {
                    (<$ty>::try_from(self.len())? + end) as usize
                } else {
                    end as usize
                };

                self.get_or_err(start_us..end_us)
            }
        }
    };
}

macro_rules! impl_get_with_int_single {
    ($ty:ty) => {
        impl<T> SliceExt2<RangeFrom<$ty>> for [T] {
            type Output = [T];

            /// Index the slice from `index` to the end.
            fn get_with_int(&self, index: RangeFrom<$ty>) -> Result<&Self::Output, Error> {
                let start = index.start;

                let start_us = if start < 0 {
                    (<$ty>::try_from(self.len())? + start) as usize
                } else {
                    start as usize
                };

                self.get_or_err(start_us..)
            }
        }
    };
}

impl_get_with_int_range!(i32);
impl_get_with_int_range!(i64);
impl_get_with_int_single!(i32);
impl_get_with_int_single!(i64);

pub(crate) trait BytesExt {
    fn try_as_str(&self) -> Result<&str, Utf8Error>;
}

impl BytesExt for &[u8] {
    fn try_as_str(&self) -> Result<&str, Utf8Error> {
        str::from_utf8(self)
    }
}

#[cfg(test)]
mod tests {
    use std::sync::{Arc, atomic::AtomicUsize};

    use super::*;

    #[test]
    fn test_positive_range_i32() {
        let data = [10, 20, 30, 40, 50];
        // Standard slice: 1 to 4 -> [20, 30, 40]
        let res = data.get_with_int(1..4i32).unwrap();
        assert_eq!(res, &[20, 30, 40]);
    }

    #[test]
    fn test_negative_range_i32() {
        let data = [10, 20, 30, 40, 50];
        // Python: data[-3:-1] -> [30, 40]
        let res = data.get_with_int(-3..-1i32).unwrap();
        assert_eq!(res, &[30, 40]);
    }

    #[test]
    fn test_mixed_range_i32() {
        let data = [10, 20, 30, 40, 50];
        // Python: data[1:-1] -> [20, 30, 40]
        let res = data.get_with_int(1..-1i32).unwrap();
        assert_eq!(res, &[20, 30, 40]);
    }

    #[test]
    fn test_range_from_positive_i32() {
        let data = [10, 20, 30, 40, 50];
        // Python: data[2:] -> [30, 40, 50]
        let res = data.get_with_int(2i32..).unwrap();
        assert_eq!(res, &[30, 40, 50]);
    }

    #[test]
    fn test_range_from_negative_i32() {
        let data = [10, 20, 30, 40, 50];
        // Python: data[-2:] -> [40, 50]
        let res = data.get_with_int(-2i32..).unwrap();
        assert_eq!(res, &[40, 50]);
    }

    #[test]
    fn test_i64_indexing() {
        let data = [1, 2, 3, 4, 5];
        // Check that i64 works the same way
        let res = data.get_with_int(-2i64..).unwrap();
        assert_eq!(res, &[4, 5]);
    }

    #[test]
    fn test_out_of_bounds_errors() {
        let data = [1, 2, 3];

        // Positive OOB
        assert!(data.get_with_int(0..5i32).is_err());

        // Negative OOB (Start before 0)
        // len=3. -4 implies index -1.
        // -4 + 3 = -1, which is < 0. Your current logic (start_us) handles this via casting?
        // Let's trace your code: (-4 + 3) as usize -> -1 as usize -> HUGE number.
        // get_or_err(HUGE..) will return Err. So this is safe.
        assert!(data.get_with_int(-4i32..).is_err());
    }

    #[test]
    fn test_empty_slice() {
        let data: [i32; 0] = [];
        // 0..0 is valid for empty slice
        assert_eq!(data.get_with_int(0..0i32).unwrap(), &[]);

        // 0..1 is OOB
        assert!(data.get_with_int(0..1i32).is_err());

        // -1.. is OOB (len 0 + -1 = -1)
        assert!(data.get_with_int(-1i32..).is_err());
    }

    #[test]
    fn test_full_range_logic() {
        let data = [1, 2, 3];
        // 0..-0 (Python style) -> 0..3
        // In your logic: end < 0 check?
        // Wait, -0 is just 0. So end=0.
        // 0..0 returns empty slice.
        // Note: In Python list[0:-0] is empty. Correct.
        assert_eq!(data.get_with_int(0..0i32).unwrap(), &[]);
    }

    #[test]
    fn test_usize_as() {
        let data = [1, 2, 3];
        let l = 2 as usize;
        assert_eq!(data.get_with_int(-(l as i32)..).unwrap(), &[2, 3]);
    }

    #[test]
    fn test_print_exception_and_continue_threshold() {
        let mut conf = crate::conf::Configuration::default();
        conf.exception_counter = Arc::new(AtomicUsize::new(0));
        let region =
            crate::data::region::Region::new("chr1".to_string(), 1, 10, "GENE".to_string());

        for _ in 0..crate::conf::Configuration::MAX_EXCEPTION_COUNT {
            print_exception_and_continue(
                &anyhow!("test error"),
                "record",
                "read1",
                Some(&region),
                &conf,
            )
            .expect("should continue before threshold is exceeded");
        }

        let result = print_exception_and_continue(
            &anyhow!("test error"),
            "record",
            "read1",
            Some(&region),
            &conf,
        );

        assert!(result.is_err());
    }
}
