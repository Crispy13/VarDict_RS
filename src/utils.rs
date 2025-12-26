use anyhow::{Error, anyhow};
use std::{ops::Range, slice::SliceIndex};

pub mod aligner;

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

impl<T> SliceExt<Range<usize>> for [T] {
    fn get_or_err(
        &self,
        index: Range<usize>,
    ) -> Result<&<Range<usize> as SliceIndex<[T]>>::Output, Error> {
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
