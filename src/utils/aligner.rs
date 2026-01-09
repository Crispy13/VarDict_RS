pub(crate) enum Aligner {
    BWA,
    Star,
}

impl Default for Aligner {
    fn default() -> Self {
        Self::BWA
    }
}

impl Aligner {
    pub(crate) fn nm_tag(&self) -> &'static [u8; 2]  {
        match self {
            Aligner::BWA => b"NM",
            Aligner::Star => b"nM",
        }
    }
}
