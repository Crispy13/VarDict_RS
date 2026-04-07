#[test]
#[ignore = "Requires somatic test data (tumor/normal BAM pair) - see tests/parity/somatic_amplicon_README.md"]
fn test_somatic_simple_tumor_normal_parity() {
    // TODO: Implement once somatic test BAMs are available.
    // Expected flow:
    // 1. Run Java VarDict with -b tumor.bam -b2 normal.bam.
    // 2. Run Rust VarDict with the same flags.
    // 3. Compare output byte-for-byte.
    todo!("Somatic parity test not yet implemented");
}

#[test]
#[ignore = "Requires somatic fixture coverage for low-VAF and filtering cases - see tests/parity/somatic_amplicon_README.md"]
fn test_somatic_low_vaf_filter_parity() {
    // TODO: Implement once low-VAF somatic fixture data is available.
    // Expected flow:
    // 1. Use a tumor/normal region with a known low-frequency somatic event.
    // 2. Run Java and Rust with identical filtering flags.
    // 3. Compare output byte-for-byte.
    todo!("Somatic low-VAF parity test not yet implemented");
}

#[test]
#[ignore = "Requires somatic structural-variant fixtures - see tests/parity/somatic_amplicon_README.md"]
fn test_somatic_structural_variant_parity() {
    // TODO: Implement once somatic structural-variant fixture data is available.
    // Expected flow:
    // 1. Select a tumor/normal region that exercises somatic SV behavior.
    // 2. Run Java and Rust with matching flags.
    // 3. Compare output byte-for-byte.
    todo!("Somatic structural-variant parity test not yet implemented");
}
