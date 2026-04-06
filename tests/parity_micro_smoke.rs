use std::path::PathBuf;
use std::process::Command;

fn normalize_newlines(text: &str) -> String {
    text.replace("\r\n", "\n")
}

fn assert_output_matches(
    config: &str,
    chrom: &str,
    start: u32,
    end: u32,
    expected_fixture: &str,
    actual_stdout: &str,
) {
    let expected = normalize_newlines(expected_fixture);
    let actual = normalize_newlines(actual_stdout);
    let expected_lines: Vec<&str> = expected.lines().collect();
    let actual_lines: Vec<&str> = actual.lines().collect();

    for line_index in 0..expected_lines.len().max(actual_lines.len()) {
        let expected_line = expected_lines.get(line_index).copied();
        let actual_line = actual_lines.get(line_index).copied();

        match (expected_line, actual_line) {
            (Some(expected_line), Some(actual_line)) if expected_line == actual_line => continue,
            (Some(expected_line), Some(actual_line)) => {
                panic!(
                    "fixture mismatch for {config} {chrom}:{start}-{end} at line {}\nexpected: {}\nactual:   {}",
                    line_index + 1,
                    expected_line,
                    actual_line,
                );
            }
            (Some(expected_line), None) => {
                panic!(
                    "fixture mismatch for {config} {chrom}:{start}-{end}: actual output ended at line {}\nexpected extra line: {}",
                    line_index + 1,
                    expected_line,
                );
            }
            (None, Some(actual_line)) => {
                panic!(
                    "fixture mismatch for {config} {chrom}:{start}-{end}: actual output has extra line {}\nactual extra line: {}",
                    line_index + 1,
                    actual_line,
                );
            }
            (None, None) => break,
        }
    }
}

fn run_smoke(config: &str, chrom: &str, start: u32, end: u32, expected_fixture: &str) {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let ref_path = manifest_dir.join("testdata/hs37d5.fa");
    let bam_path =
        manifest_dir.join("testdata/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam");

    assert!(ref_path.exists(), "missing reference: {}", ref_path.display());
    assert!(bam_path.exists(), "missing BAM: {}", bam_path.display());
    let region = format!("{chrom}:{start}-{end}");

    let mut command = Command::new(env!("CARGO_BIN_EXE_vardict"));
    command.current_dir(&manifest_dir);
    command.arg("-G").arg(&ref_path);
    command.arg("-b").arg(&bam_path);
    command.arg("-N").arg("smoke_sample");
    command.arg("-R").arg(&region);

    match config {
        "default" => {}
        "nosv" => {
            command.arg("-U");
        }
        other => panic!("unsupported smoke config: {other}"),
    }

    let output = command
        .output()
        .expect("failed to run vardict binary");

    assert!(
        output.status.success(),
        "vardict failed for {config} {chrom}:{start}-{end}\nstatus: {}\nstderr:\n{}\nstdout:\n{}",
        output.status,
        String::from_utf8_lossy(&output.stderr),
        String::from_utf8_lossy(&output.stdout),
    );

    let actual_stdout = String::from_utf8(output.stdout).expect("vardict stdout was not UTF-8");
    assert_output_matches(config, chrom, start, end, expected_fixture, &actual_stdout);
}

macro_rules! parity_micro_smoke_tests {
    ($( $name:ident => ($config:literal, $chrom:literal, $start:literal, $end:literal, $fixture:literal) ),+ $(,)?) => {
        $(
            #[test]
            fn $name() {
                run_smoke($config, $chrom, $start, $end, include_str!($fixture));
            }
        )+
    };
}

parity_micro_smoke_tests!(
    parity_micro_smoke_default_chr20_126269_126333 => (
        "default",
        "20",
        126269,
        126333,
        "parity/fixtures/default_chr20_126269_126333.expected.tsv"
    ),
    parity_micro_smoke_default_chr20_168500_168800 => (
        "default",
        "20",
        168500,
        168800,
        "parity/fixtures/default_chr20_168500_168800.expected.tsv"
    ),
    parity_micro_smoke_default_chr20_168600_168800 => (
        "default",
        "20",
        168600,
        168800,
        "parity/fixtures/default_chr20_168600_168800.expected.tsv"
    ),
    parity_micro_smoke_default_chr20_25456879_25457078 => (
        "default",
        "20",
        25456879,
        25457078,
        "parity/fixtures/default_chr20_25456879_25457078.expected.tsv"
    ),
    parity_micro_smoke_default_chr20_30000000_30000300 => (
        "default",
        "20",
        30000000,
        30000300,
        "parity/fixtures/default_chr20_30000000_30000300.expected.tsv"
    ),
    parity_micro_smoke_default_chr22_40000000_40000300 => (
        "default",
        "22",
        40000000,
        40000300,
        "parity/fixtures/default_chr22_40000000_40000300.expected.tsv"
    ),
    parity_micro_smoke_default_chr22_42522500_42522800 => (
        "default",
        "22",
        42522500,
        42522800,
        "parity/fixtures/default_chr22_42522500_42522800.expected.tsv"
    ),
    parity_micro_smoke_default_chr_mt_1_300 => (
        "default",
        "MT",
        1,
        300,
        "parity/fixtures/default_chrMT_1_300.expected.tsv"
    ),
    parity_micro_smoke_default_chr_mt_300_600 => (
        "default",
        "MT",
        300,
        600,
        "parity/fixtures/default_chrMT_300_600.expected.tsv"
    ),
    parity_micro_smoke_nosv_chr20_126269_126333 => (
        "nosv",
        "20",
        126269,
        126333,
        "parity/fixtures/nosv_chr20_126269_126333.expected.tsv"
    ),
    parity_micro_smoke_nosv_chr20_168500_168800 => (
        "nosv",
        "20",
        168500,
        168800,
        "parity/fixtures/nosv_chr20_168500_168800.expected.tsv"
    ),
    parity_micro_smoke_nosv_chr20_168600_168800 => (
        "nosv",
        "20",
        168600,
        168800,
        "parity/fixtures/nosv_chr20_168600_168800.expected.tsv"
    ),
    parity_micro_smoke_nosv_chr20_25456879_25457078 => (
        "nosv",
        "20",
        25456879,
        25457078,
        "parity/fixtures/nosv_chr20_25456879_25457078.expected.tsv"
    ),
    parity_micro_smoke_nosv_chr20_30000000_30000300 => (
        "nosv",
        "20",
        30000000,
        30000300,
        "parity/fixtures/nosv_chr20_30000000_30000300.expected.tsv"
    ),
    parity_micro_smoke_nosv_chr22_40000000_40000300 => (
        "nosv",
        "22",
        40000000,
        40000300,
        "parity/fixtures/nosv_chr22_40000000_40000300.expected.tsv"
    ),
    parity_micro_smoke_nosv_chr22_42522500_42522800 => (
        "nosv",
        "22",
        42522500,
        42522800,
        "parity/fixtures/nosv_chr22_42522500_42522800.expected.tsv"
    ),
    parity_micro_smoke_nosv_chr_mt_1_300 => (
        "nosv",
        "MT",
        1,
        300,
        "parity/fixtures/nosv_chrMT_1_300.expected.tsv"
    ),
    parity_micro_smoke_nosv_chr_mt_300_600 => (
        "nosv",
        "MT",
        300,
        600,
        "parity/fixtures/nosv_chrMT_300_600.expected.tsv"
    )
);