//! HG002 exome BAM parity regression tests.
//! Each test captures a known M12 sweep parity failure against cached Java shards.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

const HG002_BAM: &str = "testdata/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam";
const REFERENCE: &str = "testdata/hs37d5.fa";
const SAMPLE: &str = "Sample_Diag-excap51-HG002-EEogPU";

#[derive(Debug, Clone, Copy)]
struct Hg002ParityCase {
    label: &'static str,
    chrom: &'static str,
    shard: &'static str,
    region: &'static str,
    extra_args: &'static [&'static str],
    expected_failure: &'static str,
}

#[derive(Debug)]
struct RawMismatchDiagnostic {
    line_index: usize,
    reason: String,
    java_line: Option<String>,
    rust_line: Option<String>,
    java_line_count: usize,
    rust_line_count: usize,
}

const RC1_FREQ_LOW_CHR4: Hg002ParityCase = Hg002ParityCase {
    label: "freq-low",
    chrom: "4",
    shard: "146",
    region: "4:145000001-146000000",
    extra_args: &["-f", "0.001"],
    expected_failure: "RC1 SV DUP depth inflation: Rust depth=288 vs Java depth=86 around 4:144801519-145041686",
};

const RC1_NOREALIGN_CHR4: Hg002ParityCase = Hg002ParityCase {
    label: "no-realign",
    chrom: "4",
    shard: "146",
    region: "4:145000001-146000000",
    extra_args: &["-k", "0"],
    expected_failure: "RC1 SV DUP depth inflation persists in the no-realign config",
};

const RC1_FISHER_CHR4: Hg002ParityCase = Hg002ParityCase {
    label: "fisher",
    chrom: "4",
    shard: "146",
    region: "4:145000001-146000000",
    extra_args: &["--fisher"],
    expected_failure: "RC1 SV DUP depth inflation persists when Fisher columns are enabled",
};

const RC2_MT_NOREALIGN: Hg002ParityCase = Hg002ParityCase {
    label: "no-realign",
    chrom: "MT",
    shard: "001",
    region: "MT:1-16569",
    extra_args: &["-k", "0"],
    expected_failure: "RC2 chrMT SV alt-depth mismatch: AltDepth 416->308 and AltFwdReads 207->99",
};

const RC3_FREQ_LOW_CHR15: Hg002ParityCase = Hg002ParityCase {
    label: "freq-low",
    chrom: "15",
    shard: "067",
    region: "15:66000001-67000000",
    extra_args: &["-f", "0.001"],
    expected_failure: "RC3 off-by-one depth mismatch: Rust depth=354 vs Java depth=355 on the complex variant",
};

const RC4_NOREALIGN_CHR14: Hg002ParityCase = Hg002ParityCase {
    label: "no-realign",
    chrom: "14",
    shard: "107",
    region: "14:106000001-107000000",
    extra_args: &["-k", "0"],
    expected_failure: "RC4 extra DEL materialization: Java has 1549 lines but Rust emits 1550",
};

const RC5A_NOREALIGN_CHR13: Hg002ParityCase = Hg002ParityCase {
    label: "no-realign",
    chrom: "13",
    shard: "103",
    region: "13:102000001-103000000",
    extra_args: &["-k", "0"],
    expected_failure: "RC5a no-realign row ordering mismatch: Java SNV at 102334267 where Rust emits DEL at 102305594",
};

const RC5B_NOREALIGN_CHR15: Hg002ParityCase = Hg002ParityCase {
    label: "no-realign",
    chrom: "15",
    shard: "042",
    region: "15:41000001-42000000",
    extra_args: &["-k", "0"],
    expected_failure: "RC5b INV genotype inflation: Rust genotype string contains extra sequence",
};

fn manifest_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn normalize_newlines(text: &str) -> String {
    text.replace("\r\n", "\n")
}

fn java_output_path(case: Hg002ParityCase) -> PathBuf {
    manifest_dir()
        .join("tmp")
        .join("hg002_parity")
        .join(case.label)
        .join(case.chrom)
        .join("java")
        .join(format!("shard_{}.tsv.gz", case.shard))
}

fn create_hard_link_if_missing(source: &Path, target: &Path) -> Result<(), String> {
    if target.exists() {
        return Ok(());
    }

    match fs::hard_link(source, target) {
        Ok(()) => Ok(()),
        Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => Ok(()),
        Err(error) => Err(format!(
            "Failed to hardlink {} -> {}: {}",
            source.display(),
            target.display(),
            error
        )),
    }
}

fn prepare_hg002_bam_with_index() -> Result<PathBuf, String> {
    let manifest_dir = manifest_dir();
    let source_bam = manifest_dir.join(HG002_BAM);
    let source_bai = source_bam.with_extension("bai");

    if !source_bam.exists() {
        return Err(format!("Missing HG002 BAM: {}", source_bam.display()));
    }
    if !source_bai.exists() {
        return Err(format!("Missing HG002 BAI: {}", source_bai.display()));
    }

    let target_dir = manifest_dir.join("tmp").join("hg002_parity_testdata");
    fs::create_dir_all(&target_dir)
        .map_err(|error| format!("Failed to create {}: {}", target_dir.display(), error))?;

    let bam_name = source_bam
        .file_name()
        .ok_or_else(|| format!("HG002 BAM has no file name: {}", source_bam.display()))?;
    let target_bam = target_dir.join(bam_name);
    let target_bam_bai = PathBuf::from(format!("{}.bai", target_bam.display()));

    create_hard_link_if_missing(&source_bam, &target_bam)?;
    create_hard_link_if_missing(&source_bai, &target_bam_bai)?;

    Ok(target_bam)
}

fn run_rust_case(case: Hg002ParityCase) -> Result<Vec<String>, String> {
    let manifest_dir = manifest_dir();
    let bam_path = prepare_hg002_bam_with_index()?;
    let reference_path = manifest_dir.join(REFERENCE);

    if !reference_path.exists() {
        return Err(format!("Missing reference FASTA: {}", reference_path.display()));
    }

    let mut command = Command::new(env!("CARGO_BIN_EXE_vardict"));
    command
        .current_dir(&manifest_dir)
        .arg("-G")
        .arg(&reference_path)
        .arg("-b")
        .arg(&bam_path)
        .arg("-N")
        .arg(SAMPLE)
        .arg("-R")
        .arg(case.region);

    for arg in case.extra_args {
        command.arg(arg);
    }

    let output = command
        .output()
        .map_err(|error| format!("Failed to run vardict for {}: {}", case.region, error))?;

    if !output.status.success() {
        return Err(format!(
            "vardict failed for {} [{} {} shard {}]\nstatus: {}\nstderr:\n{}\nstdout:\n{}",
            case.region,
            case.label,
            case.chrom,
            case.shard,
            output.status,
            String::from_utf8_lossy(&output.stderr),
            String::from_utf8_lossy(&output.stdout),
        ));
    }

    let stdout = String::from_utf8(output.stdout)
        .map_err(|error| format!("vardict stdout was not UTF-8 for {}: {}", case.region, error))?;
    Ok(normalize_newlines(&stdout)
        .lines()
        .map(|line| line.to_string())
        .collect())
}

fn read_gzip_lines(path: &Path) -> Result<Vec<String>, String> {
    let output = Command::new("gzip")
        .arg("-dc")
        .arg(path)
        .output()
        .map_err(|error| format!("Failed to run gzip for {}: {}", path.display(), error))?;

    if !output.status.success() {
        return Err(format!(
            "gzip failed for {}\nstatus: {}\nstderr:\n{}",
            path.display(),
            output.status,
            String::from_utf8_lossy(&output.stderr),
        ));
    }

    let text = String::from_utf8(output.stdout)
        .map_err(|error| format!("Java shard was not UTF-8 for {}: {}", path.display(), error))?;
    Ok(normalize_newlines(&text)
        .lines()
        .map(|line| line.to_string())
        .collect())
}

fn first_raw_mismatch(
    expected_lines: &[String],
    rust_output: &[String],
) -> Option<RawMismatchDiagnostic> {
    let shared = expected_lines.len().min(rust_output.len());
    for index in 0..shared {
        if expected_lines[index] != rust_output[index] {
            return Some(RawMismatchDiagnostic {
                line_index: index,
                reason: "line_content_mismatch".to_string(),
                java_line: Some(expected_lines[index].clone()),
                rust_line: Some(rust_output[index].clone()),
                java_line_count: expected_lines.len(),
                rust_line_count: rust_output.len(),
            });
        }
    }

    if expected_lines.len() != rust_output.len() {
        return Some(RawMismatchDiagnostic {
            line_index: shared,
            reason: format!(
                "line_count_mismatch(java={},rust={})",
                expected_lines.len(),
                rust_output.len()
            ),
            java_line: expected_lines.get(shared).cloned(),
            rust_line: rust_output.get(shared).cloned(),
            java_line_count: expected_lines.len(),
            rust_line_count: rust_output.len(),
        });
    }

    None
}

fn case_mismatch_message(case: Hg002ParityCase) -> Option<String> {
    let java_path = java_output_path(case);
    assert!(
        java_path.exists(),
        "Missing cached Java shard for {} [{} {} shard {}]: {}",
        case.region,
        case.label,
        case.chrom,
        case.shard,
        java_path.display()
    );

    let expected_lines = match read_gzip_lines(&java_path) {
        Ok(lines) => lines,
        Err(error) => return Some(format!("Failed to read cached Java shard: {}", error)),
    };
    let rust_output = match run_rust_case(case) {
        Ok(lines) => lines,
        Err(error) => return Some(format!("Failed to run Rust parity case: {}", error)),
    };

    if let Some(diag) = first_raw_mismatch(&expected_lines, &rust_output) {
        return Some(format!(
            "HG002 parity mismatch for {} [{} {} shard {}]\nexpected failure: {}\nreason: {}\nline: {}\njava_count: {}\nrust_count: {}\nJAVA: {}\nRUST: {}",
            case.region,
            case.label,
            case.chrom,
            case.shard,
            case.expected_failure,
            diag.reason,
            diag.line_index + 1,
            diag.java_line_count,
            diag.rust_line_count,
            diag.java_line.unwrap_or_else(|| "<none>".to_string()),
            diag.rust_line.unwrap_or_else(|| "<none>".to_string()),
        ));
    }

    None
}

fn assert_case_matches_java(case: Hg002ParityCase) {
    if let Some(message) = case_mismatch_message(case) {
        panic!("{}", message);
    }
}

#[test]
#[ignore = "known HG002 parity regression; requires full HG002 BAM + cached tmp/hg002_parity Java shards"]
fn test_hg002_rc1_sv_dup_depth_chr4_freq_low() {
    // Expected failure until the chr4 SV DUP depth inflation bug is fixed across all three configs.
    let failures = [RC1_FREQ_LOW_CHR4, RC1_NOREALIGN_CHR4, RC1_FISHER_CHR4]
        .into_iter()
        .filter_map(case_mismatch_message)
        .collect::<Vec<_>>();

    if !failures.is_empty() {
        panic!("{}", failures.join("\n\n"));
    }
}

#[test]
#[ignore = "known HG002 parity regression; requires full HG002 BAM + cached tmp/hg002_parity Java shards"]
fn test_hg002_rc2_sv_alt_depth_chrmt_norealign() {
    // Expected failure until the chrMT SV alt-depth accounting matches Java.
    assert_case_matches_java(RC2_MT_NOREALIGN);
}

#[test]
#[ignore = "known HG002 parity regression; requires full HG002 BAM + cached tmp/hg002_parity Java shards"]
fn test_hg002_rc3_depth_offbyone_chr15_freq_low() {
    // Expected failure until the chr15 complex-variant depth matches Java exactly.
    assert_case_matches_java(RC3_FREQ_LOW_CHR15);
}

#[test]
#[ignore = "known HG002 parity regression; requires full HG002 BAM + cached tmp/hg002_parity Java shards"]
fn test_hg002_rc4_extra_del_chr14_norealign() {
    // Expected failure until the extra chr14 DEL row is no longer materialized in Rust.
    assert_case_matches_java(RC4_NOREALIGN_CHR14);
}

#[test]
#[ignore = "known HG002 parity regression; requires full HG002 BAM + cached tmp/hg002_parity Java shards"]
fn test_hg002_rc5a_variant_order_chr13_norealign() {
    // Expected failure until Rust emits chr13 rows in the same order as Java.
    assert_case_matches_java(RC5A_NOREALIGN_CHR13);
}

#[test]
#[ignore = "known HG002 parity regression; requires full HG002 BAM + cached tmp/hg002_parity Java shards"]
fn test_hg002_rc5b_inv_genotype_chr15_norealign() {
    // Expected failure until the chr15 INV genotype string stops inflating in Rust.
    assert_case_matches_java(RC5B_NOREALIGN_CHR15);
}