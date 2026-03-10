use std::env;
use std::fs;
use std::path::PathBuf;
use std::sync::Arc;
use std::time::{SystemTime, UNIX_EPOCH};

use crackle_kit::tracing::level_filters::LevelFilter;
use vardict_rs::conf::Configuration;
use vardict_rs::data::bam_reader::BamReader;
use vardict_rs::data::region::Region;
use vardict_rs::data::shared_reference::load_shared_reference_chroms;
use vardict_rs::mods::vardict_pipeline::VarDictPipeline;
use vardict_rs::scopedata::global_read_only_scope::{GlobalReadOnlyScope, INSTANCE};

fn test_log_level() -> LevelFilter {
    env::var("VARDICT_TEST_LOG")
        .ok()
        .and_then(|value| value.to_ascii_lowercase().parse::<LevelFilter>().ok())
        .unwrap_or(LevelFilter::WARN)
}

struct EnvGuard {
    key: &'static str,
    previous: Option<String>,
}

impl EnvGuard {
    fn new(key: &'static str, value: &str) -> Self {
        let previous = env::var(key).ok();
        unsafe {
            env::set_var(key, value);
        }
        Self { key, previous }
    }
}

impl Drop for EnvGuard {
    fn drop(&mut self) {
        match &self.previous {
            Some(value) => unsafe {
                env::set_var(self.key, value);
            },
            None => unsafe {
                env::remove_var(self.key);
            },
        }
    }
}

struct TempFileGuard {
    path: PathBuf,
}

impl TempFileGuard {
    fn new(path: PathBuf) -> Self {
        Self { path }
    }
}

impl Drop for TempFileGuard {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.path);
    }
}

#[test]
fn test_structural_variants_snapshot_fixture() {
    let _ = crackle_kit::tracing_kit::setup_logging_stderr_only_verbose(test_log_level());

    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let fixture_path =
        manifest_dir.join("tests/fixtures/structural_variants_chr20_168600_168800.jsonl");
    let bam_path = manifest_dir.join("testdata/test_168714.bam");
    let ref_path = manifest_dir.join("VarDictJava/tests/integration/reference/hs37d5.fa");

    assert!(fixture_path.exists(), "Missing fixture: {:?}", fixture_path);
    assert!(bam_path.exists(), "Missing BAM: {:?}", bam_path);
    assert!(ref_path.exists(), "Missing reference: {:?}", ref_path);

    let temp_name = format!(
        "vardict_structural_variants_snapshot_{}_{}.jsonl",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("time")
            .as_nanos()
    );
    let output_path = env::temp_dir().join(temp_name);

    let _env_guard = EnvGuard::new(
        "VARDICT_STRUCTURAL_VARIANTS_JSONL",
        output_path.to_str().expect("output path"),
    );
    let _temp_guard = TempFileGuard::new(output_path.clone());

    let shared_reference =
        load_shared_reference_chroms(ref_path.to_str().expect("reference path"), &["20"])
            .expect("failed to load reference");

    let mut scope = if let Some(existing) = INSTANCE.get() {
        existing.clone()
    } else {
        let mut conf = Configuration::default();
        conf.disable_sv = true;
        conf.perform_local_realignment = true;
        let mut scope = GlobalReadOnlyScope::default();
        scope.conf = conf;
        scope
    };
    scope.chr_lens = shared_reference.get_chromosome_lengths();
    let _ = INSTANCE.set(scope.clone());

    let instance = Arc::new(scope);

    let region = Region::new("20".to_string(), 168600, 168800, String::new());

    let mut bam_reader =
        BamReader::open(bam_path.to_str().expect("bam path")).expect("failed to open BAM");

    let pipeline = VarDictPipeline::new("test_sample")
        .with_min_frequency(0.01)
        .with_min_base_quality(22.5)
        .with_min_mapping_quality(0)
        .with_pileup(false);

    let _ = pipeline
        .process_region_from_bam(&region, &shared_reference, &mut bam_reader, instance)
        .expect("pipeline failed");

    let expected = fs::read_to_string(&fixture_path).expect("read fixture");
    let actual = fs::read_to_string(&output_path).expect("read snapshot");

    assert_eq!(actual, expected, "StructuralVariants snapshot mismatch");
}
