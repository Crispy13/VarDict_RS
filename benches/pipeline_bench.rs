use std::path::{Path, PathBuf};
use std::sync::{Arc, OnceLock};

use anyhow::{Result, anyhow};
use criterion::{BatchSize, Criterion, black_box, criterion_group, criterion_main};
use rust_htslib::bam::Record;
use vardict_rs::conf::Configuration;
use vardict_rs::data::bam_reader::BamReader;
use vardict_rs::data::region::Region;
use vardict_rs::data::shared_reference::{SharedReferenceHandle, load_shared_reference_chroms};
use vardict_rs::mods::vardict_pipeline::VarDictPipeline;
use vardict_rs::scopedata::global_read_only_scope::{GlobalReadOnlyScope, INSTANCE};

struct PipelineBenchInput {
    region: Region,
    shared_reference: SharedReferenceHandle,
    records: Vec<Record>,
    target_names: Vec<String>,
    pipeline: VarDictPipeline,
    instance: Arc<GlobalReadOnlyScope>,
}

static PIPELINE_BENCH_INPUT: OnceLock<PipelineBenchInput> = OnceLock::new();

fn benchmark_input() -> &'static PipelineBenchInput {
    PIPELINE_BENCH_INPUT.get_or_init(|| {
        build_pipeline_bench_input().unwrap_or_else(|error| {
            panic!("failed to initialize pipeline benchmark input: {error}")
        })
    })
}

fn build_pipeline_bench_input() -> Result<PipelineBenchInput> {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let reference_path = manifest_dir.join("testdata/hs37d5.fa");
    let bam_path = manifest_dir.join("testdata/test_168714.bam");

    ensure_exists(&reference_path)?;
    ensure_exists(&bam_path)?;

    let region = Region::new("20".to_string(), 168600, 168800, String::new());
    let shared_reference = load_shared_reference_chroms(reference_path.as_path(), &["20"])?;
    let scope = initialize_scope(&shared_reference, &bam_path, false);
    let (records, target_names) = load_cached_records(&bam_path, &region)?;
    let pipeline = VarDictPipeline::new("bench_sample")
        .with_min_frequency(0.01)
        .with_min_base_quality(22.5)
        .with_min_mapping_quality(0)
        .with_pileup(false);

    Ok(PipelineBenchInput {
        region,
        shared_reference,
        records,
        target_names,
        pipeline,
        instance: Arc::new(scope),
    })
}

fn initialize_scope(
    shared_reference: &SharedReferenceHandle,
    bam_path: &Path,
    disable_sv: bool,
) -> GlobalReadOnlyScope {
    let mut scope = if let Some(existing) = INSTANCE.get() {
        existing.clone()
    } else {
        let mut conf = Configuration::default();
        conf.disable_sv = disable_sv;
        conf.perform_local_realignment = true;

        let mut new_scope = GlobalReadOnlyScope::default();
        new_scope.conf = conf;
        new_scope
    };

    scope.conf.disable_sv = disable_sv;
    scope.conf.perform_local_realignment = true;
    scope.chr_lens = shared_reference.get_chromosome_lengths();
    scope.bam_paths = vec![bam_path.display().to_string()];

    let _ = INSTANCE.set(scope.clone());
    scope
}

fn ensure_exists(path: &Path) -> Result<()> {
    if path.exists() {
        Ok(())
    } else {
        Err(anyhow!(
            "required benchmark input is missing: {}",
            path.display()
        ))
    }
}

fn load_cached_records(bam_path: &Path, region: &Region) -> Result<(Vec<Record>, Vec<String>)> {
    let mut bam_reader = BamReader::open(bam_path)?;
    let target_names = bam_reader.target_names();
    bam_reader.fetch(region.chr(), region.start(), region.end())?;

    let mut records = Vec::new();
    let mut record = Record::new();
    while bam_reader.read(&mut record)? {
        records.push(record.clone());
    }

    Ok((records, target_names))
}

fn bench_pipeline(c: &mut Criterion) {
    let input = benchmark_input();

    c.bench_function(
        "vardict_pipeline/process_region_from_cached_records_chr20_168600_168800",
        |b| {
            b.iter_batched(
                || input.records.clone(),
                |records| {
                    let output_lines = input
                        .pipeline
                        .process_region_from_cached_records(
                            &input.region,
                            &input.shared_reference,
                            records,
                            &input.target_names,
                            Arc::clone(&input.instance),
                        )
                        .unwrap_or_else(|error| panic!("pipeline benchmark failed: {error}"));
                    black_box(output_lines.len());
                },
                BatchSize::LargeInput,
            );
        },
    );
}

criterion_group!(pipeline_benches, bench_pipeline);
criterion_main!(pipeline_benches);
