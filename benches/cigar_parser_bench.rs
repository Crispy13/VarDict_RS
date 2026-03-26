use std::path::{Path, PathBuf};
use std::sync::{Arc, OnceLock};

use anyhow::{Result, anyhow};
use criterion::{BatchSize, Criterion, black_box, criterion_group, criterion_main};
use rust_htslib::bam::Record;
use vardict_rs::conf::Configuration;
use vardict_rs::data::bam_reader::BamReader;
use vardict_rs::data::reference::Reference;
use vardict_rs::data::region::Region;
use vardict_rs::data::shared_reference::{SharedReferenceHandle, load_shared_reference_chroms};
use vardict_rs::mods::cigar_parser::CigarParser;
use vardict_rs::scopedata::global_read_only_scope::{GlobalReadOnlyScope, INSTANCE};

struct CigarParserBenchInput {
    region: Region,
    reference: Reference,
    records: Vec<Record>,
    instance: Arc<GlobalReadOnlyScope>,
}

static CIGAR_PARSER_BENCH_INPUT: OnceLock<CigarParserBenchInput> = OnceLock::new();

fn benchmark_input() -> &'static CigarParserBenchInput {
    CIGAR_PARSER_BENCH_INPUT.get_or_init(|| {
        build_cigar_parser_bench_input().unwrap_or_else(|error| {
            panic!("failed to initialize cigar parser benchmark input: {error}")
        })
    })
}

fn build_cigar_parser_bench_input() -> Result<CigarParserBenchInput> {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let reference_path = manifest_dir.join("testdata/hs37d5.fa");
    let bam_path = manifest_dir.join("testdata/test_168714.bam");

    ensure_exists(&reference_path)?;
    ensure_exists(&bam_path)?;

    let region = Region::new("20".to_string(), 168600, 168800, String::new());
    let shared_reference =
        load_shared_reference_chroms(reference_path.as_path(), &["20"])?;
    let scope = initialize_scope(&shared_reference, &bam_path, true);
    let reference = build_reference(&region, &shared_reference, &scope)?;
    let (records, _) = load_cached_records(&bam_path, &region)?;

    Ok(CigarParserBenchInput {
        region,
        reference,
        records,
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
        Err(anyhow!("required benchmark input is missing: {}", path.display()))
    }
}

fn build_reference(
    region: &Region,
    shared_reference: &SharedReferenceHandle,
    scope: &GlobalReadOnlyScope,
) -> Result<Reference> {
    let extend =
        (scope.conf.number_nucleotide_to_extend + scope.conf.reference_extension).max(0) as usize;
    let extended_start = if region.start() > extend {
        region.start() - extend
    } else {
        1
    };
    let mut extended_end = region.end() + extend;
    if let Some(&chr_len) = scope.chr_lens.get(region.chr()) {
        if extended_end > chr_len {
            extended_end = chr_len;
        }
    }

    let ref_seq = shared_reference
        .get_subseq(region.chr(), extended_start, extended_end)
        .ok_or_else(|| {
            anyhow!(
                "failed to get reference for {}:{}-{}",
                region.chr(),
                extended_start,
                extended_end
            )
        })?
        .to_vec();

    let mut reference = Reference::new_with_start(ref_seq, extended_start as i64);
    let chr_len = scope.chr_lens.get(region.chr()).copied();
    reference.build_seed_map(extended_end as i64, chr_len);
    Ok(reference)
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

fn bench_cigar_parser(c: &mut Criterion) {
    let input = benchmark_input();

    c.bench_function("cigar_parser/process_records_chr20_168600_168800", |b| {
        b.iter_batched(
            || input.records.clone(),
            |mut records| {
                let mut parser = CigarParser::new(
                    input.region.clone(),
                    input.reference.clone(),
                    Arc::clone(&input.instance),
                );
                parser.process_records(records.iter_mut()).unwrap_or_else(|error| {
                    panic!("cigar parser benchmark failed: {error}")
                });
                black_box(parser.get_non_insertion_vars().len());
            },
            BatchSize::LargeInput,
        );
    });
}

criterion_group!(cigar_parser_benches, bench_cigar_parser);
criterion_main!(cigar_parser_benches);