//! Benchmark: VecMap vs HashMap at various entry counts.
//!
//! Tests the actual types used in the pipeline:
//!   RawVarMap = VecMap<VarDesc, RawVariant>  (current)
//!   vs HashMap<VarDesc, RawVariant, FxBuildHasher> (java-like original)
//!
//! Run: cargo bench --bench vecmap_vs_hashmap_bench

use std::collections::HashMap;
use std::hint::black_box;

use criterion::{BenchmarkId, Criterion, criterion_group, criterion_main};
use rustc_hash::FxBuildHasher;
use smallvec::SmallVec;
use vardict_rs::utils::vec_map::VecMap;
use vardict_rs::variants::variants::{VarDesc, Variant as RawVariant};

/// Entry counts to benchmark — covers typical (1-4) through pathological (32+)
const SIZES: &[usize] = &[1, 2, 4, 8, 16, 32];

// ---------------------------------------------------------------------------
// Helpers to generate realistic keys and values
// ---------------------------------------------------------------------------

fn make_key(i: usize) -> VarDesc {
    match i % 5 {
        0 => VarDesc::snv_key(b"ACGT"[i % 4]),
        1 => VarDesc::insertion(&format!("INS{i}").into_bytes()),
        2 => VarDesc::deletion((i as u32) + 1),
        3 => VarDesc::complex(b"ACG", &format!("T{i}").into_bytes()),
        _ => VarDesc::Raw {
            desc: SmallVec::from_slice(&format!("RAW{i}").into_bytes()),
        },
    }
}

fn make_value(i: usize) -> RawVariant {
    RawVariant {
        alt_depth: i + 1,
        alt_depth_fwd: i,
        alt_depth_rev: 1,
        mean_pos: (i as f64) * 10.0,
        mean_qual: 35.0 + (i as f64),
        mean_mapq: 40.0,
        nm: 1.0,
        low_qual_read_cnt: 0,
        high_qual_read_cnt: i + 1,
        pstd: i > 0,
        qstd: i > 1,
        pp: i * 50,
        pq: 30.0,
        extra_cnt: 0,
    }
}

fn make_keys(n: usize) -> Vec<VarDesc> {
    (0..n).map(make_key).collect()
}

// ---------------------------------------------------------------------------
// Pre-filled map builders
// ---------------------------------------------------------------------------

fn prefill_vecmap(keys: &[VarDesc]) -> VecMap<VarDesc, RawVariant> {
    let mut m = VecMap::with_capacity(keys.len());
    for (i, k) in keys.iter().enumerate() {
        m.insert(k.clone(), make_value(i));
    }
    m
}

fn prefill_hashmap(keys: &[VarDesc]) -> HashMap<VarDesc, RawVariant, FxBuildHasher> {
    let mut m = HashMap::with_capacity_and_hasher(keys.len(), FxBuildHasher);
    for (i, k) in keys.iter().enumerate() {
        m.insert(k.clone(), make_value(i));
    }
    m
}

// ---------------------------------------------------------------------------
// Benchmarks
// ---------------------------------------------------------------------------

/// Benchmark: insert N entries into an empty map
fn bench_insert(c: &mut Criterion) {
    let mut group = c.benchmark_group("insert");
    for &n in SIZES {
        let keys = make_keys(n);

        group.bench_with_input(BenchmarkId::new("VecMap", n), &n, |b, &n| {
            b.iter(|| {
                let mut m = VecMap::with_capacity(n);
                for (i, k) in keys.iter().enumerate() {
                    m.insert(k.clone(), make_value(i));
                }
                black_box(&m);
            });
        });

        group.bench_with_input(BenchmarkId::new("HashMap", n), &n, |b, &n| {
            b.iter(|| {
                let mut m = HashMap::with_capacity_and_hasher(n, FxBuildHasher);
                for (i, k) in keys.iter().enumerate() {
                    m.insert(k.clone(), make_value(i));
                }
                black_box(&m);
            });
        });
    }
    group.finish();
}

/// Benchmark: look up every key in a pre-filled map (hit)
fn bench_get_hit(c: &mut Criterion) {
    let mut group = c.benchmark_group("get_hit");
    for &n in SIZES {
        let keys = make_keys(n);
        let vm = prefill_vecmap(&keys);
        let hm = prefill_hashmap(&keys);

        group.bench_with_input(BenchmarkId::new("VecMap", n), &n, |b, _| {
            b.iter(|| {
                for k in &keys {
                    black_box(vm.get(k));
                }
            });
        });

        group.bench_with_input(BenchmarkId::new("HashMap", n), &n, |b, _| {
            b.iter(|| {
                for k in &keys {
                    black_box(hm.get(k));
                }
            });
        });
    }
    group.finish();
}

/// Benchmark: look up a key that does NOT exist (miss)
fn bench_get_miss(c: &mut Criterion) {
    let mut group = c.benchmark_group("get_miss");
    let miss_key = VarDesc::Raw {
        desc: SmallVec::from_slice(b"NONEXISTENT"),
    };

    for &n in SIZES {
        let keys = make_keys(n);
        let vm = prefill_vecmap(&keys);
        let hm = prefill_hashmap(&keys);

        group.bench_with_input(BenchmarkId::new("VecMap", n), &n, |b, _| {
            b.iter(|| {
                black_box(vm.get(&miss_key));
            });
        });

        group.bench_with_input(BenchmarkId::new("HashMap", n), &n, |b, _| {
            b.iter(|| {
                black_box(hm.get(&miss_key));
            });
        });
    }
    group.finish();
}

/// Benchmark: entry().or_default() pattern (the hot path in CigarParser)
fn bench_entry_or_default(c: &mut Criterion) {
    let mut group = c.benchmark_group("entry_or_default");
    for &n in SIZES {
        let keys = make_keys(n);

        group.bench_with_input(BenchmarkId::new("VecMap", n), &n, |b, _| {
            b.iter_batched(
                || prefill_vecmap(&keys),
                |mut m| {
                    // Access existing entries + one new entry
                    for k in &keys {
                        let v = m.entry(k.clone()).or_default();
                        v.alt_depth += 1;
                    }
                    let new_key = VarDesc::Raw {
                        desc: SmallVec::from_slice(b"NEW"),
                    };
                    m.entry(new_key).or_default();
                    black_box(&m);
                },
                criterion::BatchSize::SmallInput,
            );
        });

        group.bench_with_input(BenchmarkId::new("HashMap", n), &n, |b, _| {
            b.iter_batched(
                || prefill_hashmap(&keys),
                |mut m| {
                    for k in &keys {
                        let v = m.entry(k.clone()).or_default();
                        v.alt_depth += 1;
                    }
                    let new_key = VarDesc::Raw {
                        desc: SmallVec::from_slice(b"NEW"),
                    };
                    m.entry(new_key).or_default();
                    black_box(&m);
                },
                criterion::BatchSize::SmallInput,
            );
        });
    }
    group.finish();
}

/// Benchmark: iterate all entries and sum a field
fn bench_iterate(c: &mut Criterion) {
    let mut group = c.benchmark_group("iterate");
    for &n in SIZES {
        let keys = make_keys(n);
        let vm = prefill_vecmap(&keys);
        let hm = prefill_hashmap(&keys);

        group.bench_with_input(BenchmarkId::new("VecMap", n), &n, |b, _| {
            b.iter(|| {
                let sum: usize = vm.iter().map(|(_, v)| v.alt_depth).sum();
                black_box(sum);
            });
        });

        group.bench_with_input(BenchmarkId::new("HashMap", n), &n, |b, _| {
            b.iter(|| {
                let sum: usize = hm.iter().map(|(_, v)| v.alt_depth).sum();
                black_box(sum);
            });
        });
    }
    group.finish();
}

/// Benchmark: remove an entry from the middle
fn bench_remove(c: &mut Criterion) {
    let mut group = c.benchmark_group("remove");
    for &n in SIZES {
        let keys = make_keys(n);
        let mid = n / 2;

        group.bench_with_input(BenchmarkId::new("VecMap", n), &n, |b, _| {
            b.iter_batched(
                || prefill_vecmap(&keys),
                |mut m| {
                    black_box(m.remove(&keys[mid]));
                },
                criterion::BatchSize::SmallInput,
            );
        });

        group.bench_with_input(BenchmarkId::new("HashMap", n), &n, |b, _| {
            b.iter_batched(
                || prefill_hashmap(&keys),
                |mut m| {
                    black_box(m.remove(&keys[mid]));
                },
                criterion::BatchSize::SmallInput,
            );
        });
    }
    group.finish();
}

/// Benchmark: memory size of map instance (reported via custom measurement)
fn bench_mem_size(c: &mut Criterion) {
    let mut group = c.benchmark_group("mem_size_report");
    // This isn't a timing benchmark — just prints sizes for reference.
    // We use a trivial iter to let criterion run it.
    for &n in SIZES {
        let keys = make_keys(n);
        let vm = prefill_vecmap(&keys);
        let hm = prefill_hashmap(&keys);

        let vm_size =
            std::mem::size_of_val(&vm) + vm.len() * std::mem::size_of::<(VarDesc, RawVariant)>();
        let hm_size = std::mem::size_of_val(&hm)
            + hm.capacity() * (std::mem::size_of::<(VarDesc, RawVariant)>() + 1); // +1 for control byte

        // Print sizes as a side effect so they show up in bench output
        eprintln!(
            "[mem @{n} entries] VecMap: {vm_size} bytes, HashMap: {hm_size} bytes, ratio: {:.2}x",
            hm_size as f64 / vm_size as f64
        );

        group.bench_with_input(BenchmarkId::new("VecMap", n), &n, |b, _| {
            b.iter(|| black_box(vm_size));
        });
        group.bench_with_input(BenchmarkId::new("HashMap", n), &n, |b, _| {
            b.iter(|| black_box(hm_size));
        });
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_insert,
    bench_get_hit,
    bench_get_miss,
    bench_entry_or_default,
    bench_iterate,
    bench_remove,
    bench_mem_size,
);
criterion_main!(benches);
