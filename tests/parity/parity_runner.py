#!/usr/bin/env python3
"""Single-config parity runner for one chromosome or an FAI-driven chromosome sweep."""

from __future__ import annotations

import argparse
import difflib
from itertools import islice
import json
import os
from concurrent.futures import Future, ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
import gzip
from pathlib import Path
import shutil
import sys
from typing import Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tests.parity.lib import comparator, config, reporter, resources, runner, shard


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_REF_FASTA = Path("testdata/hs37d5.fa")
DEFAULT_BAM_PATH = Path("testdata/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam")
DEFAULT_SAMPLE_NAME = "NA12878"
JAVA_JAR = Path("VarDictJava/build/libs/VarDict-1.8.3.jar")
DEFAULT_RUST_BIN = Path("target/debug-release/vardict")
DEFAULT_BASE_DIR = Path("tmp/na12878_parity")
DISK_ABORT_GB = 5.0

REF_FASTA = DEFAULT_REF_FASTA
REF_FASTA_FAI = Path(f"{DEFAULT_REF_FASTA}.fai")
BAM_PATH = DEFAULT_BAM_PATH
SAMPLE_NAME = DEFAULT_SAMPLE_NAME
BASE_DIR = DEFAULT_BASE_DIR


@dataclass(frozen=True, slots=True)
class CompareArtifact:
    status: str
    reason: str | None
    first_diff_line: int | None
    summary: str


@dataclass(frozen=True, slots=True)
class ChromosomeRunOutcome:
    chrom: str
    result: reporter.ChromResult
    elapsed_seconds: float


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run Java vs Rust parity for one option config on one chromosome.",
    )
    parser.add_argument("--bam", type=Path, default=DEFAULT_BAM_PATH, help=f"Path to BAM file (default: {DEFAULT_BAM_PATH})")
    parser.add_argument("--ref", type=Path, default=DEFAULT_REF_FASTA, help=f"Path to reference FASTA (default: {DEFAULT_REF_FASTA})")
    parser.add_argument("--sample", default=DEFAULT_SAMPLE_NAME, help=f"Sample name for -N (default: {DEFAULT_SAMPLE_NAME})")
    parser.add_argument("--base-dir", type=Path, default=DEFAULT_BASE_DIR, help=f"Output directory (default: {DEFAULT_BASE_DIR})")
    parser.add_argument("--chr", default="20", help="Chromosome name to run (default: 20)")
    parser.add_argument(
        "--all-chr",
        action="store_true",
        help="Run all chromosomes from the reference FAI sequentially",
    )
    parser.add_argument(
        "--chr-len",
        type=int,
        help="Override chromosome length instead of looking it up in the FAI",
    )
    parser.add_argument(
        "--config",
        help="Config ID from tests.parity.lib.config.TEST_MATRIX, for example T1-01",
    )
    parser.add_argument(
        "--config-label",
        help="Config label from tests.parity.lib.config.TEST_MATRIX, for example pileup",
    )
    parser.add_argument(
        "--opts",
        default="",
        help="Custom CLI flags passed to both Java and Rust; cannot be combined with --config or --config-label",
    )
    parser.add_argument(
        "--opts-label",
        default="default",
        help="Artifact directory label when using --opts (default: default)",
    )
    parser.add_argument("--freq", type=float, default=0.01, help="Default frequency threshold when flags omit -f")
    parser.add_argument("--parallel", type=int, default=10, help="Rust worker count (default: 10)")
    parser.add_argument("--java-parallel", type=int, default=10, help="Java worker count (default: 10)")
    parser.add_argument("--java-heap", default="8g", help="Java heap size per worker (default: 8g)")
    parser.add_argument("--timeout", type=int, default=600, help="Per-shard timeout in seconds (default: 600)")
    parser.add_argument("--retry", type=int, default=0, help="Retry transient shard failures this many times (default: 0)")
    parser.add_argument("--rust-only", action="store_true", help="Skip Java generation and require cached Java shards")
    parser.add_argument("--clean-rust", action="store_true", help="Remove only Rust and diff outputs before running")
    parser.add_argument("--clean-all", action="store_true", help="Delete Java cache as well during cleanup")
    parser.add_argument("--no-cleanup", action="store_true", help="Keep all shard outputs after comparison")
    parser.add_argument("--no-build", action="store_true", help="Skip cargo build and staleness checks")
    parser.add_argument(
        "--release",
        action="store_true",
        help="Use the repository debug-release binary path and build profile",
    )
    parser.add_argument("--rust-bin", help="Use a specific Rust binary path")
    parser.add_argument(
        "--no-stop",
        action="store_true",
        help="Continue through an all-chromosome sweep after failures",
    )
    parser.add_argument("--shard-size", type=int, default=1_000_000, help="Shard size in bases (default: 1000000)")
    parser.add_argument("--mem-warn-gb", type=float, default=4.0, help="Warn below this MemAvailable threshold")
    parser.add_argument("--mem-abort-gb", type=float, default=2.0, help="Abort below this MemAvailable threshold")
    parser.add_argument("--disk-warn-gb", type=float, default=10.0, help="Warn below this disk threshold")
    return parser


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def fail(message: str) -> int:
    print(f"ERROR: {message}", file=sys.stderr)
    return 1


def resolve_rust_bin(args: argparse.Namespace) -> Path:
    if args.release and args.rust_bin:
        raise ValueError("--release and --rust-bin cannot be used together")
    if args.rust_bin:
        return Path(args.rust_bin)
    return DEFAULT_RUST_BIN


def validate_args(args: argparse.Namespace) -> None:
    if args.parallel <= 0:
        raise ValueError("--parallel must be greater than zero")
    if args.java_parallel <= 0:
        raise ValueError("--java-parallel must be greater than zero")
    if args.timeout <= 0:
        raise ValueError("--timeout must be greater than zero")
    if args.retry < 0:
        raise ValueError("--retry must be zero or greater")
    if args.shard_size <= 0:
        raise ValueError("--shard-size must be greater than zero")
    if args.mem_warn_gb <= args.mem_abort_gb:
        raise ValueError("--mem-warn-gb must be greater than --mem-abort-gb")
    if args.disk_warn_gb <= DISK_ABORT_GB:
        raise ValueError(f"--disk-warn-gb must be greater than {DISK_ABORT_GB:.0f}")
    if args.chr_len is not None and args.chr_len <= 0:
        raise ValueError("--chr-len must be greater than zero")
    if (args.config or args.config_label) and args.opts:
        raise ValueError("--opts cannot be combined with --config or --config-label")


def build_env() -> dict[str, str]:
    env = os.environ.copy()
    for key in ("CFLAGS", "CXXFLAGS", "CPPFLAGS", "LDFLAGS"):
        env.pop(key, None)
    return env


def auto_build_if_needed(rust_bin: Path, *, no_build: bool) -> None:
    if no_build:
        return
    resolved_bin = project_path(rust_bin)
    if resolved_bin.exists():
        newer_source = any(path.stat().st_mtime > resolved_bin.stat().st_mtime for path in PROJECT_ROOT.glob("src/**/*.rs"))
        if newer_source:
            print(f"WARNING: {rust_bin} may be stale; rebuilding to keep parity artifacts current", file=sys.stderr)
        else:
            return

    import subprocess

    command = ["cargo", "build", "--profile", "debug-release"]
    print("Building Rust binary via: cargo build --profile debug-release")
    subprocess.run(command, cwd=PROJECT_ROOT, env=build_env(), check=True)


def project_path(path: Path) -> Path:
    return path if path.is_absolute() else PROJECT_ROOT / path


def reference_fai_path(ref_fasta: Path) -> Path:
    return Path(f"{ref_fasta}.fai")


def configure_runtime_paths(*, bam_path: Path, ref_fasta: Path, sample_name: str, base_dir: Path) -> None:
    global BAM_PATH, REF_FASTA, REF_FASTA_FAI, SAMPLE_NAME, BASE_DIR

    BAM_PATH = bam_path
    REF_FASTA = ref_fasta
    REF_FASTA_FAI = reference_fai_path(ref_fasta)
    SAMPLE_NAME = sample_name
    BASE_DIR = base_dir


def resolved_base_dir(base_dir: Path | None = None) -> Path:
    return BASE_DIR if base_dir is None else base_dir


def resolve_test_config(args: argparse.Namespace) -> config.TestConfig:
    if args.config:
        try:
            return config.TEST_MATRIX_BY_ID[args.config]
        except KeyError as exc:
            raise ValueError(f"Unknown config ID: {args.config}") from exc

    if args.config_label:
        matches = [entry for entry in config.TEST_MATRIX if entry.label == args.config_label]
        if not matches:
            raise ValueError(f"Unknown config label: {args.config_label}")
        if len(matches) > 1:
            raise ValueError(f"Config label is ambiguous: {args.config_label}")
        return matches[0]

    if args.opts:
        return config.TestConfig(config_id="CUSTOM", label=args.opts_label, flags=args.opts, tier="custom")

    return config.TestConfig(config_id="DEFAULT", label=args.opts_label, flags="", tier="custom")


def load_chrom_lengths() -> dict[str, int]:
    fai_path = project_path(REF_FASTA_FAI)
    if not fai_path.exists():
        raise FileNotFoundError(f"Reference FAI not found: {REF_FASTA_FAI}")
    return shard.load_fai(fai_path)


def rust_binary_mtime(rust_bin: Path) -> str | None:
    resolved = project_path(rust_bin)
    if not resolved.exists():
        return None
    return str(int(resolved.stat().st_mtime))


def chrom_dir(test_config: config.TestConfig, chrom: str, base_dir: Path | None = None) -> Path:
    return resolved_base_dir(base_dir) / test_config.label / chrom


def diff_dir(test_config: config.TestConfig, chrom: str, base_dir: Path | None = None) -> Path:
    return chrom_dir(test_config, chrom, base_dir) / "diff"


def shard_status_path(test_config: config.TestConfig, shard_info: shard.Shard, base_dir: Path | None = None) -> Path:
    return diff_dir(test_config, shard_info.chrom, base_dir) / f"shard_{shard_info.label}.status"


def shard_meta_path(test_config: config.TestConfig, shard_info: shard.Shard, base_dir: Path | None = None) -> Path:
    return diff_dir(test_config, shard_info.chrom, base_dir) / f"shard_{shard_info.label}.meta.json"


def shard_diff_path(test_config: config.TestConfig, shard_info: shard.Shard, base_dir: Path | None = None) -> Path:
    return diff_dir(test_config, shard_info.chrom, base_dir) / f"shard_{shard_info.label}.diff"


def shard_verified_path(test_config: config.TestConfig, shard_info: shard.Shard, base_dir: Path | None = None) -> Path:
    return diff_dir(test_config, shard_info.chrom, base_dir) / f"shard_{shard_info.label}.verified"


def chrom_verified_path(test_config: config.TestConfig, chrom: str, base_dir: Path | None = None) -> Path:
    return chrom_dir(test_config, chrom, base_dir) / ".verified"


def results_json_path(test_config: config.TestConfig, chrom: str, base_dir: Path | None = None) -> Path:
    return chrom_dir(test_config, chrom, base_dir) / "results.json"


def is_marker_current(marker_path: Path, rust_mtime: str | None) -> bool:
    if not marker_path.exists() or rust_mtime is None:
        return False
    try:
        payload = json.loads(marker_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return False
    return payload.get("rust_binary_mtime") == rust_mtime


def write_marker(marker_path: Path, *, rust_mtime: str | None, payload: dict[str, object]) -> None:
    marker_path.parent.mkdir(parents=True, exist_ok=True)
    body = {"rust_binary_mtime": rust_mtime, "timestamp": utc_now(), **payload}
    marker_path.write_text(json.dumps(body, indent=2, sort_keys=True), encoding="utf-8")


def remove_if_exists(path: Path) -> None:
    path.unlink(missing_ok=True)


def prepare_dirs(test_config: config.TestConfig, chrom: str, *, clean_rust: bool) -> None:
    root = chrom_dir(test_config, chrom)
    if clean_rust:
        shutil.rmtree(root / "rust", ignore_errors=True)
        shutil.rmtree(root / "diff", ignore_errors=True)
    for child in (root / "java", root / "rust", root / "diff"):
        child.mkdir(parents=True, exist_ok=True)


def check_resource_thresholds(args: argparse.Namespace, *, phase_label: str) -> None:
    messages = resources.check_resources(
        args.mem_warn_gb,
        args.mem_abort_gb,
        args.disk_warn_gb,
        path=PROJECT_ROOT / "tmp",
        abort_disk_gb=DISK_ABORT_GB,
    )
    aborts = [message for message in messages if message.startswith("ABORT:")]
    warnings = [message for message in messages if message.startswith("WARN:")]
    for message in warnings:
        print(f"WARNING: {phase_label}: {message}", file=sys.stderr)
    if aborts:
        raise RuntimeError(f"{phase_label}: {'; '.join(aborts)}")


def list_chromosomes(args: argparse.Namespace, lengths: dict[str, int]) -> list[str]:
    if args.all_chr:
        return list(lengths)
    return [args.chr]


def chrom_length(chrom: str, args: argparse.Namespace, lengths: dict[str, int]) -> int:
    if args.chr_len is not None and not args.all_chr:
        return args.chr_len
    try:
        return lengths[chrom]
    except KeyError as exc:
        raise ValueError(f"Chromosome '{chrom}' not found in {REF_FASTA_FAI}") from exc


def ensure_inputs_exist() -> None:
    for path in (REF_FASTA, BAM_PATH):
        if not project_path(path).exists():
            raise FileNotFoundError(f"Required file not found: {path}")


def phase_indices_for_java(
    parity: runner.ParityRunner,
    shards: Sequence[shard.Shard],
    test_config: config.TestConfig,
    rust_mtime: str | None,
) -> list[shard.Shard]:
    needed: list[shard.Shard] = []
    for shard_info in shards:
        if is_marker_current(shard_verified_path(test_config, shard_info), rust_mtime):
            continue
        if not parity.is_java_cached(shard_info, test_config):
            needed.append(shard_info)
    return needed


def phase_indices_for_rust(
    parity: runner.ParityRunner,
    shards: Sequence[shard.Shard],
    test_config: config.TestConfig,
    rust_mtime: str | None,
) -> list[shard.Shard]:
    needed: list[shard.Shard] = []
    for shard_info in shards:
        if is_marker_current(shard_verified_path(test_config, shard_info), rust_mtime):
            continue
        if parity.is_rust_stale(shard_info, test_config):
            needed.append(shard_info)
    return needed


def validate_java_cache(parity: runner.ParityRunner, shards: Sequence[shard.Shard], test_config: config.TestConfig) -> None:
    missing = [shard_info.label for shard_info in shards if not parity.is_java_cached(shard_info, test_config)]
    if missing:
        raise FileNotFoundError(
            f"--rust-only: {len(missing)} Java shard(s) missing for {shards[0].chrom if shards else 'unknown'}"
        )


def run_parallel_phase(
    phase_label: str,
    shards_to_run: Sequence[shard.Shard],
    *,
    max_workers: int,
    args: argparse.Namespace,
    executor_fn,
) -> dict[str, runner.RunResult]:
    results_by_label: dict[str, runner.RunResult] = {}
    if not shards_to_run:
        return results_by_label

    print(f"--- {phase_label}: {len(shards_to_run)} shards, {max_workers} workers ---")
    check_resource_thresholds(args, phase_label=phase_label)

    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        future_map: dict[Future[runner.RunResult], shard.Shard] = {}
        for shard_info in shards_to_run:
            check_resource_thresholds(args, phase_label=f"{phase_label} launch")
            future = pool.submit(executor_fn, shard_info)
            future_map[future] = shard_info

        completed = 0
        for future in as_completed(future_map):
            shard_info = future_map[future]
            result = future.result()
            results_by_label[shard_info.label] = result
            completed += 1
            if completed == len(shards_to_run) or completed % 10 == 0:
                print(f"  {phase_label}: {completed}/{len(shards_to_run)} shards complete")
            if not result.success:
                status = "timeout" if result.timed_out else f"exit {result.returncode}"
                print(
                    f"  {phase_label} failed for shard_{shard_info.label} ({shard.shard_region_string(shard_info)}): {status}",
                    file=sys.stderr,
                )
    return results_by_label


def compare_summary(
    compare_result: comparator.CompareResult,
    *,
    java_run: runner.RunResult | None,
    rust_run: runner.RunResult | None,
    timeout_seconds: int,
) -> CompareArtifact:
    if compare_result.status == "PASS":
        return CompareArtifact("PASS", None, None, "pass")
    if compare_result.status == "EMPTY":
        return CompareArtifact("EMPTY", None, None, "both outputs empty")
    if compare_result.reason == "missing_output":
        if java_run is not None and not java_run.success:
            summary = f"java timeout ({timeout_seconds}s)" if java_run.timed_out else f"java command failed (exit {java_run.returncode})"
            return CompareArtifact("FAIL", "java_missing", None, summary)
        if rust_run is not None and not rust_run.success:
            summary = f"rust timeout ({timeout_seconds}s)" if rust_run.timed_out else f"rust command failed (exit {rust_run.returncode})"
            return CompareArtifact("FAIL", "rust_missing", None, summary)
        return CompareArtifact("FAIL", "missing_output", None, "missing output")
    if compare_result.reason == "java_empty_suspect":
        summary = (
            "SUSPECT: java output empty but rust has "
            f"{compare_result.rust_line_count} lines (likely transient Java failure; delete java cache and re-run)"
        )
        return CompareArtifact("FAIL", compare_result.reason, compare_result.first_diff_line, summary)
    if compare_result.reason == "rust_empty_suspect":
        summary = (
            "SUSPECT: rust output empty but java has "
            f"{compare_result.java_line_count} lines (likely transient Rust failure; delete rust cache and re-run)"
        )
        return CompareArtifact("FAIL", compare_result.reason, compare_result.first_diff_line, summary)
    if compare_result.reason == "line_count_mismatch":
        summary = f"line_count_mismatch (java={compare_result.java_line_count}, rust={compare_result.rust_line_count})"
        return CompareArtifact("FAIL", compare_result.reason, compare_result.first_diff_line, summary)
    if compare_result.first_diff_line is not None:
        return CompareArtifact("FAIL", compare_result.reason, compare_result.first_diff_line, f"first difference at line {compare_result.first_diff_line}")
    return CompareArtifact("FAIL", compare_result.reason, compare_result.first_diff_line, compare_result.reason or "byte_mismatch")


def write_meta(
    meta_path: Path,
    *,
    shard_info: shard.Shard,
    java_path: Path,
    rust_path: Path,
    diff_path: Path,
    compare: comparator.CompareResult,
    artifact: CompareArtifact,
) -> None:
    payload = {
        "region": shard.shard_region_string(shard_info),
        "shard": shard_info.label,
        "java_file": str(java_path),
        "rust_file": str(rust_path),
        "diff_file": str(diff_path),
        "java_lines": compare.java_line_count,
        "rust_lines": compare.rust_line_count,
        "first_diff_line": artifact.first_diff_line,
        "summary": artifact.summary,
        "reason": artifact.reason,
    }
    meta_path.parent.mkdir(parents=True, exist_ok=True)
    meta_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def write_diff_preview(java_path: Path, rust_path: Path, diff_path: Path) -> None:
    if not java_path.exists() or not rust_path.exists():
        diff_path.write_text("missing output\n", encoding="utf-8")
        return
    with gzip.open(java_path, "rt", encoding="utf-8", errors="replace") as left_handle:
        java_lines = left_handle.readlines()
    with gzip.open(rust_path, "rt", encoding="utf-8", errors="replace") as right_handle:
        rust_lines = right_handle.readlines()
    diff_lines = list(islice(difflib.unified_diff(java_lines, rust_lines, fromfile=str(java_path), tofile=str(rust_path)), 500))
    diff_path.write_text("".join(diff_lines) if diff_lines else "byte_mismatch\n", encoding="utf-8")


def compare_shards(
    parity: runner.ParityRunner,
    shards: Sequence[shard.Shard],
    test_config: config.TestConfig,
    *,
    rust_mtime: str | None,
    java_runs: dict[str, runner.RunResult],
    rust_runs: dict[str, runner.RunResult],
    timeout_seconds: int,
) -> tuple[reporter.ChromResult, list[dict[str, object]]]:
    failures: list[reporter.ShardFailure] = []
    failure_rows: list[dict[str, object]] = []
    pass_count = 0
    fail_count = 0
    empty_count = 0

    print("--- Phase 3: Compare ---")
    for shard_info in shards:
        status_path = shard_status_path(test_config, shard_info)
        meta_path = shard_meta_path(test_config, shard_info)
        diff_path = shard_diff_path(test_config, shard_info)
        marker_path = shard_verified_path(test_config, shard_info)
        java_path = parity._output_path("java", shard_info, test_config)
        rust_path = parity._output_path("rust", shard_info, test_config)

        remove_if_exists(status_path)
        remove_if_exists(meta_path)
        remove_if_exists(diff_path)

        if is_marker_current(marker_path, rust_mtime):
            status_path.parent.mkdir(parents=True, exist_ok=True)
            status_path.write_text("PASS\n", encoding="utf-8")
            pass_count += 1
            continue

        compare_result = comparator.compare_shard_gz(java_path, rust_path)
        artifact = compare_summary(
            compare_result,
            java_run=java_runs.get(shard_info.label),
            rust_run=rust_runs.get(shard_info.label),
            timeout_seconds=timeout_seconds,
        )

        status_path.parent.mkdir(parents=True, exist_ok=True)
        status_path.write_text(f"{artifact.status}\n", encoding="utf-8")
        if artifact.status == "PASS":
            write_marker(marker_path, rust_mtime=rust_mtime, payload={"region": shard.shard_region_string(shard_info)})
            pass_count += 1
            continue
        if artifact.status == "EMPTY":
            write_marker(marker_path, rust_mtime=rust_mtime, payload={"region": shard.shard_region_string(shard_info), "empty": True})
            pass_count += 1
            empty_count += 1
            continue

        remove_if_exists(marker_path)
        fail_count += 1
        if artifact.reason != "missing_output":
            write_diff_preview(java_path, rust_path, diff_path)
        else:
            diff_path.write_text(f"{artifact.summary}\n", encoding="utf-8")
        write_meta(
            meta_path,
            shard_info=shard_info,
            java_path=java_path,
            rust_path=rust_path,
            diff_path=diff_path,
            compare=compare_result,
            artifact=artifact,
        )
        failures.append(
            reporter.ShardFailure(
                shard_label=shard_info.label,
                region=shard.shard_region_string(shard_info),
                reason=artifact.reason or "mismatch",
                first_diff_line=artifact.first_diff_line,
                summary=artifact.summary,
            )
        )
        failure_rows.append(
            {
                "shard": shard_info.label,
                "region": shard.shard_region_string(shard_info),
                "reason": artifact.reason,
                "first_diff_line": artifact.first_diff_line,
                "summary": artifact.summary,
            }
        )

    status = "FAIL" if fail_count else ("EMPTY" if empty_count == len(shards) else "PASS")
    chrom_result = reporter.ChromResult(
        status=status,
        total_shards=len(shards),
        pass_count=pass_count - empty_count,
        fail_count=fail_count,
        empty_count=empty_count,
        failures=failures,
    )
    return chrom_result, failure_rows


def write_results_json(
    test_config: config.TestConfig,
    chrom: str,
    *,
    elapsed_seconds: float,
    chrom_result: reporter.ChromResult,
    failure_rows: Sequence[dict[str, object]],
) -> None:
    output_path = results_json_path(test_config, chrom)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "chromosome": chrom,
        "config_id": test_config.config_id,
        "opts_label": test_config.label,
        "total_shards": chrom_result.total_shards,
        "pass": chrom_result.pass_count,
        "fail": chrom_result.fail_count,
        "empty": chrom_result.empty_count,
        "elapsed_seconds": elapsed_seconds,
        "status": chrom_result.status,
        "failures": list(failure_rows),
    }
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def write_chrom_verified(test_config: config.TestConfig, chrom: str, *, chrom_result: reporter.ChromResult, rust_mtime: str | None) -> None:
    marker_path = chrom_verified_path(test_config, chrom)
    write_marker(
        marker_path,
        rust_mtime=rust_mtime,
        payload={
            "chromosome": chrom,
            "total_shards": chrom_result.total_shards,
            "pass": chrom_result.pass_count,
            "empty": chrom_result.empty_count,
        },
    )


def cleanup_outputs(test_config: config.TestConfig, chrom: str, *, clean_all: bool, chrom_result: reporter.ChromResult) -> None:
    root = chrom_dir(test_config, chrom)
    java_dir = root / "java"
    rust_dir = root / "rust"
    diff_root = root / "diff"
    if chrom_result.fail_count == 0:
        shutil.rmtree(rust_dir, ignore_errors=True)
        shutil.rmtree(diff_root, ignore_errors=True)
        if clean_all:
            shutil.rmtree(java_dir, ignore_errors=True)
        return

    failing = {failure.shard_label for failure in chrom_result.failures}
    if not clean_all:
        for artifact in rust_dir.glob("*.tsv.gz"):
            label = artifact.stem.removeprefix("shard_").removesuffix(".tsv")
            if label not in failing:
                remove_if_exists(artifact)
        for artifact in rust_dir.glob("*.log"):
            label = artifact.stem.removeprefix("shard_")
            if label not in failing:
                remove_if_exists(artifact)
        return

    for mode in ("java", "rust"):
        mode_dir = root / mode
        for artifact in mode_dir.glob("*"):
            label = artifact.stem.removeprefix("shard_").removesuffix(".tsv")
            if label not in failing:
                remove_if_exists(artifact)


def run_retry_attempts(
    parity: runner.ParityRunner,
    shards: Sequence[shard.Shard],
    test_config: config.TestConfig,
    *,
    args: argparse.Namespace,
    rust_mtime: str | None,
    java_runs: dict[str, runner.RunResult],
    rust_runs: dict[str, runner.RunResult],
) -> tuple[reporter.ChromResult, list[dict[str, object]]]:
    chrom_result, failure_rows = compare_shards(
        parity,
        shards,
        test_config,
        rust_mtime=rust_mtime,
        java_runs=java_runs,
        rust_runs=rust_runs,
        timeout_seconds=args.timeout,
    )
    if args.retry <= 0:
        return chrom_result, failure_rows

    for attempt in range(1, args.retry + 1):
        retryable = [
            failure
            for failure in chrom_result.failures
            if failure.reason in {"java_empty_suspect", "rust_empty_suspect", "java_missing", "rust_missing", "missing_output"}
        ]
        if not retryable:
            break
        labels = {failure.shard_label for failure in retryable}
        print(f"--- Retry {attempt}/{args.retry}: {len(labels)} shards ---")
        for shard_info in shards:
            if shard_info.label not in labels:
                continue
            remove_if_exists(shard_verified_path(test_config, shard_info))
            remove_if_exists(shard_status_path(test_config, shard_info))
            remove_if_exists(shard_meta_path(test_config, shard_info))
            remove_if_exists(shard_diff_path(test_config, shard_info))
            java_path = parity._output_path("java", shard_info, test_config)
            rust_path = parity._output_path("rust", shard_info, test_config)
            remove_if_exists(rust_path)
            remove_if_exists(parity._log_path("rust", shard_info, test_config))
            if not args.rust_only:
                remove_if_exists(java_path)
                remove_if_exists(parity._log_path("java", shard_info, test_config))
                java_runs[shard_info.label] = parity.run_java_shard(shard_info, test_config)
            rust_runs[shard_info.label] = parity.run_rust_shard(shard_info, test_config)

        chrom_result, failure_rows = compare_shards(
            parity,
            shards,
            test_config,
            rust_mtime=rust_mtime,
            java_runs=java_runs,
            rust_runs=rust_runs,
            timeout_seconds=args.timeout,
        )
        if chrom_result.fail_count == 0:
            break
    return chrom_result, failure_rows


def run_chromosome(
    parity: runner.ParityRunner,
    test_config: config.TestConfig,
    chrom: str,
    *,
    chrom_len: int,
    args: argparse.Namespace,
    rust_mtime: str | None,
) -> ChromosomeRunOutcome:
    shards = shard.generate_shards(chrom, chrom_len, args.shard_size)
    if not shards:
        result = reporter.ChromResult(status="EMPTY", total_shards=0, pass_count=0, fail_count=0, empty_count=0)
        write_results_json(test_config, chrom, elapsed_seconds=0.0, chrom_result=result, failure_rows=[])
        return ChromosomeRunOutcome(chrom=chrom, result=result, elapsed_seconds=0.0)

    verified_path = chrom_verified_path(test_config, chrom)
    if is_marker_current(verified_path, rust_mtime):
        result = reporter.ChromResult(
            status="PASS",
            total_shards=len(shards),
            pass_count=len(shards),
            fail_count=0,
            empty_count=0,
        )
        write_results_json(test_config, chrom, elapsed_seconds=0.0, chrom_result=result, failure_rows=[])
        print(f"=== {chrom}: already verified ({len(shards)} shards) ===")
        return ChromosomeRunOutcome(chrom=chrom, result=result, elapsed_seconds=0.0)

    prepare_dirs(test_config, chrom, clean_rust=args.clean_rust)
    java_runs: dict[str, runner.RunResult] = {}
    rust_runs: dict[str, runner.RunResult] = {}
    java_needed = phase_indices_for_java(parity, shards, test_config, rust_mtime)
    rust_needed = phase_indices_for_rust(parity, shards, test_config, rust_mtime)

    print(f"=== Processing {chrom} ({len(shards)} shards, rust workers={args.parallel}) ===")
    import time

    start_time = time.monotonic()
    if args.rust_only:
        validate_java_cache(parity, shards, test_config)
        print("--- Phase 1: Java cache validated (--rust-only) ---")
    elif java_needed:
        java_runs = run_parallel_phase(
            "Phase 1 (Java)",
            java_needed,
            max_workers=args.java_parallel,
            args=args,
            executor_fn=lambda shard_info: parity.run_java_shard(shard_info, test_config),
        )
    else:
        print("--- Phase 1: Java (all cached) ---")

    if rust_needed:
        rust_runs = run_parallel_phase(
            "Phase 2 (Rust)",
            rust_needed,
            max_workers=args.parallel,
            args=args,
            executor_fn=lambda shard_info: parity.run_rust_shard(shard_info, test_config),
        )
    else:
        print("--- Phase 2: Rust (all cached) ---")

    chrom_result, failure_rows = run_retry_attempts(
        parity,
        shards,
        test_config,
        args=args,
        rust_mtime=rust_mtime,
        java_runs=java_runs,
        rust_runs=rust_runs,
    )
    elapsed_seconds = time.monotonic() - start_time
    write_results_json(test_config, chrom, elapsed_seconds=elapsed_seconds, chrom_result=chrom_result, failure_rows=failure_rows)
    if chrom_result.fail_count == 0:
        write_chrom_verified(test_config, chrom, chrom_result=chrom_result, rust_mtime=rust_mtime)
    else:
        remove_if_exists(verified_path)

    if not args.no_cleanup:
        print("--- Phase 4: Cleanup ---")
        cleanup_outputs(test_config, chrom, clean_all=args.clean_all, chrom_result=chrom_result)

    print_chrom_summary(chrom, chrom_result, elapsed_seconds)
    return ChromosomeRunOutcome(chrom=chrom, result=chrom_result, elapsed_seconds=elapsed_seconds)


def print_chrom_summary(chrom: str, chrom_result: reporter.ChromResult, elapsed_seconds: float) -> None:
    print()
    print(f"=== {SAMPLE_NAME} Parity Summary ({chrom}) ===")
    print(
        "Shards: "
        f"{chrom_result.total_shards} (pass: {chrom_result.pass_count}, fail: {chrom_result.fail_count}, empty: {chrom_result.empty_count})"
    )
    print(f"Time: {elapsed_seconds:.1f}s")
    if chrom_result.failures:
        print()
        print("Failing shards:")
        for failure in chrom_result.failures:
            print(f"  shard_{failure.shard_label} ({failure.region}): {failure.summary or failure.reason}")


def print_global_summary(outcomes: Sequence[ChromosomeRunOutcome]) -> None:
    total_chroms = len(outcomes)
    chrom_pass = sum(1 for outcome in outcomes if outcome.result.fail_count == 0)
    chrom_fail = total_chroms - chrom_pass
    total_shards = sum(outcome.result.total_shards for outcome in outcomes)
    total_pass = sum(outcome.result.pass_count for outcome in outcomes)
    total_fail = sum(outcome.result.fail_count for outcome in outcomes)
    total_empty = sum(outcome.result.empty_count for outcome in outcomes)

    print()
    print("=== ALL CHROMOSOMES SUMMARY ===")
    print(f"Chromosomes: {total_chroms} (pass: {chrom_pass}, fail: {chrom_fail})")
    print(f"Shards: {total_shards} (pass: {total_pass}, fail: {total_fail}, empty: {total_empty})")
    if chrom_fail:
        print()
        print("Failing chromosomes:")
        for outcome in outcomes:
            if outcome.result.fail_count:
                print(f"  {outcome.chrom}: {outcome.result.fail_count}/{outcome.result.total_shards} failed")


def create_runner_config(args: argparse.Namespace, rust_bin: Path) -> runner.RunnerConfig:
    return runner.RunnerConfig(
        ref_fasta=REF_FASTA,
        bam_path=BAM_PATH,
        java_jar=JAVA_JAR,
        rust_bin=rust_bin,
        base_dir=BASE_DIR,
        parallel=args.parallel,
        java_parallel=args.java_parallel,
        java_heap=args.java_heap,
        timeout=args.timeout,
        retry=0,
        freq=args.freq,
        compress=True,
        rust_only=args.rust_only,
        sample_name=SAMPLE_NAME,
    )


def run(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        configure_runtime_paths(
            bam_path=args.bam,
            ref_fasta=args.ref,
            sample_name=args.sample,
            base_dir=args.base_dir,
        )
        validate_args(args)
        ensure_inputs_exist()
        rust_bin = resolve_rust_bin(args)
        auto_build_if_needed(rust_bin, no_build=args.no_build)
        resolved_config = resolve_test_config(args)
        lengths = load_chrom_lengths()
        chromosomes = list_chromosomes(args, lengths)
        rust_mtime = rust_binary_mtime(rust_bin)
        parity = runner.ParityRunner(PROJECT_ROOT, create_runner_config(args, rust_bin))
    except Exception as exc:  # pragma: no cover - CLI error handling
        return fail(str(exc))

    outcomes: list[ChromosomeRunOutcome] = []
    for chrom in chromosomes:
        try:
            outcome = run_chromosome(
                parity,
                resolved_config,
                chrom,
                chrom_len=chrom_length(chrom, args, lengths),
                args=args,
                rust_mtime=rust_mtime,
            )
            outcomes.append(outcome)
        except Exception as exc:  # pragma: no cover - CLI error handling
            print(f"ERROR: {chrom}: {exc}", file=sys.stderr)
            if not args.no_stop:
                return 1

    if args.all_chr:
        print_global_summary(outcomes)

    return 1 if any(outcome.result.fail_count for outcome in outcomes) else 0


def main(argv: Sequence[str] | None = None) -> int:
    return run(argv)


if __name__ == "__main__":
    raise SystemExit(main())