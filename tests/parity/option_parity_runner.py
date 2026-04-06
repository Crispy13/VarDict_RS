#!/usr/bin/env python3
"""Config sweep orchestrator for the Python parity harness."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tests.parity import parity_runner
from tests.parity.lib import checkpoint, config, reporter, shard


PROJECT_ROOT = Path(__file__).resolve().parents[2]
BASE_DIR = parity_runner.DEFAULT_BASE_DIR
REPORT_FILE = BASE_DIR / "option_parity_report.json"
REF_FASTA_FAI = parity_runner.reference_fai_path(parity_runner.DEFAULT_REF_FASTA)


def checkpoint_path(base_dir: Path | None = None) -> Path:
    resolved_base_dir = BASE_DIR if base_dir is None else base_dir
    return resolved_base_dir / ".sweep_checkpoint.json"


def report_file_path(base_dir: Path | None = None) -> Path:
    resolved_base_dir = BASE_DIR if base_dir is None else base_dir
    return resolved_base_dir / "option_parity_report.json"


def configure_runtime_paths(*, bam_path: Path, ref_fasta: Path, sample_name: str, base_dir: Path) -> None:
    global BASE_DIR, REPORT_FILE, REF_FASTA_FAI

    BASE_DIR = base_dir
    REPORT_FILE = report_file_path(base_dir)
    REF_FASTA_FAI = parity_runner.reference_fai_path(ref_fasta)
    parity_runner.configure_runtime_paths(
        bam_path=bam_path,
        ref_fasta=ref_fasta,
        sample_name=sample_name,
        base_dir=base_dir,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run parity across a config x chromosome matrix.")
    parser.add_argument("--bam", type=Path, default=parity_runner.DEFAULT_BAM_PATH, help=f"Path to BAM file (default: {parity_runner.DEFAULT_BAM_PATH})")
    parser.add_argument("--ref", type=Path, default=parity_runner.DEFAULT_REF_FASTA, help=f"Path to reference FASTA (default: {parity_runner.DEFAULT_REF_FASTA})")
    parser.add_argument("--sample", default=parity_runner.DEFAULT_SAMPLE_NAME, help=f"Sample name for -N (default: {parity_runner.DEFAULT_SAMPLE_NAME})")
    parser.add_argument("--base-dir", type=Path, default=parity_runner.DEFAULT_BASE_DIR, help=f"Output directory (default: {parity_runner.DEFAULT_BASE_DIR})")
    parser.add_argument(
        "--preset",
        choices=tuple(config.PRESETS),
        help="Preset selection: smoke, dev, tier1, config-spread, core-wide, pairwise, full-gate, release",
    )
    parser.add_argument("--tier", action="append", help="Filter by tier; repeatable")
    parser.add_argument("--config-id", action="append", help="Run only specific config IDs; repeatable")
    parser.add_argument("--chr", dest="chromosomes", action="append", help="Override chromosomes; repeatable")
    parser.add_argument("--all-chr", action="store_true", help="Use all chromosomes from the reference FAI")
    parser.add_argument("--no-stop", action="store_true", help="Continue after cell failures")
    parser.add_argument("--no-cleanup", action="store_true", help="Keep shard artifacts after each cell")
    parser.add_argument("--rust-only", action="store_true", help="Skip Java generation and require cached Java shards")
    parser.add_argument("--clean-all", action="store_true", help="Delete Java cache during cleanup as well")
    parser.add_argument("--release", action="store_true", help="Use the repository debug-release build path")
    parser.add_argument("--rust-bin", help="Use a specific Rust binary path")
    parser.add_argument("--no-build", action="store_true", help="Skip cargo build in both orchestrator and inner runner")
    parser.add_argument("--parallel", type=int, help="Rust worker count override")
    parser.add_argument("--java-parallel", type=int, help="Java worker count override")
    parser.add_argument("--java-heap", help="Java heap override")
    parser.add_argument("--timeout", type=int, help="Per-shard timeout override")
    parser.add_argument("--retry", type=int, default=0, help="Per-shard retry count")
    parser.add_argument("--freq", type=float, help="Default frequency when the config flags omit -f")
    parser.add_argument("--shard-size", type=int, help="Shard size override")
    parser.add_argument("--mem-warn-gb", type=float, help="Memory warning threshold override")
    parser.add_argument("--mem-abort-gb", type=float, help="Memory abort threshold override")
    parser.add_argument("--disk-warn-gb", type=float, help="Disk warning threshold override")
    parser.add_argument("--dry-run", action="store_true", help="Print the matrix that would run")
    parser.add_argument("--status", action="store_true", help="Print the current parity status matrix")
    parser.add_argument("--fresh", action="store_true", help="Delete the checkpoint before running")
    return parser


def fail(message: str) -> int:
    print(f"ERROR: {message}", file=sys.stderr)
    return 1


def load_real_chromosomes() -> list[str]:
    fai = PROJECT_ROOT / REF_FASTA_FAI
    if not fai.exists():
        raise FileNotFoundError(f"Reference FAI not found: {REF_FASTA_FAI}")
    chrom_lengths = shard.load_fai(fai)
    return [chrom for chrom in chrom_lengths if chrom in {str(number) for number in range(1, 23)} | {"X", "Y", "MT"}]


def resolve_matrix(args: argparse.Namespace) -> tuple[list[config.TestConfig], list[str], list[str]]:
    if args.preset and args.tier:
        raise ValueError("--preset cannot be combined with --tier")

    real_chromosomes = load_real_chromosomes()
    preset_name = args.preset or "config-spread"

    if preset_name == "full-gate":
        return [], [], list(config.PRESETS["full-gate"].stages)

    configs, chromosomes = config.resolve_preset(
        preset_name,
        real_chromosomes=real_chromosomes,
    )

    if args.tier:
        allowed_tiers = set(args.tier)
        configs = [entry for entry in config.TEST_MATRIX if str(entry.tier) in allowed_tiers]
    if args.config_id:
        requested = set(args.config_id)
        configs = [entry for entry in configs if entry.config_id in requested]
    if args.all_chr:
        chromosomes = real_chromosomes
    elif args.chromosomes:
        chromosomes = list(args.chromosomes)

    if not configs:
        raise ValueError("No configs matched the selected filters")
    if not chromosomes:
        raise ValueError("No chromosomes matched the selected filters")
    return configs, chromosomes, []


def current_rust_mtime(args: argparse.Namespace) -> str | None:
    rust_bin = parity_runner.resolve_rust_bin(args)
    return parity_runner.rust_binary_mtime(rust_bin)


def cell_status(test_config: config.TestConfig, chrom: str, *, rust_mtime: str | None) -> str:
    verified_path = parity_runner.chrom_verified_path(test_config, chrom)
    if parity_runner.is_marker_current(verified_path, rust_mtime):
        return "PASS"
    if verified_path.exists():
        return "STALE"

    results_path = parity_runner.results_json_path(test_config, chrom)
    if results_path.exists():
        payload = json.loads(results_path.read_text(encoding="utf-8"))
        if int(payload.get("fail", 0)) > 0:
            return "FAIL"
        if int(payload.get("pass", 0)) > 0 or int(payload.get("empty", 0)) > 0:
            return "PEND"

    chrom_root = parity_runner.chrom_dir(test_config, chrom)
    if (chrom_root / "java").exists() or (chrom_root / "rust").exists() or (chrom_root / "diff").exists():
        return "PEND"
    return "PEND"


def dry_run_status_text(status: str) -> str:
    return {
        "PASS": "verified",
        "STALE": "stale",
        "FAIL": "failed",
        "PEND": "unverified",
    }.get(status, status.lower())


def print_dry_run(configs: Sequence[config.TestConfig], chromosomes: Sequence[str], *, rust_mtime: str | None) -> None:
    cell_count = len(configs) * len(chromosomes)
    print(f"DRY-RUN: {len(configs)} configs x {len(chromosomes)} chromosomes = {cell_count} cells")
    chrom_lengths = shard.load_fai(PROJECT_ROOT / REF_FASTA_FAI)

    for entry in configs:
        joined_chroms = " ".join(chromosomes)
        print(f'DRY-RUN: {entry.config_id} ({entry.label}) tier={entry.tier} flags="{entry.flags}" chrs="{joined_chroms}"')
        for chrom in chromosomes:
            count = len(shard.generate_shards(chrom, chrom_lengths[chrom]))
            status = dry_run_status_text(cell_status(entry, chrom, rust_mtime=rust_mtime))
            display_chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
            print(f"  {display_chrom}: {count} shards, status={status}")


def print_status(configs: Sequence[config.TestConfig], chromosomes: Sequence[str], *, rust_mtime: str | None) -> None:
    matrix: dict[str, dict[str, str]] = {}
    for entry in configs:
        matrix[entry.config_id] = {}
        for chrom in chromosomes:
            matrix[entry.config_id][chrom] = cell_status(entry, chrom, rust_mtime=rust_mtime)
    reporter.print_status_matrix(matrix)


def passthrough_args(args: argparse.Namespace) -> list[str]:
    forwarded: list[str] = [
        "--bam",
        str(args.bam),
        "--ref",
        str(args.ref),
        "--sample",
        args.sample,
        "--base-dir",
        str(args.base_dir),
    ]
    if args.no_stop:
        forwarded.append("--no-stop")
    if args.no_cleanup:
        forwarded.append("--no-cleanup")
    if args.rust_only:
        forwarded.append("--rust-only")
    if args.clean_all:
        forwarded.append("--clean-all")
    if args.release:
        forwarded.append("--release")
    if args.no_build:
        forwarded.append("--no-build")
    if args.rust_bin:
        forwarded.extend(["--rust-bin", args.rust_bin])
    if args.parallel is not None:
        forwarded.extend(["--parallel", str(args.parallel)])
    if args.java_parallel is not None:
        forwarded.extend(["--java-parallel", str(args.java_parallel)])
    if args.java_heap is not None:
        forwarded.extend(["--java-heap", str(args.java_heap)])
    if args.timeout is not None:
        forwarded.extend(["--timeout", str(args.timeout)])
    if args.retry:
        forwarded.extend(["--retry", str(args.retry)])
    if args.freq is not None:
        forwarded.extend(["--freq", str(args.freq)])
    if args.shard_size is not None:
        forwarded.extend(["--shard-size", str(args.shard_size)])
    if args.mem_warn_gb is not None:
        forwarded.extend(["--mem-warn-gb", str(args.mem_warn_gb)])
    if args.mem_abort_gb is not None:
        forwarded.extend(["--mem-abort-gb", str(args.mem_abort_gb)])
    if args.disk_warn_gb is not None:
        forwarded.extend(["--disk-warn-gb", str(args.disk_warn_gb)])
    return forwarded


def apply_resource_overrides(entry: config.TestConfig, args: argparse.Namespace, forwarded: list[str]) -> list[str]:
    overrides = config.resource_overrides(entry)
    option_specs = (
        ("parallel", args.parallel, "--parallel"),
        ("java_parallel", args.java_parallel, "--java-parallel"),
        ("java_heap", args.java_heap, "--java-heap"),
    )
    for override_key, explicit_value, option_flag in option_specs:
        if explicit_value is not None:
            continue
        override_value = overrides.get(override_key)
        if override_value is not None:
            forwarded.extend([option_flag, str(override_value)])
    return forwarded


def load_cell_result(entry: config.TestConfig, chrom: str, *, elapsed_seconds: float, status_override: str | None = None) -> reporter.ChromResult:
    results_path = parity_runner.results_json_path(entry, chrom)
    if not results_path.exists():
        status = status_override or "FAIL"
        return reporter.ChromResult(status=status, total_shards=0, pass_count=0, fail_count=1 if status == "FAIL" else 0, empty_count=0)

    payload = json.loads(results_path.read_text(encoding="utf-8"))
    failures = [
        reporter.ShardFailure(
            shard_label=str(item.get("shard", "")),
            region=str(item.get("region", "")),
            reason=str(item.get("reason", "unknown")),
            first_diff_line=item.get("first_diff_line"),
            summary=str(item.get("summary", "")),
        )
        for item in payload.get("failures", [])
    ]
    return reporter.ChromResult(
        status=status_override or str(payload.get("status", "PASS")),
        total_shards=int(payload.get("total_shards", 0)),
        pass_count=int(payload.get("pass", 0)),
        fail_count=int(payload.get("fail", 0)),
        empty_count=int(payload.get("empty", 0)),
        failures=failures,
    )


def run_cell(entry: config.TestConfig, chrom: str, args: argparse.Namespace) -> tuple[int, reporter.ChromResult, float]:
    forwarded = apply_resource_overrides(entry, args, passthrough_args(args))
    if "--no-build" not in forwarded:
        forwarded.append("--no-build")
    argv = [
        "--chr",
        chrom,
        "--config",
        entry.config_id,
        *forwarded,
    ]
    import time

    start = time.monotonic()
    rc = parity_runner.main(argv)
    elapsed = time.monotonic() - start
    result = load_cell_result(entry, chrom, elapsed_seconds=elapsed, status_override="FAIL" if rc else None)
    return rc, result, elapsed


def build_once_if_needed(args: argparse.Namespace) -> None:
    rust_bin = parity_runner.resolve_rust_bin(args)
    parity_runner.auto_build_if_needed(rust_bin, no_build=args.no_build)


def run_stage(configs: Sequence[config.TestConfig], chromosomes: Sequence[str], args: argparse.Namespace) -> int:
    if args.fresh:
        checkpoint.clear_checkpoint(checkpoint_path(args.base_dir))

    state = checkpoint.load_checkpoint(checkpoint_path(args.base_dir)) or checkpoint.SweepCheckpoint()
    completed_cells = {tuple(item) for item in state.completed}
    config_results: dict[str, reporter.ConfigResult] = {
        entry.config_id: reporter.ConfigResult(
            config_id=entry.config_id,
            label=entry.label,
            tier=entry.tier,
            chromosomes={},
            elapsed_seconds=0.0,
        )
        for entry in configs
    }

    failed = False
    for entry in configs:
        for chrom in chromosomes:
            cell_key = (entry.config_id, chrom)
            if cell_key in completed_cells:
                status = cell_status(entry, chrom, rust_mtime=current_rust_mtime(args))
                config_results[entry.config_id].chromosomes[chrom] = load_cell_result(entry, chrom, elapsed_seconds=0.0, status_override=status)
                continue

            print()
            print(f"--- {entry.config_id}: chromosome {chrom} ---")
            rc, chrom_result, elapsed = run_cell(entry, chrom, args)
            config_results[entry.config_id].chromosomes[chrom] = chrom_result
            config_results[entry.config_id].elapsed_seconds += elapsed

            if rc == 0:
                state.failed = [item for item in state.failed if tuple(item[:2]) != cell_key]
                state.completed.append(cell_key)
                checkpoint.save_checkpoint(state, checkpoint_path(args.base_dir))
            else:
                state.failed = [item for item in state.failed if tuple(item[:2]) != cell_key]
                state.failed.append((entry.config_id, chrom, chrom_result.status))
                checkpoint.save_checkpoint(state, checkpoint_path(args.base_dir))
                failed = True
                if not args.no_stop:
                    results = [config_results[entry_id] for entry_id in config_results if config_results[entry_id].chromosomes]
                    reporter.write_json_report(results, report_file_path(args.base_dir))
                    reporter.print_summary_table(results)
                    return 1

    results = list(config_results.values())
    reporter.write_json_report(results, report_file_path(args.base_dir))
    reporter.print_summary_table(results)
    if not failed:
        checkpoint.clear_checkpoint(checkpoint_path(args.base_dir))
    return 1 if failed else 0


def run_full_gate(args: argparse.Namespace) -> int:
    for stage_name in config.PRESETS["full-gate"].stages:
        print()
        print("========================================")
        print(f"FULL-GATE: Running {stage_name}")
        print("========================================")
        stage_args = argparse.Namespace(**vars(args))
        stage_args.preset = stage_name
        stage_args.fresh = True
        rc = run_from_args(stage_args)
        if rc != 0:
            print(f"FULL-GATE: {stage_name} FAILED")
            return rc
        print(f"FULL-GATE: {stage_name} PASSED")
    print()
    print("FULL-GATE: All stages passed.")
    return 0


def run_from_args(args: argparse.Namespace) -> int:
    if args.dry_run and args.status:
        return fail("--dry-run and --status cannot be used together")

    try:
        configs, chromosomes, stages = resolve_matrix(args)
    except Exception as exc:  # pragma: no cover - CLI error handling
        return fail(str(exc))

    if stages:
        return run_full_gate(args)

    rust_mtime = current_rust_mtime(args)
    if args.status:
        print_status(configs, chromosomes, rust_mtime=rust_mtime)
        return 0
    if args.dry_run:
        print_dry_run(configs, chromosomes, rust_mtime=rust_mtime)
        return 0

    try:
        build_once_if_needed(args)
    except Exception as exc:  # pragma: no cover - CLI error handling
        return fail(str(exc))

    return run_stage(configs, chromosomes, args)


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    configure_runtime_paths(
        bam_path=args.bam,
        ref_fasta=args.ref,
        sample_name=args.sample,
        base_dir=args.base_dir,
    )
    return run_from_args(args)


if __name__ == "__main__":
    raise SystemExit(main())