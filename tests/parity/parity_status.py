#!/usr/bin/env python3
"""Status matrix display for parity cache state."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tests.parity import parity_runner
from tests.parity.lib import config, reporter


PROJECT_ROOT = Path(__file__).resolve().parents[2]
BASE_DIR = PROJECT_ROOT / parity_runner.BASE_DIR
DEFAULT_CHROMOSOMES = ("20", "22", "MT")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Display parity status for config x chromosome cells.")
    parser.add_argument("--preset", choices=tuple(config.PRESETS), help="Filter configs/chromosomes by preset")
    parser.add_argument("--config-id", action="append", help="Limit output to specific config IDs; repeatable")
    parser.add_argument("--chr", dest="chromosomes", action="append", help="Limit output to specific chromosomes; repeatable")
    parser.add_argument("--tier", action="append", help="Limit output to specific tiers; repeatable")
    parser.add_argument("--json", dest="json_path", type=Path, help="Write the status payload to a JSON file")
    return parser


def fail(message: str) -> int:
    print(f"ERROR: {message}", file=sys.stderr)
    return 1


def normalize_chr(chromosome: str) -> str:
    return chromosome[3:] if chromosome.startswith("chr") else chromosome


def load_real_chromosomes() -> list[str]:
    fai_path = PROJECT_ROOT / parity_runner.REF_FASTA_FAI
    if not fai_path.exists():
        return list(DEFAULT_CHROMOSOMES)
    chrom_lengths = parity_runner.shard.load_fai(fai_path)
    allowed = {str(number) for number in range(1, 23)} | {"X", "Y", "MT"}
    return [chrom for chrom in chrom_lengths if chrom in allowed]


def resolve_selection(args: argparse.Namespace) -> tuple[list[config.TestConfig], list[str]]:
    if args.preset and args.tier:
        raise ValueError("--preset cannot be combined with --tier")

    real_chromosomes = load_real_chromosomes()
    if args.preset:
        preset = config.PRESETS[args.preset]
        if preset.stages:
            seen_config_ids: set[str] = set()
            seen_chromosomes: set[str] = set()
            configs = []
            chromosomes = []
            for stage_configs, stage_chromosomes in config.resolve_preset_sequence(args.preset, real_chromosomes=real_chromosomes):
                for entry in stage_configs:
                    if entry.config_id in seen_config_ids:
                        continue
                    seen_config_ids.add(entry.config_id)
                    configs.append(entry)
                for chrom in stage_chromosomes:
                    if chrom in seen_chromosomes:
                        continue
                    seen_chromosomes.add(chrom)
                    chromosomes.append(chrom)
        else:
            configs, chromosomes = config.resolve_preset(args.preset, real_chromosomes=real_chromosomes)
    else:
        configs = list(config.TEST_MATRIX)
        chromosomes = list(DEFAULT_CHROMOSOMES)

    if args.tier:
        allowed_tiers = set(args.tier)
        configs = [entry for entry in configs if str(entry.tier) in allowed_tiers]
    if args.config_id:
        requested = set(args.config_id)
        configs = [entry for entry in configs if entry.config_id in requested]
    if args.chromosomes:
        chromosomes = [normalize_chr(chrom) for chrom in args.chromosomes]

    if not configs:
        raise ValueError("No configs matched the selected filters")
    if not chromosomes:
        raise ValueError("No chromosomes matched the selected filters")
    return configs, chromosomes


def status_for_cell(test_config: config.TestConfig, chrom: str, *, rust_mtime: str | None) -> str:
    verified_path = parity_runner.chrom_verified_path(test_config, chrom)
    if parity_runner.is_marker_current(verified_path, rust_mtime):
        return "PASS"
    if verified_path.exists():
        return "STALE"

    chrom_root = parity_runner.chrom_dir(test_config, chrom)
    diff_dir = chrom_root / "diff"
    if diff_dir.is_dir():
        for status_path in sorted(diff_dir.glob("*.status")):
            if status_path.read_text(encoding="utf-8").strip() == "FAIL":
                return "FAIL"

    results_path = parity_runner.results_json_path(test_config, chrom)
    if results_path.exists():
        payload = json.loads(results_path.read_text(encoding="utf-8"))
        if int(payload.get("fail", 0)) > 0:
            return "FAIL"

    return "PEND"


def selected_has_any_data(configs: list[config.TestConfig], chromosomes: list[str]) -> bool:
    for entry in configs:
        for chrom in chromosomes:
            if parity_runner.chrom_dir(entry, chrom).exists():
                return True
    return False


def build_payload(configs: list[config.TestConfig], chromosomes: list[str], *, rust_mtime: str | None) -> dict[str, object]:
    rows: list[dict[str, object]] = []
    matrix: dict[str, dict[str, str]] = {}
    for entry in configs:
        chrom_statuses: dict[str, str] = {}
        for chrom in chromosomes:
            chrom_statuses[chrom] = status_for_cell(entry, chrom, rust_mtime=rust_mtime)
        matrix[entry.config_id] = chrom_statuses
        rows.append(
            {
                "config_id": entry.config_id,
                "label": entry.label,
                "tier": entry.tier,
                "chromosomes": chrom_statuses,
            }
        )
    return {
        "base_dir": str(parity_runner.BASE_DIR),
        "rust_binary": str(parity_runner.DEFAULT_RUST_BIN),
        "rust_binary_mtime": rust_mtime,
        "chromosomes": chromosomes,
        "rows": rows,
        "matrix": matrix,
    }


def write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        configs, chromosomes = resolve_selection(args)
    except ValueError as exc:
        return fail(str(exc))

    rust_mtime = parity_runner.rust_binary_mtime(parity_runner.DEFAULT_RUST_BIN)
    payload = build_payload(configs, chromosomes, rust_mtime=rust_mtime)

    if args.json_path is not None:
        write_json(args.json_path, payload)

    if not BASE_DIR.exists() or not selected_has_any_data(configs, chromosomes):
        print("No parity data found")
        return 0

    reporter.print_status_matrix(payload["matrix"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())