#!/usr/bin/env python3
"""Cleanup utility for parity cache artifacts."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import sys

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tests.parity import parity_runner
from tests.parity.lib import config


PROJECT_ROOT = Path(__file__).resolve().parents[2]
BASE_DIR = PROJECT_ROOT / parity_runner.BASE_DIR
JAVA_JAR = PROJECT_ROOT / parity_runner.JAVA_JAR
RUST_BIN = PROJECT_ROOT / parity_runner.DEFAULT_RUST_BIN
DIFF_SUBDIR = "diff"
JAVA_SUBDIR = "java"
RUST_SUBDIR = "rust"


@dataclass(frozen=True, slots=True)
class CleanupTarget:
    path: Path
    group: str
    size_bytes: int


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Remove stale or scoped parity cache artifacts.")
    parser.add_argument(
        "--scope",
        choices=("stale", "all", "java", "rust", "failed", "verified"),
        required=True,
        help="Cleanup scope: stale, all, java, rust, failed, or verified",
    )
    parser.add_argument("--config-id", action="append", help="Limit cleanup to specific config IDs; repeatable")
    parser.add_argument("--chr", dest="chromosomes", action="append", help="Limit cleanup to specific chromosomes; repeatable")
    parser.add_argument("--dry-run", action="store_true", help="List what would be deleted without deleting")
    return parser


def fail(message: str) -> int:
    print(f"ERROR: {message}", file=sys.stderr)
    return 1


def normalize_chr(chromosome: str) -> str:
    return chromosome[3:] if chromosome.startswith("chr") else chromosome


def resolve_config_labels(config_ids: list[str] | None) -> set[str] | None:
    if not config_ids:
        return None
    labels: set[str] = set()
    missing: list[str] = []
    for config_id in config_ids:
        entry = config.TEST_MATRIX_BY_ID.get(config_id)
        if entry is None:
            missing.append(config_id)
            continue
        labels.add(entry.label)
    if missing:
        raise ValueError(f"Unknown config ID(s): {', '.join(sorted(missing))}")
    return labels


def iter_chr_dirs(base_dir: Path) -> list[Path]:
    if not base_dir.exists():
        return []
    return sorted(
        path
        for config_dir in base_dir.iterdir()
        if config_dir.is_dir()
        for path in config_dir.iterdir()
        if path.is_dir()
    )


def matches_scope(path: Path, *, labels: set[str] | None, chromosomes: set[str] | None) -> bool:
    label = path.parent.name
    chrom = normalize_chr(path.name)
    if labels is not None and label not in labels:
        return False
    if chromosomes is not None and chrom not in chromosomes:
        return False
    return True


def path_size(path: Path) -> int:
    try:
        return path.lstat().st_size
    except OSError:
        return 0


def human_bytes(size_bytes: int) -> str:
    units = ("B", "K", "M", "G", "T", "P")
    value = float(size_bytes)
    unit_index = 0
    while value >= 1024.0 and unit_index < len(units) - 1:
        value /= 1024.0
        unit_index += 1
    if unit_index == 0:
        return f"{int(value)}{units[unit_index]}"
    if unit_index == 1 and value >= 10:
        return f"{value:.0f}{units[unit_index]}"
    return f"{value:.1f}{units[unit_index]}"


def verified_marker_current(marker_path: Path, rust_mtime: str | None) -> bool:
    return parity_runner.is_marker_current(marker_path, rust_mtime)


def iter_shard_ids(*directories: Path) -> list[str]:
    shard_ids: set[str] = set()
    for directory in directories:
        if not directory.is_dir():
            continue
        for path in directory.iterdir():
            if not path.is_file() and not path.is_symlink():
                continue
            name = path.name
            if not name.startswith("shard_"):
                continue
            remainder = name[len("shard_") :]
            shard_id = remainder.split(".", 1)[0]
            if shard_id:
                shard_ids.add(shard_id)
    return sorted(shard_ids)


def shard_entries(directory: Path, shard_id: str) -> list[Path]:
    if not directory.is_dir():
        return []
    prefix = f"shard_{shard_id}."
    return sorted(
        path
        for path in directory.iterdir()
        if (path.is_file() or path.is_symlink()) and path.name.startswith(prefix)
    )


def shard_output_path(directory: Path, shard_id: str) -> Path | None:
    for suffix in (".tsv.gz", ".tsv"):
        candidate = directory / f"shard_{shard_id}{suffix}"
        if candidate.exists():
            return candidate
    return None


def tree_entries(directory: Path) -> list[Path]:
    if not directory.is_dir():
        return []
    return sorted(path for path in directory.rglob("*") if path.is_file() or path.is_symlink())


def is_within_base(path: Path, base_dir: Path) -> bool:
    try:
        path.resolve(strict=False).relative_to(base_dir.resolve(strict=False))
    except ValueError:
        return False
    return True


def add_target(targets: dict[Path, CleanupTarget], path: Path, *, base_dir: Path) -> None:
    if not path.exists() and not path.is_symlink():
        return
    if not is_within_base(path, base_dir):
        raise ValueError(f"Refusing to delete outside {base_dir}: {path}")
    relative = path.relative_to(base_dir)
    group_path = relative.parent
    group = "./" if group_path == Path(".") else f"{group_path.as_posix()}/"
    targets.setdefault(path, CleanupTarget(path=path, group=group, size_bytes=path_size(path)))


def collect_stale_java(base_dir: Path, *, labels: set[str] | None, chromosomes: set[str] | None) -> dict[Path, CleanupTarget]:
    if not JAVA_JAR.exists():
        raise FileNotFoundError(f"Java JAR not found: {JAVA_JAR}")
    jar_mtime = JAVA_JAR.stat().st_mtime
    targets: dict[Path, CleanupTarget] = {}
    for chr_dir in iter_chr_dirs(base_dir):
        if not matches_scope(chr_dir, labels=labels, chromosomes=chromosomes):
            continue
        java_dir = chr_dir / JAVA_SUBDIR
        for shard_id in iter_shard_ids(java_dir):
            output_path = shard_output_path(java_dir, shard_id)
            if output_path is None:
                for candidate in shard_entries(java_dir, shard_id):
                    add_target(targets, candidate, base_dir=base_dir)
                continue
            if output_path.stat().st_mtime < jar_mtime:
                for candidate in shard_entries(java_dir, shard_id):
                    add_target(targets, candidate, base_dir=base_dir)
    return targets


def collect_stale_rust(base_dir: Path, *, labels: set[str] | None, chromosomes: set[str] | None) -> dict[Path, CleanupTarget]:
    if not RUST_BIN.exists():
        raise FileNotFoundError(f"Rust binary not found: {RUST_BIN}")
    rust_mtime = RUST_BIN.stat().st_mtime
    targets: dict[Path, CleanupTarget] = {}
    for chr_dir in iter_chr_dirs(base_dir):
        if not matches_scope(chr_dir, labels=labels, chromosomes=chromosomes):
            continue
        rust_dir = chr_dir / RUST_SUBDIR
        diff_dir = chr_dir / DIFF_SUBDIR
        for shard_id in iter_shard_ids(rust_dir, diff_dir):
            output_path = shard_output_path(rust_dir, shard_id)
            if output_path is None or output_path.stat().st_mtime < rust_mtime:
                for directory in (rust_dir, diff_dir):
                    for candidate in shard_entries(directory, shard_id):
                        add_target(targets, candidate, base_dir=base_dir)
    return targets


def collect_failed(base_dir: Path, *, labels: set[str] | None, chromosomes: set[str] | None) -> dict[Path, CleanupTarget]:
    targets: dict[Path, CleanupTarget] = {}
    for chr_dir in iter_chr_dirs(base_dir):
        if not matches_scope(chr_dir, labels=labels, chromosomes=chromosomes):
            continue
        java_dir = chr_dir / JAVA_SUBDIR
        rust_dir = chr_dir / RUST_SUBDIR
        diff_dir = chr_dir / DIFF_SUBDIR
        for shard_id in iter_shard_ids(java_dir, rust_dir, diff_dir):
            status_path = diff_dir / f"shard_{shard_id}.status"
            if status_path.read_text(encoding="utf-8").strip() != "FAIL" if status_path.exists() else True:
                continue
            for directory in (java_dir, rust_dir, diff_dir):
                for candidate in shard_entries(directory, shard_id):
                    add_target(targets, candidate, base_dir=base_dir)
    return targets


def collect_verified(base_dir: Path, *, labels: set[str] | None, chromosomes: set[str] | None) -> dict[Path, CleanupTarget]:
    if not RUST_BIN.exists():
        return {}
    rust_mtime = parity_runner.rust_binary_mtime(parity_runner.DEFAULT_RUST_BIN)
    targets: dict[Path, CleanupTarget] = {}
    for chr_dir in iter_chr_dirs(base_dir):
        if not matches_scope(chr_dir, labels=labels, chromosomes=chromosomes):
            continue
        if verified_marker_current(chr_dir / ".verified", rust_mtime):
            for candidate in tree_entries(chr_dir):
                add_target(targets, candidate, base_dir=base_dir)
    return targets


def collect_all(base_dir: Path, *, labels: set[str] | None, chromosomes: set[str] | None) -> dict[Path, CleanupTarget]:
    targets: dict[Path, CleanupTarget] = {}
    if labels is None and chromosomes is None:
        for candidate in tree_entries(base_dir):
            add_target(targets, candidate, base_dir=base_dir)
        return targets
    for chr_dir in iter_chr_dirs(base_dir):
        if not matches_scope(chr_dir, labels=labels, chromosomes=chromosomes):
            continue
        for candidate in tree_entries(chr_dir):
            add_target(targets, candidate, base_dir=base_dir)
    return targets


def collect_targets(scope: str, base_dir: Path, *, labels: set[str] | None, chromosomes: set[str] | None) -> dict[Path, CleanupTarget]:
    if scope == "stale":
        combined = collect_stale_java(base_dir, labels=labels, chromosomes=chromosomes)
        combined.update(collect_stale_rust(base_dir, labels=labels, chromosomes=chromosomes))
        return combined
    if scope == "java":
        return collect_stale_java(base_dir, labels=labels, chromosomes=chromosomes)
    if scope == "rust":
        return collect_stale_rust(base_dir, labels=labels, chromosomes=chromosomes)
    if scope == "failed":
        return collect_failed(base_dir, labels=labels, chromosomes=chromosomes)
    if scope == "verified":
        return collect_verified(base_dir, labels=labels, chromosomes=chromosomes)
    if scope == "all":
        return collect_all(base_dir, labels=labels, chromosomes=chromosomes)
    raise ValueError(f"Unsupported scope: {scope}")


def print_summary(targets: dict[Path, CleanupTarget], *, dry_run: bool) -> None:
    total_bytes = sum(target.size_bytes for target in targets.values())
    total_files = len(targets)
    prefix = "DRY-RUN: Would remove" if dry_run else "DELETED:"
    if dry_run:
        print(f"{prefix} {human_bytes(total_bytes)} across {total_files} files")
    else:
        print(f"{prefix} {human_bytes(total_bytes)} across {total_files} files")

    grouped: dict[str, tuple[int, int]] = {}
    for target in targets.values():
        bytes_total, count_total = grouped.get(target.group, (0, 0))
        grouped[target.group] = (bytes_total + target.size_bytes, count_total + 1)

    for group in sorted(grouped):
        size_bytes, count = grouped[group]
        print(f"  {group}: {human_bytes(size_bytes)} ({count} files)")


def delete_targets(targets: dict[Path, CleanupTarget], *, base_dir: Path) -> None:
    for path in sorted(targets):
        if is_within_base(path, base_dir):
            path.unlink(missing_ok=True)
    empty_dirs = sorted((path for path in base_dir.rglob("*") if path.is_dir()), key=lambda path: len(path.parts), reverse=True)
    for directory in empty_dirs:
        try:
            directory.rmdir()
        except OSError:
            continue


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if not BASE_DIR.exists():
        print("No parity data found")
        return 0

    try:
        labels = resolve_config_labels(args.config_id)
        chromosomes = {normalize_chr(chrom) for chrom in args.chromosomes} if args.chromosomes else None
        targets = collect_targets(args.scope, BASE_DIR, labels=labels, chromosomes=chromosomes)
    except (FileNotFoundError, ValueError) as exc:
        return fail(str(exc))

    if not targets:
        print("Nothing to clean")
        return 0

    print_summary(targets, dry_run=args.dry_run)
    if not args.dry_run:
        delete_targets(targets, base_dir=BASE_DIR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())