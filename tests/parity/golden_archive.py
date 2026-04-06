#!/usr/bin/env python3
"""Archive and restore Java parity cache directories."""

from __future__ import annotations

import argparse
from datetime import datetime
import hashlib
from pathlib import Path
import sys
import tarfile

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tests.parity import parity_runner
from tests.parity.lib import config


PROJECT_ROOT = Path(__file__).resolve().parents[2]
BASE_DIR = PROJECT_ROOT / parity_runner.BASE_DIR


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Create or restore gzip-compressed Java parity archives.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    create_parser = subparsers.add_parser("create", help="Create an archive from Java cache directories")
    create_parser.add_argument("--output", type=Path, help="Archive path (default: tmp/golden_java_YYYYMMDD.tar.gz)")
    create_parser.add_argument("--config-id", action="append", help="Limit to specific config IDs; repeatable")
    create_parser.add_argument("--chr", dest="chromosomes", action="append", help="Limit to specific chromosomes; repeatable")
    create_parser.add_argument("--dry-run", action="store_true", help="List archive members without writing the archive")

    restore_parser = subparsers.add_parser("restore", help="Restore an archive into the Java cache")
    restore_parser.add_argument("archive", type=Path, help="Archive to restore")
    restore_parser.add_argument("--chr", dest="chromosomes", action="append", help="Limit restore to specific chromosomes; repeatable")
    restore_parser.add_argument("--config-id", action="append", help="Limit restore to specific config IDs; repeatable")
    restore_parser.add_argument("--dry-run", action="store_true", help="List archive members without extracting")

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


def matches_scope(path: Path, *, labels: set[str] | None, chromosomes: set[str] | None) -> bool:
    label = path.parts[0]
    chrom = normalize_chr(path.parts[1])
    if labels is not None and label not in labels:
        return False
    if chromosomes is not None and chrom not in chromosomes:
        return False
    return True


def default_output_path() -> Path:
    stamp = datetime.now().strftime("%Y%m%d")
    return PROJECT_ROOT / "tmp" / f"golden_java_{stamp}.tar.gz"


def archive_members(*, labels: set[str] | None, chromosomes: set[str] | None) -> list[Path]:
    if not BASE_DIR.exists():
        return []
    members: list[Path] = []
    for java_dir in sorted(BASE_DIR.glob("*/*/java")):
        relative = java_dir.relative_to(BASE_DIR)
        if matches_scope(relative, labels=labels, chromosomes=chromosomes):
            members.append(java_dir)
    return members


def file_count(paths: list[Path]) -> int:
    total = 0
    for directory in paths:
        total += sum(1 for path in directory.rglob("*") if path.is_file())
    return total


def total_size(paths: list[Path]) -> int:
    size_bytes = 0
    for directory in paths:
        for path in directory.rglob("*"):
            if path.is_file():
                size_bytes += path.stat().st_size
    return size_bytes


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


def checksum_path(archive_path: Path) -> Path:
    return archive_path.with_name(f"{archive_path.name}.sha256")


def write_checksum(archive_path: Path) -> None:
    digest = hashlib.sha256()
    with archive_path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    checksum_path(archive_path).write_text(f"{digest.hexdigest()}  {archive_path.name}\n", encoding="utf-8")


def verify_checksum(archive_path: Path) -> None:
    checksum_file = checksum_path(archive_path)
    if not checksum_file.exists():
        print(f"WARNING: No checksum file found ({checksum_file}), skipping verification")
        return
    expected_line = checksum_file.read_text(encoding="utf-8").strip()
    expected = expected_line.split()[0] if expected_line else ""
    digest = hashlib.sha256()
    with archive_path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    if digest.hexdigest() != expected:
        raise ValueError("Checksum verification failed")
    print("Checksum OK")


def safe_members(handle: tarfile.TarFile, *, labels: set[str] | None, chromosomes: set[str] | None) -> list[tarfile.TarInfo]:
    selected: list[tarfile.TarInfo] = []
    for member in handle.getmembers():
        member_path = Path(member.name)
        if member_path.is_absolute() or ".." in member_path.parts:
            raise ValueError(f"Unsafe archive member: {member.name}")
        if member.issym() or member.islnk():
            raise ValueError(f"Refusing symlink/hardlink in archive: {member.name}")
        try:
            relative = member_path.relative_to(parity_runner.BASE_DIR)
        except ValueError:
            raise ValueError(f"Archive member is outside parity cache: {member.name}") from None
        if member.isdir():
            continue
        if member_path.parts[-2:] and member_path.parts[-1] == "java":
            continue
        if matches_scope(relative, labels=labels, chromosomes=chromosomes):
            selected.append(member)
    return selected


def create_archive(args: argparse.Namespace) -> int:
    try:
        labels = resolve_config_labels(args.config_id)
    except ValueError as exc:
        return fail(str(exc))
    chromosomes = {normalize_chr(chrom) for chrom in args.chromosomes} if args.chromosomes else None
    members = archive_members(labels=labels, chromosomes=chromosomes)
    if not members:
        return fail(f"No Java cache directories found in {BASE_DIR}")

    archive_path = args.output or default_output_path()
    archive_path = archive_path if archive_path.is_absolute() else PROJECT_ROOT / archive_path

    print(f"Found {len(members)} Java cache directories")
    print(f"Total files: {file_count(members)}")
    print(f"Uncompressed size: {human_bytes(total_size(members))}")

    if args.dry_run:
        print(f"DRY-RUN: Would create archive {archive_path}")
        for directory in members:
            print(f"  {directory.relative_to(PROJECT_ROOT).as_posix()}")
        return 0

    archive_path.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive_path, mode="w:gz") as handle:
        for directory in members:
            handle.add(directory, arcname=directory.relative_to(PROJECT_ROOT))
    write_checksum(archive_path)
    print("Archive created successfully:")
    print(f"  Path: {archive_path.relative_to(PROJECT_ROOT).as_posix()}")
    print(f"  Size: {human_bytes(archive_path.stat().st_size)}")
    print(f"  Checksum: {checksum_path(archive_path).relative_to(PROJECT_ROOT).as_posix()}")
    return 0


def restore_archive(args: argparse.Namespace) -> int:
    archive_path = args.archive if args.archive.is_absolute() else PROJECT_ROOT / args.archive
    if not archive_path.exists():
        return fail(f"Archive not found: {archive_path}")
    try:
        labels = resolve_config_labels(args.config_id)
        verify_checksum(archive_path)
    except ValueError as exc:
        return fail(str(exc))
    chromosomes = {normalize_chr(chrom) for chrom in args.chromosomes} if args.chromosomes else None

    with tarfile.open(archive_path, mode="r:gz") as handle:
        try:
            members = safe_members(handle, labels=labels, chromosomes=chromosomes)
        except ValueError as exc:
            return fail(str(exc))
        if not members:
            return fail("No archive members matched the requested filters")
        if args.dry_run:
            print(f"DRY-RUN: Would restore {len(members)} files from {archive_path.relative_to(PROJECT_ROOT).as_posix()}")
            for member in members[:50]:
                print(f"  {member.name}")
            if len(members) > 50:
                print(f"  ... {len(members) - 50} more")
            return 0
        handle.extractall(PROJECT_ROOT, members=members)

    restored_dirs = sorted(BASE_DIR.glob("*/*/java"))
    restored_files = sum(1 for directory in restored_dirs for path in directory.rglob("*") if path.is_file())
    print("Archive restored successfully:")
    print(f"  Directories: {len(restored_dirs)}")
    print(f"  Files: {restored_files}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "create":
        return create_archive(args)
    if args.command == "restore":
        return restore_archive(args)
    return fail(f"Unknown command: {args.command}")


if __name__ == "__main__":
    raise SystemExit(main())