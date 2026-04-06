"""Shard comparison helpers for the parity harness."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import gzip
from typing import TextIO


@dataclass(frozen=True, slots=True)
class CompareResult:
    status: str
    reason: str | None
    first_diff_line: int | None
    java_line: str | None
    rust_line: str | None
    java_line_count: int = 0
    rust_line_count: int = 0


def compare_shard(java_path: Path, rust_path: Path) -> CompareResult:
    if not java_path.exists() or not rust_path.exists():
        return CompareResult("FAIL", "missing_output", None, None, None)
    with java_path.open("r", encoding="utf-8", errors="replace") as java_handle, rust_path.open(
        "r", encoding="utf-8", errors="replace"
    ) as rust_handle:
        return _compare_streams(java_path, rust_path, java_handle, rust_handle)


def compare_shard_gz(java_path: Path, rust_path: Path) -> CompareResult:
    if not java_path.exists() or not rust_path.exists():
        return CompareResult("FAIL", "missing_output", None, None, None)
    with gzip.open(java_path, "rt", encoding="utf-8", errors="replace") as java_handle, gzip.open(
        rust_path, "rt", encoding="utf-8", errors="replace"
    ) as rust_handle:
        return _compare_streams(java_path, rust_path, java_handle, rust_handle)


def _compare_streams(
    java_path: Path,
    rust_path: Path,
    java_handle: TextIO,
    rust_handle: TextIO,
) -> CompareResult:
    java_line_count = 0
    rust_line_count = 0
    first_diff_line: int | None = None
    first_java_line: str | None = None
    first_rust_line: str | None = None

    while True:
        java_line = java_handle.readline()
        rust_line = rust_handle.readline()

        if java_line:
            java_line_count += 1
        if rust_line:
            rust_line_count += 1

        if not java_line and not rust_line:
            break

        if java_line != rust_line and first_diff_line is None:
            first_diff_line = max(java_line_count, rust_line_count)
            first_java_line = _clean_line(java_line)
            first_rust_line = _clean_line(rust_line)

    if java_line_count == 0 and rust_line_count == 0:
        return CompareResult("EMPTY", None, None, None, None)
    if java_line_count == 0:
        return CompareResult(
            "FAIL",
            "java_empty_suspect",
            1 if rust_line_count else None,
            None,
            first_rust_line,
            java_line_count,
            rust_line_count,
        )
    if rust_line_count == 0:
        return CompareResult(
            "FAIL",
            "rust_empty_suspect",
            1 if java_line_count else None,
            first_java_line,
            None,
            java_line_count,
            rust_line_count,
        )
    if first_diff_line is None:
        return CompareResult("PASS", None, None, None, None, java_line_count, rust_line_count)
    if java_line_count != rust_line_count:
        return CompareResult(
            "FAIL",
            "line_count_mismatch",
            first_diff_line,
            first_java_line,
            first_rust_line,
            java_line_count,
            rust_line_count,
        )
    return CompareResult(
        "FAIL",
        "content_mismatch",
        first_diff_line,
        first_java_line,
        first_rust_line,
        java_line_count,
        rust_line_count,
    )


def _clean_line(value: str) -> str | None:
    if not value:
        return None
    return value.rstrip("\n")
